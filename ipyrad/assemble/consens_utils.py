#!/usr/bin/env python

"""
Utilities for s5_consensus
"""

import os
import sys
import glob
import gzip
import subprocess as sps
from collections import Counter
import h5py
import numpy as np
import pandas as pd
import scipy.stats
from ipyrad.assemble.utils import (
    clustdealer, IPyradError, DCONS, TRANS, CIGARDICT
)


BIN_SAMTOOLS = os.path.join(sys.prefix, "bin", "samtools")


class Processor:
    """
    The consensus calling process that calls alleles and genotypes and
    applies filters to consens seqs. It writes results to two files
    for each chunk, tmpcons (SAM) and tmpcats (HDF5). 
    """
    def __init__(self, data, sample, chunkfile, is_ref):

        # input params
        self.data = data
        self.sample = sample
        self.chunkfile = chunkfile
        self.is_ref = is_ref

        # set max limits and params from data
        self.tmpnum = int(self.chunkfile.split(".")[-1])
        self.optim = int(self.chunkfile.split(".")[-2])

        self.maxlen = self.data.max_frag
        self.maxhet = self.data.params.max_h_consens
        self.maxalleles = self.data.params.max_alleles_consens
        self.maxn = (
            self.data.params.max_n_consens if not self.is_ref else int(1e6))

        # target attrs to be filled
        self.counters = {}
        self.filters = {}

        # tmp attrs to be filled on each cluster
        self.names = None
        self.seqs = None
        self.ref_position = None
        self.consens = None
        self.hidx = None
        self.nheteros = None
        self.nalleles = 1

        # local copies to use to fill the arrays
        self.catarr = np.zeros((self.optim, self.data.max_frag, 4), dtype=np.uint32)
        self.nallel = np.zeros((self.optim, ), dtype=np.uint8)
        self.refarr = np.zeros((self.optim, 3), dtype=np.int64)

        # store data for each chunk with a stats counters.
        self.counters = {
            "name": self.tmpnum,
            "heteros": 0,
            "nsites": 0,
            "nconsens": 0,
        }

        # store filters applied to each chunk
        self.filters = {
            "depth": 0,
            "maxh": 0,
            "maxn": 0,
            "maxalleles": 0,
        }

        # store data for writing
        self.storeseq = {}

        # if reference-mapped then parse the fai (TSV) to get all scaffs.
        # Names and order of scaffs is not predictable so we create a 
        # dict to map {int: scaffname} and {scaffname: int}
        if self.is_ref:
            fai = pd.read_csv(
                self.data.params.reference_sequence + ".fai",
                names=['scaffold', 'size', 'sumsize', 'a', 'b'],
                sep="\t",
                dtype=object)
            self.faidict = {j: i for i, j in enumerate(fai.scaffold)}
            self.revdict = {j: i for i, j in self.faidict.items()}


    def run(self):
        """
        The main function to run the processor
        """
        self.process_chunk()
        self.write_chunk()


    def process_chunk(self):
        """
        iterate over clusters processing each sequentially.
        """
        # load 2 lines at a time to get <name> <sequence>
        inclust = open(self.chunkfile, 'rt')
        pairdealer = zip(*[inclust] * 2)
        done = 0

        # stream through the clusters
        while not done:

            # repeat pairdealer until an entire cluster is pulled in.
            # returns 1 if it reaches the end of the file.
            done, chunk = clustdealer(pairdealer, 1)

            # if chunk was returned then process it.
            if not chunk:
                continue

            # fills .name and .seqs attributes
            # seqs is uint8 with missing as 78 or 45                
            self.parse_cluster(chunk)

            # denovo only: mask repeats (drops lowcov sites w/ dashes) and
            # converts remaining dashes to Ns
            if not self.is_ref:
                # return 1 if entire seq was not masked
                if self.mask_repeats():
                    continue

            # return 1 if cov > mindepth_mj (or >max), filter[depth] + 1
            if self.filter_mindepth():
                continue

            # fills .consens with genotype calls and trims .seqs
            # consens is uint8 with all missing as Ns (78)
            triallele, mindepth2 = self.new_build_consens()

            # all sites are below mindepth even though read cov was not
            # filter[depth] + 1
            if self.filter_mindepth2(mindepth2):
                continue

            # simple 3rd-allele filter, filter[maxallele] + 1
            # Returns 1 if allele not allowed. See also allele filter.
            if self.filter_triallele(triallele):
                continue

            # fills .hidx and .nheteros
            self.get_heteros()

            # return 1 if too many heterozygote sites, filter[maxh] + 1
            if self.filter_maxhetero():
                continue

            # return 1 if too many N or too short, filter[maxn] + 1
            # maxN set arbitrarily high for REF to allow contig-building.
            if self.filter_max_n_min_len():
                continue

            # Infer haplotypes from reads at variable geno sites
            # filter[maxalleles] + 1
            if self.nheteros:
                if self.filter_alleles():
                    continue

            # store result
            self.store_data()

        # cleanup close handle
        inclust.close()



    def parse_cluster(self, chunk):
        """
        Read in cluster chunk to get .names & .seqs and ref position
        """
        # get names and seqs
        piece = chunk[0].strip().split("\n")
        self.names = piece[0::2]
        seqs = piece[1::2]

        if self.data.hackers.declone_PCR_duplicates:
            reps = [1 for i in self.names]
        else:
            reps = [int(n.split(";")[-2][5:]) for n in self.names]

        # fill array with copies of each seq unless decloning is turned
        # on in the hackersdict, in which case we do not count reps.
        sarr = np.zeros((sum(reps), len(seqs[0])), dtype=np.uint8)
        idx = 0
        for ireps, iseq in zip(reps, seqs):
            for _ in range(ireps):
                sarr[idx] = np.array(list(iseq)).astype(bytes).view(np.uint8)
                idx += 1

        # trim array to allow maximum length
        self.seqs = sarr[:, :self.maxlen]

        # ref positions (chromint, pos_start, pos_end)
        self.ref_position = (-1, 0, 0)
        if self.is_ref:
            # parse position from name string
            rname = self.names[0].rsplit(";", 2)[0]
            chrom, posish = rname.rsplit(":")
            pos0, pos1 = posish.split("-")
            # drop `tag=ACGTACGT` if 3rad and declone_PCR_duplicates is set
            pos1 = pos1.split(";")[0]
            # pull idx from .fai reference dict
            chromint = self.faidict[chrom] + 1
            self.ref_position = (int(chromint), int(pos0), int(pos1))


    def mask_repeats(self):
        """
        Removes mask columns with low depth repeats from denovo clusters.
        """
        # get column counts of dashes
        idepths = np.sum(self.seqs == 45, axis=0).astype(float)

        # get proportion of bases that are dashes at each site
        props = idepths / self.seqs.shape[0] 

        # if proportion of dashes sites is more than 0.8?
        keep = props < 0.9

        # report to logger
        # if np.any(props >= 0.9):
        #     logger.debug(
        #         "repeat masked {}:\n{}"
        #         .format(np.where(~keep)[0], self.seqs[:, ~keep])
        #     )

        # drop sites from seqs
        self.seqs = self.seqs[:, keep]
        self.seqs[self.seqs == 45] = 78

        # apply filter in case ALL sites were dropped
        if self.seqs.size:
            return 0          
        self.filters['depth'] += 1
        return 1


    def filter_mindepth(self):
        """
        return 1 if READ depth < minimum
        """
        sumreps = self.seqs.shape[0]
        bool1 = sumreps >= self.data.params.min_depth_majrule
        bool2 = sumreps <= self.data.params.max_depth
        if bool1 & bool2:       
            return 0

        # return that this cluster was filtered out
        self.filters['depth'] += 1
        return 1


    def new_build_consens(self):
        """
        Replace function below using int arrays. 
        Returns 1 if evidence of a triallele.
        """
        # get unphased genotype call of this sequence
        consens, triallele = new_base_caller(
            self.seqs, 
            self.data.params.min_depth_majrule, 
            self.data.params.min_depth_statistical, 
            self.data.hetero_est,
            self.data.error_est,
        )

        # TODO: for denovo we could consider trimming Ns from internal
        # split near pair separator, applied here.
        # logger.debug(consens)

        # trim Ns (uncalled genos) from the left and right
        trim = np.where(consens != 78)[0]

        # if everything was trimmed then return empty 
        if not trim.any():
            return 0, 1

        # otherwise trim edges
        ltrim = trim.min()
        rtrim = trim.max() + 1
        self.consens = consens[ltrim: rtrim]
        self.seqs = self.seqs[:, ltrim: rtrim]

        # update position for trimming
        self.ref_position = (
            self.ref_position[0], 
            self.ref_position[1] + ltrim,
            self.ref_position[1] + ltrim + rtrim,
        )
        # return triallele filter, mindepth2 filter
        return triallele, 0


    def filter_triallele(self, triallele):
        """
        A simple filter on 3rd alleles if only allowing max 2
        """
        if triallele:
            if self.data.params.max_alleles_consens < 3:
                self.filters['maxalleles'] += 1
                return 1
        return 0


    def filter_mindepth2(self, mindepth2):
        """
        A simple filter on 3rd alleles if only allowing max 2
        """
        if mindepth2:
            self.filters['depth'] += 1
            return 1
        return 0


    def get_heteros(self):
        """
        Record the indices of heterozygous sites which will be used
        to identify the number of alleles at this locus. Easier to 
        redo this now after the seqs have been trimmed.

        # ambiguous sites
               R   K   S   Y   W   M
        array([82, 75, 83, 89, 87, 77], dtype=uint8)
        """
        hsites = np.any([
            self.consens == 82,
            self.consens == 75,
            self.consens == 83,
            self.consens == 89,
            self.consens == 87,
            self.consens == 77,                                                            
        ], axis=0)

        self.hidx = np.where(hsites)[0]
        self.nheteros = self.hidx.size


    def filter_maxhetero(self):
        """
        Return 1 if it PASSED the filter, else 0
        """
        if self.nheteros > (self.consens.size * self.maxhet):
            self.filters['maxh'] += 1
            # print(bytes(self.consens).decode())
            return 1
        return 0


    def filter_max_n_min_len(self):
        """
        Return 1 if it PASSED the filter, else 0. The minLen filter
        here is treated the same as maxN since Ns compose the adjacent
        edges of the locus that were uncalled.
        """
        # trimmed too short in size then count as maxN
        if self.consens.size < self.data.params.filter_min_trim_len:
            self.filters['maxn'] += 1
            return 1
        if (self.consens == 78).sum() > (self.consens.size * self.maxn):
            self.filters['maxn'] += 1
            return 1
        return 0


    def filter_alleles(self):
        """
        Infer the number of alleles from haplotypes.
        
        Easy case:
        AAAAAAAAAATAAAAAAAAAACAAAAAAAA
        AAAAAAAAAACAAAAAAAAAATAAAAAAAA

        T,C
        C,T

        Challenging case: more likely paralogs slip by.
        AAAAAAAAAATAAAAAAAA-----------
        AAAAAAAAAACAAAAAAAA-----------
        -------------AAAAAAAACAAAAAAAA
        -------------AAAAAAAATAAAAAAAA

        T -
        C -
        - C
        - T
        ---------
        """
        # if only one hetero site then there are only two bi-alleles
        if self.nheteros == 1:
            self.nalleles = 2
            return 0

        # array of hetero sites (denovo: this can include Ns)
        harray = self.seqs[:, self.hidx]

        # remove varsites containing N
        harray = harray[~np.any(harray == 78, axis=1)]

        # remove alleles containing a site not called in bi-allele genos
        calls0 = np.array([DCONS[i][0] for i in self.consens[self.hidx]])
        calls1 = np.array([DCONS[i][1] for i in self.consens[self.hidx]])
        bicalls = (harray == calls0) | (harray == calls1)
        harray = harray[np.all(bicalls, axis=1), :]

        # get counts of each allele (e.g., AT:2, CG:2)
        ccx = Counter([tuple(i) for i in harray])

        # PASSED PARALOG FILTERING
        # there are fewer alleles than the limit (default=2)
        if len(ccx) <= self.maxalleles:
            self.nalleles = len(ccx)
            return 0

        # MAYBE CAN PASS ---------------------------------------
        # below here tries to correct alleles in case of errors
        # ------------------------------------------------------ 

        # Try dropping alleles with very low coverage (< 10%)
        alleles = []
        for allele in ccx:
            if (ccx[allele] / self.seqs.shape[0]) >= 0.1:
                alleles.append(allele)  # ccx[allele]
        
        # in case an allele was dropped, now check again.
        # apply filter if all were dropped as lowcov
        ccx = Counter(alleles)
        if len(ccx) <= self.maxalleles:
            self.nalleles = len(ccx)
            return 0

        # DID NOT PASS FILTER ----------------------------------
        self.filters['maxalleles'] += 1            
        self.nalleles = len(ccx)

        # debugging to logger: this cluster was filtered
        # print(f"filtered by max-alleles\n{ccx}")
        return 1


    def store_data(self):
        """
        Store the site count data for writing to HDF5, and stats.
        """
        # current counter
        cidx = self.counters["nconsens"]
        self.nallel[cidx] = self.nalleles
        self.refarr[cidx] = self.ref_position

        # store a reduced array with only CATG
        catg = np.array(
            [np.sum(self.seqs == i, axis=0) for i in (67, 65, 84, 71)],
            dtype=np.uint16,
        ).T

        # base counts: do not allow ints larger than 65535 (uint16)
        # print("catg:\n", catg)
        self.catarr[cidx, :catg.shape[0], :] = catg

        # store the seqdata and advance counters (ints->bytes->str)
        self.storeseq[cidx] = bytes(self.consens).decode()
        # print('storeseq: ', self.storeseq[cidx])
        self.counters["name"] += 1
        self.counters["nconsens"] += 1
        self.counters["heteros"] += self.nheteros


    def write_chunk(self):
        """
        Writes chunk of consens reads to disk, stores depths, alleles, and 
        chroms, and stores stats. For denovo data it writes consens chunk
        as a fasta file. For reference data it writes as a SAM file that is
        compatible to be converted to BAM (very stringent about cigars).
        """
        # the reference arr is bigger than the actual seqarr, so find the
        # size of the actually filled data array (end).
        end = np.where(np.all(self.refarr == 0, axis=1))[0]
        if np.any(end):
            end = end.min()
        else:
            end = self.refarr.shape[0]

        # write final consens string chunk
        consenshandle = os.path.join(
            self.data.tmpdir,
            f"{self.sample.name}_tmpcons.{end}.{self.tmpnum}"
        )

        # write chunk. 
        if self.storeseq:
            outfile = open(consenshandle, 'wt')

            # denovo just write the consens simple
            if not self.is_ref:
                seqlist = []
                for key in self.storeseq:
                    seq = self.storeseq[key]
                    cstring = ">{}_{}\n{}".format(self.sample.name, key, seq)
                    seqlist.append(cstring)
                outfile.write("\n".join(seqlist))

            # reference needs to store if the read is revcomp to reference
            else:
                constring = "\n".join(
                    ["{}:{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                    .format(
                        self.sample.name,
                        self.refarr[i][0],
                        self.refarr[i][1],
                        self.refarr[i][1] + len(self.storeseq[i]),
                        #self.refarr[i][2],
                        0,
                        self.revdict[self.refarr[i][0] - 1],
                        self.refarr[i][1],
                        0,
                        "{}M".format(len(self.storeseq[i])),
                        "*",
                        0,
                        #self.refarr[i][2] - self.refarr[i][1],
                        len(self.storeseq[i]),
                        self.storeseq[i],
                        "*",
                    ) for i in self.storeseq]
                )
                outfile.write(constring)
            outfile.close()

        # store reduced size arrays with indexes matching to keep indexes
        tmp5 = consenshandle.replace("_tmpcons.", "_tmpcats.")
        with h5py.File(tmp5, 'w') as io5:
            # reduce catarr to uint16 size for easier storage
            self.catarr[self.catarr > 65535] = 65535
            self.catarr = self.catarr.astype(np.uint16)
            io5.create_dataset(name="cats", data=self.catarr[:end])
            io5.create_dataset(name="alls", data=self.nallel[:end])
            io5.create_dataset(name="chroms", data=self.refarr[:end])
        del self.catarr
        del self.nallel
        del self.refarr

        # return stats and skip sites that are Ns (78)
        #self.counters['nsites'] = sum([len(i) for i in self.storeseq.values()])
        self.counters['nsites'] = sum(
            sum(1 if j != 78 else 0 for j in i) 
            for i in self.storeseq.values()
        )
        del self.storeseq



def new_base_caller(sarr, mindepth_mj, mindepth_st, est_het, est_err):
    """
    Call site genotypes from site counts. Uses scipy binomial.
    """
    # start with an array of Ns
    cons = np.zeros(sarr.shape[1], dtype=np.uint8)
    cons.fill(78)

    # record if evidence of a tri-allele
    triallele = 0

    # iterate over columns
    for cidx in range(sarr.shape[1]):

        # get col, nmask, dashmask, and bothmask
        col = sarr[:, cidx]
        nmask = col == 78
        dmask = col == 45
        bmask = nmask | dmask

        # if non-masked site depth is below mindepth majrule fill with N (78)
        if np.sum(~bmask) < mindepth_mj:
            cons[cidx] = 78
            continue
                
        # if not variable (this will include the 'nnnn' pair separator (110))
        # in denovo data sets
        masked = col[~bmask]
        if np.all(masked == masked[0]):
            cons[cidx] = masked[0]
            continue

        # make statistical base calls on allele frequencies
        counts = np.bincount(col, minlength=79)
        counts[78] = 0
        counts[45] = 0

        # get allele freqs (first-most, second, third = p, q, r)
        pbase = np.argmax(counts)
        nump = counts[pbase]
        counts[pbase] = 0

        qbase = np.argmax(counts)
        numq = counts[qbase]
        counts[qbase] = 0

        rbase = np.argmax(counts)
        numr = counts[rbase]
        counts[rbase] = 0

        # if third allele occurs >X% then fill N and mark as paralog
        # 1/6 as the cutoff
        if (numr / (nump + numq + numr)) >= 0.15:
            triallele = 1

        # based on biallelic depth
        bidepth = nump + numq
        if bidepth < mindepth_mj:
            cons[cidx] = 78
            continue

        # if depth is too high, reduce to sampled int
        if bidepth > 500:
            nump = int(500 * (nump / float(bidepth)))
            numq = int(500 * (numq / float(bidepth)))

        # make majority-rule call
        if bidepth < mindepth_st:
            if nump == numq:
                cons[cidx] = TRANS[(pbase, qbase)]
            else:
                cons[cidx] = pbase        

        # make statistical base call
        ishet, prob = get_binom(nump, numq, est_err, est_het)
        if prob < 0.95:
            cons[cidx] = 78
        else:
            if ishet:
                cons[cidx] = TRANS[(pbase, qbase)]
                # print(f'cidx {cidx}; cons[cidx] {cons[cidx]}')                
            else:
                cons[cidx] = pbase
    return cons, triallele



def get_binom(base1, base2, est_err, est_het):
    """
    return probability of base call
    """
    prior_homo = (1. - est_het) / 2.
    prior_hete = est_het

    ## calculate probs
    bsum = base1 + base2
    hetprob = scipy.special.comb(bsum, base1) / (2. ** (bsum))
    homoa = scipy.stats.binom.pmf(base2, bsum, est_err)
    homob = scipy.stats.binom.pmf(base1, bsum, est_err)

    ## calculate probs
    hetprob *= prior_hete
    homoa *= prior_homo
    homob *= prior_homo

    ## final
    probabilities = [homoa, homob, hetprob]
    bestprob = max(probabilities) / float(sum(probabilities))

    ## return
    if hetprob > homoa:
        return True, bestprob
    return False, bestprob



def make_cigar(arr):
    """
    Writes a cigar string with locations of indels and lower case ambigs
    """
    # simplify data
    arr[np.char.islower(arr)] = '.'
    indel = np.bool_(arr == "-")
    ambig = np.bool_(arr == ".")
    arr[~(indel + ambig)] = "A"

    # counters
    cigar = ""
    mcount = 0
    tcount = 0
    lastbit = arr[0]
    for _, j in enumerate(arr):

        # write to cigarstring when state change
        if j != lastbit:
            if mcount:
                cigar += "{}{}".format(mcount, "M")
                mcount = 0
            else:
                cigar += "{}{}".format(tcount, CIGARDICT.get(lastbit))
                tcount = 0
            mcount = 0
            
        # increase counters
        if j in ('.', '-'):
            tcount += 1
        else:
            mcount += 1
        lastbit = j
        
    # write final block
    if mcount:
        cigar += "{}{}".format(mcount, "M")
    if tcount:
        cigar += '{}{}'.format(tcount, CIGARDICT.get(lastbit))
    return cigar



def concat_catgs(data, sample, isref):
    """
    Concat catgs into a single sample catg and remove tmp files
    """
    # collect tmpcat files written by write_chunks()
    tmpcats = glob.glob(os.path.join(data.tmpdir, f"{sample.name}_tmpcats.*"))
    tmpcats.sort(key=lambda x: int(x.split(".")[-1]))

    # get full nrows of the new h5 from the tmpcat filenames
    nrows = sum([int(i.rsplit(".", 2)[-2]) for i in tmpcats])

    # get shape info from the first cat, (optim, maxlen, 4)
    optim = min(nrows, 5000)
    with h5py.File(tmpcats[0], 'r') as io5:
        _, maxlen, _ = io5['cats'].shape

    # Check values of nrows and maxlen are > 0
    # This literally shouldn't happen, but it does, or has at least twice.
    # Related to issue #369
    if not all([nrows, maxlen]):
        raise IPyradError(
            "Error in concat_catgs both nrows and maxlen must be positive."
            "\nsample: {}\tnrows: {}\tmaxlen: {}"
            .format(sample.name, nrows, maxlen))

    # fill in the chunk array
    with h5py.File(sample.files.depths, 'w') as ioh5:
        dcat = ioh5.create_dataset(
            name="catg",
            shape=(nrows, maxlen, 4),
            dtype=np.uint32,
            chunks=(optim, maxlen, 4),
            compression="gzip")
        dall = ioh5.create_dataset(
            name="nalleles", 
            shape=(nrows, ),
            dtype=np.uint8,
            chunks=(optim, ),
            compression="gzip")
        
        # only create chrom for reference-aligned data
        if isref:
            dchrom = ioh5.create_dataset(
                name="chroms",
                shape=(nrows, 3),
                dtype=np.int64,
                chunks=(optim, 3),
                compression="gzip")

        # Combine all those tmp cats into the big cat
        start = 0
        for icat in tmpcats:
            addon = int(icat.rsplit(".", 2)[-2])
            end = start + addon
            io5 = h5py.File(icat, 'r')
            dcat[start:end] = io5['cats'][:]
            dall[start:end] = io5['alls'][:]
            if isref:
                dchrom[start:end] = io5['chroms'][:]
            start = end
            io5.close()
            os.remove(icat)



def concat_denovo_consens(data, sample):
    """
    Concatenate consens bits into fasta file for denovo assemblies
    """
    # collect consens chunk files
    combs1 = glob.glob(os.path.join(data.tmpdir, f"{sample.name}_tmpcons.*"))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    # stream through file renaming consens reads and write out
    idx = 0    
    with gzip.open(sample.files.consens, 'wt') as out:
        for fname in combs1:
            cstack = []
            #starter = int(fname.rsplit(".")[-1])
            with open(fname) as infile:
                for line in infile:
                    if line.startswith(">"):
                        name, _ = line.rsplit("_", 1)
                        line = name + "_{}".format(idx)
                        idx += 1
                    cstack.append(line.strip())
                out.write("\n".join(cstack) + "\n")
            os.remove(fname)



def concat_reference_consens(data, sample):
    """
    Concatenates consens bits into SAM for reference assemblies
    """
    # write sequences to SAM file
    chunks = glob.glob(os.path.join(data.tmpdir, f"{sample.name}_tmpcons.*"))
    chunks.sort(key=lambda x: int(x.split(".")[-1]))

    # open sam handle for writing to bam
    bamfile = os.path.join(data.stepdir, f"{sample.name}.bam")
    samfile = os.path.join(data.tmpdir, f"{sample.name}.sam")
    with open(samfile, 'wt') as outf:

        # parse fai file for writing headers
        fai = "{}.fai".format(data.params.reference_sequence)
        fad = pd.read_csv(fai, sep="\t", names=["SN", "LN", "POS", "N1", "N2"])
        headers = ["@HD\tVN:1.0\tSO:coordinate"]
        headers += [
            "@SQ\tSN:{}\tLN:{}".format(i, j)
            for (i, j) in zip(fad["SN"], fad["LN"])
        ]
        outf.write("\n".join(headers) + "\n")

        # write to file with sample names imputed to line up with catg array
        counter = 0
        for fname in chunks:
            with open(fname) as infile:
                # impute catg ordered seqnames 
                data = infile.readlines()
                fdata = []
                for line in data:
                    name, chrom, rest = line.rsplit(":", 2)
                    fdata.append(
                        "{}_{}:{}:{}".format(name, counter, chrom, rest)
                        )
                    counter += 1
                outf.write("".join(fdata) + "\n")

    # convert to bam
    cmd = [BIN_SAMTOOLS, "view", "-Sb", samfile]
    with open(bamfile, 'w') as outbam:
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=outbam)
        comm = proc.communicate()[1]
    if proc.returncode:
        raise IPyradError("error in samtools: {}".format(comm))



# def encode_alleles(consens, hidx, alleles):
#     """ 
#     Store phased allele data for diploids 
#     """
#     ## find the first hetero site and choose the priority base
#     ## example, if W: then priority base in A and not T. PRIORITY=(order: CATG)
#     bigbase = PRIORITY[consens[hidx[0]]]

#     ## find which allele has priority based on bigbase
#     bigallele = [i for i in alleles if i[0] == bigbase][0]

#     ## uplow other bases relative to this one and the priority list
#     ## e.g., if there are two hetero sites (WY) and the two alleles are
#     ## AT and TC, then since bigbase of (W) is A second hetero site should
#     ## be stored as y, since the ordering is swapped in this case; the priority
#     ## base (C versus T) is C, but C goes with the minor base at h site 1.
#     #consens = list(consens)
#     for hsite, pbase in zip(hidx[1:], bigallele[1:]):
#         if PRIORITY[consens[hsite]] != pbase:
#             consens[hsite] = consens[hsite].lower()

#     ## return consens
#     return consens

