#!/usr/bin/env python

"""
Helper classes for step7 (write_outputs.py) for filtering loci for the .loci output.
"""


import numpy as np
from numba import njit



class ChunkProcessor(object):
    """
    Takes a chunk of aligned loci and (1) applies filters to it; 
    (2) gets edges, (3) builds snpstring, (4) returns chunk and stats.
    (5) writes 
    """
    def __init__(self, data, chunksize, chunkfile):
        # init data
        self.data = data
        self.chunksize = chunksize
        self.chunkfile = chunkfile
        self.isref = self.data.isref
        self.ispair = self.data.ispair
        self.minsamp = self.data.params.min_samples_locus

        # Minsamp is calculated _before_ the reference sequence is removed
        # and so if we want the minsamp param to be honored as it is written
        # in the params file we need to _add_ 1 to the value, so that when
        # the ref is excluded the minsamp value will be accurate.
        # If the ref is _included_ then it counts toward minsample and no
        # adjustment is necessary.
        if self.isref:
            if self.data.hackersonly.exclude_reference:
                self.minsamp += 1

        # filters (dups, minsamp, maxind, maxall, maxvar, maxshared)
        self.filters = np.zeros((self.chunksize, 5), dtype=np.bool_)
        self.filterlabels = (
            'dups', 
            'maxind',  
            'maxvar', 
            'maxshared',
            'minsamp',
        )
        # (R1>, <R1, R2>, <R2)
        self.edges = np.zeros((self.chunksize, 4), dtype=np.uint16)

        # check filter settings
        self.fmaxsnps = self.data.params.max_SNPs_locus
        if isinstance(self.fmaxsnps, tuple):
            self.fmaxsnps = self.fmaxsnps[0]
        if isinstance(self.fmaxsnps, int):
            self.fmaxsnps = 0.10  # backwards compatibility make as a float

        self.fmaxhet = self.data.params.max_shared_Hs_locus
        if isinstance(self.fmaxhet, tuple):
            self.fmaxhet = self.fmaxhet[0]
        # TODO: This backwards compatibility is hard coded. Maybe better to 
        # just raise an error here, or really during parsing of the params
        # file is best.
        if isinstance(self.fmaxhet, int):
            self.fmaxhet = 0.5  # backwards compatibility make as a float

        self.maxinds = self.data.params.max_Indels_locus
        if isinstance(self.maxinds, tuple):
            self.maxinds = self.maxinds[0]  # backwards compatibility

        # store stats on sample coverage and locus coverage
        self.scov = {i: 0 for i in self.data.snames}
        self.lcov = {i: 0 for i in range(1, len(self.data.snames) + 1)}
        self.var = {i: 0 for i in range(5000)}
        self.pis = {i: 0 for i in range(5000)}
        self.nbases = 0

        # tmp outfile list and filename
        self.outlist = []
        self.outfile = self.chunkfile + '.loci'
        self.outpickle = self.chunkfile + '.p'
        self.outarr = self.chunkfile + '.npy'

        # open a generator to the chunks
        self.chunkio = open(self.chunkfile, 'rb')
        self.loci = enumerate(iter(self.chunkio.read().split(b"//\n//\n")))

        # filled in each chunk
        self.names = []
        self.nidxs = []
        self.aseqs = []
        self.useqs = []


    def next_locus(self):
        self.names = []
        self.nidxs = []
        self.aseqs = []
        self.useqs = []

        # advance locus to next, parse names and seqs
        self.iloc, lines = next(self.loci)
        lines = lines.decode().strip().split("\n")
        for line in lines:
            if line[0] == ">":
                name, nidx = line[1:].rsplit("_", 1)
                self.names.append(name)
                self.nidxs.append(nidx)
            else:
                self.aseqs.append(list(bytes(line.encode())))
                self.useqs.append(list(bytes(line.upper().encode())))

        # filter to include only samples in this assembly
        mask = np.array([i in self.data.snames for i in self.names])
        self.names = np.array(self.names)[mask].tolist()

        if not self.filter_dups():
            # [ref] store consens read start position as mapped to ref
            self.nidxs = np.array(self.nidxs)[mask].tolist()
            self.useqs = np.array(self.useqs)[mask, :].astype(np.uint8)
            self.aseqs = np.array(self.aseqs)[mask, :].astype(np.uint8)



    def run(self):

        # iterate through loci in the chunk
        while 1:
            try:
                self.next_locus()
            except StopIteration:
                break

            # fill filter 0
            if self.filter_dups():
                continue

            # apply filters 
            edges = Edges(self.data, self.useqs)
            edges.get_edges()
            self.edges[self.iloc] = edges.edges

            # fill filter 4
            self.filter_minsamp_pops()
            self.filters[self.iloc, 4] += int(edges.bad)

            # trim edges, need to use uppered seqs for maxvar & maxshared
            edg = self.edges[self.iloc]
            ublock = self.useqs[:, edg[0]:edg[3]]
            ablock = self.aseqs[:, edg[0]:edg[3]]

            # filter if are any empty samples after trimming
            self.filters[self.iloc, 4] += np.sum(np.all(ublock == 45, axis=1))

            # bail out of locus now if it is already bad...
            if self.filters[self.iloc].sum():
                continue

            # [denovo]: store shift of left edge start position from 
            # alignment, this position is needed for pulling depths in VCF.
            # [ref]: nidx string will be updated in to_locus() with edg
            self.masked = None
            if not self.isref:

                # what is the leftmost consens edge (not -)
                ishift = [
                    np.where(self.aseqs[i] != 45)[0].min() 
                    for i in range(self.aseqs.shape[0])
                ]

                # fill nidxs with nidxs and shift info
                inidxs = []
                for idx, (i, j) in enumerate(zip(self.nidxs, ishift)):

                    # add to ishift if trimmed region contains indels
                    indshift = (self.aseqs[idx, j:edges.edges[0]] == 45).size
                    inidxs.append("{}-{}".format(i, j + indshift))
                self.nidxs = inidxs

                # mask insert in denovo data
                self.aseqs[:, edges.edges[1]:edges.edges[2]] = 110  # n
                self.useqs[:, edges.edges[1]:edges.edges[2]] = 78   # N

            # for is-ref we need to mask the insert between pairs
            else:
                if self.ispair and self.data.params.min_samples_locus > 1:
                    inserts = np.all(ublock[1:, :] == 78, axis=0)
                    self.masked = ublock[:, np.invert(inserts)]

            # apply filters on edge trimmed reads
            self.filter_maxindels(ublock)

            # get snpstring on trimmed reads
            if self.isref and self.data.hackersonly.exclude_reference:
                snparr = self.get_snpsarrs(ublock, True)
            else:
                snparr = self.get_snpsarrs(ublock)                
            self.filter_maxvars(ublock, snparr)

            # apply filters on edge trimmed reads
            self.filter_maxshared(ublock)

            # store stats for the locus that passed filtering
            if not self.filters[self.iloc, :].sum():
                # do sample and locus counters
                for name in self.names:
                    self.scov[name] += 1

                # advance locus counter
                if self.isref and self.data.hackersonly.exclude_reference:
                    self.lcov[self.useqs.shape[0] - 1] += 1
                else:
                    self.lcov[self.useqs.shape[0]] += 1

                # do SNP distribution counter
                if self.masked is None:
                    self.nbases += ublock.shape[1]
                else:
                    self.nbases += self.masked.shape[1]
                self.var[snparr[:, :].sum()] += 1
                self.pis[snparr[:, 1].sum()] += 1                   

                # write to .loci string
                locus = self.to_locus(ablock, snparr, edg)
                self.outlist.append(locus)

        # If no loci survive filtering then don't write the files
        if np.fromiter(self.lcov.values(), dtype=int).sum() > 0:
            # write the chunk to tmpdir
            with open(self.outfile, 'w') as outchunk:
                outchunk.write("\n".join(self.outlist) + "\n")

            # thin edgelist to filtered loci and write to array
            mask = np.invert(self.filters.sum(axis=1).astype(np.bool_))
            np.save(self.outarr, self.edges[mask, 0])

        # close file handle
        self.chunkio.close()


    def to_locus(self, block, snparr, edg):
        """
        write chunk to a loci string
        """
        # store as a list 
        locus = []

        # convert snparrs to snpstrings
        snpstring = "".join([
            "-" if snparr[i, 0] else "*" if snparr[i, 1] else " " 
            for i in range(len(snparr))
            ])

        # get nidx string for getting vcf depths to match SNPs
        if self.isref:
            # get ref position from nidxs 
            refpos = ":".join(self.nidxs[0].rsplit(":", 2)[-2:])

            # trim ref position string for edge trims
            chrom, pos = refpos.split(":")
            ostart, end = pos.split("-")
            start = int(ostart) + edg[0]
            end = start + (edg[3] - edg[0])

            # get consens hit indexes and start positions
            nidbits = []
            for bit in self.nidxs[1:]:
                # handle multiple consens merged
                bkey = []
                for cbit in bit.split(";"):
                    cidx, _, pos = cbit.split(":")

                    # start pos of sample is its consens start pos + ostart
                    # where ostart is the ref position start after trim. So 
                    # how far ahead of ref start does the consens read start.
                    posplus = int(pos.split("-")[0]) - int(ostart)
                    bkey.append("{}:{}".format(cidx, posplus))
                nidbits.append("-".join(bkey))                       

            # put ref back into string and append consens hits
            refpos = "{}:{}-{}".format(chrom, start, end)
            nidbits = [refpos] + nidbits
            nidxstring = ",".join(nidbits)

        # denovo stores start read start position in the nidx string
        else:
            nidxstring = ",".join(self.nidxs)

        # if not paired data (with an insert)
        for idx, name in enumerate(self.names):
            locus.append(
                "{}{}".format(
                    self.data.pnames[name],
                    block[idx, :].tostring().decode())
            )
        locus.append("{}{}|{}|".format(
            self.data.snppad, snpstring, nidxstring))
        return "\n".join(locus)


    def filter_dups(self):
        """
        Filter a locus because it contains duplicate names. This happens
        in denovo data if consens reads are slighly more similar to each 
        other than the raw data were.
        """
        if len(set(self.names)) < len(self.names):
            self.filters[self.iloc, 0] = 1
            return True
        return False


    def filter_minsamp_pops(self):
        """
        filter by minsamp or by minsamp x populations.
        """
        # default: no population information
        if not self.data.populations:
            if len(self.names) < self.minsamp:  # data.params.min_samples_locus:
                # store locus filter
                self.filters[self.iloc, 4] = 1
                # return True
            # return False

        # use populations 
        else:
            minfilters = []
            for pop in self.data.populations:
                samps = self.data.populations[pop][1]
                minsamp = self.data.populations[pop][0]
                if len(set(samps).intersection(set(self.names))) < minsamp:
                    minfilters.append(pop)
            if any(minfilters):
                self.filters[self.iloc, 4] = 1
                # return True
            # return False


    def filter_maxindels(self, ublock):
        """
        Max size of internal indels. Denovo vs. Ref, single versus paired.
        """
        # get max indels for read1, read2
        inds = maxind_numba(ublock)        
        if inds > self.maxinds:
            self.filters[self.iloc, 1] = 1
            # return True
        # return False


    def filter_maxvars(self, ublock, snpstring):
        """
        Filter on the max allowed variants 
        """
        # mask insert area
        if self.masked is not None:
            if snpstring.sum() > (self.masked.shape[1] * self.fmaxsnps):
                self.filters[self.iloc, 2] = 1
                # return True

        # use full locus
        else:
            if snpstring.sum() > (ublock.shape[1] * self.fmaxsnps):
                self.filters[self.iloc, 2] = 1
                # return True
        # return False


    def filter_maxshared(self, ublock):
        """
        Paralog detection based on shared heterozygosity across samples.
        """
        nhs = count_maxhet_numba(ublock)
        if nhs > (self.fmaxhet * ublock.shape[0]):
            self.filters[self.iloc, 3] = 1
            # return True
        # return False


    def get_snpsarrs(self, block, exclude_ref=False):
        """
        Count nsnps with option to exclude reference sample from count
        """
        snpsarr = np.zeros((block.shape[1], 2), dtype=np.bool_)
        return snpcount_numba(block, snpsarr, int(bool(exclude_ref)))







class Edges:
    """
    Trims edges of overhanging sequences, cutsites, and pair inserts
    """
    def __init__(self, data, seqs):
        self.data = data
        self.seqs = seqs

        # params
        self.bad = False
        self.exclude_ref = self.data.hackersonly.exclude_reference
        self.edges = np.array([0, 0, 0, self.seqs.shape[1]])
        self.trims = np.array([0, 0, 0, 0])  # self.seqs.shape[1]])
        self.minlen = self.data.params.filter_min_trim_len

        # to be filled
        self.trimseq = None


    def get_edges(self):
        """

        """
        # -1 to site coverage if ref is excluded from the count
        minsites_left = self.data.hackersonly.trim_loci_min_sites
        minsites_right = self.data.hackersonly.trim_loci_min_sites
        if "reference" in self.data.params.assembly_method:
            if self.exclude_ref:
                minsites_left -= 1
                minsites_right -= 1

        # get .edges of good locus or .bad
        self.trim_for_coverage(
            minsite_left=minsites_left,
            minsite_right=minsites_right,
        )

        # fill trimseq with the trimmed sequence array
        self.trimseq = self.seqs[:, self.edges[0]:self.edges[3]]

        # apply edge filtering to locus
        try:
            if not self.bad:
                self.trim_overhangs()
                self.trim_param_trim_loci()
        except Exception:  # TypeError
            self.bad = True
            # TODO: logger here for errors

        # check that locus has anything left
        self.trim_check()


    def trim_for_coverage(self, minsite_left=4, minsite_right=4):
        "trim edges to where data is not N or -"

        # what is the limit of site coverage for trimming?
        minsamp_left = min(minsite_left, self.seqs.shape[0])
        minsamp_right = min(minsite_right, self.seqs.shape[0])        

        # how much cov is there at each site?
        mincovs = np.sum((self.seqs != 78) & (self.seqs != 45), axis=0)

        # locus left trim
        self.edges[0] = locus_left_trim(self.seqs, minsamp_left, mincovs)
        self.edges[3] = locus_right_trim(self.seqs, minsamp_right, mincovs)
        if self.edges[3] <= self.edges[0]:
            self.bad = True

        # find insert region for paired data to mask it...
        self.edges[1] = 0
        self.edges[2] = 0


    def trim_overhangs(self):
        "fuzzy match to trim the restriction_overhangs from r1 and r2"

        # trim left side for overhang
        for cutter in self.data.params.restriction_overhang:

            # skip if None
            if not cutter:
                continue

            # will be ints for py2/3
            cutter = np.array(list(bytes(cutter.encode())))

            # compare match over cut size skipping Ns and allow .25 diffs
            slx = slice(0, cutter.shape[0])
            matching = self.trimseq[:, slx] == cutter
            mask = np.where(
                (self.trimseq[:, slx] != 78) & (self.trimseq[:, slx] != 45))
            matchset = matching[mask]
            if float(matchset.sum()) / matchset.size >= 0.75:
                self.trims[0] = len(cutter)

            # trim right side for overhang
            if self.data.params.restriction_overhang[1]:
                # revcomp the cutter (string not array)
                # cutter = np.array(list(bcomp(cutter.encode())[::-1]))
                slx = slice(
                    self.trimseq.shape[1] - cutter.shape[0], self.trimseq.shape[1])
                matching = self.trimseq[:, slx] == cutter
                mask = np.where(
                    (self.trimseq[:, slx] != 78) & (self.trimseq[:, slx] != 45))
                matchset = matching[mask]
                if float(matchset.sum()) / matchset.size >= 0.75:
                    self.trims[3] = len(cutter)


    def trim_param_trim_loci(self):
        "user entered hard trims"
        self.trims[0] = max([self.trims[0], self.data.params.trim_loci[0]])
        self.trims[1] = (self.trims[1] - self.data.params.trim_loci[1]
            if self.trims[1] else 0)
        self.trims[2] = (self.trims[2] + self.data.params.trim_loci[2]
            if self.trims[2] else 0)
        self.trims[3] = max([self.trims[3], self.data.params.trim_loci[3]])


    def trim_check(self):
        self.edges[0] += self.trims[0]
        self.edges[1] -= self.trims[1]
        self.edges[2] += self.trims[2]
        self.edges[3] -= self.trims[3]

        # checks
        if any(self.edges < 0):
            self.bad = True
        if self.edges[3] <= self.edges[0]:
            self.bad = True
        if self.edges[1] > self.edges[2]:
            self.bad = True
        # check total length including insert
        if (self.edges[3] - self.edges[0]) < self.minlen:
            self.bad = True




# -------------------------------------------------------------
# jitted Edges functions (njit = nopython mode)
# -------------------------------------------------------------

@njit
def locus_left_trim(seqs, minsamp, mincovs):
    leftmost = np.where(mincovs >= minsamp)[0]
    if leftmost.size:
        return leftmost.min()
    return 0

@njit
def locus_right_trim(seqs, minsamp, mincovs):
    rightmost = np.where(mincovs >= minsamp)[0]
    if rightmost.size:
        return rightmost.max() + 1
    return 0




# -------------------------------------------------------------
# jitted ChunkProcessor functions (njit = nopython mode)
# -------------------------------------------------------------

@njit
def maxind_numba(block):
    "count the max size of internal indels"
    inds = 0
    for row in range(block.shape[0]):
        where = np.where(block[row] != 45)[0]
        left = np.min(where)
        right = np.max(where)
        obs = np.sum(block[row, left:right] == 45)
        if obs > inds:
            inds = obs
    return inds


@njit
def snpcount_numba(block, snpsarr, rowstart):
    "Used to count the number of unique bases in a site for snpstring."  

    # iterate over all loci
    for site in range(block.shape[1]):

        # make new array
        catg = np.zeros(4, dtype=np.int16)

        # a list for only catgs
        ncol = block[rowstart:, site]
        for idx in range(ncol.shape[0]):
            if ncol[idx] == 67:    # C
                catg[0] += 1
            elif ncol[idx] == 65:  # A
                catg[1] += 1
            elif ncol[idx] == 84:  # T
                catg[2] += 1
            elif ncol[idx] == 71:  # G
                catg[3] += 1
            elif ncol[idx] == 82:  # R
                catg[1] += 1       # A
                catg[3] += 1       # G
            elif ncol[idx] == 75:  # K
                catg[2] += 1       # T
                catg[3] += 1       # G
            elif ncol[idx] == 83:  # S
                catg[0] += 1       # C
                catg[3] += 1       # G
            elif ncol[idx] == 89:  # Y
                catg[0] += 1       # C
                catg[2] += 1       # T
            elif ncol[idx] == 87:  # W
                catg[1] += 1       # A
                catg[2] += 1       # T
            elif ncol[idx] == 77:  # M
                catg[0] += 1       # C
                catg[1] += 1       # A

        # get second most common site
        catg.sort()

        # if invariant e.g., [0, 0, 0, 9], then nothing (" ")
        if not catg[2]:
            pass
        # store that site is variant as synapomorphy or autapomorphy
        else:           
            if catg[2] > 1:
                snpsarr[site, 1] = True
            else:
                snpsarr[site, 0] = True
    return snpsarr


@njit
def count_maxhet_numba(block):
    """
    Fast jit'd counter of the max number of hetero sites in a locus
    """
    counts = np.zeros(block.shape[1], dtype=np.int16)
    for fidx in range(block.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(block[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts.max()



AMBIGARR = np.array(list(bytes(b"RSKYWM"))).astype(np.uint8)