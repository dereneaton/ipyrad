#!/usr/bin/env python

# py2/3 compatibility
from __future__ import print_function
try:
    from builtins import range
    from itertools import izip  # , chain
except ImportError:
    # from itertools import chain
    izip = zip

# standard lib imports
import os
import glob
import time
import shutil

# third party imports
import numpy as np
import pandas as pd
from numba import njit
from .utils import IPyradError, clustdealer  #, splitalleles, BTS, AMBIGS, TRANSFULL

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

# globals
AMBIGARR = np.array(list(b"RSKYWM"))

# classes
class Step7:
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.samples = self.get_subsamples()
        self.setup_dirs()
        self.get_chunksize()


    def run(self):
        # split clusters into bits.
        self.split_clusters()

        # get filter and snp info on edge trimmed data.
        # write to chunks for building output files from.
        # write stats file.
        self.remote_process_chunks()

        # iterate over new processed chunks to fill the h5 array.
        # build loci and stats while counting nsnps and nbases.
        # build arrays and fill in second iteration over chunks.
        self.write_loci_and_fill_arrays()


        # 
        # and fill the seq and snp arrays.
        
        self.build_loci()
        self.build_alleles()

        # build conversions 
        self.remote_fill_seqarrs()
        self.remote_build_conversions()

        # send jobs to build vcf
        self.remote_fill_depths()
        self.remote_build_vcf()

    ## init functions -----------------------------------
    def get_subsamples(self):
        "get subsamples for this assembly. All must have been in step6"
        # filter samples by state
        state5 = self.data.stats.index[self.data.stats.state < 6]
        state6 = self.data.stats.index[self.data.stats.state == 6]
        state7 = self.data.stats.index[self.data.stats.state > 6]

        # tell user which samples are not ready for step5
        if state5.any():
            print("skipping samples not in state==6:\n{}"
                  .format(state5.tolist()))

        if self.force:
            # run all samples above state 5
            subs = self.data.stats.index[self.data.stats.state > 5]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 6
            if state7.any():
                raise IPyradError(
                    "some samples are already in state==7. If you wish to \n" \
                  + "create new outfiles for this assembly use the force \n" \
                  + "argument.")
            # run all samples in state 6
            subsamples = [self.data.samples[i] for i in state6]

        # check that kept samples were in the step6 database
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_consens:
                checked_samples.append(sample)
            else:
                print("skipping {}; no consensus reads found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 7")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_consens,
            reverse=True,
        )
        return checked_samples


    def setup_dirs(self):
        "Create temp h5 db for storing filters and depth variants"

        # make new output directory
        self.data.dirs.outfiles = os.path.join(
            self.data.paramsdict["project_dir"],
            "{}_outfiles".format(self.data.name),
            )
        if os.path.exists(self.data.dirs.outfiles):
            shutil.rmtree(self.data.dirs.outfiles)
        if not os.path.exists(self.data.dirs.outfiles):
            os.makedirs(self.data.dirs.outfiles)

        # make tmpdir directory
        self.data.tmpdir = os.path.join(
            self.data.paramsdict["project_dir"],
            "{}_outfiles".format(self.data.name),
            "tmpdir",
            )
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.makedirs(self.data.tmpdir)

        # make new database file
        self.data.database = os.path.join(
            self.data.dirs.outfiles,
            self.data.name + ".hdf5",
            )


    def get_chunksize(self):
        "get nloci and ncpus to chunk and distribute work across processors"
        # count number of loci
        self.rawloci = os.path.join(
            self.data.dirs.across, 
            self.data.name + "_raw_loci.fa")
        with open(self.rawloci, 'r') as inloci:
            self.nraws = sum(1 for i in inloci if i == "//\n") // 2

        # chunk to approximately 2 chunks per core
        self.ncpus = len(self.ipyclient.ids)
        self.chunks = ((self.nraws // (self.ncpus * 2)) + \
                       (self.nraws % (self.ncpus * 2)))


    def write_loci_and_fill_arrays(self):

        # open new database file handle
        with h5py.File(self.data.database, 'w') as io5:

            # attributes
            io5.attrs["samples"] = [i.name.encode() for i in self.samples]
            io5.attrs["filters"] = [
                b"duplicates",
                b"minsamples",
                b"maxindels", 
                b"maxalleles",
                b"maxvars",
                b"maxsharedhets",
                ]

            # arrays
            # io5.create_dataset(
            #     name="edges", 
            #     shape=(self.nraws, 5),
            #     dtype=np.uint16,
            #     chunks=(self.chunks, 5),
            #     compression='gzip')


    def split_clusters(self):
        with open(self.rawloci, 'rb') as clusters:
            pairdealer = izip(*[iter(clusters)] * 2)

            # grab a chunk of clusters
            idx = 0
            while 1:

                # if an engine is available pull off a chunk
                try:
                    done, chunk = clustdealer(pairdealer, self.chunks)
                except IndexError:
                    raise IPyradError("rawloci formatting error in %s", chunk)

                # write to tmpdir and increment counter
                if chunk:
                    chunkpath = os.path.join(
                        self.data.tmpdir, 
                        "chunk-{}".format(idx),
                        )
                    with open(chunkpath, 'wb') as outfile:
                        outfile.write(b"//\n//\n".join(chunk))
                    idx += 1

                # break on final chunk
                if done:
                    break

    def remote_process_chunks(self):
        "Calls process chunks in parallel."
        start = time.time()
        printstr = ("applying filters    ", "s7")
        rasyncs = {}

        jobs = glob.glob(os.path.join(self.data.tmpdir, "chunk-*"))
        for jobfile in jobs:
            args = (self.data, jobfile)
            rasyncs[jobfile] = self.lbview.apply(process_chunks, *args)
        
        # iterate until all chunks are processed
        while 1:
            # get and enter results into hdf5 as they come in
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                break          

        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].exception())



class Processor:
    def __init__(self, data, chunks, chunkfile):
        """
        Takes a chunk of aligned loci and (1) applies filters to it; 
        (2) gets edges, (3) builds snpstring, (4) returns chunk and stats.
        """
        # init data
        self.data = data
        self.chunks = chunks
        self.chunkfile = chunkfile

        # filters (dups, minsamp, maxind, maxall, maxvar, maxshared)
        self.filters = np.zeros((self.chunks, 5), dtype=np.bool_)
        self.filterlabels = (
            'dups', 'minsamp', 'maxind',  # 'maxall', 
            'maxvar', 'maxshared')
        self.edges = np.zeros((self.chunks, 4), dtype=np.uint16)

        # dict mapping of samples to padded names for loci file aligning.
        self.snames = sorted(list(self.data.samples.keys()))
        self.pnames, self.snppad = self.get_padded_names()      

        # store stats on sample coverage and locus coverage
        self.scov = {i: 0 for i in self.snames}
        self.lcov = {i: 0 for i in range(1, len(self.snames) + 1)}

        # tmp outfile list and filename
        self.outlist = []
        self.outfile = self.chunkfile + '.loci'

    def run(self):
        # todo: this could be an iterator...
        with open(self.chunkfile, 'rb') as infile:
            loci = infile.read().split(b"//\n//\n")
    
            # iterate over loci
            for iloc, loc in enumerate(loci):                              
                # load names and seqs 
                lines = loc.decode().strip().split("\n")
                names = []
                nidxs = []
                seqs = []
                for line in lines:
                    if line[0] == ">":
                        name, nidx = line[1:].rsplit("_", 1)
                        names.append(name)
                        nidxs.append(nidx)
                    else:
                        seqs.append(list(line))

                #names = [i[1:].rsplit("_", 1)[0] for i in lines[::2]]
                #name_indices = [i[1:].rsplit("_", 1)[1] for i in lines[::2]]
                #seqs = np.array([list(i) for i in lines[1::2]])
                #seqs = np.array([list(i.upper()) for i in lines[1::2]])
                #seqs = seqs.astype(bytes).view(np.uint8)

                # filter to only include only samples in this assembly
                mask = [i in self.snames for i in names]
                names = np.array(names)[mask].tolist()
                nidxs = np.array(nidxs)[mask].tolist()
                seqs = np.array(seqs)[mask, :].astype(bytes).view(np.uint8)
                
                # apply filters
                efilter, edges = self.get_edges(seqs)
                self.edges[iloc] = edges
                self.filters[iloc, 0] = self.filter_dups(names)
                self.filters[iloc, 1] = self.filter_minsamp_pops(names)
                self.filters[iloc, 1] += efilter

                # should we fill terminal indels as N's here?
                #...

                # trim edges
                edg = self.edges[iloc]
                block1 = seqs[:, edg[0]:edg[1]]
                block2 = seqs[:, edg[2]:edg[3]]

                # apply filters on edge trimmed reads
                self.filters[iloc, 2] += self.filter_maxindels(block1, block2)

                # get snpstring on trimmed reads
                snparr1, snparr2 = self.get_snpstrings(block1, block2)
                self.filters[iloc, 3] = self.filter_maxvars(snparr1, snparr2)

                # apply filters on edge trimmed reads
                self.filters[iloc, 4] = self.filter_maxshared(block1, block2)
                #self.filters[iloc, 5] = self.filter_maxalleles()

                # store stats for the locus that passed filtering
                for name in names:
                    self.scov[name] += 1
                self.lcov[seqs.shape[0]] += 1

                # write to .loci string
                locus = self.to_locus(names, nidxs, block1, block2, snparr1, snparr2)
                self.outlist.append(locus)

        # write the chunk to tmpdir
        with open(self.outfile, 'w') as outchunk:
            outchunk.write("\n".join(self.outlist))

    def to_locus(self, names, nidxs, block1, block2, snparr1, snparr2):
        "write chunk to a loci string"

        # store as a list 
        locus = []

        # convert snparrs to snpstrings
        snpstring1 = "".join([
            "-" if snparr1[i, 0] else \
            "*" if snparr1[i, 1] else \
            " " for i in range(len(snparr1))
        ])
        snpstring2 = "".join([
            "-" if snparr2[i, 0] else \
            "*" if snparr2[i, 1] else \
            " " for i in range(len(snparr2))
        ])
        nidxstring = ",".join(nidxs)

        # if not paired data (with an insert)
        if not block2.size:
            for idx, name in enumerate(names):
                locus.append(
                    "{}{}".format(
                        self.pnames[name],
                        block1[idx, :].tostring().decode())
                )
            locus.append("{}{}|{}".format(self.snppad, snpstring1, nidxstring))
        else:
            raise NotImplementedError("see here")

        return "\n".join(locus)

    ## filters based on names -----
    def filter_dups(self, names):
        unames = [i.rsplit("_", 1)[0] for i in names]
        if len(set(unames)) < len(names):
            return True
        return False

    # TODO: support pops
    def filter_minsamp_pops(self, names):
        if not self.data.populations:
            mins = self.data.paramsdict["min_samples_locus"]
            if len(names) < mins:
                return True
            return False

        else:
            # TODO:...
            raise NotImplementedError(
                "please contact the developers about this issue.")
            for pop in self.data._populations:
                return True

    ## filters based on seqs ------
    def filter_maxindels(self, block1, block2):
        # get max indels for read1, read2
        maxi = self.data.paramsdict["max_Indels_locus"]
        maxi = np.array(maxi).astype(np.int64)

        # single-end
        if "pair" not in self.data.paramsdict["datatype"]:
            inds = maxind_numba(block1)
            if inds > maxi[0]:
                return True
            return False

        # todo: if paired then worry about splits
        else:
            inds1 = maxind_numba(block1)
            inds2 = maxind_numba(block2)            
            if (inds1 > maxi[0]) or (inds2 > maxi[1]):
                return True
            return False

    def filter_maxvars(self, snpstring1, snpstring2):

        # get max indels for read1, read2
        maxs = self.data.paramsdict["max_SNPs_locus"]
        maxs = np.array(maxs).astype(np.int64)

        # single-end
        if "pair" not in self.data.paramsdict["datatype"]:
            if np.sum(snpstring1, axis=1).sum() > maxs[0]:
                return True
            return False

        else:
            if (np.sum(snpstring1, axis=1).sum() > maxs[0]) or \
               (np.sum(snpstring2, axis=1).sum() > maxs[1]):
                return True
            return False

    def filter_maxshared(self, block1, block2):
        blocks = np.concatenate([block1, block2], axis=1)
        # calculate as a proportion
        maxhet = self.data.paramsdict["max_shared_Hs_locus"]       
        if isinstance(maxhet, float):
            maxhet = np.floor(maxhet * blocks.shape[0]).astype(np.int16)
        elif isinstance(maxhet, int):
            maxhet = np.int16(maxhet)

        # get max from combined block
        if maxhet_numba(blocks, maxhet).max() > maxhet:
            return True
        return False

    def get_edges(self, seqs):
        """
        Trim based on three criteria and take the max for each edge.
        1. user entered hard trimming.
        2. removing cutsite overhangs.
        3. trimming singleton-like overhangs from seqs of diff lengths.
        """
        # 1. hard trim edges
        trim1 = np.array(self.data.paramsdict["trim_loci"])

        # 2. fuzzy match for trimming restriction site where it's expected.
        trim2 = np.array([0, 0, 0, 0])
        overhangs = np.array([
            i.encode() for i in self.data.paramsdict["restriction_overhang"]
            ])
        for pos, overhang in enumerate(overhangs):
            if overhang:
                cutter = np.array(list(overhang))
                trim2[pos] = check_cutter(seqs, pos, cutter, 0.75)

        # 3. find where the edge is not indel marked (really unknown ("N"))
        trim3 = np.array([0, 0, 0, 0])
        for pos in range(4):
            trim3[pos] = check_minsamp(seqs, pos, 4)
        
        # get max edges
        trim = np.max([trim1, trim2, trim3], axis=0)

        # return edges as slice indices
        r1left = trim[0]
        r1right = seqs.shape[1] - trim[1]

        # TODO: resolve r2
        r2left = r2right = r1right
        edges = (r1left, r1right, r2left, r2right)
        
        # get filter
        bad = False
        if (r1right < r1left) or (r2left < r1right) or (r2right < r2left):
            bad = True
        return bad, edges

    def get_snpstrings(self, block1, block2):
        if "pair" not in self.data.paramsdict["datatype"]:
            snpstring1 = np.zeros((block1.shape[1], 2), dtype=np.bool_)
            snpstring1 = snpcount_numba(block1, snpstring1)
            snpstring2 = np.array([])
        else:
            snpstring1 = np.zeros((block1.shape[1], 2), dtype=np.bool_)
            snpstring1 = snpcount_numba(block1, snpstring1)
            snpstring2 = np.zeros((block2.shape[1], 2), dtype=np.bool_)
            snpstring2 = snpcount_numba(block2, snpstring2)
        return snpstring1, snpstring2

    def get_padded_names(self):
        # get longest name
        longlen = max(len(i) for i in self.snames)
        # Padding distance between name and seq.
        padding = 5
        # add pad to names
        pnames = {
            name: "{}{}".format(name, " " * (longlen - len(name) + padding))
            for name in self.snames
        }
        snppad = "//" + " " * (longlen - 2 + padding)
        return pnames, snppad


# ----------------------------------------
# Step7 external functions
# ----------------------------------------
def process_chunks(step, chunkfile):
    # process chunk returns writes to files and returns proc with features.
    proc = Processor(step.data, step.chunks, chunkfile)
    proc.run()

    # organize stats into dataframes for statsfile or arrays for h5 database
    ftable = pd.DataFrame(
        columns=["total", "applied_order", "retained_loci"],
        index=[
            "total_prefiltered_loci",
            "filtered_by_rm_duplicates",
            "filtered_by_max_indels",
            "filtered_by_max_snps",
            "filtered_by_max_shared_het",
            "filtered_by_min_sample",  # "filtered_by_max_alleles",
            "total_filtered_loci"],
    )
    ftable.loc[0, :] = step.nraws

    # filter rm dups
    ftable.iloc[1, 0:2] = proc.filters[:, 0].sum()
    ftable.iloc[1, 2] = ftable.iloc[0, 2] - ftable.iloc[1, 1]
    mask = proc.filters[:, 0]

    # filter max indels
    ftable.iloc[2, 0] = proc.filters[:, 1].sum()
    ftable.iloc[2, 1] = proc.filters[mask, 1].sum()
    ftable.iloc[2, 2] = ftable.iloc[2, 0] - ftable.iloc[2, 1]
    mask = proc.filters[:, 0:2].sum(axis=1).astype(np.bool)

    # filter max snps
    ftable.iloc[3, 0] = proc.filters[:, 2].sum()
    ftable.iloc[3, 1] = proc.filters[mask, 2].sum()
    ftable.iloc[3, 2] = ftable.iloc[3, 0] - ftable.iloc[3, 1]


    #proc.filters
    #proc.lcov
    #proc.scov

    # stats dfs
    ltable = pd.Series(proc.lcov)
    stable = pd.Series(proc.scov)
    dataframes = {'ftable': ftable, 'stable': stable, 'ltable': ltable}

    # arrays for outputs
    arrays = {}

    return dataframes



# ----------------------------------------
# Processor external functions
# ----------------------------------------
def check_minsamp(seq, position, minsamp):
    "used in Processor.get_edges() for trimming edges of - or N sites."
    minsamp = min(minsamp, seq.shape[0])
    if position == 0:       
        mincovs = np.sum((seq != 78) & (seq != 45), axis=0)
        return np.where(mincovs >= minsamp)[0].min()
    
    elif position == 1:
        mincovs = np.sum((seq != 78) & (seq != 45), axis=0)
        maxhit = np.where(mincovs >= minsamp)[0].max()
        return seq.shape[1] - (maxhit + 1)

    ## TODO...    
    elif position == 2:
        return 0
    
    else:
        return 0


def check_cutter(seq, position, cutter, cutoff):
    "used in Processor.get_edges() for trimming edges for cutters"    
    # Set slice position
    if position == 0:
        check = slice(0, cutter.shape[0])
    elif position == 1:
        check = slice(seq.shape[0] - cutter.shape[0], seq.shape[0])
    ## TODO...
    else:
        raise TypeError("position argument must be an int in 0, 1, 2, 3")
        
    # Find matching bases
    matching = seq[:, check] == cutter
    
    # check for cutter at 5' read1
    mask = np.where(
        (seq[:, check] != 78) & (seq[:, check] != 45) 
    )
    
    # subset matching bases to include only unmasked
    matchset = matching[mask]
    
    # allow for N percent matching
    if matchset.sum() / matchset.size <= cutoff:
        return 0
    return cutter.shape[0]


# ---------------------------------------
# jitted Processor functions
# ---------------------------------------
@njit
def maxind_numba(block):
    "count the max number of internal indels"
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
def snpcount_numba(block, snpstring):
    "Used to count the number of unique bases in a site for snpstring."  

    # iterate over all loci
    for site in range(block.shape[1]):

        # make new array
        catg = np.zeros(4, dtype=np.int16)

        # a list for only catgs
        ncol = block[:, site]
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
                snpstring[site, 1] = True
            else:
                snpstring[site, 0] = True
    return snpstring


@njit
def maxhet_numba(block, maxhet):
    counts = np.zeros(block.shape[1], dtype=np.int16)
    for fidx in range(block.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(block[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts
