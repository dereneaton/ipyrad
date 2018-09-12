#!/usr/bin/env python

# py2/3 compatibility
from __future__ import print_function
try:
    from builtins import range
    from itertools import izip
except ImportError:
    izip = zip

# standard lib imports
import os
import glob
import time
import shutil
import pickle
from collections import Counter

# third party imports
import numpy as np
import pandas as pd
from numba import njit
from .utils import IPyradError, clustdealer, splitalleles

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


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

        # dict mapping of samples to padded names for loci file aligning.
        self.data.snames = sorted(list(self.data.samples.keys()))
        self.data.pnames, self.data.snppad = self.get_padded_names()      

    def run(self):
        # split clusters into bits.
        self.split_clusters()

        # get filter and snp info on edge trimmed data.
        # write to chunks for building output files and save dimensions.
        self.remote_process_chunks()

        # write stats file while counting nsnps and nbases.
        # runs pretty fast
        self.collect_stats()

        # can be sent to remote.
        self.write_loci_and_alleles()

        # iterate over new processed chunks to fill the h5 array.
        # single-threaded
        self.fill_seq_array()

        # write outputs -- all parallelizable after phy array is built
        #self.remote_build_conversions()


        #self.remote_write_outfiles()

        #self.write_phylip()
        #self.write_nexus()

        # write outputs -- all parallelizable after snp array is built
        #self.write_snps()
        #self.write_structure()


        # build loci and alleles files stats 
        # build arrays and fill in second iteration over chunks.
        #self.write_loci()

        # and fill the seq and snp arrays.       
        #self.build_loci()
        #self.build_alleles()

        # build conversions 
        #self.remote_fill_seqarrs()
        #self.remote_build_conversions()

        # send jobs to build vcf
        #self.remote_fill_depths()
        #self.remote_build_vcf()

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

    def write_loci_and_alleles(self):

        # write alleles file
        allel = 'a' in self.data.paramsdict["output_formats"]

        # store output handle
        self.data.outfiles.loci = os.path.join(
            self.data.dirs.outfiles, self.data.name + ".loci")
        if allel:
            self.data.outfiles.alleles = os.path.join(
                self.data.dirs.outfiles, self.data.name + ".alleles.loci")

        # gather all loci bits
        locibits = glob.glob(os.path.join(self.data.tmpdir, "*.loci"))
        sortbits = sorted(locibits, 
            key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

        # what is the length of the name padding?
        with open(sortbits[0], 'r') as test:
            pad = np.where(np.array(list(test.readline())) == " ")[0].max()

        # write to file while adding counters to the ordered loci
        # TODO: this will need to grab chrom info on reference.
        outloci = open(self.data.outfiles.loci, 'w')
        if allel:
            outalleles = open(self.data.outfiles.alleles, 'w')

        for bit in sortbits:
            # store until writing
            lchunk = []
            achunk = []
            idx = 0

            # LOCI ONLY: iterate through chunk files
            if not allel:
                for bit in sortbits:
                    for line in iter(open(bit, 'r')):
                        if "|\n" not in line:
                            lchunk.append(line[:pad] + line[pad:].upper())
                        else:
                            lchunk.append(
                                "{}|{}|\n".format(line.rsplit("|", 2)[0], idx))
                            idx += 1

            # ALLELES: iterate through chunk files
            else:
                for bit in sortbits:
                    for line in iter(open(bit, 'r')):
                        if "|\n" not in line:
                            name = line[:pad]
                            seq = line[pad:]
                            lchunk.append(name + seq.upper())

                            all1, all2 = splitalleles(seq)
                            aname, spacer = name.split(" ", 1)
                            achunk.append(aname + "_0 " + spacer + all1)
                            achunk.append(aname + "_1 " + spacer + all2)
                        else:
                            lchunk.append(
                                "{}|{}|\n".format(line.rsplit("|", 2)[0], idx))
                            achunk.append(
                                "{}|{}|\n".format(line.rsplit("|", 2)[0], idx))
                            idx += 1
                outalleles.write("".join(achunk))
            outloci.write("".join(lchunk))
        outloci.close()
        if allel:
            outalleles.close()

    def fill_seq_array(self):
       
        # init/reset hdf5 database
        with h5py.File(self.data.database, 'w') as io5:
            io5.attrs["samples"] = [i.name.encode() for i in self.samples]

            # temporary array data sets 
            io5.create_dataset(
                name="phy",
                shape=(self.ntaxa, self.nbases), 
                dtype=np.uint8,
            )
            # temporary array data sets 
            io5.create_dataset(
                name="phymap",
                shape=(self.nloci,),
                dtype=np.int64,
            )

            # gather all loci bits
            locibits = glob.glob(os.path.join(self.data.tmpdir, "*.loci"))
            sortbits = sorted(locibits, 
                key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

            # name order for entry in array
            sidxs = {sample: i for (i, sample) in enumerate(self.data.samples)}

            # iterate through file
            gstart = 0
            start = end = 0
            maxsize = 100000
            tmploc = {}
            maplist = []
            mapstart = mapend = 0
            locidx = 0

            # array to store until writing
            tmparr = np.zeros((self.ntaxa, maxsize + 5000), dtype=np.uint8)
            
            # iterate over chunkfiles
            for bit in sortbits:
                # iterate lines of file until locus endings
                for line in iter(open(bit, 'r')):
                    
                    # still filling locus until |\n
                    if "|\n" not in line:
                        name, seq = line.split()
                        tmploc[name] = seq

                    # locus is full, dump it
                    else:
                        # convert seqs to an array
                        locidx += 1
                        loc = np.array([list(i) for i in tmploc.values()]).astype(bytes).view(np.uint8)
                        
                        # drop the site that are all N or -
                        mask = np.all((loc == 0) | (loc == 71), axis=0)
                        loc = loc[:, ~mask]
                        
                        # store end position of locus for map
                        end = start + loc.shape[1]
                        for idx, name in enumerate(tmploc):
                            tmparr[sidxs[name], start:end] = loc[idx]
                        maplist.append(gstart + end)
                        
                        # reset locus
                        start = end
                        tmploc = {}
                        
                    # dump tmparr when it gets large
                    if end > maxsize:
                        
                        # trim right overflow from tmparr
                        trim = np.where(tmparr != 0)[1]
                        if trim.size:
                            trim = trim.max() + 1
                        else:
                            trim = tmparr.shape[1]

                        # fill missing with 78 (N)
                        tmparr[tmparr == 0] = 78
                        
                        # dump tmparr to hdf5
                        io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]                       
                        io5['phymap'][mapstart:locidx] = np.array(maplist, dtype=np.int64)
                        mapstart = locidx
                        maplist = []
                        
                        # reset
                        tmparr = np.zeros((self.ntaxa, maxsize + 5000), dtype=np.uint8)
                        gstart += trim
                        start = end = 0
                        
            # write final chunk
            trim = np.where(tmparr != 0)[1]
            if trim.size:
                trim = trim.max() + 1
            else:
                trim = tmparr.shape[1]

            # fill missing with 78 (N)
            tmparr[tmparr == 0] = 78

            # dump tmparr to hdf5
            io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]           
            mapend = mapstart + len(maplist)
            io5['phymap'][mapstart:mapend] = np.array(maplist, dtype=np.int64)

    def fill_snp_array(self):
       
        # open new database file handle
        with h5py.File(self.data.database, 'a') as io5:

            # temporary array data sets 
            io5.create_dataset(
                name="snps",
                shape=(self.ntaxa, self.nsnps),
                dtype=np.uint8,
            )
            # temporary array data sets 
            io5.create_dataset(
                name="snpsmap",
                shape=(self.nsnps, 4),
                dtype=np.uint32,
            )

            # gather all loci bits
            locibits = glob.glob(os.path.join(self.data.tmpdir, "*.loci"))
            sortbits = sorted(locibits, 
                key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

            # name order for entry in array
            sidxs = {sample: i for (i, sample) in enumerate(self.data.samples)}

            # iterate through file
            gstart = 0
            start = end = 0
            maxsize = 100000
            tmploc = {}
            maplist = []
            mapstart = mapend = 0

            # array to store until writing
            tmparr = np.zeros((self.ntaxa, maxsize + 5000), dtype=np.uint8)

            # iterate over chunkfiles
            for bit in sortbits:
                # iterate lines of file until locus endings
                for line in iter(open(bit, 'r')):
                    
                    # still filling locus until |\n
                    if "|\n" not in line:
                        name, seq = line.split()
                        tmploc[name] = seq

                    # locus is full, dump it
                    else:
                        # convert seqs to an array
                        loc = np.array([list(i) for i in tmploc.values()]).astype(bytes).view(np.uint8)
                        snps = line[len(self.data.snppad):].rsplit("|", 1)[0]
                        snpsarr = np.array(list(snps)) != ""

                        # select only the SNP sites
                        snpsites = loc[snpsarr]
                        print(snpsites)
                        
                        # store end position of locus for map
                        end = start + loc.shape[1]
                        for idx, name in enumerate(tmploc):
                            tmparr[sidxs[name], start:end] = loc[idx]
                        maplist.append(end)
                        
                        # reset locus
                        start = end
                        tmploc = {}
                        
                    # dump tmparr when it gets large
                    if end > maxsize:
                        
                        # trim right overflow from tmparr
                        trim = np.where(tmparr != 0)[1]
                        if trim.size:
                            trim = trim.max() + 1
                        else:
                            trim = tmparr.shape[1]

                        # fill missing with 78 (N)
                        tmparr[tmparr == 0] = 78
                        
                        # dump tmparr to hdf5
                        io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]
                        
                        mapend = mapstart + len(maplist)
                        io5['phymap'][mapstart:mapend] = np.array(maplist, dtype=np.int64)
                        mapstart += mapend
                        maplist = []
                        
                        # reset
                        tmparr = np.zeros((self.ntaxa, maxsize + 5000), dtype=np.uint8)
                        gstart += trim
                        start = end = 0
                        
            # write final chunk
            trim = np.where(tmparr != 0)[1]
            if trim.size:
                trim = trim.max() + 1
            else:
                trim = tmparr.shape[1]

            # fill missing with 78 (N)
            tmparr[tmparr == 0] = 78

            # dump tmparr to hdf5
            io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]
            #print(end, gstart, gstart + trim)
            
            mapend = mapstart + len(maplist)
            io5['phymap'][mapstart:mapend] = np.array(maplist, dtype=np.int64)

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
            args = (self.data, self.chunks, jobfile)
            rasyncs[jobfile] = self.lbview.apply(process_chunks, *args)
        
        # iterate until all chunks are processed
        while 1:
            # get and enter results into hdf5 as they come in
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                break          

        # write stats
        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].exception())

    def collect_stats(self):
        "Collect results from Processor and write stats file."

        # organize stats into dataframes 
        ftable = pd.DataFrame(
            columns=["total_filters", "applied_order", "retained_loci"],
            index=[
                "total_prefiltered_loci",
                "filtered_by_rm_duplicates",
                "filtered_by_max_indels",
                "filtered_by_max_SNPs",
                "filtered_by_max_shared_het",
                "filtered_by_min_sample",  # "filtered_by_max_alleles",
                "total_filtered_loci"],
        )

        # load pickled dictionaries into a dict
        pickles = glob.glob(os.path.join(self.data.tmpdir, "*.p"))
        pdicts = {}
        for pkl in pickles:
            with open(pkl, 'rb') as inp:
                pdicts[pkl.rsplit("-", 1)[-1][:-2]] = pickle.load(inp)

        # join dictionaries into global stats
        afilts = np.concatenate([i['filters'] for i in pdicts.values()])
        lcovs = Counter({})
        scovs = Counter({})
        cvar = Counter({})
        cpis = Counter({})
        nbases = 0
        for lcov in [i['lcov'] for i in pdicts.values()]:
            lcovs.update(lcov)
        for scov in [i['scov'] for i in pdicts.values()]:
            scovs.update(scov)
        for var in [i['var'] for i in pdicts.values()]:
            cvar.update(var)
        for pis in [i['pis'] for i in pdicts.values()]:
            cpis.update(pis)
        for count in [i['nbases'] for i in pdicts.values()]:
            nbases += count

        # make into nice DataFrames
        ftable.iloc[0, :] = (0, 0, self.nraws)

        # filter rm dups
        ftable.iloc[1, 0:2] = afilts[:, 0].sum()
        ftable.iloc[1, 2] = ftable.iloc[0, 2] - ftable.iloc[1, 1]
        mask = afilts[:, 0]

        # filter max indels
        ftable.iloc[2, 0] = afilts[:, 1].sum()
        ftable.iloc[2, 1] = afilts[~mask, 1].sum()
        ftable.iloc[2, 2] = ftable.iloc[1, 2] - ftable.iloc[2, 1]
        mask = afilts[:, 0:2].sum(axis=1).astype(np.bool)

        # filter max snps
        ftable.iloc[3, 0] = afilts[:, 2].sum()
        ftable.iloc[3, 1] = afilts[~mask, 2].sum()
        ftable.iloc[3, 2] = ftable.iloc[2, 2] - ftable.iloc[3, 1]
        mask = afilts[:, 0:3].sum(axis=1).astype(np.bool)

        # filter max shared H
        ftable.iloc[4, 0] = afilts[:, 3].sum()
        ftable.iloc[4, 1] = afilts[~mask, 3].sum()
        ftable.iloc[4, 2] = ftable.iloc[3, 2] - ftable.iloc[4, 1]
        mask = afilts[:, 0:4].sum(axis=1).astype(np.bool)

        # filter minsamp
        ftable.iloc[5, 0] = afilts[:, 4].sum()
        ftable.iloc[5, 1] = afilts[~mask, 4].sum()
        ftable.iloc[5, 2] = ftable.iloc[4, 2] - ftable.iloc[5, 1]
        mask = afilts[:, 0:4].sum(axis=1).astype(np.bool)

        ftable.iloc[6, 0] = ftable.iloc[:, 0].sum()
        ftable.iloc[6, 1] = ftable.iloc[:, 1].sum()
        ftable.iloc[6, 2] = ftable.iloc[5, 2]

        # save stats to the data object
        self.data.stats_dfs.s7_filters = ftable
        self.data.stats_dfs.s7_samples = pd.DataFrame(
            pd.Series(scovs, name="sample_coverage"))

        ## get locus cov and sums 
        lrange = range(1, len(self.samples) + 1)
        covs = pd.Series(lcovs, name="locus_coverage", index=lrange)
        start = self.data.paramsdict["min_samples_locus"] - 1
        sums = pd.Series(
            {i: np.sum(covs[start:i]) for i in lrange},            
            name="sum_coverage", 
            index=lrange)
        self.data.stats_dfs.s7_loci = pd.concat([sums, covs], axis=1)

        ## get SNP distribution       
        maxsnps = sum(self.data.paramsdict['max_SNPs_locus'])
        sumd = {}
        sump = {}
        for i in range(maxsnps):
            sumd[i] = np.sum([i * cvar[i] for i in range(i + 1)])
            sump[i] = np.sum([i * cpis[i] for i in range(i + 1)])        
        self.data.stats_dfs.s7_snps = pd.concat([
            pd.Series(cvar, name="var"),
            pd.Series(sumd, name="sum_var"),
            pd.Series(cpis, name="pis"),
            pd.Series(sump, name="sum_pis"),           
            ],          
            axis=1
        )
        ## trim SNP distribution to exclude unobserved endpoints
        varmin = (self.data.stats_dfs.s7_snps['var'] != 0).idxmin()
        pismin = (self.data.stats_dfs.s7_snps['pis'] != 0).idxmin()
        amin = max([varmin, pismin])
        self.data.stats_dfs.s7_snps = self.data.stats_dfs.s7_snps.iloc[:amin]

        ## store dimensions for array building 
        self.nloci = ftable.iloc[6, 2]
        self.nbases = nbases
        self.nsnps = self.data.stats_dfs.s7_snps["sum_var"].max()
        self.ntaxa = len(self.samples)

        # write to file
        self.data.stats_files.s7 = os.path.join(
            self.data.dirs.outfiles, "{}_stats.txt".format(self.data.name))
        with open(self.data.stats_files.s7, 'w') as outstats:
            print(STATS_HEADER_1, file=outstats)
            self.data.stats_dfs.s7_filters.to_string(buf=outstats)

            print(STATS_HEADER_2, file=outstats)
            self.data.stats_dfs.s7_samples.to_string(buf=outstats)

            print(STATS_HEADER_3, file=outstats)
            self.data.stats_dfs.s7_loci.to_string(buf=outstats)

            print(STATS_HEADER_4, file=outstats)
            self.data.stats_dfs.s7_snps.to_string(buf=outstats)

            print("\n\n\nFinal Sample stats summary", file=outstats)
            statcopy = self.data.stats.copy()
            statcopy.state = 7
            statcopy['loci_in_assembly'] = self.data.stats_dfs.s7_samples
            statcopy.to_string(buf=outstats)

    def get_padded_names(self):
        # get longest name
        longlen = max(len(i) for i in self.data.snames)
        # Padding distance between name and seq.
        padding = 5
        # add pad to names
        pnames = {
            name: "{}{}".format(name, " " * (longlen - len(name) + padding))
            for name in self.data.snames
        }
        snppad = "//" + " " * (longlen - 2 + padding)
        return pnames, snppad

    def remote_write_outfiles(self):

        oformats = ["phy", "nex"]
        for oformat in oformats:

            # store handle to data object
            self.data.outfiles[oformat] = os.path.join(
                self.data.dirs.outfiles,
                self.data.name + oformat)

            # run conversion on remote engine
            #args = ()
            #rasync = self.lbview.apply()
            convert_output(self.data, oformat)


# ------------------------------------------------------------
# Classes initialized and run on remote engines.
# ------------------------------------------------------------

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
            'dups', 
            'maxind',  
            'maxvar', 
            # 'maxall',             
            'maxshared',
            'minsamp', 
            )
        self.edges = np.zeros((self.chunks, 4), dtype=np.uint16)

        # store stats on sample coverage and locus coverage
        self.maxsnps = sum(self.data.paramsdict['max_SNPs_locus'])
        self.scov = {i: 0 for i in self.data.snames}
        self.lcov = {i: 0 for i in range(1, len(self.data.snames) + 1)}
        self.var = {i: 0 for i in range(self.maxsnps)}
        self.pis = {i: 0 for i in range(self.maxsnps)}
        self.nbases = 0

        # tmp outfile list and filename
        self.outlist = []
        self.outfile = self.chunkfile + '.loci'
        self.outarr = self.chunkfile + '.npy'
        self.outpickle = self.chunkfile + '.p'

    def run(self):

        # store list of edge trims for VCF building
        edgelist = []
        ns = len(self.data.snames)

        # todo: this could be an iterator...
        with open(self.chunkfile, 'rb') as infile:
            loci = infile.read().split(b"//\n//\n")
    
            # iterate over loci
            for iloc, loc in enumerate(loci):                              
                # load names and seqs 
                lines = loc.decode().strip().split("\n")
                names = []
                nidxs = []
                aseqs = []
                useqs = []
                for line in lines:
                    if line[0] == ">":
                        name, nidx = line[1:].rsplit("_", 1)
                        names.append(name)
                        nidxs.append(nidx)
                    else:
                        aseqs.append(list(line))
                        useqs.append(list(line.upper()))

                #names = [i[1:].rsplit("_", 1)[0] for i in lines[::2]]
                #name_indices = [i[1:].rsplit("_", 1)[1] for i in lines[::2]]
                #seqs = np.array([list(i) for i in lines[1::2]])
                #seqs = np.array([list(i.upper()) for i in lines[1::2]])
                #seqs = seqs.astype(bytes).view(np.uint8)

                # filter to only include only samples in this assembly
                mask = [i in self.data.snames for i in names]
                names = np.array(names)[mask].tolist()
                nidxs = np.array(nidxs)[mask].tolist()
                useqs = np.array(useqs)[mask, :].astype(bytes).view(np.uint8)
                aseqs = np.array(aseqs)[mask, :].astype(bytes).view(np.uint8)
                #aseqs = np.array(seqs)[mask, :].astype(bytes).view(np.uint8)                
                
                # apply filters
                efilter, edges = self.get_edges(useqs)
                self.edges[iloc] = edges
                self.filters[iloc, 0] = self.filter_dups(names)
                self.filters[iloc, 1] = self.filter_minsamp_pops(names)
                self.filters[iloc, 1] += efilter

                # should we fill terminal indels as N's here?
                #...

                # trim edges, need to use uppered seqs for maxvar & maxshared
                edg = self.edges[iloc]
                ublock1 = useqs[:, edg[0]:edg[1]]
                ublock2 = useqs[:, edg[2]:edg[3]]
                block1 = aseqs[:, edg[0]:edg[1]]
                block2 = aseqs[:, edg[2]:edg[3]]

                # apply filters on edge trimmed reads
                self.filters[iloc, 2] += self.filter_maxindels(block1, block2)

                # get snpstring on trimmed reads
                snparr1, snparr2 = self.get_snpsarrs(ublock1, ublock2)
                self.filters[iloc, 3] = self.filter_maxvars(snparr1, snparr2)

                # apply filters on edge trimmed reads
                self.filters[iloc, 4] = self.filter_maxshared(ublock1, ublock2)
                #self.filters[iloc, 5] = self.filter_maxalleles()

                # store stats for the locus that passed filtering
                if not self.filters[iloc, :].sum():                   
                    # do sample and locus counters
                    for name in names:
                        self.scov[name] += 1
                    self.lcov[useqs.shape[0]] += 1             

                    # do SNP distribution counter
                    if snparr2.size:
                        snps = np.concatenate([snparr1, snparr2])
                        self.nbases += ublock1.shape[1] + ublock1.shape[1]
                    else:
                        snps = snparr1
                        self.nbases += ublock1.shape[1]
                    self.var[snps[:, 0].sum()] += 1
                    self.pis[snps[:, 1].sum()] += 1                   

                    # write to .loci string
                    locus = self.to_locus(
                        names, nidxs, block1, block2, snparr1, snparr2,
                    )
                    self.outlist.append(locus)

                    # if VCF: store information on edge trimming so we can line
                    # up depth information from catgs.
                    #if "V" in ...
                    edgearr = np.zeros(ns, dtype=np.uint8)
                    for name, row in zip(names, aseqs):
                        sidx = self.data.snames.index(name)
                        trim = np.where(row != 45)[0]
                        if trim.size:
                            edgearr[sidx] = trim.min()
                        edgelist.append(edgearr)

        # write the chunk to tmpdir
        with open(self.outfile, 'w') as outchunk:
            outchunk.write("\n".join(self.outlist) + "\n")

        # convert edgelist to an array and save as .npy
        np.save(self.outarr, np.array(edgelist))


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
                        self.data.pnames[name],
                        block1[idx, :].tostring().decode())
                )
            locus.append("{}{}|{}|".format(
                self.data.snppad, snpstring1, nidxstring))
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
        "max internal indels. Denovo vs. Ref, single versus paired."
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

    def get_snpsarrs(self, block1, block2):
        if "pair" not in self.data.paramsdict["datatype"]:
            snpsarr1 = np.zeros((block1.shape[1], 2), dtype=np.bool_)
            snpsarr1 = snpcount_numba(block1, snpsarr1)
            snpsarr2 = np.array([])
        else:
            snpsarr1 = np.zeros((block1.shape[1], 2), dtype=np.bool_)
            snpsarr1 = snpcount_numba(block1, snpsarr1)
            snpsarr2 = np.zeros((block2.shape[1], 2), dtype=np.bool_)
            snpsarr2 = snpcount_numba(block2, snpsarr2)
        return snpsarr1, snpsarr2


class Converter:
    "functions for converting hdf5 arrays into output files"
    def __init__(self, data):
        self.data = data
        self.output_formats = self.data.paramsdict["output_formats"]
        self.database = self.data.database

    def run(self, oformat):
        pass


    def write_phy(self):
        # write from hdf5 array
        with open(self.data.outfiles.phy, 'w') as out:
            with h5py.File(self.data.database, 'r') as io5:
                # load seqarray
                seqarr = io5['phy'][:]
                arrsize = io5['phymap'][-1]

                # write dims
                out.write("{} {}\n".format(self.ntaxa, arrsize))

                # write to disk
                for idx in range(io5['phy'].shape[0]):
                    seq = seqarr[idx, :].view("S1")
                    out.write(
                        "{}{}".format(
                            self.data.pnames[self.data.snames[idx]],
                            b"".join(seq).decode().upper() + "\n",
                        )
                    )
                    
                    
    def write_nex(self):
        # write from hdf5 array
        with open(self.data.outfiles.nex, 'w') as out:
            with h5py.File(self.data.database, 'r') as io5:
                # load seqarray
                seqarr = io5['phy'][:]
                arrsize = io5['phymap'][-1]

                ## write nexus seq header
                out.write(NEXHEADER.format(seqarr.shape[0], arrsize))

                ## grab a big block of data
                chunksize = 100000  # this should be a multiple of 100
                for bidx in range(0, arrsize, chunksize):
                    bigblock = seqarr[:, bidx:bidx + chunksize]
                    lend = arrsize - bidx

                    ## write interleaved seqs 100 chars with longname+2 before
                    tmpout = []            
                    for block in range(0, min(chunksize, lend), 100):
                        stop = min(block + 100, arrsize)

                        for idx, name in enumerate(self.data.pnames):
                            seqdat = bigblock[idx, block:stop]
                            tmpout.append("  {}{}\n".format(
                                self.data.pnames[name], 
                                b"".join(seqdat).decode().upper()))
                        tmpout.append("\n")

                    ## print intermediate result and clear
                    if any(tmpout):
                        out.write("".join(tmpout))
                ## closer
                out.write(NEXCLOSER)
                
                ## add partition information from maparr
                maparr = io5["phymap"][:]
                charsetblock = []
                charsetblock.append("BEGIN SETS;")
                for idx in range(0, maparr.shape[0] - 1):
                    charsetblock.append("charset {} = {}-{};".format(
                        idx, maparr[idx], maparr[idx + 1]
                        )
                    )
                charsetblock.append("END;")
                out.write("\n".join(charsetblock))


    def write_snps_map(self):
        """ write a map file with linkage information for SNPs file"""

        ## grab map data from tmparr
        tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
        with h5py.File(tmparrs, 'r') as io5:
            maparr = io5["maparr"][:]

            ## get last data 
            end = np.where(np.all(maparr[:] == 0, axis=1))[0]
            if np.any(end):
                end = end.min()
            else:
                end = maparr.shape[0]

            ## write to map file (this is too slow...)
            outchunk = []
            with open(data.outfiles.snpsmap, 'w') as out:
                for idx in range(end):
                    ## build to list
                    line = maparr[idx, :]
                    #print(line)
                    outchunk.append(
                        "{}\trad{}_snp{}\t{}\t{}\n"
                        .format(line[0], line[1], line[2], 0, line[3]))
                    ## clear list
                    if not idx % 10000:
                        out.write("".join(outchunk))
                        outchunk = []
                ## write remaining
                out.write("".join(outchunk))


class Vcfbuilder:
    def __init__(self):
        pass


# -----------------------------------------------------------
# Step7 external functions that are run on engines
# -----------------------------------------------------------
def process_chunks(data, chunks, chunkfile):
    # process chunk writes to files and returns proc with features.
    proc = Processor(data, chunks, chunkfile)
    proc.run()

    # write process stats to a pickle file for collating later.
    out = {
        "filters": proc.filters, 
        "lcov": proc.lcov, 
        "scov": proc.scov,
        "var": proc.var,
        "pis": proc.pis,
        "nbases": proc.nbases
    }
    with open(proc.outpickle, 'wb') as outpickle:
        pickle.dump(out, outpickle)


def convert_outputs(data, oformat):
    Converter(data).run(oformat)


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


# -------------------------------------------------------------
# jitted Processor functions
# -------------------------------------------------------------
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
def snpcount_numba(block, snpsarr):
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
                snpsarr[site, 1] = True
            else:
                snpsarr[site, 0] = True
    return snpsarr


@njit
def maxhet_numba(block, maxhet):
    counts = np.zeros(block.shape[1], dtype=np.int16)
    for fidx in range(block.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(block[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts


# -----------------------------------------------------------
# GLOBALS
# -----------------------------------------------------------
AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)
STATS_HEADER_1 = """
## The number of loci caught by each filter.
## ipyrad API location: [assembly].stats_dfs.s7_filters
"""
STATS_HEADER_2 = """\n\n
## The number of loci recovered for each Sample.
## ipyrad API location: [assembly].stats_dfs.s7_samples
"""
STATS_HEADER_3 = """\n\n
## The number of loci for which N taxa have data.
## ipyrad API location: [assembly].stats_dfs.s7_loci
"""
STATS_HEADER_4 = """\n\n
The distribution of SNPs (var and pis) per locus.
## var = Number of loci with n variable sites (pis + autapomorphies)
## pis = Number of loci with n parsimony informative site (minor allele in >1 sample)
## ipyrad API location: [assembly].stats_dfs.s7_snps
"""
NEXHEADER = """#nexus
begin data;
  dimensions ntax={} nchar={};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""
NEXCLOSER = """  ;
end;
"""
