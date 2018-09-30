#!/usr/bin/env python

# py2/3 compatibility
from __future__ import print_function
try:
    from builtins import range
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

# standard lib imports
import os
import sys
import glob
import time
import shutil
import pickle
from collections import Counter

# third party imports
import numpy as np
import pandas as pd
import ipyrad
from numba import njit
from .utils import IPyradError, clustdealer, splitalleles
from .utils import BTS, GETCONS, DCONS


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
        self.data.isref = bool("ref" in self.data.paramsdict["assembly_method"])

        # returns samples in the order we want them in the outputs
        self.samples = self.get_subsamples()
        self.setup_dirs()
        self.get_chunksize()
        #self.chroms2ints()

        # dict mapping of samples to padded names for loci file aligning.
        self.data.snames = [i.name for i in self.samples]
        self.data.pnames, self.data.snppad = self.get_padded_names()      

        # output file formats to produce
        self.formats = set(['l']).union(
            set(self.data.paramsdict["output_formats"]))


    def run(self):
        # split clusters into bits.
        self.split_clusters()

        # get filter and snp info on edge trimmed data.
        # write to chunks for building output files and save dimensions.
        self.remote_process_chunks()

        # write stats file while counting nsnps and nbases.
        self.collect_stats()
        self.store_file_handles()

        # write loci and alleles outputs (parallelized on 3 engines)
        self.remote_build_arrays_and_write_loci()

        # send conversion jobs from array files to engines
        self.remote_write_outfiles()

        # send jobs to build vcf
        if 'v' in self.formats:
            self.remote_fill_depths()
            self.remote_build_vcf()    

    ## init functions ------
    def get_subsamples(self):
        "get subsamples for this assembly. All must have been in step6"

        # get samples from the database file
        if not os.path.exists(self.data.clust_database):
            raise IPyradError("You must first complete step6.")
        with open(self.data.clust_database, 'r') as inloci:
            dbsamples = inloci.readline()[1:].strip().split(",@")

        # samples are in this assembly but not database (raise error)
        nodb = set(self.data.samples).difference(set(dbsamples))
        if nodb:
            raise IPyradError(MISSING_SAMPLE_IN_DB.format(nodb))

        # samples in database not in this assembly, that's OK, you probably
        # branched to drop some samples. 

        # samples in populations file that are not in this assembly. Raise 
        # an error, it's probably a typo and should be corrected. 
        poplists = [i[1] for i in self.data.populations.values()]
        popset = set(chain(*poplists))
        badpop = popset.difference(set(self.data.samples))
        if badpop:
            raise IPyradError(BADPOP_SAMPLES.format(badpop))

        # output files already exist for this assembly (check stats). Raise
        # error unless using the force flag to prevent overwriting. 
        self.data.stats_files.s7 = os.path.abspath(
            os.path.join(
                self.data.dirs.outfiles, 
                "{}_stats.txt".format(self.data.name),
            )
        )
        if not self.force:
            if os.path.exists(self.data.stats_files.s7):
                raise IPyradError(
        "Step 7 results already exist for this Assembly. Use force to overwrite.")

        # if ref init a new sample for reference if including
        if self.data.paramsdict['assembly_method'] == 'reference':
            ref = ipyrad.Sample("reference")
            snames = [ref] + sorted(
                list(set(self.data.samples.values())), 
                key=lambda x: x.name)
            return snames

        snames = sorted(
            list(set(self.data.samples.values())),
            key=lambda x: x.name)
        return snames

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

        # make new database files
        self.data.seqs_database = os.path.join(
            self.data.dirs.outfiles,
            self.data.name + ".seqs.hdf5",
            )
        self.data.snps_database = os.path.join(
            self.data.dirs.outfiles,
            self.data.name + ".snps.hdf5",
            )
        for dbase in [self.data.snps_database, self.data.seqs_database]:
            if os.path.exists(dbase):
                os.remove(dbase)

    def get_chunksize(self):
        "get nloci and ncpus to chunk and distribute work across processors"
        # this file is inherited from step 6 to allow step7 branching.
        with open(self.data.clust_database, 'r') as inloci:
            # skip header
            inloci.readline()
            # get nraw loci
            self.nraws = sum(1 for i in inloci if i == "//\n") // 2

        # chunk to approximately 2 chunks per core
        self.ncpus = len(self.ipyclient.ids)
        self.chunks = ((self.nraws // (self.ncpus * 2)) + \
                       (self.nraws % (self.ncpus * 2)))

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

    ## core functions ------
    def store_file_handles(self):
        # always produce a .loci file + whatever they ask for.
        testformats = list(self.formats)
        for outf in testformats:

            # if it requires a pop file and they don't have one then skip
            # and print a warning:
            if (outf in ("t", "m")) and (not self.data.populations):
                print(POPULATION_REQUIRED.format(outf), file=sys.stderr)

                # remove format from the set
                self.formats.discard(outf)
                continue

            else:                
                # store handle to data object
                for ending in OUT_SUFFIX[outf]:
                    
                    # store 
                    self.data.outfiles[ending[1:]] = os.path.join(
                        self.data.dirs.outfiles,
                        self.data.name + ending)           

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
        self.data.stats_dfs.s7_loci = pd.concat([covs, sums], axis=1)

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
        with open(self.data.stats_files.s7, 'w') as outstats:
            print(STATS_HEADER_1, file=outstats)
            self.data.stats_dfs.s7_filters.to_string(buf=outstats)

            print(STATS_HEADER_2, file=outstats)
            self.data.stats_dfs.s7_samples.to_string(buf=outstats)

            print(STATS_HEADER_3, file=outstats)
            self.data.stats_dfs.s7_loci.to_string(buf=outstats)

            print(STATS_HEADER_4, file=outstats)
            self.data.stats_dfs.s7_snps.to_string(buf=outstats)

            print("\n\n\n## Final Sample stats summary", file=outstats)
            statcopy = self.data.stats.copy()
            statcopy.state = 7
            statcopy['loci_in_assembly'] = self.data.stats_dfs.s7_samples
            statcopy.to_string(buf=outstats)
            print("\n\n\n## Alignment matrix statistics:", file=outstats)

    def split_clusters(self):
        with open(self.data.clust_database, 'rb') as clusters:
            # skip header
            clusters.readline()

            # build iterator
            pairdealer = izip(*[iter(clusters)] * 2)

            # grab a chunk of clusters
            idx = 0
            while 1:

                # if an engine is available pull off a chunk
                try:
                    done, chunk = clustdealer(pairdealer, self.chunks)
                except IndexError:
                    raise IPyradError(
                        "clust_database formatting error in %s", chunk)

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
                raise IPyradError(rasyncs[job].get())

    def remote_build_arrays_and_write_loci(self):
        # start loci concatenating job on a remote
        start = time.time()
        printstr = ("building arrays     ", "s7")
        rasyncs = {}
        args0 = (self.data,)
        args1 = (self.data, self.ntaxa, self.nbases, self.nloci)
        args2 = (self.data, self.ntaxa, self.nsnps)
        rasyncs[0] = self.lbview.apply(write_loci_and_alleles, *args0)
        rasyncs[1] = self.lbview.apply(fill_seq_array, *args1)
        rasyncs[2] = self.lbview.apply(fill_snp_array, *args2)
        # track progress.
        while 1:
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                break
        # check for errors
        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].get())

    def remote_write_outfiles(self):
        "Calls Converter object funcs in parallel."
        start = time.time()
        printstr = ("writing conversions ", "s7")        
        rasyncs = {}

        for outf in self.formats:
            rasyncs[outf] = self.lbview.apply(
                convert_outputs, *(self.data, outf))

        # iterate until all chunks are processed
        while 1:
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                break          

        # write stats
        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].get())
                
    def remote_fill_depths(self):
        "send each sample to build depth arrays"
        start = time.time()
        printstr = ("indexing vcf depths ", "s7")        
        rasyncs = {}

        for sample in self.data.samples.values():
            if not sample.name == "reference":
                rasyncs[sample.name] = self.lbview.apply(
                    fill_vcf_depths, *(self.data, self.nsnps, sample))

        # iterate until all chunks are processed
        while 1:
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                break          

        # write stats
        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].get())

    def remote_build_vcf(self):
        "write VCF file"
        start = time.time()
        printstr = ("writing vcf output  ", "s7")        
        rasync = self.lbview.apply(build_vcf, self.data)          

        # iterate until all chunks are processed
        while 1:
            ready = rasync.ready()
            self.data._progressbar(1, ready, start, printstr)
            time.sleep(0.5)
            if ready:
                break          

        # write stats
        print("")
        if not rasync.successful():
            raise IPyradError(rasync.get())


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
        self.isref = data.paramsdict["assembly_method"] == "reference"

        # filters (dups, minsamp, maxind, maxall, maxvar, maxshared)
        self.filters = np.zeros((self.chunks, 5), dtype=np.bool_)
        self.filterlabels = (
            'dups', 
            'maxind',  
            'maxvar', 
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

                # filter to only include only samples in this assembly
                mask = [i in self.data.snames for i in names]
                names = np.array(names)[mask].tolist()
                nidxs = np.array(nidxs)[mask].tolist()
                useqs = np.array(useqs)[mask, :].astype(bytes).view(np.uint8)
                aseqs = np.array(aseqs)[mask, :].astype(bytes).view(np.uint8)
                
                # apply filters
                efilter, edges = self.get_edges(useqs)
                self.edges[iloc] = edges
                self.filters[iloc, 0] = self.filter_dups(names)
                self.filters[iloc, 4] = self.filter_minsamp_pops(names)
                self.filters[iloc, 4] += efilter

                # should we fill terminal indels as N's here?
                #...

                # trim edges, need to use uppered seqs for maxvar & maxshared
                edg = self.edges[iloc]
                ublock1 = useqs[:, edg[0]:edg[1]]
                ublock2 = useqs[:, edg[2]:edg[3]]
                block1 = aseqs[:, edg[0]:edg[1]]
                block2 = aseqs[:, edg[2]:edg[3]]

                # apply filters on edge trimmed reads
                self.filters[iloc, 1] += self.filter_maxindels(block1, block2)

                # get snpstring on trimmed reads
                snparr1, snparr2 = self.get_snpsarrs(ublock1, ublock2)
                self.filters[iloc, 2] = self.filter_maxvars(snparr1, snparr2)

                # apply filters on edge trimmed reads
                self.filters[iloc, 3] = self.filter_maxshared(ublock1, ublock2)

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
                    self.var[snps[:, :].sum()] += 1
                    self.pis[snps[:, 1].sum()] += 1                   

                    # write to .loci string
                    locus = self.to_locus(
                        names, nidxs, block1, block2, snparr1, snparr2, edg,
                    )
                    self.outlist.append(locus)

                    # if VCF: store edge trim amount so we can line
                    edgelist.append(edges[0])

        # write the chunk to tmpdir
        with open(self.outfile, 'w') as outchunk:
            outchunk.write("\n".join(self.outlist) + "\n")

        # convert edgelist to an array and save as .npy
        np.save(self.outarr, np.array(edgelist))


    def to_locus(self, names, nidxs, block1, block2, snparr1, snparr2, edg):
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
        if self.isref:
            # get ref position from nidx 
            nidbits = [i.rsplit(":", 2)[0] for i in nidxs[1:]]
            refpos = ":".join(nidxs[0].rsplit(":", 2)[-2:])

            # trim ref position for edge trims
            chrom, pos = refpos.split(":")
            start, end = pos.split("-")
            start = int(start) + edg[0]
            end = start + (edg[3] - edg[0])

            # put back into string
            refpos = "{}:{}-{}".format(chrom, start, end)
            nidbits = [refpos] + nidbits
            nidxstring = ",".join(nidbits)

        else:
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
   
    def filter_minsamp_pops(self, names):
        if not self.data.populations:
            mins = self.data.paramsdict["min_samples_locus"]
            if len(names) < mins:
                return True
            return False

        else:
            minfilters = []
            for pop in self.data.populations:
                samps = self.data.populations[pop][1]
                minsamp = self.data.populations[pop][0]
                if len(set(samps).intersection(set(names))) < minsamp:
                    minfilters.append(pop)
            if any(minfilters):
                return True
            return False

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

        # get max snps for read1, read2
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

        # unless maxhet was set to 0, use 1 as a min setting
        if self.data.paramsdict["max_shared_Hs_locus"] != 0:
            maxhet = max(1, maxhet)

        # DEBUGGING
        # with open("/home/deren/test-indels.txt", 'a') as test:
        #     print("{} {}".format(maxhet_numba(blocks, maxhet).max(), maxhet),
        # file=test)
        #     for row in range(blocks.shape[0]):
        #         print(blocks[row, :].tostring().decode(), file=test)
        #     print("\n\n", file=test)

        # get max from combined block
        if maxhet_numba(blocks, maxhet).max() > maxhet:
            return True
        return False

    def get_edges(self, seqs):
        """
        Trim terminal edges or mask internal edges based on three criteria and
        take the max for each edge.
        1. user entered hard trimming.
        2. removing cutsite overhangs.
        3. trimming singleton-like overhangs from seqs of diff lengths.
        """

        # record whether to filter this locus based on sample coverage
        bad = False
        
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
        try:            
            for pos in range(4):
                trim3[pos] = check_minsamp(seqs, pos, 4)
        except ValueError:
            bad = True
        
        # get max edges
        trim = np.max([trim1, trim2, trim3], axis=0)

        # return edges as slice indices
        r1left = trim[0]
        r1right = seqs.shape[1] - trim[1]

        # TODO: resolve r2
        r2left = r2right = r1right
        edges = (r1left, r1right, r2left, r2right)
        
        # get filter
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
        self.seqs_database = self.data.seqs_database
        self.snps_database = self.data.snps_database        

    def run(self, oformat):

        # phy array outputs
        if oformat == "p":
            self.write_phy()

        # phy array + phymap outputs
        if oformat == "n":
            self.write_nex()
        
        if oformat == "G":
            self.write_gphocs()

        # phy array + phymap + populations outputs
        if oformat == "m":
            pass

        # snps array + snpsmap outputs
        if oformat == "s":
            self.write_snps()
            self.write_snps_map()

        # recommended to use analysis tools for unlinked sampling.
        #if oformat == "u":
        #    pass

        if oformat == "k":
            self.write_str()

        if oformat == "g":
            self.write_geno()


    def write_phy(self):
        # write from hdf5 array
        with open(self.data.outfiles.phy, 'w') as out:
            with h5py.File(self.seqs_database, 'r') as io5:
                # load seqarray
                seqarr = io5['phy']
                arrsize = io5['phymap'][-1]

                # write dims
                out.write("{} {}\n".format(len(self.data.snames), arrsize))

                # write to disk
                for idx in range(io5['phy'].shape[0]):
                    seq = seqarr[idx, :arrsize].view("S1")
                    out.write(
                        "{}{}".format(
                            self.data.pnames[self.data.snames[idx]],
                            b"".join(seq).decode().upper() + "\n",
                        )
                    )


    def write_nex(self):
        # write from hdf5 array
        with open(self.data.outfiles.nex, 'w') as out:
            with h5py.File(self.seqs_database, 'r') as io5:
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


    def write_snps(self):
        # write from hdf5 array
        with open(self.data.outfiles.snps, 'w') as out:
            with h5py.File(self.snps_database, 'r') as io5:
                # load seqarray
                seqarr = io5['snps']

                # write dims
                out.write("{} {}\n".format(
                    len(self.data.snames), seqarr.shape[1]))

                # write to disk one row at a time
                # (todo: chunk optimize for this.)
                for idx in range(io5['snps'].shape[0]):
                    seq = seqarr[idx, :].view("S1")
                    out.write(
                        "{}{}".format(
                            self.data.pnames[self.data.snames[idx]],
                            b"".join(seq).decode().upper() + "\n",
                        )
                    )


    def write_snps_map(self):
        "write a map file with linkage information for SNPs file"
        with open(self.data.outfiles.snpsmap, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                # access array of data
                maparr = io5["snpsmap"]

                ## write to map file in chunks of 10000
                for start in range(0, maparr.shape[1], 10000):
                    outchunk = []

                    # grab chunk
                    rdat = maparr[start:start + 10000, :]

                    # get chroms
                    if self.data.isref:
                        revdict = chroms2ints(self.data, 1)
                        for i in rdat:
                            outchunk.append(
                                "{}\t{}:{}\t{}\t{}\n"
                                .format(
                                    i[0], 
                                    revdict[i[3]], 
                                    i[4], 
                                    i[3], 
                                    i[2] + 1,
                                )
                            )
                    else:    
                        # convert to text for writing
                        for i in rdat:
                            outchunk.append(
                                "{}\trad{}_snp{}\t{}\t{}\t{}\n"
                                .format(
                                    i[0], 
                                    i[0] - 1, 
                                    i[1], 
                                    revdict[i[3]], 
                                    i[2] + 1, i[4],
                                )
                            )

                    # write chunk to file
                    out.write("".join(outchunk))
                    outchunk = []
                    

    def write_str(self):
        # write data from snps database, resolve ambiguous bases and numeric.
        with open(self.data.outfiles.str, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                snparr = io5["snps"]

                if self.data.paramsdict["max_alleles_consens"] > 1:
                    for idx, name in enumerate(self.data.pnames):
                        # get row of data
                        snps = snparr[idx, :].view("S1")
                        # expand for ambiguous bases
                        snps = [BTS[i.upper()] for i in snps]
                        # convert to numbers and write row for each resolution
                        sequence = "\t".join([STRDICT[i[0]] for i in snps])
                        out.write(
                            "{}\t\t\t\t\t{}\n"
                            .format(self.data.pnames[name], sequence))
                        sequence = "\t".join([STRDICT[i[1]] for i in snps])                            
                        out.write(
                            "{}\t\t\t\t\t{}\n"
                            .format(self.data.pnames[name], sequence))

                else:
                    for idx, name in enumerate(self.data.pnames):
                        # get row of array data
                        snps = snparr[idx, :].view("S1")
                        # expand for ambiguous bases
                        snps = [BTS[i.upper()] for i in snps]
                        # convert to numbers and write row for each resolution
                        sequence = "\t".join([STRDICT[i[0]] for i in snps])
                        out.write(
                            "{}\t\t\t\t\t{}\n"
                            .format(self.data.pnames[name], sequence))


    def write_gphocs(self):
        "b/c it is similar to .loci we just parse .loci and modify it."
        with open(self.data.outfiles.gphocs, 'w') as out:
            indat = iter(open(self.data.outfiles.loci, 'r'))

            # write nloci header
            out.write("{}\n".format(
                self.data.stats_dfs.s7_loci["sum_coverage"].max()))

            # read in each locus at a time
            idx = 0
            loci = []
            locus = []
            while 1:
                try:
                    line = next(indat)
                except StopIteration:
                    indat.close()
                    break

                # end of locus
                if line.endswith("|\n"):
                    
                    # write stats and locus to string and store
                    nsamp = len(locus)
                    slen = len(locus[0].split()[-1])
                    locstr = ["locus{} {} {}\n".format(idx, nsamp, slen)]
                    loci.append("".join(locstr + locus))

                    # reset locus
                    idx += 1
                    locus = []

                else:
                    locus.append(line)

                if not idx % 10000:
                    out.write("\n".join(loci))
                    loci = []
                    
            # write to file
            if loci:
                out.write("\n".join(loci))


    def write_geno(self):
        with open(self.data.outfiles.geno, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                genos = io5["genos"][:]
                snpgenos = np.zeros(genos.shape[:2], dtype=np.uint8)
                snpgenos.fill(9)
                
                # fill (0, 0)
                snpgenos[np.all(genos == 0, axis=2)] = 2
                
                # fill (0, 1) and (1, 0)
                snpgenos[np.sum(genos, axis=2) == 1] = 1
                
                # fill (1, 1)
                snpgenos[np.all(genos == 1, axis=2)] = 0
                
                # write to file
                np.savetxt(out, snpgenos, delimiter="", fmt="%d")


# -----------------------------------------------------------
# Step7 external functions that are run on engines. These do not require
# The step class object but only the data class object, which avoids
# problems with sending ipyclient (open file) to an engine.
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


# ------------------------------------------------------------
# funcs parallelized on remote engines 
# -------------------------------------------------------------
def write_loci_and_alleles(data):

    # write alleles file
    allel = 'a' in data.paramsdict["output_formats"]

    # gather all loci bits
    locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
    sortbits = sorted(locibits, 
        key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

    # what is the length of the name padding?
    with open(sortbits[0], 'r') as test:
        pad = np.where(np.array(list(test.readline())) == " ")[0].max()

    # write to file while adding counters to the ordered loci
    outloci = open(data.outfiles.loci, 'w')
    if allel:
        outalleles = open(data.outfiles.alleles, 'w')

    for bit in sortbits:
        # store until writing
        lchunk = []
        achunk = []
        idx = 0

        # LOCI ONLY: iterate through chunk files
        if not allel:
            for line in iter(open(bit, 'r')):
                if "|\n" not in line:
                    lchunk.append(line[:pad] + line[pad:].upper())
                else:
                    snpstring, nidxs = line.rsplit("|", 2)[:2]
                    if data.paramsdict["assembly_method"] == 'reference':
                        refpos = nidxs.split(",")[0]
                        lchunk.append(
                            "{}|{}|\n".format(snpstring, refpos))                   
                    else:
                        lchunk.append(
                            "{}|{}|\n".format(snpstring, idx))
                    idx += 1

        # ALLELES: iterate through chunk files to write LOCI AND ALLELES
        else:
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
                    snpstring, nidxs = line.rsplit("|", 2)[:2]
                    asnpstring = "//  " + snpstring[2:]
                    if data.paramsdict["assembly_method"] == 'reference':
                        refpos = nidxs.split(",")[0]
                        lchunk.append(
                            "{}|{}|\n".format(snpstring, refpos))
                        achunk.append(
                            "{}|{}|\n".format(asnpstring, refpos))
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


#VCF can use revdict and snpsmap to build everything, I think. BUT, merged
#consens reads still need to be supported in fill_vcf_depths... 
# somehow (sums?)

def pseudoref2ref(pseudoref, ref):
    """
    Reorder psuedoref (observed bases at snps sites) to have the ref allele
    listed first. On rare occasions when ref is 'N' then 
    """
    # create new empty array
    npseudo = np.zeros(pseudoref.shape, dtype=np.uint8)

    # at all sites where pseudo 0 matches reference, leave it
    matched = np.where(pseudoref[:, 0] == ref)[0]
    npseudo[matched] = pseudoref[matched, :]

    # at other sites, shift order so ref is first
    notmatched = np.where(pseudoref[:, 0] != ref)[0]
    for row in notmatched:
        dat = list(pseudoref[row])

        # skips if ref allele is missing (N)
        try:
            # pop ref and insert it
            new = dat.pop(dat.index(ref[row]))    
            dat.insert(0, new)
            npseudo[row] = dat
        except ValueError:
            npseudo[row] = pseudoref[row]

    return npseudo


def chroms2ints(data, which):
    # if reference-mapped then parse the fai to get index number of chroms
    fai = pd.read_csv(
        data.paramsdict["reference_sequence"] + ".fai",
        names=['scaffold', 'size', 'sumsize', 'a', 'b'],
        sep="\t",
    )
    faidict = {j: i for i, j in enumerate(fai.scaffold)}
    if not which:
        return faidict
    else:
        revdict = {j: i for i, j in faidict.items()}
        return revdict


def fill_seq_array(data, ntaxa, nbases, nloci):
   
    # init/reset hdf5 database
    with h5py.File(data.seqs_database, 'w') as io5:

        # temporary array data sets 
        io5.create_dataset(
            name="phy",
            shape=(ntaxa, nbases), 
            dtype=np.uint8,
        )
        # temporary array data sets 
        io5.create_dataset(
            name="phymap",
            shape=(nloci,),
            dtype=np.int64,
        )

        # gather all loci bits
        locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
        sortbits = sorted(locibits, 
            key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

        # name order for entry in array
        snames = data.snames
        sidxs = {sample: i for (i, sample) in enumerate(snames)}

        # iterate through file
        gstart = 0
        start = end = 0
        maxsize = 100000
        tmploc = {}
        maplist = []
        mapstart = mapend = 0
        locidx = 0

        # array to store until writing
        tmparr = np.zeros((ntaxa, maxsize + 5000), dtype=np.uint8)
        
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
                    loc = (np.array([list(i) for i in tmploc.values()])
                        .astype(bytes).view(np.uint8))
                    
                    # drop the site that are all N or -
                    mask = np.all((loc == 45) | (loc == 78), axis=0)
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
                    
                    # trim right overflow from tmparr (end filled as 0s)
                    trim = np.where(tmparr != 0)[1]
                    if trim.size:
                        trim = trim.max() + 1
                    else:
                        trim = tmparr.shape[1]

                    # fill missing with 78 (N)
                    tmparr[tmparr == 0] = 78
                    
                    # dump tmparr to hdf5
                    io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]                       
                    io5['phymap'][mapstart:locidx] = (
                        np.array(maplist, dtype=np.int64))
                    mapstart = locidx
                    maplist = []
                    
                    # reset
                    tmparr = np.zeros((ntaxa, maxsize + 5000), dtype=np.uint8)
                    gstart += trim
                    start = end = 0
                    
        # trim final chunk tmparr to size
        trim = np.where(tmparr != 0)[1]
        if trim.size:
            trim = trim.max() + 1
        else:
            trim = tmparr.shape[1]

        # fill missing with 78 (N)
        tmparr[tmparr == 0] = 78

        # dump tmparr and maplist to hdf5
        io5['phy'][:, gstart:gstart + trim] = tmparr[:, :trim]           
        mapend = mapstart + len(maplist)
        io5['phymap'][mapstart:mapend] = np.array(maplist, dtype=np.int64)

        # write stats to the output file
        with open(data.stats_files.s7, 'a') as outstats:
            trim = io5["phymap"][locidx - 1]
            missmask = io5["phy"][:trim] == 78
            missmask += io5["phy"][:trim] == 45
            missing = 100 * (missmask.sum() / io5["phy"][:trim].size)
            print("sequence matrix size: ({}, {}), {:.2f}% missing sites."
                .format(
                    len(data.samples), 
                    trim, 
                    max(0, missing),
                ),
                file=outstats,
            )


def fill_snp_array(data, ntaxa, nsnps):

    # get faidict to convert chroms to ints
    faidict = chroms2ints(data, 0)    

    # open new database file handle
    with h5py.File(data.snps_database, 'w') as io5:

        # Database files for storing arrays on disk. 
        # Should optimize for slicing by rows if we run into slow writing, or 
        # it uses too much mem. For now letting h5py to auto-chunking.
        io5.create_dataset(
            name="snps",
            shape=(ntaxa, nsnps),
            dtype=np.uint8,
        )
        # store snp locations:
        # (loc-counter, loc-snp-counter, loc-snp-pos, chrom, chrom-snp-pos)
        io5.create_dataset(
            name="snpsmap",
            shape=(nsnps, 5),
            dtype=np.uint32,
        )
        # store snp locations
        io5.create_dataset(
            name="pseudoref",
            shape=(nsnps, 4),
            dtype=np.uint8,
        )
        # store genotype calls (0/0, 0/1, 0/2, etc.)
        io5.create_dataset(
            name="genos",
            shape=(nsnps, ntaxa, 2),
            dtype=np.uint8,
        )

        # gather all loci bits
        locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
        sortbits = sorted(locibits, 
            key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

        # name order for entry in array
        sidxs = {sample: i for (i, sample) in enumerate(data.snames)}

        # iterate through file
        start = end = 0
        tmploc = {}
        locidx = 1
        snpidx = 1            

        # array to store until writing
        tmparr = np.zeros((ntaxa, nsnps), dtype=np.uint8)
        tmpmap = np.zeros((nsnps, 5), dtype=np.uint32)

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
                    loc = (
                        np.array([list(i) for i in tmploc.values()])
                        .astype(bytes).view(np.uint8)
                        )
                    snps, idxs, _ = line[len(data.snppad):].rsplit("|", 2)
                    snpsmask = np.array(list(snps)) != " "
                    snpsidx = np.where(snpsmask)[0]

                    # select only the SNP sites
                    snpsites = loc[:, snpsmask]

                    # store end position of locus for map
                    end = start + snpsites.shape[1]
                    for idx, name in enumerate(tmploc):
                        tmparr[sidxs[name], start:end] = snpsites[idx, :]
                        
                    # store snpsmap data 1-indexed with chroms info
                    if data.isref:
                        chrom, pos = idxs.split(",")[0].split(":")
                        start = int(pos.split("-")[0])
                        chromidx = faidict[chrom]
                        for isnp in range(snpsites.shape[1]):
                            isnpx = snpsidx[isnp]
                            tmpmap[snpidx - 1] = (
                                locidx, isnp, isnpx, chromidx, isnpx + start,
                            )
                            snpidx += 1
                    # store snpsmap data (snpidx is 1-indexed)
                    else:
                        for isnp in range(snpsites.shape[1]):
                            tmpmap[snpidx - 1] = (locidx, isnp, snpsidx[isnp])
                            snpidx += 1
                    locidx += 1

                    # reset locus
                    start = end
                    tmploc = {}
                        
        # fill missing with 78 (N)
        tmparr[tmparr == 0] = 78

        # dump tmparr and maplist to hdf5       
        io5['snps'][:] = tmparr[:]
        io5['snpsmap'][:] = tmpmap
        del tmparr

        # write stats output
        with open(data.stats_files.s7, 'a') as outstats:
            missmask = io5["snps"][:] == 78
            missmask += io5["snps"][:] == 45
            missing = 100 * (missmask.sum() / io5["snps"][:nsnps].size)
            print("snps matrix size: ({}, {}), {:.2f}% missing sites."
                .format(
                    len(data.samples),
                    nsnps,
                    missing,
                ),
                file=outstats,
            )


        # fill in the reference and geno arrays
        # convert snps to numeric.
        snparr = io5["snps"][:].view("S1")
        snparr = np.char.upper(snparr).view(np.uint8)

        # store pseudo-ref (most common base)
        # with ambiguous bases resolved: (87, 78, 0, 0).
        if data.paramsdict['assembly_method'] != 'reference':
            io5['pseudoref'][:] = reftrick(snparr, GETCONS)
    
        else:
            ref = snparr[data.snames.index('reference')]   
            pseudoref = reftrick(snparr, GETCONS)
            io5['pseudoref'][:] = pseudoref2ref(pseudoref, ref)

        # fill for each taxon
        for sidx in range(ntaxa):
            resos = [DCONS[i] for i in snparr[sidx, :]]

            # pseudoref version
            io5['genos'][:, sidx, :] = get_genos(
                np.array([i[0] for i in resos]), 
                np.array([i[1] for i in resos]),
                io5['pseudoref'][:]
                )


def fill_vcf_depths(data, nsnps, sample):
    "fills an array with SNP depths for VCF file."
    
    # array to store vcfdepths for this taxon
    vcfd = np.zeros((nsnps, 4), dtype=np.uint8)
    
    # load snpsmap with locations of SNPS on trimmed loci
    with h5py.File(data.snps_database, 'r') as io5:
        snpsmap = io5['snpsmap'][:, [0, 2]]   

    # load catgs for this sample
    with h5py.File(sample.files.database, 'r') as io5:
        catgs = io5['catg'][:]
        
    # loci and indel chunks with sidx labels
    locibits = glob.glob(os.path.join(data.tmpdir, "chunk*.loci"))
    sortbits = sorted(locibits, key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))
    indels = glob.glob(os.path.join(data.tmpdir, "chunk*.npy"))
    sindels = sorted(indels, key=lambda x: int(x.rsplit("-", 1)[-1][:-4]))           

    # counters
    locidx = 0
    snpidx = 0
    names = []
    seqs = []

    # iterate through loci to impute indels into catgs
    for locbit, indbit in zip(sortbits, sindels):
        localidx = 0

        # load trim array with left trim values of loci
        trim = np.load(indbit)
        with open(locbit, 'r') as inloci:
            for line in iter(inloci):

                # continue filling locus
                if "|\n" not in line:
                    name, seq = line.split()
                    names.append(name)
                    seqs.append(seq)

                # locus filled, process it.
                else:
                    # get snps for this locus
                    locsnps = snpsmap[snpsmap[:, 0] == locidx + 1]
                        
                    # is this sample in the locus?
                    if sample.name in names:
                        
                        # store catg idx
                        sidxs = [i for i in line.rsplit("|", 2)[1].split(',')]
                        if data.isref:
                            sidxs = [0] + [int(i) for i in sidxs[1:]]
                        else:
                            sidxs = [int(i) for i in sidxs]
                        nidx = names.index(sample.name)
                        sidx = sidxs[nidx]
                        seq = seqs[nidx]

                        # get trim for this locus
                        itrim = trim[localidx]

                        # enter snps into vcfd
                        if locsnps.size:
                            seqarr = np.array(list(seq))
                            inds = np.where(seqarr == "-")[0]
                            for snp in locsnps[:, 1]:
                                # how many indels come before this position?
                                shift = np.sum(inds < snp)

                                # if this site itself is an indel then offset
                                # would confuse it, and it should be zero.
                                try:
                                    if seqarr[snp] == "-":
                                        vcfd[snpidx] = (0, 0, 0, 0)
                                    else:                                                                       
                                        # If IndexError on right end this means
                                        # indels are shifted so that ---- 
                                        # filled the end of alignment
                                        vcfd[snpidx] = (
                                            catgs[sidx, (snp + itrim) - shift]
                                        )
                                except IndexError:
                                    vcfd[snpidx] = (0, 0, 0, 0)
                                snpidx += 1

                    # sample not in this locus still advance SNP counter
                    else:
                        snpidx += locsnps.shape[0]

                    # advance counter and reset dict
                    locidx += 1
                    localidx += 1
                    names = []
                    seqs = []

    # write vcfd to file
    vcfout = os.path.join(data.tmpdir, sample.name + ".depths.hdf5")
    with h5py.File(vcfout, 'w') as io5:
        io5.create_dataset(
            name="depths",
            data=vcfd,
        )


def build_vcf(data, chunksize=1000):
    "pull in data from several hdf5 databases to construct vcf string"

    # dictionary to translate locus numbers to chroms
    if data.isref:
        revdict = chroms2ints(data, 1)

    # pull locus numbers and positions from snps database
    with h5py.File(data.snps_database, 'r') as io5:

        # iterate over chunks
        for chunk in range(0, io5['genos'].shape[0], chunksize):

            # load array chunks
            if data.isref:
                genos = io5['genos'][chunk:chunk + chunksize, :, :]
            else:
                genos = io5['genos'][chunk:chunk + chunksize, 1:, :]

            # if reference then psuedo ref is already ordered with REF first.
            pref = io5['pseudoref'][chunk:chunk + chunksize]
            snpmap = io5['snpsmap'][chunk:chunk + chunksize]

            # get alt genotype calls
            alts = [
                b",".join(i).decode().strip(",")
                for i in pref[:, 1:].view("S1") 
            ]

            # get chrom names
            if data.isref:
                chroms = [revdict[i] for i in snpmap[:, 3]]
            else:
                chroms = ["locus_{}".format(i - 1) for i in snpmap[:, 0]]
            
            # build df label cols
            df_pos = pd.DataFrame({
                '#CHROM': chroms,
                'POS': snpmap[:, -1] + 1,  # 1-indexed
                'ID': ['.'] * genos.shape[0],
                'REF': [i.decode() for i in pref[:, 0].view("S1")],
                'ALT': alts,
                'QUAL': [13] * genos.shape[0],
                'FILTER': ['PASS'] * genos.shape[0],
            })

            # get number of samples for each site
            offset = 0
            if data.isref:
                offset = 1

            nsamps = (
                genos.shape[1] - offset - np.any(genos == 9, axis=2).sum(axis=1)
                )

            # store sum of coverage at each site
            asums = []
            
            # build depth columns for each sample
            df_depth = pd.DataFrame({})
            for sname in data.snames:
                if sname != "reference":
                
                    # build geno strings
                    genostrs = [
                        b"/".join(i).replace(b"9", b".").decode()
                        for i in genos[:, data.snames.index(sname)]
                        .astype(bytes)
                    ]
                    
                    # build depth and depthsum strings
                    dpth = os.path.join(data.tmpdir, sname + ".depths.hdf5")
                    with h5py.File(dpth, 'r') as s5:
                        dpt = s5['depths'][chunk:chunk + chunksize]
                        sums = [sum(i) for i in dpt]
                        strs = [
                            ",".join([str(k) for k in i.tolist()])
                            for i in dpt
                        ]
                        
                        # save concat string to name
                        df_depth[sname] = [
                            "{}:{}:{}".format(i, j, k) for (i, j, k) in 
                            zip(genostrs, sums, strs)]

                        # add sums to global list
                        asums.append(np.array(sums))

            # make final columns
            nsums = sum(asums)
            colinfo = pd.Series(
                name="INFO",
                data=[
                    "NS={};DP={}".format(i, j) for (i, j) in zip(nsamps, nsums)
                ])
            colform = pd.Series(
                name="FORMAT",
                data=["GT:DP:CATG"] * genos.shape[0],
                )
                    
            # debugging           
            arr = pd.concat([df_pos, colinfo, colform, df_depth], axis=1)

            ## PRINTING VCF TO FILE
            ## choose reference string
            if data.paramsdict["reference_sequence"]:
                reference = data.paramsdict["reference_sequence"]
            else:
                reference = "pseudo-reference (most common base at site)"

            header = VCFHEADER.format(
                date=time.strftime("%Y/%m/%d"),
                version=ipyrad.__version__,
                reference=os.path.basename(reference)
                ) 

            with open(data.outfiles.vcf, 'a') as out:
                if chunk == 0:
                    out.write(header)
                    arr.to_csv(out, sep='\t', index=False)
                else:
                    arr.to_csv(out, sep='\t', index=False, header=False)


# ----------------------------------------
# Processor external functions
# ----------------------------------------
def check_minsamp(seq, position, minsamp):
    "used in Processor.get_edges() for trimming edges of - or N sites."
    minsamp = min(minsamp, seq.shape[0])
    if position == 0:       
        mincovs = np.sum((seq != 78) & (seq != 45), axis=0)
        # find leftmost edge with minsamp coverage
        leftmost = np.where(mincovs >= minsamp)[0]
        if leftmost.size:
            return leftmost[0].min()
        # no sites actually have minsamp coverage although there are minsamp
        # rows of data, this can happen when reads only partially overlap. Loc
        # should be excluded for minsamp filter.
        else:
            raise ValueError("no sites above minsamp coverage in edge trim")
    
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
# jitted Processor functions (njit = nopython mode)
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


@njit
def reftrick(iseq, consdict):
    "Returns the most common base at each site in order."

    altrefs = np.zeros((iseq.shape[1], 4), dtype=np.uint8)
    altrefs[:, 1] = 46

    for col in range(iseq.shape[1]):
        ## expand colums with ambigs and remove N-
        fcounts = np.zeros(111, dtype=np.int64)
        counts = np.bincount(iseq[:, col])  #, minlength=90)
        fcounts[:counts.shape[0]] = counts
        ## set N and - to zero, wish numba supported minlen arg
        fcounts[78] = 0
        fcounts[45] = 0
        ## add ambig counts to true bases
        for aidx in range(consdict.shape[0]):
            nbases = fcounts[consdict[aidx, 0]]
            for _ in range(nbases):
                fcounts[consdict[aidx, 1]] += 1
                fcounts[consdict[aidx, 2]] += 1
            fcounts[consdict[aidx, 0]] = 0

        ## now get counts from the modified counts arr
        who = np.argmax(fcounts)
        altrefs[col, 0] = who
        fcounts[who] = 0

        ## if an alt allele fill over the "." placeholder
        who = np.argmax(fcounts)
        if who:
            altrefs[col, 1] = who
            fcounts[who] = 0

            ## if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 2] = who
            fcounts[who] = 0

            ## if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 3] = who

    return altrefs


@njit
def get_genos(f10, f01, pseudoref):
    res = np.zeros((f10.size, 2), dtype=np.uint8)
    
    for i in range(f10.size):
        match = np.where(f10[i] == pseudoref[i])[0]
        if match.size:
            res[i, 0] = match[0]
        else:
            res[i, 0] = 9
        match = np.where(f01[i] == pseudoref[i])[0]
        if match.size:
            res[i, 1] = match[0]
        else:
            res[i, 1] = 9
    return res


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
MISSING_SAMPLE_IN_DB = """
There are samples in this assembly that were not present in step 6. This is 
likely due to branching or merging. The following samples are not in the step6
database:
{}
"""
BADPOP_SAMPLES = """
There are sample names in the populations assignments that are not present in 
this assembly. This is likely due to a typo and should be corrected. The 
following sample names are in the pop assignments but not in this Assembly:
{}
"""
POPULATION_REQUIRED = """\
Warning: Skipping output format '{}'. Requires population assignments.\
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
OUT_SUFFIX = {
    'l': ('.loci',),
    'p': ('.phy', '.phymap',),
    's': ('.snps', '.snpsmap',),
    'n': ('.nex',),
    'k': ('.str',),
    'a': ('.alleles',),
    'g': ('.geno',),
    'G': ('.gphocs',),
    'u': ('.usnps',),
    'v': ('.vcf',),
    't': ('.treemix',),
    'm': ('.migrate',),
    }
STRDICT = {
    'A': '0', 
    'T': '1', 
    'G': '2', 
    'C': '3', 
    'N': '-9', 
    '-': '-9',
}

VCFHEADER = """\
##fileformat=VCFv4.0
##fileDate={date}
##source=ipyrad_v.{version}
##reference={reference}
##phasing=unphased
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
"""

