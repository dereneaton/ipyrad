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
from .utils import IPyradError, clustdealer, splitalleles, chroms2ints
from .utils import BTS, GETCONS, DCONS  # , bcomp


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
        self.data.isref = bool("ref" in self.data.params.assembly_method)
        self.data.ispair = bool("pair" in self.data.params.datatype)

        # returns samples in the order we want them in the outputs
        self.print_headers()
        self.samples = self.get_subsamples()
        self.setup_dirs()
        self.get_chunksize()

        # dict mapping of samples to padded names for loci file aligning.
        self.data.snames = [i.name for i in self.samples]
        self.data.pnames, self.data.snppad = self.get_padded_names()      

        # output file formats to produce ('l' is required).
        self.formats = set(['l']).union(
            set(self.data.params.output_formats))


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
            # throttle job to avoid memory errors based on catg size
            self.remote_fill_depths()
            self.remote_build_vcf()

        # cleanup
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)


    def print_headers(self):
        if self.data._cli:
            self.data._print(
                "\n{}Step 7: Filtering and formatting output files "
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        "get subsamples for this assembly. All must have been in step6"

        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 7")

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

        # output files already exist for this assembly. Raise
        # error unless using the force flag to prevent overwriting. 
        if not self.force:
            _outdir = os.path.join(
                self.data.params.project_dir,
                "{}_outfiles".format(self.data.name),
                )
            _outdir = os.path.realpath(_outdir)
            if os.path.exists(os.path.join(_outdir, 
                "{}.loci".format(self.data.name),
                )):
                raise IPyradError(
        "Step 7 results exist for this Assembly. Use force to overwrite.")

        # if ref init a new sample for reference if including
        if self.data.params.assembly_method == 'reference':
            ref = ipyrad.core.sample.Sample("reference")
            snames = [ref] + sorted(
                list(set(self.data.samples.values())), 
                key=lambda x: x.name)
            return snames
        else:
            snames = sorted(
                list(set(self.data.samples.values())),
                key=lambda x: x.name)
            return snames


    def setup_dirs(self):
        "Create temp h5 db for storing filters and depth variants"

        # reset outfiles paths
        for key in self.data.outfiles:
            self.data.outfiles[key] = ""

        # make new output directory
        self.data.dirs.outfiles = os.path.join(
            self.data.params.project_dir,
            "{}_outfiles".format(self.data.name),
            )
        self.data.dirs.outfiles = os.path.realpath(self.data.dirs.outfiles)
        if os.path.exists(self.data.dirs.outfiles):
            shutil.rmtree(self.data.dirs.outfiles)
        if not os.path.exists(self.data.dirs.outfiles):
            os.makedirs(self.data.dirs.outfiles)

        # stats output handle
        self.data.stats_files.s7 = os.path.abspath(
            os.path.join(
                self.data.dirs.outfiles, 
                "{}_stats.txt".format(self.data.name),
            )
        )

        # make tmpdir directory
        self.data.tmpdir = os.path.join(
            self.data.dirs.outfiles, 
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

        # chunk to approximately 4 chunks per core
        self.ncpus = len(self.ipyclient.ids)
        self.chunksize = sum([
            (self.nraws // (self.ncpus * 4)),
            (self.nraws % (self.ncpus * 4)),
        ])


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


    def store_file_handles(self):
        # always produce a .loci file + whatever they ask for.
        testformats = list(self.formats)
        for outf in testformats:

            # if it requires a pop file and they don't have one then skip
            # and print a warning:
            if (outf in ("t", "m")) and (not self.data.populations):
                #print(POPULATION_REQUIRED.format(outf), file=sys.stderr)

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
        start = self.data.params.min_samples_locus - 1
        sums = pd.Series(
            {i: np.sum(covs[start:i]) for i in lrange},            
            name="sum_coverage", 
            index=lrange)
        self.data.stats_dfs.s7_loci = pd.concat([covs, sums], axis=1)

        # fill pis to match var
        for i in cvar:
            if not cpis.get(i):
                cpis[i] = 0

        ## get SNP distribution       
        sumd = {}
        sump = {}
        for i in range(max(cvar.keys()) + 1):
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

        # trim SNP distribution to exclude unobserved endpoints
        snpmax = np.where(
            np.any(
                self.data.stats_dfs.s7_snps.loc[:, ["var", "pis"]] != 0, axis=1
                )
            )[0]
        if snpmax.size:
            snpmax = snpmax.max()
            self.data.stats_dfs.s7_snps = (
                self.data.stats_dfs.s7_snps.loc[:snpmax])

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
                    done, chunk = clustdealer(pairdealer, self.chunksize)
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
        jobs = sorted(jobs, key=lambda x: int(x.rsplit("-")[-1]))        
        for jobfile in jobs:
            args = (self.data, self.chunksize, jobfile)
            rasyncs[jobfile] = self.lbview.apply(process_chunks, *args)
        
        # iterate until all chunks are processed
        while 1:
            # get and enter results into hdf5 as they come in
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.5)
            if len(ready) == sum(ready):
                self.data._print("")
                break          

        # write stats       
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
                self.data._print("")
                break
        # check for errors
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
                self.data._print("")
                break          

        # write stats
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
                self.data._print("")
                break

        # write stats
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
                self.data._print("")
                break          

        # write stats
        if not rasync.successful():
            raise IPyradError(rasync.get())


# ------------------------------------------------------------
# Classes initialized and run on remote engines.
# ------------------------------------------------------------
def process_chunks(data, chunksize, chunkfile):
    # process chunk writes to files and returns proc with features.
    proc = Processor(data, chunksize, chunkfile)
    proc.run()

    # shorten dictionaries
    iii = max([i for i in proc.var if proc.var[i]])
    proc.var = {i: j for (i, j) in proc.var.items() if i <= iii}
    iii = max([i for i in proc.pis if proc.pis[i]])
    proc.pis = {i: j for (i, j) in proc.pis.items() if i <= iii}

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


##############################################################

class Processor:
    def __init__(self, data, chunksize, chunkfile):
        """
        Takes a chunk of aligned loci and (1) applies filters to it; 
        (2) gets edges, (3) builds snpstring, (4) returns chunk and stats.
        (5) writes 
        """
        # init data
        self.data = data
        self.chunksize = chunksize
        self.chunkfile = chunkfile
        self.isref = self.data.isref
        self.ispair = self.data.ispair

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
        self.io = open(self.chunkfile, 'rb')
        self.loci = enumerate(iter(self.io.read().split(b"//\n//\n")))

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

        self.iloc, lines = next(self.loci)
        lines = lines.decode().strip().split("\n")
        for line in lines:
            if line[0] == ">":
                name, nidx = line[1:].rsplit("_", 1)
                self.names.append(name)
                self.nidxs.append(nidx)
            else:
                self.aseqs.append(list(line))
                self.useqs.append(list(line.upper()))

        # filter to include only samples in this assembly
        mask = [i in self.data.snames for i in self.names]
        self.names = np.array(self.names)[mask].tolist()

        # [ref] store consens read start position as mapped to ref
        self.nidxs = np.array(self.nidxs)[mask].tolist()
        self.useqs = np.array(self.useqs)[mask, :].astype(bytes).view(np.uint8)
        self.aseqs = np.array(self.aseqs)[mask, :].astype(bytes).view(np.uint8)



    def run(self):

        while 1:
            try:
                self.next_locus()
            except StopIteration:
                break

            # apply filters 
            edges = Edges(self.data, self.useqs)
            self.edges[self.iloc] = edges.edges

            # fill filter 0
            self.filter_dups()

            # fill filter 4
            self.filter_minsamp_pops()
            self.filters[self.iloc, 4] += int(edges.bad)

            # bail out of locus now if it is already bad...
            if self.filters[self.iloc].sum():
                continue

            # trim edges, need to use uppered seqs for maxvar & maxshared
            edg = self.edges[self.iloc]
            ublock = self.useqs[:, edg[0]:edg[3]]
            ablock = self.aseqs[:, edg[0]:edg[3]]

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
            snparr = self.get_snpsarrs(ublock)
            self.filter_maxvars(ublock, snparr)

            # apply filters on edge trimmed reads
            self.filter_maxshared(ublock)

            # store stats for the locus that passed filtering
            if not self.filters[self.iloc, :].sum():
                # do sample and locus counters
                for name in self.names:
                    self.scov[name] += 1
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

        # write the chunk to tmpdir
        with open(self.outfile, 'w') as outchunk:
            outchunk.write("\n".join(self.outlist) + "\n")

        # thin edgelist to filtered loci and write to array
        mask = np.invert(self.filters.sum(axis=1).astype(np.bool_))
        np.save(self.outarr, self.edges[mask, 0])

        # close file handle
        self.io.close()


    def to_locus(self, block, snparr, edg):
        "write chunk to a loci string"

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
        if len(set(self.names)) < len(self.names):
            self.filters[self.iloc, 0] = 1
        return False


    def filter_minsamp_pops(self):
        if not self.data.populations:
            if len(self.names) < self.data.params.min_samples_locus:
                self.filters[self.iloc, 4] = 1
                return True
            return False

        else:
            minfilters = []
            for pop in self.data.populations:
                samps = self.data.populations[pop][1]
                minsamp = self.data.populations[pop][0]
                if len(set(samps).intersection(set(self.names))) < minsamp:
                    minfilters.append(pop)
            if any(minfilters):
                self.filters[self.iloc, 4] = 1
                return True
            return False


    def filter_maxindels(self, ublock):
        "max size of internal indels. Denovo vs. Ref, single versus paired."
        # get max indels for read1, read2
        inds = maxind_numba(ublock)        
        if inds > self.maxinds:
            self.filters[self.iloc, 1] = 1
            return True
        return False


    def filter_maxvars(self, ublock, snpstring):
        # mask insert area
        if self.masked is not None:
            if snpstring.sum() > (self.masked.shape[1] * self.fmaxsnps):
                self.filters[self.iloc, 2] = 1
                return True

        # use full locus
        else:
            if snpstring.sum() > (ublock.shape[1] * self.fmaxsnps):
                self.filters[self.iloc, 2] = 1
                return True
        return False


    def filter_maxshared(self, ublock):
        nhs = count_maxhet_numba(ublock)
        if nhs > (self.fmaxhet * ublock.shape[0]):
            self.filters[self.iloc, 3] = 1
            return True
        return False


    def get_snpsarrs(self, block):
        snpsarr = np.zeros((block.shape[1], 2), dtype=np.bool_)
        return snpcount_numba(block, snpsarr)


##############################################################

class Edges:
    "Trims edges of overhanging sequences, cutsites, and pair inserts"
    def __init__(self, data, seqs):
        self.data = data
        self.seqs = seqs
        self.trimseq = None
        self.bad = False
        self.edges = np.array([0, 0, 0, self.seqs.shape[1]])
        self.trims = np.array([0, 0, 0, 0])  # self.seqs.shape[1]])
        self.minlen = self.data.params.filter_min_trim_len

        # get edges of good locus
        self.trim_for_coverage()
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


    def trim_for_coverage(self):
        "trim edges to where data is not N or -"

        # at least this much cov at each site
        minsamp = min(4, self.seqs.shape[0])

        # how much cov at each site?
        mincovs = np.sum((self.seqs != 78) & (self.seqs != 45), axis=0)

        # locus left trim
        self.edges[0] = locus_left_trim(self.seqs, minsamp, mincovs)
        self.edges[3] = locus_right_trim(self.seqs, minsamp, mincovs)
        if self.edges[3] <= self.edges[0]:
            self.bad = True

        # find insert region for paired data to mask it
        self.edges[1] = 0
        self.edges[2] = 0



    def trim_overhangs(self):
        "fuzzy match to trim the restriction_overhangs from r1 and r2"

        # trim left side for overhang
        for cutter in self.data.params.restriction_overhang:
            cutter = np.array(list(cutter.encode()))
            slx = slice(0, cutter.shape[0])
            matching = self.trimseq[:, slx] == cutter
            mask = np.where(
                (self.trimseq[:, slx] != 78) & (self.trimseq[:, slx] != 45))
            matchset = matching[mask]
            if matchset.sum() / matchset.size >= 0.75:
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
                if matchset.sum() / matchset.size >= 0.75:
                    self.trims[3] = len(cutter)


    def trim_param_trim_loci(self):
        "user entered hard trims"
        self.trims[0] = max([self.trims[0], self.data.params.trim_loci[0]])
        self.trims[1] = (self.trims[1] - self.data.params.trim_loci[0]
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


# jit
def locus_left_trim(seqs, minsamp, mincovs):
    leftmost = np.where(mincovs >= minsamp)[0]
    if leftmost.size:
        return leftmost.min()
    return 0

# jit
def locus_right_trim(seqs, minsamp, mincovs):
    rightmost = np.where(mincovs >= minsamp)[0]
    if rightmost.size:
        return rightmost.max() + 1
    return 0

###############################################################

def convert_outputs(data, oformat):
    Converter(data).run(oformat)


###############################################################

class Converter:
    "functions for converting hdf5 arrays into output files"
    def __init__(self, data):
        self.data = data
        self.output_formats = self.data.params.output_formats
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
                arrsize = io5['phymap'][-1, 2]

                # if reference then inserts are not trimmed from phy
                # ...

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
                # load seqarray (this could be chunked, this can be >50Gb)
                seqarr = io5['phy'][:]
                arrsize = io5['phymap'][-1, 2]

                ## write nexus seq header
                out.write(NEXHEADER.format(seqarr.shape[0], arrsize))

                ## grab a big block of data
                chunksize = 100000  # this should be a multiple of 100
                for bidx in range(0, arrsize, chunksize):
                    bigblock = seqarr[:, bidx:bidx + chunksize]
                    lend = int(arrsize - bidx)

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
                maparr = io5["phymap"][:, 2]
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
        counter = 1
        with open(self.data.outfiles.snpsmap, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                # access array of data
                maparr = io5["snpsmap"]

                ## write to map file in chunks of 10000
                for start in range(0, maparr.shape[0], 10000):
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
                                    # 1-index to 0-index fix (1/6/19)
                                    revdict[i[3] - 1], i[4],
                                    i[2] + 1,
                                    counter,
                                    #i[4], 
                                )
                            )
                            counter += 1
                    else:    
                        # convert to text for writing
                        for i in rdat:
                            outchunk.append(
                                "{}\tloc{}_snp{}_pos{}\t{}\t{}\n"
                                .format(
                                    i[0], 
                                    i[0] - 1, i[4] - 1, i[2],
                                    i[2] + 1, 
                                    counter,
                                    #i[4],
                                )
                            )
                            counter += 1

                    # write chunk to file
                    out.write("".join(outchunk))
                    outchunk = []
                    

    def write_str(self):
        # write data from snps database, resolve ambiguous bases and numeric.
        with open(self.data.outfiles.str, 'w') as out:
            with h5py.File(self.data.snps_database, 'r') as io5:
                snparr = io5["snps"]

                if self.data.params.max_alleles_consens > 1:
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


# ------------------------------------------------------------
# funcs parallelized on remote engines 
# -------------------------------------------------------------
def write_loci_and_alleles(data):

    # get faidict to convert chroms to ints
    if data.isref:
        faidict = chroms2ints(data, True)

    # write alleles file
    allel = 'a' in data.params.output_formats

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

    idx = 0
    for bit in sortbits:
        # store until writing
        lchunk = []
        achunk = []

        # LOCI ONLY: iterate through chunk files
        if not allel:
            indata = open(bit, 'r')
            for line in iter(indata):

                # write name, seq pairs
                if "|\n" not in line:
                    lchunk.append(line[:pad] + line[pad:].upper())
                
                # write snpstring and info
                else:
                    snpstring, nidxs = line.rsplit("|", 2)[:2]
                    if data.params.assembly_method == 'reference':
                        refpos = nidxs.split(",")[0]

                        # translate refpos chrom idx (1-indexed) to chrom name
                        cid, rid = refpos.split(":")
                        cid = faidict[int(cid) - 1]
                        lchunk.append(
                            "{}|{}:{}:{}|\n".format(snpstring, idx, cid, rid))
                    else:
                        lchunk.append(
                            "{}|{}|\n".format(snpstring, idx))
                    idx += 1
            # close bit handle            
            indata.close()

        # ALLELES: iterate through chunk files to write LOCI AND ALLELES
        else:
            indata = open(bit, 'r')
            for line in iter(indata):
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
                    if data.params.assembly_method == 'reference':
                        refpos = nidxs.split(",")[0]

                        # translate refpos chrom idx (1-indexed) to chrom name
                        cid, rid = refpos.split(":")
                        cid = faidict[int(cid) - 1]

                        lchunk.append(
                            "{}|{}:{}:{}|\n".format(snpstring, idx, cid, rid))
                        achunk.append(
                            "{}|{}:{}:{}|\n".format(asnpstring, idx, cid, rid))
                    else:
                        lchunk.append(
                            "{}|{}|\n".format(line.rsplit("|", 2)[0], idx))
                        achunk.append(
                            "{}|{}|\n".format(line.rsplit("|", 2)[0], idx))
                    idx += 1
            indata.close()
            outalleles.write("".join(achunk))
        outloci.write("".join(lchunk))
    outloci.close()
    if allel:
        outalleles.close()


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


def fill_seq_array(data, ntaxa, nbases, nloci):
   
    # init/reset hdf5 database
    with h5py.File(data.seqs_database, 'w') as io5:

        # temporary array data sets 
        phy = io5.create_dataset(
            name="phy",
            shape=(ntaxa, nbases), 
            dtype=np.uint8,
        )
        # temporary array data sets 
        phymap = io5.create_dataset(
            name="phymap",
            shape=(nloci, 5),
            dtype=np.uint64,
        )

        # store attrs of the reference genome to the phymap
        if data.params.assembly_method == 'reference':
            io5["scaffold_lengths"] = get_fai_values(data, "length")
            io5["scaffold_names"] = get_fai_values(data, "scaffold").astype("S")
            phymap.attrs["reference"] = data.params.reference_sequence
            phymap.attrs["phynames"] = [i.encode() for i in data.pnames]
            phymap.attrs["columns"] = [
                b"chroms", b"phy0", b"phy1", b"pos0", b"pos1", 
            ]

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
        mapends = []
        mapchroms = []
        mappos0 = []
        mappos1 = []
        mapstart = mapend = 0
        locidx = 0

        # array to store until writing
        tmparr = np.zeros((ntaxa, maxsize + 5000), dtype=np.uint8)
        
        # iterate over chunkfiles
        for bit in sortbits:
            # iterate lines of file until locus endings
            indata = open(bit, 'r')
            for line in iter(indata):
                
                # still filling locus until |\n
                if "|\n" not in line:

                    # if empty skip
                    try:
                        name, seq = line.split()
                        tmploc[name] = seq
                    except ValueError:
                        continue

                # locus is full, dump it
                else:
                    # convert seqs to an array
                    locidx += 1

                    # parse chrom:pos-pos                   
                    if data.isref:
                        lineend = line.split("|")[1]
                        chrom = int(lineend.split(":")[0])
                        pos0, pos1 = 0, 0                    
                        pos0, pos1 = (
                            int(i) for i in lineend
                            .split(":")[1]
                            .split(",")[0]
                            .split("-")
                        )

                    # seq into an array
                    loc = (np.array([list(i) for i in tmploc.values()])
                        .astype(bytes).view(np.uint8))
                    
                    # drop the site that are all N or - (e.g., pair inserts)
                    if (data.isref and data.ispair):
                        mask = np.all(loc[1:, :] == 78, axis=0)
                    else:
                        mask = np.all((loc == 45) | (loc == 78), axis=0)
                    loc = loc[:, np.invert(mask)]
                    
                    # store end position of locus for map
                    end = start + loc.shape[1]
                    for idx, name in enumerate(tmploc):
                        tmparr[sidxs[name], start:end] = loc[idx]
                    mapends.append(gstart + end)

                    if data.isref:
                        mapchroms.append(chrom)
                        mappos0.append(pos0)
                        mappos1.append(pos1)
                    else:
                        mapchroms.append(locidx - 1)
                    
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
                    phy[:, gstart:gstart + trim] = tmparr[:, :trim]                       
                    phymap[mapstart:locidx, 0] = mapchroms
                    phymap[mapstart:locidx, 2] = mapends
                    if data.isref:
                        phymap[mapstart:locidx, 3] = mappos0
                        phymap[mapstart:locidx, 4] = mappos1                       
                    mapstart = locidx
                    mapends = []
                    mapchroms = []
                    mappos0 = []
                    mappos1 = []
                    
                    # reset
                    tmparr = np.zeros((ntaxa, maxsize + 5000), dtype=np.uint8)
                    gstart += trim
                    start = end = 0

            # close bit handle
            indata.close()
                    
        # trim final chunk tmparr to size
        trim = np.where(tmparr != 0)[1]
        if trim.size:
            trim = trim.max() + 1
        else:
            trim = tmparr.shape[1]

        # fill missing with 78 (N)
        tmparr[tmparr == 0] = 78

        # dump tmparr and maplist to hdf5. Because we dropped sites that are 
        # all N or - the length of phy data can be less than nbases and so 
        # there can be 000 at the end of phy. This is ok, we trim it when
        # writing phylip file, but good to be aware it's there for other things
        phy[:, gstart:gstart + trim] = tmparr[:, :trim]           
        mapend = mapstart + len(mapends) + 1
        phymap[mapstart:mapend, 0] = mapchroms
        phymap[mapstart:mapend, 2] = mapends
        if data.isref:
            phymap[mapstart:mapend, 3] = mappos0
            phymap[mapstart:mapend, 4] = mappos1
        phymap[1:, 1] = phymap[:-1, 2]

        # write stats to the output file
        with open(data.stats_files.s7, 'a') as outstats:
            trim = phymap[-1, 2]  # locidx - 1]
            missmask = phy[:trim] == 78
            missmask += phy[:trim] == 45
            missing = 100 * (missmask.sum() / phy[:trim].size)
            print("sequence matrix size: ({}, {}), {:.2f}% missing sites."
                .format(
                    len(data.samples), 
                    trim, 
                    max(0, missing),
                ),
                file=outstats,
            )


def fill_snp_array(data, ntaxa, nsnps):

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
            indata = open(bit, 'r')
            for line in iter(indata):
                
                # while still filling locus until |\n store name,seq in dict
                if "|\n" not in line:
                    try:
                        name, seq = line.split()
                        tmploc[name] = seq
                    except ValueError:
                        continue

                # locus is full, dump it
                else:
                    # convert seqs to an array
                    loc = np.array(
                        [list(i) for i in tmploc.values()]
                    ).astype(bytes).view(np.uint8)                        
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
                        #chromidx = faidict[chrom]
                        chromidx = int(chrom)
                        for isnp in range(snpsites.shape[1]):
                            isnpx = snpsidx[isnp]
                            tmpmap[snpidx - 1] = (
                                locidx, isnp, isnpx, chromidx, isnpx + start,
                            )
                            snpidx += 1

                    # store snpsmap data (snpidx is 1-indexed)
                    else:
                        for isnp in range(snpsites.shape[1]):
                            tmpmap[snpidx - 1] = (
                                locidx, isnp, snpsidx[isnp], 0, snpidx,
                            )
                            snpidx += 1
                    locidx += 1

                    # reset locus
                    start = end
                    tmploc = {}

            # close file handle
            indata.close()

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
        if data.params.assembly_method != 'reference':
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


###############################################################

class VCF_filler:
    """
    Incorporate indels and trim amounts when grabbing depths from CATG arrays
    (depth arrays from step 5). Indels are only releveant to denovo data.
    """
    def __init__(self, data, nsnps, sample):

        # input locus bits
        self.locbits = glob.glob(os.path.join(data.tmpdir, "chunk*.loci"))
        self.locbits = sorted(
            self.locbits, key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))
        self.loclines = None

        # input arrays of indels arrays
        self.indbits = glob.glob(os.path.join(data.tmpdir, "chunk*.indels*"))
        if not self.indbits:
            self.indbits = [None] * len(self.locbits)

        # input trim arrays
        self.trimbits = glob.glob(os.path.join(data.tmpdir, "chunk*.npy"))
        self.trimbits = sorted(
            self.trimbits, key=lambda x: int(x.rsplit("-", 1)[-1][:-4]))

        # array to store vcfdepths for this taxon
        self.vcfd = np.zeros((nsnps, 4), dtype=np.uint32)

        # the sample for this comp
        self.sname = sample.name
        self.isref = bool(data.isref)
        
        # snpsmap has locations of SNPs on trimmed loci, e.g., 
        # no SNPs are on loc 1 and 2, first is on 3 at post-trim pos 11
        # [    3     0    11     1 41935]
        # [    4     0    57     1 56150]
        with h5py.File(data.snps_database, 'r') as io5:
            self.snpsmap = io5['snpsmap'][:, [0, 2]]   

        # TODO: scaffs should be ordered (right?) so no need to load it all!
        # All catgs for this sample (this could be done more mem efficient...)
        with h5py.File(sample.files.database, 'r') as io5:
            self.catgs = io5['catg'][:]
            self.maxlen = self.catgs.shape[1]

        # Sample-level counters
        self.locidx = 0
        self.snpidx = 0


    def run(self):
        "loops over chunked files streaming through all loci for this sample"
        for idx in range(len(self.locbits)):
            self.localidx = 0
            self.locfill(idx)


    def locfill(self, idx):
        "iterates over loci in chunkfile to get and enter catgs for snps"
        # load the arrays for this bit
        edges = np.load(self.trimbits[idx])
        inds = self.indbits[idx]
        if inds:
            inds = np.load(inds)

        # iterate over the chunk of trimmed loci
        self.loclines = iter(open(self.locbits[idx], 'r'))
        while 1:

            # yield increments locidx by 1
            try:
                self.yield_loc()
            except StopIteration:
                break

            # get snps for this locus (1-indexed locus idxs)
            self.locsnps = self.snpsmap[self.snpsmap[:, 0] == self.locidx]

            # get global trim for this locus (0-indexed edge arr)
            self.gtrim = edges[self.localidx - 1]

            # if SNPs and data for this sample enter catgs
            if (self.locsnps.size) and (self.sname in self.names):
                if self.isref:
                    self.ref_enter_catgs()
                else:
                    self.denovo_enter_catgs()
            else:
                # advance SNP counter even though this sample wasn't in SNP
                self.snpidx += self.locsnps.shape[0]


    def ref_enter_catgs(self):

        # map SNP position to pre-trim locus position
        nidx = self.names.index(self.sname)
        sidx = self.sidxs[nidx]
        tups = [[int(j) for j in i.split(":")] for i in sidx.split("-")]

        # SNP is in samples, so get and store catg data for locidx
        # [0] post-trim chrom:start-end of locus
        # [1:] how far ahead of start does this sample start
        # FOR DEBUGGING 
        # seq = seqs[nidx]
        # seqarr = np.array(list(seq))

        # enter each SNP 
        for snp in self.locsnps[:, 1]:
            # in case multiple consens were merged in step 6 of this sample
            for tup in tups:
                cidx, coffset = tup
                pos = snp + (self.gtrim - coffset)
                if (pos >= 0) & (pos < self.maxlen):
                    self.vcfd[self.snpidx] += self.catgs[cidx, pos]
            self.snpidx += 1


    def denovo_enter_catgs(self):
        """
        Grab catg depths for each SNP position -- needs to take into account
        trim from left end, and impution of indels.
        """
        nidx = self.names.index(self.sname)
        sidx = self.sidxs[nidx]
        tups = [[int(j) for j in i.split("-")] for i in sidx.split(":")]

        # SNP is in samples, so get and store catg data for locidx
        # [0] post-trim chrom:start-end of locus
        # [1:] how far ahead of start does this sample start
        # FOR DEBUGGING 
        seq = self.seqs[nidx]

        # enter each SNP 
        for snp in self.locsnps[:, 1]:
            # indels before this SNP
            ishift = seq[:snp].count("-")

            # in case multiple consens were merged in step 6 of this sample
            for tup in tups:
                cidx, coffset = tup
                # pos = snp + (self.gtrim - coffset) - ishift
                pos = snp + coffset - ishift                
                if (pos >= 0) & (pos < self.maxlen):
                    self.vcfd[self.snpidx] += self.catgs[cidx, pos]
            self.snpidx += 1


    def yield_loc(self):
        self.names = []
        self.seqs = []
        while 1:
            line = next(self.loclines)
            if "|\n" not in line:
                # skip if .loci chunk is empty
                try:
                    name, seq = line.split()
                    self.names.append(name)
                    self.seqs.append(seq)
                except ValueError:
                    continue
            else:
                self.locidx += 1
                self.localidx += 1                
                self.sidxs = [i for i in line.rsplit("|", 2)[1].split(',')]
                break



def fill_vcf_depths(data, nsnps, sample):
    "get catg depths for this sample."
    filler = VCF_filler(data, nsnps, sample)
    filler.run()

    # write vcfd to file and cleanup
    vcfout = os.path.join(data.tmpdir, sample.name + ".depths.hdf5")
    with h5py.File(vcfout, 'w') as io5:
        io5.create_dataset(
            name="depths",
            data=filler.vcfd,
            )
    del filler


def build_vcf(data, chunksize=1000):
    # removed at init of Step function anyway.
    if os.path.exists(data.outfiles.vcf):
        os.remove(data.outfiles.vcf)
        
    # dictionary to translate locus numbers to chroms
    if data.isref:
        revdict = chroms2ints(data, True)

    # pull locus numbers and positions from snps database
    with h5py.File(data.snps_database, 'r') as io5:

        # iterate over chunks
        for chunk in range(0, io5['genos'].shape[0], chunksize):

            # if reference then psuedo ref is already ordered with REF first.
            pref = io5['pseudoref'][chunk:chunk + chunksize]
            snpmap = io5['snpsmap'][chunk:chunk + chunksize]

            # load array chunks
            if data.isref:
                genos = io5['genos'][chunk:chunk + chunksize, 1:, :]
                snames = data.snames[1:]

                # 1-indexed to 0-indexed (1/9/2019)
                chroms = [revdict[i - 1] for i in snpmap[:, 3]]
                ids = [
                    "loc{}_pos{}".format(i - 1, j) for (i, j) 
                    in snpmap[:, [0, 2]]
                ]

                # reference based positions: pos on scaffold: 4, yes. tested.
                pos = snpmap[:, 4]
                #offset = 1
                
            else:
                genos = io5['genos'][chunk:chunk + chunksize, :, :]
                snames = data.snames
                chroms = ["RAD_{}".format(i - 1) for i in snpmap[:, 0]]
                ids = [
                    "loc{}_pos{}".format(i - 1, j) for (i, j) 
                    in snpmap[:, [0, 2]]
                ]
                # denovo based positions: pos on locus. tested. works. right.
                pos = snpmap[:, 2]
                # offset = 0

            # get alt genotype calls
            alts = [
                b",".join(i).decode().strip(",")
                for i in pref[:, 1:].view("S1") 
            ]

            # build df label cols
            df_pos = pd.DataFrame({
                '#CHROM': chroms,
                'POS': pos,            # 1-indexed
                'ID': ids,             # 0-indexed
                'REF': [i.decode() for i in pref[:, 0].view("S1")],
                'ALT': alts,
                'QUAL': [13] * genos.shape[0],
                'FILTER': ['PASS'] * genos.shape[0],
            })

            # get sample coverage at each site
            nsamps = (
                genos.shape[1] - np.any(genos == 9, axis=2).sum(axis=1)
            )

            # store sum of coverage at each site
            asums = []

            # build depth columns for each sample
            df_depth = pd.DataFrame({})
            for sname in snames:

                # build geno strings
                genostrs = [
                    b"/".join(i).replace(b"9", b".").decode()
                    for i in genos[:, snames.index(sname)]
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

            # concat and order columns correctly
            infocols = pd.concat([df_pos, colinfo, colform], axis=1)
            infocols = infocols[["#CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT"]]
            arr = pd.concat([infocols, df_depth], axis=1)

            # debugging                       
            #print(arr.head())
            ## PRINTING VCF TO FILE
            ## choose reference string
            if data.isref:
                reference = data.params.reference_sequence
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


# -------------------------------------------------------------
# jitted Processor functions (njit = nopython mode)
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
def count_maxhet_numba(block):
    counts = np.zeros(block.shape[1], dtype=np.int16)
    for fidx in range(block.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(block[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts.max()


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


def get_fai_values(data, value):
    reference_file = data.params.reference_sequence
    fai = pd.read_csv(
        reference_file + ".fai",   
        names=['scaffold', 'length', 'sumsize', 'a', 'b'],
        sep="\t",
    )
    return fai[value].values  



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

