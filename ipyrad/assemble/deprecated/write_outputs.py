#!/usr/bin/env python

"""
Filter loci and generate output files.
"""

# py2/3 compatibility
from __future__ import print_function
try:
    from builtins import range, bytes
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

# standard lib imports
import os
import glob
import shutil
import pickle
from collections import Counter

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

# third party imports
from loguru import logger
import numpy as np
import pandas as pd
from numba import njit

# internal imports
import ipyrad
from ipyrad.assemble.utils import IPyradError, clustdealer, splitalleles
from ipyrad.assemble.utils import GETCONS, DCONS, chroms2ints
from ipyrad.assemble.utils import AssemblyProgressBar

# helper classes ported to separate files.
from ipyrad.assemble.write_outputs_helpers import ChunkProcessor
from ipyrad.assemble.write_outputs_converter import Converter
from ipyrad.assemble.write_outputs_vcf import FillVCF, build_vcf



class Step7:
    """
    Organization for step7 funcs.
    """
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

        # dimensions filled by collect_stats and used for array building
        self.nloci = 0
        self.nbases = 0
        self.nsnps = 0
        self.ntaxa = 0

        # dict mapping of samples to padded names for loci file aligning.
        self.data.snames = [i.name for i in self.samples]
        self.data.pnames, self.data.snppad = self.get_padded_names()

        # output file formats to produce ('l' is required).
        self.formats = set(['l']).union(
            set(self.data.params.output_formats))


    def run(self):
        """
        All steps to complete step7 assembly
        """
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
        # throttle job to avoid memory errors based on catg size
        # if 'v' in self.formats:
            # self.remote_fill_depths()
            # self.remote_build_vcf()

        # cleanup
        # if os.path.exists(self.data.tmpdir):
            # shutil.rmtree(self.data.tmpdir)


    def print_headers(self):
        """
        print the CLI header
        """
        if self.data._cli:
            self.data._print(
                "\n{}Step 7: Filtering and formatting output files "
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        """
        get subsamples for this assembly. All must have been in step6
        """

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
            samples = [ref] + sorted(
                list(set(self.data.samples.values())), 
                key=lambda x: x.name)
            return samples
        
        samples = sorted(
            list(set(self.data.samples.values())),
            key=lambda x: x.name)
        return samples


    def setup_dirs(self):
        """
        Create temp h5 db for storing filters and depth variants
        """
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
        """
        get nloci and ncpus to chunk and distribute work across processors
        """
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
        """
        Get padded names for print .loci file
        """
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
        """
        Fills self.data.outfiles with names of the files to be produced.
        """
        # always produce a .loci file + whatever they ask for.
        testformats = list(self.formats)
        for outf in testformats:

            # if it requires a pop file and they don't have one then skip
            # and write the warning to the expected file, to prevent an
            # annoying message every time if you don't have a pops file, but
            # still to be transparent about skipping some files. This caused
            # me some real amount of pain, like "why isnt' the treemix file
            # being created, fudkckkk!!!1" And then like 10 minutes later, oh
            # yeah, no pops file, fml. 3/2020 iao.
            if (outf in ("t", "m")) and (not self.data.populations):
                outfile = os.path.join(
                    self.data.dirs.outfiles,
                    self.data.name + OUT_SUFFIX[outf][0],
                )
                with open(outfile, 'w') as out:
                    out.write(POPULATION_REQUIRED.format(outf))

                # remove format from the set
                self.formats.discard(outf)
                continue

            # store handle to data object
            for ending in OUT_SUFFIX[outf]:

                # store 
                self.data.outfiles[ending[1:]] = os.path.join(
                    self.data.dirs.outfiles,
                    self.data.name + ending)           


    def collect_stats(self):
        """
        Collect results from ChunkProcessor and write stats file.
        """
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

        # store dimensions for array building 
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

        # bail out here if no loci were found
        if not self.nloci:
            raise IPyradError("No loci passed filters.")


    def split_clusters(self):
        """
        Splits the step6 clust_database into chunks to be processed
        in parallel by ChunkProcessor to apply filters.
        """
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
                except IndexError as err:
                    msg = "clust_database formatting error in {}".format(chunk)
                    logger.exception(msg)
                    raise IPyradError from err

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
        """
        Calls process_chunk() function in parallel.
        """
        printstr = ("applying filters    ", "s7")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        jobs = glob.glob(os.path.join(self.data.tmpdir, "chunk-*"))
        jobs = sorted(jobs, key=lambda x: int(x.rsplit("-")[-1]))        
        rasyncs = {}
        for jobfile in jobs:
            args = (self.data, self.chunksize, jobfile)
            rasyncs[jobfile] = self.lbview.apply(process_chunk, *args)

        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_build_arrays_and_write_loci(self):
        """
        Calls write_loci_and_alleles(), fill_seq_array() and fill_snp_array().
        """
        # start loci concatenating job on a remote
        printstr = ("building arrays     ", "s7")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        args1 = (self.data, self.ntaxa, self.nbases, self.nloci)
        args2 = (self.data, self.ntaxa, self.nsnps)

        # fill with filtered loci chunks from ChunkProcessor
        rasyncs = {}
        rasyncs[0] = self.lbview.apply(write_loci_and_alleles, self.data)
        rasyncs[1] = self.lbview.apply(fill_seq_array, *args1)
        rasyncs[2] = self.lbview.apply(fill_snp_array, *args2)

        # track progress.
        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_write_outfiles(self):
        """
        Calls convert_outputs() in parallel.
        """
        printstr = ("writing conversions ", "s7")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasyncs = {}
        for outf in self.formats:
            rasyncs[outf] = self.lbview.apply(
                convert_outputs, *(self.data, outf))

        # iterate until all chunks are processed
        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_fill_depths(self):
        """
        Call fill_vcf_depths() in parallel.
        """
        printstr = ("indexing vcf depths ", "s7")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasyncs = {}
        for sample in self.data.samples.values():
            if not sample.name == "reference":
                rasyncs[sample.name] = self.lbview.apply(
                    fill_vcf_depths, *(self.data, self.nsnps, sample))
        # iterate until all chunks are processed
        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_build_vcf(self):
        """
        Calls build_vcf() in parallel. 
        """
        printstr = ("writing vcf output  ", "s7")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasync = self.lbview.apply(build_vcf, self.data)          
        prog.jobs = {0: rasync}
        prog.block()
        prog.check()



# ------------------------------------------------------------
# Classes initialized and run on remote engines.
# ------------------------------------------------------------
def process_chunk(data, chunksize, chunkfile):
    """
    init a ChunkProcessor, run it and collect results.
    """
    # process chunk writes to files and returns proc with features.
    proc = ChunkProcessor(data, chunksize, chunkfile)
    proc.run()

    # check for variants or set max to 0
    try:
        mvar = max([i for i in proc.var if proc.var[i]])
    except ValueError:
        mvar = 0
    try:
        mpis = max([i for i in proc.pis if proc.pis[i]])
    except ValueError:
        mpis = 0

    # shorten dictionaries   
    proc.var = {i: j for (i, j) in proc.var.items() if i <= mvar}
    proc.pis = {i: j for (i, j) in proc.pis.items() if i <= mpis}

    # write process stats to a pickle file for collating later.
    # We have to write stats for each process, even if it returns
    # no loci in order for the filtering stats to make sense.
    # https://github.com/dereneaton/ipyrad/issues/358
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

<<<<<<< HEAD
class Processor(object):
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
            return True
        return False


    def filter_minsamp_pops(self):
        "filter by minsamp or by minsamp x populations"

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
        "max size of internal indels. Denovo vs. Ref, single versus paired."
        # get max indels for read1, read2
        inds = maxind_numba(ublock)        
        if inds > self.maxinds:
            self.filters[self.iloc, 1] = 1
            # return True
        # return False


    def filter_maxvars(self, ublock, snpstring):
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
        nhs = count_maxhet_numba(ublock)
        if nhs > (self.fmaxhet * ublock.shape[0]):
            self.filters[self.iloc, 3] = 1
            # return True
        # return False


    def get_snpsarrs(self, block, exclude_ref=False):
        "count nsnps with option to exclude reference sample from count"
        snpsarr = np.zeros((block.shape[1], 2), dtype=np.bool_)
        return snpcount_numba(block, snpsarr, int(bool(exclude_ref)))


##############################################################

class Edges:
    "Trims edges of overhanging sequences, cutsites, and pair inserts"
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
=======
>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7

###############################################################

def convert_outputs(data, oformat):
    """
    Call the Converter class functions to write formatted output files
    from the HDF5 database inputs.
    """
    try:
        Converter(data).run(oformat)
    except Exception as inst:
        # Allow one file to fail without breaking all step 7
        msg = ("Error creating outfile: {}\n{}\t{}"
            .format(OUT_SUFFIX[oformat], type(inst).__name__, inst))
        logger.exception(msg)
        raise IPyradError(msg)


###############################################################



def fill_vcf_depths(data, nsnps, sample):
    """
    Get catg depths for this sample.
    """
    filler = FillVCF(data, nsnps, sample)
    filler.run()

    # write vcfd to file and cleanup
    vcfout = os.path.join(data.tmpdir, sample.name + ".depths.hdf5")
    with h5py.File(vcfout, 'w') as io5:
        io5.create_dataset(
            name="depths",
            data=filler.vcfd,
            )
    del filler



# ------------------------------------------------------------
# funcs parallelized on remote engines 
# -------------------------------------------------------------
def write_loci_and_alleles(data):
    """
    Write the .loci file from processed loci chunks. Tries to write
    allele files with phased diploid calls if present.
    """
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

                # skip reference lines if excluding
                if data.hackersonly.exclude_reference:
                    if "reference     " in line:
                        continue

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

                # skip reference lines if excluding
                if data.hackersonly.exclude_reference:
                    if "reference     " in line:
                        continue

                if "|\n" not in line:
                    name = line[:pad]
                    seq = line[pad:]
                    lchunk.append(name + seq.upper())

                    all1, all2 = splitalleles(seq)
                    aname, spacer = name.split(" ", 1)
                    # adjust seqnames for proper buffering of the snpstring
                    achunk.append(aname + "_0 " + spacer[2:] + all1)
                    achunk.append(aname + "_1 " + spacer[2:] + all2)
                else:
                    snpstring, nidxs = line.rsplit("|", 2)[:2]
                    # adjust length of snpstring so it lines up for refseq
                    asnpstring = "//" + snpstring[2:]
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
    Called in fill_snps_array.
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
    """
    Fills the HDF5 seqs array from loci chunks and stores phymap.
    This contains the full sequence data for all sites >mincov. 
    """
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
        else:
            phymap.attrs["reference"] = "pseudoref"

        # store names and 
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

        # array to store until writing; TODO: Accomodate large files...
        tmparr = np.zeros((ntaxa, maxsize + 50000), dtype=np.uint8)
        
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

                    # seq ordered into array by snames as int8 (py2/3 checked)
                    loc = np.array([
                        list(bytes(tmploc[i].encode())) for i in snames
                        if i in tmploc
                        ]).astype(np.int8)
                    # loc = (np.array([list(i) for i in tmploc.values()])
                    # .astype(bytes).view(np.uint8))
                    
                    # TODO: check code here for reference excluded...
                    # drop the site that are all N or - (e.g., pair inserts)
                    if (data.isref and data.ispair):
                        mask = np.all(loc[1:, :] == 78, axis=0)
                    else:
                        mask = np.all((loc == 45) | (loc == 78), axis=0)
                    loc = loc[:, np.invert(mask)]
                    
                    # store end position of locus for map
                    end = start + loc.shape[1]

                    # checked for py2/3 (keeping name order straight important)
                    lidx = 0
                    for name in snames:
                        if name in tmploc:
                            sidx = sidxs[name]
                            tmparr[sidx, start:end] = loc[lidx]
                            lidx += 1
                    # tnames = sorted(tmploc.keys())
                    # for idx, name in enumerate(snames):
                        # if name in tmploc
                        # sidx = sidxs[name]
                        # tmparr[sidx, start:end] = loc[idx]
                    # for idx, name in enumerate(tmploc):
                        # tmparr[sidxs[name], start:end] = loc[idx]
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
                    tmparr = np.zeros((ntaxa, maxsize + 50000), dtype=np.uint8)
                    gstart += trim
                    start = end = 0

            # close bit handle
            indata.close()

        if start == 0 and end == 0:
            # The last chunk fell exactly on the maxsize boundary so it has
            # already been trimmed and dumped to the phy. In this case the for
            # loop on the line iterator has fallen through (no more data) so
            # there is no final chunk to trim.
            pass
        else:
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
            phymap[mapstart:locidx, 0] = mapchroms
            phymap[mapstart:locidx, 2] = mapends
            if data.isref:
                phymap[mapstart:locidx, 3] = mappos0
                phymap[mapstart:locidx, 4] = mappos1
        phymap[1:, 1] = phymap[:-1, 2]

        # fill 'scaffold' information for denovo data sets from the data
        if "reference" not in data.params.assembly_method:
            # 1-index the "chromosomes" and store lengths as loclens
            phymap[:, 0] += 1  # <- does not cause problems for outfiles...
            io5["scaffold_names"] = (io5["phymap"][:, 0]).astype("S")
            io5["scaffold_lengths"] = io5["phymap"][:, 2] - io5["phymap"][:, 1]

        # write stats to the output file
        with open(data.stats_files.s7, 'a') as outstats:
            trim = phymap[-1, 2]  # locidx - 1]
            missmask = phy[:trim] == 78
            missmask += phy[:trim] == 45
            missing = 100 * (missmask.sum() / float(phy[:trim].size))
            print("sequence matrix size: ({}, {}), {:.2f}% missing sites."
                .format(
                    len(snames), 
                    trim, 
                    max(0, missing),
                ),
                file=outstats,
            )


def fill_snp_array(data, ntaxa, nsnps):
    """
    Fills the SNPS HDF5 database with variables sites.
    """
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

        # store sample names and snpmap columns names as attributes
        io5["snps"].attrs["names"] = [i.encode() for i in data.snames]
        io5["snpsmap"].attrs["columns"] = [
            b"locus", b"locidx", b"locpos", b"scaf", b"scafpos",
            # b"arrpos",
        ]

        # gather all loci bits
        locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
        sortbits = sorted(
            locibits, 
            key=lambda x: int(x.rsplit("-", 1)[-1][:-5])
        )

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
                    # convert seqs to an np.int8 array, checked py2/3
                    loc = np.array(
                        [list(bytes(tmploc[i].encode())) for i in data.snames 
                         if i in tmploc]
                        ).astype(np.int8)
                    # loc = np.array(
                    # [list(i) for i in tmploc.values()]
                    # ).astype(bytes).view(np.uint8)                        
                    snps, idxs, _ = line[len(data.snppad):].rsplit("|", 2)
                    snpsmask = np.array(list(snps)) != " "
                    snpsidx = np.where(snpsmask)[0]

                    # select only the SNP sites
                    snpsites = loc[:, snpsmask]

                    # store end position of locus for map
                    end = start + snpsites.shape[1]

                    # checked for py2/3 (keeping name order straight important)
                    lidx = 0
                    for name in data.snames:
                        if name in tmploc:
                            sidx = sidxs[name]
                            tmparr[sidx, start:end] = snpsites[lidx, :]
                            lidx += 1
                    # for idx, name in enumerate(tmploc):
                        # tmparr[sidxs[name], start:end] = snpsites[idx, :]

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
            missing = 100 * (missmask.sum() / float(io5["snps"][:nsnps].size))
            print(
                "snps matrix size: ({}, {}), {:.2f}% missing sites."
                .format(
                    len(data.snames),
                    nsnps,
                    missing,
                ),
                file=outstats,
            )

        # fill in the reference and geno arrays
        # convert snps to characters uppered to get most common as pseudoref
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


@njit
def reftrick(iseq, consdict):
    """
    Returns the most common base at each site in order.
    Called from fill_snp_array().
    """
    altrefs = np.zeros((iseq.shape[1], 4), dtype=np.uint8)
    altrefs[:, 1] = 46

    for col in range(iseq.shape[1]):
        # expand colums with ambigs and remove N-
        fcounts = np.zeros(111, dtype=np.int64)
        counts = np.bincount(iseq[:, col])  #, minlength=90)
        fcounts[:counts.shape[0]] = counts
        # set N and - to zero, wish numba supported minlen arg
        fcounts[78] = 0
        fcounts[45] = 0
        # add ambig counts to true bases
        for aidx in range(consdict.shape[0]):
            nbases = fcounts[consdict[aidx, 0]]
            for _ in range(nbases):
                fcounts[consdict[aidx, 1]] += 1
                fcounts[consdict[aidx, 2]] += 1
            fcounts[consdict[aidx, 0]] = 0

        # now get counts from the modified counts arr
        who = np.argmax(fcounts)
        altrefs[col, 0] = who
        fcounts[who] = 0

        # if an alt allele fill over the "." placeholder
        who = np.argmax(fcounts)
        if who:
            altrefs[col, 1] = who
            fcounts[who] = 0

            # if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 2] = who
            fcounts[who] = 0

            # if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 3] = who

    return altrefs


@njit
def get_genos(f10, f01, pseudoref):
    """
    Returns genotype as 0/1/2/3 or 9 for missing.
    Called from fill_snp_array
    """
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
    """
    Returns ... from the indexed reference genome as an array
    """
    reference_file = data.params.reference_sequence
    fai = pd.read_csv(
        reference_file + ".fai",   
        names=['scaffold', 'length', 'sumsize', 'a', 'b'],
        sep="\t",
    )
    return fai[value].values  







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
## The "reference" sample is included if present unless 'exclude_reference=True'
"""
MISSING_SAMPLE_IN_DB = """
There are samples in this assembly that were not present in step 6. This is 
likely due to failed samples retained in the assembly from prior to step 5, or
branching/merging. Either branch and remove these samples, or run them through
step 6. The following samples are not in the step6 database:
{}
Simplest solution is to branch and remove these from the assembly.
"""
BADPOP_SAMPLES = """
There are sample names in the populations assignments that are not present in 
this assembly. This is likely due to a typo and should be corrected. The 
following sample names are in the pop assignments but not in this Assembly:
{}
"""
POPULATION_REQUIRED = """\
Warning: Skipping output format '{}'. Requires population assignments.
You can alternatively create this type of file using ipyrad-analysis 
after assembling your data.
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
    'u': ('.usnps', '.ustr', '.ugeno'),
    'v': ('.vcf',),
    't': ('.treemix',),
    'm': ('.migrate',),
    }



if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    # tdata = ip.load_json("/tmp/test-simpairddrad.json")
    # tdata.params.output_formats = "lpsnkaguvtm"
    # tdata.run("7", auto=True, force=True)
    # logger.info(tdata.stats.T)

    tdata = ip.load_json("/tmp/test-amaranth.json")
    tdata.run("7", auto=True, force=True)
    print(tdata.stats)
    # print(tdata.stats_dfs.s5)

    # self.data.hackersonly.declone_PCR_duplicates:
    # tdata = ip.load_json("/tmp/test-amaranth-denovo.json")
    # tdata.ipcluster['cores'] = 4
    # tdata.run("7", auto=True, force=True)
    # logger.info(tdata.stats.T)
