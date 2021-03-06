#!/usr/bin/env python

"""
Call consensus base calls on paired or single-end stacks/contigs
"""

# py2/3 compatible
from __future__ import print_function
try:
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

import os
import time
import gzip
import glob
import shutil
import warnings
import subprocess as sps

from loguru import logger
import numpy as np
import pandas as pd

import ipyrad as ip
from ipyrad.assemble.jointestimate import recal_hidepth
from ipyrad.assemble.utils import IPyradError, clustdealer
from ipyrad.assemble.utils import AssemblyProgressBar
from ipyrad.assemble.consens_process import Processor

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class Step5:
    """
    Organized Step 5 functions for all datatype and methods
    """
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.print_headers()
        self.samples = self.get_subsamples()
        self.isref = bool("reference" in data.params.assembly_method)
        self.ipyclient = ipyclient
        self.lbview = ipyclient.load_balanced_view()
        self.setup_dirs()


    def print_headers(self):
        """
        print headers for the CLI
        """
        if self.data._cli:
            self.data._print(
                "\n{}Step 5: Consensus base/allele calling "
                .format(self.data._spacer)
            )



    def get_subsamples(self):
        """
        Apply state, ncluster, and force filters to select samples
        """
        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 5")

        # filter samples by state
        state3 = self.data.stats.index[self.data.stats.state < 4]
        state4 = self.data.stats.index[self.data.stats.state == 4]
        state5 = self.data.stats.index[self.data.stats.state > 4]

        # tell user which samples are not ready for step5
        if state3.any():
            print("skipping samples not in state==4:\n{}"
                  .format(state3.tolist()))

        if self.force:
            # run all samples above state 3
            subs = self.data.stats.index[self.data.stats.state > 3]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 5
            if state5.any():
                print("skipping samples already finished step 5:\n{}"
                      .format(state5.tolist()))

            # run all samples in state 4
            subsamples = [self.data.samples[i] for i in state4]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.clusters_hidepth:
                checked_samples.append(sample)
            else:
                print("skipping {}; no clusters found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 5")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.clusters_hidepth,
            reverse=True,
        )

        # if sample is already done skip
        if "hetero_est" not in self.data.stats:
            for sample in checked_samples:
                sample.stats.hetero_est = 0.001
                sample.stats.error_est = 0.0001

        if self.data._cli:
            print(u"  Mean error  [{:.5f} sd={:.5f}]".format(
                self.data.stats.error_est.mean(),
                self.data.stats.error_est.std(),
            ))
            print(u"  Mean hetero [{:.5f} sd={:.5f}]".format(
                self.data.stats.hetero_est.mean(),
                self.data.stats.hetero_est.std(),
            ))
        return checked_samples



    def setup_dirs(self):
        """
        setup directories, remove old tmp files
        """
        # final results dir
        self.data.dirs.consens = os.path.join(
            self.data.dirs.project, 
            "{}_consens".format(self.data.name))
        if not os.path.exists(self.data.dirs.consens):
            os.mkdir(self.data.dirs.consens)

        # tmpfile dir (zap it if it exists)
        self.data.tmpdir = os.path.join(
            self.data.dirs.project, 
            "{}-tmpdir".format(self.data.name))
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)

        # assign output file handles for s6
        for sample in self.samples:
            if not self.isref:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.gz".format(sample.name))
            else:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.bam".format(sample.name))

        # zap existing consens files for selected samples
        for sample in self.samples:
            if os.path.exists(sample.files.consens):
                os.remove(sample.files.consens)
            dfile = sample.files.consens.rsplit(".", 2)[0]
            sample.files.database = dfile + ".catg.hdf5"
            if os.path.exists(dfile):
                os.remove(dfile)

        # set up parallel client: allow user to throttle cpus
        self.lbview = self.ipyclient.load_balanced_view()
        if self.data.ipcluster["cores"]:
            self.ncpus = self.data.ipcluster["cores"]
        else:
            self.ncpus = len(self.ipyclient.ids)



    def run(self):
        """
        Run the main functions on the parallel client
        """
        self.remote_calculate_depths()
        self.remote_make_chunks()
        statsdicts = self.remote_process_chunks()
        self.remote_concatenate_chunks()
        self.data_store(statsdicts)
        self.data.save()



    def remote_calculate_depths(self):
        """
        Checks whether mindepth has changed and calc nclusters and maxlen
        """
        # send jobs to be processed on engines
        start = time.time()
        printstr = ("calculating depths  ", "s5")
        prog = AssemblyProgressBar({}, start, printstr, self.data)
        prog.update()

        jobs = {}
        maxlens = []
        for sample in self.samples:
            rasync = self.lbview.apply(recal_hidepth, *(self.data, sample))
            jobs[sample.name] = rasync

        # block until finished
        prog.jobs = jobs
        prog.block()
        prog.check()

        # check for failures and collect results
        for sample in self.samples:
            hidepth, maxlen, _, _ = jobs[sample.name].get()
            sample.stats["clusters_hidepth"] = hidepth
            sample.stats_dfs.s3["clusters_hidepth"] = hidepth
            maxlens.append(maxlen)
        
        # update hackersdict with max fragement length
        self.data.hackersonly.max_fragment_length = max(maxlens)
        logger.debug("max_fragment_length set to {}".format(max(maxlens)))



    def remote_make_chunks(self):
        """
        split clusters into chunks for parallel processing
        """
        printstr = ("chunking clusters   ", "s5")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        # send off samples to be chunked
        jobs = {}
        for sample in self.samples:
            args = (self.data, sample, len(self.ipyclient))
            rasync = self.lbview.apply(make_chunks, *args)
            jobs[sample.name] = rasync

        # block until finished
        prog.jobs = jobs
        prog.block()
        prog.check()



    def remote_process_chunks(self):
        """
        Process the cluster chunks into arrays and consens or bam files
        """
        # send chunks to be processed
        printstr = ("consens calling     ", "s5")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        # submit jobs (10 per sample === can be hundreds of jobs...)
        jobs = {sample.name: [] for sample in self.samples}
        for sample in self.samples:
            chunks = glob.glob(os.path.join(
                self.data.tmpdir,
                "{}.chunk.*".format(sample.name)))
            chunks.sort(key=lambda x: int(x.split('.')[-1]))

            # submit jobs
            for chunk in chunks:
                args = (self.data, sample, chunk, self.isref)
                rasync = self.lbview.apply(process_chunks, *args)
                jobs[sample.name].append(rasync)
               
        # track progress - just wait for all to finish before concat'ing
        allsyncs = chain(*[jobs[i] for i in jobs])
        prog.jobs = dict(enumerate(allsyncs))
        prog.update()
        prog.block()
        prog.check()

        # collect all results for a sample and store stats 
        statsdicts = {}
        for sample in self.samples:
            statsdicts[sample.name] = [i.get() for i in jobs[sample.name]]
        return statsdicts



    def remote_concatenate_chunks(self):
        """
        Concatenate chunks and relabel for joined chunks. This spends
        most of its time storing CATG data that will probably not be used,
        but is important for saving SNP depths info.
        """
        # concatenate and store catgs
        printstr = ("indexing alleles    ", "s5")
        prog = AssemblyProgressBar({1:1}, None, printstr, self.data)
        prog.update()

        # concat catgs for each sample
        asyncs1 = {}
        for sample in self.samples:
            args = (self.data, sample, self.isref)
            asyncs1[sample.name] = self.lbview.apply(concat_catgs, *args)

        # select the next job
        if self.isref:
            concat_job = concat_reference_consens
        else:
            concat_job = concat_denovo_consens

        # collect all results for a sample and store stats 
        asyncs2 = {}
        for sample in self.samples:
            args = (self.data, sample)
            rasync = self.lbview.apply(concat_job, *args)
            asyncs2[sample.name] = rasync
            
        # track progress of stats storage
        alljobs = list(asyncs1.values()) + list(asyncs2.values())
        prog.jobs = dict(enumerate(alljobs))
        prog.update()
        prog.block()
        prog.check()



    def data_store(self, statsdicts):
        """
        Store assembly object stats
        """       
        # store sample stats
        for sample in self.samples:
            logger.info("store sample stats: {}".format(sample.name))
            store_sample_stats(sample, statsdicts[sample.name])

        # build Assembly stats
        self.data.stats_dfs.s5 = self.data._build_stat("s5")

        # write stats file
        self.data.stats_files.s5 = os.path.join(
            self.data.dirs.consens, 
            's5_consens_stats.txt')

        with open(self.data.stats_files.s5, 'w') as out:
            self.data.stats_dfs.s5.to_string(
                buf=out,
                formatters={
                    'clusters_total': '{:.0f}'.format,
                    'filtered_by_depth': '{:.0f}'.format,
                    'filtered_by_maxH': '{:.0f}'.format,
                    'filtered_by_maxN': '{:.0f}'.format,
                    'filtered_by_max_alleles': '{:.0f}'.format,
                    'reads_consens': '{:.0f}'.format,
                    'nsites': '{:.0f}'.format,
                    'nhetero': '{:.0f}'.format,
                    'heterozygosity': '{:.5f}'.format
                })



def make_chunks(data, sample, ncpus):
    """
    Split job into bits and pass to the client
    """
    # counter for split job submission
    num = 0

    # set optim size for chunks in N clusters. The first few chunks take longer
    # because they contain larger clusters, so we create 4X as many chunks as
    # processors so that they are split more evenly.
    optim = int(
        (sample.stats.clusters_total // ncpus) + \
        (sample.stats.clusters_total % ncpus))

    # open to clusters
    with gzip.open(sample.files.clusters, 'rb') as clusters:
        # create iterator to sample 2 lines at a time
        pairdealer = izip(*[iter(clusters)] * 2)

        # Use iterator to sample til end of cluster
        done = 0
        while not done:
            # grab optim clusters and write to file.
            done, chunk = clustdealer(pairdealer, optim)
            chunk = [i.decode() for i in chunk]

            # make file handle
            chunkhandle = os.path.join(
                data.tmpdir,
                "{}.chunk.{}.{}".format(sample.name, optim, num * optim))

            # write to file
            if chunk:
                with open(chunkhandle, 'wt') as outchunk:
                    outchunk.write("//\n//\n".join(chunk) + "//\n//\n")
                num += 1
    logger.debug("chunked {} to {} files".format(sample.files.clusters, num))



def process_chunks(data, sample, chunkfile, isref):
    """
    Remote callable function to use Processor class instance
    """
    proc = Processor(data, sample, chunkfile, isref)
    proc.run()
    return proc.counters, proc.filters     



def concat_catgs(data, sample, isref):
    """
    Concat catgs into a single sample catg and remove tmp files
    """
    # collect tmpcat files written by write_chunks()
    tmpcats = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcats.*".format(sample.name))
    )
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
    with h5py.File(sample.files.database, 'w') as ioh5:
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
    combs1 = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcons.*".format(sample.name)))
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
                        name, num = line.rsplit("_", 1)
                        line = name + "_{}".format(idx)
                        idx += 1
                        #int(num) + starter)
                    cstack.append(line.strip())
                out.write("\n".join(cstack) + "\n")
            os.remove(fname)


def concat_reference_consens(data, sample):
    """
    Concatenates consens bits into SAM for reference assemblies
    """
    # write sequences to SAM file
    combs1 = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcons.*".format(sample.name)))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    # open sam handle for writing to bam
    samfile = sample.files.consens.replace(".bam", ".sam")    
    with open(samfile, 'w') as outf:

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
        for fname in combs1:
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
    cmd = [ip.bins.samtools, "view", "-Sb", samfile]
    with open(sample.files.consens, 'w') as outbam:
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=outbam)
        comm = proc.communicate()[1]
    if proc.returncode:
        raise IPyradError("error in samtools: {}".format(comm))
    else:
        pass
        #os.remove(samfile)
        #for fname in combs1:
        #    os.remove(fname)


def store_sample_stats(sample, statsdicts):
    """
    Not parallel, store the sample objects stats
    """
    # record results
    xcounters = {
        "nconsens": 0,
        "heteros": 0,
        "nsites": 0,
    }
    xfilters = {
        "depth": 0,
        "maxh": 0,
        "maxn": 0,
        "maxalleles": 0,
    }

    # merge finished consens stats
    for counters, filters in statsdicts:
        # sum individual counters
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]

    # set Sample stats_dfs values
    if int(xcounters['nsites']):
        prop = int(xcounters["heteros"]) / float(xcounters['nsites'])
    else:
        prop = 0

    # store stats attributes to the sample
    sample.stats_dfs.s5.nsites = int(xcounters["nsites"])
    sample.stats_dfs.s5.nhetero = int(xcounters["heteros"])
    sample.stats_dfs.s5.filtered_by_depth = xfilters['depth']
    sample.stats_dfs.s5.filtered_by_maxH = xfilters['maxh']
    sample.stats_dfs.s5.filtered_by_maxN = xfilters['maxn']
    sample.stats_dfs.s5.filtered_by_max_alleles = int(xfilters['maxalleles'])
    sample.stats_dfs.s5.reads_consens = int(xcounters["nconsens"])
    sample.stats_dfs.s5.clusters_total = sample.stats_dfs.s3.clusters_total
    sample.stats_dfs.s5.heterozygosity = float(prop)

    # set the Sample stats summary value
    sample.stats.reads_consens = int(xcounters["nconsens"])

    # save state to Sample if successful
    if sample.stats.reads_consens:
        sample.stats.state = 5
    else:
        print("No clusters passed filtering in Sample: {}".format(sample.name))



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG")


    tdata = ip.load_json("/tmp/test-amaranth-ref.json")
    tdata.run("5", auto=True, force=True)
    print(tdata.stats)
    print(tdata.stats_dfs.s5)

    # tdata = ip.load_json("/tmp/test-amaranth.json")
    # tdata.run("5", auto=True, force=True)
    # print(tdata.stats)
    # print(tdata.stats_dfs.s5)

    # tdata = ip.load_json("/tmp/test-simpairddrad.json")
    # tdata.run("5", auto=True, force=True)
    # logger.info(tdata.stats_dfs.s5.T)

    # self.data.hackersonly.declone_PCR_duplicates:
    # tdata = ip.load_json("/tmp/test-amaranth-denovo.json")
    # tdata.ipcluster['cores'] = 4
    # tdata.run("5", auto=True, force=True)
    # logger.info(tdata.stats.T)
