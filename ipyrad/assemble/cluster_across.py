#!/usr/bin/env python

"""
cluster across samples using vsearch with options for
hierarchical clustering.
"""

# py2/3 compatible
from __future__ import print_function
try:
    from builtins import range
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

import os
import pty
import gzip
import glob
import time
import numba
import shutil
import random
import select
import socket
import logging
import warnings
import subprocess as sps

import numpy as np
import pandas as pd
import dask.array as da
import ipyrad
from .util import IPyradWarningExit, IPyradError, clustdealer, fullcomp


with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

LOGGER = logging.getLogger(__name__)

TODO = """

1. provide method to split/ contigs in step3.
2. fix progress bar for cluster2 in step3.
3. finish databasing here.
4. fix progress bar for cluster here.
5. write concat bams code here.
6. 
"""

class Step6:
    def __init__(self, data, samples, ipyclient, randomseed=0, force=False):
        self.data = data
        self.samples = samples
        self.ipyclient = ipyclient
        self.randomseed = randomseed
        self.setup_dirs(force)

        # groups/threading information
        self.cgroups = {}
        self.assign_groups()
        self.hostd = {}
        self.tune_hierarchical_threading()


    def setup_dirs(self, force=False):
        "set up across and tmpalign dirs and init h5 database file"
        self.data.dirs.across = os.path.realpath(
            os.path.join(
                self.data.paramsdict["project_dir"],
                "{}_across".format(self.data.name)))
        self.data.tmpdir = os.path.join(
            self.data.dirs.across,
            "{}-tmpalign".format(self.data.name))
        self.data.clust_database = os.path.join(
            self.data.dirs.across,
            "{}.clust.hdf5".format(self.data.name))

        # clear out
        if force:
            odir = self.data.dirs.across
            if os.path.exists(odir):
                shutil.rmtree(odir)

        # make dirs
        if not os.path.exists(self.data.dirs.across):
            os.mkdir(self.data.dirs.across)
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)


    def assign_groups(self):
        "assign samples to groups if not user provided for hierarchical clust"
        # use population info to split samples into groups; or assign random
        if self.data.populations:
            self.cgroups = {}
            for idx, val in enumerate(self.data.populations.values()):
                self.cgroups[idx] = val[1]

        # by default let's split taxa into groups of 20-50 samples at a time
        else:
            # calculate the number of cluster1 jobs to perform:
            if len(self.samples) <= 100:
                groupsize = 20
            elif len(self.samples) <= 500:
                groupsize = 50
            else:
                groupsize = 100

            alls = self.samples
            self.cgroups = {}
            idx = 0
            for samps in range(0, len(alls), groupsize):
                self.cgroups[idx] = alls[samps: samps + groupsize]
                idx += 1


    def tune_hierarchical_threading(self):
        "tune threads for across-sample clustering"

        # get engine data, skips busy engines.
        hosts = {}
        for eid in self.ipyclient.ids:
            engine = self.ipyclient[eid]
            if not engine.outstanding:
                hosts[eid] = engine.apply(socket.gethostname)

        # get targets on each hostname for spreading jobs out.
        self.ipyclient.wait()
        hosts = [(eid, i.get()) for (eid, i) in hosts.items()]
        hostnames = set([i[1] for i in hosts])
        self.hostd = {x: [i[0] for i in hosts if i[1] in x] for x in hostnames}

        # calculate the theading of cluster1 jobs:
        self.data.ncpus = len(self.ipyclient.ids)
        njobs = len(self.cgroups)
        nnodes = len(self.hostd)

        # how to load-balance cluster2 jobs
        # maxthreads = 8 cuz vsearch isn't v efficient above that.
        ## e.g., 24 cpus; do 2 12-threaded jobs
        ## e.g., 2 nodes; 40 cpus; do 2 20-threaded jobs or 4 10-threaded jobs
        ## e.g., 4 nodes; 80 cpus; do 8 10-threaded jobs
        if nnodes == 1:
            thr = np.floor(self.data.ncpus / njobs).astype(int)
            eids = max(1, thr)
            eids = max(eids, len(list(self.hostd.values())[0]))

        else:
            eids = []
            for node in self.hostd:
                sids = self.hostd[node]
                nids = len(sids)
                thr = np.floor(nids / (njobs / nnodes)).astype(int)
                thr = max(1, thr)
                thr = min(thr, nids)
                eids.extend(self.hostd[node][::thr])

        # set nthreads based on _ipcluster dict (default is 2)        
        #if "threads" in self.data._ipcluster.keys():
        #    self.nthreads = int(self.data._ipcluster["threads"])
        self.nthreads = 2
        if self.data.ncpus > 4:
            self.nthreads = 4
        eids = self.ipyclient.ids[::self.nthreads]

        # create load-balancers
        self.lbview = self.ipyclient.load_balanced_view()
        self.thview = self.ipyclient.load_balanced_view(targets=eids)


    def run(self):

        # DENOVO
        if self.data.paramsdict["assembly_method"] == "denovo":

            # prepare clustering inputs for hierarchical clustering
            self.remote_build_concats_tier1()

            # send initial clustering jobs (track finished jobs)
            self.remote_cluster1()

            # prepare second tier inputs
            self.remote_build_concats_tier2()

            # send cluster2 job (track actual progress)
            self.remote_cluster2()

            # build clusters
            self.remote_build_denovo_clusters()

            # align denovo clusters
            self.remote_align_clusters()

            # enter database values
            self.build_database()

        elif self.data.paramsdict["assembly_method"] == "reference":

            # prepare bamfiles
            self.remote_concat_bams()

            # build clusters from bedtools merge
            self.remote_build_ref_clusters()

            # enter database values
            self.build_database()


    def remote_run(self, printstr, function, args, threaded=False):
        # submit job
        start = time.time()
        rasyncs = {}
        for sample in self.samples:
            fargs = [self.data, sample] + list(args)
            if threaded:
                rasyncs[sample.name] = self.thview.apply(function, *fargs)
            else:
                rasyncs[sample.name] = self.lbview.apply(function, *fargs)

        # track job
        while 1:
            ready = [rasyncs[i].ready() for i in rasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        # check for errors
        print("")
        for job in rasyncs:
            if not rasyncs[job].successful():
                raise IPyradError(rasyncs[job].exception())


    def remote_build_concats_tier1(self):
        "prepares concatenated consens input files for each clust1 group"

        start = time.time()
        printstr = ("concatenating inputs", "s6")
        rasyncs = {}
        for jobid, group in self.cgroups.items():
            samples = [i for i in self.samples if i.name in group]
            args = (self.data, jobid, samples, self.randomseed)
            rasyncs[jobid] = self.lbview.apply(build_concat_files, *args)
        
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
                raise IPyradError(rasyncs[job].exception())


    def remote_cluster1(self):
        "send threaded jobs to remote engines"
        start = time.time()
        printstr = ("clustering across 1 ", "s6")
        rasyncs = {}
        for jobid in self.cgroups:
            args = (self.data, jobid, self.nthreads)
            rasyncs[jobid] = self.thview.apply(cluster, *args)
        
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
                raise IPyradError(rasyncs[job].exception())


    def remote_build_concats_tier2(self):
        start = time.time()
        printstr = ("concatenating inputs", "s6")
        args = (self.data, list(self.cgroups.keys()), self.randomseed)
        rasync = self.lbview.apply(build_concat_two, *args)
        
        while 1:
            ready = rasync.ready()
            self.data._progressbar(int(ready), 1, start, printstr)
            time.sleep(0.5)
            if ready:
                break

        # check for errors
        print("")
        if not rasync.successful():
            raise IPyradError(rasync.exception())        


    def remote_cluster2(self):
        start = time.time()
        printstr = ("clustering across 2 ", "s6")
        args = (self.data, 'x', 0)
        rasync = self.thview.apply(cluster, *args)
        
        while 1:
            time.sleep(0.2)
            self.data._progressbar(1, int(rasync.ready()), start, printstr)
            if rasync.ready():
                break

        # check for errors
        print("")
        if not rasync.successful():
            raise IPyradError(rasync.exception())          


    def remote_build_denovo_clusters(self):
       
        # filehandles
        uhandle = os.path.join(
            self.data.dirs.across, 
            "{}-x.utemp".format(self.data.name))
        usort = os.path.join(
            self.data.dirs.across, 
            "{}-x.utemp.sort".format(self.data.name))

        # sort utemp files, count seeds, ...
        start = time.time()
        printstr = ("building clusters   ", "s6")
        async1 = self.lbview.apply(sort_seeds, *(uhandle, usort))
        while 1:
            ready = [async1.ready()]
            self.data._progressbar(3, sum(ready), start, printstr)
            time.sleep(0.1)
            if all(ready):
                break

        async2 = self.lbview.apply(count_seeds, usort)
        while 1:
            ready = [async1.ready(), async2.ready()]
            self.data._progressbar(3, sum(ready), start, printstr)
            time.sleep(0.1)
            if all(ready):
                break
        nseeds = async2.result()

        # send the clust bit building job to work and track progress
        async3 = self.lbview.apply(
            build_denovo_clusters, 
            *(self.data, usort, nseeds, list(self.cgroups.keys())))
        while 1:
            ready = [async1.ready(), async2.ready(), async3.ready()]
            self.data._progressbar(3, sum(ready), start, printstr)
            time.sleep(0.1)
            if all(ready):
                break
        print("")

        # check for errors
        for job in [async1, async2, async3]:
            if not job.successful():
                raise IPyradWarningExit(job.result())


    def remote_align_denovo_clusters(self):
    
        # get files
        globpath = os.path.join(self.data.tmpdir, self.data.name + ".chunk_*")
        clustbits = glob.glob(globpath)

        # submit jobs to engines
        start = time.time()
        printstr = ("aligning clusters   ", "s6")
        jobs = {}
        for idx, _ in enumerate(clustbits):
            args = [self.data, self.samples, clustbits[idx]]
            jobs[idx] = self.lbview.apply(align_to_array, *args)
        allwait = len(jobs)

        # print progress while bits are aligning
        while 1:
            finished = [i.ready() for i in jobs.values()]
            fwait = sum(finished)
            self.data._progressbar(allwait, fwait, start, printstr)
            time.sleep(0.1)
            if all(finished):
                break

        # check for errors in muscle_align_across
        keys = list(jobs.keys())
        for idx in keys:
            if not jobs[idx].successful():
                raise IPyradWarningExit(
                    "error in step 6 {}".format(jobs[idx].exception()))
            del jobs[idx]
        print("")


    def remote_concat_bams(self):
        pass


    def remote_build_concat_bams(self):
        pass


def build_concat_bams(data):
    "create .sam files from consens data and make a concatenated bam file"
    pass


def build_concat_two(data, jobids, randomseed):
    seeds = [
        os.path.join(
            data.dirs.across, 
            "{}-{}.htemp".format(data.name, jobid)) for jobid in jobids
        ]
    allseeds = os.path.join(
        data.dirs.across, 
        "{}-x-catshuf.fa".format(data.name))
    cmd1 = ['cat'] + seeds
    cmd2 = [ipyrad.bins.vsearch, '--sortbylength', '-', '--output', allseeds]
    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, close_fds=True)
    proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, close_fds=True)
    proc2.communicate()
    proc1.stdout.close()


def build_concat_files(data, jobid, samples, randomseed):
    """
    [This is run on an ipengine]
    Make a concatenated consens file with sampled alleles (no RSWYMK/rswymk).
    Orders reads by length and shuffles randomly within length classes
    """
    conshandles = [
        sample.files.consens[0] for sample in samples if 
        sample.stats.reads_consens]
    conshandles.sort()
    assert conshandles, "no consensus files found"

    ## concatenate all of the gzipped consens files
    cmd = ['cat'] + conshandles
    groupcons = os.path.join(
        data.dirs.across, 
        "{}-{}-catcons.gz".format(data.name, jobid))
    LOGGER.debug(" ".join(cmd))
    with open(groupcons, 'w') as output:
        call = sps.Popen(cmd, stdout=output, close_fds=True)
        call.communicate()

    ## a string of sed substitutions for temporarily replacing hetero sites
    ## skips lines with '>', so it doesn't affect taxon names
    subs = ["/>/!s/W/A/g", "/>/!s/w/A/g", "/>/!s/R/A/g", "/>/!s/r/A/g",
            "/>/!s/M/A/g", "/>/!s/m/A/g", "/>/!s/K/T/g", "/>/!s/k/T/g",
            "/>/!s/S/C/g", "/>/!s/s/C/g", "/>/!s/Y/C/g", "/>/!s/y/C/g"]
    subs = ";".join(subs)

    ## impute pseudo-haplo information to avoid mismatch at hetero sites
    ## the read data with hetero sites is put back into clustered data later.
    ## pipe passed data from gunzip to sed.
    cmd1 = ["gunzip", "-c", groupcons]
    cmd2 = ["sed", subs]
    LOGGER.debug(" ".join(cmd1))
    LOGGER.debug(" ".join(cmd2))

    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, close_fds=True)
    allhaps = groupcons.replace("-catcons.gz", "-cathaps.fa")
    with open(allhaps, 'w') as output:
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=output, close_fds=True)
        proc2.communicate()
    proc1.stdout.close()

    ## now sort the file using vsearch
    allsort = groupcons.replace("-catcons.gz", "-catsort.fa")
    cmd1 = [ipyrad.bins.vsearch,
            "--sortbylength", allhaps,
            "--fasta_width", "0",
            "--output", allsort]
    LOGGER.debug(" ".join(cmd1))
    proc1 = sps.Popen(cmd1, close_fds=True)
    proc1.communicate()

    ## shuffle sequences within size classes. Tested seed (8/31/2016)
    ## shuffling works repeatably with seed.
    random.seed(randomseed)

    ## open an iterator to lengthsorted file and grab two lines at at time
    allshuf = groupcons.replace("-catcons.gz", "-catshuf.fa")
    outdat = open(allshuf, 'wt')
    indat = open(allsort, 'r')
    idat = izip(iter(indat), iter(indat))
    done = 0

    chunk = [next(idat)]
    while not done:
        ## grab 2-lines until they become shorter (unless there's only one)
        oldlen = len(chunk[-1][-1])
        while 1:
            try:
                dat = next(idat)
            except StopIteration:
                done = 1
                break
            if len(dat[-1]) == oldlen:
                chunk.append(dat)
            else:
                ## send the last chunk off to be processed
                random.shuffle(chunk)
                outdat.write("".join(chain(*chunk)))
                ## start new chunk
                chunk = [dat]
                break

    ## do the last chunk
    random.shuffle(chunk)
    outdat.write("".join(chain(*chunk)))

    indat.close()
    outdat.close()


def cluster(data, jobid, nthreads):

    # get files for this jobid
    catshuf = os.path.join(
        data.dirs.across, 
        "{}-{}-catshuf.fa".format(data.name, jobid))
    uhaplos = os.path.join(
        data.dirs.across, 
        "{}-{}.utemp".format(data.name, jobid))
    hhaplos = os.path.join(
        data.dirs.across, 
        "{}-{}.htemp".format(data.name, jobid))

    ## parameters that vary by datatype
    ## (too low of cov values yield too many poor alignments)
    strand = "plus"
    cov = 0.75         # 0.90
    if data.paramsdict["datatype"] in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.60
    elif data.paramsdict["datatype"] == "pairgbs":
        strand = "both"
        cov = 0.75     # 0.90

    cmd = [ipyrad.bins.vsearch,
           "-cluster_smallmem", catshuf,
           "-strand", strand,
           "-query_cov", str(cov),
           "-minsl", str(0.5),
           "-id", str(data.paramsdict["clust_threshold"]),
           "-userout", uhaplos,
           "-notmatched", hhaplos,
           "-userfields", "query+target+qstrand",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-fasta_width", "0",
           "-threads", str(nthreads),  # "0",
           "-fulldp",
           "-usersort",
           ]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    out = proc.communicate()
    if proc.returncode:
        raise IPyradError(out)


def count_seeds(usort):
    """
    uses bash commands to quickly count N seeds from utemp file
    """
    with open(usort, 'r') as insort:
        cmd1 = ["cut", "-f", "2"]
        cmd2 = ["uniq"]
        cmd3 = ["wc"]
        proc1 = sps.Popen(cmd1, stdin=insort, stdout=sps.PIPE, close_fds=True)
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, close_fds=True)
        proc3 = sps.Popen(cmd3, stdin=proc2.stdout, stdout=sps.PIPE, close_fds=True)
        res = proc3.communicate()
        nseeds = int(res[0].split()[0])
        proc1.stdout.close()
        proc2.stdout.close()
        proc3.stdout.close()
    return nseeds


def sort_seeds(uhandle, usort):
    """ sort seeds from cluster results"""
    cmd = ["sort", "-k", "2", uhandle, "-o", usort]
    proc = sps.Popen(cmd, close_fds=True)
    proc.communicate()


def build_denovo_clusters(data, usort, nseeds, jobids):

    # load all concat fasta files into a dictionary (memory concerns here...)
    allcons = {}
    conshandles = [
        os.path.join(
            data.dirs.across, "{}-{}-catcons.gz".format(data.name, jobid))
        for jobid in jobids]
    for conshandle in conshandles:
        subcons = {}
        with gzip.open(conshandle, 'rt') as iocons:
            cons = izip(*[iter(iocons)] * 2)
            for namestr, seq in cons:
                nnn, sss = [i.strip() for i in (namestr, seq)]
                subcons[nnn[1:]] = sss
        allcons.update(subcons)
        del subcons

    # load all utemp files into a dictionary
    subdict = {}
    usortfiles = [
        os.path.join(data.dirs.across, "{}-{}.utemp".format(data.name, jobid))
        for jobid in jobids
    ]
    for ufile in usortfiles:
        with open(ufile, 'r') as inhits:
            for line in inhits:
                hit, seed, ori = line.strip().split()
                if seed not in subdict:
                    subdict[seed] = [(hit, ori)]
                else:
                    subdict[seed].append((hit, ori))

    # set optim to approximately 4 chunks per core. Smaller allows for a bit
    # cleaner looking progress bar. 40 cores will make 160 files.
    optim = ((nseeds // (data.ncpus * 4)) + (nseeds % (data.ncpus * 4)))

    # iterate through usort grabbing seeds and matches
    with open(usort, 'rt') as insort:
        # iterator, seed null, and seqlist null
        isort = iter(insort)
        loci = 0
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0

        while 1:
            try:
                hit, seed, ori = next(isort).strip().split()
            except StopIteration:
                break
        
            # if same seed append match
            if seed != lastseed:
                # store the last fseq, count it, and clear it
                if fseqs:
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                # occasionally write to file
                if seqsize >= optim:
                    if seqlist:
                        loci += seqsize
                        pathname = os.path.join(
                            data.tmpdir, 
                            "{}.chunk_{}".format(data.name, loci))
                        with open(pathname, 'wt') as clustout:
                            clustout.write(
                                "\n//\n//\n".join(seqlist) + "\n//\n//\n")
                        # reset counter and list
                        seqlist = []
                        seqsize = 0

                # store the new seed on top of fseqs
                fseqs.append(">{}\n{}".format(seed, allcons[seed]))
                lastseed = seed

                # expand subhits to seed
                uhits = subdict.get(seed)
                if uhits:
                    for ahit, ori in uhits:
                        if ori == "-":
                            seq = fullcomp(allcons[ahit])[::-1]
                        else:
                            seq = allcons[ahit]
                        fseqs.append(">{}\n{}".format(ahit, seq))

            # expand matches with subdict
            hitseqs = [(hit, allcons[hit], ori)]
            uhits = subdict.get(hit)
            if uhits:
                for hit in uhits:
                    hitseqs.append((hit[0], allcons[hit[0]], hit[1]))

            # revcomp if orientation is reversed
            for sname, seq, ori in hitseqs:
                if ori == "-":
                    seq = fullcomp(seq)[::-1]
                fseqs.append(">{}\n{}".format(sname, seq))

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
        seqsize += 1
        loci += seqsize
    if seqlist:
        pathname = os.path.join(data.tmpdir, 
            data.name + ".chunk_{}".format(loci))
        with open(pathname, 'wt') as clustsout:
            clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")

    ## final progress and cleanup
    del allcons


def align_to_array(data, samples, chunk):

    # data are already chunked, read in the whole thing
    with open(chunk, 'rt') as infile:
        clusts = infile.read().split("//\n//\n")[:-1]    

    # snames to ensure sorted order
    samples.sort(key=lambda x: x.name)
    snames = [sample.name for sample in samples]

    # make a tmparr to store metadata (this can get huge, consider using h5)
    maxvar = 25
    maxlen = data._hackersonly["max_fragment_length"] + 20
    maxinds = data.paramsdict["max_Indels_locus"]
    ifilter = np.zeros(len(clusts), dtype=np.bool_)
    dfilter = np.zeros(len(clusts), dtype=np.bool_)
    varpos = np.zeros((len(clusts), maxvar), dtype='u1')
    indels = np.zeros((len(clusts), len(samples), sum(maxinds)), dtype='u1')
    edges = np.zeros((len(clusts), len(samples), 4), dtype='u1')
    consens = np.zeros((len(clusts), maxlen), dtype="S1")

    # create a persistent shell for running muscle in. 
    proc = sps.Popen(["bash"], stdin=sps.PIPE, stdout=sps.PIPE, bufsize=0)

    # iterate over clusters until finished
    allstack = []
    for ldx in range(len(clusts)):
        istack = []
        lines = clusts[ldx].strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]
        seqs = [i[:90] for i in seqs]    
        align1 = []
        
        # find duplicates and skip aligning but keep it for downstream.
        unames = set([i.rsplit("_", 1)[0] for i in names])
        if len(unames) < len(names):
            dfilter[ldx] = True
            istack = ["{}\n{}".format(i[1:], j) for i, j in zip(names, seqs)]

        else:
            # append counter to names because muscle doesn't retain order
            nnames = [">{};*{}".format(j[1:], i) for i, j in enumerate(names)]

            # make back into strings
            cl1 = "\n".join(["\n".join(i) for i in zip(nnames, seqs)])                

            # store allele (lowercase) info
            amask, abool = store_alleles(seqs)

            # send align1 to the bash shell (TODO: check for pipe-overflow)
            cmd1 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
                    .format(cl1, ipyrad.bins.muscle, "//\n"))
            proc.stdin.write(cmd1.encode())

            # read the stdout by line until splitter is reached
            for line in iter(proc.stdout.readline, b'//\n'):
                align1.append(line.decode())

            # reorder b/c muscle doesn't keep order
            lines = "".join(align1)[1:].split("\n>")
            dalign1 = dict([i.split("\n", 1) for i in lines])
            keys = sorted(
                dalign1.keys(), 
                key=lambda x: int(x.rsplit("*")[-1])
            )
            seqarr = np.zeros(
                (len(nnames), len(dalign1[keys[0]].replace("\n", ""))),
                dtype='S1',
                )
            for kidx, key in enumerate(keys):
                concatseq = dalign1[key].replace("\n", "")
                seqarr[kidx] = list(concatseq)
                    
                # fill in edges for each seq
                edg = (
                    len(concatseq) - len(concatseq.lstrip('-')),
                    len(concatseq.rstrip('-')),
                    len(concatseq.rstrip('-')),
                    len(concatseq.rstrip('-')),
                    )
                sidx = snames.index(key.rsplit('_', 1)[0])
                edges[ldx, sidx] = edg

            # get ref/variant counts in an array
            varmat = refbuild(seqarr.view('u1'), PSEUDO_REF)
            constmp = varmat[:, 0].view('S1')
            consens[ldx, :constmp.size] = constmp
            varstmp = np.nonzero(varmat[:, 1])[0]
            varpos[ldx, :varstmp.size] = varstmp 

            # impute allele information back into the read b/c muscle
            wkeys = np.argsort([i.rsplit("_", 1)[0] for i in keys])

            # sort in sname (alphanumeric) order for indels storage
            for idx in wkeys:
                # seqarr and amask are in input order here
                args = (seqarr, idx, amask)
                seqarr[idx], indidx = retrieve_indels_and_alleles(*args)
                sidx = snames.index(names[idx].rsplit("_", 1)[0][1:])
                if indidx.size:
                    indsize = min(indidx.size, sum(maxinds))
                    indels[ldx, sidx, :indsize] = indidx[:indsize]
                wname = names[idx]

                # store (only for visually checking alignments)
                istack.append(
                    "{}\n{}".format(wname, b"".join(seqarr[idx]).decode()))

            # indel filter
            if indels.sum(axis=1).max() >= sum(maxinds):
                ifilter[ldx] = True

        # store the stack (only for visually checking alignments)
        if istack:
            allstack.append("\n".join(istack))

    # cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    # write to file after (only for visually checking alignments)
    odx = chunk.rsplit("_")[-1]
    alignfile = os.path.join(data.tmpdir, "aligned_{}.fa".format(odx))
    with open(alignfile, 'wt') as outfile:
        outfile.write("\n//\n//\n".join(allstack) + "\n")
    #os.remove(chunk)

    # save indels array to tmp dir
    np.save(
        os.path.join(data.tmpdir, "ifilter_{}.tmp.npy".format(odx)),
        ifilter)
    np.save(
        os.path.join(data.tmpdir, "dfilter_{}.tmp.npy".format(odx)),
        dfilter)


def store_alleles(seqs):
    shape = (len(seqs), max([len(i) for i in seqs]))
    arrseqs = np.zeros(shape, dtype=np.bytes_)
    for row in range(arrseqs.shape[0]):
        seqsrow = seqs[row]
        arrseqs[row, :len(seqsrow)] = list(seqsrow)
    amask = np.char.islower(arrseqs)
    if np.any(amask):
        return amask, True
    else:
        return amask, False


def retrieve_indels_and_alleles(seqarr, idx, amask):
    concatarr = seqarr[idx]
    newmask = np.zeros(len(concatarr), dtype=np.bool_)                        

    # check for indels and impute to amask
    indidx = np.where(concatarr == b"-")[0]
    if np.sum(amask):
        if indidx.size:

            # impute for position of variants
            allrows = np.arange(amask.shape[1])
            mask = np.ones(allrows.shape[0], dtype=np.bool_)
            for idx in indidx:
                if idx < mask.shape[0]:
                    mask[idx] = False
            not_idx = allrows[mask == 1]

            # fill in new data into all other spots
            newmask[not_idx] = amask[idx, :not_idx.shape[0]]

        else:
            newmask = amask[idx]
                
        # lower the alleles
        concatarr[newmask] = np.char.lower(concatarr[newmask])

    return concatarr, indidx


def init_database(data, samples, nloci):
    """
    generate empty h5 database to be filled and setup chunk sizes
    """   
    # sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    # get maxlen dim
    maxlen = data._hackersonly["max_fragment_length"] + 20
    maxvar = data.paramsdict[""]

    # open file handle
    with h5py.File(data.clust_database, 'w') as io5:

        # chunk to approximately 2 chunks per core
        chunks = ((nloci // (data.cpus * 2)) + (nloci % (data.cpus * 2)))
        chunklen = chunks * len(samples) * maxlen * 4
        while chunklen > int(500e6):
            chunks = (chunks // 2) + (chunks % 2)
            chunklen = chunks * len(samples) * maxlen * 4
        data.chunks = chunks

        # INIT FULL CATG ARRAY
        io5.attrs["samples"] = [i.name.encode() for i in samples]
        io5.create_dataset(
            name="catg", 
            shape=(nloci, len(samples), maxvar, 4),
            dtype=np.uint32,
            chunks=(chunks, 1, maxvar, 4),
            compression="gzip")
        io5.create_dataset(
            name="cons", 
            shape=(nloci, len(samples), maxlen),
            dtype="S1",
            chunks=(chunks, len(samples), maxlen),
            compression='gzip')
        io5.create_dataset(
            name="edge", 
            shape=(nloci, len(samples), 4),
            dtype="u1",
            chunks=(chunks, len(samples), 4),
            compression='gzip')
        io5.create_dataset(
            name="chro", 
            shape=(nloci, 3), 
            dtype=np.int64, 
            chunks=(chunks, 3),
            compression="gzip")

        # filters
        io5.create_dataset("dfilter", (nloci, ), dtype=np.bool_)
        io5.create_dataset("ifilter", (nloci, ), dtype=np.bool_)
        io5.create_dataset("afilter", (nloci, ), dtype=np.bool_)


def fill_database(data, samples, nloci):
    """
    fill h5 database with npy tmp files from align steps.
    """
    # init datasets in the h5 array in appropriate size and chunking
    init_database(data, samples, nloci)

    # collect all arrays
    glob.glob("")

    # enter consens and filters into h5 datasets
    # ...

    # construct subarrays to collect depths of variants with indels imputed.
    # ...

    # enter variants into the variant h5 dataset
    # ...




##########################################

@numba.jit(nopython=True)
def refbuild(iseq, consdict):
    """ Returns the most common base at each site in order. """

    altrefs = np.zeros((iseq.shape[1], 4), dtype=np.uint8)

    for col in range(iseq.shape[1]):
        ## expand colums with ambigs and remove N-
        fcounts = np.zeros(111, dtype=np.int64)
        counts = np.bincount(iseq[:, col])
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


PSEUDO_REF = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    [78, 78, 78],
    [45, 45, 45],
    ], dtype=np.uint8)
 
#############################################
















def persistent_popen_align3(data, samples, chunk):
    """  notes """

    ## data are already chunked, read in the whole thing
    with open(chunk, 'rt') as infile:
        clusts = infile.read().split("//\n//\n")[:-1]    

    ## snames to ensure sorted order
    samples.sort(key=lambda x: x.name)
    snames = [sample.name for sample in samples]

    ## make a tmparr to store metadata (this can get huge, consider using h5)
    maxlen = data._hackersonly["max_fragment_length"] + 20
    indels = np.zeros((len(samples), len(clusts), maxlen), dtype=np.bool_)
    duples = np.zeros(len(clusts), dtype=np.bool_)

    ## create a persistent shell for running muscle in. 
    proc = sps.Popen(["bash"], stdin=sps.PIPE, stdout=sps.PIPE, bufsize=0)

    ## iterate over clusters until finished
    allstack = []
    for ldx in range(len(clusts)):
        ## new alignment string for read1s and read2s
        aligned = []
        istack = []
        lines = clusts[ldx].strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]        
        align1 = []
        align2 = ""
        
        # find duplicates and skip aligning but keep it for downstream.
        if len(names) == len(set([i.rsplit("_", 1)[0] for i in names])):
            duples[ldx] = True
            istack = ["{}\n{}".format(i[1:], j) for i, j in zip(names, seqs)]

        else:
            # append counter to names because muscle doesn't retain order
            names = [">{};*{}".format(j[1:], i) for i, j in enumerate(names)]

            try:
                ## try to split names on nnnn splitter
                clust1, clust2 = zip(*[i.split("nnnn") for i in seqs])

                ## make back into strings
                cl1 = "\n".join(chain(*zip(names, clust1)))
                cl2 = "\n".join(chain(*zip(names, clust2)))

                ## store allele (lowercase) info
                shape = (len(seqs), max([len(i) for i in seqs]))
                arrseqs = np.zeros(shape, dtype="S1")
                for row in range(arrseqs.shape[0]):
                    seqsrow = seqs[row]
                    arrseqs[row, :len(seqsrow)] = list(seqsrow)
                amask = np.char.islower(arrseqs)
                save_alleles = np.any(amask)

                ## send align1 to the bash shell
                ## TODO: check for pipe-overflow here and use files for i/o                
                cmd1 = "echo -e '{}' | {} -quiet -in - ; echo {}"\
                       .format(cl1, ipyrad.bins.muscle, "//")
                print(cmd1, file=proc.stdin)

                ## read the stdout by line until splitter is reached
                for line in iter(proc.stdout.readline, "//\n"):
                    align1 += line

                ## send align2 to the bash shell
                ## TODO: check for pipe-overflow here and use files for i/o                
                cmd2 = "echo -e '{}' | {} -quiet -in - ; echo {}"\
                       .format(cl2, ipyrad.bins.muscle, "//")
                print(cmd2, file=proc.stdin)

                ## read the stdout by line until splitter is reached
                for line in iter(proc.stdout.readline, "//\n"):
                    align2 += line

                ## join the aligned read1 and read2 and ensure name order match
                la1 = align1[1:].split("\n>")
                la2 = align2[1:].split("\n>")
                dalign1 = dict([i.split("\n", 1) for i in la1])
                dalign2 = dict([i.split("\n", 1) for i in la2])
                keys = sorted(dalign1.keys(), key=get_derep_num)
                keys2 = sorted(dalign2.keys(), key=get_derep_num)

                ## Make sure R1 and R2 actually exist for each sample. If not
                ## bail out of this cluster.
                if not len(keys) == len(keys2):
                    LOGGER.error("R1 and R2 results differ in length: "\
                               + "\nR1 - {}\nR2 - {}".format(keys, keys2))
                    continue

                ## impute allele (lowercase) info back into alignments
                for kidx, key in enumerate(keys):
                    concatseq = dalign1[key].replace("\n", "") + \
                        "nnnn" + dalign2[key].replace("\n", "")

                    ## impute alleles
                    if save_alleles:
                        newmask = np.zeros(len(concatseq), dtype=np.bool_)                        
                        ## check for indels and impute to amask
                        indidx = np.where(np.array(list(concatseq)) == "-")[0]
                        if indidx.size:
                            allrows = np.arange(amask.shape[1])
                            mask = np.ones(allrows.shape[0], dtype=np.bool_)
                            for idx in indidx:
                                if idx < mask.shape[0]:
                                    mask[idx] = False
                            not_idx = allrows[mask == 1]
                            ## fill in new data into all other spots
                            newmask[not_idx] = amask[kidx, :not_idx.shape[0]]
                        else:
                            newmask = amask[kidx]
                        
                        ## lower the alleles
                        concatarr = np.array(list(concatseq))
                        concatarr[newmask] = np.char.lower(concatarr[newmask])
                        concatseq = concatarr.tostring()
                        #LOGGER.info(concatseq)
                        
                    ## fill list with aligned data
                    aligned.append("{}\n{}".format(key, concatseq))

                ## put into a dict for writing to file
                #aligned = []
                #for key in keys:
                #    aligned.append("\n".join(
                #        [key, 
                #         dalign1[key].replace("\n", "")+"nnnn"+\
                #         dalign2[key].replace("\n", "")]))
            except IndexError as inst:
                LOGGER.debug("Error in PE - ldx: {}".format())
                LOGGER.debug("Vars: {}".format(dict(globals(), **locals())))
                raise

            except ValueError:
                ## make back into strings
                cl1 = "\n".join(["\n".join(i) for i in zip(names, seqs)])                

                ## store allele (lowercase) info
                shape = (len(seqs), max([len(i) for i in seqs]))
                arrseqs = np.zeros(shape, dtype=np.bytes_)
                for row in range(arrseqs.shape[0]):
                    seqsrow = seqs[row]
                    arrseqs[row, :len(seqsrow)] = list(seqsrow)
                amask = np.char.islower(arrseqs)
                save_alleles = np.any(amask)

                # send align1 to the bash shell (TODO: check for pipe-overflow)
                cmd1 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
                        .format(cl1, ipyrad.bins.muscle, "//\n"))
                proc.stdin.write(cmd1.encode())

                ## read the stdout by line until splitter is reached
                for line in iter(proc.stdout.readline, b'//\n'):
                    align1.append(line.decode())

                ## ensure name order match
                lines = "".join(align1)[1:].split("\n>")
                dalign1 = dict([i.split("\n", 1) for i in lines])
                keys = sorted(
                    dalign1.keys(), 
                    key=lambda x: int(x.rsplit("*")[-1])
                    )

                ## put into dict for writing to file
                for kidx, key in enumerate(keys):
                    concatseq = dalign1[key].replace("\n", "")
                    
                    ## impute alleles
                    if save_alleles:
                        newmask = np.zeros(len(concatseq), dtype=np.bool_)                        
                    
                        ## check for indels and impute to amask
                        indidx = np.where(np.array(list(concatseq)) == "-")[0]
                        if indidx.size:
                            allrows = np.arange(amask.shape[1])
                            mask = np.ones(allrows.shape[0], dtype=np.bool_)
                            for idx in indidx:
                                if idx < mask.shape[0]:
                                    mask[idx] = False
                            not_idx = allrows[mask == 1]
                    
                            ## fill in new data into all other spots
                            newmask[not_idx] = amask[kidx, :not_idx.shape[0]]

                        else:
                            newmask = amask[kidx]
                        
                        ## lower the alleles
                        concatarr = np.array(list(concatseq))
                        concatarr[newmask] = np.char.lower(concatarr[newmask])
                        concatseq = concatarr.tostring()

                    ## fill list with aligned data
                    aligned.append("{}\n{}".format(key, concatseq))

            ## enforce maxlen on aligned seqs
            aseqs = np.vstack([list(i.split("\n")[1]) for i in aligned])
            LOGGER.info("\naseqs here: %s", aseqs)

            ## index names by snames order
            sidxs = [snames.index(key.rsplit("_", 1)[0]) for key in keys]
            thislen = min(maxlen, aseqs.shape[1])
            for idx in range(aseqs.shape[0]):
                ## enter into stack
                newn = aligned[idx].split(";", 1)[0]
                #newn = key[idx].split(";", 1)[0]
                istack.append(
                    "{}\n{}".format(newn, "".join(aseqs[idx, :thislen])))               
                #  "{}\n{}".format(newn, aseqs[idx, :thislen].tostring()))
                ## name index in sorted list (indels order)
                sidx = sidxs[idx]
                indels[sidx, ldx, :thislen] = aseqs[idx, :thislen] == "-"

        if istack:
            allstack.append("\n".join(istack))
            #LOGGER.debug("\n\nSTACK (%s)\n%s\n", duples[ldx], "\n".join(istack))

    ## cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    #LOGGER.info("\n\nALLSTACK %s\n", "\n".join(i) for i in allstack[:5]])

    ## write to file after
    odx = chunk.rsplit("_")[-1]
    alignfile = os.path.join(data.tmpdir, "align_{}.fa".format(odx))
    with open(alignfile, 'wt') as outfile:
        outfile.write("\n//\n//\n".join(allstack) + "\n")
        os.remove(chunk)

    ## save indels array to tmp dir
    ifile = os.path.join(data.tmpdir, "indels_{}.tmp.npy".format(odx))
    np.save(ifile, indels)
    dfile = os.path.join(data.tmpdir, "duples_{}.tmp.npy".format(odx))
    np.save(dfile, duples)



def multi_muscle_align(data, samples, ipyclient):
    """
    Sends the cluster bits to nprocessors for muscle alignment. They return
    with indel.h5 handles to be concatenated into a joint h5.
    """
    LOGGER.info("starting alignments")

    ## get client
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    printstr = ("aligning clusters   ", "s6")
    data._progressbar(20, 0, start, printstr)

    # submit clustbits as jobs to engines. The chunkfiles are removed when they
    # are finished so this job can even be restarted if it was half finished, 
    # though that is probably rare. 
    path = os.path.join(data.tmpdir, data.name + ".chunk_*")
    clustbits = glob.glob(path)
    jobs = {}
    for idx in range(len(clustbits)):
        args = [data, samples, clustbits[idx]]
        jobs[idx] = lbview.apply(persistent_popen_align3, *args)
    allwait = len(jobs)
    data._progressbar(20, 0, start, printstr)

    ## print progress while bits are aligning
    while 1:
        finished = [i.ready() for i in jobs.values()]
        fwait = sum(finished)
        data._progressbar(allwait, fwait, start, printstr)
        time.sleep(0.1)
        if all(finished):
            break

    ## check for errors in muscle_align_across
    keys = list(jobs.keys())
    for idx in keys:
        if not jobs[idx].successful():
            LOGGER.error(
                "error in persistent_popen_align %s", jobs[idx].exception())
            raise IPyradWarningExit(
                "error in step 6 {}".format(jobs[idx].exception()))
        del jobs[idx]
    print("")



def concatclusts(outhandle, alignbits):
    """ concatenates sorted aligned cluster tmpfiles and removes them."""
    with gzip.open(outhandle, 'wt') as out:
        for fname in alignbits:
            with open(fname) as infile:
                out.write(infile.read() + "//\n//\n")
            #os.remove(fname)



def build_indels(data, samples, ipyclient):
    """
    Builds the indels array and catclust.gz file from the aligned clusters.
    Building catclust is very fast. Entering indels into h5 array is a bit
    slow but can probably be sped up. (todo). NOT currently parallelized.
    """

    ## progress bars
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    printstr = ("database indels     ", "s6")
    njobs = len(glob.glob(os.path.join(data.tmpdir, "align_*.fa"))) + 1

    ## build tmparrs
    rasync = lbview.apply(build_tmp_h5, *(data, samples))

    ## track progress
    while 1:
        ready = bool(rasync.ready())
        data._progressbar(njobs, ready, start, printstr)
        if ready:
            break
        else:
            time.sleep(0.1)

    ## check for errors
    if not rasync.successful():
        raise IPyradWarningExit(rasync.result())

    ## start subfunc
    rasync = lbview.apply(sub_build_indels, *(data, samples))
    
    prog = 1
    while 1:
        if rasync.stdout:
            prog = int(rasync.stdout.split()[-1]) + 1
        data._progressbar(njobs, prog, start, printstr)
        if rasync.ready():
            break
        else:
            time.sleep(0.1)

    ## check for errors
    if not rasync.successful():
        raise IPyradWarningExit(rasync.result())
    print("")

    ## prepare for next substep by removing the singlecat result files if 
    ## they exist. 
    snames = [i.name for i in samples]
    snames.sort()
    smpios = [os.path.join(data.dirs.across, i + '.tmp.h5') for i in snames]
    for smpio in smpios:
        if os.path.exists(smpio):
            os.remove(smpio)



def sub_build_indels(data, samples):
    """ sub func in `build_indels()`. """

    ## get file handles
    indelfiles = glob.glob(os.path.join(data.tmpdir, "indels_*.tmp.npy"))
    alignbits = glob.glob(os.path.join(data.tmpdir, "align_*.fa"))

    ## sort into input order by chunk names
    indelfiles.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))
    alignbits.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-3]))
    LOGGER.info("indelfiles %s", indelfiles)
    LOGGER.info("alignbits %s", alignbits)
    chunksize = int(indelfiles[0].rsplit("_", 1)[-1][:-8])

    ## concatenate finished seq clusters into a tmp file
    outhandle = os.path.join(data.dirs.across, data.name + "_catclust.gz")
    concatclusts(outhandle, alignbits)

    ## get dims for full indel array
    maxlen = data._hackersonly["max_fragment_length"] + 20
    nloci = get_nloci(data)
    LOGGER.info("maxlen inside build is %s", maxlen)
    LOGGER.info("nloci for indels %s", nloci)

    ## INIT TEMP INDEL ARRAY
    ## build an indel array for ALL loci in cat.clust.gz,
    ## chunked so that individual samples can be pulled out
    ipath = os.path.join(data.dirs.across, data.name + ".tmp.indels.hdf5")
    with h5py.File(ipath, 'w') as io5:
        iset = io5.create_dataset(
            "indels",
            shape=(len(samples), nloci, maxlen),
            dtype=np.bool_,
            chunks=(1, chunksize, maxlen))

        ## again make sure names are ordered right
        samples.sort(key=lambda x: x.name)

        #iset.attrs["chunksize"] = (1, data.nloci, maxlen)
        iset.attrs["samples"] = [i.name.encode() for i in samples]

        ## enter all tmpindel arrays into full indel array
        done = 0
        init = 0
        for indf in indelfiles:
            end = int(indf.rsplit("_", 1)[-1][:-8])
            inarr = np.load(indf)
            LOGGER.info('inarr shape %s', inarr.shape)
            LOGGER.info('iset shape %s', iset.shape)
            iset[:, init:end, :] = inarr[:, :end - init]
            init += end - init
            done += 1
            print(done)



def call_cluster(data, noreverse, ipyclient):
    """
    distributes 'cluster()' function to an ipyclient to make sure it runs
    on a high memory node. 
    """
    ## Find host with the most engines, for now just using first.
    # lbview = ipyclient.load_balanced_view()

    ## request engine data, skips busy engines.    
    asyncs = {}
    for eid in ipyclient.ids:
        engine = ipyclient[eid]
        if not engine.outstanding:
            asyncs[eid] = engine.apply(socket.gethostname)
    ## get results
    hosts = {}
    for key in asyncs:
        hosts[key] = asyncs[key].get()
    ## count them
    results = {}
    for eid, hostname in hosts.items():
        if hostname in results:
            results[hostname].append(eid)
        else:
            results[hostname] = [eid] 

    ## which is largest
    hosts = sorted(results.items(), key=lambda x: len(x[1]), reverse=True)
    _, eids = hosts[0]
    bighost = ipyclient[eids[0]]

    # nthreads is len eids, or ipcluster.threads, unless ipcluster.threads 
    # is really small, then we assume threads should not apply here.
    #    ipyrad -p params.txt -s 6 -c 20 would give:
    #    min(20, max(2, 10)) = 8
    # while 
    #    ipyrad -p params.txt -s 6 -c 20 -t 4 would give:
    #    min(20, max(4, 10)) = 10
    # and 
    #    ipyrad -p params.txt -s 6 -c 20 -t 15 would give:
    #    min(20, max(15, 10)) = 15
    # and
    #    ipyrad -p params.txt -s 6 -c 16 --MPI (on 2 X 8-core nodes) gives:
    #    min(8, max(2, 10)) = 8
    nthreads = min(len(eids), max(data._ipcluster["threads"], 10))

    ## submit job to the host with the most
    rasync = bighost.apply(cluster, *(data, noreverse, nthreads))
    #async = lbview.apply(cluster, *(data, noreverse, nthreads))
    
    ## track progress
    prog = 0
    start = time.time()
    printstr = ("clustering across   ", "s6")
    
    while 1:
        if rasync.stdout:
            prog = int(rasync.stdout.split()[-1])
        data._progressbar(100, prog, start, printstr)
        if rasync.ready():
            data._progressbar(100, prog, start, printstr)
            print("")
            break
        else:
            time.sleep(0.5)

    ## store log result
    ipyclient.wait()
    data.stats_files.s6 = os.path.join(
        data.dirs.across, "s6_cluster_stats.txt")



# def cluster(data, noreverse, nthreads):
#     """
#     Calls vsearch for clustering across samples.
#     """

#     ## input and output file handles
#     cathaplos = os.path.join(data.dirs.across, data.name + "_catshuf.tmp")
#     uhaplos = os.path.join(data.dirs.across, data.name + ".utemp")
#     hhaplos = os.path.join(data.dirs.across, data.name + ".htemp")
#     logfile = os.path.join(data.dirs.across, "s6_cluster_stats.txt")

#     ## parameters that vary by datatype
#     ## (too low of cov values yield too many poor alignments)
#     strand = "plus"
#     cov = 0.75    # 0.90
#     if data.paramsdict["datatype"] in ["gbs", "2brad"]:
#         strand = "both"
#         cov = 0.60
#     elif data.paramsdict["datatype"] == "pairgbs":
#         strand = "both"
#         cov = 0.75   # 0.90

#     ## nthreads is calculated in 'call_cluster()'
#     cmd = [ipyrad.bins.vsearch,
#            "-cluster_smallmem", cathaplos,
#            "-strand", strand,
#            "-query_cov", str(cov),
#            "-minsl", str(0.5),
#            "-id", str(data.paramsdict["clust_threshold"]),
#            "-userout", uhaplos,
#            "-notmatched", hhaplos,
#            "-userfields", "query+target+qstrand",
#            "-maxaccepts", "1",
#            "-maxrejects", "0",
#            "-fasta_width", "0",
#            "-threads", str(nthreads),  # "0",
#            "-fulldp",
#            "-usersort",
#            "-log", logfile]

#     ## override reverse clustering option
#     if noreverse:
#         strand = "plus"  # -leftjust "

#     try:
#         ## this seems to start vsearch on a different pid than the engine
#         ## and so it's hard to kill... 
#         LOGGER.info(cmd)
#         (dog, owner) = pty.openpty()
#         proc = sps.Popen(cmd, stdout=owner, stderr=owner, close_fds=True)
                                     
#         prog = 0
#         newprog = 0
#         while 1:
#             isdat = select.select([dog], [], [], 0)
#             if isdat[0]:
#                 dat = os.read(dog, 80192).decode()
#             else:
#                 dat = ""
#             if "Clustering" in dat:
#                 try:
#                     newprog = int(dat.split()[-1][:-1])
#                 ## may raise value error when it gets to the end
#                 except ValueError:
#                     pass

#             ## break if done
#             ## catches end chunk of printing if clustering went really fast
#             elif "Clusters:" in dat:
#                 LOGGER.info("ended vsearch tracking loop")
#                 break
#             else:
#                 time.sleep(0.1)
#             ## print progress
#             if newprog != prog:
#                 print(newprog)
#                 prog = newprog

#         ## another catcher to let vsearch cleanup after clustering is done
#         proc.wait()
#         print(100)


#     except KeyboardInterrupt:
#         LOGGER.info("interrupted vsearch here: %s", proc.pid)
#         os.kill(proc.pid, 2)
#         raise KeyboardInterrupt()
#     except sps.CalledProcessError as inst:
#         raise IPyradWarningExit("""
#         Error in vsearch: \n{}\n{}""".format(inst, sps.STDOUT))
#     except OSError as inst:
#         raise IPyradWarningExit("""
#         Failed to allocate pty: \n{}""".format(inst))

#     finally:
#         data.stats_files.s6 = logfile
#         ## cleanup processes
#         #del proc, dog, owner



def build_h5_array(data, samples, nloci):
    """
    Sets up all of the h5 arrays that we will fill. 
    The catg array of prefiltered loci  is 4-dimensional (Big), so one big 
    array would overload memory, we need to fill it in slices. 
    This will be done in multicat (singlecat) and fill_superseqs.
    """

    # sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    # get maxlen dim
    maxlen = data._hackersonly["max_fragment_length"] + 20
    LOGGER.info("maxlen inside build_h5_array is %s", maxlen)
    LOGGER.info("nloci inside build_h5_array is %s", nloci)

    # open new h5 handle
    data.clust_database = os.path.join(
        data.dirs.across, data.name + ".clust.hdf5")
    io5 = h5py.File(data.clust_database, 'w')

    # chunk to approximately 2 chunks per core
    chunks = ((nloci // (data.cpus * 2)) + (nloci % (data.cpus * 2)))

    # Number of elements in hdf5 chunk may not exceed 500MB
    # This is probably not actually optimal, to have such
    # enormous chunk sizes, could probably explore efficiency
    # of smaller chunk sizes on very very large datasets
    chunklen = chunks * len(samples) * maxlen * 4
    while chunklen > int(500e6):
        chunks = (chunks // 2) + (chunks % 2)
        chunklen = chunks * len(samples) * maxlen * 4
    LOGGER.info("chunks in build_h5_array: %s", chunks)

    data.chunks = chunks
    LOGGER.info("nloci is %s", nloci)
    LOGGER.info("chunks is %s", data.chunks)

    ## INIT FULL CATG ARRAY
    ## store catgs with a .10 loci chunk size. For huuuge speed improvements
    # we should change this to only store catgs for variable sites...
    supercatg = io5.create_dataset(
        "catgs", 
        (nloci, len(samples), maxlen, 4),
        dtype=np.uint32,
        chunks=(chunks, 1, maxlen, 4),
        compression="gzip")
    superseqs = io5.create_dataset(
        "seqs", 
        (nloci, len(samples), maxlen),
        dtype="S1",
        #dtype=np.uint8,
        chunks=(chunks, len(samples), maxlen),
        compression='gzip')
    superalls = io5.create_dataset(
        "nalleles",
        (nloci, len(samples)),
        dtype=np.uint8,
        chunks=(chunks, len(samples)),
        compression="gzip")
    superchroms = io5.create_dataset(
        "chroms", (nloci, 3), 
        dtype=np.int64, 
        chunks=(chunks, 3),
        compression="gzip")

    ## allele count storage
    supercatg.attrs["chunksize"] = (chunks, 1, maxlen, 4)
    supercatg.attrs["samples"] = [i.name.encode() for i in samples]
    superseqs.attrs["chunksize"] = (chunks, len(samples), maxlen)
    superseqs.attrs["samples"] = [i.name.encode() for i in samples]
    superalls.attrs["chunksize"] = (chunks, len(samples))
    superalls.attrs["samples"] = [i.name.encode() for i in samples]
    superchroms.attrs["chunksize"] = (chunks, len(samples))
    superchroms.attrs["samples"] = [i.name.encode() for i in samples]

    ## array for pair splits locations, dup and ind filters
    io5.create_dataset("splits", (nloci, ), dtype=np.uint16)
    io5.create_dataset("duplicates", (nloci, ), dtype=np.bool_)

    ## close the big boy
    io5.close()



def fill_dups_arr(data):
    """
    fills the duplicates array from the multi_muscle_align tmp files
    """
    ## build the duplicates array
    duplefiles = glob.glob(os.path.join(data.tmpdir, "duples_*.tmp.npy"))
    duplefiles.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))

    ## enter the duplicates filter into super h5 array
    io5 = h5py.File(data.clust_database, 'r+')
    dfilter = io5["duplicates"]

    ## enter all duple arrays into full duplicates array
    init = 0
    for dupf in duplefiles:
        end = int(dupf.rsplit("_", 1)[-1][:-8])
        inarr = np.load(dupf)
        dfilter[init:end] = inarr
        init += end - init
        #os.remove(dupf)
    #del inarr

    ## continued progress bar
    LOGGER.info("all duplicates: %s", dfilter[:].sum())
    io5.close()



def build_tmp_h5(data, samples):
    """ build tmp h5 arrays that can return quick access for nloci"""
    ## get samples and names, sorted
    snames = [i.name for i in samples]
    snames.sort()

    ## Build an array for quickly indexing consens reads from catg files.
    ## save as a npy int binary file.
    uhandle = os.path.join(data.dirs.across, data.name + ".utemp.sort")
    bseeds = os.path.join(data.dirs.across, data.name + ".tmparrs.h5")

    ## send as first async1 job
    get_seeds_and_hits(uhandle, bseeds, snames)



def get_nloci(data):
    """ return nloci from the tmp h5 arr"""
    bseeds = os.path.join(data.dirs.across, data.name + ".tmparrs.h5")
    with h5py.File(bseeds) as io5:
        return io5["seedsarr"].shape[0]



def get_seeds_and_hits(uhandle, bseeds, snames):
    """
    builds a seeds and hits (uarr) array of ints from the utemp.sort file.
    Saves outputs to files ...
    """
    ## read in the utemp.sort file
    updf = np.loadtxt(uhandle, dtype=np.bytes_)

    ## Get seeds for all matches from usort
    seeds = np.unique(updf[:, 1])
    seedsarr = np.column_stack([
        np.arange(len(seeds)),
        [i.rsplit(b"_", 1)[0] for i in seeds],
        [i.rsplit(b"_", 1)[1] for i in seeds],
        ])
    seedsarr[:, 1] = [snames.index(i.decode()) for i in seedsarr[:, 1]]
    seedsarr = seedsarr.astype(np.int64)
    LOGGER.info("got a seedsarr %s", seedsarr.shape)

    ## Get matches from usort and create an array for fast entry
    uarr = np.zeros((updf.shape[0], 3), dtype=np.int64)
    idx = -1
    lastloc = None
    for ldx in range(updf.shape[0]):
        tloc = updf[ldx, 1]
        if tloc != lastloc:
            idx += 1
        uarr[ldx, 0] = idx
        lastloc = tloc

    ## create a column with sample index
    uarr[:, 1] = [int(snames.index(i.rsplit(b"_", 1)[0].decode()))
                  for i in updf[:, 0]]
    ## create a column with only consens index for sample
    uarr[:, 2] = [int(i.rsplit(b"_", 1)[1]) for i in updf[:, 0]]
    uarr = uarr.astype(np.int64)
    LOGGER.info("got a uarr %s", uarr.shape)

    ## save as h5 to we can grab by sample slices
    with h5py.File(bseeds, 'w') as io5:
        io5.create_dataset("seedsarr", data=seedsarr, dtype=np.int64)
        io5.create_dataset("uarr", data=uarr, dtype=np.int64)



def new_multicat(data, samples, ipyclient):
    """
    Calls 'singlecat()' for all samples to build index files.
    """

    ## track progress
    LOGGER.info("in the multicat")
    start = time.time()
    printstr = ("indexing clusters   ", "s6")

    ## Build the large h5 array. This will write a new HDF5 file and overwrite
    ## existing data. 
    nloci = get_nloci(data)
    build_h5_array(data, samples, nloci)

    ## parallel client (reserve engine 0 for data entry), if/else here in case
    ## user has only one engine.
    if len(ipyclient) > 1:
        filler = ipyclient.load_balanced_view(targets=[0])
        smallview = ipyclient.load_balanced_view(targets=ipyclient.ids[1::2])
    else:
        filler = ipyclient.load_balanced_view(targets=[0])
        smallview = ipyclient.load_balanced_view(targets=[0])                

    ## First submit a sleeper job as temp_flag for cleanups
    last_sample = 0
    cleanups = {}
    cleanups[last_sample] = filler.apply(time.sleep, 0.0)

    ## fill the duplicates filter array
    rasync = smallview.apply(fill_dups_arr, data)
    while 1:
        data._progressbar(20, 0, start, printstr)
        time.sleep(0.1)
        if rasync.ready():
            break
    if not rasync.successful():
        raise IPyradWarningExit(rasync.result())

    ## Get all existing .tmp.h5 files. If files exist then assume that we are
    ## restarting an interrupted job. We need to check for each one whether it 
    ## has it finished being built, and whether it has been written to the 
    ## large array yet.
    snames = [i.name for i in samples]
    snames.sort()
    smpios = {i: os.path.join(data.dirs.across, i + '.tmp.h5') for i in snames}

    ## send 'singlecat()' jobs to engines
    bseeds = os.path.join(data.dirs.across, data.name + ".tmparrs.h5")
    jobs = {}
    for sample in samples:
        sidx = snames.index(sample.name)
        args = (data, sample, bseeds, sidx, nloci)

        ## Only build it if it doesn't already exist. Singlecat removes
        ## unfinished files if interrupted, so .tmp.h5 should not exist
        ## unless the file is ready to be entered. 
        if not os.path.exists(smpios[sample.name]):
            jobs[sample.name] = smallview.apply(singlecat, *args)

    ## track progress of singlecat jobs and submit writing jobs for finished
    ## singlecat files (.tmp.h5).
    alljobs = len(jobs)
    while 1:
        ## check for finished jobs
        curkeys = list(jobs.keys())
        for key in curkeys:
            rasync = jobs[key]
            if rasync.ready():
                if rasync.successful():
                    ## submit cleanup for finished job
                    args = (data, data.samples[key], snames.index(key))
                    with filler.temp_flags(after=cleanups[last_sample]):
                        cleanups[key] = filler.apply(write_to_fullarr, *args)
                        last_sample = key
                        del jobs[key]
                else:
                    if not rasync.successful():
                        raise IPyradWarningExit(rasync.result())

        ## print progress or break
        data._progressbar(alljobs, alljobs - len(jobs), start, printstr)
        time.sleep(0.1)
        if not jobs:
            break        

    ## add the dask_chroms func for reference data
    if 'reference' in data.paramsdict["assembly_method"]:
        with filler.temp_flags(after=cleanups.values()):
            cleanups['ref'] = filler.apply(dask_chroms, *(data, samples))

    ## ------- print breakline between indexing and writing database ---------
    print("")

    ## track progress of databaseing
    start = time.time()
    printstr = ("building database   ", "s6")
    while 1:
        finished = [i for i in cleanups.values() if i.ready()]
        data._progressbar(len(cleanups), len(finished), start, printstr)
        time.sleep(0.1)
        ## break if one failed, or if finished
        if not all([i.successful() for i in finished]):
            break
        if len(cleanups) == len(finished):
            break

    ## check for errors
    for job in cleanups:
        if cleanups[job].ready():
            if not cleanups[job].successful():
                raise IPyradWarningExit((job, cleanups[job].result()))


## This is where indels are imputed
def singlecat(data, sample, bseeds, sidx, nloci):
    """
    Orders catg data for each sample into the final locus order. This allows
    all of the individual catgs to simply be combined later. They are also in
    the same order as the indels array, so indels are inserted from the indel
    array that is passed in.
    """

    LOGGER.info("in single cat here")
    ## enter ref data?
    isref = 'reference' in data.paramsdict["assembly_method"]

    ## grab seeds and hits info for this sample
    with h5py.File(bseeds, 'r') as io5:
        ## get hits just for this sample and sort them by sample order index
        hits = io5["uarr"][:]
        hits = hits[hits[:, 1] == sidx, :]
        #hits = hits[hits[:, 2].argsort()]
        ## get seeds just for this sample and sort them by sample order index
        seeds = io5["seedsarr"][:]
        seeds = seeds[seeds[:, 1] == sidx, :]
        #seeds = seeds[seeds[:, 2].argsort()]
        full = np.concatenate((seeds, hits))
        full = full[full[:, 0].argsort()]

    ## still using max+20 len limit, rare longer merged reads get trimmed
    ## we need to allow room for indels to be added too
    maxlen = data._hackersonly["max_fragment_length"] + 20

    ## we'll fill a new catg and alleles arr for this sample in locus order,
    ## which is known from seeds and hits
    ocatg = np.zeros((nloci, maxlen, 4), dtype=np.uint32)
    onall = np.zeros(nloci, dtype=np.uint8)
    ochrom = np.zeros((nloci, 3), dtype=np.int64)
    
    ## grab the sample's data and write to ocatg and onall
    if not sample.files.database:
        raise IPyradWarningExit("missing catg file - {}".format(sample.name))

    with h5py.File(sample.files.database, 'r') as io5:
        ## get it and delete it
        catarr = io5["catg"][:]
        tmp = catarr[full[:, 2], :maxlen, :]
        del catarr
        ocatg[full[:, 0], :tmp.shape[1], :] = tmp
        del tmp

        ## get it and delete it
        nall = io5["nalleles"][:]
        onall[full[:, 0]] = nall[full[:, 2]]
        del nall

        ## fill the reference data
        if isref:
            chrom = io5["chroms"][:]
            ochrom[full[:, 0]] = chrom[full[:, 2]]
            del chrom

    ## get indel locations for this sample
    ipath = os.path.join(data.dirs.across, data.name + ".tmp.indels.hdf5")
    with h5py.File(ipath, 'r') as ih5:
        indels = ih5["indels"][sidx, :, :maxlen]

    ## insert indels into ocatg
    newcatg = inserted_indels(indels, ocatg)
    del ocatg, indels
    
    ## save individual tmp h5 data
    smpio = os.path.join(data.dirs.across, sample.name + '.tmp.h5')
    with h5py.File(smpio, 'w') as oh5:
        oh5.create_dataset("icatg", data=newcatg, dtype=np.uint32)
        oh5.create_dataset("inall", data=onall, dtype=np.uint8)
        if isref:
            oh5.create_dataset("ichrom", data=ochrom, dtype=np.int64)



## This func could potentially be replaced entirely by making 
## a dask array made up concatenating all of the individual 
## .tmp.h5 arrays. Reading from that might be a bit slower, tho.
## Something to consider since it would be parallel, while this 
## step is not. If this step proves to be too slow in uber
## big assemblies then we should try the other approach.
def write_to_fullarr(data, sample, sidx):
    """ writes arrays to h5 disk """

    ## enter ref data?
    #isref = 'reference' in data.paramsdict["assembly_method"]
    LOGGER.info("writing fullarr %s %s", sample.name, sidx)

    ## save big arrays to disk temporarily
    with h5py.File(data.clust_database, 'r+') as io5:
        ## open views into the arrays we plan to fill
        chunk = io5["catgs"].attrs["chunksize"][0]
        catg = io5["catgs"]
        nall = io5["nalleles"]

        ## adding an axis to newcatg makes it write about 1000X faster.
        smpio = os.path.join(data.dirs.across, sample.name + '.tmp.h5')
        with h5py.File(smpio) as indat:

            ## grab all of the data from this sample's arrays
            newcatg = indat["icatg"]  # [:]
            onall = indat["inall"]    # [:]

            ## enter it into the full array one chunk at a time
            for cidx in range(0, catg.shape[0], chunk):
                end = cidx + chunk
                catg[cidx:end, sidx:sidx + 1, :] = np.expand_dims(newcatg[cidx:end, :], axis=1)
                nall[:, sidx:sidx + 1] = np.expand_dims(onall, axis=1)



def dask_chroms(data, samples):
    """
    A dask relay function to fill chroms for all samples
    """
    
    ## example concatenating with dask
    h5s = [os.path.join(data.dirs.across, s.name + ".tmp.h5") for s in samples]
    handles = [h5py.File(i) for i in h5s]
    dsets = [i['/ichrom'] for i in handles]
    arrays = [da.from_array(dset, chunks=(10000, 3)) for dset in dsets]
    stack = da.stack(arrays, axis=2)

    ## max chrom (should we check for variable hits? if so, things can get wonk
    maxchrom = da.max(stack, axis=2)[:, 0]

    ## max pos
    maxpos = da.max(stack, axis=2)[:, 2]

    ## min pos
    mask = stack == 0
    stack[mask] = 9223372036854775807  # max int64 value
    minpos = da.min(stack, axis=2)[:, 1]
    final = da.stack([maxchrom, minpos, maxpos], axis=1)
    final.to_hdf5(data.clust_database, "/chroms")

    ## close the h5 handles
    [i.close() for i in handles]

#max int64 = 9223372036854775807
#max uint64 = 18446744073709551615
#max32 = 4294967295


@numba.jit(nopython=True)
def inserted_indels(indels, ocatg):
    """
    inserts indels into the catg array
    """
    ## return copy with indels inserted
    newcatg = np.zeros(ocatg.shape, dtype=np.uint32)

    ## iterate over loci and make extensions for indels
    for iloc in range(ocatg.shape[0]):
        ## get indels indices
        indidx = np.where(indels[iloc, :])[0]
        if np.any(indidx):
            ## which new (empty) rows will be added
            allrows = np.arange(ocatg.shape[1])
            mask = np.ones(allrows.shape[0], dtype=np.bool_)
            for idx in indidx:
                mask[idx] = False
            not_idx = allrows[mask == 1]

            ## fill in new data into all other spots
            newcatg[iloc][not_idx] = ocatg[iloc, :not_idx.shape[0]]
        else:
            newcatg[iloc] = ocatg[iloc]
    return newcatg



def fill_superseqs(data, samples):
    """
    Fills the superseqs array with seq data from cat.clust
    and fill the edges array with information about paired split locations.
    """

    ## load super to get edges
    io5 = h5py.File(data.clust_database, 'r+')
    superseqs = io5["seqs"]
    splits = io5["splits"]

    ## samples are already sorted
    snames = [i.name for i in samples]
    LOGGER.info("snames %s", snames)

    ## get maxlen again
    maxlen = data._hackersonly["max_fragment_length"] + 20
    LOGGER.info("maxlen inside fill_superseqs is %s", maxlen)

    ## data has to be entered in blocks
    infile = os.path.join(data.dirs.across, data.name + "_catclust.gz")
    clusters = gzip.open(infile, 'r')
    pairdealer = izip(*[iter(clusters)] * 2)

    ## iterate over clusters
    chunks = superseqs.attrs["chunksize"]
    chunksize = chunks[0]
    done = 0
    iloc = 0
    cloc = 0
    chunkseqs = np.zeros(chunks, dtype="S1")
    chunkedge = np.zeros(chunksize, dtype=np.uint16)

    while 1:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradWarningExit("clustfile formatting error in %s", chunk)

        ## if chunk is full put into superseqs and reset counter
        if cloc == chunksize:
            LOGGER.info("cloc chunk writing %s", cloc)
            superseqs[iloc - cloc:iloc] = chunkseqs
            splits[iloc - cloc:iloc] = chunkedge
            ## reset chunkseqs, chunkedge, cloc
            cloc = 0
            chunkseqs = np.zeros((chunksize, len(samples), maxlen), dtype="S1")
            chunkedge = np.zeros((chunksize), dtype=np.uint16)

        ## get seq and split it
        if chunk:
            try:
                fill = np.zeros((len(samples), maxlen), dtype="S1")
                fill.fill("N")
                piece = chunk[0].decode().strip().split("\n")
                names = piece[0::2]
                seqs = np.array([list(i) for i in piece[1::2]])
                
                ## fill in the separator if it exists
                separator = np.where(np.all(seqs == 'n', axis=0))[0]
                if np.any(separator):
                    chunkedge[cloc] = separator.min()

                # fill in the hits
                # seqs will be (5,) IF the seqs are variable lengths, which 
                # can happen if it had duplicaes AND there were indels, and 
                # so the indels did not get aligned
                try:
                    shlen = seqs.shape[1]
                except IndexError as inst:
                    shlen = min([len(x) for x in seqs])

                for name, seq in zip(names, seqs):
                    sidx = snames.index(name.rsplit("_", 1)[0])
                    #fill[sidx, :shlen] = seq[:maxlen]
                    fill[sidx, :shlen] = seq[:shlen]

                ## PUT seqs INTO local ARRAY
                chunkseqs[cloc] = fill

            except Exception as inst:
                LOGGER.info(inst)
                LOGGER.info("\nfill: %s\nshlen %s\nmaxlen %s", fill.shape, shlen, maxlen)
                LOGGER.info("dupe chunk \n{}".format("\n".join(chunk)))

            ## increase counters if there was a chunk
            cloc += 1
            iloc += 1
        if done:
            break

    ## write final leftover chunk
    superseqs[iloc - cloc:, ] = chunkseqs[:cloc]
    splits[iloc - cloc:] = chunkedge[:cloc]

    ## close super
    io5.close()
    clusters.close()

    ## edges is filled with splits for paired data.
    LOGGER.info("done filling superseqs")

    ## close handle
    #os.remove(infile)




def sub_build_clustbits(data, usort, nseeds):
    """
    A subfunction of build_clustbits to allow progress tracking. This func
    splits the unaligned clusters into bits for aligning on separate cores.
    """

    ## load FULL concat fasta file into a dict. This could cause RAM issues.
    ## this file has iupac codes in it, not ambigs resolved, and is gzipped.
    LOGGER.info("loading full _catcons file into memory")
    allcons = {}
    conshandle = os.path.join(data.dirs.across, data.name + "_catcons.tmp")
    with gzip.open(conshandle, 'rt') as iocons:
        cons = izip(*[iter(iocons)] * 2)
        for namestr, seq in cons:
            nnn, sss = [i.strip() for i in (namestr, seq)]
            allcons[nnn[1:]] = sss

    ## set optim to approximately 4 chunks per core. Smaller allows for a bit
    ## cleaner looking progress bar. 40 cores will make 160 files.
    optim = ((nseeds // (data.cpus * 4)) + (nseeds % (data.cpus * 4)))
    LOGGER.info(
        "building clustbits, optim=%s, nseeds=%s, cpus=%s",
        optim, nseeds, data.cpus)

    ## iterate through usort grabbing seeds and matches
    with open(usort, 'rt') as insort:
        ## iterator, seed null, and seqlist null
        isort = iter(insort)
        loci = 0
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0

        while 1:
            ## grab the next line
            try:
                hit, seed, ori = next(isort).strip().split()
            except StopIteration:
                break

            try:
                ## if same seed, append match
                if seed != lastseed:
                    ## store the last fseq, count it, and clear it
                    if fseqs:
                        seqlist.append("\n".join(fseqs))
                        seqsize += 1
                        fseqs = []
                    ## occasionally write to file
                    if seqsize >= optim:
                        if seqlist:
                            loci += seqsize
                            pathname = os.path.join(data.tmpdir,
                                data.name + ".chunk_{}".format(loci))
                            with open(pathname, 'wt') as clustsout:
                                clustsout.write(
                                    "\n//\n//\n".join(seqlist) + "\n//\n//\n")
                                LOGGER.debug(
                                    "writing chunk - seqsize {} loci {} {}"
                                    .format(seqsize, loci, clustsout.name))

                            ## reset list and counter
                            seqlist = []
                            seqsize = 0
    
                    ## store the new seed on top of fseq
                    fseqs.append(">{}\n{}".format(seed, allcons[seed]))
                    lastseed = seed
    
                ## add match to the seed
                seq = allcons[hit]
                ## revcomp if orientation is reversed
                if ori == "-":
                    seq = fullcomp(seq)[::-1]
                fseqs.append(">{}\n{}".format(hit, seq))
            except KeyError as inst:
                ## Caught bad seed or hit? Log and continue.
                LOGGER.error(
                    "Bad Seed/Hit: seqsize {}\tloci {}\tseed {}\thit {}"
                    .format(seqsize, loci, seed, hit))

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
        seqsize += 1
        loci += seqsize
    if seqlist:
        pathname = os.path.join(data.tmpdir, 
            data.name + ".chunk_{}".format(loci))
        with open(pathname, 'wt') as clustsout:
            clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")

    ## final progress and cleanup
    del allcons
    clustbits = glob.glob(os.path.join(data.tmpdir, data.name + ".chunk_*"))

    ## return stuff
    return clustbits, loci



def build_input_file(data, samples, randomseed):
    """
    [This is run on an ipengine]
    Make a concatenated consens file with sampled alleles (no RSWYMK/rswymk).
    Orders reads by length and shuffles randomly within length classes
    """

    ## get all of the consens handles for samples that have consens reads
    ## this is better than using sample.files.consens for selecting files
    ## b/c if they were moved we only have to edit data.dirs.consens

    ## scratch the statement above, people shouldn't be moving files, 
    ## they should be using merge/branch, and so sample.files.consens
    ## is needed to keep track of samples from different dirs if they
    ## are later merged into the same assembly.
    #conshandles = [os.path.join(data.dirs.consens, sample.name+".consens.gz") 
    #              for sample in samples if \
    #              sample.stats.reads_consens]
    conshandles = [
        sample.files.consens[0] for sample in samples if 
        sample.stats.reads_consens]
    conshandles.sort()
    assert conshandles, "no consensus files found"

    ## concatenate all of the gzipped consens files
    cmd = ['cat'] + conshandles
    #allcons = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    allcons = os.path.join(data.dirs.across, data.name + "_catcons.tmp")
    LOGGER.debug(" ".join(cmd))
    with open(allcons, 'w') as output:
        call = sps.Popen(cmd, stdout=output, close_fds=True)
        call.communicate()

    ## a string of sed substitutions for temporarily replacing hetero sites
    ## skips lines with '>', so it doesn't affect taxon names
    subs = ["/>/!s/W/A/g", "/>/!s/w/A/g", "/>/!s/R/A/g", "/>/!s/r/A/g",
            "/>/!s/M/A/g", "/>/!s/m/A/g", "/>/!s/K/T/g", "/>/!s/k/T/g",
            "/>/!s/S/C/g", "/>/!s/s/C/g", "/>/!s/Y/C/g", "/>/!s/y/C/g"]
    subs = ";".join(subs)

    ## impute pseudo-haplo information to avoid mismatch at hetero sites
    ## the read data with hetero sites is put back into clustered data later.
    ## pipe passed data from gunzip to sed.
    cmd1 = ["gunzip", "-c", allcons]
    cmd2 = ["sed", subs]
    LOGGER.debug(" ".join(cmd1))
    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, close_fds=True)
    allhaps = allcons.replace("_catcons.tmp", "_cathaps.tmp")
    with open(allhaps, 'w') as output:
        LOGGER.debug(" ".join(cmd2))
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=output, close_fds=True)
        proc2.communicate()
    proc1.stdout.close()

    ## now sort the file using vsearch
    allsort = allcons.replace("_catcons.tmp", "_catsort.tmp")
    cmd1 = [ipyrad.bins.vsearch,
            "--sortbylength", allhaps,
            "--fasta_width", "0",
            "--output", allsort]
    LOGGER.debug(" ".join(cmd1))
    proc1 = sps.Popen(cmd1, close_fds=True)
    proc1.communicate()

    ## shuffle sequences within size classes. Tested seed (8/31/2016)
    ## shuffling works repeatably with seed.
    random.seed(randomseed)

    ## open an iterator to lengthsorted file and grab two lines at at time
    allshuf = allcons.replace("_catcons.tmp", "_catshuf.tmp")
    outdat = open(allshuf, 'wt')
    indat = open(allsort, 'r')
    idat = izip(iter(indat), iter(indat))
    done = 0

    chunk = [next(idat)]
    while not done:
        ## grab 2-lines until they become shorter (unless there's only one)
        oldlen = len(chunk[-1][-1])
        while 1:
            try:
                dat = next(idat)
            except StopIteration:
                done = 1
                break
            if len(dat[-1]) == oldlen:
                chunk.append(dat)
            else:
                ## send the last chunk off to be processed
                random.shuffle(chunk)
                outdat.write("".join(chain(*chunk)))
                ## start new chunk
                chunk = [dat]
                break

    ## do the last chunk
    random.shuffle(chunk)
    outdat.write("".join(chain(*chunk)))

    indat.close()
    outdat.close()



def clean_and_build_concat(data, samples, randomseed, ipyclient):
    """ 
    STEP 6-1:
    Clears dirs and databases and calls 'build_input_file()'
    """
    ## but check for new clust database name if this is a new branch
    cleanup_tempfiles(data)
    catclust = os.path.join(data.dirs.across, data.name + "_catclust.gz")
    if os.path.exists(catclust):
        os.remove(catclust)
    if os.path.exists(data.clust_database):
        os.remove(data.clust_database)

    ## get parallel view
    start = time.time()
    printstr = ("concat/shuffle input", "s6")

    ## make a vsearch input fasta file with all samples reads concat
    rasync = ipyclient[0].apply(build_input_file, *[data, samples, randomseed])
    while 1:
        ready = int(rasync.ready())
        data._progressbar(1, ready, start, printstr)
        if ready:
            break
        else:
            time.sleep(0.1)
    print("")

    ## store that this step was successful
    if not rasync.successful():
        raise IPyradWarningExit(rasync.result())



def run(data, samples, noreverse, force, randomseed, ipyclient, **kwargs):
    """
    For step 6 the run function is sub divided a bit so that users with really
    difficult assemblies can possibly interrupt and restart the step from a 
    checkpoint. 

    Substeps that are run:
     1. build concat consens file, 
     2. cluster all consens,
     3. split clusters into bits,
     4. align bits, 
     5. build indel array
     6. build h5 array.
     7. Enter seq data & cleanup
    """

    ## if force then set checkpoint to zero and run all substeps for just
    ## the user specified steps. 
    if force:
        data._checkpoint = 0
        if kwargs.get('substeps'):
            substeps = kwargs.get('substeps')
        else:  
            substeps = range(1, 8)

    ## if {data}._checkpoint attribute exists then find the checkpoint where
    ## this assembly left off (unless force) and build step list from there.
    else:
        if kwargs.get('substeps'):
            substeps = kwargs.get('substeps')
        else:  
            if hasattr(data, '_checkpoint'):
                substeps = range(max(1, data._checkpoint), 8)
            else:
                data._checkpoint = 0
                substeps = range(1, 8)

    ## build substeps list to subset which funtions need to be run
    if isinstance(substeps, (int, float, str)):
        substeps = [substeps]
        substeps = [int(i) for i in substeps]

    ## print continuation message
    if substeps[0] != 1:
        print("{}Continuing from checkpoint 6.{}"
              .format(data._spacer, substeps[0]))
    LOGGER.info("checkpoint = %s", data._checkpoint)
    LOGGER.info("substeps = %s", substeps)

    ## Set variables on data that are needed for all steps;
    data.dirs.across = os.path.realpath(
        os.path.join(data.paramsdict["project_dir"], data.name + "_across"))
    data.tmpdir = os.path.join(data.dirs.across, data.name + "-tmpalign")
    data.clust_database = os.path.join(
        data.dirs.across, data.name + ".clust.hdf5")
    if not os.path.exists(data.dirs.across):
        os.mkdir(data.dirs.across)
    if not os.path.exists(data.tmpdir):
        os.mkdir(data.tmpdir)
    data.cpus = data._ipcluster["cores"]
    if not data.cpus:
        data.cpus = len(ipyclient)

    ## STEP 6-1: Clean database and build input concat file for clustering
    if 1 in substeps:
        clean_and_build_concat(data, samples, randomseed, ipyclient)
        data._checkpoint = 1
        data.save()

    ## STEP 6-2: Cluster across w/ vsearch; uses all threads on largest host 
    if 2 in substeps:
        call_cluster(data, noreverse, ipyclient)
        data._checkpoint = 2
        data.save()

    ## builds consens cluster bits and writes them to the tmp directory. These
    ## will not be deleted until either step 6-6 is complete, or the force flag
    ## is used. This will clear the tmpdir if it is run.
    if 3 in substeps:
        build_clustbits(data, ipyclient, force)
        data._checkpoint = 3
        data.save()

    ## muscle align the cluster bits and create tmp hdf5 indel arrays for the
    ## next step. These will not be deleted until...
    if 4 in substeps:
        multi_muscle_align(data, samples, ipyclient)
        data._checkpoint = 4
        data.save()

    ## fill the indel array with the indel tmp arrays from aligning step.
    if 5 in substeps:
        build_indels(data, samples, ipyclient)
        data._checkpoint = 5
        data.save()

    if 6 in substeps:
        ## builds the final HDF5 array which includes three main keys
        ## /catg -- contains all indiv catgs and has indels inserted
        ##   .attr['samples'] = [samples]
        ## /filters -- filled for dups, left empty for others until step 7.
        ##   .attr['filters'] = [f1, f2, f3, f4, f5]
        ## /seqs -- contains the clustered sequence data as string arrays
        ##   .attr['samples'] = [samples]
        ## /edges -- gets the paired split locations for now.
        ## /snps  -- left empty for now

        ## FILL SUPERCATG and fills dupfilter, indfilter, and nalleles
        ## this function calls singlecat() on each sample and enters their
        ## resulting arrays into the superarray. If all singlecats are built
        ## then it will continue to enter them into the database. 
        LOGGER.info("multicat -- building full database")
        new_multicat(data, samples, ipyclient)
        data._checkpoint = 6
        data.save()

    if 7 in substeps:
        ## FILL SUPERSEQS and fills edges(splits) for paired-end data
        fill_superseqs(data, samples)
        data._checkpoint = 7
        data.save()

        ## remove files but not dir (used in step 1 too)
        cleanup_tempfiles(data)
        ## remove the tmpdir
        if os.path.exists(data.tmpdir):
            shutil.rmtree(data.tmpdir)

        ## set sample states
        for sample in samples:
            sample.stats.state = 6
        print("")



def cleanup_tempfiles(data):
    """ 
    Function to remove older files. This is called either in substep 1 or after
    the final substep so that tempfiles are retained for restarting interrupted
    jobs until we're sure they're no longer needed. 
    """

    ## remove align-related tmp files
    tmps1 = glob.glob(os.path.join(data.tmpdir, "*.fa"))
    tmps2 = glob.glob(os.path.join(data.tmpdir, "*.npy"))
    for tmp in tmps1 + tmps2:
        if os.path.exists(tmp):
            os.remove(tmp)

    ## remove cluster related files
    removal = [
        os.path.join(data.dirs.across, data.name + ".utemp"),
        os.path.join(data.dirs.across, data.name + ".htemp"),
        os.path.join(data.dirs.across, data.name + "_catcons.tmp"),
        os.path.join(data.dirs.across, data.name + "_cathaps.tmp"),
        os.path.join(data.dirs.across, data.name + "_catshuf.tmp"),
        os.path.join(data.dirs.across, data.name + "_catsort.tmp"),
        os.path.join(data.dirs.across, data.name + ".tmparrs.h5"),
        os.path.join(data.dirs.across, data.name + ".tmp.indels.hdf5"),
        ]
    for rfile in removal:
        if os.path.exists(rfile):
            os.remove(rfile)

    ## remove singlecat related h5 files
    smpios = glob.glob(os.path.join(data.dirs.across, '*.tmp.h5'))
    for smpio in smpios:
        if os.path.exists(smpio):
            os.remove(smpio)            