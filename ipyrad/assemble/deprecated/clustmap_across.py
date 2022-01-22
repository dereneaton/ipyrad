n#!/usr/bin/env python

"""
Denovo: Cluster across samples w/ vsearch and build clusters.
Reference: Extract mapping from bam files and bedtools to build clusters
  and apply filters to identify bad mapping.
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
import shutil
import random
import select
import subprocess as sps

from loguru import logger
import numpy as np
from pysam import AlignmentFile, FastaFile
import ipyrad
from ipyrad.assemble.utils import IPyradError, fullcomp, chroms2ints
from ipyrad.assemble.utils import AssemblyProgressBar, VsearchProgressBar



class Step6:
    """
    Group loci across samples by orthology using clustering in denovo
    or by pulling from map positions in reference.
    """
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.randomseed = int(self.data.hackersonly.random_seed)
        self.isref = bool('ref' in self.data.params.assembly_method)
        self.force = force
        self.ipyclient = ipyclient
        self.print_headers()
        self.samples = self.get_subsamples()
        self.setup_dirs(force)

        # groups/threading information for hierarchical clustering
        # ----- DISABLED FOR NOW -------------
        # self.cgroups = {}
        # self.assign_groups()
        # self.hostd = {}
        # self.tune_hierarchical_threading()

        # NEW CODE TO OVERRIDE HIERARCH CLUSTERING
        self.cgroups = {
            0: self.samples,
        }
        self.data.ncpus = len(self.ipyclient.ids)
        self.nthreads = len(self.ipyclient.ids)
        self.lbview = self.ipyclient.load_balanced_view()
        self.thview = self.ipyclient.load_balanced_view()

        # beds 
        self.regions = []


    def print_headers(self):
        """
        Print CLI header 
        """
        if self.data._cli:
            self.data._print(
                "\n{}Step 6: Clustering/Mapping across samples "
                .format(self.data._spacer)
            )


    def setup_dirs(self, force=False):
        """
        set up across and tmpalign dirs and init h5 database file
        """
        self.data.dirs.across = os.path.realpath(os.path.join(
            self.data.params.project_dir,
            "{}_across".format(self.data.name)))
        self.data.tmpdir = os.path.join(
            self.data.dirs.across,
            "{}-tmpalign".format(self.data.name))
        self.data.clust_database = os.path.join(
            self.data.dirs.across, 
            self.data.name + "_clust_database.fa")

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


    def get_subsamples(self):
        """
        Apply state, ncluster, and force filters to select samples
        """
        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 6")

        # filter samples by state
        state4 = self.data.stats.index[self.data.stats.state < 5]
        state5 = self.data.stats.index[self.data.stats.state == 5]
        state6 = self.data.stats.index[self.data.stats.state > 5]

        # tell user which samples are not ready for step5
        if state4.any():
            print("skipping samples not in state==5:\n{}"
                  .format(state4.tolist()))

        if self.force:
            # run all samples above state 4
            subs = self.data.stats.index[self.data.stats.state > 4]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 6
            if state6.any():
                raise IPyradError(
                    "Some samples are already in state==6. If you wish to "
                    "create a new database for across sample comparisons "
                    "use the force=True (-f) argument.")
            # run all samples in state 5
            subsamples = [self.data.samples[i] for i in state5]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_consens:
                checked_samples.append(sample)
            else:
                print("skipping {}; no consensus reads found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 6")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_consens,
            reverse=True,
        )
        return checked_samples        


    # def assign_groups(self):
    #     "assign samples to groups if not user provided for hierarchical clust"

    #     # to hold group ints mapping to list of sample objects
    #     # {0: [a, b, c], 1: [d, e, f]}
    #     self.cgroups = {}

    #     # use population info to split samples into groups; or assign random
    #     if self.data.populations:
    #         self.cgroups = {}
    #         for idx, val in enumerate(self.data.populations.values()):
    #             self.cgroups[idx] = [self.data.samples[x] for x in val[1]]

    #     # by default let's split taxa into groups of 20-50 samples at a time
    #     else:
    #         # calculate the number of cluster1 jobs to perform:
    #         if len(self.samples) <= 100:
    #             groupsize = 20
    #         elif len(self.samples) <= 500:
    #             groupsize = 50
    #         else:
    #             groupsize = 100

    #         # split samples evenly into groups
    #         alls = self.samples
    #         nalls = len(self.samples)
    #         ngroups = int(np.ceil(nalls / groupsize))
    #         gsize = int(np.ceil(nalls / ngroups))

    #         idx = 0
    #         for samps in range(0, nalls, gsize):
    #             self.cgroups[idx] = alls[samps: samps + gsize]
    #             idx += 1


    # def tune_hierarchical_threading(self):
    #     "tune threads for across-sample clustering used in denovo assemblies"

    #     # get engine data, skips busy engines.
    #     hosts = {}
    #     for eid in self.ipyclient.ids:
    #         engine = self.ipyclient[eid]
    #         if not engine.outstanding:
    #             hosts[eid] = engine.apply(socket.gethostname)

    #     # get targets on each hostname for spreading jobs out.
    #     self.ipyclient.wait()
    #     hosts = [(eid, i.get()) for (eid, i) in hosts.items()]
    #     hostnames = set([i[1] for i in hosts])
    #     self.hostd = {x: [i[0] for i in hosts if i[1] in x] for x in hostnames}

    #     # calculate the theading of cluster1 jobs:
    #     self.data.ncpus = len(self.ipyclient.ids)
    #     njobs = len(self.cgroups)
    #     nnodes = len(self.hostd)

    #     # how to load-balance cluster2 jobs
    #     # maxthreads = 8 cuz vsearch isn't v efficient above that.
    #     ## e.g., 24 cpus; do 2 12-threaded jobs
    #     ## e.g., 2 nodes; 40 cpus; do 2 20-threaded jobs or 4 10-threaded jobs
    #     ## e.g., 4 nodes; 80 cpus; do 8 10-threaded jobs
    #     if nnodes == 1:
    #         thr = np.floor(self.data.ncpus / njobs).astype(int)
    #         eids = max(1, thr)
    #         eids = max(eids, len(list(self.hostd.values())[0]))

    #     else:
    #         eids = []
    #         for node in self.hostd:
    #             sids = self.hostd[node]
    #             nids = len(sids)
    #             thr = np.floor(nids / (njobs / nnodes)).astype(int)
    #             thr = max(1, thr)
    #             thr = min(thr, nids)
    #             eids.extend(self.hostd[node][::thr])

    #     # set nthreads based on ipcluster dict (default is 2)        
    #     #if "threads" in self.data.ipcluster.keys():
    #     #    self.nthreads = int(self.data.ipcluster["threads"])
    #     self.nthreads = 2
    #     if self.data.ncpus > 4:
    #         self.nthreads = int(np.floor(
    #             self.data.ncpus) / len(self.cgroups))
    #     eids = self.ipyclient.ids[::self.nthreads]

    #     # create load-balancers
    #     self.lbview = self.ipyclient.load_balanced_view()
    #     self.thview = self.ipyclient.load_balanced_view(targets=eids)


    def run(self):
        """
        Runs the set of methods for denovo or reference method
        """
        # DENOVO
        if self.data.params.assembly_method == "denovo":

            # prepare clustering inputs for hierarchical clustering
            self.remote_build_concats_tier1()

            # if multiple clusters:
            if len(self.cgroups.keys()) == 1:
                self.remote_cluster_tiers(0)

            else:
                # send initial clustering jobs (track finished jobs)
                self.remote_cluster1()

                # prepare second tier inputs
                self.remote_build_concats_tier2()

                # send cluster2 job (track actual progress)
                self.remote_cluster_tiers('x')

            # build clusters
            self.remote_build_denovo_clusters()

            # align denovo clusters
            self.remote_align_denovo_clusters()

            # concat aligned files
            self.concat_alignments()


        elif self.data.params.assembly_method == "reference":

            # prepare bamfiles (merge and sort)
            self.remote_concat_bams()

            # get extents of regions using bedtools merge
            self.remote_build_ref_regions()

            # build clusters from regions
            self.remote_build_ref_clusters()

            # concat aligned files (This is not necessary, chunk again in s7)
            self.concat_alignments()

        # clean up step here...
        self.data.stats_files.s6 = self.data.clust_database

        # set sample states
        for sample in self.samples:
            sample.stats.state = 6


    def remote_build_concats_tier1(self):
        """
        Prepares concatenated consens input files for each clust1 group
        """
        printstr = ("concatenating inputs", "s6")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasyncs = {}
        for jobid, group in self.cgroups.items():
            # should we use sample objects or sample names in cgroups?
            # Well you gotta choose one! W/o pops file it uses sample objects
            # so I made it use sample objects if pop_assign_file is set iao
            samples = [i for i in self.samples if i in group]
            args = (self.data, jobid, samples, self.randomseed)
            rasyncs[jobid] = self.lbview.apply(build_concat_files, *args)

        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_cluster1(self):
        """
        send threaded jobs to remote engines
        """
        printstr = ("clustering tier 1   ", "s6")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasyncs = {}
        for jobid in self.cgroups:
            args = (self.data, jobid, self.nthreads)
            rasyncs[jobid] = self.thview.apply(cluster, *args)

        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_build_concats_tier2(self):
        """
        Hierarchical clustering build
        """
        printstr = ("concatenating inputs", "s6")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        args = (self.data, list(self.cgroups.keys()), self.randomseed)
        rasync = self.lbview.apply(build_concat_two, *args)

        prog.jobs = rasync
        prog.block()
        prog.check()


    def remote_cluster_tiers(self, jobid):
        """
        Hierarchical cluster. Uses vsearch output to track progress.
        """
        printstr = ("clustering across   ", "s6")

        # nthreads=0 defaults to using all cores
        args = (self.data, jobid, 0, True)
        rasync = self.thview.apply(cluster, *args)
        prog = VsearchProgressBar({0: rasync}, None, printstr, self.data)
        prog.block()
        prog.check()



    def remote_build_denovo_clusters(self):
        """
        build denovo clusters from vsearch clustered seeds
        """
        # filehandles; if not multiple tiers then 'x' is jobid 0
        uhandle = os.path.join(
            self.data.dirs.across, 
            "{}-x.utemp".format(self.data.name))
        buildfunc = build_hierarchical_denovo_clusters
        if not os.path.exists(uhandle):
            uhandle = uhandle.replace("-x.utemp", "-0.utemp")
            buildfunc = build_single_denovo_clusters
        usort = uhandle + ".sort"

        # sort utemp files, count seeds.
        start = time.time()
        printstr = ("building clusters   ", "s6")
        async1 = self.lbview.apply(sort_seeds, uhandle)
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
        nseeds = async2.get()

        # send the clust bit building job to work and track progress
        async3 = self.lbview.apply(
            buildfunc, *(self.data, usort, nseeds, list(self.cgroups.keys())))
        while 1:
            ready = [async1.ready(), async2.ready(), async3.ready()]
            self.data._progressbar(3, sum(ready), start, printstr)
            time.sleep(0.1)
            if all(ready):
                break
        self.data._print("")

        # check for errors
        for job in [async1, async2, async3]:
            if not job.successful():
                job.get()



    def remote_align_denovo_clusters(self):
        """
        Distributes parallel jobs to align_to_array() function. 
        """
        # submit jobs to engines
        printstr = ("aligning clusters   ", "s6")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        # get files
        globpath = os.path.join(self.data.tmpdir, self.data.name + ".chunk_*")
        clustbits = glob.glob(globpath)
        jobs = {}
        for idx, _ in enumerate(clustbits):
            args = [self.data, self.samples, clustbits[idx]]
            jobs[idx] = self.lbview.apply(align_to_array, *args)

        # print progress while bits are aligning
        prog.jobs = jobs
        prog.block()
        prog.check()



    def concat_alignments(self):
        """
        This step is not necessary... we just chunk it up again in step 7...
        it's nice having a file as a product, but why bother...
        It creates a header with names of all samples that were present when
        step 6 was completed. 
        """
        # get files
        globlist = glob.glob(os.path.join(self.data.tmpdir, "aligned_*.fa"))
        clustbits = sorted(
            globlist,
            key=lambda x: int(x.rsplit("_", 1)[1].split(".")[0]),
        )

        # store path to clust database 
        self.data.clust_database = os.path.join(
            self.data.dirs.across, 
            self.data.name + "_clust_database.fa")

        # TODO: count nsnps and save it to the JSON for step 7


        # TODO: use cat to concatenate chunks


        # TODO: with cat be sure empty chunks don't cause problems.


        # write clusters to file with a header that has all samples in db        
        snames = sorted([i.name for i in self.samples])
        with open(self.data.clust_database, 'wt') as out:
            out.write("#{}\n".format(",@".join(snames)))
            for clustfile in clustbits:
                with open(clustfile, 'r') as indata:
                    dat = indata.read()
                    if dat:
                        out.write(dat)  # + "//\n//\n")

        # final cleanup
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)


    ## REFERENCE BASED FUNCTIONS ---------------------------------


    def remote_concat_bams(self):
        """
        Merge bam files into a single large sorted indexed bam
        """
        start = time.time()
        printstr = ("concatenating bams  ", "s6")

        catbam = os.path.join(
            self.data.dirs.across, 
            "{}.cat.bam".format(self.data.name)
        )

        # concatenate consens bamfiles for all samples in this assembly
        cmd1 = [
            ipyrad.bins.samtools,
            "merge", 
            "-f", 
            catbam,
        ]

        # Use the sample.files.consens info, rather than data.dirs to allow
        # for merging assemblies after step 5 where data.dirs is invalid/empty.
        for sample in self.samples:
            cmd1.append(sample.files.consens)

        proc = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)

        # progress bar
        while proc.poll() != 0:
            self.data._progressbar(3, 0, start, printstr)
            time.sleep(0.1)

        # parse result
        err = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(
                "error in: {}: {}".format(" ".join(cmd1), err))

        # sort the bam file
        cmd2 = [
            ipyrad.bins.samtools,
            "sort",
            "-T",
            catbam + '.tmp',
            "-o", 
            os.path.join(
                self.data.dirs.across, 
                "{}.cat.sorted.bam".format(self.data.name)
            ),
            catbam,
        ]
        proc = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)

        # progress bar
        while not proc.poll() == 0:
            self.data._progressbar(3, 1, start, printstr)
            time.sleep(0.1)

        # parse result
        err = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(
                "error in: {}: {}".format(" ".join(cmd2), err))
        os.remove(catbam)

        try:
            # index the bam file
            cmd3 = [
                ipyrad.bins.samtools,
                "index",
                "-c",
                os.path.join(
                    self.data.dirs.across,
                    "{}.cat.sorted.bam".format(self.data.name)
                ),
            ]
            proc = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE)

            # progress bar
            while not proc.poll() == 0:
                self.data._progressbar(3, 2, start, printstr)
                time.sleep(0.1)

            # parse result
            err = proc.communicate()[0].decode()
            if proc.returncode:
                raise IPyradError(
                    "error in: {}: {}".format(" ".join(cmd3), err))
        except IPyradError as ipyerror:
            # For bam files with large chromosomes (>~500Mb) the .bai indexing
            # will fail with this exit message. Try again with .csi indexing.
            # https://github.com/dereneaton/ipyrad/issues/435
            # Will keep bai as default because this has never come up before
            # but it doesn't hurt as a fallback.
            if not "hts_idx_check_range" in str(ipyerror):
                raise ipyerror

            # index the bam file
            cmd3 = [
                ipyrad.bins.samtools,
                "index", "-c",
                os.path.join(
                    self.data.dirs.across,
                    "{}.cat.sorted.bam".format(self.data.name)
                ),
            ]
            proc = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE)

            # progress bar
            while not proc.poll() == 0:
                self.data._progressbar(3, 2, start, printstr)
                time.sleep(0.1)

            # parse result
            err = proc.communicate()[0].decode()
            if proc.returncode:
                raise IPyradError(
                    "error in: {}: {}".format(" ".join(cmd3), err))

        self.data._progressbar(3, 3, start, printstr)
        self.data._print("")



    def remote_build_ref_regions(self):
        """
        call bedtools remotely and track progress
        """
        printstr = ("fetching regions    ", "s6")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()
        rasync = self.ipyclient[0].apply(build_ref_regions, self.data)
        prog.jobs = {1: rasync}        
        prog.block()
        prog.check()
        self.regions = rasync.get()
        logger.debug('regions: {}...'.format(self.regions[:10]))


    def remote_build_ref_clusters(self):
        """
        build clusters and find variants/indels to store
        """       
        # send N jobs each taking chunk of regions
        ncpus = self.data.ncpus
        nloci = len(self.regions)
        optim = int((nloci // ncpus) + (nloci % ncpus))
        optim = int(np.ceil(optim / 2))

        # send jobs to func
        printstr = ("building database   ", "s6")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()
        prog.jobs = {}
        for idx, chunk in enumerate(range(0, nloci, optim)):
            region = self.regions[chunk: chunk + optim]
            if region:
                args = (self.data, idx, region)
                prog.jobs[idx] = self.lbview.apply(build_ref_clusters, *args)

        # print progress while bits are aligning
        prog.block()
        prog.check()



def resolve_duplicates(keys, arr):
    """
    Tries to join together duplicate consens reads that were not previously
    collapsed, likely because there was no overlap of the sequences for one 
    or more samples, but there was for others. Joins two consens reads if the
    """
    newkeys = []
    snames = np.array([i.rsplit(":", 2)[0].rsplit("_", 1)[0] for i in keys])
    newarr = np.zeros((len(set(snames)) + 1, arr.shape[1]), dtype="S1")
    
    # put reference into arr
    newarr[0] = arr[0]
    
    # fill rest while merging dups
    nidx = 1
    seen = set()
    for sidx, key in enumerate(keys):
        sname = snames[sidx]
        if sname not in seen:
            # add to list of seen names
            seen.add(sname)

            # get all rows of data for this sname (+1 b/c ref)
            didx = np.where(snames == sname)[0] + 1
            if didx.size > 1:
                iarr = arr[didx, :].view(np.uint8)
                iarr[iarr == 78] = 0
                iarr[iarr == 45] = 0
                if np.all(np.any(iarr == 0, axis=0)):
                    newarr[nidx] = iarr.max(axis=0).view("S1")
                else:
                    raise IPyradError("duplicate could not be resolved")
                # store key with reference to all dups
                ikeys = [keys[i - 1] for i in didx]
                fidxs = ";".join([i.rsplit("_", 1)[-1] for i in ikeys])
                newkeys.append("{}_{}".format(sname, fidxs))
                
            else:
                # store array data and orig key
                newarr[nidx] = arr[didx]
                newkeys.append(keys[sidx])
                
            nidx += 1
    
    # fill terminal edges with N again since array can increase
    newarr[newarr == b""] = b"N"
    return newkeys, newarr



def build_ref_regions(data):
    """
    Use bedtools to pull in consens reads overlapping some region of ref
    """
    cmd1 = [
        ipyrad.bins.bedtools,
        "bamtobed",
        "-i", 
        os.path.join(
            data.dirs.across,
            "{}.cat.sorted.bam".format(data.name)
        )
    ]

    cmd2 = [
        ipyrad.bins.bedtools, 
        "merge", 
        "-d", "0",
        "-i", "-",
    ]

    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc2 = sps.Popen(
        cmd2, 
        stdin=proc1.stdout,
        stderr=sps.STDOUT,
        stdout=sps.PIPE,
    )
    result = proc2.communicate()[0].decode()
    if proc2.returncode:
        raise IPyradError(
            "error in {}: {}".format(" ".join(cmd2), result))
    regs = [i.split("\t") for i in result.strip().split("\n")]
    return [(i, int(j), int(k)) for i, j, k in regs]



def build_ref_clusters(data, idx, iregion):
    """
    Given a chunk of regions this will pull in the reference for each region
    and then pull in all consens reads matching to that region. It uses cigar
    info to align the consens reads with the ref. This also merges consens
    from the same sample that were not merged earlier, which is why we expect
    no duplicate samples in the output of reference assemblies.
    """

    # prepare i/o for bamfile with mapped reads
    bamfile = AlignmentFile(
        os.path.join(
            data.dirs.across,
            "{}.cat.sorted.bam".format(data.name)),
        'rb')

    # dict to map chromosome names to integers
    faidict = chroms2ints(data, False)

    # prepare i/o for pysam reference indexed
    reffai = FastaFile(data.params.reference_sequence)

    # store path to cluster bit
    outbit = os.path.join(data.tmpdir, "aligned_{}.fa".format(idx))

    # get clusters
    iregions = iter(iregion)
    clusts = []

    while 1:
        # pull in all consens reads mapping to a bed region
        try:
            region = next(iregions)
            reads = bamfile.fetch(*region)
        except StopIteration:
            break

        # build a dict to reference seqs and cigars by name
        mstart = 9e12
        mend = 0
        rdict = {}
        for read in reads:
            rstart = read.reference_start
            rend = rstart + read.qlen
            mstart = min(mstart, rstart)
            mend = max(mend, rend)
            rdict[read.qname] = (read.seq, read.cigar, rstart, rend)
        keys = sorted(rdict.keys(), key=lambda x: x.rsplit(":", 2)[0])

        # pull in the reference for this region (1-indexed)
        refs = reffai.fetch(region[0], mstart, mend)

        # make empty array
        rlen = mend - mstart
        arr = np.zeros((len(keys) + 1, rlen), dtype=bytes)
        arr[0] = list(refs.upper())

        # fill arr with remaining samples
        for idx, key in enumerate(keys):
            seq, cigar, start, end = rdict[key]

            # how far ahead of ref start and short of ref end is this read
            fidx = start - mstart
            eidx = arr.shape[1] - (mend - end)

            # enter into the array, trim end if longer than pulled ref
            arr[idx + 1, fidx:eidx] = list(seq)[:eidx - fidx]

            # mod sequence according to cigar for indels and ambigs
            # csums is the location of impute on the seq, so it must be 
            # incremented by fidx and not extend past eidx
            for cidx, cig in enumerate(cigar):
                if cig[0] == 4:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[idx + 1, csums] = arr[idx + 1, csums].lower()
                if cig[0] == 1:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[idx + 1, csums] = b"-"

        # fill terminal edges with N
        arr[arr == b""] = b"N"

        # duplicates merge here (only perfect merge on all Ns) and reshape
        # the array to match. This will need to be resolved in catgs...
        # if it does not merge then
        try:
            keys, arr = resolve_duplicates(keys, arr)
        except IPyradError:
            pass

        # get consens seq and variant site index 
        clust = [">reference_{}:{}:{}-{}\n{}".format(
            0, 
            faidict[region[0]] + 1, mstart + 1, mend + 1,   # 1-indexed
            b"".join(arr[0]).decode()
        )]
        for idx, key in enumerate(keys):    
            clust.append(
                ">{}\n{}".format(key, b"".join(arr[idx + 1]).decode())
            )
        clusts.append("\n".join(clust))

    # dump to temp file until concat in next step.
    with open(outbit, 'w') as outfile:
        if clusts:
            outfile.write("\n//\n//\n".join(clusts) + "\n//\n//\n")



def build_concat_two(data, jobids, randomseed):
    """
    
    """
    seeds = [
        os.path.join(
            data.dirs.across, 
            "{}-{}.htemp".format(data.name, jobid)) for jobid in jobids
    ]
    allseeds = os.path.join(
        data.dirs.across, 
        "{}-x-catshuf.fa".format(data.name))
    cmd1 = ['cat'] + seeds
    cmd2 = [
        ipyrad.bins.vsearch, 
        '--sortbylength', '-', 
        '--fasta_width', '0',
        '--output', allseeds,
    ]
    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, close_fds=True)
    proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, close_fds=True)
    proc2.communicate()
    proc1.stdout.close()



def build_concat_files(data, jobid, samples, randomseed):
    """
    [This is returnn on an ipengine]
    Make a concatenated consens file with sampled alleles (no RSWYMK/rswymk).
    Orders reads by length and shuffles randomly within length classes
    """
    conshandles = [
        sample.files.consens for sample in samples if 
        sample.stats.reads_consens]
    conshandles.sort()
    assert conshandles, "no consensus files found"

    ## concatenate all of the gzipped consens files
    cmd = ['cat'] + conshandles
    groupcons = os.path.join(
        data.dirs.across, 
        "{}-{}-catcons.gz".format(data.name, jobid))
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



def cluster(data, jobid, nthreads, print_progress=False):
    """

    """
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
    cov = 0.5         # 0.90
    if data.params.datatype in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.60
    elif data.params.datatype == "pairgbs":
        strand = "both"
        cov = 0.75     # 0.90

    cmd = [ipyrad.bins.vsearch,
           "-cluster_smallmem", catshuf,
           "-strand", strand,
           "-query_cov", str(cov),
           "-minsl", str(0.5),
           "-id", str(data.params.clust_threshold),
           "-userout", uhaplos,
           "-notmatched", hhaplos,
           "-userfields", "query+target+qstrand",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-fasta_width", "0",
           "--minseqlength", str(data.params.filter_min_trim_len),
           "-threads", str(nthreads),  # "0",
           "-fulldp",
           "-usersort",
           ]

    if not print_progress:
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            raise IPyradError(out)

    else:
        (worker, boss) = pty.openpty()
        proc = sps.Popen(cmd, stdout=boss, stderr=boss, close_fds=True)
        prog = 0
        newprog = 0
        while 1:
            isdat = select.select([worker], [], [], 0)
            if isdat[0]:
                dat = os.read(worker, 80192).decode()
            else:
                dat = ""
            if "Clustering" in dat:
                try:
                    newprog = int(dat.split()[-1][:-1])
                # may raise value error when it gets to the end
                except ValueError:
                    pass
    
            # print progress (do not remove print statement, stdout is parsed)
            if newprog != prog:
                print(int(newprog))
                prog = newprog
                time.sleep(0.1)
        
            # break if done
            # catches end chunk of printing if clustering went really fast
            if "Clusters:" in dat:
                break

        # another catcher to let vsearch cleanup after clustering is done
        proc.wait()
        print(100)



def count_seeds(uhandle):
    """
    uses bash commands to quickly count N seeds from utemp file
    """
    with open(uhandle, 'r') as insort:
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



def sort_seeds(uhandle):
    "sort seeds from cluster results"
    cmd = ["sort", "-k", "2", uhandle, "-o", uhandle + ".sort"]
    proc = sps.Popen(cmd, close_fds=True)
    proc.communicate()



def build_single_denovo_clusters(data, usort, nseeds, *args):
    "use this function when not hierarchical clustering"
    # load all concat fasta files into a dictionary (memory concerns here...)
    conshandle = os.path.join(
        data.dirs.across, 
        "{}-0-catcons.gz".format(data.name),
    )
    allcons = {}
    with gzip.open(conshandle, 'rt') as iocons:
        cons = izip(*[iter(iocons)] * 2)
        for namestr, seq in cons:
            nnn, sss = [i.strip() for i in (namestr, seq)]
            allcons[nnn[1:]] = sss

    # load all utemp files into a dictionary
    usortfile = os.path.join(
        data.dirs.across,
        "{}-0.utemp.sort".format(data.name)
    )

    # set optim to approximately 4 chunks per core. Smaller allows for a bit
    # cleaner looking progress bar. 40 cores will make 160 files.
    # This often does not work as intended. iao 10/26/19
    # optim = ((nseeds // (data.ncpus * 4)) + (nseeds % (data.ncpus * 4)))
    optim = np.ceil(nseeds / (data.ncpus * 4))

    # iterate through usort grabbing seeds and matches
    with open(usortfile, 'rt') as insort:
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
        
            # store hit if still matching to same seed
            if seed == lastseed:
                if ori == "-":
                    seq = fullcomp(allcons[hit])[::-1]
                else:
                    seq = allcons[hit]
                fseqs.append(">{}\n{}".format(hit, seq))

            # store seed and hit (to a new cluster) if new seed.
            else:  
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

                # store the first hit to the seed
                seq = allcons[hit]
                if ori == "-":
                    seq = fullcomp(seq)[::-1]
                fseqs.append(">{}\n{}".format(hit, seq))

    # write whatever is left over to the clusts file
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


def build_hierarchical_denovo_clusters(data, usort, nseeds, jobids):
    "use this function when building clusters from hierarchical clusters"
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
    insort = open(usort, 'rt')
    isort = iter(insort)

    # seed null, and seqlist null
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

    # close handle
    insort.close()

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
        seqsize += 1
        loci += seqsize
    if seqlist:
        pathname = os.path.join(
            data.tmpdir, 
            data.name + ".chunk_{}".format(loci))
        with open(pathname, 'wt') as clustsout:
            clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")

    ## final progress and cleanup
    del allcons


def align_to_array(data, samples, chunk):
    """
    Opens a tmp clust chunk and iterates over align jobs.
    """
    # data are already chunked, read in the whole thing
    with open(chunk, 'rt') as infile:
        clusts = infile.read().split("//\n//\n")[:-1]    

    # snames to ensure sorted order
    samples.sort(key=lambda x: x.name)

    # create a persistent shell for running muscle in. 
    proc = sps.Popen(["bash"], stdin=sps.PIPE, stdout=sps.PIPE, bufsize=0)

    # iterate over clusters until finished
    allstack = []
    for ldx in range(len(clusts)):
        istack = []
        lines = clusts[ldx].strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]

        # skip aligning and continue if duplicates present (locus too big)
        # but reshape locs to be same lengths by adding --- to end, this 
        # simplifies handling them in step7 (they're still always filtered)
        unames = set([i.rsplit("_", 1)[0] for i in names])
        if len(unames) < len(names):
            longname = max([len(i) for i in seqs])
            seqs = [i.ljust(longname, "-") for i in seqs]
            istack = [">{}\n{}".format(i[1:], j) for i, j in zip(names, seqs)]
            allstack.append("\n".join(istack))
            continue

        # else locus looks good, align it.
        # is there a paired-insert in any samples in the locus?
        try:

            # try to split cluster list at nnnn separator for each read
            left = [i.split("nnnn")[0] for i in seqs]
            right = [i.split("nnnn")[1] for i in seqs]

            # align separately
            istack1 = muscle_it(proc, names, left)
            istack2 = muscle_it(proc, names, right)

            # combine in order
            for sdx in range(len(istack1)):
                n1, s1 = istack1[sdx].split("\n")
                s2 = istack2[sdx].split("\n")[-1]
                istack.append(n1 + "\n" + s1 + "nnnn" + s2)

        # no insert just align a single locus
        except IndexError:
            istack = muscle_it(proc, names, seqs)

        # store the locus
        if istack:
            allstack.append("\n".join(istack))

    # cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    # write to file when chunk is finished
    odx = chunk.rsplit("_")[-1]
    alignfile = os.path.join(data.tmpdir, "aligned_{}.fa".format(odx))
    with open(alignfile, 'wt') as outfile:
        outfile.write("\n//\n//\n".join(allstack) + "\n//\n//\n")



def muscle_it(proc, names, seqs):
    """
    Align with muscle, ensure name order, and return as string
    """  
    istack = []

    # append counter to names because muscle doesn't retain order
    nnames = [">{};*{}".format(j[1:], i) for i, j in enumerate(names)]
    
    # make back into strings
    cl1 = "\n".join(["\n".join(i) for i in zip(nnames, seqs)])

    # store allele (lowercase) info, returns mask with lowercases
    amask, abool = store_alleles(seqs)

    # send align1 to the bash shell (TODO: check for pipe-overflow)
    cmd1 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
            .format(cl1, ipyrad.bins.muscle, "//\n"))
    proc.stdin.write(cmd1.encode())

    # read the stdout by line until splitter is reached
    align1 = []
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
        (len(names), len(dalign1[keys[0]].replace("\n", ""))),
        dtype='S1',
    )
    for kidx, key in enumerate(keys):
        concatseq = dalign1[key].replace("\n", "")
        seqarr[kidx] = list(concatseq)

    # get alleles back using fast jit'd function.
    if np.sum(amask):
        intarr = seqarr.view(np.uint8)
        iamask = retrieve_alleles_after_aligning(intarr, amask)
        seqarr[iamask] = np.char.lower(seqarr[iamask])

    # sort in sname (alphanumeric) order. 
    istack = []    
    wkeys = np.argsort([i.rsplit("_", 1)[0] for i in keys])
    for widx in wkeys:
        wname = names[widx]
        istack.append(
            "{}\n{}".format(wname, b"".join(seqarr[widx]).decode()))
    return istack


def store_alleles(seqs):
    """
    Returns a mask selecting columns with lower case calls, and 
    a boolean of whether or not any exist. This is used to put them 
    back into alignments after muscle destroys all this info during
    alignment.
    """
    # get shape of the array and new empty array
    shape = (len(seqs), max([len(i) for i in seqs]))
    arrseqs = np.zeros(shape, dtype=np.bytes_)

    # iterate over rows 
    for row in range(arrseqs.shape[0]):
        seqsrow = seqs[row]
        arrseqs[row, :len(seqsrow)] = list(seqsrow)

    # mask the lo...
    amask = np.char.islower(arrseqs)
    if np.any(amask):
        return amask, True
    else:
        return amask, False


def retrieve_alleles_after_aligning(intarr, amask):
    """
    Imputes lower case allele calls back into alignments 
    while taking account for spacing caused by insertions.
    """
    newmask = np.zeros(intarr.shape, dtype=np.bool_)
    
    for ridx in range(intarr.shape[0]):
        iarr = intarr[ridx]
        indidx = np.where(iarr == 45)[0]
        
        # if no indels then simply use the existing mask
        if not indidx.size:
            newmask[ridx] = amask[ridx]
            
        # if indels that impute 
        else:
            allrows = np.arange(amask.shape[1])
            mask = np.ones(allrows.shape[0], dtype=np.bool_)
            for idx in indidx:
                if idx < mask.shape[0]:
                    mask[idx] = False
            not_idx = allrows[mask == 1]
            
            # fill in new data into all other spots
            newmask[ridx, not_idx] = amask[ridx, :not_idx.shape[0]]
    return newmask



if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    # tdata = ip.load_json("/tmp/test-simpairddrad.json")
    # tdata.run("6", auto=True, force=True)
    # logger.info(tdata.stats.T)


    tdata = ip.load_json("/tmp/test-amaranth.json")
    tdata.run("6", auto=True, force=True)
    print(tdata.stats)
    # print(tdata.stats_dfs.s5)

    # self.data.hackersonly.declone_PCR_duplicates:
    # tdata = ip.load_json("/tmp/test-amaranth-denovo.json")
    # tdata.ipcluster['cores'] = 4
    # tdata.run("6", auto=True, force=True)
    # logger.info(tdata.stats.T)
