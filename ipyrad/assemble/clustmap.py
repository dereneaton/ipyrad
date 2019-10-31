#!/usr/bin/env python

"""
de-replicates edit files and clusters de-replciated reads
by sequence similarity using vsearch
"""

from __future__ import print_function
try:
    from itertools import izip, islice, chain
except ImportError:
    from itertools import islice, chain
    izip = zip

import os
import gzip
import glob
import time
import shutil
import warnings
import subprocess as sps

import numpy as np
import pysam
import ipyrad as ip
from .utils import IPyradError, bcomp, comp


class Step3:
    "Class for organizing step functions across datatypes and read formats"

    def __init__(self, data, force, ipyclient):

        # store attributes
        self.data = data
        # self.noreverse = noreverse
        self.maxindels = 8
        self.force = force
        self.ipyclient = ipyclient
        self.gbs = bool("gbs" in self.data.params.datatype)
        self.print_headers()
        self.samples = self.get_subsamples()

        # init funcs
        self.setup_dirs()
        self.tune_threads()


    def run(self):
        "Run the assembly functions for this step"

        # paired-end data methods ------------------------------
        if "pair" in self.data.params.datatype:

            # Filter out by mapping to alt reference
            # if self.data.params.reference_as_filter:
            #     self.remote_index_refs("alt")
            #     self.remote_run(
            #         function=map_alt_reference,
            #         printstr=("filter map to alt ref", "s3"),
            #         args=(),
            #         threaded=True,
            #     )

            # DENOVO ----
            if self.data.params.assembly_method == "denovo":
                self.remote_run(
                    function=concat_multiple_edits,
                    printstr=("concatenating       ", "s3"),
                    args=(),
                )
                self.remote_run(
                    function=merge_pairs_with_vsearch,
                    printstr=("join merged pairs   ", "s3"),
                    args=(True,),
                )
                self.remote_run(
                    function=merge_end_to_end,
                    printstr=("join unmerged pairs ", "s3"),
                    args=(True, True,),
                )
                self.remote_run(
                    function=dereplicate,
                    printstr=("dereplicating       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run_cluster_build()
                self.remote_run_align_cleanup()


            # REFERENCE ----
            elif self.data.params.assembly_method == "reference":
                self.remote_index_refs()
                self.remote_run(
                    function=concat_multiple_edits,
                    printstr=("concatenating       ", "s3"),
                    args=(),
                    threaded=True,
                )
                self.remote_run(
                    function=merge_end_to_end,
                    printstr=("join unmerged pairs ", "s3"),
                    args=(False, False,),
                    threaded=True,
                )
                if self.data.hackersonly.declone_PCR_duplicates:
                    self.remote_run(
                        function=declone_3rad,
                        printstr=("declone 3RAD        ", "s3"),
                        args=(self.nthreads, self.force),
                    )
                self.remote_run(
                    function=dereplicate,
                    printstr=("dereplicating       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run(
                    function=split_endtoend_reads,
                    printstr=("splitting dereps    ", "s3"),
                    args=(),
                )
                self.remote_run(
                    function=mapping_reads,
                    printstr=("mapping reads       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run(
                    function=build_clusters_from_cigars,
                    printstr=("building clusters   ", "s3"),
                    args=(),
                )

            # DENOVO MINUS
            elif self.data.params.assembly_method == "denovo-reference":
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

            elif self.data.params.assembly_method == "denovo+reference":
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

            else:
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

        ## single-end methods ------------------------------------
        else:
            # DENOVO
            if self.data.params.assembly_method == "denovo":
                self.remote_run(
                    function=concat_multiple_edits,
                    printstr=("concatenating       ", "s3"),
                    args=(),
                )
                self.remote_run(
                    function=dereplicate,
                    printstr=("dereplicating       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run_cluster_build()
                self.remote_run_align_cleanup()

            # REFERENCE
            elif self.data.params.assembly_method == "reference":
                self.remote_index_refs()
                self.remote_run(
                    function=concat_multiple_edits,
                    printstr=("concatenating       ", "s3"),
                    args=(),
                )
                self.remote_run(
                    function=dereplicate,
                    printstr=("dereplicating       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run(
                    function=mapping_reads,
                    printstr=("mapping reads       ", "s3"),
                    args=(self.nthreads,),
                    threaded=True,
                )
                self.remote_run(
                    function=build_clusters_from_cigars,
                    printstr=("building clusters   ", "s3"),
                    args=(),
                )

            # DENOVO MINUS
            elif self.data.params.assembly_method == "denovo-reference":
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

            elif self.data.params.assembly_method == "denovo+reference":
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

            else:
                raise NotImplementedError(
                    "datatype + assembly_method combo not yet supported")

        self.remote_run_sample_cleanup()
        self.cleanup()


    def print_headers(self):
        # print headers
        if self.data._cli:
            self.data._print(
                "\n{}Step 3: Clustering/Mapping reads within samples"
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        "Apply state, ncluster, and force filters to select samples"

        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 3")

        # filter samples by state
        state1 = self.data.stats.index[self.data.stats.state < 2]
        state2 = self.data.stats.index[self.data.stats.state == 2]
        state3 = self.data.stats.index[self.data.stats.state > 2]

        # tell user which samples are not ready for step 3
        if state1.any():
            self.data._print("skipping samples not yet in state==2:\n{}"
                  .format(state1.tolist()))

        if self.force:
            # run all samples above state 1
            subs = self.data.stats.index[self.data.stats.state > 1]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already cmopleted step 3
            if state3.any():
                self.data._print("skipping samples already finished step 3:\n{}"
                      .format(state3.tolist()))

            # run all samples in state 2
            subsamples = [self.data.samples[i] for i in state2]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_passed_filter:
                checked_samples.append(sample)
            else:
                self.data._print("skipping {}; no reads found.")
        if not any(checked_samples):
            raise IPyradError("No samples ready for step 3.")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_passed_filter,
            reverse=True,
        )
        return checked_samples


    def setup_dirs(self):
        "setup directories for the tmp files and cluster/ref outputs"
        # make output folder for clusters
        pdir = os.path.realpath(self.data.params.project_dir)
        self.data.dirs.clusts = os.path.join(
            pdir, "{}_clust_{}"
            .format(self.data.name, self.data.params.clust_threshold))

        if not os.path.exists(self.data.dirs.clusts):
            os.mkdir(self.data.dirs.clusts)

        # make a tmpdir for align files
        self.data.tmpdir = os.path.abspath(os.path.expanduser(
            os.path.join(pdir, self.data.name + '-tmpalign')))
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)

        # If ref mapping, init samples and make refmapping output directory
        if self.data.params.assembly_method != "denovo":
            self.data.dirs.refmapping = os.path.join(
                pdir, "{}_refmapping".format(self.data.name))
            if not os.path.exists(self.data.dirs.refmapping):
                os.mkdir(self.data.dirs.refmapping)

        # set a filepath for stored cluster results
        for sname in self.data.samples:
            self.data.samples[sname].files.clusters = os.path.join(
                self.data.dirs.clusts,
                "{}.clustS.gz".format(sname))


    def tune_threads(self):
        "setup threading to efficiently run clust/ref across HPC"
        # set nthreads based on _ipcluster dict (default is 2)
        if "threads" in self.data.ipcluster.keys():
            self.nthreads = int(self.data.ipcluster["threads"])

        # create standard load-balancers
        self.lbview = self.ipyclient.load_balanced_view()
        self.thview = self.ipyclient.load_balanced_view()  # to be threaded

        # if nthreads then scale thview to use threads
        eids = self.ipyclient.ids
        if self.nthreads:
            if self.nthreads <= len(self.ipyclient.ids):
                self.thview = self.ipyclient.load_balanced_view(
                    targets=eids[::self.nthreads])

        # else try auto-tuning to 2 or 4 threaded
        else:
            if len(self.ipyclient) >= 40:
                self.thview = self.ipyclient.load_balanced_view(
                    targets=eids[::4])
            else:
                self.thview = self.ipyclient.load_balanced_view(
                    targets=eids[::2])


    def cleanup(self):
        "cleanup / statswriting function for Assembly obj"
        self.data.stats_dfs.s3 = self.data._build_stat("s3")
        self.data.stats_files.s3 = os.path.join(
            self.data.dirs.clusts, "s3_cluster_stats.txt")
        with open(self.data.stats_files.s3, 'w') as outfile:
            self.data.stats_dfs.s3.to_string(
                buf=outfile,
                formatters={
                    'merged_pairs': '{:.0f}'.format,
                    'clusters_total': '{:.0f}'.format,
                    'clusters_hidepth': '{:.0f}'.format,
                    'filtered_bad_align': '{:.0f}'.format,
                    'avg_depth_stat': '{:.2f}'.format,
                    'avg_depth_mj': '{:.2f}'.format,
                    'avg_depth_total': '{:.2f}'.format,
                    'sd_depth_stat': '{:.2f}'.format,
                    'sd_depth_mj': '{:.2f}'.format,
                    'sd_depth_total': '{:.2f}'.format
                })

        # remove temporary alignment chunk and derep files
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)


    def remote_index_refs(self):
        "index the reference seq for bwa and samtools (pysam)"
        start = time.time()
        printstr = ("indexing reference  ", "s3")
        rasync1 = self.lbview.apply(
            index_ref_with_bwa, self.data)
        rasync2 = self.lbview.apply(
            index_ref_with_sam, self.data)

        # track job
        while 1:
            ready = [rasync1.ready(), rasync2.ready()]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        # check for errors
        self.data._print("")
        for job in [rasync1, rasync2]:
            if not job.successful():
                job.get()


    def remote_run_cluster_build(self):
        # submit clustering/mapping job
        start = time.time()
        casyncs = {}
        for sample in self.samples:
            casyncs[sample.name] = self.thview.apply(
                cluster,
                *(self.data, sample, self.nthreads, self.force)
            )

        # submit cluster building job
        basyncs = {}
        for sample in self.samples:
            with self.lbview.temp_flags(after=casyncs[sample.name]):
                basyncs[sample.name] = self.lbview.apply(
                    build_clusters,
                    *(self.data, sample, self.maxindels)
                )

        # submit cluster chunking job
        hasyncs = {}
        for sample in self.samples:
            with self.lbview.temp_flags(after=basyncs[sample.name]):
                hasyncs[sample.name] = self.lbview.apply(
                    muscle_chunker,
                    *(self.data, sample)
                )

        # track job progress
        printstr = ("clustering/mapping  ", "s3")
        while 1:
            ready = [casyncs[i].ready() for i in casyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        self.data._print("")
        for job in casyncs:            
            if not casyncs[job].successful():
                casyncs[job].get()

        # track job progress
        start = time.time()
        printstr = ("building clusters   ", "s3")
        while 1:
            ready = [basyncs[i].ready() for i in basyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        self.data._print("")
        for job in basyncs:
            if not basyncs[job].successful():
                basyncs[job].get()

        # track job progress
        start = time.time()
        printstr = ("chunking clusters   ", "s3")
        while 1:
            ready = [hasyncs[i].ready() for i in hasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        self.data._print("")
        for job in hasyncs:
            if not hasyncs[job].successful():
                hasyncs[job].get()


    def remote_run_align_cleanup(self):

        # submit ten aligning jobs for each sample
        start = time.time()
        aasyncs = {}
        for sample in self.samples:
            aasyncs[sample.name] = []
            for idx in range(10):
                handle = os.path.join(
                    self.data.tmpdir,
                    "{}_chunk_{}.ali".format(sample.name, idx))

                rasync = self.lbview.apply(
                    align_and_parse,
                    *(handle, self.maxindels, self.gbs)
                )
                aasyncs[sample.name].append(rasync)

        # a list with all aasyncs concatenated
        allasyncs = list(chain(*[aasyncs[i] for i in aasyncs]))

        # submit cluster building job for each sample *after* all align jobs
        basyncs = {}
        for sample in self.samples:
            with self.lbview.temp_flags(after=aasyncs[sample.name]):
                basyncs[sample.name] = self.lbview.apply(
                    reconcat,
                    *(self.data, sample)
                )

        # track job 1 progress
        printstr = ("aligning clusters   ", "s3")
        while 1:
            ready = [i.ready() for i in allasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        self.data._print("")
        for job in allasyncs:
            if not job.successful():
                job.get()

        # track job 2 progress
        start = time.time()
        printstr = ("concat clusters     ", "s3")
        while 1:
            ready = [basyncs[i].ready() for i in basyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        self.data._print("")
        for job in basyncs:
            if not basyncs[job].successful():
                basyncs[job].get()


    def remote_run_sample_cleanup(self):
        # submit job
        printstr = ("calc cluster stats  ", "s3")
        start = time.time()
        rasyncs = {}
        njobs = len(self.samples)
        for sample in self.samples:
            args = [self.data, sample]
            rasyncs[sample.name] = self.lbview.apply(get_quick_depths, *args)

        # enter result stats as the jobs finish
        finished = 0       
        while 1:
            samplelist = list(rasyncs.keys())
            for sname in samplelist:
                if rasyncs[sname].ready():

                    # enter results to sample object and checks for errors
                    maxlens, depths = rasyncs[sname].get()
                    store_sample_stats(
                        self.data, self.data.samples[sname], maxlens, depths)
                    finished += 1

                    # remove sample from todo list, and del from rasyncs mem
                    rasyncs.pop(sname)

            # progress bar
            self.data._progressbar(njobs, finished, start, printstr)
            time.sleep(0.1)
            if finished == njobs:
                break
        self.data._print("")


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

        # check for errors, will raise ipp.RemoteError
        self.data._print("")
        for job in rasyncs:
            rasyncs[job].get()

        # clean up to free any RAM
        self.ipyclient.purge_everything()


def dereplicate(data, sample, nthreads):
    """
    Dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file. Paired
    reads are dereplicated as one concatenated read and later split again.
    Updated this function to take infile and outfile to support the double
    dereplication that we need for 3rad (5/29/15 iao).
    """
    # find input file with following precedence:
    # .trimmed.fastq.gz, .concatedit.fq.gz, ._merged.fastq, ._declone.fastq
    infiles = [
        os.path.join(data.dirs.edits,
            "{}.trimmed_R1_.fastq.gz".format(sample.name)),
        os.path.join(data.dirs.edits, 
            "{}_R1_concatedit.fq.gz".format(sample.name)),
        os.path.join(data.tmpdir, 
            "{}_merged.fastq".format(sample.name)),
        os.path.join(data.tmpdir, 
            "{}_declone.fastq".format(sample.name)),
    ]
    infiles = [i for i in infiles if os.path.exists(i)]
    infile = infiles[-1]

    # datatypes options
    strand = "plus"
    if data.params.datatype is ('gbs' or '2brad'):
        strand = "both"

    # do dereplication with vsearch
    cmd = [
        ip.bins.vsearch,
        "--derep_fulllength", infile,
        "--strand", strand,
        "--output", os.path.join(data.tmpdir, sample.name + "_derep.fastq"),
        "--threads", str(nthreads),
        "--fasta_width", str(0),
        "--sizeout", 
        "--relabel_md5",
        "--quiet",
        #"--fastq_qmax", "1000",        
    ]

    # decompress argument
    if infile.endswith(".gz"):
        cmd.append("--gzip_decompress")

    # build PIPEd job   
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    errmsg = proc.communicate()[0]
    if proc.returncode:
        raise IPyradError(errmsg.decode())


def concat_multiple_edits(data, sample):

    # define output files
    concat1 = os.path.join(
        data.dirs.edits,
        "{}_R1_concatedit.fq.gz".format(sample.name))
    concat2 = os.path.join(
        data.dirs.edits,
        "{}_R2_concatedit.fq.gz".format(sample.name))

    if len(sample.files.edits) > 1:
        # cat all inputs; index 0 b/c files are in tuples for r1, r2
        cmd1 = ["cat"] + [i[0] for i in sample.files.edits]

        # write to new concat handle
        with open(concat1, 'w') as cout1:
            proc1 = sps.Popen(
                cmd1, stderr=sps.STDOUT, stdout=cout1, close_fds=True)
            res1 = proc1.communicate()[0]
            if proc1.returncode:
                raise IPyradError("error in: %s, %s", cmd1, res1)

        ## Only set conc2 if R2 actually exists
        if os.path.exists(str(sample.files.edits[0][1])):
            cmd2 = ["cat"] + [i[1] for i in sample.files.edits]
            with gzip.open(concat2, 'w') as cout2:
                proc2 = sps.Popen(
                    cmd2, stderr=sps.STDOUT, stdout=cout2, close_fds=True)
                res2 = proc2.communicate()[0]
                if proc2.returncode:
                    raise IPyradError("error in: %s, %s", cmd2, res2)


def merge_pairs_with_vsearch(data, sample, revcomp):
    "Merge PE reads using vsearch to find overlap."

    # define input files
    concat1 = os.path.join(
        data.dirs.edits, 
        "{}_R1_concatedit.fq.gz".format(sample.name))
    concat2 = os.path.join(
        data.dirs.edits, 
        "{}_R2_concatedit.fq.gz".format(sample.name))
    if os.path.exists(concat1):
        infile1 = concat1
    else:
        infile1 = sample.files.edits[0][0]      
    if os.path.exists(concat2):
        infile2 = concat2
    else:
        infile2 = sample.files.edits[0][1]

    # define output files
    mergedfile = os.path.join(
        data.tmpdir,
        "{}_merged.fastq".format(sample.name))
    nonmerged1 = os.path.join(
        data.tmpdir,
        "{}_nonmerged_R1_.fastq".format(sample.name))
    nonmerged2 = os.path.join(
        data.tmpdir,
        "{}_nonmerged_R2_.fastq".format(sample.name))

    ## get the maxn and minlen values
    try:
        maxn = sum(data.params.max_low_qual_bases)
    except TypeError:
        maxn = data.params.max_low_qual_bases
    minlen = str(max(32, data.params.filter_min_trim_len))

    # vsearch merge can now take gzipped files (v.2.8)
    cmd = [
        ip.bins.vsearch,
        "--fastq_mergepairs", infile1,
        "--reverse", infile2,
        "--fastqout", mergedfile,
        "--fastqout_notmerged_fwd", nonmerged1,
        "--fastqout_notmerged_rev", nonmerged2,
        "--fasta_width", "0",
        "--fastq_minmergelen", minlen,
        "--fastq_maxns", str(maxn),
        "--fastq_minovlen", "20",
        "--fastq_maxdiffs", "4",
        "--label_suffix", "_m1",
        #"--fastq_qmax", "100",
        "--threads", "2",
        "--fastq_allowmergestagger",
    ]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    res = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError("Error merge pairs:\n {}\n{}".format(cmd, res))


def merge_end_to_end(data, sample, revcomp, append):
    """ 
    Combines read1 and read2 with a 'nnnn' separator. If the data are going 
    to be refmapped then do not revcomp the read2. 
    """
    # input files; 
    mergedfile = os.path.join(
        data.tmpdir, 
        "{}_merged.fastq".format(sample.name))
    nonmerged1 = os.path.join(
        data.tmpdir, 
        "{}_nonmerged_R1_.fastq".format(sample.name))
    nonmerged2 = os.path.join(
        data.tmpdir, 
        "{}_nonmerged_R2_.fastq".format(sample.name))
    concat1 = os.path.join(
        data.tmpdir,
        "{}_R1_concatedit.fq.gz".format(sample.name))
    concat2 = os.path.join(
        data.tmpdir, 
        "{}_R2_concatedit.fq.gz".format(sample.name))
    edits1 = os.path.join(
        data.dirs.edits,
        "{}.trimmed_R1_.fastq.gz".format(sample.name))
    edits2 = os.path.join(
        data.dirs.edits, 
        "{}.trimmed_R2_.fastq.gz".format(sample.name))

    # file precedence
    nonm1 = [i for i in (edits1, concat1, nonmerged1) if os.path.exists(i)][-1]
    nonm2 = [i for i in (edits2, concat2, nonmerged2) if os.path.exists(i)][-1]

    # Combine the unmerged pairs and append to the merge file
    if append:
        combout = open(mergedfile, 'a')
    else:
        combout = open(mergedfile, 'w')

    # read in paired end read files 4 lines at a time
    if nonm1.endswith(".gz"):
        fr1 = gzip.open(nonm1, 'rb')
    else:
        fr1 = open(nonm1, 'rb')
    quart1 = izip(*[iter(fr1)] * 4)

    if nonm2.endswith(".gz"):
        fr2 = gzip.open(nonm2, 'rb')
    else:
        fr2 = open(nonm2, 'rb')
    quart2 = izip(*[iter(fr2)] * 4)
    quarts = izip(quart1, quart2)

    ## a list to store until writing
    writing = []
    counts = 0

    ## iterate until done
    while 1:
        try:
            read1s, read2s = next(quarts)
        except StopIteration:
            break
        if revcomp:
            writing.append(b"".join([
                read1s[0],
                read1s[1].strip() + b"nnnn" + (
                    bcomp(read2s[1].strip()[::-1]) + b"\n"),
                read1s[2],
                read1s[3].strip() + b"nnnn" + (
                    read2s[3].strip()[::-1] + b"\n"),
                ]))
        else:
            writing.append(b"".join([
                read1s[0],
                read1s[1].strip() + b"nnnn" + (
                    read2s[1]),
                read1s[2],
                read1s[3].strip() + b"nnnn" + (
                    read2s[3]),
                ]))

        counts += 1
        if not counts % 5000:
            combout.write(b"".join(writing).decode())
            writing = []

    if writing:
        combout.write(b"".join(writing).decode())

    ## close handles
    fr1.close()
    fr2.close()
    combout.close()


def count_merged_reads(data, sample):
    # record how many read pairs were merged
    mergedfile = os.path.join(
        data.tmpdir, 
        "{}_merged.fastq".format(sample.name))  
    with open(mergedfile, 'r') as tmpf:
        nmerged = sum(1 for i in tmpf.readlines()) // 4
    return nmerged


def cluster(data, sample, nthreads, force):
    """
    Calls vsearch for clustering. cov varies by data type, values were chosen
    based on experience, but could be edited by users
    """
    # get dereplicated reads for denovo+reference or denovo-reference
    derephandle = os.path.join(
        data.tmpdir, "{}_derep.fastq".format(sample.name))
    assert os.path.exists(derephandle), "bad derep handle"

    # create handles for the outfiles
    uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")
                       
    ## datatype specific optimization
    ## minsl: the percentage of the seed that must be matched
    ##    smaller values for RAD/ddRAD where we might want to combine, say 50bp
    ##    reads and 100bp reads in the same analysis.
    ## query_cov: the percentage of the query sequence that must match seed
    ##    smaller values are needed for gbs where only the tips might overlap
    ##    larger values for pairgbs where they should overlap near completely
    ##    small minsl and high query cov allows trimmed reads to match to untrim
    ##    seed for rad/ddrad/pairddrad.
    strand = "plus"
    cov = 0.5
    minsl = 0.5
    if data.params.datatype in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.5
        minsl = 0.5
    elif data.params.datatype == 'pairgbs':
        strand = "both"
        cov = 0.75
        minsl = 0.75

    ## If this value is not null (which is the default) then override query cov
    if data.hackersonly.query_cov:
        cov = str(data.hackersonly.query_cov)
        assert float(cov) <= 1, "query_cov must be <= 1.0"

    ## get call string
    cmd = [ip.bins.vsearch,
           "-cluster_smallmem", derephandle,
           "-strand", strand,
           "-query_cov", str(cov),
           "-id", str(data.params.clust_threshold),
           "-minsl", str(minsl),
           "-userout", uhandle,
           "-userfields", "query+target+id+gaps+qstrand+qcov",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-threads", str(nthreads),
           "-notmatched", temphandle,
           "-fasta_width", "0",
           # "-fastq_qmax", "100",
           "-fulldp",
           "-usersort"]

    ## run vsearch
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    res = proc.communicate()[0]

    # check for errors
    if proc.returncode:
        raise IPyradError("cmd {}: {}".format(cmd, res))


def build_clusters(data, sample, maxindels):
    """
    Combines information from .utemp and .htemp files to create .clust files,
    which contain un-aligned clusters. Hits to seeds are only kept in the
    cluster if the number of internal indels is less than 'maxindels'.
    By default, we set maxindels=6 for this step (within-sample clustering).
    """

    ## If reference assembly then here we're clustering the unmapped reads
    if "reference" in data.params.assembly_method:
        derepfile = os.path.join(
            data.tmpdir, sample.name + "-refmap_derep.fastq")
    else:
        derepfile = os.path.join(data.tmpdir, sample.name + "_derep.fastq")

    ## i/o vsearch files
    uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
    usort = os.path.join(data.dirs.clusts, sample.name + ".utemp.sort")
    hhandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")
    clustsout = open(
        os.path.join(
            data.dirs.clusts, 
            "{}.clust.txt".format(sample.name)), 
        'w')

    ## Sort the uhandle file so we can read through matches efficiently
    cmd = ["sort", "-k", "2", uhandle, "-o", usort]
    proc = sps.Popen(cmd, close_fds=True)
    proc.communicate()[0]

    ## load ALL derep reads into a dictionary (this can be a few GB of RAM)
    ## and is larger if names are larger. We are grabbing two lines at a time.
    alldereps = {}
    with open(derepfile, 'rt') as ioderep:
        dereps = izip(*[iter(ioderep)] * 2)
        for namestr, seq in dereps:
            nnn, sss = [i.strip() for i in (namestr, seq)]  
            alldereps[nnn[1:]] = sss

    ## store observed seeds (this could count up to >million in bad data sets)
    seedsseen = set()

    ## Iterate through the usort file grabbing matches to build clusters
    with open(usort, 'rt') as insort:
        ## iterator, seed null, seqlist null
        isort = iter(insort)
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0
        while 1:
            ## grab the next line
            try:
                hit, seed, _, ind, ori, _ = next(isort).strip().split()
            except StopIteration:
                break

            ## same seed, append match
            if seed != lastseed:
                seedsseen.add(seed)
                ## store the last cluster (fseq), count it, and clear fseq
                if fseqs:
                    ## sort fseqs by derep after pulling out the seed
                    fseqs = [fseqs[0]] + sorted(fseqs[1:], 
                        key=lambda x: 
                            int(x.split(";size=")[-1].split("\n")[0][:-2]), 
                        reverse=True)                    
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                # occasionally write/dump stored clusters to file and clear mem
                if not seqsize % 10000:
                    if seqlist:
                        clustsout.write(
                            "\n//\n//\n".join(seqlist) + "\n//\n//\n")
                        ## reset list and counter
                        seqlist = []

                ## store the new seed on top of fseq list
                fseqs.append(">{};*\n{}".format(seed, alldereps[seed]))
                lastseed = seed

            ## add match to the seed
            ## revcomp if orientation is reversed (comp preserves nnnn)
            if ori == "-":
                seq = comp(alldereps[hit])[::-1]
            else:
                seq = alldereps[hit]
            # only save if not too many indels
            if int(ind) <= maxindels:
                fseqs.append(">{};{}\n{}".format(hit, ori, seq))

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
    if seqlist:
        clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")

    ## now write the seeds that had no hits. Make dict from htemp
    with open(hhandle, 'rt') as iotemp:
        nohits = izip(*[iter(iotemp)] * 2)
        seqlist = []
        seqsize = 0
        while 1:
            try:
                nnn, _ = [i.strip() for i in next(nohits)]
            except StopIteration:
                break

            ## occasionally write to file
            if not seqsize % 10000:
                if seqlist:
                    clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")
                    ## reset list and counter
                    seqlist = []

            ## append to list if new seed
            if nnn[1:] not in seedsseen:
                seqlist.append("{};*\n{}".format(nnn, alldereps[nnn[1:]]))
                seqsize += 1

    ## write whatever is left over to the clusts file
    if seqlist:
        clustsout.write("\n//\n//\n".join(seqlist))

    ## close the file handle
    clustsout.close()
    del alldereps


def muscle_chunker(data, sample):
    """
    Splits the muscle alignment into chunks. Each chunk is run on a separate
    computing core. Because the largest clusters are at the beginning of the 
    clusters file, assigning equal clusters to each file would put all of the 
    large cluster, that take longer to align, near the top. So instead we 
    randomly distribute the clusters among the files. If assembly method is
    reference then this step is just a placeholder and nothing happens. 
    """

    ## only chunk up denovo data, refdata has its own chunking method which 
    ## makes equal size chunks, instead of uneven chunks like in denovo
    if data.params.assembly_method != "reference":
        ## get the number of clusters
        clustfile = os.path.join(data.dirs.clusts, sample.name + ".clust.txt")
        with iter(open(clustfile, 'rt')) as clustio:
            nloci = sum(1 for i in clustio if "//" in i) // 2
            optim = (nloci // 20) + (nloci % 20)
            #io
            #optim = int(nloci/10)

        ## write optim clusters to each tmp file
        clustio = open(clustfile, 'rt')
        inclusts = iter(clustio.read().strip().split("//\n//\n"))
        
        ## splitting loci so first file is smaller and last file is bigger
        inc = optim // 10
        for idx in range(10):
            ## how big is this chunk?
            this = optim + (idx * inc)
            left = nloci - this
            ## Distribute all loci equally among chunks
            ## io
            ##left = nloci - optim
            ##this = optim
            if idx == 9:
                ## grab everything left
                grabchunk = list(islice(inclusts, int(1e9)))
            else:
                ## grab next chunks-worth of data
                grabchunk = list(islice(inclusts, this))
                nloci = left

            ## write the chunk to file
            tmpfile = os.path.join(
                data.tmpdir, sample.name + "_chunk_{}.ali".format(idx))
            with open(tmpfile, 'wb') as out:
                out.write(str.encode("//\n//\n".join(grabchunk)))
        clustio.close()


def align_and_parse(handle, max_internal_indels=5, is_gbs=False):
    """ much faster implementation for aligning chunks """

    # CHECK: data are already chunked, read in the whole thing. bail if no data.
    clusts = []
    try:
        with open(handle, 'rb') as infile:
            clusts = infile.read().decode().split("//\n//\n")
            # remove any empty spots
            clusts = [i for i in clusts if i]
            # Skip entirely empty chunks; return 0 if no clusters in file
            # Allows some chunks to be empty without raising an error.
            if not clusts:
                return 0

    # return 0 if file not read for some reason...
    except IOError:
        return 0

    ## count discarded clusters for printing to stats later
    highindels = 0

    ## iterate over clusters sending each to muscle, splits and aligns pairs
    aligned = persistent_popen_align3(clusts, 200, is_gbs)

    ## store good alignments to be written to file
    refined = []

    ## filter and trim alignments
    for clust in aligned:
        # check for too many internal indels
        if not aligned_indel_filter(clust, max_internal_indels):
            refined.append(clust)
        else:
            highindels += 1

    ## write to file after
    if refined:
        outhandle = handle.rsplit(".", 1)[0] + ".aligned"
        with open(outhandle, 'wb') as outfile:
            try:
                outfile.write("\n//\n//\n".join(refined) + "\n")
            except TypeError:
                outfile.write(("\n//\n//\n".join(refined) + "\n").encode())
    return highindels


def reconcat(data, sample):
    """ takes aligned chunks (usually 10) and concatenates them """

    ## get chunks
    chunks = glob.glob(os.path.join(data.tmpdir,
             sample.name + "_chunk_[0-9].aligned"))

    ## sort by chunk number, cuts off last 8 =(aligned)
    chunks.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))

    ## concatenate finished reads
    sample.files.clusters = os.path.join(
        data.dirs.clusts, sample.name + ".clustS.gz")

    ## reconcats aligned clusters
    with gzip.open(sample.files.clusters, 'wb') as out:
        for fname in chunks:
            with open(fname) as infile:
                dat = infile.read().strip()
                dat += "\n//\n//\n"
                try:
                    out.write(dat)
                except TypeError:
                    out.write(dat.encode())
            os.remove(fname)


def persistent_popen_align3(clusts, maxseqs=200, is_gbs=False):
    "keeps a persistent bash shell open and feeds it muscle alignments"

    # create a separate shell for running muscle in, this is much faster
    # than spawning a separate subprocess for each muscle call
    proc = sps.Popen(
        ["bash"],
        stdin=sps.PIPE,
        stdout=sps.PIPE,
        bufsize=0,
    )

    ## iterate over clusters in this file until finished
    aligned = []
    for clust in clusts:

        ## new alignment string for read1s and read2s
        align1 = []
        align2 = []

        ## don't bother aligning if only one seq
        if clust.count(">") == 1:
            aligned.append(clust.replace(">", "").strip())
        else:

            # do we need to split the alignment? (is there a PE insert?)
            try:
                # make into list (only read maxseqs lines, 2X cuz names)
                lclust = clust.split()[:maxseqs * 2]

                # try to split cluster list at nnnn separator for each read
                lclust1 = list(chain(*zip(
                    lclust[::2], [i.split("nnnn")[0] for i in lclust[1::2]])))
                lclust2 = list(chain(*zip(
                    lclust[::2], [i.split("nnnn")[1] for i in lclust[1::2]])))

                # put back into strings
                clust1 = "\n".join(lclust1)
                clust2 = "\n".join(lclust2)

                # Align the first reads.
                # The muscle command with alignment as stdin and // as split
                cmd1 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
                        .format(clust1, ip.bins.muscle, "//\n"))

                # send cmd1 to the bash shell
                proc.stdin.write(cmd1.encode())

                # read the stdout by line until splitter is reached
                # meaning that the alignment is finished.
                for line in iter(proc.stdout.readline, b'//\n'):
                    align1.append(line.decode())

                # Align the second reads.
                # The muscle command with alignment as stdin and // as split
                cmd2 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
                        .format(clust2, ip.bins.muscle, "//\n"))

                # send cmd2 to the bash shell
                proc.stdin.write(cmd2.encode())

                # read the stdout by line until splitter is reached
                # meaning that the alignment is finished.
                for line in iter(proc.stdout.readline, b'//\n'):
                    align2.append(line.decode())

                # join up aligned read1 and read2 and ensure names order match
                lines1 = "".join(align1)[1:].split("\n>")
                lines2 = "".join(align2)[1:].split("\n>")
                dalign1 = dict([i.split("\n", 1) for i in lines1])
                dalign2 = dict([i.split("\n", 1) for i in lines2])

                # sort the first reads
                keys = list(dalign1.keys())
                seed = [i for i in keys if i[-1] == "*"][0]
                keys.pop(keys.index(seed))
                order = [seed] + sorted(
                    keys, key=get_derep_num, reverse=True)                

                # combine in order
                alignpe = []                
                for key in order:
                    alignpe.append("\n".join([
                        key, 
                        dalign1[key].replace("\n", "") + "nnnn" + \
                        dalign2[key].replace("\n", "")]))

                ## append aligned cluster string
                aligned.append("\n".join(alignpe).strip())

            # Malformed clust. Dictionary creation with only 1 element 
            except ValueError as inst:
                print("Bad PE cluster - {}\nla1 - {}\nla2 - {}"
                      .format(clust, lines1, lines2)
                      )

            ## Either reads are SE, or at least some pairs are merged.
            except IndexError:
                    
                # limit the number of input seqs
                # use lclust already built before checking pairs
                lclust = "\n".join(clust.split()[:maxseqs * 2])

                # the muscle command with alignment as stdin and // as splitter
                cmd = ("echo -e '{}' | {} -quiet -in - ; echo {}"
                       .format(lclust, ip.bins.muscle, "//\n"))

                ## send cmd to the bash shell (TODO: PIPE could overflow here!)
                proc.stdin.write(cmd.encode())

                ## read the stdout by line until // is reached. This BLOCKS.
                for line in iter(proc.stdout.readline, b'//\n'):
                    align1.append(line.decode())

                ## remove '>' from names, and '\n' from inside long seqs                
                lines = "".join(align1)[1:].split("\n>")

                ## find seed of the cluster and put it on top.
                #seed = [i for i in lines if i.split(";")[-1][0] == "*"][0]
                seed = [i for i in lines if i.split('\n')[0][-1] == "*"][0]
                lines.pop(lines.index(seed))
                lines = [seed] + sorted(
                    lines, key=get_derep_num, reverse=True)

                ## format remove extra newlines from muscle
                aa = [i.split("\n", 1) for i in lines]
                align1 = [i[0] + '\n' + "".join([j.replace("\n", "")
                          for j in i[1:]]) for i in aa]

                # trim edges in sloppy gbs/ezrad data.
                # Maybe relevant to other types too...
                if is_gbs:
                    align1 = gbs_trim(align1)

                ## append to aligned
                aligned.append("\n".join(align1))

    # cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    ## return the aligned clusters
    return aligned


def aligned_indel_filter(clust, max_internal_indels):
    """ checks for too many internal indels in muscle aligned clusters """

    ## make into list
    lclust = clust.split()
    
    ## paired or not
    try:
        seq1 = [i.split("nnnn")[0] for i in lclust[1::2]]
        seq2 = [i.split("nnnn")[1] for i in lclust[1::2]]
        intindels1 = [i.rstrip("-").lstrip("-").count("-") for i in seq1]
        intindels2 = [i.rstrip("-").lstrip("-").count("-") for i in seq2]
        intindels = intindels1 + intindels2
        if max(intindels) > max_internal_indels:
            return 1
    except IndexError:
        seq1 = lclust[1::2]
        intindels = [i.rstrip("-").lstrip("-").count("-") for i in seq1]
        if max(intindels) > max_internal_indels:
            return 1     
    return 0


def get_derep_num(read):
    "return the number of replicates in a derep read"
    return int(read.split("=")[-1].split("\n")[0][:-2])


def gbs_trim(align1):
    """
    No reads can go past the left of the seed, or right of the least extended
    reverse complement match. Example below. m is a match. u is an area where 
    lots of mismatches typically occur. The cut sites are shown.
    
    Original locus*
    Seed           TGCAG************************************-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm------------------------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuuuuuuuu
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu

    Trimmed locus*
    Seed           TGCAG************************************---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmm
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    """
    leftmost = rightmost = None
    dd = {k: v for k, v in [j.rsplit("\n", 1) for j in align1]}
    seed = [i for i in dd.keys() if i.rsplit(";")[-1][0] == "*"][0]
    leftmost = [i != "-" for i in dd[seed]].index(True)
    revs = [i for i in dd.keys() if i.rsplit(";")[-1][0] == "-"]
    if revs:
        subright = max([[i != "-" for i in seq[::-1]].index(True) 
                        for seq in [dd[i] for i in revs]])
    else:
        subright = 0
    rightmost = len(dd[seed]) - subright

    ## if locus got clobbered then print place-holder NNN
    names, seqs = zip(*[i.rsplit("\n", 1) for i in align1])
    if rightmost > leftmost:
        newalign1 = [n + "\n" + i[leftmost:rightmost] 
                     for i, n in zip(seqs, names)]
    else:
        newalign1 = [n + "\nNNN" for i, n in zip(seqs, names)]
    return newalign1


def index_ref_with_bwa(data):
    "Index the reference sequence, unless it already exists"
    # get ref file from params
    refseq_file = data.params.reference_sequence
    index_files = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if not os.path.exists(refseq_file):
        raise IPyradError(
            REQUIRE_REFERENCE_PATH.format(data.params.assembly_method))

    # If reference sequence already exists then bail out of this func
    if all([os.path.isfile(refseq_file + i) for i in index_files]):
        return

    # bwa index <reference_file>
    cmd = [ip.bins.bwa, "index", refseq_file]
    proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=None)
    error = proc.communicate()[1].decode()

    # error handling for one type of error on stderr
    if proc.returncode:
        if "please use bgzip" in error:
            raise IPyradError(NO_ZIP_BINS.format(refseq_file))
        else:
            raise IPyradError(error)


def index_ref_with_sam(data):
    "Index ref for building scaffolds w/ index numbers in steps 5-6"
    # get ref file from params
    refseq_file = data.params.reference_sequence
    if not os.path.exists(refseq_file):
        raise IPyradError(
            REQUIRE_REFERENCE_PATH.format(data.params.assembly_method))

    # If reference index exists then bail out unless force
    if os.path.exists(refseq_file + ".fai"):
        return

    # complain if file is bzipped
    if refseq_file.endswith(".gz"):
        raise IPyradError("You must decompress your genome file.") 

    # index the file
    pysam.faidx(refseq_file)


def mapping_reads(data, sample, nthreads):
    """
    Map reads to reference sequence. This reads in the fasta files
    (samples.files.edits), and maps each read to the reference. Unmapped reads
    are dropped right back in the de novo pipeline.
    Mapped reads end up in a sam file.
    """

    # outfiles
    samout = os.path.join(
        data.dirs.refmapping,
        "{}.sam".format(sample.name))

    bamout = os.path.join(
        data.dirs.refmapping,
        "{}-mapped-sorted.bam".format(sample.name))

    ubamout = os.path.join(
        data.dirs.refmapping,
        "{}-unmapped.bam".format(sample.name))

    ufastqout = os.path.join(
        data.dirs.refmapping,
        "{}-unmapped.fastq".format(sample.name))

    # infiles: dereplicated read1 and read2 files
    if "pair" not in data.params.datatype:
        infiles = [
            os.path.join(data.tmpdir, "{}_derep.fastq".format(sample.name))]
    else:
        inread1 = os.path.join(
            data.tmpdir, "{}_R1-tmp.fastq".format(sample.name))
        inread2 = os.path.join(
            data.tmpdir, "{}_R2-tmp.fastq".format(sample.name))
        infiles = []
        for infile in [inread1, inread2]:
            if os.path.exists(infile):
                infiles.append(infile)

    # (cmd1) bwa mem [OPTIONS] <index_name> <file_name_A> [<file_name_B>]
    #  -t #         : Number of threads
    #  -M           : Mark split alignments as secondary.

    # (cmd2) samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
    #  -b = write to .bam
    #  -q = Only keep reads with mapq score >= 30 (seems to be pretty standard)
    #  -F = Select all reads that DON'T have these flags.
    #        0x4 (segment unmapped)
    #        0x100 (Secondary alignment) 
    #        0x800 (supplementary alignment) (chimeric-like, not likely)
    #        0x71
    #        0xb1 
    #  -U = Write out all reads that don't pass the -F filter
    #       (all unmapped reads go to this file).

    # (cmd3) samtools sort [options...] [in.bam]
    #  -T = Temporary file name, this is required by samtools, ignore it
    #       Here we hack it to be samhandle.tmp cuz samtools cleans it up
    #  -O = Output file format, in this case bam
    #  -o = Output file name

    # (cmd5) samtools bam2fq -v 45 [in.bam]
    #   -v45 set the default qscore arbirtrarily high
    #
    cmd1 = [
        ip.bins.bwa, "mem",
        "-t", str(max(1, nthreads)),
        "-M",
        data.params.reference_sequence,
    ]
    cmd1 += infiles

    # Insert optional flags for bwa
    bwa_args = data.hackersonly.bwa_args.split()
    bwa_args.reverse()
    for arg in bwa_args:
        cmd1.insert(2, arg)

    with open(samout, 'wb') as outfile:
        proc1 = sps.Popen(cmd1, stderr=None, stdout=outfile)
        error1 = proc1.communicate()[0]
        if proc1.returncode:
            raise IPyradError("bwa error: {}".format(error1))

    # sends unmapped reads to a files and will PIPE mapped reads to cmd3
    cmd2 = [
        ip.bins.samtools, "view",
        "-b",
        "-F", "0x904",
        "-U", ubamout,
        samout,
    ]

    # this is gonna catch mapped bam output from cmd2 and write to file
    cmd3 = [
        ip.bins.samtools, "sort",
        "-T", os.path.join(data.dirs.refmapping, sample.name + ".sam.tmp"),
        "-O", "bam",
        "-o", bamout]

    # Later we're gonna use samtools to grab out regions using 'view'
    cmd4 = [ip.bins.samtools, "index", bamout]

    # convert unmapped reads to fastq
    cmd5 = [
        ip.bins.samtools, "bam2fq",
        "-v 45",
        ubamout,
    ]

    # Insert additional arguments for paired data to the commands.
    # We assume Illumina paired end reads for the orientation
    # of mate pairs (orientation: ---> <----).
    if 'pair' in data.params.datatype:
        # add samtools filter for only keep if both pairs hit
        # 0x1 - Read is paired
        # 0x2 - Each read properly aligned
        cmd2.insert(2, "0x3")
        cmd2.insert(2, "-f")

        # tell bam2fq that there are output files for each read pair
        cmd5.insert(2, os.path.join(
            data.tmpdir, sample.name + "-tmp-umap1.fastq"))
        cmd5.insert(2, "-1")
        cmd5.insert(2, os.path.join(
            data.tmpdir, sample.name + "-tmp-umap2.fastq"))
        cmd5.insert(2, "-2")
    else:
        cmd5.insert(2, ufastqout)
        cmd5.insert(2, "-0")


    # cmd2 writes to sname.unmapped.bam and fills pipe with mapped BAM data
    proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)

    # cmd3 pulls mapped BAM from pipe and writes to sname.mapped-sorted.bam
    proc3 = sps.Popen(
        cmd3, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc2.stdout)
    error3 = proc3.communicate()[0]
    if proc3.returncode:
        raise IPyradError(error3)
    proc2.stdout.close()

    # cmd4 indexes the bam file
    proc4 = sps.Popen(cmd4, stderr=sps.STDOUT, stdout=sps.PIPE)
    error4 = proc4.communicate()[0]
    if proc4.returncode:
        raise IPyradError(error4)

    # Running cmd5 writes to either edits/sname-refmap_derep.fastq for SE
    # or it makes edits/sname-tmp-umap{12}.fastq for paired data, which
    # will then need to be merged.
    proc5 = sps.Popen(cmd5, stderr=sps.STDOUT, stdout=sps.PIPE)
    error5 = proc5.communicate()[0]
    if proc5.returncode:
        raise IPyradError(error5)


def check_insert_size(data, sample):
    """
    check mean insert size for this sample and update 
    hackersonly.max_inner_mate_distance if need be. This value controls how 
    far apart mate pairs can be to still be considered for bedtools merging 
    downstream.
    """
    sbam = os.path.join(
        data.dirs.refmapping, 
        "{}-mapped-sorted.bam".format(sample.name)
    )
    stats = pysam.stats(sbam)
    statslines = [i for i in stats.split("\n") if i.startswith("SN")]

    ## starting vals
    avg_insert = 0
    stdv_insert = 0
    avg_len = 0

    ## iterate over results
    for line in statslines:
        if "insert size average" in line:
            avg_insert = float(line.split(":")[-1].strip())

        elif "insert size standard deviation" in line:
            ## hack to fix sim data when stdv is 0.0. Shouldn't
            ## impact real data bcz stdv gets rounded up below
            stdv_insert = float(line.split(":")[-1].strip()) + 0.1
       
        elif "average length" in line:
            avg_len = float(line.split(":")[-1].strip())

    ## If all values return successfully set the max inner mate distance.
    ## This is tricky. avg_insert is the average length of R1+R2+inner mate
    ## distance. avg_len is the average length of a read. If there are lots
    ## of reads that overlap then avg_insert will be close to but bigger than
    ## avg_len. We are looking for the right value for `bedtools merge -d`
    ## which wants to know the max distance between reads. 
    if all([avg_insert, stdv_insert, avg_len]):
        ## If 2 * the average length of a read is less than the average
        ## insert size then most reads DO NOT overlap
        if stdv_insert < 5:
            stdv_insert = 5.
        if (2 * avg_len) < avg_insert:
            hack = avg_insert + (3 * np.ceil(stdv_insert)) - (2 * avg_len)

        ## If it is > than the average insert size then most reads DO
        ## overlap, so we have to calculate inner mate distance a little 
        ## differently.
        else:
            hack = (avg_insert - avg_len) + (3 * np.ceil(stdv_insert))           

        ## set the hackerdict value
        data.hackersonly.max_inner_mate_distance = int(np.ceil(hack))

    else:
        ## If something fsck then set a relatively conservative distance
        data.hackersonly.max_inner_mate_distance = 300


def bedtools_merge(data, sample):
    """
    Get all contiguous genomic regions with one or more overlapping
    reads. This is the shell command we'll eventually run

    bedtools bamtobed -i 1A_0.sorted.bam | bedtools merge [-d 100]
        -i <input_bam>  :   specifies the input file to bed'ize
        -d <int>        :   For PE set max distance between reads
    """
    mappedreads = os.path.join(
        data.dirs.refmapping,
        "{}-mapped-sorted.bam".format(sample.name))

    # command to call `bedtools bamtobed`, and pipe output to stdout
    # Usage:   bedtools bamtobed [OPTIONS] -i <bam>
    # Usage:   bedtools merge [OPTIONS] -i <bam>
    cmd1 = [ip.bins.bedtools, "bamtobed", "-i", mappedreads]
    cmd2 = [ip.bins.bedtools, "merge", "-i", "-"]

    # If PE the -d flag to tell bedtools how far apart to allow mate pairs.
    # If SE the -d flag is negative, specifying that SE reads need to
    # overlap by at least a specific number of bp. This prevents the
    # stairstep syndrome when a + and - read are both extending from
    # the same cutsite. Passing a negative number to `merge -d` gets this done.
    # +++ scrath the above, we now deal with step ladder data
    if 'pair' in data.params.datatype:
        check_insert_size(data, sample)
        cmd2.insert(2, str(data.hackersonly.max_inner_mate_distance))
        cmd2.insert(2, "-d")
    #else:
    #    cmd2.insert(2, str(-1 * data._hackersonly["min_SE_refmap_overlap"]))
    #    cmd2.insert(2, "-d")

    ## pipe output from bamtobed into merge
    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc1.stdout)
    result = proc2.communicate()[0].decode()
    proc1.stdout.close()

    # check for errors and do cleanup
    if proc2.returncode:
        raise IPyradError("error in %s: %s", cmd2, result)

    # Report the number of regions we're returning
    # nregions = len(result.strip().split("\n"))
    return result


def build_clusters_from_cigars(data, sample):
    """
    Directly building clusters relative to reference. Uses the function 
    cigared() to impute indels relative to reference. This means add - for 
    insertion and skip* deletions. Skipping is not a good final solution.
    """
    # get all regions with reads. Generator to yield (str, int, int)
    fullregions = bedtools_merge(data, sample).strip().split("\n")
    regions = (i.split("\t") for i in fullregions)
    regions = ((i, int(j), int(k)) for (i, j, k) in regions)

    # access reads from bam file using pysam
    bamfile = pysam.AlignmentFile(
        os.path.join(
            data.dirs.refmapping,
            "{}-mapped-sorted.bam".format(sample.name)),
        'rb')

    # output path 
    opath = os.path.join(
        data.dirs.clusts, "{}.clustS.gz".format(sample.name))
    out = gzip.open(opath, 'wt')
    idx = 0

    # iterate over all regions to build clusters
    clusters = []
    for reg in regions:
        # uncomment and compare against ref sequence when testing
        # ref = get_ref_region(data.paramsdict["reference_sequence"], *reg)
        reads = bamfile.fetch(*reg)

        # store reads in a dict
        rdict = {}

        # paired-end data cluster building
        if "pair" in data.params.datatype:

            # match paired reads together in a dictionary
            for read in reads:
                if read.qname not in rdict:
                    rdict[read.qname] = [read, None]
                else:
                    rdict[read.qname][1] = read

            # sort keys by derep number
            keys = sorted(
                rdict.keys(),
                key=lambda x: int(x.split("=")[-1]), reverse=True)

            # build the cluster based on map positions, orientation, cigar
            clust = []
            for key in keys:
                r1, r2 = rdict[key]
                if r1 and r2:

                    #lref = len(ref[1])
                    lref = reg[2] - reg[1]
                    arr1 = np.zeros(lref, dtype="U1")
                    arr2 = np.zeros(lref, dtype="U1")
                    arr1.fill("-")
                    arr2.fill("-")

                    # how far ahead of the start does this read begin
                    seq = cigared(r1.seq, r1.cigar)
                    start = r1.reference_start - reg[1] 
                    arr1[start:start + len(seq)] = list(seq)
                    
                    seq = cigared(r2.seq, r2.cigar)
                    start = r2.reference_start - reg[1] 
                    arr2[start:start + len(seq)] = list(seq)
                    
                    arr3 = join_arrays(arr1, arr2)
                    pairseq = "".join(arr3)

                    ori = "+"
                    if r1.is_reverse:
                        ori = "-"
                    derep = r1.qname.split("=")[-1]
                    rname = "{}:{}-{};size={};{}".format(
                        reg[0], reg[1], reg[2], derep, ori)
                    # *reg, derep, ori)
                    clust.append("{}\n{}".format(rname, pairseq))

        # single-end data cluster building
        else:   
            mstart = int(9e12)
            mend = 0

            for read in reads:
                rdict[read.qname] = read
                mstart = min(mstart, read.aend - read.alen)
                mend = max(mend, read.aend)
        
            # sort keys by derep number
            keys = sorted(
                rdict.keys(),
                key=lambda x: int(x.split("=")[-1]), reverse=True)

            # build the cluster based on map positions, orientation, cigar
            clust = []
            for key in keys:
                r1 = rdict[key]

                #aref = np.array(list(ref[1]))
                lref = mend - mstart
                arr1 = np.zeros(lref, dtype="U1")
                arr1.fill("-")

                # how far ahead of the start does this read begin
                seq = cigared(r1.seq, r1.cigar)
                rstart = (r1.aend - r1.alen) - mstart
                arr1[rstart:rstart + len(seq)] = list(seq)
                aseq = "".join(arr1)

                ori = "+"
                if r1.is_reverse:
                    ori = "-"
                derep = r1.qname.split("=")[-1]
                rname = "{}:{}-{};size={};{}".format(
                    reg[0], mstart, mend, derep, ori)
                clust.append("{}\n{}".format(rname, aseq))

        # store this cluster
        if clust:
            clusters.append("\n".join(clust))
            idx += 1

        # if 1000 clusters stored then write to disk
        if not idx % 10000:
            if clusters:
                out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
                clusters = []
 
    # write final remaining clusters to disk
    if clusters:
        out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
    out.close()


def split_endtoend_reads(data, sample):
    """
    Takes R1nnnnR2 derep reads from paired data and splits it back into
    separate R1 and R2 parts for read mapping.
    """
    inp = os.path.join(data.tmpdir, "{}_derep.fastq".format(sample.name))
    out1 = os.path.join(data.tmpdir, "{}_R1-tmp.fastq".format(sample.name))
    out2 = os.path.join(data.tmpdir, "{}_R2-tmp.fastq".format(sample.name))

    splitderep1 = open(out1, 'w')
    splitderep2 = open(out2, 'w')

    with open(inp, 'r') as infile:
        # Read in the infile two lines at a time: (seqname, sequence)
        duo = izip(*[iter(infile)] * 2)

        ## lists for storing results until ready to write
        split1s = []
        split2s = []

        ## iterate over input splitting, saving, and writing.
        idx = 0
        while 1:
            try:
                itera = next(duo)
            except StopIteration:
                break
            ## split the duo into separate parts and inc counter
            part1, part2 = itera[1].split("nnnn")
            idx += 1

            ## R1 needs a newline, but R2 inherits it from the original file
            ## store parts in lists until ready to write
            split1s.append("{}{}\n".format(itera[0], part1))
            split2s.append("{}{}".format(itera[0], part2))

            ## if large enough then write to file
            if not idx % 10000:
                splitderep1.write("".join(split1s))
                splitderep2.write("".join(split2s))
                split1s = []
                split2s = []

    ## write final chunk if there is any
    if any(split1s):
        splitderep1.write("".join(split1s))
        splitderep2.write("".join(split2s))

    ## close handles
    splitderep1.close()
    splitderep2.close()


# DEPRECATED: SLOW
# def get_ref_region(reference, contig, rstart, rend):
#     "returns the reference sequence over a given region"
#     cmd = [
#         ip.bins.samtools, 'faidx',
#         reference,
#         "{}:{}-{}".format(contig, rstart + 1, rend),
#     ]
#     stdout = sps.Popen(cmd, stdout=sps.PIPE).communicate()[0]
#     name, seq = stdout.decode().split("\n", 1)
#     listseq = [name, seq.replace("\n", "")]
#     return listseq


def join_arrays(arr1, arr2):
    "join read1 and read2 arrays and resolve overlaps"
    arr3 = np.zeros(arr1.size, dtype="U1")
    for i in range(arr1.size):

        if arr1[i] == arr2[i]:
            arr3[i] = arr1[i]

        elif arr1[i] == "N":
            if arr2[i] == "-":
                arr3[i] = "N"
            else:
                arr3[i] = arr2[i]

        elif arr2[i] == "N":
            if arr1[i] == "-":
                arr3[i] = "N"
            else:
                arr3[i] = arr1[i]

        elif arr1[i] == "-":
            if arr2[i] == "N":
                arr3[i] = "N"
            else:
                arr3[i] = arr2[i]

        elif arr2[i] == "-":
            if arr1[i] == "N":
                arr3[i] = "N"
            else:
                arr3[i] = arr1[i]

        else:
            arr3[i] = "N"
    return arr3


def cigared(sequence, cigartups):
    "modify sequence based on its cigar string"
    start = 0
    seq = ""
    for tup in cigartups:
        flag, add = tup
        if flag == 0:
            seq += sequence[start:start + add]
        if flag == 1:
            pass
        if flag == 2:
            seq += "-" * add
            start -= add
        if flag == 4:
            pass
        start += add
    return seq


def get_quick_depths(data, sample):
    "iterate over clustS files to get data"

    ## use existing sample cluster path if it exists, since this
    ## func can be used in step 4 and that can occur after merging
    ## assemblies after step3, and if we then referenced by data.dirs.clusts
    ## the path would be broken.
    if not sample.files.clusters:
        sample.files.clusters = os.path.join(
            data.dirs.clusts,
            "{}.clustS.gz".format(sample.name))

    try:
        # get new clustered loci
        with gzip.open(sample.files.clusters, 'rt') as infile:
            pairdealer = izip(*[iter(infile)] * 2)

            ## storage
            depths = []
            maxlen = []

            ## start with cluster 0
            tdepth = 0
            tlen = 0

            ## iterate until empty
            while 1:
                ## grab next
                try:
                    name, seq = next(pairdealer)
                except StopIteration:
                    break

                # if not the end of a cluster
                if name.strip() == seq.strip():
                    depths.append(tdepth)
                    maxlen.append(tlen)
                    tlen = 0
                    tdepth = 0

                else:
                    tdepth += int(name.strip().split("=")[-1][:-2])
                    tlen = len(seq)
    except TypeError:
        raise IPyradError(
            "error in get_quick_depths(): {}".format(sample.files.clusters))

    # return
    return np.array(maxlen), np.array(depths)


def store_sample_stats(data, sample, maxlens, depths):
    "stats, cleanup, and link to samples"

    # Test if depths is non-empty, but just full of zeros.
    if not depths.max():
        print("    no clusters found for {}".format(sample.name))
        return

    else:
        # store which min was used to calculate hidepth here
        sample.stats_dfs.s3["hidepth_min"] = data.params.mindepth_majrule

        # If our longest sequence is longer than the current max_fragment_len
        # then update max_fragment_length. For assurance we require that
        # max len is 4 greater than maxlen, to allow for pair separators.
        hidepths = depths >= data.params.mindepth_majrule
        maxlens = maxlens[hidepths]

        # Handle the case where there are no hidepth clusters
        if maxlens.any():
            maxlen = int(maxlens.mean() + (2. * maxlens.std()))
        else:
            maxlen = 0
        if maxlen > data.hackersonly.max_fragment_length:
            data.hackersonly.max_fragment_length = int(maxlen + 4)

        # make sense of stats
        keepmj = depths[depths >= data.params.mindepth_majrule]
        keepstat = depths[depths >= data.params.mindepth_statistical]

        # sample summary stat assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = int(depths.shape[0])
        sample.stats["clusters_hidepth"] = int(keepmj.shape[0])

        # store depths histogram as a dict. Limit to first 25 bins
        bars, bins = np.histogram(depths, bins=range(1, 26))
        sample.depths = {int(i): int(v) for i, v in zip(bins, bars) if v}

        # sample stat assignments
        # Trap numpy warnings ("mean of empty slice") for samps w/ few reads
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
            sample.stats_dfs.s3["clusters_total"] = int(depths.shape[0])
            try:
                sample.stats_dfs.s3["clusters_hidepth"] = (
                    int(sample.stats["clusters_hidepth"]))
            except ValueError:
                # Handle clusters_hidepth == NaN
                sample.stats_dfs.s3["clusters_hidepth"] = 0
            sample.stats_dfs.s3["avg_depth_total"] = float(depths.mean())
            sample.stats_dfs.s3["avg_depth_mj"] = float(keepmj.mean())
            sample.stats_dfs.s3["avg_depth_stat"] = float(keepstat.mean())
            sample.stats_dfs.s3["sd_depth_total"] = float(depths.std())
            sample.stats_dfs.s3["sd_depth_mj"] = float(keepmj.std())
            sample.stats_dfs.s3["sd_depth_stat"] = float(keepstat.std())

    # store results
    # If PE, samtools reports the _actual_ number of reads mapped, both
    # R1 and R2, so here if PE divide the results by 2 to stay consistent
    # with how we've been reporting R1 and R2 as one "read pair"
    if "pair" in data.params.datatype:
        sample.stats["refseq_mapped_reads"] = sum(depths)
        sample.stats["refseq_unmapped_reads"] = int(
            sample.stats.reads_passed_filter - \
            sample.stats["refseq_mapped_reads"])

    # cleanup
    if not data.params.assembly_method == "denovo":    
        unmapped = os.path.join(
            data.dirs.refmapping, 
            sample.name + "-unmapped.bam")
        samplesam = os.path.join(
            data.dirs.refmapping, 
            sample.name + ".sam")
        for rfile in [unmapped, samplesam]:
            if os.path.exists(rfile):
                os.remove(rfile)

    # if loglevel==DEBUG
    # log_level = ip.logger.getEffectiveLevel()
    # if log_level != 10:
        ## Clean up loose files only if not in DEBUG
        ##- edits/*derep, utemp, *utemp.sort, *htemp, *clust.gz
    derepfile = os.path.join(data.dirs.edits, sample.name + "_derep.fastq")
    mergefile = os.path.join(data.dirs.edits, sample.name + "_merged_.fastq")
    uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
    usort = os.path.join(data.dirs.clusts, sample.name + ".utemp.sort")
    hhandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")
    clusters = os.path.join(data.dirs.clusts, sample.name + ".clust.txt")

    for rfile in [derepfile, mergefile, uhandle, usort, hhandle, clusters]:
        if os.path.exists(rfile):
            os.remove(rfile)


def declone_3rad(data, sample):
    """
    3rad uses random adapters to identify pcr duplicates. We will
    remove pcr dupes here. Basically append the radom adapter to
    each sequence, do a regular old vsearch derep, then trim
    off the adapter, and push it down the pipeline. This will
    remove all identical seqs with identical random i5 adapters.
    """

    # tmp_outfile = tempfile.NamedTemporaryFile(
    #     mode='wb',
    #     delete=False,
    #     dir=data.dirs.edits,
    #     suffix="_append_adapters_.fastq",
    # )

    try:
        ## Remove adapters from head of sequence and write out
        ## tmp_outfile is now the input file for the next step
        ## first vsearch derep discards the qscore so we iterate
        ## by pairs
        with open(tmp_outfile.name) as infile:
            with open(os.path.join(
                data.dirs.edits, 
                sample.name + "_declone.fastq"), 'wb') as outfile:
                duo = izip(*[iter(infile)] * 2)

                ## a list to store until writing
                writing = []
                counts2 = 0

                while 1:
                    try:
                        read = next(duo)
                    except StopIteration:
                        break

                    ## Peel off the adapters. There's probably a faster
                    ## way of doing this.
                    writing.append("\n".join([
                        read[0].strip(),
                        read[1].strip()[8:]]
                    ))

                    ## Write the data in chunks
                    counts2 += 1
                    if not counts2 % 1000:
                        outfile.write("\n".join(writing) + "\n")
                        writing = []
                if writing:
                    outfile.write("\n".join(writing))
                    outfile.close()

    except Exception as inst:
        raise IPyradError(
            "    Caught error while decloning 3rad data - {}".format(inst))

    finally:
        ## failed samples will cause tmp file removal to raise.
        ## just ignore it.
        try:
            ## Clean up temp files
            if os.path.exists(tmp_outfile.name):
                os.remove(tmp_outfile.name)
            if os.path.exists(tmp_outfile.name):
                os.remove(tmp_outfile.name)
        except Exception as inst:
            pass


# globals
NO_ZIP_BINS = """
  Reference sequence must be de-compressed fasta or bgzip compressed,
  your file is probably gzip compressed. The simplest fix is to gunzip
  your reference sequence by running this command:
 
      gunzip {}

  Then edit your params file to remove the `.gz` from the end of the
  path to your reference sequence file and rerun step 3 with the `-f` flag.

      error {}
  """
REQUIRE_REFERENCE_PATH = """\
  Assembly method {} requires that you enter a 'reference_sequence_path'.
"""
