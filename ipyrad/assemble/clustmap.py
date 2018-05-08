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
import tempfile
import warnings
import numpy as np
import ipyrad as ip
import subprocess as sps

from .refmap import refmap_init, refmap_stats
from .refmap import mapreads, ref_muscle_chunker, ref_build_and_muscle_chunk
from .util import IPyradError, IPyradWarningExit
from .util import bcomp, comp

import logging
LOGGER = logging.getLogger(__name__)


class Step3:
    def __init__(self, data, samples, noreverse, maxindels, force, ipyclient):

        # store attributes
        self.data = data
        self.noreverse = noreverse
        self.maxindels = maxindels
        self.force = force
        self.ipyclient = ipyclient
        self.gbs = bool("gbs" in self.data.paramsdict["datatype"])
        self.samples = self.check_samples(samples)

        # init funcs
        self.setup_dirs()
        self.refmap_init()
        self.tune_threads()
        self.tune_load_balancer()


    def run(self):
        # jobs to be run on the cluster
        if self.data.paramsdict["assembly_method"] == "denovo":
            self.remote_run_dereps()
            self.remote_run_cluster_map_build()
            self.remote_run_align_cleanup()

        elif self.data.paramsdict["assembly_method"] == "denovo+reference":
            self.remote_run_dereps()
            self.remote_run_cluster_map_build()
            self.remote_run_align_cleanup()

        elif self.data.paramsdict["assembly_method"] == "denovo-reference":
            self.remote_run_dereps()
            self.remote_run_cluster_map_build()
            self.remote_run_align_cleanup()

        elif self.data.paramsdict["assembly_method"] == "reference":
            self.remote_run_dereps()
            self.remote_run_cluster_map_build()
            self.remote_run_align_cleanup()
        self.cleanup()  

    ## init functions ------------------------------------
    def check_samples(self, samples):
        ## list of samples to submit to queue
        subsamples = []

        ## if sample is already done skip
        for sample in samples:
            if sample.stats.state < 2:
                print("Sample not ready for clustering. First run step2 on: {}"
                      .format(sample.name))
                continue

            if not self.force:
                if sample.stats.state >= 3:
                    print("Skipping {}; aleady clustered. Use force to re-cluster"
                          .format(sample.name))
                else:
                    if sample.stats.reads_passed_filter:
                        subsamples.append(sample)
            else:
                ## force to overwrite
                if sample.stats.reads_passed_filter:
                    subsamples.append(sample)

        ## run subsamples
        if not subsamples:
            raise IPyradError(
                "No Samples ready to be clustered. First run step 2.")
        return subsamples


    def setup_dirs(self):
        # make output folder for clusters
        pdir = os.path.realpath(self.data.paramsdict["project_dir"])
        self.data.dirs.clusts = os.path.join(
            pdir, "{}_clust_{}"
            .format(self.data.name, self.data.paramsdict["clust_threshold"]))
        if not os.path.exists(self.data.dirs.clusts):
            os.mkdir(self.data.dirs.clusts)

        # make a tmpdir for align files
        self.data.tmpdir = os.path.abspath(os.path.expanduser(
            os.path.join(pdir, self.data.name + '-tmpalign')))
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)

        # If ref mapping, init samples and make the refmapping output directory.
        if not self.data.paramsdict["assembly_method"] == "denovo":
            # make output directory for read mapping process
            self.data.dirs.refmapping = os.path.join(
                pdir, "{}_refmapping".format(self.data.name))
            if not os.path.exists(self.data.dirs.refmapping):
                os.mkdir(self.data.dirs.refmapping)


    def refmap_init(self):
        # if refmapping make filehandles that will be persistent
        if not self.data.paramsdict["assembly_method"] == "denovo":
            for sample in self.samples:
                refmap_init(self.data, sample, self.force)

            # set thread-count to 2 for paired-data
            self.nthreads = 2
        
        # set thread-count to 1 for single-end data          
        else:
            self.nthreads = 1    


    def tune_threads(self):
        # overwrite nthreads if value in _ipcluster dict
        if "threads" in self.data._ipcluster.keys():
            self.nthreads = int(self.data._ipcluster["threads"])
            # if more CPUs than there are samples then increase threads
            _ncpus = len(self.ipyclient)
            if _ncpus > 2 * len(self.data.samples):
                self.nthreads *= 2


    def tune_load_balancer(self):
        # TODO: for HPC this should make sure targets are spread on diff nodes.
        eids = self.ipyclient.ids
        if self.nthreads:
            if self.nthreads < len(self.ipyclient.ids):
                thview = self.ipyclient.load_balanced_view(
                    targets=eids[::self.nthreads])
            elif self.nthreads == 1:
                thview = self.ipyclient.load_balanced_view()
            else:
                if len(self.ipyclient) > 40:
                    thview = self.ipyclient.load_balanced_view(
                        targets=eids[::4])
                else:
                    thview = self.ipyclient.load_balanced_view(
                        targets=eids[::2])
        self.lbview = thview


    def cleanup(self):
        """ cleanup / statswriting function for Assembly obj """
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

    ## run functions ---------------------------------------
    def remote_run_dereps(self):
        # submit job
        start = time.time()
        printstr = ("dereplicating     ", "s3")
        rasyncs = {}
        for sample in self.samples:
            rasyncs[sample.name] = self.lbview.apply(
                derep_sort_map, 
                *(self.data, sample, self.force, self.nthreads)
                )

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


    def remote_run_cluster_map_build(self):
        # submit clustering/mapping job
        start = time.time()
        casyncs = {}
        for sample in self.samples:
            casyncs[sample.name] = self.lbview.apply(
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
        printstr = ("clustering/mapping", "s3")
        while 1:
            ready = [casyncs[i].ready() for i in casyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        print("")
        for job in casyncs:
            if not casyncs[job].successful():
                raise IPyradError(casyncs[job].exception())

        # track job progress
        start = time.time()
        printstr = ("building clusters ", "s3")
        while 1:
            ready = [basyncs[i].ready() for i in basyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        print("")
        for job in basyncs:
            if not basyncs[job].successful():
                raise IPyradError(basyncs[job].exception())

        # track job progress
        start = time.time()        
        printstr = ("chunking clusters ", "s3")
        while 1:
            ready = [hasyncs[i].ready() for i in hasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        print("")
        for job in hasyncs:
            if not hasyncs[job].successful():
                raise IPyradError(hasyncs[job].exception())


    def remote_run_align_cleanup(self):

        # submit ten aligning jobs for each sample
        start = time.time()
        aasyncs = {}
        for sample in self.samples:
            aasyncs[sample.name] = []
            for idx in range(10):
                handle = os.path.join(self.data.tmpdir, 
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
        printstr = ("aligning clusters ", "s3")
        while 1:
            ready = [i.ready() for i in allasyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        print("")
        for job in allasyncs:
            if not job.successful():
                raise IPyradError(job.exception())

        # track job 2 progress
        start = time.time()
        printstr = ("concat clusters   ", "s3")
        while 1:
            ready = [basyncs[i].ready() for i in basyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break
        print("")
        for job in basyncs:
            if not basyncs[job].successful():
                raise IPyradError(basyncs[job].exception())

        # store stats results and cleanup
        self.data.stats_dfs.s3 = self.data._build_stat("s3")
        for sample in self.samples:

            # get bad align stats
            bad_aligns = sum([int(i.get()) for i in aasyncs[sample.name]])
            self.data.stats_dfs.s3.filtered_bad_align = bad_aligns

            # purge temp files
            _sample_cleanup(self.data, sample)


### step 3.1 funcs
def derep_sort_map(data, sample, force, nthreads):
    """
    Read carefully, this functions acts very differently depending on 
    datatype. Refmaps, then merges, then dereplicates, then denovo 
    clusters reads.
    """

    ## CONCAT FILES FOR MERGED ASSEMBLIES   
    # returns sample.files.edits if nothing to concat
    sample.concatfiles = concat_multiple_edits(data, sample)

    ## PAIRED DATA ONLY:
    ## Denovo: merge and concat fastq pairs [sample.files.pairs]
    ## Reference: only concat fastq pairs  []
    ## Denovo + Reference: merge & concat non mapped reads?
    # reads in concatfiles and writes to mergefiles
    sample.mergedfile = merge_pairs(data, sample, 1, 1)  

    ## 3rad uses random adapters to identify pcr duplicates. We will
    ## remove pcr dupes here. Basically append the random adapter to
    ## each sequence, do a regular old vsearch derep, then trim
    ## off the adapter, and push it down the pipeline. This will
    ## remove all identical seqs with identical random i5 adapters.
    if "3rad" not in data.paramsdict["datatype"]:
        new_derep_and_sort(
            data,
            sample.mergedfile,
            os.path.join(data.tmpdir, sample.name + "_derep.fastq"),
            nthreads)
    else:
        declone_3rad(data, sample)
        new_derep_and_sort(data,
            os.path.join(data.dirs.edits, sample.name + "_declone.fastq"),
            os.path.join(data.tmpdir, sample.name + "_derep.fastq"),
            nthreads)


def new_derep_and_sort(data, infile, outfile, nthreads):
    """
    Dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file. Paired
    reads are dereplicated as one concatenated read and later split again.
    Updated this function to take infile and outfile to support the double
    dereplication that we need for 3rad (5/29/15 iao).
    """
    ## datatypes options
    strand = "plus"
    if data.paramsdict["datatype"] is ('gbs' or '2brad'):
        strand = "both"

    ## do dereplication with vsearch
    cmd = [
        ip.bins.vsearch,
        "--derep_fulllength", infile,
        "--strand", strand,
        "--output", outfile,
        "--threads", str(nthreads),
        "--fasta_width", str(0),
        "--fastq_qmax", "1000",
        "--sizeout", 
        "--relabel_md5",
        ]
    ip.logger.info("derep cmd %s", cmd)

    ## build PIPEd job
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    errmsg = proc.communicate()[0]
    if proc.returncode:
        ip.logger.error("error inside derep_and_sort %s", errmsg)
        raise IPyradWarningExit(errmsg)


def declone_3rad(data, sample):
    """
    3rad uses random adapters to identify pcr duplicates. We will
    remove pcr dupes here. Basically append the radom adapter to
    each sequence, do a regular old vsearch derep, then trim
    off the adapter, and push it down the pipeline. This will
    remove all identical seqs with identical random i5 adapters.
    """
    ip.logger.info("Entering declone_3rad - {}".format(sample.name))

    ## Append i5 adapter to the head of each read. Merged file is input, and
    ## still has fq qual score so also have to append several qscores for the
    ## adapter bases. Open the merge file, get quarts, go through each read
    ## and append the necessary stuff.
    adapter_seqs_file = tempfile.NamedTemporaryFile(
        mode='wb',
        delete=False,
        dir=data.dirs.edits,
        suffix="_append_adapters_.fastq"
        )

    try:
        with open(sample.files.edits[0][0]) as infile:
            quarts = izip(*[iter(infile)] * 4)

            ## a list to store until writing
            writing = []
            counts = 0

            while 1:
                try:
                    read = next(quarts)
                except StopIteration:
                    break

                ## Split on +, get [1], split on "_" (can be either _r1 or
                ## _m1 if merged reads) and get [0] for the i5
                ## prepend "EEEEEEEE" as qscore for the adapters
                i5 = read[0].split("+")[1].split("_")[0]

                ## If any non ACGT in the i5 then drop this sequence
                if 'N' in i5:
                    continue
                writing.append("\n".join([
                    read[0].strip(),
                    i5 + read[1].strip(),
                    read[2].strip(),
                    "E" * 8 + read[3].strip()]
                ))

                ## Write the data in chunks
                counts += 1
                if not counts % 1000:
                    adapter_seqs_file.write("\n".join(writing)+"\n")
                    writing = []
            if writing:
                adapter_seqs_file.write("\n".join(writing))
                adapter_seqs_file.close()

        tmp_outfile = tempfile.NamedTemporaryFile(mode='wb',
                                        delete=False,
                                        dir=data.dirs.edits,
                                        suffix="_decloned_w_adapters_.fastq")

        ## Close the tmp file bcz vsearch will write to it by name, then
        ## we will want to reopen it to read from it.
        tmp_outfile.close()
        ## Derep the data (adapters+seq)
        derep_and_sort(data, adapter_seqs_file.name,
                       os.path.join(data.dirs.edits, tmp_outfile.name), 2)

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

        ip.logger.info("Removed pcr duplicates from {} - {}"
                    .format(sample.name, counts - counts2))

    except Exception as inst:
        raise IPyradError(
            "    Caught error while decloning 3rad data - {}".format(inst))

    finally:
        ## failed samples will cause tmp file removal to raise.
        ## just ignore it.
        try:
            ## Clean up temp files
            if os.path.exists(adapter_seqs_file.name):
                os.remove(adapter_seqs_file.name)
            if os.path.exists(tmp_outfile.name):
                os.remove(tmp_outfile.name)
        except Exception as inst:
            pass


def concat_multiple_edits(data, sample):
    if len(sample.files.edits) < 2:
        return sample.files.edits
    else:
        # cat all inputs; index 0 b/c files are in tuples for r1, r2
        cmd1 = ["cat"] + [i[0] for i in sample.files.edits]

        # write to new concat handle
        conc1 = os.path.join(
            data.dirs.edits, sample.name + "_R1_concatedit.fq.gz")
        with open(conc1, 'w') as cout1:
            proc1 = sps.Popen(cmd1, 
                stderr=sps.STDOUT, stdout=cout1, close_fds=True)
            res1 = proc1.communicate()[0]
        if proc1.returncode:
            raise IPyradWarningExit("error in: %s, %s", cmd1, res1)

        ## Only set conc2 if R2 actually exists
        conc2 = 0
        if os.path.exists(str(sample.files.edits[0][1])):
            cmd2 = ["cat"] + [i[1] for i in sample.files.edits]
            conc2 = os.path.join(
                data.dirs.edits, sample.name + "_R2_concatedit.fq.gz")
            with gzip.open(conc2, 'w') as cout2:
                proc2 = sps.Popen(cmd2, 
                    stderr=sps.STDOUT, stdout=cout2, close_fds=True)
                res2 = proc2.communicate()[0]
            if proc2.returncode:
                raise IPyradWarningExit("error in: %s, %s", cmd2, res2)

        ## store new file handles
        return [(conc1, conc2)]


### step 3.1 subfuncs
def merge_pairs(data, sample, revcomp, vsearch_merge):
    """
    Merge PE reads for denovo paired data. This uses vsearch to check for 
    merged pairs, and for those that do not merge it concatenates the pairs
    together with a 'nnnn' separator. 
    """

    # returns the file that should be derep'd
    if 'pair' not in data.paramsdict['datatype']:
        return sample.concatfiles[0][0]

    if "reference" in data.paramsdict["assembly_method"]:
        return sample.concatfiles

    sample.mergedfile = os.path.join(
        data.tmpdir, sample.name + "_merged.fastq")
    
    # if vsearch_merge then catch nonmerged in a separate file
    if not vsearch_merge:
        return sample.concatfiles
    else:       
        nonmerged1 = os.path.join(
            data.tmpdir, sample.name + "_nonmerged_R1_.fastq")
        nonmerged2 = os.path.join(
            data.tmpdir, sample.name + "_nonmerged_R2_.fastq")

    ## get the maxn and minlen values
    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']
    minlen = str(max(32, data.paramsdict["filter_min_trim_len"]))

    # vsearch merge can now take gzipped files (v.2.8)
    if vsearch_merge:
        cmd = [
            ip.bins.vsearch,
            "--fastq_mergepairs", sample.concatfiles[0][0],
            "--reverse", sample.concatfiles[0][1],
            "--fastqout", sample.mergedfile,
            "--fastqout_notmerged_fwd", nonmerged1,
            "--fastqout_notmerged_rev", nonmerged2,
            "--fasta_width", "0",
            "--fastq_minmergelen", minlen,
            "--fastq_maxns", str(maxn),
            "--fastq_minovlen", "20",
            "--fastq_maxdiffs", "4",
            "--label_suffix", "_m1",
            "--fastq_qmax", "1000",
            "--threads", "2",
            "--fastq_allowmergestagger",
            ]

        LOGGER.debug("merge cmd: %s", " ".join(cmd))
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        res = proc.communicate()[0]
        LOGGER.info(res.decode())
        if proc.returncode:
            LOGGER.error("Error: %s %s", cmd, res)
            raise IPyradWarningExit("Error merge pairs:\n %s\n%s", cmd, res)

    # record how many read pairs were merged
    with open(sample.mergedfile, 'r') as tmpf:
        nmerged = sum(1 for i in tmpf.readlines()) // 4
        sample.stats.reads_merged = nmerged

    # Combine the unmerged pairs and append to the merge file
    with open(sample.mergedfile, 'a') as combout:

        # read in paired end read files 4 lines at a time
        if nonmerged1.endswith(".gz"):
            fr1 = gzip.open(nonmerged1, 'rb')
        else:
            fr1 = open(nonmerged1, 'rb')
        quart1 = izip(*[iter(fr1)] * 4)
        if nonmerged2.endswith(".gz"):
            fr2 = gzip.open(nonmerged2, 'rb')
        else:
            fr2 = open(nonmerged2, 'rb')
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
            if not counts % 50:
                combout.write(b"".join(writing).decode())
                writing = []

        if writing:
            combout.write(b"".join(writing).decode())

    ## close handles
    fr1.close()
    fr2.close()
    combout.close()
    return sample.mergedfile


### step 3.2 funcs
def cluster(data, sample, nthreads, force):
    """
    Calls vsearch for clustering. cov varies by data type, values were chosen
    based on experience, but could be edited by users
    """

    # get the dereplicated reads
    if "reference" not in data.paramsdict["assembly_method"]:
        derephandle = os.path.join(data.tmpdir, sample.name + "_derep.fastq")
    else:
        derephandle = os.path.join(
            data.tmpdir, sample.name + "-refmap_derep.fastq")
        
        # In the event all reads for all samples map successfully then 
        # clustering the unmapped reads makes no sense, so just bail out.
        if not os.stat(derephandle).st_size:
            # In this case you do have to create empty, dummy vsearch output
            # files so building_clusters will not fail.
            uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
            usort = os.path.join(data.dirs.clusts, sample.name + ".utemp.sort")
            hhandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")
            for f in [uhandle, usort, hhandle]:
                open(f, 'a').close()
            return
    assert os.path.exists(derephandle), "bad derep handle: {}".format(derephandle)

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
    cov = 0.75
    minsl = 0.5
    if data.paramsdict["datatype"] in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.5
        minsl = 0.5
    elif data.paramsdict["datatype"] == 'pairgbs':
        strand = "both"
        cov = 0.75
        minsl = 0.75

    ## If this value is not null (which is the default) then override query cov
    if data._hackersonly["query_cov"]:
        cov = str(data._hackersonly["query_cov"])
        assert float(cov) <= 1, "query_cov must be <= 1.0"

    ## get call string
    cmd = [ip.bins.vsearch,
           "-cluster_smallmem", derephandle,
           "-strand", strand,
           "-query_cov", str(cov),
           "-id", str(data.paramsdict["clust_threshold"]),
           "-minsl", str(minsl),
           "-userout", uhandle,
           "-userfields", "query+target+id+gaps+qstrand+qcov",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-threads", str(nthreads),
           "-notmatched", temphandle,
           "-fasta_width", "0",
           "-fastq_qmax", "100",
           "-fulldp",
           "-usersort"]

    ## run vsearch
    ip.logger.debug("%s", cmd)
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    res = proc.communicate()[0]

    # check for errors
    if proc.returncode:
        ip.logger.error("error %s: %s", cmd, res)
        raise IPyradWarningExit("cmd {}: {}".format(cmd, res))


def build_clusters(data, sample, maxindels):
    """
    Combines information from .utemp and .htemp files to create .clust files,
    which contain un-aligned clusters. Hits to seeds are only kept in the
    cluster if the number of internal indels is less than 'maxindels'.
    By default, we set maxindels=6 for this step (within-sample clustering).
    """

    ## If reference assembly then here we're clustering the unmapped reads
    if "reference" in data.paramsdict["assembly_method"]:
        derepfile = os.path.join(
            data.tmpdir, sample.name + "-refmap_derep.fastq")
    else:
        derepfile = os.path.join(data.tmpdir, sample.name + "_derep.fastq")

    ## i/o vsearch files
    uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
    usort = os.path.join(data.dirs.clusts, sample.name + ".utemp.sort")
    hhandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")

    ## create an output file to write clusters to
    sample.files.clusters = os.path.join(
        data.dirs.clusts, sample.name + ".clust.gz")
    clustsout = gzip.open(sample.files.clusters, 'wt')

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
                            int(x.split(";size=")[-1].split("\n")[0][:-1]), 
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
                fseqs.append(">{}*\n{}".format(seed, alldereps[seed]))
                lastseed = seed

            ## add match to the seed
            ## revcomp if orientation is reversed (comp preserves nnnn)
            if ori == "-":
                seq = comp(alldereps[hit])[::-1]
            else:
                seq = alldereps[hit]
            ## only save if not too many indels
            if int(ind) <= maxindels:
                fseqs.append(">{}{}\n{}".format(hit, ori, seq))
            else:
                ip.logger.info("filtered by maxindels: %s %s", ind, seq)

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
                seqlist.append("{}*\n{}".format(nnn, alldereps[nnn[1:]]))
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
    ## log our location for debugging
    ip.logger.info("inside muscle_chunker")

    ## only chunk up denovo data, refdata has its own chunking method which 
    ## makes equal size chunks, instead of uneven chunks like in denovo
    if data.paramsdict["assembly_method"] != "reference":
        ## get the number of clusters
        clustfile = os.path.join(data.dirs.clusts, sample.name + ".clust.gz")
        with iter(gzip.open(clustfile, 'rt')) as clustio:
            nloci = sum(1 for i in clustio if "//" in i) // 2
            #tclust = clustio.read().count("//")//2
            optim = (nloci // 20) + (nloci % 20)
            ip.logger.info("optim for align chunks: %s", optim)

        ## write optim clusters to each tmp file
        clustio = gzip.open(clustfile, 'rt')
        inclusts = iter(clustio.read().strip().split("//\n//\n"))
        
        ## splitting loci so first file is smaller and last file is bigger
        inc = optim // 10
        for idx in range(10):
            ## how big is this chunk?
            this = optim + (idx * inc)
            left = nloci - this
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

    ## data are already chunked, read in the whole thing. bail if no data.
    try:
        with open(handle, 'rb') as infile:
            clusts = infile.read().decode().split("//\n//\n")
            # remove any empty spots
            clusts = [i for i in clusts if i]
            # Skip entirely empty chunks
            if not clusts:
                raise IPyradError("no clusters in file: {}".format(handle))

    except (IOError, IPyradError):
        LOGGER.debug("skipping empty chunk - {}".format(handle))
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
        if not _aligned_indel_filter(clust, max_internal_indels):
            refined.append(clust)
        else:
            highindels += 1

    ## write to file after
    if refined:
        outhandle = handle.rsplit(".", 1)[0] + ".aligned"
        with open(outhandle, 'wb') as outfile:
            outfile.write(str.encode("\n//\n//\n".join(refined) + "\n"))

    ## remove the old tmp file
    if not LOGGER.getEffectiveLevel() == 10:
        os.remove(handle)
    return highindels


def reconcat(data, sample):
    """ takes aligned chunks (usually 10) and concatenates them """

    ## get chunks
    chunks = glob.glob(os.path.join(data.tmpdir,
             sample.name + "_chunk_[0-9].aligned"))

    ## sort by chunk number, cuts off last 8 =(aligned)
    chunks.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))
    ip.logger.info("chunk %s", chunks)

    ## concatenate finished reads
    sample.files.clusters = os.path.join(
        data.dirs.clusts, sample.name + ".clustS.gz")

    ## reconcats aligned clusters
    with gzip.open(sample.files.clusters, 'wb') as out:
        for fname in chunks:
            with open(fname) as infile:
                dat = infile.read().strip()
                out.write(str.encode(dat + "\n//\n//\n"))                   
            os.remove(fname)


### step 3.2 subfuncs
def persistent_popen_align3(clusts, maxseqs=200, is_gbs=False):
    """ keeps a persistent bash shell open and feeds it muscle alignments """

    ## create a separate shell for running muscle in, this is much faster
    ## than spawning a separate subprocess for each muscle call
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
                    keys, key=_get_derep_num, reverse=True)                

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
                ip.logger.debug(
                    "Bad PE cluster - {}\nla1 - {}\nla2 - {}"
                    .format(clust, lines1, lines2))

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
                    lines, key=_get_derep_num, reverse=True)

                ## format remove extra newlines from muscle
                aa = [i.split("\n", 1) for i in lines]
                align1 = [i[0] + '\n' + "".join([j.replace("\n", "") 
                          for j in i[1:]]) for i in aa]
                
                # trim edges in sloppy gbs/ezrad data. 
                # Maybe relevant to other types too...
                if is_gbs:
                    align1 = _gbs_trim(align1)

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


def get_quick_depths(data, sample):
    """ iterate over clustS files to get data """

    ## use existing sample cluster path if it exists, since this
    ## func can be used in step 4 and that can occur after merging
    ## assemblies after step3, and if we then referenced by data.dirs.clusts
    ## the path would be broken.
    if sample.files.clusters:
        pass
    else:
        ## set cluster file handles
        sample.files.clusters = os.path.join(
            data.dirs.clusts, sample.name + ".clustS.gz")

    ## get new clustered loci
    fclust = data.samples[sample.name].files.clusters
    clusters = gzip.open(fclust, 'rt')
    pairdealer = izip(*[iter(clusters)] * 2)

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

        ## if not the end of a cluster
        #print name.strip(), seq.strip()
        if name.strip() == seq.strip():
            depths.append(tdepth)
            maxlen.append(tlen)
            tlen = 0
            tdepth = 0

        else:
            tdepth += int(name.strip().split("=")[-1][:-1])
            tlen = len(seq)

    ## return
    clusters.close()
    return np.array(maxlen), np.array(depths)


def _sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """

    # get maxlen and depths array from clusters
    maxlens, depths = get_quick_depths(data, sample)

    try:
        depths.max()
    except ValueError:
        ## If depths is an empty array max() will raise
        print("    no clusters found for {}".format(sample.name))
        return

    ## Test if depths is non-empty, but just full of zeros.
    if depths.max():
        ## store which min was used to calculate hidepth here
        sample.stats_dfs.s3["hidepth_min"] = (
            data.paramsdict["mindepth_majrule"])

        # If our longest sequence is longer than the current max_fragment_len
        # then update max_fragment_length. For assurance we require that
        # max len is 4 greater than maxlen, to allow for pair separators.
        hidepths = depths >= data.paramsdict["mindepth_majrule"]
        maxlens = maxlens[hidepths]

        ## Handle the case where there are no hidepth clusters
        if maxlens.any():
            maxlen = int(maxlens.mean() + (2. * maxlens.std()))
        else:
            maxlen = 0
        if maxlen > data._hackersonly["max_fragment_length"]:
            data._hackersonly["max_fragment_length"] = maxlen + 4

        ## make sense of stats
        keepmj = depths[depths >= data.paramsdict["mindepth_majrule"]]
        keepstat = depths[depths >= data.paramsdict["mindepth_statistical"]]

        ## sample summary stat assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = depths.shape[0]
        sample.stats["clusters_hidepth"] = keepmj.shape[0]

        ## store depths histogram as a dict. Limit to first 25 bins
        bars, bins = np.histogram(depths, bins=range(1, 26))
        sample.depths = {int(i): int(v) for i, v in zip(bins, bars) if v}

        ## sample stat assignments
        ## Trap numpy warnings ("mean of empty slice") printed by samples
        ## with few reads.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
            sample.stats_dfs.s3["clusters_total"] = depths.shape[0]
            try:
                sample.stats_dfs.s3["clusters_hidepth"] = (
                    int(sample.stats["clusters_hidepth"]))
            except ValueError:
                ## Handle clusters_hidepth == NaN
                sample.stats_dfs.s3["clusters_hidepth"] = 0
            sample.stats_dfs.s3["avg_depth_total"] = depths.mean()
            sample.stats_dfs.s3["avg_depth_mj"] = keepmj.mean()
            sample.stats_dfs.s3["avg_depth_stat"] = keepstat.mean()
            sample.stats_dfs.s3["sd_depth_total"] = depths.std()
            sample.stats_dfs.s3["sd_depth_mj"] = keepmj.std()
            sample.stats_dfs.s3["sd_depth_stat"] = keepstat.std()

    else:
        print("    no clusters found for {}".format(sample.name))

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["assembly_method"] == "denovo":
        refmap_stats(data, sample)

    # if loglevel==DEBUG
    log_level = LOGGER.getEffectiveLevel()
    if not log_level == 10:
        # Clean up loose files only if not in DEBUG
        # edits/*derep, utemp, *utemp.sort, *htemp, *clust.gz
        derepfile = os.path.join(data.dirs.edits, sample.name + "_derep.fastq")
        mergefile = os.path.join(data.dirs.edits, sample.name + "_merged.fastq")
        uhandle = os.path.join(data.dirs.clusts, sample.name + ".utemp")
        usort = os.path.join(data.dirs.clusts, sample.name + ".utemp.sort")
        hhandle = os.path.join(data.dirs.clusts, sample.name + ".htemp")
        clusters = os.path.join(data.dirs.clusts, sample.name + ".clust.gz")

        for f in [derepfile, mergefile, uhandle, usort, hhandle, clusters]:
            if os.path.exists(f):
                os.remove(f)


def _aligned_indel_filter(clust, max_internal_indels):
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


def _get_derep_num(read):
    "return the number of replicates in a derep read"
    return int(read.split("=")[-1].split("\n")[0][:-1])


def _gbs_trim(align1):
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

