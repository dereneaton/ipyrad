#!/usr/bin/env python

""" 
Modifies and/or trims reads based on quality scores, presence of adapters, 
and user entered options to trim ends of reads. Uses the probabilistic trimming
methods implemented in the 'cutadapt' software.
"""

from __future__ import print_function

import os
import time
import numpy as np
import subprocess as sps
from .utils import IPyradWarningExit, IPyradError, fullcomp


class Step2(object):
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::2])
        self.print_headers()
        self.samples = self.get_subsamples()
        self.check_binaries()
        self.setup_dirs()
        self.check_adapters()


    def print_headers(self):
        if self.data._cli:
            self.data._print(
                "\n{}Step 2: Filtering and trimming reads"
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        "Apply state, ncluster, and force filters to select samples"

        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 2")

        # filter samples by state
        state1 = self.data.stats.index[self.data.stats.state == 1]
        statex = self.data.stats.index[self.data.stats.state > 1]

        # build list to run for samples being forced
        if self.force:
            subsamples = list(self.data.samples.values())
        else:
            # tell user which samples have already completed step 2
            if statex.any():
                print("skipping samples already finished step 2:\n{}"
                      .format(statex.tolist()))
            # run all samples in state 1
            subsamples = [self.data.samples[i] for i in state1]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_raw:
                checked_samples.append(sample)
            else:
                print("skipping {}; no reads found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 3")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_raw,
            reverse=True,
        )
        return checked_samples


    def setup_dirs(self):
        self.data.dirs.edits = os.path.join(
            self.data.params.project_dir,
            "{}_edits".format(self.data.name))
        if not os.path.exists(self.data.dirs.edits):
            os.makedirs(self.data.dirs.edits)        


    def check_binaries(self):
        cmd = ['which', 'cutadapt']
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if not comm:
            raise IPyradError("program 'cutadapt' not found.")


    def check_adapters(self):
        """
        Allow extra adapters if filters=3, and add poly repeats if not 
        in list of adapters. 
        """
        if int(self.data.params.filter_adapters) == 3:
            if not self.data.hackersonly.p3_adapters_extra:
                for poly in ["A" * 8, "T" * 8, "C" * 8, "G" * 8]:
                    self.data.hackersonly.p3_adapters_extra = (
                        self.data.hackersonly.p3_adapters_extra + [poly])

            if not self.data.hackersonly.p5_adapters_extra:
                for poly in ["A" * 8, "T" * 8, "C" * 8, "G" * 8]:
                    self.data.hackersonly.p5_adapters_extra = (
                        self.data.hackersonly.p5_adapters_extra + [poly])
        else:
            self.data.hackersonly.p5_adapters_extra = []
            self.data.hackersonly.p3_adapters_extra = []        


    def run(self):
        self.remote_concat_multiple_raws()
        self.remote_run_cutadapt()
        self.assembly_cleanup()
        self.data.save()


    def remote_concat_multiple_raws(self):
        "concatenate multiple raw files into a single file."

        # if no samples have multiple then just move on
        if not any([len(i.files.fastqs) > 1 for i in self.samples]):
            for sample in self.samples:
                sample.files.concat = sample.files.fastqs

        # otherwise concatenate them
        else:
            # run on single engine due to i/o limits
            start = time.time()
            printstr = ("concatenating inputs", "s2")
            finished = 0
            catjobs = {}
            for sample in self.samples:
                if len(sample.files.fastqs) > 1:
                    catjobs[sample.name] = self.lbview.apply(
                        concat_multiple_inputs, *(self.data, sample))
                else:
                    sample.files.concat = sample.files.fastqs

            # wait for all to finish
            while 1:
                finished = sum([i.ready() for i in catjobs.values()])
                self.data._progressbar(len(catjobs), finished, start, printstr)
                time.sleep(0.1)
                if finished == len(catjobs):
                    break

            # collect results, which are concat file handles.
            self.data._print("")
            for rasync in catjobs:
                if catjobs[rasync].successful():
                    self.data.samples[rasync].files.concat = (
                        catjobs[rasync].result())
                else:
                    error = catjobs[rasync].result()  # exception()
                    # ip.logger.error("error in step2 concat %s", error)
                    raise IPyradError(
                        "error in step2 concat: {}".format(error))


    def remote_run_cutadapt(self):
        # choose cutadapt function based on datatype
        start = time.time()
        printstr = ("processing reads    ", "s2")
        finished = 0
        rawedits = {}

        # send samples to cutadapt filtering
        for sample in self.samples:
            if "pair" in self.data.params.datatype:
                rasync = self.lbview.apply(
                    cutadaptit_pairs, *(self.data, sample))
            else:
                rasync = self.lbview.apply(
                    cutadaptit_single, *(self.data, sample))
            rawedits[sample.name] = rasync

        ## wait for all to finish
        while 1:
            finished = sum([i.ready() for i in rawedits.values()])
            self.data._progressbar(len(rawedits), finished, start, printstr)
            time.sleep(0.1)
            if finished == len(rawedits):
                print("")
                break

        ## collect results, report failures, store stats. async = sample.name
        for rasync in rawedits:
            if rawedits[rasync].successful():
                res = rawedits[rasync].result()

                ## if single cleanup is easy
                if "pair" not in self.data.params.datatype:
                    parse_single_results(
                        self.data, self.data.samples[rasync], res)
                else:
                    parse_pair_results(
                        self.data, self.data.samples[rasync], res)
            else:
                raise IPyradError(
                    "error in {}: {}"
                    .format(rasync, rawedits[rasync].get())
                )


    def assembly_cleanup(self):
        # build s2 results data frame
        self.data.stats_dfs.s2 = self.data._build_stat("s2")
        self.data.stats_files.s2 = os.path.join(
            self.data.dirs.edits, 
            's2_rawedit_stats.txt')

        # write stats for all samples
        with open(self.data.stats_files.s2, 'w') as outfile:
            (
                self.data.stats_dfs.s2.fillna(value=0)
                .astype(np.int)
                .to_string(outfile)
            )



def assembly_cleanup(data):
    "cleanup for assembly object"

    ## build s2 results data frame
    data.stats_dfs.s2 = data._build_stat("s2")
    data.stats_files.s2 = os.path.join(data.dirs.edits, 's2_rawedit_stats.txt')

    ## write stats for all samples
    with open(data.stats_files.s2, 'w') as outfile:
        data.stats_dfs.s2.fillna(value=0).astype(np.int).to_string(outfile)



def parse_single_results(data, sample, res1):
    """ parse results from cutadapt into sample data"""

    ## set default values 
    #sample.stats_dfs.s2["reads_raw"] = 0
    sample.stats_dfs.s2["trim_adapter_bp_read1"] = 0
    sample.stats_dfs.s2["trim_quality_bp_read1"] = 0
    sample.stats_dfs.s2["reads_filtered_by_Ns"] = 0
    sample.stats_dfs.s2["reads_filtered_by_minlen"] = 0
    sample.stats_dfs.s2["reads_passed_filter"] = 0

    ## parse new values from cutadapt results output
    lines = res1.decode().strip().split("\n")
    for line in lines:

        if "Total reads processed:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["reads_raw"] = value

        if "Reads with adapters:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["trim_adapter_bp_read1"] = value

        if "Quality-trimmed" in line:
            value = int(line.split()[1].replace(",", ""))
            sample.stats_dfs.s2["trim_quality_bp_read1"] = value

        if "Reads that were too short" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_minlen"] = value

        if "Reads with too many N" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_Ns"] = value
   
        if "Reads written (passing filters):" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["reads_passed_filter"] = value

    ## save to stats summary
    if sample.stats_dfs.s2.reads_passed_filter:
        sample.stats.state = 2
        sample.stats.reads_passed_filter = (
            sample.stats_dfs.s2.reads_passed_filter)
        sample.files.edits = [
            (os.path.join(
                data.dirs.edits, sample.name + ".trimmed_R1_.fastq.gz"), 0)]
        ## write the long form output to the log file.
        # ip.logger.info(res1)

    else:
        print("{}No reads passed filtering in Sample: {}"
              .format(data._spacer, sample.name))



def parse_pair_results(data, sample, res):
    """ parse results from cutadapt for paired data"""
    # ip.logger.info("in parse pair mod results\n%s", res)   
    ## set default values
    sample.stats_dfs.s2["trim_adapter_bp_read1"] = 0
    sample.stats_dfs.s2["trim_adapter_bp_read2"] = 0    
    sample.stats_dfs.s2["trim_quality_bp_read1"] = 0
    sample.stats_dfs.s2["trim_quality_bp_read2"] = 0    
    sample.stats_dfs.s2["reads_filtered_by_Ns"] = 0
    sample.stats_dfs.s2["reads_filtered_by_minlen"] = 0
    sample.stats_dfs.s2["reads_passed_filter"] = 0

    lines = res.decode().strip().split("\n")
    qprimed = 0
    for line in lines:
        ## set primer to catch next line
        if "Quality-trimmed" in line:
            qprimed = 1

        ## grab read1 and read2 lines when qprimed
        if "Read 1:" in line:
            if qprimed:
                value = int(line.split()[2].replace(",", ""))
                sample.stats_dfs.s2["trim_quality_bp_read1"] = value

        if "Read 2:" in line:
            if qprimed:
                value = int(line.split()[2].replace(",", ""))
                sample.stats_dfs.s2["trim_quality_bp_read2"] = value
                qprimed = 0

        if "Read 1 with adapter:" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["trim_adapter_bp_read1"] = value

        if "Read 2 with adapter:" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["trim_adapter_bp_read2"] = value

        if "Total read pairs processed:" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["reads_raw"] = value

        if "Pairs that were too short" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_minlen"] = value

        if "Pairs with too many N" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_Ns"] = value

        if "Pairs written (passing filters):" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["reads_passed_filter"] = value

    ## save to stats summary
    if sample.stats_dfs.s2.reads_passed_filter:
        sample.stats.state = 2
        sample.stats.reads_passed_filter = (
            sample.stats_dfs.s2.reads_passed_filter)
        sample.files.edits = [(
            os.path.join(
                data.dirs.edits, sample.name + ".trimmed_R1_.fastq.gz"), 
            os.path.join(
                data.dirs.edits, sample.name + ".trimmed_R2_.fastq.gz")
            )]
    else:
        print("No reads passed filtering in Sample: {}".format(sample.name))


# CALLED BY STEP
def cutadaptit_single(data, sample):
    """ 
    Applies quality and adapter filters to reads using cutadapt. If the ipyrad
    filter param is set to 0 then it only filters to hard trim edges and uses
    mintrimlen. If filter=1, we add quality filters. If filter=2 we add
    adapter filters. 
    """
    sname = sample.name
    # if (GBS, ddRAD) we look for the second cut site + adapter. For SE
    # data we don't bother trying to remove the second barcode since it's not
    # as critical as with PE data.
    if data.params.datatype == "rad":
        adapter = data.hackersonly.p3_adapter
    
    else:
        # if GBS then the barcode can also be on the other side. 
        if data.param.datatype == "gbs":
           

            if data.barcodes:

                # make full adapter (-revcompcut-revcompbarcode-adapter)                
                adapter = "".join([
                    fullcomp(data.params.restriction_overhang[1])[::-1], 
                    fullcomp(data.barcodes[sample.name])[::-1], 
                    data.hackersonly.p3_adapter, 
                ])
            
                # and add adapter without revcompbarcode (incomplete)
                incomplete_adapter = "".join([
                    fullcomp(data.params.restriction_overhang[1])[::-1], 
                    data.hackersonly.p3_adapter
                ])                

                # append incomplete adapter to extras (-recompcut-adapter)
                data.hackersonly.p3_adapters_extra = [
                    data.hackersonly.p3_adapters_extra, 
                    incomplete_adapter,
                ]

            else:
                # else no search for barcodes on 3'
                adapter = "".join([
                    fullcomp(data.params.restriction_overhang[1])[::-1], 
                    data.hackersonly.p3_adapter, 
                ])

        # not GBS, simple
        else:
            adapter = "".join([
                fullcomp(data.params.restriction_overhang[1])[::-1], 
                data.hackersonly.p3_adapter, 
            ])

    # get length trim parameter from new or older version of ipyrad params
    trim5r1 = trim3r1 = []
    trimlen = data.params.trim_reads
        
    # trim 5' end
    if trimlen[0]:
        trim5r1 = ["-u", str(trimlen[0])]
    if trimlen[1] < 0:
        trim3r1 = ["-u", str(trimlen[1])]
    if trimlen[1] > 0:
        trim3r1 = ["--length", str(trimlen[1])]
    # else:
        # trimlen = data.paramsdict.get("edit_cutsites")
        # trim5r1 = ["--cut", str(trimlen[0])]

    ## testing new 'trim_reads' setting
    cmdf1 = ["cutadapt"]
    if trim5r1:
        cmdf1 += trim5r1
    if trim3r1:
        cmdf1 += trim3r1
    cmdf1 += [
        "--minimum-length", str(data.params.filter_min_trim_len),
        "--max-n", str(data.params.max_low_qual_bases),
        "--trim-n", 
        "--output", os.path.join(
            data.dirs.edits, sname + ".trimmed_R1_.fastq.gz"),
        sample.files.concat[0][0],
        ]

    if int(data.params.filter_adapters):
        ## NEW: only quality trim the 3' end for SE data.
        cmdf1.insert(1, "20")
        cmdf1.insert(1, "-q")
        cmdf1.insert(1, str(data.params.phred_Qscore_offset))
        cmdf1.insert(1, "--quality-base")

    ## if filter_adapters==3 then p3_adapters_extra will already have extra
    ## poly adapters added to its list. 
    if int(data.params.filter_adapters) > 1:
        ## first enter extra cuts (order of input is reversed)
        for extracut in list(set(data.hackersonly.p3_adapters_extra))[::-1]:
            cmdf1.insert(1, extracut)
            cmdf1.insert(1, "-a")
        ## then put the main cut so it appears first in command
        cmdf1.insert(1, adapter)
        cmdf1.insert(1, "-a")

    ## do modifications to read1 and write to tmp file
    proc1 = sps.Popen(cmdf1, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    try:
        res1 = proc1.communicate()[0]
    except KeyboardInterrupt:
        proc1.kill()
        raise KeyboardInterrupt

    ## raise errors if found
    if proc1.returncode:
        raise IPyradWarningExit(" error in {}\n {}".format(" ".join(cmdf1), res1))

    ## return result string to be parsed outside of engine
    return res1


# CALLED BY STEP
## BEING MODIFIED FOR MULTIPLE BARCODES (i.e., merged samples. NOT PERFECT YET)
def cutadaptit_pairs(data, sample):
    """
    Applies trim & filters to pairs, including adapter detection. If we have
    barcode information then we use it to trim reversecut+bcode+adapter from 
    reverse read, if not then we have to apply a more general cut to make sure 
    we remove the barcode, this uses wildcards and so will have more false 
    positives that trim a little extra from the ends of reads. Should we add
    a warning about this when filter_adapters=2 and no barcodes?
    """
    sname = sample.name

    ## applied to read pairs
    finput_r1 = sample.files.concat[0][0]
    finput_r2 = sample.files.concat[0][1]

    ## Get adapter sequences. This is very important. For the forward adapter
    ## we don't care all that much about getting the sequence just before the 
    ## Illumina adapter, b/c it will either be random (in RAD), or the reverse
    ## cut site of cut1 or cut2 (gbs or ddrad). Either way, we can still trim it 
    ## off later in step7 with trim overhang if we want. And it should be invar-
    ## iable unless the cut site has an ambiguous char. The reverse adapter is 
    ## super important, however b/c it can contain the inline barcode and 
    ## revcomp cut site. We def want to trim out the barcode, and ideally the 
    ## cut site too to be safe. Problem is we don't always know the barcode if 
    ## users demultiplexed their data elsewhere. So, if barcode is missing we 
    ## do a very fuzzy match before the adapter and trim it out. 

    ## this just got more complicated now that we allow merging technical
    ## replicates in step 1 since a single sample might have multiple barcodes
    ## associated with it and so we need to search for multiple adapter+barcode
    ## combinations.
    ## We will assume that if they are 'linking_barcodes()' here then there are
    ## no technical replicates in the barcodes file. If there ARE technical
    ## replicates, then they should run step1 so they are merged, in which case
    ## the sample specific barcodes will be saved to each Sample under its
    ## .barcode attribute as a list. 

    # try linking barcodes again in case user just added a barcodes path
    if not data.barcodes:
        try:
            data._link_barcodes()
        except IPyradError as inst:
            pass

    # only effects this engine copy of barcodes (not lost from real data)
    if data.params.datatype == "pair3rad":
        data.barcodes = {}       

    # barcodes are present meaning they were parsed to the samples in step 1.
    if data.barcodes:
        try:
            adapter1 = "".join([
                fullcomp(data.params.restriction_overhang[1])[::-1], 
                data.hackersonly.p3_adapter, 
            ])
            
            # which barcode
            if isinstance(sample.barcode, list):
                bcode = fullcomp(sample.barcode[0])[::-1]
            elif isinstance(data.barcodes[sample.name], list):
                bcode = fullcomp(data.barcodes[sample.name][0][::-1])
            else:
                bcode = fullcomp(data.barcodes[sample.name])[::-1]

            # add full adapter (-revcompcut-revcompbcode-adapter)
            adapter2 = "".join([
                fullcomp(data.params.restriction_overhang[0])[::-1], 
                bcode,
                data.hackersonly.p5_adapter, 
            ])

        except KeyError as inst:
            msg = """
    Sample name does not exist in the barcode file. The name in the barcode file
    for each sample must exactly equal the raw file name for the sample minus
    `_R1`. So for example a sample called WatDo_PipPrep_R1_100.fq.gz must
    be referenced in the barcode file as WatDo_PipPrep_100. The name in your
    barcode file for this sample must match: {}
    """.format(sample.name)
            raise IPyradWarningExit(msg)

    else:
        if data.params.datatype != "pair3rad":
            print(NO_BARS_GBS_WARNING)
        adapter1 = data.hackersonly.p3_adapter
        adapter2 = fullcomp(data.hackersonly.p5_adapter)

    # parse trim_reads
    trim5r1 = trim5r2 = trim3r1 = trim3r2 = []
    trimlen = data.params.trim_reads
        
    # trim 5' end
    if trimlen[0]:
        trim5r1 = ["-u", str(trimlen[0])]
    if trimlen[1] < 0:
        trim3r1 = ["-u", str(trimlen[1])]
    if trimlen[1] > 0:
        trim3r1 = ["--length", str(trimlen[1])]

    # legacy support for trimlen = 0,0 default
    if len(trimlen) > 2:
        if trimlen[2]:
            trim5r2 = ["-U", str(trimlen[2])]

    if len(trimlen) > 3:
        if trimlen[3]:
            if trimlen[3] < 0:
                trim3r2 = ["-U", str(trimlen[3])]
            if trimlen[3] > 0:            
                trim3r2 = ["--length", str(trimlen[3])]

    # testing new 'trim_reads' setting
    cmdf1 = ["cutadapt"]
    if trim5r1:
        cmdf1 += trim5r1
    if trim3r1:
        cmdf1 += trim3r1
    if trim5r2:
        cmdf1 += trim5r2
    if trim3r2:
        cmdf1 += trim3r2

    cmdf1 += [
        "--trim-n",
        "--max-n", str(data.params.max_low_qual_bases),
        "--minimum-length", str(data.params.filter_min_trim_len),
        "-o", os.path.join(
            data.dirs.edits, sname + ".trimmed_R1_.fastq.gz"),
        "-p", os.path.join(
            data.dirs.edits, sname + ".trimmed_R2_.fastq.gz"),
        finput_r1,
        finput_r2,
        ]

    ## additional args
    if int(data.params.filter_adapters) < 2:
        # add a dummy adapter to let cutadapt know we are not using legacy-mode
        cmdf1.insert(1, "XXX")
        cmdf1.insert(1, "-A")

    if int(data.params.filter_adapters):
        cmdf1.insert(1, "20,20")
        cmdf1.insert(1, "-q")
        cmdf1.insert(1, str(data.params.phred_Qscore_offset))
        cmdf1.insert(1, "--quality-base")

    if int(data.params.filter_adapters) > 1:
        ## if technical replicates then add other copies
        if isinstance(sample.barcode, list):
            for extrabar in sample.barcode[1:]:
                
                data.hackersonly.p5_adapters_extra.append(
                    "".join([
                        fullcomp(data.params.restriction_overhang[0])[::-1], 
                        fullcomp(extrabar)[::-1], 
                        data.hackersonly.p5_adapter
                    ])
                )

                data.hackersonly.p5_adapters_extra.append(
                    "".join([
                        fullcomp(data.params.restriction_overhang[1])[::-1], 
                        data.hackersonly.p3_adapter,
                    ])
                )

        ## first enter extra cuts
        zcut1 = list(set(data.hackersonly.p3_adapters_extra))[::-1]
        zcut2 = list(set(data.hackersonly.p5_adapters_extra))[::-1]
        for ecut1, ecut2 in zip(zcut1, zcut2):
            cmdf1.insert(1, ecut1)
            cmdf1.insert(1, "-a")
            cmdf1.insert(1, ecut2)
            cmdf1.insert(1, "-A")
        ## then put the main cut first
        cmdf1.insert(1, adapter1)
        cmdf1.insert(1, '-a')        
        cmdf1.insert(1, adapter2)
        cmdf1.insert(1, '-A')         

    ## do modifications to read1 and write to tmp file
    proc1 = sps.Popen(cmdf1, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    res1 = proc1.communicate()[0]
    if proc1.returncode:
        raise IPyradWarningExit("error in cutadapt: {}".format(res1.decode()))
    return res1



# CALLED BY STEP
def concat_multiple_inputs(data, sample):
    """ 
    If multiple fastq files were appended into the list of fastqs for samples
    then we merge them here before proceeding. 
    """
    ## if more than one tuple in fastq list
    if len(sample.files.fastqs) > 1:
        ## create a cat command to append them all (doesn't matter if they 
        ## are gzipped, cat still works). Grab index 0 of tuples for R1s.
        cmd1 = ["cat"] + [i[0] for i in sample.files.fastqs]

        isgzip = ".gz"
        if not sample.files.fastqs[0][0].endswith(".gz"):
            isgzip = ""

        ## write to new concat handle
        conc1 = os.path.join(
            data.dirs.edits, sample.name + "_R1_concat.fq{}".format(isgzip))
        with open(conc1, 'w') as cout1:
            proc1 = sps.Popen(
                cmd1, stderr=sps.STDOUT, stdout=cout1, close_fds=True)
            res1 = proc1.communicate()[0]
        if proc1.returncode:
            raise IPyradWarningExit("error in: {}, {}".format(cmd1, res1))

        ## Only set conc2 if R2 actually exists
        conc2 = 0
        if "pair" in data.params.datatype:
            cmd2 = ["cat"] + [i[1] for i in sample.files.fastqs]
            conc2 = os.path.join(
                data.dirs.edits, sample.name + "_R2_concat.fq{}".format(isgzip))
            with open(conc2, 'w') as cout2:
                proc2 = sps.Popen(
                    cmd2, stderr=sps.STDOUT, stdout=cout2, close_fds=True)
                proc2.communicate()[0]
            if proc2.returncode:
                raise IPyradWarningExit(
                    "Error concatenating fastq files. Make sure all "\
                  + "these files exist: {}\nError message: {}"
                    .format(cmd2, proc2.returncode))

        ## store new file handles
        sample.files.concat = [(conc1, conc2)]
    return sample.files.concat


## GLOBALS
NO_BARS_GBS_WARNING = """\
    This is a just a warning: 
    You set 'filter_adapters' to 2 (stringent), however, b/c your data
    is paired-end gbs or ddrad the second read is likely to contain the 
    barcode in addition to the adapter sequence. Actually, it would be:
                  [sequence][cut overhang][barcode][adapter]
    If you add a barcode table to your params it will improve the accuracy of 
    adapter trimming, especially if your barcodes are of varying lengths, and 
    ensure proper removal of barcodes. If you do not have this information 
    that's OK, but we will apply a slightly more rigorous trimming of 3' edges 
    on R2 that results in more false positives (more bp trimmed off of R2). 
    """
