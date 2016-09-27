#!/usr/bin/env ipython2

""" 
Modifies and/or trims reads based on quality scores, presence of adapters, 
and user entered options to trim ends of reads. Uses the probabilistic trimming
methods implemented in the 'cutadapt' software.
"""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0301

import os
import time
import datetime
import subprocess as sps
import numpy as np
import glob
import shutil
from .util import *

import logging
LOGGER = logging.getLogger(__name__)

## shortcut name for os.path.join
OPJ = os.path.join




def assembly_cleanup(data):
    """ cleanup for assembly object """

    ## build s2 results data frame
    data.stats_dfs.s2 = data.build_stat("s2")
    data.stats_files.s2 = os.path.join(data.dirs.edits, 's2_rawedit_stats.txt')

    ## write stats for all samples
    with open(data.stats_files.s2, 'w') as outfile:
        data.stats_dfs.s2.astype(np.int).to_string(outfile)



def parse_single_cut_results(data, sample, res1):
    """ parse results from cutadapt into sample data"""

    ## set default values 
    #sample.stats_dfs.s2["reads_raw"] = 0
    sample.stats_dfs.s2["trimmed_bp_adapter"] = 0
    sample.stats_dfs.s2["trimmed_bp_quality"] = 0
    sample.stats_dfs.s2["reads_filtered_by_trimlen"] = 0
    sample.stats_dfs.s2["reads_filtered_by_Ns"] = 0
    sample.stats_dfs.s2["reads_passed_filter"] = 0

    ## parse new values from cutadapt results output
    lines = res1.strip().split("\n")
    for line in lines:

        if "Total reads processed:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["reads_raw"] = value

        if "Reads with adapters:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["trimmed_bp_adapter"] = value

        if "Quality-trimmed" in line:
            value = int(line.split()[1].replace(",", ""))
            sample.stats_dfs.s2["trimmed_bp_quality"] = value

        if "Reads that were too short" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_trimlen"] = value

        if "Reads with too many N" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_Ns"] = value
   
        if "Reads written (passing filters):" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["reads_passed_filter"] = value

    ## save to stats summary
    if sample.stats_dfs.s2.reads_passed_filter:
        sample.stats.state = 2
        sample.stats.reads_passed_filter = sample.stats_dfs.s2.reads_passed_filter
        sample.files.edits = [
            (OPJ(data.dirs.edits, sample.name+".trimmed_R1_.fastq.gz"), 0)]
        ## write the long form output to the log file.
        LOGGER.info(res1)

    else:
        print("No reads passed filtering in Sample: {}".format(sample.name))



def parse_pair_adapt_results(data, sample, res):
    """ get stats from cutadapt on pairs """
    LOGGER.info(res)
    lines = res.strip().split("\n")
    for line in lines:

        if "Reads with adapters:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["trimmed_bp_adapter"] = value

        if "Total reads processed:" in line:
            value = int(line.split()[3].replace(",", ""))
            sample.stats_dfs.s2["reads_raw"] = value

        if "Reads that were too short" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_trimlen"] = value

        if "Reads with too many N" in line:
            value = int(line.split()[5].replace(",", ""))
            sample.stats_dfs.s2["reads_filtered_by_Ns"] = value

        if "Reads written (passing filters):" in line:
            value = int(line.split()[4].replace(",", ""))
            sample.stats_dfs.s2["reads_passed_filter"] = value

    ## save to stats summary
    if sample.stats_dfs.s2.reads_passed_filter:
        sample.stats.state = 2
        sample.stats.reads_passed_filter = sample.stats_dfs.s2.reads_passed_filter
        sample.files.edits = [(
             OPJ(data.dirs.edits, sample.name+".trimmed_R1_.fastq.gz"), 
             OPJ(data.dirs.edits, sample.name+".trimmed_R2_.fastq.gz")
             )]
        ## write the long form output to the log file.
        LOGGER.info(res)

    else:
        print("No reads passed filtering in Sample: {}".format(sample.name))



def parse_pair_mod_results(sample, res):
    """ parse results from cutadapt for paired data"""

    lines = res.strip().split("\n")
    for line in lines:
        if "Quality-trimmed" in line:
            value = int(line.split()[1].replace(",", ""))
            if np.isnan(sample.stats_dfs.s2["trimmed_bp_quality"]):
                sample.stats_dfs.s2["trimmed_bp_quality"] = value
            else:
                sample.stats_dfs.s2["trimmed_bp_quality"] += value
   




def cutadaptit_single(data, sample):
    """ 
    Applies quality and adapter filters to reads using cutadapt. If the ipyrad
    filter param is set to 0 then it only filters to hard trim edges and uses
    mintrimlen. If filter=1, we add quality filters. If filter=2 we add
    adapter filters. 
    """
    sname = sample.name

    ## if (GBS, ddRAD) then use second cut site with adapter.
    if data.paramsdict["datatype"] != "rad":
        adapter = fullcomp(data.paramsdict["restriction_overhang"][1])[::-1]+\
                  data._hackersonly["p3_adapter"]
    else:
        adapter = data._hackersonly["p3_adapter"]

    ## get command line argument
    if not data.paramsdict["filter_adapters"]:
        cmdf1 = ["cutadapt", 
                 "--cut", str(data.paramsdict["edit_cutsites"][0]),
                 "--minimum-length", str(data.paramsdict["filter_min_trim_len"]),
                 "--trim-n", 
                 "--max-n", str(data.paramsdict["max_low_qual_bases"]),
                 "--output", OPJ(data.dirs.edits, sname+".trimmed_R1_.fastq.gz"),
                 sample.files.concat[0][0]]

    elif data.paramsdict["filter_adapters"] == 1:
        cmdf1 = ["cutadapt", 
                 "--cut", str(data.paramsdict["edit_cutsites"][0]),
                 "--minimum-length", str(data.paramsdict["filter_min_trim_len"]),
                 "--trim-n", 
                 "--max-n", str(data.paramsdict["max_low_qual_bases"]),
                 "--quality-cutoff", "20,20", 
                 "--quality-base", str(data.paramsdict["phred_Qscore_offset"]),
                 "--output", OPJ(data.dirs.edits, sname+".trimmed_R1_.fastq.gz"),
                 sample.files.concat[0][0]]

    else:
        cmdf1 = ["cutadapt", 
                 "-a", adapter,
                 "--cut", str(data.paramsdict["edit_cutsites"][0]),
                 "--minimum-length", str(data.paramsdict["filter_min_trim_len"]),
                 "--trim-n", 
                 "--max-n", str(data.paramsdict["max_low_qual_bases"]),
                 "--quality-cutoff", "20,20", 
                 "--quality-base", str(data.paramsdict["phred_Qscore_offset"]),
                 "--output", OPJ(data.dirs.edits, sname+".trimmed_R1_.fastq.gz"),
                 sample.files.concat[0][0]]

    ## do modifications to read1 and write to tmp file
    LOGGER.info(cmdf1)
    proc1 = sps.Popen(cmdf1, stderr=sps.STDOUT, stdout=sps.PIPE)
    try:
        res1 = proc1.communicate()[0]
    except KeyboardInterrupt:
        res1.kill()

    ## raise errors if found
    if proc1.returncode:
        raise IPyradWarningExit(" error in %s, %s", cmdf1, res1)

    ## return result string to be parsed outside of engine
    return res1



def cutadaptit_pair_mod(data, sample, readpair):
    """ 
    Applies quality and edge trimming to reads using cutadapt but does not yet
    apply filters because we do can't filter R1 without R2. 
    """

    sname = sample.name

    ## FIRST round of filtering:
    ## apply filter for trim_edges and quality trimmer, but no filters yet 
    ## b/c we wouldn't want to filter R1 without R2
    if readpair == 1:
        trimedge = str(data.paramsdict["edit_cutsites"][0])
        foutput = OPJ(data.dirs.edits, sname+".tmp_R1_.fastq.gz")
        finput = sample.files.concat[0][0]
    else:
        trimedge = str(data.paramsdict["edit_cutsites"][1])
        foutput = OPJ(data.dirs.edits, sname+".tmp_R2_.fastq.gz")
        finput = sample.files.concat[0][1]

    if not data.paramsdict["filter_adapters"]:
        cmdf1 = ["cutadapt", "--cut", trimedge,
                             "--trim-n",
                             "--output", foutput,
                             finput]

    else:
        cmdf1 = ["cutadapt", "--cut", trimedge,
                             "--trim-n",
                             "--quality-cutoff", "20,20", 
                             "--quality-base", str(data.paramsdict["phred_Qscore_offset"]),         
                             "--output", foutput,
                             finput]

    ## do modifications to read1 and write to tmp file
    proc1 = sps.Popen(cmdf1, stderr=sps.STDOUT, stdout=sps.PIPE)
    try:
        res1 = proc1.communicate()[0]
    except KeyboardInterrupt:
        res1.kill()
    ## raise errors if found
    if proc1.returncode:
        raise IPyradWarningExit(" error in %s, %s", cmdf1, res1)

    return res1



def cutadaptit_pair_adapt(data, sample):
    """
    Applies filters to readpairs, including adapter detection. If we have
    barcode information then we can use it to trim reversecut+bcode+adapter from 
    reverse read, if not then we have to apply a more general cut to make sure 
    we remove the barcode. 
    """

    sname = sample.name

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

    if data.barcodes:
        adapter1 = fullcomp(data.paramsdict["restriction_overhang"][1])[::-1]+\
                   data._hackersonly["p3_adapter"]
        adapter2 = fullcomp(data.paramsdict["restriction_overhang"][0])[::-1]+\
                   data.barcodes[sample.name]+\
                   data._hackersonly["p5_adapter"]
    else:
        adapter1 = fullcomp(data.paramsdict["restriction_overhang"][1])[::-1]+\
                   data._hackersonly["p3_adapter"]
        adapter2 = "N"*len(data.paramsdict["restriction_overhang"][0])+\
                   "N"*6+\
                   data._hackersonly["p5_adapter"]


    ## SECOND round of filtering:
    ## apply filter for adapter detection and filter pairs
    ## filter reads1
    if data.paramsdict["filter_adapters"] != 2:    
        cmdf3 = ["cutadapt", 
                  "--max-n", str(data.paramsdict["max_low_qual_bases"]), 
                  "--minimum-length", str(data.paramsdict["filter_min_trim_len"]),
                  "-o", OPJ(data.dirs.edits, sname+".trimmed_R1_.fastq.gz"), 
                  "-p", OPJ(data.dirs.edits, sname+".trimmed_R2_.fastq.gz"),
                  OPJ(data.dirs.edits, sname+".tmp_R1_.fastq.gz"), 
                  OPJ(data.dirs.edits, sname+".tmp_R2_.fastq.gz")]
    else:
        cmdf3 = ["cutadapt", 
                  "-a", adapter1, 
                  "-A", adapter2,
                  "--max-n", str(data.paramsdict["max_low_qual_bases"]), 
                  "--minimum-length", str(data.paramsdict["filter_min_trim_len"]),
                  "-o", OPJ(data.dirs.edits, sname+".trimmed_R1_.fastq.gz"), 
                  "-p", OPJ(data.dirs.edits, sname+".trimmed_R2_.fastq.gz"),
                  OPJ(data.dirs.edits, sname+".tmp_R1_.fastq.gz"), 
                  OPJ(data.dirs.edits, sname+".tmp_R2_.fastq.gz")]

    ## do modifications to read1 and write to tmp file
    proc3 = sps.Popen(cmdf3, stderr=sps.STDOUT, stdout=sps.PIPE)
    try:
        res3 = proc3.communicate()[0]
    except KeyboardInterrupt:
        res3.kill()
    ## raise errors if found
    if proc3.returncode:
        raise IPyradWarningExit(" error in %s, %s", cmdf3, res3)

    ## return results string to be parsed outside of engine
    return res3



def run2(data, samples, force, ipyclient):
    """ 
    Filter for samples that are already finished with this step, allow others
    to run, pass them to parallel client function to filter with cutadapt. 
    """
    ## check for cutadapt (at least until we get it into our own conda channel)
    try:
        import cutadapt
    except ImportError:
        raise IPyradWarningExit("""\
    Required dependency 'cutadapt' is missing. Run the following command to 
    install it: `conda install -c bioconda cutadapt`
    """)

    ## create output directories 
    data.dirs.edits = os.path.join(os.path.realpath(
                                   data.paramsdict["project_dir"]), 
                                   data.name+"_edits")
    if not os.path.exists(data.dirs.edits):
        os.makedirs(data.dirs.edits)

    ## get samples
    subsamples = choose_samples(samples, force)

    ## parallel client
    lbview = ipyclient.load_balanced_view()

    ## concatenate reads if they come from merged assemblies.
    if any([len(i.files.fastqs) > 1 for i in subsamples]):
        ## run on single engine for now
        start = time.time()
        finished = 0
        catjobs = {}
        for sample in subsamples:
            catjobs[sample.name] = ipyclient[0].apply(concat_muliple_inputs, *(data, sample))

        ## wait for all to finish
        while 1:
            finished = sum([i.ready() for i in catjobs.values()])
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(len(subsamples), finished, 
                       " concatenating inputs  | {}".format(elapsed))
            time.sleep(0.1)
            if finished == len(subsamples):
                break

        ## collect results, which are concat file handles.
        for async in catjobs:
            if catjobs[async].successful():
                data.samples[async].files.concat = catjobs[async].result
            else:
                raise IPyradWarningExit(catjobs[async].exception())
    else:
        for sample in subsamples:
            ## just copy fastqs handles to concat attribute
            sample.files.concat = sample.files.fastqs

    ## choose cutadapt function based on datatype
    start = time.time()
    finished = 0
    rawedits = {}

    ## send samples to cutadapt filtering
    if "pair" in data.paramsdict["datatype"]:
        for sample in subsamples:
            key1 = sample.name+"_read1"
            key2 = sample.name+"_read2"
            key3 = sample.name+"_pair"
            rawedits[key1] = lbview.apply(cutadaptit_pair_mod, *(data, sample, 1))
            rawedits[key2] = lbview.apply(cutadaptit_pair_mod, *(data, sample, 2))
            deps = [rawedits[key1], rawedits[key2]]
            with lbview.temp_flags(after=deps):
                rawedits[key3] = lbview.apply(cutadaptit_pair_adapt, *(data, sample))
    else:
        for sample in subsamples:
            rawedits[sample.name] = lbview.apply(cutadaptit_single, *(data, sample))

    ## wait for all to finish
    while 1:
        finished = sum([i.ready() for i in rawedits.values()])
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(rawedits), finished, 
                    " processing reads      | {}".format(elapsed))
        time.sleep(0.1)
        if finished == len(rawedits):
            break

    ## collect results, report failures, and store stats. async = sample.name
    for async in rawedits:
        if rawedits[async].successful():
            res = rawedits[async].result()

            ## if single cleanup is easy
            if "pair" not in data.paramsdict["datatype"]:
                parse_single_cut_results(data, data.samples[async], res)
            else:
                ## is it pairmod or pairadapt?
                name, key = async.rsplit("_", 1)
                LOGGER.info("name, key %s %s", name, key)
                if key == "read1":
                    parse_pair_mod_results(data.samples[name], res)
                elif key == "read2":
                    parse_pair_mod_results(data.samples[name], res)
                else:
                    parse_pair_adapt_results(data, data.samples[name], res)
        else:
            raise IPyradWarningExit(rawedits[async].exception())

    ## store sample results in data stats
    assembly_cleanup(data)



def choose_samples(samples, force):
    """ filter out samples that are already done with this step, unless force"""

    ## hold samples that pass
    subsamples = []

    ## filter the samples again
    if not force:
        for sample in samples:
            if sample.stats.state >= 2:
                print("""\
    Skipping Sample {}; Already filtered. Use force argument to overwrite.\
    """.format(sample.name))
            elif not sample.stats.reads_raw:
                print("""\
    Skipping Sample {}; No reads found in file {}\
    """.format(sample.name, sample.files.fastqs))
            else:
                subsamples.append(sample)

    else:
        for sample in samples:
            if not sample.stats.reads_raw:
                print("""\
    Skipping Sample {}; No reads found in file {}\
    """.format(sample.name, sample.files.fastqs))
            else:
                subsamples.append(sample)
    return subsamples



def concat_muliple_inputs(data, sample):
    """ 
    if multiple fastq files were appended into the list of fastqs for samples
    then we merge them here before proceeding. 
    """

    ## if more than one tuple in fastq list
    if len(sample.files.fastqs) > 1:
        ## create a cat command to append them all (doesn't matter if they 
        ## are gzipped, cat still works). Grab index 0 of tuples for R1s.
        cmd1 = ["cat"] + [i[0] for i in sample.files.fastqs]

        ## write to new concat handle
        conc1 = os.path.join(data.dirs.edits, sample.name+"_R1_concat.fq")
        with open(conc1, 'w') as cout1:
            proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=cout1)
            res1 = proc1.communicate()
        if proc1.returncode:
            raise IPyradWarningExit("error in: %s, %s", cmd1, res1)

        ## Only set conc2 if R2 actually exists
        conc2 = 0
        if os.path.exists(sample.files.fastqs[0][1]):
            cmd2 = ["cat"] + [i[1] for i in sample.files.fastqs]
            conc2 = os.path.join(data.dirs.edits, sample.name+"_R2_concat.fq")
            with open(conc2, 'w') as cout2:
                proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=cout2)
                res2 = proc2.communicate()
            if proc2.returncode:
                raise IPyradWarningExit("error in: %s, %s", cmd2, res2)

        ## store new file handles
        sample.files.concat = [(conc1, conc2)]
    return sample.files.concat



# def run(data, samples, nreplace, force, preview, ipyclient):
#     """ run the major functions for editing raw reads """

#     ## hold samples that pass
#     subsamples = []
#     asyncs = []

#     ## filter the samples again
#     if not force:
#         for sample in samples:
#             if sample.stats.state >= 2:
#                 print("""\
#     Skipping Sample {}; Already filtered. Use force argument to overwrite.\
#     """.format(sample.name))
#             elif not sample.stats.reads_raw:
#                 print("""\
#     Skipping Sample {}; No reads found in file {}\
#     """.format(sample.name, sample.files.fastqs))
#             else:
#                 subsamples.append(sample)

#     else:
#         for sample in samples:
#             if not sample.stats.reads_raw:
#                 print("""\
#     Skipping Sample {}; No reads found in file {}\
#     """.format(sample.name, sample.files.fastqs))
#             else:
#                 subsamples.append(sample)

#     ## get optim which is used to slice a sample into 10 equal sized chunks
#     ## if preview-mode then 10 chunks are made out of 'preview' reads.
#     optims = {}
#     for sample in subsamples:
#         if preview:
#             tots = data._hackersonly["preview_step2"]
#             assert not tots % 4, \
#             "_hackersonly preview_step2 value (nlines) must be divisible by 4."
#             optims[sample.name] = (tots // 10) + (tots % 10)
#         else:
#             tots = int(sample.stats.reads_raw)
#             optims[sample.name] = (tots // 10) + (tots % 10)

#     if preview:
#         if data._headers:
#             print("""\
#     Running preview mode: subselecting maximum of {} reads per sample\
#     """.format(data._hackersonly["preview_step2"]))

#     ## link output directories 
#     data.dirs.edits = os.path.join(os.path.realpath(
#                                    data.paramsdict["project_dir"]), 
#                                    data.name+"_edits")

#     ## create output directory if doesn't exist
#     if not os.path.exists(data.dirs.edits):
#         os.makedirs(data.dirs.edits)

#     ## client
#     lbview = ipyclient.load_balanced_view()

#     ## start progress
#     start = time.time()
#     elapsed = datetime.timedelta(seconds=int(0))
#     progressbar(len(subsamples), 0, 
#                " processing reads      | {}".format(elapsed))            

#     ## save sliced asyncs
#     sliced = {i.name:[] for i in subsamples}    

#     ## send jobs to queue to get slices and process them
#     for sample in subsamples:
#         ## get optim slice size for this sample
#         optim = optims[sample.name]
#         ## if multiple fastq files for this sample, create a tmp concat file
#         ## TODO: This could be parallelized, for merged real data, but only if 
#         ## disk-writing is parallelized...
#         if len(sample.files.fastqs) > 1:
#             ## just copy as a temporary placeholder for fastqs
#             conc1 = os.path.join(data.dirs.edits, sample.name+"_R1_concat.fq")
#             with open(conc1, 'w') as cout1:
#                 for tups in sample.files.fastqs:
#                     cout1.write(gzip.open(tups[0]).read())

#             ## Only set conc2 if R2 actually exists
#             conc2 = 0
#             try:
#                 if os.path.exists(sample.files.fastqs[0][1]):
#                     conc2 = os.path.join(data.dirs.edits, sample.name+"_R2_concat.fq")
#                     with open(conc2, 'w') as cout2:
#                         for tups in sample.files.fastqs:
#                             cout2.write(gzip.open(tups[1]).read())
#             except IndexError as _:
#                 ## if SE then no R2 files
#                 pass

#             sample.files.concat = [(conc1, conc2)]
#         else:
#             ## just copy as a temporary placeholder for fastqs
#             sample.files.concat = sample.files.fastqs

#         for job in range(10):
#             args = [data, sample, job, nreplace, optim]
#             async = lbview.apply(rawedit, args)
#             sliced[sample.name].append(async)

#     ## print progress
#     try:
#         while 1:
#             asyncs = list(itertools.chain(*sliced.values()))
#             tots = len(asyncs)
#             done = sum([i.ready() for i in asyncs])
#             elapsed = datetime.timedelta(seconds=int(time.time()-start))
#             if tots != done:
#                 progressbar(tots, done,
#                            " processing reads      | {}".format(elapsed))            
#                 time.sleep(1)
#             else:
#                 progressbar(tots, done,
#                            " processing reads      | {}".format(elapsed))
#                 #if data._headers:
#                 print("")
#                 break

#     except (KeyboardInterrupt, SystemExit):
#         print('\n  Interrupted! Cleaning up... ')
#         raise 

#     ## enforced cleanup
#     finally:
#         ## if all jobs were successful in a sample then cleanup
#         if asyncs:
#             for sample in subsamples:
#                 ## if finished
#                 if all([i.ready() for i in sliced[sample.name]]):
#                     ## if no errors
#                     if all([i.successful() for i in sliced[sample.name]]):                
#                         results = [i.get() for i in sliced[sample.name]]
#                         sample_cleanup(data, sample, results)
#                     ## print errors if they occurred
#                     else:
#                         for async in sliced[sample.name]:
#                             if not async.successful():
#                                 print("Error: %s", async.metadata.error)
#         if all([i.completed for i in asyncs]):
#             ## do final stats and cleanup
#             assembly_cleanup(data)
#         ## clean up concat files (if any)
#         concats = glob.glob(os.path.join(data.dirs.edits, "*_concat.fq"))
#         for concat in concats:
#             os.remove(concat)



if __name__ == "__main__":

    import ipyrad as ip

    ## get path to root (ipyrad) dir/ 
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__))))

    ## run tests
    TESTDIRS = ["test_rad", "test_pairgbs"]

    for tdir in TESTDIRS:
        TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", tdir, "data1"))
        TEST.step2(force=True)
        print(TEST.stats)

