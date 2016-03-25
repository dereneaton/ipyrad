#!/usr/bin/env ipython2

""" edits reads based on quality scores. Can be used
to check for adapters and primers, but is not optimized 
for all types of cutters """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=W0142

import os
import glob
import gzip
import math
import itertools
import numpy as np
import shutil
from .util import *

import logging
LOGGER = logging.getLogger(__name__)



## define globals for afilter

def afilter(data, sample, bases, cuts1, cuts2, read):
    """ 
    Applies filter for primers & adapters. Does three checks in order from 
    left to right of seq. Checks cut+adapter, adapter, and part of adapter
    """
    ## set empty
    where1 = where2 = where3 = check1 = check2 = None
    adapter1 = "AGATCGG"
    adapter2 = "AGAGCGT"

    ## if ddrad we're looking for the second cutter on read1
    rvcuts = [comp(i)[::-1] for i in cuts1]
    if "ddrad" in data.paramsdict["datatype"]:
        rvcuts = [comp(i)[::-1] for i in cuts2]

    ## if not ambiguous cuts then make second null (Zs)
    if not rvcuts[1]:
        rvcuts[1] = "Z"*4

    ## only look for cut sites if there are cut sites
    if any(rvcuts):
        ## Look for adapters -- rvcut+[adapter]
        if read == 1:
            lookfor1 = rvcuts[0]+"AGA"
            lookfor2 = rvcuts[1]+"AGA"
        else:
            ## read 2 will have the barcode between: rvcut+[barcode]+[adapter]
            ## and so should only be checked for if barcode info is available
            if sample.name in data.barcodes:
                barcode = data.barcodes[sample.name]
                lookfor1 = rvcuts[0]+comp(barcode)[::-1][:3]
                lookfor2 = rvcuts[1]+comp(barcode)[::-1][:3]
            else:
                ## this gets cancelled out if there's no barcodes info
                ## unless setting is strict (2)
                lookfor1 = rvcuts[0]+"NNN"
                lookfor2 = rvcuts[1]+"NNN"

        ## if strict then shorter lookfor
        if data.paramsdict["filter_adapters"] == 2:
            lookfor1 = lookfor1[:-3]
            lookfor2 = lookfor2[:-3]
        ## if cutter has ambiguous base (2 in cuts) look for both 
        ## otherwise just one
        try: 
            check1 = max(0, bases[1].tostring().find(lookfor1))
        except Exception as inst:
            LOGGER.error("EXCEPTION %s", [bases, lookfor1, inst])
        
        ### looks for second resolution
        if sum([1 for i in rvcuts]) == 2:
            check2 = max(0, bases[1].tostring().find(lookfor2))
        ### CHECK FOR FIRST FILTER
        if check1 or check2:
            where1 = min([i for i in [check1, check2] if i])

    #LOGGER.debug("where1:%s, ch1:%s, ch2:%s, read:%s", 
    #              where1, check1, check2, read)

    ## look for adapter sequence directly in two parts: "AGATCGGA.AGAGCGTC"
    ## if strict then shorten the lookfor string
    dist = None
    if data.paramsdict["filter_adapters"] == 2:
        dist = -2    
    lookfor1 = adapter1[:dist]
    lookfor2 = adapter2[:dist]

    check1 = max(0, bases[1].tostring().find(lookfor1))
    check2 = max(0, bases[1].tostring().find(lookfor2))
    ### CHECK FOR SECOND FILTER
    mincheck = min(check1, check2)

    ## How far back from adapter to trim to remove cutsite and barcodes
    if mincheck:
        ## trim further for check2
        backup = 1
        if not check1:
            backup = 9
        ## find where to trim
        if read == 1:
            where2 = max(0, mincheck - (len(cuts2[0]) + backup))
        else:
            if "ddrad" in data.paramsdict["datatype"]:
                trim = mincheck - ((len(cuts1[0]) + len(cuts2[0])) + backup)
                where2 = max(0, trim)
            else:
                trim = mincheck - ((len(cuts1[0]) + len(cuts1[0])) + backup)
                where2 = max(0, trim)
    else:
        where2 = 0

    #LOGGER.debug("where2:%s, ch1:%s, ch2:%s, read:%s",
    #              where2, check1, check2, read)

    ## CHECK FOR THIRD FILTER
    ## if strict filter, do additional search for cut site near edges
    if any(rvcuts):
        where3 = 0
        if data.paramsdict["filter_adapters"] == 2:
            if not (where1 or where2):
                if read == 1:
                    cutback = -15
                    tail = bases[1].tostring()[cutback:]
                    if any([i in tail for i in rvcuts]):
                        where3 = len(bases[1]) + cutback
                else:
                    cutback = -10
                    tail = bases[1].tostring()[cutback:]
                    if any([i in tail for i in rvcuts]):
                        where3 = len(bases[1]) + cutback
    
    #LOGGER.debug("where3:%s, ......, ......, read:%s", where3, read)
    try:
        where = min([i for i in [where1, where2, where3] if i])
    except (TypeError, ValueError):
        where = 0
    return where 



def adapterfilter(args):
    """ filters for adapters """

    ## parse args
    data, sample, read1, read2, cuts1, cuts2, write1, write2, point = args

    ## look for adapters in sequences
    wheretocut1 = wheretocut2 = cutter = None

    ## keep counter
    kept = 1

    ## check for adapter
    if data.paramsdict["filter_adapters"]:
        ## get trim length for read1
        wheretocut1 = afilter(data, sample, read1, cuts1, cuts2, 1)
        if 'pair' in data.paramsdict["datatype"]:
            wheretocut2 = afilter(data, sample, read2, cuts1, cuts2, 2)

        ## cut one or both reads depending on detection of adapters
        cutter = 0
        if wheretocut1 and wheretocut2:
            cutter = min(wheretocut1, wheretocut2)
        elif wheretocut1 or wheretocut2:
            cutter = max(wheretocut1, wheretocut2)

    ## pairgbs need to be trimmed to the same length. This means it will
    ## always have _c in read names
    if data.paramsdict['datatype'] == 'pairgbs':
        readlens = [len(i) for i in (read1[1], read2[1])]
        cutter = min([i for i in readlens+[cutter] if i])
    #LOGGER.debug("@cutter w1:%s w2:%s ", wheretocut1, wheretocut2)

    if cutter:
        ## if trimmed frag is still long enough
        if cutter >= max(32, data.paramsdict["filter_min_trim_len"]):
            ## write fastq format
            sseq1 = "\n".join(["@"+sample.name+"_"+str(point)+"_c1", 
                               read1[1].tostring()[:cutter],
                               "+",
                               #read1[2].tostring(),
                               read1[3].tostring()[:cutter]])
            write1.append(sseq1)
            #LOGGER.debug("%s", read1[1].tostring())
            #LOGGER.debug("%s -trim r1", read1[1].tostring()[:cutter])

            if len(read2):
                sseq2 = "\n".join(["@"+sample.name+"_"+str(point)+"_c2",
                                   read2[1].tostring()[:cutter],
                                   "+",
                                   #read2[2].tostring(),
                                   read2[3].tostring()[:cutter]])
                                   #bases1[:cutter]+"SS"+bases2[:cutter]+"\n"
                write2.append(sseq2)
                #LOGGER.debug("%s", read2[1].tostring())
                #LOGGER.debug("%s -trim r2", read2[1].tostring()[:cutter])                
        else:
            kept = 0

    ## not trimmed
    else:
        sseq1 = "\n".join(["@"+sample.name+"_"+str(point)+"_r1", 
                           read1[1].tostring(),
                           "+", 
                           #read1[2].tostring(),
                           read1[3].tostring()])
        write1.append(sseq1)
        if len(read2):
            sseq2 = "\n".join(["@"+sample.name+"_"+str(point)+"_r2", 
                               read2[1].tostring(),
                               "+",
                               #read2[2].tostring(),
                               read2[3].tostring()])
            write2.append(sseq2)
    return write1, write2, read1, read2, point, kept



def rawedit(args):
    """ 
    Applies two filters: Nfilter and adapterfilter.
    """
    ## parse args
    data, sample, tmptuple, nreplace, point = args

    ## get cut sites 
    cuts1, cuts2 = [ambigcutters(i) for i in \
                    data.paramsdict["restriction_overhang"]]
    #LOGGER.debug("cutsites %s %s", cuts1, cuts2)
    #LOGGER.debug([i for i in data.paramsdict["restriction_overhang"]])

    ## the read1 demultiplexed reads file
    if tmptuple[0].endswith(".gz"):
        fr1 = gzip.open(os.path.realpath(tmptuple[0]), 'rb')
    else:
        fr1 = open(os.path.realpath(tmptuple[0]), 'rb')

    ## the read2 demultiplexed reads file, if paired
    if "pair" in data.paramsdict["datatype"]:
        if tmptuple[1].endswith(".gz"):
            fr2 = gzip.open(os.path.realpath(tmptuple[1]), 'rb')
        else:
            fr2 = open(os.path.realpath(tmptuple[1]), 'rb')

    ## create iterators to sample 4 lines at a time 
    quart1 = itertools.izip(*[iter(fr1)]*4)
    if "pair" in data.paramsdict["datatype"]:
        ## pair second read files, quarts samples both
        quart2 = itertools.izip(*[iter(fr2)]*4)
        quarts = itertools.izip(quart1, quart2)
    else:
        ## read in single end read file, quarts samples 1 for other
        quarts = itertools.izip(quart1, iter(int, 1))

    ## counters 
    counts = {}
    counts["orig"] = 0
    counts["quality"] = 0
    counts["adapter"] = 0    
    counts["keep"] = 0    

    ## in memory storage to write to file
    write1 = []
    write2 = []

    ## apply filters to each read (pair)
    while 1:
        try: 
            quart = quarts.next()
        except StopIteration: 
            break

        ## strip reads
        counts['orig'] += 1 
        read1 = np.array([np.array(list(i.strip())) for i in quart[0]])
        read2 = []        
        if "pair" in data.paramsdict["datatype"]:
            read2 = np.array([np.array(list(i.strip())) for i in quart[1]])

        ## filter for Ns
        ## base1 is only returned if both passed for paired-end reads        
        read1, read2 = nfilter(data, read1, read2, nreplace)

        if len(read1):
            ## replace cut sites unless cuts already removed
            if any(data.paramsdict["edit_cutsites"]):
                read1, read2 = modify_cuts(data, read1, read2)

            ## filter for adapters, and trim to size if pairgbs
            args = [data, sample, read1, read2, cuts1, cuts2, 
                    write1, write2, point]
            write1, write2, read1, read2, point, kept = adapterfilter(args)
            #LOGGER.debug("down HERE")

            ## counters
            if kept:
                counts["keep"] += 1
            else:
                counts["adapter"] += 1
        else:
            counts["quality"] += 1

        ## advance number for read names
        point += 1

    ## close file handles
    fr1.close()
    ## write to file
    handle = os.path.join(data.dirs.edits, 
                          "tmp1_"+sample.name+"_"+str(point)+".fq")
    with open(handle, 'wb') as out1:
        out1.write("\n".join(write1))
        out1.write("\n")

    if "pair" in data.paramsdict["datatype"]:    
        fr2.close()
        handle = os.path.join(data.dirs.edits, 
                          "tmp2_"+sample.name+"_"+str(point)+".fq")
        with open(handle, 'wb') as out2:
            out2.write("\n".join(write2))
            out2.write("\n")            

    ## the number of ...
    return counts



def modify_cuts(data, read1, read2):
    """ fix cut sites to be error free and not carry ambiguities """
    ## shorthand
    cutsmod = data.paramsdict["edit_cutsites"]

    ## replace cut region with a single resolution of cut w/o ambiguities
    if len(read1) and cutsmod[0]:
        if isinstance(cutsmod[0], int):
            read1[1] = read1[1][abs(cutsmod[0]):]      
            read1[3] = read1[3][abs(cutsmod[0]):]        
        elif isinstance(cutsmod[0], str):
            read1[1][:len(cutsmod[0])] = list(cutsmod[0])
            read1[3][:len(cutsmod[0])] = ["B"]*len(cutsmod[0])

    ## same for cut2 and end of second read
    if len(read2) and cutsmod[1]:
        ## fix cut sites to be error free before counting Ns
        #LOGGER.debug("before R2: %s", "".join(read2[1]))
        if isinstance(cutsmod[1], int):
            read2[1] = read2[1][abs(cutsmod[1]):]   
            read2[3] = read2[3][abs(cutsmod[1]):]
        elif isinstance(cutsmod[1], str):
            read2[1][:len(cutsmod[1])] = list(cutsmod[1])
            read2[3][:len(cutsmod[1])] = ["B"]*len(cutsmod[1])
        #LOGGER.debug("after_R2: %s", "".join(read2[1]))

    return read1, read2



def nfilter(data, read1, read2, nreplace=True):
    """ 
    Filters based on max allowed Ns, and either replaces the Ns, or leaves 
    them. Right now replaces. Also fixes cut site to be error free.
    """
    ## get qscores
    qscore = np.array([ord(i) for i in read1[3]])
    phred1 = qscore - data.paramsdict["phred_Qscore_offset"]

    ## replace low qual with N
    if nreplace:
        read1[1][phred1 < 20] = "N"

    ## maxN filter only applies to first reads if not paired
    if not "pair" in data.paramsdict["datatype"]:
        if sum(phred1 < 20) < data.paramsdict["max_low_qual_bases"]:
            return read1, np.array([])
        else:
            return np.array([]), np.array([])

    else:
        ## get qscores
        qscore = np.array([ord(i) for i in read2[3]])
        phred2 = qscore - data.paramsdict["phred_Qscore_offset"]

        ## replace low qual with N
        if nreplace:
            read2[1][phred2 < 20] = "N"

        ## return passed 
        if max([sum(phred1 < 20), sum(phred2 < 20)]) \
               < data.paramsdict["max_low_qual_bases"]:
            return read1, read2
        else:
            return np.array([]), np.array([])



def prechecks(data):
    """ checks before starting analysis """
    ## link output directories 
    data.dirs.edits = os.path.join(os.path.realpath(
                                   data.paramsdict["project_dir"]), 
                                   data.name+"_edits")
    ## create output directory if doesn't exist
    if not os.path.exists(data.dirs.edits):
        os.makedirs(data.dirs.edits)



def roundup(num):
    """ round to nearest hundred """
    return int(math.ceil(num / 100.0)) * 100



def run_sample(data, sample, nreplace, preview, ipyclient):
    """ 
    Splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them .
    """
    ## before starting
    prechecks(data)
    num = 0

    ## set optim size, can only be optimized if reads_raw 
    optim = 10000
    if sample.stats.reads_raw:
        ### split reads among N processors
        inputreads = sample.stats.reads_raw
        ncpus = len(ipyclient)
        ## The goal is to get the optimum number of lines per chunk so
        ## as to split into 2 * ncpus chunks. Devide number of input reads
        ## by 2*ncpus, then multipy by 4 to get the number of lines
        ## Cast to an int, since numreads should be a whole number and
        ## zcat_make_temps doesn't like the '100.0' float part.
        ## Here optim is measured in nreads.
        optim = int(inputreads // (ncpus * 2))
        ## multiply by 4 to ensure fastq quartet sampling
        ## The unit for optim is now nlines now
        optim *= 4
        LOGGER.info("optim=%s", optim)

    ## truncate fastq if preview-mode. For step 2 preview mode we slice off
    ## 4 optim size chunks so there is one chunk for each of the default
    ## 4 ipyclients to work on.
    if preview:
        nlines = data._hackersonly["preview_step2"]
        optim = nlines // len(ipyclient)
        LOGGER.warn("""\
    Running preview mode: subselecting {} reads from {}
    """.format(nlines, sample.files.fastqs))
        sample_fastq = preview_truncate_fq(data, sample.files.fastqs, nlines)
    else:
        sample_fastq = sample.files.fastqs

    LOGGER.debug("sample_fastq back here is %s", sample_fastq)
    ## break up the file into smaller tmp files for each processor
    tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    chunkslist = []
    for fastqtuple in sample_fastq:
        args = [data, fastqtuple, num, tmpdir, optim]
        _, achunk = zcat_make_temps(args)
        chunkslist += achunk

    LOGGER.info("Executing %s file, in %s chunks, across %s cpus", \
                 len(sample_fastq), len(chunkslist), len(ipyclient))

    ## send chunks across processors, will delete if fail
    try:
        submitted_args = []
        for tmptuple in chunkslist:
            ## point is used to increment names across processors
            point = num*(optim // 4)  
            submitted_args.append([data, sample, tmptuple, nreplace, point])
            num += 1

        ## create view to auto spread chunks across engines
        lbview = ipyclient.load_balanced_view()
        results = lbview.map_async(rawedit, submitted_args)
        results.get()

    finally:
        ## if process failed at any point delete temp files
        tmpdirs = glob.glob(os.path.join(data.dirs.project, 
                                         data.name+"-tmpchunks"))
        if tmpdirs:
            for tmpdir in tmpdirs:
                try:
                    shutil.rmtree(tmpdir)   
                except OSError:
                    ## block again nfs error
                    LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

        ## Cleans up tmp files in the edits directory if rawedit passes,
        ## but gets killed before sample_cleanup
        tmpedit1 = glob.glob(os.path.join(data.dirs.edits,
                    "tmp1_"+sample.name+"_*.fq"))
        tmpedit2 = glob.glob(os.path.join(data.dirs.edits,
                    "tmp2_"+sample.name+"_*.fq"))
        try:
            for i in tmpedit1:
                os.remove(i)
            for i in tmpedit2:
                os.remove(i)
        except:
            pass

        if preview:
            ## cleans up tmp files generated by preview
            for tmpfile in sample_fastq[0]:
                os.remove(tmpfile)

    return num, results



def sample_cleanup(data, sample, submitted, results):
    """ cleaning up for one sample """

    ## rejoin chunks
    combs1 = glob.glob(os.path.join(
                        data.dirs.edits,
                        "tmp1_"+sample.name+"_*.fq"))
    combs1.sort(key=lambda x: int(x.split("_")[-1].replace(".fq", "")))
    handle1 = os.path.join(data.dirs.edits, sample.name+"_R1_.fastq")
    handle2 = ""

    ## same for pairs
    if "pair" in data.paramsdict["datatype"]:
        combs2 = glob.glob(os.path.join(
                            data.dirs.edits,
                            "tmp2_"+sample.name+"_*.fq"))
        combs2.sort(key=lambda x: int(x.split("_")[-1].replace(".fq", "")))    
        handle2 = os.path.join(data.dirs.edits, sample.name+"_R2_.fastq")        
        assert len(combs1) == len(combs2), \
            "mismatched number of paired read files - {} {}"\
            .format(len(combs1), len(combs2))

    ## Clean up tmp files
    with open(handle1, 'wb') as out:
        for fname in combs1:
            with open(fname) as infile:
                out.write(infile.read())
            os.remove(fname)

    if "pair" in data.paramsdict["datatype"]:
        with open(handle2, 'wb') as out:
            for fname in combs2:
                with open(fname) as infile:
                    out.write(infile.read())
                os.remove(fname)

    ## record results
    fcounts = {"orig": 0,
               "quality": 0,
               "adapter": 0, 
               "keep": 0}

    ## merge finished edits
    for i in range(submitted):
        counts = results[i]
        fcounts["orig"] += counts["orig"]
        fcounts["quality"] += counts["quality"]
        fcounts["adapter"] += counts["adapter"]
        fcounts["keep"] += counts["keep"]

    ## store full s2 stats
    sample.stats_dfs.s2["reads_raw"] = fcounts["orig"]
    sample.stats_dfs.s2["filtered_by_qscore"] = fcounts["quality"]
    sample.stats_dfs.s2["filtered_by_adapter"] = fcounts["adapter"]
    sample.stats_dfs.s2["reads_passed"] = fcounts["keep"]

    ## store summary stats
    sample.stats.reads_filtered = fcounts["keep"]

    ## save stats to Sample if successful
    if sample.stats.reads_filtered:
        sample.stats.state = 2
        sample.files.edits = [(handle1, handle2)]
    else:
        print("No reads passed filtering in Sample: {}".format(sample.name))



def assembly_cleanup(data):
    """ cleanup for assembly object """

    ## remove tmpchunks dir
    tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    if os.path.exists(tmpdir):
        try:
            shutil.rmtree(tmpdir)   
        except OSError:
            ## In some instances nfs creates hidden dot files in directories
            ## that claim to be "busy" when you try to remove them. Don't
            ## kill the run if you can't remove this directory.
            LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

    ## build s2 results data frame
    data.stats_dfs.s2 = data.build_stat("s2")
    data.stats_files.s2 = os.path.join(data.dirs.edits, 's2_rawedit_stats.txt')

    ## write stats for all samples
    with open(data.stats_files.s2, 'w') as outfile:
        data.stats_dfs.s2.to_string(outfile)



def run(data, samples, nreplace, force, preview, ipyclient):
    """ run the major functions for editing raw reads """

    ## hold samples that pass
    subsamples = []

    ## filter the samples again
    if not force:
        for sample in samples:
            if sample.stats.state >= 2:
                print("""
    Skipping Sample {}; Already filtered. Use force=True to overwrite.""".\
    format(sample.name))

            elif not sample.stats.reads_raw:
                print("""
    Skipping Sample {}; Too few reads ({}). Use force=True to run anyways.""".\
    format(sample.name, sample.stats.reads_raw))

            else:
                subsamples.append(sample)

    else:
        for sample in samples:
            if not sample.stats.reads_raw:
                print("""
    Skipping Sample {}; No reads found in file {}""".\
    format(sample.name, sample.files.fastqs))

            else:
                subsamples.append(sample)

    ## print to screen if preview method
    if preview:
        nlines = data._hackersonly["preview_step2"]
        if data._headers:
            print("""\
    Running preview mode: subselecting maximum of {} reads per sample\
    """.format(nlines))

    ## RUN THE SAMPLES
    try:
        for njob, sample in enumerate(subsamples):
            sub, res = run_sample(data, sample, nreplace, preview, ipyclient)
            sample_cleanup(data, sample, sub, res)
            progressbar(len(subsamples), njob)
        progressbar(len(subsamples), len(subsamples))
        print("")
    ## ensure the assembly is cleaned up even if stopped
    finally:
        assembly_cleanup(data)



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

