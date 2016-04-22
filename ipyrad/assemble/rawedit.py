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
import time
import datetime
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
    """ 
    Filters one read for Illumina adapters and appends to a list for writing
    """

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
    return write1, write2, kept



def rawedit(args):
    """ 
    Applies two filters: Nfilter and adapterfilter.
    """
    ## parse args
    data, sample, num, nreplace, optim = args
    point = num*optim

    ## get cut sites 
    cuts1, cuts2 = [ambigcutters(i) for i in \
                    data.paramsdict["restriction_overhang"]]

    #LOGGER.info("cutsites %s %s", cuts1, cuts2)
    #LOGGER.info([i for i in data.paramsdict["restriction_overhang"]])

    ## get data slices as iterators and open file handles
    tups = sample.files.fastqs[0]
    fr1, fr2, io1, io2 = get_slice(tups, optim, num)

    ## make quarts iterator to sample (r1,r2) or (r1, 1)
    #tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    #read1 = os.path.join(tmpdir, "tmp_{}_{}_R1.fastq".format(sample.name, num))
    #read2 = os.path.join(tmpdir, "tmp_{}_{}_R2.fastq".format(sample.name, num))

    # fr1 = open(read1, 'rb')
    # quart1 = itertools.izip(*[iter(fr1)]*4)
    quart1 = itertools.izip(fr1, fr1, fr1, fr1)
    if "pair" in data.paramsdict["datatype"]:
        #fr2 = open(read2, 'rb')
        quart2 = itertools.izip(fr2, fr2, fr2, fr2)
        quarts = itertools.izip(quart1, quart2)
    else:
        quarts = itertools.izip(quart1, iter(int, 1))

    ## counters 
    counts = {}
    for key in ["orig", "quality", "adapter", "keep"]:
        counts[key] = 0

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
            args = [data, sample, read1, read2, 
                    cuts1, cuts2, write1, write2, point]
            write1, write2, kept = adapterfilter(args)

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
    #fr1.close()
    io1.close()

    ## write to file
    handle = os.path.join(data.dirs.edits, 
                          "tmp1_"+sample.name+"_"+str(point)+".fq")
    with open(handle, 'wb') as out1:
        out1.write("\n".join(write1))
        out1.write("\n")

    if "pair" in data.paramsdict["datatype"]:    
        #fr2.close()
        io2.close()
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



def roundup(num):
    """ round to nearest hundred """
    return int(math.ceil(num / 100.0)) * 100



    # ## set optim size, can only be optimized if reads_raw 
    # optim = 10000
    # ncpus = detect_cpus()
    # inputreads = int(data.stats.reads_raw.max())
    # if inputreads:
    #     ## Here optim is measured in reads.
    #     optim = (inputreads // (ncpus)) + (inputreads % (ncpus)) 
    #     ## multiply by 4 to get the number of lines
    #     optim = int(optim*4)
    #     LOGGER.info("optim=%s lines", optim)

    # ## truncate fastq if preview-mode. For step 2 preview mode we slice off
    # ## 4 optim size chunks so there is one chunk for each of the default
    # ## 4 ipyclients to work on.
    # if preview:
    #     nlines = data._hackersonly["preview_step2"]
    #     optim = int(nlines // ncpus())
    #     LOGGER.warn("""\
    # Running preview mode: subselecting only the first {} reads"""\
    # .format(nlines))

    ## tmpdir for slices
    # tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    # if not os.path.exists(tmpdir):
    #     os.mkdir(tmpdir)

    #return optim



# def send_slices(args):
#     """ 
#     Splits fastq file into smaller chunks.
#     """
#     ## parse args
#     data, sample, optim = args

#     ### location for slices
#     tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")

#     ## iterate over samples
#     num = 0
#     ## open file handles
#     tups = sample.files.fastqs[0]

#     if tups[0].endswith(".gz"):
#         io1 = gzip.open(tups[0])
#         rawr1 = iter(io1)
#         if tups[1]:
#             io2 = gzip.open(tups[1])
#             rawr2 = iter(io2)
#     else:
#         io1 = open(tups[0])
#         rawr1 = iter(io1)
#         if tups[1]:
#             io2 = open(tups[1])
#             rawr2 = iter(io2)

#     ## creating 10 slices per sample
#     for num in range(10):
#         dat1 = "".join(itertools.islice(rawr1, int(optim*4)))
#         if tups[1]:
#             dat2 = "".join(itertools.islice(rawr2, int(optim*4)))                

#         ## send slices
#         if dat1:
#             ## get outhandle
#             out1 = os.path.join(tmpdir, "tmp_{}_{}_R1.fastq"\
#                             .format(sample.name, num))
#             ## write slice
#             with open(out1, 'w') as tmpout:
#                 tmpout.write(dat1)
                
#             ## if paired-end
#             if tups[1]:
#                 ## get outhandle
#                 out2 = os.path.join(tmpdir, "tmp_{}_{}_R2.fastq"\
#                                 .format(sample.name, num))
#                 ## write slice
#                 with open(out2, 'w') as tmpout:
#                     tmpout.write(dat2)
#                 ## advance counter
#                 num += 1

#     io1.close()
#     if tups[1]:
#         io2.close()



def get_slice(tups, optim, jnum):
    """ 
    Slices a chunk from a fastq file and returns it as a list
    """
    ## open file handles
    if tups[0].endswith(".gz"):
        io1 = gzip.open(tups[0])
        rawr1 = iter(io1)
        if tups[1]:
            io2 = gzip.open(tups[1])
            rawr2 = iter(io2)
        else:
            io2 = 0
    else:
        io1 = open(tups[0])
        rawr1 = iter(io1)
        if tups[1]:
            io2 = open(tups[1])
            rawr2 = iter(io2)
        else:
            io2 = 0

    ## skip until jnum slice 
    ## sum takes about 5 seconds for 2M gzipped reads.
    skip = 0
    for _ in range(jnum):
        skip += sum(1 for i in itertools.islice(rawr1, int(optim*4)))
        if tups[1]:
            _ = sum(1 for i in itertools.islice(rawr2, int(optim*4)))
    LOGGER.info("%s. skipped %s lines", jnum, skip)

    ## now return the correct slice as a generator
    dat1 = itertools.islice(rawr1, int(optim*4))
    if tups[1]:
        dat2 = itertools.islice(rawr2, int(optim*4))
    else:
        dat2 = iter(int, 1)

    return dat1, dat2, io1, io2




# def fake():
#     ## Cleans up tmp files in the edits directory if rawedit passes,
#         ## but gets killed before sample_cleanup
#         tmpedit1 = glob.glob(os.path.join(data.dirs.edits,
#                     "tmp1_"+sample.name+"_*.fq"))
#         tmpedit2 = glob.glob(os.path.join(data.dirs.edits,
#                     "tmp2_"+sample.name+"_*.fq"))
#         try:
#             for i in tmpedit1:
#                 os.remove(i)
#             for i in tmpedit2:
#                 os.remove(i)
#         except:
#             pass
#         raise inst

#     finally:
#         ## if process failed at any point delete temp files
#         tmpdirs = glob.glob(os.path.join(data.dirs.project, 
#                                          data.name+"-tmpchunks"))
#         if tmpdirs:
#             for tmpdir in tmpdirs:
#                 try:
#                     shutil.rmtree(tmpdir)   
#                 except OSError:
#                     ## block again nfs error
#                     LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

#         if preview:
#             ## cleans up tmp files generated by preview
#             for tmpfile in sample_fastq[0]:
#                 os.remove(tmpfile)

#     return num, results



def sample_cleanup(data, sample, results):
    """ cleaning up for one sample """

    ## get tmp files and sort them 
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

    ## Concatenate tmp files
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
    for i in range(10):
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
    # tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    # if os.path.exists(tmpdir):
    #     try:
    #         shutil.rmtree(tmpdir)   
    #     except OSError:
    #         ## In some instances nfs creates hidden dot files in directories
    #         ## that claim to be "busy" when you try to remove them. Don't
    #         ## kill the run if you can't remove this directory.
    #         LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

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
    Skipping Sample {}; Already filtered. Use force argument to overwrite.""".\
    format(sample.name))
            elif not sample.stats.reads_raw:
                print("""
    Skipping Sample {}; No reads found in file {}"""\
    .format(sample.name, sample.files.fastqs))
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

    ## get optim which is used to slice a sample into 10 equal sized chunks
    ## if preview-mode then 10 chunks are made out of 'preview' reads.
    optims = {}
    for sample in subsamples:
        if preview:
            tots = data._hackersonly["preview_step2"]
            assert not tots % 4, \
            "_hackersonly preview_step2 value (nlines) must be divisible by 4."
            optims[sample.name] = (tots // 10) + (tots % 10)
        else:
            tots = int(sample.stats.reads_raw)
            optims[sample.name] = (tots // 10) + (tots % 10)

    if preview:
        if data._headers:
            print("""\
    Running preview mode: subselecting maximum of {} reads per sample\
    """.format(data._hackersonly["preview_step2"]))

    ## link output directories 
    data.dirs.edits = os.path.join(os.path.realpath(
                                   data.paramsdict["project_dir"]), 
                                   data.name+"_edits")

    ## create output directory if doesn't exist
    if not os.path.exists(data.dirs.edits):
        os.makedirs(data.dirs.edits)

    ## client
    lbview = ipyclient.load_balanced_view()

    ## start progress
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(0))
    progressbar(len(subsamples), 0, 
               " processing reads  | {}".format(elapsed))            

    ## save sliced asyncs
    sliced = {i.name:[] for i in subsamples}    
    ## send jobs to queue to get slices and process them
    for sample in subsamples:
        ## get optim slice size for this sample
        optim = optims[sample.name]        
        for job in range(10):
            args = [data, sample, job, nreplace, optim]
            async = lbview.apply(rawedit, args)
            sliced[sample.name].append(async)

    ## print progress
    try:
        while 1:
            asyncs = list(itertools.chain(*sliced.values()))
            tots = len(asyncs)
            done = sum([i.ready() for i in asyncs])
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            if tots != done:
                progressbar(tots, done,
                           " processing reads  | {}".format(elapsed))            
                time.sleep(1)
            else:
                progressbar(tots, done,
                           " processing reads  | {}".format(elapsed))
                if data._headers:
                    print("")
                break

    except (KeyboardInterrupt, SystemExit):
        print('\n  Interrupted! Cleaning up... ')
        raise #raise IPyradWarningExit("Keyboard Interrupt")

    ## enforced cleanup
    finally:
        ## if all jobs were successful in a sample then cleanup
        for sample in subsamples:
            ## if finished
            if all([i.ready() for i in sliced[sample.name]]):
                ## if no errors
                if all([i.successful() for i in sliced[sample.name]]):                
                    results = [i.get() for i in sliced[sample.name]]
                    sample_cleanup(data, sample, results)
                ## print errors if they occurred
                else:
                    for async in sliced[sample.name]:
                        if not async.successful():
                            print("Error: %s", async.metadata.error)

        ## do final stats and cleanup
        #assembly_cleanup(data)



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

