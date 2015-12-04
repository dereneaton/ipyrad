#!/usr/bin/env ipython2

""" edits reads based on quality scores. Can be used
to check for adapters and primers, but is not optimized 
for all types of cutters """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212

import os
import glob
import gzip
import math
import time
import itertools
import numpy as np
import ipyparallel as ipp
from .demultiplex import ambigcutters
from .demultiplex import zcat_make_temps

import logging
LOGGER = logging.getLogger(__name__)


def afilter(data, sample, bases, cuts1, cuts2, read):
    """ 
    Applies filter for primers & adapters. Does three checks in order from 
    left to right of seq. Checks cut+adapter, adapter, and part of adapter
    """
    ## set empty
    where1 = where2 = where3 = check1 = check2 = None

    ## if ddrad we're looking for the second cutter on read1
    rvcuts = [comp(i)[::-1] for i in cuts1]
    if "ddrad" in data.paramsdict["datatype"]:
        rvcuts = [comp(i)[::-1] for i in cuts2]

    ## if not ambiguous cuts then make second Zs
    if not rvcuts[1]:
        rvcuts[1] = "Z"*10

    ## Look for adapters -- rvcut+[adapter]
    if read == 1:
        lookfor1 = rvcuts[0]+"AGA"
        lookfor2 = rvcuts[1]+"AGA"
    else:
        ## read 2 will have the barcode between: rvcut+[barcode]+[adapter]
        if sample.barcode:
            lookfor1 = rvcuts[0]+comp(sample.barcode)[::-1][:3]
            lookfor2 = rvcuts[1]+comp(sample.barcode)[::-1][:3]
        else:
            lookfor1 = rvcuts[0]+"NN"
            lookfor2 = rvcuts[1]+"NN"

    ## if strict then shorter lookfor
    if data.paramsdict["filter_adapters"] == 2:
        lookfor1 = lookfor1[:-2]
        lookfor2 = lookfor2[:-2]

    ## if cutter has ambiguous base (2 in cuts) look for both otherwise just one
    try: 
        check1 = max(0, bases[1].tostring().rfind(lookfor1))
    except Exception as inst:
        LOGGER.debug([bases, lookfor1, inst])
    if sum([1 for i in rvcuts]) == 2:
        check2 = max(0, bases[1].tostring().rfind(lookfor2))

    if check1 or check2:
        where1 = min([i for i in [check1, check2] if i])

    LOGGER.debug("where1:%s, ch1:%s, ch2:%s, read:%s", 
                 where1, check1, check2, read)

    ## look for adapter sequence directly in two parts: "AGATCGGA.AGAGCGTC"
    lookfor1 = "AGATCGGA"
    lookfor2 = "AGAGCGTC"

    ## if strict then shorter lookfor
    if data.paramsdict["filter_adapters"] == 2:        
        lookfor1 = lookfor1[:-2]
        lookfor2 = lookfor2[:-2]

    check1 = max(0, bases[1].tostring().rfind(lookfor1))
    check2 = max(0, bases[1].tostring().rfind(lookfor2))
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

    LOGGER.debug("where2:%s, ch1:%s, ch2:%s, read:%s",
                  where2, check1, check2, read)

    ## if strict filter, do additional search for cut site near edges
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
    
    LOGGER.debug("where3:%s, ......, ......, read:%s", where3, read)
    try:
        where = min([i for i in [where1, where2, where3] if i])
    except (TypeError, ValueError):
        where = 0
    return where 


def comp(seq):
    """ returns a seq with small complement"""
    return seq.replace("A", 't')\
           .replace('T', 'a')\
           .replace('C', 'g')\
           .replace('G', 'c')\
           .replace('n', 'Z')\
           .upper().replace("Z", "n")



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

    ## if the read needs to be trimmed
    LOGGER.debug("@cutter w1:%s w2:%s ", wheretocut1, wheretocut2)

    if cutter:
        ## if trimmed frag is still long enough
        if cutter > max(32, data.paramsdict["filter_min_trim_len"]):
            ## write fastq format
            sseq1 = "\n".join(["@"+sample.name+"_"+str(point)+"_c1", 
                               read1[1].tostring()[:cutter],
                               read1[2].tostring(),
                               read1[3].tostring()[:cutter]])
            write1.append(sseq1)
            LOGGER.debug("%s", read1[1].tostring())
            LOGGER.debug("%s -trim r1", read1[1].tostring()[:cutter])

            if len(read2):
                sseq2 = "\n".join(["@"+sample.name+"_"+str(point)+"_c2",
                                   read2[1].tostring()[:cutter],
                                   read2[2].tostring(),
                                   read2[3].tostring()[:cutter]])
                                   #bases1[:cutter]+"SS"+bases2[:cutter]+"\n"
                write2.append(sseq2)
                LOGGER.debug("%s", read2[1].tostring())
                LOGGER.debug("%s -trim r2", read2[1].tostring()[:cutter])                
        else:
            kept = 0

    ## not trimmed
    else:
        sseq1 = "\n".join(["@"+sample.name+"_"+str(point)+"_r1", 
                           read1[1].tostring(),
                           read1[2].tostring(),
                           read1[3].tostring()])
        write1.append(sseq1)
        if len(read2):
            sseq2 = "\n".join(["@"+sample.name+"_"+str(point)+"_r2", 
                               read2[1].tostring(),
                               read2[2].tostring(),
                               read2[3].tostring()])
            write2.append(sseq2)

    return write1, write2, read1, read2, point, kept



def rawedit(args):
    """ 
    Applies two filters: Nfilter and adapterfilter.
    """

    ## get args
    data, sample, tmptuple, nreplace, point = args

    ## get cut sites
    LOGGER.debug(data.paramsdict["restriction_overhang"])
    LOGGER.debug([i for i in data.paramsdict["restriction_overhang"]])
    cuts1, cuts2 = [ambigcutters(i) for i in \
                    data.paramsdict["restriction_overhang"]]
    LOGGER.info("cutsites %s %s", cuts1, cuts2)

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
            ## replace cut sites            
            read1, read2 = replace_cuts(data, read1, read2, cuts1, cuts2)

            ## filter for adapters 
            args = [data, sample, read1, read2, cuts1, cuts2, 
                    write1, write2, point]
            write1, write2, read1, read2, point, kept = adapterfilter(args)
            if kept:
                counts["keep"] += 1
            else:
                counts["adapter"] += 1
        else:
            counts["quality"] += 1

        ## TODO: print to file after N to retain low memory load


        ## advance number for read names
        point += 1

    ## close file handles
    fr1.close()
    ## write to file
    handle = os.path.join(data.dirs.edits, 
                          "tmp1_"+sample.name+"_"+str(point)+".gz")
    with gzip.open(handle, 'wb') as out1:
        out1.write("\n".join(write1))
        out1.write("\n")

    if "pair" in data.paramsdict["datatype"]:    
        fr2.close()
        handle = os.path.join(data.dirs.edits, 
                          "tmp2_"+sample.name+"_"+str(point)+".gz")
        with gzip.open(handle, 'wb') as out2:
            out2.write("\n".join(write2))
            out2.write("\n")            

    ## the number of ...
    return counts



def replace_cuts(data, read1, read2, cuts1, cuts2):
    """ fix cut sites to be error free and not carry ambiguities """

    ## replace cut region with a single resolution of cut w/o ambiguities
    read1[1][:len(cuts1[0])] = list(cuts1[0])
    read1[3][:len(cuts1[0])] = ["B"]*len(cuts1[0])
    #read1[1] = read1[1][:len(cut1[0])] 

    ## same for cut2 and end of second read
    if len(read2):
        ## fix cut sites to be error free before counting Ns
        if 'gbs' in data.paramsdict["datatype"]:
            ## single digest use overhang for enzyme 1
            read2[1][:len(cuts1[0])] = list(cuts1[0])
            read2[3][:len(cuts1[0])] = ["B"]*len(cuts1[0])
            #read2[1] = read2[1][len(cut1[0]):] # = list(cut1[0])            
        else:
            ## double digest use overhang for enzyme 2
            read2[1][:len(cuts2[0])] = list(cuts2[0])
            read2[3][:len(cuts2[0])] = ["B"]*len(cuts2[0])
            #read2[1] = read2[1][len(cut2[0]):] # = list(cut2[0])            
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

        ## reverse complement as string
        ## bases2 = comp(bases2.tostring())[::-1]

        ## return passed 
        if max([sum(phred1 < 20), sum(phred2 < 20)]) \
               < data.paramsdict["max_low_qual_bases"]:
            return read1, read2
        else:
            return np.array([]), np.array([])



def prechecks(data, preview):
    """ checks before starting analysis """
    ## link output directories 
    data.dirs.edits = os.path.join(os.path.realpath(
                                   data.paramsdict["working_directory"]), 
                                   data.name+"_edits")
    ## create output directory if doesn't exist
    if not os.path.exists(data.dirs.edits):
        os.makedirs(data.dirs.edits)
    ## preview
    if preview:
        print("created new directory: {}".format(data.dirs.edits))



def roundup(num):
    """ round to nearest hundred """
    return int(math.ceil(num / 100.0)) * 100



def run_full(data, sample, ipyclient, nreplace):
    """ 
    Splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them .
    """
    ## before starting
    preview = 0
    prechecks(data, preview)

    ## load up work queue
    submitted = 0
    num = 0

    ## set optim size, can only be optimized if reads_raw 
    optim = 10000
    if sample.stats.reads_raw:
        if sample.stats.reads_raw > 1e5:
            optim = roundup(sample.stats.reads_raw/len(ipyclient.ids))*4

    ## break up the file into smaller tmp files for each processor
    chunkslist = []
    for fastqtuple in sample.files.fastqs:
        args = [data, fastqtuple, num, optim]
        _, achunk = zcat_make_temps(args)
        chunkslist += achunk

    ## send chunks across processors, will delete if fail
    try: 
        submitted_args = []
        for tmptuple in chunkslist:
            ## used to increment names across processors
            point = num*optim   #10000 #num*(chunksize/2)
            args = [data, sample, tmptuple, nreplace, point]
            submitted_args.append(args)
            submitted += 1
            num += 1

        ## first launch of ipyclient
        lbview = ipyclient.load_balanced_view()
        results = lbview.map_async(rawedit, submitted_args)

        ## return errors if an engine fails
        try:
            results.get()
        except (AttributeError, TypeError):
            for key in ipyclient.history:
                if ipyclient.metadata[key].error:
                    LOGGER.error("step2 error: %s", 
                        ipyclient.metadata[key].error)
                    raise SystemExit
                if ipyclient.metadata[key].stdout:
                    LOGGER.error("step2 stdout:%s", 
                        ipyclient.metadata[key].stdout)
                    raise SystemExit            
        del lbview
    
    finally:
        ## if process failed at any point delete temp files
        for tmptuple in chunkslist:
            os.remove(tmptuple[0])
            if "pair" in data.paramsdict["datatype"]:
                os.remove(tmptuple[1]) 

    return submitted, results



def cleanup(data, sample, submitted, results):
    """ cleaning up """

    ## rejoin chunks
    combs1 = glob.glob(os.path.join(
                        data.dirs.edits,
                        "tmp1_"+sample.name+"_*.gz"))
    combs1.sort(key=lambda x: int(x.split("_")[-1].replace(".gz", "")))
    handle1 = os.path.join(data.dirs.edits, sample.name+"_R1_.fastq.gz")
    handle2 = ""

    with gzip.open(handle1, 'wb') as out:
        for fname in combs1:
            with gzip.open(fname) as infile:
                out.write(infile.read())
            os.remove(fname)

    if "pair" in data.paramsdict["datatype"]:
        combs2 = glob.glob(os.path.join(
                            data.dirs.edits,
                            "tmp2_"+sample.name+"_*.gz"))
        combs2.sort(key=lambda x: int(x.split("_")[-1].replace(".gz", "")))    
        handle2 = os.path.join(data.dirs.edits, sample.name+"_R2_.fastq.gz")        
        assert len(combs1) == len(combs2), \
               "mismatched number of paired read files"

        with gzip.open(handle2, 'wb') as out:
            for fname in combs2:
                with gzip.open(fname) as infile:
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

    data.statsfiles.s2 = os.path.join(data.dirs.edits, 's2_rawedit_stats.txt')

    ## find longest name to make printing code block
    longestname = max([len(i) for i in data.samples.keys()])
    printblock = "{:<%d}  {:>13} {:>13} {:>13} {:>13}\n" % (longestname + 4)

    if not os.path.exists(data.statsfiles.s2):
        with open(data.statsfiles.s2, 'w') as outfile:
            outfile.write(printblock.format("sample", "Nreads_orig", "-qscore", 
                                            "-adapters", "Nreads_kept"))

    ## append stats to file
    outfile = open(data.statsfiles.s2, 'a+')
    outfile.write(printblock.format(sample.name, 
                                    str(fcounts["orig"]),
                                    str(fcounts["quality"]),
                                    str(fcounts["adapter"]), 
                                    str(fcounts["keep"])))
    outfile.close()

    ## always overwrite stats b/c even if multiple fastq
    ## files it concatenates them before running.
    sample.stats.reads_filtered = fcounts["keep"]

    ## save stats to Sample if successful
    if sample.stats.reads_filtered:
        sample.stats.state = 2
        sample.files.edits.append((handle1, handle2))
        sample.files.edits = list(set(sample.files.edits))

        ## save stats to the sample??
        data._stamp("s2 rawediting on "+sample.name)
    else:
        print("No reads passed filtering in Sample: {}".format(sample.name))



def run(data, samples, ipyclient, force=False, nrep=True):
    """ run the major functions for editing raw reads """
    ## print warning if skipping all samples
    if not force:
        if all([i.stats.state >= 2 for i in samples]):
            print("Skipping step2: All {} ".format(len(data.samples))\
                 +"Samples already filtered in `{}`".format(data.name))
            
        else:
            for sample in samples:
                if sample.stats.state >= 2:
                    print("Skipping Sample {}; ".format(sample.name)
                         +"Already filtered. Use force=True to overwrite.")
                elif sample.stats.reads_raw < 100:
                    print("Skipping Sample {}; ".format(sample.name)
                         +"Too few reads ({}). Use force=True to overwrite.".\
                           format(sample.stats.reads_raw))
                else:
                    submitted, results = run_full(data, sample, ipyclient, nrep)
                    cleanup(data, sample, submitted, results)
    else:
        for sample in samples:
            if not sample.stats.reads_raw:
                print("Skipping Sample {}; ".format(sample.name)
                     +"No reads found in file {}".\
                     format(sample.files.fastqs))
            else:
                submitted, results = run_full(data, sample, ipyclient, nrep)
                cleanup(data, sample, submitted, results)



if __name__ == "__main__":
    ## run test
    import ipyrad as ip

    ## test rad
    TEST = ip.load_assembly("testrad")
    TEST.step2(force=True)

    ## test gbs
    TEST = ip.load_assembly("testgbs")
    TEST.step2(force=True)

    ## test pairgbs
    TEST = ip.load_assembly("testpairgbs")
    TEST.step2(force=True)

