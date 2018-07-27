#!/usr/bin/env python2

""" demultiplex raw sequence data given a barcode map."""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0301

import os
import io
import gzip
import glob
import time
import shutil
import datetime
import itertools
import cPickle as pickle
import numpy as np
import subprocess as sps
from ipyrad.core.sample import Sample
from ipyrad.assemble.util import *
from collections import defaultdict, Counter

import logging
LOGGER = logging.getLogger(__name__)



def combinefiles(filepath):
    """ Joins first and second read file names """
    ## unpack seq files in filepath
    fastqs = glob.glob(filepath)
    firsts = [i for i in fastqs if "_R1_" in i]

    ## check names
    if not firsts:
        raise IPyradWarningExit("First read files names must contain '_R1_'.")

    ## get paired reads
    seconds = [ff.replace("_R1_", "_R2_") for ff in firsts]
    return zip(firsts, seconds)



def findbcode(cutters, longbar, read1):
    """ find barcode sequence in the beginning of read """
    ## default barcode string
    for cutter in cutters[0]:
        ## If the cutter is unambiguous there will only be one.
        if not cutter:
            continue
        search = read1[1][:int(longbar[0]+len(cutter)+1)]
        barcode = search.rsplit(cutter, 1)
        if len(barcode) > 1:
            return barcode[0]
    ## No cutter found
    return barcode[0] 



def find3radbcode(cutters, longbar, read1):
    """ find barcode sequence in the beginning of read """
    ## default barcode string
    for ambigcuts in cutters:
        for cutter in ambigcuts:
            ## If the cutter is unambiguous there will only be one.
            if not cutter:
                continue
            search = read1[1][:int(longbar[0]+len(cutter)+1)]
            splitsearch = search.rsplit(cutter, 1)
            if len(splitsearch) > 1:
                return splitsearch[0]
    ## No cutter found
    return splitsearch[0] 



def make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars):
    """
    Write stats and stores to Assembly object.
    """

    ## out file
    outhandle = os.path.join(data.dirs.fastqs, 's1_demultiplex_stats.txt')
    outfile = open(outhandle, 'w')

    ## write the header for file stats ------------------------------------
    outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                  format("raw_file", "total_reads", "cut_found", "bar_matched"))

    ## write the file stats
    r1names = sorted(perfile)
    for fname in r1names:
        dat = perfile[fname]
        #dat = [perfile[fname][i] for i in ["ftotal", "fcutfound", "fmatched"]]
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format(fname, dat[0], dat[1], dat[2]))
        ## repeat for pairfile
        if 'pair' in data.paramsdict["datatype"]:
            fname = fname.replace("_R1_", "_R2_")
            outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                format(fname, dat[0], dat[1], dat[2]))

    ## spacer, how many records for each sample --------------------------
    outfile.write('\n{:<35}  {:>13}\n'.format("sample_name", "total_reads"))

    ## names alphabetical. Write to file. Will save again below to Samples.
    snames = set()
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        snames.add(sname)
        
    for sname in sorted(list(snames)):
        outfile.write("{:<35}  {:>13}\n".format(sname, fsamplehits[sname]))

    ## spacer, which barcodes were found -----------------------------------
    outfile.write('\n{:<35}  {:>13} {:>13} {:>13}\n'.\
                  format("sample_name", "true_bar", "obs_bar", "N_records"))

    ## write sample results
    for sname in sorted(data.barcodes):
        if "-technical-replicate-" in sname:
            fname = sname.rsplit("-technical-replicate", 1)[0]  
        else:
            fname = sname
            
        ## write perfect hit
        hit = data.barcodes[sname]
        offhitstring = ""
    
        ## write off-n hits
        ## sort list of off-n hits  
        if fname in fdbars:
            offkeys = list(fdbars.get(fname))
            for offhit in offkeys[::-1]:
                ## exclude perfect hit
                if offhit not in data.barcodes.values():
                    offhitstring += '{:<35}  {:>13} {:>13} {:>13}\n'.\
                        format(sname, hit, offhit, fbarhits[offhit]/2)
                    #sumoffhits += fbarhits[offhit]
        
            ## write string to file
            outfile.write('{:<35}  {:>13} {:>13} {:>13}\n'.\
                #format(sname, hit, hit, fsamplehits[fname]-sumoffhits))
                format(sname, hit, hit, fbarhits[hit]/2))
            outfile.write(offhitstring)
        
    ## write misses
    misskeys = list(fmisses.keys())
    misskeys.sort(key=fmisses.get)
    for key in misskeys[::-1]:
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format("no_match", "_", key, fmisses[key]))
    outfile.close()        


    ## Link Sample with this data file to the Assembly object
    for sname in snames:

        ## make the sample
        sample = Sample()
        sample.name = sname

        ## allow multiple barcodes if its a replicate. 
        barcodes = []
        for n in xrange(500):
            fname = sname+"-technical-replicate-{}".format(n)
            fbar = data.barcodes.get(fname)
            if fbar:
                barcodes.append(fbar)
        if barcodes:
            sample.barcode = barcodes
        else:
            sample.barcode = data.barcodes[sname]

        ## file names        
        if 'pair' in data.paramsdict["datatype"]:
            sample.files.fastqs = [(os.path.join(data.dirs.fastqs,
                                                  sname+"_R1_.fastq.gz"),
                                     os.path.join(data.dirs.fastqs,
                                                  sname+"_R2_.fastq.gz"))]
        else:
            sample.files.fastqs = [(os.path.join(data.dirs.fastqs,
                                                  sname+"_R1_.fastq.gz"), "")]
        ## fill in the summary stats
        sample.stats["reads_raw"] = int(fsamplehits[sname])
        ## fill in the full df stats value
        sample.stats_dfs.s1["reads_raw"] = int(fsamplehits[sname])

        ## Only link Sample if it has data
        if sample.stats["reads_raw"]:
            sample.stats.state = 1
            data.samples[sample.name] = sample
        else:
            print("Excluded sample: no data found for", sname)

    ## initiate s1 key for data object
    data.stats_dfs.s1 = data._build_stat("s1")
    data.stats_files.s1 = outhandle



## EXPERIMENTAL; not yet implemented
def barmatch2(data, tups, cutters, longbar, matchdict, fnum):
    """
    cleaner barmatch func...
    """

    ## how many reads to store before writing to disk
    waitchunk = int(1e6)
    ## pid name for this engine
    epid = os.getpid()

    ## counters for total reads, those with cutsite, and those that matched
    filestat = np.zeros(3, dtype=np.int)
    ## store observed sample matches
    samplehits = {}
    ## dictionaries to store first and second reads until writing to file
    dsort1 = {} 
    dsort2 = {} 
    ## dictionary for all bars matched in sample
    dbars = {} 

    ## fill for sample names
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        samplehits[sname] = 0
        dsort1[sname] = []
        dsort2[sname] = []
        dbars[sname] = set()

    ## store observed bars
    barhits = {}
    for barc in matchdict:
        barhits[barc] = 0

    ## store others
    misses = {}
    misses['_'] = 0

    ## build func for finding barcode
    getbarcode = get_barcode_func(data, longbar)

    ## get quart iterator of reads
    if tups[0].endswith(".gz"):
        ofunc = gzip.open
    else:
        ofunc = open

    ## create iterators 
    ofile1 = ofunc(tups[0], 'r')
    fr1 = iter(ofile1) 
    quart1 = itertools.izip(fr1, fr1, fr1, fr1)
    if tups[1]:
        ofile2 = ofunc(tups[1], 'r')
        fr2 = iter(ofile2)  
        quart2 = itertools.izip(fr2, fr2, fr2, fr2)
        quarts = itertools.izip(quart1, quart2)
    else:
        quarts = itertools.izip(quart1, iter(int, 1))

    ## go until end of the file
    while 1:
        try:
            read1, read2 = quarts.next()
            read1 = list(read1)
            filestat[0] += 1
        except StopIteration:
            break
    
        barcode = ""
        ## Get barcode_R2 and check for matching sample name
        if '3rad' in data.paramsdict["datatype"]:
            ## Here we're just reusing the findbcode function
            ## for R2, and reconfiguring the longbar tuple to have the
            ## maxlen for the R2 barcode
            ## Parse barcode. Use the parsing function selected above.
            barcode1 = find3radbcode(cutters=cutters, 
                                longbar=longbar, read1=read1)
            barcode2 = find3radbcode(cutters=cutters, 
                                longbar=(longbar[2], longbar[1]), read1=read2)
            barcode = barcode1 + "+" + barcode2
        else:
            ## Parse barcode. Uses the parsing function selected above.
            barcode = getbarcode(cutters, read1, longbar)
   
        ## find if it matches 
        sname_match = matchdict.get(barcode)

        if sname_match:
            #sample_index[filestat[0]-1] = snames.index(sname_match) + 1
            ## record who matched
            dbars[sname_match].add(barcode)
            filestat[1] += 1
            filestat[2] += 1
            samplehits[sname_match] += 1
            barhits[barcode] += 1
            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1
    
            ## trim off barcode
            lenbar = len(barcode)
            if '3rad' in data.paramsdict["datatype"]:
                ## Iff 3rad trim the len of the first barcode
                lenbar = len(barcode1)
    
            if data.paramsdict["datatype"] == '2brad':
                overlen = len(cutters[0][0]) + lenbar + 1
                read1[1] = read1[1][:-overlen] + "\n"
                read1[3] = read1[3][:-overlen] + "\n"
            else:
                read1[1] = read1[1][lenbar:]
                read1[3] = read1[3][lenbar:]
    
            ## Trim barcode off R2 and append. Only 3rad datatype
            ## pays the cpu cost of splitting R2
            if '3rad' in data.paramsdict["datatype"]:
                read2 = list(read2)
                read2[1] = read2[1][len(barcode2):]
                read2[3] = read2[3][len(barcode2):]
    
            ## append to dsort
            dsort1[sname_match].append("".join(read1))
            if 'pair' in data.paramsdict["datatype"]:
                dsort2[sname_match].append("".join(read2))

        else:
            misses["_"] += 1
            if barcode:
                filestat[1] += 1

        ## how can we make it so all of the engines aren't trying to write to
        ## ~100-200 files all at the same time? This is the I/O limit we hit..
        ## write out at 100K to keep memory low. It is fine on HPC which can 
        ## write parallel, but regular systems might crash
        if not filestat[0] % waitchunk:
            ## write the remaining reads to file"
            writetofile(data, dsort1, 1, epid)
            if 'pair' in data.paramsdict["datatype"]:
                writetofile(data, dsort2, 2, epid)
            ## clear out dsorts
            for sample in data.barcodes:
                if "-technical-replicate-" in sname:
                    sname = sname.rsplit("-technical-replicate", 1)[0]
                dsort1[sname] = []
                dsort2[sname] = []
            ## reset longlist
            #longlist = np.zeros(waitchunk, dtype=np.uint32)                

    ## close open files
    ofile1.close()
    if tups[1]:
        ofile2.close()

    ## write the remaining reads to file
    writetofile(data, dsort1, 1, epid)
    if 'pair' in data.paramsdict["datatype"]:
        writetofile(data, dsort2, 2, epid)

    ## return stats in saved pickle b/c return_queue is too small
    ## and the size of the match dictionary can become quite large
    samplestats = [samplehits, barhits, misses, dbars]
    outname = os.path.join(data.dirs.fastqs, "tmp_{}_{}.p".format(epid, fnum))
    with open(outname, 'w') as wout:
        pickle.dump([filestat, samplestats], wout)

    return outname



def get_barcode_func(data, longbar):
    """ returns the fastest func given data & longbar"""
    ## build func for finding barcode
    if longbar[1] == 'same':
        if data.paramsdict["datatype"] == '2brad':
            def getbarcode(cutters, read1, longbar):
                """ find barcode for 2bRAD data """
                return read1[1][:-(len(cutters[0][0]) + 1)][-longbar[0]:]

        else:
            def getbarcode(_, read1, longbar):
                """ finds barcode for invariable length barcode data """
                return read1[1][:longbar[0]]
    else:
        def getbarcode(cutters, read1, longbar):
            """ finds barcode for variable barcode lengths"""
            return findbcode(cutters, longbar, read1)
    return getbarcode



def get_quart_iter(tups):
    """ returns an iterator to grab four lines at a time """

    if tups[0].endswith(".gz"):
        ofunc = gzip.open
    else:
        ofunc = open

    ## create iterators 
    ofile1 = ofunc(tups[0], 'r')
    fr1 = iter(ofile1) 
    quart1 = itertools.izip(fr1, fr1, fr1, fr1)
    if tups[1]:
        ofile2 = ofunc(tups[1], 'r')
        fr2 = iter(ofile2)  
        quart2 = itertools.izip(fr2, fr2, fr2, fr2)
        quarts = itertools.izip(quart1, quart2)
    else:
        ofile2 = 0
        quarts = itertools.izip(quart1, iter(int, 1))

    ## make a generator
    def feedme(quarts):
        for quart in quarts:
            yield quart
    genquarts = feedme(quarts)

    ## return generator and handles
    return genquarts, ofile1, ofile2



## called by demux2()
def barmatch(data, tups, cutters, longbar, matchdict, fnum):
    """
    Matches reads to barcodes in barcode file and writes to individual temp 
    files, after all read files have been split, temp files are collated into 
    .fastq files
    """

    ## how many reads to store before writing to disk
    waitchunk = int(1e6)

    ## pid name for this engine
    epid = os.getpid()

    ## counters for total reads, those with cutsite, and those that matched
    filestat = np.zeros(3, dtype=np.int)
    
    ## dictionary to store barcode hits for each sample
    samplehits = {}
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        samplehits[sname] = 0

    ## dict to record all barcodes
    barhits = {}
    for barc in matchdict:
        barhits[barc] = 0

    ## record everything else that is found
    misses = {}
    misses['_'] = 0

    ## dictionaries to store first and second reads until writing to file
    dsort1 = {} 
    dsort2 = {} 

    ## dictionary for all bars matched in sample
    dbars = {} 
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        dsort1[sname] = []
        dsort2[sname] = []
        dbars[sname] = set()
    
    ## get func for finding barcode
    if longbar[1] == 'same':
        if data.paramsdict["datatype"] == '2brad':
            def getbarcode(cutters, read1, longbar):
                """ find barcode for 2bRAD data """
                ## +1 is for the \n at the end of the sequence line
                lencut = len(cutters[0][0]) + 1
                return read1[1][:-lencut][-longbar[0]:]
        else:
            def getbarcode(_, read1, longbar):
                """ finds barcode for invariable length barcode data """
                return read1[1][:longbar[0]]
    else:
        def getbarcode(cutters, read1, longbar):
            """ finds barcode for variable barcode lengths"""
            return findbcode(cutters, longbar, read1)

    ## get quart iterator of reads
    if tups[0].endswith(".gz"):
        ofunc = gzip.open
    else:
        ofunc = open

    ## create iterators 
    ofile1 = ofunc(tups[0], 'r')
    fr1 = iter(ofile1) 
    quart1 = itertools.izip(fr1, fr1, fr1, fr1)
    if tups[1]:
        ofile2 = ofunc(tups[1], 'r')
        fr2 = iter(ofile2)  
        quart2 = itertools.izip(fr2, fr2, fr2, fr2)
        quarts = itertools.izip(quart1, quart2)
    else:
        quarts = itertools.izip(quart1, iter(int, 1))

    LOGGER.debug("Doing chunk %s", tups[0])

    ## go until end of the file
    while 1:
        try:
            read1, read2 = quarts.next()
            read1 = list(read1)
            filestat[0] += 1
        except StopIteration:
            break
    
        barcode = ""
        ## Get barcode_R2 and check for matching sample name
        if '3rad' in data.paramsdict["datatype"]:
            ## Here we're just reusing the findbcode function
            ## for R2, and reconfiguring the longbar tuple to have the
            ## maxlen for the R2 barcode
            ## Parse barcode. Use the parsing function selected above.
            barcode1 = find3radbcode(cutters=cutters, 
                                longbar=longbar, read1=read1)
            barcode2 = find3radbcode(cutters=cutters, 
                                longbar=(longbar[2], longbar[1]), read1=read2)
            barcode = barcode1 + "+" + barcode2
        else:
            ## Parse barcode. Uses the parsing function selected above.
            barcode = getbarcode(cutters, read1, longbar)
   
        ## find if it matches 
        sname_match = matchdict.get(barcode)

        if sname_match:
            #sample_index[filestat[0]-1] = snames.index(sname_match) + 1
            ## record who matched
            dbars[sname_match].add(barcode)
            filestat[1] += 1
            filestat[2] += 1
            samplehits[sname_match] += 1
            barhits[barcode] += 1
            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1
    
            ## trim off barcode
            lenbar = len(barcode)
            if '3rad' in data.paramsdict["datatype"]:
                ## Iff 3rad trim the len of the first barcode
                lenbar = len(barcode1)
    
            ## for 2brad we trim the barcode AND the synthetic overhang
            ## The `+1` is because it trims the newline
            if data.paramsdict["datatype"] == '2brad':
                overlen = len(cutters[0][0]) + lenbar + 1
                read1[1] = read1[1][:-overlen] + "\n"
                read1[3] = read1[3][:-overlen] + "\n"
            else:
                read1[1] = read1[1][lenbar:]
                read1[3] = read1[3][lenbar:]
    
            ## Trim barcode off R2 and append. Only 3rad datatype
            ## pays the cpu cost of splitting R2
            if '3rad' in data.paramsdict["datatype"]:
                read2 = list(read2)
                read2[1] = read2[1][len(barcode2):]
                read2[3] = read2[3][len(barcode2):]
    
            ## append to dsort
            dsort1[sname_match].append("".join(read1))
            if 'pair' in data.paramsdict["datatype"]:
                dsort2[sname_match].append("".join(read2))

        else:
            misses["_"] += 1
            if barcode:
                filestat[1] += 1

        ## how can we make it so all of the engines aren't trying to write to
        ## ~100-200 files all at the same time? This is the I/O limit we hit..
        ## write out at 100K to keep memory low. It is fine on HPC which can 
        ## write parallel, but regular systems might crash
        if not filestat[0] % waitchunk:
            ## write the remaining reads to file"
            writetofile(data, dsort1, 1, epid)
            if 'pair' in data.paramsdict["datatype"]:
                writetofile(data, dsort2, 2, epid)
            ## clear out dsorts
            for sname in data.barcodes:
                if "-technical-replicate-" in sname:
                    sname = sname.rsplit("-technical-replicate", 1)[0]
                dsort1[sname] = []
                dsort2[sname] = []             

    ## close open files
    ofile1.close()
    if tups[1]:
        ofile2.close()

    ## write the remaining reads to file
    writetofile(data, dsort1, 1, epid)
    if 'pair' in data.paramsdict["datatype"]:
        writetofile(data, dsort2, 2, epid)

    ## return stats in saved pickle b/c return_queue is too small
    ## and the size of the match dictionary can become quite large
    samplestats = [samplehits, barhits, misses, dbars]
    outname = os.path.join(data.dirs.fastqs, "tmp_{}_{}.p".format(epid, fnum))
    with open(outname, 'w') as wout:
        pickle.dump([filestat, samplestats], wout)

    return outname



def writetofastq(data, dsort, read):
    """ 
    Writes sorted data 'dsort dict' to a tmp files
    """
    if read == 1:
        rrr = "R1"
    else:
        rrr = "R2"

    for sname in dsort:
        ## skip writing if empty. Write to tmpname
        handle = os.path.join(data.dirs.fastqs, 
                "{}_{}_.fastq".format(sname, rrr))
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))



def writetofile(data, dsort, read, pid):
    """ 
    Writes sorted data 'dsort dict' to a tmp files
    """
    if read == 1:
        rrr = "R1"
    else:
        rrr = "R2"

    for sname in dsort:
        ## skip writing if empty. Write to tmpname
        handle = os.path.join(data.dirs.fastqs, 
                "tmp_{}_{}_{}.fastq".format(sname, rrr, pid))
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))



def new_collate_files(data, sname, tmp1s, tmp2s):
    
    ## outfile handle (before gzip)
    out1 = os.path.join(data.dirs.fastqs, "{}_R1_.fastq".format(sname))

    ## build command
    cmd1 = ["cat"]
    for tmpfile in tmp1s:
        cmd1 += tmpfile

    ## run concat command
    proc = sps.Popen(cmd1, stderr=sps.PIPE, stdout=out)
    proc.communicate()

    ## run compression on the file
    ## first look to see if pigz is available
    proc = sps.Popen(['which', 'pigz'], stderr=sps.PIPE, stdout=sps.PIPE).communicate()
    if proc[0].strip():
        compress = ["pigz", "-p", "4"]
    else:
        compress = ["gzip"]



def collate_files(data, sname, tmp1s, tmp2s):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    ## out handle
    out1 = os.path.join(data.dirs.fastqs, "{}_R1_.fastq.gz".format(sname))
    out = io.BufferedWriter(gzip.open(out1, 'w'))

    ## build cmd
    cmd1 = ['cat']
    for tmpfile in tmp1s:
        cmd1 += [tmpfile]

    ## compression function
    proc = sps.Popen(['which', 'pigz'], stderr=sps.PIPE, stdout=sps.PIPE).communicate()
    if proc[0].strip():
        compress = ["pigz"]
    else:
        compress = ["gzip"]

    ## call cmd
    proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
    proc2 = sps.Popen(compress, stdin=proc1.stdout, stderr=sps.PIPE, stdout=out)
    err = proc2.communicate()
    if proc2.returncode:
        raise IPyradWarningExit("error in collate_files R1 %s", err)
    proc1.stdout.close()
    out.close()

    ## then cleanup
    for tmpfile in tmp1s:
        os.remove(tmpfile)

    if 'pair' in data.paramsdict["datatype"]:
        ## out handle
        out2 = os.path.join(data.dirs.fastqs, "{}_R2_.fastq.gz".format(sname))
        out = io.BufferedWriter(gzip.open(out2, 'w'))

        ## build cmd
        cmd1 = ['cat']
        for tmpfile in tmp2s:
            cmd1 += [tmpfile]

        ## call cmd
        proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
        proc2 = sps.Popen(compress, stdin=proc1.stdout, stderr=sps.PIPE, stdout=out)
        err = proc2.communicate()
        if proc2.returncode:
            raise IPyradWarningExit("error in collate_files R2 %s", err)
        proc1.stdout.close()
        out.close()

        ## then cleanup
        for tmpfile in tmp2s:
            os.remove(tmpfile)


def prechecks2(data, force):
    """
    A new simplified version of prechecks func before demux
    Checks before starting analysis. 
    -----------------------------------
    1) Is there data in raw_fastq_path
    2) Is there a barcode file
    3) Is there a workdir and fastqdir
    4) remove old fastq/tmp_sample_R*_ dirs/
    5) return file names as pairs (r1, r2) or fakepairs (r1, 1)
    6) get ambiguous cutter resolutions
    7) get optim size
    """

    ## check for data using glob for fuzzy matching
    if not glob.glob(data.paramsdict["raw_fastq_path"]):
        raise IPyradWarningExit(NO_RAWS.format(data.paramsdict["raw_fastq_path"]))

    ## find longest barcode
    try:
        ## Handle 3rad multi-barcodes. Gets len of the first one. 
        ## Should be harmless for single barcode data
        barlens = [len(i.split("+")[0]) for i in data.barcodes.values()]
        if len(set(barlens)) == 1:
            longbar = (barlens[0], 'same')
        else:
            longbar = (max(barlens), 'diff')

        ## For 3rad we need to add the length info for barcodes_R2
        if "3rad" in data.paramsdict["datatype"]:
            barlens = [len(i.split("+")[1]) for i in data.barcodes.values()]
            longbar = (longbar[0], longbar[1], max(barlens))
    except ValueError:
        raise IPyradWarningExit(NO_BARS.format(data.paramsdict["barcodes_path"]))

    ## setup dirs: [workdir] and a [workdir/name_fastqs]
    opj = os.path.join

    ## create project dir
    pdir = os.path.realpath(data.paramsdict["project_dir"])
    if not os.path.exists(pdir):
        os.mkdir(pdir)

    ## create fastq dir
    data.dirs.fastqs = opj(pdir, data.name+"_fastqs")
    if os.path.exists(data.dirs.fastqs) and force:
        print(OVERWRITING_FASTQS.format(**{"spacer":data._spacer}))
        shutil.rmtree(data.dirs.fastqs)
    if not os.path.exists(data.dirs.fastqs):
        os.mkdir(data.dirs.fastqs)

    ## insure no leftover tmp files from a previous run (there shouldn't be)
    oldtmps = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R1_"))
    oldtmps += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R2_"))    
    for oldtmp in oldtmps:
        os.remove(oldtmp)

    ## gather raw sequence filenames (people want this to be flexible ...)
    if 'pair' in data.paramsdict["datatype"]:
        raws = combinefiles(data.paramsdict["raw_fastq_path"])
    else:
        raws = zip(glob.glob(data.paramsdict["raw_fastq_path"]), iter(int, 1))

    ## returns a list of both resolutions of cut site 1
    ## (TGCAG, ) ==> [TGCAG, ]
    ## (TWGC, ) ==> [TAGC, TTGC]
    ## (TWGC, AATT) ==> [TAGC, TTGC]
    cutters = [ambigcutters(i) for i in data.paramsdict["restriction_overhang"]]
    assert cutters, "Must enter a `restriction_overhang` for demultiplexing."

    ## get matchdict
    matchdict = inverse_barcodes(data)

    ## return all
    return raws, longbar, cutters, matchdict



def inverse_barcodes(data):
    """ Build full inverse barcodes dictionary """

    matchdict = {}
    bases = set("CATGN")
    poss = set()

    ## do perfect matches
    for sname, barc in data.barcodes.items():
        ## remove -technical-replicate-N if present
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        matchdict[barc] = sname
        poss.add(barc)

        if data.paramsdict["max_barcode_mismatch"] > 0:
            ## get 1-base diffs
            for idx1, base in enumerate(barc):
                diffs = bases.difference(base)
                for diff in diffs:
                    lbar = list(barc)
                    lbar[idx1] = diff
                    tbar1 = "".join(lbar)
                    if tbar1 not in poss:
                        matchdict[tbar1] = sname                    
                        poss.add(tbar1)
                    else:
                        if matchdict.get(tbar1) != sname:
                            print("""\
        Note: barcodes {}:{} and {}:{} are within {} base change of each other
            Ambiguous barcodes that match to both samples will arbitrarily
            be assigned to the first sample. If you do not like this idea 
            then lower the value of max_barcode_mismatch and rerun step 1\n"""\
        .format(sname, barc, 
                matchdict[tbar1], data.barcodes[matchdict[tbar1]],
                data.paramsdict["max_barcode_mismatch"]))

                ## if allowing two base difference things get big
                ## for each modified bar, allow one modification to other bases
                if data.paramsdict["max_barcode_mismatch"] > 1:
                    for idx2, _ in enumerate(tbar1):
                        ## skip the base that is already modified
                        if idx2 != idx1:
                            for diff in bases.difference(tbar1[idx2]):
                                ltbar = list(tbar1)
                                ltbar[idx2] = diff
                                tbar2 = "".join(ltbar)
                                if tbar2 not in poss:
                                    matchdict[tbar2] = sname                    
                                    poss.add(tbar2)
                                else:
                                    if matchdict.get(tbar2) != sname:
                                        print("""\
        Note: barcodes {}:{} and {}:{} are within {} base change of each other\
             Ambiguous barcodes that match to both samples will arbitrarily
             be assigned to the first sample. If you do not like this idea 
             then lower the value of max_barcode_mismatch and rerun step 1\n"""\
        .format(sname, barc, 
                            matchdict[tbar2], data.barcodes[matchdict[tbar2]],
                            data.paramsdict["max_barcode_mismatch"]))
    return matchdict


def estimate_optim(data, testfile, ipyclient):
    """ 
    Estimate a reasonable optim value by grabbing a chunk of sequences, 
    decompressing and counting them, to estimate the full file size.
    """
    ## count the len of one file and assume all others are similar len
    insize = os.path.getsize(testfile)
    tmp_file_name = os.path.join(data.paramsdict["project_dir"], "tmp-step1-count.fq")
    if testfile.endswith(".gz"):
        infile = gzip.open(testfile)
        outfile = gzip.open(tmp_file_name, 'wb', compresslevel=5)
    else:
        infile = open(testfile)
        outfile = open(tmp_file_name, 'w')
        
    ## We'll take the average of the size of a file based on the
    ## first 10000 reads to approximate number of reads in the main file
    outfile.write("".join(itertools.islice(infile, 40000)))
    outfile.close()
    infile.close()

    ## Get the size of the tmp file
    tmp_size = os.path.getsize(tmp_file_name)

    ## divide by the tmp file size and multiply by 10000 to approximate
    ## the size of the input .fq files
    inputreads = int(insize / tmp_size) * 10000
    os.remove(tmp_file_name)

    return inputreads


def run2(data, ipyclient, force):
    """
    One input file (or pair) is run on two processors, one for reading 
    and decompressing the data, and the other for demuxing it.
    """

    ## get file handles, name-lens, cutters, and matchdict
    raws, longbar, cutters, matchdict = prechecks2(data, force)

    ## wrap funcs to ensure we can kill tmpfiles
    kbd = 0
    try:
        ## if splitting files, split files into smaller chunks for demuxing
        chunkfiles = splitfiles(data, raws, ipyclient)

        ## send chunks to be demux'd
        statdicts = demux2(data, chunkfiles, cutters, longbar, matchdict, ipyclient)

        ## concat tmp files
        concat_chunks(data, ipyclient)

        ## build stats from dictionaries
        perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts    
        make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)


    except KeyboardInterrupt:
        print("\n  ...interrupted, just a second while we ensure proper cleanup")
        kbd = 1

    ## cleanup
    finally:
        ## cleaning up the tmpdir is safe from ipyclient
        tmpdir = os.path.join(data.paramsdict["project_dir"], "tmp-chunks-"+data.name)
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)

        if kbd:
            raise KeyboardInterrupt("s1")
        else:
            _cleanup_and_die(data)


def _cleanup_and_die(data):
    """ cleanup func for step 1 """
    tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*.fastq"))
    tmpfiles += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.p"))
    for tmpf in tmpfiles:            
        os.remove(tmpf)


## EXPERIMENTAL; not yet implemented. Tries to skip chunking big files. Maybe 
## faster for some cases where most time is spent chunking.
def run3(data, ipyclient, force):
    """
    One input file (or pair) is run on two processors, one for reading 
    and decompressing the data, and the other for demuxing it.
    """

    start = time.time()
    ## get file handles, name-lens, cutters, and matchdict, 
    ## and remove any existing files if a previous run failed.
    raws, longbar, cutters, matchdict = prechecks2(data, force)

    ## wrap funcs to ensure we can kill tmpfiles
    kbd = 0
    try:
        ## send chunks to be demux'd, nothing is parallelized yet.
        lbview = ipyclient.load_balanced_view()
        args = (data, raws, cutters, longbar, matchdict)
        async = lbview.apply(demux3, *args)

        ## track progress
        while 1:
            ## how many of this func have finished so far
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            printstr = ' writing/compressing   | {} | s1 |'
            progressbar(len(ready), sum(ready), printstr, spacer=spacer)
            time.sleep(0.1)
            if async.ready():
                print("")
                break
        
        if async.successful():
            statdicts = async.get()
        else:
            raise IPyradWarningExit(async.get())

        ## build stats from dictionaries
        perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts
        make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)


    except KeyboardInterrupt:
        print("\n  ...interrupted, just a second while we ensure proper cleanup")
        kbd = 1

    ## cleanup
    finally:
        ## cleaning up the tmpdir is safe from ipyclient
        tmpdir = os.path.join(data.paramsdict["project_dir"], "tmp-chunks-"+data.name)
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)

        tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*.fastq"))
        tmpfiles += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.p"))
        for tmpf in tmpfiles:
            if os.path.exists(tmpf):
                os.remove(tmpf)

        if kbd:
            raise 



def splitfiles(data, raws, ipyclient):
    """ sends raws to be chunked"""

    ## create a tmpdir for chunked_files and a chunk optimizer 
    tmpdir = os.path.join(data.paramsdict["project_dir"], "tmp-chunks-"+data.name)
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    ## chunk into 8M reads
    totalreads = estimate_optim(data, raws[0][0], ipyclient)
    optim = int(8e6)
    njobs = int(totalreads/(optim/4.)) * len(raws)

    ## if more files than cpus: no chunking
    nosplit = 0
    if (len(raws) > len(ipyclient)) or (totalreads < optim):
        nosplit = 1

    ## send slices N at a time. The dict chunkfiles stores a tuple of rawpairs
    ## dictionary to store asyncresults for sorting jobs
    start = time.time()
    chunkfiles = {}
    for fidx, tups in enumerate(raws):
        handle = os.path.splitext(os.path.basename(tups[0]))[0]
        ## if number of lines is > 20M then just submit it
        if nosplit:
            chunkfiles[handle] = [tups]
        else:
            ## chunk the file using zcat_make_temps
            chunklist = zcat_make_temps(data, tups, fidx, tmpdir, optim, njobs, start)
            chunkfiles[handle] = chunklist
    if not nosplit:
        print("")

    return chunkfiles



def concat_chunks(data, ipyclient):
    """ 
    Concatenate chunks. If multiple chunk files match to the same sample name
    but with different barcodes (i.e., they are technical replicates) then this
    will assign all the files to the same sample name file.
    """

    ## collate files progress bar
    start = time.time()
    printstr = ' writing/compressing   | {} | s1 |'
    lbview = ipyclient.load_balanced_view()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(10, 0, printstr.format(elapsed), spacer=data._spacer) 
    ## get all the files
    ftmps = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.fastq"))

    ## a dict to assign tmp files to names/reads
    r1dict = {}
    r2dict = {}
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        r1dict[sname] = []
        r2dict[sname] = []

    ## assign to name keys
    for ftmp in ftmps:
        base, orient, _ = ftmp.rsplit("_", 2)
        sname = base.rsplit("/", 1)[-1].split("tmp_", 1)[1]
        if orient == "R1":
            r1dict[sname].append(ftmp)
        else:
            r2dict[sname].append(ftmp)

    ## concatenate files
    snames = []
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        snames.append(sname)

    writers = []
    for sname in set(snames):
        tmp1s = sorted(r1dict[sname])
        tmp2s = sorted(r2dict[sname])
        writers.append(lbview.apply(collate_files, *[data, sname, tmp1s, tmp2s]))

    total = len(writers)
    while 1:
        ready = [i.ready() for i in writers]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(total, sum(ready), printstr.format(elapsed), spacer=data._spacer) 
        time.sleep(0.1)
        if all(ready):
            print("")
            break



## EXPERIMENTAL; NOT YET IMPLEMENTED
def demux3(data, raws, cutters, longbar, matchdict):

    ## store counters
    chunksize = int(1e6)
    perfile = {}
    #start = time.time()

    ## choose func for finding barcode
    if longbar[1] == 'same':
        if data.paramsdict["datatype"] == '2brad':
            def getbarcode(_, read1, longbar):
                """ find barcode for 2bRAD data """
                ## +1 is for the \n at the end of the sequence line
                lencut = len(cutters[0][0]) + 1
                return read1[1][:-lencut][-longbar[0]:]
        else:
            def getbarcode(_, read1, longbar):
                """ finds barcode for invariable length barcode data """
                return read1[1][:longbar[0]]
    else:
        def getbarcode(cutters, read1, longbar):
            """ finds barcode for variable barcode lengths"""
            return findbcode(cutters, longbar, read1)    

    ## start streaming data in 
    for rawtup in raws:

        ## get the read generator
        handle = rawtup[0].rsplit("/", 1)[1]
        filestat = np.zeros(3)

        ## this returns a generator (gen) to grab pairs of sequences 
        ## 4-lines at a time, and it returns the two open file objects
        ## so we can close them later.
        gen, handle1, handle2 = get_quart_iter(rawtup)

        ## dictionaries to store hits samples and barcodes, and misses.
        samplehits = {}
        barhits = {}
        misses = {"_": 0}
        dsort1 = {} 
        dsort2 = {} 
        dbars = {} 
        for sname in data.barcodes:
            samplehits[sname] = 0
            dsort1[sname] = []
            dsort2[sname] = []
            dbars[sname] = set()
        for barc in matchdict:
            barhits[barc] = 0

        ## iterate over the file 
        for read1, read2 in gen:
            filestat[0] += 1

            ## Parse barcode. Use the parsing function selected above.
            if '3rad' in data.paramsdict["datatype"]:
                barcode1 = find3radbcode(cutters=cutters, longbar=longbar, read1=read1)
                barcode2 = find3radbcode(cutters=cutters, longbar=(longbar[2], longbar[1]), read1=read2)
                barcode = barcode1 + "+" + barcode2
            else:
                barcode = getbarcode(cutters, read1, longbar)
            sname_match = matchdict.get(barcode)

            if sname_match:
                ## need to be lists to modify them
                read1 = list(read1)
                if read2:
                    read2 = list(read2)

                ## record the match
                filestat[1] += 1
                filestat[2] += 1
                dbars[sname_match].add(barcode)
                samplehits[sname_match] += 1
                barhits[barcode] += 1
                if barhits.get(barcode):
                    barhits[barcode] += 1
                else:
                    barhits[barcode] = 1

                ## trim off barcode
                lenbar = len(barcode)
                if data.paramsdict["datatype"] == "3rad":
                    lenbar = len(barcode1)
                if data.paramsdict["datatype"] == "2brad":
                    overlen = len(cutters[0][0]) + lenbar + 1
                    read1[1] = read1[1][:-overlen] + "\n"
                    read1[3] = read1[3][:-overlen] + "\n"
                else:
                    read1[1] = read1[1][lenbar:]
                    read1[3] = read1[3][lenbar:]

                ## trim barcode from R2 and more for 3rad
                if data.paramsdict["datatype"] == "3rad":
                    read2 = list(read2)
                    read2[1] = read2[1][len(barcode2):]
                    read2[3] = read2[3][len(barcode2):]

                ## append to dsort
                dsort1[sname_match].append("".join(read1))
                if "pair" in data.paramsdict["datatype"]:
                    dsort2[sname_match].append("".join(read2))

            else:
                ## did not match to a sample
                misses["_"] += 1
                ## but store that a barcode was observed
                if barcode:
                    filestat[1] += 1

            ## write to file and clear dictionaries
            if not filestat[0] % chunksize:
                ## dump reads to file
                writetofastq(data, dsort1, 1)
                if 'pair' in data.paramsdict["datatype"]:
                    writetofastq(data, dsort2, 2)
                ## clear out dsorts
                for sample in data.barcodes:
                    dsort1[sample] = []
                    dsort2[sample] = []

        ## close open files
        handle1.close()
        if rawtup[1]:
            handle2.close()
        perfile[handle] = filestat

        ## write the remaining reads to file
        writetofastq(data, dsort1, 1)
        if 'pair' in data.paramsdict["datatype"]:
            writetofastq(data, dsort2, 2)

    ## return stats in saved pickle b/c return_queue is too small
    ## and the size of the match dictionary can become quite large
    #samplestats = [samplehits, barhits, misses, dbars]
    #outname = os.path.join(data.dirs.fastqs, "tmp_{}_{}.p".format(epid, fnum))
    #with open(outname, 'w') as wout:
    #    pickle.dump([filestat, samplestats], wout)
    #return outname
    return perfile, samplehits, barhits, misses, dbars
                     


def demux2(data, chunkfiles, cutters, longbar, matchdict, ipyclient):
    """ 
    Submit chunks to be sorted by the barmatch() function then 
    calls putstats().
    """

    ## parallel stuff, limit to 1/4 of available cores for RAM limits.
    start = time.time()
    printstr = ' sorting reads         | {} | s1 |'
    lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::4])

    ## store statcounters and async results in dicts
    perfile = {}
    filesort = {}
    total = 0
    done = 0 

    ## chunkfiles is a dict with {handle: chunkslist, ...}. The func barmatch
    ## writes results to samplename files with PID number, and also writes a 
    ## pickle for chunk specific results with fidx suffix, which it returns.
    for handle, rawtuplist in chunkfiles.items():
        ## get args for job
        for fidx, rawtuple in enumerate(rawtuplist):
            #handle = os.path.splitext(os.path.basename(rawtuple[0]))[0]
            args = (data, rawtuple, cutters, longbar, matchdict, fidx)

            ## submit the job
            async = lbview.apply(barmatch, *args)
            filesort[total] = (handle, async)
            total += 1

            ## get ready to receive stats: 'total', 'cutfound', 'matched'
            perfile[handle] = np.zeros(3, dtype=np.int)

    ## stats for each sample
    fdbars = {}
    fsamplehits = Counter()
    fbarhits = Counter()
    fmisses = Counter()
    ## a tuple to hold my dictionaries
    statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars

    ## wait for jobs to finish
    while 1:
        fin = [i for i, j in filesort.items() if j[1].ready()]
        #fin = [i for i in jobs if i[1].ready()]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(total, done, printstr.format(elapsed), spacer=data._spacer)
        time.sleep(0.1)

        ## should we break?
        if total == done:
            print("")
            break

        ## cleanup
        for key in fin:
            tup = filesort[key]
            if tup[1].successful():
                pfile = tup[1].result()
                handle = tup[0]
                if pfile:
                    ## check if this needs to return data
                    putstats(pfile, handle, statdicts)
                    ## purge to conserve memory
                    del filesort[key]
                    done += 1

    return statdicts



## DEPRECATED FOR DEMUX2()
def demux(data, chunkfiles, cutters, longbar, matchdict, ipyclient):
    """ submit chunks to be sorted """

    ## parallel stuff
    start = time.time()
    printstr = ' sorting reads         | {} | s1 |'
    lbview = ipyclient.load_balanced_view()

    ## store statcounters and async results in dicts
    perfile = {}
    filesort = {}
    for handle, rawtuplist in chunkfiles.items():
        ## get args for job
        for fidx, rawtuple in enumerate(rawtuplist):
            #handle = os.path.splitext(os.path.basename(rawtuple[0]))[0]
            args = (data, rawtuple, cutters, longbar, matchdict, fidx)

            ## submit the job
            filesort[handle] = lbview.apply(barmatch, *args)

            ## get ready to receive stats: 'total', 'cutfound', 'matched'
            perfile[handle] = np.zeros(3, dtype=np.int)

    ## stats for each sample
    fdbars = {}
    fsamplehits = Counter()
    fbarhits = Counter()
    fmisses = Counter()
    ## a tuple to hold my dictionaries
    statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars

    try:
        kbd = 0
        total = len(chunkfiles)
        done = 0
        ## wait for jobs to finish
        while 1:
            fin = [i for i, j in filesort.items() if j.ready()]
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(total, done, printstr.format(elapsed), spacer=data._spacer)
            time.sleep(0.1)

            ## should we break?
            if total == done:
                print("")
                break

            ## cleanup
            for job in fin:
                if filesort[job].successful():
                    pfile = filesort[job].result()
                    #if result:
                    if pfile:
                        ## check if this needs to return data
                        putstats(pfile, handle, statdicts)
                        
                        ## purge to conserve memory
                        del filesort[job]
                        done += 1

        ## keep tacking progreess during writing stage
        start = time.time()
        printstr = ' writing/compressing   | {} | s1 |'
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(10, 0, printstr.format(elapsed), spacer=data._spacer)


    except KeyboardInterrupt:
        ## wait to cleanup
        kbd = 1
        raise


    ## only proceed here if barmatch jobs were not interrupted
    else:
        ## collate files and do progress bar
        ftmps = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.fastq"))

        ## a dict to assign tmp files to names/reads
        r1dict = {}
        r2dict = {}
        for sname in data.barcodes:
            r1dict[sname] = []
            r2dict[sname] = []

        ## assign to name keys
        for ftmp in ftmps:
            ## split names
            base, orient, _ = ftmp.rsplit("_", 2)
            sname = base.rsplit("/", 1)[-1].split("tmp_", 1)[1]
            ## put into dicts
            if orient == "R1":
                r1dict[sname].append(ftmp)
            else:
                r2dict[sname].append(ftmp)

        ## concatenate files
        total = len(data.barcodes)
        done = 0

        ## store asyncs of collate jobs
        writers = []
        for sname in data.barcodes:
            tmp1s = sorted(r1dict[sname])
            tmp2s = sorted(r2dict[sname])
            writers.append(lbview.apply(collate_files, 
                           *[data, sname, tmp1s, tmp2s]))
        
        ## track progress of collate jobs
        while 1:
            ready = [i.ready() for i in writers]
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(total, sum(ready), printstr.format(elapsed), spacer=data._spacer)  
            time.sleep(0.1)
            if all(ready):
                print("")
                break

    finally:
        ## clean up junk files
        tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*.fastq"))
        tmpfiles += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.p"))
        for tmpf in tmpfiles:
            os.remove(tmpf)

        if kbd:
            raise KeyboardInterrupt()
        else:
            ## build stats from dictionaries
            perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts    
            make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)


def putstats(pfile, handle, statdicts):
    """ puts stats from pickles into a dictionary """

    ## load in stats
    with open(pfile, 'r') as infile:
        filestats, samplestats = pickle.load(infile)

    ## get dicts from statdicts tuple
    perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts

    ## pull new stats
    #handle = os.path.splitext(os.path.basename(handle))[0]
    perfile[handle] += filestats

    ## update sample stats
    samplehits, barhits, misses, dbars = samplestats
    fsamplehits.update(samplehits)
    fbarhits.update(barhits)
    fmisses.update(misses)
    fdbars.update(dbars)

    ## repack the tuple and return
    statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars
    return statdicts


## used by splitfiles()
def zcat_make_temps(data, raws, num, tmpdir, optim, njobs, start):
    """ 
    Call bash command 'cat' and 'split' to split large files. The goal
    is to create N splitfiles where N is a multiple of the number of processors
    so that each processor can work on a file in parallel.
    """

    printstr = ' chunking large files  | {} | s1 |'

    ## split args
    tmpdir = os.path.realpath(tmpdir)
    LOGGER.info("zcat is using optim = %s", optim)

    ## read it, is it gzipped?
    catcmd = ["cat"]
    if raws[0].endswith(".gz"):
        catcmd = ["gunzip", "-c"]

    ## get reading commands for r1s, r2s
    cmd1 = catcmd + [raws[0]]
    cmd2 = catcmd + [raws[1]]

    ## second command splits and writes with name prefix
    cmd3 = ["split", "-a", "4", "-l", str(int(optim)), "-", 
            os.path.join(tmpdir, "chunk1_"+str(num)+"_")]
    cmd4 = ["split", "-a", "4", "-l", str(int(optim)), "-", 
            os.path.join(tmpdir, "chunk2_"+str(num)+"_")]

    ### run splitter
    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc3 = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc1.stdout)

    ## wrap the actual call so we can kill it if anything goes awry
    while 1:
        try:
            if not isinstance(proc3.poll(), int):
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                done = len(glob.glob(os.path.join(tmpdir, 'chunk1_*')))
                progressbar(njobs, min(njobs, done), printstr.format(elapsed), spacer=data._spacer)
                time.sleep(0.1)
            else:
                res = proc3.communicate()[0]
                proc1.stdout.close()
                break

        except KeyboardInterrupt:
            proc1.kill()
            proc3.kill()
            raise KeyboardInterrupt()

    if proc3.returncode:
        raise IPyradWarningExit(" error in %s: %s", cmd3, res)

    ## grab output handles
    chunks1 = glob.glob(os.path.join(tmpdir, "chunk1_"+str(num)+"_*"))
    chunks1.sort()

    if "pair" in data.paramsdict["datatype"]:
        proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc4 = sps.Popen(cmd4, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc2.stdout)

        ## wrap the actual call so we can kill it if anything goes awry
        while 1:
            try:
                if not isinstance(proc4.poll(), int):
                    elapsed = datetime.timedelta(seconds=int(time.time()-start))
                    done = len(glob.glob(os.path.join(tmpdir, 'chunk1_*')))
                    progressbar(njobs, min(njobs, done), printstr.format(elapsed), data._spacer)
                    time.sleep(0.1)
                else:
                    res = proc4.communicate()[0]
                    proc2.stdout.close()
                    break

            except KeyboardInterrupt:
                proc2.kill()
                proc4.kill()
                raise KeyboardInterrupt()

        if proc4.returncode:
            raise IPyradWarningExit(" error in %s: %s", cmd4, res)

        ## grab output handles
        chunks2 = glob.glob(os.path.join(tmpdir, "chunk2_"+str(num)+"_*"))
        chunks2.sort()
    
    else:
        chunks2 = [0]*len(chunks1)

    assert len(chunks1) == len(chunks2), \
        "R1 and R2 files are not the same length."

    ## ensure full progress bar b/c estimates njobs could be off
    progressbar(10, 10, printstr.format(elapsed), spacer=data._spacer)
    return zip(chunks1, chunks2)



## GLOBALS
NO_BARS = """\
    Barcodes file not found. You entered: '{}'
    """

NO_RAWS = """\
    No data found in {}. Fix path to data files.
    """

OVERWRITING_FASTQS = """\
{spacer}[force] overwriting fastq files previously created by ipyrad.
{spacer}This _does not_ affect your original/raw data files."""


if __name__ == "__main__":

    ## run test
    import ipyrad as ip

    ## get current location
    #PATH = os.path.abspath(os.path.dirname(__file__))
    ## get location of test files
    #IPATH = os.path.dirname(os.path.dirname(PATH))
    #DATA = os.path.join(IPATH, "tests", "test_rad")

    ## create an Assembly
    TEST = ip.Assembly("test")

    ## this expects that you have an ipcluster running...
    TEST.set_params(1, "./maintest")
    TEST.set_params(2, "./ipsimdata/rad_example_R1_.fastq.gz")
    TEST.set_params(3, "./ipsimdata/rad_example_barcodes.txt")

    ## run demultiplexing
    TEST.run(1)
