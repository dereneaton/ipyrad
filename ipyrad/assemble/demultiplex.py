#!/usr/bin/env ipython2

""" demultiplex raw sequence data given a barcode map."""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212

import os
import gzip
import glob
import time
import shutil
import datetime
import itertools
import numpy as np
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



def findbcode(cut, longbar, read1):
    """ find barcode sequence in the beginning of read """
    ## default barcode string
    search = read1[1][:int(longbar[0]+len(cut)+1)]
    barcode = search.rsplit(cut, 1)[0]
    return barcode 



def make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars):
    """
    Write stats and stores to Assembly object.
    """

    ## out file
    outhandle = os.path.join(data.dirs.fastqs, 's1_demultiplex_stats.txt')
    outfile = open(outhandle, 'w')

    ## how many from each rawfile
    outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                  format("raw_file", "total_reads", "cut_found", "bar_matched"))

    ## sort rawfile names
    rawfilenames = sorted(perfile)
    for rawstat in rawfilenames:
        dat = [perfile[rawstat][i] for i in ["ftotal", "fcutfound", "fmatched"]]
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format(*[rawstat]+[str(i) for i in dat]))
        if "pair" in data.paramsdict["datatype"]:
            rawstat2 = rawstat.replace("_R1_", "_R2_")
            outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                format(*[rawstat2]+[str(i) for i in dat]))

    ## spacer, how many records for each sample
    outfile.write('\n{:<35}  {:>13}\n'.\
                  format("sample_name", "total_reads"))

    ## names alphabetical. Write to file. Will save again below to Samples.
    names_sorted = sorted(data.barcodes)
    for name in names_sorted:
        outfile.write("{:<35}  {:>13}\n".format(name, fsamplehits[name]))

    ## spacer, which barcodes were found
    outfile.write('\n{:<35}  {:>13}{:>13}{:>13}\n'.\
                  format("sample_name", "true_bar", "obs_bar", "N_records"))

    ## write sample results
    for name in names_sorted:
        ## write perfect hit
        hit = data.barcodes[name]
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format(name, hit, hit, fsamplehits[name]))

        ## write off-n hits
        ## sort list of off-n hits
        if name in fdbars:
            offkeys = list(fdbars.get(name))
            offkeys.sort(key=fbarhits.get)
            for offhit in offkeys[::-1]:
                ## exclude perfect hit
                if offhit not in data.barcodes.values():
                    outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                        format(name, hit, offhit, fbarhits[offhit]))

    ## write misses
    misskeys = list(fmisses.keys())
    misskeys.sort(key=fmisses.get)
    for key in misskeys[::-1]:
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format("no_match", "_", key, fmisses[key]))
    outfile.close()

    ## Link Sample with this data file to the Assembly object
    for name in data.barcodes:
        sample = Sample()
        sample.name = name
        sample.barcode = data.barcodes[name]
        if "pair" in data.paramsdict["datatype"]:
            sample.files.fastqs = [(os.path.join(data.dirs.fastqs,
                                                  name+"_R1_.fastq.gz"),
                                     os.path.join(data.dirs.fastqs,
                                                  name+"_R2_.fastq.gz"))]
        else:
            sample.files.fastqs = [(os.path.join(data.dirs.fastqs,
                                                  name+"_R1_.fastq.gz"), "")]
        ## fill in the summary stats
        sample.stats["reads_raw"] = int(fsamplehits[name])
        ## fill in the full df stats value
        sample.stats_dfs.s1["reads_raw"] = int(fsamplehits[name])

        ## Only link Sample if it has data
        if sample.stats["reads_raw"]:
            sample.stats.state = 1
            data.samples[sample.name] = sample
        else:
            print("Excluded sample: no data found for", name)

    ## initiate s1 key for data object
    data.stats_dfs.s1 = data.build_stat("s1")
    data.stats_files.s1 = outhandle



def barmatch(args):
    """
    Matches reads to barcodes in barcode file and writes to individual temp 
    files, after all read files have been split, temp files are collated into 
    .fastq files
    """

    ## read in list of args
    data, filenum, subnum, cutters, longbar, matchdict, fr1, fr2 = args
    #LOGGER.info("Entering barmatch - {}.{}".format(filenum, subnum))

    ## counters for total reads, those with cutsite, and those that matched
    total = 0
    cutfound = 0
    matched = 0 

    ## dictionary to sample hits
    samplehits = {}
    for sname in data.barcodes:
        samplehits[sname] = 0
    ## dict to record matchable barcodes that hit
    barhits = {}
    for barc in matchdict:
        barhits[barc] = 0
    ## record everything else that is found
    misses = {}
    misses['_'] = 0
    
    ## get slice of data as a generator
    # tups = rawtuple
    # fr1, fr2, io1, io2 = get_slice(tups, optim, subnum)

    fr1 = iter(fr1)
    fr2 = iter(fr2)
    # ## create iterators to sample 4 lines at a time
    quart1 = itertools.izip(fr1, fr1, fr1, fr1)
    if 'pair' in data.paramsdict["datatype"]:
        quart2 = itertools.izip(fr2, fr2, fr2, fr2)
        quarts = itertools.izip(quart1, quart2)
    else:
        quarts = itertools.izip(quart1, iter(int, 1))


    ## dictionaries to store first and second reads
    dsort1 = {} 
    dsort2 = {} 
    ## dictionary for all bars matched in sample
    dbars = {} 
    for sample in data.barcodes:
        dsort1[sample] = []
        dsort2[sample] = []
        dbars[sample] = set()

    ## get func for finding barcode
    if longbar[1] == 'same':
        if data.paramsdict["datatype"] == '2brad':
            def getbarcode(_, read1, longbar):
                """ find barcode for 2bRAD data """
                return read1[1][-longbar[0]:]
        else:
            def getbarcode(_, read1, longbar):
                """ finds barcode for invariable length barcode data """
                return read1[1][:longbar[0]]
    else:
        def getbarcode(cutters, read1, longbar):
            """ finds barcode for variable barcode lengths"""
            return findbcode(cutters, longbar, read1)

    ## go until end of the file
    while 1:
        try:
            read1, read2 = quarts.next()
            read1 = list(read1)
            total += 1
        except StopIteration:
            break

        ## Parse barcode. Uses the parsing function selected above.
        barcode = getbarcode(cutters, read1, longbar)

        ## find if it matches 
        sname_match = matchdict.get(barcode)

        if sname_match:
            ## record who matched
            dbars[sname_match].add(barcode)
            matched += 1
            cutfound += 1
            samplehits[sname_match] += 1
            barhits[barcode] += 1
            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1

            ## trim off barcode
            lenbar = len(barcode)
            if data.paramsdict["datatype"] == '2brad':
                read1[1] = read1[1][:-lenbar]
                read1[3] = read1[3][:-lenbar]
            else:
                read1[1] = read1[1][lenbar:]
                read1[3] = read1[3][lenbar:]

            ## append to dsort
            dsort1[sname_match].append("".join(read1))
            if 'pair' in data.paramsdict["datatype"]:
                dsort2[sname_match].append("".join(read2))

        else:
            ## record whether cut found                
            if barcode:
                cutfound += 1
                if barcode in misses:
                    misses[barcode] += 1
                else:
                    misses[barcode] = 1
            else:
                misses["_"] += 1

        ## write out at 10K to keep memory low
        if not total % 10000:
            ## write the remaining reads to file"
            writetofile(data, dsort1, 1, filenum, subnum)
            if 'pair' in data.paramsdict["datatype"]:
                writetofile(data, dsort2, 2, filenum, subnum) 
            ## clear out dsorts
            for sample in data.barcodes:
                dsort1[sample] = []
                dsort2[sample] = []

    ## write the remaining reads to file"
    writetofile(data, dsort1, 1, filenum, subnum)
    if 'pair' in data.paramsdict["datatype"]:
        writetofile(data, dsort2, 2, filenum, subnum)        

    ## return stats in saved pickle b/c return_queue is too tiny
    filestats = [subnum, total, cutfound, matched]
    samplestats = [samplehits, barhits, misses, dbars]

    return filestats, samplestats



def writetofile(data, dsort, read, filenum, subnum):
    """ writes dsort dict to a tmp file. Used in barmatch. """
    if read == 1:
        rrr = "R1"
    else:
        rrr = "R2"

    for sname in dsort:
        ## skip writing if empty. Write to tmpname
        #if dsort[sname]:
        tmpdir = os.path.join(data.dirs.fastqs, "tmp_{}_{}_"\
                              .format(sname, rrr))
        handle = os.path.join(tmpdir, "tmp_{}_{}_{}_{}.fastq"\
                              .format(sname, rrr, filenum, subnum))
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))#+"\n")



def collate_subs(args):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    ## parse args
    data, filenum, sublist = args
    #LOGGER.info("in collatesubs %s: %s", filenum, sublist)

    ## output R1 collated file handle
    checkpoint = 0
    while 1:
        out1 = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(data.barcodes.keys()[0]),
                        "coll_{}_R1_{}.fastq".format(filenum, checkpoint))
        if not os.path.exists(out1):
            break
        else:
            checkpoint += 1

    ## sample names are in data.barcodes
    for sname in data.barcodes:

        ## get chunks
        chunks = []
        for subnum in sublist:
            inchunk = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(sname), 
                        "tmp_{}_R1_{}_{}.fastq".format(sname, filenum, subnum))
            chunks.append(inchunk)
        ## one outfile to write to, 
        out1 = os.path.join(data.dirs.fastqs, 
                    "tmp_{}_R1_".format(sname),             
                    "coll_{}_R1_{}.fastq".format(filenum, checkpoint))

        with open(out1, 'a') as tmpout:
            chunks.sort(key=lambda x: int(x.split("_")[-1][:-6]))
            for inchunk in chunks:
                with open(inchunk, 'r') as tmpin:
                    tmpout.write(tmpin.read())
                ## clean up
                os.remove(inchunk)


        ## do second reads
        if "pair" in data.paramsdict["datatype"]:
            ## get chunks
            chunks = []
            for subnum in sublist:
                inchunk = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname), 
                        "tmp_{}_R2_{}_{}.fastq".format(sname, filenum, subnum))
                chunks.append(inchunk)

            out2 = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname),             
                        "coll_{}_R2_{}.fastq".format(filenum, checkpoint))
            with open(out2, 'a') as tmpout:
                chunks.sort(key=lambda x: int(x.split("_")[-1][:-6]))
                for inchunk in chunks:
                    with open(inchunk, 'r') as tmpin:
                        tmpout.write(tmpin.read())
                    ## clean up
                    os.remove(inchunk)



def collate_files(data, filenum):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """

    ## sample names are in data.barcodes
    for sname in data.barcodes:

        ## get chunks
        incols = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(sname), 
                        "coll_*_R1_{}.fastq".format(filenum))
        out1 = os.path.join(data.dirs.fastqs, 
                        "{}_R1_.fastq.gz".format(sname))

        with gzip.open(out1, 'a') as tmpout:
            chunks = glob.glob(incols)
            chunks.sort(key=lambda x: int(x.split("_")[-3]))
            for incol in chunks:
                with open(incol, 'r') as tmpin:
                    print(incol, out1)
                    tmpout.write(tmpin.read())

        if "pair" in data.paramsdict["datatype"]:
            ## get chunks
            incols = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname), 
                        "coll_*_R2_{}.fastq".format(filenum))
            out2 = os.path.join(data.dirs.fastqs, 
                        "{}_R2_.fastq.gz".format(sname))
            with gzip.open(out2, 'a') as tmpout:
                chunks = glob.glob(incols)
                chunks.sort(key=lambda x: int(x.split("_")[-3]))
                for incol in chunks:
                    with open(incol, 'r') as tmpin:
                        tmpout.write(tmpin.read())



def prechecks(data, preview):
    """ 
    Checks before starting analysis. 
    -----------------------------------
    1) Is there data in raw_fastq_path
    2) Is there a barcode file
    3) Is there a workdir and fastqdir
    4) remove old fastq/tmp_sample_R*_ dirs/
    5) make tmpdirs for each barcode name
    6) return file names as pairs (r1, r2) or fakepairs (r1, 1)
    7) get ambiguous cutter resolutions
    8) get optim size
    """

    ## check for data, do glob for fuzzy matching
    if not glob.glob(data.paramsdict["raw_fastq_path"]):
        raise IPyradWarningExit("""\
    No data found in {}. Fix path to data files.""".
    format(data.paramsdict["raw_fastq_path"]))

    ## find longest barcode
    try:
        barlens = [len(i) for i in data.barcodes.values()]
        if len(set(barlens)) == 1:
            longbar = (barlens[0], 'same')
        else:
            longbar = (max(barlens), 'diff')
    except ValueError:
        raise IPyradWarningExit("    Barcodes file not found.")

    ## make sure there is a [workdir] and a [workdir/name_fastqs]
    data.dirs.fastqs = os.path.join(
                        data.paramsdict["project_dir"], data.name+"_fastqs")
    if not os.path.exists(data.paramsdict["project_dir"]):
        os.mkdir(data.paramsdict["project_dir"])
    if not os.path.exists(data.dirs.fastqs):
        os.mkdir(data.dirs.fastqs)

    ## insure no leftover tmp files from a previous run (there shouldn't be)
    oldtmps = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R1_"))
    oldtmps += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R2_"))    
    for oldtmp in oldtmps:
        try:
            shutil.rmtree(oldtmp)
        except OSError as inst:
            ## In some instances nfs creates hidden dot files in directories
            ## that claim to be "busy" when you try to remove them. Don't
            ## kill the run if you can't remove this directory.
            LOGGER.warn("Failed to remove tmpdir {}".format(inst))

    ## make fresh tmpdirs for each sample
    for sample in data.barcodes:
        tmpname = os.path.join(data.dirs.fastqs, "tmp_"+sample+"_R1_")
        os.mkdir(tmpname)
        tmpname = os.path.join(data.dirs.fastqs, "tmp_"+sample+"_R2_")
        os.mkdir(tmpname)

    ## gather raw sequence filenames
    if "pair" in data.paramsdict["datatype"]:
        raws = combinefiles(data.paramsdict["raw_fastq_path"])
    else:
        raws = zip(glob.glob(data.paramsdict["raw_fastq_path"]), iter(int, 1))

    ## returns a list of both resolutions of cut site 1
    ## (TGCAG, ) ==> [TGCAG, ]
    ## (TWGC, ) ==> [TAGC, TTGC]
    ## (TWGC, AATT) ==> [TAGC, TTGC]
    cutters = ambigcutters(data.paramsdict["restriction_overhang"][0])
    assert cutters, "Must enter a `restriction_overhang` for demultiplexing."

    ## set optim chunk size
    if preview:
        optim = ((data._hackersonly["preview_step1"]) // (data.cpus))
    else:
        optim = estimate_optim(raws[0][0], data.cpus)

    ## Build full inverse barcodes dictionary
    matchdict = {}
    bases = set("CATGN")
    poss = set()
    ## do perfect matches
    for sname, barc in data.barcodes.items():
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
                        print("""\
        warning: barcodes {}:{} and {}:{} are within {} base change of each other"""\
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
                                    print("""\
        warning: barcodes {}:{} and {}:{} are within {} base change of each other"""\
        .format(sname, barc, 
                            matchdict[tbar2], data.barcodes[matchdict[tbar2]],
                            data.paramsdict["max_barcode_mismatch"]))

    return raws, longbar, cutters, optim, matchdict



def estimate_optim(testfile, ncpus):
    """ 
    Estimate a reasonable optim value by grabbing a chunk of sequences, 
    decompressing and counting them, to estimate the full file size.
    """
    ## count the len of one file and assume all others are similar len
    insize = os.path.getsize(testfile)
    tmp_file_name = "./tmp-step1-count.fq"
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

    ## break into ncpu chunks for sending to optim lines of data to each cpu
    optim = (inputreads // (ncpus)) + (inputreads % (ncpus))
    optim = int(optim*4)
    LOGGER.info("total reads: %s, total lines: %s", inputreads, inputreads*4)
    LOGGER.info("optim=%s, optim*ncpus=%s", optim, optim*ncpus)

    return optim



def run(data, preview, ipyclient):
    """ 
    Demultiplexes raw fastq files given a barcodes file.
    """

    ## checks on data before starting
    raws, longbar, cutters, optim, matchdict = prechecks(data, preview)
    LOGGER.info('ncpus %s', data.cpus)
    LOGGER.info('optim %s', optim)

    ## Truncate the input fq so it'll run faster. 
    if preview:
        warning = """\
    Running preview mode: subselecting maximum of {} reads\
    """.format(data._hackersonly["preview_step1"], raws)
        ## print the preview message    
        if data._headers:
            print(warning)
        LOGGER.warn(warning)

    ## stat counters for each rawdata file
    perfile = {}
    for rawtuple in raws:
        handle = os.path.splitext(os.path.basename(rawtuple[0]))[0]
        perfile[handle] = {}
        perfile[handle]["ftotal"] = 0
        perfile[handle]["fcutfound"] = 0
        perfile[handle]["fmatched"] = 0
    ## stats for each sample
    fdbars = {}
    fsamplehits = Counter()
    fbarhits = Counter()
    fmisses = Counter()

    ## a dictionary to hold my dictionaries
    statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars

    ## initial progress bar
    start = time.time()

    ## set up parallel client
    lbview = ipyclient.load_balanced_view()

    ## count expected total number of jobs, set total a bit less than 100%
    ## optim is the number of reads per cpu, if it is too large then we 
    ## load too much into memory. So if optim > 1M then we subsample it
    ## which makes things run just a bit slower
    multi = 1
    while optim > int(4e6):
        optim //= 10
        multi *= 10

    nfiles = len(raws)
    totaljobs = data.cpus * nfiles * multi
    done = 0

    ### progress
    LOGGER.info("doing {} chunks [size = {} lines] on {} cpus"\
                .format(totaljobs, optim, data.cpus))

    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(totaljobs, min(totaljobs, done), 
               ' sorting reads  | {}'.format(elapsed))

    ## dictionary to store asyncresults for barmatch jobs
    filesort = {}
    collatesubs = {}

    ## send slices N at a time
    filenum = 0
    for tups in raws:
        ## submit some fake jobs to filesort
        filesort.setdefault(filenum, [])                
        collatesubs.setdefault(filenum, [])        
        for i in range(data.cpus):
            filesort[filenum].append(lbview.apply(time.sleep, 0.1))

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

        ###############################
        subnum = 0
        sublist = []        
        while 1:
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(totaljobs, done, 
               ' sorting reads  | {}'.format(elapsed))

            ## if engines are available add jobs to fill up engines
            ready = [i.ready() for i in filesort[filenum]]

            if any(ready):
                ## grab a 4K slice
                dat1 = list(itertools.islice(rawr1, optim))
                if tups[1]:
                    dat2 = list(itertools.islice(rawr2, optim))
                else:
                    dat2 = [""]

                ## done collate and quit
                if not dat1:
                    if sublist:
                        ## wait for last jobs to finish (I know, not ideal)
                        ipyclient.wait()
                        lbview.apply(collate_subs, [data, filenum, sublist])
                    ## close files
                    io1.close()
                    if tups[1]:
                        io2.close()
                    break

                ## if too many subs, collate
                elif len(sublist) >= 20:
                    lbview.apply(collate_subs, [data, filenum, sublist])
                    sublist = []

                ## updates stat dicts for finished jobs
                asyncs = [i for (i, j) in zip(filesort[filenum], ready) if j]
                for async in asyncs:
                    result = async.result()
                    if result:
                        statdicts = putstats(
                                        result, raws[filenum][0], statdicts)
                    ## purge from memory
                    ipyclient.purge_local_results(async)

                ## update dictionary to know about missing asyncs
                filesort[filenum] = [i for (i, j) in \
                                     zip(filesort[filenum], ready) if not j]

                ## build a list of args
                args = [data, filenum, subnum, cutters, longbar,
                        matchdict, dat1, dat2]

                ## send the list to N engines
                async = lbview.apply_async(barmatch, args)
                filesort[filenum].append(async)

                ## advance subs
                sublist.append(subnum)
                subnum += 1
                done += 1
        filenum += 1
        ## collate colls from across files
        ipyclient.wait()


    for filenum in range(filenum):
        for async in filesort[filenum]:
            result = async.get()
            if result:
                statdicts = putstats(result, raws[filenum][0], statdicts)
        ## collate across files
        collate_files(data, filenum)

    ## clean up junk files
    tmpdirs = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*_"))
    for tmpdir in tmpdirs:
        try:
            shutil.rmtree(tmpdir)
        except Exception:
            pass

    ## final progress
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 20,
               ' sorting reads  | {}'.format(elapsed))

    ## build stats from dictionaries
    perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts    
    make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)



def putstats(result, handle, statdicts):
    """ puts stats from async into dictionaries """
    filestats, samplestats = result
    _, total, cutfound, matched = filestats
    samplehits, barhits, misses, dbars = samplestats

    ## get dicts
    perfile, fsamplehits, fbarhits, fmisses, fdbars = statdicts

    ## update file stats
    handle = os.path.splitext(os.path.basename(handle))[0]
    perfile[handle]["ftotal"] += total
    perfile[handle]["fcutfound"] += cutfound
    perfile[handle]["fmatched"] += matched    

    ## update sample stats
    fsamplehits.update(samplehits)
    fbarhits.update(barhits)        
    fmisses.update(misses)
    fdbars.update(dbars)

    statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars
    return statdicts



        ## send new slices from iterator as others finish
        # done = 0
        # sublist = []
        # while 1:
        #     time.sleep(0.1)
        #     ## progress bar is not perfect in this case
        #     #elapsed = datetime.timedelta(seconds=int(time.time()-start))
        #     #progressbar(totaljobs, min(totaljobs, chunkdone+done), 
        #     #    ' sorting reads  | {}'.format(elapsed))

        #     ## check if any engine jobs finished
        #     subsyncs = filesort[filenum]

        #     ## for each individual job, if done, send new, store dones
        #     for ridx, rasync in enumerate(subsyncs):
        #         if rasync.ready():
        #             ## if done then get stats, bounces if not done
        #             results = rasync.get(0)
        #             filestats, samplestats = results
        #             rsubnum, total, cutfound, matched = filestats
        #             samplehits, barhits, misses, dbars = samplestats

        #             ## save that this sub finished so we can add it to collate
        #             sublist.append(rsubnum)

        #             ## update file stats
        #             handle = os.path.splitext(
        #                         os.path.basename(raws[filenum][0]))[0]
        #             perfile[handle]["ftotal"] += total
        #             perfile[handle]["fcutfound"] += cutfound
        #             perfile[handle]["fmatched"] += matched    

        #             ## update sample stats
        #             fsamplehits.update(samplehits)
        #             fbarhits.update(barhits)        
        #             fmisses.update(misses)
        #             fdbars.update(dbars) 

        #             ## store as done
        #             del filesort[filenum][ridx]
        #             done += 1

        #             ## grab a new 40K slice and send to open Engine
        #             dat1 = list(itertools.islice(rawr1, chunksize))
        #             if tups[1]:
        #                 dat2 = list(itertools.islice(rawr2, chunksize))
        #             else:
        #                 dat2 = [""]

        #             if dat1:
        #                 args = [data, filenum, subnum, cutters, longbar,
        #                         matchdict, dat1, dat2]
        #                 async = lbview.apply(barmatch, args)
        #                 filesort[filenum].append(async)
        #                 subnum += 1
        #             else:
        #                 continue

        #     ## if all jobs in key are done then collate subs and break
        #     if done == subnum:
        #         if sublist:
        #             casync = lbview.apply(collate_subs, 
        #                                   [data, filenum, sublist])
        #             collatesubs[filenum].append(casync)
        #         ## close files
        #         io1.close()
        #         if tups[1]:
        #             io2.close()
        #         break

        #     ## if 20 chunks are done then collate them
        #     elif len(sublist) >= 20:
        #         casync = lbview.apply(collate_subs, [data, filenum, sublist])
        #         collatesubs[filenum].append(casync)
        #         sublist = []
            
        # filenum += 1
        # chunkdone += done
    
    ## set wait job until all finished. THen collate files
    # tmpids = []
    # for fnum in collatesubs:
    #     tmpl = list(itertools.chain(*collatesubs[fnum]))
    #     tmpids.append(tmpl)
    # tmpids = list(itertools.chain(*collatesubs[tmpids]))

    # print(tmpids)
    # with lbview.temp_flags(after=tmpids):
    #     res = lbview.apply(collate_files, data)

    # ## print final progress
    # elapsed = datetime.timedelta(seconds=int(time.time()-start))
    # progressbar(totaljobs, totaljobs, 
    #             ' sorting reads  | {}'.format(elapsed))
    # if data._headers:
    #     print("")
    # ## wait on final collate files job. progress bar?
    # res.wait()
    ipyclient.wait()

    # except (KeyboardInterrupt, SystemExit):
    #     print('\n  Interrupted! Cleaning up... ')
    #     raise

    ## clean up junk files
    tmpdirs = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*_"))
    for tmpdir in tmpdirs:
        try:
            shutil.rmtree(tmpdir)
        except Exception:
            pass

    #     ## if success: make stats
    ## build stats from dictionaries
    make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)


    #     for job in range(filenum):
    #         for async in sorting[job]:
    #             if async.ready():
    #                 if async.successful():
    #                     filestats, samplestats = async.get()
    #                     handle, total, cutfound, matched = filestats
    #                     samplehits, barhits, misses, dbars = samplestats

    #                     ## update file stats
    #                     perfile[handle]["ftotal"] += total
    #                     perfile[handle]["fcutfound"] += cutfound
    #                     perfile[handle]["fmatched"] += matched    

    #                     ## update sample stats
    #                     fsamplehits.update(samplehits)
    #                     fbarhits.update(barhits)        
    #                     fmisses.update(misses)
    #                     fdbars.update(dbars)

    #                 else:
    #                     ## grab metadata and complain
    #                     meta = async.metadata
    #                     ## if not done do nothing, if failure print error
    #                     if meta.error:
    #                         LOGGER.error("error: %s", meta.error)

                         



if __name__ == "__main__":

    ## run test
    import ipyrad as ip
    #from ipyrad.core.assembly import Assembly

    ## get current location
    #PATH = os.path.abspath(os.path.dirname(__file__))
    ## get location of test files
    #IPATH = os.path.dirname(os.path.dirname(PATH))
    #DATA = os.path.join(IPATH, "tests", "test_rad")

    TEST = ip.Assembly("profile_s1")
    TEST.set_params(1, "./maintest")
    TEST.set_params(2, "./ipsimdata/sim_rad_test_R1_.fastq.gz")
    TEST.set_params(3, "./ipsimdata/sim_rad_test_barcodes.txt")
    print(TEST.cpus)
    TEST.cpus = 4
    TEST.step1()
