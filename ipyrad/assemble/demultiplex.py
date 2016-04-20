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
import cPickle as pickle
from ipyrad.core.sample import Sample
from .util import *
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
        raise AssertionError("First read files names must contain '_R1_'.")

    ## get paired reads
    seconds = [ff.replace("_R1_", "_R2_") for ff in firsts]
    return zip(firsts, seconds)



def matching(barcode, data):
    "allows for N base difference between barcodes"
    for name, realbar in data.barcodes.items():
        if len(barcode) == len(realbar):
            sames = sum([i == j for (i, j) in zip(barcode, realbar)])
            diffs = len(barcode) - sames
            if diffs <= data.paramsdict["max_barcode_mismatch"]:
                return name
    return 0



def findbcode(cut, longbar, read1):
    """ find barcode sequence in the beginning of read """
    ## default barcode string
    search = read1[1][:int(longbar[0]+len(cut)+1)]
    countcuts = search.count(cut)
    if countcuts == 1:
        barcode = search.split(cut, 1)[0]
    elif countcuts == 2:
        barcode = search.rsplit(cut, 2)[0]
    else:
        barcode = ""
    return barcode



def barmatch(args):
    """
    Matches reads to barcodes in barcode file and writes to individual temp 
    files, after all read files have been split, temp files are collated into 
    .fastq files
    """

    ## read in list of args
    data, rawtuple, filenum, subnum, cut, longbar = args
    LOGGER.info("Entering barmatch - {} | {}.{}"\
                 .format(rawtuple, filenum, subnum))

    ## counters for total reads, those with cutsite, and those that matched
    total = 0
    cutfound = 0
    matched = 0 

    ## dictionary to record barcode hits & misses
    samplehits = {}
    barhits = {}
    misses = {}
    misses['_'] = 0
    
    ## read in paired reads files
    tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    chunk = (os.path.join(tmpdir, "tmp_R1_{}_{}.fastq".format(filenum, subnum)),
             os.path.join(tmpdir, "tmp_R2_{}_{}.fastq".format(filenum, subnum)))

    fr1 = open(chunk[0], 'rb')
    ## create iterators to sample 4 lines at a time
    quart1 = itertools.izip(*[iter(fr1)]*4)
    if 'pair' in data.paramsdict["datatype"]:
        fr2 = open(chunk[1], 'rb')
        quart2 = itertools.izip(*[iter(fr2)]*4)
        quarts = itertools.izip(quart1, quart2)
    else:
        ## read in single end read file"
        quarts = itertools.izip(quart1, iter(int, 1))

    ## dictionaries to store first and second reads
    dsort1 = defaultdict(list)
    dsort2 = defaultdict(list)
    dbars = defaultdict(list)    ## all bars matched in sample
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
        def getbarcode(cut, read1, longbar):
            """ finds barcode for variable barcode lengths"""
            return findbcode(cut, longbar, read1)

    ## go until end of the file
    while 1:
        try: 
            read1s, read2s = quarts.next()
        except StopIteration: 
            break
        total += 1

        ## strip
        read1 = [i.strip() for i in read1s]
        if 'pair' in data.paramsdict["datatype"]:
            read2 = [i.strip() for i in read2s]

        ## Parse barcode. Use the parsing function selected above.
        barcode = getbarcode(cut, read1, longbar)

        ## find if it matches 
        didmatch = matching(barcode, data)
        if didmatch:

            ## record who matched
            dbars[didmatch].add(barcode)
            matched += 1
            cutfound += 1
            if didmatch in samplehits:
                samplehits[didmatch] += 1
            else:
                samplehits[didmatch] = 1

            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1

            ## trim off barcode
            if data.paramsdict["datatype"] == '2brad':
                read1[1] = read1[1][:-len(barcode)]
                read1[3] = read1[3][:-len(barcode)]
            else:
                read1[1] = read1[1][len(barcode):]
                read1[3] = read1[3][len(barcode):]

            ## append to dsort
            dsort1[didmatch].append("\n".join(read1).strip())
            if 'pair' in data.paramsdict["datatype"]:
                dsort2[didmatch].append("\n".join(read2).strip())

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

    ## close file handles
    fr1.close()
    if 'pair' in data.paramsdict["datatype"]:
        fr2.close()

    ## return stats in saved pickle b/c return_queue is too tiny
    handle = os.path.splitext(os.path.basename(rawtuple[0]))[0]
    filestats = [handle, total, cutfound, matched]
    samplestats = [samplehits, barhits, misses, dbars]

    return filestats, samplestats



def writetofile(data, dsort, read, filenum, chunknum):
    """ writes dsort dict to a tmp file. Used in barmatch. """
    if read == 1:
        rname = "_R1_"
    else:
        rname = "_R2_"

    for sample in dsort:
        ## skip writing if empty. Write to tmpname
        if dsort[sample]:
            tmpdir = os.path.join(data.dirs.fastqs, "tmp_"+sample+rname)
            handle = os.path.join(tmpdir, "tmp_"+sample+rname+\
                                  str(filenum)+"_"+str(chunknum))
            with open(handle, 'a+') as out:
                out.write("\n".join(dsort[sample])+"\n")



def collate_subs(args):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    ## parse args
    data, filenum = args

    ## sample names are in data.barcodes
    for sname in data.barcodes:

        ## get chunks
        inchunks = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(sname), 
                        "tmp_{}_R1_{}_*.fastq".format(sname, filenum))

        ## one outfile to write to, 
        out1 = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(sname),             
                        "coll_{}_R1_.fastq.gz".format(filenum))
        with gzip.open(out1, 'a') as tmpout:
            for inchunk in inchunks:
                with open(inchunk, 'r') as tmpin:
                    tmpout.write(tmpin.read())

        if "pair" in data.paramsdict["datatype"]:
            ## get chunks
            inchunks = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname), 
                        "tmp_{}_R2_{}_*.fastq".format(sname, filenum))
            ## one outfile to write to, 
            out2 = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname),             
                        "coll_{}_R2_.fastq.gz".format(filenum))
            with gzip.open(out2, 'a') as tmpout:
                for inchunk in inchunks:
                    with open(inchunk, 'r') as tmpin:
                        tmpout.write(tmpin.read())



def collate_files(data):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    ## sample names are in data.barcodes
    for sname in data.barcodes:

        ## get chunks
        incols = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R1_".format(sname), 
                        "coll_{}_R1_*.fastq.gz".format(sname))
        out1 = os.path.join(data.dirs.fastqs, 
                        "{}_R1_.fastq.gz".format(sname))

        with gzip.open(out1, 'w') as tmpout:
            for incol in incols:
                with gzip.open(incol, 'r') as tmpin:
                    tmpout.write(tmpin.read())

        if "pair" in data.paramsdict["datatype"]:
            ## get chunks
            incols = os.path.join(data.dirs.fastqs, 
                        "tmp_{}_R2_".format(sname), 
                        "coll_{}_R2_*.fastq.gz".format(sname))
            out2 = os.path.join(data.dirs.fastqs, 
                        "{}_R2_.fastq.gz".format(sname))
            with gzip.open(out2, 'w') as tmpout:
                for incol in incols:
                    with gzip.open(incol, 'r') as tmpin:
                        tmpout.write(tmpin.read())




def prechecks(data, ipyclient, preview):
    """ 
    Checks before starting analysis. 
    -----------------------------------
    1) Is there data in raw_fastq_path
    2) Is there a barcode file
    3) Is there a workdir and fastqdir
    4) remove old tmpdirs
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

    ## gather raw sequence data
    if "pair" in data.paramsdict["datatype"]:
        raws = combinefiles(data.paramsdict["raw_fastq_path"])
    else:
        raws = zip(glob.glob(data.paramsdict["raw_fastq_path"]), iter(int, 1))

    ## returns a list of both resolutions of cut site 1
    ## (TGCAG, ) ==> [TGCAG, ]
    ## (TWGC, ) ==> [TAGC, TTGC]
    ## (TWGC, AATT) ==> [TAGC, TTGC]
    cutters = [ambigcutters(i) for i in \
               data.paramsdict["restriction_overhang"]][0]
    assert cutters, "Must enter a `restriction_overhang` for demultiplexing."

    ## set optim chunk size
    ncpus = detect_cpus()
    if preview:
        optim = ((data._hackersonly["preview_step1"]) // ncpus) 

    else:
        ## count the len of one file and assume all others are similar len
        testfile = raws[0][0]
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
        for i in xrange(40000):
            outfile.write(next(infile))
        outfile.close()
        infile.close()

        ## Get the size of the tmp file
        tmp_size = os.path.getsize(tmp_file_name)

        ## divide by the tmp file size and multiply by 10000 to approximate
        ## the size of the input .fq files
        inputreads = int(insize / tmp_size) * 10000
        #LOGGER.info("inputreads estimate - {}, lines estimate - {}"\
        #            .format(inputreads, inputreads*4))
        os.remove(tmp_file_name)

        ## it's going to be multiplied by 4 to ensure its divisible
        ## and then again by 4 if inside preview truncate, so x32 here.
        ## should result in 2X as many chunk files as cpus, but files are 
        ## split by reads not by lines so divide this by 4 to get optim
        ## num reads to split per tmp file. 
        ## send 4 chunks to each cpu
        optim = (inputreads // (ncpus * 2)) + \
                (inputreads % (ncpus * 2))

    LOGGER.info("optim=%s, optim*4=%s", optim, optim*4)
    return raws, longbar, cutters, optim



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
        sample.stats["reads_raw"] = fsamplehits[name]
        ## fill in the full df stats value
        sample.stats_dfs.s1["reads_raw"] = fsamplehits[name]

        ## Only link Sample if it has data
        if sample.stats["reads_raw"]:
            sample.stats.state = 1
            data.samples[sample.name] = sample
        else:
            print("Excluded sample: no data found for", name)

    ## initiate s1 key for data object
    data.stats_dfs.s1 = data.build_stat("s1")
    data.stats_files.s1 = outhandle



def run(data, preview, ipyclient):
    """ 
    Demultiplexes raw fastq files given a barcodes file.
    """

    ## checks on data before starting
    raws, longbar, cutters, optim = prechecks(data, ipyclient, preview)

    ## Truncate the input fq so it'll run faster. 
    if preview:
        warning = """\
    Running preview mode: subselecting maximum of {} reads\
    """.format(data._hackersonly["preview_step1"], raws)
        ## print the preview message    
        if data._headers:
            print(warning)
        LOGGER.warn(warning)

        ## enter step1 truncate length, returns one truncated fq (file,)
        nlines = data._hackersonly["preview_step1"]
        raws = preview_truncate_fq(data, raws, nlines)

    ## initial progress bar
    start = time.time()

    ## set up parallel client
    lbview = ipyclient.load_balanced_view()

    ## make tmpdir to hold chunks
    tmpdir = os.path.join(data.dirs.project, data.name+"-tmpchunks")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    ## dictionary to store asyncresults for barmatch jobs
    sorting = {}

    ## set up a job to send slices to be sorted on the queue
    filenum = 0
    for rawtuple in raws:
        sorting[filenum] = []

        ## open file handles
        if rawtuple[0].endswith(".gz"):
            rawr1 = iter(gzip.open(rawtuple[0]))
            if rawtuple[1]:
                rawr2 = iter(gzip.open(rawtuple[1]))
        else:
            rawr1 = iter(open(rawtuple[0]))
            if rawtuple[1]:
                rawr2 = iter(open(rawtuple[1]))

        subnum = 0
        while 1:
            ## slice out R1 chunk
            dat1 = "".join(itertools.islice(rawr1, int(optim)*4))
            if dat1:
                out1 = os.path.join(tmpdir, "tmp_R1_{}_{}.fastq"\
                                      .format(filenum, subnum))
                with open(out1, 'w') as tmpout:
                    tmpout.write(dat1)

                ## if paired write R2s
                if "pair" in data.paramsdict["datatype"]:
                    dat2 = "".join(itertools.islice(rawr2, int(optim)*4)) 
                    out2 = os.path.join(tmpdir, "tmp_R1_{}_{}.fastq"\
                                      .format(filenum, subnum))
                    with open(out2, 'w') as tmpout:
                        tmpout.write(dat2)

                ## apply job to be sorted
                args = [data, rawtuple, filenum, subnum, cutters, longbar]
                sorting[filenum].append(lbview.apply(barmatch, args))
                subnum += 1
            else:
                break

        ## set up a job to collate sorted slices WITHIN each rawtuple
        collate_s = {}
        for subnum in range(subnum):
            tmpids = list(itertools.chain(\
                         *[i.msg_ids for i in sorting[filenum]]))
            #LOGGER.info("len %s, tmpids; %s", len(tmpids), tmpids)

            ## do not let multiple attempts collate at once!!
            ## this requires that all samples within filenum are done
            with lbview.temp_flags(after=tmpids):
                collate_s[filenum] = lbview.apply(collate_subs, 
                                           [data, filenum, subnum])
        filenum += 1

    ## set wait job until all finished. 
    tmpids = list(itertools.chain(*[i.msg_ids for i in collate_s.values()]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(collate_files)

    try:
        allsorts = list(itertools.chain(*[sorting[i] for i in sorting]))
        allwait = len(allsorts)
        while 1:
            if not res.ready():
                fwait = sum([i.ready() for i in allsorts])
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(allwait, fwait, 
                            " demultiplexing reads  | {}".format(elapsed))
                ## got to next print row when done
                sys.stdout.flush()
                time.sleep(1)
            else:
                ## print final statement
                elapsed = datetime.timedelta(seconds=int(time.time()-start))                
                progressbar(20, 20, 
                            " demultiplexing reads  | {}".format(elapsed))
                break
        if data._headers:
            print("")

    except (KeyboardInterrupt, SystemExit):
        print('\n  Interrupted! Cleaning up... ')
        raise

    finally:
        ## clean up junk files
        tmpdirs = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*_"))
        tmpdirs += glob.glob(os.path.join(data.dirs.project, 
                                          data.name+"-tmpchunks"))
        for tmpdir in tmpdirs:
            try:
                shutil.rmtree(tmpdir)
            except Exception:
                pass

        ## if success: make stats
        ## stats for each rawdata file
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

        for job in range(filenum):
            for async in sorting[job]:
                if async.successful():
                    filestats, samplestats = async.get()
                    handle, total, cutfound, matched = filestats
                    samplehits, barhits, misses, dbars = samplestats

                    ## update file stats
                    perfile[handle]["ftotal"] += total
                    perfile[handle]["fcutfound"] += cutfound
                    perfile[handle]["fmatched"] += matched    

                    ## update sample stats
                    fsamplehits.update(samplehits)
                    fbarhits.update(barhits)        
                    fmisses.update(misses)
                    fdbars.update(dbars)

                else:
                    ## grab metadata and complain
                    meta = async.metadata
                    ## if not done do nothing, if failure print error
                    LOGGER.error('  error occurred')
                    if meta.error:
                        LOGGER.error("""\
                stdout: %s
                stderr: %s 
                error: %s""", meta.stdout, meta.stderr, meta.error)

        ## build stats from dictionaries
        make_stats(data, perfile, fsamplehits, fbarhits, fmisses, fdbars)
                         


    # finally:
    #     try:
    #         ## cleans up chunk files and stats pickles
    #                 except OSError:
    #                     ## encountered a strange error once where nfs creates 
    #                     ## hidden dot files in directories that claim to be 
    #                     ## "busy" when you try to remove them. Don't
    #                     ## kill the run if you can't remove this directory.
    #                     LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

    #         ## cleans up pickle files and tmp files generated by preview
    #         tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "*.pickle"))
    #         tmpfiles += glob.glob(os.path.join(data.dirs.fastqs, "*.preview_t"))            
    #         for tmpfile in tmpfiles:
    #             if os.path.exists(tmpfile):
    #                 os.remove(tmpfile)

    #     except AttributeError as inst:
    #         ## If barcodes file is fsck, then finally frags because 
    #         ## data.dirs.fastqs doesn't exist
    #         LOGGER.warn(inst)



if __name__ == "__main__":

    ## run test
    import ipyrad as ip
    #from ipyrad.core.assembly import Assembly

    ## get current location
    PATH = os.path.abspath(os.path.dirname(__file__))
    ## get location of test files
    IPATH = os.path.dirname(os.path.dirname(PATH))
    
    DATA = os.path.join(IPATH, "tests", "test_rad")
    TEST = ip.Assembly("test-demultiplex")
    #TEST = ip.load_assembly(os.path.join(DATA, "testrad"))
    TEST.set_params(1, "./")
    TEST.set_params(2, "./tests/data/sim_rad_test_R1_.fastq.gz")
    TEST.set_params(3, "./tests/data/sim_rad_test_barcodes.txt")
    TEST.step1()
