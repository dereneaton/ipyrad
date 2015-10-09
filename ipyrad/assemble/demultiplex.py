#!/usr/bin/env python2

""" demultiplex raw sequence data given a barcode map """

from __future__ import print_function
import os
import sys
import gzip
import glob
import tempfile
import itertools
#import multiprocessing
import cPickle as pickle
import ipyparallel as ipp
from ipyrad.core.sample import Sample
from collections import defaultdict, Counter


def combinefiles(filepath):
    """ Joins first and second read file names """
    ## unpack seq files in filepath
    fastqs = glob.glob(filepath)
    firsts = [i for i in fastqs if "_R1_" in i]

    ## check names
    if not firsts:
        sys.exit("\n\tFirst read files names must contain '_R1_'.")

    ## get paired reads
    seconds = [ff.replace("_R1_", "_R2_") for ff in firsts]

    return zip(firsts, seconds)



def revcomp(sequence):
    "returns reverse complement of a string"
    sequence = sequence[::-1].strip()\
                             .replace("A", "t")\
                             .replace("T", "a")\
                             .replace("C", "g")\
                             .replace("G", "c").upper()
    return sequence


def matching(barcode, data):
    "allows for N base difference between barcodes"
    match = 0
    for name, realbar in data.barcodes.items():
        if len(barcode) == len(realbar):
            sames = sum([i == j for (i, j) in zip(barcode, realbar)])
            diffs = len(barcode) - sames
            if diffs <= data.paramsdict["max_barcode_mismatch"]:
                match = name
    return match



def ambigcutters(seq):
    """ returns both resolutions of a cut site
    that has an ambiguous base in it, else the 
    single cut site """
    resos = []
    ambigs = {"R":("G", "A"),
              "K":("G", "T"),
              "S":("G", "C"),
              "Y":("T", "C"),
              "W":("T", "A"),
              "M":("C", "A")}
    if any([i in list("RKSYWM") for i in seq]):
        for base in list("RKSYWM"):
            if base in seq:
                resos.append(seq.replace(base, ambigs[base][0]))
                resos.append(seq.replace(base, ambigs[base][1]))
        return resos
    else:
        return [seq, ""]



def findbcode(cut, longbar, read1):
    """ find barcode sequence in the beginning of read """
    ## default barcode string
    #barcode = 'N'*20
    search = read1[1][:longbar+len(cut)]
    countcuts = search.count(cut)
    if countcuts == 1:
        barcode = search.split(cut, 1)[0]
    elif countcuts == 2:
        barcode = search.rsplit(cut, 2)[0]
    else:
        barcode = ""
    return barcode



def barmatch(args):
    """matches reads to barcodes in barcode file
    and writes to individual temp files, after all
    read files have been split, temp files are collated
    into .fq files"""

    ## read in list of args
    data, rawfile, chunk, cut, longbar, chunknum, filenum = args

    ## is paired?
    paired = bool('pair' in data.paramsdict["datatype"])

    ## counters for stats output
    total = 0
    cutfound = 0      ## cut site found
    matched = 0       ## bar matches 

    ## dictionary to record barcode hits & misses
    samplehits = {}
    barhits = {}
    misses = {}
    misses['_'] = 0
    
    ## read in paired end read files"
    fr1 = open(chunk[0], 'rb')
    ## create iterators to sample 4 lines at a time
    quart1 = itertools.izip(*[iter(fr1)]*4)
    if paired:
        fr2 = open(chunk[1], 'rb')
        quart2 = itertools.izip(*[iter(fr2)]*4)
        quarts = itertools.izip(quart1, quart2)
    else:
        ## read in single end read file"
        quarts = itertools.izip(quart1, iter(int, 1))

    ## dictionaries to store first and second reads
    dsort1 = defaultdict(list)
    dsort2 = defaultdict(list)
    dbars = defaultdict(list)   ## all bars matched in sample
    for sample in data.barcodes:
        dsort1[sample] = []
        dsort2[sample] = []
        dbars[sample] = set()

    ## go until end of the file
    while 1:
        try: 
            read1s, read2s = quarts.next()
        except StopIteration: 
            break
        total += 1

        ## strip
        read1 = [i.strip() for i in read1s]
        if paired:
            read2 = [i.strip() for i in read2s]

        ## Parse barcode from sequence 
        ## very simple sorting if barcodes are invariable
        if longbar[1] == 'same':
            if data.paramsdict["datatype"] == '2brad':
                barcode = read1[1][-longbar[0]:]
            else:
                barcode = read1[1][:longbar[0]]
        else:
            barcode = findbcode(cut, longbar[0], read1)

        ## find if it matches 
        didmatch = matching(barcode, data)
        if didmatch:
            dbars[didmatch].add(barcode)
            matched += 1
            cutfound += 1
            ## trim off barcode
            if data.paramsdict["datatype"] == '2brad':
                read1[1] = read1[1][:-len(barcode)]
                read1[3] = read1[3][:-len(barcode)]
            else:
                read1[1] = read1[1][len(barcode):]
                read1[3] = read1[3][len(barcode):]
            ## record who matched
            if didmatch in samplehits:
                samplehits[didmatch] += 1
            else:
                samplehits[didmatch] = 1

            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1                
            ## append to dsort
            dsort1[didmatch].append("\n".join(read1).strip())
            if paired:
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

    ## write the remaining reads to file"
    for sample in dsort1:
        handle = os.path.join(data.dirs.fastqs,
                              "tmp_"+sample+"_R1_"+\
                              str(filenum)+"_"+str(chunknum)+".gz")
        with gzip.open(handle, 'w') as out:
            out.write("\n".join(dsort1[sample])+"\n")
    if paired:
        for sample in dsort2:
            handle = os.path.join(data.dirs.fastqs,
                                  "tmp_"+sample+"_R2_"+\
                                  str(filenum)+"_"+str(chunknum)+".gz")
            with gzip.open(handle, 'w') as out:
                out.write("\n".join(dsort2[sample])+"\n")

    ## close file handles
    fr1.close()
    if paired:
        fr2.close()

    ## return stats in saved pickle b/c return_queue is too tiny
    handle = os.path.splitext(os.path.basename(rawfile))[0]
    filestats = [handle, total, cutfound, matched]
    samplestats = [samplehits, barhits, misses, dbars]

    pickout = open(os.path.join(
                      data.dirs.fastqs,
                      handle+"_"+str(chunknum)+"_"+\
                      str(filenum)+".pickle"), "wb")
    pickle.dump([filestats, samplestats], pickout)
    pickout.close()



def blocks(files, size=50120000):
    """ reads in block 5MB in size"""
    while 1:
        block = files.read(size)
        if not block:
            break
        yield block



def chunker(args):
    """ splits fastq file into nchunks, preserving that
    the data are in 4-line groups. Split across N processors. """
    print("inside this chunk")
    data, fastq, paired, num, optim, pickleout = args
    ## is gzipped?
    gzipped = bool(fastq[0].endswith(".gz"))

    ## count nlines in read1 file
    if gzipped:
        with gzip.open(fastq[0], 'rb') as indata:
            totlen = sum([ind.count("\n") for ind in blocks(indata)])
    else:
        with open(fastq[0], 'rb') as indata:
            totlen = sum([ind.count("\n") for ind in blocks(indata)])
    added = 0

    ## data in at optim lines at a time
    if gzipped:
        indata1 = gzip.open(fastq[0])
    else:
        indata1 = open(fastq[0])
    ## 
    if paired:
        if gzipped:
            indata2 = gzip.open(fastq[1])
        else:
            indata2 = open(fastq[1])            

    ## list for tempname tuples
    chunks1 = []
    chunks2 = []
    while added < totlen:
        with tempfile.NamedTemporaryFile('w+b', delete=False) as out:
            out.write("".join(list(itertools.islice(indata1, optim))))
        chunks1.append(out.name)
        if paired:
            with tempfile.NamedTemporaryFile('w+b', delete=False) as out:
                out.write("".join(list(itertools.islice(indata2, optim))))
            chunks2.append(out.name)
        else:
            chunks2.append(tempfile.NamedTemporaryFile(delete=False).name)
        added += optim

    ## cleanup
    indata1.close()
    if paired:
        indata2.close()

    ## returns a list with tuple of (tmp1, "") or (tmp1, tmp2) and chunk size
    ## GEEEEEz, return queue is too small again, causes freezing
    ## use pickles instead, again.
    ## return [fastq[0], zip(chunks1, chunks2)]
    if pickleout:
        pickout = open(os.path.join(
                          data.dirs.fastqs, 
                          "chunks_"+str(num)+".pickle"), "wb")
        pickle.dump([fastq[0], zip(chunks1, chunks2)], pickout)
        pickout.close()
    else:
        return [fastq[0], zip(chunks1, chunks2)]



def parallel_chunker(data, raws, paired):
    """ iterate over raw data files and split into N pieces """
    ## count how many rawfiles have been done
    #work_queue = multiprocessing.Queue()    
    submitted_args = []
    #chunk_queue = multiprocessing.Queue()        
    num = 0
    for rawtuple in list(raws):
        ## for each (R1, R2) put on work queue, w/ fixblock=1
        #work_queue.put([data, rawtuple, 
        #                paired, num, 400000, 1])
        submitted_args.append([data, rawtuple, paired, num, 400000, 1])
        num += 1

    ## call to ipp
    ipyclient = ipp.Client()
    workers = ipyclient.load_balanced_view()
    res = workers.map(chunker, submitted_args)
    res.get() #ipyclient.wait()

    #apply_async(chunker, submitted_args)
    ## spawn workers, run barmatch function"
    # jobs = []        
    # for _ in xrange(data.paramsdict["N_processors"]):
    #     work = worker.Worker(work_queue, chunk_queue, chunker)
    #     work.start()
    #     jobs.append(work)
    # for job in jobs:
    #     job.join()



def parallel_sorter(data, rawfilename, chunks, cutter, longbar, filenum):
    """ takes list of chunk files and runs barmatch function
    on them across N processors and outputs temp file results.
    """
    ## send file to multiprocess queue"
    chunknum = 0
    submitted_args = []
    for tmptuple in chunks:
        submitted_args.append([data, rawfilename, tmptuple, cutter,
                               longbar, chunknum, filenum])
        chunknum += 1

    ## uses all available processors
    ipyclient = ipp.Client()
    workers = ipyclient.load_balanced_view()
    res = workers.map_async(barmatch, submitted_args)
    res.get()  

 

def collate_tmps(data, paired):
    """ collate temp files back into 1 sample """
    for name in data.barcodes:
        ## nproc len list of chunks
        combs = glob.glob(os.path.join(
                          data.dirs.fastqs, "tmp_"+name)+"_R1_*.gz")
        combs.sort(key=lambda x: int(x.split("_")[-1].replace(".gz", "")[0]))

        ## one outfile to write to
        handle_r1 = os.path.join(data.dirs.fastqs, name+"_R1_.gz")
        with gzip.open(handle_r1, 'w') as out:
            for fname in combs:
                with gzip.open(fname) as infile:
                    out.write(infile.read())
        if paired:
            ## nproc len list of chunks
            combs = glob.glob(os.path.join(
                              data.dirs.fastqs, "tmp_"+name)+"_R2_*.gz")
            combs.sort()                        
            ## one outfile to write to
            handle_r2 = os.path.join(data.dirs.fastqs, name+"_R2_.gz")
            with gzip.open(handle_r2, 'w') as out:
                for fname in combs:
                    with gzip.open(fname) as infile:
                        out.write(infile.read())



def prechecks(data, preview):
    """ todo before starting analysis """
    ## check for data
    if not glob.glob(data.paramsdict["raw_fastq_path"]):
        sys.exit("\tNo data found in "+data.paramsdict["raw_fastq_path"]+\
                 ". Fix path to the data files\n")

    ## find longest barcode
    barlens = [len(i) for i in data.barcodes.values()]
    if len(set(barlens)) == 1:
        longbar = (barlens[0], 'same')
    else:
        longbar = (max(barlens), 'diff')

    ## make sure there is an out directory
    data.dirs.fastqs = os.path.join(data.paramsdict["working_directory"],
                                    'fastq')
    if not os.path.exists(data.paramsdict["working_directory"]):
        os.mkdir(data.paramsdict["working_directory"])
    if not os.path.exists(data.dirs.fastqs):
        os.mkdir(data.dirs.fastqs)

    ## if leftover tmp files, remove
    oldtmps = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.gz"))
    if oldtmps:
        for oldtmp in oldtmps:
            os.remove(oldtmp)

    ## gather raw sequence data
    paired = bool("pair" in data.paramsdict["datatype"])
    if paired:
        raws = combinefiles(data.paramsdict["raw_fastq_path"])
    else:
        raws = zip(glob.glob(data.paramsdict["raw_fastq_path"]), iter(int, 1))

    ## preview rawdata files
    if preview:
        for sample in data.barcodes:
            print("preview:", sample, data.barcodes[sample])
        for rawfile in raws:
            print("preview:", rawfile)

    ## make list of all perfect matching cut sites
    cut1, _ = [ambigcutters(i) for i in \
               data.paramsdict["restriction_overhang"]]

    return raws, longbar, cut1, paired




def make_stats(data, raws, paired):
    """ reads in pickled stats, collates, and writes to file """
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

    ## get stats from each file pickle
    pickles = glob.glob(os.path.join(data.dirs.fastqs, "*.pickle"))
    for picfile in pickles:
        with open(picfile, "rb") as pickin:
            filestats, samplestats = pickle.load(pickin)

        #counts = [total, cutfound, matched]
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


    data.statsfiles.s1 = os.path.join(data.dirs.fastqs, 
                                      's1_demultiplex_stats.txt')
    outfile = open(data.statsfiles.s1, 'w')

    ## how many from each rawfile
    outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                  format("raw_file", "total_reads", 
                         "cut_found", "bar_matched"))
    ## sort rawfile names
    rawfilenames = sorted(perfile)
    for rawstat in rawfilenames:
        dat = [perfile[rawstat][i] for i in ["ftotal", "fcutfound", "fmatched"]]
        outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
            format(*[rawstat]+[str(i) for i in dat]))
        if paired:
            rawstat2 = rawstat.replace("_R1_", "_R2_")
            outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'.\
                format(*[rawstat2]+[str(i) for i in dat]))

    ## spacer, how many records for each sample
    outfile.write('\n{:<35}  {:>13}\n'.\
                  format("sample_name", "total_R1_reads"))

    ## names alphabetical
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
        if paired:
            sample.files["fastq"] = (os.path.join(data.dirs.fastqs,
                                                  name+"_R1_.gz"),
                                     os.path.join(data.dirs.fastqs,
                                                  name+"_R2_.gz"))
        else:
            sample.files["fastq"] = (os.path.join(data.dirs.fastqs,
                                                  name+"_R1_.gz"),)
        sample.stats["reads_raw"] = fsamplehits[name]
        if sample.stats["reads_raw"]:
            sample.stats.state = 1
            data.samples[sample.name] = sample
        else:
            print("Excluded sample: no data found for", name)




def run(data, preview):
    """ demultiplexes raw fastq files given a barcodes file"""

    ## checks on data before starting
    raws, longbar, cut1, paired = prechecks(data, preview)

    if preview:
        print('raws', raws)
    ## nested structure to prevent abandoned temp files
    try: 
        ## splits up all files into chunks, returns list of list
        ## of chunks names in tuples
        parallel_chunker(data, raws, paired)
        chunkspickles = glob.glob(os.path.join(data.dirs.fastqs,
                                              "chunks_*.pickle"))
        chunkslist = []
        for picklefile in chunkspickles:
            ## put chunk names in list
            with open(picklefile) as pickdata:
                rawfilename, chunks = pickle.load(pickdata)
            chunkslist.append([rawfilename, chunks])
        for picklefile in chunkspickles:            
            os.remove(picklefile)
        if preview:
            print("chunkslist", chunkslist)
        try:
            ## sorts all chunk files into tmp files by barcode
            filenum = 0            
            for rawfilename, chunks in chunkslist:
                for cutter in cut1:
                    if cutter:     
                        ## sort chunks for this list     
                        parallel_sorter(data, rawfilename, chunks, 
                                        cutter, longbar, filenum)
                filenum += 1
                ## combine tmps for ambiguous cuts
                ## somefunc()
            ## collate tmps back into one file
            collate_tmps(data, paired)
            make_stats(data, raws, paired)
        finally:
            ## cleans up tmp files
            tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.gz"))
            for tmpfile in tmpfiles:
                os.remove(tmpfile)
    finally:
        ## cleans up chunk files and stats pickles
        pickles = glob.glob(os.path.join(data.dirs.fastqs, "*.pickle"))
        for pickfile in pickles:
            os.remove(pickfile)
        tmpfiles = glob.glob(os.path.join(
                             tempfile.gettempdir(), 
                             "tmp*"))
        for tmpfile in tmpfiles:
            if os.path.isfile(tmpfile):
                os.remove(tmpfile)



if __name__ == "__main__":
    pass
    #PARAMS = {}
    #run(PARAMS)
