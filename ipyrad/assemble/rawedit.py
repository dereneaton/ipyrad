#!/usr/bin/env ipython2

""" edits reads based on quality scores. Can be used
to check for adapters and primers, but is not optimized 
for all types of cutters """

from __future__ import print_function
# pylint: disable=E1101
import os
import glob
import gzip
import itertools
import numpy as np
from .demultiplex import ambigcutters
from .demultiplex import chunker  #, blocks




def afilter(data, sample, bases, cut1, cut2, strict, preview, read):
    """ applies filter for primers & adapters """

    ## if pear merged, then read2 is revcomp again
    if sample.pear:
        if read == 2:
            bases = comp(bases)[::-1]
            #print('pear:', sample.pear)

    ## which cutter will be next to adapters?
    if "ddrad" in data.paramsdict["datatype"]:
        cuts = [comp(i)[::-1] for i in cut2]
    else:
        cuts = [comp(i)[::-1] for i in cut1]

    # ## if second read
    # if read == 2:
    #     cuta, cutb = [comp(i)[::-1] for i in cuts]
    # else:
    cuta, cutb = cuts

    ## start looking for cut+[beginning of adapter]
    check1 = check2 = where = None
    if cuts[1]:
        ## if cut sites are ambiguous (what a pain...)        
        if strict == 2:
            lookfor1 = cuta+"A"
            lookfor2 = cutb+"A"
        else:
            lookfor1 = cuta+"AGA"
            lookfor2 = cutb+"AGA"
        if lookfor1 in bases:
            check1 = bases.rindex(lookfor1)
        if lookfor2 in bases:
            check2 = bases.rindex(lookfor2)
        if check1 or check2:
            where = min([i for i in [check1, check2] if i])
    else:
        if strict == 2:
            lookfor = cuta+"A"
        else:
            lookfor = cuta+"AGA"
        if lookfor in bases:
            where = bases.rindex(lookfor)
                
    if not where:
        ## look for adapter sequence directly"
        if strict == 2:
            lookfor1 = "AGATCG"
        else:
            lookfor1 = "AGATCGGA"
        if lookfor1 in bases:
            if read == 1:
                where = bases.rindex(lookfor1) - (len(cuta) + 1)
            else:
                where = bases.rindex(lookfor1) - (len(cuta) + 6)

    ## look for CUT at end of seq
    if not where:
        if cuta in bases[-len(cuta)-5:]:
            where = bases.rindex(cuta)

    ## if preview: show what you've been doing
    if preview:
        print(lookfor, where)
    return where



def comp(seq):
    """ returns a seq with small complement"""
    return seq.replace("A", 't')\
           .replace('T', 'a')\
           .replace('C', 'g')\
           .replace('G', 'c')\
           .replace('n', 'Z')\
           .upper().replace("Z", "n")



def rawedit(args):
    """ three functions:
    (1) replaces low quality base calls with Ns,
    (2) checks for adapter sequence if strict set to 1 or 2 """

    ## get args
    data, sample, tmptuple, paired, preview, point = args

    ## get cut sites
    cut1, cut2 = [ambigcutters(i) for i in \
                  data.paramsdict["restriction_overhang"]]

    ## the read1 demultiplexed reads file
    if tmptuple[0].endswith(".gz"):
        fr1 = gzip.open(tmptuple[0], 'rb')
    else:
        fr1 = open(tmptuple[0], 'rb')
    ## the read2 demultiplexed reads file, if paired
    if paired:
        if tmptuple[1].endswith(".gz"):
            fr2 = gzip.open(tmptuple[1], 'rb')
        else:
            fr2 = open(tmptuple[1], 'rb')

    ## create iterators to sample 4 lines at a time 
    quart1 = itertools.izip(*[iter(fr1)]*4)
    if paired:
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
    writing = []

    ## apply filters to each read (pair)
    while 1:
        try: 
            quart = quarts.next()
        except StopIteration: 
            break

        ## strip reads
        counts['orig'] += 1 
        read1 = [i.strip() for i in quart[0]]
        if paired:
            read2 = [i.strip() for i in quart[1]]
        else:
            read2 = []

        ## filter for Ns
        bases1, bases2 = Nfilter(read1, read2, data, paired, cut1, cut2)

        ## if passed the quality filter, base1 is only returned if both 
        ## pass for paired-end reads
        if bases1:
            ## trim edges if merge data
            if "merge" in data.paramsdict["datatype"]:
                bases1 = trim_merge_edge(bases1, cut1, cut2)
            ## for revcomp 
            if data.paramsdict["datatype"] == "pairgbs":
                #bases1 += "N"*len(sample.barcode)
                bases2 = bases2[-len(bases1):]
                #bases2 = bases2[:len(bases1)]
                #bases2 = "N"*(max([len(i) for i in data.barcodes]) -\
                #               len(sample.barcode))+bases2

            ## look for adapters in sequences
            filternum = data.paramsdict["filter_adapters"]
            wheretocut1 = wheretocut2 = None

            if filternum:
                ## get trim length for read1
                wheretocut1 = afilter(data, sample, bases1, cut1, cut2,
                                      filternum, preview, read=1)
                wheretocut2 = afilter(data, sample, bases2, cut1, cut2,
                                      filternum, preview, read=2)

            ## cut one or both reads depending on detection of adapters
            if wheretocut1 and wheretocut2:
                cutter = min(wheretocut1, wheretocut2)
            elif wheretocut1 or wheretocut2:
                cutter = max(wheretocut1, wheretocut2)
            else:
                cutter = 0

            ## if strict filter, do additional search for cut site near edges
            if filternum == 2:
                if not cutter:
                    if (cut1 in bases1[-15:]) or \
                       (cut2 in bases2[-10:]):
                        cutter = len(bases1)-15

            ## if the read needs to be trimmed
            if cutter:
                if preview:
                    print('cutter', cutter)
                ## if trimmed frag is still long enough
                if cutter > max(32, data.paramsdict["filter_min_trim_len"]):
                    ## data trimmed

                    ## write fasta format
                    if 'pair' in data.paramsdict["datatype"]:
                        sseq = ">"+sample.name+"_"+str(point)+\
                               "_trimpair\n"+bases1[:cutter]+\
                               "SS"+bases2[:cutter]+"\n"
                    else:
                        sseq = ">"+sample.name+"_"+str(point)+\
                               "_c1\n"+bases1[:cutter]+"\n"
                    writing.append(sseq)
                else:
                    ## adapter filtered
                    counts["keep"] += 1

            else:
                ## no cutter
                counts["keep"] += 1
                if 'pair' in data.paramsdict["datatype"]:
                    sseq = ">"+sample.name+"_"+str(point)+\
                           "_pair\n"+bases1+"SS"+\
                           bases2+"\n"
                else:
                    sseq = ">"+sample.name+"_"+str(point)+\
                           "_r1\n"+bases1+"\n"
                writing.append(sseq)
        else:
            ## N filtered out
            counts["quality"] += 1
        ## advance number for read names
        point += 1

    ## write to file
    outdir = os.path.join(data.paramsdict["working_directory"], "edits")
    handle = os.path.join(outdir, "tmp_"+sample.name+"_"+str(point)+".gz")

    ## close file handles
    fr1.close()
    if paired:
        fr2.close()

    if preview:
        print("".join(writing[:20]))

    with gzip.open(handle, 'wb') as out:
        out.write("".join(writing))
    return counts



def Nfilter(read1, read2, data, paired, cut1, cut2):
    """ Filters based on max allowed Ns, and either replaces the Ns, 
        or leaves them. Right now replaces. Also fixes cut site to
        be error free. """

    ## max number of low quality bases allowed
    max_n = data.paramsdict["max_low_qual_bases"]

    ## get qscores
    qscore = np.array([ord(i) for i in read1[3].strip()])
    phred = qscore - data.paramsdict["phred_Qscore_offset"]

    ## replace low qual with N
    bases1 = np.array([i for i in read1[1]])
    bases1[phred < 20] = "N"

    ## fix cut sites to be error free and not carry ambiguities
    bases1[:len(cut1[0])] = list(cut1[0])

    ## count Ns
    bases1 = bases1.tostring()
    bases1_ns = bases1.count("N")

    ## return passed 
    if not paired:
        if bases1_ns < max_n:
            return bases1, False
        else:
            return False, False

    else:
        ## get qscores
        qscore = np.array([ord(i) for i in read2[3].strip('\n')])
        phred = qscore - data.paramsdict["phred_Qscore_offset"]

        ## replace low qual with N
        bases2 = np.array([i for i in read2[1]])
        bases2[phred < 20] = "N"

        ## fix cut sites to be error free before counting Ns
        if 'gbs' in data.paramsdict["datatype"]:
            ## single digest use overhang for enzyme 1
            bases2[:len(cut1[0])] = list(cut1[0])
        else:
            ## double digest use overhang for enzyme 2
            bases2[:len(cut2[0])] = list(cut2[0])

        ## reverse complement as string
        bases2 = comp(bases2.tostring())[::-1]

        ## count Ns
        #bases2 = bases2.tostring()    ## comp return string
        bases2_ns = bases2.count("N")

        ## return passed 
        if (bases1_ns < max_n) and (bases2_ns < max_n):
            return bases1, bases2
        else:
            return False, False



def trim_merge_edge(bases, cut1, cut2):
    """ trim off farside cutsite in merged reads --
    to replace would require knowing how the cut site is cut """

    ## check if two cutters
    if cut2[0]:
        bases = bases[:-(len(cut2[0])+1)]
    else:
        bases = bases[:-(len(cut1[0])+1)]
    return bases



def prechecks(data, preview):
    """ checks before starting analysis """
    ## create output directories 
    #statsdir = os.path.join(data.paramsdict["working_directory"], 'stats')
    #if not os.path.exists(statsdir):
    #    os.makedirs(statsdir)
    data.dirs.edits = os.path.join(data.paramsdict["working_directory"], 
                                  'edits')
    if not os.path.exists(data.dirs.edits):
        os.makedirs(data.dirs.edits)
    ## preview
    if preview:
        print("preview")



def run_full(data, sample, ipyclient, preview):
    """ splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them """

    ## before starting
    prechecks(data, preview)

    ## load up work queue
    submitted = 0
    num = 0

    ## is paired?
    paired = bool("pair" in data.paramsdict["datatype"])

    ## set optim size
    optim = 4000
    if sample.stats.reads_raw:
        if sample.stats.reads_raw > 1e5:
            optim = int(4e4)
        if sample.stats.reads_raw > 1e6:
            optim = int(1e5)
        if sample.stats.reads_raw > 5e6:
            optim = int(4e5)

    ## break up the file into smaller tmp files for each processor
    args = [data, sample.files["fastq"], paired, num, optim, 0]
    _, chunkslist = chunker(args)

    ## 
    try: 
        ## send file to multiprocess queue, will delete if fail
        submitted_args = []
        for tmptuple in chunkslist:
            ## used to increment names across processors
            point = num*optim #10000 #num*(chunksize/2)
            submitted_args.append([data, sample, tmptuple, 
                                    paired, preview, point])
            #work_queue.put([data, sample, tmptuple, paired, preview, point])
            submitted += 1
            num += 1

        ## call to ipp
        dview = ipyclient.load_balanced_view()
        results = dview.map_async(rawedit, submitted_args)
        results.get()
        del dview
    
    finally:
        ## if process failed at any point delete temp files
        for tmptuple in chunkslist:
            os.remove(tmptuple[0])
            os.remove(tmptuple[1]) 
    return submitted, results



def cleanup(data, sample, submitted, results):
    """ cleaning up """

    ## rejoin chunks
    combs = glob.glob(os.path.join(
                        data.dirs.edits,
                        "tmp_"+sample.name+"_*.gz"))

    ## one outfile to write to
    editout = os.path.join(data.dirs.edits,
                           sample.name+".fasta")
    combs.sort(key=lambda x: int(x.split("_")[-1].replace(".gz", "")))
    with open(editout, 'wb') as out:
        for fname in combs:
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
    if not os.path.exists(data.statsfiles.s2):
        with open(data.statsfiles.s2, 'w') as outfile:
            outfile.write('{:<35}  {:>13} {:>13} {:>13} {:>13}\n'.\
                format("sample", "Nreads_orig", "-qscore", 
                       "-adapters", "Nreads_kept"))

    ## append stats to file
    outfile = open(data.statsfiles.s2, 'a+')
    outfile.write('{:<35}  {:>13} {:>13} {:>13} {:>13}\n'.\
                  format(sample.name, 
                         str(fcounts["orig"]),
                         str(fcounts["quality"]),
                         str(fcounts["adapter"]), 
                         str(fcounts["keep"])))
    outfile.close()

    ## save stats to Sample if successful
    sample.stats.state = 2
    sample.files.edits = editout
    sample.stats.reads_filtered = fcounts["keep"]
    ## save stats to the sample??
    data.stamp("s2 rawediting on "+sample.name)        



def run(data, sample, ipyclient, preview=0, force=False):
    """ run the major functions for editing raw reads """
    ## if sample is already done skip
    if not force:
        if sample.stats.state >= 2:
            print("skipping {}. Already edited. Use force=True to overwrite"\
                  .format(sample.name))
        elif sample.stats.reads_raw < 1000:
            print("skipping {}. Too few reads ({}). Use force=True \
                  to override".format(sample.name, sample.stats.reads_raw))
        else:
            submitted, results = run_full(data, sample, ipyclient, preview)
            cleanup(data, sample, submitted, results)
    else:
        submitted, results = run_full(data, sample, ipyclient, preview)
        cleanup(data, sample, submitted, results)


if __name__ == "__main__":
    pass
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #main(PARAMS, FASTQS, QUIET)
