#!/usr/bin/env python2.7

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

from __future__ import print_function
import multiprocessing
import os
import sys
import itertools
import subprocess
import gzip
import numpy as np
from .rawedit import comp
from ipyrad.assemble import worker

# pylint: disable=E1101


def cleanup(data, sample):
    """ stats, cleanup, and sample """
    
    ## grab file
    clustdir = os.path.join(data.paramsdict["working_directory"],
        "clust_"+str(data.paramsdict["clust_threshold"]))
    clusthandle = os.path.join(clustdir, sample.name+".clustS.gz")

    ## get depth stats
    infile = gzip.open(clusthandle)
    duo = itertools.izip(*[iter(infile)]*2)
    depth = []
    thisdepth = []
    while 1:
        try:
            itera = duo.next()[0]
        except StopIteration:
            break
        if itera != "//\n":
            thisdepth.append(int(itera.split(";")[1][5:]))
        else:
            ## append and reset
            depth.append(sum(thisdepth))
            thisdepth = []
    infile.close()

    if depth:
        ## make sense of stats
        depth = np.array(depth)
        keepmj = depth[depth >= data.paramsdict["mindepth_majrule"]]    
        keepstat = depth[depth >= data.paramsdict["mindepth_statistical"]]

        ## sample assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = len(depth)
        sample.stats["clusters_kept"] = max([len(i) for i in \
                                             (keepmj, keepstat)])
        sample.files["clusters"] = str(clusthandle)
        sample.depths["total"] = depth
        sample.depths["mjmin"] = keepmj
        sample.depths["statmin"] = keepstat

        data.stamp("s3 clustering on "+sample.name)        
    else:
        print("no clusters found for {}".format(sample.name))



def muscle_align(data, sample):
    """ multiple sequence alignment of cluster with muscle. 
    If paired, first and second reads are split and aligned 
    separately """
    ## iterator to read clust file 2 lines at a time
    clustdir = os.path.join(data.paramsdict["working_directory"],
        "clust_"+str(data.paramsdict["clust_threshold"]))
    clusthandle = os.path.join(clustdir, sample.name+".clust.gz")
    infile = gzip.open(clusthandle)
    duo = itertools.izip(*[iter(infile)]*2)

    ## lists for storing first and second aligned loci
    out = []
    cnts = 0

    ## remove clustS file if it already exists
    #clustshandle = clusthandle.replace(".clust.gz", ".clustS.gz")
    #if os.path.exists(clustshandle):
    #    os.remove(clustshandle)

    ## iterate over loci
    while 1:
        ## reset lists
        names = []
        seqs = []
        stack = []
        nameiter = 0
        try: 
            itera = duo.next()
        except StopIteration:
            break
        ## read in all data for this stack
        while itera[0] != "//\n":
            names.append(itera[0].strip()+str(nameiter))
            seqs.append(itera[1].strip())
            itera = duo.next()
            nameiter += 1
        ## if not singleton then align locus
        if len(names) > 1:
            ## separate first and second reads
            stringnames = alignfast(data, names[:200], seqs[:200])
            anames, aseqs = sortalign(stringnames)
            ## a dict for names2seqs post alignment and indel check
            somedic = {}
            leftlimit = 0
            for i in range(len(anames)):
                ## if not too many indels
                if aseqs[i].rstrip('-').lstrip('-').count('-') \
                        <= data.paramsdict["max_Indels_locus"][0]:
                    somedic[anames[i]] = aseqs[i]

                ## do not allow sequence to the left of the seed
                ## to ensure exclusion of adapter/barcodes in gbs
                if anames[i].split(";")[-1][0] == "*":
                    leftlimit = min([aseqs[i].index(j) for j in \
                        aseqs[i] if j != "-"])

            ## reorder keys by derep number
            keys = somedic.keys()
            keys.sort(key=lambda x: int(x.split(";")[-1][1:]))
            for key in keys:
                if key[-1] == '-':
                    if set(somedic[key][:4]) == {"-"}:
                        stack.append(key+"\n"+somedic[key][leftlimit:])
                else:
                    stack.append(key+"\n"+somedic[key][leftlimit:])
        else:
            if names:
                stack = [names[0]+"\n"+seqs[0]]

        if stack:
            out.append("\n".join(stack))

        cnts += 1
        ## only write to file after 5000 aligned loci
        if not cnts % 5000:
            if out:
                outfile = gzip.open(
                    clusthandle.replace(".clust", ".clustS"), 'a')
                outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
                outfile.close()
            out = []

    outfile = gzip.open(clusthandle.replace(".clust", ".clustS"), 'a')
    if out:
        outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
    outfile.close()



def sortalign(stringnames):
    """ parses muscle output from a string to two lists """
    objs = stringnames.split("\n>")
    seqs = [i.split("\n")[0].replace(">", "")+"\n"+\
              "".join(i.split('\n')[1:]) for i in objs]
    aligned = [i.split("\n") for i in seqs]
    newnames = [">"+i[0] for i in aligned]
    seqs = [i[1] for i in aligned]     
    return newnames, seqs



def alignfast(data, names, seqs):
    """ makes subprocess call to muscle """
    inputstring = "\n".join(">"+i+"\n"+j for i, j in zip(names, seqs))
    cmd = "/bin/echo '"+inputstring+"' | "+data.muscle+" -quiet -in -"
    piped = subprocess.Popen(cmd, shell=True, 
                       stdin=subprocess.PIPE,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT,
                       close_fds=True)
    _, fout = piped.stdin, piped.stdout
    return fout.read()



def build_clusters(data, sample):
    """ combines information from .u and ._temp files 
    to create .clust files, which contain un-aligned clusters """

    ## find files for this sample
    clustdir = os.path.join(
        data.paramsdict["working_directory"], 
        "clust_"+str(data.paramsdict["clust_threshold"]))

    ## derepfile 
    derepfile = sample.files["edits"].replace(".fasta", ".derep")

    ## vsearch results files
    ufile = os.path.join(clustdir, sample.name+".utemp")
    tempfile = os.path.join(clustdir, sample.name+".htemp")

    ## create an output file to write clusters to        
    clustfile = gzip.open(os.path.join(clustdir, 
                          sample.name+".clust.gz"), 'wb')
    sample.files["clusts"] = clustfile

    ## if .u files are present read them in as userout
    if os.path.exists(ufile):
        userout = open(ufile, 'rb').readlines()
    else:
        userout = []
        print("\n\tSkipping: no '.utemp' file. No clust matches")

    ## load reads into a dictionary"
    hits = {}  

    ## iter 2 lines at a time
    ioderep = open(derepfile, 'rb')
    dereps = itertools.izip(*[iter(ioderep)]*2)
    for namestr, seq in dereps:
        _, count = namestr.strip()[:-1].split(";size=")
        hits[namestr] = [int(count), seq.strip()]
    ioderep.close()

    ## create dictionary of .utemp cluster hits
    udic = {} 
    for uline in userout:
        uitems = uline.strip().split("\t")
        ## if the seed is in udic
        if ">"+uitems[1]+"\n" in udic:
            ## append hit to udict
            udic[">"+uitems[1]+'\n'].append([">"+uitems[0]+"\n", uitems[4],
                                            uitems[5].strip(),
                                            uitems[3]])
        else:
            udic[">"+uitems[1]+'\n'] = [[">"+uitems[0]+"\n", uitems[4],
                                       uitems[5].strip(),
                                       uitems[3]]]

    ## map s equences to clust file in order
    seqslist = [] 
    for key, values in udic.iteritems():
        seq = [key.strip()+"*\n"+hits[key][1]]

        ## allow only 6 internal indels in hits to seed for within-sample clust
        if not any([int(i[3]) > 6 for i in values]):
            for i in xrange(len(values)):
                if values[i][1] == "+":
                    ## only match forward reads if high Cov
                    #if cov[i] >= 90:
                    seq.append(values[i][0].strip()+"+\n"+\
                               hits[values[i][0]][1])
                else:
                    ## name change for reverse hits
                    ## allow low Cov for reverse hits
                    ## .replace("_r1;", "_c1;")+\
                    ## reverse matches should have left terminal gaps
                    ## if not it is a bad alignment and exclude
                    revseq = comp(hits[values[i][0]][1][::-1])
                    seq.append(values[i][0].strip()+"-\n"+revseq)
        seqslist.append("\n".join(seq))
    clustfile.write("\n//\n//\n".join(seqslist)+"\n")

    ## make Dict. from seeds (_temp files) 
    iotemp = open(tempfile, 'rb')
    invars = itertools.izip(*[iter(iotemp)]*2)
    seedsdic = {k:v for (k, v) in invars}  
    iotemp.close()

    ## create a set for keys in I not in seedsdic
    set1 = set(seedsdic.keys())   ## htemp file (no hits) seeds
    set2 = set(udic.keys())       ## utemp file (with hits) seeds
    diff = set1.difference(set2)  ## seeds in 'temp not matched to in 'u
    if diff:
        for i in list(diff):
            clustfile.write("//\n//\n"+i.strip()+"*\n"+hits[i][1]+'\n')
    clustfile.write("//\n//\n")
    clustfile.close()
    del dereps
    del userout
    del udic
    del seedsdic



def split_among_processors(data, samples, preview, noreverse, nthreads):
    """ pass the samples across N_processors to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param preview: toggle preview printing to stdout
    :param override_reverse: toggle revcomp clustering despite datatype default
    :param nthreads: default 4 per parallel job. Set to override default. 

    :returns: None
    """
    ## make output folder for clusters  
    clustdir = os.path.join(
                  data.paramsdict["working_directory"],
                  "clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(clustdir):
        os.makedirs(clustdir)

    ## check for stats dir
    statsdir = os.path.join(
                  data.paramsdict["working_directory"], "stats")
    if not os.path.exists(statsdir):
        os.makedirs(statsdir)

    ## queue items and counters for multiprocessing
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    ## sort samples by size so biggest go on first
    samples = sorted(samples, key=lambda x: os.stat(x.files["edits"]).st_size)

    ## submit files and args to queue, for func run_all
    for sample in samples:
        work_queue.put([data, sample, preview, noreverse, nthreads])

    ## make a job list and start workers on jobs
    jobs = []
    for _ in xrange(data.paramsdict["N_processors"]):
        work = worker.Worker(work_queue, result_queue, runprocess)
        jobs.append(work)
        work.start()
    for j in jobs:
        j.join()

    ## 
    if preview:
        print("finished clustering")
    
    ## write stats to samples
    for sample in samples:
        cleanup(data, sample)

    ## summarize results to stats file
    statsdir = os.path.join(data.paramsdict["working_directory"], "stats")
    if not os.path.exists(statsdir):
        os.mkdir(statsdir)
    statsfile = os.path.join(statsdir, 's3_cluster_stats.txt')
    if not os.path.exists(statsfile):
        with open(statsfile, 'w') as outfile:
            outfile.write(""+\
            "{:<20}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}\n""".\
                format("sample", "N_reads", "clusts_tot", 
                       "clusts_kept", "avg.depth.tot", 
                       "avg.depth>mj", "avg.depth>stat"))

    ## append stats to file
    outfile = open(statsfile, 'a+')
    samples.sort(key=lambda x: x.name)
    for sample in samples:
        outfile.write(""+\
    "{:<20}   {:>9}   {:>10}   {:>11}   {:>13.2f}   {:>11.2f}   {:>13.2f}\n".\
                      format(sample.name, 
                             int(sample.stats["reads_filtered"]),
                             int(sample.stats["clusters_total"]),
                             int(sample.stats["clusters_kept"]),
                             np.mean(sample.depths["total"]),
                             np.mean(sample.depths["mjmin"]),
                             np.mean(sample.depths["statmin"]),
                             ))
    outfile.close()
    ## link statsfile to dataobj
    data.statsfiles["s3"] = statsfile




def derep_and_sort(data, sample, preview):
    """ dereplicates reads and write to .step file
    ...sorts dereplicated file (.step) so reads that were highly
    replicated are at the top, and singletons at bottom, writes
    output to .derep file """

    ## select file
    handle = sample.files["edits"]

    ## reverse complement clustering for some types    
    if data.paramsdict["datatype"] in ['pairgbs', 'gbs', 'merged']:
        reverse = " -strand both "
    else:
        reverse = " "

    ## only make derep if .derep doesn't already exist
    if not os.path.exists(handle.replace(".fasta", ".derep")):
        if preview:
            print("dereplicating...")
        ## do dereplication with vsearch
        cmd = data.vsearch+\
            " -derep_fulllength "+handle+\
            reverse+\
            " -output "+handle.replace(".fasta", ".derep")+\
            " -sizeout "+\
            " -threads 4"+\
            " -fasta_width 0"
        subprocess.call(cmd, shell=True, 
                             stderr=subprocess.STDOUT, 
                             stdout=subprocess.PIPE)

        ## do sorting with vsearch
        #cmd = data.paramsdict["vsearch_path"]+\
        #      " -sortbysize "+handle.replace(".fasta", ".step")+\
        #      " -output "+handle.replace(".fasta", ".derep")
        #subprocess.call(cmd, shell=True, 
        #                     stderr=subprocess.STDOUT,
        #                     stdout=subprocess.PIPE)
        ## remove step files
        #if os.path.exists(handle.replace(".fasta", ".step")):
        #    os.remove(handle.replace(".fasta", ".step"))
    else:
        ## pass
        if preview:
            print('skipping dereplication, derep file already exists')



def cluster(data, sample, preview, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    derephandle = sample.files["edits"].replace(".fasta", ".derep")
    clustdir = os.path.join(data.paramsdict["working_directory"],
        "clust_"+str(data.paramsdict["clust_threshold"]))
    uhandle = os.path.join(clustdir, sample.name+".utemp")
    temphandle = os.path.join(clustdir, sample.name+".htemp")

    ## datatype variables
    if data.paramsdict["datatype"] in ['gbs']:
        reverse = " -strand both "
        cov = " -query_cov .35 " 
    elif data.paramsdict["datatype"] in ['pairgbs', 'merged']:
        reverse = "  -strand both "
        cov = " -query_cov .60 " 
    else:  ## rad, ddrad, ddradmerge
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    ## get call string
    cmd = data.vsearch+\
        " -cluster_smallmem "+derephandle+\
        reverse+\
        cov+\
        " -id "+str(data.paramsdict["clust_threshold"])+\
        " -userout "+uhandle+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -minsl 0.5"+\
        " -fulldp"+\
        " -threads "+str(nthreads)+\
        " -usersort "+\
        " -notmatched "+temphandle+\
        " -fasta_width 0"

    #print(cmd)

    ## run vsearch
    if preview:
        ## make this some kind of wait command that kills after a few mins
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    else:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)



def multi_muscle_align(data, sample, nthreads):
    """ Splits the muscle alignment across four processors, since it 
    isn't able to run above ~.25 due to all the I/O. This is a kludge 
    until I find how to write a better (faster) wrapper for muscle...
    Split parts are unequal in size since the first contains most of the 
    large alignments and the later chunks include mostly singletons that 
    do not need aligning. 
    """
    muscle_align(data, sample)



def runprocess(data, sample, preview, noreverse, nthreads):
    """ splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them """
    ## preview
    if preview:
        print("preview: in run_full")

    ## derep and sort reads by their size
    derep_and_sort(data, sample, preview)

    ## cluster_vsearch
    cluster(data, sample, preview, noreverse, nthreads)

    ## cluster_rebuild
    build_clusters(data, sample)

    ## align clusters on single thread
    muscle_align(data, sample)

    ## align clusters across multiple threads
    #multi_muscle_align(data, sample, nthreads)


def run(data, samples, preview=0, noreverse=0, nthreads=4):
    """ run the major functions for editing raw reads """
    ## print to screen?
    # if not quiet:
    # 	toscreen = sys.stdout
    # else:
    # 	toscreen = "dev/null"
    #print(samples)

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for _, sample in samples:
        if sample.stats['state'] >= 3:
            print("{} aleady clustered. Sample.stats['state'] == {}".\
                  format(sample.name, int(sample.stats['state'])))
        else:
            if sample.stats['state'] == 2:
                subsamples.append(sample)
            else:
                print("skipping,", sample.name,
                      "not yet edited. Sample.stats['state'] ==", 
                      str(sample.stats['state']))

    ## run subsamples 
    if subsamples:
        split_among_processors(data, subsamples, preview, noreverse, nthreads)
    else:
        print(
"""
No samples found in state 2. To reset states do:

for samp in data.samples.values():
    samp.stats['state'] = 2
""")




if __name__ == "__main__":
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
    pass
    
