#!/usr/bin/env python2.7

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=F0401

import os
import sys
import gzip
import tempfile
import itertools
import subprocess
import numpy as np
import ipyparallel as ipp
from .demultiplex import blocks
from .rawedit import comp



def cleanup(data, sample):
    """ stats, cleanup, and sample """
    
    ## get clustfile
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")

    ## get depth stats
    infile = gzip.open(sample.files.clusters)
    duo = itertools.izip(*[iter(infile)]*2)
    depth = []
    thisdepth = []
    while 1:
        try:
            itera = duo.next()[0]
        except StopIteration:
            break
        if itera != "//\n":
            thisdepth.append(int(itera.split(";")[-2][5:]))
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
        sample.depths.total = depth
        sample.depths.mjmin = keepmj
        sample.depths.statmin = keepstat

        data.stamp("s3 clustering on "+sample.name)        
    else:
        print("no clusters found for {}".format(sample.name))



def muscle_align2(data, sample, chunk, num):
    """ new """
    ## 
    clusts = gzip.open(chunk).read().split("//\n//\n")
    out = []

    ##
    for clust in clusts:
        stack = []        
        lines = clust.split("\n")
        names = lines[::2]
        seqs = lines[1::2]
        ## append counter to end of names
        names = [j+str(i) for i, j in enumerate(names)]        

        ## don't bother aligning singletons
        if len(names) > 1:
            stringnames = alignfast(data, names[:200], seqs[:200])
            anames, aseqs = sortalign(stringnames)
            ## a dict for names2seqs post alignment and indel check
            somedic = {}
            leftlimit = 0
            for i in range(len(anames)):
                ## filter for max indels
                nindels = aseqs[i].rstrip('-').lstrip('-').count('-')
                if nindels <= data.paramsdict["max_Indels_locus"][0]:
                    somedic[anames[i]] = aseqs[i]

                ## do not allow sequence to the left of the seed
                ## to ensure exclusion of adapter/barcodes in gbs
                if anames[i][-1] == "0":
                    leftlimit = min([i for (i, j) in enumerate(aseqs[i]) \
                                     if j != "-"])

            ## reorder keys by derep number
            keys = somedic.keys()
            keys.sort(key=lambda x: int(x[-1]))
            for key in keys:
                if key[-1] == '-':
                    ## reverse matches should have --- overhang
                    if set(somedic[key][:4]) == {"-"}:
                        stack.append(key+"\n"+somedic[key][leftlimit:])
                    else:
                        pass ## excluded
                else:
                    stack.append(key+"\n"+somedic[key][leftlimit:])
        else:
            if names:
                stack = [names[0]+"\n"+seqs[0]]

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    sys.stderr.write(str(len(out)))
    outfile = gzip.open(chunk)
    outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
    outfile.close()



def muscle_align(data, sample, chunk):
    """ multiple sequence alignment of clusters with muscle. """

    ## iterator to read clust file 2 lines at a time
    infile = gzip.open(chunk).read().strip().split("//\n//\n")
    duo = itertools.izip(*[iter(infile)]*2)

    ## list for storing finished loci and a counter
    out = []

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
                if anames[i][-1] == "*":
                    leftlimit = min([i for (i, j) in enumerate(aseqs[i]) \
                                     if j != "-"])

            ## reorder keys by derep number
            keys = somedic.keys()
            keys.sort(key=lambda x: int(x.split(";")[-1][1:]))
            for key in keys:
                if key[-1] == '-':
                    if set(somedic[key][:4]) == {"-"}:
                        stack.append(key+"\n"+somedic[key][leftlimit:])
                    else:
                        pass ## excluded
                else:
                    stack.append(key+"\n"+somedic[key][leftlimit:])
        else:
            if names:
                stack = [names[0]+"\n"+seqs[0]]

        if stack:
            out.append("\n".join(stack))

        ## write to file after 1000 aligned loci
        outfile = gzip.open(os.path.join(
                             data.dirs.clusts,
                             "tmp_"+sample.name+"_"+str(point)+".gz"))
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

    ## derepfile 
    derepfile = sample.files["edits"].replace(".fasta", ".derep")

    ## vsearch results files
    ufile = os.path.join(data.dirs.clusts, sample.name+".utemp")
    tempfile = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## create an output file to write clusters to        
    clustfile = gzip.open(os.path.join(
                          data.dirs.clusts,
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



def split_among_processors(data, samples, preview, noreverse, force):
    """ pass the samples across N_processors to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param preview: toggle preview printing to stdout
    :param noreverse: toggle revcomp clustering despite datatype default
    :param threaded_view: ipyparallel load_balanced_view client

    :returns: None
    """
    ## start a parallel client
    ipyclient = ipp.Client()

    ## nthreads per job for clustering
    threaded_view = ipyclient.load_balanced_view(
                      targets=ipyclient.ids[::data.paramsdict["N_processors"]])
    tpp = len(threaded_view)

    ## make output folder for clusters  
    data.dirs.clusts = os.path.join(
                          data.dirs.edits,
                         "clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.makedirs(data.dirs.clusts)

    ## sort samples by size so biggest go on first
    samples = sorted(samples, key=lambda x: os.stat(x.files["edits"]).st_size)

    ## submit files and args to queue, for func clustall
    submitted_args = []
    for sample in samples:
        if force:
            submitted_args.append([data, sample, preview, noreverse, tpp])
        else:
            ## if not already clustered 
            if sample.stats.state != 3.5:
                submitted_args.append([data, sample, preview, noreverse, tpp])
            else:
                ## clustered but not aligned
                pass

    ## call to ipp for clustering
    results = threaded_view.map(clustall, submitted_args)
    results.get()    

    ## call to ipp for aligning
    lbview = ipyclient.load_balanced_view()
    for sample in samples:
        multi_muscle_align(data, sample, lbview)

    ## write stats to samples
    for sample in samples:
        cleanup(data, sample)

    ## summarize results to stats file
    data.statsfiles.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    if not os.path.exists(data.statsfiles.s3):
        with open(data.statsfiles.s3, 'w') as outfile:
            outfile.write(""+\
            "{:<20}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}\n""".\
                format("sample", "N_reads", "clusts_tot", 
                       "clusts_kept", "avg.depth.tot", 
                       "avg.depth>mj", "avg.depth>stat"))

    ## append stats to file
    outfile = open(data.statsfiles.s3, 'a+')
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



def derep_and_sort(data, sample, preview, nthreads):
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
            " -threads "+str(nthreads)+\
            " -fasta_width 0"
        subprocess.call(cmd, shell=True, 
                             stderr=subprocess.STDOUT, 
                             stdout=subprocess.PIPE)
    else:
        ## pass
        if preview:
            print('skipping dereplication, derep file already exists')



def cluster(data, sample, preview, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    derephandle = sample.files["edits"].replace(".fasta", ".derep")
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

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



def multi_muscle_align(data, sample, lbview):
    """ Splits the muscle alignment across nthreads processors, each runs on 
    1000 clusters at a time. This is a kludge until I find how to write a 
    better wrapper for muscle. 
    """
    ## split clust.gz file into nthreads*10 bits cluster bits
    tmpnames = []

    try: 
        ## get the number of clusters
        clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
        with gzip.open(clustfile, 'rb') as indata:
            totlen = sum([ind.count(";*\n") for ind in blocks(indata)])
        added = 0

        ## calculate optim
        if totlen < 5000:
            ## split into 10 bits
            optim = 250
        elif totlen < 50000:
            ## split into >50 bits
            optim = 1000
        else:
            ## split into 100 bits        
            optim = 2000

        ## write optim clusters to each tmp file
        inclusts = iter(gzip.open(clustfile, 'rb').read().\
                                    strip().split("//\n//\n"))
        grabchunk = itertools.izip(*[inclusts]*optim)        
        while added < totlen:
            with tempfile.NamedTemporaryFile('w+b', delete=False, 
                                             prefix=sample.name, 
                                             suffix='.ali') as out:
                out.write("//\n//\n".join(grabchunk.next())+"//\n//\n")
            tmpnames.append(out.name)
            print(added, added+optim, totlen)
            added += optim

        submitted_args = []
        for num, fname in enumerate(tmpnames):
            submitted_args.append([data, sample, fname, num])

        ## run muscle on all tmp files            
        lbview.map(muscle_align2, submitted_args)

        ## concatenate finished reads
        sample.files.clusters = clustfile.replace("clust.gz", "clustS.gz")
        with gzip.open(sample.files.clusters, 'wb') as out:
            for fname in tmpnames:
                with open(fname) as infile:
                    out.write(infile.read())
                    #os.remove(fname)

    finally:
        ## still delete tmpfiles if job was interrupted
        pass
        #for fname in tmpnames:
        #    if os.path.exists(fname):
        #        os.remove(fname)




def clustall(args):
    """ splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them """

    ## get args
    data, sample, preview, noreverse, nthreads = args

    ## preview
    if preview:
        print("preview: in run_full, using", nthreads)

    ## derep and sort reads by their size
    derep_and_sort(data, sample, preview, nthreads)

    ## cluster_vsearch
    cluster(data, sample, preview, noreverse, nthreads)

    ## cluster_rebuild
    build_clusters(data, sample)



def run(data, samples, preview, noreverse, force):
    """ run the major functions for clustering within samples """

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for _, sample in samples:
        if not force:
            if not sample.files.isnull().clusters:
                print("{} aleady clustered. Sample.stats.state == {}".\
                      format(sample.name, int(sample.stats['state'])))
            else:
                subsamples.append(sample)
        
        else:
            ## clean up existing files from this sample
            if not sample.files.isnull().clusters:
                if os.path.exists(sample.files.clusters):
                    os.remove(sample.files.clusters)
            subsamples.append(sample)

    ## run subsamples 
    if subsamples:
        args = [data, subsamples, preview, noreverse, force]
        split_among_processors(*args)
    else:
        print("\nNo samples found in state 2. To rewrite existing data use\n"+\
              "force=True, or change Sample states to 2.\n")




if __name__ == "__main__":
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
    pass
    
