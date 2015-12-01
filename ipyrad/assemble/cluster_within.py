#!/usr/bin/env python2.7

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142

import os
import sys
import gzip
import time
import tempfile
import itertools
import subprocess
import numpy as np
import ipyparallel as ipp
from .rawedit import comp

import logging
LOGGER = logging.getLogger(__name__)



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

        data._stamp("s3 clustering on "+sample.name)        
    else:
        print("no clusters found for {}".format(sample.name))

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["assembly_method"] == "denovo":
        cmd = data.samtools+\
            " flagstat "+sample.files.unmapped_reads
        result = subprocess.check_output(cmd, shell=True,
                                              stderr=subprocess.STDOUT)
        sample.stats["refseq_unmapped_reads"] = int(result.split()[0])

        cmd = data.samtools+\
            " flagstat "+sample.files.mapped_reads
        result = subprocess.check_output(cmd, shell=True,
                                              stderr=subprocess.STDOUT)
        sample.stats["refseq_mapped_reads"] = int(result.split()[0])



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
    #piped = subprocess.Popen(cmd, shell=True, 
    #                   stdout=subprocess.PIPE,
    #                   stderr=subprocess.STDOUT,
    #                   close_fds=True)
    #_, fout = piped.stdin, piped.stdout
    return fout.read()



def build_clusters(data, sample):
    """ combines information from .u and ._temp files 
    to create .clust files, which contain un-aligned clusters """

    ## derepfile 
    derepfile = os.path.join(data.dirs.edits, sample.name+".derep")

    ## vsearch results files
    ufile = os.path.join(data.dirs.clusts, sample.name+".utemp")
    htempfile = os.path.join(data.dirs.clusts, sample.name+".htemp")

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
    iotemp = open(htempfile, 'rb')
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




def split_among_processors(data, samples, ipyclient, preview, noreverse, force):
    """ pass the samples to N engines to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param preview: toggle preview printing to stdout
    :param noreverse: toggle revcomp clustering despite datatype default
    :param threaded_view: ipyparallel load_balanced_view client

    :returns: None
    """
    ## nthreads per job for clustering
    threaded_view = ipyclient.load_balanced_view(
                    targets=ipyclient.ids[::data.paramsdict["engines_per_job"]])
    tpp = len(threaded_view)

    ## make output folder for clusters  
    data.dirs.clusts = os.path.join(
                        os.path.realpath(data.paramsdict["working_directory"]),
                        data.name+"_"+
                       "clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.makedirs(data.dirs.clusts)

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

    # If reference sequence is specified then try read mapping, else pass.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                        os.path.realpath(data.paramsdict["working_directory"]),
                        data.name+"_refmapping")
      
        if not os.path.exists(data.dirs.refmapping):
            os.makedirs(data.dirs.refmapping)
        ## Set the mapped and unmapped reads files for this sample
        for sample in samples:
            sorted_unmapped_bamhandle = os.path.join(data.dirs.refmapping,
                                            sample.name+"-sorted-unmapped.bam")
            sorted_mapped_bamhandle = os.path.join(data.dirs.refmapping, 
                                            sample.name+"-sorted-mapped.bam")
            sample.files.unmapped_reads = sorted_unmapped_bamhandle
            sample.files.mapped_reads = sorted_mapped_bamhandle

        ## call to ipp for read mapping
        results = threaded_view.map(mapreads, submitted_args)
        try:
            results.get()
        except (AttributeError, TypeError):
            for key in ipyclient.history:
                if ipyclient.metadata[key].error:
                    LOGGER.error("step3 readmapping error: %s", 
                        ipyclient.metadata[key].error)
                    raise SystemExit("step3 readmapping error.\n({})."\
                                     .format(ipyclient.metadata[key].error))
                if ipyclient.metadata[key].stdout:
                    LOGGER.error("step3 readmapping stdout:%s", 
                        ipyclient.metadata[key].stdout)
                    raise SystemExit("step3 readmapping error.")

            LOGGER.error(ipyclient.metadata)
            sys.exit("")

    ## call to ipp for clustering
    results = threaded_view.map(clustall, submitted_args)
    try:
        results.get()
    except (AttributeError, TypeError):
        for key in ipyclient.history:
            if ipyclient.metadata[key].error:
                LOGGER.error("step3 clustering error: %s", 
                    ipyclient.metadata[key].error)
                raise SystemExit("step3 clustering error.\n({})."\
                                 .format(ipyclient.metadata[key].error))
            if ipyclient.metadata[key].stdout:
                LOGGER.error("step3 clustering stdout:%s", 
                    ipyclient.metadata[key].stdout)
                raise SystemExit("step3 clustering error.")
    del threaded_view 

    ## call to ipp for aligning
    lbview = ipyclient.load_balanced_view()
    for sample in samples:
        multi_muscle_align(data, sample, lbview)
    del lbview

    ## If reference sequence is specified then pull in alignments from 
    ## mapped bam files and write them out to the clustS files to fold
    ## them back into the pipeline.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                        os.path.realpath(data.paramsdict["working_directory"]),
                        data.name+"_refmapping")
        for sample in samples:
            results = getalignedreads(data, sample)

    ## write stats to samples
    for sample in samples:
        cleanup(data, sample)

    ## summarize results to stats file
    ## TODO: update statsfile with mapped and unmapped reads for reference mapping
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



def combine_pairs(data, sample):
    """ in prep """
    LOGGER.info(sample.files.edits)

    ## open file for writing to
    combined = os.path.join(data.dirs.edits, sample.name+"_pairs.fastq")
    combout = open(combined, 'wb')

    ## get nonmerged pairs
    nonmerged1 = os.path.join(data.dirs.edits, 
                              sample.name+"_nonmerged_R1_.fastq")
    nonmerged2 = os.path.join(data.dirs.edits, 
                              sample.name+"_nonmerged_R2_.fastq")

    ## read in paired end read files"
    ## create iterators to sample 4 lines at a time
    fr1 = open(nonmerged1, 'rb')
    quart1 = itertools.izip(*[iter(fr1)]*4)
    fr2 = open(nonmerged2, 'rb')
    quart2 = itertools.izip(*[iter(fr2)]*4)
    quarts = itertools.izip(quart1, quart2)

    ## a list to store until writing
    writing = []
    counts = 0

    ## iterate until done
    while 1:
        try:
            read1s, read2s = quarts.next()
        except StopIteration:
            break
        writing.append("\n".join([
                        read1s[0].strip(),
                        read1s[1].strip()+\
                            "SSSS"+comp(read2s[1].strip())[::-1],
                        read1s[2].strip(),
                        read1s[3].strip()+\
                            "SSSS"+read2s[3].strip()[::-1]]
                        ))
        counts += 1
        if not counts % 1000:
            combout.write("\n".join(writing)+"\n")

    combout.write("\n".join(writing))
    combout.close()

    sample.files.edits = [(combined, )]
    return sample



def merge_fastq_pairs(data, sample):
    """ Merge paired fastq reads. """

    ## tempnames for merge files
    merged = os.path.join(data.dirs.edits,
                          sample.name+"_merged_.fastq")
    nonmerged1 = os.path.join(data.dirs.edits, 
                              sample.name+"_nonmerged_R1_.fastq")
    nonmerged2 = os.path.join(data.dirs.edits, 
                              sample.name+"_nonmerged_R2_.fastq")

    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']

    assert os.path.exists(sample.files.edits[0][1]), \
           "No paired read file (_R2_ file) found."

    ## vsearch merging
    cmd = data.vsearch \
      +" --fastq_mergepairs "+sample.files.edits[0][0] \
      +" --reverse "+sample.files.edits[0][1] \
      +" --fastqout "+merged \
      +" --fastqout_notmerged_fwd "+nonmerged1 \
      +" --fastqout_notmerged_rev "+nonmerged2 \
      +" --fasta_width 0 " \
      +" --fastq_allowmergestagger " \
      +" --fastq_minmergelen 32 " \
      +" --fastq_maxns "+str(maxn)

    try:
        subprocess.check_call(cmd, shell=True,   
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(subprocess.STDOUT)
        LOGGER.error(cmd)        
        sys.exit("Error in merging pairs: \n({}).".format(inst))
    finally:
        sample.merged = 1
    return sample



def concat_edits(data, sample):
    """ concatenate if multiple edits files for a sample """
    ## if more than one tuple in the list
    if len(sample.files.edits) > 1:
        ## create temporary concat file
        tmphandle1 = os.path.join(data.dirs.edits,
                              "tmp1_"+sample.name+".concat")       
        with open(tmphandle1, 'wb') as tmp:
            for editstuple in sample.files.edits:
                with open(editstuple[0]) as inedit:
                    tmp.write(inedit)

        if 'pair' in data.paramsdict['datatype']:
            tmphandle2 = os.path.join(data.dirs.edits,
                                "tmp2_"+sample.name+".concat")       
            with open(tmphandle2, 'wb') as tmp:
                for editstuple in sample.files.edits:
                    with open(editstuple[1]) as inedit:
                        tmp.write(inedit)
            sample.files.edits = [(tmphandle1, tmphandle2)]
        else:
            sample.files.edits = [(tmphandle1, )]
    return sample



def derep_and_sort(data, sample, preview, nthreads):
    """ dereplicates reads and sorts so reads that were highly
    replicated are at the top, and singletons at bottom, writes
    output to .derep file """

    ## reverse complement clustering for some types    
    #if data.paramsdict["datatype"] in ['pairgbs', 'gbs', 'merged']:
    if sample.merged:
        reverse = " -strand both "
    else:
        reverse = " "

    LOGGER.debug("derep FILE %s", sample.files.edits[0][0])

    ## do dereplication with vsearch
    cmd = data.vsearch+\
          " -derep_fulllength "+sample.files.edits[0][0]+\
          reverse+\
          " -output "+os.path.join(data.dirs.edits, sample.name+".derep")+\
          " -sizeout "+\
          " -threads "+str(nthreads)+\
          " -fasta_width 0"

    ## run vsearch
    try:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.info(cmd)
        LOGGER.info(inst)
        sys.exit("Error in vsearch: \n{}\n{}\n{}."\
                 .format(inst, subprocess.STDOUT, cmd))



def cluster(data, sample, preview, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    derephandle = os.path.join(data.dirs.edits, sample.name+".derep")
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
    try:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



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
        optim = 1000

        ## write optim clusters to each tmp file
        inclusts = iter(gzip.open(clustfile, 'rb').read().\
                                  strip().split("//\n//\n"))
        grabchunk = list(itertools.islice(inclusts, optim))
        while grabchunk:
            with tempfile.NamedTemporaryFile('w+b', delete=False, 
                                             dir=data.dirs.clusts,
                                             prefix=sample.name+"_", 
                                             suffix='.ali') as out:
                out.write("//\n//\n".join(grabchunk)+"//\n//\n")
            tmpnames.append(out.name)
            grabchunk = list(itertools.islice(inclusts, optim))
    
        ## create job queue
        submitted_args = []
        for num, fname in enumerate(tmpnames):
            submitted_args.append([data, sample, fname, num])

        ## run muscle on all tmp files            
        lbview.map(muscle_align2, submitted_args)

        ## concatenate finished reads
        sample.files.clusters = os.path.join(data.dirs.clusts,
                                             sample.name+".clustS.gz")
        with gzip.open(sample.files.clusters, 'wb') as out:
            for fname in tmpnames:
                with open(fname) as infile:
                    out.write(infile.read())

    finally:
        ## still delete tmpfiles if job was interrupted
        for fname in tmpnames:
            if os.path.exists(fname):
                os.remove(fname)



def clustall(args):
    """ splits fastq file into smaller chunks and distributes them across
    multiple processors, and runs the rawedit func on them """

    ## get args
    data, sample, preview, noreverse, nthreads = args

    ## preview
    if preview:
        print("preview: in clustall, using", nthreads)

    ## concatenate edits files in case a Sample has multiple, and 
    ## returns a new Sample.files.edits with the concat file
    sample = concat_edits(data, sample)

    ## if reference do ref align here...
    ## ...

    ## merge fastq pairs
    if 'pair' in data.paramsdict['datatype']:
        ## merge pairs that overlap
        sample = merge_fastq_pairs(data, sample)
        ## combine end-to-end remaining pairs for denovo clustering
        sample = combine_pairs(data, sample)

    ## convert fastq to fasta, then derep and sort reads by their size
    derep_and_sort(data, sample, preview, nthreads)

    ## cluster derep fasta files in vsearch 
    cluster(data, sample, preview, noreverse, nthreads)

    ## cluster_rebuild
    build_clusters(data, sample)



def mapreads(args):
    """ Attempt to map reads to reference sequence. This reads in the 
    samples.files.edits .fasta files, and attempts to map each read to the 
    reference sequences. Unmapped reads are dropped right back in the de novo 
    pipeline. Reads that map successfully are processed and pushed downstream 
    and joined with the rest of the data post musle_align. The read mapping 
    produces a sam file, unmapped reads are pulled out of this and dropped back 
    in place of the edits (.fasta) file. The raw edits .fasta file is moved to 
    .<sample>.fasta to hide it in case the mapping screws up and we need to roll-back.
    Mapped reads stay in the sam file and are pulled out of the pileup later."""

    ## TODO: Right now this is coded to read in the edits which are .fasta format
    ## It would be ideal to read in the fastq instead, since this contains info
    ## about base quality scores, improves the mapping confidence.

    ## get args
    data, sample, preview, noreverse, nthreads = args

    ## preview
    if preview:
        print("preview: in run_full, using", nthreads)

    ## Files we'll use during reference sequence mapping
    ##
    ## samhandle - Raw output from smalt reference mapping
    ## unmapped_bamhandle - bam file containing only unmapped reads
    ## mapped_bamhandle - bam file with only successfully mapped reads
    ##    In the case of paired end, we only consider pairs where both
    ##    reads map successfully
    ## sorted_*_bamhandle - Sorted bam files for mapped and unmapped
    ##    This is a precursor to bam2fq which requires sorted bams
    ## unmapped_fastq_handle - The final output file of this function
    ##    writes out unmapped reads to the 'edits' directory as .fq
    ##    which is what downstream analysis expects
    samhandle = os.path.join(data.dirs.refmapping, sample.name+".sam")
    unmapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-unmapped.bam")
    mapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-mapped.bam")
    sorted_unmapped_bamhandle = sample.files["unmapped_reads"]
    sorted_mapped_bamhandle = sample.files["mapped_reads"]

    ## TODO, figure out why edits is a tuple? There could be multiple edits files, yes?
    ## but by the time we get to step three they are all collapsed in to one big file
    ## which is the first element of the tuple, yes?
    ##
    ## TODO: This is hackish, we are overwriting the fastq that contains all reads
    ## might want to preserve the full .fq file in case of a later 'force' command
    unmapped_fastq_handle = sample.files.edits[0][0] #os.path.join(data.dirs.edits, sample.name+".fastq")

################
# TODO: Paired end isn't handled yet, but it needs to be.
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
################


    ## get call string
    cmd = data.smalt+\
        " map -f sam -n " + str(nthreads) +\
        " -o " + samhandle +\
        " " + data.paramsdict['reference_sequence'] +\
        " " + sample.files.edits[0][0]

    ## run smalt
    if preview:
        ## make this some kind of wait command that kills after a few mins
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    else:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)

    ## Get the reads that map successfully. For PE both reads must map
    ## successfully in order to qualify.
    ##   1) Get mapped reads and convert to bam
    ##   2) Sort them and save the path to the bam to sample.files.mapped_reads
    ## The mapped reads are synced back into the pipeline downstream, after
    ## muscle aligns the umapped reads.
    ##
    ## samtools view arguments
    ##   -b = write to .bam
    ##   -F = Select all reads that DON'T have this flag (0x4 means unmapped)
    ##        <TODO>: This is deeply hackish right now, it will need some
    ##                serious thinking to make this work for PE, etc.
    cmd = data.samtools+\
        " view -b -F 0x4 "+samhandle+\
            " > " + mapped_bamhandle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)
    ## Step 2 sort mapped bam
    ##   -T = Temporary file name, this is required by samtools, you can ignore it
    ##        Here we just hack it to be samhandle.tmp cuz samtools will clean it up
    ##   -O = Output file format, in this case bam
    ##   -o = Output file name
    cmd = data.samtools+\
        " sort -T "+samhandle+".tmp" +\
        " -O bam "+mapped_bamhandle+\
        " -o "+sorted_mapped_bamhandle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ##############################################
    ## Do unmapped
    ##############################################

    ## Fetch up the unmapped reads and write them to the edits .fasta file
    ## In order to do this we have to go through a little process:
    ##   1) Get the unmapped reads and convert to bam
    ##   2) Sort them 
    ##   3) Dump bam to fastq

    ## Step 1 get unmapped reads with samtools view
    ##   -b = write out to .bam
    ##   -f = Select only reads with associated flags set (0x4 means unmapped)
    cmd = data.samtools+\
        " view -b -f 0x4 "+samhandle+\
            " > " + unmapped_bamhandle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## Step 2 sort unmapped bam
    ##   -T = Temporary file name, this is required by samtools, you can ignore it
    ##        Here we just hack it to be samhandle.tmp cuz samtools will clean it up
    ##   -O = Output file format, in this case bam
    ##   -o = Output file name
    cmd = data.samtools+\
        " sort -T "+samhandle+".tmp" +\
        " -O bam "+unmapped_bamhandle+\
        " -o "+sorted_unmapped_bamhandle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## Step 3 Dump sorted bam to fastq
    ## No args for this one. Pipe it through gzip, because the downstream
    ## analysis expects gzipped fq
    cmd = data.samtools+\
        " bam2fq "+sorted_unmapped_bamhandle+\
        " | gzip > "+unmapped_fastq_handle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## bs for experimentation. is samtools/smalt fscking the output?
    ## TODO: get rid of this when mapped reads get returned to the 
    ## pipeline right. or whenever, it's not useful.
    cmd = data.samtools+\
        " bam2fq "+sorted_mapped_bamhandle+\
        " | gzip >> "+unmapped_fastq_handle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## This is the end of processing for each sample. Stats
    ## are appended to the sample for mapped and unmapped reads 
    ## during cluster_within.cleanup()

    ## This is hax to get fastq to fasta to get this off the ground.
    ## samtools bam2fq natively returns fastq, you just delete this code
    ## when fastq pipleline is working
#    writing = []
#    with open(os.path.realpath(unmapped_fastq_handle), 'rb') as fq:
#        quart1 = itertools.izip(*[iter(fq)]*4)
#        quarts = itertools.izip(quart1, iter(int, 1))
#        writing = []
#        j=0
#        while 1:
#            try:
#                quart = quarts.next()
#            except StopIteration:
#                break
#            read1 = [i.strip() for i in quart[0]]
#            sseq = ">"+sample.name+"_"+str(j)+\
#                           "_r1\n"+read1[1]+"\n"
#            writing.append(sseq)
#            j = j+1
#
#    with open( sample.files.edits[0], 'w' ) as out:
#        out.write("".join(writing))
    #####################################################################
    ## End block of dummy code, delete this when fastq works


def getalignedreads(data, sample):
    """Pull aligned reads out of sorted mapped bam files and
    append them to the clustS.gz file so the fall into downstream analysis """

    mapped_fastq_handle = "/tmp/wat.fq"
    ## Build the samtools bam2fq command to push bam out
    cmd = data.samtools+\
        " bam2fq "+sample.files["mapped_reads"]+\
        " > "+mapped_fastq_handle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)



def run(data, samples, ipyclient, preview, noreverse, force):
    """ run the major functions for clustering within samples """

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for _, sample in samples:
        if not force:
            if sample.files.clusters:
                print("Skipping {}; aleady clustered.".\
                      format(sample.name))
            else:
                subsamples.append(sample)
        else:
            ## force to overwrite
            subsamples.append(sample)

    ## run subsamples 
    assert subsamples, "No Samples ready to be clustered. To rewrite existing" \
                      +"data use force=True."
    args = [data, subsamples, ipyclient, preview, noreverse, force]
    split_among_processors(*args)




if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.set_params(28, '/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa')
    DATA.get_params()
    print(DATA.log)
    DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
    pass
    
