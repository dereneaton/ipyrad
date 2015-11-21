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
#from .demultiplex import blocks
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

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["reference_sequence"] == "":
        cmd = data.samtools+\
            " flagstat "+sample.files.unmapped_reads
        result = subprocess.check_output( cmd, shell=True,
                                           stderr=subprocess.STDOUT )
        sample.stats["refseq_unmapped_reads"] = int(result.split()[0])

        cmd = data.samtools+\
            " flagstat "+sample.files.mapped_reads
        result = subprocess.check_output( cmd, shell=True,
                                           stderr=subprocess.STDOUT )
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
    """ pass the samples across N_processors to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param preview: toggle preview printing to stdout
    :param noreverse: toggle revcomp clustering despite datatype default
    :param threaded_view: ipyparallel load_balanced_view client

    :returns: None
    """
    ## nthreads per job for clustering
    threaded_view = ipyclient.load_balanced_view(
                    targets=ipyclient.ids[::data.paramsdict["N_processors"]])
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
    if not data.paramsdict["reference_sequence"] == "":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                                os.path.realpath(data.paramsdict["working_directory"]),
                                data.name+"_"+
                                "refmapping")
        if not os.path.exists(data.dirs.refmapping):
            os.makedirs(data.dirs.refmapping)
        ## Set the mapped and unmapped reads files for this sample
        for sample in samples:
            sorted_unmapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-sorted-unmapped.bam")
            sorted_mapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-sorted-mapped.bam")
            sample.files.unmapped_reads = sorted_unmapped_bamhandle
    	    sample.files.mapped_reads = sorted_mapped_bamhandle
            print(sample.files.mapped_reads)

        ## call to ipp for read mapping
        results = threaded_view.map(mapreads, submitted_args)
        results.get()

    ## call to ipp for clustering
    results = threaded_view.map(clustall, submitted_args)
    results.get()   
    del threaded_view 

    ## call to ipp for aligning
    lbview = ipyclient.load_balanced_view()
    for sample in samples:
        multi_muscle_align(data, sample, lbview)
    del lbview

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

    ## concatenate if multiple edits files for a sample
    if len(sample.files.edits) > 1:
        ## create temporary concat file
        tmphandle = os.path.join(data.dirs.edits,
                              "tmp_"+sample.name+".concat")       
        with open(tmphandle, 'wb') as tmp:
            for editfile in sample.files.edits:
                with open(editfile) as inedit:
                    tmp.write(inedit)
        handle = tmphandle
    else:
        handle = sample.files.edits[0]

    ## reverse complement clustering for some types    
    if data.paramsdict["datatype"] in ['pairgbs', 'gbs', 'merged']:
        reverse = " -strand both "
    else:
        reverse = " "

    ## do dereplication with vsearch
    cmd = data.vsearch+\
        " -derep_fulllength "+handle+\
        reverse+\
        " -output "+os.path.join(data.dirs.edits, sample.name+".derep")+\
        " -sizeout "+\
        " -threads "+str(nthreads)+\
        " -fasta_width 0"
    subprocess.call(cmd, shell=True, 
                         stderr=subprocess.STDOUT, 
                         stdout=subprocess.PIPE)

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
        print("preview: in run_full, using", nthreads)

    ## derep and sort reads by their size
    derep_and_sort(data, sample, preview, nthreads)

    ## cluster_vsearch
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
    unmapped_fastq_handle = os.path.join(data.dirs.edits, sample.name+".fastq")

    ## Stash this one for later, because we'll need access to the successfully
    ## mapped reads downstream
    sample.files.mapped_reads = sorted_mapped_bamhandle

################
# Paired end isn't handled yet, but it needs to be.
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
        " " + data.get_params(27) +\
        " " + sample.files.fastq[0]

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
    ## No args for this one
    ## TODO: output is actually fastq, so fix that
    ## TODO: output file name is fsck, and it also needs to be .gz to kosher with
    ##       the other files for de novo
    cmd = data.samtools+\
        " bam2fq "+sorted_unmapped_bamhandle+\
        " > "+unmapped_fastq_handle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    cmd = data.samtools+\
        " flagstat "+sorted_unmapped_bamhandle
    result = subprocess.check_output( cmd, shell=True,
                                           stderr=subprocess.STDOUT )
    sample.stats.refseq_unmapped_reads = int(result.split()[0])

    cmd = data.samtools+\
        " flagstat "+sorted_mapped_bamhandle
    result = subprocess.check_output( cmd, shell=True,
                                           stderr=subprocess.STDOUT )
    sample.stats.refseq_mapped_reads = int(result.split()[0])

    with open('/tmp/wat.do', 'a') as out:
        out.write(result.split()[0]+" "+sample.name+ " " + sample.files.mapped_reads)

    ## This is hax to get fastq to fasta to get this off the ground.
    ## samtools bam2fq natively returns fastq, you just delete this code
    ## when fastq pipleline is working
    writing = []
    with open(os.path.realpath(unmapped_fastq_handle), 'rb') as fq:
        quart1 = itertools.izip(*[iter(fq)]*4)
        quarts = itertools.izip(quart1, iter(int, 1))
        writing = []
        while 1:
            try:
                quart = quarts.next()
            except StopIteration:
                break
            read1 = [i.strip() for i in quart[0]]
            sseq = ">"+sample.name+"_"+str(0)+\
                           "_c1\n"+read1[1]+"\n"
            writing.append(sseq)

    with open( sample.files.edits[0], 'w' ) as out:
        out.write("".join(writing))
    #####################################################################
    ## End block of dummy code, delete this when fastq works

def run(data, samples, ipyclient, preview, noreverse, force):
    """ run the major functions for clustering within samples """

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for _, sample in samples:
        if not force:
            if sample.files.clusters:
                print("{} aleady clustered. Sample.stats.state == {}".\
                      format(sample.name, int(sample.stats['state'])))
            else:
                subsamples.append(sample)
        
        else:
            ## clean up existing files from this sample and overwrite
            #if sample.files.clusters:
            #    if os.path.exists(sample.files.clusters):
            #        os.remove(sample.files.clusters)
            subsamples.append(sample)

    ## run subsamples 
    if subsamples:
        args = [data, subsamples, ipyclient, preview, noreverse, force]
        split_among_processors(*args)
    else:
        print("\nNo samples found in state 2. To rewrite existing data use\n"+\
              "force=True, or change Sample states to 2.\n")




if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.set_params(27, '/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa')
    DATA.get_params()
    print(DATA.log)
    DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
    pass
    
