#!/usr/bin/env python2.7

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142
# pylint: disable=W0212

import os
import sys
import gzip
#import shlex
import socket
import tempfile
import itertools
import subprocess
import numpy as np
from collections import OrderedDict, defaultdict
from .refmap import *
from .util import * #comp, merge_pairs

import logging
LOGGER = logging.getLogger(__name__)


def sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """
    
    ## get clustfile
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")

    if "pair" in data.paramsdict["datatype"]:
        ## record merge file name temporarily
        sample.files.merged = os.path.join(data.dirs.edits,
                                           sample.name+"_merged_.fastq")
        ## record how many read pairs were merged
        with open(sample.files.merged, 'r') as tmpf:
            ## divide by four b/c its fastq
            sample.stats.reads_merged = len(tmpf.readlines())/4

    ## get depth stats
    infile = gzip.open(sample.files.clusters)
    duo = itertools.izip(*[iter(infile)]*2)
    depth = []
    thisdepth = 0
    while 1:
        try:
            itera = duo.next()[0]
        except StopIteration:
            break
        if itera != "//\n":
            try:
                thisdepth += int(itera.split(";")[-2][5:])
            except IndexError:
                ## TODO: if no clusts pass align filter this will raise
                LOGGER.debug("Here %s, %s", sample.name, itera)
                raise IPyradError("bad cluster file: %s", sample.name)
        else:
            ## append and reset
            depth.append(thisdepth)
            thisdepth = 0
    infile.close()

    if depth:
        ## make sense of stats
        depth = np.array(depth)
        keepmj = depth[depth >= data.paramsdict["mindepth_majrule"]]    
        keepstat = depth[depth >= data.paramsdict["mindepth_statistical"]]

        ## sample assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = len(depth)
        sample.stats["clusters_hidepth"] = max([len(i) for i in \
                                             (keepmj, keepstat)])
        sample.depths = depth

        data._stamp("s3 clustering on "+sample.name)        
    else:
        print("no clusters found for {}".format(sample.name))

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["assembly_method"] == "denovo":
        refmap_stats(data, sample)



def muscle_align(args):
    """ aligns reads, does split then aligning for paired reads """
    ## parse args
    data, chunk = args

    ## data are already chunked, read in the whole thing
    infile = open(chunk, 'rb')
    clusts = infile.read().split("//\n//\n")
    out = []

    ## iterate over clusters and align
    for clust in clusts:
        stack = []
        lines = clust.split("\n")
        names = lines[::2]
        seqs = lines[1::2]

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+str(i) for i, j in enumerate(names)]        

        ## don't bother aligning singletons
        if len(names) <= 1:
            if names:
                stack = [names[0]+"\n"+seqs[0]]
        else:
            ## split seqs if paired end seqs
            if all(["nnnn" in i for i in seqs]):
                seqs1 = [i.split("nnnn")[0] for i in seqs] 
                seqs2 = [i.split("nnnn")[1] for i in seqs]
                ## muscle align
                string1 = muscle_call(data, names[:200], seqs1[:200])
                string2 = muscle_call(data, names[:200], seqs2[:200])
                ## resort so they're in same order as original
                anames, aseqs1 = parsemuscle(string1)
                anames, aseqs2 = parsemuscle(string2)
                ## get leftlimit of seed, no hits can go left of this 
                ## this can save pairgbs from garbage
                idxs = [i for i, j in enumerate(aseqs1[0]) if j != "-"]
                leftlimit = min(0, idxs)
                aseqs1 = [i[leftlimit:] for i in aseqs1]
                ## preserve order in ordereddict
                aseqs = zip(aseqs1, aseqs2) 
                somedic = OrderedDict()

                ## append to somedic if not bad align
                for i in range(len(anames)):
                    ## filter for max internal indels
                    intindels1 = aseqs[i][0].rstrip('-').lstrip('-').count('-')
                    intindels2 = aseqs[i][1].rstrip('-').lstrip('-').count('-')

                    if intindels1 <= data.paramsdict["max_Indels_locus"][0] & \
                       intindels2 <= data.paramsdict["max_Indels_locus"][1]:
                        somedic[anames[i]] = aseqs[i][0]+"nnnn"+aseqs[i][1]
                    else:
                        LOGGER.info("high indels: %s", aseqs[i])

            else:
                string1 = muscle_call(data, names[:200], seqs[:200])
                anames, aseqs = parsemuscle(string1)
                ## Get left and right limits, no hits can go outside of this. 
                ## This can save gbs overlap data significantly. 
                if 'gbs' in data.paramsdict['datatype']:
                    ## left side filter is the left limit of the seed
                    idxs = [i for i, j in enumerate(aseqs[0]) if j != "-"]
                    if idxs:
                        leftlimit = max(0, min(idxs))
                        aseqs = [i[leftlimit:] for i in aseqs]
                        #LOGGER.info('leftlimit %s', leftlimit)

                    ## right side filter is the reverse seq that goes the least
                    ## far to the right.
                    revs = [aseqs[i] for i in range(len(aseqs)) if \
                            anames[i].split(";")[-1][0] == '-']
                    idxs = []
                    for rseq in revs:
                        idxs.append(max(\
                            [i for i, j in enumerate(rseq) if j != "-"]))
                    if idxs:
                        rightlimit = min(idxs)
                        aseqs = [i[:rightlimit] for i in aseqs]
                        #LOGGER.info('rightlimit %s', leftlimit)                    

                somedic = OrderedDict()
                for i in range(len(anames)):
                    ## filter for max internal indels and max indel inserts
                    intind = aseqs[i].rstrip('-').lstrip('-')
                    ind1 = intind.count('-') <= \
                                data.paramsdict["max_Indels_locus"][0]
                    #ind2 = len([i.split("-") for i in intind if i]) < 3
                    if ind1:
                        somedic[anames[i]] = aseqs[i]
                    else:
                        LOGGER.info("high indels: %s", aseqs[i])

            ## save dict into a list ready for printing
            for key in somedic.keys():
                stack.append(key+"\n"+somedic[key])

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    infile.close()
    outfile = open(chunk, 'wb')#+"_tmpout_"+str(num))
    outfile.write("\n//\n//\n".join(out)+"\n")#//\n//\n")
    outfile.close()



def parsemuscle(out):
    """ parse muscle string output into two sorted lists. Sorts them first."""
    lines = out[1:].split("\n>")
    names = [line.split("\n", 1)[0] for line in lines]
    seqs = [line.split("\n", 1)[1].replace("\n", "") for line in lines]
    tups = zip(names, seqs)
    ## who knew, zip(*) is the inverse of zip
    anames, aseqs = zip(*sorted(tups, 
                        key=lambda x: int(x[0].split(";")[-1][1:])))
    return anames, aseqs



def muscle_call(data, names, seqs):
    """ Makes subprocess call to muscle. A little faster than before. 
    TODO: Need to make sure this works on super large strings and does not 
    overload the PIPE buffer.
     """
    inputstring = "\n".join(">"+i+"\n"+j for i, j in zip(names, seqs))
    return subprocess.Popen(data.bins.muscle, 
                            stdin=subprocess.PIPE, 
                            stdout=subprocess.PIPE)\
                            .communicate(inputstring)[0]



def build_clusters(data, sample):
    """ combines information from .utemp and .htemp files 
    to create .clust files, which contain un-aligned clusters """

    ## derepfile 
    derepfile = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")

    ## vsearch results files
    ufile = os.path.join(data.dirs.clusts, sample.name+".utemp")
    htempfile = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## create an output file to write clusters to        
    clustfile = gzip.open(os.path.join(data.dirs.clusts,
                          sample.name+".clust.gz"), 'wb')
    sample.files["clusts"] = clustfile

    ## if .u files are present read them in as userout
    try:
        userout = open(ufile, 'rb').readlines()
    except IOError:
        inst = """
    No clusters (.utemp hits) found for {}. If you are running preview mode and
    the size of the truncated input file isn't big enough try increasing the 
    size of <your_assembly>._hackersonly[\"preview_truncate_length\"
    """.format(sample.name)
        LOGGER.warn(inst)
        raise IPyradError(inst)

    ## load derep reads into a dictionary
    hits = {}  
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
            udic[">"+uitems[1]+'\n'].append([">"+uitems[0]+"\n", 
                                            uitems[4],
                                            uitems[5].strip(),
                                            uitems[3]])
        else:
            ## write as seed
            udic[">"+uitems[1]+'\n'] = [[">"+uitems[0]+"\n", 
                                         uitems[4],
                                         uitems[5].strip(),
                                         uitems[3]]]

    ## map sequences to clust file in order
    seqslist = [] 
    for key, values in udic.iteritems():
        ## this is the seed. Store the left most non indel base for the seed. 
        ## Do not allow any other hits to go left of this (for pairddrad)
        seedhit = hits[key][1]
        seq = [key.strip()+"*\n"+seedhit]

        ## allow only N internal indels in hits to seed for within-sample clust
        ## prior to alignment. This improves alignments. 
        ## could be written a little more cleanly but this allows better debug.
        for i in xrange(len(values)):
            inserts = int(values[i][3])
            if values[i][1] == "+":
                if inserts < 6:
                    fwdseq = hits[values[i][0]][1]
                    seq.append(values[i][0].strip()+"+\n"+fwdseq)
                else:
                    LOGGER.debug("exc indbld: %s %s", inserts, fwdseq)
            ## flip to the right orientation 
            else:
                if inserts < 6:
                    revseq = comp(hits[values[i][0]][1][::-1])
                    seq.append(values[i][0].strip()+"-\n"+revseq)
                else:
                    LOGGER.debug("exc indbld: %s %s", inserts, revseq)

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
    #clustfile.write("//\n//\n\n")
    clustfile.close()
    del dereps
    del userout
    del udic
    del seedsdic



def get_threaded_view(ipyclient, split=True):
    """ gets optimum threaded view of ids given the host setup """
    ## engine ids
    ## e.g., [0, 1, 2, 3, 4, 5, 6, 7, 8]
    eids = ipyclient.ids  

    ## get host names    
    ## e.g., ['a', 'a', 'b', 'b', 'a', 'c', 'c', 'c', 'c']
    dview = ipyclient.direct_view()
    hosts = dview.apply_sync(socket.gethostname)
    del dview

    ## group ids into a dict by their hostnames
    ## e.g., {a: [0, 1, 4], b: [2, 3], c: [5, 6, 7, 8]}
    hostdict = defaultdict(list)
    for host, eid in zip(hosts, eids):
        hostdict[host].append(eid)

    ## Now split threads on the same host into separate proc if there are many
    hostdictkeys = hostdict.keys()
    for key in hostdictkeys:
        gids = hostdict[key]
        maxt = len(gids)
        if len(gids) >= 4:
            maxt = 2
        if len(gids) >= 6:
            maxt = 3
        if len(gids) >= 16:
            maxt = 4
        ## split ids into groups of maxt
        threaded = [gids[i:i+maxt] for i in xrange(0, len(gids), maxt)]
        lth = len(threaded)
        ## if anything was split (lth>1) update hostdict with new proc
        if lth > 1:
            hostdict.pop(key)
            for hostid in range(lth):
                hostdict[str(key)+"_"+str(hostid)] = threaded[hostid]

    ## logprint for viewing
    toview = "\n".join([str(key)+":"+str(val) for key, val in hostdict.items()])
    LOGGER.info("threaded_view %s", toview)

    ## make sure split numbering is correct
    threaded = hostdict.values()                
    #assert len(ipyclient.ids) <= len(list(itertools.chain(*threaded)))
    print("threaded_view:", dict(hostdict))
    return threaded



def split_among_processors(data, samples, ipyclient, noreverse, force, preview):
    """ pass the samples to N engines to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param ipyclient: ipyparallel load_balanced_view client
    :param noreverse: toggle revcomp clustering despite datatype default
    :param force: force
    :param preview: run preview
    :param align_only: skips clustering/mapping and aligns existing files

    :returns: None
    """
    ## make output folder for clusters  
    data.dirs.clusts = os.path.join(
                os.path.realpath(data.paramsdict["working_directory"]),
                data.name+"_"+"clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.makedirs(data.dirs.clusts)

    ## If ref mapping, init samples and make the refmapping output directory.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                        os.path.realpath(data.paramsdict["working_directory"]),
                        data.name+"_refmapping")
        if not os.path.exists(data.dirs.refmapping):
            os.makedirs(data.dirs.refmapping)

    ## Create threaded_view of engines by grouping only ids that are threaded
    ## load_balanced is not the right thing to use here, since we need to know
    ## which engines are being used for each job

    ## multi-node threading
    #threaded = get_threaded_view(ipyclient)
    #threaded_view = ipyclient.load_balanced_view(
    #                    targets=[i[0] for i in threaded])
    #tpproc = [len(i) for i in threaded][0]

    try:
        if preview:
            ## Truncate the input files and temporarily swap out the values for
            ## sample.files.edits

            for sample in samples:
                sample.files["edits_preview_bak"] = sample.files.edits
                sample.files.edits = preview_truncate_fq(data, sample.files.edits_preview_bak)

        ## single node threading
        tpp = 3
        threaded_view = ipyclient.load_balanced_view(targets=ipyclient.ids[::tpp])
        tpproc = len(threaded_view)

        ## Initialize the mapped and unmapped file paths per sample
        if not data.paramsdict["assembly_method"] == "denovo":
            for sample in samples:
                sample = refmap_init(data, sample)

        ## LIST ARGS for funcs mapreads and clustall
        submitted_args = []
        for sample in samples:
            #if not align_only:
            submitted_args.append([data, sample, noreverse, tpproc])
    
        ## MAP READS if reference sequence is specified
        if not data.paramsdict["assembly_method"] == "denovo":
            results = threaded_view.map(mapreads, submitted_args)
            results.get()
    
        ## DENOVO CLUSTER returns 0/1 of whether clust.gz built without fail
        ## for samples in the order in which they were submitted
        results = threaded_view.map(clustall, submitted_args)
        results = results.get()
    
        ## record that sample is clustered but not yet aligned
        for success, sample in zip(results, samples):
            if success:
                sample.stats.state = 2.5
    
        ## TODO: should it still look for REFSEQ reads if it had no utemp hits?
    
        ## REFSEQ reads (not denovo or denovo_only) pull in alignments from mapped 
        ## bam files and write them to the clust.gz files to fold them back into 
        ## the pipeline.
        if "denovo" not in data.paramsdict["assembly_method"]: 
            for sample in samples:
                finalize_aligned_reads(data, sample, ipyclient)
    
        ## call to ipp for muscle aligning only if the Sample passed clust/mapping
        for sample in samples:
            if sample.stats.state == 2.5:
                LOGGER.info("muscle aligning")            
                multi_muscle_align(data, sample, ipyclient)
                ## run sample cleanup 
                sample_cleanup(data, sample)
    
        ## run data cleanup
        data_cleanup(data, samples)

    except Exception as e:
        LOGGER.warn(e)
        raise
    finally:
        ## For preview/refmap restore original sample.files.edits paths and clean up
        ## the tmp files.

        ## If we did refmapping return the samples.files.edits to their original
        ## condition. Have to do this before restoring preview files because
        ## the refmap edits backup will be a backup of the preview truncated files.
        ## The order of these two functions matters.
        if "denovo" not in data.paramsdict["assembly_method"]: 
            for sample in samples:
                refmap_cleanup(data, sample)

        if preview:
            for sample in samples:
                try:
                    ## If edits and edits_preview_bak are the same then something borked
                    ## so don't delete any files
                    if sample.files.edits == sample.files.edits_preview_bak:
                        sample.files.pop("edits_preview_bak", None)
                        continue

                    for f in sample.files.edits[0]:
                        if(os.path.exists(f)):
                            os.remove( f )

                    ## Restore original paths to full fastq files
                    sample.files.edits = sample.files.edits_preview_bak
                    ## Remove the tmp file reference. The second arg defines what to return
                    ## if the key doesn't exist.
                    sample.files.pop("edits_preview_bak", None)
                except:
                    pass


def data_cleanup(data, samples):
    """ cleanup / statswriting function for Assembly obj """
    ## TODO: update statsfile with mapped and unmapped reads for refmapping
    data.statsfiles.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    if not os.path.exists(data.statsfiles.s3):
        with open(data.statsfiles.s3, 'w') as outfile:
            outfile.write(""+\
            "{:<20}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}   {:>9}\n""".\
                format("sample", "N_reads", "clusts_tot", 
                       "clusts_hidepth", "avg.depth.tot", 
                       "avg.depth>mj", "avg.depth>stat"))

    ## append stats to file
    outfile = open(data.statsfiles.s3, 'a+')
    samples.sort(key=lambda x: x.name)
    for sample in samples:
        ## only do for finished Samples
        if sample.stats.state == 3:
            outfile.write(""+\
    "{:<20}   {:>9}   {:>10}   {:>11}   {:>13.2f}   {:>11.2f}   {:>13.2f}\n".\
                      format(sample.name, 
                             int(sample.stats["reads_filtered"]),
                             int(sample.stats["clusters_total"]),
                             int(sample.stats["clusters_hidepth"]),
                             np.mean(sample.depths),
                             np.mean(sample.depths[sample.depths >= \
                                     data.paramsdict["mindepth_majrule"]]),
                             np.mean(sample.depths[sample.depths >= \
                                     data.paramsdict["mindepth_statistical"]])
                             ))
    outfile.close()



def concat_edits(data, sample):
    """ concatenate if multiple edits files for a sample """
    LOGGER.debug("Entering concat_edits: %s", sample.name)
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



def derep_and_sort(data, sample, nthreads):
    """ dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file """

    LOGGER.debug("Entering derep_and_sort: %s", sample.name)

    ## reverse complement clustering for some types    
    if "gbs" in data.paramsdict["datatype"]:
        reverse = " -strand both "
    else:
        reverse = " "

    LOGGER.debug("derep FILE %s", sample.files.edits[0][0])

    ## do dereplication with vsearch
    cmd = data.bins.vsearch+\
          " -derep_fulllength "+sample.files.edits[0][0]\
         +reverse \
         +" -output "+os.path.join(data.dirs.edits, sample.name+"_derep.fastq")\
         +" -sizeout " \
         +" -threads "+str(nthreads) \
         +" -fasta_width 0"

    ## run vsearch
    LOGGER.debug("%s", cmd)
    try:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(cmd)
        LOGGER.error(inst)
        raise IPyradError("Error in vsearch: \n{}\n{}\n{}."\
                          .format(inst, subprocess.STDOUT, cmd))



def cluster(data, sample, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    derephandle = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## datatype specific optimization
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .30 " 
    elif data.paramsdict["datatype"] == 'pairgbs':
        reverse = "  -strand both "
        cov = " -query_cov .90 " 
    else:  ## rad, ddrad
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override query cov
    if data._hackersonly["query_cov"]:
        cov = " -query_cov "+str(data._hackersonly["query_cov"])
        assert data._hackersonly["query_cov"] < 1, "query_cov must be < 1.0"

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        LOGGER.warn(noreverse, "not performing reverse complement clustering")

    ## get call string
    cmd = data.bins.vsearch+\
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
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



def multi_muscle_align(data, sample, ipyclient):
    """ Splits the muscle alignment across nthreads processors, each runs on 
    1000 clusters at a time. This is a kludge until I find how to write a 
    better wrapper for muscle. 
    """
    ## load balanced is fine for this b/c its on chunks and no threads
    lbview = ipyclient.load_balanced_view()

    ## split clust.gz file into nthreads*10 bits cluster bits
    tmpnames = []
    tmpdir = os.path.join(data.dirs.working, 'tmpalign')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    try: 
        ## get the number of clusters
        clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
        clustio = gzip.open(clustfile, 'rb')
        optim = 100
        if sample.stats.clusters_total > 2000:
            optim = int(sample.stats.clusters_total/len(ipyclient))/2
            print("optim: ", optim)
        ## write optim clusters to each tmp file
        inclusts = iter(clustio.read().strip().split("//\n//\n"))
        grabchunk = list(itertools.islice(inclusts, optim))
        while grabchunk:
            with tempfile.NamedTemporaryFile('w+b', delete=False, dir=tmpdir,
                                             prefix=sample.name+"_", 
                                             suffix='.ali') as out:
                out.write("//\n//\n".join(grabchunk))
            tmpnames.append(out.name)
            grabchunk = list(itertools.islice(inclusts, optim))
        clustio.close()

        ## create job queue
        submitted_args = []
        for fname in tmpnames:
            submitted_args.append([data, fname])

        ## run muscle on all tmp files. 'map' returns in same order          
        results = lbview.map_async(muscle_align, submitted_args)
        results.get()

        ## concatenate finished reads
        sample.files.clusters = os.path.join(data.dirs.clusts,
                                             sample.name+".clustS.gz")
        with gzip.open(sample.files.clusters, 'wb') as out:
            for fname in tmpnames:
                with open(fname) as infile:
                    out.write(infile.read()+"//\n//\n")

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        ## still delete tmpfiles if job was interrupted
        for fname in tmpnames:
            if os.path.exists(fname):
                os.remove(fname)
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)

        del lbview



def clustall(args):
    """ Running on remote Engine. Refmaps, then merges, then dereplicates, 
    then denovo clusters reads. """

    ## get args
    data, sample, noreverse, nthreads = args
    LOGGER.debug("clustall() %s", sample.name)

    ## concatenate edits files in case a Sample has multiple, and 
    ## returns a new Sample.files.edits with the concat file
    sample = concat_edits(data, sample)

    ## merge fastq pairs
    if 'pair' in data.paramsdict['datatype']:
        ## merge pairs that overlap and combine non-overlapping
        ## pairs into one merged file. merge_pairs takes the unmerged
        ## files list as an argument because we're reusing this code 
        ## in the refmap pipeline, trying to generalize.
        LOGGER.debug("Merging pairs - %s", sample.files)
        #unmerged_files = [sample.files.edits[0][0], sample.files.edits[0][1]]
        mergefile, nmerged = merge_pairs(data, sample) #, unmerged_files)
        sample.files.edits = [(mergefile, )]
        sample.files.pairs = mergefile
        sample.stats.reads_merged = nmerged
        sample.merged = 1
        LOGGER.debug("Merged file - {}".format(mergefile))

        ## OLD DEREN CODE w/ combine_pairs (keeping for now)
        ## merge pairs that overlap into a merge file
        #sample = merge_fastq_pairs(data, sample)
        ## pair together remaining pairs into a pairs file
        #sample = combine_pairs(data, sample)
        ## if merge file, append to pairs file
        #with open(sample.files.pairs, 'a') as tmpout:
        #    with open(sample.files.merged, 'r') as tmpin:
        #        tmpout.write(tmpin.read().replace("_c1\n", "_m4\n"))
        #LOGGER.info(sample.files.edits)

    ## convert fastq to fasta, then derep and sort reads by their size
    #LOGGER.debug(data, sample, nthreads)
    derep_and_sort(data, sample, nthreads)
    #LOGGER.debug(data, sample, nthreads)
    
    ## cluster derep fasta files in vsearch 
    cluster(data, sample, noreverse, nthreads)

    ## cluster_rebuild. Stop and print warning if no .utemp hits
    try:
        build_clusters(data, sample)
        ## record that it passed the clustfile build
        return 1
    except IPyradError as inst:
        print(inst)
        return 0



def run(data, samples, noreverse, force, preview, ipyclient):
    """ run the major functions for clustering within samples """

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for sample in samples:
        ## If sample not in state 2 don't try to cluster it.
        if sample.stats.state < 2:
            print("  Sample not ready for clustering. Please run "
                    + "step 2 on sample: {}".format(sample.name))
            continue

        if not force:
            if sample.stats.state >= 3:
                print("  Skipping {}; aleady clustered. Use force to recluster".\
                      format(sample.name))
            else:
                if sample.stats.reads_filtered:
                    subsamples.append(sample)
        else:
            ## force to overwrite
            if sample.stats.reads_filtered:            
                subsamples.append(sample)

    ## run subsamples 
    if not subsamples:
        print("  No Samples ready to be clustered. First run step2().")
    else:
        args = [data, subsamples, ipyclient, noreverse, force, preview]
        ## Because this includes long running jobs we want to make sure
        ## that if any Engine fails then all Engines stop and new jobs
        ## are not submitted. This also means we need good Exception 
        ## handling within 'split_among_processors' so that exceptions 
        ## we know can arise (e.g., no utemp file) are handled without fails.
        ## This will catch other exceptions such as KeyboardInterrupt, shut
        ## everything down, and then re-raise the Exception. 
        #try:
            ## get pids in case you need to kill them later            
            #pids = ipyclient[:].apply_async(os.getpid).get()
            ## submit jobs
        split_among_processors(*args)

        #except Exception as inst:
            ## prevents waiting jobs from being added to queue
        #    ipyclient.abort()
            ## kill Engines left executing (e.g., vsearch, smalt). 
            ## currently it kills your ipcluster instance tho... not ideal.
            ## maybe this should be in the Assembly wrapper...?
            ##ipyclient.shutdown()
            ##ipyclient = ipp.Client()
        #    raise inst



if __name__ == "__main__":
    ## test...
    import ipyrad as ip

    ## reload autosaved data. In case you quit and came back 
    data1 = ip.load.load_assembly("test_rad/data1.assembly")

    ## run step 6
    data1.step3(force=True)

    # DATA = Assembly("test")
    # DATA.get_params()
    # DATA.set_params(1, "./")
    # DATA.set_params(28, '/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa')
    # DATA.get_params()
    # print(DATA.log)
    # DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
    pass
    
