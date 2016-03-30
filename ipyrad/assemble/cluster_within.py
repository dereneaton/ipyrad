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
import time
import tempfile
import itertools
import subprocess
import numpy as np
import ipyrad

from collections import OrderedDict
from refmap import *
from util import * 

import logging
LOGGER = logging.getLogger(__name__)


def sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """
    
    ## set cluster file handles
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")

    if 'pair' in data.paramsdict["datatype"]:
        sample.files.merged = os.path.join(data.dirs.edits,
                                           sample.name+"_merged_.fastq")
        ## record how many read pairs were merged
        with open(sample.files.merged, 'r') as tmpf:
            sample.stats.reads_merged = len(tmpf.readlines()) // 4

    ## get depth stats
    infile = gzip.open(sample.files.clusters)
    duo = itertools.izip(*[iter(infile)]*2)
    depth = []
    thisdepth = 0
    maxlen = 0
    while 1:
        try:
            itera = duo.next()
            seqlen = len(itera[1])
            itera = itera[0]
            if seqlen > maxlen:
                maxlen = seqlen
        except StopIteration:
            break
        if itera != "//\n":
            try:
                tdepth = int(itera.split(";")[-2][5:])
                ## double depth for merged reads
                if "_m1;s" in itera:
                    tdepth *= 2 
                thisdepth += tdepth

            except IndexError:
                ## TODO: if no clusts pass align filter this will raise
                LOGGER.debug("Here %s, %s", sample.name, itera)
                raise IPyradError("ERROR 63: bad cluster file: %s", sample.name)
        else:
            ## append and reset
            depth.append(thisdepth)
            thisdepth = 0
    infile.close()

    ## If our longest sequence is longer than the current max_fragment_length
    ## then update max_fragment_length
    ## Padding with 4 extra bases to account for pair + sequence separator.
    if maxlen > data._hackersonly["max_fragment_length"]:
        data._hackersonly["max_fragment_length"] = maxlen + 4

    if depth:
        ## make sense of stats
        depth = np.array(depth)
        keepmj = depth[depth >= data.paramsdict["mindepth_majrule"]]    
        keepstat = depth[depth >= data.paramsdict["mindepth_statistical"]]

        ## sample summary stat assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = len(depth)
        sample.stats["clusters_hidepth"] = \
                                    max([len(i) for i in (keepmj, keepstat)])
        ## store large list of depths. Maybe turn this into a Counter
        ## dict to save space. We'll see...
        sample.depths = depth
        mmj = np.mean(sample.depths[sample.depths >= \
                                    data.paramsdict["mindepth_majrule"]])
        mms = np.mean(sample.depths[sample.depths >= \
                                    data.paramsdict["mindepth_statistical"]])
       
        ## sample stat assignments
        sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
        sample.stats_dfs.s3["clusters_total"] = len(depth)
        sample.stats_dfs.s3["clusters_hidepth"] = \
                                        int(sample.stats["clusters_hidepth"])
        sample.stats_dfs.s3["avg_depth_total"] = np.mean(sample.depths)
        sample.stats_dfs.s3["avg_depth_mj"] = mmj
        sample.stats_dfs.s3["avg_depth_stat"] = mms

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

    LOGGER.debug("Doing chunk %s", chunk)

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
            try:
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

                ## append to somedic if not bad align
                for i in range(len(anames)):
                    ## filter for max internal indels
                    intindels1 = aseqs[i][0].rstrip('-').lstrip('-').count('-')
                    intindels2 = aseqs[i][1].rstrip('-').lstrip('-').count('-')
                    #if intindels1 <= data.paramsdict["max_Indels_locus"][0] & \
                    #   intindels2 <= data.paramsdict["max_Indels_locus"][1]:
                    if (intindels1 <= 5) and (intindels2 <= 5):
                        stack.append("\n".join([anames[i], 
                                        aseqs[i][0]+"nnnn"+aseqs[i][1]]))
                        #somedic[anames[i]] = aseqs[i]]))
                    else:
                        LOGGER.info("""
                high indels: %s
                1, 2: (%s, %s)
                """, aseqs[i], intindels1, intindels2)

            except IndexError:
                string1 = muscle_call(data, names[:200], seqs[:200])
                anames, aseqs = parsemuscle(string1)
                ## Get left and right limits, no hits can go outside of this. 
                ## This can save gbs overlap data significantly. 
                if 'gbs' in data.paramsdict['datatype']:
                    ## left side filter is the left limit of the seed, unless 
                    ## there is a mixture of merged and non-merged reads. If 
                    ## any reads are merged then the left limit is the left
                    ## limit of a merged read.

                    # ## if merged
                    # if any(["_m1;s" in nam for nam in anames]):
                    #     LOGGER.info("""
                    #                 THIS IS A MERGE/NON_MERGE MATCH: 
                    #                 %s""", aseqs)
                    #     idxs = []
                    #     for nam, seq in zip(anames, aseqs):
                    #         if "_m1;s" in nam:
                    #             idxs.append(\
                    #           min([i for i, j in enumerate(seq) if j != "-"]))

                    ## much simpler if no merged reads
                    ## else:
                    idxs = [i for i, j in enumerate(aseqs[0]) if j != "-"]

                    ## apply left side filter
                    if idxs:
                        leftlimit = max(0, min(idxs))
                        aseqs = [i[leftlimit:] for i in aseqs]
                        LOGGER.info('leftlimit %s', leftlimit)

                    ## right side filter is the reverse seq that goes the least
                    ## far to the right. 
                    ## Get reverse seqs and their index (index, rseq)
                    revs = [(i, aseqs[i]) for i in range(len(aseqs)) if \
                            anames[i].split(";")[-1][0] == '-']
                    ## get right side filter and remove if there's a bad match.
                    idxs = []
                    for ridx, rseq in revs:
                        try:
                            idxs.append(max(\
                                [i for i, j in enumerate(rseq) if j != "-"]))
                        except ValueError as _:
                            LOGGER.debug("\
                            Found chunk that contains a locus that's all "\
                          +"indels. Throw it out and count it as filtered.")
                            ## Remove the seq name from the names list, and 
                            ## continue with the next iteration of the for loop,
                            ## effectively drops the rseq. Use list comprehension
                            ## to drop the idx'th element and then convert back to tuple
                            anames = tuple([x for i,x in enumerate(anames) if i!=ridx])
                            aseqs = tuple([x for i,x in enumerate(aseqs) if i!=ridx])
                            continue
                    if idxs:
                        rightlimit = min(idxs)
                        aseqs = [i[:rightlimit] for i in aseqs]
                        #LOGGER.info('rightlimit %s', rightlimit)                    

                for i in range(len(anames)):
                    ## filter for max internal indels 
                    intind = aseqs[i].rstrip('-').lstrip('-')
                    ind1 = intind.count('-') <= \
                                data.paramsdict["max_Indels_locus"][0]
                    ## could also filter by max indel inserts                                
                    #ind2 = len([i.split("-") for i in intind if i]) < 3
                    if ind1:
                        stack.append("\n".join([anames[i], aseqs[i]]))
                        #somedic[anames[i]] = aseqs[i]
                    else:
                        LOGGER.info("high indels: %s", aseqs[i])

            ## save dict into a list ready for printing
            #for key in somedic.keys():
            #    stack.append(key+"\n"+somedic[key])

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
    return subprocess.Popen(ipyrad.bins.muscle, 
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
                fwdseq = hits[values[i][0]][1]
                if inserts < 6:
                    seq.append(values[i][0].strip()+"+\n"+fwdseq)
                else:
                    fwdseq = hits[values[i][0]][1]                    
                    LOGGER.debug("exc indbld: %s %s", inserts, fwdseq)
            ## flip to the right orientation 
            else:
                revseq = comp(hits[values[i][0]][1][::-1])
                if inserts < 6:
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



def setup_dirs(data):
    """ sets up directories for step3 data """
    ## make output folder for clusters  
    data.dirs.clusts = os.path.join(
                os.path.realpath(data.paramsdict["project_dir"]),
                data.name+"_"+"clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.makedirs(data.dirs.clusts)

    ## If ref mapping, init samples and make the refmapping output directory.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                        os.path.realpath(data.paramsdict["project_dir"]),
                        data.name+"_refmapping")
        if not os.path.exists(data.dirs.refmapping):
            os.makedirs(data.dirs.refmapping)



def null_func(arglist):
    """ 
    Takes a list of args and prints them to the logger. Used as a null func
    to skip steps that are only used in some assembly steps. 
    """
    LOGGER.info("null func skipping %s", str(arglist))



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
    ## make directories
    setup_dirs(data)

    ## Create threaded_view of engines by grouping only ids that are threaded
    ## load_balanced is not the right thing to use here, since we need to know
    ## which engines are being used for each job
    hostdict = get_threaded_view(ipyclient)
    threaded_views = {}
    for host in hostdict:
        ## e.g., client.load_balanced_view([1,3])
        threaded_views[host] = ipyclient.load_balanced_view(host.values())
    ## A single load-balanced view for FUNCs 3-4
    lbview = ipyclient.load_balanced_view()

    ## FUNC 1: refmap
    ## A dictionary to store finished 
    res_ref = {sample:[] for sample in samples}

    ## default null func skips if no reference
    mcfunc = null_func
    if "reference" in data.paramsdict["assembly_method"]:
        ## set the func to mapreads
        mcfunc = mapreads
        ## Initialize the mapped and unmapped file paths per sample
        for sample in samples:
            sample = refmap_init(data, sample)

    ## Submit initial func for all samples (mapreads and/or clustall)    
    for sample, thv in zip(samples, itertools.cycle(threaded_views)):
        ## make args list with nthreads for this client
        args = [data, sample, noreverse, len(thv.values())]
        ## submit job to client w/ args
        res_ref[sample] = thv.apply_async(mcfunc, args)

    ## FUNC 2: clustall
    ## A dictionary to store whether sample has finished clustering
    res_clust = {sample:[] for sample in samples}

    ## default null func skips if no denovo
    if data.paramsdict["assembly_method"] in ["denovo", "denovo+reference"]: 
        mcfunc = clustall

    ## require that the sample successfully finished previous step
    for sample, thv in zip(samples, itertools.cycle(threaded_views)):
        with thv.temp_flags(after=res_ref[sample]):
            res_clust[sample] = thv.apply_async(mcfunc, args)

    ## FUNC 3: reference cleanup            
    ## A dictionary to store whether sample has finished cleanup
    res_clean = {sample:[] for sample in samples}

    ## Pull in alignments from mapped bam files and write them to the clust.gz 
    ## to fold them back into the pipeline. If we are doing "denovo" then 
    ## don't call this, but less obvious, "denovo-reference" intentionally 
    ## doesn't call this to effectively discard reference mapped reads.
    mcfunc = null_func
    if data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]: 
        mcfunc = finalize_aligned_reads

    ## requires sample to have finished the previous step before running
    for sample in samples:
        with lbview.temp_flags(after=res_clust[sample]):
            res_clean[sample] = lbview.apply_async(mcfunc, [data, sample])
            #TODO : finalize_aligned_reads(data, sample, ipyclient)

    ## FUNC 4: split up and submit chunks of sample file for aligning
    for sample in samples:
        with lbview.temp_flags(after=res_clean[sample]):
            lbview.apply(multi_muscle_align, [data, sample])


    ## FUNC 5: split up and submit chunks of sample file for aligning
    for sample in samples:
        with lbview.temp_flags(after=res_clean[sample]):
            lbview.apply(muscle_align, [data, fname])


    ## FUNC 6: split up and submit chunks of sample file for aligning
    for sample in samples:
        with lbview.temp_flags(after=res_clean[sample]):
            lbview.apply(cleanup, [data, sample])


    ## RUN ALL FUNCTIONS:
    lbview.wait_interactive()

    ## do sample cleanup
    for sample in samples:
        sample_cleanup(data, sample)
    LOGGER.debug("Finished sample cleanup. max_fragment_length = "\
                + "{}".format(data._hackersonly["max_fragment_length"]))
 
    ## run data cleanup
    data_cleanup(data, samples)



def split_among_cleanup(data, samples, preview):
    """
        ## For preview/refmap restore original sample.files.edits paths and 
        ## clean up the tmp files.

        ## If we did refmapping return the samples.files.edits to their original
        ## condition. Have to do this before restoring preview files because
        ## the refmap edits backup will be a backup of the preview truncated 
        ## files. The order of these two functions matters.
    """
    if not "denovo" == data.paramsdict["assembly_method"]: 
        for sample in samples:
            refmap_cleanup(data, sample)

    if preview:
        for sample in samples:
            try:
                ## If edits and edits_preview_bak are the same then 
                ## something borked so don't delete any files
                if sample.files.edits == sample.files.edits_preview_bak:
                    sample.files.pop("edits_preview_bak", None)
                    continue

                for tmpf in sample.files.edits[0]:
                    if os.path.exists(tmpf):
                        os.remove(tmpf)

                ## Restore original paths to full fastq files
                sample.files.edits = sample.files.edits_preview_bak
                ## Remove the tmp file reference. The second arg defines 
                ## what to return if the key doesn't exist.
                sample.files.pop("edits_preview_bak", None)
            except Exception:
                pass



def data_cleanup(data, samples):
    """ cleanup / statswriting function for Assembly obj """
    ## TODO: update stats_dfs with mapped and unmapped reads for refmapping
    data.stats_dfs.s3 = data.build_stat("s3")
    data.stats_files.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    with open(data.stats_files.s3, 'w') as outfile:
        data.stats_dfs.s3.to_string(outfile)



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
    cmd = ipyrad.bins.vsearch+\
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
    ## minsl: the percentage of the seed that must be matched
    ##    smaller values for RAD/ddRAD where we might want to combine, say 50bp 
    ##    reads and 100bp reads in the same analysis. 
    ## query_cov: the percentage of the query sequence that must match seed
    ##    smaller values are needed for gbs where only the tips might overlap
    ##    larger values for pairgbs where they should overlap near completely
    ##    small minsl and high query cov allows trimmed reads to match to untrim
    ##    seed for rad/ddrad/pairddrad. 
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .33 " 
        minsl = " 0.33"
    elif data.paramsdict["datatype"] == 'pairgbs':
        reverse = "  -strand both "
        cov = " -query_cov .75 " 
        minsl = " 0.75"  
    else:  ## rad, ddrad
        reverse = " -leftjust "
        cov = " -query_cov .90 "
        minsl = " 0.5"

    ## override query cov
    if data._hackersonly["query_cov"]:
        cov = " -query_cov "+str(data._hackersonly["query_cov"])
        assert data._hackersonly["query_cov"] <= 1, "query_cov must be <= 1.0"

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        LOGGER.warn(noreverse, "not performing reverse complement clustering")

    ## get call string
    cmd = ipyrad.bins.vsearch+\
        " -cluster_smallmem "+derephandle+\
        reverse+\
        cov+\
        " -id "+str(data.paramsdict["clust_threshold"])+\
        " -userout "+uhandle+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -minsl "+str(minsl)+\
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



def muscle_chunker(data, sample):
    """ 
    Splits the muscle alignment into chunks. Each runs on N clusters at a time.
    """
    ## split clust.gz file into nthreads*10 bits cluster bits
    tmpnames = []
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    ## get the number of clusters
    clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
    clustio = gzip.open(clustfile, 'rb')
    tclust = clustio.read().count("//")//2
    optim = int(tclust/10)
    LOGGER.debug("optim for align chunks: %s", optim)

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

    return tmpnames



def multi_muscle_align(data, sample, tmpname):
    """ aligns chunked seqs using muscle """



        ## concatenate finished reads
        sample.files.clusters = os.path.join(data.dirs.clusts,
                                             sample.name+".clustS.gz")
        with gzip.open(sample.files.clusters, 'wb') as out:
            for fname in tmpnames:
                with open(fname) as infile:
                    out.write(infile.read()+"//\n//\n")



def alignment_cleanup(tmpnames):
    ## still delete tmpfiles if job was interrupted
    for fname in tmpnames:
        if os.path.exists(fname):
            os.remove(fname)

    if os.path.exists(tmpdir):
        try:
            shutil.rmtree(tmpdir)
        except OSError as inst:
            ## In some instances nfs creates hidden dot files in directories
            ## that claim to be "busy" when you try to remove them. Don't
            ## kill the run if you can't remove this directory.
            LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))



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
        #mergefile, nmerged 
        sample = merge_pairs(data, sample) #, unmerged_files)
        sample.files.edits = [(sample.files.merged, )]
        sample.files.pairs = sample.files.merged
        #sample.stats.reads_merged = nmerged
        sample.merged = 1
        LOGGER.debug("Merged file - {}".format(sample.files.merged))

    ## convert fastq to fasta, then derep and sort reads by their size
    derep_and_sort(data, sample, nthreads)
    
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
            print("""\
    Sample not ready for clustering. First run step2 on sample: {}""".\
    format(sample.name))
            continue

        if not force:
            if sample.stats.state >= 3:
                print("""\
    Skipping {}; aleady clustered. Use force to re-cluster""".\
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
        LOGGER.debug("splitting among processors")
        try:
            split_among_processors(*args)
        except Exception as inst:
            LOGGER.warn(inst)
            raise
        finally:
            split_among_cleanup(data, samples, preview)





if __name__ == "__main__":
    ## test...

    ## reload autosaved data. In case you quit and came back 
    data1 = ip.load_json("test_rad/data1.assembly")

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

    
