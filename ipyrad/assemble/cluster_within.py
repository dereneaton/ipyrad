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
# pylint: disable=C0301

import os
import sys
import gzip
import itertools
import subprocess
import numpy as np
import ipyrad
import time
import datetime

from collections import Counter
from ipyparallel import Dependency, RemoteError
from refmap import *
from util import * 

import logging
LOGGER = logging.getLogger(__name__)



def sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """
    
    ## set cluster file handles
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")

    ## number of merged reads is updated dynamically
    ## TODO: this won't capture merged reads that are merged during refmap
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
    maxlen = 0
    thisdepth = []
    while 1:
        try:
            itera = duo.next()
        except StopIteration:
            break
        ## keep going until the end of cluster
        if itera[0] != "//\n":
            try:
                ## get longest seqlen 
                seqlen = len(itera[1])            
                ## get nreps in read
                tdepth = int(itera[0].split(";")[-2][5:])
                ## double depth for merged reads
                if "_m1;s" in itera[0]:
                    tdepth *= 2 
                thisdepth.append(tdepth)
            except IndexError:
                ## if no clusts pass align filter this will raise
                LOGGER.info("No aligned clusters passed for %s", sample.name)
                raise IPyradError("No aligned clusters passed: {}".format(sample.name))
        else:
            ## update maxlen
            if seqlen > maxlen:
                maxlen = seqlen
            ## append and reset
            depth.append(sum(thisdepth))
            thisdepth = []
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
                                max([i.shape[0] for i in (keepmj, keepstat)])
        ## store large list of depths as a counter dict
        sample.depths = dict(Counter(depth))

        ## get stats for depths
        mmj = depth[depth >= data.paramsdict["mindepth_majrule"]]
        mms = depth[depth >= data.paramsdict["mindepth_statistical"]]
       
        ## sample stat assignments
        sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
        sample.stats_dfs.s3["clusters_total"] = depth.shape[0]
        sample.stats_dfs.s3["clusters_hidepth"] = int(sample.stats["clusters_hidepth"])
        sample.stats_dfs.s3["avg_depth_total"] = depth.mean()
        sample.stats_dfs.s3["avg_depth_mj"] = mmj.mean()
        sample.stats_dfs.s3["avg_depth_stat"] = mms.mean()

        sample.stats_dfs.s3["sd_depth_total"] = depth.std()
        sample.stats_dfs.s3["sd_depth_mj"] = mmj.std()
        sample.stats_dfs.s3["sd_depth_stat"] = mms.std()

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

    LOGGER.debug("aligning chunk %s", chunk)

    ## data are already chunked, read in the whole thing. 
    ## bail out if there is not a chunk file, means there were few clusts
    try:
        infile = open(chunk, 'rb')
    except IOError:
        return 0
    clusts = infile.read().split("//\n//\n")
    out = []

    ## a counter for discarded clusters due to poor alignment
    highindels = 0

    ## iterate over clusters and align
    for clust in clusts:
        stack = []
        lines = clust.split("\n")
        names = lines[::2]
        seqs = lines[1::2]
        badalign = 0

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+str(i) for i, j in enumerate(names)]

        ## 0 length names indicates back to back //\n//\n sequences, so pass
        ## TODO: This should probably be fixed upstream so it doesn't happen
        ## but it's protective regardless... but still, it's slowing us down...
        if len(names) == 0:
            pass
        ## don't bother aligning singletons
        elif len(names) <= 1:
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
                    else:
                        highindels += 1
                        badalign = 1
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
                            anames = tuple([x for i, x in enumerate(anames) if i != ridx])
                            aseqs = tuple([x for i, x in enumerate(aseqs) if i != ridx])
                            continue
                    if idxs:
                        rightlimit = min(idxs)
                        aseqs = [i[:rightlimit] for i in aseqs]
                        #LOGGER.info('rightlimit %s', rightlimit)                    

                badalign = 0
                for i in range(len(anames)):
                    ## filter for max internal indels 
                    intind = aseqs[i].rstrip('-').lstrip('-')
                    ind1 = intind.count('-') <= \
                                data.paramsdict["max_Indels_locus"][0]
                    if ind1:
                        stack.append("\n".join([anames[i], aseqs[i]]))
                    else:
                        highindels += 1
                        badalign = 1
                        LOGGER.info("high indels: %s", aseqs[i])

        ## finally, add to outstack if alignment is good
        if stack:
            if not badalign:
                out.append("\n".join(stack))

    ## write to file after
    infile.close()
    outfile = open(chunk, 'wb')
    outfile.write("\n//\n//\n".join(out)+"\n")
    outfile.close()
    return highindels



def parsemuscle(out):
    """ parse muscle string output into two sorted lists. Sorts them first."""
    lines = out[1:].split("\n>")
    names = [line.split("\n", 1)[0] for line in lines]
    seqs = [line.split("\n", 1)[1].replace("\n", "") for line in lines]
    tups = zip(names, seqs)
    ## who knew, zip(*) is the inverse of zip
    anames, aseqs = zip(*sorted(tups, 
                        key=lambda x: int(x[0].split(";")[-1][1:])))
    ## remove reference if present
    if "_REF;+0" in anames[0]:
        return anames[1:], aseqs[1:]
    else:
        return anames, aseqs





def muscle_call(data, names, seqs):
    """ 
    Makes subprocess call to muscle. A little faster than before. 
    TODO: Need to make sure this works on super large strings and does not 
    overload the PIPE buffer.
    """

    ## get input string
    inputstr = "\n".join([">{}\n{}".format(i, j) for i, j in zip(names, seqs)])
    args = [ipyrad.bins.muscle, "-quiet"]

    ## increase gap penalty if reference region is included
    if "_REF;+0" in names:
        args += ["-gapopen", "-1200"]

    ## make a call arg
    proc1 = subprocess.Popen(args,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE)
    ## return result
    return proc1.communicate(inputstr)[0]




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
        ## this is the seed. 
        seedhit = hits[key][1]
        seq = [key.strip()+"*\n"+seedhit]

        ## allow only N internal indels in hits to seed for within-sample clust
        ## prior to alignment. Exclude read that match poorly. This improves 
        ## alignments. Whole stack will be excluded after alignment if poor. 
        for i in xrange(len(values)):
            inserts = int(values[i][3])
            if values[i][1] == "+":
                fwdseq = hits[values[i][0]][1]
                if inserts < 6:
                    seq.append(values[i][0].strip()+"+\n"+fwdseq)
                else:
                    LOGGER.info("exc indbld: %s %s", inserts, fwdseq)
            ## flip to the right orientation 
            else:
                revseq = comp(hits[values[i][0]][1][::-1])
                if inserts < 6:
                    seq.append(values[i][0].strip()+"-\n"+revseq)
                else:
                    LOGGER.info("exc indbld: %s %s", inserts, revseq)

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

    ## clust directory
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    return tmpdir


def null_func(args):
    """ 
    Takes a list of args and prints them to the logger. Used as a null func
    to skip steps that are only used in some assembly steps. 
    """
    LOGGER.info("null func skipping")#  %s", str(arglist))



def apply_jobs(data, samples, ipyclient, noreverse, force, preview):
    """ pass the samples to N engines to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param ipyclient: ipyparallel load_balanced_view client
    :param noreverse: toggle revcomp clustering despite datatype default
    :param force: force
    :param preview: run preview
    :param align_only: skips clustering/mapping and aligns existing files

    ## This is all the internal documentation from the prior version of
    ## this function which I copied here but haven't updated much.
    Step 1: Derep and concat
    Step 2: mapreads ----------------------------------------------------
     Submit samples to reference map else null
    Step 3: clustering ---------------------------------------------------
     Cluster reads for all assembly methods except 'reference' since
     the refmapping clusters them for us and we don't care about the reads
     that don't map.
    Step 4: reference cleanup -------------------------------------------
     Pull in alignments from mapped bam files and write them to the clust.gz 
     to fold them back into the pipeline. If we are doing "denovo" then 
     don't call this, but less obvious, "denovo-reference" intentionally 
     doesn't call this to effectively discard reference mapped reads.
    Step 5: split up clusters into chunks -------------------------------
    Step 6: align chunks -------------------------------------------------
    Step 7: concat chunks -------------------------------------------------

    :returns: None
    """
    ## make directories. tmpdir is returned because muscle_align
    ## and muscle_chunk both rely on it
    tmpdir = setup_dirs(data)

    ## Create an ordered dict to aggregate info for each step. This dict must be ordered
    ## because the order determins both the order of the jobs that get spawned as well
    ## as the order of printing of progress bars. The contents of this dict will 
    ## eventually include:
    ##  printstr - the pretty printing string to make the progress bar look nice.
    ##  function - the function to all inside apply()
    ##  extra_args - Any extra args for the above function beyond [data, sample]
    ##               which get passed to all functions
    ##  read_flag - The asyncResult object for the monitor task that is dependent
    ##              on all the other tasks for this step completing. Use this
    ##              to determine when all the tasks for this step are done. derp.
    ##  async_results - A dict of samples and their corresponding asyncresults
    ##              objects for the call to the function at this step.
    ## "muscle_align" as has an additional record called `all_aligns` that is just
    ##              a giant list of all the asyncresults for all the samples.
    ##
    ## I made the spacing goofy to make it easier to see that all the printstr
    ## values are all the same length.
    ##
    ## Also, in order to control which steps run you simply don't add the steps
    ## you don't want to run to this dict. Since all dependencies are created
    ## relative to the members of the extant dict, the steps in the rocess are
    ## agnostic about which step preceded it, so adding/removing steps is super easy.
    ## All the logic for controlling which steps run for the various assembly
    ## methods is right here!
    steps = collections.OrderedDict()
    steps["derep_concat_split"] =   {"printstr": "dereplicating    ",\
                                    "function":derep_concat_split, "extra_args":[]}
    if "reference" in data.paramsdict["assembly_method"]:
        steps["mapreads"] =         {"printstr": "mapping          ",\
                                    "function":mapreads, "extra_args":[noreverse, 1]}
    if data.paramsdict["assembly_method"] != "reference":
        steps["clust_and_build"] =  {"printstr": "clustering       ",\
                                    "function":clust_and_build, "extra_args":[noreverse, 1]}
    if data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]: 
        steps["ref_muscle_chunker"] = {"printstr": "finalize mapping ",\
                                    "function":ref_muscle_chunker, "extra_args":[]}
    steps["muscle_chunker"] =       {"printstr": "chunking         ",\
                                    "function":muscle_chunker, "extra_args":[tmpdir]}
    steps["muscle_align"] =         {"printstr": "aligning         ",\
                                    "function":muscle_align, "extra_args":[]}
    steps["reconcat"] =             {"printstr": "concatenating    ",\
                                    "function":reconcat, "extra_args":[]}

    ## Create threaded_view of engines by grouping only ids that are threaded
    hostdict = get_threaded_view(ipyclient)
    threaded_views = {}
    for key, val in hostdict.items():
        ## e.g., client.load_balanced_view([1,3])
        threaded_views[key] = ipyclient.load_balanced_view(val)

    ## A single load-balanced view for FUNCs 3-4
    lbview = ipyclient.load_balanced_view()
    first_dependency = lbview.apply(os.getpid)
    while not first_dependency.ready():
        time.sleep(1)

    ## If doing reference sequence mapping in any fashion then init
    ## the samples. This is a very fast operation.
    if "reference" in data.paramsdict["assembly_method"]:
        samples = [refmap_init([data, s]) for s in samples]

    ## Now iterate through the ordered dictionary of steps and create the tasks.
    ## Roughly:
    ##  - Create the dependencies. In all cases except the first step, the task 
    ##      for the current sample is dependent on the state of the task for that
    ##      sample from the prefious step. For the first step you just use a dummy
    ##      dependency.
    ##  - Apply the function for the current step for each sample, pulling in any
    ##      extra args from the steps dictionary. Store the async_result for each
    ##      sample in the steps dict. In the case of the "muscle_align" step we
    ##      also store one giant list of all the async results because it's
    ##      convenient not to have to keep recreating it downstream. 
    ##
    ## There are only a very few # of difference in the process for each step.
    ## "muscle_align" is significantly different so it gets its own block, as well
    ## "reconcat" is a little different because the dependencies are weird.
    for i, step in enumerate(steps):
        stepdict = steps[step]

        stepdict["async_results"] = {}
        for sample in samples:
            ######################
            ## Set up dependencies
            ######################
            if i == 0:
                check_deps = first_dependency
            else:
                laststep = steps.items()[i-1][1]
                if step == "reconcat":
                    tmpids = list(itertools.chain(*[j.msg_ids for j in laststep["async_results"][sample]]))
                else:
                    tmpids = laststep["async_results"][sample]
                check_deps = Dependency(tmpids, failure=False, success=True)
            with lbview.temp_flags(after=check_deps):
                #######################################################
                ## Call apply with our functions and args for this step
                #######################################################
                if step == "muscle_align":
                    stepdict["async_results"][sample] = []
                    for j in range(10):
                        chunk = os.path.join(tmpdir, sample.name+"_chunk_{}.ali".format(j))
                        stepdict["async_results"][sample].append(lbview.apply(muscle_align, [data, chunk]))
                    all_aligns = list(itertools.chain(*stepdict["async_results"].values()))
                    stepdict["all_aligns"] = all_aligns
                else:
                        args = [data, sample]
                        args.extend(stepdict["extra_args"])
                        stepdict["async_results"][sample] = lbview.apply(stepdict["function"], args)

        ## Each step has one job that gets created after all other jobs are created
        ## This job is dependent on all the other child jobs for this step and we'll use it
        ## as the `ready_flag` to indicate when this step has entirely completed.
        if step == "muscle_align":
            tmpids = stepdict["all_aligns"]
        else:
            tmpids = [i for i in stepdict["async_results"].values()]
        with lbview.temp_flags(after=tmpids):
            stepdict["ready_flag"] = lbview.apply(time.sleep, 1)    

    ## Create one final task that sits at the back of the queue. Basically this is
    ## now only being used to keep track of elapsed time, we do not rely
    ## on the state of this job for anything else.
    tmpids = list(itertools.chain(*[i.msg_ids for i in steps["reconcat"]["async_results"].values()]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1) 

        ## If you want to test out any individual step synchronously then uncomment
        ## this block, it will wait for all samples to finish that step before
        ## proceding. This can be helpful for debugging.
#       try:
#           [step["async_results"][i].get() for i in step["async_results"]]
#       except Exception as inst:
#           print(inst)
#       for i in step["async_results"]:
#           print(step["async_results"][i].metadata.status)
#       return 1

    ## All ipclient jobs have now been created, so now we watch the results pile up
    ## Here we go through each step in order, check the ready_flag. If our current
    ## state is not ready then we update the progress bar. If the current state
    ## _is_ ready, then we finalize the progress bar and check the results.
    elapsed = "0:00:00"
    for step,vals in steps.iteritems():
        async_results = vals["async_results"]
        while not vals["ready_flag"].ready():
            ## prints a progress bar
            if step == "muscle_align":
                ## This is a little hackish. You have to divide by 10
                ## because there are ten muscle calls per sample.
                finished = sum([i.ready() for i in vals["all_aligns"]])/10
            else:
                finished = sum([async_results[i].ready() for i in async_results])
            elapsed = datetime.timedelta(seconds=int(res.elapsed))
            #print(finished)
            progressbar(len(vals["async_results"]), finished,
                " {}     | {}".format(vals["printstr"], elapsed))
            sys.stdout.flush()
            time.sleep(1)

        ## Print the finished progress bar for this step and go to the next line
        progressbar(100, 100,
            " {}     | {}".format(vals["printstr"], elapsed))
        print("")

        ## Check the results to see if any/all samples failed.
        try:
            if step == "muscle_align":
                failed_samples = check_results_alignment(async_results)
            else:
                failed_samples = check_results(async_results)
            if failed_samples:
                print(failed_samples)
        except IPyradError as inst:
            print("All samples failed {} - {}".format(step, inst))
            sys.exit()

    ## Cleanup -------------------------------------------------------
    for sample in samples:
        sample_cleanup(data, sample)
        ## Helpful debug for seeing which samples had bad alignments
        #for i in steps["muscle_align"][sample]:
        #    print(sample, i.get())
        badaligns = sum([i.get() for i in steps["muscle_align"]["async_results"][sample]])
        sample.stats_dfs.s3.filtered_bad_align = badaligns

    data_cleanup(data)


def check_results_alignment(async_results):
    """
    Since the alignment step forks a bunch more processes the results
    are a little different. The input for this function is a dict containing:
    {sample_name:[ipyparallel.asyncResult()]}. Just iterate through the list
    and call check_results() to do the real work for us.
    """
    failed_chunks = {}
    failed_samples = []
    errors = {}
    for sample, result in async_results.iteritems():
        failed_chunks[sample] = []
        errors[sample] = []
        for r in result:
            try:
                failed_chunks[sample].append(check_results({sample:r}))
            except Exception as inst:
                errors[sample].append(inst)

        ## If the length of the errors list for this sample == the length
        ## of the expected results then all chunks failed
        if len(errors[sample]) == len(result):
            failed_samples.append(sample)

    if len(failed_samples) == len(async_results):
        raise IPyradError("All samples failed alignment\n{}".format(errors.values()))

    return failed_samples

def check_results(async_results):
    """
    Check the output from each step. If all samples fail then raise an
    IPyradError. If only some samples fail try to print helpful messages.
    The argument for this function is a dictionary containing:
    {sample_name:ipyparallel.asyncResult()}
    """
#    for i in async_results:
#        print(i, async_results[i].ready(), async_results[i].error)
    errors = {}
    for sample, result in async_results.iteritems():
        try:
            result.get()
        except Exception as inst:
            errors[sample] = str(inst)
            LOGGER.warn(inst)

    ## if all samples fail then raise an error
    if len(errors) == len(async_results):
        if all([x == "None" for x in errors.values()]):
            msg = "All samples failed this stage because all samples"\
                    + " failed some previous stage."
        else:
            ## Fetch the actual error message
            msg = "\n".join([x for x in errors.values()])
        raise IPyradError("All samples failed:\n"+msg)

    if errors:
        print("  Samples failed this step:"+" ".join(errors.keys()))
    return errors.keys()


def jobs_cleanup(data, samples, preview):
    """
    For preview/refmap restore original sample.files.edits paths and 
    clean up the tmp files.

    If we did refmapping return the samples.files.edits to their original
    condition. Have to do this before restoring preview files because
    the refmap edits backup will be a backup of the preview truncated 
    files. The order of these two functions matters.
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



def derep_and_sort(data, sample):
    """ 
    Dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file. Paired
    reads are dereplicated as one concatenated read and later split again. 
    """

    LOGGER.info("in the real derep; %s", sample.name)
    ## reverse complement clustering for some types    
    if "gbs" in data.paramsdict["datatype"]:
        reverse = " -strand both "
    else:
        reverse = " "

    ## do dereplication with vsearch
    cmd = ipyrad.bins.vsearch\
         +" -derep_fulllength "+sample.files.edits[0][0]\
         +reverse \
         +" -output "+os.path.join(data.dirs.edits, sample.name+"_derep.fastq")\
         +" -sizeout " \
         +" -threads 1 "\
         +" -fasta_width 0"
    LOGGER.info(cmd)

    ## run vsearch
    try:
        LOGGER.info(cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(inst)
        raise IPyradError(inst)



def data_cleanup(data):
    """ cleanup / statswriting function for Assembly obj """
    data.stats_dfs.s3 = data.build_stat("s3")
    data.stats_files.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    with open(data.stats_files.s3, 'w') as outfile:
        data.stats_dfs.s3.to_string(outfile)



def concat_edits(data, sample):
    """ 
    Concatenate if multiple edits files for a sample
    """
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



def cluster(data, sample, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    if "reference" in data.paramsdict["assembly_method"]:
        derephandle = os.path.join(data.dirs.edits, sample.name+"-refmap_derep.fastq")
    else:
        derephandle = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## If derep file doesn't exist then bail out
    if not os.path.isfile(derephandle):
        LOGGER.warn("Bad derephandle - {}".format(derephandle))
        raise IPyradError("Input file for clustering doesn't exist - {}"\
                        .format(derephandle))    

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



def muscle_chunker(args):
    """ 
    Splits the muscle alignment into chunks. Each runs on N clusters at a time.
    """
    ## parse args
    LOGGER.info("inside chunker")
    data, sample, tmpdir = args

    ## get the number of clusters
    clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
    with gzip.open(clustfile, 'rb') as clustio:
        tclust = clustio.read().count("//")//2
        optim = (tclust//10) + (tclust%10)
        LOGGER.info("optim for align chunks: %s", optim)

    ## write optim clusters to each tmp file
    clustio = gzip.open(clustfile, 'rb')
    inclusts = iter(clustio.read().strip().split("//\n//\n"))
    grabchunk = list(itertools.islice(inclusts, optim))

    idx = 0
    while grabchunk:
        tmpfile = os.path.join(tmpdir, sample.name+"_chunk_{}.ali".format(idx)) 
        with open(tmpfile, 'wb') as out:
            out.write("//\n//\n".join(grabchunk))
        idx += 1
        grabchunk = list(itertools.islice(inclusts, optim))
    clustio.close()



def reconcat(args):
    """ aligns chunked seqs using muscle """

    ## parse args
    data, sample = args

    ## get chunks
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    chunks = glob.glob(os.path.join(tmpdir, sample.name+"_chunk_*"))
   
    ## sort by chunk number
    chunks.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-4]))

    LOGGER.info("chunk %s", chunks)
    ## concatenate finished reads
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")
    ## reconcats aligned clusters
    with gzip.open(sample.files.clusters, 'wb') as out:
        for fname in chunks:
            with open(fname) as infile:
                dat = infile.read()
                ## avoids mess if last chunk was empty
                if dat.endswith("\n"):
                    out.write(dat+"//\n//\n")
                else:
                    out.write(dat+"\n//\n//\n")
            os.remove(fname)



def alignment_cleanup(data):
    """ delete tmp-align dir and ..."""
    ## still delete tmpfiles if job was interrupted
    # for fname in tmpnames:
    #     if os.path.exists(fname):
    #         os.remove(fname)

    ## remove aligns dir
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    if os.path.exists(tmpdir):
        try:
            shutil.rmtree(tmpdir)
        except OSError as inst:
            ## In some instances nfs creates hidden dot files in directories
            ## that claim to be "busy" when you try to remove them. Don't
            ## kill the run if you can't remove this directory.
            LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))



def derep_concat_split(args):
    """ 
    Running on remote Engine. Refmaps, then merges, then dereplicates, 
    then denovo clusters reads. 
    """

    ## get args
    data, sample = args
    LOGGER.info("INSIDE derep %s", sample.name)

    ## concatenate edits files in case a Sample has multiple, and 
    ## returns a new Sample.files.edits with the concat file
    sample = concat_edits(data, sample)

    ## Denovo: merge or concat fastq pairs [sample.files.pairs]
    ## Reference: only concat fastq pairs  []
    ## Denovo + Reference: 
    if 'pair' in data.paramsdict['datatype']:
        ## merge pairs that overlap and combine non-overlapping
        ## pairs into one merged file. merge_pairs takes the unmerged
        ## files list as an argument because we're reusing this code 
        ## in the refmap pipeline, trying to generalize.
        LOGGER.debug("Merging pairs - %s", sample.files.edits)
        merge = 1
        revcomp = 1
        ## If doing any kind of reference mapping do not merge
        ## only concatenate so the reads can be split later and mapped
        ## separately. 
        if "reference" in data.paramsdict["assembly_method"]:
            merge = 0
            revcomp = 0
        sample.files.merged = os.path.join(data.dirs.edits,
                                        sample.name+"_merged_.fastq")
        sample.stats.reads_merged = merge_pairs(data, sample.files.edits, 
                                        sample.files.merged, revcomp, merge)
        LOGGER.info("Merged pairs - {} - {}".format(sample.name, \
                                        sample.stats.reads_merged))
        sample.files.edits = [(sample.files.merged, )]
        LOGGER.debug("Merged file - {}".format(sample.files.merged))

    ## convert fastq to fasta, then derep and sort reads by their size
    derep_and_sort(data, sample)
    


def clust_and_build(args):
    """ 
    cluster and build clusters
    """
    ## parse args
    data, sample, noreverse, nthreads = args
    LOGGER.debug("Entering clust_and_build - {}".format(sample))

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



def cleanup_and_die(async_results):
    LOGGER.debug("Entering cleanup_and_die")
    
    LOGGER.debug(async_results)
    res = []
    ## Sort through msg_ids and check metadata stats to figure out what went wrong
    for i in async_results:
        if not i == None:
            res.extend(i)
            LOGGER.warn("Got error - {}".format(i))
    return res



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
        ## arguments to apply_jobs, inst catches exceptions
        args = [data, subsamples, ipyclient, noreverse, force, preview]
        inst = ""

        ## wrap job in try/finally to ensure cleanup
        try:
            apply_jobs(*args)
        except Exception as inst:
            LOGGER.warn(inst)
            raise
        finally:
            alignment_cleanup(data)



if __name__ == "__main__":
    ## test...

    ## reload autosaved data. In case you quit and came back 
    DATA = ipyrad.load_json("cli/cli.json")

    ## run step 6
    DATA.step3(force=True)

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

    
