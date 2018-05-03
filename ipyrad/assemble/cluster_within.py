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
# pylint: disable=R0915
# pylint: disable=R0914
# pylint: disable=R0912

import os
import io
import gzip
import glob
import itertools

import numpy as np
import ipyrad
import time
import datetime
import warnings
import networkx as nx
import ipyparallel as ipp

from refmap import *
from util import *

## Python3 subprocess is faster for muscle-align
try:
    import subprocess32 as sps
except ImportError:
    import subprocess as sps

import logging
LOGGER = logging.getLogger(__name__)



def get_quick_depths(data, sample):
    """ iterate over clustS files to get data """

    ## use existing sample cluster path if it exists, since this
    ## func can be used in step 4 and that can occur after merging
    ## assemblies after step3, and if we then referenced by data.dirs.clusts
    ## the path would be broken.
    if sample.files.clusters:
        pass
    else:
        ## set cluster file handles
        sample.files.clusters = os.path.join(
            data.dirs.clusts, sample.name+".clustS.gz")

    ## get new clustered loci
    fclust = data.samples[sample.name].files.clusters
    clusters = gzip.open(fclust, 'r')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## storage
    depths = []
    maxlen = []

    ## start with cluster 0
    tdepth = 0
    tlen = 0

    ## iterate until empty
    while 1:
        ## grab next
        try:
            name, seq = pairdealer.next()
        except StopIteration:
            break

        ## if not the end of a cluster
        #print name.strip(), seq.strip()
        if name.strip() == seq.strip():
            depths.append(tdepth)
            maxlen.append(tlen)
            tlen = 0
            tdepth = 0

        else:
            tdepth += int(name.split(";")[-2][5:])
            tlen = len(seq)

    ## return
    clusters.close()
    return np.array(maxlen), np.array(depths)



def sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """

    ## get maxlen and depths array from clusters
    maxlens, depths = get_quick_depths(data, sample)

    try:
        depths.max()
    except ValueError:
        ## If depths is an empty array max() will raise
        print("    no clusters found for {}".format(sample.name))
        return

    ## Test if depths is non-empty, but just full of zeros.
    if depths.max():
        ## store which min was used to calculate hidepth here
        sample.stats_dfs.s3["hidepth_min"] = data.paramsdict["mindepth_majrule"]

        ## If our longest sequence is longer than the current max_fragment_length
        ## then update max_fragment_length. For assurance we require that
        ## max len is 4 greater than maxlen, to allow for pair separators.
        hidepths = depths >= data.paramsdict["mindepth_majrule"]
        maxlens = maxlens[hidepths]

        ## Handle the case where there are no hidepth clusters
        if maxlens.any():
            maxlen = int(maxlens.mean() + (2.*maxlens.std()))
        else:
            maxlen = 0
        if maxlen > data._hackersonly["max_fragment_length"]:
            data._hackersonly["max_fragment_length"] = maxlen + 4

        ## make sense of stats
        keepmj = depths[depths >= data.paramsdict["mindepth_majrule"]]
        keepstat = depths[depths >= data.paramsdict["mindepth_statistical"]]

        ## sample summary stat assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = depths.shape[0]
        sample.stats["clusters_hidepth"] = keepmj.shape[0]

        ## store depths histogram as a dict. Limit to first 25 bins
        bars, bins = np.histogram(depths, bins=range(1, 26))
        sample.depths = {int(i):v for i, v in zip(bins, bars) if v}

        ## sample stat assignments
        ## Trap numpy warnings ("mean of empty slice") printed by samples
        ## with few reads.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
            sample.stats_dfs.s3["clusters_total"] = depths.shape[0]
            try:
                sample.stats_dfs.s3["clusters_hidepth"] = int(sample.stats["clusters_hidepth"])
            except ValueError:
                ## Handle clusters_hidepth == NaN
                sample.stats_dfs.s3["clusters_hidepth"] = 0
            sample.stats_dfs.s3["avg_depth_total"] = depths.mean()
            sample.stats_dfs.s3["avg_depth_mj"] = keepmj.mean()
            sample.stats_dfs.s3["avg_depth_stat"] = keepstat.mean()
            sample.stats_dfs.s3["sd_depth_total"] = depths.std()
            sample.stats_dfs.s3["sd_depth_mj"] = keepmj.std()
            sample.stats_dfs.s3["sd_depth_stat"] = keepstat.std()

    else:
        print("    no clusters found for {}".format(sample.name))

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["assembly_method"] == "denovo":
        refmap_stats(data, sample)

    log_level = logging.getLevelName(LOGGER.getEffectiveLevel())

    if not log_level == "DEBUG":
        ## Clean up loose files only if not in DEBUG
        ##- edits/*derep, utemp, *utemp.sort, *htemp, *clust.gz
        derepfile = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")
        mergefile = os.path.join(data.dirs.edits, sample.name+"_merged_.fastq")
        uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
        usort = os.path.join(data.dirs.clusts, sample.name+".utemp.sort")
        hhandle = os.path.join(data.dirs.clusts, sample.name+".htemp")
        clusters = os.path.join(data.dirs.clusts, sample.name+".clust.gz")

        for f in [derepfile, mergefile, uhandle, usort, hhandle, clusters]:
            try:
                os.remove(f)
            except:
                pass



## winner, rigorously testing in sequential and parallel against other funcs
def persistent_popen_align3(clusts, maxseqs=200, is_gbs=False):
    """ keeps a persistent bash shell open and feeds it muscle alignments """

    ## create a separate shell for running muscle in, this is much faster
    ## than spawning a separate subprocess for each muscle call
    proc = sps.Popen(["bash"], 
                     stdin=sps.PIPE, 
                     stdout=sps.PIPE, 
                     universal_newlines=True)

    ## iterate over clusters in this file until finished
    aligned = []
    for clust in clusts:

        ## new alignment string for read1s and read2s
        align1 = ""
        align2 = ""

        ## don't bother aligning if only one seq
        if clust.count(">") == 1:
            aligned.append(clust.replace(">", "").strip())
        else:

            ## do we need to split the alignment? (is there a PE insert?)
            try:
                ## make into list (only read maxseqs lines, 2X cuz names)
                lclust = clust.split()[:maxseqs*2]

                ## try to split cluster list at nnnn separator for each read
                lclust1 = list(itertools.chain(*zip(\
                     lclust[::2], [i.split("nnnn")[0] for i in lclust[1::2]])))
                lclust2 = list(itertools.chain(*zip(\
                     lclust[::2], [i.split("nnnn")[1] for i in lclust[1::2]])))

                ## put back into strings
                clust1 = "\n".join(lclust1)
                clust2 = "\n".join(lclust2)

                ## Align the first reads.
                ## The muscle command with alignment as stdin and // as splitter
                cmd1 = "echo -e '{}' | {} -quiet -in - ; echo {}"\
                        .format(clust1, ipyrad.bins.muscle, "//")

                ## send cmd1 to the bash shell
                print(cmd1, file=proc.stdin)

                ## read the stdout by line until splitter is reached
                ## meaning that the alignment is finished.
                for line in iter(proc.stdout.readline, '//\n'):
                    align1 += line

                ## Align the second reads.
                ## The muscle command with alignment as stdin and // as splitter
                cmd2 = "echo -e '{}' | {} -quiet -in - ; echo {}"\
                        .format(clust2, ipyrad.bins.muscle, "//")

                ## send cmd2 to the bash shell
                print(cmd2, file=proc.stdin)

                ## read the stdout by line until splitter is reached
                ## meaning that the alignment is finished.
                for line in iter(proc.stdout.readline, '//\n'):
                    align2 += line

                ## join up aligned read1 and read2 and ensure names order matches
                la1 = align1[1:].split("\n>")
                la2 = align2[1:].split("\n>")
                dalign1 = dict([i.split("\n", 1) for i in la1])
                dalign2 = dict([i.split("\n", 1) for i in la2])
                align1 = []
                try:
                    keys = sorted(dalign1.keys(), key=DEREP, reverse=True)
                except ValueError as inst:
                    ## Lines is empty. This means the call to muscle alignment failed.
                    ## Not sure how to handle this, but it happens only very rarely.
                    LOGGER.error("Muscle alignment failed: Bad clust - {}\nBad lines - {}"\
                                .format(clust, lines))
                    continue

                ## put seed at top of alignment
                seed = [i for i in keys if i.split(";")[-1][0]=="*"][0]
                keys.pop(keys.index(seed))
                keys = [seed] + keys
                for key in keys:
                    align1.append("\n".join([key, 
                                    dalign1[key].replace("\n", "")+"nnnn"+\
                                    dalign2[key].replace("\n", "")]))

                ## append aligned cluster string
                aligned.append("\n".join(align1).strip())

            ## Malformed clust. Dictionary creation with only 1 element will raise.
            except ValueError as inst:
                LOGGER.debug("Bad PE cluster - {}\nla1 - {}\nla2 - {}".format(\
                                clust, la1, la2))

            ## Either reads are SE, or at least some pairs are merged.
            except IndexError:
                    
                ## limit the number of input seqs
                lclust = "\n".join(clust.split()[:maxseqs*2])

                ## the muscle command with alignment as stdin and // as splitter
                cmd = "echo -e '{}' | {} -quiet -in - ; echo {}"\
                            .format(lclust, ipyrad.bins.muscle, "//")

                ## send cmd to the bash shell (TODO: PIPE could overflow here!)
                print(cmd, file=proc.stdin)

                ## read the stdout by line until // is reached. This BLOCKS.
                for line in iter(proc.stdout.readline, '//\n'):
                    align1 += line

                ## remove '>' from names, and '\n' from inside long seqs                
                lines = align1[1:].split("\n>")

                try:
                    ## find seed of the cluster and put it on top.
                    seed = [i for i in lines if i.split(";")[-1][0]=="*"][0]
                    lines.pop(lines.index(seed))
                    lines = [seed] + sorted(lines, key=DEREP, reverse=True)
                except ValueError as inst:
                    ## Lines is empty. This means the call to muscle alignment failed.
                    ## Not sure how to handle this, but it happens only very rarely.
                    LOGGER.error("Muscle alignment failed: Bad clust - {}\nBad lines - {}"\
                                .format(clust, lines))
                    continue

                ## format remove extra newlines from muscle
                aa = [i.split("\n", 1) for i in lines]
                align1 = [i[0]+'\n'+"".join([j.replace("\n", "") for j in i[1:]]) for i in aa]
                
                ## trim edges in sloppy gbs/ezrad data. Maybe relevant to other types too...
                if is_gbs:
                    align1 = gbs_trim(align1)

                ## append to aligned
                aligned.append("\n".join(align1).strip())
               
    # cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    ## return the aligned clusters
    return aligned   


def gbs_trim(align1):
    """
    No reads can go past the left of the seed, or right of the least extended
    reverse complement match. Example below. m is a match. u is an area where 
    lots of mismatches typically occur. The cut sites are shown.
    
    Original locus*
    Seed           TGCAG************************************-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm------------------------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuuuuuuuu
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu

    Trimmed locus*
    Seed           TGCAG************************************---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmm
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    """
    leftmost = rightmost = None
    dd = {k:v for k,v in [j.rsplit("\n", 1) for j in align1]}
    seed = [i for i in dd.keys() if i.rsplit(";")[-1][0] == "*"][0]
    leftmost = [i != "-" for i in dd[seed]].index(True)
    revs = [i for i in dd.keys() if i.rsplit(";")[-1][0] == "-"]
    if revs:
        subright = max([[i!="-" for i in seq[::-1]].index(True) \
            for seq in [dd[i] for i in revs]])
    else:
        subright = 0
    rightmost = len(dd[seed]) - subright

    ## if locus got clobbered then print place-holder NNN
    names, seqs = zip(*[i.rsplit("\n", 1) for i in align1])
    if rightmost > leftmost:
        newalign1 = [n+"\n"+i[leftmost:rightmost] for i,n in zip(seqs, names)]
    else:
        newalign1 = [n+"\nNNN" for i,n in zip(seqs, names)]
    return newalign1


## quick lambda func to get derep number from reads
DEREP = lambda x: int(x.split("=")[-1].split(";")[0])


## max-internal-indels could be modified if we add it to hackerz dict.
def align_and_parse(handle, max_internal_indels=5, is_gbs=False):
    """ much faster implementation for aligning chunks """

    ## data are already chunked, read in the whole thing. bail if no data.
    try:
        with open(handle, 'rb') as infile:
            clusts = infile.read().split("//\n//\n")
            ## remove any empty spots
            clusts = [i for i in clusts if i]
            ## Skip entirely empty chunks
            if not clusts:
                raise IPyradError
    except (IOError, IPyradError):
        LOGGER.debug("skipping empty chunk - {}".format(handle))
        return 0

    ## count discarded clusters for printing to stats later
    highindels = 0

    ## iterate over clusters sending each to muscle, splits and aligns pairs
    try:
        aligned = persistent_popen_align3(clusts, 200, is_gbs)
    except Exception as inst:
        LOGGER.debug("Error in handle - {} - {}".format(handle, inst))
        #raise IPyradWarningExit("error hrere {}".format(inst))
        aligned = []        

    ## store good alignments to be written to file
    refined = []

    ## filter and trim alignments
    for clust in aligned:

        ## check for too many internal indels
        filtered = aligned_indel_filter(clust, max_internal_indels)

        ## reverse complement matches. No longer implemented.
        #filtered = overshoot_filter(clust)

        ## finally, add to outstack if alignment is good
        if not filtered:
            refined.append(clust)#"\n".join(stack))
        else:
            highindels += 1

    ## write to file after
    if refined:
        outhandle = handle.rsplit(".", 1)[0]+".aligned"
        with open(outhandle, 'wb') as outfile:
            outfile.write("\n//\n//\n".join(refined)+"\n")

    ## remove the old tmp file
    log_level = logging.getLevelName(LOGGER.getEffectiveLevel())
    if not log_level == "DEBUG":
        os.remove(handle)
    return highindels



def aligned_indel_filter(clust, max_internal_indels):
    """ checks for too many internal indels in muscle aligned clusters """

    ## make into list
    lclust = clust.split()
    
    ## paired or not
    try:
        seq1 = [i.split("nnnn")[0] for i in lclust[1::2]]
        seq2 = [i.split("nnnn")[1] for i in lclust[1::2]]
        intindels1 = [i.rstrip("-").lstrip("-").count("-") for i in seq1]
        intindels2 = [i.rstrip("-").lstrip("-").count("-") for i in seq2]
        intindels = intindels1 + intindels2
        if max(intindels) > max_internal_indels:
            return 1
       
    except IndexError:
        seq1 = lclust[1::2]
        intindels = [i.rstrip("-").lstrip("-").count("-") for i in seq1]
        if max(intindels) > max_internal_indels:
            return 1 
    
    return 0



def build_clusters(data, sample, maxindels):
    """
    Combines information from .utemp and .htemp files to create .clust files,
    which contain un-aligned clusters. Hits to seeds are only kept in the
    cluster if the number of internal indels is less than 'maxindels'.
    By default, we set maxindels=6 for this step (within-sample clustering).
    """

    ## If reference assembly then here we're clustering the unmapped reads
    if "reference" in data.paramsdict["assembly_method"]:
        derepfile = os.path.join(data.dirs.edits, sample.name+"-refmap_derep.fastq")
    else:
        derepfile = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")
    ## i/o vsearch files
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    usort = os.path.join(data.dirs.clusts, sample.name+".utemp.sort")
    hhandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## create an output file to write clusters to
    sample.files.clusters = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
    clustsout = gzip.open(sample.files.clusters, 'wb')

    ## Sort the uhandle file so we can read through matches efficiently
    cmd = ["sort", "-k", "2", uhandle, "-o", usort]
    proc = sps.Popen(cmd, close_fds=True)
    _ = proc.communicate()[0]

    ## load ALL derep reads into a dictionary (this can be a few GB of RAM)
    ## and is larger if names are larger. We are grabbing two lines at a time.
    alldereps = {}
    with open(derepfile, 'rb') as ioderep:
        dereps = itertools.izip(*[iter(ioderep)]*2)
        for namestr, seq in dereps:
            nnn, sss = [i.strip() for i in namestr, seq]
            alldereps[nnn[1:]] = sss

    ## store observed seeds (this could count up to >million in bad data sets)
    seedsseen = set()

    ## Iterate through the usort file grabbing matches to build clusters
    with open(usort, 'rb') as insort:
        ## iterator, seed null, seqlist null
        isort = iter(insort)
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0
        while 1:
            ## grab the next line
            try:
                hit, seed, _, ind, ori, _ = isort.next().strip().split()
            except StopIteration:
                break

            ## same seed, append match
            if seed != lastseed:
                seedsseen.add(seed)
                ## store the last cluster (fseq), count it, and clear fseq
                if fseqs:
                    ## sort fseqs by derep after pulling out the seed
                    fseqs = [fseqs[0]] + sorted(fseqs[1:], key=lambda x: \
                        int(x.split(";size=")[1].split(";")[0]), reverse=True)                    
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                ## occasionally write/dump stored clusters to file and clear mem
                if not seqsize % 10000:
                    if seqlist:
                        clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")
                        ## reset list and counter
                        seqlist = []

                ## store the new seed on top of fseq list
                fseqs.append(">{}*\n{}".format(seed, alldereps[seed]))
                lastseed = seed

            ## add match to the seed
            ## revcomp if orientation is reversed (comp preserves nnnn)
            if ori == "-":
                seq = comp(alldereps[hit])[::-1]
            else:
                seq = alldereps[hit]
            ## only save if not too many indels
            if int(ind) <= maxindels:
                fseqs.append(">{}{}\n{}".format(hit, ori, seq))
            else:
                LOGGER.info("filtered by maxindels: %s %s", ind, seq)

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
    if seqlist:
        clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")

    ## now write the seeds that had no hits. Make dict from htemp
    with open(hhandle, 'rb') as iotemp:
        nohits = itertools.izip(*[iter(iotemp)]*2)
        seqlist = []
        seqsize = 0
        while 1:
            try:
                nnn, _ = [i.strip() for i in nohits.next()]
            except StopIteration:
                break

            ## occasionally write to file
            if not seqsize % 10000:
                if seqlist:
                    clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")
                    ## reset list and counter
                    seqlist = []

            ## append to list if new seed
            if nnn[1:] not in seedsseen:
                seqlist.append("{}*\n{}".format(nnn, alldereps[nnn[1:]]))
                seqsize += 1

    ## write whatever is left over to the clusts file
    if seqlist:
        clustsout.write("\n//\n//\n".join(seqlist))#+"\n//\n//\n")

    ## close the file handle
    clustsout.close()
    del alldereps



def setup_dirs(data):
    """ sets up directories for step3 data """
    ## make output folder for clusters
    pdir = os.path.realpath(data.paramsdict["project_dir"])
    data.dirs.clusts = os.path.join(pdir, "{}_clust_{}"\
                       .format(data.name, data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.mkdir(data.dirs.clusts)

    ## make a tmpdir for align files
    data.tmpdir = os.path.abspath(os.path.expanduser(
        os.path.join(pdir, data.name+'-tmpalign')))
    if not os.path.exists(data.tmpdir):
        os.mkdir(data.tmpdir)

    ## If ref mapping, init samples and make the refmapping output directory.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(pdir, "{}_refmapping".format(data.name))
        if not os.path.exists(data.dirs.refmapping):
            os.mkdir(data.dirs.refmapping)



def new_apply_jobs(data, samples, ipyclient, nthreads, maxindels, force):
    """
    Create a DAG of prealign jobs to be run in order for each sample. Track
    Progress, report errors. Each assembly method has a slightly different
    DAG setup, calling different functions.
    """

    ## is datatype gbs? used in alignment-trimming by align_and_parse()
    is_gbs = bool("gbs" in data.paramsdict["datatype"])

    ## Two view objects, threaded and unthreaded
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    firstfunc = "derep_concat_split"
    printstr = " {}    | {} | s3 |".format(PRINTSTR[firstfunc], elapsed)
    #printstr = " {}      | {} | s3 |".format(PRINTSTR[], elapsed)
    progressbar(10, 0, printstr, spacer=data._spacer)

    ## TODO: for HPC systems this should be done to make sure targets are spread
    ## among different nodes.
    if nthreads:
        if nthreads < len(ipyclient.ids):
            thview = ipyclient.load_balanced_view(targets=ipyclient.ids[::nthreads])
        elif nthreads == 1:
            thview = ipyclient.load_balanced_view()
        else:
            if len(ipyclient) > 40:
                thview = ipyclient.load_balanced_view(targets=ipyclient.ids[::4])
            else:
                thview = ipyclient.load_balanced_view(targets=ipyclient.ids[::2])


    ## get list of jobs/dependencies as a DAG for all pre-align funcs.
    dag, joborder = build_dag(data, samples)

    ## dicts for storing submitted jobs and results
    results = {}

    ## submit jobs to the engines in single or threaded views. The topological
    ## sort makes sure jobs are input with all dependencies found.
    for node in nx.topological_sort(dag):
        ## get list of async results leading to this job
        deps = [results.get(n) for n in dag.predecessors(node)]
        deps = ipp.Dependency(dependencies=deps, failure=True)

        ## get func, sample, and args for this func (including [data, sample])
        funcstr, chunk, sname = node.split("-", 2)
        func = FUNCDICT[funcstr]
        sample = data.samples[sname]

        ## args vary depending on the function
        if funcstr in ["derep_concat_split", "cluster"]:
            args = [data, sample, nthreads, force]
        elif funcstr in ["mapreads"]:
            args = [data, sample, nthreads, force]
        elif funcstr in ["build_clusters"]:
            args = [data, sample, maxindels]
        elif funcstr in ["muscle_align"]:
            handle = os.path.join(data.tmpdir, 
                        "{}_chunk_{}.ali".format(sample.name, chunk))
            args = [handle, maxindels, is_gbs]
        else:
            args = [data, sample]

        # submit and store AsyncResult object. Some jobs are threaded.
        if nthreads and (funcstr in THREADED_FUNCS):
            #LOGGER.info('submitting %s to %s-threaded view', funcstr, nthreads)
            with thview.temp_flags(after=deps, block=False):
                results[node] = thview.apply(func, *args)
        else:
            #LOGGER.info('submitting %s to single-threaded view', funcstr)
            with lbview.temp_flags(after=deps, block=False):
                results[node] = lbview.apply(func, *args)

    ## track jobs as they finish, abort if someone fails. This blocks here
    ## until all jobs are done. Keep track of which samples have failed so
    ## we only print the first error message.
    sfailed = set()
    for funcstr in joborder + ["muscle_align", "reconcat"]:
        errfunc, sfails, msgs = trackjobs(funcstr, results, spacer=data._spacer)
        LOGGER.info("{}-{}-{}".format(errfunc, sfails, msgs))
        if errfunc:
            for sidx in xrange(len(sfails)):
                sname = sfails[sidx]
                errmsg = msgs[sidx]
                if sname not in sfailed:
                    print("  sample [{}] failed. See error in ./ipyrad_log.txt"\
                          .format(sname))
                    LOGGER.error("sample [%s] failed in step [%s]; error: %s",
                                  sname, errfunc, errmsg)
                    sfailed.add(sname)

    ## Cleanup of successful samples, skip over failed samples
    badaligns = {}
    for sample in samples:
        ## The muscle_align step returns the number of excluded bad alignments
        for async in results:
            func, chunk, sname = async.split("-", 2)
            if (func == "muscle_align") and (sname == sample.name):
                if results[async].successful():
                    badaligns[sample] = int(results[async].get())

    ## for the samples that were successful:
    for sample in badaligns:
        ## store the result
        sample.stats_dfs.s3.filtered_bad_align = badaligns[sample]
        ## store all results
        try:
            sample_cleanup(data, sample)
        except Exception as inst:
            msg = "  Sample {} failed this step. See ipyrad_log.txt.\
                  ".format(sample.name)
            print(msg)
            LOGGER.error("%s - %s", sample.name, inst)

    ## store the results to data
    data_cleanup(data)

    ## uncomment to plot the dag
    #_plot_dag(dag, results, snames)



def build_dag(data, samples):
    """
    build a directed acyclic graph describing jobs to be run in order.
    """

    ## Create DAGs for the assembly method being used, store jobs in nodes
    snames = [i.name for i in samples]
    dag = nx.DiGraph()

    ## get list of pre-align jobs from globals based on assembly method
    joborder = JOBORDER[data.paramsdict["assembly_method"]]

    ## WHICH JOBS TO RUN: iterate over the sample names
    for sname in snames:
        ## append pre-align job for each sample to nodes list
        for func in joborder:
            dag.add_node("{}-{}-{}".format(func, 0, sname))

        ## append align func jobs, each will have max 10
        for chunk in xrange(10):
            dag.add_node("{}-{}-{}".format("muscle_align", chunk, sname))

        ## append final reconcat jobs
        dag.add_node("{}-{}-{}".format("reconcat", 0, sname))

    ## ORDER OF JOBS: add edges/dependency between jobs: (first-this, then-that)
    for sname in snames:
        for sname2 in snames:
            ## enforce that clust/map cannot start until derep is done for ALL
            ## samples. This is b/c...
            dag.add_edge("{}-{}-{}".format(joborder[0], 0, sname2),
                         "{}-{}-{}".format(joborder[1], 0, sname))

        ## add remaining pre-align jobs 
        for idx in xrange(2, len(joborder)):
            dag.add_edge("{}-{}-{}".format(joborder[idx-1], 0, sname),
                         "{}-{}-{}".format(joborder[idx], 0, sname))

        ## Add 10 align jobs, none of which can start until all chunker jobs
        ## are finished. Similarly, reconcat jobs cannot start until all align
        ## jobs are finished.
        for sname2 in snames:
            for chunk in range(10):
                dag.add_edge("{}-{}-{}".format("muscle_chunker", 0, sname2),
                             "{}-{}-{}".format("muscle_align", chunk, sname))
                ## add that the final reconcat job can't start until after
                ## each chunk of its own sample has finished aligning.
                dag.add_edge("{}-{}-{}".format("muscle_align", chunk, sname),
                             "{}-{}-{}".format("reconcat", 0, sname))
    ## return the dag
    return dag, joborder



def _plot_dag(dag, results, snames):
    """
    makes plot to help visualize the DAG setup. For developers only.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.dates import date2num
        from matplotlib.cm import gist_rainbow

        ## first figure is dag layout
        plt.figure("dag_layout", figsize=(10, 10))
        nx.draw(dag,
                pos=nx.spring_layout(dag),
                node_color='pink',
                with_labels=True)
        plt.savefig("./dag_layout.png", bbox_inches='tight', dpi=200)

        ## second figure is times for steps
        pos = {}
        colors = {}

        for node in dag:
            #jobkey = "{}-{}".format(node, sample)
            mtd = results[node].metadata
            start = date2num(mtd.started)
            #runtime = date2num(md.completed)# - start
            ## sample id to separate samples on x-axis
            _, _, sname = node.split("-", 2)
            sid = snames.index(sname)
            ## 1e6 to separate on y-axis
            pos[node] = (start+sid, start*1e6)
            colors[node] = mtd.engine_id

        ## x just spaces out samples;
        ## y is start time of each job with edge leading to next job
        ## color is the engine that ran the job
        ## all jobs were submitted as 3 second wait times
        plt.figure("dag_starttimes", figsize=(10, 16))
        nx.draw(dag, pos,
                node_list=colors.keys(),
                node_color=colors.values(),
                cmap=gist_rainbow,
                with_labels=True)
        plt.savefig("./dag_starttimes.png", bbox_inches='tight', dpi=200)

    except Exception as inst:
        LOGGER.warning(inst)



def trackjobs(func, results, spacer):
    """
    Blocks and prints progress for just the func being requested from a list
    of submitted engine jobs. Returns whether any of the jobs failed.

    func = str
    results = dict of asyncs
    """

    ## TODO: try to insert a better way to break on KBD here.
    LOGGER.info("inside trackjobs of %s", func)

    ## get just the jobs from results that are relevant to this func
    asyncs = [(i, results[i]) for i in results if i.split("-", 2)[0] == func]

    ## progress bar
    start = time.time()
    while 1:
        ## how many of this func have finished so far
        ready = [i[1].ready() for i in asyncs]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        printstr = " {}    | {} | s3 |".format(PRINTSTR[func], elapsed)
        progressbar(len(ready), sum(ready), printstr, spacer=spacer)
        time.sleep(0.1)
        if len(ready) == sum(ready):
            print("")
            break

    sfails = []
    errmsgs = []
    for job in asyncs:
        if not job[1].successful():
            sfails.append(job[0])
            errmsgs.append(job[1].result())

    return func, sfails, errmsgs

    ## did any samples fail?
    #success = [i[1].successful() for i in asyncs]

    ## return functionstring and error message on failure
    #if not all(success):
    #    ## get error messages
    #    errmsgs = [i[1].exception() for i in asyncs if not i[1].successful()]
    #    ## get samlpes that failed
    #    sfails = [i[0].split("-", 2)[-1] for i in asyncs if not i[1].successful()]
    #    return func, sfails, errmsgs
    #else:
    #    return 0, [], []



def declone_3rad(data, sample):
    """
    3rad uses random adapters to identify pcr duplicates. We will
    remove pcr dupes here. Basically append the radom adapter to
    each sequence, do a regular old vsearch derep, then trim
    off the adapter, and push it down the pipeline. This will
    remove all identical seqs with identical random i5 adapters.
    """

    LOGGER.info("Entering declone_3rad - {}".format(sample.name))

    ## Append i5 adapter to the head of each read. Merged file is input, and
    ## still has fq qual score so also have to append several qscores for the
    ## adapter bases. Open the merge file, get quarts, go through each read
    ## and append the necessary stuff.

    adapter_seqs_file = tempfile.NamedTemporaryFile(mode='wb',
                                        delete=False,
                                        dir=data.dirs.edits,
                                        suffix="_append_adapters_.fastq")

    try:
        with open(sample.files.edits[0][0]) as infile:
            quarts = itertools.izip(*[iter(infile)]*4)

            ## a list to store until writing
            writing = []
            counts = 0

            while 1:
                try:
                    read = quarts.next()
                except StopIteration:
                    break

                ## Split on +, get [1], split on "_" (can be either _r1 or
                ## _m1 if merged reads) and get [0] for the i5
                ## prepend "EEEEEEEE" as qscore for the adapters
                i5 = read[0].split("+")[1].split("_")[0]

                ## If any non ACGT in the i5 then drop this sequence
                if 'N' in i5:
                    continue
                writing.append("\n".join([
                                read[0].strip(),
                                i5 + read[1].strip(),
                                read[2].strip(),
                                "E"*8 + read[3].strip()]
                            ))

                ## Write the data in chunks
                counts += 1
                if not counts % 1000:
                    adapter_seqs_file.write("\n".join(writing)+"\n")
                    writing = []
            if writing:
                adapter_seqs_file.write("\n".join(writing))
                adapter_seqs_file.close()

        tmp_outfile = tempfile.NamedTemporaryFile(mode='wb',
                                        delete=False,
                                        dir=data.dirs.edits,
                                        suffix="_decloned_w_adapters_.fastq")

        ## Close the tmp file bcz vsearch will write to it by name, then
        ## we will want to reopen it to read from it.
        tmp_outfile.close()
        ## Derep the data (adapters+seq)
        derep_and_sort(data, adapter_seqs_file.name,
                       os.path.join(data.dirs.edits, tmp_outfile.name), 2)

        ## Remove adapters from head of sequence and write out
        ## tmp_outfile is now the input file for the next step
        ## first vsearch derep discards the qscore so we iterate
        ## by pairs
        with open(tmp_outfile.name) as infile:
            with open(os.path.join(data.dirs.edits, sample.name+"_declone.fastq"),\
                                'wb') as outfile:
                duo = itertools.izip(*[iter(infile)]*2)

                ## a list to store until writing
                writing = []
                counts2 = 0

                while 1:
                    try:
                        read = duo.next()
                    except StopIteration:
                        break

                    ## Peel off the adapters. There's probably a faster
                    ## way of doing this.
                    writing.append("\n".join([
                                    read[0].strip(),
                                    read[1].strip()[8:]]
                                ))

                    ## Write the data in chunks
                    counts2 += 1
                    if not counts2 % 1000:
                        outfile.write("\n".join(writing)+"\n")
                        writing = []
                if writing:
                    outfile.write("\n".join(writing))
                    outfile.close()

        LOGGER.info("Removed pcr duplicates from {} - {}".format(sample.name, counts-counts2))

    except Exception as inst:
        raise IPyradError("    Caught error while decloning "\
                                + "3rad data - {}".format(inst))

    finally:
        ## failed samples will cause tmp file removal to raise.
        ## just ignore it.
        try:
            ## Clean up temp files
            if os.path.exists(adapter_seqs_file.name):
                os.remove(adapter_seqs_file.name)
            if os.path.exists(tmp_outfile.name):
                os.remove(tmp_outfile.name)
        except Exception as inst:
            pass



def derep_and_sort(data, infile, outfile, nthreads):
    """
    Dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file. Paired
    reads are dereplicated as one concatenated read and later split again.
    Updated this function to take infile and outfile to support the double
    dereplication that we need for 3rad (5/29/15 iao).
    """

    ## datatypes options
    strand = "plus"
    if "gbs" in data.paramsdict["datatype"]\
        or "2brad" in data.paramsdict["datatype"]:
        strand = "both"

    ## pipe in a gzipped file
    if infile.endswith(".gz"):
        catcmd = ["gunzip", "-c", infile]
    else:
        catcmd = ["cat", infile]

    ## do dereplication with vsearch
    cmd = [ipyrad.bins.vsearch,
            "--derep_fulllength", "-",
            "--strand", strand,
            "--output", outfile,
            "--threads", str(nthreads),
            "--fasta_width", str(0),
            "--fastq_qmax", "1000",
            "--sizeout", 
            "--relabel_md5",
            ]
    LOGGER.info("derep cmd %s", cmd)

    ## run vsearch
    proc1 = sps.Popen(catcmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    proc2 = sps.Popen(cmd, stdin=proc1.stdout, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)

    try:
        errmsg = proc2.communicate()[0]
    except KeyboardInterrupt:
        LOGGER.info("interrupted during dereplication")
        raise KeyboardInterrupt()

    if proc2.returncode:
        LOGGER.error("error inside derep_and_sort %s", errmsg)
        raise IPyradWarningExit(errmsg)



def data_cleanup(data):
    """ cleanup / statswriting function for Assembly obj """
    data.stats_dfs.s3 = data._build_stat("s3")
    data.stats_files.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    with io.open(data.stats_files.s3, 'w') as outfile:
        data.stats_dfs.s3.to_string(
            buf=outfile,
            formatters={
                'merged_pairs':'{:.0f}'.format,
                'clusters_total':'{:.0f}'.format,
                'clusters_hidepth':'{:.0f}'.format,
                'filtered_bad_align':'{:.0f}'.format,
                'avg_depth_stat':'{:.2f}'.format,
                'avg_depth_mj':'{:.2f}'.format,
                'avg_depth_total':'{:.2f}'.format,
                'sd_depth_stat':'{:.2f}'.format,
                'sd_depth_mj':'{:.2f}'.format,
                'sd_depth_total':'{:.2f}'.format
            })



def concat_multiple_edits(data, sample):
    """
    if multiple fastq files were appended into the list of fastqs for samples
    then we merge them here before proceeding.
    """

    ## if more than one tuple in fastq list
    if len(sample.files.edits) > 1:
        ## create a cat command to append them all (doesn't matter if they
        ## are gzipped, cat still works). Grab index 0 of tuples for R1s.
        cmd1 = ["cat"] + [i[0] for i in sample.files.edits]

        ## write to new concat handle
        conc1 = os.path.join(data.dirs.edits, sample.name+"_R1_concatedit.fq.gz")
        with open(conc1, 'w') as cout1:
            proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=cout1, close_fds=True)
            res1 = proc1.communicate()[0]
        if proc1.returncode:
            raise IPyradWarningExit("error in: %s, %s", cmd1, res1)

        ## Only set conc2 if R2 actually exists
        conc2 = 0
        if os.path.exists(str(sample.files.edits[0][1])):
            cmd2 = ["cat"] + [i[1] for i in sample.files.edits]
            conc2 = os.path.join(data.dirs.edits, sample.name+"_R2_concatedit.fq.gz")
            with gzip.open(conc2, 'w') as cout2:
                proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=cout2, close_fds=True)
                res2 = proc2.communicate()[0]
            if proc2.returncode:
                raise IPyradWarningExit("error in: %s, %s", cmd2, res2)

        ## store new file handles
        sample.files.edits = [(conc1, conc2)]
    return sample.files.edits



def cluster(data, sample, nthreads, force):
    """
    Calls vsearch for clustering. cov varies by data type, values were chosen
    based on experience, but could be edited by users
    """

    ## get the dereplicated reads
    if "reference" in data.paramsdict["assembly_method"]:
        derephandle = os.path.join(data.dirs.edits, sample.name+"-refmap_derep.fastq")
        ## In the event all reads for all samples map successfully then clustering
        ## the unmapped reads makes no sense, so just bail out.
        if not os.stat(derephandle).st_size:
            ## In this case you do have to create empty, dummy vsearch output
            ## files so building_clusters will not fail.
            uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
            usort = os.path.join(data.dirs.clusts, sample.name+".utemp.sort")
            hhandle = os.path.join(data.dirs.clusts, sample.name+".htemp")
            for f in [uhandle, usort, hhandle]:
                open(f, 'a').close()
            return
    else:
        derephandle = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")

    ## create handles for the outfiles
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## If derep file doesn't exist then bail out
    if not os.path.isfile(derephandle):
        LOGGER.warn("Bad derephandle - {}".format(derephandle))
        raise IPyradError("Input file for clustering doesn't exist - {}"\
                        .format(derephandle))

    ## testing one sample fail
    #if sample.name == "1C_0":
    #    x

    ## datatype specific optimization
    ## minsl: the percentage of the seed that must be matched
    ##    smaller values for RAD/ddRAD where we might want to combine, say 50bp
    ##    reads and 100bp reads in the same analysis.
    ## query_cov: the percentage of the query sequence that must match seed
    ##    smaller values are needed for gbs where only the tips might overlap
    ##    larger values for pairgbs where they should overlap near completely
    ##    small minsl and high query cov allows trimmed reads to match to untrim
    ##    seed for rad/ddrad/pairddrad.
    strand = "plus"
    cov = 0.75
    minsl = 0.5
    if data.paramsdict["datatype"] in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.5
        minsl = 0.5
    elif data.paramsdict["datatype"] == 'pairgbs':
        strand = "both"
        cov = 0.75
        minsl = 0.75

    ## If this value is not null (which is the default) then override query cov
    if data._hackersonly["query_cov"]:
        cov = str(data._hackersonly["query_cov"])
        assert float(cov) <= 1, "query_cov must be <= 1.0"

    ## get call string
    cmd = [ipyrad.bins.vsearch,
           "-cluster_smallmem", derephandle,
           "-strand", strand,
           "-query_cov", str(cov),
           "-id", str(data.paramsdict["clust_threshold"]),
           "-minsl", str(minsl),
           "-userout", uhandle,
           "-userfields", "query+target+id+gaps+qstrand+qcov",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-threads", str(nthreads),
           "-notmatched", temphandle,
           "-fasta_width", "0",
           "-fastq_qmax", "100",
           "-fulldp",
           "-usersort"]

    ## not sure what the benefit of this option is exactly, needs testing,
    ## might improve indel detection on left side, but we don't want to enforce
    ## aligning on left side if not necessarily, since quality trimmed reads
    ## might lose bases on left side in step2 and no longer align.
    #if data.paramsdict["datatype"] in ["rad", "ddrad", "pairddrad"]:
    #    cmd += ["-leftjust"]

    ## run vsearch
    LOGGER.debug("%s", cmd)
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)

    ## This is long running so we wrap it to make sure we can kill it
    try:
        res = proc.communicate()[0]
    except KeyboardInterrupt:
        proc.kill()
        raise KeyboardInterrupt

    ## check for errors
    if proc.returncode:
        LOGGER.error("error %s: %s", cmd, res)
        raise IPyradWarningExit("cmd {}: {}".format(cmd, res))



def muscle_chunker(data, sample):
    """
    Splits the muscle alignment into chunks. Each chunk is run on a separate
    computing core. Because the largest clusters are at the beginning of the 
    clusters file, assigning equal clusters to each file would put all of the 
    large cluster, that take longer to align, near the top. So instead we 
    randomly distribute the clusters among the files. If assembly method is
    reference then this step is just a placeholder and nothing happens. 
    """
    ## log our location for debugging
    LOGGER.info("inside muscle_chunker")

    ## only chunk up denovo data, refdata has its own chunking method which 
    ## makes equal size chunks, instead of uneven chunks like in denovo
    if data.paramsdict["assembly_method"] != "reference":
        ## get the number of clusters
        clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
        with iter(gzip.open(clustfile, 'rb')) as clustio:
            nloci = sum(1 for i in clustio if "//" in i) // 2
            #tclust = clustio.read().count("//")//2
            optim = (nloci//20) + (nloci%20)
            LOGGER.info("optim for align chunks: %s", optim)

        ## write optim clusters to each tmp file
        clustio = gzip.open(clustfile, 'rb')
        inclusts = iter(clustio.read().strip().split("//\n//\n"))
        
        ## splitting loci so first file is smaller and last file is bigger
        inc = optim // 10
        for idx in range(10):
            ## how big is this chunk?
            this = optim + (idx * inc)
            left = nloci-this
            if idx == 9:
                ## grab everything left
                grabchunk = list(itertools.islice(inclusts, int(1e9)))
            else:
                ## grab next chunks-worth of data
                grabchunk = list(itertools.islice(inclusts, this))
                nloci = left

            ## write the chunk to file
            tmpfile = os.path.join(data.tmpdir, sample.name+"_chunk_{}.ali".format(idx))
            with open(tmpfile, 'wb') as out:
                out.write("//\n//\n".join(grabchunk))

        ## write the chunk to file
        #grabchunk = list(itertools.islice(inclusts, left))
        #if grabchunk:
        #    tmpfile = os.path.join(data.tmpdir, sample.name+"_chunk_9.ali")
        #    with open(tmpfile, 'a') as out:
        #        out.write("\n//\n//\n".join(grabchunk))
        clustio.close()



def reconcat(data, sample):
    """ takes aligned chunks (usually 10) and concatenates them """

    try:
        ## get chunks
        chunks = glob.glob(os.path.join(data.tmpdir,
                 sample.name+"_chunk_[0-9].aligned"))

        ## sort by chunk number, cuts off last 8 =(aligned)
        chunks.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))
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
    except Exception as inst:
        LOGGER.error("Error in reconcat {}".format(inst))
        raise



def derep_concat_split(data, sample, nthreads, force):
    """
    Running on remote Engine. Refmaps, then merges, then dereplicates,
    then denovo clusters reads.
    """

    ## report location for debugging
    LOGGER.info("INSIDE derep %s", sample.name)

    ## MERGED ASSEMBIES ONLY:
    ## concatenate edits files within Samples. Returns a new sample.files.edits 
    ## with the concat file. No change if not merged Assembly.
    mergefile = os.path.join(data.dirs.edits, sample.name+"_merged_.fastq")
    if not force:
        if not os.path.exists(mergefile):
            sample.files.edits = concat_multiple_edits(data, sample)
        else:
            LOGGER.info("skipped concat_multiple_edits: {} exists"\
                        .format(mergefile))
    else:
        sample.files.edits = concat_multiple_edits(data, sample)

    ## PAIRED DATA ONLY:
    ## Denovo: merge or concat fastq pairs [sample.files.pairs]
    ## Reference: only concat fastq pairs  []
    ## Denovo + Reference: ...
    if 'pair' in data.paramsdict['datatype']:
        ## the output file handle for merged reads
        

        ## modify behavior of merging vs concating if reference
        if "reference" in data.paramsdict["assembly_method"]:
            nmerged = merge_pairs(data, sample.files.edits, mergefile, 0, 0)
        else:
            nmerged = merge_pairs(data, sample.files.edits, mergefile, 1, 1)

        ## store results
        sample.files.edits = [(mergefile, )]
        sample.stats.reads_merged = nmerged

    ## 3rad uses random adapters to identify pcr duplicates. We will
    ## remove pcr dupes here. Basically append the radom adapter to
    ## each sequence, do a regular old vsearch derep, then trim
    ## off the adapter, and push it down the pipeline. This will
    ## remove all identical seqs with identical random i5 adapters.
    if "3rad" in data.paramsdict["datatype"]:
        declone_3rad(data, sample)
        derep_and_sort(data,
                os.path.join(data.dirs.edits, sample.name+"_declone.fastq"),
                os.path.join(data.dirs.edits, sample.name+"_derep.fastq"),
                nthreads)
    else:
        ## convert fastq to fasta, then derep and sort reads by their size.
        ## we pass in only one file b/c paired should be merged by now.
        derep_and_sort(data,
                sample.files.edits[0][0],
                os.path.join(data.dirs.edits, sample.name+"_derep.fastq"),
                nthreads)



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



def run(data, samples, noreverse, maxindels, force, preview, ipyclient):
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
                if sample.stats.reads_passed_filter:
                    subsamples.append(sample)
        else:
            ## force to overwrite
            if sample.stats.reads_passed_filter:
                subsamples.append(sample)

    ## run subsamples
    if not subsamples:
        print("  No Samples ready to be clustered. First run step2().")

    else:
        ## arguments to apply_jobs, inst catches exceptions
        try:
            ## make dirs that are needed including tmpdir
            setup_dirs(data)

            ## if refmapping make filehandles that will be persistent
            if not data.paramsdict["assembly_method"] == "denovo":
                for sample in subsamples:
                    refmap_init(data, sample, force)

                    ## set thread-count to 2 for paired-data
                    nthreads = 2
            ## set thread-count to 1 for single-end data          
            else:
                nthreads = 1

            ## overwrite nthreads if value in _ipcluster dict
            if "threads" in data._ipcluster.keys():
                nthreads = int(data._ipcluster["threads"])

                ## if more CPUs than there are samples then increase threads
                _ncpus = len(ipyclient)
                if _ncpus > 2*len(data.samples):
                    nthreads *= 2

            ## submit jobs to be run on cluster
            args = [data, subsamples, ipyclient, nthreads, maxindels, force]
            new_apply_jobs(*args)


        finally:
            ## this can fail if jobs were not stopped properly and are still
            ## writing to tmpdir. don't cleanup if debug is on.
            try:
                log_level = logging.getLevelName(LOGGER.getEffectiveLevel())
                if not log_level == "DEBUG":

                    if os.path.exists(data.tmpdir):
                        shutil.rmtree(data.tmpdir)
                    ## get all refmap_derep.fastqs
                    rdereps = glob.glob(os.path.join(data.dirs.edits, "*-refmap_derep.fastq"))
                    ## Remove the unmapped fastq files
                    for rmfile in rdereps:
                        os.remove(rmfile)

            except Exception as _:
                LOGGER.warning("failed to cleanup files/dirs")




### GLOBALS

THREADED_FUNCS = ["derep_concat_split", "cluster", "mapreads"]

PRINTSTR = {
    #"derep_concat_split" : "concat+dereplicate",
    "derep_concat_split" : "dereplicating     ",
    "mapreads" :           "mapping           ",
    "cluster" :            "clustering        ",
    "build_clusters" :     "building clusters ",
    "ref_muscle_chunker" : "finalize mapping  ",
    "muscle_chunker" :     "chunking          ",
    "muscle_align" :       "aligning          ",
    "reconcat" :           "concatenating     ",
    "ref_build_and_muscle_chunk" : 
                           "fetch mapped reads",
    }


FUNCDICT = {
    "derep_concat_split" : derep_concat_split,
    "mapreads" :           mapreads,
    "cluster" :            cluster,
    "build_clusters" :     build_clusters,
    "ref_muscle_chunker" : ref_muscle_chunker,
    "ref_build_and_muscle_chunk" : ref_build_and_muscle_chunk,    
    "muscle_chunker" :     muscle_chunker,
    "muscle_align" :       align_and_parse, #muscle_align,
    "reconcat" :           reconcat
    }


## Pre-align funcs for the four assembly methods
JOBORDER = {
    "denovo" : [
        "derep_concat_split",
        "cluster",
        "build_clusters",
        "muscle_chunker"
        ], 
    "reference" : [
        "derep_concat_split",
        "mapreads",
        "ref_build_and_muscle_chunk",
        "muscle_chunker",  ## <- doesn't do anything but hold back aligning jobs
        #"ref_muscle_chunker",
        ], 
    "denovo+reference" : [
        "derep_concat_split",
        "mapreads",
        "cluster",
        "build_clusters",
        "ref_build_and_muscle_chunk",
        "muscle_chunker"
        ], 
    "denovo-reference" : [
        "derep_concat_split",
        "mapreads",
        "cluster",
        "build_clusters",
        "muscle_chunker"
        ],
    }   


NO_UHITS_ERROR = """\
    No clusters (.utemp hits) found for {}. If you are running preview mode and
    the size of the truncated input file isn't big enough try increasing the
    size of <your_assembly>._hackersonly[\"preview_truncate_length\"
    """



if __name__ == "__main__":
    ## test...

    ## reload autosaved data. In case you quit and came back
    JSONPATH = "/home/deren/Documents/ipyrad/tests/cli/cli.json"
    DATA = ipyrad.load_json(JSONPATH)
    DATA.run('3', force=True)

    ## reload autosaved data. In case you quit and came back
    JSONPATH = "/home/deren/Documents/ipyrad/tests/pairtest/pairtest.json"
    DATA = ipyrad.load_json(JSONPATH)
    DATA.run('3', force=True)

