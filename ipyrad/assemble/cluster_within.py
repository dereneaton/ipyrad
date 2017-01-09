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
import gzip
import glob
import itertools

import numpy as np
import ipyrad
import time
import datetime
import warnings
import networkx as nx

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

    ## Clean up loose files
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


def muscle_align(data, sample, chunk, maxindels):
    """ aligns reads, does split then aligning for paired reads """

    ## data are already chunked, read in the whole thing. bail if no data.
    try:
        handle = os.path.join(data.tmpdir, "{}_chunk_{}.ali"\
                              .format(sample.name, chunk))
        with open(handle, 'rb') as infile:
            clusts = infile.read().split("//\n//\n")
        if not clusts[0]:
            return 0
    except IOError:
        LOGGER.debug("skipping empty chunk - %s", handle)
        return 0

    LOGGER.debug("aligning chunk %s", handle)
    ## storage and a counter for discarded clusters due to poor alignment
    out = []
    highindels = 0

    ## iterate over clusters and align
    for clust in clusts:
        stack = []
        lines = clust.strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]
        badalign = 0

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+str(i) for i, j in enumerate(names)]

        ## don't bother aligning singletons
        if len(names) == 1:
            stack = ["{}\n{}".format(names[0], seqs[0])]
        else:
            ## split seqs if paired end seqs
            try:
                seqs1 = [i.split("nnnn")[0] for i in seqs]
                seqs2 = [i.split("nnnn")[1] for i in seqs]
                ## muscle align
                string1 = muscle_call(data, names[:200], seqs1[:200])
                string2 = muscle_call(data, names[:200], seqs2[:200])
                ## resort so they're in same order as original
                anames, aseqs1 = parsemuscle(data, string1)
                anames, aseqs2 = parsemuscle(data, string2)
                ## get leftlimit of seed, no hits can go left of this
                ## this can save pairgbs from garbage
                idxs = [i for i, j in enumerate(aseqs1[0]) if j != "-"]
                leftlimit = min(0, idxs)
                aseqs1 = [i[leftlimit:] for i in aseqs1]
                ## preserve order in ordereddict
                aseqs = zip(aseqs1, aseqs2)

                ## append to somedic if not bad align
                for i in xrange(len(anames)):
                    ## filter for max internal indels
                    intindels1 = aseqs[i][0].rstrip('-').lstrip('-').count('-')
                    intindels2 = aseqs[i][1].rstrip('-').lstrip('-').count('-')

                    if (intindels1 <= maxindels) and (intindels2 <= maxindels):
                        stack.append("{}\n{}nnnn{}".format(anames[i], aseqs[i][0], aseqs[i][1]))

                    else:
                        highindels += 1
                        badalign = 1
                        LOGGER.info("""
                              high indels: %s
                              1, 2, max: (%s, %s, %s)
                              """, aseqs[i], intindels1, intindels2, maxindels)

            except IndexError:
                ## if nnnn is not present in ALL seqs, then we end up here, and
                ## so we should remove "nnnn" spacer from any seqs that had it
                seqs = [i.replace('nnnn', '') for i in seqs]
                string1 = muscle_call(data, names[:200], seqs[:200])
                anames, aseqs = parsemuscle(data, string1)

                ## Get left and right limits, no hits can go outside of this.
                ## This can save gbs overlap data significantly.
                if 'gbs' in data.paramsdict['datatype']:

                    ## do not allow any indels left of seed's left side
                    aseqs = np.array([list(i) for i in aseqs])
                    leftlimit = max(0, np.min(np.where(aseqs[0] != "-")[0]))

                    ## right filter is the revcomped seq that goes least right
                    isrev = np.array([i.split(";")[-1][0] == "-" for i in anames])

                    if np.sum(isrev):
                        revd = aseqs[isrev, :]
                        maxrev = [np.max(np.where(i != "-")[0]) for i in revd]
                        rightlimit = np.min(maxrev)
                    else:
                        rightlimit = aseqs.shape[1]

                    ## trim all ready down to the trimmed edge length
                    aseqs = aseqs[:, leftlimit:rightlimit]
                    allindels = np.all(aseqs == "-", axis=1)
                    aseqs = aseqs[~allindels, :]
                    anames = list(np.array(anames)[~allindels])

                    ## if all seqs were trimmed then skip this loc
                    if not anames:
                        continue
                    ## make seqs from array back into lists
                    aseqs = [i.tostring() for i in aseqs]

                badalign = 0
                for i in range(len(anames)):
                    ## filter for max internal indels
                    intind = aseqs[i].rstrip('-').lstrip('-').count('-')
                    if intind <= maxindels:
                        stack.append("{}\n{}".format(anames[i], aseqs[i]))
                    else:
                        highindels += 1
                        badalign = 1
                        LOGGER.info("high indels: %s", aseqs[i])

        ## finally, add to outstack if alignment is good
        if stack:
            if not badalign:
                out.append("\n".join(stack))

    ## write to file after
    outhandle = handle.rsplit(".", 1)[0]+".aligned"
    with open(outhandle, 'wb') as outfile:
        outfile.write("\n//\n//\n".join(out)+"\n")

    ## remove the old tmp file
    os.remove(handle)
    return highindels



def parsemuscle(data, out):
    """
    parse muscle string output into two *sorted* lists.
    """
    ## remove '>' and split lines
    lines = out[1:].split("\n>")
    ## grab the names
    names = [line.split("\n", 1)[0] for line in lines]
    ## grab the seqs
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

    ## make input string
    inputstr = "\n".join(["{}\n{}".format(i, j) for i, j in zip(names, seqs)])
    cmd = [ipyrad.bins.muscle, "-quiet"]

    ## increase gap penalty if reference region is included
    ## This could use more testing/refining!
    if "_REF;+0" in names:
        cmd += ["-gapopen", "-1200"]

    ## make a call arg
    proc1 = sps.Popen(cmd, stdin=sps.PIPE, stdout=sps.PIPE, close_fds=True)
    ## return result
    return proc1.communicate(inputstr)[0]



def build_clusters(data, sample, maxindels):
    """
    Combines information from .utemp and .htemp files to create .clust files,
    which contain un-aligned clusters. Hits to seeds are only kept in the
    cluster if the number of internal indels is less than 'maxindels'.
    By default, we set maxindels=6 for this step (within-sample clustering).
    """

    ## derepfile
    derepfile = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")

    ## vsearch results files
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
                ## store the last fseq, count it, and clear fseq
                if fseqs:
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                ## occasionally write to file
                if not seqsize % 10000:
                    if seqlist:
                        clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")
                        ## reset list and counter
                        seqlist = []

                ## store the new seed on top of fseq
                fseqs.append(">{}*\n{}".format(seed, alldereps[seed]))
                lastseed = seed

            ## add match to the seed
            seq = alldereps[hit]
            ## revcomp if orientation is reversed (comp preserves nnnn)
            if ori == "-":
                seq = comp(seq)[::-1]
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
    data.tmpdir = os.path.join(pdir, data.name+'-tmpalign')
    if not os.path.exists(data.tmpdir):
        os.mkdir(data.tmpdir)

    ## If ref mapping, init samples and make the refmapping output directory.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(pdir, "{}_refmapping".format(data.name))
        if not os.path.exists(data.dirs.refmapping):
            os.mkdir(data.dirs.refmapping)




def new_apply_jobs(data, samples, ipyclient, nthreads, maxindels):
    """
    Create a DAG of prealign jobs to be run in order for each sample. Track
    Progress, report errors. Each assembly method has a slightly different
    DAG setup, calling different functions.
    """

    ## Two view objects, threaded and unthreaded
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(10, 0, " {}     | {} | s3 |"\
                .format(PRINTSTR["derep_concat_split"], elapsed))

    ## for HPC systems this should be done to make sure targets are spread
    ## among different nodes.
    if nthreads:
        if nthreads < len(ipyclient.ids):
            thview = ipyclient.load_balanced_view(targets=ipyclient.ids[::nthreads])
        elif nthreads == 1:
            thview = ipyclient.load_balanced_view()
        else:
            thview = ipyclient.load_balanced_view(targets=ipyclient.ids[::2])
    snames = [i.name for i in samples]

    ## Create DAGs for the assembly method being used, store jobs in nodes
    dag = nx.DiGraph()
    nodes = []

    ## iterate over the sample names
    for sname in snames:
        ## get list of pre-align jobs from globals based on assembly method
        method = data.paramsdict["assembly_method"]
        if method == "denovo":
            joborder = DENOVO
        elif method == "reference":
            joborder = REFERENCE
        elif method == "denovo+reference":
            joborder = DENOVO_PLUS
        else:
            joborder = DENOVO_MINUS

        ## append pre-align job for each sample to nodes list
        for func in joborder:
            nodes.append("{}-{}-{}".format(func, 0, sname))

        ## add jobs for the align funcs, each will have max 10
        for chunk in range(10):
            nodes.append("{}-{}-{}".format("muscle_align", chunk, sname))

        ## add final reconcat jobs
        nodes.append("{}-{}-{}".format("reconcat", 0, sname))

    ## add all nodes to the DAG
    for node in nodes:
        dag.add_node(node)

    ## add edges/dependencies bewteen jobs to enforce an order of operations
    ## the pattern is (first-this, then-that).
    for sname in snames:
        ## set dependencies on all samples having finished something
        for sname2 in snames:
            ## enforce clust/map func doesn't run until all derep jobs are done
            dag.add_edge("{}-{}-{}".format(joborder[0], 0, sname2),
                         "{}-{}-{}".format(joborder[1], 0, sname))
            ## enforce that job after clust/map doesn't run until all clust/map
            ## jobs are finished, this protect threaded view from being diluted
            #dag.add_edge("{}-{}-{}".format(joborder[1], 0, sname2),
            #             "{}-{}-{}".format(joborder[2], 0, sname))

        ## add more order restrictions to pre-align jobs. jobs 1,2 are above.
        ## 3 and 4 (ref and muscle chunk) are single threaded, so
        ## we'll just restrict that within sample the previous has to finish.
        for idx in xrange(2, len(joborder)):
            ## remaining steps of prealign jobs
            dag.add_edge("{}-{}-{}".format(joborder[idx-1], 0, sname),
                         "{}-{}-{}".format(joborder[idx], 0, sname))

        ## add the first align job for each sample such that it cannot start
        ## until all samples have finished chunking, and then add remaining
        ## 9 align jobs that can't start until after the first chunk is
        ## finished. This is nice because it makes the largest chunks go first,
        ## and it *greatly* simplifies the dag. However, it slows performance
        ## if there are fewer than four samples, since it will wait on the
        ## first chunk, but such small jobs are uncommon.
        for sname2 in snames:
            for chunk in range(10):
                dag.add_edge("{}-{}-{}".format("muscle_chunker", 0, sname2),
                             "{}-{}-{}".format("muscle_align", chunk, sname))   
                ## add that the final reconcat job can't start until after
                ## each chunk of its own sample has finished aligning.
                dag.add_edge("{}-{}-{}".format("muscle_align", chunk, sname),
                             "{}-{}-{}".format("reconcat", 0, sname))

        ### old style in which first muscle-align job was prioritized...
        # for sname2 in snames:
        #     dag.add_edge("{}-{}-{}".format("muscle_chunker", 0, sname2),
        #                  "{}-{}-{}".format("muscle_align", 0, sname))

        # ## add remaining 9 align jobs dependent on first (0) being finished
        # for chunk in range(1, 10):
        #     dag.add_edge("{}-{}-{}".format("muscle_align", 0, sname),
        #                  "{}-{}-{}".format("muscle_align", chunk, sname))

        #     ## add that the final reconcat job can't start until after
        #     ## each chunk of its own sample has finished aligning.
        #     dag.add_edge("{}-{}-{}".format("muscle_align", chunk, sname),
        #                  "{}-{}-{}".format("reconcat", 0, sname))

    ## dicts for storing submitted jobs and results
    results = {}

    ## submit jobs to the engines in single or threaded views. The topological
    ## sort makes sure jobs are input with all dependencies found.
    for node in nx.topological_sort(dag):
        ## get list of async results leading to this job
        deps = [results[n] for n in dag.predecessors(node)]

        ## get func, sample, and args for this func (including [data, sample])
        #print(node, deps)
        funcstr, chunk, sname = node.split("-", 2)
        func = FUNCDICT[funcstr]
        sample = data.samples[sname]

        ## args vary depending on the function
        if funcstr in ["derep_concat_split", "mapreads", "cluster"]:
            args = [data, sample, nthreads]
        elif funcstr in ["build_clusters"]:
            args = [data, sample, maxindels]
        elif funcstr in ["muscle_align"]:
            args = [data, sample, chunk, maxindels]
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
        errfunc, sfails, msgs = trackjobs(funcstr, results)
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
            msg = """
  Sample failed this step. See ipyrad_log.txt for details - {}
""".format(sample.name)
            print(msg)
            LOGGER.error("{} - {}".format(sample.name, inst))

    ## store the results to data
    data_cleanup(data)

    ## uncomment to plot the dag
    #_plot_dag(dag, results, snames)



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



def trackjobs(func, results):
    """
    Blocks and prints progress for just the func being requested from a list
    of submitted engine jobs. Returns whether any of the jobs failed.
    """
    LOGGER.info("inside trackjobs of %s", func)

    ## get just the jobs from results that are relevant to this func
    asyncs = [(i, results[i]) for i in results if i.split("-", 2)[0] == func]

    ## progress bar
    start = time.time()
    while 1:
        ## how many of this func have finished so far
        ready = [i[1].ready() for i in asyncs]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(ready), sum(ready),
                    " {}     | {} | s3 |".format(PRINTSTR[func], elapsed))
        time.sleep(0.1)
        if len(ready) == sum(ready):
            print("")
            break

    ## did any samples fail?
    success = [i[1].successful() for i in asyncs]

    ## return functionstring and error message on failure
    if not all(success):
        ## get error messages
        errmsgs = [i[1].exception() for i in asyncs if not i[1].successful()]
        ## get samlpes that failed
        sfails = [i[0].split("-", 2)[-1] for i in asyncs if not i[1].successful()]
        return func, sfails, errmsgs
    else:
        return 0, [], []



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

    strand = "plus"
    if "gbs" in data.paramsdict["datatype"]:
        strand = "both"

    ## testing bailouts (comment out)
    #if "1A_0" in infile:
    #    infile = 'xxx'

    ## pipe in a gzipped file
    if infile.endswith(".gz"):
        catcmd = ["gunzip", "-c", infile]
    else:
        catcmd = ["cat", infile]
    LOGGER.info("catcmd %s", catcmd)

    ## do dereplication with vsearch
    cmd = [ipyrad.bins.vsearch,
            "-derep_fulllength", "-",
            "-strand", strand,
            "-output", outfile,
            "-threads", str(nthreads),
            "-fasta_width", str(0),
            "-fastq_qmax", "1000",
            "-sizeout"]
    LOGGER.info("derep cmd %s", cmd)

    ## run vsearch
    proc1 = sps.Popen(catcmd, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    proc2 = sps.Popen(cmd, stdin=proc1.stdout, stderr=sps.STDOUT, stdout=sps.PIPE, close_fds=True)
    errmsg = proc2.communicate()[0]
    if proc2.returncode:
        LOGGER.error("error inside derep_and_sort %s", errmsg)
        raise IPyradWarningExit(errmsg)



def data_cleanup(data):
    """ cleanup / statswriting function for Assembly obj """
    data.stats_dfs.s3 = data._build_stat("s3")
    data.stats_files.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    with open(data.stats_files.s3, 'w') as outfile:
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



def cluster(data, sample, nthreads):
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
    cov = 0.90
    minsl = 0.5
    if data.paramsdict["datatype"] == "gbs":
        strand = "both"
        cov = 0.33
        minsl = 0.33
    elif data.paramsdict["datatype"] == 'pairgbs':
        strand = "both"
        cov = 0.75
        minsl = 0.75

    ## If this value is not null (which is the default) then override query cov
    if data._hackersonly["query_cov"]:
        cov = " -query_cov "+str(data._hackersonly["query_cov"])
        assert data._hackersonly["query_cov"] <= 1, "query_cov must be <= 1.0"

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

    ## check for errors
    if proc.returncode:
        LOGGER.error("error %s: %s", cmd, res)
        raise IPyradWarningExit("cmd %s: %s", res)



def muscle_chunker(data, sample):
    """
    Splits the muscle alignment into chunks. Each chunk is run on a separate
    computing core. Because the largest clusters are at the beginning of the 
    clusters file, assigning equal clusters to each file would put all of the 
    large cluster, that take longer to align, near the top. So instead we 
    randomly distribute the clusters among the files. 
    """
    ## log our location for debugging
    LOGGER.info("inside muscle_chunker")

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



def derep_concat_split(data, sample, nthreads):
    """
    Running on remote Engine. Refmaps, then merges, then dereplicates,
    then denovo clusters reads.
    """

    ## report location for debugging
    LOGGER.info("INSIDE derep %s", sample.name)

    ## concatenate edits files within Samples if an Assembly was formed from
    ## merging several assemblies. This returns a new sample.files.edits with
    ## the concat file. No change if not merged Assembly.
    #sample = concat_edits(data, sample)
    sample.files.edits = concat_multiple_edits(data, sample)
    LOGGER.info("Passed concat edits")

    ## Denovo: merge or concat fastq pairs [sample.files.pairs]
    ## Reference: only concat fastq pairs  []
    ## Denovo + Reference:
    if 'pair' in data.paramsdict['datatype']:
        ## merge pairs that overlap and concatenate non-overlapping pairs with
        ## a "nnnn" separator. merge_pairs takes the unmerged files list as an
        ## argument because we're reusing this code in the refmap pipeline.
        LOGGER.debug("Merging pairs - %s", sample.files.edits)
        ## If doing any reference mapping do not merge only concatenate so the
        ## reads can be split later and mapped separately.
        merge = rcomp = 1
        if "reference" in data.paramsdict["assembly_method"]:
            merge = rcomp = 0

        ## merge R1 and R2 before we derep
        mergefile = os.path.join(data.dirs.edits, sample.name+"_merged_.fastq")
        nmerged = merge_pairs(data, sample.files.edits, mergefile, rcomp, merge)
        sample.files.edits = [(mergefile, )]
        sample.stats.reads_merged = nmerged
        LOGGER.info("Merged pairs %s %s", sample.name, sample.stats.reads_merged)
        LOGGER.debug("Merged file - {}".format(sample.files.edits))

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
            ## make dirs that are needed
            setup_dirs(data)
            ## if refmapping make filehandles that will be persistent
            if not data.paramsdict["assembly_method"] == "denovo":
                for sample in subsamples:
                    refmap_init(data, sample)

            ## use saved value if present, else use hard-coded thread count
            nthreads = 2
            if "threads" in data._ipcluster.keys():
                nthreads = int(data._ipcluster["threads"])

                ## TODO: choose threading based on ncores and nsamples...
                _ncpus = len(ipyclient)
                if _ncpus > len(data.samples):
                    pass

            maxindels = 8
            args = [data, subsamples, ipyclient, nthreads, maxindels]
            new_apply_jobs(*args)

        finally:
            ## this can fail if jobs were not stopped properly and are still
            ## writing to tmpdir.
            try:
                if os.path.exists(data.tmpdir):
                    shutil.rmtree(data.tmpdir)
                ## get all refmap_derep.fastqs
                rdereps = glob.glob(os.path.join(data.dirs.edits,
                                     "*-refmap_derep.fastq"))
                ## Remove the unmapped fastq files
                for rmfile in rdereps:
                    os.remove(rmfile)

            except Exception as _:
                LOGGER.warning("failed to cleanup files/dirs")







### GLOBALS

THREADED_FUNCS = ["derep_concat_split", "cluster", "mapreads"]

PRINTSTR = {
    "derep_concat_split" : "dereplicating    ",
    "mapreads" :           "mapping          ",
    "cluster" :            "clustering       ",
    "build_clusters" :     "building clusters",
    "ref_muscle_chunker" : "finalize mapping ",
    "muscle_chunker" :     "chunking         ",
    "muscle_align" :       "aligning         ",
    "reconcat" :           "concatenating    "
    }


FUNCDICT = {
    "derep_concat_split" : derep_concat_split,
    "mapreads" :           mapreads,
    "cluster" :            cluster,
    "build_clusters" :     build_clusters,
    "ref_muscle_chunker" : ref_muscle_chunker,
    "muscle_chunker" :     muscle_chunker,
    "muscle_align" :       muscle_align,
    "reconcat" :           reconcat
    }


## Pre-align funcs for the four assembly methods
DENOVO = ["derep_concat_split",
          "cluster",
          "build_clusters",
          "muscle_chunker"]

REFERENCE = ["derep_concat_split",
             "mapreads",
             "ref_muscle_chunker",
             "muscle_chunker"]

DENOVO_PLUS = ["derep_concat_split",
               "mapreads",
               "cluster",
               "build_clusters",
               "ref_muscle_chunker",
               "muscle_chunker"]

DENOVO_MINUS = ["derep_concat_split",
                "mapreads",
                "cluster",
                "build_clusters",
                "muscle_chunker"]

ALIGNFUNCS = ["muscle_align",
              "reconcat"]


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

    ## run step 6
    DATA.run('3', force=True)
