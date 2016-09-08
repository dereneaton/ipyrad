#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0301

import os
import sys
import pty
import gzip
import h5py
import time
import numba
import shutil
import random
import datetime
import itertools
import subprocess
import numpy as np
import pandas as pd
import ipyrad as ip
from ipyparallel.error import TimeoutError
from ipyrad.assemble.util import *
from ipyrad.assemble.cluster_within import muscle_call, parsemuscle


import logging
LOGGER = logging.getLogger(__name__)



def muscle_align_across(args):
    """ 
    Reads in a chunk of the rebuilt clusters. Aligns clusters and writes 
    tmpaligned chunk of seqs. Also fills a tmparray with indels locations. 
    tmpchunk and tmparr will be compiled by multi_muscle_align func. 
    """
    ## parse args, names are used to order arrays by taxon names
    data, samples, chunk = args
    ## ensure sorted order
    samples.sort(key=lambda x: x.name)
    snames = [sample.name for sample in samples]

    try:
        ## data are already chunked, read in the whole thing
        infile = open(chunk, 'rb')
        clusts = infile.read().split("//\n//\n")[:-1]
        out = []

        ## tmparray to store indel information 
        maxlen = data._hackersonly["max_fragment_length"] + 30
        indels = np.zeros((len(samples), len(clusts), maxlen), dtype=np.bool)

        ## iterate over clusters and align
        loc = 0
        for loc, clust in enumerate(clusts):
            stack = []
            lines = clust.strip().split("\n")
            names = [i.split()[0][1:] for i in lines]
            seqs = [i.split()[1] for i in lines]

            LOGGER.info("seqs first %s %s", len(seqs), seqs)

            ## append counter to end of names b/c muscle doesn't retain order
            names = [j+";*"+str(i) for i, j in enumerate(names)]

            ## don't bother aligning singletons
            if len(names) <= 1:
                if names:
                    stack = [names[0]+"\n"+seqs[0]]
            else:
                ## split seqs before align if PE. If 'nnnn' not found (single end 
                ## or merged reads) then `except` will pass it to SE alignment. 
                paired = 1
                try:
                    seqs1 = [i.split("nnnn")[0] for i in seqs] 
                    seqs2 = [i.split("nnnn")[1] for i in seqs]
                except IndexError:
                    paired = 0

                if paired:
                    string1 = muscle_call(data, names, seqs1)
                    string2 = muscle_call(data, names, seqs2)
                    anames, aseqs1 = parsemuscle(data, string1)
                    anames, aseqs2 = parsemuscle(data, string2)

                    ## resort so they're in same order
                    aseqs = [i+"nnnn"+j for i, j in zip(aseqs1, aseqs2)]
                    for i in xrange(len(anames)):
                        stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])
                        sidx = [snames.index(anames[i].rsplit("_", 1)[0])]
                        ## store the indels and separator regions as indels
                        locinds = np.zeros(maxlen, dtype=np.bool)
                        for idx in range(min(maxlen, len(aseqs[i]))):
                            if aseqs[i][idx] == "-":
                                locinds[idx] = True
                        indels[sidx, loc, :] = locinds

                else:
                    string1 = muscle_call(data, names, seqs)
                    anames, aseqs = parsemuscle(data, string1)
                    for i in xrange(len(anames)):
                        stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])
                        sidx = snames.index(anames[i].rsplit("_", 1)[0])
                        ## store the indels
                        locinds = np.zeros(maxlen, dtype=np.bool)
                        for idx in range(min(maxlen, len(aseqs[i]))):
                            if aseqs[i][idx] == "-":
                                locinds[idx] = True
                        indels[sidx, loc, :] = locinds

            if stack:
                out.append("\n".join(stack))

        ## write to file after
        infile.close()
        with open(chunk, 'wb') as outfile:
            outfile.write("\n//\n//\n".join(out)+"\n")

        ## save indels array to tmp dir
        tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
        ifile = os.path.join(tmpdir, chunk+".h5")
        iofile = h5py.File(ifile, 'w')
        iofile.create_dataset('indels', data=indels)
        iofile.close()

    except Exception as inst:
        LOGGER.debug("Caught exception in muscle_align_across - {}".format(inst))
        raise inst

    return ifile, loc+1



def multi_muscle_align(data, samples, clustbits, ipyclient):
    """ 
    Sends the cluster bits to nprocessors for muscle alignment. 
    """

    ## get client
    LOGGER.info("starting alignments")
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, " aligning clusters     | {}".format(elapsed))

    ## create job queue for clustbits
    jobs = {}
    for idx, fname in enumerate(clustbits):
        jobs[idx] = lbview.apply(muscle_align_across, [data, samples, fname])
        #elapsed = datetime.timedelta(seconds=int(time.time()-start))
        #progressbar(20, 0, " aligning clusters     | {}".format(elapsed))

    LOGGER.info("clustbits %s", len(clustbits))
    LOGGER.info("submitted %s jobs to muscle_align_across", len(jobs))
    allwait = len(jobs)

    try:
        ## align clusters
        while 1:
            finished = [i.ready() for i in jobs.values()]
            if not all(finished):
                fwait = sum(finished)
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(allwait, fwait, 
                            " aligning clusters     | {}".format(elapsed))
                time.sleep(0.1)

            else:
                fwait = sum(finished)
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(allwait, fwait, 
                            " aligning clusters     | {}".format(elapsed))
                break

        ## build indels array                
        indeltups = []
        keys = jobs.keys()
        for idx in keys:
            if jobs[idx].completed and jobs[idx].successful():
                indeltups.append(jobs[idx].get())
            del jobs[idx]

        ## submit to indel entry
        if indeltups:
            build(data, samples, indeltups, clustbits, start)
        else:
            msg = "\n\n  All samples failed alignment step. Check the log files."
            raise IPyradError(msg)

    finally:
        ## Do tmp file cleanup
        for path in ["_cathaps.tmp", "_catcons.tmp", ".utemp", ".htemp"]:
            fname = os.path.join(data.dirs.consens, data.name+path)
            #if os.path.exists(fname):
            #    os.remove(fname)

        tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
        if os.path.exists(tmpdir):
            try:
                shutil.rmtree(tmpdir)
            except OSError as _:
                ## In some instances nfs creates hidden dot files in directories
                ## that claim to be "busy" when you try to remove them. Don't
                ## kill the run if you can't remove this directory.
                LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

    #if data._headers:
    print("")



def build(data, samples, indeltups, clustbits, init):
    """ 
    Builds the indels array and catclust.gz file from the aligned clusters.
    Both are tmpfiles used for filling the supercatg and superseqs arrays.
    NOT currently parallelizable
    """

    LOGGER.info(indeltups)
    ## sort into input order by chunk names
    indeltups.sort(key=lambda x: int(x[0].rsplit("_", 1)[-1][:-3]))

    ## get dims for full indel array
    maxlen = data._hackersonly["max_fragment_length"] + 30
    LOGGER.info("maxlen inside build is %s", maxlen)

    ## INIT INDEL ARRAY
    ## build an indel array for ALL loci in cat.clust.gz, 
    ## chunked so that individual samples can be pulled out
    ipath = os.path.join(data.dirs.consens, data.name+".indels")
    io5 = h5py.File(ipath, 'w')
    iset = io5.create_dataset("indels", (len(samples), data.nloci, maxlen),
                              dtype=np.bool)
                              #chunks=(len(samples), data.nloci/10, maxlen),
                              #compression="gzip")

    ## again make sure names are ordered right
    samples.sort(key=lambda x: x.name)

    ## enter all tmpindel arrays into full indel array
    ## TODO: This could be sped up with dask or numba
    for tup in indeltups:
        LOGGER.info('indeltups: %s, %s', tup[0], tup[1])
        start = int(tup[0].rsplit("_", 1)[-1][:-3])
        ioinds = h5py.File(tup[0], 'r')
        iset[:, start:start+tup[1], :] += ioinds['indels'][:] 
        ioinds.close()

        ## continued progress bar from multi_muscle_align
        elapsed = datetime.timedelta(seconds=int(time.time()-init))
        progressbar(100, 99, " aligning clusters     | {}".format(elapsed))

    io5.close()
    elapsed = datetime.timedelta(seconds=int(time.time()-init))
    progressbar(100, 100, " aligning clusters     | {}".format(elapsed))

    ## concatenate finished seq clusters into a tmp file 
    outhandle = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    with gzip.open(outhandle, 'wb') as out:
        for fname in clustbits:
            with open(fname) as infile:
                out.write(infile.read()+"//\n//\n")



def cluster(data, noreverse, ipyclient):
    """ 
    Calls vsearch for clustering across samples. 
    """

    ## input and output file handles
    cathaplos = os.path.join(data.dirs.consens, data.name+"_catshuf.tmp")
    uhaplos = os.path.join(data.dirs.consens, data.name+".utemp")
    hhaplos = os.path.join(data.dirs.consens, data.name+".htemp")
    logfile = os.path.join(data.dirs.consens, "s6_cluster_stats.txt")

    ## parameters that vary by datatype 
    ## (too low of cov values yield too many poor alignments)
    strand = "plus"
    cov = ".90"
    if data.paramsdict["datatype"] == "gbs":
        strand = "both"
        cov = ".60"
    elif data.paramsdict["datatype"] == "pairgbs":
        strand = "both"
        cov = ".90"

    ## get call string. Thread=0 means all. 
    ## old userfield: -userfields query+target+id+gaps+qstrand+qcov" \
    cmd = [ip.bins.vsearch, 
           "-cluster_smallmem", cathaplos, 
           "-strand", strand, 
           "-query_cov", cov, 
           "-id", str(data.paramsdict["clust_threshold"]), 
           "-userout", uhaplos, 
           "-notmatched", hhaplos, 
           "-userfields", "query+target+qstrand", 
           "-maxaccepts", "1", 
           "-maxrejects", "0", 
           "-minsl", "0.5", 
           "-fasta_width", "0", 
           "-threads", "0", 
           "-fulldp", 
           "-usersort", 
           "-log", logfile]

    ## override reverse clustering option
    if noreverse:
        ##strand = " -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    try:
        LOGGER.info(cmd)
        start = time.time()
        
        (dog, owner) = pty.openpty()
        proc = subprocess.Popen(cmd, stdout=owner, stderr=owner, 
                                     close_fds=True)
        done = 0
        while 1:
            dat = os.read(dog, 80192)
            if "Clustering" in dat:
                try:
                    done = int(dat.split()[-1][:-1])
                ## may raise value error when it gets to the end
                except ValueError:
                    pass
            ## break if done
            ## catches end chunk of printing if clustering went really fast
            elif "Clusters:" in dat:
                LOGGER.info("ended vsearch tracking loop")
                break
            else:
                time.sleep(0.1)
            ## print progress
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(100, done, 
                " clustering across     | {}".format(elapsed))

        ## another catcher to let vsearch cleanup after clustering is done
        proc.wait()
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 100, 
                    " clustering across     | {}".format(elapsed))

    except subprocess.CalledProcessError as inst:
        raise IPyradWarningExit("""
        Error in vsearch: \n{}\n{}""".format(inst, subprocess.STDOUT))

    finally:
        ## progress bar
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 100, " clustering across     | {}".format(elapsed))
        #if data._headers:
        print("")
        data.stats_files.s6 = logfile
        ## cleanup processes
        del proc, dog, owner



def build_h5_array(data, samples, ipyclient):
    """ build full catgs file with imputed indels. """
    ## catg array of prefiltered loci (4-dimensional) aye-aye!
    ## this can be multiprocessed!! just sum arrays at the end
    ## but one big array will overload memory, so need to be done in slices

    ## sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    #############################tempororary
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)
    ## add name and index columns to dataframe
    updf.loc[:, 3] = [i.rsplit("_", 1)[0] for i in updf[0]]
    ## create a column with only consens index
    updf.loc[:, 4] = [i.rsplit("_", 1)[1] for i in updf[0]]
    ## group by seed
    udic = updf.groupby(by=1, sort=True)
    ## get number of clusters
    data.nloci = sum(1 for _ in udic.groups.iterkeys())
    del updf, udic
    ################################


    ## get maxlen dim
    maxlen = data._hackersonly["max_fragment_length"] + 30    
    LOGGER.info("maxlen inside build_h5_array is %s", maxlen)

    ## open new h5 handle
    data.clust_database = os.path.join(
                                    data.dirs.consens, data.name+".clust.hdf5")
    io5 = h5py.File(data.clust_database, 'w')

    ## choose chunk optim size
    chunks = 100
    if data.nloci < 100:
        chunks = data.nloci
    if data.nloci > 10000:
        chunks = 500
    if data.nloci > 50000:
        chunks = 1000
    if data.nloci > 200000:
        chunks = 2000

    ## Number of elements in hdf5 chunk may not exceed 4GB
    ## This is probably not actually optimal, to have such
    ## enormous chunk sizes, could probably explore efficiency
    ## of smaller chunk sizes on very very large datasets
    chunklen = chunks * len(samples) * maxlen * 4
    if chunklen > 4000000000:
        chunks = int(round(4000000000/(len(samples) * maxlen * 4)))

    data.chunks = chunks
    LOGGER.info("data.nloci is %s", data.nloci)
    LOGGER.info("chunks is %s", data.chunks)

    # ## very big data set
    # if (data.nloci > 100000) and len(data.samples.keys()) > 100:
    #     chunks = 200
    # ## very very big data set
    # if (data.nloci > 100000) and len(data.samples.keys()) > 200:
    #     chunks = 100
    ### Warning: for some reason increasing the chunk size to 5000 
    ### caused enormous problems in step7. Not sure why. hdf5 inflate error. 
    #if data.nloci > 100000:
    #    chunks = 5000        

    ## INIT FULL CATG ARRAY
    ## store catgs with a .10 loci chunk size
    supercatg = io5.create_dataset("catgs", (data.nloci, len(samples), maxlen, 4),
                                    dtype=np.uint32,
                                    chunks=(chunks, len(samples), maxlen, 4), 
                                    compression="gzip")
    supercatg.attrs["chunksize"] = chunks
    supercatg.attrs["samples"] = [i.name for i in samples]


    ## INIT FULL SEQS ARRAY
    ## array for clusters of consens seqs
    superseqs = io5.create_dataset("seqs", (data.nloci, len(samples), maxlen),
                                    dtype="|S1",
                                    chunks=(chunks, len(samples), maxlen), 
                                    compression='gzip')
    superseqs.attrs["chunksize"] = chunks
    superseqs.attrs["samples"] = [i.name for i in samples]

    ## allele count storage
    io5.create_dataset("nalleles", (data.nloci, len(samples)), 
                                 dtype=np.uint8,
                                 chunks=(chunks, len(samples)),
                                 compression="gzip")

    ## array for filter that will be applied in step7
    io5.create_dataset("duplicates", (data.nloci, ), dtype=np.bool)
    ## array for filter that will be applied in step7
    io5.create_dataset("indels", (data.nloci, ), dtype=np.uint16)
    ## array for pair splits locations
    io5.create_dataset("splits", (data.nloci, ), dtype=np.uint16)

    ## close the big boy
    io5.close()

    ## FILL SUPERCATG and fills dupfilter, indfilter, and nalleles
    multicat(data, samples, ipyclient)

    ## FILL SUPERSEQS and fills edges(splits) for paired-end data
    fill_superseqs(data, samples)

    ## set sample states
    for sample in samples:
        sample.stats.state = 6



## TODO: if we wanted to skip depths (and singlecats) can we if there's indels?
def multicat(data, samples, ipyclient):
    """
    Runs singlecat for each sample.
    For each sample this fills its own hdf5 array with catg data & indels. 
    ## maybe this can be parallelized. Can't right now since we pass it 
    ## an open file object (indels). Room for speed improvements, tho.
    """

    ## create parallel client
    ## it might be a good idea to limit the number of engines here based on 
    ## some estimate of the individual file sizes, to avoid memory limits.

    ## read in the biggest individual file (e.g., 4GB), and limit the number
    ## of concurrent engines so that if N files that big were loaded into
    ## memory they would not exceed 75% memory. 

    ## but what if we're on multiple machines and each has its own memory?...
    ## then we need to ask each machine and use the min value

    lbview = ipyclient.load_balanced_view()

    ## make a list of jobs
    jobs = {}
    for sidx, sample in enumerate(samples):
        jobs[sample.name] = lbview.apply(singlecat, [data, sample, sidx])
        LOGGER.info("submitting %s to singlecat", sample.name)

    ## set wait job until all finished. 
    start = time.time()
    
    ## create an empty job to queue up cleaning
    #async = lbview.apply(time.sleep, 0.01)
    cleaning = {}
    last_sample = 0
    cleaning[last_sample] = lbview.apply(time.sleep, 0.1)

    ## counters for progress bar
    allwait = len(jobs)
    fwait = 0
    cwait = 0

    ## loop until finished
    while 1:
        ## if any jobs have finished then do a process
        finish_sc = [i.ready() for i in jobs.values()]
        finish_cl = [i.ready() for i in cleaning.values()]

        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        progressbar(allwait, fwait, 
                    " indexing clusters     | {}".format(elapsed))

        if any(finish_sc):
            ## clear job memory of finished jobs
            snames = jobs.keys()
            ## iterate over remaining samples/keys
            for sname in snames:
                ## if async not already finished/deleted
                if jobs[sname].completed:
                    if jobs[sname].successful():
                        LOGGER.info("cleanup here %s", cleaning)
                        ## track finished
                        fwait += 1
                        ## purge memory of the old one
                        del jobs[sname]                    

                        ## Get the async result for the last sample in the cleanup
                        ## queue. If the last sample in the cleanup queue
                        ## gets cleaned up and removed from the cleaning dict
                        ## before the current sample completes then this will
                        ## throw a KeyError. In this case just set the current
                        ## sample to clean up after the 0th job (the dummy job).
                        ## Not having this try/except creates a nasty race condition
                        ## where usually everything works fine, but sometimes
                        ## on some datasets it'll throw, and choke the whole run.
                        try:
                            last_sample_async = cleaning[last_sample]
                        except KeyError:
                            LOGGER.debug("Last sample already cleaned up. Cleaning "\
                                + "dict should be empty here - {}".format(cleaning))
                            last_sample_async = cleaning[0]

                        ## submit a clean up job. Can't start til after last one
                        with lbview.temp_flags(after=last_sample_async):
                            cleaning[sname] = lbview.apply(insert_and_cleanup, 
                                                           *[data, sname])
                        LOGGER.info("submitting %s to cleanup", sname)
                        last_sample = sname

                    else:
                        ## print error if something went wrong
                        meta = jobs[sname].metadata
                        if meta.error:
                            LOGGER.error('  sample %s did not finish', sname)
                            LOGGER.error("""\
                stdout: %s
                stderr: %s 
                error: %s""", meta.stdout, meta.stderr, meta.error)
                        del jobs[sname]

        ## but, while in singlecats loop, remove cleanups as they finish
        ## to avoid memory from climbing
        if any(finish_cl):
            snames = cleaning.keys()
            ## iterate over remaining samples/keys
            for sname in snames:
                ## Don't clean up the dummy async object because if it gets
                ## cleaned up before any real samples finish singlecat
                ## then you get a KeyError. It creates a nasty race condition.
                if cleaning[sname].completed and cleaning[sname].successful()\
                        and not sname == 0:
                    cwait += 1
                    del cleaning[sname]

        ## if finished with singlecats, move on to next progress bar. 
        if not jobs.keys():
            break

        ## print progress
        time.sleep(0.1)

    ## print final progress
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    progressbar(20, 20,
        " indexing clusters     | {}".format(elapsed))
    print("")

    ## check each sample for success and enter into full DB
    start = time.time()

    while 1:
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(allwait, max(0, cwait-1),
            " building database     | {}".format(elapsed))

        ## get finished
        finish_cl = [i.ready() for i in cleaning.values()]
        if any(finish_cl):
            snames = cleaning.keys()
            ## iterate over remaining samples/keys
            for sname in snames:
                if cleaning[sname].completed and cleaning[sname].successful():
                    cwait += 1
                    del cleaning[sname]

        ## if all finished then quit
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(allwait, max(0, cwait-1),
            " building database     | {}".format(elapsed))

        if not cleaning.keys():
            # print final progress
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(20, 20,
                " building database     | {}".format(elapsed))
            break
        else:
            time.sleep(0.1)
    print("")

    ## remove indels array
    ifile = os.path.join(data.dirs.consens, data.name+".indels")
    if os.path.exists(ifile):
        os.remove(ifile)



def insert_and_cleanup(data, sname):
    """ 
    Result is a sample name returned by singlecat. This function loads 
    three tmp arrays saved by singlecat: dupfilter, indfilter, and indels, 
    and inserts them into the superfilters array and the catg array. 
    It is NOT run in parallel, but rather, sequentially, on purpose.
    If this is changed then the way the filters are filled needs to change.
    """
    ## grab supercatg from super and get index of this sample
    io5 = h5py.File(data.clust_database, 'r+')
    catg = io5["catgs"]
    nall = io5["nalleles"]
    sidx = list(catg.attrs["samples"]).index(sname)
    LOGGER.info("insert & cleanup : %s sidx: %s", sname, sidx)

    ## get individual h5 saved in locus order and with filters
    smp5 = os.path.join(data.dirs.consens, sname+".tmp.h5")
    ind5 = h5py.File(smp5, 'r')
    icatg = ind5["icatg"]
    inall = ind5["inall"]

    ## FILL SUPERCATG -- catg is chunked by nchunk loci
    chunk = catg.attrs["chunksize"]
    for chu in xrange(0, catg.shape[0], chunk):
        catg[chu:chu+chunk, sidx, ] += icatg[chu:chu+chunk]

    ## FILL allelic information here.
    for chu in xrange(0, catg.shape[0], chunk):
        nall[chu:chu+chunk, sidx] += inall[chu:chu+chunk]

    ## put in loci that were filtered b/c of dups [0] or inds [1]
    ## this is not chunked in any way, and so needs to be changed if we 
    ## later decide to have multiple samples write to the disk at once.
    dfilters = io5["duplicates"]
    dfilters[:] += ind5["dfilter"][:]

    ## read in the existing indels counter, line it up with this samples 
    ## indels per loc, and get the max at each locus.
    ifilters = io5["indels"]
    ifilters[:] = np.array([ifilters[:], ind5["ifilter"][:]]).max(axis=0)

    ## close h5s
    io5.close()
    ind5.close()

    ## clean up / remove individual h5 catg file
    if os.path.exists(smp5):
        os.remove(smp5)



## This is where indels are imputed
def singlecat(args):
    """ 
    Orders catg data for each sample into the final locus order. This allows
    all of the individual catgs to simply be combined later. They are also in 
    the same order as the indels array, so indels are inserted from the indel
    array that is passed in. 
    """
    ## parse args
    data, sample, sidx = args

    ## load utemp cluster hits as pandas data frame
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)
    ## add name and index columns to dataframe
    updf.loc[:, 3] = [i.rsplit("_", 1)[0] for i in updf[0]]
    ## create a column with only consens index
    updf.loc[:, 4] = [i.rsplit("_", 1)[1] for i in updf[0]]
    ## group by seed
    udic = updf.groupby(by=1, sort=True)
    ## get number of clusters
    nloci = sum(1 for _ in udic.groups.iterkeys())

    ## set maximum allowable length of clusters. For single end RAD this only
    ## has to allow for a few extra indels, so 30 should be fine. For GBS, we 
    ## have to allow for the fact that two clusters could now partially overlap
    ## which previously did not within samples. Again, 30 is pretty good, but
    ## we have to set this as a max cut off so anything over this will get 
    ## trimmed. Fine, it's probably quite rare. 
    maxlen = data._hackersonly["max_fragment_length"] + 30
    LOGGER.info("maxlen inside singlecat is %s", maxlen)

    ## INIT SINGLE CATG ARRAY
    ## create an h5 array for storing catg tmp for this sample
    ## size has nloci from utemp, maxlen from hackersdict
    smpio = os.path.join(data.dirs.consens, sample.name+".tmp.h5")
    if os.path.exists(smpio):
        os.remove(smpio)
    smp5 = h5py.File(smpio, 'w')
    icatg = smp5.create_dataset('icatg', (nloci, maxlen, 4), dtype=np.uint32)
    inall = smp5.create_dataset('inall', (nloci,), dtype=np.uint8)

    ## get catg and nalleles from step5. the shape of catg: (nconsens, maxlen)
    old_h5 = h5py.File(sample.files.database, 'r')
    catarr = old_h5["catg"][:]
    nall = old_h5["nalleles"][:]
    chunksize = data.chunks

    ## local filters to fill while filling catg array
    dfilter = np.zeros(nloci, dtype=np.bool)
    ifilter = np.zeros(nloci, dtype=np.uint16)

    ## 
    LOGGER.info("filling catg")

    ## use numba compiled function to fill arrays
    ikeys = udic.groups.keys()
    dfilter = fillcats(sample.name, udic, chunksize, catarr, nall, icatg, inall, dfilter, maxlen)

    # step = iloc = 0
    # for iloc, seed in enumerate(udic.groups.iterkeys()):
    #     ipdf = udic.get_group(seed)
    #     ask = ipdf.where(ipdf[3] == sample.name).dropna()

    #     ## write to disk after 1000 writes
    #     if not iloc % chunksize:
    #         icatg[step:iloc] = locatg[:]
    #         inall[step:iloc] = lonall[:]
    #         ## reset
    #         step = iloc
    #         locatg = np.zeros((chunksize, maxlen, 4), dtype=np.uint32)
    #         lonall = np.zeros((chunksize,), dtype=np.uint8)

    #     ## fill local array one at a time
    #     if ask.shape[0] == 1: 
    #         ## if multiple hits of a sample to a locus then it is not added
    #         ## to the locus, and instead the locus is masked for exclusion
    #         ## using the filters array. Using (+=) instead of just (=) allows
    #         ## for catarr to be smaller minlen than locatg w/o error
    #         locatg[iloc-step, :catarr.shape[1]] = catarr[int(ask[4]), :icatg.shape[1], :]
    #         lonall[iloc-step] = nall[int(ask[4]),]
    #         #icatg[iloc] = catarr[int(ask[4]), :icatg.shape[1], :]
    #     elif ask.shape[0] > 1:
    #         ## store that this cluster failed b/c it had duplicate samples. 
    #         dfilter[iloc] = True

    ## store dfilter and clear memory
    smp5.create_dataset("dfilter", data=dfilter)


    LOGGER.info("filling seeds")
    ## for each locus in which Sample was the seed. 
    mask = np.array([sample.name in i for i in ikeys])
    aran = np.arange(len(mask))
    masked = aran[mask]
    for iloc in masked:
        sfill = int(ikeys[iloc].split("_")[-1])
        icatg[iloc, :catarr.shape[1]] = catarr[sfill, :icatg.shape[1]]
        inall[iloc] = nall[sfill]
    LOGGER.info("done seedmatching")

    ## close the old hdf5 connections
    old_h5.close()

    ## clear memory for hit map
    del udic
    del updf

    ## get indel locations for this sample
    ipath = os.path.join(data.dirs.consens, data.name+".indels")
    indh5 = h5py.File(ipath, 'r')
    indels = indh5["indels"][sidx, :, :]

    tmpstart = time.time()
    LOGGER.info("loading indels for %s %s, %s", 
                sample.name, sidx, np.sum(indels[:10]))

    ## insert indels into new_h5 (icatg array) which now has loci in the same
    ## order as the final clusters from the utemp file
    # for iloc in xrange(icatg.shape[0]):
    #     ## indels locations
    #     indidx = np.where(indels[iloc, :])[0]
    #     #LOGGER.info("indidx %s, len %s", indidx.shape, len(indidx))
    #     if np.any(indidx):
    #         ## store number of indels for this sample at this locus
    #         ifilter[iloc] = len(indidx)

    #         ## insert indels into catg array
    #         newrows = icatg[iloc].shape[0] + len(indidx)
    #         not_idx = np.array([k for k in range(newrows) if k not in indidx])
    #         ## create an empty new array with the right dims
    #         newfill = np.zeros((newrows, 4), dtype=icatg.dtype)
    #         ## fill with the old vals leaving the indels rows blank
    #         newfill[not_idx, :] = icatg[iloc]
    #         ## store new data into icatg
    #         icatg[iloc] = newfill[:icatg.shape[1]]

    for iloc in xrange(icatg.shape[0]):
        indidx, found = where_numba(indels, np.uint32(iloc))
        if found:
            ifilter[iloc] = indidx.shape[0]
            newfill = getfill_numba(np.uint32(indidx), icatg[iloc])
            icatg[iloc] = newfill[:icatg.shape[1]]

    elapsed = datetime.timedelta(seconds=int(time.time()- tmpstart))
    LOGGER.info("how long to do indels for %s %s", 
                sample.name, elapsed)

    ## put filteres into h5 object
    smp5.create_dataset("ifilter", data=ifilter)

    ## close the new h5 that was written to
    smp5.close()
    indh5.close()

    ## clear memory
    del indels

    ## return name for
    return sample.name



@numba.jit()
def fillcats(name, udic, chunksize, catarr, nall, icatg, inall, dfilter, maxlen):
    ## create a local array to fill until writing to disk for write efficiency
    locatg = np.zeros((chunksize, maxlen, 4), dtype=np.uint32)
    lonall = np.zeros((chunksize,), dtype=np.uint8)

    ## go go 
    step = iloc = 0
    ikeys = udic.groups.keys()
    for iloc, seed in enumerate(ikeys):
        ipdf = udic.get_group(seed)
        ask = ipdf.where(ipdf[3] == name).dropna()

        ## write to disk after 1000 writes
        if not iloc % chunksize:
            icatg[step:iloc] = locatg[:]
            inall[step:iloc] = lonall[:]
            ## reset
            step = iloc
            locatg = np.zeros((chunksize, maxlen, 4), dtype=np.uint32)
            lonall = np.zeros((chunksize,), dtype=np.uint8)

        ## fill local array one at a time
        if ask.shape[0] == 1: 
            ## if multiple hits of a sample to a locus then it is not added
            ## to the locus, and instead the locus is masked for exclusion
            ## using the filters array. Using (+=) instead of just (=) allows
            ## for catarr to be smaller minlen than locatg w/o error
            locatg[iloc-step, :catarr.shape[1]] = catarr[int(ask[4]), :icatg.shape[1], :]
            lonall[iloc-step] = nall[int(ask[4]),]
            #icatg[iloc] = catarr[int(ask[4]), :icatg.shape[1], :]
        elif ask.shape[0] > 1:
            ## store that this cluster failed b/c it had duplicate samples. 
            dfilter[iloc] = True

    ## write the leftover chunk 
    icatg[step:iloc] = locatg[:icatg[step:iloc].shape[0]]
    inall[step:iloc] = lonall[:inall[step:iloc].shape[0]]    
    return dfilter



@numba.jit(nopython=True)
def where_numba(indels, iloc):
    """ 
    numba compiled function for grabbing indel locations. Doesn't really
    speed up np.where, but is abut 6X faster for np.any. 
    """
    indidx = np.where(indels[iloc, :])[0]        
    found = np.any(indidx)
    return indidx, found



## TODO: simplify these funcs / subfuncs
def getfill_numba(indidx, tmpcatg):
    """ 
    Fast array filling function for singlecat. 
    Cannot be fully numba compiled, but has a compiled subfunction
    """
    newrows = np.uint32(tmpcatg.shape[0] + indidx.shape[0])
    newfill = np.zeros(shape=(newrows, 4), dtype=np.uint32)
    res = getfill_numba_sub(tmpcatg, indidx, newfill, newrows)
    return res



@numba.jit("u4[:,:](u4[:,:],u4[:],u4[:,:],u4)", nopython=True)
def getfill_numba_sub(tmpcatg, indidx, newfill, newrows):
    """ compiled subfunction of getfill_numba """
    allrows = np.arange(newrows)
    mask = np.ones(shape=newrows, dtype=np.uint32)
    for idx in indidx:
        mask[idx] = 0
    not_idx = allrows[mask == 1]
    newfill[not_idx, :] = tmpcatg
    return tmpcatg



def fill_superseqs(data, samples):
    """ 
    Fills the superseqs array with seq data from cat.clust 
    and fill the edges array with information about paired split locations.
    """

    ## load super to get edges
    io5 = h5py.File(data.clust_database, 'r+')
    superseqs = io5["seqs"]
    splits = io5["splits"]

    ## samples are already sorted
    snames = [i.name for i in samples]
    LOGGER.info("snames %s", snames)

    ## get maxlen again
    maxlen = data._hackersonly["max_fragment_length"] + 30
    LOGGER.info("maxlen inside fill_superseqs is %s", maxlen)

    ## data has to be entered in blocks
    infile = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    clusters = gzip.open(infile, 'r')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## iterate over clusters
    chunksize = superseqs.attrs["chunksize"]
    done = 0
    iloc = 0
    cloc = 0
    chunkseqs = np.zeros((chunksize, len(samples), maxlen), dtype="|S1")
    chunkedge = np.zeros((chunksize), dtype=np.uint16)

    while 1:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)    
        if done:
            break

        ## if chunk is full put into superseqs and reset counter
        if cloc == chunksize:
            superseqs[iloc-cloc:iloc] = chunkseqs
            splits[iloc-cloc:iloc] = chunkedge
            ## reset chunkseqs, chunkedge, cloc
            cloc = 0
            chunkseqs = np.zeros((chunksize, len(samples), maxlen), dtype="|S1")
            chunkedge = np.zeros((chunksize), dtype=np.uint16)

        ## get seq and split it
        if chunk:
            fill = np.zeros((len(samples), maxlen), dtype="|S1")
            fill.fill("N")
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = np.array([list(i) for i in piece[1::2]])
            ## fill in the separator if it exists
            separator = np.where(np.all(seqs == "n", axis=0))[0]
            if np.any(separator):
                chunkedge[cloc] = separator.min()

            ## fill in the hits
            LOGGER.info("seqs : %s", seqs)
            LOGGER.info("seqs.shape : %s", seqs.shape) 
            shlen = seqs.shape[1]
            for name, seq in zip(names, seqs):
                sidx = snames.index(name.rsplit("_", 1)[0])
                fill[sidx, :shlen] = seq[:maxlen]

            ## PUT seqs INTO local ARRAY 
            chunkseqs[cloc] = fill

        ## increase counters
        cloc += 1
        iloc += 1

    ## write final leftover chunk
    superseqs[iloc-cloc:,] = chunkseqs[:cloc]
    splits[iloc-cloc:] = chunkedge[:cloc]

    ## close super
    io5.close()

    ## edges is filled with splits for paired data.
    LOGGER.info("done filling superseqs")

    ## close handle
    clusters.close()



## TODO: This could use to be parallelized...
def build_reads_file(data, ipyclient):
    """ 
    Reconstitutes clusters from .utemp and htemp files and writes them 
    to chunked files for aligning in muscle. Return a dictionary with 
    seed:hits info from utemp file.
    """
    ## parallel client
    #lbview = ipyclient.load_balanced_view()
    start = time.time()

    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0,
        " building clusters     | {}".format(elapsed))

    LOGGER.info("building reads file -- loading utemp file into mem")
    ## read in cluster hits as pandas data frame
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)

    ## load full fasta file into a Dic... this really doesn't need to be done
    ## with a pandas DF, and right now it unnecessarily loads in all the data
    ## at once, which is not memory efficient. however, it seems fast enough
    ## and hasn't overloaded memory yet, but keep an eye on it.
    LOGGER.info("loading full _catcons file into memory")
    conshandle = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    consdf = pd.read_table(conshandle, delim_whitespace=1, 
                           header=None, compression='gzip')
    printstring = "{:<%s}    {}" % str(4+max([len(i) for i in set(consdf[::2][0])]))
    consdic = {i:j for i, j in \
                    zip(\
                        itertools.chain(*consdf[::2].values.tolist()), 
                        itertools.chain(*consdf[1::2].values.tolist())
                    )
               }

    ## make an tmpout directory and a printstring for writing to file
    tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    ## groupby index 1 (seeds) 
    groups = updf.groupby(by=1, sort=False)

    ## a chunker for writing every N
    optim = 100
    if len(groups) < 100:
        optim = len(groups)
    if len(groups) > 10000:
        optim = 1000
    if len(groups) > 200000:
        optim = 2000
    if len(groups) > 500000:
        optim = 5000
    #if len(groups) > 2000:
    #    optim = len(groups) // 10

    ## get seqs back from consdic
    clustbits = []
    locilist = []
    loci = 0

    LOGGER.info("building reads file -- loading building loci")
    tots = len(set(updf[1].values))

    ## iterate over seeds and add in hits seqs
    for seed in set(updf[1].values):
        ## get dataframe for this locus/group (gdf) 
        gdf = groups.get_group(seed)
        ## set seed name and sequence
        seedseq = consdic.get(">"+seed)
        names = [">"+seed]
        seqs = [seedseq]
        ## iterate over group dataframe (gdf) and add hits names and seqs
        ## revcomp the hit if not '+' in the df.
        for i in gdf.index:
            hit = gdf[0][i]
            names.append(">"+hit)
            if gdf[2][i] == "+":
                seqs.append(consdic[">"+hit])
            else:
                seqs.append(fullcomp(consdic[">"+hit][::-1]))

        ## append the newly created locus to the locus list
        locilist.append("\n".join([printstring.format(i, j) \
                        for i, j in zip(names, seqs)]))
        loci += 1

        ## print progress
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(tots, loci,
            " building clusters     | {}".format(elapsed))

        ## if enough loci have finished write to file to clear the mem
        if not loci % optim:
            ## a file to write results to
            handle = os.path.join(tmpdir, "tmp_"+str(loci - optim))
            with open(handle, 'w') as tmpout:
                tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
                locilist = []
                clustbits.append(handle)

    ## write the final remaining to file
    if locilist:
        handle = os.path.join(tmpdir, "tmp_"+str(loci - len(locilist)))
        with open(handle, 'w') as tmpout:
            tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
            clustbits.append(handle)

    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 100, 
        " building clusters     | {}".format(elapsed))
    print("")

    ## cleanup
    del consdic, consdf, updf

    ## return stuff
    return clustbits, loci



## This is slow currently, high memory, 
def build_input_file(data, samples, randomseed):
    """ 
    Make a concatenated consens file with sampled alleles (no RSWYMK/rswymk). 
    Orders reads by length and shuffles randomly within length classes
    """

    ## get all of the consens handles for samples that have consens reads
    ## this is better than using sample.files.consens for selecting files 
    ## b/c if they were moved we only have to edit data.dirs.consens 
    conshandles = [os.path.join(data.dirs.consens, sample.name+".consens.gz") \
                  for sample in samples if \
                  sample.stats.reads_consens]
    conshandles.sort()
    assert conshandles, "no consensus files found"
    
    ## concatenate all of the consens files
    cmd = ['cat'] + conshandles
    #cmd = ['cat'] + glob.glob(os.path.join(data.dirs.consens, '*.consens.gz'))
    allcons = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    LOGGER.debug(" ".join(cmd))
    with open(allcons, 'w') as output:
        call = subprocess.Popen(cmd, stdout=output)
        call.communicate()

    ## a string of sed substitutions for temporarily replacing hetero sites    
    subs = ["/>/!s/W/A/g", "/>/!s/w/A/g", "/>/!s/R/A/g", "/>/!s/r/A/g", 
            "/>/!s/M/A/g", "/>/!s/m/A/g", "/>/!s/K/T/g", "/>/!s/k/T/g", 
            "/>/!s/S/C/g", "/>/!s/s/C/g", "/>/!s/Y/C/g", "/>/!s/y/C/g"]
    subs = ";".join(subs)

    ## impute pseudo-haplo information to avoid mismatch at hetero sites
    ## the read data with hetero sites is put back into clustered data later
    cmd1 = ["gunzip", "-c", allcons]
    cmd2 = ["sed", subs]

    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)

    #allhaps = open(allcons.replace("_catcons.tmp", "_cathaps.tmp"), 'w')
    LOGGER.debug(" ".join(cmd1))
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    allhaps = allcons.replace("_catcons.tmp", "_cathaps.tmp")
    with open(allhaps, 'w') as output:
        LOGGER.debug(" ".join(cmd2))
        proc2 = proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=output) 
        proc2.communicate()
    proc1.stdout.close()

    ## now sort the file using vsearch
    allsort = allcons.replace("_catcons.tmp", "_catsort.tmp")  
    cmd1 = [ipyrad.bins.vsearch, "--sortbylength", allhaps, 
            "--fasta_width", "0", "--output", allsort]
    LOGGER.debug(" ".join(cmd1))
    proc1 = subprocess.Popen(cmd1)
    proc1.communicate()

    ## shuffle sequences within size classes
    random.seed(randomseed)

    allshuf = allcons.replace("_catcons.tmp", "_catshuf.tmp")      
    outdat = open(allshuf, 'w')
    indat = open(allsort, 'r')
    idat = itertools.izip(iter(indat), iter(indat))
    done = 0

    chunk = [idat.next()]
    while not done:
        ## grab 2-lines until they become shorter (unless there's only one)
        oldlen = len(chunk[-1][-1])

        while 1:
            try:
                dat = idat.next()
            except StopIteration:
                done = 1
                break
            if len(dat[-1]) == oldlen:
                chunk.append(dat)
            else:
                ## send the last chunk off to be processed
                random.shuffle(chunk)
                outdat.write("".join(itertools.chain(*chunk)))
                ## start new chunk
                chunk = [dat]
                break

    ## do the last chunk
    random.shuffle(chunk)    
    outdat.write("".join(itertools.chain(*chunk)))

    indat.close()
    outdat.close()



def run(data, samples, noreverse, force, randomseed, ipyclient):
    """ 
    Master function to run :
     1. build input file, 2. cluster, 3. split clusters into bits, 
     4. align bits, 5. build h5 array. 
    """

    ## clean the slate
    if os.path.exists(data.clust_database):
        os.remove(data.clust_database)

    ## parallel
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    
    ## make file with all samples reads    
    binput = lbview.apply(build_input_file, *[data, samples, randomseed])
    while not binput.ready():
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 0, 
            " concat/shuffle input  | {}".format(elapsed))
        time.sleep(0.5)
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 100, 
        " concat/shuffle input  | {}".format(elapsed))
    print("")

    ## call vsearch
    cluster(data, noreverse, ipyclient)

    # # # build consens clusters and returns chunk handles to be aligned
    clustbits, nloci = build_reads_file(data, ipyclient)
    data.nloci = nloci

    # ## muscle align the consens reads and creates hdf5 indel array
    # LOGGER.info("muscle alignment & building indel database")
    multi_muscle_align(data, samples, clustbits, ipyclient)

    ## builds the final HDF5 array which includes three main keys
    ## /catg -- contains all indiv catgs and has indels inserted
    ##   .attr['samples'] = [samples]
    ## /filters -- filled for dups, left empty for others until step 7.
    ##   .attr['filters'] = [f1, f2, f3, f4, f5]
    ## /seqs -- contains the clustered sequence data as string arrays
    ##   .attr['samples'] = [samples]
    ## /edges -- gets the paired split locations for now.
    ## /snps  -- left empty for now
    LOGGER.info("building full database")    
    ## calls singlecat func inside
    build_h5_array(data, samples, ipyclient)



if __name__ == "__main__":

    ## get path to test dir/ 
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__)
               )
           )
       )

    ## run test on pairgbs data1
    # TEST = ip.load.load_assembly(os.path.join(\
    #                      ROOT, "tests", "Ron", "Ron"))
    # TEST.step6(force=True)
    # print(TEST.stats)

    ## run test on pairgbs data1
    # TEST = ip.load.load_assembly(os.path.join(\
    #                      ROOT, "tests", "test_pairgbs", "test_pairgbs"))
    # TEST.step6(force=True)
    # print(TEST.stats)

    ## run test on rad data1
    #TEST = ip.load.load_assembly(os.path.join(\
    #                     ROOT, "tests", "test_rad", "data1"))
    #TEST.step6(force=True)
    #print(TEST.stats)

    ## load test data (pairgbs)
    DATA = ip.load_json(
         "/home/deren/Documents/RADmissing/rad1/half_min4.json")
    #SAMPLES = DATA.samples.values()

    # ## run step 6
    DATA.step6(force=True)


