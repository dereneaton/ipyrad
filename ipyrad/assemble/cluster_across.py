#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0212
# pylint: disable=C0301

import os
import sys
import gzip
import h5py
import time
import random
import shutil
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
    snames = [sample.name for sample in samples]

    ## data are already chunked, read in the whole thing
    infile = open(chunk, 'rb')
    clusts = infile.read().split("//\n//\n")[:-1]
    out = []

    ## tmparray to store indel information 
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2
    indels = np.zeros((len(samples), len(clusts), maxlen), dtype=np.bool)

    ## iterate over clusters and align
    loc = 0
    for loc, clust in enumerate(clusts):
        stack = []
        lines = clust.strip().split("\n")
        names = [i.split()[0][1:] for i in lines]
        seqs = [i.split()[1] for i in lines]

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+";*"+str(i) for i, j in enumerate(names)]

        ## don't bother aligning singletons
        if len(names) <= 1:
            if names:
                stack = [names[0]+"\n"+seqs[0]]
        else:
            ## split seqs before align if PE. If 'nnnn' not found (single end 
            ## or merged reads) then `except` will pass it to SE alignment. 
            try:
                seqs1 = [i.split("nnnn")[0] for i in seqs] 
                seqs2 = [i.split("nnnn")[1] for i in seqs]

                string1 = muscle_call(data, names, seqs1)
                string2 = muscle_call(data, names, seqs2)
                anames, aseqs1 = parsemuscle(string1)
                anames, aseqs2 = parsemuscle(string2)

                ## resort so they're in same order
                aseqs = [i+"nnnn"+j for i, j in zip(aseqs1, aseqs2)]
                for i in xrange(len(anames)):
                    stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])
                    ## store the indels and separator regions as indels
                    locinds = np.array(list(aseqs[i])) == "-"
                    sidx = [snames.index(anames[i].rsplit("_", 1)[0])]
                    indels[sidx, loc, :locinds.shape[0]] = locinds

            except IndexError:
                string1 = muscle_call(data, names, seqs)
                anames, aseqs = parsemuscle(string1)
                for i in xrange(len(anames)):
                    stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])                    
                    ## store the indels
                    locinds = np.array(list(aseqs[i])) == "-"
                    sidx = snames.index(anames[i].rsplit("_", 1)[0])
                    indels[sidx, loc, :locinds.shape[0]] = locinds

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    infile.close()
    with open(chunk, 'wb') as outfile:
        outfile.write("\n//\n//\n".join(out)+"\n")

    return chunk, indels, loc+1



def multi_muscle_align(data, samples, clustbits, ipyclient):
    """ 
    Sends the cluster bits to nprocessors for muscle alignment. 
    """

    ## get client
    LOGGER.info("starting alignments")
    lbview = ipyclient.load_balanced_view()

    ## create job queue for clustbits
    submitted_args = []
    for fname in clustbits:
        submitted_args.append([data, samples, fname])

    ## submit queued jobs
    jobs = {}
    for idx, job in enumerate(submitted_args):
        jobs[idx] = lbview.apply(muscle_align_across, job)

    ## set wait job until all finished. 
    tmpids = list(itertools.chain(*[i.msg_ids for i in jobs.values()]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1)

    ## wait on results inside try statement. If interrupted we ensure clean up
    ## and save of finished results
    try:
        allwait = len(jobs)
        while 1:
            if not res.ready():
                fwait = sum([jobs[i].ready() for i in jobs])
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(allwait, fwait, 
                            " aligning clusters 2/3  | {}".format(elapsed))
                ## got to next print row when done
                sys.stdout.flush()
                time.sleep(1)
            else:
                ## print final statement
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(20, 20, " aligning clusters 2/3  | {}".format(elapsed))
                break

    except (KeyboardInterrupt, SystemExit):
        print('\n  Interrupted! Cleaning up... ')
        raise
    
    finally:
        ## save result tuples into a list
        indeltups = []

        ## clean up jobs
        for idx in jobs:
            ## grab metadata
            meta = jobs[idx].metadata            

            ## do this if success
            if meta.completed:
                ## get the results
                indeltups.append(jobs[idx].get())

            ## if not done do nothing, if failure print error
            else:
                LOGGER.error('  chunk %s did not finish', idx)
                if meta.error:
                    LOGGER.error("""\
                stdout: %s
                stderr: %s 
                error: %s""", meta.stdout, meta.stderr, meta.error)

        ## build indels array and seqs if res finished succesful
        if res.metadata.completed:
            build(data, samples, indeltups, clustbits)

        ## Do tmp file cleanup
        for path in ["_cathaps.tmp", "_catcons.tmp", ".utemp", ".htemp"]:
            fname = os.path.join(data.dirs.consens, data.name+path)
            if os.path.exists(path):
                os.remove(path)

        tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
        if os.path.exists(tmpdir):
            try:
                shutil.rmtree(tmpdir)
            except OSError as _:
                ## In some instances nfs creates hidden dot files in directories
                ## that claim to be "busy" when you try to remove them. Don't
                ## kill the run if you can't remove this directory.
                LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))

    if data._headers:
        print("")



def build(data, samples, indeltups, clustbits):
    """ 
    Builds the indels array and catclust.gz file from the aligned clusters.
    """
    ## sort into input order by chunk names
    indeltups.sort(key=lambda x: int(x[0].rsplit("_", 1)[1]))

    ## get dims for full indel array
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2

    ## INIT INDEL ARRAY
    ## build an indel array for ALL loci in cat.clust.gz
    ipath = os.path.join(data.dirs.consens, data.name+".indels")
    io5 = h5py.File(ipath, 'w')
    iset = io5.create_dataset("indels", (len(samples), data.nloci, maxlen),
                              dtype=np.bool, 
                              chunks=(1, data.nloci, maxlen))
                              #, compression="gzip")

    ## enter all tmpindel arrays into full indel array
    for tup in indeltups:
        start = int(tup[0].rsplit("_", 1)[1])
        iset[:, start:start+tup[2], :] = tup[1]
    io5.close()

    ## concatenate finished seq clusters into a tmp file 
    outhandle = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    with gzip.open(outhandle, 'wb') as out:
        for fname in clustbits:
            with open(fname) as infile:
                out.write(infile.read()+"//\n//\n")




def cluster(data, noreverse):
    """ 
    Calls vsearch for clustering across samples. 
    """

    ## input and output file handles
    cathaplos = os.path.join(data.dirs.consens, data.name+"_cathaps.tmp")
    uhaplos = os.path.join(data.dirs.consens, data.name+".utemp")
    hhaplos = os.path.join(data.dirs.consens, data.name+".htemp")

    ## parameters that vary by datatype 
    ## (too low of cov values yield too many poor alignments)
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .60 "
    elif data.paramsdict["datatype"] == "pairgbs":
        reverse = " -strand both "
        cov = " -query_cov .90 "
    else:
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    ## get call string. Thread=0 means all. 
    ## old userfield: -userfields query+target+id+gaps+qstrand+qcov" \
    cmd = ip.bins.vsearch+\
        " -cluster_smallmem "+cathaplos \
       +reverse \
       +cov \
       +" -id "+str(data.paramsdict["clust_threshold"]) \
       +" -userout "+uhaplos \
       +" -notmatched "+hhaplos \
       +" -userfields query+target+qstrand" \
       +" -maxaccepts 1" \
       +" -maxrejects 0" \
       +" -minsl 0.5" \
       +" -fasta_width 0" \
       +" -threads 0" \
       +" -fulldp " \
       +" -usersort " \
       +" -quiet"

    try:
        LOGGER.debug(cmd)
        #proc = subprocess.Popen(cmd, stderr=subprocess.STDOUT,            
        #                             stdout=subprocess.PIPE)
        ## fork to non-blocking stream so we can read the output
        start = time.time()

        # nbsr = NonBlockingStreamReader(proc.stdout)
        # #proc.stdin.write(cmd)
        # ## parse stdout as it comes in
        # while 1:
        #     output = nbsr.readline(0.1)
        #     print(output)
        #     if not output:
        #         break
        #     print(output)
        #     elapsed = datetime.timedelta(seconds=int(time.time()-start))
        #     progressbar(100, 5, " clustering across 1/3  | {}".format(elapsed))
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 1, " clustering across 1/3  | {}".format(elapsed))
        subprocess.check_call(cmd, shell=True, 
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
        ## change to subprocess.call and parse the PIPE while it streams in
        ## to create a progress bar. 
        ## ...

    except subprocess.CalledProcessError as inst:
        raise IPyradWarningExit("""
        Error in vsearch: \n{}\n{}""".format(inst, subprocess.STDOUT))

    ## progress bar
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 100, " clustering across 1/3  | {}".format(elapsed))
    print("")




def build_h5_array(data, samples, ipyclient):
    """ build full catgs file """
    ## catg array of prefiltered loci (4-dimensional) aye-aye!
    ## this can be multiprocessed!! just sum arrays at the end
    ## but one big array will overload memory, so need to be done in slices

    ## sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    ## initialize an hdf5 array of the super catg matrix with dims
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2

    ## open h5 handle
    data.database = os.path.join(data.dirs.consens, data.name+".hdf5")
    io5 = h5py.File(data.database, 'w')

    ## choose chunk optim size
    chunks = data.nloci
    if data.nloci > 10000:
        chunks = 1000

    ## INIT FULL CATG ARRAY
    ## store catgs with a .10 loci chunk size
    supercatg = io5.create_dataset("catgs", (data.nloci, len(samples), maxlen, 4),
                                    dtype=np.uint32,
                                    chunks=(chunks, len(samples), maxlen, 4))
                                    #compression="gzip")
    supercatg.attrs["chunksize"] = chunks
    supercatg.attrs["samples"] = [i.name for i in samples]


    ## INIT FULL SEQS ARRAY
    ## array for clusters of consens seqs
    superseqs = io5.create_dataset("seqs", (data.nloci, len(samples), maxlen),
                                    dtype="|S1",
                                    chunks=(chunks, len(samples), maxlen))
    superseqs.attrs["chunksize"] = chunks
    superseqs.attrs["samples"] = [i.name for i in samples]


    ## INIT FULL SNPS ARRAY
    ## array for snp string, 2 cols, - and *
    snps = io5.create_dataset("snps", (data.nloci, maxlen, 2), dtype=np.bool,
                              chunks=(chunks, maxlen, 2))
    snps.attrs["chunksize"] = chunks
    snps.attrs["names"] = ["-", "*"]


    ## INIT FULL FILTERS ARRAY
    ## array for filters that will be applied in step7
    filters = io5.create_dataset("filters", (data.nloci, 6), dtype=np.bool)#,
    filters.attrs["filters"] = ["duplicates", "max_indels", "max_snps", 
                                "max_hets", "min_samps", "bad_edges"]


    ## INIT FULL EDGE ARRAY
    ## array for edgetrimming 
    edges = io5.create_dataset("edges", (data.nloci, 5), dtype=np.uint16,
                               chunks=(chunks, 5))
    edges.attrs["chunksize"] = chunks
    edges.attrs["names"] = ["R1_L", "R1_R", "R2_L", "R2_R", "sep"]

    ## close the big boy
    io5.close()

    ## FILL SUPERCATG and fills dupfilter and indfilter
    multicat(data, samples, ipyclient)

    ## FILL SUPERSEQS and fills edges(splits) for paired-end data
    fill_superseqs(data, samples)

    ## set sample states
    for sample in samples:
        sample.stats.state = 6

    ## write database as stats output?
    data.stats_files.s6 = "[this step has no stats ouput]"



def multicat(data, samples, ipyclient):
    """
    Runs singlecat for each sample.
    For each sample this fills its own hdf5 array with catg data & indels. 
    ## maybe this can be parallelized. Can't right now since we pass it 
    ## an open file object (indels). Room for speed improvements, tho.
    """

    ## make a list of jobs
    submitted_args = []
    for sidx, sample in enumerate(samples):
        submitted_args.append([data, sample, sidx])

    ## create parallel client
    lbview = ipyclient.load_balanced_view()

    ## submit initial nproc jobs
    LOGGER.info("submitting singlecat jobs")
    jobs = {}
    for args in submitted_args:
        jobs[args[1].name] = lbview.apply(singlecat, args)

    ## set wait job until all finished. 
    tmpids = list(itertools.chain(*[i.msg_ids for i in jobs.values()]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1)

    try:
        allwait = len(jobs)
        while 1:
            if not res.ready():
                fwait = sum([i.ready() for i in jobs.values()])
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(allwait, fwait, 
                            " building database 3/3  | {}".format(elapsed))
                ## got to next print row when done
                sys.stdout.flush()
                time.sleep(1)
            else:
                ## print final statement
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(20, 20, " building database 3/3  | {}".format(elapsed))
                break

    except (KeyboardInterrupt, SystemExit):
        print('\n  Interrupted! Cleaning up... ')
    
    finally:
        for sname in jobs:
            ## grab metadata
            meta = jobs[sname].metadata

            ## do this if success
            if meta.completed:
                ## get the results
                insert_and_cleanup(data, sname)

            ## if not done do nothing, if failure print error
            else:
                LOGGER.error('  sample %s did not finish', sname)
                if meta.error:
                    LOGGER.error("""\
                stdout: %s
                stderr: %s 
                error: %s""", meta.stdout, meta.stderr, meta.error)

    if data._headers:
        print("")

    ## remove indels array
    os.remove(os.path.join(data.dirs.consens, data.name+".indels"))



def insert_and_cleanup(data, sname):
    """ 
    Result is a sample name returned by singlecat. This function loads 
    three tmp arrays saved by singlecat: dupfilter, indfilter, and indels, 
    and inserts them into the superfilters array and the catg array. 
    """
    ## grab supercatg from super and get index of this sample
    io5 = h5py.File(data.database, 'r+')
    catg = io5["catgs"]
    sidx = list(catg.attrs["samples"]).index(sname)

    ## get individual h5 saved in locus order and with filters
    smp5 = os.path.join(data.dirs.consens, sname+".tmp.h5")
    ind5 = h5py.File(smp5, 'r')
    icatg = ind5["icatg"]

    ## FILL SUPERCATG -- catg is chunked by nchunk loci
    chunk = catg.attrs["chunksize"]
    for chu in xrange(0, catg.shape[0], chunk):
        catg[chu:chu+chunk, sidx, ...] = icatg[chu:chu+chunk]

    ## put in loci that were filtered b/c of dups [0] or inds [1]
    filters = io5["filters"]
    filters[:, 0] = ind5["dupfilter"][:]
    filters[:, 1] = ind5["indfilter"][:]

    ## close h5s
    io5.close()
    ind5.close()

    ## clean up / remove individual h5 catg file
    if os.path.exists(smp5):
        os.remove(smp5)

    # ## FILL SUPERCATG -- combine individual hdf5 arrays into supercatg
    # h5handles = []
    # for sidx, sample in enumerate(samples):
    #     h5handle = os.path.join(data.dirs.consens, sample.name+".hdf5")
    #     h5handles.append(h5handle)
    #     ## open, read, and close individual data base
    #     with h5py.File(h5handle, 'r') as singleh5:
    #         icatg = singleh5[sample.name]
    #         supercatg[:, sidx, :, :] = icatg




def singlecat(args):
    """ 
    Orders catg data for each sample into the final locus order. This allows
    all of the individual catgs to simply be combined later. They are also in 
    the same order as the indels array, so indels are inserted from the indel
    array that is passed in. 
    """
    #LOGGER.debug("singlecat: %s", sample.name)
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

    ## get indel locations for this sample
    ipath = os.path.join(data.dirs.consens, data.name+".indels")
    with h5py.File(ipath, 'r') as indh5:
        indels = indh5["indels"][sidx, :, :]

    ## get maxlen
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2

    ## INIT SINGLE CATG ARRAY
    ## create an h5 array for storing catg tmp for this sample
    ## size has nloci from utemp, maxlen from hackersdict
    smpio = os.path.join(data.dirs.consens, sample.name+".tmp.h5")
    ## if h5 for this samples exists delete it.
    if os.path.exists(smpio):
        os.remove(smpio)
    smp5 = h5py.File(smpio, 'w')
    icatg = smp5.create_dataset('icatg', (nloci, maxlen, 4), dtype=np.uint32)

    ## LOAD IN STEP5 CATG ARRAY
    ## get catg from step5 for this sample, the shape is (nconsens, maxlen)
    old_h5 = h5py.File(sample.files.database, 'r')
    catarr = old_h5["catg"][:]

    ## filters to fill while filling catg array
    dupfilter = np.zeros(nloci, dtype=np.bool)
    indfilter = np.zeros(nloci, dtype=np.bool)

    ## create a local array to fill until writing to disk for mem efficiency
    local = np.zeros((10000, maxlen, 4), dtype=np.uint32)
    start = iloc = 0
    for iloc, seed in enumerate(udic.groups.iterkeys()):
        ipdf = udic.get_group(seed)
        ask = ipdf.where(ipdf[3] == sample.name).dropna()

        ## write to disk after 10000 writes
        if not iloc % 10000:
            icatg[start:iloc] = local[:]
            ## reset
            start = iloc
            local = np.zeros((10000, maxlen, 4), dtype=np.uint32)

        ## fill local array one at a time
        if ask.shape[0] == 1: 
            ## if multiple hits of a sample to a locus then it is not added
            ## to the locus, and instead the locus is masked for exclusion
            ## using the filters array.
            local[iloc-start] = catarr[int(ask[4]), :icatg.shape[1], :]
            #icatg[iloc] = catarr[int(ask[4]), :icatg.shape[1], :]
        elif ask.shape[0] > 1:
            ## store that this cluster failed b/c it had duplicate samples. 
            dupfilter[iloc] = True

    ## write the leftover 
    icatg[start:iloc] = local[:icatg[start:iloc].shape[0]]

    ## for each locus in which Sample was the seed. This writes quite a bit
    ## slower than the local setting
    seedmatch1 = (sample.name in i for i in udic.groups.keys())
    seedmatch2 = (i for i, j in enumerate(seedmatch1) if j)
    LOGGER.debug("starting seedmatching")
    for iloc in seedmatch2:
        try:
            sfill = udic.groups.keys()[iloc].split("_")[-1]
            icatg[iloc] = catarr[sfill, :, :]
        except IndexError:
            LOGGER.error("""
                sample %s, 
                no seedmatch at iloc %s, 
                """)
    ## close the old hdf5 connections
    old_h5.close()

    ## insert indels into new_h5 (icatg array) which now has loci in the same
    ## order as the final clusters from the utemp file
    for iloc in xrange(icatg.shape[0]):
        ## indels locations
        indidx = np.where(indels[iloc, :])[0]
        if np.any(indidx):
            ### apply indel filter 
            if "pair" in data.paramsdict["datatype"]:
                ## get lenght of locus 
                if len(indidx) > sum(data.paramsdict["max_Indels_locus"]):
                    indfilter[iloc] = True
            else:
                if len(indidx) > data.paramsdict["max_Indels_locus"][0]:
                    indfilter[iloc] = True

            ## insert indels into catg array
            newrows = icatg[iloc].shape[0] + len(indidx)
            not_idx = np.array([k for k in range(newrows) if k not in indidx])
            ## create an empty new array with the right dims
            newfill = np.zeros((newrows, 4), dtype=icatg.dtype)
            ## fill with the old vals leaving the indels rows blank
            newfill[not_idx, :] = icatg[iloc]
            ## store new data into icatg
            icatg[iloc] = newfill[:icatg[iloc].shape[0]]

    ## put filteres into h5 object
    smp5.create_dataset("indfilter", data=indfilter)
    smp5.create_dataset("dupfilter", data=dupfilter)
    ## close the new h5 that was written to
    smp5.close()

    ## return name for
    return sample.name



def fill_superseqs(data, samples):
    """ fills the superseqs array with seq data from cat.clust """

    ## load super
    io5 = h5py.File(data.database, 'r+')
    superseqs = io5["seqs"]
    edges = io5["edges"]

    ## samples are already sorted
    snames = [i.name for i in samples]
    ## get maxlen again
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2

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
    #LOGGER.info("chunkseqs is %s", chunkseqs.shape)
    while 1:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)    
        if done:
            break
        ## input array must be this shape
        if chunk:
            fill = np.zeros((len(samples), maxlen), dtype="|S1")
            fill.fill("N")
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = np.array([list(i) for i in piece[1::2]])
            ## fill in the separator if it exists
            separator = np.where(np.all(seqs == "n", axis=0))[0]
            if np.any(separator):
                chunkedge[iloc] = separator.min()

            ## fill in the hits
            #indices = range(len(snames))
            shlen = seqs.shape[1]
            for name, seq in zip(names, seqs):
                sidx = snames.index(name.rsplit("_", 1)[0])
                fill[sidx, :shlen] = seq

        ## if chunk is full put into superseqs
        if cloc == chunksize:
            superseqs[iloc-cloc:iloc] = chunkseqs
            edges[iloc-cloc:iloc, 4] = chunkedge
            ## reset chunkseqs
            cloc = 0
            chunkseqs = np.zeros((chunksize, len(samples), maxlen), dtype="|S1")
            chunkedge = np.zeros((chunksize), dtype=np.uint16)

        ## PUT seqs INTO local ARRAY 
        chunkseqs[cloc] = fill

        ## increase counters
        cloc += 1
        iloc += 1

    ## write final leftover chunk
    superseqs[iloc-cloc:,] = chunkseqs[:cloc]
    edges[iloc-cloc:, 4] = chunkedge[:cloc]

    ## close super
    io5.close()

    ## edges is filled with splits for paired data.
    LOGGER.info("done filling superseqs")
    #LOGGER.info('edges preview %s', edges[:10])
    #LOGGER.info('edges preview %s', edges[-10:])

    ## close handle
    clusters.close()



def build_reads_file(data):
    """ reconstitutes clusters from .utemp and htemp files and writes them 
    to chunked files for aligning in muscle. Return a dictionary with 
    seed:hits info from utemp file.
    """
    LOGGER.info("building reads file -- loading utemp file into mem")
    ## read in cluster hits as pandas data frame
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)

    ## load full fasta file into a Dic
    LOGGER.info("loading full _catcons file into memory")
    conshandle = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    consdf = pd.read_table(conshandle, delim_whitespace=1, header=None)
    printstring = "{:<%s}    {}" % max([len(i) for i in set(consdf[0])])
    consdic = consdf.set_index(0)[1].to_dict()
    del consdf

    ## make an tmpout directory and a printstring for writing to file
    tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    ## groupby index 1 (seeds) 
    groups = updf.groupby(by=1, sort=False)

    ## a chunker for writing every N
    optim = 100
    if len(groups) > 2000:
        optim = len(groups) // 10

    ## get seqs back from consdic
    clustbits = []
    locilist = []
    loci = 0

    LOGGER.info("building reads file -- loading building loci")
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

    ## store nloci 
    data.nloci = loci

    ## return stuff
    return clustbits



def build_input_file(data, samples, outgroups, randomseed):
    """ 
    Make a concatenated consens files to input to vsearch. 
    Orders reads by length and shuffles randomly within length classes
    """

    ##  make a copy list that will not have outgroups excluded
    conshandles = list(itertools.chain(*[samp.files.consens \
                                         for samp in samples]))
    ## remove outgroup sequences, add back in later to bottom after shuffling
    ## outgroups could be put to end of sorted list
    conshandles.sort()
    assert conshandles, "no consensus files found"
                
    ## output file for consens seqs from all taxa in consfiles list
    allcons = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    allhaps = open(allcons.replace("_catcons.tmp", "_cathaps.tmp"), 'w')

    LOGGER.info("concatenating sequences into _catcons.tmp file")
    ## combine cons files and write as pandas readable format to all file
    ## this is the file that will be read in later to build clusters
    printstr = "{:<%s}    {}" % 100   ## max len allowable name
    with open(allcons, 'wb') as consout:
        for qhandle in conshandles:
            with gzip.open(qhandle, 'r') as tmpin:
                data = tmpin.readlines()
                names = iter(data[::2])
                seqs = iter(data[1::2])
                consout.write("".join(printstr.format(i.strip(), j) \
                              for i, j in zip(names, seqs)))

    start = time.time()
    LOGGER.info("sorting sequences into len classes")
    ## created version with haplos that is also shuffled within seqlen classes
    random.seed(randomseed)
    ## read back in cons as a data frame
    consdat = pd.read_table(allcons, delim_whitespace=1, header=None)
    ## make a new column with seqlens and then groupby seqlens
    consdat[2] = pd.Series([len(i) for i in consdat[1]], index=consdat.index)
    lengroups = consdat.groupby(by=2)
    lenclasses = sorted(set(consdat[2]), reverse=True)

    LOGGER.info("shuffling sequences within len classes & sampling alleles")
    ## write all cons in pandas readable format
    for lenc in lenclasses:
        group = lengroups.get_group(lenc)
        ## shuffle the subgroup
        shuf = group.reindex(np.random.permutation(group.index))
        ## convert ambiguity codes into a sampled haplotype for any sample
        ## to use for clustering, but ambiguities are still saved in allcons
        writinghaplos = []
        for idx, ind in enumerate(shuf.index):
            ## append to list
            writinghaplos.append("\n".join([shuf[0][ind], 
                                            splitalleles(shuf[1][ind])[0]]))
            ## write and reset
            if not idx % 1000:
                allhaps.write("\n".join(writinghaplos)+"\n")
                writinghaplos = []
        ## write final chunk
        allhaps.write("\n".join(writinghaplos)+"\n")
    allhaps.close()
    LOGGER.info("sort/shuf/samp took %s seconds", int(time.time()-start))



def run(data, samples, noreverse, force, randomseed, ipyclient):
    """ subselect and pass args for across-sample clustering """

    ## clean the slate
    if os.path.exists(data.database):
        os.remove(data.database)

    ## make file with all samples reads, allow priority to shunt outgroups
    ## to the end of the file
    outgroups = []
    LOGGER.info("creating input files")
    build_input_file(data, samples, outgroups, randomseed)

    ## call vsearch
    LOGGER.info("clustering")    
    cluster(data, noreverse)

    ## build consens clusters and returns chunk handles to be aligned
    LOGGER.info("building consens clusters")        
    clustbits = build_reads_file(data)

    ## muscle align the consens reads and creates hdf5 indel array
    LOGGER.info("muscle alignment & building indel database")
    multi_muscle_align(data, samples, clustbits, ipyclient)

    ## builds the final HDF5 array which includes three main keys
    ## /catg -- contains all indiv catgs and has indels inserted
    ##   .attr['samples'] = [samples]
    ## /filters -- empty for now, will be filled in step 7
    ##   .attr['filters'] = [f1, f2, f3, f4, f5]
    ## /seqs -- contains the clustered sequence data as string arrays
    ##   .attr['samples'] = [samples]
    LOGGER.info("building full database")    
    ## calls singlecat func inside
    build_h5_array(data, samples, ipyclient)

    ## do we need to load in singleton clusters?
    ## invarcats()
    ## invarcats()


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


