#!/usr/bin/env python2

""" cluster across samples using vsearch with options for
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0301
# pylint: disable=R0914
# pylint: disable=R0915

import os
import pty
import gzip
import glob
import h5py
import time
import numba
import shutil
import random
import select
import datetime
import itertools
import subprocess
import numpy as np
import ipyrad as ip
from ipyrad.assemble.util import *
from ipyrad.assemble.cluster_within import muscle_call, parsemuscle


import logging
LOGGER = logging.getLogger(__name__)



def muscle_align_across(data, samples, chunk):
    """
    Reads in a chunk of the rebuilt clusters. Aligns clusters and writes
    tmpaligned chunk of seqs. Also fills a tmparray with indels locations.
    tmpchunk and tmparr will be compiled by multi_muscle_align func.
    """

    ## snames to ensure sorted order
    samples.sort(key=lambda x: x.name)
    snames = [sample.name for sample in samples]

    ## data are already chunked, read in the whole thing
    with open(chunk, 'rb') as infile:
        clusts = infile.read().split("//\n//\n")[:-1]

    ## tmparray to store indel information, for super huge data it may
    ## be necessary to build h5 here instead of below. Quick test shows that
    ## (300, 2M, 300) will give a memory error.
    maxlen = data._hackersonly["max_fragment_length"] + 20
    indels = np.zeros((len(samples), len(clusts), maxlen), dtype=np.bool_)
    duples = np.zeros(len(clusts), dtype=np.bool_)

    ## iterate over clusters and align, record excluded loci containing dups
    out = []
    for ldx in xrange(len(clusts)):
        stack = []
        lines = clusts[ldx].strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]

        ## append counter to end of names b/c muscle doesn't retain order
        names = [">{};*{}".format(j[1:], i) for i, j in enumerate(names)]

        ## don't bother aligning singletons
        if len(names) <= 1:
            pass
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
                aseqs = ["{}nnnn{}".format(i, j) for i, j in zip(aseqs1, aseqs2)]
                aseqs = np.array([list(i) for i in aseqs])
                thislen = min(maxlen, aseqs.shape[1])
                for idx in xrange(aseqs.shape[0]):
                    newn = anames[idx].rsplit(';', 1)[0]
                    ## save to aligned cluster
                    stack.append("{}\n{}".format(newn, aseqs[idx, :thislen].tostring()))
                    ## name index from sorted list (indels order)
                    sidx = snames.index(anames[idx].rsplit("_", 1)[0])
                    ## store the indels
                    LOGGER.info(np.where(aseqs[idx, :thislen] == "-")[0])
                    indels[sidx, ldx, :thislen] = aseqs[idx, :thislen] == "-"

            else:
                seqs = [i.replace('nnnn','') for i in seqs]
                string1 = muscle_call(data, names, seqs)
                anames, aseqs = parsemuscle(data, string1)
                ## aseqs is the length of the data
                aseqs = np.array([list(i) for i in aseqs])
                ## this len is at most maxlen
                thislen = min(maxlen, aseqs.shape[1])
                ## get sname indexes of samples, also used to check for dups
                sidxs = [snames.index(anames[idx].rsplit("_", 1)[0]) for idx \
                         in xrange(aseqs.shape[0])]
                ## if no dups
                if len(set(sidxs)) != aseqs.shape[0]:
                    duples[ldx] = 1
                ## arrange data into the correct order to input to indels arr
                for idx in xrange(aseqs.shape[0]):
                    newn = anames[idx].rsplit(';', 1)[0]
                    ## name index from sorted list (indels order)
                    sidx = sidxs[idx] #snames.index(anames[idx].rsplit("_", 1)[0])
                    ## save to aligned cluster
                    stack.append("{}\n{}".format(newn, aseqs[idx, :thislen].tostring()))
                    ## store the indels
                    indels[sidx, ldx, :thislen] = aseqs[idx, :thislen] == "-"

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    odx = chunk.rsplit("_")[-1]
    alignfile = os.path.join(data.tmpdir, "align_{}.fa".format(odx))
    with open(alignfile, 'wb') as outfile:
        outfile.write("\n//\n//\n".join(out)+"\n")
        os.remove(chunk)

    ## save indels array to tmp dir
    ifile = os.path.join(data.tmpdir, "indels_{}.tmp.npy".format(odx))
    np.save(ifile, indels)
    dfile = os.path.join(data.tmpdir, "duples_{}.tmp.npy".format(odx))
    np.save(dfile, duples)
    #with h5py.File(ifile, 'w') as iofile:
    #    iofile.create_dataset('indels', data=indels)




def multi_muscle_align(data, samples, clustbits, ipyclient):
    """
    Sends the cluster bits to nprocessors for muscle alignment. They return
    with indel.h5 handles to be concatenated into a joint h5.
    """
    LOGGER.info("starting alignments")

    ## get client
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, " aligning clusters     | {} | s6 |".format(elapsed))

    ## submit clustbits as jobs to engines
    jobs = {}
    for idx in xrange(len(clustbits)):
        args = [data, samples, clustbits[idx]]
        jobs[idx] = lbview.apply(muscle_align_across, *args)
    allwait = len(jobs)
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, " aligning clusters     | {} | s6 |".format(elapsed))

    ## print progress while bits are aligning
    while 1:
        finished = [i.ready() for i in jobs.values()]
        fwait = sum(finished)
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(allwait, fwait,
                    " aligning clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)
        if all(finished):
            break

    ## check for errors in muscle_align_across
    keys = jobs.keys()
    for idx in keys:
        if not jobs[idx].successful():
            LOGGER.error("error in muscle_align_across %s", jobs[idx].exception())
            raise IPyradWarningExit("error in step 6 %s", jobs[idx].exception())
        del jobs[idx]
    print("")



def concatclusts(outhandle, alignbits):
    """ concatenates sorted aligned cluster tmpfiles and removes them."""
    with gzip.open(outhandle, 'wb') as out:
        for fname in alignbits:
            with open(fname) as infile:
                out.write(infile.read()+"//\n//\n")
            #os.remove(fname)



def build_indels(data, samples, ipyclient):
    """
    Builds the indels array and catclust.gz file from the aligned clusters.
    Building catclust is very fast. Entering indels into h5 array is a bit
    slow but can probably be sped up. (todo).
    NOT currently parallelized.
    """

    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 0, " database indels       | {} | s6 |".format(elapsed))

    ## get file handles
    indelfiles = glob.glob(os.path.join(data.tmpdir, "indels_*.tmp.npy"))
    alignbits = glob.glob(os.path.join(data.tmpdir, "align_*.fa"))

    ## sort into input order by chunk names
    indelfiles.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))
    alignbits.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-3]))
    LOGGER.info("indelfiles %s", indelfiles)
    LOGGER.info("alignbits %s", alignbits)
    chunksize = int(indelfiles[0].rsplit("_", 1)[-1][:-8])

    ## concatenate finished seq clusters into a tmp file
    outhandle = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    async = ipyclient[0].apply(concatclusts, *(outhandle, alignbits))
    while not async.ready():
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(alignbits), 0, " database indels       | {} | s6 |".format(elapsed))
        time.sleep(0.1)

    ## get dims for full indel array
    maxlen = data._hackersonly["max_fragment_length"] + 20
    LOGGER.info("maxlen inside build is %s", maxlen)

    ## INIT TEMP INDEL ARRAY
    ## build an indel array for ALL loci in cat.clust.gz,
    ## chunked so that individual samples can be pulled out
    ipath = os.path.join(data.dirs.consens, data.name+".tmp.indels")
    io5 = h5py.File(ipath, 'w')
    iset = io5.create_dataset("indels",
                              shape=(len(samples), data.nloci, maxlen),
                              dtype=np.bool_,
                              chunks=(1, chunksize, maxlen))
                              #compression="gzip") this is temp why bother zip?
                              #chunks=(len(samples), min(1000, data.nloci), maxlen),
    ## again make sure names are ordered right
    samples.sort(key=lambda x: x.name)

    #iset.attrs["chunksize"] = (1, data.nloci, maxlen)
    iset.attrs["samples"] = [i.name for i in samples]

    ## enter all tmpindel arrays into full indel array
    done = 0
    init = 0
    for indf in indelfiles:
        end = int(indf.rsplit("_", 1)[-1][:-8])
        inarr = np.load(indf)
        LOGGER.info('inarr shape %s', inarr.shape)
        LOGGER.info('iset shape %s', iset[:].shape)
        iset[:, init:end, :] = inarr[:, :end-init]
        init += end-init
        done += 1
        #os.remove(indf)

        ## continued progress bar from multi_muscle_align
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(alignbits), done, " database indels       | {} | s6 |".format(elapsed))

    io5.close()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 100, " database indels       | {} | s6 |".format(elapsed))
    print("")



def cluster(data, noreverse):
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
    cov = 0.90
    if data.paramsdict["datatype"] == "gbs":
        strand = "both"
        cov = 0.60
    elif data.paramsdict["datatype"] == "pairgbs":
        strand = "both"
        cov = 0.90

    ## get call string. Thread=0 means all (default) but we may want to set
    ## an upper limit otherwise threads=60 could clobber RAM on large dsets.
    ## old userfield: -userfields query+target+id+gaps+qstrand+qcov" \
    cmd = [ip.bins.vsearch,
           "-cluster_smallmem", cathaplos,
           "-strand", strand,
           "-query_cov", str(cov),
           "-minsl", str(0.5),
           "-id", str(data.paramsdict["clust_threshold"]),
           "-userout", uhaplos,
           "-notmatched", hhaplos,
           "-userfields", "query+target+qstrand",
           "-maxaccepts", "1",
           "-maxrejects", "0",
           "-fasta_width", "0",
           "-threads", "0",
           "-fulldp",
           "-usersort",
           "-log", logfile]

    ## override reverse clustering option
    if noreverse:
        strand = "plus"  # -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    try:
        LOGGER.info(cmd)
        start = time.time()

        (dog, owner) = pty.openpty()
        proc = subprocess.Popen(cmd, stdout=owner, stderr=owner,
                                     close_fds=True)
        done = 0
        while 1:
            #dat = os.read(dog, 80192)
            isdat, _, _ = select.select([dog], [], [], 0)
            if isdat:
                dat = os.read(dog, 80192)
            else:
                dat = ""

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
                " clustering across     | {} | s6 |".format(elapsed))

        ## another catcher to let vsearch cleanup after clustering is done
        proc.wait()
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 100,
                    " clustering across     | {} | s6 |".format(elapsed))

    except subprocess.CalledProcessError as inst:
        raise IPyradWarningExit("""
        Error in vsearch: \n{}\n{}""".format(inst, subprocess.STDOUT))
    except OSError as inst:
        raise IPyradWarningExit("""
        Failed to allocate pty: \n{}""".format(inst))

    finally:
        ## progress bar
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 100, " clustering across     | {} | s6 |".format(elapsed))
        #if data._headers:
        print("")
        data.stats_files.s6 = logfile
        ## cleanup processes
        del proc, dog, owner



def build_h5_array(data, samples, ipyclient):
    """
    Build full catgs file with imputed indels. catg array of prefiltered loci
    is 4-dimensional (Big), so one big array would overload memory, we need
    to it in slices. Calls multicat (singlecat), and fill_superseqs.
    """

    ## sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    ## get maxlen dim
    maxlen = data._hackersonly["max_fragment_length"] + 20
    LOGGER.info("maxlen inside build_h5_array is %s", maxlen)

    ## open new h5 handle
    data.clust_database = os.path.join(data.dirs.consens, data.name+".clust.hdf5")
    io5 = h5py.File(data.clust_database, 'w')

    ## chunk to approximately 2 chunks per core
    chunks = ((data.nloci // (data.cpus*2)) + (data.nloci % (data.cpus*2)))

    ## Number of elements in hdf5 chunk may not exceed 500MB
    ## This is probably not actually optimal, to have such
    ## enormous chunk sizes, could probably explore efficiency
    ## of smaller chunk sizes on very very large datasets
    chunklen = chunks * len(samples) * maxlen * 4
    while chunklen > int(500e6):
        chunks = (chunks // 2) + (chunks % 2)
        chunklen = chunks * len(samples) * maxlen * 4

    data.chunks = chunks
    LOGGER.info("data.nloci is %s", data.nloci)
    LOGGER.info("chunks is %s", data.chunks)

    ## INIT FULL CATG ARRAY
    ## store catgs with a .10 loci chunk size
    supercatg = io5.create_dataset("catgs", (data.nloci, len(samples), maxlen, 4),
                                    dtype=np.uint32,
                                    chunks=(chunks, 1, maxlen, 4),
                                    compression="gzip")
    superseqs = io5.create_dataset("seqs", (data.nloci, len(samples), maxlen),
                                    dtype="|S1",
                                    chunks=(chunks, len(samples), maxlen),
                                    compression='gzip')
    superalls = io5.create_dataset("nalleles", (data.nloci, len(samples)),
                                    dtype=np.uint8,
                                    chunks=(chunks, len(samples)),
                                    compression="gzip")
    ## allele count storage
    supercatg.attrs["chunksize"] = (chunks, 1, maxlen, 4)
    supercatg.attrs["samples"] = [i.name for i in samples]
    superseqs.attrs["chunksize"] = (chunks, len(samples), maxlen)
    superseqs.attrs["samples"] = [i.name for i in samples]
    superalls.attrs["chunksize"] = (chunks, len(samples))
    superalls.attrs["samples"] = [i.name for i in samples]

    ## array for pair splits locations, dup and ind filters
    io5.create_dataset("splits", (data.nloci, ), dtype=np.uint16)
    io5.create_dataset("duplicates", (data.nloci, ), dtype=np.bool_)

    ## close the big boy
    io5.close()

    ## FILL SUPERCATG and fills dupfilter, indfilter, and nalleles
    multicat(data, samples, ipyclient)

    ## FILL SUPERSEQS and fills edges(splits) for paired-end data
    fill_superseqs(data, samples)

    ## set sample states
    for sample in samples:
        sample.stats.state = 6



def fill_dups_arr(data):
    """
    fills the duplicates array from the multi_muscle_align tmp files
    """
    ## build the duplicates array
    duplefiles = glob.glob(os.path.join(data.tmpdir, "duples_*.tmp.npy"))
    duplefiles.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-8]))

    ## enter the duplicates filter into super h5 array
    io5 = h5py.File(data.clust_database, 'r+')
    dfilter = io5["duplicates"]

    ## enter all duple arrays into full duplicates array
    init = 0
    for dupf in duplefiles:
        end = int(dupf.rsplit("_", 1)[-1][:-8])
        inarr = np.load(dupf)
        dfilter[init:end] = inarr
        init += end-init
        os.remove(dupf)
    #del inarr

    ## continued progress bar
    LOGGER.info("all duplicates: %s", dfilter[:].sum())
    io5.close()



def get_seeds_and_hits(uhandle, bseeds, snames):
    """
    builds a seeds and hits (uarr) array of ints from the utemp.sort file.
    Saves outputs to files ...
    """
    ## read in the utemp.sort file
    updf = np.loadtxt(uhandle, dtype="S")

    ## Get seeds for all matches from usort
    seeds = np.unique(updf[:, 1])
    seedsarr = np.column_stack([
                    np.arange(len(seeds)),
                   [i.rsplit("_", 1)[0] for i in seeds],
                   [i.rsplit("_", 1)[1] for i in seeds]])
    seedsarr[:, 1] = [snames.index(i) for i in seedsarr[:, 1]]
    seedsarr = seedsarr.astype(np.int64)
    LOGGER.info("got a seedsarr %s", seedsarr.shape)

    ## Get matches from usort and create an array for fast entry
    uarr = np.zeros((updf.shape[0], 3), dtype=np.int64)
    idx = -1
    lastloc = None
    for ldx in xrange(updf.shape[0]):
        tloc = updf[ldx, 1]
        if tloc != lastloc:
            idx += 1
        uarr[ldx, 0] = idx
        lastloc = tloc
    ## create a column with sample index
    uarr[:, 1] = [int(snames.index(i.rsplit("_", 1)[0])) for i in updf[:, 0]]
    ## create a column with only consens index for sample
    uarr[:, 2] = [int(i.rsplit("_", 1)[1]) for i in updf[:, 0]]
    uarr = uarr.astype(np.int64)
    LOGGER.info("got a uarr %s", uarr.shape)

    ## save as h5 to we can grab by sample slices
    with h5py.File(bseeds, 'w') as io5:
        io5.create_dataset("seedsarr", data=seedsarr, dtype=np.int64)
        io5.create_dataset("uarr", data=uarr, dtype=np.int64)



def multicat(data, samples, ipyclient):
    """
    Runs singlecat and cleanup jobs for each sample.
    For each sample this fills its own hdf5 array with catg data & indels.
    This is messy, could use simplifiying.
    """

    ## notes for improvements:
    LOGGER.info("in the multicat")
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    progressbar(20, 0, " indexing clusters     | {} | s6 |".format(elapsed))

    ## parallel client
    lbview = ipyclient.load_balanced_view()
    ## First submit a sleeper job as temp_flag for cleanups
    last_sample = 0
    cleanups = {}
    cleanups[last_sample] = lbview.apply(time.sleep, 0.0)

    ## get samples and names, sorted
    snames = [i.name for i in samples]
    snames.sort()

    ## Build an array for quickly indexing consens reads from catg files.
    ## save as a npy int binary file.
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp.sort")
    bseeds = os.path.join(data.dirs.consens, data.name+".tmparrs.h5")

    ## send as first async1 job
    async1 = lbview.apply(get_seeds_and_hits, *(uhandle, bseeds, snames))
    async2 = lbview.apply(fill_dups_arr, data)

    ## progress
    while not (async1.ready() and async2.ready()):
        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        progressbar(20, 0, " indexing clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)
    if not async1.successful():
        raise IPyradWarningExit("error in get_seeds: %s", async1.exception())
    if not async2.successful():
        raise IPyradWarningExit("error in fill_dups: %s", async2.exception())

    ## make a limited njobs view based on mem limits
    smallview = ipyclient.load_balanced_view(targets=ipyclient.ids[::4])

    ## make a list of jobs
    jobs = {}
    for sample in samples:
        ## grab just the hits for this sample
        sidx = snames.index(sample.name)
        jobs[sample.name] = smallview.apply(singlecat, *(data, sample, bseeds, sidx))
        #jobs[sample.name] = lbview.apply(singlecat, *(data, sample, bseeds, sidx))

    ## check for finished and submit disk-writing job when finished
    allwait = len(jobs)
    fwait = 0
    while 1:
        ## check for new jobs to submit
        jkeys = jobs.keys()
        for name in jkeys:
            if jobs[name].ready():
                if jobs[name].successful():
                    ## get args for writing
                    args = (data, data.samples[name], snames.index(name))
                    ## store that job has finished
                    fwait += 1
                    # submit writing job
                    LOGGER.info("sending %s:%s after %s", snames.index(name), name, last_sample)
                    with lbview.temp_flags(after=cleanups[last_sample]):
                        cleanups[name] = lbview.apply(write_to_fullarr, *args)
                    last_sample = name
                    del jobs[name]
                else:
                    LOGGER.error(" error in singlecat (%s) %s",
                                 name, jobs[name].exception())
                    raise IPyradWarningExit(" error in singlecat ({}) {}"\
                                        .format(name), jobs[name].exception())

        ## print progress bar
        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        progressbar(allwait, fwait, " indexing clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)
        if fwait == allwait:
            break

    ## wait for writing jobs to finish
    print("")
    start = time.time()
    while 1:
        finished = [i for i in cleanups.values() if i.ready()]
        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        progressbar(len(cleanups), len(finished),
                    " building database     | {} | s6 |".format(elapsed))
        time.sleep(0.1)
        ## break if one failed, or if finished
        if not all([i.successful() for i in finished]):
            break
        if len(cleanups) == len(finished):
            break

    ## check for errors
    for job in cleanups:
        if cleanups[job].ready():
            if not cleanups[job].successful():
                LOGGER.error(" error in write_to_fullarr (%s) %s",
                             job, cleanups[job].exception())
                raise IPyradWarningExit(" error in write_to_fullarr ({}) {}"\
                                    .format(job, cleanups[job].exception()))

    ## remove large indels array file
    ifile = os.path.join(data.dirs.consens, data.name+".tmp.indels")
    if os.path.exists(ifile):
        os.remove(ifile)

    ## print final progress
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    progressbar(10, 10, " building database     | {} | s6 |".format(elapsed))
    print("")



## This is where indels are imputed
def singlecat(data, sample, bseeds, sidx):
    """
    Orders catg data for each sample into the final locus order. This allows
    all of the individual catgs to simply be combined later. They are also in
    the same order as the indels array, so indels are inserted from the indel
    array that is passed in.
    """

    ## grab seeds and hits info for this sample
    with h5py.File(bseeds, 'r') as io5:
        ## get hits just for this sample and sort them by sample order index
        hits = io5["uarr"][:]
        hits = hits[hits[:, 1] == sidx, :]
        #hits = hits[hits[:, 2].argsort()]
        ## get seeds just for this sample and sort them by sample order index
        seeds = io5["seedsarr"][:]
        seeds = seeds[seeds[:, 1] == sidx, :]
        #seeds = seeds[seeds[:, 2].argsort()]
        full = np.concatenate((seeds, hits))
        full = full[full[:, 0].argsort()]

    ## still using max+30 len limit, rare longer merged reads get trimmed
    ## we need to allow room for indels to be added too
    maxlen = data._hackersonly["max_fragment_length"] + 20

    ## we'll fill a new catg and alleles arr for this sample in locus order,
    ## which is known from seeds and hits
    ocatg = np.zeros((data.nloci, maxlen, 4), dtype=np.uint32)
    onall = np.zeros(data.nloci, dtype=np.uint8)

    LOGGER.info("single cat here")
    ## grab the sample's data and write to ocatg and onall
    with h5py.File(sample.files.database, 'r') as io5:
        ## get it and delete it
        catarr = io5["catg"][:]
        #LOGGER.info("catarr shape %s", catarr.shape)
        #LOGGER.info("ocatg shape %s", ocatg.shape)
        tmp = catarr[full[:, 2], :maxlen, :]
        del catarr
        ocatg[full[:, 0], :tmp.shape[1], :] = tmp
        del tmp

        ## get it and delete it
        nall = io5["nalleles"][:]
        onall[full[:, 0]] = nall[full[:, 2]]
        del nall

    ## get indel locations for this sample
    ipath = os.path.join(data.dirs.consens, data.name+".tmp.indels")
    with h5py.File(ipath, 'r') as ih5:
        indels = ih5["indels"][sidx, :, :maxlen]
        #LOGGER.info("h5 indels shape %s", ih5["indels"].shape)
        #LOGGER.info("indels shape %s", indels.shape)

    ## insert indels into ocatg
    newcatg = inserted_indels(indels, ocatg)
    del ocatg, indels

    ## save big arrays to disk temporarily
    smpio = os.path.join(data.dirs.consens, sample.name+'.tmp.h5')
    if os.path.exists(smpio):
        os.remove(smpio)
    with h5py.File(smpio, 'w') as oh5:
        oh5.create_dataset("icatg", data=newcatg, dtype=np.uint32)
                           #chunks=(chunksize, maxlen, 4), gzip=True)
        oh5.create_dataset("inall", data=onall, dtype=np.uint8)
                           #chunks=(chunksize, maxlen, 4), gzip=True)
    # LOGGER.info("submitting %s %s to singlecat", sidx, hits.shape)
    #del ocatg
    #del onall

    #return newcatg, onall, sidx



#def write_to_fullarr(data, newcatg, onall, sidx):
def write_to_fullarr(data, sample, sidx):
    """ writes arrays to h5 disk """
    ## save big arrays to disk temporarily


    with h5py.File(data.clust_database, 'r+') as io5:
        chunk = io5["catgs"].attrs["chunksize"][0]
        catg = io5["catgs"]
        nall = io5["nalleles"]

        ## adding an axis to newcatg makes it write about 1000X faster.
        ## so instead of e.g.,
        ##        newcatg (1000, 150, 4) -> catg (1000, 150, 4)
        ## we do:
        ##        newcatg (1000, 1, 150, 4) -> catg (1000, 1, 150, 4)
        smpio = os.path.join(data.dirs.consens, sample.name+'.tmp.h5')
        with h5py.File(smpio) as indat:
            newcatg = indat["icatg"][:]
            onall = indat["inall"][:]
            for cidx in xrange(0, catg.shape[0], chunk):
                #LOGGER.info(catg.shape, nall.shape, newcatg.shape, onall.shape, sidx)
                end = cidx + chunk
                catg[cidx:end, sidx:sidx+1, :] = np.expand_dims(newcatg[cidx:end, :], axis=1)
                nall[:, sidx:sidx+1] = np.expand_dims(onall, axis=1)
        os.remove(smpio)


            # for cidx in xrange(0, catg.shape[0], chunks):
            #     #LOGGER.info(catg.shape, nall.shape, newcatg.shape, onall.shape, sidx)
            #     catg[init:cidx, sidx:sidx+1, :] = np.expand_dims(newcatg[init:cidx, :], axis=1)
            #     nall[:, sidx:sidx+1] = np.expand_dims(onall, axis=1)
            #     init = cidx
    # smpio = os.path.join(data.dirs.consens, sample.name+'.tmp.h5')
    # if os.path.exists(smpio):
    #     os.remove(smpio)
    # with h5py.File(smpio, 'w') as oh5:
    #     oh5.create_dataset("icatg", data=newcatg, dtype=np.uint32)
    #     oh5.create_dataset("inall", data=onall, dtype=np.uint8)

    # with h5py.File(data.clust_database, 'r+') as io5:
    #     catg = io5["catgs"]
    #     nall = io5["nalleles"]

    #     #LOGGER.info(catg.shape, nall.shape, newcatg.shape, onall.shape, sidx)
    #     catg[:, sidx, :] = newcatg[:]
    #     nall[:, sidx] = onall[:]


# def insert_and_cleanup(data, sname, sidx, newcatg, onall):
#     """
#     Enter results from singlecat into the super h5 arrays. This is
#     not parallelized. Why isn't it faster?
#     """

#     ## grab supercatg from super and get index of this sample
#     with h5py.File(data.clust_database, 'r+') as io5:
#         catg = io5["catgs"]
#         nall = io5["nalleles"]

#         #smpio = os.path.join(data.dirs.consens, sname+".tmp.h5")
#         #with h5py.File(smpio, 'r') as smp5:
#         #    newcatg = smp5["icatg"]
#         #    onall = smp5["inall"]

#             ## get this samples index in the h5
#             sidx = list(catg.attrs["samples"]).index(sname)
#             LOGGER.info("insert & cleanup : %s sidx: %s", sname, sidx)

#             ## single fill, does this crush memory? should we chunk?
#             ## I think it's slow because it's entering into a column
#             ## but catg is chunked by rows... C vs. F, amiright?
#             catg[:, sidx, :] = newcatg[:]
#             nall[:, sidx] = onall[:]

            ## FILL SUPERCATG -- catg is chunked by nchunk loci
            # chunk = catg.attrs["chunksize"]
            # for chu in xrange(0, catg.shape[0], chunk):
            #     catg[chu:chu+chunk, sidx, :] = newcatg[chu:chu+chunk]

            # ## FILL allelic information here.
            # for chu in xrange(0, catg.shape[0], chunk):
            #     nall[chu:chu+chunk, sidx] = onall[chu:chu+chunk]



# @numba.jit(nopython=True)
# def fill_duplicates_parallel(base, block, seedsarr, uarr, ntax):
#     """
#     Applies filter to block of loci. Requires numba v.0.28.
#     """
#     ## filter for duplicates
#     size = min(block, seedsarr.shape[0] - (base))
#     dfilter = np.zeros(size, dtype=np.bool_)

#     ## fill dfilter
#     for bidx in xrange(size):
#         idx = bidx + base
#         ## get idx of seed and matches
#         sidxs = uarr[uarr[:, 0] == idx, 1]
#         seedx = seedsarr[seedsarr[:, 0] == idx, 1]
#         bins = np.concatenate((sidxs, seedx))
#         counts = np.bincount(bins)
#         if counts.max() > 1:
#             dfilter[bidx] = True
#     return dfilter



# @numba.jit(nopython=True)
# def fill_cats_and_alleles(maxlen, seeds, hits, catarr, nall):
#     """
#     Fill catg and nalleles data in locus sorted order and returns
#     the filled arrays. Uses the seeds and hits arrays. Needs the entire
#     catarr and nall arrays loaded into memory, unfortunately.
#     """
#     ocatg = np.zeros((10000, maxlen, 4), dtype=np.uint32)
#     onall = np.zeros(10000, dtype=np.uint8)

#     ## fill the locus data where sample was a seed
#     #for idx in xrange(seeds.shape[0]):
#     for idx in xrange(hslice, hslice+10000):
#         ## this is the samples data
#         cidx = seeds[idx, 2]

#         ## we can fill it into the array if its in this 10K chunk
#         if cidx in seeds[:, ]

#         ## is this sample the seed?
#         if idx in seeds[:, 0]:

#         loc = seeds[idx, 0]
#         cidx = seeds[idx, 2]

#         ## set locus with sample data
#         ocatg[loc, :catarr.shape[1]] = catarr[cidx]
#         onall[loc] = nall[cidx]

#     ## fill the locus data where sample is a hit
#     for idx in xrange(hits.shape[0]):
#         ## get the locus number
#         loc = hits[idx, 0]
#         cidx = hits[idx, 2]

#         ## set locus with sample data
#         ocatg[loc, :catarr.shape[1]] = catarr[cidx]
#         onall[loc] = nall[cidx]

#     return ocatg, onall



@numba.jit(nopython=True)
def inserted_indels(indels, ocatg):
    """
    inserts indels into the catg array
    """
    ## return copy with indels inserted
    newcatg = np.zeros(ocatg.shape, dtype=np.uint32)

    ## iterate over loci and make extensions for indels
    for iloc in xrange(ocatg.shape[0]):
        ## get indels indices
        indidx = np.where(indels[iloc, :])[0]
        if np.any(indidx):
            ## which new (empty) rows will be added
            allrows = np.arange(ocatg.shape[1])
            mask = np.ones(allrows.shape[0], dtype=np.bool_)
            for idx in indidx:
                mask[idx] = False
            not_idx = allrows[mask == 1]

            ## fill in new data into all other spots
            newcatg[iloc][not_idx] = ocatg[iloc, :not_idx.shape[0]]
        else:
            newcatg[iloc] = ocatg[iloc]
    return newcatg



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
    maxlen = data._hackersonly["max_fragment_length"] + 20
    LOGGER.info("maxlen inside fill_superseqs is %s", maxlen)

    ## data has to be entered in blocks
    infile = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    clusters = gzip.open(infile, 'r')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## iterate over clusters
    chunks = superseqs.attrs["chunksize"]
    chunksize = chunks[0]
    done = 0
    iloc = 0
    cloc = 0
    chunkseqs = np.zeros(chunks, dtype="|S1")
    chunkedge = np.zeros(chunksize, dtype=np.uint16)

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
            separator = np.where(np.all(seqs == 'n', axis=0))[0]
            if np.any(separator):
                chunkedge[cloc] = separator.min()

            ## fill in the hits
            #LOGGER.info("seqs : %s", seqs)
            #LOGGER.info("seqs.shape : %s", seqs.shape)
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
    clusters.close()

    ## edges is filled with splits for paired data.
    LOGGER.info("done filling superseqs")

    ## close handle
    os.remove(infile)



def count_seeds(usort):
    """
    uses bash commands to quickly count N seeds from utemp file
    """
    with open(usort, 'r') as insort:
        cmd1 = ["cut", "-f", "2"]
        cmd2 = ["uniq"]
        cmd3 = ["wc"]
        proc1 = subprocess.Popen(cmd1, stdin=insort, stdout=subprocess.PIPE)
        proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)
        proc3 = subprocess.Popen(cmd3, stdin=proc2.stdout, stdout=subprocess.PIPE)
        res = proc3.communicate()
        nseeds = int(res[0].split()[0])
        proc1.stdout.close()
        proc2.stdout.close()
        proc3.stdout.close()
    return nseeds



def build_clustbits(data, ipyclient):
    """
    Reconstitutes clusters from .utemp and htemp files and writes them
    to chunked files for aligning in muscle.
    """

    ## parallel client
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(3, 0, " building clusters     | {} | s6 |".format(elapsed))
    LOGGER.info("building reads file -- loading utemp file into mem")

    ## send sort job to engines. Sorted seeds allows us to work through
    ## the utemp file one locus at a time instead of reading all into mem.
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    usort = os.path.join(data.dirs.consens, data.name+".utemp.sort")
    cmd = ["sort", "-k", "2", uhandle, "-o", usort]
    async1 = subprocess.Popen(cmd)
    while async1.poll() == None:
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(3, 0, " building clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)

    ## send count seeds job to engines.
    async2 = ipyclient[0].apply(count_seeds, usort)
    while not async2.ready():
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(3, 1, " building clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)

    ## wait for both to finish while printing progress timer
    nseeds = async2.result()

    ## send the clust bit building job to work and track progress
    async3 = ipyclient[0].apply(sub_build_clustbits, *(data, usort, nseeds))
    while not async3.ready():
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(3, 2, " building clusters     | {} | s6 |".format(elapsed))
        time.sleep(0.1)
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(3, 3, " building clusters     | {} | s6 |".format(elapsed))
    print("")

    ## return the nloci and clustbits
    clustbits, nloci = async3.result()
    return clustbits, nloci



def sub_build_clustbits(data, usort, nseeds):
    """
    A subfunction of build_clustbits to allow progress tracking. This func
    splits the unaligned clusters into bits for aligning on separate cores.
    """

    ## load FULL concat fasta file into a dict. This could cause RAM issues.
    ## this file has iupac codes in it, not ambigs resolved, and is gzipped.
    LOGGER.info("loading full _catcons file into memory")
    allcons = {}
    conshandle = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    with gzip.open(conshandle, 'rb') as iocons:
        cons = itertools.izip(*[iter(iocons)]*2)
        for namestr, seq in cons:
            nnn, sss = [i.strip() for i in namestr, seq]
            allcons[nnn[1:]] = sss

    ## set optim to approximately 4 chunks per core. Smaller allows for a bit
    ## cleaner looking progress bar. 40 cores will make 160 files.
    optim = ((nseeds // (data.cpus*4)) + (nseeds % (data.cpus*4)))
    LOGGER.info("building clustbits, optim=%s, nseeds=%s, cpus=%s",
                optim, nseeds, data.cpus)

    ## iterate through usort grabbing seeds and matches
    with open(usort, 'rb') as insort:
        ## iterator, seed null, and seqlist null
        isort = iter(insort)
        loci = 0
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0

        while 1:
            ## grab the next line
            try:
                hit, seed, ori = isort.next().strip().split()
            except StopIteration:
                break

            ## if same seed, append match
            if seed != lastseed:

                ## store the last fseq, count it, and clear it
                if fseqs:
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                ## occasionally write to file
                if seqsize >= optim:
                    if seqlist:
                        loci += seqsize
                        with open(os.path.join(data.tmpdir,
                            data.name+".chunk_{}".format(loci)), 'w') as clustsout:
                            clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")
                        ## reset list and counter
                        seqlist = []
                        seqsize = 0

                ## store the new seed on top of fseq
                fseqs.append(">{}\n{}".format(seed, allcons[seed]))
                lastseed = seed

            ## add match to the seed
            seq = allcons[hit]
            ## revcomp if orientation is reversed
            if ori == "-":
                seq = fullcomp(seq)[::-1]
            fseqs.append(">{}\n{}".format(hit, seq))

    ## write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
        seqsize += 1
        loci += seqsize
    if seqlist:
        with open(os.path.join(data.tmpdir,
            data.name+".chunk_{}".format(loci)), 'w') as clustsout:
            clustsout.write("\n//\n//\n".join(seqlist)+"\n//\n//\n")

    ## final progress and cleanup
    del allcons
    clustbits = glob.glob(os.path.join(data.tmpdir, data.name+".chunk_*"))

    ## return stuff
    return clustbits, loci



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

    ## concatenate all of the gzipped consens files
    cmd = ['cat'] + conshandles
    allcons = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    LOGGER.debug(" ".join(cmd))
    with open(allcons, 'w') as output:
        call = subprocess.Popen(cmd, stdout=output)
        call.communicate()

    ## a string of sed substitutions for temporarily replacing hetero sites
    ## skips lines with '>', so it doesn't affect taxon names
    subs = ["/>/!s/W/A/g", "/>/!s/w/A/g", "/>/!s/R/A/g", "/>/!s/r/A/g",
            "/>/!s/M/A/g", "/>/!s/m/A/g", "/>/!s/K/T/g", "/>/!s/k/T/g",
            "/>/!s/S/C/g", "/>/!s/s/C/g", "/>/!s/Y/C/g", "/>/!s/y/C/g"]
    subs = ";".join(subs)

    ## impute pseudo-haplo information to avoid mismatch at hetero sites
    ## the read data with hetero sites is put back into clustered data later.
    ## pipe passed data from gunzip to sed.
    cmd1 = ["gunzip", "-c", allcons]
    cmd2 = ["sed", subs]
    LOGGER.debug(" ".join(cmd1))
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    allhaps = allcons.replace("_catcons.tmp", "_cathaps.tmp")
    with open(allhaps, 'w') as output:
        LOGGER.debug(" ".join(cmd2))
        proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=output)
        proc2.communicate()
    proc1.stdout.close()

    ## now sort the file using vsearch
    allsort = allcons.replace("_catcons.tmp", "_catsort.tmp")
    cmd1 = [ipyrad.bins.vsearch,
            "--sortbylength", allhaps,
            "--fasta_width", "0",
            "--output", allsort]
    LOGGER.debug(" ".join(cmd1))
    proc1 = subprocess.Popen(cmd1)
    proc1.communicate()

    ## shuffle sequences within size classes. Tested seed (8/31/2016)
    ## shuffling works repeatably with seed.
    random.seed(randomseed)

    ## open an iterator to lengthsorted file and grab two lines at at time
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

    ## get parallel view
    start = time.time()

    ## a chunker for writing every N loci. This uses core info, meaning that
    ## if users do not supply -n arg then it might be poorly estimated.
    ## if no info we use detect_cpus to get info for this node.
    data.cpus = data._ipcluster["cores"]
    if not data.cpus:
        ## should this use detect_cpus? more appropriate would be len ipyclient
        #data.cpus = detect_cpus()
        data.cpus = len(ipyclient)

    ## make a vsearch input fasta file with all samples reads concat
    binput = ipyclient[0].apply(build_input_file, *[data, samples, randomseed])
    while not binput.ready():
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(100, 0, " concat/shuffle input  | {} | s6 |".format(elapsed))
        time.sleep(0.1)
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(100, 100, " concat/shuffle input  | {} | s6 |".format(elapsed))
    print("")

    ## calls vsearch, uses all threads available to head node
    cluster(data, noreverse)

    # ## make an tmpout directory
    data.tmpdir = os.path.join(data.dirs.consens, data.name+"-tmpaligns")
    if not os.path.exists(data.tmpdir):
        os.mkdir(data.tmpdir)

    ## wrap everything involving tmpdir to make sure we delete it on failure
    try:
        ## build consens clusters and returns chunk handles to be aligned
        clustbits, nloci = build_clustbits(data, ipyclient)
        data.nloci = nloci

        ## muscle align the consens reads and creates hdf5 indel array
        multi_muscle_align(data, samples, clustbits, ipyclient)

        ## submit to indel entry
        build_indels(data, samples, ipyclient)

        ## builds the final HDF5 array which includes three main keys
        ## /catg -- contains all indiv catgs and has indels inserted
        ##   .attr['samples'] = [samples]
        ## /filters -- filled for dups, left empty for others until step 7.
        ##   .attr['filters'] = [f1, f2, f3, f4, f5]
        ## /seqs -- contains the clustered sequence data as string arrays
        ##   .attr['samples'] = [samples]
        ## /edges -- gets the paired split locations for now.
        ## /snps  -- left empty for now

        ## calls singlecat func inside
        LOGGER.info("building full database")
        build_h5_array(data, samples, ipyclient)

    #except Exception as inst:
    #    LOGGER.error(inst)
    #    raise IPyradWarningExit(inst)

    finally:
        ## delete the tmpdir
        shutil.rmtree(data.tmpdir)




if __name__ == "__main__":

    ## get path to test dir/
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__)
               )
           )
       )


    ## load test data (pairgbs)
    DATA = ip.load_json("/home/deren/Documents/RADmissing/rad1/half_min4.json")

    # ## run step 6
    DATA.run("6", force=True)
