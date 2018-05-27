#!/usr/bin/env python

""" call consensus base calls on single-end data """

# py2/3 compatible
from __future__ import print_function
try:
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

import os
import time
import gzip
import glob
import shutil
import warnings
from collections import Counter

import numpy as np
import pandas as pd
import scipy.stats
import scipy.misc

import ipyrad as ip
from .jointestimate import recal_hidepth
from .util import IPyradError, clustdealer, PRIORITY

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class Step5:
    "Organized Step 5 functions for all datatype and methods"

    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.samples = self.get_subsamples()
        self.ref = bool("reference" in data.paramsdict["assembly_method"])
        self.ipyclient = ipyclient
        self.lbview = ipyclient.load_balanced_view()
        self.setup_dirs()


    def setup_dirs(self):
        "setup directories, remove old tmp files"

        # final results dir
        self.data.dirs.consens = os.path.join(
            self.data.dirs.project, 
            "{}_consens".format(self.data.name))
        if not os.path.exists(self.data.dirs.consens):
            os.mkdir(self.data.dirs.consens)

        # tmpfile dir (zap it if it exists)
        self.data.tmpdir = os.path.join(
            self.data.dirs.project, 
            "{}-tmpdir".format(self.data.name))
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)

        # assign output file handles for s6
        for sample in self.samples:
            if not self.ref:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.gz".format(sample.name))
            else:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.sam".format(sample.name))

        # set up parallel client: allow user to throttle cpus
        self.lbview = self.ipyclient.load_balanced_view()
        if self.data._ipcluster["cores"]:
            self.ncpus = self.data._ipcluster["cores"]
        else:
            self.ncpus = len(self.ipyclient.ids)


    def get_subsamples(self):
        "Apply state, ncluster, and force filters to select samples"

        # filter samples by state
        state3 = self.data.stats.index[self.data.stats.state < 4]
        state4 = self.data.stats.index[self.data.stats.state == 4]
        state5 = self.data.stats.index[self.data.stats.state > 4]

        # tell user which samples are not ready for step5
        if state3.any():
            print("skipping samples not in state==4:\n{}"
                  .format(state3.tolist()))

        if self.force:
            # run all samples above state 3
            subs = self.data.stats.index[self.data.stats.state > 3]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 5
            if state5.any():
                print("skipping samples already finished step 5:\n{}"
                      .format(state5.tolist()))

            # run all samples in state 4
            subsamples = [self.data.samples[i] for i in state4]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.clusters_hidepth:
                checked_samples.append(sample)
            else:
                print("skipping {}; no clusters found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 5")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.clusters_hidepth,
            reverse=True,
        )

        # if sample is already done skip
        if "hetero_est" not in self.data.stats:
            for sample in checked_samples:
                sample.stats.hetero_est = 0.001
                sample.stats.error_est = 0.0001

        if self.data._headers:
            print(u"  Mean error  [{:.5f} sd={:.5f}]".format(
                self.data.stats.error_est.mean(),
                self.data.stats.error_est.std(),
            ))
            print(u"  Mean hetero [{:.5f} sd={:.5f}]".format(
                self.data.stats.hetero_est.mean(),
                self.data.stats.hetero_est.std(),
            ))
        return checked_samples


    def run(self):
        "run the main functions on the parallel client"
        # this isn't setup yet to allow restarting if interrupted mid run
        try:
            self.remote_calculate_depths()
            self.remote_make_chunks()
            statsdicts = self.remote_process_chunks()
            self.remote_concatenate_chunks()
            self.data_store(statsdicts)
        finally:
            shutil.rmtree(self.data.tmpdir)
            self.data.save()


    def remote_calculate_depths(self):
        "checks whether mindepth has changed and calc nclusters and maxlen"
        # send jobs to be processed on engines
        start = time.time()
        printstr = ("calculating depths  ", "s5")
        jobs = {}
        maxlens = []
        for sample in self.samples:
            jobs[sample.name] = self.lbview.apply(
                recal_hidepth,
                *(self.data, sample))

        # block until finished
        while 1:
            ready = [i.ready() for i in jobs.values()]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        # check for failures and collect results
        print("")
        for sample in self.samples:
            if not jobs[sample.name].successful():
                raise IPyradError(jobs[sample.name].exception())
            else:
                hidepth, maxlen, _, _ = recal_hidepth(self.data, sample)
                # (not saved) stat values are for majrule min
                sample.stats["clusters_hidepth"] = hidepth
                sample.stats_dfs.s3["clusters_hidepth"] = hidepth
                maxlens.append(maxlen)
        
        # update hackersdict with max fragement length
        self.data._hackersonly["max_fragment_length"] = max(maxlens)


    def remote_make_chunks(self):
        "split clusters into chunks for parallel processing"

        # first progress bar
        start = time.time()
        printstr = ("chunking clusters   ", "s5")

        # send off samples to be chunked
        jobs = {}
        for sample in self.samples:
            jobs[sample.name] = self.lbview.apply(
                make_chunks,
                *(self.data, sample, len(self.ipyclient)))

        ## block until finished
        while 1:
            ready = [i.ready() for i in jobs.values()]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        ## check for failures
        print("")
        for sample in self.samples:
            if not jobs[sample.name].successful():
                raise IPyradError(jobs[sample.name].exception())


    def remote_process_chunks(self):
        "process the cluster chunks into arrays and consens or bam files"

        # send chunks to be processed
        start = time.time()
        jobs = {sample.name: [] for sample in self.samples}
        printstr = ("consens calling     ", "s5")

        # submit jobs
        for sample in self.samples:
            chunks = glob.glob(os.path.join(
                self.data.tmpdir,
                "{}.chunk.*".format(sample.name)))
            chunks.sort(key=lambda x: int(x.split('.')[-1]))

            # submit jobs
            for chunk in chunks:
                jobs[sample.name].append(
                    self.lbview.apply(
                        process_chunks,
                        *(self.data, sample, chunk, self.ref)))
                
        # track progress - just wait for all to finish before concat'ing
        allsyncs = list(chain(*[jobs[i] for i in jobs]))
        while 1:
            ready = [i.ready() for i in allsyncs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        # check for failures
        for job in allsyncs:
            if not job.successful():
                ip.logger.error("error in step 5: %s", job.exception())
                raise IPyradError("error in consensus_calls():\n{}"
                                  .format(job.exception()))

        # collect all results for a sample and store stats 
        statsdicts = {}
        for sample in self.samples:
            statsdicts[sample.name] = [i.result() for i in jobs[sample.name]]
        return statsdicts


    def remote_concatenate_chunks(self):
        # concatenate and store catgs
        print("")
        start = time.time()
        printstr = ("indexing alleles    ", "s5")

        # concat catgs for each sample
        asyncs1 = {}
        for sample in self.samples:
            asyncs1[sample.name] = self.lbview.apply(
                concat_catgs,                
                *(self.data, sample, self.ref))           

        # collect all results for a sample and store stats 
        if self.ref:
            concat_job = concat_reference_consens
        else:
            concat_job = concat_denovo_consens
        asyncs2 = {}
        for sample in self.samples:
            asyncs2[sample.name] = self.lbview.apply(
                concat_job,
                *(self.data, sample))
            
        # track progress of stats storage
        alljobs = list(asyncs1.values()) + list(asyncs2.values())
        while 1:
            ready = [i.ready() for i in alljobs]
            self.data._progressbar(len(ready), sum(ready), start, printstr)
            time.sleep(0.1)
            if len(ready) == sum(ready):
                break

        # check for failures:
        print("")
        for job in alljobs:
            if not job.successful():
                ip.logger.error("error in step 5: %s", job.exception())
                raise IPyradError(job.exception())


    def data_store(self, statsdicts):
        "store assembly object stats"
        
        # store sample stats
        for sample in self.samples:
            store_sample_stats(self.data, sample, statsdicts[sample.name])

        # build Assembly stats
        self.data.stats_dfs.s5 = self.data._build_stat("s5")

        # write stats file
        self.data.stats_files.s5 = os.path.join(
            self.data.dirs.consens, 
            's5_consens_stats.txt')

        with open(self.data.stats_files.s5, 'w') as out:
            self.data.stats_dfs.s5.to_string(
                buf=out,
                formatters={
                    'clusters_total': '{:.0f}'.format,
                    'filtered_by_depth': '{:.0f}'.format,
                    'filtered_by_maxH': '{:.0f}'.format,
                    'filtered_by_maxN': '{:.0f}'.format,
                    'reads_consens': '{:.0f}'.format,
                    'nsites': '{:.0f}'.format,
                    'nhetero': '{:.0f}'.format,
                    'heterozygosity': '{:.5f}'.format
                })



def concat_catgs(data, sample, isref):
    "concat catgs into a single sample catg and remove tmp files"

    # collect tmpcat files
    tmpcats = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcats.*".format(sample.name)))
    tmpcats.sort(key=lambda x: int(x.split(".")[-1]))

    # get full nrows of the new h5 from the tmpcat filenames
    nrows = sum([int(i.rsplit(".", 2)[-2]) for i in tmpcats])

    # get shape info from the first cat, (optim, maxlen, 4)
    optim = min(nrows, 5000)
    with h5py.File(tmpcats[0], 'r') as io5:
        _, maxlen, _ = io5['cats'].shape

    # file handle for concat'd catg array
    outcatg = os.path.join(
        data.dirs.consens,
        "{}.catg".format(sample.name))

    # fill in the chunk array
    with h5py.File(outcatg, 'w') as ioh5:
        dcat = ioh5.create_dataset(
            name="catg",
            shape=(nrows, maxlen, 4),
            dtype=np.uint32,
            chunks=(optim, maxlen, 4),
            compression="gzip")
        dall = ioh5.create_dataset(
            name="nalleles", 
            shape=(nrows, ),
            dtype=np.uint8,
            chunks=(optim, ),
            compression="gzip")
        
        # only create chrom for reference-aligned data
        if isref:
            dchrom = ioh5.create_dataset(
                name="chroms",
                shape=(nrows, 3),
                dtype=np.int64,
                chunks=(optim, 3),
                compression="gzip")

        # Combine all those tmp cats into the big cat
        start = 0
        for icat in tmpcats:
            addon = int(icat.rsplit(".", 2)[-2])
            end = start + addon
            io5 = h5py.File(icat, 'r')
            dcat[start:end] = io5['cats'][:]
            dall[start:end] = io5['alls'][:]
            if isref:
                dchrom[start:end] = io5['chroms'][:]
            start = end
            io5.close()
            os.remove(icat)


def concat_denovo_consens(data, sample):
    "concatenate consens bits into fasta file for denovo assemblies"

    # collect consens chunk files
    combs1 = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcons.*".format(sample.name)))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    # write to the file
    with gzip.open(sample.files.consens, 'wt') as out:
        for fname in combs1:
            with open(fname) as infile:
                out.write(infile.read() + "\n")
            os.remove(fname)


def concat_reference_consens(data, sample):
    "concatenates consens bits into SAM for reference assemblies"
    
    with open(sample.files.consens, 'w') as outf:

        # parse fai file for writing headers
        fai = "{}.fai".format(data.paramsdict["reference_sequence"])
        fad = pd.read_csv(fai, sep="\t", names=["SN", "LN", "POS", "N1", "N2"])
        headers = [
            "@SQ     SN:{} LN:{}".format(i, j)
            for (i, j) in zip(fad["SN"], fad["LN"])
        ]
        outf.write("\n".join(headers) + "\n")

        # write sequences to SAM file
        combs1 = glob.glob(os.path.join(
            data.tmpdir,
            "{}_tmpcons.*".format(sample.name)))
        combs1.sort(key=lambda x: int(x.split(".")[-1]))

        # write to file with sample names imputed to line up with catg array
        counter = 0
        for fname in combs1:
            with open(fname) as infile:
                # impute catg ordered seqnames 
                data = infile.readlines()
                fdata = []
                for line in data:
                    name, chrom, rest = line.rsplit(":", 2)
                    fdata.append(
                        "{}_{}:{}:{}".format(name, counter, chrom, rest)
                        )
                    counter += 1
                outf.write("".join(fdata) + "\n")
            os.remove(fname)


def store_sample_stats(data, sample, statsdicts):
    "not parallel, store the sample objects stats"

    # record results
    xcounters = {
        "nconsens": 0,
        "heteros": 0,
        "nsites": 0,
    }
    xfilters = {
        "depth": 0,
        "maxh": 0,
        "maxn": 0,
    }

    # merge finished consens stats
    for counters, filters in statsdicts:
        # sum individual counters
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]

    # set Sample stats_dfs values
    if int(xcounters['nsites']):
        prop = int(xcounters["heteros"]) / float(xcounters['nsites'])
    else:
        prop = 0

    # store stats attributes to the sample
    sample.stats_dfs.s5.nsites = int(xcounters["nsites"])
    sample.stats_dfs.s5.nhetero = int(xcounters["heteros"])
    sample.stats_dfs.s5.filtered_by_depth = xfilters['depth']
    sample.stats_dfs.s5.filtered_by_maxH = xfilters['maxh']
    sample.stats_dfs.s5.filtered_by_maxN = xfilters['maxn']
    sample.stats_dfs.s5.reads_consens = int(xcounters["nconsens"])
    sample.stats_dfs.s5.clusters_total = sample.stats_dfs.s3.clusters_total
    sample.stats_dfs.s5.heterozygosity = float(prop)

    # set the Sample stats summary value
    sample.stats.reads_consens = int(xcounters["nconsens"])

    # save state to Sample if successful
    if sample.stats.reads_consens:
        sample.stats.state = 5
    else:
        print("No clusters passed filtering in Sample: {}".format(sample.name))


def make_chunks(data, sample, ncpus):
    "split job into bits and pass to the client"

    # counter for split job submission
    num = 0

    # set optim size for chunks in N clusters. The first few chunks take longer
    # because they contain larger clusters, so we create 4X as many chunks as
    # processors so that they are split more evenly.
    optim = int(
        (sample.stats.clusters_total // ncpus) + \
        (sample.stats.clusters_total % ncpus))

    # open to clusters
    with gzip.open(sample.files.clusters, 'rb') as clusters:
        # create iterator to sample 2 lines at a time
        pairdealer = izip(*[iter(clusters)] * 2)

        # Use iterator to sample til end of cluster
        done = 0
        while not done:
            # grab optim clusters and write to file.
            done, chunk = clustdealer(pairdealer, optim)
            chunk = [i.decode() for i in chunk]

            # make file handle
            chunkhandle = os.path.join(
                data.tmpdir,
                "{}.chunk.{}.{}".format(sample.name, optim, num * optim))

            # write to file
            if chunk:
                with open(chunkhandle, 'wt') as outchunk:
                    outchunk.write("//\n//\n".join(chunk) + "//\n//\n")
                num += 1


def process_chunks(data, sample, chunk, isref):
    "make consensus base and allele calls"

    # temporarily store the mean estimates to Assembly
    este = data.stats.error_est.mean()
    esth = data.stats.hetero_est.mean()

    # get index number from tmp file name
    tmpnum = int(chunk.split(".")[-1])
    optim = int(chunk.split(".")[-2])

    # prepare data for reading and use global maxfraglen
    clusters = open(chunk, 'rb')
    pairdealer = izip(*[iter(clusters)] * 2)
    maxlen = data._hackersonly["max_fragment_length"]

    # local copies to use to fill the arrays
    catarr = np.zeros((optim, maxlen, 4), dtype=np.uint32)
    nallel = np.zeros((optim, ), dtype=np.uint8)
    refarr = np.zeros((optim, 3), dtype=np.int64)

    # if reference-mapped then parse the fai to get index number of chroms
    if isref:
        fai = pd.read_csv(
            data.paramsdict["reference_sequence"] + ".fai",
            names=['scaffold', 'size', 'sumsize', 'a', 'b'],
            sep="\t")
        faidict = {j: i for i, j in enumerate(fai.scaffold)}

    # store data for stats counters.
    counters = {
        "name": tmpnum,
        "heteros": 0,
        "nsites": 0,
        "nconsens": 0,
    }

    # store data for what got filtered
    filters = {
        "depth": 0,
        "maxh": 0,
        "maxn": 0,
    }

    # store data for writing
    storeseq = {}

    ## set max limits
    if 'pair' in data.paramsdict["datatype"]:
        maxhet = sum(data.paramsdict["max_Hs_consens"])
        maxn = sum(data.paramsdict["max_Ns_consens"])
    else:
        maxhet = data.paramsdict["max_Hs_consens"][0]
        maxn = data.paramsdict["max_Ns_consens"][0]

    # load the refmap dictionary if refmapping
    done = 0
    while not done:
        done, chunk = clustdealer(pairdealer, 1)
        if chunk:            
            # get names and seqs
            piece = chunk[0].decode().strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]

            # pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]

            # IF this is a reference mapped read store the chrom and pos info
            # -1 defaults to indicating an anonymous locus, since we are using
            # the faidict as 0 indexed. If chrompos fails it defaults to -1
            # this is used for the 'denovo + reference' method.
            ref_position = (-1, 0, 0)
            if isref:
                # parse position from name string
                rname = names[0].rsplit(";", 2)[0]
                chrom, posish = rname.rsplit(":")
                pos0, pos1 = posish.split("-")
                
                # pull idx from .fai reference dict
                chromint = faidict[chrom] + 1
                ref_position = (int(chromint), int(pos0), int(pos1))

            # apply filters and fill arrays
            if nfilter1(data, reps):

                # get stacks of base counts
                sseqs = [list(seq) for seq in seqs]
                arrayed = np.concatenate(
                    [[seq] * rep for seq, rep in zip(sseqs, reps)]
                ).astype(bytes)
                arrayed = arrayed[:, :maxlen]
                
                # get consens call for each site, applies paralog-x-site filter
                consens = base_caller(
                    arrayed, 
                    data.paramsdict["mindepth_majrule"], 
                    data.paramsdict["mindepth_statistical"],
                    esth, 
                    este,
                )

                # apply a filter to remove low coverage sites/Ns that
                # are likely sequence repeat errors. This is only applied to
                # clusters that already passed the read-depth filter (1)
                # and is not applied to ref aligned data which is expected to
                # have internal spacers that we leave in place here. 
                if b"N" in consens:
                    if not isref:
                        try:
                            consens, arrayed = removerepeats(consens, arrayed)

                        except ValueError:
                            ip.logger.info(
                                "Caught bad chunk w/ all Ns. Skip it.")
                            continue

                # get hetero sites
                hidx = [
                    i for (i, j) in enumerate(consens) if j.decode() 
                    in list("RKSYWM")]
                nheteros = len(hidx)

                # filter for max number of hetero sites
                if nfilter2(nheteros, maxhet):
                    # filter for maxN, & minlen
                    if nfilter3(consens, maxn):
                        # counter right now
                        current = counters["nconsens"]
                        # get N alleles and get lower case in consens
                        consens, nhaps = nfilter4(consens, hidx, arrayed)
                        # store the number of alleles observed
                        nallel[current] = nhaps
                        # store a reduced array with only CATG
                        catg = np.array(
                            [np.sum(arrayed == i, axis=0)
                                for i in [b'C', b'A', b'T', b'G']],
                            dtype='uint32').T
                        catarr[current, :catg.shape[0], :] = catg
                        refarr[current] = ref_position

                        # store the seqdata for tmpchunk
                        # storeseq[counters["name"]] = b"".join(list(consens))
                        storeseq[current] = b"".join(list(consens))
                        counters["name"] += 1
                        counters["nconsens"] += 1
                        counters["heteros"] += nheteros
                    else:
                        filters['maxn'] += 1
                else:
                    filters['maxh'] += 1
            else:
                filters['depth'] += 1

    # close infile io
    clusters.close()

    # find last empty size
    end = np.where(np.all(refarr == 0, axis=1))[0]
    if np.any(end):
        end = end.min()
    else:
        end = refarr.shape[0]

    # write final consens string chunk
    consenshandle = os.path.join(
        data.tmpdir,
        "{}_tmpcons.{}.{}".format(sample.name, end, tmpnum))

    # write chunk
    if storeseq:
        with open(consenshandle, 'wt') as outfile:
            if not isref:
                outfile.write(
                    "\n".join([">" + sample.name + "_" + str(key) + \
                    "\n" + storeseq[key].decode() for key in storeseq]))
            else:
                outfile.write(
                    "\n".join(
                        ["{}:{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                        .format(
                            sample.name,
                            refarr[i][0],
                            refarr[i][1],
                            refarr[i][2],
                            0,
                            refarr[i][0],
                            refarr[i][1],
                            0,
                            "*",
                            "*",
                            0,
                            refarr[i][2] - refarr[i][1],
                            storeseq[i].decode(),
                            "*",
                        ) for i in storeseq.keys()]
                        ))

    # store reduced size arrays with indexes matching to keep indexes
    tmp5 = consenshandle.replace("_tmpcons.", "_tmpcats.")
    with h5py.File(tmp5, 'w') as io5:
        io5.create_dataset(name="cats", data=catarr[:end])
        io5.create_dataset(name="alls", data=nallel[:end])
        io5.create_dataset(name="chroms", data=refarr[:end])
    del catarr
    del nallel
    del refarr

    # return stats
    counters['nsites'] = sum([len(i) for i in storeseq.values()])
    del storeseq
    return counters, filters


def apply_filters_and_fill_arrs():
    pass


def base_caller(arrayed, mindepth_majrule, mindepth_statistical, estH, estE):
    "call all sites in a locus array. Can't be jit'd yet b/c scipy"

    # an array to fill with consensus site calls
    cons = np.zeros(arrayed.shape[1], dtype=np.uint8)
    cons.fill(78)
    arr = arrayed.view(np.uint8)

    # iterate over columns
    for col in range(arr.shape[1]):
        # the site of focus
        carr = arr[:, col]

        # if site is all dash then fill it
        if np.all(carr == 45):
            cons[col] = 45
            
        # else mask all N and - sites for base call
        else:
            mask = carr == 45
            mask += carr == 78
            marr = carr[~mask]
            
            # skip if only empties
            if not marr.shape[0]:
                cons[col] = 78
                
            # skip if not variable
            elif np.all(marr == marr[0]):
                cons[col] = marr[0]

            # estimate variable site call
            else:
                # get allele freqs (first-most, second, third = p, q, r)
                counts = np.bincount(marr)

                pbase = np.argmax(counts)
                nump = counts[pbase]
                counts[pbase] = 0

                qbase = np.argmax(counts)
                numq = counts[qbase]
                counts[qbase] = 0

                ## based on biallelic depth
                bidepth = nump + numq
                if bidepth < mindepth_majrule:
                    cons[col] = 78

                else:
                    # if depth is too high, reduce to sampled int
                    if bidepth > 500:
                        base1 = int(500 * (nump / float(bidepth)))
                        base2 = int(500 * (numq / float(bidepth)))
                    else:
                        base1 = nump
                        base2 = numq

                    # make statistical base call
                    if bidepth >= mindepth_statistical:
                        ishet, prob = get_binom(base1, base2, estE, estH)
                        if prob < 0.95:
                            cons[col] = 78
                        else:
                            if ishet:
                                cons[col] = TRANS[(pbase, qbase)]
                            else:
                                cons[col] = pbase

                    # make majrule base call
                    else:
                        if nump == numq:
                            cons[col] = TRANS[(pbase, qbase)]
                        else:
                            cons[col] = pbase
    return cons.view("S1")



TRANS = {
    (71, 65): 82,
    (71, 84): 75,
    (71, 67): 83,
    (84, 67): 89,
    (84, 65): 87,
    (67, 65): 77,
    (65, 67): 77,
    (65, 84): 87,
    (67, 84): 89,
    (67, 71): 83,
    (84, 71): 75,
    (65, 71): 82,
}


def get_binom(base1, base2, estE, estH):
    "return probability of base call"
    prior_homo = (1. - estH) / 2.
    prior_hete = estH

    ## calculate probs
    bsum = base1 + base2
    hetprob = scipy.misc.comb(bsum, base1) / (2. ** (bsum))
    homoa = scipy.stats.binom.pmf(base2, bsum, estE)
    homob = scipy.stats.binom.pmf(base1, bsum, estE)

    ## calculate probs
    hetprob *= prior_hete
    homoa *= prior_homo
    homob *= prior_homo

    ## final
    probabilities = [homoa, homob, hetprob]
    bestprob = max(probabilities) / float(sum(probabilities))

    ## return
    if hetprob > homoa:
        return True, bestprob
    return False, bestprob



def removerepeats(consens, arrayed):
    """
    Checks for interior Ns in consensus seqs and removes those that are at
    low depth, here defined as less than 1/3 of the average depth. The prop 1/3
    is chosen so that mindepth=6 requires 2 base calls that are not in [N,-].

    Python3 notes:
    consens and arrayed are both in bytes in entry. Consens is converted to
    unicode for operations, and both are returned as bytes.
    """

    ## default trim no edges
    consens[consens == b"-"] = b"N"
    consens = b"".join(consens)

    ## split for pairs
    try:
        cons1, cons2 = consens.split(b"nnnn")
        split = consens.index(b"nnnn")
        arr1 = arrayed[:, :split]
        arr2 = arrayed[:, split + 4:]
    except ValueError:
        cons1 = consens
        cons2 = ""
        arr1 = arrayed

    ## trim from left and right of cons1
    edges = [None, None]
    lcons = len(cons1)
    cons1 = cons1.lstrip(b"N")
    edges[0] = lcons - len(cons1)

    ## trim from right if nonzero
    lcons = len(cons1)
    cons1 = cons1.rstrip(b"N")
    if lcons - len(cons1):
        edges[1] = -1 * (lcons - len(cons1))

    ## trim same from arrayed
    arr1 = arr1[:, edges[0]:edges[1]]

    ## trim from left and right of cons2 if present
    if cons2:
        ## trim from left and right of cons1
        edges = [None, None]
        lcons = len(cons2)
        cons2 = cons2.lstrip(b"N")
        edges[0] = lcons - len(cons2)

        ## trim from right if nonzero
        lcons = len(cons2)
        cons2 = cons2.rstrip(b"N")
        if lcons - len(cons2):
            edges[1] = -1 * (lcons - len(cons2))

        ## trim same from arrayed
        arr2 = arr2[:, edges[0]:edges[1]]

        ## reconstitute pairs
        consens = cons1 + b"nnnn" + cons2
        consens = np.array(list(consens), dtype=bytes)
        sep = np.array(arr1.shape[0] * [list(b"nnnn")])
        arrayed = np.hstack([arr1, sep, arr2])

    ## if single-end...
    else:
        consens = np.array(list(cons1), dtype=bytes)
        arrayed = arr1

    ## get column counts of Ns and -s
    ndepths = np.sum(arrayed == b'N', axis=0)
    idepths = np.sum(arrayed == b'-', axis=0)

    ## get proportion of bases that are N- at each site
    nons = ((ndepths + idepths) / float(arrayed.shape[0])) >= 0.75
    ## boolean of whether base was called N
    isn = consens == b"N"
    ## make ridx
    ridx = nons * isn

    ## apply filter
    consens = consens[~ridx]
    arrayed = arrayed[:, ~ridx]

    return consens, arrayed



# def newconsensus(data, sample, tmpchunk, optim):
#     """
#     new faster replacement to consensus
#     """
#     ## do reference map funcs?
#     isref = bool("reference" in data.paramsdict["assembly_method"])

#     ## temporarily store the mean estimates to Assembly
#     data._este = data.stats.error_est.mean()
#     data._esth = data.stats.hetero_est.mean()

#     ## get number relative to tmp file
#     tmpnum = int(tmpchunk.split(".")[-1])

#     ## prepare data for reading
#     clusters = open(tmpchunk, 'rb')
#     pairdealer = izip(*[iter(clusters)] * 2)
#     maxlen = data._hackersonly["max_fragment_length"]

#     ## write to tmp cons to file to be combined later
#     consenshandle = os.path.join(
#         data.dirs.consens, sample.name + "_tmpcons." + str(tmpnum))
#     tmp5 = consenshandle.replace("_tmpcons.", "_tmpcats.")
#     with h5py.File(tmp5, 'w') as io5:
#         io5.create_dataset("cats", (optim, maxlen, 4), dtype=np.uint32)
#         io5.create_dataset("alls", (optim, ), dtype=np.uint8)
#         io5.create_dataset("chroms", (optim, 3), dtype=np.int64)

#         ## local copies to use to fill the arrays
#         catarr = io5["cats"][:]
#         nallel = io5["alls"][:]
#         refarr = io5["chroms"][:]

#     ## if reference-mapped then parse the fai to get index number of chroms
#     if isref:
#         fai = pd.read_csv(
#             data.paramsdict["reference_sequence"] + ".fai",
#             names=['scaffold', 'size', 'sumsize', 'a', 'b'],
#             sep="\t")
#         faidict = {j: i for i, j in enumerate(fai.scaffold)}

#     ## store data for stats counters
#     counters = {"name": tmpnum,
#                 "heteros": 0,
#                 "nsites": 0,
#                 "nconsens": 0}

#     ## store data for what got filtered
#     filters = {"depth": 0,
#                "maxh": 0,
#                "maxn": 0}

#     ## store data for writing
#     storeseq = {}

#     ## set max limits
#     if 'pair' in data.paramsdict["datatype"]:
#         maxhet = sum(data.paramsdict["max_Hs_consens"])
#         maxn = sum(data.paramsdict["max_Ns_consens"])
#     else:
#         maxhet = data.paramsdict["max_Hs_consens"][0]
#         maxn = data.paramsdict["max_Ns_consens"][0]

#     ## load the refmap dictionary if refmapping
#     done = 0
#     while not done:
#         try:
#             done, chunk = clustdealer(pairdealer, 1)
#         except IndexError:
#             raise IPyradError("clustfile formatting error in {}".format(chunk))

#         if chunk:
#             ## get names and seqs
#             piece = chunk[0].decode().strip().split("\n")
#             names = piece[0::2]
#             seqs = piece[1::2]

#             ## pull replicate read info from seqs
#             reps = [int(sname.split(";")[-2][5:]) for sname in names]

#             ## IF this is a reference mapped read store the chrom and pos info
#             ## -1 defaults to indicating an anonymous locus, since we are using
#             ## the faidict as 0 indexed. If chrompos fails it defaults to -1
#             ref_position = (-1, 0, 0)
#             if isref:
#                 try:
#                     ## parse position from name string
#                     name, _, _ = names[0].rsplit(";", 2)
#                     chrom, pos0, pos1 = name.rsplit(":", 2)

#                     ## pull idx from .fai reference dict
#                     chromint = faidict[chrom] + 1
#                     ref_position = (int(chromint), int(pos0), int(pos1))

#                 except Exception as inst:
#                     ip.logger.debug(
#                         "Reference sequence chrom/pos failed for {}"
#                         .format(names[0]))
#                     ip.logger.debug(inst)
                    
#             ## apply read depth filter
#             if nfilter1(data, reps):

#                 ## get stacks of base counts
#                 sseqs = [list(seq) for seq in seqs]
#                 arrayed = np.concatenate(
#                     [[seq] * rep for seq, rep in zip(sseqs, reps)]
#                 ).astype(bytes)
#                 arrayed = arrayed[:, :maxlen]
                
#                 # get consens call for each site, applies paralog-x-site filter
#                 # consens = np.apply_along_axis(basecall, 0, arrayed, data)
#                 consens = basecaller(
#                     arrayed, 
#                     data.paramsdict["mindepth_majrule"], 
#                     data.paramsdict["mindepth_statistical"],
#                     data._esth, 
#                     data._este,
#                 )

#                 ## apply a filter to remove low coverage sites/Ns that
#                 ## are likely sequence repeat errors. This is only applied to
#                 ## clusters that already passed the read-depth filter (1)
#                 if "N" in consens:
#                     try:
#                         consens, arrayed = removerepeats(consens, arrayed)

#                     except ValueError:
#                         ip.logger.info("Caught bad chunk w/ all Ns. Skip it.")
#                         continue

#                 ## get hetero sites
#                 hidx = [i for (i, j) in enumerate(consens)
#                         if j in list("RKSYWM")]
#                 nheteros = len(hidx)

#                 ## filter for max number of hetero sites
#                 if nfilter2(nheteros, maxhet):
#                     ## filter for maxN, & minlen
#                     if nfilter3(consens, maxn):
#                         ## counter right now
#                         current = counters["nconsens"]
#                         ## get N alleles and get lower case in consens
#                         consens, nhaps = nfilter4(consens, hidx, arrayed)
#                         ## store the number of alleles observed
#                         nallel[current] = nhaps

#                         ## store a reduced array with only CATG
#                         catg = np.array(
#                             [np.sum(arrayed == i, axis=0)
#                                 for i in list("CATG")],
#                             dtype='uint32').T
#                         catarr[current, :catg.shape[0], :] = catg
#                         refarr[current] = ref_position

#                         ## store the seqdata for tmpchunk
#                         storeseq[counters["name"]] = b"".join(list(consens))
#                         counters["name"] += 1
#                         counters["nconsens"] += 1
#                         counters["heteros"] += nheteros
#                     else:
#                         #ip.logger.debug("@haplo")
#                         filters['maxn'] += 1
#                 else:
#                     #ip.logger.debug("@hetero")
#                     filters['maxh'] += 1
#             else:
#                 #ip.logger.debug("@depth")
#                 filters['depth'] += 1

#     ## close infile io
#     clusters.close()

#     ## write final consens string chunk
#     if storeseq:
#         with open(consenshandle, 'wt') as outfile:
#             outfile.write(
#                 "\n".join([">" + sample.name + "_" + str(key) + \
#                 "\n" + storeseq[key].decode() for key in storeseq]))

#     ## write to h5 array, this can be a bit slow on big data sets and is not
#     ## currently convered by progressbar movement.
#     with h5py.File(tmp5, 'a') as io5:
#         io5["cats"][:] = catarr
#         io5["alls"][:] = nallel
#         io5["chroms"][:] = refarr
#     del catarr
#     del nallel
#     del refarr

#     ## return stats
#     counters['nsites'] = sum([len(i) for i in storeseq.values()])
#     return counters, filters



# def basecaller(arrayed, mindepth_majrule, mindepth_statistical, estH, estE):
#     "call all sites in a locus array. Can be jit'd..."

#     ## an array to fill with consensus site calls
#     cons = np.zeros(arrayed.shape[1], dtype=np.uint8)
#     cons.fill(78)
#     arr = arrayed.view(np.uint8)

#     ## iterate over columns
#     for col in range(arr.shape[1]):
#         ## the site of focus
#         carr = arr[:, col]

#         ## make mask of N and - sites
#         mask = carr == 45
#         mask += carr == 78
#         marr = carr[~mask]

#         ## skip if only empties (e.g., N-)
#         if not marr.shape[0]:
#             cons[col] = 78

#         ## skip if not variable
#         elif np.all(marr == marr[0]):
#             cons[col] = marr[0]

#         ## estimate variable site call
#         else:
#             ## get allele freqs (first-most, second, third = p, q, r)
#             counts = np.bincount(marr)

#             pbase = np.argmax(counts)
#             nump = counts[pbase]
#             counts[pbase] = 0

#             qbase = np.argmax(counts)
#             numq = counts[qbase]
#             counts[qbase] = 0

#             #rbase = np.argmax(counts)
#             #numr = counts[rbase]          # not used

#             ## based on biallelic depth
#             bidepth = nump + numq
#             if bidepth < mindepth_majrule:
#                 cons[col] = 78

#             else:
#                 ## if depth is too high, reduce to sampled int
#                 if bidepth > 500:
#                     base1 = int(500 * (nump / float(bidepth)))
#                     base2 = int(500 * (numq / float(bidepth)))
#                 else:
#                     base1 = nump
#                     base2 = numq

#                 ## make statistical base call
#                 if bidepth >= mindepth_statistical:
#                     ishet, prob = get_binom(base1, base2, estE, estH)
#                     #ip.logger.info("ishet, prob, b1, b2: %s %s %s %s", ishet, prob, base1, base2)
#                     if prob < 0.95:
#                         cons[col] = 78
#                     else:
#                         if ishet:
#                             cons[col] = TRANS[(pbase, qbase)]
#                         else:
#                             cons[col] = pbase

#                 ## make majrule base call
#                 else:  # if bidepth >= mindepth_majrule:
#                     if nump == numq:
#                         cons[col] = TRANS[(pbase, qbase)]
#                     else:
#                         cons[col] = pbase
#     return cons.view("S1")



def nfilter1(data, reps):
    "applies read depths filter"
    if sum(reps) >= data.paramsdict["mindepth_majrule"] and \
        sum(reps) <= data.paramsdict["maxdepth"]:
        return 1
    return 0



def nfilter2(nheteros, maxhet):
    "applies max heteros in a seq filter"
    if nheteros <= maxhet:
        return 1
    return 0



def nfilter3(consens, maxn):
    "applies filter for maxN and hard minlen (32)"
    # minimum length for clustering in vsearch
    if consens.size >= 32:
        if consens[consens == b"N"].size <= maxn:
            return 1
        return 0
    return 0



def nfilter4(consens, hidx, arrayed):
    "applies max haplotypes filter returns pass and consens"

    # if less than two Hs then there is only one allele
    if len(hidx) < 2:
        return consens, 1

    # store base calls for hetero sites
    harray = arrayed[:, hidx]

    # remove any reads that have N or - base calls at hetero sites
    # these cannot be used when calling alleles currently.
    harray = harray[~np.any(harray == b"-", axis=1)]
    harray = harray[~np.any(harray == b"N", axis=1)]

    # get counts of each allele (e.g., AT:2, CG:2)
    ccx = Counter([tuple(i) for i in harray])

    ## Two possibilities we would like to distinguish, but we can't. Therefore,
    ## we just throw away low depth third alleles that are within seq. error.
    ## 1) a third base came up as a sequencing error but is not a unique allele
    ## 2) a third or more unique allele is there but at low frequency

    ## remove low freq alleles if more than 2, since they may reflect
    ## sequencing errors at hetero sites, making a third allele, or a new
    ## allelic combination that is not real.
    if len(ccx) > 2:
        totdepth = harray.shape[0]
        cutoff = max(1, totdepth // 10)
        alleles = [i for i in ccx if ccx[i] > cutoff]
    else:
        alleles = ccx.keys()

    ## how many high depth alleles?
    nalleles = len(alleles)

    ## if 2 alleles then save the phase using lowercase coding
    if nalleles == 2:
        try:
            consens = storealleles(consens, hidx, alleles)
        except (IndexError, KeyError):
            ## the H sites do not form good alleles
            ip.logger.info("failed at phasing loc, skipping")
            ip.logger.info("""
    consens %s
    hidx %s
    alleles %s
                """, consens, hidx, alleles)
        return consens, nalleles

    ## just return the info for later filtering
    else:
        return consens, nalleles



def storealleles(consens, hidx, alleles):
    """ store phased allele data for diploids """
    ## find the first hetero site and choose the priority base
    ## example, if W: then priority base in A and not T. PRIORITY=(order: CATG)
    bigbase = PRIORITY[consens[hidx[0]]]

    ## find which allele has priority based on bigbase
    bigallele = [i for i in alleles if i[0] == bigbase][0]

    ## uplow other bases relative to this one and the priority list
    ## e.g., if there are two hetero sites (WY) and the two alleles are
    ## AT and TC, then since bigbase of (W) is A second hetero site should
    ## be stored as y, since the ordering is swapped in this case; the priority
    ## base (C versus T) is C, but C goes with the minor base at h site 1.
    #consens = list(consens)
    for hsite, pbase in zip(hidx[1:], bigallele[1:]):
        if PRIORITY[consens[hsite]] != pbase:
            consens[hsite] = consens[hsite].lower()

    ## return consens
    return consens
