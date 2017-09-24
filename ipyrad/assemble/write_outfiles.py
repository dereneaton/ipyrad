#!/usr/bin/env python2.7

""" 
Apply filters and write output files. The basic body plan of this
code is as follows:
  * Read in the final aligned clusters file
  * Make sizeable chunks of loci
  * Distribute these chunks to parallel filters
  * Combine output of parallel filters into final loci file
  * Write out the output in full vcf format
"""

# pylint: disable=W0142
# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=C0301
# pylint: disable=C0103

from __future__ import print_function

import pandas as pd
import numpy as np
import datetime
import shutil
import numba
import copy
import time
import glob
import gzip
import h5py
import re
import os
from collections import Counter, OrderedDict
from ipyrad import __version__
from util import *

try:
    import subprocess32 as sps
except ImportError:
    import subprocess as sps

import logging
LOGGER = logging.getLogger(__name__)


## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized. LOCI and VCF are default. Some others are created as
## dependencies of others.

OUTPUT_FORMATS = {'l': 'loci',
                  'p': 'phy',
                  's': 'snps',
                  'n': 'nex',
                  'k': 'struct',
                  'a': 'alleles',
                  'g': 'geno',
                  'G': "gphocs",
                  'u': 'usnps',
                  'v': 'vcf',
                  't': 'treemix',
                  'm': 'migrate-n'}
                  #'V': 'vcfFull',   ## currently hidden


def run(data, samples, force, ipyclient):
    """
    Check all samples requested have been clustered (state=6), make output
    directory, then create the requested outfiles. Excluded samples are already
    removed from samples.
    """

    ## prepare dirs
    data.dirs.outfiles = os.path.join(data.dirs.project, data.name+"_outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    ## make the snps/filters data base, fills the dups and inds filters
    ## and fills the splits locations
    data.database = os.path.join(data.dirs.outfiles, data.name+".hdf5")
    init_arrays(data)

    ## Apply filters to supercatg and superhdf5 with selected samples
    ## and fill the filters and edge arrays.
    filter_all_clusters(data, samples, ipyclient)

    ## Everything needed is in the now filled h5 database. Filters were applied
    ## with 'samples' taken into account. Now we create the loci file (default)
    ## output and build a stats file.
    data.outfiles.loci = os.path.join(data.dirs.outfiles, data.name+".loci")
    data.outfiles.alleles = os.path.join(data.dirs.outfiles, data.name+".alleles.loci")
    make_loci_and_stats(data, samples, ipyclient)

    ## OPTIONAL OUTPUTS:
    output_formats = data.paramsdict["output_formats"]

    ## held separate from *output_formats cuz it's big and parallelized
    if any([x in output_formats for x in ["v", "V"]]):
        full = "V" in output_formats
        try:
            make_vcf(data, samples, ipyclient, full=full)
        except IPyradWarningExit as inst:
            ## Something fsck vcf build. Sometimes this is simply a memory
            ## issue, so trap the exception and allow it to try building
            ## the other output formats.
            print("  Error building vcf. See ipyrad_log.txt for details.")
            LOGGER.error(inst)

    ## make other array-based formats, recalcs keeps and arrays
    make_outfiles(data, samples, output_formats, ipyclient)

    ## print friendly message
    shortpath = data.dirs.outfiles.replace(os.path.expanduser("~"), "~")
    print("{}Outfiles written to: {}\n".format(data._spacer, shortpath))



def make_stats(data, samples, samplecounts, locuscounts):
    """ write the output stats file and save to Assembly obj."""

    ## get meta info
    with h5py.File(data.clust_database, 'r') as io5:
        anames = io5["seqs"].attrs["samples"]
        nloci = io5["seqs"].shape[0]
        optim = io5["seqs"].attrs["chunksize"][0]

    ## open the out handle. This will have three data frames saved to it.
    ## locus_filtering, sample_coverages, and snp_distributions
    data.stats_files.s7 = os.path.join(data.dirs.outfiles,
                                       data.name+"_stats.txt")
    outstats = open(data.stats_files.s7, 'w')

    ########################################################################
    ## get stats for locus_filtering, use chunking.
    filters = np.zeros(6, dtype=int)
    passed = 0
    start = 0

    piscounts = Counter()
    varcounts = Counter()
    for i in range(200):
        piscounts[i] = 0
        varcounts[i] = 0


    applied = pd.Series([0]*8,
        name="applied_order",
        index=[
        "total_prefiltered_loci",
        "filtered_by_rm_duplicates",
        "filtered_by_max_indels",
        "filtered_by_max_snps",
        "filtered_by_max_shared_het",
        "filtered_by_min_sample",
        "filtered_by_max_alleles",
        "total_filtered_loci"])

    ## load the h5 database
    co5 = h5py.File(data.database, 'r')

    while start < nloci:
        hslice = [start, start+optim]
        ## load each array
        afilt = co5["filters"][hslice[0]:hslice[1], ]
        asnps = co5["snps"][hslice[0]:hslice[1], ]

        ## get subarray results from filter array
        # max_indels, max_snps, max_hets, min_samps, bad_edges, max_alleles
        filters += afilt.sum(axis=0)
        applied["filtered_by_rm_duplicates"] += afilt[:, 0].sum()
        mask = afilt[:, 0].astype(np.bool)
        applied["filtered_by_max_indels"] += afilt[~mask, 1].sum()
        mask = afilt[:, 0:2].sum(axis=1).astype(np.bool)
        applied["filtered_by_max_snps"] += afilt[~mask, 2].sum()
        mask = afilt[:, 0:3].sum(axis=1).astype(np.bool)
        applied["filtered_by_max_shared_het"] += afilt[~mask, 3].sum()
        mask = afilt[:, 0:4].sum(axis=1).astype(np.bool)
        applied["filtered_by_min_sample"] += afilt[~mask, 4].sum()
        mask = afilt[:, 0:5].sum(axis=1).astype(np.bool)
        applied["filtered_by_max_alleles"] += afilt[~mask, 5].sum()

        passed += np.sum(afilt.sum(axis=1) == 0)

        ## get filter to count snps for only passed loci
        ## should we filter by all vars, or just by pis? doing all var now.
        apply_filter = afilt.sum(axis=1).astype(np.bool)

        ## get snps counts
        snplocs = asnps[~apply_filter, :].sum(axis=1)
        varlocs = snplocs.sum(axis=1)
        varcounts.update(Counter(varlocs))
        #snpcounts.update(Counter(snplocs[:, 0]))
        piscounts.update(Counter(snplocs[:, 1]))

        ## increase counter to advance through h5 database
        start += optim

    ## record filtering of loci from total to final
    filtdat = pd.Series(np.concatenate([[nloci], filters, [passed]]),
        name="total_filters",
        index=[
        "total_prefiltered_loci",
        "filtered_by_rm_duplicates",
        "filtered_by_max_indels",
        "filtered_by_max_snps",
        "filtered_by_max_shared_het",
        "filtered_by_min_sample",
        "filtered_by_max_alleles",
        "total_filtered_loci"])


    retained = pd.Series([0]*8,
        name="retained_loci",
        index=[
        "total_prefiltered_loci",
        "filtered_by_rm_duplicates",
        "filtered_by_max_indels",
        "filtered_by_max_snps",
        "filtered_by_max_shared_het",
        "filtered_by_min_sample",
        "filtered_by_max_alleles",
        "total_filtered_loci"])
    retained["total_prefiltered_loci"] = nloci
    retained["filtered_by_rm_duplicates"] = nloci - applied["filtered_by_rm_duplicates"]
    retained["filtered_by_max_indels"] = retained["filtered_by_rm_duplicates"] - applied["filtered_by_max_indels"]
    retained["filtered_by_max_snps"] = retained["filtered_by_max_indels"] - applied["filtered_by_max_snps"]
    retained["filtered_by_max_shared_het"] = retained["filtered_by_max_snps"] - applied["filtered_by_max_shared_het"]
    retained["filtered_by_min_sample"] = retained["filtered_by_max_shared_het"] - applied["filtered_by_min_sample"]
    retained["filtered_by_max_alleles"] = retained["filtered_by_min_sample"] - applied["filtered_by_max_alleles"]
    retained["total_filtered_loci"] = passed


    print("\n\n## The number of loci caught by each filter."+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_filters\n",
          file=outstats)
    data.stats_dfs.s7_filters = pd.DataFrame([filtdat, applied, retained]).T
    data.stats_dfs.s7_filters.to_string(buf=outstats)


    ########################################################################
    ## make dataframe of sample_coverages
    ## samplecounts is len of anames from db. Save only samples in samples.
    #print(samplecounts)
    #samples = [i.name for i in samples]
    ## get sample names in the order of anames
    #sids = [list(anames).index(i) for i in samples]
    #covdict = {name: val for name, val in zip(np.array(samples)[sidx], samplecounts)}
    #covdict = {name: val for name, val in zip(samples, samplecounts[sidx])}
    covdict = pd.Series(samplecounts, name="sample_coverage", index=anames)
    covdict = covdict[covdict != 0]
    print("\n\n\n## The number of loci recovered for each Sample."+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_samples\n",
          file=outstats)
    data.stats_dfs.s7_samples = pd.DataFrame(covdict)
    data.stats_dfs.s7_samples.to_string(buf=outstats)


    ########################################################################
    ## get stats for locus coverage
    lrange = range(1, len(samples)+1)
    locdat = pd.Series(locuscounts, name="locus_coverage", index=lrange)
    start = data.paramsdict["min_samples_locus"]-1
    locsums = pd.Series({i: np.sum(locdat.values[start:i]) for i in lrange},
                        name="sum_coverage", index=lrange)
    print("\n\n\n## The number of loci for which N taxa have data."+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_loci\n",
          file=outstats)
    data.stats_dfs.s7_loci = pd.concat([locdat, locsums], axis=1)
    data.stats_dfs.s7_loci.to_string(buf=outstats)


    #########################################################################
    ## get stats for SNP_distribution
    try:
        smax = max([i+1 for i in varcounts if varcounts[i]])
    except Exception as inst:
        raise IPyradWarningExit("""
    Exception: empty varcounts array. This could be because no samples 
    passed filtering, or it could be because you have overzealous filtering.
    Check the values for `trim_loci` and make sure you are not trimming the
    edge too far
    """)

    vardat = pd.Series(varcounts, name="var", index=range(smax)).fillna(0)
    sumd = {}
    for i in range(smax):
        sumd[i] = np.sum([i*vardat.values[i] for i in range(i+1)])
    varsums = pd.Series(sumd, name="sum_var", index=range(smax))

    pisdat = pd.Series(piscounts, name="pis", index=range(smax)).fillna(0)
    sumd = {}
    for i in range(smax):
        sumd[i] = np.sum([i*pisdat.values[i] for i in range(i+1)])
    pissums = pd.Series(sumd, name="sum_pis", index=range(smax))

    print("\n\n\n## The distribution of SNPs (var and pis) per locus."+\
          "\n## var = Number of loci with n variable sites (pis + autapomorphies)"+\
          "\n## pis = Number of loci with n parsimony informative site (minor allele in >1 sample)"+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_snps\n",
          file=outstats)
    data.stats_dfs.s7_snps = pd.concat([vardat, varsums, pisdat, pissums],
                                        axis=1)
    data.stats_dfs.s7_snps.to_string(buf=outstats)

    ##########################################################################
    ## print the stats summary (-r summary) with final sample loci data.
    fullstat = data.stats
    fullstat['state'] = 7
    fullstat["loci_in_assembly"] = data.stats_dfs.s7_samples

    print("\n\n\n## Final Sample stats summary\n", file=outstats)
    fullstat.to_string(buf=outstats)

    ## close it
    outstats.close()
    co5.close()



def select_samples(dbsamples, samples):
    """
    Get the row index of samples that are included. If samples are in the
    'excluded' they were already filtered out of 'samples' during _get_samples.
    """
    ## get index from dbsamples
    samples = [i.name for i in samples]
    sidx = [list(dbsamples).index(i) for i in samples]
    sidx.sort()
    return sidx



def filter_all_clusters(data, samples, ipyclient):
    """
    Open the clust_database HDF5 array with seqs, catg, and filter data. 
    Fill the remaining filters.
    """
    ## create loadbalanced ipyclient
    lbview = ipyclient.load_balanced_view()

    ## get chunk size from the HD5 array and close
    with h5py.File(data.clust_database, 'r') as io5:
        ## the size of chunks for reading/writing
        optim = io5["seqs"].attrs["chunksize"][0]
        ## the samples in the database in their locus order
        dbsamples = io5["seqs"].attrs["samples"]
        ## the total number of loci
        nloci = io5["seqs"].shape[0]

    ## make a tmp directory for saving chunked arrays to
    chunkdir = os.path.join(data.dirs.outfiles, data.name+"_tmpchunks")
    if not os.path.exists(chunkdir):
        os.mkdir(chunkdir)

    ## get the indices of the samples that we are going to include
    sidx = select_samples(dbsamples, samples)
    ## do the same for the populations samples
    if data.populations:
        data._populations = {}
        for pop in data.populations:
            _samps = [data.samples[i] for i in data.populations[pop][1]]
            data._populations[pop] = (data.populations[pop][0],
                                      select_samples(dbsamples, _samps))
                                                
    LOGGER.info("samples %s \n, dbsamples %s \n, sidx %s \n",
                samples, dbsamples, sidx)

    ## Put inside a try statement so we can delete tmpchunks
    try:
        ## load a list of args to send to Engines. Each arg contains the index
        ## to sample optim loci from catg, seqs, filters &or edges, which will
        ## be loaded on the remote Engine.

        ## create job queue
        start = time.time()
        printstr = " filtering loci        | {} | s7 |"
        fasyncs = {}
        submitted = 0
        while submitted < nloci:
            hslice = np.array([submitted, submitted+optim])
            fasyncs[hslice[0]] = lbview.apply(filter_stacks, *(data, sidx, hslice))
            submitted += optim

        ## run filter_stacks on all chunks
        while 1:
            readies = [i.ready() for i in fasyncs.values()]
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(len(readies), sum(readies),
                printstr.format(elapsed), spacer=data._spacer)
            time.sleep(0.1)
            if sum(readies) == len(readies):
                print("")
                break

        ## raise error if any jobs failed
        for async in fasyncs:
            if not fasyncs[async].successful():
                LOGGER.error("error in filter_stacks on chunk %s: %s",
                             async, fasyncs[async].exception())
                raise IPyradWarningExit("error in filter_stacks on chunk {}: {}"\
                             .format(async, fasyncs[async].exception()))
        ipyclient.purge_everything()

        ## get all the saved tmp arrays for each slice
        tmpsnp = glob.glob(os.path.join(chunkdir, "snpf.*.npy"))
        tmphet = glob.glob(os.path.join(chunkdir, "hetf.*.npy"))
        tmpmin = glob.glob(os.path.join(chunkdir, "minf.*.npy"))
        tmpedg = glob.glob(os.path.join(chunkdir, "edgf.*.npy"))
        tmppld = glob.glob(os.path.join(chunkdir, "pldf.*.npy"))
        tmpind = glob.glob(os.path.join(chunkdir, "indf.*.npy"))

        ## sort array files within each group
        arrdict = OrderedDict([('ind', tmpind),
                               ('snp', tmpsnp), ('het', tmphet),
                               ('min', tmpmin), ('edg', tmpedg),
                               ('pld', tmppld)])
        for arrglob in arrdict.values():
            arrglob.sort(key=lambda x: int(x.rsplit(".")[-2]))

        ## re-load the full filter array who's order is
        ## ["duplicates", "max_indels", "max_snps", "max_hets", "min_samps", "max_alleles"]
        io5 = h5py.File(data.database, 'r+')
        superfilter = np.zeros(io5["filters"].shape, io5["filters"].dtype)

        ## iterate across filter types (dups is already filled)
        ## we have [4,4] b/c minf and edgf both write to minf
        for fidx, ftype in zip([1, 2, 3, 4, 4, 5], arrdict.keys()):
            ## fill in the edgefilters
            for ffile in arrdict[ftype]:
                ## grab a file and get it's slice
                hslice = int(ffile.split(".")[-2])
                ## load in the array
                arr = np.load(ffile)
                ## store slice into full array (we use += here because the minf
                ## and edgf arrays both write to the same filter).
                superfilter[hslice:hslice+optim, fidx] += arr

        ## store to DB
        io5["filters"][:] += superfilter
        del arr, superfilter

        ## store the other arrayed values (edges, snps)
        edgarrs = glob.glob(os.path.join(chunkdir, "edgearr.*.npy"))
        snparrs = glob.glob(os.path.join(chunkdir, "snpsarr.*.npy"))

        ## sort array files within each group
        arrdict = OrderedDict([('edges', edgarrs), ('snps', snparrs)])
        for arrglob in arrdict.values():
            arrglob.sort(key=lambda x: int(x.rsplit(".")[-2]))

        ## fill the edge array, splits are already in there.
        superedge = np.zeros(io5['edges'].shape, io5['edges'].dtype)
        for ffile in arrdict['edges']:
            ## grab a file and get it's slice
            hslice = int(ffile.split(".")[-2])
            ## load in the array w/ shape (hslice, 5)
            arr = np.load(ffile)
            ## store slice into full array
            superedge[hslice:hslice+optim, :] = arr
        io5["edges"][:, :] = superedge
        del arr, superedge

        ## fill the snps array. shape= (nloci, maxlen, 2)
        supersnps = np.zeros(io5['snps'].shape, io5['snps'].dtype)
        for ffile in arrdict['snps']:
            ## grab a file and get it's slice
            hslice = int(ffile.split(".")[-2])
            ## load in the array w/ shape (hslice, maxlen, 2)
            arr = np.load(ffile)
            ## store slice into full array
            LOGGER.info("shapes, %s %s", supersnps.shape, arr.shape)
            supersnps[hslice:hslice+optim, :, :] = arr
        io5["snps"][:] = supersnps
        del arr
        io5.close()

    finally:
        ## clean up the tmp files/dirs even if we failed.
        try:
            LOGGER.info("finished filtering")
            shutil.rmtree(chunkdir)
        except (IOError, OSError):
            pass



def padnames(names):
    """ pads names for loci output """

    ## get longest name
    longname_len = max(len(i) for i in names)
    ## Padding distance between name and seq.
    padding = 5
    ## add pad to names
    pnames = [name + " " * (longname_len - len(name)+ padding) \
              for name in names]
    snppad = "//" + " " * (longname_len - 2 + padding)
    return np.array(pnames), snppad



## incorporating samples...
def make_loci_and_stats(data, samples, ipyclient):
    """
    Makes the .loci file from h5 data base. Iterates by optim loci at a
    time and write to file. Also makes alleles file if requested.
    """
    ## start vcf progress bar
    start = time.time()
    printstr = " building loci/stats   | {} | s7 |"
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, printstr.format(elapsed), spacer=data._spacer)

    ## get some db info
    with h5py.File(data.clust_database, 'r') as io5:
        ## will iterate optim loci at a time
        optim = io5["seqs"].attrs["chunksize"][0]
        nloci = io5["seqs"].shape[0]
        anames = io5["seqs"].attrs["samples"]

    ## get name and snp padding
    pnames, snppad = padnames(anames)
    snames = [i.name for i in samples]
    smask = np.array([i not in snames for i in anames])

    ## keep track of how many loci from each sample pass all filters
    samplecov = np.zeros(len(anames), dtype=np.int32)

    ## set initial value to zero for all values above min_samples_locus
    #for cov in range(data.paramsdict["min_samples_locus"], len(anames)+1):
    locuscov = Counter()
    for cov in range(len(anames)+1):
        locuscov[cov] = 0

    ## client for sending jobs to parallel engines
    lbview = ipyclient.load_balanced_view()

    ## send jobs in chunks
    loci_asyncs = {}
    for istart in xrange(0, nloci, optim):
        args = [data, optim, pnames, snppad, smask, istart, samplecov, locuscov, 1]
        loci_asyncs[istart] = lbview.apply(locichunk, args)

    while 1:
        done = [i.ready() for i in loci_asyncs.values()]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(done), sum(done), printstr.format(elapsed), spacer=data._spacer)
        time.sleep(0.1)
        if len(done) == sum(done):
            print("")
            break

    ## check for errors
    for job in loci_asyncs:
        if loci_asyncs[job].ready() and not loci_asyncs[job].successful():
            LOGGER.error("error in building loci [%s]: %s",
                         job, loci_asyncs[job].exception())
            raise IPyradWarningExit(loci_asyncs[job].exception())

    ## concat and cleanup
    results = [i.get() for i in loci_asyncs.values()]
    ## update dictionaries
    for chunk in results:
        samplecov += chunk[0]
        locuscov.update(chunk[1])

    ## get all chunk files
    tmploci = glob.glob(data.outfiles.loci+".[0-9]*")
    ## sort by start value
    tmploci.sort(key=lambda x: int(x.split(".")[-1]))

    ## write tmpchunks to locus file
    locifile = open(data.outfiles.loci, 'w')
    for tmploc in tmploci:
        with open(tmploc, 'r') as inloc:
            locdat = inloc.read()
            locifile.write(locdat)
            os.remove(tmploc)
    locifile.close()

    ## make stats file from data
    make_stats(data, samples, samplecov, locuscov)

    ## repeat for alleles output
    if "a" in data.paramsdict["output_formats"]:

        loci_asyncs = {}
        for istart in xrange(0, nloci, optim):
            args = [data, optim, pnames, snppad, smask, istart, samplecov, locuscov, 0]
            loci_asyncs[istart] = lbview.apply(locichunk, args)

        while 1:
            done = [i.ready() for i in loci_asyncs.values()]
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(len(done), sum(done),
                " building alleles      | {} | s7 |".format(elapsed), 
                spacer=data._spacer)
            time.sleep(0.1)
            if len(done) == sum(done):
                print("")
                break

        ## check for errors
        for job in loci_asyncs:
            if loci_asyncs[job].ready() and not loci_asyncs[job].successful():
                LOGGER.error("error in building alleles [%s]: %s",
                             job, loci_asyncs[job].exception())
                raise IPyradWarningExit(loci_asyncs[job].exception())

        ## concat and cleanup
        #results = [i.get() for i in loci_asyncs.values()]

        ## get all chunk files
        tmploci = glob.glob(data.outfiles.loci+".[0-9]*")
        ## sort by start value
        tmploci.sort(key=lambda x: int(x.split(".")[-1]))

        ## write tmpchunks to locus file
        locifile = open(data.outfiles.alleles, 'w')
        for tmploc in tmploci:
            with open(tmploc, 'r') as inloc:
                locdat = inloc.read()
                inalleles = get_alleles(locdat)
                locifile.write(inalleles)
                os.remove(tmploc)
        locifile.close()



def get_alleles(locdat):
    locs = []
    for loc in locdat.split("|\n"):
        lines = loc.split("\n")
        inloc = []
        for line in lines[:-1]:
            try:
                #LOGGER.info("line %s...", line)
                firstspace = line.index(" ")
                lastspace = line.rindex(" ")
                spacer = lastspace - firstspace + 1
                name, seq = line.split()
                seq1, seq2 = splitalleles(seq)
                LOGGER.info("""
                    seqx %s
                    seq1 %s
                    seq2 %s
                    """, seq, seq1, seq2)
                inloc.append(name+"_0" + spacer * " " + seq1)
                inloc.append(name+"_1" + spacer * " " + seq2)
            except ValueError as inst:
                LOGGER.debug("Found empty locus, all samples filtered.")
        if inloc:
            inloc.append("//  "+lines[-1][2:]+"|")
            locs.append("\n".join(inloc))

    return "\n".join(locs) + "\n"



def locichunk(args):
    """
    Function from make_loci to apply to chunks. smask is sample mask.
    """
    ## parse args
    data, optim, pnames, snppad, smask, start, samplecov, locuscov, upper = args

    ## this slice
    hslice = [start, start+optim]

    ## get filter db info
    co5 = h5py.File(data.database, 'r')
    afilt = co5["filters"][hslice[0]:hslice[1], ]
    aedge = co5["edges"][hslice[0]:hslice[1], ]
    asnps = co5["snps"][hslice[0]:hslice[1], ]

    ## get seqs db
    io5 = h5py.File(data.clust_database, 'r')
    if upper:
        aseqs = np.char.upper(io5["seqs"][hslice[0]:hslice[1], ])
    else:
        aseqs = io5["seqs"][hslice[0]:hslice[1], ]

    ## which loci passed all filters
    keep = np.where(np.sum(afilt, axis=1) == 0)[0]
    store = []

    ## write loci that passed after trimming edges, then write snp string
    for iloc in keep:
        edg = aedge[iloc]
        #LOGGER.info("!!!!!! iloc edg %s, %s", iloc, edg)
        args = [iloc, pnames, snppad, edg, aseqs, asnps, smask, samplecov, locuscov, start]
        if edg[4]:
            outstr, samplecov, locuscov = enter_pairs(*args)
            store.append(outstr)
        else:
            outstr, samplecov, locuscov = enter_singles(*args)
            store.append(outstr)

    ## write to file and clear store
    tmpo = os.path.join(data.dirs.outfiles, data.name+".loci.{}".format(start))
    with open(tmpo, 'w') as tmpout:
        tmpout.write("\n".join(store) + "\n")

    ## close handles
    io5.close()
    co5.close()

    ## return sample counter
    return samplecov, locuscov, start



def enter_pairs(iloc, pnames, snppad, edg, aseqs, asnps, smask, samplecov, locuscov, start):
    """ enters funcs for pairs """

    ## snps was created using only the selected samples.
    LOGGER.info("edges in enter_pairs %s", edg)
    seq1 = aseqs[iloc, :, edg[0]:edg[1]+1]
    snp1 = asnps[iloc, edg[0]:edg[1]+1, ]

    ## the 2nd read edges are +5 for the spacer
    seq2 = aseqs[iloc, :, edg[2]:edg[3]+1]
    snp2 = asnps[iloc, edg[2]:edg[3]+1, ]

    ## remove rows with all Ns, seq has only selected samples
    nalln = np.all(seq1 == "N", axis=1)

    ## make mask of removed rows and excluded samples. Use the inverse
    ## of this to save the coverage for samples
    nsidx = nalln + smask
    LOGGER.info("nsidx %s, nalln %s, smask %s", nsidx, nalln, smask)
    samplecov = samplecov + np.invert(nsidx).astype(np.int32)
    LOGGER.info("samplecov %s", samplecov)
    idx = np.sum(np.invert(nsidx).astype(np.int32))
    LOGGER.info("idx %s", idx)
    locuscov[idx] += 1

    ## select the remaining names in order
    seq1 = seq1[~nsidx, ]
    seq2 = seq2[~nsidx, ]
    names = pnames[~nsidx]

    ## save string for printing, excluding names not in samples
    outstr = "\n".join(\
        [name + s1.tostring()+"nnnn"+s2.tostring() for name, s1, s2 in \
         zip(names, seq1, seq2)])

    #LOGGER.info("s1 %s", s1.tostring())
    #LOGGER.info("s2 %s", s2.tostring())

    ## get snp string and add to store
    snpstring1 = ["-" if snp1[i, 0] else \
                 "*" if snp1[i, 1] else \
                 " " for i in range(len(snp1))]
    snpstring2 = ["-" if snp2[i, 0] else \
                 "*" if snp2[i, 1] else \
                 " " for i in range(len(snp2))]

    #npis = str(snpstring1+snpstring2).count("*")
    #nvars = str(snpstring1+snpstring2).count("-") + npis
    outstr += "\n" + snppad + "".join(snpstring1)+\
              "    "+"".join(snpstring2)+"|{}|".format(iloc+start)
              #"|LOCID={},DBID={},NVAR={},NPIS={}|"\
              #.format(1+iloc+start, iloc, nvars, npis)

    return outstr, samplecov, locuscov



def enter_singles(iloc, pnames, snppad, edg, aseqs, asnps, smask, samplecov, locuscov, start):
    """ enter funcs for SE or merged data """

    ## grab all seqs between edges
    seq = aseqs[iloc, :, edg[0]:edg[1]+1]
    ## snps was created using only the selected samples, and is edge masked.
    ## The mask is for counting snps quickly, but trimming is still needed here
    ## to make the snps line up with the seqs in the snp string.
    snp = asnps[iloc, edg[0]:edg[1]+1, ]

    ## remove rows with all Ns, seq has only selected samples
    nalln = np.all(seq == "N", axis=1)

    ## make mask of removed rows and excluded samples. Use the inverse
    ## of this to save the coverage for samples
    nsidx = nalln + smask
    samplecov = samplecov + np.invert(nsidx).astype(np.int32)
    idx = np.sum(np.invert(nsidx).astype(np.int32))
    locuscov[idx] += 1

    ## select the remaining names in order
    seq = seq[~nsidx, ]
    names = pnames[~nsidx]

    ## save string for printing, excluding names not in samples
    outstr = "\n".join(\
        [name + s.tostring() for name, s in zip(names, seq)])

    ## get snp string and add to store
    snpstring = ["-" if snp[i, 0] else \
                 "*" if snp[i, 1] else \
                 " " for i in range(len(snp))]
    outstr += "\n" + snppad + "".join(snpstring) + "|{}|".format(iloc+start)
    #LOGGER.info("outstr %s", outstr)
    return outstr, samplecov, locuscov



def init_arrays(data):
    """
    Create database file for storing final filtered snps data as hdf5 array.
    Copies splits and duplicates info from clust_database to database.
    """

    ## get stats from step6 h5 and create new h5
    co5 = h5py.File(data.clust_database, 'r')
    io5 = h5py.File(data.database, 'w')

    ## get maxlen and chunk len
    maxlen = data._hackersonly["max_fragment_length"] + 20
    chunks = co5["seqs"].attrs["chunksize"][0]
    nloci = co5["seqs"].shape[0]

    ## make array for snp string, 2 cols, - and *
    snps = io5.create_dataset("snps", (nloci, maxlen, 2),
                              dtype=np.bool,
                              chunks=(chunks, maxlen, 2),
                              compression='gzip')
    snps.attrs["chunksize"] = chunks
    snps.attrs["names"] = ["-", "*"]

    ## array for filters that will be applied in step7
    filters = io5.create_dataset("filters", (nloci, 6), dtype=np.bool)
    filters.attrs["filters"] = ["duplicates", "max_indels",
                                "max_snps", "max_shared_hets",
                                "min_samps", "max_alleles"]

    ## array for edgetrimming
    edges = io5.create_dataset("edges", (nloci, 5),
                               dtype=np.uint16,
                               chunks=(chunks, 5),
                               compression="gzip")
    edges.attrs["chunksize"] = chunks
    edges.attrs["names"] = ["R1_L", "R1_R", "R2_L", "R2_R", "sep"]

    ## xfer data from clustdb to finaldb
    edges[:, 4] = co5["splits"][:]
    filters[:, 0] = co5["duplicates"][:]

    ## close h5s
    io5.close()
    co5.close()



def filter_stacks(data, sidx, hslice):
    """
    Grab a chunk of loci from the HDF5 database. Apply filters and fill the
    the filters boolean array.

    The design of the filtering steps intentionally sacrifices some performance
    for an increase in readability, and extensibility. Calling multiple filter
    functions ends up running through the sequences per stack several times,
    but I felt this design made more sense, and also will easily allow us to
    add more filters in the future.
    """
    LOGGER.info("Entering filter_stacks")

    ## open h5 handles
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')
    ## get a chunk (hslice) of loci for the selected samples (sidx)
    #superseqs = io5["seqs"][hslice[0]:hslice[1], sidx,]
    ## get an int view of the seq array
    #superints = io5["seqs"][hslice[0]:hslice[1], sidx, :].view(np.int8)

    ## we need to use upper to skip lowercase allele storage
    ## this slows down the rate of loading in data by a ton.
    superints = np.char.upper(io5["seqs"][hslice[0]:hslice[1], sidx,]).view(np.int8)
    LOGGER.info("superints shape %s", superints)

    ## fill edge filter
    ## get edges of superseqs and supercats, since edges need to be trimmed
    ## before counting hets, snps, inds. Technically, this could edge trim
    ## clusters to the point that they are below the minlen, and so this
    ## also constitutes a filter, though one that is uncommon. For this
    ## reason we have another filter called edgfilter.
    splits = co5["edges"][hslice[0]:hslice[1], 4]
    edgfilter, edgearr = get_edges(data, superints, splits)
    del splits
    LOGGER.info('passed edges %s', hslice[0])

    ## minsamp coverages filtered from superseqs
    minfilter = filter_minsamp(data, superints)
    LOGGER.info('passed minfilt %s', hslice[0])

    ## maxhets per site column from superseqs after trimming edges
    hetfilter = filter_maxhet(data, superints, edgearr)
    LOGGER.info('passed minhet %s', hslice[0])

    ## ploidy filter
    pldfilter = io5["nalleles"][hslice[0]:hslice[1]].max(axis=1) > \
                                         data.paramsdict["max_alleles_consens"]

    ## indel filter, needs a fresh superints b/c get_edges does (-)->(N)
    indfilter = filter_indels(data, superints, edgearr)
    LOGGER.info('passed minind %s', hslice[0])

    ## Build the .loci snpstring as an array (snps)
    ## shape = (chunk, 1) dtype=S1, or should it be (chunk, 2) for [-,*] ?
    snpfilter, snpsarr = filter_maxsnp(data, superints, edgearr)

    LOGGER.info("edg %s", edgfilter.sum())
    LOGGER.info("min %s", minfilter.sum())
    LOGGER.info("het %s", hetfilter.sum())
    LOGGER.info("pld %s", pldfilter.sum())
    LOGGER.info("snp %s", snpfilter.sum())
    LOGGER.info("ind %s", indfilter.sum())

    ## SAVE FILTERS AND INFO TO DISK BY SLICE NUMBER (.0.tmp.h5)
    chunkdir = os.path.join(data.dirs.outfiles, data.name+"_tmpchunks")

    handle = os.path.join(chunkdir, "edgf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, edgfilter)

    handle = os.path.join(chunkdir, "minf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, minfilter)

    handle = os.path.join(chunkdir, "hetf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, hetfilter)

    handle = os.path.join(chunkdir, "snpf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, snpfilter)

    handle = os.path.join(chunkdir, "pldf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, pldfilter)

    handle = os.path.join(chunkdir, "indf.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, indfilter)

    handle = os.path.join(chunkdir, "snpsarr.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, snpsarr)

    handle = os.path.join(chunkdir, "edgearr.{}.npy".format(hslice[0]))
    with open(handle, 'w') as out:
        np.save(out, edgearr)

    io5.close()
    co5.close()



def get_edges(data, superints, splits):
    """
    Gets edge trimming based on the overlap of sequences at the edges of
    alignments and the tuple arg passed in for edge_trimming. Trims as
    (R1 left, R1 right, R2 left, R2 right). We also trim off the restriction
    site if it present. This modifies superints, and so should be run on an
    engine so it doesn't affect local copy. If this is changed to run locally
    for some reason make sure we copy the superints instead.
    """
    ## the filtering arg and parse it into minsamp numbers
    if "trim_overhang" in data.paramsdict:
        edgetrims = np.array(data.paramsdict["trim_overhang"]).astype(np.int16)
    else:
        edgetrims = np.array(data.paramsdict["trim_loci"]).astype(np.int16)

    ## Cuts 3 and 4 are only for 3rad/radcap
    ## TODO: This is moderately hackish, it's not using cut3/4
    ## correctly, just assuming the length is the same as cut1/2
    try:
        cut1, cut2, _, _ = data.paramsdict["restriction_overhang"]
        LOGGER.debug("Found 3Rad cut sites")
    except ValueError:
        cut1, cut2 = data.paramsdict["restriction_overhang"]

    cuts = np.array([len(cut1), len(cut2)], dtype=np.int16)

    ## a local array for storing edge trims
    edges = np.zeros((superints.shape[0], 5), dtype=np.int16)

    ## a local array for storing edge filtered loci, these are stored
    ## eventually as minsamp excludes.
    edgefilter = np.zeros((superints.shape[0],), dtype=np.bool)

    ## TRIM GUIDE. The cut site lengths are always trimmed. In addition,
    ## edge overhangs are trimmed to min(4, minsamp), and then additional
    ## number of columns is trimmed based on edgetrims values.
    ## A special case, -1 value means no trim at all.
    if data.paramsdict["min_samples_locus"] <= 4:
        minedge = np.int16(data.paramsdict["min_samples_locus"])
    else:
        minedge = np.int16(max(4, data.paramsdict["min_samples_locus"]))

    ## convert all - to N to make this easier
    nodashints = copy.deepcopy(superints)#.copy()
    nodashints[nodashints == 45] = 78

    ## trim overhanging edges
    ## get the number not Ns in each site,
    #ccx = np.sum(superseqs != "N", axis=1)
    ccx = np.sum(nodashints != 78, axis=1, dtype=np.uint16)
    efi, edg = edgetrim_numba(splits, ccx, edges, edgefilter, edgetrims, cuts, minedge)
    return efi, edg



@numba.jit(nopython=True)
def edgetrim_numba(splits, ccx, edges, edgefilter, edgetrims, cuts, minedge):
    ## get splits
    for idx in xrange(splits.shape[0]):

        ## get the data, single or paired
        if splits[idx]:
            r1s = ccx[idx, :splits[idx]]
            r2s = ccx[idx, splits[idx]+4:]
        else:
            r1s = ccx[idx, :]

        ## if 0 or positive, get edge 0
        x = np.where(r1s >= minedge)[0]
        if edgetrims[0] >= 0:
            if x.shape[0] > 20:
                edges[idx][0] = np.min(x[cuts[0]+edgetrims[0]:])
            else:
                edges[idx][0] = np.uint16(0)
                edgefilter[idx] = True

        ## fill in edge 1
        if edgetrims[1] >= 0:
            if x.shape[0] > 20:
                edges[idx][1] = np.max(x-edgetrims[1])
            else:
                edges[idx][1] = np.uint16(1)
                edgefilter[idx] = True

        ## If paired, do second read edges
        if splits[idx]:
            x = np.where(r2s >= minedge)[0]
            if edgetrims[2] >= 0:
                if x.shape[0] > 20:
                    edges[idx][2] = splits[idx] + np.int16(4) + np.min(x) + edgetrims[2]
                else:
                    edges[idx][2] = edges[idx][1] + 4
                    edgefilter[idx] = True

            if edgetrims[3] >= 0:
                ## get farthest site that is not all Ns
                if x.shape[0] > 20:
                    ## return index w/ spacers
                    edges[idx][3] = (splits[idx] + np.int16(4) + np.max(x)) - edgetrims[3]
                else:
                    edges[idx][3] = edges[idx][2] + 1
                    edgefilter[idx] = True

            ## enter the pair splitter
            edges[idx][4] = splits[idx]

        ## sanity check filter on edges
        if (edges[idx][1] < edges[idx][0]) or (edges[idx][3] < edges[idx][2]):
            edgefilter[idx] = True

        ## minlen check on edges (set arbitrarily to 20 right now)
        #if ((edges[idx][1] - edges[idx][0]) < 20) or \
        #   ((edges[idx][3] - edges[idx][2]) < 20):
        #    edgefilter[idx] = True

    return edgefilter, edges



def filter_minsamp(data, superints):
    """
    Filter minimum # of samples per locus from superseqs[chunk]. The shape
    of superseqs is [chunk, sum(sidx), maxlen]
    """
    ## global minsamp
    minsamp = data.paramsdict["min_samples_locus"]

    ## use population minsamps
    if data.populations:
        ## data._populations will look like this:
        ## {'a': (3, [0, 1, 2, 3],
        ##  'b': (3, [4, 5, 6, 7],
        ##  'c': (3, [8, 9, 10, 11]}
        LOGGER.info("POPULATIONS %s", data.populations)
        
        ## superints has already been subsampled by sidx
        ## get the sidx values for each pop
        minfilters = []
        for pop in data._populations:
            samps = data._populations[pop][1]
            minsamp = data._populations[pop][0]
            mini = np.sum(~np.all(superints[:, samps, :] == 78, axis=2), axis=1) < minsamp
            minfilters.append(mini)
        ## get sum across all pops for each locus
        minfilt = np.any(minfilters, axis=0)

    else:
        ## if not pop-file use global minsamp filter
        minfilt = np.sum(~np.all(superints == 78, axis=2), axis=1) < minsamp
        #LOGGER.info("Filtered by min_samples_locus - {}".format(minfilt.sum()))
    return minfilt



def ucount(sitecol):
    """
    Used to count the number of unique bases in a site for snpstring.
    returns as a spstring with * and -
    """

    ## a list for only catgs
    catg = [i for i in sitecol if i in "CATG"]

    ## find sites that are ambigs
    where = [sitecol[sitecol == i] for i in "RSKYWM"]

    ## for each occurrence of RSKWYM add ambig resolution to catg
    for ambig in where:
        for _ in range(ambig.size):
            catg += list(AMBIGS[ambig[0]])

    ## if invariant return " "
    if len(set(catg)) < 2:
        return " "
    else:
        ## get second most common site
        second = Counter(catg).most_common()[1][1]
        if second > 1:
            return "*"
        else:
            return "-"



def filter_maxsnp(data, superints, edgearr):
    """
    Filter max # of SNPs per locus. Do R1 and R2 separately if PE.
    Also generate the snpsite line for the .loci format and save in the snp arr
    This uses the edge filters that have been built based on trimming, and
    saves the snps array with edges filtered. **Loci are not yet filtered.**
    """

    ## an empty array to count with failed loci
    snpfilt = np.zeros(superints.shape[0], dtype=np.bool)
    snpsarr = np.zeros((superints.shape[0], superints.shape[2], 2), dtype=np.bool)
    maxsnps = np.array(data.paramsdict['max_SNPs_locus'], dtype=np.int16)

    ## get the per site snp string | shape=(chunk, maxlen)
    # snpsarr[:, :, 0] = snps == "-"
    # snpsarr[:, :, 1] = snps == "*"
    snpsarr = snpcount_numba(superints, snpsarr)
    LOGGER.info("---found the snps: %s", snpsarr.sum())
    snpfilt, snpsarr = snpfilter_numba(snpsarr, snpfilt, edgearr, maxsnps)
    LOGGER.info("---filtered snps: %s", snpfilt.sum())
    return snpfilt, snpsarr



@numba.jit(nopython=True)
def snpfilter_numba(snpsarr, snpfilt, edgearr, maxsnps):
    for idx in xrange(snpsarr.shape[0]):
        if not edgearr[idx, 4]:
            ## exclude snps that are outside of the edges
            for sidx in xrange(snpsarr.shape[1]):
                if sidx < edgearr[idx, 0]:
                    snpsarr[idx, sidx, :] = False
                elif sidx > edgearr[idx, 1]:
                    snpsarr[idx, sidx, :] = False

            nvar = snpsarr[idx, :].sum()
            if nvar > maxsnps[0]:
                snpfilt[idx] = True

        else:
            for sidx in xrange(snpsarr.shape[1]):
                if sidx < edgearr[idx, 0]:
                    snpsarr[idx, sidx, :] = False
                else:
                    if sidx > edgearr[idx, 1]:
                        if sidx < edgearr[idx, 2]:
                            snpsarr[idx, sidx, :] = False
                    if sidx > edgearr[idx, 3]:
                        snpsarr[idx, sidx, :] = False

            nvar1 = snpsarr[idx, :][:edgearr[idx, 4]].sum()
            if nvar1 > maxsnps[0]:
                snpfilt[idx] = True
            nvar2 = snpsarr[idx, :][edgearr[idx, 4]:].sum()
            if nvar2 > maxsnps[1]:
                snpfilt[idx] = True

    return snpfilt, snpsarr



@numba.jit(nopython=True)
def snpcount_numba(superints, snpsarr):
    """
    Used to count the number of unique bases in a site for snpstring.
    """
    ## iterate over all loci
    for iloc in xrange(superints.shape[0]):
        for site in xrange(superints.shape[2]):

            ## make new array
            catg = np.zeros(4, dtype=np.int16)

            ## a list for only catgs
            ncol = superints[iloc, :, site]
            for idx in range(ncol.shape[0]):
                if ncol[idx] == 67: #C
                    catg[0] += 1
                elif ncol[idx] == 65: #A
                    catg[1] += 1
                elif ncol[idx] == 84: #T
                    catg[2] += 1
                elif ncol[idx] == 71: #G
                    catg[3] += 1
                elif ncol[idx] == 82: #R
                    catg[1] += 1        #A
                    catg[3] += 1        #G
                elif ncol[idx] == 75: #K
                    catg[2] += 1        #T
                    catg[3] += 1        #G
                elif ncol[idx] == 83: #S
                    catg[0] += 1        #C
                    catg[3] += 1        #G
                elif ncol[idx] == 89: #Y
                    catg[0] += 1        #C
                    catg[2] += 1        #T
                elif ncol[idx] == 87: #W
                    catg[1] += 1        #A
                    catg[2] += 1        #T
                elif ncol[idx] == 77: #M
                    catg[0] += 1        #C
                    catg[1] += 1        #A


            ## get second most common site
            catg.sort()
            ## if invariant e.g., [0, 0, 0, 9], then nothing (" ")
            if not catg[2]:
                pass
            else:
                if catg[2] > 1:
                    snpsarr[iloc, site, 1] = True
                else:
                    snpsarr[iloc, site, 0] = True
    return snpsarr



def filter_maxhet(data, superints, edgearr):
    """
    Filter max shared heterozygosity per locus. The dimensions of superseqs
    are (chunk, sum(sidx), maxlen). Don't need split info since it applies to
    entire loci based on site patterns (i.e., location along the seq doesn't
    matter.) Current implementation does ints, but does not apply float diff
    to every loc based on coverage...
    """
    ## the filter max
    ## The type of max_shared_Hs_locus is determined and the cast to either
    ## int or float is made at assembly load time
    maxhet = data.paramsdict["max_shared_Hs_locus"]
    if isinstance(maxhet, float):
        ## get an array with maxhet fraction * ntaxa with data for each locus
        #maxhet = np.array(superints.shape[1]*maxhet, dtype=np.int16)
        maxhet = np.floor(
            maxhet * (superints.shape[1] - 
                np.all(superints == 78, axis=2).sum(axis=1))).astype(np.int16)
    elif isinstance(maxhet, int):
        maxhet = np.zeros(superints.shape[0], dtype=np.int16)
        maxhet.fill(data.paramsdict["max_shared_Hs_locus"])

    ## an empty array to fill with failed loci
    LOGGER.info("--------------maxhet mins %s", maxhet)
    hetfilt = np.zeros(superints.shape[0], dtype=np.bool)
    hetfilt = maxhet_numba(superints, edgearr, maxhet, hetfilt)
    LOGGER.info("--------------maxhet sums %s", hetfilt.sum())
    return hetfilt



def filter_indels(data, superints, edgearr):
    """
    Filter max indels. Needs to split to apply to each read separately.
    The dimensions of superseqs are (chunk, sum(sidx), maxlen).
    """

    maxinds = np.array(data.paramsdict["max_Indels_locus"]).astype(np.int64)

    ## an empty array to fill with failed loci
    ifilter = np.zeros(superints.shape[0], dtype=np.bool_)

    ## if paired then worry about splits
    if "pair" in data.paramsdict["datatype"]:
        for idx in xrange(superints.shape[0]):
            block1 = superints[idx, :, edgearr[idx, 0]:edgearr[idx, 1]]
            block2 = superints[idx, :, edgearr[idx, 2]:edgearr[idx, 3]]

            sums1 = maxind_numba(block1)
            sums2 = maxind_numba(block2)

            if (sums1 > maxinds[0]) or (sums2 > maxinds[1]):
                ifilter[idx] = True

    else:
        for idx in xrange(superints.shape[0]):
            ## get block based on edge filters
            block = superints[idx, :, edgearr[idx, 0]:edgearr[idx, 1]]
            ## shorten block to exclude terminal indels
            ## if data at this locus (not already filtered by edges/minsamp)
            if block.shape[1] > 1:
                sums = maxind_numba(block)
                #LOGGER.info("maxind numba %s %s", idx, sums)
                #LOGGER.info("sums, maxinds[0], compare: %s %s %s",
                #             sums, maxinds[0], sums > maxinds[0])
                if sums > maxinds[0]:
                    ifilter[idx] = True

    LOGGER.info("--------------maxIndels sums %s", ifilter.sum())
    return ifilter



@numba.jit(nopython=True)
def maxind_numba(block):
    """ filter for indels """
    ## remove terminal edges
    inds = 0
    for row in xrange(block.shape[0]):
        where = np.where(block[row] != 45)[0]
        left = np.min(where)
        right = np.max(where)
        obs = np.sum(block[row, left:right] == 45)
        if obs > inds:
            inds = obs
    return inds



@numba.jit(nopython=True)
def maxhet_numba(superints, edgearr, maxhet, hetfilt):
    for idx in xrange(edgearr.shape[0]):
        ## do read1s
        fcolumns = superints[idx, :, edgearr[idx, 0]:edgearr[idx, 1]]
        count1s = np.zeros(fcolumns.shape[1], dtype=np.int16)
        for fidx in xrange(fcolumns.shape[1]):
            subcount = 0
            for ambig in AMBIGARR:
                subcount += np.sum(fcolumns[:, fidx] == ambig)
            count1s[fidx] = subcount

        ## do read2s
        fcolumns = superints[idx, :, edgearr[idx, 2]:edgearr[idx, 3]]
        count2s = np.zeros(fcolumns.shape[1], dtype=np.int16)
        for fidx in xrange(fcolumns.shape[1]):
            subcount = 0
            for ambig in AMBIGARR:
                subcount += np.sum(fcolumns[:, fidx] == ambig)
            count2s[fidx] = subcount

        ## check against max
        if count1s.max() > maxhet[idx]:
            hetfilt[idx] = True
        elif count2s.size:
            if count2s.max() > maxhet[idx]:
                hetfilt[idx] = True
    return hetfilt

## MAKE GLOBAL
AMBIGARR = np.array(list("RSKYWM")).view(np.int8)




def make_outfiles(data, samples, output_formats, ipyclient):
    """
    Get desired formats from paramsdict and write files to outfiles
    directory.
    """

    ## will iterate optim loci at a time
    with h5py.File(data.clust_database, 'r') as io5:
        optim = io5["seqs"].attrs["chunksize"][0]
        nloci = io5["seqs"].shape[0]

        ## get name and snp padding
        anames = io5["seqs"].attrs["samples"]
        snames = [i.name for i in samples]
        ## get only snames in this data set sorted in the order they are in io5
        names = [i for i in anames if i in snames]
        pnames, _ = padnames(names)


    ## get names boolean
    sidx = np.array([i in snames for i in anames])
    assert len(pnames) == sum(sidx)

    ## get names index in order of pnames
    #sindx = [list(anames).index(i) for i in snames]

    ## send off outputs as parallel jobs
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    results = {}

    ## build arrays and outputs from arrays.
    ## these arrays are keys in the tmp h5 array: seqarr, snparr, bisarr, maparr
    boss_make_arrays(data, sidx, optim, nloci, ipyclient)

    start = time.time()
    ## phy and partitions are a default output ({}.phy, {}.phy.partitions)
    if "p" in output_formats:
        data.outfiles.phy = os.path.join(data.dirs.outfiles, data.name+".phy")
        async = lbview.apply(write_phy, *[data, sidx, pnames])
        results['phy'] = async

    ## nexus format includes ... additional information ({}.nex)
    if "n" in output_formats:
        data.outfiles.nexus = os.path.join(data.dirs.outfiles, data.name+".nex")
        async = lbview.apply(write_nex, *[data, sidx, pnames])
        results['nexus'] = async

    ## snps is actually all snps written in phylip format ({}.snps.phy)
    if "s" in output_formats:
        data.outfiles.snpsmap = os.path.join(data.dirs.outfiles, data.name+".snps.map")
        data.outfiles.snpsphy = os.path.join(data.dirs.outfiles, data.name+".snps.phy")
        async = lbview.apply(write_snps, *[data, sidx, pnames])
        results['snps'] = async
        async = lbview.apply(write_snps_map, data)
        results['snpsmap'] = async

    ## usnps is one randomly sampled snp from each locus ({}.u.snps.phy)
    if "u" in output_formats:
        data.outfiles.usnpsphy = os.path.join(data.dirs.outfiles, data.name+".u.snps.phy")
        async = lbview.apply(write_usnps, *[data, sidx, pnames])
        results['usnps'] = async

    ## str and ustr are for structure analyses. A fairly outdated format, six
    ## columns of empty space. Full and subsample included ({}.str, {}.u.str)
    if "k" in output_formats:
        data.outfiles.str = os.path.join(data.dirs.outfiles, data.name+".str")
        data.outfiles.ustr = os.path.join(data.dirs.outfiles, data.name+".ustr")        
        async = lbview.apply(write_str, *[data, sidx, pnames])
        results['structure'] = async

    ## geno output is for admixture and other software. We include all SNPs,
    ## but also a .map file which has "distances" between SNPs.
    if 'g' in output_formats:
        data.outfiles.geno = os.path.join(data.dirs.outfiles, data.name+".geno")
        data.outfiles.ugeno = os.path.join(data.dirs.outfiles, data.name+".u.geno")
        async = lbview.apply(write_geno, *[data, sidx])
        results['geno'] = async

    ## G-PhoCS output. Have to use cap G here cuz little g is already taken, lol.
    if 'G' in output_formats:
        data.outfiles.gphocs = os.path.join(data.dirs.outfiles, data.name+".gphocs")
        async = lbview.apply(write_gphocs, *[data, sidx])
        results['gphocs'] = async

    ## wait for finished outfiles
    while 1:
        readies = [i.ready() for i in results.values()]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(readies), sum(readies),
            " writing outfiles      | {} | s7 |".format(elapsed), 
            spacer=data._spacer)
        time.sleep(0.1)
        if all(readies):
            break
    print("")

    ## check for errors
    for suff, async in results.items():
        if not async.successful():
            print("  Warning: error encountered while writing {} outfile: {}"\
                  .format(suff, async.exception()))
            LOGGER.error("  Warning: error in writing %s outfile: %s", \
                         suff, async.exception())

    ## remove the tmparrays
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name))
    os.remove(tmparrs)



def boss_make_arrays(data, sidx, optim, nloci, ipyclient):
    
    ## make a list of slices to distribute in parallel
    hslices = [start for start in range(0, nloci, optim)]
    
    ## parallel client setup
    lbview = ipyclient.load_balanced_view()
        
    ## load the h5 database and grab some needed info
    maxlen = data._hackersonly["max_fragment_length"] + 20

    ## distribute jobs in chunks of optim
    asyncs = []
    for hslice in hslices:
        args = (data, sidx, hslice, optim, maxlen)
        async = lbview.apply(worker_make_arrays, *args)
        asyncs.append(async)
               
    ## shape of arrays is sidx, we will subsample h5 w/ sidx to match. Seqarr
    ## will actually be a fair bit smaller since its being loaded with edge
    with h5py.File(data.database, 'r') as co5:
        maxsnp = co5["snps"][:].sum()
        afilt = co5["filters"][:]
        nkeeps = np.sum(np.sum(afilt, axis=1) == 0)

    ## a tmp h5 to hold working arrays (the seq array is not filtered and trimmed
    ## It is the phylip output, essentially.
    h5name = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name))
    with h5py.File(h5name, 'w') as tmp5:

        ## ensure chunksize is not greater than array size
        shape1 = maxlen*nkeeps
        tmp5.create_dataset("seqarr", (sum(sidx), shape1), dtype="S1", 
                            chunks=(sum(sidx), min(shape1, maxlen*optim)))
        tmp5.create_dataset("snparr", (sum(sidx), maxsnp), dtype="S1", 
                            chunks=(sum(sidx), min(maxsnp, optim)))
        tmp5.create_dataset('bisarr', (sum(sidx), nkeeps), dtype="S1", 
                            chunks=(sum(sidx), min(nkeeps, optim)))
        tmp5.create_dataset('maparr', (maxsnp, 4), dtype=np.uint32)

        ## track progress, catch errors, and enter results into h5 as it arrive
        start = time.time()
        njobs = len(asyncs)
        ## axis1 counters
        seqidx = snpidx = bisidx = mapidx = locidx = 0
        
        while 1:
            ## we need to collect results in order!
            if asyncs[0].ready():
                if asyncs[0].successful():
                    ## enter results and del async
                    seqarr, snparr, bisarr, maparr = asyncs[0].result()
                    tmp5["seqarr"][:, seqidx:seqidx+seqarr.shape[1]] = seqarr
                    seqidx += seqarr.shape[1]
                    tmp5["snparr"][:, snpidx:snpidx+snparr.shape[1]] = snparr
                    snpidx += snparr.shape[1]
                    tmp5["bisarr"][:, bisidx:bisidx+bisarr.shape[1]] = bisarr
                    bisidx += bisarr.shape[1]

                    ## mapfile needs idxs summed, only bother if there is data
                    ## that passed filtering for this chunk
                    if maparr.shape[0]:
                        mapcopy = maparr.copy()
                        mapcopy[:, 0] += locidx
                        mapcopy[:, 3] += mapidx
                        tmp5["maparr"][mapidx:mapidx+maparr.shape[0], :] = mapcopy
                        locidx = mapcopy[-1, 0]
                        mapidx += mapcopy.shape[0]
                    del asyncs[0]
                        
                else:
                    print(asyncs[0].exception())

            ## print progress            
            time.sleep(0.1)
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(njobs, njobs-len(asyncs), 
                " building arrays       | {} | s7 |".format(elapsed), 
                spacer=data._spacer)
            
            ## are we done?
            if not asyncs:
                print("")
                break

    #return 
    

    
def worker_make_arrays(data, sidx, hslice, optim, maxlen):
    """
    Parallelized worker to build array chunks for output files. One main 
    goal here is to keep seqarr to less than ~1GB RAM.
    """
    
    ## big data arrays
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')
    
    ## temporary storage until writing to h5 array    
    maxsnp = co5["snps"][hslice:hslice+optim].sum()         ## concat later
    maparr = np.zeros((maxsnp, 4), dtype=np.uint32)
    snparr = np.zeros((sum(sidx), maxsnp), dtype="S1")
    bisarr = np.zeros((sum(sidx), maxsnp), dtype="S1")
    seqarr = np.zeros((sum(sidx), maxlen*optim), dtype="S1")

    ## apply all filters and write loci data
    seqleft = 0
    snpleft = 0
    bis = 0                     

    ## edge filter has already been applied to snps, but has not yet been
    ## applied to seqs. The locus filters have not been applied to either yet.
    mapsnp = 0
    totloc = 0

    afilt = co5["filters"][hslice:hslice+optim, :]
    aedge = co5["edges"][hslice:hslice+optim, :]
    asnps = co5["snps"][hslice:hslice+optim, :]
    #aseqs = io5["seqs"][hslice:hslice+optim, sidx, :]
    ## have to run upper on seqs b/c they have lowercase storage of alleles
    aseqs = np.char.upper(io5["seqs"][hslice:hslice+optim, sidx, :])

    ## which loci passed all filters
    keep = np.where(np.sum(afilt, axis=1) == 0)[0]

    ## write loci that passed after trimming edges, then write snp string
    for iloc in keep:
        ## grab r1 seqs between edges
        edg = aedge[iloc]

        ## grab SNPs from seqs already sidx subsampled and edg masked.
        ## needs to be done here before seqs are edgetrimmed.
        getsnps = asnps[iloc].sum(axis=1).astype(np.bool)
        snps = aseqs[iloc, :, getsnps].T

        ## trim edges and split from seqs and concatenate for pairs.
        ## this seq array will be the phy output.
        if not "pair" in data.paramsdict["datatype"]:
            seq = aseqs[iloc, :, edg[0]:edg[1]+1]
        else:
            seq = np.concatenate([aseqs[iloc, :, edg[0]:edg[1]+1],
                                  aseqs[iloc, :, edg[2]:edg[3]+1]], axis=1)

        ## remove cols from seq (phy) array that are all N-
        lcopy = seq
        lcopy[lcopy == "-"] = "N"
        bcols = np.all(lcopy == "N", axis=0)
        seq = seq[:, ~bcols]

        ## put into large array (could put right into h5?)
        seqarr[:, seqleft:seqleft+seq.shape[1]] = seq
        seqleft += seq.shape[1]

        ## subsample all SNPs into an array
        snparr[:, snpleft:snpleft+snps.shape[1]] = snps
        snpleft += snps.shape[1]

        ## Enter each snp into the map file
        for i in xrange(snps.shape[1]):
            ## 1-indexed loci in first column
            ## actual locus number in second column
            ## counter for this locus in third column
            ## snp counter total in fourth column
            maparr[mapsnp, :] = [totloc+1, hslice+iloc, i, mapsnp+1]
            mapsnp += 1

        ## subsample one SNP into an array
        if snps.shape[1]:
            samp = np.random.randint(snps.shape[1])
            bisarr[:, bis] = snps[:, samp]
            bis += 1
            totloc += 1
            
    ## clean up
    io5.close()
    co5.close()
    
    ## trim trailing edges b/c we made the array bigger than needed.
    ridx = np.all(seqarr == "", axis=0)
    seqarr = seqarr[:, ~ridx]
    ridx = np.all(snparr == "", axis=0)
    snparr = snparr[:, ~ridx]
    ridx = np.all(bisarr == "", axis=0)
    bisarr = bisarr[:, ~ridx]
    ridx = np.all(maparr == 0, axis=1)
    maparr = maparr[~ridx, :]

    ## return these three arrays which are pretty small
    ## catg array gets to be pretty huge, so we return only
    return seqarr, snparr, bisarr, maparr
  
  

def write_phy(data, sidx, pnames):
    """ 
    write the phylip output file from the tmparr[seqarray] 
    """

    ## grab seq data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        seqarr = io5["seqarr"]

        ## trim to size b/c it was made longer than actual
        end = np.where(np.all(seqarr[:] == "", axis=0))[0]
        if np.any(end):
            end = end.min()
        else:
            end = seqarr.shape[1] 

        ## write to phylip 
        with open(data.outfiles.phy, 'w') as out:
            ## write header
            out.write("{} {}\n".format(seqarr.shape[0], end))

            ## write data rows
            for idx, name in enumerate(pnames):
                out.write("{}{}\n".format(name, "".join(seqarr[idx, :end])))
    LOGGER.debug("finished writing phy in: %s", time.time() - start)



def write_nex(data, sidx, pnames):
    """ 
    write the nexus output file from the tmparr[seqarray] and tmparr[maparr]
    """    

    ## grab seq data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        seqarr = io5["seqarr"]

        ## trim to size b/c it was made longer than actual
        end = np.where(np.all(seqarr[:] == "", axis=0))[0]
        if np.any(end):
            end = end.min()
        else:
            end = seqarr.shape[1] 

        ## write to nexus
        data.outfiles.nex = os.path.join(data.dirs.outfiles, data.name+".nex")

        with open(data.outfiles.nex, 'w') as out:

            ## write nexus seq header
            out.write(NEXHEADER.format(seqarr.shape[0], end))

            ## grab a big block of data
            chunksize = 100000  # this should be a multiple of 100
            for bidx in xrange(0, end, chunksize):
                bigblock = seqarr[:, bidx:bidx+chunksize]
                lend = end-bidx
                #LOGGER.info("BIG: %s %s %s %s", bigblock.shape, bidx, lend, end)

                ## write interleaved seqs 100 chars with longname+2 before
                tmpout = []            
                for block in xrange(0, min(chunksize, lend), 100):
                    stop = min(block+100, end)

                    for idx, name in enumerate(pnames):
                        seqdat = bigblock[idx, block:stop]
                        tmpout.append("  {}{}\n".format(name, "".join(seqdat)))
                    tmpout.append("\n")

                ## print intermediate result and clear
                if any(tmpout):
                    out.write("".join(tmpout))
            ## closer
            out.write(NEXCLOSER)
    LOGGER.debug("finished writing nex in: %s", time.time() - start)



## TODO: this could have much more information for reference aligned data
def write_snps_map(data):
    """ write a map file with linkage information for SNPs file"""

    ## grab map data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        maparr = io5["maparr"][:]

        ## get last data 
        end = np.where(np.all(maparr[:] == 0, axis=1))[0]
        if np.any(end):
            end = end.min()
        else:
            end = maparr.shape[0]

        ## write to map file (this is too slow...)
        outchunk = []
        with open(data.outfiles.snpsmap, 'w') as out:
            for idx in xrange(end):
                ## build to list
                line = maparr[idx, :]
                #print(line)
                outchunk.append(\
                    "{}\trad{}_snp{}\t{}\t{}\n"\
                    .format(line[0], line[1], line[2], 0, line[3]))
                ## clear list
                if not idx % 10000:
                    out.write("".join(outchunk))
                    outchunk = []
            ## write remaining
            out.write("".join(outchunk))
    LOGGER.debug("finished writing snps_map in: %s", time.time() - start)



def write_snps(data, sidx, pnames):
    """ write the snp string """

    ## grab snp data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        snparr = io5["snparr"]

        ## trim to size b/c it was made longer than actual
        end = np.where(np.all(snparr[:] == "", axis=0))[0]
        if np.any(end):
            end = end.min()
        else:
            end = snparr.shape[1]        

        ## write to snps file
        with open(data.outfiles.snpsphy, 'w') as out:
            out.write("{} {}\n".format(snparr.shape[0], end))
            for idx, name in enumerate(pnames):
                out.write("{}{}\n".format(name, "".join(snparr[idx, :end])))
    LOGGER.debug("finished writing snps in: %s", time.time() - start)



def write_usnps(data, sidx, pnames):
    """ write the bisnp string """

    ## grab bis data from tmparr
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        bisarr = io5["bisarr"]

        ## trim to size b/c it was made longer than actual
        end = np.where(np.all(bisarr[:] == "", axis=0))[0]
        if np.any(end):
            end = end.min()
        else:
            end = bisarr.shape[1]        

        ## write to usnps file
        with open(data.outfiles.usnpsphy, 'w') as out:
            out.write("{} {}\n".format(bisarr.shape[0], end))
            for idx, name in enumerate(pnames):
                out.write("{}{}\n".format(name, "".join(bisarr[idx, :end])))



def write_str(data, sidx, pnames):
    """ Write STRUCTURE format for all SNPs and unlinked SNPs """

    ## grab snp and bis data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        snparr = io5["snparr"]
        bisarr = io5["bisarr"]

        ## trim to size b/c it was made longer than actual
        bend = np.where(np.all(bisarr[:] == "", axis=0))[0]
        if np.any(bend):
            bend = bend.min()
        else:
            bend = bisarr.shape[1]        

        send = np.where(np.all(snparr[:] == "", axis=0))[0]       
        if np.any(send):
            send = send.min()
        else:
            send = snparr.shape[1]        

        ## write to str and ustr
        out1 = open(data.outfiles.str, 'w')
        out2 = open(data.outfiles.ustr, 'w')
        numdict = {'A': '0', 'T': '1', 'G': '2', 'C': '3', 'N': '-9', '-': '-9'}
        if data.paramsdict["max_alleles_consens"] > 1:
            for idx, name in enumerate(pnames):
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in snparr[idx, :send]])))
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][1]] for i in snparr[idx, :send]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in bisarr[idx, :bend]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][1]] for i in bisarr[idx, :bend]])))
        else:
            ## haploid output
            for idx, name in enumerate(pnames):
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in snparr[idx, :send]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in bisarr[idx, :bend]])))
        out1.close()
        out2.close()
    LOGGER.debug("finished writing str in: %s", time.time() - start)



def write_geno(data, sidx):
    """
    write the geno output formerly used by admixture, still supported by 
    adegenet, perhaps. Also, sNMF still likes .geno.
    """

    ## grab snp and bis data from tmparr
    start = time.time()
    tmparrs = os.path.join(data.dirs.outfiles, "tmp-{}.h5".format(data.name)) 
    with h5py.File(tmparrs, 'r') as io5:
        snparr = io5["snparr"]
        bisarr = io5["bisarr"]

        ## trim to size b/c it was made longer than actual
        bend = np.where(np.all(bisarr[:] == "", axis=0))[0]
        if np.any(bend):
            bend = bend.min()
        else:
            bend = bisarr.shape[1]        

        send = np.where(np.all(snparr[:] == "", axis=0))[0]       
        if np.any(send):
            send = send.min()
        else:
            send = snparr.shape[1]        

        ## get most common base at each SNP as a pseudo-reference
        ## and record 0,1,2 or missing=9 for counts of the ref allele
        snpref = reftrick(snparr[:, :send].view(np.int8), GETCONS).view("S1")
        bisref = reftrick(bisarr[:, :bend].view(np.int8), GETCONS).view("S1")

        ## geno matrix to fill (9 is empty)
        snpgeno = np.zeros((snparr.shape[0], send), dtype=np.uint8)
        snpgeno.fill(9)
        bisgeno = np.zeros((bisarr.shape[0], bend), dtype=np.uint8)
        bisgeno.fill(9)

        ##--------------------------------------------------------------------
        ## fill in complete hits (match to first column ref base)
        mask2 = np.array(snparr[:, :send] == snpref[:, 0])
        snpgeno[mask2] = 2

        ## fill in single hits (heteros) match to hetero of first+second column
        ambref = np.apply_along_axis(lambda x: TRANSFULL[tuple(x)], 1, snpref[:, :2])
        mask1 = np.array(snparr[:, :send] == ambref)
        snpgeno[mask1] = 1

        ## fill in zero hits, meaning a perfect match to the second column base
        ## anything else is left at 9 (missing), b/c it's either missing or it
        ## is not bi-allelic. 
        mask0 = np.array(snparr[:, :send] == snpref[:, 1])
        snpgeno[mask0] = 0

        ##--------------------------------------------------------------------

        ## fill in complete hits
        mask2 = np.array(bisarr[:, :bend] == bisref[:, 0])
        bisgeno[mask2] = 2

        ## fill in single hits (heteros)
        ambref = np.apply_along_axis(lambda x: TRANSFULL[tuple(x)], 1, bisref[:, :2])
        mask1 = np.array(bisarr[:, :bend] == ambref)
        bisgeno[mask1] = 1

        ## fill in zero hits (match to second base)
        mask0 = np.array(bisarr[:, :bend] == bisref[:, 1])
        bisgeno[mask0] = 0

        ##---------------------------------------------------------------------
        ## print to files
        np.savetxt(data.outfiles.geno, snpgeno.T, delimiter="", fmt="%d")
        np.savetxt(data.outfiles.ugeno, bisgeno.T, delimiter="", fmt="%d")
    LOGGER.debug("finished writing geno in: %s", time.time() - start)



def write_gphocs(data, sidx):
    """
    write the g-phocs output. This code is hella ugly bcz it's copy/pasted
    directly from the old loci2gphocs script from pyrad. I figure having it
    get done the stupid way is better than not having it done at all, at
    least for the time being. This could probably be sped up significantly.
    """

    outfile = data.outfiles.gphocs
    infile = data.outfiles.loci

    infile = open(infile)
    outfile = open(outfile, 'w')

    ## parse the loci
    ## Each set of reads at a locus is appended with a line
    ## beginning with // and ending with |x, where x in the locus id.
    ## so after this call 'loci' will contain an array
    ## of sets of each read per locus.
    loci = re.compile("\|[0-9]+\|").split(infile.read())[:-1]

    # Print the header, the number of loci in this file
    outfile.write(str(len(loci)) + "\n\n")

    # iterate through each locus, print out the header for each locus:
    # <locus_name> <n_samples> <locus_length>
    # Then print the data for each sample in this format:
    # <individual_name> <sequence>
    for i, loc in enumerate(loci):
        ## Get rid of the line that contains the snp info
        loc = loc.rsplit("\n", 1)[0]

        # Separate out each sequence within the loc block. 'sequences'
        # will now be a list strings containing name/sequence pairs.
        # We select each line in the locus string that starts with ">"
        names = [line.split()[0] for line in loc.strip().split("\n")]
        try:
            sequences = [line.split()[1] for line in loc.strip().split("\n")]
        except:
            pass
        # Strips off 'nnnn' separator for paired data
        # replaces '-' with 'N'
        editsequences = [seq.replace("n","").replace('-','N') for seq in sequences]
        sequence_length = len(editsequences[0])

        # get length of longest name and add 4 spaces
        longname = max(map(len,names))+4

        # Print out the header for this locus
        outfile.write('locus{} {} {}\n'.format(str(i), len(sequences), sequence_length))

        # Iterate through each sequence read at this locus and write it to the file.
        for name,sequence in zip(names,sequences):
            # Clean up the sequence data to make gphocs happy. Only accepts UPPER
            # case chars for bases, and only accepts 'N' for missing data.
            outfile.write(name+" "*(longname-len(name))+sequence + "\n")
        ## Separate loci with so it's prettier
        outfile.write("\n")



def make_vcf(data, samples, ipyclient, full=0):
    """
    Write the full VCF for loci passing filtering. Other vcf formats are
    possible, like SNPs-only, or with filtered loci included but the filter
    explicitly labeled. These are not yet supported, however.
    """
    ## start vcf progress bar
    start = time.time()
    printstr = " building vcf file     | {} | s7 |"
    LOGGER.info("Writing .vcf file")
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, printstr.format(elapsed), spacer=data._spacer)

    ## create outputs for v and V, gzip V to be friendly
    data.outfiles.vcf = os.path.join(data.dirs.outfiles, data.name+".vcf")
    if full:
        data.outfiles.VCF = os.path.join(data.dirs.outfiles, data.name+".vcf.gz")

    ## get some db info
    with h5py.File(data.clust_database, 'r') as io5:
        ## will iterate optim loci at a time
        optim = io5["seqs"].attrs["chunksize"][0]
        nloci = io5["seqs"].shape[0]
        ## get name and snp padding
        anames = io5["seqs"].attrs["samples"]
        snames = [i.name for i in samples]
        names = [i for i in anames if i in snames]

    ## get names index
    sidx = np.array([i in snames for i in anames])

    ## client for sending jobs to parallel engines, for this step we'll limit
    ## to half of the available cpus if
    lbview = ipyclient.load_balanced_view()

    ## send jobs in chunks
    vasyncs = {}
    total = 0
    for chunk in xrange(0, nloci, optim):
        vasyncs[chunk] = lbview.apply(vcfchunk, *(data, optim, sidx, chunk, full))
        total += 1

    ## tmp files get left behind and intensive processes are left running when a
    ## a job is killed/interrupted during vcf build, so we try/except wrap.
    try:
        while 1:
            keys = [i for (i, j) in vasyncs.items() if j.ready()]
            ## check for failures
            for job in keys:
                if not vasyncs[job].successful():
                    ## raise exception
                    err = " error in vcf build chunk {}: {}"\
                          .format(job, vasyncs[job].result())
                    LOGGER.error(err)
                    raise IPyradWarningExit(err)
                else:
                    ## free up memory
                    del vasyncs[job]

            finished = total - len(vasyncs) #sum([i.ready() for i in vasyncs.values()])
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(total, finished, printstr.format(elapsed), spacer=data._spacer)
            time.sleep(0.5)
            if not vasyncs:
                break
        print("")

    except Exception as inst:
        ## make sure all future jobs are aborted
        keys = [i for (i, j) in vasyncs.items() if not j.ready()]
        try:
            for job in keys:
                #vasyncs[job].abort()
                vasyncs[job].cancel()
        except Exception:
            pass
        ## make sure all tmp files are destroyed
        vcfchunks = glob.glob(os.path.join(data.dirs.outfiles, "*.vcf.[0-9]*"))
        h5chunks = glob.glob(os.path.join(data.dirs.outfiles, ".tmp.[0-9]*.h5"))
        for dfile in vcfchunks+h5chunks:
            os.remove(dfile)
        ## reraise the error
        raise inst


    ## writing full vcf file to disk
    start = time.time()
    printstr = " writing vcf file      | {} | s7 |"
    res = lbview.apply(concat_vcf, *(data, names, full))
    ogchunks = len(glob.glob(data.outfiles.vcf+".*"))
    while 1:
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        curchunks = len(glob.glob(data.outfiles.vcf+".*"))
        progressbar(ogchunks, ogchunks-curchunks, printstr.format(elapsed), spacer=data._spacer)
        time.sleep(0.1)
        if res.ready():
            break
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(1, 1, printstr.format(elapsed), spacer=data._spacer)
    print("")



def concat_vcf(data, names, full):
    """
    Sorts, concatenates, and gzips VCF chunks. Also cleans up chunks.
    """
    ## open handle and write headers
    if not full:
        writer = open(data.outfiles.vcf, 'w')
    else:
        writer = gzip.open(data.outfiles.VCF, 'w')
    vcfheader(data, names, writer)
    writer.close()

    ## get vcf chunks
    vcfchunks = glob.glob(data.outfiles.vcf+".*")
    vcfchunks.sort(key=lambda x: int(x.rsplit(".")[-1]))

    ## concatenate
    if not full:
        writer = open(data.outfiles.vcf, 'a')
    else:
        writer = gzip.open(data.outfiles.VCF, 'a')

    ## what order do users want? The order in the original ref file?
    ## Sorted by the size of chroms? that is the order in faidx.
    ## If reference mapping then it's nice to sort the vcf data by
    ## CHROM and POS. This is doing a very naive sort right now, so the
    ## CHROM will be ordered, but not the pos within each chrom.
    if data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]:
        ## Some unix sorting magic to get POS sorted within CHROM
        ## First you sort by POS (-k 2,2), then you do a `stable` sort 
        ## by CHROM. You end up with POS ordered and grouped correctly by CHROM
        ## but relatively unordered CHROMs (locus105 will be before locus11).
        cmd = ["cat"] + vcfchunks + [" | sort -k 2,2 -n | sort -k 1,1 -s"]
        cmd = " ".join(cmd)
        proc = sps.Popen(cmd, shell=True, stderr=sps.STDOUT, stdout=writer, close_fds=True)
    else:
        proc = sps.Popen(["cat"] + vcfchunks, stderr=sps.STDOUT, stdout=writer, close_fds=True)

    err = proc.communicate()[0]
    if proc.returncode:
        raise IPyradWarningExit("err in concat_vcf: %s", err)
    writer.close()

    for chunk in vcfchunks:
        os.remove(chunk)



def vcfchunk(data, optim, sidx, chunk, full):
    """
    Function called within make_vcf to run chunks on separate engines.
    """
    ## empty array to be filled before writing
    ## will not actually be optim*maxlen, extra needs to be trimmed
    maxlen = data._hackersonly["max_fragment_length"] + 20

    ## get data sliced (optim chunks at a time)
    hslice = [chunk, chunk+optim]

    ## read all taxa from disk (faster), then subsample taxa with sidx and
    ## keepmask to greatly reduce the memory load
    with h5py.File(data.database, 'r') as co5:
        afilt = co5["filters"][hslice[0]:hslice[1], :]
        keepmask = afilt.sum(axis=1) == 0
        ## apply mask to edges
        aedge = co5["edges"][hslice[0]:hslice[1], :]
        aedge = aedge[keepmask, :]
    del afilt
    ## same memory subsampling.
    with h5py.File(data.clust_database, 'r') as io5:
        ## apply mask to edges to aseqs and acatg
        #aseqs = io5["seqs"][hslice[0]:hslice[1], :, :].view(np.uint8)
        ## need to read in seqs with upper b/c lowercase allele info
        aseqs = np.char.upper(io5["seqs"][hslice[0]:hslice[1], :, :]).view(np.uint8)
        aseqs = aseqs[keepmask, :]
        aseqs = aseqs[:, sidx, :]
        acatg = io5["catgs"][hslice[0]:hslice[1], :, :, :]
        acatg = acatg[keepmask, :]
        acatg = acatg[:, sidx, :, :]
        achrom = io5["chroms"][hslice[0]:hslice[1]]
        achrom = achrom[keepmask, :]        

    LOGGER.info('acatg.shape %s', acatg.shape)
    ## to save memory some columns are stored in diff dtypes until printing
    if not full:
        with h5py.File(data.database, 'r') as co5:
            snps = co5["snps"][hslice[0]:hslice[1], :]
            snps = snps[keepmask, :]
            snps = snps.sum(axis=2)
        snpidxs = snps > 0
        maxsnplen = snps.sum()

    ## vcf info to fill, this is bigger than the actual array
    nrows = maxsnplen
    cols0 = np.zeros(nrows, dtype=np.int64) #h5py.special_dtype(vlen=bytes))
    cols1 = np.zeros(nrows, dtype=np.uint32)
    cols34 = np.zeros((nrows, 2), dtype="S5")
    cols7 = np.zeros((nrows, 1), dtype="S20")

    ## when nsamples is high this blows up memory (e.g., dim=(5M x 500))
    ## so we'll instead create a list of arrays with 10 samples at a time.
    ## maybe later replace this with a h5 array
    tmph = os.path.join(data.dirs.outfiles, ".tmp.{}.h5".format(hslice[0]))
    htmp = h5py.File(tmph, 'w')
    htmp.create_dataset("vcf", shape=(nrows, sum(sidx)), dtype="S20")

    ## which loci passed all filters
    init = 0

    ## write loci that passed after trimming edges, then write snp string
    locindex = np.where(keepmask)[0]
    for iloc in xrange(aseqs.shape[0]):
        edg = aedge[iloc]
        ## grab all seqs between edges
        if not 'pair' in data.paramsdict["datatype"]:
            seq = aseqs[iloc, :, edg[0]:edg[1]+1]
            catg = acatg[iloc, :, edg[0]:edg[1]+1]
            if not full:
                snpidx = snpidxs[iloc, edg[0]:edg[1]+1]
                seq = seq[:, snpidx]
                catg = catg[:, snpidx]
        else:
            seq = np.hstack([aseqs[iloc, :, edg[0]:edg[1]+1],
                             aseqs[iloc, :, edg[2]:edg[3]+1]])
            catg = np.hstack([acatg[iloc, :, edg[0]:edg[1]+1],
                              acatg[iloc, :, edg[2]:edg[3]+1]])
            if not full:
                snpidx = np.hstack([snpidxs[iloc, edg[0]:edg[1]+1],
                                    snpidxs[iloc, edg[2]:edg[3]+1]])
                seq = seq[:, snpidx]
                catg = catg[:, snpidx]

        ## empty arrs to fill
        alleles = np.zeros((nrows, 4), dtype=np.uint8)
        genos = np.zeros((seq.shape[1], sum(sidx)), dtype="S4")
        genos[:] = "./.:"

        ## ----  build string array ----
        pos = 0
        ## If any < 0 this indicates an anonymous locus in denovo+ref assembly
        if achrom[iloc][0] > 0:
            pos = achrom[iloc][1]
            cols0[init:init+seq.shape[1]] = achrom[iloc][0]
            cols1[init:init+seq.shape[1]] = pos + np.where(snpidx)[0] + 1
        else:
            if full:
                cols1[init:init+seq.shape[1]] = pos + np.arange(seq.shape[1]) + 1
            else:
                cols1[init:init+seq.shape[1]] = pos + np.where(snpidx)[0] + 1
                cols0[init:init+seq.shape[1]] = (chunk + locindex[iloc] + 1) * -1

        ## fill reference base
        alleles = reftrick(seq, GETCONS)

        ## get the info string column
        tmp0 = np.sum(catg, axis=2)
        tmp1 = tmp0 != 0
        tmp2 = tmp1.sum(axis=1) > 0
        nsamp = np.sum(tmp1, axis=0)
        depth = np.sum(tmp0, axis=0)
        list7 = [["NS={};DP={}".format(i, j)] for i, j in zip(nsamp, depth)]
        if list7:
            cols7[init:init+seq.shape[1]] = list7

        ## default fill cons sites where no variants
        genos[tmp1.T] = "0/0:"

        ## fill cons genotypes for sites with alt alleles for taxa in order
        mask = alleles[:, 1] == 46
        mask += alleles[:, 1] == 45

        obs = alleles[~mask, :]
        alts = seq[:, ~mask]
        who = np.where(mask == False)[0]
        ## fill variable sites
        for site in xrange(alts.shape[1]):
            bases = alts[:, site]
            #LOGGER.info("bases %s", bases)
            ohere = obs[site][obs[site] != 0]
            #LOGGER.info("ohere %s", ohere)
            alls = np.array([DCONS[i] for i in bases], dtype=np.uint32)
            #LOGGER.info("all %s", alls)
            for jdx in xrange(ohere.shape[0]):
                alls[alls == ohere[jdx]] = jdx

            #LOGGER.info("all2 %s", alls)
            ## fill into array
            for cidx in xrange(catg.shape[0]):
                if tmp2[cidx]:
                    if alls[cidx][0] < 5:
                        genos[who[site], cidx] = "/".join(alls[cidx].astype("S1").tolist())+":"
                    else:
                        genos[who[site], cidx] = "./.:"
                    #LOGGER.info("genos filled: %s %s %s", who[site], cidx, genos)

        ## build geno+depth strings
        ## for each taxon enter 4 catg values
        fulltmp = np.zeros((seq.shape[1], catg.shape[0]), dtype="S20")
        for cidx in xrange(catg.shape[0]):
            ## fill catgs from catgs
            tmp0 = [str(i.sum()) for i in catg[cidx]]
            tmp1 = [",".join(i) for i in catg[cidx].astype("S4").tolist()]
            tmp2 = ["".join(i+j+":"+k) for i, j, k in zip(genos[:, cidx], tmp0, tmp1)]
            ## fill tmp allcidx
            fulltmp[:, cidx] = tmp2

        ## write to h5 for this locus
        htmp["vcf"][init:init+seq.shape[1], :] = fulltmp

        cols34[init:init+seq.shape[1], 0] = alleles[:, 0].view("S1")
        cols34[init:init+seq.shape[1], 1] = [",".join([j for j in i if j]) \
                                    for i in alleles[:, 1:].view("S1").tolist()]

        ## advance counter
        init += seq.shape[1]

    ## trim off empty rows if they exist
    withdat = cols0 != 0
    tot = withdat.sum()

    ## get scaffold names
    faidict = {}
    if (data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]) and \
       (os.path.exists(data.paramsdict["reference_sequence"])):
        fai = pd.read_csv(data.paramsdict["reference_sequence"] + ".fai", 
                names=['scaffold', 'size', 'sumsize', 'a', 'b'],
                sep="\t")
        faidict = {i+1:j for i,j in enumerate(fai.scaffold)}
    try:
        ## This is hax, but it's the only way it will work. The faidict uses positive numbers
        ## for reference sequence mapped loci for the CHROM/POS info, and it uses negative
        ## numbers for anonymous loci. Both are 1 indexed, which is where that last `+ 2` comes from.
        faidict.update({-i:"locus_{}".format(i-1) for i in xrange(chunk+1, chunk + optim + 2)})
        chroms = [faidict[i] for i in cols0]
    except Exception as inst:
        LOGGER.error("Invalid chromosome dictionary indexwat: {}".format(inst))
        LOGGER.debug("faidict {}".format([str(k)+"/"+str(v) for k, v in faidict.items() if "locus" in v]))
        LOGGER.debug("chroms {}".format([x for x in cols0 if x < 0]))
        raise
    cols0 = np.array(chroms)
    #else:
    #    cols0 = np.array(["locus_{}".format(i) for i in cols0-1])

    ## Only write if there is some data that passed filtering
    if tot:
        LOGGER.debug("Writing data to vcf")
        if not full:
            writer = open(data.outfiles.vcf+".{}".format(chunk), 'w')
        else:
            writer = gzip.open(data.outfiles.vcf+".{}".format(chunk), 'w')

        try:
            ## write in iterations b/c it can be freakin huge.
            ## for cols0 and cols1 the 'newaxis' slice and the transpose
            ## are for turning the 1d arrays into column vectors.
            np.savetxt(writer,
                np.concatenate(
                   (cols0[:tot][np.newaxis].T,
                    cols1[:tot][np.newaxis].T,
                    np.array([["."]]*tot, dtype="S1"),
                    cols34[:tot, :],
                    np.array([["13", "PASS"]]*tot, dtype="S4"),
                    cols7[:tot, :],
                    np.array([["GT:DP:CATG"]]*tot, dtype="S10"),
                    htmp["vcf"][:tot, :],
                    ),
                    axis=1),
                delimiter="\t", fmt="%s")
        except Exception as inst:
            LOGGER.error("Error building vcf file - ".format(inst))
            raise
        writer.close()

    ## close and remove tmp h5
    htmp.close()
    os.remove(tmph)



@numba.jit(nopython=True)
def reftrick(iseq, consdict):
    """ Returns the most common base at each site in order. """

    altrefs = np.zeros((iseq.shape[1], 4), dtype=np.uint8)
    altrefs[:, 1] = 46

    for col in xrange(iseq.shape[1]):
        ## expand colums with ambigs and remove N-
        fcounts = np.zeros(90, dtype=np.int64)
        counts = np.bincount(iseq[:, col])#, minlength=90)
        fcounts[:counts.shape[0]] = counts
        ## set N and - to zero, wish numba supported minlen arg
        fcounts[78] = 0
        fcounts[45] = 0
        ## add ambig counts to true bases
        for aidx in xrange(consdict.shape[0]):
            nbases = fcounts[consdict[aidx, 0]]
            for _ in xrange(nbases):
                fcounts[consdict[aidx, 1]] += 1
                fcounts[consdict[aidx, 2]] += 1
            fcounts[consdict[aidx, 0]] = 0

        ## now get counts from the modified counts arr
        who = np.argmax(fcounts)
        altrefs[col, 0] = who
        fcounts[who] = 0

        ## if an alt allele fill over the "." placeholder
        who = np.argmax(fcounts)
        if who:
            altrefs[col, 1] = who
            fcounts[who] = 0

            ## if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 2] = who
            fcounts[who] = 0

            ## if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 3] = who

    return altrefs



GETCONS = np.array([[82, 71, 65],
                    [75, 71, 84],
                    [83, 71, 67],
                    [89, 84, 67],
                    [87, 84, 65],
                    [77, 67, 65]], dtype=np.uint8)

## with N and - masked to 9
GETCONS2 = np.array([[82, 71, 65],
                     [75, 71, 84],
                     [83, 71, 67],
                     [89, 84, 67],
                     [87, 84, 65],
                     [77, 67, 65],
                     [78, 9, 9],
                     [45, 9, 9]], dtype=np.uint8)



DCONS = {67: [67, 67],
         65: [65, 65],
         84: [84, 84],
         71: [71, 71],
         82: [71, 65],
         75: [71, 84],
         83: [71, 67],
         89: [84, 67],
         87: [84, 65],
         77: [67, 65],
         78: [46, 46],
         45: [46, 46]}

# GETCONS = np.array([["C", "C", "C"],
#                     ["A", "A", "A"],
#                     ["T", "T", "T"],
#                     ["G", "G", "G"],
#                     ["R", "G", "A"],
#                     ["K", "G", "T"],
#                     ["S", "G", "C"],
#                     ["Y", "T", "C"],
#                     ["W", "T", "A"],
#                     ["M", "C", "A"],
#                     ["N", ".", "."],
#                     ["-", ".", "."]])


def viewgeno(site, reso):
    """
    Get resolution. Only 1 or 0 allowed. Used for geno
    """
    return VIEW[site][reso]


## vectorize the viewgeno func.
## basically just makes it run as a for loop
vecviewgeno = np.vectorize(viewgeno)



def vcfheader(data, names, ofile):
    """
    Prints header for vcf files
    """
    ## choose reference string
    if data.paramsdict["reference_sequence"]:
        reference = data.paramsdict["reference_sequence"]
    else:
        reference = "pseudo-reference (most common base at site)"


    ##FILTER=<ID=minCov,Description="Data shared across <{mincov} samples">
    ##FILTER=<ID=maxSH,Description="Heterozygosous site shared across >{maxsh} samples">
    header = """\
##fileformat=VCFv4.0
##fileDate={date}
##source=ipyrad_v.{version}
##reference={reference}
##phasing=unphased
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{names}
""".format(date=time.strftime("%Y/%m/%d"),
           version=__version__,
           reference=os.path.basename(reference),
           mincov=data.paramsdict["min_samples_locus"],
           maxsh=data.paramsdict["max_shared_Hs_locus"],
           names="\t".join(names))
    ## WRITE
    ofile.write(header)


NEXHEADER = \
"""#nexus
begin data;
  dimensions ntax={} nchar={};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""

NEXCLOSER = """\
  ;
end;
"""



if __name__ == "__main__":
    import ipyrad as ip

    ## get path to test dir/
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__)
               )
           )
       )

    ## run test on pairgbs data1
    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "test_rad", "data1"))
    TEST.step7(force=True)
    print(TEST.stats)
