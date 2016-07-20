#!/usr/bin/env python2.7

""" Apply filters and write output files. The basic body plan of this
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
import itertools
import datetime
import shutil
import numba
import time
import glob
import gzip
import h5py
import os
from collections import Counter, OrderedDict
#from ipyrad.file_conversion import *
from ipyrad import __version__
from util import *

import logging
LOGGER = logging.getLogger(__name__)


## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized. LOCI and VCF are default. Some others are created as 
## dependencies of others.

#OUTPUT_FORMATS = ['alleles', 'phy', 'nex', 'snps', 'usnps', 'vcf',
#                  'str', 'geno', 'treemix', 'migrate', 'gphocs']

OUTPUT_FORMATS = ['phy', 'nex', 'snps', 'usnps', 'str', 'geno', 'vcf']



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
    make_loci_and_stats(data, samples, ipyclient)

    ## OPTIONAL OUTPUTS -- grab from params as a string, if commas then split 
    ## to a list, e.g., [phy, nex, str]. 
    output_formats = data.paramsdict["output_formats"]
    if "," in output_formats:
        output_formats = [i.strip() for i in output_formats.split(",")]
    if "*" in output_formats:
        output_formats = OUTPUT_FORMATS

    ## held separate from *output_formats cuz it's big and parallelized 
    if 'vcf' in output_formats:
        make_vcf(data, samples, ipyclient)

    ## make other array-based formats, recalcs keeps and arrays
    make_outfiles(data, samples, output_formats, ipyclient)

    ## print friendly message
    shortpath = data.dirs.outfiles.replace(os.path.expanduser("~"), "~")
    print("  Outfiles written to: {}".format(shortpath))



def make_stats(data, samples, samplecounts, locuscounts):
    """ write the output stats file and save to Assembly obj."""
    ## load the h5 database
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')    
    ## get meta info 
    anames = io5["seqs"].attrs["samples"]
    nloci = io5["seqs"].shape[0]
    optim = io5["seqs"].attrs["chunksize"]

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
    #snpcounts = Counter()
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

    #bisnps = Counter()
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
          "\n## ipyrad API location: [assembly].statsfiles.s7_filters\n", 
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
    smax = max([i+1 for i in varcounts if varcounts[i]])

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

    print("\n\n\n## The distribution of SNPs (var and pis) across loci."+\
          "\n## var = all variable sites (pis + autapomorphies)"+\
          "\n## pis = parsimony informative site (minor allele in >1 sample)"+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_snps\n",
          file=outstats)
    data.stats_dfs.s7_snps = pd.concat([vardat, varsums, pisdat, pissums], 
                                        axis=1)    
    data.stats_dfs.s7_snps.to_string(buf=outstats)

    ## close it
    outstats.close()
    io5.close()
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
    Open the HDF5 array with seqs, catg, and filter data. Fill the remaining
    filters. 
    """
    ## create loadbalanced ipyclient
    lbview = ipyclient.load_balanced_view()

    ## get chunk size from the HD5 array and close
    io5 = h5py.File(data.clust_database, 'r')
    ## the size of chunks for reading/writing
    optim = io5["seqs"].attrs["chunksize"]
    ## the samples in the database in their locus order
    dbsamples = io5["seqs"].attrs["samples"]
    ## the total number of loci        
    nloci = io5["seqs"].shape[0]
    io5.close()

    ## make a tmp directory for saving chunked arrays to
    chunkdir = os.path.join(data.dirs.outfiles, data.name+"_tmpchunks")
    if not os.path.exists(chunkdir):
        os.mkdir(chunkdir)

    ## get the indices of the samples that we are going to include
    sidx = select_samples(dbsamples, samples)
    LOGGER.info("samples %s \n, dbsamples %s \n, sidx %s \n", 
        samples, dbsamples, sidx)

    ## Put inside a try statement so we can delete tmpchunks 
    try:
        ## load a list of args to send to Engines. Each arg contains the index
        ## to sample optim loci from catg, seqs, filters &or edges, which will
        ## be loaded on the remote Engine. 

        ## create job queue
        start = time.time()
        results = [] 
        submitted = 0
        while submitted < nloci:
            hslice = np.array([submitted, submitted+optim])
            async = lbview.apply(filter_stacks, [data, sidx, hslice])
            results.append(async)
            submitted += optim

        ## run filter_stacks on all chunks
        while 1:
            readies = [i.ready() for i in results]
            if not all(readies):
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(len(readies), sum(readies), 
                    " filtering loci        | {}".format(elapsed))
                time.sleep(1)
            else:
                break
        progressbar(20, 20, " filtering loci        | {}".format(elapsed))        
        print("")

        ## make sure that all are done
        [i.get() for i in results]

        ## get all the saved tmp arrays for each slice
        tmpsnp = glob.glob(os.path.join(chunkdir, "snpf.*.npy"))
        tmphet = glob.glob(os.path.join(chunkdir, "hetf.*.npy"))
        tmpmin = glob.glob(os.path.join(chunkdir, "minf.*.npy"))
        tmpedg = glob.glob(os.path.join(chunkdir, "edgf.*.npy"))
        tmppld = glob.glob(os.path.join(chunkdir, "pldf.*.npy"))        

        ## sort array files within each group
        arrdict = OrderedDict([
            ('snp', tmpsnp), ('het', tmphet),
            ('min', tmpmin), ('edg', tmpedg), 
            ('pld', tmppld)])
        for arrglob in arrdict.values():
            arrglob.sort(key=lambda x: int(x.rsplit(".")[-2]))

        ## re-load the full filter array who's order is
        ## ["duplicates", "max_indels", "max_snps", "max_hets", "min_samps", "max_alleles"]
        io5 = h5py.File(data.database, 'r+')
        superfilter = np.zeros(io5["filters"].shape, io5["filters"].dtype)

        ## iterate across filter types (dups & indels are already filled)
        ## we have [4,4] b/c minf and edgf both write to minf
        for fidx, ftype in zip([2, 3, 4, 4, 5], arrdict.keys()):
            ## fill in the edgefilters
            for ffile in arrdict[ftype]:
                ## grab a file and get it's slice            
                hslice = int(ffile.split(".")[-2])
                ## load in the array
                arr = np.load(ffile)
                ## store slice into full array (we use += here because the minf
                ## and edgf arrays both write to the same filter). 
                superfilter[hslice:hslice+optim, fidx] += arr
        io5["filters"][:] += superfilter
        del arr

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
            LOGGER.info("arr is %s", arr.shape)
            LOGGER.info("superedge is %s", superedge.shape)            
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
            supersnps[hslice:hslice+optim, :, :] = arr
        io5["snps"][:] = supersnps            
        del arr

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
    time and write to file. 
    """
    ## start vcf progress bar
    start = time.time()
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, 
        " building loci/stats   | {}".format(elapsed))

    ## get some db info
    io5 = h5py.File(data.clust_database, 'r')
    ## will iterate optim loci at a time
    optim = io5["seqs"].attrs["chunksize"]
    nloci = io5["seqs"].shape[0]
    anames = io5["seqs"].attrs["samples"]

    ## get name and snp padding
    pnames, snppad = padnames(anames)
    snames = [i.name for i in samples]
    smask = np.array([i not in snames for i in anames])

    ## keep track of how many loci from each sample pass all filters
    samplecov = np.zeros(len(anames), dtype=np.int32)
    locuscov = Counter()

    ## set initial value to zero for all values above min_samples_locus
    #for cov in range(data.paramsdict["min_samples_locus"], len(anames)+1):
    for cov in range(len(anames)+1):        
        locuscov[cov] = 0

    ## client for sending jobs to parallel engines
    lbview = ipyclient.load_balanced_view()

    ## send jobs in chunks
    loci_asyncs = []
    for istart in xrange(0, nloci, optim):
        args = [data, optim, pnames, snppad, smask, istart, samplecov, locuscov]
        loci_asyncs.append(lbview.apply(locichunk, args))

    ### FOR DEBUGGING
    # ipyclient.wait()
    # for job in loci_asyncs:
    #     if not job.successful():
    #         print(job.metadata)

    ## just a waiting function for chunks to finish
    tmpids = list(itertools.chain(*[i.msg_ids for i in loci_asyncs]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1)

    while 1:
        done = [i.ready() for i in loci_asyncs]
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        progressbar(len(done), sum(done),
            " building loci/stats   | {}".format(elapsed))
        if not res.completed:
            time.sleep(1)
        else:
            print("")   
            break            

    ## check for errors
    for job in loci_asyncs:
        if job.ready() and not job.successful():
            print(job.metadata)

    ## concat and cleanup
    results = [i.get() for i in loci_asyncs]
    ## update dictionaries
    for chunk in results:
        samplecov += chunk[0]
        locuscov.update(chunk[1])

    ## get all chunk files
    tmploci = glob.glob(data.outfiles.loci+".[0-9]*")
    ## sort by start value
    tmploci.sort(key=lambda x: int(x.split(".")[-1]))
    ## write tmpchunks to locus file
    with open(data.outfiles.loci, 'w') as locifile:
        for tmploc in tmploci:
            with open(tmploc, 'r') as inloc:
                locifile.write(inloc.read())
            os.remove(tmploc)

    ## make stats file from data
    make_stats(data, samples, samplecov, locuscov)
    io5.close()



def locichunk(args):
    """ 
    Function from make_loci to apply to chunks. smask is sample mask. 
    """
    ## parse args
    data, optim, pnames, snppad, smask, start, samplecov, locuscov = args

    ## this slice 
    hslice = [start, start+optim]

    ## get filter db info
    co5 = h5py.File(data.database, 'r')
    afilt = co5["filters"][hslice[0]:hslice[1], ]
    aedge = co5["edges"][hslice[0]:hslice[1], ]
    asnps = co5["snps"][hslice[0]:hslice[1], ]

    ## get seqs db 
    io5 = h5py.File(data.clust_database, 'r')
    aseqs = io5["seqs"][hslice[0]:hslice[1], ]

    ## which loci passed all filters
    #LOGGER.info("edgfilt is %s", aedge.sum())
    keep = np.where(np.sum(afilt, axis=1) == 0)[0]
    #LOGGER.info("keep is %s", keep.sum())
    ## store until printing
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
    """

    ## get stats from step6 h5 and create new h5
    co5 = h5py.File(data.clust_database, 'r')
    io5 = h5py.File(data.database, 'w')

    ## get maxlen and chunk len
    maxlen = data._hackersonly["max_fragment_length"] + 30
    chunks = co5["seqs"].attrs["chunksize"]
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
    
    ## apply indel filter 
    if "pair" in data.paramsdict["datatype"]:
        ## get lenght of locus 
        maxinds = sum(data.paramsdict["max_Indels_locus"])
    else:
        maxinds = data.paramsdict["max_Indels_locus"][0]
    filters[:, 1] = co5["indels"][:] > maxinds

    io5.close()
    co5.close()



def filter_stacks(args):
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
    data, sidx, hslice = args

    ## open h5 handles
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r+')
    ## get a chunk (hslice) of loci for the selected samples (sidx)
    superseqs = io5["seqs"][hslice[0]:hslice[1], sidx,]

    ## duplicate filter is already filled
    ## indels filter is already filled

    ## get an int view of the seq array
    superints = superseqs.view(np.int8)

    ## fill edge filter
    ## get edges of superseqs and supercats, since edges need to be trimmed 
    ## before counting hets and snps. Technically, this could edge trim 
    ## clusters to the point that they are below the minlen, and so this 
    ## also constitutes a filter, though one that is uncommon. For this 
    ## reason we have another filter called edgfilter.
    splits = co5["edges"][hslice[0]:hslice[1], 4]
    edgfilter, edgearr = get_edges(data, superints, splits)
    del splits
    LOGGER.info('passed edges %s', hslice[0])

    ## minsamp coverages filtered from superseqs
    minfilter = filter_minsamp(data.paramsdict["min_samples_locus"], superints)
    LOGGER.info('passed minfilt %s', hslice[0])

    ## maxhets per site column from superseqs after trimming edges
    hetfilter = filter_maxhet(data, superints, edgearr)
    LOGGER.info('passed minhet %s', hslice[0])

    ## ploidy filter
    pldfilter = io5["nalleles"][hslice[0]:hslice[1]].max(axis=1) > \
                                         data.paramsdict["max_alleles_consens"]

    ## Build the .loci snpstring as an array (snps) 
    ## shape = (chunk, 1) dtype=S1, or should it be (chunk, 2) for [-,*] ?
    snpfilter, snpsarr = filter_maxsnp(data, superints, edgearr)

    LOGGER.info("edg %s", edgfilter.sum())
    LOGGER.info("min %s", minfilter.sum())
    LOGGER.info("het %s", hetfilter.sum())
    LOGGER.info("pld %s", pldfilter.sum())
    LOGGER.info("snp %s", snpfilter.sum())

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
    site if it present.
    """
    ## the filtering arg and parse it into minsamp numbers
    edgetrims = np.array(data.paramsdict["trim_overhang"])

    ## Cuts 3 and 4 are only for 3rad/radcap
    ## TODO: This is moderately hackish, it's not using cut3/4
    ## correctly, just assuming the length is the same as cut1/2
    if "3rad" in data.paramsdict["datatype"]:
        cut1, cut2, _, _ = data.paramsdict["restriction_overhang"] 
    else:
        cut1, cut2 = data.paramsdict["restriction_overhang"]
    cuts = np.int16([len(cut1), len(cut2)])

    ## a local array for storing edge trims
    edges = np.zeros((superints.shape[0], 5), dtype=np.int16)

    ## a local array for storing edge filtered loci, these are stored 
    ## eventually as minsamp excludes.
    edgefilter = np.zeros((superints.shape[0],), dtype=np.bool)

    ## TRIM GUIDE. The cut site is always trimmed. Number if what else gets it.
    ## TODO: all of these are not yet fully implemented. 
    ## If negative, then N bases are trimmed from that edge.
    ## If zero: no additional trimming.
    ## If positive: N is the minimum number of samples with data to trim to.
    ## default is (4, *, *, 4)
    ## A * value means no missing data (overhang) at edges.
    #edgetrims[edgetrims == "*"] = minsamp
    edgetrims = edgetrims.astype(np.int16)

    ## convert all - to N to make this easier, (only affects local copy)
    #superseqs[superseqs == "-"] = "N"
    superints[superints == 45] = 78

    ## trim overhanging edges
    ## get the number not Ns in each site, 
    #ccx = np.sum(superseqs != "N", axis=1)
    ccx = np.sum(superints != 78, axis=1, dtype=np.int16)
    efi, edg = edgetrim_numba(splits, ccx, edges, edgefilter, edgetrims, cuts)
    return efi, edg


@numba.jit(nopython=True)
def edgetrim_numba(splits, ccx, edges, edgefilter, edgetrims, cuts):
    ## get splits
    for idx in xrange(splits.shape[0]):
        
        ## get the data
        if splits[idx]:
            r1s = ccx[idx, :splits[idx]]
            r2s = ccx[idx, splits[idx]+4:]
        else:
            r1s = ccx[idx, :]
            
        ## fill in edge 0
        if edgetrims[0] > 0:
            x = np.where(r1s >= edgetrims[0])[0]
            if np.any(x):
                edges[idx][0] = np.min(x[cuts[0]:])
            else:
                edges[idx][0] = np.uint16(0)
                edgefilter[idx] = True
            
        elif edgetrims[0] < 0:
            x = np.where(r1s >= 0)[0]
            if np.any(x):
                edges[idx][0] = np.max(x[cuts[0]:]) - edgetrims[0]
            else:
                edges[idx][0] = np.uint16(0)
                edgefilter[idx] = True
                
        ## fill in edge 1
        if edgetrims[1] > 0:
            x = np.where(r1s >= edgetrims[1])[0]
            if np.any(x):
                edges[idx][1] = np.max(x)
            else:
                edges[idx][1] = np.uint16(1)
                edgefilter[idx] = True
            
        elif edgetrims[1] < 0:
            x = np.where(r1s > 0)[0]
            if np.any(x):
                edges[idx][1] = np.max(x) - edgetrims[1]
            else:
                edges[idx][1] = np.uint16(1)
                edgefilter[idx] = True       
           
        ## If paired, do second read edges    
        if splits[idx]:
            if edgetrims[2] > 0:
                x = np.where(r2s >= edgetrims[2])[0]
                if np.any(x):
                    edges[idx][2] = splits[idx] + np.int16(4) + np.min(x)
                else:
                    edges[idx][2] = edges[idx][1] + 4
                    edgefilter[idx] = True
            elif edgetrims[2] < 0:
                x = np.where(r2s >= 0)[0]
                if np.any(x):
                    edges[idx][2] = splits[idx] + np.int16(4) + np.max(x) + edgetrims[2]
                else:
                    edges[idx][2] = edges[idx][1] + 4
                    edgefilter[idx] = True
            
            if edgetrims[3] > 0:
                ## get farthest site that is not all Ns
                x = np.where(r2s > 0)[0]
                if np.any(x):
                    ## remove the cut site
                    cutx = x[:np.max(x) - cuts[1]]
                    ## find farthest with high cov
                    covmask = r2s[:np.max(cutx)+np.int16(1)] >= edgetrims[3]
                    ## return index w/ spacers
                    edges[idx][3] = splits[idx] + np.int16(4) + np.max(np.where(covmask)[0])
                else:
                    edges[idx][3] = edges[idx][2] + 1
                    edgefilter[idx] = True   
            elif edgetrims[3] < 0:
                edges[idx][3] = r2s[edges[idx][3]] + edgetrims[3]
                
            ## enter the pair splitter
            edges[idx][4] = splits[idx]
            
    return edgefilter, edges
        


def filter_minsamp(minsamp, superints):
    """ 
    Filter minimum # of samples per locus from superseqs[chunk]. The shape
    of superseqs is [chunk, sum(sidx), maxlen]
    """
    ## the minimum filter
    #minsamp = data.paramsdict["min_samples_locus"]
    ## ask which rows are not all N along seq dimension, then sum along sample 
    ## dimension to get the number of samples that are not all Ns.
    ## minfilt is a boolean array where True means it failed the filter.
    #minfilt = np.sum(~np.all(superseqs == "N", axis=2), axis=1) < minsamp
    minfilt = np.sum(~np.all(superints == 78, axis=2), axis=1) < minsamp
    #LOGGER.debug("minfilt %s", minfilt)

    ## print the info
    #LOGGER.info("Filtered by min_samples_locus - {}".format(minfilt.sum()))
    return minfilt



def fakeref(sitecol, idx=0):
    """
    Used to find the most frequent base at each column for making a 
    stand-in reference sequence for denovo loci that have no reference. 
    Simply a way for representing the results in geno and VCF outputs.
    If idx, returns others alleles that are non-REF. 
    """
    ## a list for only catgs
    catg = [i for i in sitecol if i in "CATG"]

    ## find sites that are ambigs
    where = [sitecol[sitecol == i] for i in "RSKYWM"]

    ## for each occurrence of RSKWYM add ambig resolution to catg
    for ambig in where:
        for _ in range(ambig.size):
            catg += list(AMBIGS[ambig[0]])

    ## return the most common element, if tied, it doesn't matter who wins
    setg = set(catg)
    base = max(setg or ["."], key=catg.count)
    if not idx:
        return base
    else:
        setg.discard(base)
        base = max(setg or ["."], key=catg.count)
        if idx == 1:
            return base
        else:
            setg.discard(base)
            base = max(setg or [""], key=catg.count)
            if idx == 2:
                return base
            else:
                setg.discard(base)
                base = max(setg or [""], key=catg.count)
                return base



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
                ## if one var, e.g., [0, 0, 1, 8] then (-), else (*)
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
        maxhet = np.array(superints.shape[1]*maxhet, dtype=np.int16)
    elif isinstance(maxhet, int):
        maxhet = np.array(maxhet, dtype=np.int16)

    ## an empty array to fill with failed loci
    hetfilt = np.zeros(superints.shape[0], dtype=np.bool)
    hetfilt = maxhet_numba(superints, edgearr, maxhet, hetfilt)
    LOGGER.info("--------------maxhet sums %s", hetfilt.sum())
    return hetfilt



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
        if (count1s.max() > maxhet) or (count2s.max() > maxhet):
            hetfilt[idx] = True 
    return hetfilt    

## MAKE GLOBAL
AMBIGARR = np.array(list("RSKYWM")).view(np.int8)




def make_outfiles(data, samples, output_formats, ipyclient):
    """
    Get desired formats from paramsdict and write files to outfiles 
    directory.
    """

    ## load the h5 database
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')    

    ## will iterate optim loci at a time
    optim = io5["seqs"].attrs["chunksize"]
    nloci = io5["seqs"].shape[0]

    ## get name and snp padding
    anames = io5["seqs"].attrs["samples"]
    snames = [i.name for i in samples]
    names = [i for i in anames if i in snames]
    pnames, _ = padnames(names)
    pnames.sort()

    ## get names boolean
    sidx = np.array([i in snames for i in anames])
    assert len(pnames) == sum(sidx)

    ## get names index in order of pnames
    #sindx = [list(anames).index(i) for i in snames]

    ## build arrays and outputs from arrays. 
    ## TODO: make faster with numba, dask, and/or parallel
    arrs = make_arrays(data, sidx, optim, nloci, io5, co5)
    seqarr, snparr, bisarr, maparr = arrs

    ## send off outputs as parallel jobs
    lbview = ipyclient.load_balanced_view()
    start = time.time()
    results = []

    ## phy and partitions are a default output ({}.phy, {}.phy.partitions)
    if "phy" in output_formats:
        async = lbview.apply(write_phy, *[data, seqarr, sidx, pnames])
        results.append(async)
        #write_phy(data, seqarr, sidx, pnames)

    ## other outputs
    ## nexus format includes ... additional information ({}.nex)
    if 'nexus' in output_formats:
        async = lbview.apply(write_nex, *[data, seqarr, sidx, pnames])
        results.append(async)
        #write_nex(data, seqarr, sidx, pnames)

    ## snps is actually all snps written in phylip format ({}.snps.phy)
    if 'snps' in output_formats:
        async = lbview.apply(write_snps, *[data, snparr, sidx, pnames])
        results.append(async)
        write_snps_map(data, maparr)
        #write_snps(data, snparr, sidx, pnames)

    ## usnps is one randomly sampled snp from each locus ({}.u.snps.phy)
    if 'usnps' in output_formats:
        async = lbview.apply(write_usnps, *[data, bisarr, sidx, pnames])
        results.append(async)
        #write_usnps(data, bisarr, sidx, pnames)

    ## str and ustr are for structure analyses. A fairly outdated format, six
    ## columns of empty space. Full and subsample included ({}.str, {}.u.str)
    if 'str' in output_formats:
        async = lbview.apply(write_str, *[data, snparr, bisarr, sidx, pnames])
        results.append(async)
        #write_str(data, snparr, bisarr, sidx, pnames)

    ## geno output is for admixture and other software. We include all SNPs, 
    ## but also a .map file which has "distances" between SNPs.
    ## ({}.geno, {}.map)
    if 'geno' in output_formats:
        async = lbview.apply(write_geno, *[data, snparr, bisarr, sidx, pnames])
        results.append(async)
        #write_geno(data, snparr, bisarr, sidx, pnames)

    ## run filter_stacks on all chunks
    while 1:
        readies = [i.ready() for i in results]
        if not all(readies):
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(len(readies), sum(readies), 
                " writing outfiles      | {}".format(elapsed))
            time.sleep(1)
        else:
            break

    ## final progress bar
    elapsed = datetime.timedelta(seconds=int(time.time()-start))            
    progressbar(20, 20, " writing outfiles      | {}".format(elapsed))        
    #if data._headers:
    print("")

    ## check for errors
    for async in results:
        if not async.completed:
            print(async.metadata.error)

    ## close h5 handle
    io5.close()






## TODO: use view(np.uint8) and compile loop as numba func
def make_arrays(data, sidx, optim, nloci, io5, co5):
    """ 
    Builds arrays for seq, snp, and bis data after applying locus filters, 
    and edge filteres for the seq data. These arrays are used to build outs. 
    """

    ## make empty arrays for filling
    maxlen = data._hackersonly["max_fragment_length"] + 30
    maxsnp = co5["snps"][:].sum()

    ## shape of arrays is sidx, we will subsample h5 w/ sidx to match
    seqarr = np.zeros((sum(sidx), maxlen*nloci), dtype="S1")
    snparr = np.zeros((sum(sidx), maxsnp), dtype="S1")
    bisarr = np.zeros((sum(sidx), nloci), dtype="S1")
    maparr = np.zeros((maxsnp, 4), dtype=np.uint32)

    ## apply all filters and write loci data
    start = 0
    seqleft = 0
    snpleft = 0
    bis = 0

    ## edge filter has already been applied to snps, but has not yet been 
    ## applied to seqs. The locus filters have not been applied to either yet.
    mapsnp = 0
    totloc = 0
    while start < nloci:
        hslice = [start, start+optim]
        afilt = co5["filters"][hslice[0]:hslice[1], ...]
        aedge = co5["edges"][hslice[0]:hslice[1], ...]
        asnps = co5["snps"][hslice[0]:hslice[1], ...]
        aseqs = io5["seqs"][hslice[0]:hslice[1], sidx, ...]

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]
        #LOGGER.info("edges.shape %s", aedge.shape)

        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            ## grab r1 seqs between edges
            edg = aedge[iloc]
            #LOGGER.info("edg %s", edg)

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
            #LOGGER.info("seq %s", seq.shape)

            ## remove cols from seq (phy) array that are all N-
            lcopy = seq
            lcopy[lcopy == "-"] = "N"
            bcols = np.all(lcopy == "N", axis=0)
            seq = seq[:, ~bcols]

            ## put into large array
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
                maparr[mapsnp, :] = [totloc+1, start+iloc, i, mapsnp+1]
                mapsnp += 1

            #LOGGER.info("%s -- %s", seqleft, seqleft+seq.shape[1])
            #LOGGER.info("%s -- %s", snps, snps.shape)

            ## subsample one SNP into an array
            if snps.shape[1]:
                samp = np.random.randint(snps.shape[1])
                bisarr[:, bis] = snps[:, samp]
                bis += 1
                totloc += 1
                
        ## increase the counter
        start += optim

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



def write_phy(data, seqarr, sidx, pnames):
    """ write the phylip output file"""
    data.outfiles.phy = os.path.join(data.dirs.outfiles, data.name+".phy")
    with open(data.outfiles.phy, 'w') as out:
        ## write header
        out.write("{} {}\n".format(seqarr.shape[0], seqarr.shape[1]))
        ## write data rows
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(seqarr[sidx][idx, :])))


## TODO:
def write_nex(data, seqarr, sidx, pnames):
    """ not implemented yet """
    return 0
    ## write the phylip string
    data.outfiles.nex = os.path.join(data.dirs.outfiles, data.name+".nex")
    with open(data.outfiles.nex, 'w') as out:
        ## trim down to size
        #ridx = np.all(seqarr == "", axis=0)
        out.write("{} {}\n".format(seqarr.shape[0], seqarr.shape[1]))
                                   #seqarr[:, ~ridx].shape[1]))
        for idx, name in zip(sidx, pnames):
            out.write("{}{}\n".format(name, "".join(seqarr[idx])))


## TODO: this could have much more information for reference aligned data
def write_snps_map(data, maparr):
    """ write a map file with linkage information for SNPs file"""
    data.outfiles.snpsmap = os.path.join(data.dirs.outfiles, data.name+".snps.map")
    with open(data.outfiles.snpsmap, 'w') as out:
        for idx in range(maparr.shape[0]):
            line = maparr[idx, :]
            out.write("{}\trad{}_snp{}\t{}\t{}\n".format(line[0], line[1], line[2], 0, line[3]))


def write_snps(data, snparr, sidx, pnames):
    """ write the snp string """
    data.outfiles.snps = os.path.join(data.dirs.outfiles, data.name+".snps.phy")    
    with open(data.outfiles.snps, 'w') as out:
        out.write("{} {}\n".format(snparr.shape[0], snparr.shape[1]))
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(snparr[sidx][idx, :])))



def write_usnps(data, bisarr, sidx, pnames):
    """ write the bisnp string """
    data.outfiles.usnps = os.path.join(data.dirs.outfiles, data.name+".u.snps.phy")
    with open(data.outfiles.usnps, 'w') as out:
        out.write("{} {}\n".format(bisarr.shape[0], bisarr.shape[1]))
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(bisarr[sidx][idx, :])))



def write_str(data, snparr, bisarr, sidx, pnames):
    """ Write STRUCTURE format """
    data.outfiles.str = os.path.join(data.dirs.outfiles, data.name+".str")
    data.outfiles.ustr = os.path.join(data.dirs.outfiles, data.name+".u.str")        
    out1 = open(data.outfiles.str, 'w')
    out2 = open(data.outfiles.ustr, 'w')
    numdict = {'A': '0', 'T': '1', 'G': '2', 'C': '3', 'N': '-9', '-': '-9'}
    if data.paramsdict["max_alleles_consens"] > 1:
        for idx, name in enumerate(pnames):
            out1.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][0]] for i in snparr[sidx][idx]])))
            out1.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][1]] for i in snparr[sidx][idx]])))
            out2.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][0]] for i in bisarr[sidx][idx]])))
            out2.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][1]] for i in bisarr[sidx][idx]])))
    else:
        ## haploid output
        for idx, name in enumerate(pnames):
            out1.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][0]] for i in snparr[sidx][idx]])))
            out2.write("{}\t\t\t\t\t{}\n"\
                .format(name,
                "\t".join([numdict[DUCT[i][0]] for i in bisarr[sidx][idx]])))
    out1.close()
    out2.close()



def write_geno(data, snparr, bisarr, sidx, inh5):
    ## Write GENO format
    data.outfiles.geno = os.path.join(data.dirs.outfiles, data.name+".geno")
    data.outfiles.ugeno = os.path.join(data.dirs.outfiles, data.name+".u.geno")

    ## get most common base at each SNP as a pseudo-reference 
    ## and record 0,1,2 or missing=9 for counts of the ref allele
    snpref = np.apply_along_axis(fakeref, 0, snparr)
    bisref = np.apply_along_axis(fakeref, 0, bisarr)

    ## geno is printed as a matrix where columns are individuals
    ## I order them by same order as in .loci, which is alphanumeric
    snpgeno = np.zeros(snparr.shape, dtype=np.uint8)
    ## put in missing
    snpgeno[snparr == "N"] = 9
    snpgeno[snparr == "-"] = 9
    ## put in complete hits
    snpgeno[snparr == snpref] = 2
    ## put in hetero hits where resolve base matches ref
    for reso in range(2):
        hets = vecviewgeno(snparr, reso)
        snpgeno[hets == snpref] = 1

    ## geno is printed as a matrix where columns are individuals
    ## I order them by same order as in .loci, which is alphanumeric
    bisgeno = np.zeros(bisarr.shape, dtype=np.uint8)
    ## put in missing
    bisgeno[bisarr == "N"] = 9
    bisgeno[bisarr == "-"] = 9
    ## put in complete hits
    bisgeno[bisarr == bisref] = 2
    ## put in hetero hits where resolve base matches ref
    for reso in range(2):
        hets = vecviewgeno(bisarr, reso)
        bisgeno[hets == bisref] = 1

    ## print to files
    np.savetxt(data.outfiles.geno, snpgeno.T, delimiter="", fmt="%d")
    np.savetxt(data.outfiles.ugeno, bisgeno.T, delimiter="", fmt="%d")

    ## write a map file for use in admixture with locations of SNPs
    ## for denovo data we just have evenly spaced the unlinked SNPs
    # 1  rs123456  0  1234555
    # 1  rs234567  0  1237793
    # 1  rs224534  0  -1237697        <-- exclude this SNP
    # 1  rs233556  0  1337456        




def make_vcf(data, samples, ipyclient):
    """ 
    Write the full VCF for loci passing filtering. Other vcf formats are 
    possible, like SNPs-only, or with filtered loci included but the filter
    explicitly labeled. These are not yet supported, however.
    """
    ## start vcf progress bar
    start = time.time()
    LOGGER.info("Writing .vcf file")
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, 
        " building vcf file     | {}".format(elapsed))

    ## create output, gzip it to be disk friendly
    data.outfiles.vcf = os.path.join(data.dirs.outfiles, data.name+".vcf.gz")

    ## get some db info
    io5 = h5py.File(data.clust_database, 'r')

    ## will iterate optim loci at a time
    optim = io5["seqs"].attrs["chunksize"]
    nloci = io5["seqs"].shape[0]

    ## get name and snp padding
    anames = io5["seqs"].attrs["samples"]
    snames = [i.name for i in samples]
    names = [i for i in anames if i in snames]

    ## get names index
    sidx = np.array([i in snames for i in anames])

    ## client for sending jobs to parallel engines
    lbview = ipyclient.load_balanced_view()

    ## send jobs in chunks
    vcf_asyncs = []
    for istart in xrange(0, nloci, optim):
        args = [data, optim, sidx, istart]
        vcf_asyncs.append(lbview.apply(vcfchunk, args))

    ## after all chunks are finished sort them, concatenate, and gzip
    tmpids = list(itertools.chain(*[i.msg_ids for i in vcf_asyncs]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(concat_vcf, *[data, names])

    while 1:
        if not res.completed:
            done = [i.ready() for i in vcf_asyncs]
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(len(done), sum(done),
                " building vcf file     | {}".format(elapsed))
            time.sleep(1)
        else:
            break

    [i.get() for i in vcf_asyncs]
    ## final progress bar
    progressbar(20, 20,
            " building vcf file     | {}".format(elapsed))
    print("")   

    io5.close()



def concat_vcf(data, names):
    """ 
    Sorts, concatenates, and gzips VCF chunks. Also cleans up chunks.
    """
    ## open handle and write headers
    vout = gzip.open(data.outfiles.vcf, 'w')
    vcfheader(data, names, vout)

    ## get vcf chunks
    vcfins = glob.glob(data.outfiles.vcf+".[0-9]*")

    ## sort
    vcfins.sort(key=lambda x: int(x.split(".")[-1]))

    ## concatenate
    for vcfin in vcfins:
        with open(vcfin, 'r') as vin:
            vout.write(vin.read())
        os.remove(vcfin)



def vcfchunk(args):
    """ 
    Function called within make_vcf to run chunks on separate engines. 
    """

    ## parse args
    data, optim, sidx, start = args

    ## load the h5 database
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')    

    ## empty array to be filled before writing 
    ## will not actually be optim*maxlen, extra needs to be trimmed
    maxlen = data._hackersonly["max_fragment_length"] + 30
    gstr = np.zeros((optim*maxlen, 9+sum(sidx)), dtype="S20")

    ## get data sliced (optim chunks at a time)
    hslice = [start, start+optim]
    afilt = co5["filters"][hslice[0]:hslice[1], ]
    aedge = co5["edges"][hslice[0]:hslice[1], ]
    aseqs = io5["seqs"][hslice[0]:hslice[1], sidx, ]
    acatg = io5["catgs"][hslice[0]:hslice[1], sidx, ]

    ## which loci passed all filters
    keep = np.where(np.sum(afilt, axis=1) == 0)[0]
    seqleft = 0

    ## write loci that passed after trimming edges, then write snp string
    for iloc in keep:
        edg = aedge[iloc]
        ## grab all seqs between edges
        seq = aseqs[iloc, :, edg[0]:edg[1]+1]
        catg = acatg[iloc, :, edg[0]:edg[1]+1]

        ## ----  build string array ---- 
        ## fill (CHR) chromosome/contig (reference) or RAD-locus (denovo)
        gstr[seqleft:seqleft+seq.shape[1], 0] = "RAD_{}_".format(start+iloc)
        ## fill (POS) position
        gstr[seqleft:seqleft+seq.shape[1], 1] = range(seq.shape[1])
        ## fill (ID) what is this? missing value is .
        gstr[seqleft:seqleft+seq.shape[1], 2] = "."
        ## get reference seq and put it in 
        seqref = np.apply_along_axis(fakeref, 0, seq)
        gstr[seqleft:seqleft+seq.shape[1], 3] = seqref

        ## get alt alleles and put them in 
        coma1 = np.zeros(seqref.shape, dtype="S1")
        coma2 = np.zeros(seqref.shape, dtype="S1")
        seqa = np.apply_along_axis(fakeref, 0, seq, *[1])
        seqb = np.char.array(np.apply_along_axis(fakeref, 0, seq, *[2]))
        seqc = np.char.array(np.apply_along_axis(fakeref, 0, seq, *[3]))
        coma1[seqb != ""] = ","
        coma2[seqc != ""] = ","
        seqr = np.char.array(seqa) + np.char.array(coma1) + \
               np.char.array(seqb) + np.char.array(coma2) + \
                   np.char.array(seqc)

        ## fill ALT allele
        gstr[seqleft:seqleft+seq.shape[1], 4] = seqr

        ## add qual [just 0.95 for now]
        gstr[seqleft:seqleft+seq.shape[1], 5] = '13'

        ## add filter (only PASS for now)
        gstr[seqleft:seqleft+seq.shape[1], 6] = "PASS"

        ## add info (just NS and DP)
        info = np.zeros(catg.shape[1], dtype="S20")
        nons = np.sum(np.sum(catg, axis=2) != 0, axis=0)            
        sums = np.sum(np.sum(catg, axis=2), axis=0)
        for i in xrange(catg.shape[1]):
            info[i] = "NS={};DP={}".format(nons[i], sums[i])
        gstr[seqleft:seqleft+seq.shape[1], 7] = info

        ## add format string
        gstr[seqleft:seqleft+seq.shape[1], 8] = "GT:CATG"

        ## ----  build counts array ----
        cnts = np.zeros((catg.shape[1], catg.shape[0]), dtype="S20")
        carr = catg.astype("S20")

        for site in xrange(cnts.shape[0]):
            obs = seqref[site]+seqr[site].replace(",", "")
            ords = {i:j for (i, j) in zip(obs, '0123')}
            ords.update({"N":".", "-":"."})
            for taxon in range(cnts.shape[1]):
                base1, base2 = DUCT[seq[taxon][site]]
                cnts[site][taxon] = "{}/{}:{}".format(\
                ords[base1],
                ords[base2],
                ",".join(carr[taxon][site].tolist())
                )
        gstr[seqleft:seqleft+seq.shape[1], 9:] = cnts

        ## advance counter
        seqleft += seq.shape[1]

    ## trim off empty rows if they exist
    empties = np.all(gstr == "", axis=1)
    gstr = gstr[~empties]

    ## append to vcf output file
    with open(data.outfiles.vcf+".{}".format(start), 'w') as vout:
        np.savetxt(vout, gstr, delimiter="\t", fmt="%s")

    io5.close()
    co5.close()



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
