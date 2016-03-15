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

from __future__ import print_function

import pandas as pd
import numpy as np
import itertools
import shutil
import h5py
import os
from collections import Counter, OrderedDict
from ipyrad.file_conversion import *
from util import *

import logging
LOGGER = logging.getLogger(__name__)


## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized. LOCI and VCF are default. Some others are created as 
## dependencies of others.
OUTPUT_FORMATS = ['alleles', 'phy', 'nex', 'snps', 'usnps',
                  'str', 'geno', 'treemix', 'migrate', 'gphocs']


def run(data, samples, force, ipyclient):
    """ Check all samples requested have been clustered (state=6), make output 
    directory, then create the requested outfiles. Excluded samples are already
    removed from samples.
    """
    ## prepare dirs
    data.dirs.outfiles = os.path.join(data.dirs.project, data.name+"_outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    ## Apply filters to supercatg and superhdf5 with selected samples
    ## and fill the filters and edge arrays.
    LOGGER.info("Applying filters")
    filter_all_clusters(data, samples, ipyclient)

    ## Everything needed is in the now filled h5 database. Filters were applied
    ## with 'samples' taken into account, but still need to print only included
    LOGGER.info("Writing .loci file")
    samplecounts, locuscounts, keep = make_loci(data, samples)

    LOGGER.info("Writing stats output")
    make_stats(data, samples, samplecounts, locuscounts)

    ## 'keep' has the non-filtered loci idx
    LOGGER.info("Writing .vcf file")
    #make_vcf(data)

    LOGGER.info("Writing other formats")
    make_outfiles(data, samples, keep)

    shortpath = data.dirs.outfiles.replace(os.path.expanduser("~"), "~")
    print("    Outfiles written to: {}".format(shortpath))



def make_stats(data, samples, samplecounts, locuscounts):
    """ write the output stats file and save to Assembly obj."""
    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')
    anames = inh5["seqs"].attrs["samples"]
    nloci = inh5["seqs"].shape[0]
    optim = 1000

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
    for i in range(50):
        piscounts[i] = 0
        varcounts[i] = 0
    #bisnps = Counter()
    while start < nloci:
        hslice = [start, start+optim]
        ## load each array
        afilt = inh5["filters"][hslice[0]:hslice[1], ]
        asnps = inh5["snps"][hslice[0]:hslice[1], ]

        ## get subarray results from filter array
        # "max_indels", "max_snps", "max_hets", "min_samps", "bad_edges"
        #print('afilt - ', afilt)
        filters += afilt.sum(axis=0)
        #print('filters - ', filters)
        passed += np.sum(afilt.sum(axis=1) == 0)

        ## get snps counts
        snplocs = asnps[:].sum(axis=1)
        varlocs = snplocs.sum(axis=1)
        varcounts.update(Counter(varlocs))
        #snpcounts.update(Counter(snplocs[:, 0]))
        piscounts.update(Counter(snplocs[:, 1]))
        ## bisnps TODO snpcols[2]

        ## increase counter to advance through h5 database
        start += optim

    ## record filtering of loci from total to final
    filtdat = pd.Series(np.concatenate([[nloci], filters, [passed]]),
        name="locus_filtering", 
        index=[
        "total_prefiltered_loci", 
        "filtered_by_rm_duplicates",
        "filtered_by_max_indels", 
        "filtered_by_max_snps",
        "filtered_by_max_hetero",
        "filtered_by_min_sample",
        "filtered_by_edge_trim",
        "total_filtered_loci"])

    print("\n\n## The number of loci caught by each filter."+\
          "\n## ipyrad API location: [assembly].statsfiles.s7_filters\n", 
          file=outstats)
    data.stats_dfs.s7_filters = pd.DataFrame(filtdat)
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
    lrange = range(1, len(anames)+1)
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
    srange = range(max([i for i in varcounts if varcounts[i]])+1)
    vardat = pd.Series(varcounts, name="var", index=srange).fillna(0)
    sumd = {0: 0}
    varsums = pd.Series(sumd.update(
                        {i: np.sum(vardat.values[1:i]) for i in srange[1:]}), 
                        name="sum_var", index=srange)

    pisdat = pd.Series(piscounts, name="pis", index=srange).fillna(0)
    sumd = {0: 0}    
    pissums = pd.Series(sumd.update(
                        {i: np.sum(pisdat.values[1:i]) for i in srange[1:]}), 
                        name="sum_pis", index=srange)
    print("\n\n\n## The distribution of SNPs (var and pis) across loci."+\
          "\n## pis = parsimony informative site (minor allele in >1 sample)"+\
          "\n## var = all variable sites (pis + autapomorphies)"+\
          "\n## ipyrad API location: [assembly].stats_dfs.s7_snps\n",
          file=outstats)
    data.stats_dfs.s7_snps = pd.concat([vardat, varsums, pisdat, pissums], 
                                        axis=1)    
    data.stats_dfs.s7_snps.to_string(buf=outstats)

    ## close it
    outstats.close()



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
    LOGGER.debug("in filter_all_clusters")
    ## create loadbalanced ipyclient
    lbview = ipyclient.load_balanced_view()

    ## get chunk size from the HD5 array and close
    with h5py.File(data.database, 'r') as inh5:
        ## the size of chunks for reading/writing
        optim = inh5["seqs"].attrs["chunksize"]
        ## the samples in the database in their locus order
        dbsamples = inh5["seqs"].attrs["samples"]
        ## the total number of loci        
        nloci = inh5["seqs"].shape[0]

    ## make a tmp directory for saving chunked arrays to
    chunkdir = os.path.join(data.dirs.outfiles, data.name+"_tmpchunks")
    if not os.path.exists(chunkdir):
        os.mkdir(chunkdir)

    ## get the indices of the samples that we are going to include
    sidx = select_samples(dbsamples, samples)

    ## Split dataset into chunks
    try:
        ## load a list of args to send to Engines. Each arg contains the index
        ## to sample 1000 loci from catg, seqs, filters &or edges, which will
        ## be loaded on the remote Engine. 

        ## create job queue
        submitted_args = []
        submitted = 0
        while submitted < nloci:
            hslice = np.array([submitted, submitted+optim])
            submitted_args.append([data, sidx, hslice])
            submitted += optim

        ## run filter_stacks on all chunks
        results = lbview.map_async(filter_stacks, submitted_args)
        results.get()

        ## get all the saved tmp arrays for each slice
        tmpsnp = glob.glob(os.path.join(chunkdir, "snpf.*.npy"))
        tmphet = glob.glob(os.path.join(chunkdir, "hetf.*.npy"))
        tmpmin = glob.glob(os.path.join(chunkdir, "minf.*.npy"))
        tmpedg = glob.glob(os.path.join(chunkdir, "edgf.*.npy"))

        ## sort array files within each group
        arrdict = OrderedDict([
            ('snp', tmpsnp), ('het', tmphet),
            ('min', tmpmin), ('edg', tmpedg)])
        for arrglob in arrdict.values():
            arrglob.sort(key=lambda x: int(x.rsplit(".")[-2]))

        ## load the full filter array who's order is
        ## ["duplicates", "max_indels", 
        ## "max_snps", "max_hets", "min_samps", "bad_edges"]
        inh5 = h5py.File(data.database, 'r+')
        superfilter = inh5["filters"]

        ## iterate across filter types
        for fidx, ftype in zip(range(2, 6), arrdict.keys()):
            ## fill in the edgefilters
            for ffile in arrdict[ftype]:
                ## grab a file and get it's slice            
                hslice = int(ffile.split(".")[-2])
                ## load in the array
                arr = np.load(ffile)
                ## store slice into full array
                superfilter[hslice:hslice+optim, fidx] = arr
        del arr
        ## finished filling filter array
        #LOGGER.info('superfilter[:10] : %s', superfilter[:10])        

        ## store the other arrayed values (edges, snps)
        edgarrs = glob.glob(os.path.join(chunkdir, "edgearr.*.npy"))
        snparrs = glob.glob(os.path.join(chunkdir, "snpsarr.*.npy"))
        ## sort array files within each group
        arrdict = OrderedDict([('edges', edgarrs), ('snps', snparrs)])
        for arrglob in arrdict.values():
            arrglob.sort(key=lambda x: int(x.rsplit(".")[-2]))


        ## fill the edge array 
        superedge = inh5['edges']
        for ffile in arrdict['edges']:
            ## grab a file and get it's slice            
            hslice = int(ffile.split(".")[-2])
            ## load in the array w/ shape (hslice, 5)
            arr = np.load(ffile)
            ## store slice into full array
            superedge[hslice:hslice+optim, 0:4, ] = arr
        del arr
        ## finished with superedge
        #LOGGER.info('superedge[:10] : %s', superedge[:10])

        ## fill the snps array. shape= (nloci, maxlen, 2)
        supersnps = inh5['snps']
        for ffile in arrdict['snps']:
            ## grab a file and get it's slice            
            hslice = int(ffile.split(".")[-2])
            ## load in the array w/ shape (hslice, maxlen, 2)
            arr = np.load(ffile)
            ## store slice into full array
            supersnps[hslice:hslice+optim, :, :] = arr
        del arr
        ## finished with superedge
        #LOGGER.info('supersnps[:10] : %s', supersnps[:10])


    finally:
        ## clean up the tmp files/dirs
        try:
            LOGGER.info("finished filtering")
            shutil.rmtree(chunkdir)
        except (IOError, OSError):
            pass



def padnames(anames, longname_len):
    """ pads names for loci output """
    ## Padding distance between name and seq.
    padding = 5    
    ## 
    pnames = [name + " " * (longname_len - len(name)+ padding) \
              for name in anames]
    snppad = "//" + " " * (longname_len - 2 + padding)
    return np.array(pnames), snppad



## incorportatig samples...
def make_loci(data, samples):
    """ makes the .loci file from h5 data base. Iterates by 1000 loci at a 
    time and write to file. """

    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')

    ## open the out handle
    data.outfiles.loci = os.path.join(data.dirs.outfiles, data.name+".loci")
    locifile = open(data.outfiles.loci, 'w')

    ## iterate 1000 loci at a time
    optim = 1000
    nloci = inh5["seqs"].shape[0]

    ## get sidx of samples
    anames = inh5["seqs"].attrs["samples"]
    ## get name and snp padding
    longname_len = max(len(i) for i in anames)
    pnames, snppad = padnames(anames, longname_len)
    samples = [i.name for i in samples]
    smask = np.array([i not in samples for i in anames])

    ## keep track of how many loci from each sample pass all filters
    samplecov = np.zeros(len(anames), dtype=int)
    locuscov = Counter()
    ## set initial value to zero for all values above min_samples_locus
    for cov in range(data.paramsdict["min_samples_locus"], len(anames)+1):
        locuscov[cov] = 0

    ## apply all filters and write loci data
    start = 0
    while start < nloci:
        hslice = [start, start+optim]
        afilt = inh5["filters"][hslice[0]:hslice[1], ]
        aedge = inh5["edges"][hslice[0]:hslice[1], ]
        aseqs = inh5["seqs"][hslice[0]:hslice[1], ]
        asnps = inh5["snps"][hslice[0]:hslice[1], ]

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]

        ## store until printing
        store = []

        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            edg = aedge[iloc]
            args = [iloc, pnames, snppad, edg, aseqs, asnps, 
                    smask, samplecov, locuscov]
            if edg[4]:
                outstr, samplecov, locuscov = enter_pairs(*args)
                store.append(outstr)
            else:
                outstr, samplecov, locuscov = enter_singles(*args)
                store.append(outstr)

            # ## grab all seqs between edges
            # seq = aseqs[iloc, :, edg[0]:edg[1]]
            # ## snps was created using only the selected samples.
            # snp = asnps[iloc, edg[0]:edg[1], ]
            # ## remove rows with all Ns, seq has only selected samples
            # nalln = np.all(seq == "N", axis=1)
            # ## make mask of removed rows and excluded samples. Use the inverse
            # ## of this to save the coverage for samples
            # nsidx = nalln + smask
            # samplecov += np.invert(nsidx).astype(int)
            # locuscov[np.sum(np.invert(nsidx).astype(int))] += 1
            # ## select the remaining names in order
            # seq = seq[~nsidx, ]
            # names = pnames[~nsidx]

            # ## save string for printing, excluding names not in samples
            # store.append("\n".join(\
            #     [name + s.tostring() for name, s in zip(names, seq)]))

            # ## get snp string and add to store
            # snpstring = ["-" if snp[i, 0] else \
            #              "*" if snp[i, 1] else \
            #              " " for i in range(len(snp))]
            # store.append(snppad + "".join(snpstring))

        ## write to file and clear store
        locifile.write("\n".join(store) + "\n")
        store = []
    
        ## increase the counter
        start += optim

    ## close handle
    locifile.close()

    ## return sample counter
    return samplecov, locuscov, keep


def enter_pairs(iloc, pnames, snppad, edg, aseqs, asnps, 
                smask, samplecov, locuscov):
    """ enters funcs for pairs """
    seq1 = aseqs[iloc, :, edg[0]:edg[1]]
    seq2 = aseqs[iloc, :, edg[2]+edg[4]+4:edg[3]+edg[4]]
    ## snps was created using only the selected samples.
    snp1 = asnps[iloc, edg[0]:edg[1], ]
    snp2 = asnps[iloc, edg[2]+edg[4]+4:edg[3]+edg[4], ]    
    ## remove rows with all Ns, seq has only selected samples
    nalln = np.all(seq1 == "N", axis=1)
    ## make mask of removed rows and excluded samples. Use the inverse
    ## of this to save the coverage for samples
    nsidx = nalln + smask
    samplecov += np.invert(nsidx).astype(int)
    locuscov[np.sum(np.invert(nsidx).astype(int))] += 1
    ## select the remaining names in order
    seq1 = seq1[~nsidx, ]
    seq2 = seq2[~nsidx, ]    
    names = pnames[~nsidx]
    ## save string for printing, excluding names not in samples
    outstr = "\n".join(\
        [name + s1.tostring()+"nnnn"+s2.tostring() for name, s1, s2 in \
         zip(names, seq1, seq2)])

    ## get snp string and add to store
    snpstring1 = ["-" if snp1[i, 0] else \
                 "*" if snp1[i, 1] else \
                 " " for i in range(len(snp1))]
    snpstring2 = ["-" if snp2[i, 0] else \
                 "*" if snp2[i, 1] else \
                 " " for i in range(len(snp2))]
    outstr += "\n" + snppad + "".join(snpstring1)+\
              "    "+"".join(snpstring2)+"|"

    return outstr, samplecov, locuscov



def enter_singles(iloc, pnames, snppad, edg, aseqs, asnps, 
                  smask, samplecov, locuscov):
    """ enter funcs for SE or merged data """
    ## grab all seqs between edges
    seq = aseqs[iloc, :, edg[0]:edg[1]]
    ## snps was created using only the selected samples.
    snp = asnps[iloc, edg[0]:edg[1], ]
    ## remove rows with all Ns, seq has only selected samples
    nalln = np.all(seq == "N", axis=1)
    ## make mask of removed rows and excluded samples. Use the inverse
    ## of this to save the coverage for samples
    nsidx = nalln + smask
    samplecov += np.invert(nsidx).astype(int)
    locuscov[np.sum(np.invert(nsidx).astype(int))] += 1
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
    outstr += "\n" + snppad + "".join(snpstring) + "|"

    return outstr, samplecov, locuscov



def filter_stacks(args):
    """ Grab a chunk of loci from the HDF5 database. Apply filters and fill the 
    the filters boolean array. 

    The design of the filtering steps intentionally sacrifices some performance
    for an increase in readability, and extensibility. Calling multiple filter
    functions ends up running through the sequences per stack several times, 
    but I felt this design made more sense, and also will easily allow us to
    add more filters in the future.
    """
    data, sidx, hslice = args

    try:
        ioh5 = h5py.File(data.database, 'r')
        ## get a chunk (hslice) of loci for the selected samples (sidx)
        superseqs = ioh5["seqs"][hslice[0]:hslice[1], sidx,]
        LOGGER.info('logging from engine: %s', superseqs.shape)

        ## get first index of splits for paired reads. Only need this until 
        ## get_edges gets the edges
        splits = ioh5["edges"][hslice[0]:hslice[1], 4]
        #LOGGER.info('splits preview %s', splits[:10])
        #LOGGER.info('splits preview %s', splits[-10:])        
        ## Empty arrays to be filled in used to later fill in the full 
        ## filter array. Filter dims = (nloci, 5)
        ## Filter order is: [dups, indels, maxhets, maxsnps, minsamp]
        #dupfilter = np.zeros((hslice[1],), dtype=np.bool)
        #indfilter = np.zeros((hslice[1],), dtype=np.bool)
        #hetfilter = np.zeros((hslice[1],), dtype=np.bool)
        #minfilter = np.zeros((hslice[1],), dtype=np.bool)
        #edgfilter = np.zeros((hslice[1],), dtype=np.bool)
        #snpfilter = np.zeros((hslice[1],), dtype=np.bool)

        ## get edges of superseqs, since edges should be trimmed off before
        ## counting hets and snps. Technically, this could edge trim clusters
        ## to the point that they are below the minlen, and so this also 
        ## constitutes a filter, though one that is uncommon. For this reason
        ## we have a filter also called edgfilter. (needs to be added to filter)
        LOGGER.debug("getting edges")
        edgfilter, edgearr = get_edges(data, superseqs, splits)
        del splits
        LOGGER.debug("edgfilter %s", edgfilter[:5])
        ## duplicates are already filtered during step6.

        ## minsamp coverages filtered from superseqs
        minfilter = filter_minsamp(data, superseqs)
        #LOGGER.debug("minfilter %s", minfilter[:5])

        ## maxhets per site column from superseqs after trimming edges
        hetfilter = filter_maxhet(data, superseqs, edgearr)

        ## Build the .loci snpstring as an array (snps) 
        ## shape = (chunk, 1) dtype=S1, or should it be (chunk, 2) for [-,*] ?
        snpfilter, snpsarr = filter_maxsnp(data, superseqs, edgearr)

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

        handle = os.path.join(chunkdir, "snpsarr.{}.npy".format(hslice[0]))
        with open(handle, 'w') as out:
            np.save(out, snpsarr)

        handle = os.path.join(chunkdir, "edgearr.{}.npy".format(hslice[0]))
        with open(handle, 'w') as out:
            np.save(out, edgearr)

    except Exception as inst:
        ## Do something real here
        LOGGER.warn(inst)
        raise

    finally:
        pass



def get_edges(data, superseqs, splits):
    """ 
    Gets edge trimming based on the overlap of sequences at the edges of 
    alignments and the tuple arg passed in for edge_trimming. Trims as
    (R1 left, R1 right, R2 left, R2 right). We also trim off the restriction
    site if it present.
    """
    ## the filtering arg and parse it into minsamp numbers
    edgetuple = data.paramsdict["trim_overhang"]
    minsamp = data.paramsdict["min_samples_locus"]
    allsamp = superseqs.shape[1]
    cut1, cut2 = data.paramsdict["restriction_overhang"]        

    ## a local array for storing edge trims
    edges = np.zeros((superseqs.shape[0], 4), dtype=np.int16)
    ## a local array for storing edge filtered loci
    edgefilter = np.zeros((superseqs.shape[0],), dtype=np.bool)

    ## X means trim this number bases from the edge of all reads.
    ## 4s mean no missing data (superseqs.shape[0])
    ##    -- raise warning if minsamp < 4    
    ## 3s mean trim cut sites and overhangs til there is data for >= minsamp.
    ## 2s mean trim cut sites and overhangs til there is data for >= 4 samples.
    ##    -- raise warning if minsamp < 4
    ## 1s mean trim cut sites only.
    ## 0s mean no edge trimming at all, not even cut sites.
    ## TODO: more testing on this with pairddrad.

    ## get edges for overhang trimming first, and trim cut sites at end
    ## use .get so that if not in 2,3,4 it returns None
    edgedict = {0: 1,
                1: 1,
                2: 4,
                3: minsamp,
                4: allsamp}
    edgemins = [edgedict.get(i) for i in edgetuple]

    LOGGER.info("edgemins %s", edgemins)
    LOGGER.info("splits %s", splits)

    ## convert all - to N to make this easier
    superseqs[superseqs == "-"] = "N"
    ## the background fill of superseqs should be N instead of "", then this 
    ## wouldn't be necessary...
    ## superseqs[superseqs == ""] = "N"

    ## trim overhanging edges
    ## get the number not Ns in each site, 
    ccx = np.sum(superseqs != "N", axis=1)
    #LOGGER.debug("ccx %s", ccx)
    #LOGGER.debug("splits %s", splits)
    for idx, split in enumerate(splits):
        if split:
            r1s = ccx[idx, :split]
            r2s = ccx[idx, split+4:]
        else:
            r1s = ccx[idx, :]
        ## set default values
        edge0 = edge1 = edge2 = edge3 = 0
        
        ## if edge trim fails then locus is filtered
        try:
            edge0 = max(0, np.where(r1s >= edgemins[0])[0].min())
            edge1 = np.where(r1s >= edgemins[1])[0].max()
        except ValueError:
            #edgefilter[idx] = True
            LOGGER.debug("edge filtered")
            edge1 = np.where(r1s >= 1)[0].max()
            #LOGGER.debug("Xccx %s", ccx[idx])            
            #LOGGER.debug("Xsplit %s", split)
            #LOGGER.debug("Xr1s %s", r1s)

        ## filter cut1
        if edgetuple[0]:
            edge0 = edge0+len(cut1)
        else:
            assert edge0 < edge1             

        LOGGER.debug("edges r1 (%s, %s)", edge0, edge1)

        ## if split then do the second reads separate
        if split:
            try:
                edge2 = np.where(r2s >= edgemins[2])[0].min()
                edge3 = np.where(r2s >= edgemins[3])[0].max()
            except ValueError:
                #edgefilter[idx] = True
                edge2 = edge1 + 4
                edge3 = np.where(r2s >= 1)[0].max()

            ## filter cut2
            if edgetuple[3]:
                edge3 = edge3-len(cut2)
            else:
                assert edge2 < edge3

        #LOGGER.debug("edges r2 (%s, %s)", edge2, edge3)

        ## store edges
        edges[idx] = np.array([edge0, edge1, edge2, edge3])

    return edgefilter, edges



def filter_minsamp(data, superseqs):
    """ Filter minimum # of samples per locus from superseqs[chunk]. The shape
    of superseqs is [chunk, len(sidx), maxlen]
    """
    ## the minimum filter
    minsamp = data.paramsdict["min_samples_locus"]
    ## ask which rows are not all N along seq dimension, then sum along sample 
    ## dimension to get the number of samples that are not all Ns.
    ## minfilt is a boolean array where True means it failed the filter.
    minfilt = np.sum(~np.all(superseqs == "N", axis=2), axis=1) < minsamp
    LOGGER.debug("minfilt %s", minfilt)

    ## print the info
    LOGGER.info("Filtered by min_samples_locus - {}".format(minfilt.sum()))
    return minfilt



def ucount(sitecol):
    """ 
    Used to count the number of unique bases in a site for snpstring. 
    returns as a spstring with * and - 
    """
    ## make into set
    site = set(sitecol)

    ## get resolutions of ambiguitites
    for iupac in "RSKYWM":
        if iupac in site:
            site.discard(iupac)
            site.update(AMBIGS[iupac])

    ## remove - and N
    site.discard("N")
    site.discard("-")

    ## if site is invariant return ""
    if len(site) < 2:
        return " "
    else:
        ## get how many bases come up more than once 
        ccx = np.sum(np.array([np.sum(sitecol == base) for base in site]) > 1)
        ## if another site came up more than once return *
        if ccx > 1:
            return "*"
        else:
            return "-"



def filter_maxsnp(data, superseqs, edges):
    """ 
    Filter max # of SNPs per locus. Do R1 and R2 separately if PE. 
    Also generate the snpsite line for the .loci format and save in the snp arr
    """

    maxs1, maxs2 = data.paramsdict["max_SNPs_locus"]

    ## an empty array to fill with failed loci
    snpfilt = np.zeros(superseqs.shape[0], dtype=np.bool)
    snpsarr = np.zeros((superseqs.shape[0], superseqs.shape[2], 2), 
                       dtype=np.bool)

    ## get the per site snp string (snps) shape=(chunk, maxlen)
    snps = np.apply_along_axis(ucount, 1, superseqs)
    ## fill snps array
    snpsarr[:, :, 0] = snps == "-"
    snpsarr[:, :, 1] = snps == "*"

    ## count how many snps are within the edges and fill snpfilt
    for idx, edg in enumerate(edges):
        nsnps = snpsarr[idx, edg[0]:edg[1]].sum(axis=1).sum()
        if nsnps > maxs1:
            snpfilt[idx] = True

    ## do something much slower for paired data, iterating over each locus and 
    ## splitting r1s from r2s and treating each separately.
    if "pair" in data.paramsdict["datatype"]:
        for idx, edg in enumerate(edges):
            nsnps = snpsarr[idx, edg[2]:edg[3]].sum(axis=1).sum()
            if nsnps > maxs2:
                snpfilt[idx] = True

    ## return snpsarr and snpfilter
    return snpfilt, snpsarr



def filter_maxhet(data, superseqs, edges):
    """ 
    Filter max shared heterozygosity per locus. The dimensions of superseqs
    are (chunk, len(sidx), maxlen). Don't need split info since it applies to 
    entire loci based on site patterns (i.e., location along the seq doesn't 
    matter.) Current implementation does ints, but does not apply float diff
    to every loc based on coverage... 
    """
    ## the filter max
    maxhet = data.paramsdict["max_shared_Hs_locus"]
    if isinstance(maxhet, float):
        maxhet = superseqs.shape[1]*float(maxhet)

    ## an empty array to fill with failed loci
    hetfilt = np.zeros(superseqs.shape[0], dtype=np.bool)

    ## get the per site number of bases in ambig, and get the max per site 
    ## value for each locus, then get boolean array of those > maxhet.
    for idx, edg in enumerate(edges):
        LOGGER.debug("testing %s, %s", idx, edg)
        share = np.array(\
            [np.sum(superseqs[idx, :, edg[0]:edg[1]] == ambig, axis=0)\
            .max() for ambig in "RSKYWM"]).max()
        ## fill filter
        if not share <= maxhet:
            hetfilt[idx] = True

    LOGGER.info("Filtered max_shared_heterozygosity- {}".format(hetfilt.sum()))
    return hetfilt



def make_outfiles(data, samples, keep):
    """ Get desired formats from paramsdict and write files to outfiles 
    directory 
    """
    ## Read in the input .loci file that gets transformed into other formats
    ## locifile = os.path.join(data.dirs.outfiles, data.name+".loci")
    #samples = [i.name for i in samples]

    ## select the output formats. Should we allow some kind of fuzzy matching?
    output_formats = data.paramsdict["output_formats"]
    if "*" in output_formats:
        output_formats = OUTPUT_FORMATS

    ## build arrays and outputs from arrays
    make_phynex(data, samples, keep, output_formats)

    # ## run each func
    # for filetype in output_formats:
    #     LOGGER.info("Doing - {}".format(filetype))

    #     # phy & nex come from loci2phynex
    #     if filetype in ["phy", "nex"]:
    #         filetype = "phynex"

    #     # All these file types come from loci2SNP
    #     elif filetype in ["snps", "usnps", "str", "geno"]:
    #         filetype = "SNP"

        ## Everything else has its own loci2*.py conversion file.
        ## Get the correct module for this filetype
        ## globals() here gets the module name, and then we can call
        ## the .make() function. This is a little tricky.
        #format_module = globals()["loci2"+filetype]

        ## do the call to make the new file format
        #format_module.make(data, samples)



def make_phynex(data, samples, keep, output_formats):
    """ make phylip and nexus formats. Also pulls out SNPs...? """

    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')

    ## iterate 1000 loci at a time
    optim = inh5["seqs"].attrs["chunksize"]
    nloci = inh5["seqs"].shape[0]

    ## get name and snp padding
    anames = inh5["seqs"].attrs["samples"]
    longname_len = max(len(i) for i in anames)
    pnames, _ = padnames(anames, longname_len)
    samples = [i.name for i in samples]
    smask = np.array([i not in samples for i in anames])

    ## make empty array for filling
    maxlen = data._hackersonly["max_fragment_length"]
    maxsnp = inh5["snps"][:].sum()

    seqarr = np.zeros((smask.shape[0], maxlen*nloci), dtype="S1")
    snparr = np.zeros((smask.shape[0], maxsnp), dtype="S1")
    bisarr = np.zeros((smask.shape[0], nloci), dtype="S1")        

    ## apply all filters and write loci data
    start = 0
    seqleft = 0
    snpleft = 0
    LOGGER.info("starting to build phy")
    while start < nloci:
        hslice = [start, start+optim]
        afilt = inh5["filters"][hslice[0]:hslice[1], ]
        aedge = inh5["edges"][hslice[0]:hslice[1], ]
        aseqs = inh5["seqs"][hslice[0]:hslice[1], ]
        asnps = inh5["snps"][hslice[0]:hslice[1], ]

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]

        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            edg = aedge[iloc]
            ## grab all seqs between edges
            seq = aseqs[iloc, :, edg[0]:edg[1]]
            ## grab SNPs
            getsnps = asnps[iloc].sum(axis=1).astype(np.bool)
            snps = aseqs[iloc, :, getsnps].T

            ## remove cols that are all N-
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

            ## subsample one SNP into an array
            if snps.shape[1]:
                samp = np.random.randint(snps.shape[1])
                bisarr[:, iloc] = snps[:, samp]
    
        ## increase the counter
        start += optim

    LOGGER.info("done building phy... ")
    ## trim trailing edges
    ridx = np.all(seqarr == "", axis=0)    
    seqarr = seqarr[:, ~ridx]
    ridx = np.all(snparr == "", axis=0)
    snparr = snparr[:, ~ridx]
    ridx = np.all(bisarr == "", axis=0)
    bisarr = bisarr[:, ~ridx]
    LOGGER.info("done trimming phy... ")

    ## write the phylip string
    data.outfiles.phy = os.path.join(data.dirs.outfiles, data.name+".phy")
    with open(data.outfiles.phy, 'w') as out:
        ## trim down to size
        #ridx = np.all(seqarr == "", axis=0)
        out.write("{} {}\n".format(seqarr.shape[0], seqarr.shape[1]))
                                   #seqarr[:, ~ridx].shape[1]))
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(seqarr[idx])))

    ## write the snp string
    data.outfiles.snp = os.path.join(data.dirs.outfiles, data.name+".snp")    
    with open(data.outfiles.snp, 'w') as out:
        #ridx = np.all(snparr == "", axis=0)
        out.write("{} {}\n".format(snparr.shape[0], snparr.shape[1]))
                                   #snparr[:, ~ridx].shape[1]))
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(snparr[idx])))

    ## write the bisnp string
    data.outfiles.usnp = os.path.join(data.dirs.outfiles, data.name+".usnp")
    with open(data.outfiles.usnp, 'w') as out:
        out.write("{} {}\n".format(bisarr.shape[0], bisarr.shape[1]))
                                   #bisarr.shape[1]))
        for idx, name in enumerate(pnames):
            out.write("{}{}\n".format(name, "".join(bisarr[idx])))

    ## Write STRUCTURE format
    if "str" in output_formats:
        data.outfiles.str = os.path.join(data.dirs.outfiles, data.name+".str")
        data.outfiles.ustr = os.path.join(data.dirs.outfiles, data.name+".ustr")        
        out1 = open(data.outfiles.str, 'w')
        out2 = open(data.outfiles.ustr, 'w')
        numdict = {'A': '0', 'T': '1', 'G': '2', 'C': '3', 'N': '-9', '-': '-9'}
        if data.paramsdict["max_alleles_consens"] > 1:
            for idx, name in enumerate(pnames):
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in snparr[idx]])))
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][1]] for i in snparr[idx]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in bisarr[idx]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][1]] for i in bisarr[idx]])))
        else:
            ## haploid output
            for idx, name in enumerate(pnames):
                out1.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in snparr[idx]])))
                out2.write("{}\t\t\t\t\t{}\n"\
                    .format(name,
                    "\t".join([numdict[DUCT[i][0]] for i in bisarr[idx]])))
        out1.close()
        out2.close()

    ## Write GENO format
    if "geno" in output_formats:
        data.outfiles.geno = os.path.join(data.dirs.outfiles, 
                                          data.name+".geno")
        data.outfiles.ugeno = os.path.join(data.dirs.outfiles, 
                                          data.name+".ugeno")        
        ## need to define a reference base and record 0,1,2 or missing=9
        #snparr
        #bisarr
    LOGGER.info("done writing outputs... ")


## Utility subfunctions
def count_shared_polymorphisms(seqs):
    """ 
    Count the number of shared polymorphic sites at every base in a locus. If 
    base contains too may shared polymorphisms (paramsdict[max_shared_Hs_locus])
    then we'll throw out the whole locus.
    Input is a list of sequences, output is a list of ints representing counts
    of shared polyorphisms per base
    """
    ## Transpose the list of seqs, so we now have a list of bases at each locus.
    ## Stacks of sequences should always be the same length, but if they aren't
    ## we'll protect by filling with N's, rather than truncating, which is what
    ## zip() would do.
    ## Makes a list of counters for each stack of bases of this locus
    stacks = [Counter(i) for i in itertools.izip_longest(*seqs, fillvalue="N")]

    ## This is maybe 20% too clever, but i wanted to do it with list 
    ## comprehension, so here we are. Read this as "At each base get the max 
    ## number of counts of any ambiguity code". Below y.items is a counter for 
    ## each position, and x[0] is the counter key, x[1] is the count for that 
    ## key. Another stupid thing is that max() in python 2.7 doesn't handle 
    ## empty lists, so there's this trick max(mylist or [0]), empty list 
    ## evaluates as false, so it effectively defaults max({}) = 0.
    max_counts = [max([x[1] for x in y.items() if x[0] in "RYSWKM"] or [0]) \
                  for y in stacks]

    return max_counts



def count_snps(seqs):
    """ Finds * and - snps to create snp string for .loci file """
    ## As long as we're ripping through looking for snps, keep track of the snp
    ## array to output to the .loci file
    snpsite = [" "]*len(seqs[0])

    ## Transpose to make stacks per base, same as count_shared_polymorphisms
    ## so see that function for more explicit documentation.
    stacks = [Counter(i) for i in itertools.izip_longest(*seqs, fillvalue="N")]

    ## At each stack, if the length of the counter is > 1 then its a snp
    for i, stack in enumerate(stacks):
        if len(stack) > 1:
            LOGGER.debug("Got a snp {}".format(stack))
            ## Count the number of keys that occur once = autapomorphies
            autapomorphies = stack.values().count(1)

            LOGGER.debug("Autapomorphies - {}".format(autapomorphies))
            ## If the number of autapomorphies is one less than the total len
            ## of the counter then this is still a parsimony uninformative site.                
            if autapomorphies == (len(stack) - 1):
                snpsite[i] = "-"
            ## If not then it is a parsimony informative site
            else:
                snpsite[i] = "*"

    ## When done, rip through the snpsite string and count snp markers "*" & "-"
    nsnps = sum(1 for x in snpsite if x in "*-")

    ## Just return the snpsite line
    snps = ("//", "".join(snpsite)+"|")
    #seqs.append(("//", "".join(snpsite)+"|"))

    return nsnps, snps



def make_vcfheader(data, samples, outvcf):
    LOGGER.debug("Entering make_vcfheader()")



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
