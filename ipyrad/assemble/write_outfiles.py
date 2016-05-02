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

from __future__ import print_function

import pandas as pd
import numpy as np
import datetime
import shutil
import time
import glob
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
OUTPUT_FORMATS = ['phy', 'nex', 'snps', 'usnps', 'str', 'geno']



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

    ## Apply filters to supercatg and superhdf5 with selected samples
    ## and fill the filters and edge arrays.
    LOGGER.info("Applying filters")
    filter_all_clusters(data, samples, ipyclient)

    ## Everything needed is in the now filled h5 database. Filters were applied
    ## with 'samples' taken into account. Now we create the loci file (default)
    ## output and get some data back for building the stats file
    ## 'keep' has the non-filtered loci idxs. 
    start = time.time()
    LOGGER.info("Writing .loci file")
    elapsed = datetime.timedelta(seconds=int(time.time()-start))
    progressbar(20, 0, 
        " building output files    | {}".format(elapsed))
    samplecounts, locuscounts, keep = make_loci(data, samples)

    ## Write stats file output
    LOGGER.info("Writing stats output")
    make_stats(data, samples, samplecounts, locuscounts)

    ## OPTIONAL OUTPUTS
    ## grab from params as a string, if commas then split to a list
    ## phy, nex, str
    output_formats = data.paramsdict["output_formats"]
    if "," in output_formats:
        output_formats = [i.strip() for i in output_formats.split(",")]
    if "*" in output_formats:
        output_formats = OUTPUT_FORMATS

    ## not included in *
    if 'vcf' in output_formats:
        LOGGER.info("Writing .vcf file")        
        if data._headers:
            print("  building vcf ")
        make_vcf(data, samples, keep)

    ## make other array-based formats
    LOGGER.info("Writing other formats")
    make_outfiles(data, samples, keep, output_formats, ipyclient)

    ## print friendly message
    shortpath = data.dirs.outfiles.replace(os.path.expanduser("~"), "~")
    print("  Outfiles written to: {}".format(shortpath))



def make_stats(data, samples, samplecounts, locuscounts):
    """ write the output stats file and save to Assembly obj."""
    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')
    anames = inh5["seqs"].attrs["samples"]
    nloci = inh5["seqs"].shape[0]
    optim = inh5["seqs"].attrs["chunksize"]#    optim = 1000

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
    inh5.close()



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
        start = time.time()
        results = [] #submitted_args = []
        submitted = 0
        while submitted < nloci:
            hslice = np.array([submitted, submitted+optim])
            async = lbview.apply(filter_stacks, [data, sidx, hslice])
            results.append(async)
            #submitted_args.append([data, sidx, hslice])
            submitted += optim

        ## run filter_stacks on all chunks
        while 1:
            readies = [i.ready() for i in results]
            if not all(readies):
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(len(readies), sum(readies), 
                    " filtering loci           | {}".format(elapsed))
                time.sleep(1)
            else:
                break
        progressbar(20, 20, " filtering loci           | {}".format(elapsed))        
        if data._headers:
            print("")

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



## incorportatig samples...
def make_loci(data, samples):
    """ 
    Makes the .loci file from h5 data base. Iterates by 1000 loci at a 
    time and write to file. 
    """

    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')

    ## open the out handle
    data.outfiles.loci = os.path.join(data.dirs.outfiles, data.name+".loci")
    locifile = open(data.outfiles.loci, 'w')

    ## iterate 1000 loci at a time
    #optim = 1000
    optim = inh5["seqs"].attrs["chunksize"]    
    nloci = inh5["seqs"].shape[0]

    ## get sidx of samples
    anames = inh5["seqs"].attrs["samples"]
    ## get name and snp padding
    pnames, snppad = padnames(anames)
    snames = [i.name for i in samples]
    smask = np.array([i not in snames for i in anames])

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
                    smask, samplecov, locuscov, start]
            if edg[4]:
                outstr, samplecov, locuscov = enter_pairs(*args)
                store.append(outstr)
            else:
                outstr, samplecov, locuscov = enter_singles(*args)
                store.append(outstr)

        ## write to file and clear store
        locifile.write("\n".join(store) + "\n")
        store = []
    
        ## increase the counter
        start += optim

    ## close handle
    locifile.close()
    inh5.close()

    ## return sample counter
    return samplecov, locuscov, keep



def enter_pairs(iloc, pnames, snppad, edg, aseqs, asnps, 
                smask, samplecov, locuscov, start):
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
              "    "+"".join(snpstring2)+"|{}|".format(iloc+start)

    return outstr, samplecov, locuscov



def enter_singles(iloc, pnames, snppad, edg, aseqs, asnps, 
                  smask, samplecov, locuscov, start):
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
    outstr += "\n" + snppad + "".join(snpstring) + "|{}|".format(iloc+start)

    return outstr, samplecov, locuscov



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
    data, sidx, hslice = args

    try:
        ioh5 = h5py.File(data.database, 'r')
        ## get a chunk (hslice) of loci for the selected samples (sidx)
        superseqs = ioh5["seqs"][hslice[0]:hslice[1], sidx,]
        LOGGER.info('logging from engine: %s', superseqs.shape)

        ## get first index of splits for paired reads. Only need this until 
        ## get_edges gets the edges and then these are stored together
        splits = ioh5["edges"][hslice[0]:hslice[1], 4]
        #LOGGER.info('splits preview %s', splits[:10])
        #LOGGER.info('splits preview %s', splits[-10:])        

        ## get edges of superseqs, since edges should be trimmed off before
        ## counting hets and snps. Technically, this could edge trim clusters
        ## to the point that they are below the minlen, and so this also 
        ## constitutes a filter, though one that is uncommon. For this reason
        ## we have a filter also called edgfilter.
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

    #LOGGER.info("edgemins %s", edgemins)
    #LOGGER.info("splits %s", splits)

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
    """ 
    Filter minimum # of samples per locus from superseqs[chunk]. The shape
    of superseqs is [chunk, sum(sidx), maxlen]
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
    are (chunk, sum(sidx), maxlen). Don't need split info since it applies to 
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



def make_outfiles(data, samples, keep, output_formats, ipyclient):
    """
    Get desired formats from paramsdict and write files to outfiles 
    directory.
    """

    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')

    ## will iterate optim loci at a time
    optim = inh5["seqs"].attrs["chunksize"]
    nloci = inh5["seqs"].shape[0]

    ## get name and snp padding
    anames = inh5["seqs"].attrs["samples"]
    snames = [i.name for i in samples]
    names = [i for i in anames if i in snames]
    pnames, _ = padnames(names)
    pnames.sort()

    ## get names boolean
    sidx = np.array([i in snames for i in anames])
    assert len(pnames) == sum(sidx)
    ## get names index in order of pnames
    #sindx = [list(anames).index(i) for i in snames]

    ## build arrays and outputs from arrays
    arrs = make_arrays(data, sidx, optim, nloci, keep, inh5)
    seqarr, snparr, bisarr, _ = arrs

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
                " building output files    | {}".format(elapsed))
            time.sleep(1)
        else:
            break
    ## final progress bar
    elapsed = datetime.timedelta(seconds=int(time.time()-start))            
    progressbar(20, 20, " building output files    | {}".format(elapsed))        
    if data._headers:
        print("")

    ## check for errors
    for async in results:
        if not async.completed:
            print(async.metadata.error)

    ## close h5 handle
    inh5.close()




def make_arrays(data, sidx, optim, nloci, keep, inh5):
    """ 
    Writes loci and VCF formats and builds and returns arrays for contructing
    all other output types if specified. 
    """

    ## make empty arrays for filling
    maxlen = data._hackersonly["max_fragment_length"]
    maxsnp = inh5["snps"][:].sum()

    ## shape of arrays is sidx, we will subsample h5 w/ sidx to match
    seqarr = np.zeros((sum(sidx), maxlen*nloci), dtype="S1")
    snparr = np.zeros((sum(sidx), maxsnp), dtype="S1")
    bisarr = np.zeros((sum(sidx), nloci), dtype="S1")        

    ## apply all filters and write loci data
    start = 0
    seqleft = 0
    snpleft = 0
    bis = 0

    #LOGGER.info("starting to build phy")
    while start < nloci:
        hslice = [start, start+optim]
        afilt = inh5["filters"][hslice[0]:hslice[1], ...]
        aedge = inh5["edges"][hslice[0]:hslice[1], ...]
        aseqs = inh5["seqs"][hslice[0]:hslice[1], sidx, ...]
        asnps = inh5["snps"][hslice[0]:hslice[1], ...]

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]

        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            edg = aedge[iloc]
            ## grab all seqs between edges
            seq = aseqs[iloc, :, edg[0]:edg[1]]

            ## grab SNPs from seqs which is sidx subsampled
            ## and built after filters are applied, so subsampling
            ## does not get SNPs wrong.
            getsnps = asnps[iloc].sum(axis=1).astype(np.bool)
            snps = aseqs[iloc, :, getsnps].T

            ## remove cols that are all N-
            lcopy = seq
            lcopy[lcopy == "-"] = "N"
            bcols = np.all(lcopy == "N", axis=0)
            seq = seq[:, ~bcols]

            ## put into large array
            seqarr[:, seqleft:seqleft+seq.shape[1]] = seq
            #print(seqarr[:, seqleft:seqleft+10])
            seqleft += seq.shape[1]

            ## subsample all SNPs into an array
            snparr[:, snpleft:snpleft+snps.shape[1]] = snps
            snpleft += snps.shape[1]

            ## subsample one SNP into an array
            if snps.shape[1]:
                samp = np.random.randint(snps.shape[1])
                bisarr[:, bis] = snps[:, samp]
                bis += 1

        ## increase the counter
        start += optim

    #LOGGER.info("done building phy... ")
    ## trim trailing edges b/c we made the array bigger than needed.
    ridx = np.all(seqarr == "", axis=0)    
    seqarr = seqarr[:, ~ridx]
    ridx = np.all(snparr == "", axis=0)
    snparr = snparr[:, ~ridx]
    ridx = np.all(bisarr == "", axis=0)
    bisarr = bisarr[:, ~ridx]
    #LOGGER.info("done trimming phy... ")

    ## return these three arrays which are pretty small
    ## catg array gets to be pretty huge, so we return only 
    ## the partition info TODO:
    partitions = 0
    return seqarr, snparr, bisarr, partitions



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
    ## for denovo data we just fake it and evenly space the unlinked SNPs
    # 1  rs123456  0  1234555
    # 1  rs234567  0  1237793
    # 1  rs224534  0  -1237697        <-- exclude this SNP
    # 1  rs233556  0  1337456        



def make_vcf(data, samples, inh5):
    """ 
    Write the SNPs-only VCF output 
    """

    ## load the h5 database
    inh5 = h5py.File(data.database, 'r')

    ## create outputs
    #data.outfiles.vcf = os.path.join(data.dirs.outfiles, 
    #                                    data.name+".snps.vcf")
    data.outfiles.vcf = os.path.join(data.dirs.outfiles, 
                                        data.name+".vcf")

    ## will iterate optim loci at a time
    maxlen = data._hackersonly["max_fragment_length"]
    optim = inh5["seqs"].attrs["chunksize"]
    nloci = inh5["seqs"].shape[0]

    ## get name and snp padding
    anames = inh5["seqs"].attrs["samples"]
    snames = [i.name for i in samples]
    names = [i for i in anames if i in snames]

    ## write headers
    vout = open(data.outfiles.vcf, 'w')
    vcfheader(data, names, vout)

    ## get names index
    sidx = np.array([i in snames for i in anames])

    ## apply all filters and write loci data
    start = 0
    seqleft = 0
    snpleft = 0

    ## empty array to be filled before writing 
    ## will not actually be optim*maxlen, extra needs to be trimmed
    gstr = np.zeros((optim*maxlen, 9+sum(sidx)), dtype="S20")

    #LOGGER.info("starting to build phy")
    while start < nloci:
        hslice = [start, start+optim]
        afilt = inh5["filters"][hslice[0]:hslice[1], ]
        aedge = inh5["edges"][hslice[0]:hslice[1], ]
        aseqs = inh5["seqs"][hslice[0]:hslice[1], sidx, ]
        acatg = inh5["catgs"][hslice[0]:hslice[1], sidx, ]

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]

        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            edg = aedge[iloc]
            ## grab all seqs between edges
            seq = aseqs[iloc, :, edg[0]:edg[1]]
            catg = acatg[iloc, :, edg[0]:edg[1]]

            ## ----  build string array ---- 
            ## fill (CHR) chromosome/contig (reference) or RAD-locus (denovo)
            gstr[seqleft:seqleft+seq.shape[1], 0] = "RAD_{}_".format(iloc)
            ## fill (POS) position
            gstr[seqleft:seqleft+seq.shape[1], 1] = range(seq.shape[1])
            ## fill (ID) what is this? missing value is .
            gstr[seqleft:seqleft+seq.shape[1], 2] = "."
            ## fill (REF) must be ACGTN, TODO: indels gets complicated
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

        #pd.DataFrame(gstr)to_string(buf=vout)
        np.savetxt(vout, gstr, delimiter="\t", fmt="%s")
        ## write locus chunk to file and clear
        ## ...

        ## 
        start += optim

    vout.close()        
    inh5.close()



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
