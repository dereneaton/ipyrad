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

from __future__ import print_function

import numpy as np
import itertools
import shutil
import h5py
import os
from collections import Counter
from ipyrad.file_conversion import *
from util import *

import logging
LOGGER = logging.getLogger(__name__)


## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized.
OUTPUT_FORMATS = ['alleles', 'phy', 'nex', 'snps', 'vcf', 'usnps',
                  'str', 'geno', 'treemix', 'migrate', 'gphocs']


def run(data, samples, force, ipyclient):
    """ Check all samples requested have been clustered (state=6), 
    make output directory, then create the requested outfiles.
    """

    LOGGER.info("Checking input")
    ## Make sure all samples are ready for writing output
    samples = precheck(data, samples)

    LOGGER.info("Applying filters")
    ## Apply filters to supercatg and superhdf5 
    ## and fill the filters and edge arrays
    filter_all_clusters(data, samples, ipyclient)

    LOGGER.info("Make .loci from filtered .vcf")
    ## Make .loci and .vcf from the filtered superseq & catg arrays
    #vcf2loci.make(data, samples)

    LOGGER.info("Convert .loci to all requested output file formats")
    ## Make all requested outfiles from the filtered .loci file
    #make_outfiles(data, samples, force)



def precheck(data, samples):
    """ Check all samples requested have been clustered (state=6), 
    make output directory, then create the requested outfiles.
    """
    ## Make a list of all the samples that are actually ready
    subsample = []
    for sample in samples:
        if sample.stats.state < 6:
            print("""
    Excluding Sample {}; not ready for writing output. Run step6() first.
    """.format(sample.name))
        
        else:
            subsample.append(sample)

    ## prepare dirs
    data.dirs.outfiles = os.path.join(data.dirs.project, data.name+"_outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    return subsample



def filter_all_clusters(data, samples, ipyclient):
    """ 
    Open the HDF5 array with seqs, catg, and filter data. Fill the remaining
    filters. 
    """
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
    sidx = select_samples(data, dbsamples, [i.name for i in samples])

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
        tmpedg = glob.glob(os.path.join(chunkdir, "edgf.*.npy"))
        tmphet = glob.glob(os.path.join(chunkdir, "hetf.*.npy"))
        tmpsnp = glob.glob(os.path.join(chunkdir, "snpf.*.npy"))
        tmpmin = glob.glob(os.path.join(chunkdir, "minf.*.npy"))
        ## sort array files within each group
        arrdict = {'edg':tmpedg, 'het':tmphet, 'snp':tmpsnp, 'min':tmpmin}
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
        LOGGER.info('superfilter[:10] : %s', superfilter[:10])        

        ## store the other arrayed values (edges, snps)
        edgarrs = glob.glob(os.path.join(chunkdir, "edgearr.*.npy"))
        snparrs = glob.glob(os.path.join(chunkdir, "snpsarr.*.npy"))
        ## sort array files within each group
        arrdict = {'edges':edgarrs, 'snps':snparrs}
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
        LOGGER.info('superedge[:10] : %s', superedge[:10])

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
        LOGGER.info('supersnps[:10] : %s', supersnps[:10])

        ### MAKE THE LOCI FILE
        ## h5 handle loads superseqs, superedge, supersnps, superfilter
        ## and applies all filters.
        make_loci(data, inh5)

        ## MAKE THE STATS OUTPUT

        ## MAKE OTHER OUTPUT FILES FROM CATG/VCF

        ## MAKE OTHER OUTPUT FILES FROM LOCI



    finally:
        ## clean up the tmp files/dirs
        try:
            print("finished filtering")
            shutil.rmtree(chunkdir)
        except (IOError, OSError):
            pass



def padnames(names, longname_len):
    """ pads names for loci output """
    ## Padding distance between name and seq.
    padding = 5    
    ## 
    pnames = [name + " " * (longname_len - len(name)+ padding) \
              for name in names]
    snppad = "//" + " " * (longname_len - 2 + padding)
    return np.array(pnames), snppad



def make_loci(data, inh5):
    """ makes the .loci file from h5 data base. Iterates by 1000 loci at a 
    time and write to file. """

    ## open the out handle
    locifile = open(os.path.join(data.dirs.outfiles, data.name+".loci"), 'w')

    ## iterate 1000 loci at a time
    optim = 1000
    nloci = inh5["seqs"].shape[0]

    start = 0
    while start < nloci:
        hslice = [start, start+optim]
        afilt = inh5["filters"][hslice[0]:hslice[1], ]
        aedge = inh5["edges"][hslice[0]:hslice[1], ]
        aseqs = inh5["seqs"][hslice[0]:hslice[1], ]
        asnps = inh5["snps"][hslice[0]:hslice[1], ]
        anames = inh5["seqs"].attrs["samples"]

        ## get name and snp padding
        longname_len = max(len(i) for i in anames)
        pnames, snppad = padnames(anames, longname_len)

        ## which loci passed all filters
        keep = np.where(np.sum(afilt, axis=1) == 0)[0]

        ## store until printing
        store = []
        ## write loci that passed after trimming edges, then write snp string
        for iloc in keep:
            edg = aedge[iloc]
            seq = aseqs[iloc, :, edg[0]:edg[1]]
            snp = asnps[iloc, edg[0]:edg[1], ]
            ## remove rows with all Ns
            nsidx = np.all(seq[:, edg[0]:edg[1]] == "N", axis=1)
            seq = seq[~nsidx, ]
            ## select the remaining names in order
            names = pnames[~nsidx]
            ## save string for printing
            store.append("\n".join(\
                [name + s.tostring() for name, s in zip(names, seq)]))
            ## get snp string and add to store
            snpstring = ["-" if snp[i, 0] else \
                         "*" if snp[i, 1] else \
                         " " for i in range(len(snp))]
            store.append(snppad + "".join(snpstring))

        ## write to file and clear store
        locifile.write("\n".join(store) + "\n")
        store = []
    
        ## increase the counter
        start += optim

    ## close handle
    locifile.close()



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
        print('sup', superseqs.shape)

        ## get first index of splits for paired reads. Only need this until 
        ## get_edges gets the edges
        splits = ioh5["edges"][hslice[0]:hslice[1], 4]
        
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
        edgfilter, edgearr = get_edges(data, superseqs, splits)
        del splits
        LOGGER.debug("edgfilter %s", edgfilter[:5])
        ## duplicates are already filtered during step6.

        ## indels are already filtered during step6 (kinda... todo: pair splits)

        ## minsamp coverages filtered from superseqs
        minfilter = filter_minsamp(data, superseqs)
        LOGGER.debug("minfilter %s", minfilter[:5])

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

    ## Write out .tmp loci
    #write_tmp_loci(data, loci, fname)


    ## Write out .tmp vcf
    #write_tmp_vcf(data, loci, fname)



def select_samples(data, dbsamples, samples):
    """ Get the row index of samples that are included, minus the samples that
    are excluded, from their order in the database (dbsamples)
    """
    ## Get all the samples to exclude. The 'or' conditional is a hack to account
    ## for the fact that paramsdict may either be an empty string or a list of
    ## names to exclude. List and "" don't append right, so if you do this 'or'
    ## and either of the lists is empty you'll get ["", "1B_0"]
    excludes = (data.paramsdict["excludes"] or [""]) \
             + (data.paramsdict["outgroups"] or [""])

    ## remove excludes from the list
    includes = list(set(samples) - set(excludes))
    print('inc', includes)
    excluded = set(samples).intersection(set(excludes))
    LOGGER.info("Excluded: {}".format(excluded))

    ## get index from dbsamples
    sidx = [list(dbsamples).index(i) for i in includes]
    sidx.sort()
    print("sidx: ", sidx)
    return sidx




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
    allsamp = superseqs.shape[0]
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

    ## get edges for overhang trimming first, and trim cut sites at end
    ## use .get so that if not in 2,3,4 it returns None
    edgedict = {0: None,
                1: 0,
                2: 4,
                3: minsamp,
                4: allsamp}
    edgemins = [edgedict.get(i) for i in edgetuple]
    LOGGER.debug("edgemins %s", edgemins)

    ## convert all - to N to make this easier
    superseqs[superseqs == "-"] = "N"
    ## the background fill of superseqs should be N instead of "", then this 
    ## wouldn't be necessary...
    superseqs[superseqs == ""] = "N"

    ## trim overhanging edges
    ## get the number not Ns in each site, 
    ccx = np.sum(superseqs != "N", axis=1)

    for idx, split in enumerate(splits):
        if split:
            r1s = ccx[idx, :split]
            r2s = ccx[idx, split+4:]
        else:
            r1s = ccx[idx, ]

        ## set default values
        edge0 = edge1 = edge2 = edge3 = 0
        
        ## if edge trim fails then locus is filtered
        try:
            edge0 = np.where(r1s > edgemins[0])[0].min()
            edge1 = np.where(r1s > edgemins[1])[0].max()
        except ValueError:
            edgefilter[idx] = True
        #LOGGER.debug("edge0 %s", edge0)
        #LOGGER.debug("edge1 %s", edge1)        

        ## filter cut1
        if edgetuple[0]:
            edge0 = edge0+len(cut1)

        ## if split then do the second reads separate
        if split:
            try:
                edge2 = np.where(r2s > edgemins[2])[0].min()
                edge3 = np.where(r2s > edgemins[3])[0].max()
            except ValueError:
                edgefilter[idx] = True
            #LOGGER.debug("edge0 %s", edge0)
            #LOGGER.debug("edge1 %s", edge1)        
            
            ## filter cut2
            if edgetuple[3]:
                edge3 = edge3-len(cut2)
        
        ## store edges
        edges[idx] = np.array([edge0, edge1, edge2, edge3])

        ## assertions
        assert edge0 < edge1 
        if split:
            assert edge2 < edge3

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
        maxhet = superseqs.shape[0]/float(maxhet)

    ## an empty array to fill with failed loci
    hetfilt = np.zeros(superseqs.shape[0], dtype=np.bool)

    ## get the per site number of bases in ambig, and get the max per site 
    ## value for each locus, then get boolean array of those > maxhet.
    for idx, edg in enumerate(edges):
        share = np.array(\
            [np.sum(superseqs[idx, :, edg[0]:edg[1]] == ambig, axis=1)\
            .max() for ambig in "RSKYWM"]).max()
        ## fill filter
        if not share <= maxhet:
            hetfilt[idx] = True

    LOGGER.info("Filtered max_shared_heterozygosity- {}".format(hetfilt.sum()))
    return hetfilt



def make_outfiles(data, samples, force):
    """ Get desired formats from paramsdict and write files to outfiles 
    directory 
    """
    ## Read in the input .loci file that gets transformed into other formats
    ## locifile = os.path.join(data.dirs.outfiles, data.name+".loci")

    ## filter out samples listed as excludes in params.txt
    excludes = (data.paramsdict["excludes"] or [""]) \
             + (data.paramsdict["outgroups"] or [""])
    LOGGER.warn("Excluding these individuals - {}".format(excludes))
    samples = [i for i in samples if i.name not in excludes]

    ## select the output formats. Should we allow some kind of fuzzy matching?
    output_formats = data.paramsdict["output_formats"]
    if "*" in output_formats:
        output_formats = OUTPUT_FORMATS

    ## run each func
    for filetype in output_formats:
        LOGGER.info("Doing - {}".format(filetype))

        # phy & nex come from loci2phynex
        if filetype in ["phy", "nex"]:
            filetype = "phynex"
        # All these file types come from loci2SNP
        elif filetype in ["snps", "usnps", "str", "geno"]:
            filetype = "SNP"

        ## Everything else has its own loci2*.py conversion file.
        ## Get the correct module for this filetype
        ## globals() here gets the module name, and then we can call
        ## the .make() function. This is a little tricky.
        format_module = globals()["loci2"+filetype]

        ## do the call to make the new file format
        format_module.make(data, samples)



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



## File output subfunctions
def write_tmp_loci(data, loci, fname):
    """ Write out the filtered chunk to a tmp file which will be collated by the
    top level run() thread"""

    ## Get longest sample name for pretty printing
    longname_len = max(len(x) for x in data.samples.keys())
    ## Padding distance between name and seq.
    ## This variable is used at least here and in loci2alleles. If you _ever_
    ## have to change this, consider adding it to hackesonly.
    name_padding = 5

    with open(fname.replace("chunk", "loci"), 'w') as outfile:
        for loc in loci:
            for seq in loc:

                ## Backwards compatibility with .loci format which doesn't have 
                ## leading '>' on the snpsites line
                if seq[0] == "//":
                    name = seq[0] 
                else:
                    name = ">" + seq[0]

                name += " " * (longname_len - len(name)+ name_padding)
                outfile.write(name + seq[1] + "\n")



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
