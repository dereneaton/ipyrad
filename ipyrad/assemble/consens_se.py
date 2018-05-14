#!/usr/bin/env python

""" call consensus base calls on single-end data """

# py2/3 compatible
from __future__ import print_function
try:
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

import scipy.stats
import scipy.misc
import ipyrad as ip
import pandas as pd
import numpy as np
import time
import gzip
import glob
import os
from ipyrad.assemble.jointestimate import recal_hidepth
from .util import IPyradError, IPyradWarningExit, clustdealer, PRIORITY
from collections import Counter

# why won't hdf5 just fix this...
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


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
    """
    return probability of base call
    """      
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
    else:
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
        consens = np.array(list(consens), dtype=np.bytes_)
        sep = np.array(arr1.shape[0] * [list(b"nnnn")])
        arrayed = np.hstack([arr1, sep, arr2])

    ## if single-end...
    else:
        consens = np.array(list(cons1), dtype=np.bytes_)
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



def newconsensus(data, sample, tmpchunk, optim):
    """ 
    new faster replacement to consensus 
    """
    ## do reference map funcs?
    isref = "reference" in data.paramsdict["assembly_method"]

    ## temporarily store the mean estimates to Assembly
    data._este = data.stats.error_est.mean()
    data._esth = data.stats.hetero_est.mean()

    ## get number relative to tmp file
    tmpnum = int(tmpchunk.split(".")[-1])

    ## prepare data for reading
    clusters = open(tmpchunk, 'rb')
    pairdealer = izip(*[iter(clusters)] * 2)    
    maxlen = data._hackersonly["max_fragment_length"]

    ## write to tmp cons to file to be combined later
    consenshandle = os.path.join(
        data.dirs.consens, sample.name + "_tmpcons." + str(tmpnum))
    tmp5 = consenshandle.replace("_tmpcons.", "_tmpcats.")
    with h5py.File(tmp5, 'w') as io5:
        io5.create_dataset("cats", (optim, maxlen, 4), dtype=np.uint32)
        io5.create_dataset("alls", (optim, ), dtype=np.uint8)
        io5.create_dataset("chroms", (optim, 3), dtype=np.int64)

        ## local copies to use to fill the arrays
        catarr = io5["cats"][:]
        nallel = io5["alls"][:]
        refarr = io5["chroms"][:]

    ## if reference-mapped then parse the fai to get index number of chroms
    if isref:
        fai = pd.read_csv(data.paramsdict["reference_sequence"] + ".fai", 
            names=['scaffold', 'size', 'sumsize', 'a', 'b'],
            sep="\t")
        faidict = {j: i for i, j in enumerate(fai.scaffold)}

    ## store data for stats counters
    counters = {"name": tmpnum,
                "heteros": 0,
                "nsites": 0,
                "nconsens": 0}

    ## store data for what got filtered
    filters = {"depth": 0,
               "maxh": 0,
               "maxn": 0}

    ## store data for writing
    storeseq = {}

    ## set max limits
    if 'pair' in data.paramsdict["datatype"]:
        maxhet = sum(data.paramsdict["max_Hs_consens"])
        maxn = sum(data.paramsdict["max_Ns_consens"])
    else:
        maxhet = data.paramsdict["max_Hs_consens"][0]
        maxn = data.paramsdict["max_Ns_consens"][0]

    ## load the refmap dictionary if refmapping
    done = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)

        if chunk:
            ## get names and seqs
            piece = chunk[0].decode().strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]

            ## pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]

            ## IF this is a reference mapped read store the chrom and pos info
            ## -1 defaults to indicating an anonymous locus, since we are using
            ## the faidict as 0 indexed. If chrompos fails it defaults to -1
            ref_position = (-1, 0, 0)
            if isref:
                try:
                    ## parse position from name string
                    name, _, _ = names[0].rsplit(";", 2)
                    chrom, pos0, pos1 = name.rsplit(":", 2)
                    
                    ## pull idx from .fai reference dict 
                    chromint = faidict[chrom] + 1
                    ref_position = (int(chromint), int(pos0), int(pos1))
                    
                except Exception as inst:
                    ip.logger.debug(
                        "Reference sequence chrom/pos failed for {}"
                        .format(names[0]))
                    ip.logger.debug(inst)
                    
            ## apply read depth filter
            if nfilter1(data, reps):

                ## get stacks of base counts
                sseqs = [list(seq) for seq in seqs]
                arrayed = np.concatenate(
                    [[seq] * rep for seq, rep in zip(sseqs, reps)]
                    ).astype(np.bytes_)
                arrayed = arrayed[:, :maxlen]
                
                # get consens call for each site, applies paralog-x-site filter
                # consens = np.apply_along_axis(basecall, 0, arrayed, data)
                consens = basecaller(
                    arrayed, 
                    data.paramsdict["mindepth_majrule"], 
                    data.paramsdict["mindepth_statistical"],
                    data._esth, 
                    data._este,
                    )

                ## apply a filter to remove low coverage sites/Ns that
                ## are likely sequence repeat errors. This is only applied to
                ## clusters that already passed the read-depth filter (1)
                if "N" in consens:
                    try:
                        consens, arrayed = removerepeats(consens, arrayed)

                    except ValueError:
                        ip.logger.info("Caught bad chunk w/ all Ns. Skip it.")
                        continue

                ## get hetero sites
                hidx = [i for (i, j) in enumerate(consens) 
                        if j in list("RKSYWM")]
                nheteros = len(hidx)
                
                ## filter for max number of hetero sites
                if nfilter2(nheteros, maxhet):
                    ## filter for maxN, & minlen
                    if nfilter3(consens, maxn):
                        ## counter right now
                        current = counters["nconsens"]
                        ## get N alleles and get lower case in consens
                        consens, nhaps = nfilter4(consens, hidx, arrayed)
                        ## store the number of alleles observed
                        nallel[current] = nhaps

                        ## store a reduced array with only CATG
                        catg = np.array(
                            [np.sum(arrayed == i, axis=0)
                            for i in list("CATG")],
                            dtype='uint32').T
                        catarr[current, :catg.shape[0], :] = catg
                        refarr[current] = ref_position

                        ## store the seqdata for tmpchunk
                        storeseq[counters["name"]] = b"".join(list(consens))
                        counters["name"] += 1
                        counters["nconsens"] += 1
                        counters["heteros"] += nheteros
                    else:
                        #ip.logger.debug("@haplo")
                        filters['maxn'] += 1
                else:
                    #ip.logger.debug("@hetero")
                    filters['maxh'] += 1
            else:
                #ip.logger.debug("@depth")
                filters['depth'] += 1
                
    ## close infile io
    clusters.close()

    ## write final consens string chunk
    if storeseq:
        with open(consenshandle, 'wt') as outfile:
            outfile.write(
                "\n".join([">" + sample.name + "_" + str(key) + \
                "\n" + storeseq[key].decode() for key in storeseq]))

    ## write to h5 array, this can be a bit slow on big data sets and is not 
    ## currently convered by progressbar movement.
    with h5py.File(tmp5, 'a') as io5:
        io5["cats"][:] = catarr
        io5["alls"][:] = nallel
        io5["chroms"][:] = refarr
    del catarr
    del nallel
    del refarr

    ## return stats
    counters['nsites'] = sum([len(i) for i in storeseq.values()])
    return counters, filters



def basecaller(arrayed, mindepth_majrule, mindepth_statistical, estH, estE):
    """
    call all sites in a locus array.
    """

    ## an array to fill with consensus site calls
    cons = np.zeros(arrayed.shape[1], dtype=np.uint8)
    cons.fill(78)
    arr = arrayed.view(np.uint8)
    
    ## iterate over columns
    for col in range(arr.shape[1]):
        ## the site of focus
        carr = arr[:, col]
        
        ## make mask of N and - sites
        mask = carr == 45
        mask += carr == 78        
        marr = carr[~mask]
        
        ## skip if only empties (e.g., N-)
        if not marr.shape[0]:
            cons[col] = 78
        
        ## skip if not variable
        elif np.all(marr == marr[0]):
            cons[col] = marr[0]
        
        ## estimate variable site call
        else:
            ## get allele freqs (first-most, second, third = p, q, r)
            counts = np.bincount(marr)
            
            pbase = np.argmax(counts)
            nump = counts[pbase]
            counts[pbase] = 0

            qbase = np.argmax(counts)
            numq = counts[qbase]
            counts[qbase] = 0

            #rbase = np.argmax(counts)
            #numr = counts[rbase]          # not used
            
            ## based on biallelic depth
            bidepth = nump + numq 
            if bidepth < mindepth_majrule:
                cons[col] = 78
            
            else:
                ## if depth is too high, reduce to sampled int
                if bidepth > 500:
                    base1 = int(500 * (nump / float(bidepth)))
                    base2 = int(500 * (numq / float(bidepth)))
                else:
                    base1 = nump
                    base2 = numq

                ## make statistical base call  
                if bidepth >= mindepth_statistical:
                    ishet, prob = get_binom(base1, base2, estE, estH)
                    #ip.logger.info("ishet, prob, b1, b2: %s %s %s %s", ishet, prob, base1, base2)
                    if prob < 0.95:
                        cons[col] = 78
                    else:
                        if ishet:
                            cons[col] = TRANS[(pbase, qbase)]
                        else:
                            cons[col] = pbase
                
                ## make majrule base call
                else:  # if bidepth >= mindepth_majrule:
                    if nump == numq:
                        cons[col] = TRANS[(pbase, qbase)]
                    else:
                        cons[col] = pbase
    return cons.view("S1")
            


def nfilter1(data, reps):
    """ applies read depths filter """
    if sum(reps) >= data.paramsdict["mindepth_majrule"] and \
        sum(reps) <= data.paramsdict["maxdepth"]:
        return 1
    else:
        return 0



def nfilter2(nheteros, maxhet):
    """ applies max heteros in a seq filter """
    if nheteros <= maxhet:
        return 1
    else:
        return 0



def nfilter3(consens, maxn):
    """ applies filter for maxN and hard minlen (32) """
    ## minimum length for clustering in vsearch
    if consens.size >= 32:
        if consens[consens == "N"].size <= maxn:
            return 1
        else:
            return 0
    else:
        return 0



def nfilter4(consens, hidx, arrayed):
    """ applies max haplotypes filter returns pass and consens"""

    ## if less than two Hs then there is only one allele
    if len(hidx) < 2:
        return consens, 1

    ## store base calls for hetero sites
    harray = arrayed[:, hidx]

    ## remove any reads that have N or - base calls at hetero sites
    ## these cannot be used when calling alleles currently.
    harray = harray[~np.any(harray == "-", axis=1)]
    harray = harray[~np.any(harray == "N", axis=1)]

    ## get counts of each allele (e.g., AT:2, CG:2)
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



def cleanup(data, sample, statsdicts):
    """
    cleaning up. optim is the size (nloci) of tmp arrays
    """
    ip.logger.info("in cleanup for: %s", sample.name)
    isref = 'reference' in data.paramsdict["assembly_method"]

    ## collect consens chunk files
    combs1 = glob.glob(os.path.join(
        data.dirs.consens,
        sample.name + "_tmpcons.*"))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    ## collect tmpcat files
    tmpcats = glob.glob(os.path.join(
        data.dirs.consens,
        sample.name + "_tmpcats.*"))
    tmpcats.sort(key=lambda x: int(x.split(".")[-1]))

    ## get shape info from the first cat, (optim, maxlen, 4)
    with h5py.File(tmpcats[0], 'r') as io5:
        optim, maxlen, _ = io5['cats'].shape

    ## save as a chunked compressed hdf5 array
    handle1 = os.path.join(data.dirs.consens, sample.name + ".catg")
    with h5py.File(handle1, 'w') as ioh5:
        nloci = len(tmpcats) * optim
        dcat = ioh5.create_dataset(
            "catg", 
            (nloci, maxlen, 4),
            dtype=np.uint32,
            chunks=(optim, maxlen, 4),
            compression="gzip")
        dall = ioh5.create_dataset(
            "nalleles", (nloci, ),
            dtype=np.uint8,
            chunks=(optim, ),
            compression="gzip")
        ## only create chrom for reference-aligned data
        if isref:
            dchrom = ioh5.create_dataset(
                "chroms",
                (nloci, 3), 
                dtype=np.int64, 
                chunks=(optim, 3), 
                compression="gzip")

        ## Combine all those tmp cats into the big cat
        start = 0
        for icat in tmpcats:
            io5 = h5py.File(icat, 'r')
            end = start + optim
            dcat[start:end] = io5['cats'][:]
            dall[start:end] = io5['alls'][:]
            if isref:
                dchrom[start:end] = io5['chroms'][:]
            start += optim
            io5.close()
            os.remove(icat)

    ## store the handle to the Sample
    sample.files.database = handle1

    ## record results
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

    ## merge finished consens stats
    for counters, filters in statsdicts:
        ## sum individual counters
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]

    ## merge consens read files
    handle1 = os.path.join(data.dirs.consens, sample.name + ".consens.gz")
    with gzip.open(handle1, 'wt') as out:
        for fname in combs1:
            with open(fname) as infile:
                out.write(infile.read() + "\n")
            os.remove(fname)
    sample.files.consens = [handle1]

    ## set Sample stats_dfs values
    if int(xcounters['nsites']):
        prop = int(xcounters["heteros"]) / float(xcounters['nsites'])
    else:
        prop = 0

    sample.stats_dfs.s5.nsites = int(xcounters["nsites"])
    sample.stats_dfs.s5.nhetero = int(xcounters["heteros"])
    sample.stats_dfs.s5.filtered_by_depth = xfilters['depth']
    sample.stats_dfs.s5.filtered_by_maxH = xfilters['maxh']
    sample.stats_dfs.s5.filtered_by_maxN = xfilters['maxn']
    sample.stats_dfs.s5.reads_consens = int(xcounters["nconsens"])
    sample.stats_dfs.s5.clusters_total = sample.stats_dfs.s3.clusters_total
    sample.stats_dfs.s5.heterozygosity = float(prop)

    ## set the Sample stats summary value
    sample.stats.reads_consens = int(xcounters["nconsens"])

    ## save state to Sample if successful
    if sample.stats.reads_consens:
        sample.stats.state = 5
    else:
        print("No clusters passed filtering in Sample: {}".format(sample.name))
    return sample



def chunk_clusters(data, sample):
    """ split job into bits and pass to the client """

    # counter for split job submission
    num = 0

    # set optim size for chunks in N clusters. The first few chunks take longer
    # because they contain larger clusters, so we create 4X as many chunks as
    # processors so that they are split more evenly.
    optim = int(
        (sample.stats.clusters_total // data.cpus) + \
        (sample.stats.clusters_total % data.cpus))

    ## break up the file into smaller tmp files for each engine
    ## chunking by cluster is a bit trickier than chunking by N lines
    chunkslist = []

    ## open to clusters
    with gzip.open(sample.files.clusters, 'rb') as clusters:
        ## create iterator to sample 2 lines at a time
        pairdealer = izip(*[iter(clusters)] * 2)

        ## Use iterator to sample til end of cluster
        done = 0
        while not done:
            ## grab optim clusters and write to file.
            done, chunk = clustdealer(pairdealer, optim)
            chunk = [i.decode() for i in chunk]
            chunkhandle = os.path.join(
                data.dirs.clusts,
                "tmp_" + str(sample.name) + "." + str(num * optim))
            if chunk:
                chunkslist.append((optim, chunkhandle))
                with open(chunkhandle, 'wt') as outchunk:
                    outchunk.write("//\n//\n".join(chunk) + "//\n//\n")                      
                num += 1

    return chunkslist



def get_subsamples(data, samples, force):
    """
    Apply state, ncluster, and force filters to select samples to be run.
    """

    subsamples = []
    for sample in samples:
        if not force:
            if sample.stats.state >= 5:
                print("""\
    Skipping Sample {}; Already has consens reads. Use force arg to overwrite.\
    """.format(sample.name))
            elif not sample.stats.clusters_hidepth:
                print("""\
    Skipping Sample {}; No clusters found."""
    .format(sample.name, int(sample.stats.clusters_hidepth)))
            elif sample.stats.state < 4:
                print("""\
    Skipping Sample {}; not yet finished step4 """
    .format(sample.name))
            else:
                subsamples.append(sample)

        else:
            if not sample.stats.clusters_hidepth:
                print("""\
    Skipping Sample {}; No clusters found in {}."""
    .format(sample.name, sample.files.clusters))
            elif sample.stats.state < 4:
                print("""\
    Skipping Sample {}; not yet finished step4"""
    .format(sample.name))
            else:
                subsamples.append(sample)

    if len(subsamples) == 0:
        raise IPyradWarningExit("""
    No samples to cluster, exiting.
    """)

    ## if sample is already done skip
    if "hetero_est" not in data.stats:
        print("  No estimates of heterozygosity and error rate. Using default "
              "values")
        for sample in subsamples:
            sample.stats.hetero_est = 0.001
            sample.stats.error_est = 0.0001

    if data._headers:
        print(u"""\
  Mean error  [{:.5f} sd={:.5f}]
  Mean hetero [{:.5f} sd={:.5f}]"""
  .format(data.stats.error_est.mean(), data.stats.error_est.std(),
          data.stats.hetero_est.mean(), data.stats.hetero_est.std()))

    return subsamples



def run(data, samples, force, ipyclient):
    """ checks if the sample should be run and passes the args """
    ## prepare dirs
    data.dirs.consens = os.path.join(data.dirs.project, data.name + "_consens")
    if not os.path.exists(data.dirs.consens):
        os.mkdir(data.dirs.consens)

    ## zap any tmp files that might be leftover
    tmpcons = glob.glob(os.path.join(data.dirs.consens, "*_tmpcons.*"))
    tmpcats = glob.glob(os.path.join(data.dirs.consens, "*_tmpcats.*"))
    for tmpfile in tmpcons + tmpcats:
        os.remove(tmpfile)

    ## filter through samples for those ready
    samples = get_subsamples(data, samples, force)

    ## set up parallel client: how many cores?
    lbview = ipyclient.load_balanced_view()
    data.cpus = data._ipcluster["cores"]
    if not data.cpus:
        data.cpus = len(ipyclient.ids)

    ## wrap everything to ensure destruction of temp files
    inst = ""
    try:
        ## calculate depths, if they changed.
        samples = calculate_depths(data, samples, lbview)

        ## chunk clusters into bits for parallel processing
        lasyncs = make_chunks(data, samples, lbview)

        ## process chunks and cleanup
        process_chunks(data, samples, lasyncs, lbview)

    except KeyboardInterrupt as inst:
        raise inst

    finally:
        ## if process failed at any point delete tmp files
        tmpcons = glob.glob(os.path.join(data.dirs.clusts, "tmp_*.[0-9]*"))
        tmpcons += glob.glob(os.path.join(data.dirs.consens, "*_tmpcons.*"))
        tmpcons += glob.glob(os.path.join(data.dirs.consens, "*_tmpcats.*"))
        for tmpchunk in tmpcons:
            os.remove(tmpchunk)

        ## Finished step 5. Set step 6 checkpoint to 0 to force
        ## re-running from scratch.
        data._checkpoint = 0



def calculate_depths(data, samples, lbview):
    """
    check whether mindepth has changed, and thus whether clusters_hidepth
    needs to be recalculated, and get new maxlen for new highdepth clusts.
    if mindepth not changed then nothing changes.
    """

    ## send jobs to be processed on engines
    start = time.time()
    printstr = ("calculating depths  ", "s5")
    recaljobs = {}
    maxlens = []
    for sample in samples:
        recaljobs[sample.name] = lbview.apply(recal_hidepth, *(data, sample))

    ## block until finished
    while 1:
        ready = [i.ready() for i in recaljobs.values()]
        data._progressbar(len(ready), sum(ready), start, printstr)
        time.sleep(0.1)
        if len(ready) == sum(ready):
            break

    ## check for failures and collect results
    print("")
    modsamples = []
    for sample in samples:
        if not recaljobs[sample.name].successful():
            ip.logger.error("  sample %s failed: %s", 
                sample.name, recaljobs[sample.name].exception())
        else:
            modsample, _, maxlen, _, _ = recaljobs[sample.name].result()
            modsamples.append(modsample)
            maxlens.append(maxlen)

    ## reset global maxlen if something changed
    data._hackersonly["max_fragment_length"] = int(max(maxlens)) + 4

    return samples



def make_chunks(data, samples, lbview):
    """
    calls chunk_clusters and tracks progress.
    """
    ## first progress bar
    start = time.time()
    printstr = ("chunking clusters   ", "s5")
    data._progressbar(10, 0, start, printstr)

    ## send off samples to be chunked
    lasyncs = {}
    for sample in samples:
        lasyncs[sample.name] = lbview.apply(chunk_clusters, *(data, sample))

    ## block until finished
    while 1:
        ready = [i.ready() for i in lasyncs.values()]
        data._progressbar(len(ready), sum(ready), start, printstr)
        time.sleep(0.1)
        if len(ready) == sum(ready):
            break

    ## check for failures
    print("")
    for sample in samples:
        if not lasyncs[sample.name].successful():
            ip.logger.error("  sample %s failed: %s", sample.name, 
                        lasyncs[sample.name].exception())

    return lasyncs



def process_chunks(data, samples, lasyncs, lbview):
    """
    submit chunks to consens func and ...
    """

    ## send chunks to be processed
    start = time.time()
    asyncs = {sample.name: [] for sample in samples}
    printstr = ("consens calling     ", "s5")

    ## get chunklist from results
    for sample in samples:
        clist = lasyncs[sample.name].result()
        for optim, chunkhandle in clist:
            args = (data, sample, chunkhandle, optim)
            asyncs[sample.name].append(lbview.apply_async(newconsensus, *args))
            data._progressbar(10, 0, start, printstr)

    ## track progress
    allsyncs = list(chain(*[asyncs[i.name] for i in samples]))
    while 1:
        ready = [i.ready() for i in allsyncs]
        data._progressbar(len(ready), sum(ready), start, printstr)
        time.sleep(0.1)
        if len(ready) == sum(ready):
            break

    ## get clean samples
    casyncs = {}
    for sample in samples:
        rlist = asyncs[sample.name]
        statsdicts = [i.result() for i in rlist]
        casyncs[sample.name] = lbview.apply(cleanup, *(data, sample, statsdicts))
    while 1:
        ready = [i.ready() for i in casyncs.values()]
        data._progressbar(10, 10, start, printstr)
        time.sleep(0.1)
        if len(ready) == sum(ready):
            break

    ## check for failures:
    print("")
    for key in asyncs:
        asynclist = asyncs[key]
        for rasync in asynclist:
            if not rasync.successful():
                ip.logger.error("  async error: %s \n%s", key, rasync.exception())
    for key in casyncs:
        if not casyncs[key].successful():
            ip.logger.error("  casync error: %s \n%s", key, casyncs[key].exception())

    ## get samples back
    subsamples = [i.result() for i in casyncs.values()]
    for sample in subsamples:
        data.samples[sample.name] = sample

    ## build Assembly stats
    data.stats_dfs.s5 = data._build_stat("s5")

    ## write stats file
    data.stats_files.s5 = os.path.join(data.dirs.consens, 's5_consens_stats.txt')
    with open(data.stats_files.s5, 'w') as out:
        #out.write(data.stats_dfs.s5.to_string())
        data.stats_dfs.s5.to_string(
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
