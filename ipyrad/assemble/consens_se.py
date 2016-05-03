#!/usr/bin/env python2.7

""" call consensus base calls on single-end data """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212


import scipy.stats
import scipy.misc
import itertools
import datetime
import numpy
import time
import h5py
import gzip
import glob
import os
from util import *

from collections import Counter

import logging
LOGGER = logging.getLogger(__name__)


@memoize
def binomprobr(base1, base2, error, het):
    """
    given two bases are observed at a site n1 and n2, and the error rate e, the
    probability the site is truly aa,bb,ab is calculated using binomial 
    distribution as in Li_et al 2009, 2011, and if coverage > 500, 500 
    dereplicated reads were randomly sampled.
    """
    ## major allele freq
    mjaf = base1/float(base1+base2)
    prior_homo = ((1.-het)/2.)
    prior_het = het

    ## get probabilities. Note, b/c only allow two bases, base2 == sum-base1
    try:
        hetro = scipy.misc.comb(base1+base2, base1)/(2.**(base1+base2))
    except OverflowError:
        LOGGER.error("OVERFLOW %s %s", base1, base2)
        hetro = 0.5
    homoa = scipy.stats.binom.pmf(base2, base1+base2, error)
    homob = scipy.stats.binom.pmf(base1, base1+base2, error)

    ## calculate probs
    homoa *= prior_homo
    homob *= prior_homo
    hetro *= prior_het

    ## return 
    probabilities = [homoa, homob, hetro]
    genotypes = ['aa', 'bb', 'ab']
    bestprob = max(probabilities)/float(sum(probabilities))

    return [bestprob, mjaf, genotypes[probabilities.index(max(probabilities))]]



def simpleconsensus(base1, base2):
    """
    majority consensus calling for sites with too low of coverage for
    statistical calling. Only used with 'lowcounts' option. Returns 
    the most common base. Returns consistent alphabetical order for ties.
    """
    #qQn = ['aa','bb','ab']
    maf = base1/(base1+base2)
    return [1.0, maf, 'aa']



@memoize
def hetero(base1, base2):
    """
    returns IUPAC symbol for ambiguity bases, used for polymorphic sites.
    """
    iupac = "N"
    order1 = TRANS.get((base1, base2))
    order2 = TRANS.get((base2, base1))
    if order1:
        iupac = order1
    elif order2:
        iupac = order2
    else:
        ## one or both are Ns
        if [base1, base2].count("N") == 2:
            pass
        elif base1 == "N":
            iupac = base2
        elif base2 == "N":
            iupac = base1
        ## one or both are (-)s            
        else:
            if base1 == "-":
                iupac = base2
            elif base2 == "-":
                iupac = base1
            else:
                LOGGER.error("unrecognized sequence: %s %s", base1, base2)
    return iupac



def removerepeats(consens, arrayed):
    """ 
    Checks for interior Ns in consensus seqs and removes those that are at
    low depth, here defined as less than 1/3 of the average depth. The prop 1/3
    is chosen so that mindepth=6 requires 2 base calls that are not in [N,-].
    """

    ## default trim no edges
    consens = "".join(consens).replace("-", "N")
    edges = [None, None]

    ## trim from left else index starts at zero
    lcons = len(consens)
    consens = consens.lstrip("N")
    edges[0] = lcons - len(consens)

    ## trim from right if nonzero
    lcons = len(consens)
    consens = consens.rstrip("N")
    if lcons - len(consens):
        edges[1] = -1*(lcons - len(consens))

    ## trim same from arrayed
    consens = numpy.array(list(consens))
    arrayed = arrayed[:, edges[0]:edges[1]]

    ## get column counts of Ns and -s
    ndepths = numpy.sum(arrayed == 'N', axis=0) 
    idepths = numpy.sum(arrayed == '-', axis=0)

    ## get proportion of bases that are N- at each site
    nons = ((ndepths + idepths) / float(arrayed.shape[0])) >= 0.75
    ## boolean of whether base was called N
    isn = consens == "N"
    ## make ridx
    ridx = nons * isn

    ## apply filter
    consens = consens[~ridx]
    arrayed = arrayed[:, ~ridx]
    
    return consens, arrayed



def consensus(args):
    """
    from a clust file handle, reads in all copies at a locus and sorts
    bases at each site, tests for errors at the site according to error 
    rate, calls consensus.
    """

    ## unpack args
    data, sample, tmpchunk, optim = args

    ## temporarily store the mean estimates to Assembly
    data._este = data.stats.error_est.mean()
    data._esth = data.stats.hetero_est.mean()

    ## number relative to tmp file
    tmpnum = int(tmpchunk.split(".")[-1])

    ## prepare data for reading
    clusters = open(tmpchunk, 'rb')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## array to store all the coverage data, including consens reads that are 
    ## excluded (for now). The reason we include the low cov data is that this
    ## Assembly might be branched and the new one use a lower depth filter.
    #### dimensions: nreads_in_this_chunk, max_read_length, 4 bases
    maxlen = data._hackersonly["max_fragment_length"]
    if any(x in data.paramsdict["datatype"] for x in ['pair', 'gbs']):
        maxlen *= 2
    catarr = numpy.zeros([optim, maxlen, 4], dtype='uint32')

    ## store data for stats counters
    counters = {"name" : tmpnum,
                "heteros": 0,
                "nsites" : 0,
                "nconsens" : 0}
    ## store data for what got filtered
    filters = {"depth" : 0,
               "heteros" : 0,
               "haplos" : 0,
               "maxn" : 0}

    ## store data for writing
    storeseq = {}

    ## iterate over clusters
    done = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)

        if chunk:
            ## get names and seqs
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            ## pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]

            ## apply read depth filter
            if nfilter1(data, reps):

                ## get stacks of base counts
                sseqs = [list(seq) for seq in seqs]
                arrayed = numpy.concatenate(
                          [[seq]*rep for seq, rep in zip(sseqs, reps)])

                ## get consens call for each site, applies paralog-x-site filter
                consens = numpy.apply_along_axis(basecall, 0, arrayed, data)

                ## apply a filter to remove low coverage sites/Ns that
                ## are likely sequence repeat errors. This is only applied to 
                ## clusters that already passed the read-depth filter (1)
                if "N" in consens:
                    try:
                        consens, arrayed = removerepeats(consens, arrayed)

                    except ValueError as _:
                        LOGGER.info("Caught a bad chunk w/ all Ns. Skip it.")
                        continue

                ## get hetero sites
                hidx = [i for (i, j) in enumerate(consens) \
                            if j in list("RKSYWM")]
                nheteros = len(hidx)

                ## filter for max number of hetero sites
                if nfilter2(data, nheteros):
                    ## filter for maxN, & minlen 
                    if nfilter3(data, consens):
                        ## filter for max alleles and get lower case in consens
                        consens, passed = nfilter4(data, consens, hidx, arrayed)
                        if passed:
                            ## store a reduced array with only CATG
                            catg = numpy.array(\
                       [numpy.sum(arrayed == i, axis=0) for i in list("CATG")], 
                       dtype='uint32').T
                            catarr[counters["nconsens"]][:catg.shape[0]] = catg
                            ## store data for tmpchunk
                            storeseq[counters["name"]] = "".join(list(consens))
                            counters["name"] += 1
                            counters["nconsens"] += 1
                            counters["heteros"] += nheteros                            
                        else:
                            #LOGGER.debug("@maxn")
                            filters['haplos'] += 1
                    else:
                        #LOGGER.debug("@haplo")
                        filters['maxn'] += 1
                else:
                    #LOGGER.debug("@hetero")
                    filters['heteros'] += 1
            else:
                #LOGGER.debug("@depth")
                filters['depth'] += 1
    ## close file io
    clusters.close()

    ## write to tmp cons to file to be combined later
    consenshandle = os.path.join(data.dirs.consens, 
                                 sample.name+"_tmpcons."+str(tmpnum))
    #LOGGER.info('writing %s', consenshandle)
    #LOGGER.info('passed in this chunk: %s', len(storeseq))
    #LOGGER.info('caught in this chunk: %s', filters)
    if storeseq:
        with open(consenshandle, 'wb') as outfile:
            outfile.write("\n".join([">"+sample.name+"_"+str(key)+"\n"+\
                                   str(storeseq[key]) for key in storeseq]))

    ## save tmp catg array that will be combined into hdf5 later
    with open(consenshandle.replace("_tmpcons.", "_tmpcats."), 'wb') as dumph:
        numpy.save(dumph, catarr)

    ## final counts and return
    counters['nsites'] = sum([len(i) for i in storeseq.itervalues()])        

    return counters, filters



def nfilter1(data, reps):
    """ applies read depths filter """
    if sum(reps) >= data.paramsdict["mindepth_majrule"] and \
        sum(reps) <= data.paramsdict["maxdepth"]:
        return 1
    else:
        return 0



def nfilter2(data, nheteros):
    """ applies max heteros in a seq filter """
    if nheteros <= sum(data.paramsdict["max_Hs_consens"]):
        return 1
    else:
        return 0



def nfilter3(data, consens):
    """ applies filter for maxN and hard minlen (32) """
    ## minimum length for clustering in vsearch
    if consens.size >= 32:
        if consens[consens == "N"].size <= \
                            sum(data.paramsdict["max_Ns_consens"]):
            return 1
        else:
            return 0
    else:
        return 0



def nfilter4(data, consens, hidx, arrayed):
    """ applies max haplotypes filter returns pass and consens"""

    ## only apply if >1 hetero site
    if len(hidx) < 2:
        return consens, 1

    ## only apply if user arg is to apply diploid filter
    if data.paramsdict["max_alleles_consens"] != 2:
        return consens, 1

    ## store base calls for hetero sites
    harray = arrayed[:, hidx]

    ## remove any reads that have N or - base calls at hetero sites
    ## these cannot be used when calling alleles currently.
    harray = harray[~numpy.any(harray == "-", axis=1)]    
    harray = harray[~numpy.any(harray == "N", axis=1)]

    ## get counts of each allele (e.g., AT:2, CG:2)
    ccx = Counter([tuple(i) for i in harray])

    ## Two possibilities we would like to distinguish, but we can't. Therefore, 
    ## we just throw away low depth third alleles that are within seq. error. 
    ## 1) a third base came up as a sequencing error but is not a unique allele
    ## 2) a third or more unique allele is there but at low frequency

    ## remove low freq alleles if more than 2, since they may reflect 
    ## sequencing errors instead of true heteros
    if len(ccx) > 2:
        totdepth = harray.shape[0]
        cutoff = max(1, totdepth // 10)
        alleles = [i for i in ccx if ccx[i] > cutoff]
    else:
        alleles = ccx.keys()

    ## how many high depth alleles?
    nalleles = len(alleles)

    ## if 2 alleles then go ahead to find alleles
    if nalleles > data.paramsdict["max_alleles_consens"]:
        return consens, 0
    else:
        if nalleles == 2:
            try:
                consens = storealleles(consens, hidx, alleles)
                return consens, 1 
            except (IndexError, KeyError):
                ## the H sites do not form good alleles
                return consens, 0
        else:
            return consens, 0            



def storealleles(consens, hidx, alleles):
    """ store phased allele data for diploids """
    ## find the first hetero site and choose the priority base
    ## example, if W: then priority base in T and not A. PRIORITY=(order: CATG)
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



def basecall(rsite, data):
    """ prepares stack for making base calls """
    ## count em
    site = Counter(rsite)

    ## remove Ns and (-)s
    if "N" in site:
        site.pop("N")
    if "-" in site:
        site.pop("-")

    ## get the two most common alleles
    if site:
        base1 = base2 = 0
        comms = site.most_common(2)
        base1 = comms[0][1]
        if len(comms) > 1:
            base2 = comms[1][1]

        ## if site depth after removing Ns, (-s) and third bases is below limit
        bidepth = base1 + base2
        if bidepth < data.paramsdict["mindepth_majrule"]:
            cons = "N"

        ## speedhack: make the base call using a method depending on depth
        ## if highdepth and invariable just call the only base
        elif (bidepth > 10) and (not base2):
            cons = comms[0][0]

        else:
            ## if depth > 500 reduce to <500 at same proportion to avoid 
            ## large memerror in scipy.misc.comb function
            if bidepth >= 500: 
                sbase1 = int(500 * (base1 / float(base1)))
                sbase2 = int(500 * (base2 / float(base1)))
            else:
                sbase1 = base1
                sbase2 = base2

            ## And then use basecaller
            cons = basecaller(data, sbase1, sbase2, comms)
    else:
        cons = "N"

    return cons



def basecaller(data, base1, base2, comms):
    """ inputs data to binomprobr and gets alleles correctly oriented """

    ## make statistical base call
    if base1+base2 >= data.paramsdict["mindepth_statistical"]:
        prob, _, who = binomprobr(base1, base2, data._este, data._esth)

    elif base1+base2 >= data.paramsdict["mindepth_majrule"]:
        prob, _, who = simpleconsensus(base1, base2)

    #else:
    #    LOGGER.error("gap in mindepth settings")

    ## if the base could be called with 95% probability
    if float(prob) >= 0.95:
        if who != "ab":
            ## site is homozygous
            cons = comms[0][0]
        else:
            ## site is heterozygous
            cons = hetero(*[i[0] for i in comms])
            #LOGGER.info("GIVE: %s, GET: %s", [i[0] for i in comms], cons)
    else:
        cons = "N"
    return cons



def cleanup(args):
    """ cleaning up. optim is the size (nloci) of tmp arrays """

    ## parse args list
    data, sample, statsdicts = args
    LOGGER.info("in cleanup for: %s", sample.name)

    ## collect consens chunk files
    combs1 = glob.glob(os.path.join(
                        data.dirs.consens,
                        sample.name+"_tmpcons.*"))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    ## collect tmpcat files
    tmpcats = glob.glob(os.path.join(
                        data.dirs.consens,
                        sample.name+"_tmpcats.*"))
    tmpcats.sort(key=lambda x: int(x.split(".")[-1]))

    ## get shape info from the first cat, they're all the same size arrays
    with open(tmpcats[0]) as cat:
        catg = numpy.load(cat)
    ## (optim, maxlen, 4)
    optim, maxlen, _ = catg.shape

    ## save as a chunked compressed hdf5 array
    handle1 = os.path.join(data.dirs.consens, sample.name+".catg")
    ioh5 = h5py.File(handle1, 'w')
    nloci = len(tmpcats) * optim
    dset = ioh5.create_dataset("catg", (nloci, maxlen, 4), 
                               dtype=numpy.uint32,
                               chunks=(optim, maxlen, 4),
                               compression="gzip")

    ## Combine all those tmp cats into the big cat
    start = 0
    for icat in tmpcats:
        icatg = numpy.load(icat)
        end = start + optim        
        dset[start:end] = icatg
        start += optim
        os.remove(icat)
    ioh5.close()

    ## store the handle to the Sample
    sample.files.database = handle1

    ## record results
    xcounters = {"nconsens": 0,
                 "heteros": 0, 
                 "nsites": 0}
    xfilters = {"depth": 0, 
               "heteros": 0,
               "haplos": 0,
               "maxn": 0}

    ## merge finished consens stats
    for i in range(len(combs1)):
        counters, filters = statsdicts[i]
        ## sum individual counters
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]

    ## merge consens read files
    handle1 = os.path.join(data.dirs.consens, sample.name+".consens.gz")
    with gzip.open(handle1, 'wb') as out:
        for fname in combs1:
            with open(fname) as infile:
                out.write(infile.read()+"\n")
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
    sample.stats_dfs.s5.filtered_by_maxH = xfilters['haplos']    
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



def run_full(data, sample, lbview):
    """ split job into bits and pass to the client """

    ## counter for split job submission
    num = 0

    ## set optim size for chunks in N clusters. The first few chunks take longer
    ## because they contain larger clusters, so we create 4X as many chunks as
    ## processors so that they are split more evenly.
    ncpus = detect_cpus()
    optim = int((sample.stats.clusters_total // ncpus) + \
                (sample.stats.clusters_total % ncpus))
    #LOGGER.info("ncpus %s, optim %s", ncpus, optim)

    ## break up the file into smaller tmp files for each engine
    ## chunking by cluster is a bit trickier than chunking by N lines
    chunkslist = []

    ## open to clusters
    clusters = gzip.open(sample.files.clusters, 'rb')
    ## create iterator to sample 2 lines at a time
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## Use iterator to sample til end of cluster
    done = 0
    while not done:
        ## grab optim clusters and write to file. Clustdealer breaks by clusters
        done, chunk = clustdealer(pairdealer, optim)
        chunkhandle = os.path.join(data.dirs.clusts, 
                                   "tmp_"+str(sample.name)+"."+str(num*optim))
        if chunk:
            chunkslist.append(chunkhandle)            
            with open(chunkhandle, 'wb') as outchunk:
                outchunk.write("//\n//\n".join(chunk)+"//\n//\n")
            num += 1

    ## close clusters handle
    clusters.close()

    ## send chunks across engines, will delete tmps if failed
    asyncs = []
    for chunkhandle in chunkslist:
        ## used to increment names across processors
        args = [data, sample, chunkhandle, optim]
        async = lbview.apply_async(consensus, args)
        asyncs.append(async)

    return asyncs



def run(data, samples, force, ipyclient):
    """ checks if the sample should be run and passes the args """

    ## prepare dirs
    data.dirs.consens = os.path.join(data.dirs.project, data.name+"_consens")
    if not os.path.exists(data.dirs.consens):
        os.mkdir(data.dirs.consens)

    ## zap any tmp files that might be leftover
    tmpcons = glob.glob(os.path.join(data.dirs.consens, "*_tmpcons.*"))
    tmpcats = glob.glob(os.path.join(data.dirs.consens, "*_tmpcats.*"))
    for tmpfile in tmpcons+tmpcats:
        os.remove(tmpfile)

    ## filter through samples for those ready
    subsamples = []
    for sample in samples:
        if not force:
            if sample.stats.state >= 5:
                print("""\
    Skipping Sample {}; Already has consens reads. Use force arg to overwrite.\
    """.format(sample.name))
            elif not sample.stats.clusters_hidepth:
                print("""\
    Skipping Sample {}; No clusters found."""\
    .format(sample.name, int(sample.stats.clusters_hidepth)))
            elif sample.stats.state < 4:
                print("""\
    Skipping Sample {}; not yet finished step4 """)
            else:
                subsamples.append(sample)
                
        else:
            if not sample.stats.clusters_hidepth:
                print("""\
    Skipping Sample {}; No clusters found in {}."""\
    .format(sample.name, sample.files.clusters))
            elif sample.stats.state < 4:
                print("""\
    Skipping Sample {}; not yet finished step4"""\
    .format(sample.name))
            else:
                subsamples.append(sample)

    ## if sample is already done skip
    if "hetero_est" not in data.stats:
        print("  No estimates of heterozygosity and error rate. Using default "\
              "values")
        for sample in subsamples:
            sample.stats.hetero_est = 0.001
            sample.stats.error_est = 0.0001

    if data._headers:
        if data.paramsdict["max_alleles_consens"] == 1:
            print("  Haploid calls with paralog filter [max alleles = 1]")
        elif data.paramsdict["max_alleles_consens"] == 2:
            print("  Diploid calls with paralog filter [max alleles = 2]")
        elif data.paramsdict["max_alleles_consens"] > 2:
            print("  Diploid base calls with relaxed paralog filter "\
                    "(max alleles = {})".\
                    format(data.paramsdict["max_alleles_consens"]))
        print(u"""\
  Mean error  [{:.5f} sd={:.5f}]
  Mean hetero [{:.5f} sd={:.5f}]"""\
  .format(data.stats.error_est.mean(), data.stats.error_est.std(), 
          data.stats.hetero_est.mean(), data.stats.hetero_est.std()))

    ## start a client
    start = time.time()
    lbview = ipyclient.load_balanced_view()

    ## store asyncs for consens call functions
    lasyncs = {}

    ## first progress bar 
    elapsed = datetime.timedelta(seconds=int(time.time()-start))                        
    progressbar(10, 0, " consensus calling     | {}".format(elapsed))

    ## send off jobs to be processed
    njobs = 0
    for sample in subsamples:
        ## make chunks and submit them to apply queue
        lasyncs[sample.name] = run_full(data, sample, lbview)
        njobs += len(lasyncs[sample.name])

        ## print progress post-slice
        elapsed = datetime.timedelta(seconds=int(time.time()-start))                        
        progressbar(10, 0, " consensus calling     | {}".format(elapsed))

    ## create a waiting object
    tmpids = lbview.history
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1)

    try:
        while 1:
            if not res.ready():
                time.sleep(1)

                ## count how many finished for the progress bar
                nfinished = 0
                for joblist in lasyncs.values():
                    nfinished += sum([i.ready() for i in joblist])

                ## print progress bars while we wait
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(njobs, nfinished,
                    " consensus calling     | {}".format(elapsed))

            else:
                break

        ## get clean samples
        for sample in subsamples:
            statsdicts = [i.get() for i in lasyncs[sample.name]]
            cleanup([data, data.samples[sample.name], statsdicts])

        ## build Assembly stats
        data.stats_dfs.s5 = data.build_stat("s5")

        ## write stats file
        data.stats_files.s5 = os.path.join(data.dirs.consens, 
                                               's5_consens_stats.txt')
        with open(data.stats_files.s5, 'w') as out:
            out.write(data.stats_dfs.s5.to_string())

    finally:
        ## if process failed at any point delete tmp files
        tmpcons = glob.glob(os.path.join(data.dirs.clusts, "tmp_*.[0-9]*"))
        tmpcons += glob.glob(os.path.join(data.dirs.consens, "*_tmpcons.*"))
        tmpcons += glob.glob(os.path.join(data.dirs.consens, "*_tmpcats.*"))
        for tmpchunk in tmpcons:
            os.remove(tmpchunk)

        progressbar(20, 20, " consensus calling     | {}".format(elapsed))            
        #if data._headers:
        print("")




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
                         ROOT, "tests", "test_pairgbs", "test_pairgbs"))
    TEST.step5(force=True)
    print(TEST.stats)

    ## run test on rad data1
    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "test_rad", "data1"))
    TEST.step5(force=True)
    print(TEST.stats)

    ## run test on empirical pairgbs data1
    # TEST = ip.load.load_assembly(os.path.join(\
    #         "/home/deren/Documents/longi_test_ip", "longi"))
    # TEST.step5(force=True)
    # print(TEST.stats)
