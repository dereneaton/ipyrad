#!/usr/bin/env python2.7

""" call consensus base calls on single-end data """

from __future__ import print_function
# pylint: disable=E1101

import scipy.stats
import scipy.misc
import itertools
import time
import numpy
import gzip
import glob
import os
from .util import *

from collections import Counter

import logging
LOGGER = logging.getLogger(__name__)



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
    hetro = scipy.misc.comb(base1+base2, base1)/(2.**(base1+base2))
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


def hetero(base1, base2):
    """
    returns IUPAC symbol for ambiguity bases, used for polymorphic sites.
    """
    iupac = "N"
    trans = {('G', 'A'):"R",
             ('G', 'T'):"K",
             ('G', 'C'):"S",
             ('T', 'C'):"Y",
             ('T', 'A'):"W",
             ('C', 'A'):"M"}
    order1 = trans.get((base1, base2))
    order2 = trans.get((base2, base1))
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


def removerepeat_Ns(shortcon, stacked):
    """ checks for interior Ns in consensus seqs and removes those that arise 
    next to *single repeats* of at least 3 bases on either side, which may be
    sequencing errors on deep coverage repeats """
    ### what is the most common site depth
    #comdepth = Counter([sum(i.values()) for i in stacked]).most_common(1)[0][1]
    #print comdepth, 'comdepth'

    ### below this depth it is likely a sequencing repeat
    #mindepth = comdepth/10.
    #print mindepth, 'mindepth'

    nsites = [i for i, j in enumerate(shortcon) if j == "N"]
    repeats = set()

    for nsite in nsites:
        ## grab the five bases to the left of this N site
        isvar = len(set(list(shortcon)[nsite-5:nsite]))
        if isvar < 2:
            ## could use mindepth here to increase stringency
            repeats.add(nsite)
            LOGGER.info("N repeat - left repeat")            
        
        ## grab the five bases to the right of this n site
        isvar2 = len(set(list(shortcon)[nsite+1:nsite+6]))
        if isvar2 < 2:
            ## could use mindepth here to increase stringency            
            repeats.add(nsite)
            LOGGER.info("N repeat - right repeat")

        ## check if most bases were - at this site
        if stacked[nsite].get("-"):
            if stacked[nsite].get("-")/sum(stacked[nsite].values()) > 0.5:
                repeats.add(nsite)
                LOGGER.info("N repeat - removal")

    ## remove repeat sites from shortcon and stacked
    shortcon = [j for (i, j) in enumerate(shortcon) if i not in repeats]
    stacked = stacked[:len(shortcon)]
    stacked = [j for (i, j) in enumerate(stacked) if i not in repeats]

    return "".join(shortcon), stacked



def consensus(args):
    """
    from a clust file handle, reads in all copies at a locus and sorts
    bases at each site, tests for errors at the site according to error 
    rate, calls consensus.
    """

    ## unpack args
    LOGGER.info(args)    
    data, sample, tmpchunk, optim = args

    ## number relative to tmp file
    tmpnum = int(tmpchunk.split(".")[-1])

    ## prepare data for reading
    clusters = open(tmpchunk, 'rb')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## array to store all the coverage data, including
    ## the consens reads that are excluded (for now)
    #nreads = sample.stats.clusters_kept
    ## dimensions: nreads, max_read_length, 4 bases
    ## will store counts of bases only ne
    catarr = numpy.empty([optim, 210, 4], dtype='int16')

    ## counters 
    counters = Counter()
    counters["name"] = tmpnum
    counters["heteros"] = 0
    counters["nsites"] = 0
    counters["nconsens"] = 0    

    filters = Counter()
    filters["depth"] = 0    
    filters["heteros"] = 0
    filters["haplos"] = 0
    filters["maxn"] = 0

    ## store data for writing
    storeseq = {}
    storecat = []

    ## iterate over clusters
    done = 0
    while not done:
        done, chunk = clustdealer(pairdealer, 1)
        if chunk:
            ## get names and seqs
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            reps = [int(sname.split(";")[-2][5:]) for sname in names]

            ## apply read depth filter
            if nfilter1(data, reps):

                ## get stacks of base counts
                sseqs = [list(seq) for seq in seqs]
                arrayed = numpy.concatenate(
                            [[seq]*rep for seq, rep in zip(sseqs, reps)])
                stacked = [Counter(seq) for seq in arrayed.T] 

                ## get consens call for each site, paralog site filter
                #LOGGER.info('arr %s', arrayed)
                consens = numpy.apply_along_axis(basecall, 0, arrayed, data)
                ## get hetero sites
                heteros = [i[0] for i in enumerate(consens) \
                                 if i[1] in list("RKSYWM")]
                nheteros = len(heteros)
                counters["heteros"] += nheteros
                if nheteros > 5:
                    LOGGER.info("NH: %s", nheteros)
                    LOGGER.debug("ARR: %s", arrayed)
                    LOGGER.debug("cons: %s", consens)                    
                ## filter for max number of hetero sites
                if nfilter2(data, nheteros):

                    ## filter for max number of haplotypes
                    mpl = 1
                    if nheteros > 1:
                        consens, mpl = nfilter3(data, consens, heteros, 
                                                seqs, reps)

                    ## if the locus passed paralog filtering
                    if mpl:
                        consens = "".join(consens).replace("-", "N")
                        shortcon = consens.rstrip("N")
                        ## this function which removes low coverage sites next 
                        ## to poly repeats that are likely sequencing errors 
                        ## TODO: Until this func is optimized
                        if "N" in shortcon:
                            LOGGER.debug("preN: %s", shortcon)
                        shortcon, stacked = removerepeat_Ns(shortcon, stacked)

                        if shortcon.count("N") <= \
                                       sum(data.paramsdict["max_Ns_consens"]):
                            ## minimum length for clustering in vsearch
                            if len(shortcon) >= 32:
                                ## store sequence

                                ## store a reduced array with only CATG
                                catg = numpy.array([
                                    [i["C"] for i in stacked], 
                                    [i["A"] for i in stacked], 
                                    [i["T"] for i in stacked],
                                    [i["G"] for i in stacked]], 
                                    dtype='int16').T
                                catarr[counters["nconsens"]][:catg.shape[0]] = \
                                catg
                                #storecat.append(catg)
                                storeseq[counters["name"]] = shortcon
                                counters["name"] += 1
                                counters["nconsens"] += 1                                

                            else:
                                LOGGER.info("@shortmaxn")
                                filters['maxn'] += 1
                        else:
                            LOGGER.info("@maxn")
                            filters['maxn'] += 1
                    else:
                        LOGGER.info("@haplo")
                        filters['haplos'] += 1
                else:
                    LOGGER.info("@hetero")
                    filters['hetero'] += 1
            else:
                LOGGER.info("@depth")
                filters['depth'] += 1

    clusters.close()

    consenshandle = os.path.join(data.dirs.consens, 
                                 sample.name+"_tmpcons."+str(tmpnum))
    ## write to file
    with open(consenshandle, 'wb') as outfile:
        outfile.write("\n".join([">"+sample.name+"_"+str(key)+"\n"+\
                                 storeseq[key] for key in storeseq]))
    #storecat = numpy.array(storecat)

    with open(consenshandle.replace("_tmpcons.", "_tmpcats."), 'wb') as dumph:
        numpy.save(dumph, catarr)#storecat)

    ## count the number of polymorphic sites 
    ## TODO: make a decision about the pair separator
    # if 'ddrad' in data.paramsdict["datatype"]:
    #     if 'pair' in data.paramsdict["datatype"]:
    #         sub = 4 + len(data.paramsdict["restriction_overhang"][0])
    #     else:
    #         sub = len(data.paramsdict["restriction_overhang"][0])
    # elif 'gbs' in data.paramsdict["datatype"]:
    #     if 'pair' in data.paramsdict["datatype"]:
    #         sub = 4 + len(data.paramsdict["restriction_overhang"][0])*2 
    #         #  (len(params["cut"])*2)
    #     else:
    #         sub = len(data.paramsdict["restriction_overhang"][0])
    # elif data.paramsdict["datatype"] == "merged":
    #     sub = len(data.paramsdict["restriction_overhang"][0])*2 
    # else:
    #     sub = len(data.paramsdict["restriction_overhang"][0])

    ## final counts and return
    counters['nsites'] = sum([len(i) for i in storeseq.itervalues()])
    return counters, filters



def nfilter1(data, reps):
    """ applies read depths filter """
    if sum(reps) >= data.paramsdict["mindepth_majrule"] and \
        sum(reps) <= data.paramsdict["max_stack_size"]:
        return 1
    else:
        return 0


def nfilter2(data, nheteros):
    """ applies max heteros in a seq filter """
    if nheteros <= data.paramsdict["max_Hs_consens"]:
        return 1
    else:
        return 0


def nfilter3(data, consens, heteros, seqs, reps):
    """ applies max haplotypes filter returns pass and consens"""

    ## store 
    ordered = []

    ## get each unique order of polymorphic sites
    for read, rep in zip(seqs, reps):
        orderedpolys = []
        ## exclude low depth haplos
        if rep/float(sum(reps)) > 0.10:
            for hsite in heteros:
                orderedpolys.append(read[hsite])
            ## exclude with Ns or (-)s
            if not any([i in orderedpolys for i in ["-", "N"]]):
                ordered.extend(["".join(orderedpolys)]*rep)

    counted = Counter(ordered)
    nalleles = len(counted.keys())
    LOGGER.debug("heteros locs %s", heteros)
    # LOGGER.debug("orderedpolys %s", orderedpolys)
    # LOGGER.debug("ordered %s", ordered)
    LOGGER.debug('counted %s', counted)
    # LOGGER.info('%s alleles', len(counted.keys()))

    ## how many high depth alleles?
    if nalleles > data.paramsdict["ploidy"]:
        return consens, 0
    else:
        ## if diploid try to get two alleles and return consens with
        ## lower upper and lower casing to save phased info
        try:
            if (nalleles > 1) and (data.paramsdict["ploidy"] == 2):
                consens = findalleles(consens, heteros, counted)
        except IndexError as inst:
            LOGGER.error("nfilter3 error again: %s", inst)
            LOGGER.error("\n"+"\n".join(seqs))
        return consens, 1



def findalleles(consens, heteros, counted):
    """ store phased allele data for diploids """
    ## find the first hetero site from the left
    alleles = [i[0] for i in counted.most_common(2)]
    LOGGER.debug('INSIDE alleles %s', alleles)
    firstbigbase = uplow(unhetero(consens[heteros[0]]))
    LOGGER.debug('INSIDE firstbigbase %s', firstbigbase)
    query = [i for i in alleles if i[0] == firstbigbase][0]
    LOGGER.debug('INSIDE query %s', query)

    for idx, site in enumerate(heteros[1:]):
        if query[idx] != uplow(unhetero(consens[site])):
            consens[site] = consens[site].lower()
            LOGGER.debug('newlow %s', consens[site])
    return consens            

    #     uplow(unhetero())
    #     lastbig = uplow(unhetero(heteros[idx-1]))
    #     thisbig = uplow(unhetero(heteros[idx]))
    #     if thisbase

    # for ind, hsite in enumerate(heteros[1:]):
    #     LOGGER.info('what %s, %s, %s', ind, hsite, consens[hsite])
    #     bigbase = uplow(unhetero(consens[hsite]))
    #     ## does this bigbase match the previous bigbase?
    #     ## if not iupac is small
    #     if query[ind-1] != bigbase:
    #         consens[hsite] = consens[hsite].lower()

    # return consens



def basecall(site, data):
    """ prepares stack for making base calls """

    site = Counter(site)

    ## remove Ns and (-)s
    if "N" in site:
        site.pop("N")
    if "-" in site:
        site.pop("-")

    #assert site, "site is empty, need a catch for this..."

    ## get the most common alleles
    base1 = base2 = 0
    comms = site.most_common()
    base1 = comms[0][1]
    if len(comms) > 1:
        base2 = comms[1][1]
        
    ## if site depth after removing Ns, (-s) and third bases is below limit
    bidepth = base1 + base2
    if bidepth < data.paramsdict["mindepth_majrule"]:
        cons = "N"
    else:
        ## if depth > 500 reduce to randomly sampled 500 
        if bidepth >= 500: 
            randomsample = numpy.array(tuple("A"*base1+"B"*base2))
            numpy.random.shuffle(randomsample)
            base1 = list(randomsample[:500]).count("A")
            base2 = list(randomsample[:500]).count("B")

        ## make the base call using a method depending on depth
        ## if highdepth and invariable just call the only base
        if (bidepth > 10) and (not base2):
            cons = comms[0][0]
        ## but if variable then use basecaller
        else:
            cons = basecaller(data, site, base1, base2)
    return cons



def basecaller(data, site, base1, base2):
    """ inputs data to binomprobr and gets alleles correctly oriented """

    ## make statistical base call
    if base1+base2 >= data.paramsdict["mindepth_statistical"]:
        prob, _, who = binomprobr(base1, base2, 
                                  data.stats.error_est.mean(),
                                  data.stats.hetero_est.mean())

    elif base1+base2 >= data.paramsdict["mindepth_majrule"]:
        prob, _, who = simpleconsensus(base1, base2)

    else:
        LOGGER.error("gap in mindepth settings")

    ## if the base could be called with 95% probability
    if float(prob) >= 0.95:
        if who != "ab":
            ## site is homozygous
            cons = site.most_common(1)[0][0]
        else:
            ## site is heterozygous
            cons = hetero(*[i[0] for i in site.most_common(2)])
    else:
        cons = "N"
    return cons



def clustdealer(pairdealer, optim):
    """ return optim clusters given iterators, and whether it got all or not"""
    ccnt = 0
    chunk = []
    while ccnt < optim:
        ## try refreshing taker, else quit
        try:
            taker = itertools.takewhile(lambda x: x[0] != "//\n", pairdealer)
            oneclust = ["".join(taker.next())]
        except StopIteration:
            #LOGGER.debug('last chunk %s', chunk)
            return 1, chunk

        ## load one cluster
        while 1:
            try: 
                oneclust.append("".join(taker.next()))
            except StopIteration:
                break
        chunk.append("".join(oneclust))
        ccnt += 1
    return 0, chunk



def cleanup(data, sample, statsdicts):
    """ cleaning up """

    ## rejoin chunks
    combs1 = glob.glob(os.path.join(
                        data.dirs.consens,
                        sample.name+"_tmpcons.*"))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    ## record results
    xcounters = {"nconsens": 0,
                 "heteros": 0, 
                 "nsites": 0}
    xfilters = {"depth": 0, 
               "heteros": 0,
               "haplos": 0,
               "maxn": 0}

    ## merge catg files
    cats1 = glob.glob(os.path.join(
                      data.dirs.consens,
                      sample.name+"_tmpcats.*"))
    cats1.sort(key=lambda x: int(x.split(".")[-1]))
    handle1 = os.path.join(data.dirs.consens, sample.name+".catg")
    with open(handle1, 'wb') as outcat:
        ## open first cat
        with open(cats1[0]) as cat:
            catg = numpy.load(cat)
        ## extend with other cats
        for icat in cats1[1:]:
            icatg = numpy.load(icat)
            catg = numpy.concatenate([catg, icatg])
            os.remove(icat)
        numpy.save(outcat, catg)
        os.remove(cats1[0])

    ## merge finished consens stats
    for i in range(len(combs1)):
        counters, filters = statsdicts[i]
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]
    sample.stats.reads_consens = xcounters["nconsens"]

    ## merge consens read files
    handle1 = os.path.join(data.dirs.consens, sample.name+".consens.gz")
    with gzip.open(handle1, 'wb') as out:
        for fname in combs1:
            with open(fname) as infile:
                out.write(infile.read()+"\n")
            os.remove(fname)
    sample.files.consens = [handle1]

    ## find longest name to make printing code block
    data.statsfiles.s5 = os.path.join(data.dirs.consens, 's5_consens.txt')    
    longestname = max([len(i) for i in data.samples.keys()])
    printblock = "{:<%d} {:>11} {:>11} {:>11} {:>11} " % (longestname + 4) \
                +"{:>11} {:>11} {:>11} {:>11} {:>11}\n"
    if not os.path.exists(data.statsfiles.s5):
        with open(data.statsfiles.s5, 'w') as outfile:
            outfile.write(printblock.format("sample", "nclusters", 
                "depthfilter", "maxHfilter", "haplofilter", "maxNfilter", 
                "nconsensus", "nsites", "nhetero", "hetero"))

    ## append stats to file
    outfile = open(data.statsfiles.s5, 'a+')
    try:
        prop = xcounters["heteros"]/float(xcounters['nsites'])
    except ZeroDivisionError: 
        prop = 0
    ## redefine printblock to allow for floats
    printblock = "{:<%d} {:>11} {:>11} {:>11} {:>11} " % (longestname + 4) \
                +"{:>11} {:>11} {:>11} {:>11} {:>11.5f}\n"
    outfile.write(printblock.format(
        sample.name, 
        int(sample.stats.clusters_kept),
        int(sample.stats.clusters_kept - xfilters['depth']),
        int(sample.stats.clusters_kept - xfilters['depth'] - \
            xfilters['heteros']),
        int(sample.stats.clusters_kept - xfilters['depth'] - \
            xfilters['heteros'] - xfilters['haplos']),
        int(sample.stats.clusters_kept - xfilters['depth'] - \
            xfilters['heteros'] - xfilters['haplos'] - xfilters['maxn']),
        int(sample.stats.reads_consens),
        xcounters["nsites"],
        xcounters["heteros"],
        prop)
    )

    outfile.close()

    # ## save stats to Sample if successful
    if sample.stats.reads_consens:
        sample.stats.state = 5
        ## save stats to data
        data._stamp("s5 consensus base calling on "+sample.name)

    else:
        print("No clusters passed filtering in Sample: {}".format(sample.name))



def run_full(data, sample, ipyclient):
    """ split job into bits and pass to the client """

    ## counter for split job submission
    num = 0

    ## set optim size for chunks in N clusters TODO: (make this better)
    optim = int(sample.stats.clusters_kept/
                float(len(ipyclient.ids)*2))
    optim = 100
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
        ## grab optim clusters and write to file
        done, chunk = clustdealer(pairdealer, optim)
        chunkhandle = os.path.join(data.dirs.clusts, 
                                "tmp_"+str(sample.name)+"."+str(num*optim))
        if chunk:
            chunkslist.append(chunkhandle)            
            with open(chunkhandle, 'wb') as outchunk:
                outchunk.write("//\n//\n".join(chunk)+"//\n//\n")
            num += 1
        LOGGER.debug("chunking len:%s, done:%s, num:%s", len(chunk), done, num)

    ## close clusters handle
    clusters.close()

    ## send chunks across engines, will delete tmps if failed
    try:
        submitted_args = []
        for chunkhandle in chunkslist:
            ## used to increment names across processors
            args = [data, sample, chunkhandle, optim]
            submitted_args.append(args)
            num += 1

        lbview = ipyclient.load_balanced_view()
        results = lbview.map_async(consensus, submitted_args)
        statsdicts = results.get()
        del lbview

    except Exception as inst:
        print(inst)
        for key in ipyclient.history:
            if ipyclient.metadata[key].error:
                LOGGER.error("step5 error: %s", 
                    ipyclient.metadata[key].error)
                raise SystemExit
            if ipyclient.metadata[key].stdout:
                LOGGER.error("step5 stdout:%s", 
                    ipyclient.metadata[key].stdout)
                raise SystemExit            
        #raise SystemExit

    finally:
        ## if process failed at any point delete tmp files
        for tmpchunk in chunkslist:
            os.remove(tmpchunk)

    return statsdicts

        


def run(data, samples, ipyclient, force=False):
    """ checks if the sample should be run and passes the args """

    ## message to skip all samples
    if not force:
        if all([i.stats.state >= 5 for i in samples]):
            print("Skipping step5: All {} ".format(len(data.samples))\
                 +"Samples already have consens reads ")

    ## prepare dirs
    data.dirs.consens = os.path.join(data.dirs.working, data.name+"_consens")
    if not os.path.exists(data.dirs.consens):
        os.mkdir(data.dirs.consens)

    ## zap any tmp files that might be leftover
    tmpcons = glob.glob(os.path.join(data.dirs.consens, "*_tmpcons.*"))
    tmpcats = glob.glob(os.path.join(data.dirs.consens, "*_tmpcats.*"))
    for tmpfile in tmpcons+tmpcats:
        os.remove(tmpfile)

    ## if sample is already done skip
    if data.stats.hetero_est.empty:
        print("  No estimates of heterozygosity and error rate. Using default "\
              "values")
        for sample in data.samples:
            sample.stats.hetero_est = 0.001
            sample.stats.error_est = 0.0001

    if data.paramsdict["ploidy"] == 1:
        print("  Haploid base calls and paralog filter (max haplos = 1)")
    elif data.paramsdict["ploidy"] == 2:
        print("  Diploid base calls and paralog filter (max haplos = 2)")
    elif data.paramsdict["ploidy"] == 2:
        print("  Diploid base calls and no paralog filter "\
                "(max haplos = {})".format(data.paramsdict["ploidy"]))
    print("  error rate (mean, std):  " \
             +"{:.5f}, ".format(data.stats.error_est.mean()) \
             +"{:.5f}\n".format(data.stats.error_est.std()) \
          +"  heterozyg. (mean, std):  " \
             +"{:.5f}, ".format(data.stats.hetero_est.mean()) \
             +"{:.5f}\n".format(data.stats.hetero_est.std()))

    ## Samples on queue
    for sample in samples:
        ## not force need checks
        if not force:
            if sample.stats.state >= 5:
                print("Skipping Sample {}; ".format(sample.name)
                     +"Already has consens reads. Use force=True to overwrite.")
            elif sample.stats.clusters_kept < 100:
                print("Skipping Sample {}; ".format(sample.name)
                     +"Too few clusters ({}). Use force=True to run anyway.".\
                       format(sample.stats.clusters_kept))
            else:
                statsdicts = run_full(data, sample, ipyclient)
                cleanup(data, sample, statsdicts)

        else:
            if not sample.stats.clusters_kept:
                print("Skipping Sample {}; ".format(sample.name)
                     +"No clusters found in file {}".\
                       format(sample.files.clusters_kept))
            else:
                statsdicts = run_full(data, sample, ipyclient)
                cleanup(data, sample, statsdicts)



if __name__ == "__main__":
    import ipyrad as ip
    TESTER = ip.load_dataobj("test")
    TESTER.step5.run()
