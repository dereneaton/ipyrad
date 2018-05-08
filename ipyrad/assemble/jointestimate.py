#!/usr/bin/env python

""" jointly infers heterozygosity and error rate from stacked sequences """

# py2/3 compatible
from __future__ import print_function
try:
    from itertools import izip, combinations
except ImportError:
    from itertools import combinations
    izip = zip

import scipy.optimize
import scipy.stats
import ipyrad as ip
import numpy as np
import numba
import time
import gzip
import os

from .cluster_within import get_quick_depths
from collections import Counter
from .util import IPyradWarningExit, IPyradError, clustdealer


def likelihood1(errors, bfreqs, ustacks):
    """
    Probability homozygous. """
    ## make sure base_frequencies are in the right order
    #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    #totals = np.array([ustacks.sum(axis=1)]*4).T
    totals = np.array([ustacks.sum(axis=1)] * 4).T
    prob = scipy.stats.binom.pmf(totals - ustacks, totals, errors)
    lik1 = np.sum(bfreqs * prob, axis=1)
    return lik1



@numba.jit(nopython=True)
def nblik2_build(ustacks):
    """
    JIT'd function builds large array that can be used to calc binom pmf
    """

    ## store
    #ret = np.empty(ustacks.shape[0])

    ## fill for pmf later 
    tots = np.empty((ustacks.shape[0], 1))
    twos = np.empty((ustacks.shape[0], 6))
    thrs = np.empty((ustacks.shape[0], 6, 2))

    ## fill big arrays
    for idx in range(ustacks.shape[0]):

        ust = ustacks[idx]
        tot = ust.sum()
        tots[idx] = tot

        ## fast filling of arrays
        i = 0
        for jdx in range(4):
            for kdx in range(4):
                if jdx < kdx:
                    twos[idx][i] = tot - ust[jdx] - ust[kdx]
                    thrs[idx][i] = ust[jdx], ust[jdx] + ust[kdx]
                    i += 1

    return tots, twos, thrs



def lik2_calc(err, one, tots, twos, thrs, four):
    """ 
    vectorized calc of binom pmf on large arrays 
    """

    ## calculate twos
    _twos = scipy.stats.binom.pmf(twos, tots, 0.5)

    ## calculate threes
    _thrs = thrs.reshape(thrs.shape[0] * thrs.shape[1], thrs.shape[2])
    _thrs = scipy.stats.binom.pmf(_thrs[:, 0], _thrs[:, 1], (2. * err) / 3.)
    _thrs = _thrs.reshape(thrs.shape[0], 6)

    ## calculate return sums
    return np.sum(one * _twos * (_thrs / four), axis=1)



def nlikelihood2(errors, bfreqs, ustacks):
    """ calls nblik2_build and lik2_calc for a given err """
    one = [2. * bfreqs[i] * bfreqs[j] for i, j in combinations(range(4), 2)]
    four = 1. - np.sum(bfreqs**2) 
    tots, twos, thrs = nblik2_build(ustacks)
    res2 = lik2_calc(errors, one, tots, twos, thrs, four)
    return res2



def nget_diploid_lik(pstart, bfreqs, ustacks, counts):
    """ Log likelihood score given values [H,E] """
    hetero, errors = pstart
    if (hetero <= 0.) or (errors <= 0.):
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = (1. - hetero) * likelihood1(errors, bfreqs, ustacks)
        lik2 = (hetero) * nlikelihood2(errors, bfreqs, ustacks)
        liks = lik1 + lik2
        logliks = np.log(liks[liks > 0]) * counts[liks > 0]
        score = -logliks.sum()
    return score



def get_haploid_lik(errors, bfreqs, ustacks, counts):
    """ Log likelihood score given values [E]. """
    hetero = 0.
    ## score terribly if below 0
    if errors <= 0.:
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = ((1. - hetero) * likelihood1(errors, bfreqs, ustacks)) 
        #lik2 = (hetero)*nlikelihood2(errors, bfreqs, ustacks)
        #liks = lik1+lik2
        liks = lik1
        logliks = np.log(liks[liks > 0]) * counts[liks > 0]
        score = -logliks.sum()
    return score



def recal_hidepth(data, sample):
    """
    if mindepth setting were changed then 'clusters_hidepth' needs to be 
    recalculated. Check and recalculate if necessary.
    """
    ## the minnest depth
    majrdepth = data.paramsdict["mindepth_majrule"]
    statdepth = data.paramsdict["mindepth_statistical"]    

    ## if nothing changes return existing maxlen value
    maxlen = data._hackersonly["max_fragment_length"]

    ## if coming from older version of ipyrad this attr is new
    if not hasattr(sample.stats_dfs.s3, "hidepth_min"):
        sample.stats_dfs.s3["hidepth_min"] = data.paramsdict["mindepth_majrule"]

    ## if old value not the same as current value then recalc
    if 1:  # not sample.stats_dfs.s3["hidepth_min"] == majrdepth:
        ip.logger.info(" mindepth setting changed: recalculating clusters_hidepth and maxlen")
        ## get arrays of data
        maxlens, depths = get_quick_depths(data, sample)

        ## calculate how many are hidepth
        hidepths = depths >= majrdepth
        stathidepths = depths >= statdepth

        keepmj = depths[hidepths]
        keepst = depths[stathidepths]

        try:
            ## set assembly maxlen for sample
            statlens = maxlens[stathidepths]
            statlen = int(statlens.mean() + (2. * statlens.std()))
        except:
            ## If no clusters have depth sufficient for statistical basecalling
            ## then stathidepths will be empty and all hell breaks loose, so
            ## we'll raise here and than catch the exception in optim()
            raise IPyradError("No clusts with depth sufficient for statistical basecalling.")

        ip.logger.info("%s %s %s", maxlens.shape, maxlens.mean(), maxlens.std())
        maxlens = maxlens[hidepths]
        maxlen = int(maxlens.mean() + (2.*maxlens.std()))

        ## saved stat values are for majrule min
        sample.stats["clusters_hidepth"] = keepmj.shape[0]
        sample.stats_dfs.s3["clusters_hidepth"] = keepmj.shape[0]        

    return sample, keepmj.shape[0], maxlen, keepst.shape[0], statlen


def stackarray(data, sample):
    """ 
    Stacks clusters into arrays
    """

    ## only use clusters with depth > mindepth_statistical for param estimates
    sample, _, _, nhidepth, maxlen = recal_hidepth(data, sample)

    ## get clusters file    
    clusters = gzip.open(sample.files.clusters, 'rb')
    pairdealer = izip(*[iter(clusters)] * 2)

    ## we subsample, else use first 10000 loci.
    dims = (nhidepth, maxlen, 4)
    stacked = np.zeros(dims, dtype=np.uint64)

    ## don't use sequence edges / restriction overhangs
    cutlens = [None, None]
    try:
        cutlens[0] = len(data.paramsdict["restriction_overhang"][0])
        cutlens[1] = maxlen - len(data.paramsdict["restriction_overhang"][1])
    except TypeError:
        pass
    #ip.logger.info("cutlens: %s", cutlens)

    ## fill stacked
    nclust = 0
    done = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("  clustfile formatting error in %s", chunk)

        if chunk:
            piece = chunk[0].decode().strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            ## pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]
            
            ## double reps if the read was fully merged... (TODO: Test this!)
            #merged = ["_m1;s" in sname for sname in names]
            #if any(merged):
            #    reps = [i*2 if j else i for i, j in zip(reps, merged)]

            ## get all reps
            sseqs = [list(seq) for seq in seqs]
            arrayed = np.concatenate([
                [seq] * rep for seq, rep in zip(sseqs, reps)
                ])
            
            ## enforce minimum depth for estimates
            if arrayed.shape[0] >= data.paramsdict["mindepth_statistical"]:
                ## remove edge columns and select only the first 500 
                ## derep reads, just like in step 5
                arrayed = arrayed[:500, cutlens[0]:cutlens[1]]
                ## remove cols that are pair separator
                arrayed = arrayed[:, ~np.any(arrayed == "n", axis=0)]
                ## remove cols that are all Ns after converting -s to Ns
                arrayed[arrayed == "-"] = "N"
                arrayed = arrayed[:, ~np.all(arrayed == "N", axis=0)]
                ## store in stacked dict

                catg = np.array(
                    [np.sum(arrayed == i, axis=0) for i in list("CATG")], 
                    dtype=np.uint64).T

                stacked[nclust, :catg.shape[0], :] = catg
                nclust += 1

    ## drop the empty rows in case there are fewer loci than the size of array
    newstack = stacked[stacked.sum(axis=2) > 0]
    assert not np.any(newstack.sum(axis=1) == 0), "no zero rows"
    clusters.close()

    return newstack



def optim(data, sample):
    """ func scipy optimize to find best parameters"""

    hetero = 0.01
    errors = 0.001

    ## A flag for communicating with sample_cleanup so a nice warning
    ## message can be displayed.
    success = False

    try:
        ## get array of all clusters data
        stacked = stackarray(data, sample)

        ## get base frequencies
        bfreqs = stacked.sum(axis=0) / float(stacked.sum())
        #bfreqs = bfreqs**2
        #ip.logger.debug(bfreqs)
        if np.isnan(bfreqs).any():
            raise IPyradWarningExit(" Bad stack in getfreqs; {} {}"
                   .format(sample.name, bfreqs))

        ## put into array, count array items as Byte strings
        tstack = Counter([j.tostring() for j in stacked])

        ## get keys back as arrays and store vals as separate arrays
        ustacks = np.array([np.frombuffer(i, dtype=np.uint64)
                            for i in tstack.keys()])

        ## make bi-allelic only
        #tris = np.where(np.sum(ustacks > 0, axis=1) > 2)
        #for tri in tris:
        #    minv = np.min(ustacks[tri][ustacks[tri] > 0])
        #    delv = np.where(ustacks[tri] == minv)[0][0]
        #    ustacks[tri, delv] = 0

        counts = np.array(list(tstack.values()))
        ## cleanup
        del tstack

        ## if data are haploid fix H to 0
        if int(data.paramsdict["max_alleles_consens"]) == 1:
            pstart = np.array([0.001], dtype=np.float64)
            hetero = 0.
            errors = scipy.optimize.fmin(
                get_haploid_lik, pstart,
                (bfreqs, ustacks, counts),
                maxfun=50,
                maxiter=50,
                disp=False,
                full_output=False)
        ## or do joint diploid estimates
        else:
            pstart = np.array([0.01, 0.001], dtype=np.float64)
            hetero, errors = scipy.optimize.fmin(
                nget_diploid_lik, pstart,
                (bfreqs, ustacks, counts),
                maxfun=50,
                maxiter=50,
                disp=False,
                full_output=False)
        success = True

    except IPyradError as inst:
        ## recal_hidepth raises this exception if there are no clusters that
        ## have depth sufficient for statistical basecalling. In this case
        ## we just set the default hetero and errors values to 0.01/0.001
        ## which is the default.
        ip.logger.debug("Found sample with no clusters hidepth - {}".format(sample.name))
        pass

    return hetero, errors, success



def run(data, samples, force, ipyclient):
    """ calls the main functions """

    # if haploid data
    if data.paramsdict["max_alleles_consens"] == 1:
        print("{}Applying haploid-based test (infer E with H fixed to 0)"\
              .format(data._spacer))

    subsamples = []

    ## if sample is already done skip
    for sample in samples:
        if not force:
            if sample.stats.state >= 4:
                print("    skipping {}; ".format(sample.name)+\
                      "already estimated. Use force=True to overwrite.")
            elif sample.stats.state < 3:
                print("    skipping {}; ".format(sample.name)+\
                      "not clustered yet. Run step3() first.")
            else:
                subsamples.append(sample)
        else:
            if sample.stats.state < 3:
                print("    "+sample.name+" not clustered. Run step3() first.")
            elif sample.stats.clusters_hidepth < 2:
                print("    skipping {}. Too few high depth reads ({}). "\
                      .format(sample.name, sample.stats.clusters_hidepth))
            else:
                subsamples.append(sample)

    if subsamples:
        ## submit jobs to parallel client
        submit(data, subsamples, ipyclient)



def submit(data, subsamples, ipyclient):
    """ 
    Sends jobs to engines and cleans up failures. Print progress. 
    """

    ## first sort by cluster size
    subsamples.sort(key=lambda x: x.stats.clusters_hidepth, reverse=True)
                            
    ## send all jobs to a load balanced client
    lbview = ipyclient.load_balanced_view()
    jobs = {}

    ## stores async results using sample names    
    for sample in subsamples:
        jobs[sample.name] = lbview.apply(optim, *(data, sample))

    ## wrap in a try statement so that stats are saved for finished samples.
    ## each job is submitted to cleanup as it finishes
    start = time.time() 
    printstr = ("inferring [H, E]    ", "s4")
    try:
        kbd = 0
        ## wait for jobs to finish
        while 1:
            fin = [i.ready() for i in jobs.values()]
            data._progressbar(len(fin), sum(fin), start, printstr)
            time.sleep(0.1)
            if len(fin) == sum(fin):
                break

        ## cleanup
        print("")
        for job in jobs:
            if jobs[job].successful():
                hest, eest, success = jobs[job].result()
                sample_cleanup(data.samples[job], hest, eest, success)
            else:
                ip.logger.error(
                    "  Sample %s failed with error %s", 
                    job, jobs[job].exception())
                raise IPyradWarningExit(jobs[job].result())

    except KeyboardInterrupt as kbd:
        pass

    finally:
        assembly_cleanup(data)
        if kbd:
            raise KeyboardInterrupt



def sample_cleanup(sample, hest, eest, success):
    """ 
    Stores results to the Assembly object, writes to stats file, 
    and cleans up temp files 
    """
    ## sample summary assignments
    sample.stats.state = 4
    sample.stats.hetero_est = float(hest)
    sample.stats.error_est = float(eest)

    ## sample full assigments
    sample.stats_dfs.s4.hetero_est = float(hest)
    sample.stats_dfs.s4.error_est = float(eest)

    ## In rare cases no hidepth clusters for statistical basecalling
    ## so we warn the user, but carry on with default values
    if not success:
        msg = """    Info: Sample {} - No clusters have sufficient depth for statistical
          basecalling. Setting default heterozygosity/error to 0.01/0.001.""".format(sample.name)
        print(msg)



def assembly_cleanup(data):
    """ cleanup assembly stats """
    ## Assembly assignment
    data.stats_dfs.s4 = data._build_stat("s4")  # dtype=np.float32)

    ## Update written file
    data.stats_files.s4 = os.path.join(data.dirs.clusts, 
                                       "s4_joint_estimate.txt")
    with open(data.stats_files.s4, 'w') as outfile:
        data.stats_dfs.s4.to_string(outfile)

    fails = data.stats[data.stats["state"] == 3].index.values
    if fails:
        msg = """
        These samples failed joint estimation and will be excluded from
        downstream analysis (probably very few highdepth reads):
        {}""".format(fails)
        print(msg)



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

    ## run test on RAD data1
    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "test_pairgbs", "test_pairgbs"))
    TEST.step4(force=True)
    print(TEST.stats)

    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "test_rad", "data1"))
    TEST.step4(force=True)
    print(TEST.stats)

    ## run test on messy data set
    #TEST = ip.load_assembly(os.path.join(ROOT, "tests", "radmess", "data1"))

    ## check if results are correct

    ## cleanup

