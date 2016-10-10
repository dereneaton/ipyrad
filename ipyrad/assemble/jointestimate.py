#!/usr/bin/env python2

""" jointly infers heterozygosity and error rate from stacked sequences """

from __future__ import print_function

import scipy.stats
import scipy.optimize
import numpy as np
import itertools
import datetime
import time
import gzip
import os

from collections import Counter
from util import *


# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=W0142



def likelihood1(errors, bfreqs, ustacks):
    """
    Probability homozygous. """
    ## make sure base_frequencies are in the right order
    #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    #totals = np.array([ustacks.sum(axis=1)]*4).T
    totals = np.array([ustacks.sum(axis=1)]*4).T
    prob = scipy.stats.binom.pmf(totals-ustacks, totals, errors)
    lik1 = np.sum(bfreqs*prob, axis=1)
    return lik1



def likelihood2(errors, bfreqs, ustacks):
    """probability of heterozygous"""
    returns = np.zeros([len(ustacks)])
    for idx, ustack in enumerate(ustacks):
        spair = np.array(list(itertools.combinations(ustack, 2)))
        bpair = np.array(list(itertools.combinations(bfreqs, 2)))
        one = 2.*bpair.prod(axis=1)
        tot = ustack.sum()
        atwo = tot - spair[:, 0] - spair[:, 1]
        two = scipy.stats.binom.pmf(atwo, tot, (2.*errors)/3.)
        three = scipy.stats.binom.pmf(\
                    spair[:, 0], spair.sum(axis=1), 0.5)
        four = 1.-np.sum(bfreqs**2)
        returns[idx] = np.sum(one*two*(three/four))
    return np.array(returns)


## more verbose and slow form of the function above
def liketest2(errors, bfreqs, ustack):
    """probability of heterozygous"""

    fullsum = 0
    for idx in xrange(4):
        subsum = 0
        for jdx in xrange(4):
            one = 2. * bfreqs[idx] * bfreqs[jdx]
            tot = ustack.sum()
            two = scipy.stats.binom.pmf(tot - ustack[idx] - ustack[jdx], 
                                        tot, (2.*errors)/3.)
            three = scipy.stats.binom.pmf(ustack[idx], 
                                          ustack[idx] + ustack[jdx], 0.5)
            four = 1 - np.sum(bfreqs**2)
            subsum += one * two * (three / four)
        fullsum += subsum
    return fullsum




def get_diploid_lik(pstart, bfreqs, ustacks, counts):
    """ Log likelihood score given values [H,E] """
    hetero, errors = pstart
    if (hetero <= 0.) or (errors <= 0.):
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = (1.-hetero)*likelihood1(errors, bfreqs, ustacks)
        lik2 = (hetero)*likelihood2(errors, bfreqs, ustacks)
        liks = lik1+lik2
        logliks = np.log(liks[liks > 0])*counts[liks > 0]
        score = -logliks.sum()
    return score



def get_haploid_lik(errors, bfreqs, ustacks, counts):
    """ Log likelihood score given values [E]. This can be written to run much
    faster by executing across the whole array, and/or by also in parallel """
    hetero = 0.
    ## score terribly if below 0
    if errors <= 0.:
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = ((1.-hetero)*likelihood1(errors, bfreqs, ustacks)) 
        lik2 = (hetero)*likelihood2(errors, bfreqs, ustacks)
        liks = lik1+lik2
        logliks = np.log(liks[liks > 0])*counts[liks > 0]
        score = -logliks.sum()
    return score



def stackarray(data, sample, subloci):
    """ 
    Stacks clusters into arrays
    """

    ## get clusters file    
    #LOGGER.info("Entering stackarray - {}".format(sample.name))
    clusters = gzip.open(sample.files.clusters)
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## array will be (nclusters, readlen, 4)
    maxlen = data._hackersonly["max_fragment_length"]

    ## we subsample, else use first 10000 loci.
    dims = (data.stats.clusters_hidepth.max(), maxlen, 4)
    stacked = np.zeros(dims, dtype=np.uint64)

    ## don't use sequence edges / restriction overhangs
    cutlens = [None, None]
    try:
        cutlens[0] = len(data.paramsdict["restriction_overhang"][0])
        cutlens[1] = maxlen - len(data.paramsdict["restriction_overhang"][1])
    except TypeError:
        pass
    #LOGGER.info("cutlens: %s", cutlens)

    ## fill stacked
    done = 0
    nclust = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("  clustfile formatting error in %s", chunk)

        if chunk:
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            ## pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]
            
            ## double reps if the read was fully merged... (TODO: Test this!)
            merged = ["_m1;s" in sname for sname in names]
            if any(merged):
                reps = [i*2 if j else i for i, j in zip(reps, merged)]

            ## get all reps
            sseqs = [list(seq) for seq in seqs]
            arrayed = np.concatenate(
                          [[seq]*rep for seq, rep in zip(sseqs, reps)])
            
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
                catg = np.array(\
                    [np.sum(arrayed == i, axis=0) for i in list("CATG")], 
                    dtype=np.uint64).T

                stacked[nclust, :catg.shape[0], :] = catg
                nclust += 1

        ## finish early if subsampling
        if nclust == subloci:
            done = 1

    ## drop the empty rows in case there are fewer loci than the size of array
    newstack = stacked[stacked.sum(axis=2) > 0]
    assert not np.any(newstack.sum(axis=1) == 0), "no zero rows"

    return newstack



def optim(data, sample, subloci):
    """ func scipy optimize to find best parameters"""

    ## get array of all clusters data
    stacked = stackarray(data, sample, subloci)
    #LOGGER.info('stack %s', stacked)

    ## get base frequencies
    bfreqs = stacked.sum(axis=0) / float(stacked.sum())
    #bfreqs = bfreqs**2
    
    #LOGGER.debug(bfreqs)
    if np.isnan(bfreqs).any():
        #LOGGER.error("  Caught sample with bad stack: %s %s", 
        #                sample.name, bfreqs)
        raise IPyradWarningExit(" Bad stack in getfreqs; {} {}"\
               .format(sample.name, bfreqs))

    ## put into array, count array items as Byte strings
    tstack = Counter([j.tostring() for j in stacked])
    ## get keys back as arrays and store vals as separate arrays
    ustacks = np.array([np.fromstring(i, dtype=np.uint64) \
                        for i in tstack.iterkeys()])
    ## make bi-allelic only
    #tris = np.where(np.sum(ustacks > 0, axis=1) > 2)
    #for tri in tris:
    #    minv = np.min(ustacks[tri][ustacks[tri] > 0])
    #    delv = np.where(ustacks[tri] == minv)[0][0]
    #    ustacks[tri, delv] = 0

    counts = np.array(tstack.values())
    ## cleanup    
    del tstack


    ## if data are haploid fix H to 0
    if int(data.paramsdict["max_alleles_consens"]) == 1:
        pstart = np.array([0.001], dtype=np.float64)
        hetero = 0.
        errors = scipy.optimize.fmin(get_haploid_lik, pstart,
                                    (bfreqs, ustacks, counts),
                                     disp=False,
                                     full_output=False)
    ## or do joint diploid estimates
    else:
        pstart = np.array([0.01, 0.001], dtype=np.float64)
        hetero, errors = scipy.optimize.fmin(get_diploid_lik, pstart,
                                            (bfreqs, ustacks, counts), 
                                             disp=False,
                                             full_output=False)
    return hetero, errors



def run(data, samples, subloci, force, ipyclient):
    """ calls the main functions """

    ## speed hack == use only the first 2000 high depth clusters for estimation.
    ## based on testing this appears sufficient for accurate estimates
    ## the default is set in assembly.py

    # if haploid data
    if data.paramsdict["max_alleles_consens"] == 1:
        print("  Applying haploid-based test (infer E with H fixed to 0).")

    submitted_args = []
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
                submitted_args.append([sample, subloci])
        else:
            if sample.stats.state < 3:
                print("    "+sample.name+" not clustered. Run step3() first.")
            elif sample.stats.clusters_hidepth < 2:
                print("    skipping {}. Too few high depth reads ({}). "\
                      .format(sample.name, sample.stats.clusters_hidepth))
            else:
                submitted_args.append([sample, subloci])

    if submitted_args:    
        ## submit jobs to parallel client
        submit(data, submitted_args, ipyclient)



def submit(data, submitted_args, ipyclient):
    """ 
    Sends jobs to engines and cleans up failures. Print progress. 
    """

    ## first sort by cluster size
    submitted_args.sort(key=lambda x: x[0].stats.clusters_hidepth, reverse=True)
                                           
    ## send all jobs to a load balanced client
    lbview = ipyclient.load_balanced_view()
    jobs = {}
    for sample, subloci in submitted_args:
        ## stores async results using sample names
        jobs[sample.name] = lbview.apply(optim, *(data, sample, subloci))

    ## dict to store cleanup jobs
    start = time.time()

    ## wrap in a try statement so that stats are saved for finished samples.
    ## each job is submitted to cleanup as it finishes
    try:
        ## wait for jobs to finish
        while 1:
            fin = [i.ready() for i in jobs.values()]
            elapsed = datetime.timedelta(seconds=int(time.time() - start))
            progressbar(len(fin), sum(fin), 
                " inferring [H, E]      | {} | s4 |".format(elapsed))
            time.sleep(0.1)
            if len(fin) == sum(fin):
                print("")
                break

        ## cleanup
        for job in jobs:
            if jobs[job].successful():
                hest, eest = jobs[job].result()
                sample_cleanup(data.samples[job], hest, eest)
            else:
                LOGGER.error("  Sample %s failed with error %s", 
                             job, jobs[job].exception())
                raise IPyradWarningExit("  Sample {} failed step 4"\
                                        .format(job))

    finally:
        assembly_cleanup(data)



def sample_cleanup(sample, hest, eest):
    """ 
    Stores results to the Assembly object, writes to stats file, 
    and cleans up temp files 
    """
    ## sample summary assignments
    sample.stats.state = 4
    sample.stats.hetero_est = hest
    sample.stats.error_est = eest

    ## sample full assigments
    sample.stats_dfs.s4.hetero_est = hest
    sample.stats_dfs.s4.error_est = eest



def assembly_cleanup(data):
    """ cleanup assembly stats """
    ## Assembly assignment
    data.stats_dfs.s4 = data.build_stat("s4")#, dtype=np.float32)

    ## Update written file
    data.stats_files.s4 = os.path.join(data.dirs.clusts, 
                                       "s4_joint_estimate.txt")
    with open(data.stats_files.s4, 'w') as outfile:
        data.stats_dfs.s4.to_string(outfile)




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

