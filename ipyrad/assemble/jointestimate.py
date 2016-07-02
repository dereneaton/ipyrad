#!/usr/bin/env python2

""" jointly infers heterozygosity and error rate from stacked sequences """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=W0212

import scipy.stats
import scipy.optimize
import numpy as np
import ipyparallel
import itertools
import datetime
import time
import gzip
import os

from collections import Counter
from util import *


#from numba import jit
# try: 
#     
#     NUMBA = 1
# except ImportError:
#     NUMBA = 0



def frequencies(stacked):
    """ return frequency counts """
    totals = stacked.sum(axis=1)
    totals = totals.sum(axis=0)
    freqs = totals/np.float32(totals.sum())
    return freqs


def get_frequencies(stacked):
    """ return frequency counts """
    freqs = stacked.sum(axis=0) / np.float32(stacked.sum())
    return freqs


# @jit(['float32[:,:](float32, float32, int16[:,:])'], nopython=True)
# def jlikelihood1(errors, bfreqs, ustacks):
#     """Probability homozygous. All numpy and no loop so there was 
#     no numba improvement to speed when tested. """
#     ## make sure base_frequencies are in the right order
#     #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
#     totals = np.array([ustacks.sum(axis=1)]*4).T
#     prob = scipy.stats.binom.pmf(totals-ustacks, totals, errors)
#     return np.sum(bfreqs*prob, axis=1)


def likelihood1(errors, bfreqs, ustacks):
    """Probability homozygous. All numpy and no loop so there was 
    no numba improvement to speed when tested. """
    ## make sure base_frequencies are in the right order
    #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    totals = np.array([ustacks.sum(axis=1)]*4).T
    prob = scipy.stats.binom.pmf(totals-ustacks, totals, errors)
    return np.sum(bfreqs*prob, axis=1)


# @jit(['float32[:](float32, float32[:], int16[:,:])'])
# def jlikelihood2(errors, bfreqs, ustacks):
#     """probability of heterozygous. Very minimal speedup w/ numba."""
#     returns = np.zeros(len(ustacks), dtype=np.float32)
#     for idx, ustack in enumerate(ustacks):
#         spair = np.array(list(itertools.combinations(ustack, 2)))
#         bpair = np.array(list(itertools.combinations(bfreqs, 2)))
#         one = 2.*bpair.prod(axis=1)
#         tot = ustack.sum()
#         atwo = tot - spair[:, 0] - spair[:, 1]
#         two = scipy.stats.binom.pmf(atwo, tot, (2.*errors)/3.)
#         three = scipy.stats.binom.pmf(\
#                     spair[:, 0], spair.sum(axis=1), 0.5)
#         four = 1.-np.sum(bfreqs**2)
#         returns[idx] = np.sum(one*two*(three/four))
#     return np.array(returns)


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


# @jit
# def j_diploid_lik(pstart, bfreqs, ustacks, counts):
#     """ Log likelihood score given values [H,E]. """
#     hetero, errors = pstart
#     ## tell it to score terribly if scores are negative
#     if (hetero <= 0.) or (errors <= 0.):
#         score = np.exp(100)
#     else:
#         ## get likelihood for all sites
#         lik1 = (1.-hetero)*jlikelihood1(errors, bfreqs, ustacks)
#         lik2 = (hetero)*jlikelihood2(errors, bfreqs, ustacks)
#         liks = lik1+lik2
#         logliks = np.log(liks[liks > 0])*counts[liks > 0]
#         score = -logliks.sum()
#     return score


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



def tablestack(rstack):
    """ makes a count dict of each unique array element """
    ## goes by 10% at a time to minimize memory overhead. 
    ## Check on the last chunk.
    # table = Counter()
    # for i in xrange(0, rstack.shape[0], rstack.shape[0]//10):
    #     tmp = Counter([j.tostring() for j in rstack[i:i+rstack.shape[0]//10]])
    #     table.update(tmp)
    # return table
    table = Counter([j.tostring() for j in rstack])
    #LOGGER.info('table %s', table)
    return table



def stackarray(data, sample, sub):
    """ makes a list of lists of reads at each site """
    ## get clusters file
    clusters = gzip.open(sample.files.clusters)
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## array will be (nclusters, readlen, 4)
    if "pair" in data.paramsdict["datatype"]:
        readlen = 2*data._hackersonly["max_fragment_length"]
    else:
        readlen = data._hackersonly["max_fragment_length"]

    ## TODO: make clusters_hidepth dynamic in case params change
    if sub:
        dims = (int(sub), readlen, 4)
    else:
        dims = (int(sample.stats.clusters_hidepth), readlen, 4)
    stacked = np.zeros(dims, dtype=np.uint32)
    LOGGER.info("sample %s, dims %s", sample.name, stacked.shape)

    ## don't use sequence edges / restriction overhangs
    cutlens = [None, None, None, None]
    for cidx, cut in enumerate(data.paramsdict["restriction_overhang"]):
        if cut:
            cutlens[cidx] = len(cut)
    try:
        cutlens[1] = -1*cutlens[1]
    except TypeError:
        pass
    LOGGER.info(cutlens)

    ## fill stacked
    done = 0
    nclust = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)

        if chunk:
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            ## pull replicate read info from seqs
            reps = [int(sname.split(";")[-2][5:]) for sname in names]
            ## double reps if the read was fully merged.
            merged = ["_m1;s" in sname for sname in names]
            if any(merged):
                reps = [i*2 if j else i for i, j in zip(reps, merged)]

            sseqs = [list(seq) for seq in seqs]
            try:
                arrayed = np.concatenate(
                      [[seq]*rep for seq, rep in zip(sseqs, reps)])
            except ValueError:
                LOGGER.info("sseqs %s, reps %s", "\n".join(["".join(i) for i in sseqs]), reps)
            
            ## enforce minimum depth for estimates
            if arrayed.shape[0] >= data.paramsdict["mindepth_statistical"]:
                ## remove edge columns
                arrayed = arrayed[:, cutlens[0]:cutlens[1]]
                ## remove cols that are pair separator
                arrayed = arrayed[:, ~np.any(arrayed == "n", axis=0)]
                ## remove cols that are all Ns after converting -s to Ns
                arrayed[arrayed == "-"] = "N"
                arrayed = arrayed[:, ~np.all(arrayed == "N", axis=0)]
                ## store in stacked dict
                catg = np.array(\
                    [np.sum(arrayed == i, axis=0) for i in list("CATG")], 
                    dtype=np.uint32).T
                
                try:
                    stacked[nclust, :catg.shape[0], :] = catg
                    nclust += 1
                except IndexError:
                    LOGGER.error("nclust=%s, sname=%s, ", nclust, sample.name)

                ## finish early if subsampling
                if nclust == sub:
                    done = 1
    ## drop the empty rows in case there are fewer loci than the size of array
    #tots = sum(stacked.sum(axis=1) == 0)
    #newstack = np.zeros((tots, readlen, 4), dtype=np.uint32)
    newstack = stacked[stacked.sum(axis=2) > 0]
    LOGGER.info(newstack)
    assert np.all(newstack.sum(axis=1) > 0), "not all zeros!"

    return newstack



def optim(args):
    """ func scipy optimize to find best parameters"""

    ## split args
    data, sample, sub = args

    ## get array of all clusters data
    stacked = stackarray(data, sample, sub)
    LOGGER.info('stack %s', stacked)

    ## get base frequencies
    bfreqs = get_frequencies(stacked)
    #LOGGER.debug(bfreqs)
    if np.isnan(bfreqs).any():
        raise IPyradError("Caught sample with bad stack - {} {}".\
                          format(sample.name, bfreqs))

    ## reshape to concatenate all site rows
    #rstack = stacked.reshape(stacked.shape[0]*stacked.shape[1],
    #                         stacked.shape[2])
    rstack = stacked.astype(np.uint32)
    ## put into array, count array items as Byte strings
    tstack = tablestack(rstack)
    ## drop emtpy count [0,0,0,0] <- no longer needed, since we drop 0s earlier
    #tstack.pop(np.zeros(4, dtype=np.uint32).tostring())
    ## get keys back as arrays and store vals as separate arrays
    ustacks = np.array([np.fromstring(i, dtype=np.uint32) \
                        for i in tstack.iterkeys()])
    counts = np.array(tstack.values())
    ## cleanup    
    del rstack
    del tstack
    LOGGER.info(bfreqs)
    LOGGER.info(ustacks)
    LOGGER.info(counts)    

    ## if data are haploid fix H to 0
    if data.paramsdict["max_alleles_consens"] == 1:
        pstart = np.array([0.001], dtype=np.float32)
        hetero = 0.
        errors = scipy.optimize.fmin(get_haploid_lik, pstart,
                                    (bfreqs, ustacks, counts),
                                     disp=False,
                                     full_output=False)
    ## or do joint diploid estimates
    else:
        pstart = np.array([0.01, 0.001], dtype=np.float32)
        hetero, errors = scipy.optimize.fmin(get_diploid_lik, pstart,
                                            (bfreqs, ustacks, counts), 
                                             disp=False,
                                             full_output=False)
    return hetero, errors



def run(data, samples, subsample, force, ipyclient):
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
                submitted_args.append([data, sample, subsample])
        else:
            if sample.stats.state < 3:
                print("    "+sample.name+" not clustered. Run step3() first.")
            elif sample.stats.clusters_hidepth < 2:
                print("    skipping {}. Too few high depth reads ({}). "\
                      .format(sample.name, sample.stats.clusters_hidepth))
            else:
                submitted_args.append([data, sample, subsample])

    if submitted_args:    
        ## submit jobs to parallel client
        submit(data, submitted_args, ipyclient)
    #else:
    #    print("  no samples selected")



def submit(data, submitted_args, ipyclient):
    """ 
    sends jobs to engines and cleans up failures. Print progress. 
    """

    ## first sort by cluster size
    submitted_args.sort(key=lambda x: x[1].stats.clusters_hidepth, reverse=True)
                                           
    ## send all jobs to a load balanced client
    lbview = ipyclient.load_balanced_view()
    jobs = {}
    for args in submitted_args:
        #LOGGER.info("submitting duh %s", args)
        ## stores async results using sample names
        jobs[args[1].name] = lbview.apply(optim, args)

    ## dict to store cleanup jobs
    start = time.time()
    fwait = 0
    error = 0

    ## each job is submitted to cleanup as it finishes
    allwait = len(jobs)
    try:
        while 1:
            elapsed = datetime.timedelta(seconds=int(time.time() - start))
            progressbar(allwait, fwait,
                    " inferring [H, E]      | {}".format(elapsed))
 
            ## get job keys
            keys = jobs.keys()
            for sname in keys:
                if jobs[sname].ready():
                    if jobs[sname].successful():
                        ## get the results
                        hest, eest = jobs[sname].get()
                        sample_cleanup(data.samples[sname], hest, eest)
                        del jobs[sname]
                        fwait += 1

                    else:
                        ## grab metadata                
                        meta = jobs[sname].metadata
                        ## if not done do nothing, if failure print error
                        LOGGER.error('  sample %s did not finish', sname)
                        if meta.error:
                            LOGGER.error("""\
                        stdout: %s
                        stderr: %s 
                        error: %s""", meta.stdout, meta.stderr, meta.error)
                        del jobs[sname]
                        fwait += 1

            ## wait
            if fwait == allwait:
                break
            else:
                time.sleep(0.5)

        ## print final statement
        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        progressbar(20, 20, 
                " inferring [H, E]      | {}".format(elapsed))
        print("")


    except KeyboardInterrupt as kbt:
        ## reraise kbt after cleanup
        error = 1

    finally:
        assembly_cleanup(data)
        if error:
            raise kbt



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

