#!/usr/bin/env python

""" jointly infers heterozygosity and error rate from stacked sequences """

# py2/3 compatible
from __future__ import print_function
try:
    from itertools import izip, combinations
except ImportError:
    from itertools import combinations
    izip = zip

import os
import time
import gzip
from collections import Counter

import scipy.optimize
import scipy.stats
import numpy as np
import numba

from .clustmap import get_quick_depths
from .utils import IPyradError, clustdealer


class Step4:
    "organized functions for step 4 of assembly"
    def __init__(self, data, force, ipyclient):

        self.data = data
        self.force = force
        self.ipyclient = ipyclient
        self.haploid = bool(data.params.max_alleles_consens == 1)
        if self.haploid:
            print("Running haploid inference (infer E with H fixed to 0)")
        self.print_headers()
        self.samples = self.get_subsamples()


    def print_headers(self):
        if self.data._cli:
            self.data._print(
                "\n{}Step 4: Joint estimation of error rate and heterozygosity"
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        "Apply state, ncluster, and force filters to select samples"

        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 4")

        # filter samples by state
        state2 = self.data.stats.index[self.data.stats.state < 3]
        state3 = self.data.stats.index[self.data.stats.state == 3]
        state4 = self.data.stats.index[self.data.stats.state > 3]

        # tell user which samples are not ready for step4
        if state2.any():
            print("skipping samples not in state==3:\n{}"
                  .format(state2.tolist()))

        if self.force:
            # run all samples above state 2
            subs = self.data.stats.index[self.data.stats.state > 2]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already cmopleted step 4
            if state4.any():
                print("skipping samples already finished step 4:\n{}"
                      .format(state4.tolist()))

            # run all samples in state 4
            subsamples = [self.data.samples[i] for i in state3]

        # check that kept samples have clusters
        checked_samples = []
        nodata_samples = []
        for sample in subsamples:
            if sample.stats.clusters_hidepth:
                checked_samples.append(sample)
            else:
                nodata_samples.append(sample)
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 4")
        else:
            for sample in nodata_samples:
                print("{}skipping {}; no clusters found."
                    .format(self.data._spacer, sample))

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.clusters_hidepth,
            reverse=True,
        )
        return checked_samples


    def run(self):
        "call the remote functions"
        self.remote_run_optim()
        self.cleanup()


    def remote_run_optim(self):
        "call the optim function to run in parallel"
        # headers
        start = time.time()
        printstr = ("inferring [H, E]    ", "s4")

        # send all jobs to a load balanced client
        lbview = self.ipyclient.load_balanced_view()

        # stores async results using sample names
        jobs = {}
        for sample in self.samples:
            jobs[sample.name] = lbview.apply(optim, *(self.data, sample))

        # progress bar
        while 1:
            fin = [i.ready() for i in jobs.values()]
            self.data._progressbar(len(fin), sum(fin), start, printstr)
            time.sleep(0.1)
            if len(fin) == sum(fin):
                break

        # cleanup
        self.data._print("")        
        for job in jobs:
            # collect results
            hest, eest, success = jobs[job].get()
            # store results to sample objects
            sample_cleanup(self.data.samples[job], hest, eest, success)


    def cleanup(self):
        " assembly cleanup "

        # store results
        self.data.stats_dfs.s4 = self.data._build_stat("s4")

        # make path for stats file
        self.data.stats_files.s4 = os.path.join(
            self.data.dirs.clusts,
            "s4_joint_estimate.txt")

        # write to stats file
        with open(self.data.stats_files.s4, 'w') as outfile:
            self.data.stats_dfs.s4.to_string(outfile)

        # check for failures
        if any(self.data.stats.state == 3):
            msg = ("""
            These samples failed joint estimation and will be excluded from
            downstream analysis (probably very few highdepth reads):
            {}""".format(
                self.data.stats[self.data.stats.state == 3].index.tolist()))
            print(msg)

        # clean out old stats files
        if self.force:
            self.data.stats_files.s5 = ""
            self.data.stats_files.s6 = ""
            self.data.stats_files.s7 = ""



############################################################################


def likelihood1(errors, bfreqs, ustacks):
    """
    Probability homozygous. """
    ## make sure base_frequencies are in the right order
    # print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    # totals = np.array([ustacks.sum(axis=1)]*4).T
    totals = np.array([ustacks.sum(axis=1)] * 4).T
    prob = scipy.stats.binom.pmf(totals - ustacks, totals, errors)
    lik1 = np.sum(bfreqs * prob, axis=1)
    return lik1



@numba.jit(nopython=True)
def nblik2_build(ustacks):
    "JIT'd function builds array that can be used to calc binom pmf"
    # fill for pmf later
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
    "vectorized calc of binom pmf on large arrays"
    ## calculate twos
    _twos = scipy.stats.binom.pmf(twos, tots, 0.5)

    ## calculate threes
    _thrs = thrs.reshape(thrs.shape[0] * thrs.shape[1], thrs.shape[2])
    _thrs = scipy.stats.binom.pmf(_thrs[:, 0], _thrs[:, 1], (2. * err) / 3.)
    _thrs = _thrs.reshape(thrs.shape[0], 6)

    ## calculate return sums
    return np.sum(one * _twos * (_thrs / four), axis=1)



def nlikelihood2(errors, bfreqs, ustacks):
    "calls nblik2_build and lik2_calc for a given err"
    one = [2. * bfreqs[i] * bfreqs[j] for i, j in combinations(range(4), 2)]
    four = 1. - np.sum(bfreqs**2)
    tots, twos, thrs = nblik2_build(ustacks)
    res2 = lik2_calc(errors, one, tots, twos, thrs, four)
    return res2



def nget_diploid_lik(pstart, bfreqs, ustacks, counts):
    "Log likelihood score given values [H,E]"
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
    "Log likelihood score given values [E]"
    hetero = 0.
    ## score terribly if below 0
    if errors <= 0.:
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = ((1. - hetero) * likelihood1(errors, bfreqs, ustacks))
        liks = lik1
        logliks = np.log(liks[liks > 0]) * counts[liks > 0]
        score = -logliks.sum()
    return score



def recal_hidepth(data, sample):
    """
    if mindepth setting were changed then 'clusters_hidepth' needs to be
    recalculated. Check and recalculate if necessary.
    """
    # the minnest depth
    majrdepth = data.params.mindepth_majrule
    statdepth = data.params.mindepth_statistical

    # if nothing changes return existing maxlen value
    maxlen = data.hackersonly.max_fragment_length

    # get arrays of data
    maxlens, depths = get_quick_depths(data, sample)

    # calculate how many are hidepth
    hidepths = depths >= majrdepth
    stathidepths = depths >= statdepth
    keepmj = depths[hidepths]
    keepst = depths[stathidepths]

    try:
        # set assembly maxlen for sample
        statlens = maxlens[stathidepths]
        statlen = int(statlens.mean() + (2. * statlens.std()))

    except:
        # If no clusters have depth sufficient for statistical basecalling
        # then stathidepths will be empty and all hell breaks loose, so
        # we'll raise here and than catch the exception in optim()
        raise IPyradError(
            "No clusts with depth sufficient for statistical basecalling. "
          + "I recommend you branch to drop this sample: {}"
            .format(sample.name)
            )

    maxlens = maxlens[hidepths]
    maxlen = int(maxlens.mean() + (2. * maxlens.std()))
    return keepmj.shape[0], maxlen, keepst.shape[0], statlen



def stackarray(data, sample):
    "Stacks clusters into arrays"
    # only use clusters with depth > mindepth_statistical for param estimates
    hidepth, maxlen, shidepth, smaxlen = recal_hidepth(data, sample)

    # (not saved) stat values are for majrule min
    sample.stats["clusters_hidepth"] = hidepth
    sample.stats_dfs.s3["clusters_hidepth"] = hidepth

    # get clusters file
    clusters = gzip.open(sample.files.clusters, 'rb')
    pairdealer = izip(*[iter(clusters)] * 2)

    # we subsample, else ... (could e.g., use first 10000 loci).
    # limit maxlen b/c some ref clusters can create huge contigs
    hidepth = min(10000, hidepth)
    maxlen = min(150, maxlen)
    dims = (hidepth, maxlen, 4)
    stacked = np.zeros(dims, dtype=np.uint64)

    # don't use sequence edges / restriction overhangs
    cutlens = [None, None]
    try:
        cutlens[0] = len(data.params.restriction_overhang[0])
        cutlens[1] = maxlen - len(data.params.restriction_overhang[1])
    except TypeError:
        ## Raised if either restriction_overhang parameters is an int
        pass
    except IndexError:
        ## If you don't include ro sequences then just carries on
        pass

    # fill stacked
    nclust = 0
    done = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError(
                "  clustfile formatting error in {}".format(chunk))

        if chunk:
            piece = chunk[0].decode().strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]
            # pull replicate read info from seqs
            reps = [int(sname.split("=")[-1][:-2]) for sname in names]

            ## get all reps
            sseqs = [list(seq) for seq in seqs]
            arrayed = np.concatenate([
                [seq] * rep for seq, rep in zip(sseqs, reps)
            ])

            ## enforce minimum depth for estimates
            if arrayed.shape[0] >= data.params.mindepth_statistical:
                # remove edge columns and select only the first 500
                # derep reads, just like in step 5
                arrayed = arrayed[:500, cutlens[0]:cutlens[1]]
                # remove cols that are pair separator
                arrayed = arrayed[:, ~np.any(arrayed == "n", axis=0)]
                # remove cols that are all Ns after converting -s to Ns
                arrayed[arrayed == "-"] = "N"
                arrayed = arrayed[:, ~np.all(arrayed == "N", axis=0)]
                # store in stacked dict

                catg = np.array(
                    [np.sum(arrayed == i, axis=0) for i in list("CATG")],
                    dtype=np.uint64).T

                ## Ensure catg honors the maxlen setting. If not you get a nasty
                ## broadcast error.
                stacked[nclust, :catg.shape[0], :] = catg[:maxlen, :]
                nclust += 1
        # bail out when nclusts have been done
        if nclust == hidepth:
            done = True

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
        if np.isnan(bfreqs).any():
            raise IPyradError(
                "Bad stack in getfreqs; {} {}"
                .format(sample.name, bfreqs))

        ## put into array, count array items as Byte strings
        tstack = Counter([j.tobytes() for j in stacked])

        ## get keys back as arrays and store vals as separate arrays
        ustacks = np.array(
            [np.frombuffer(i, dtype=np.uint64) for i in tstack.keys()]
        )
        counts = np.array(list(tstack.values()))
        ## cleanup
        del tstack

        ## if data are haploid fix H to 0
        if int(data.params.max_alleles_consens) == 1:
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
        # ip.logger.debug(
            # "Found sample with no clusters hidepth - {}".format(sample.name))
        pass

    return hetero, errors, success



def sample_cleanup(sample, hest, eest, success):
    "Store results to the Sample objects"
    # sample summary assignments
    sample.stats.state = 4
    sample.stats.hetero_est = float(hest)
    sample.stats.error_est = float(eest)

    # sample full assigments
    sample.stats_dfs.s4.hetero_est = float(hest)
    sample.stats_dfs.s4.error_est = float(eest)

    # In rare cases no hidepth clusters for statistical basecalling
    # so we warn the user, but carry on with default values
    if not success:
        msg = ("""
        Info: Sample {} - No clusters have sufficient depth for statistical
        basecalling. Setting default heterozygosity/error to 0.01/0.001.
        """.format(sample.name))
        print(msg)
