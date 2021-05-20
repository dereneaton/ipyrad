#!/usr/bin/env python

"""
Jointly estimate heterozygosity and error rate
"""

import os
import gzip
from collections import Counter
from itertools import combinations
from loguru import logger
from scipy.optimize import minimize
import scipy.stats
import pandas as pd
import numpy as np
import numba

from ipyrad.assemble.base_step import BaseStep
from ipyrad.core.schema import Stats4
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import IPyradError, clustdealer


class Step4(BaseStep):
    """
    Run the step4 estimation 
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 4, quiet, force)       
        self.haploid = data.params.max_alleles_consens == 1
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()

    def run(self):
        """
        Distribute optimization jobs and store results.
        """
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                optim2, 
                *(self.data, self.samples[sname], self.haploid)
            )
        msg = "inferring [H, E]"
        prog = AssemblyProgressBar(jobs, msg, 4, self.quiet)
        prog.block()
        prog.check()

        # collect updated samples and save to JSON
        for sname in prog.results:
            sample = prog.results[sname]
            self.data.samples[sname] = sample
        self.data.save_json()

        # write to stats file
        statsdf = pd.DataFrame(
            index=sorted(self.data.samples),
            columns=["hetero_est", "error_est"],
        )
        for sname in self.data.samples:
            stats = self.data.samples[sname].stats_s4.dict()
            for i in statsdf.columns:
                statsdf.loc[sname, i] = stats[i]
        logger.info("\n" + statsdf.to_string())
        outfile = os.path.join(self.stepdir, "s4_joint_estimate.txt")
        with open(outfile, 'w') as out:
            out.write(statsdf.to_string())


def optim2(data, sample, haploid):
    """
    Maximum likelihood optimization with scipy
    """

    # get array of all clusters data
    stacked = stackarray(data, sample)

    # get base frequencies
    bfreqs = stacked.sum(axis=0) / float(stacked.sum())
    if np.isnan(bfreqs).any():
        raise IPyradError(
            f"Bad stack in getfreqs; {sample.name} {bfreqs}")

    # put into array, count array items as Byte strings
    tstack = Counter([j.tostring() for j in stacked])

    # get keys back as arrays and store vals as separate arrays
    ustacks = np.array(
        [np.frombuffer(i, dtype=np.uint64) for i in tstack.keys()]
    )
    counts = np.array(list(tstack.values()))

    # cleanup
    del tstack

    if haploid:
        # fit haploid model
        fit = minimize(
            get_haploid_loglik, 
            x0=(0.001,),
            args=(bfreqs, ustacks, counts),
            method="L-BFGS-B",
            bounds=[1e-10, 1.0],
        )
        hetero = 0.0
        error = fit.x[0]
    else:
        # fit haploid model
        fit = minimize(
            nget_diploid_loglik, 
            x0=(0.01, 0.001),
            args=(bfreqs, ustacks, counts),
            method="L-BFGS-B",            
            bounds=[(1e-10, 1.0), (1e-10, 1.0)],
        )
        hetero, error = fit.x

    sample.state = 4
    sample.stats_s4 = Stats4(
        hetero_est=hetero,
        error_est=error,
        min_depth_stat_during_step4=data.params.min_depth_statistical,
    )
    return sample


def get_haploid_loglik(errors, bfreqs, ustacks, counts):
    """
    Log likelihood score given values [E]
    """
    hetero = 0.
    lik1 = ((1. - hetero) * likelihood1(errors, bfreqs, ustacks))
    liks = lik1
    logliks = np.log(liks[liks > 0]) * counts[liks > 0]
    score = -logliks.sum()
    return score


def nget_diploid_loglik(pstart, bfreqs, ustacks, counts):
    """
    Log likelihood score given values [H,E]
    """
    hetero, errors = pstart
    lik1 = (1. - hetero) * likelihood1(errors, bfreqs, ustacks)
    lik2 = (hetero) * nlikelihood2(errors, bfreqs, ustacks)
    liks = lik1 + lik2
    logliks = np.log(liks[liks > 0]) * counts[liks > 0]
    score = -logliks.sum()
    return score


def likelihood1(errors, bfreqs, ustacks):
    """
    Probability homozygous. 
    """
    ## make sure base_frequencies are in the right order
    # print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    # totals = np.array([ustacks.sum(axis=1)]*4).T
    totals = np.array([ustacks.sum(axis=1)] * 4).T
    prob = scipy.stats.binom.pmf(totals - ustacks, totals, errors)
    lik1 = np.sum(bfreqs * prob, axis=1)
    return lik1


def nlikelihood2(errors, bfreqs, ustacks):
    """
    calls nblik2_build and lik2_calc for a given err
    """
    one = [2. * bfreqs[i] * bfreqs[j] for i, j in combinations(range(4), 2)]
    four = 1. - np.sum(bfreqs**2)
    tots, twos, thrs = nblik2_build(ustacks)
    res2 = lik2_calc(errors, one, tots, twos, thrs, four)
    return res2


@numba.jit(nopython=True)
def nblik2_build(ustacks):
    """
    JIT'd function builds array that can be used to calc binom pmf
    """
    # fill for pmf later
    tots = np.empty((ustacks.shape[0], 1))
    twos = np.empty((ustacks.shape[0], 6))
    thrs = np.empty((ustacks.shape[0], 6, 2))

    # fill big arrays
    for idx in range(ustacks.shape[0]):

        ust = ustacks[idx]
        tot = ust.sum()
        tots[idx] = tot

        # fast filling of arrays
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
    # calculate twos
    _twos = scipy.stats.binom.pmf(twos, tots, 0.5)

    # calculate threes
    _thrs = thrs.reshape(thrs.shape[0] * thrs.shape[1], thrs.shape[2])
    _thrs = scipy.stats.binom.pmf(_thrs[:, 0], _thrs[:, 1], (2. * err) / 3.)
    _thrs = _thrs.reshape(thrs.shape[0], 6)

    # calculate return sums
    return np.sum(one * _twos * (_thrs / four), axis=1)


def recal_hidepth_stat(data, sample):
    """
    Returns a mask for loci above the statistical threshold and the 
    max fragment length that will be allowed.
    """
    # if nothing changed then return max fragment length
    check1 = data.params.min_depth_majrule == sample.stats_s3.min_depth_maj_during_step3
    check2 = data.params.min_depth_statistical == sample.stats_s3.min_depth_stat_during_step3
    if check1 and check2:
        maxfrag = int(
            4 + sample.stats_s3.mean_hidepth_cluster_length + 
            (2 * sample.stats_s3.std_hidepth_cluster_length)
        )
        return None, maxfrag

    # otherwise calculate depth again given the new mindepths settings.
    with gzip.open(sample.files.clusters, 'rt') as infile:
        pairdealer = zip(*[infile] * 2)

        # storage
        counts = []
        depths = []
        maxlen = []

        # start with cluster 0
        tdepth = 0
        tlen = 0
        tcount = 0

        # iterate until empty
        while 1:
            try:
                name, seq = next(pairdealer)
            except StopIteration:
                break

            # if at the end of a cluster
            if name == "//\n":
                depths.append(tdepth)
                maxlen.append(tlen)
                counts.append(tcount)
                tlen = 0
                tdepth = 0
                tcount = 0
            else:
                tdepth += int(name.strip().split("=")[-1][:-2])
                tlen = len(seq)
                tcount += 1

    # convert to arrays
    maxlens, depths = np.array(maxlen), np.array(depths)

    # get mask of clusters that are hidepth
    stat_mask = depths >= data.params.min_depth_statistical

    # get frag lenths of clusters that are hidepth
    lens_above_st = maxlens[stat_mask]

    # calculate frag length from hidepth lens
    try:       
        maxfrag = int(4 + lens_above_st.mean() + (2. * lens_above_st.std()))
    except Exception as inst:
        raise IPyradError(
            "No clusts with depth sufficient for statistical basecalling. "
            f"I recommend you branch to drop this sample: {sample.name}"
            ) from inst
    return stat_mask, maxfrag


def stackarray(data, sample):
    """
    Stacks clusters into arrays using at most 10K clusters.
    """
    # only use clusters with depth > min_depth_statistical for param estimates
    stat_mask, maxfrag = recal_hidepth_stat(data, sample)
    stat_mask = stat_mask if stat_mask is not None else np.ones(10000)

    # get clusters file
    clustio = gzip.open(sample.files.clusters, 'rt')
    pairdealer = zip(*[clustio] * 2)

    # we subsample, else ... (could e.g., use first 10000 loci).
    # limit maxlen b/c some ref clusters can create huge contigs
    maxclusts = min(10000, stat_mask.sum())
    maxfrag = min(150, maxfrag)
    dims = (maxclusts, maxfrag, 4)
    stacked = np.zeros(dims, dtype=np.uint64)

    # don't use sequence edges / restriction overhangs
    cutlens = [None, None]
    try:
        cutlens[0] = len(data.params.restriction_overhang[0])
        cutlens[1] = maxfrag - len(data.params.restriction_overhang[1])
    except TypeError:
        ## Raised if either restriction_overhang parameters is an int
        pass
    except IndexError:
        ## If you don't include ro sequences then just carries on
        pass

    # fill stacked
    nclust = 0
    done = 0
    while 1:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except StopIteration:
            break

        if chunk:
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = piece[1::2]

            # pull replicate read info from seqs
            if data.hackers.declone_PCR_duplicates:
                reps = [1 for i in names]
            else:
                reps = [int(sname.split("=")[-1][:-2]) for sname in names]

            # get all reps
            sseqs = [list(seq) for seq in seqs]
            arrayed = np.concatenate([
                [seq] * rep for seq, rep in zip(sseqs, reps)
            ])

            # enforce minimum depth for estimates
            if arrayed.shape[0] >= data.params.min_depth_statistical:
                # select at most random 500
                if arrayed.shape[0] > 500:
                    idxs = np.random.choice(
                        range(arrayed.shape[0]), 
                        size=500, 
                        replace=False,
                    )
                    arrayed = arrayed[idxs]
                
                # trim cut segments
                arrayed = arrayed[:, cutlens[0]:cutlens[1]]               
                
                # remove cols that are pair separator or all Ns
                arrayed = arrayed[:, ~np.any(arrayed == "n", axis=0)]
                arrayed[arrayed == "-"] = "N"
                arrayed = arrayed[:, ~np.all(arrayed == "N", axis=0)]

                # store in stack
                catg = np.array([
                    np.sum(arrayed == i, axis=0) for i in list("CATG")
                    ],
                    dtype=np.uint64
                ).T
                # Ensure catg honors the maxlen setting. If not you get a nasty
                # broadcast error.
                stacked[nclust, :catg.shape[0], :] = catg[:maxfrag, :]
                nclust += 1

        if done or (nclust == maxclusts):
            break

    clustio.close()
    # drop the empty rows in case there are fewer loci than the size of array
    newstack = stacked[stacked.sum(axis=2) > 0]
    assert not np.any(newstack.sum(axis=1) == 0), "no zero rows"
    return newstack



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG", stderr=False, logfile="/tmp/test.log")
   
    TEST = ip.load_json("/tmp/TEST1.json")
    TEST.run("4", force=True, quiet=False)

