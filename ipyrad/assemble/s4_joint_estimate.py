#!/usr/bin/env python

"""Jointly estimate heterozygosity and error rate.
"""

from typing import TypeVar, Tuple
from itertools import combinations
from loguru import logger
from scipy.optimize import minimize
import scipy.stats
import pandas as pd
import numpy as np
import numba

from ipyrad.assemble.base_step import BaseStep
from ipyrad.core.schema import Stats4, Stats5
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import IPyradError, NoHighDepthClustersError
from ipyrad.assemble.clustmap_within_both import iter_clusters

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")
logger = logger.bind(name="ipyrad")


class Step4(BaseStep):
    """Run the step4 estimation """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 4, quiet, force)
        self.haploid = data.params.max_alleles_consens == 1
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()

    def run(self):
        """Distribute optimization jobs and store results."""
        jobs = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample, self.haploid)
            jobs[sname] = self.lbview.apply(optim2, *args)
        msg = "inferring [H, E]"
        prog = AssemblyProgressBar(jobs, msg, 4, self.quiet)
        prog.block()
        prog.check()

        # collect updated samples and save to JSON
        for sname, result in prog.results.items():
            sample = self.data.samples[sname]
            sample.state = 4
            sample._clear_old_results()
            sample.stats_s5 = None # sample.stats_s5 = Stats5()
            sample.stats_s4 = Stats4(
                hetero_est=result[0],
                error_est=result[1],
                min_depth_stat_during_step4=self.data.params.min_depth_statistical,
            )
        self.data.save_json()

        # write to stats file
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=["hetero_est", "error_est"],
        )

        # update samples on the Assembly object
        for sname in self.samples:
            stats = self.data.samples[sname].stats_s4.dict()
            for i in statsdf.columns:
                statsdf.loc[sname, i] = stats[i]

        # log and save stats
        logger.info("\n" + statsdf.to_string())
        outfile = self.data.stepdir / "s4_joint_estimate.txt"
        with open(outfile, 'w', encoding="utf-8") as out:
            out.write(statsdf.to_string())


def optim2(data: Assembly, sample: Sample, haploid: bool) -> Tuple[float, float]:
    """Maximum likelihood optimization with scipy."""

    # get array of all clusters data: (maxclusts, maxlen, 4)
    try:
        stacked = get_stack_array(data, sample)

    # if no high depth cluster the sample can still proceed, but it
    # will not contribute to estimating the avg H,E and will receive
    # only low depth calls in step 5 unless params are changed.
    except NoHighDepthClustersError:
        return np.nan, np.nan

    # get base frequencies
    bfreqs = stacked.sum(axis=0) / float(stacked.sum())

    # count each unique site count pattern
    ustacks, counts = np.unique(stacked, axis=0, return_counts=True)

    # fit haploid or diploid model to counts
    if haploid:
        fit = minimize(
            get_haploid_loglik, 
            x0=(0.001,),
            args=(bfreqs, ustacks, counts),
            method="L-BFGS-B",
            bounds=[1e-6, 0.1],
        )
        hetero = 0.0
        error = fit.x[0]
    else:
        fit = minimize(
            nget_diploid_loglik, 
            x0=(0.01, 0.001),
            args=(bfreqs, ustacks, counts),
            method="L-BFGS-B",
            bounds=[(1e-6, 0.1), (1e-6, 0.1)],
        )
        hetero, error = fit.x
    return hetero, error


def get_haploid_loglik(errors, bfreqs, ustacks, counts):
    """Log likelihood score given values [E]."""
    hetero = 0.
    lik1 = ((1. - hetero) * likelihood1(errors, bfreqs, ustacks))
    liks = lik1
    logliks = np.log(liks[liks > 0]) * counts[liks > 0]
    score = -logliks.sum()
    return score


def nget_diploid_loglik(
    pstart: Tuple[float, float], 
    bfreqs: np.ndarray, 
    ustacks: np.ndarray, 
    counts: np.ndarray,
) -> float:
    """Return Log likelihood score given values [H,E]"""
    hetero, errors = pstart
    lik1 = (1. - hetero) * likelihood1(errors, bfreqs, ustacks)
    lik2 = (hetero) * nlikelihood2(errors, bfreqs, ustacks)
    liks = lik1 + lik2
    logliks = np.log(liks[liks > 0]) * counts[liks > 0]
    score = -logliks.sum()
    return score


def likelihood1(errors, bfreqs, ustacks):
    """Probability homozygous."""
    # make sure base_frequencies are in the right order
    # print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    # totals = np.array([ustacks.sum(axis=1)]*4).T
    totals = np.array([ustacks.sum(axis=1)] * 4).T
    prob = scipy.stats.binom.pmf(totals - ustacks, totals, errors)
    lik1 = np.sum(bfreqs * prob, axis=1)
    return lik1


def nlikelihood2(errors, bfreqs, ustacks):
    """Calls nblik2_build and lik2_calc for a given err."""
    one = [2. * bfreqs[i] * bfreqs[j] for i, j in combinations(range(4), 2)]
    four = 1. - np.sum(bfreqs**2)
    tots, twos, thrs = nblik2_build(ustacks)
    res2 = lik2_calc(errors, one, tots, twos, thrs, four)
    return res2


@numba.jit(nopython=True)
def nblik2_build(ustacks):
    """JIT'd function builds array that can be used to calc binom pmf
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


def recal_hidepth_cluster_stats(
    data: Assembly, sample: Sample, majrule: bool = False,
) -> Tuple[np.ndarray, int]:
    """Return a mask for cluster depths, and the max frag length.

    This is useful to run first to get a sense of the depths and lens
    given the current mindepth param settings.

    Note: this func is used in both steps 4 and 5.
    """
    # otherwise calculate depth again given the new mindepths settings.
    depths = []   # read depth: sum of 'sizes'
    clens = []    # lengths of clusters
    for clust in iter_clusters(sample.files.clusters, gzipped=True):
        names = clust[::2]
        sizes = [int(i.split(";")[-2][5:]) for i in names]
        depths.append(sum(sizes))
        clens.append(len(clust[1].strip()))
    clens, depths = np.array(clens), np.array(depths)

    # get mask of clusters that are hidepth
    if majrule:
        keep = depths >= data.params.min_depth_majrule
    else:
        keep = depths >= data.params.min_depth_statistical

    # get frag lenths of clusters that are hidepth
    lens_above_st = clens[keep]
    # print(f"{sample.name}, {keep.shape}, {depths} {depths >=data.params.min_depth_majrule} {data.params.min_depth_majrule} {lens_above_st}, {clens}")

    # calculate frag length from hidepth lens
    try:
        maxfrag = int(4 + lens_above_st.mean() + (2. * lens_above_st.std()))
    except Exception as inst:
        # this exception will raise in step 4 and be caught to print an
        # warning message and then will set the samples H,E estimates to
        # nan. In step 5 the nans will be caught above...
        print(
            f"sample {sample.name} has no clusters above the "
            "`min_depth_statistical` parameter setting, and thus will "
            "include only low depth base calls in step 5.")
        raise NoHighDepthClustersError(f"{sample.name}") from inst
    return keep, maxfrag


def get_stack_array(data: Assembly, sample: Sample, size: int = 10_000) -> np.ndarray:
    """Stacks clusters into arrays using at most 10K clusters.

    Uses maxlen to limit the end of arrays, and also masks the first
    and last 6 bp from each read since these are more prone to 
    alignmentn errors in denovo assemblies are will likely be 
    trimmed later.
    """
    # only use clusters with depth > min_depth_statistical for param estimates
    stat_mask, maxfrag = recal_hidepth_cluster_stats(data, sample)

    # sample many (e.g., 10_000) clusters to use for param estimation.
    maxclusts = min(size, stat_mask.sum())
    maxfrag = min(150, maxfrag)
    dims = (maxclusts, maxfrag, 4)
    stacked = np.zeros(dims, dtype=np.uint64)

    # fill stacked
    clustgen = iter_clusters(sample.files.clusters, gzipped=True)
    sidx = 0  # stored row number
    for idx, clust in enumerate(clustgen):

        # skip masked (lowdepth) clusters
        if not stat_mask[idx]:
            continue

        # if maxclusts are stored then no need to do more.
        if sidx >= maxclusts:
            continue

        # parse cluster and expand derep depths
        names = clust[0::2]
        seqs = clust[1::2]
        reps = [int(i.split(";")[-2][5:]) for i in names]
        sseqs = [list(i.strip()) for i in seqs]
        arr = np.concatenate([[seq] * rep for seq, rep in zip(sseqs, reps)])

        # select at most random 500 reads in a cluster
        if arr.shape[0] > 500:
            ridxs = range(arr.shape[0])
            ridxs = np.random.choice(ridxs, size=500, replace=False)
            arr = arr[ridxs]
                
        # mask edges, indels, and pair inserts and remove empty columns.
        arr[:, :8] = "N"
        arr[:, -8:] = "N"
        arr[arr == "-"] = "N"
        arr[:, np.any(arr == "n", axis=0)] = "N"
        arr = arr[:, ~np.all(arr == "N", axis=0)]

        # store in stack shape=(nsites, 4)
        catg = [np.sum(arr == i, axis=0) for i in list("CATG")]
        catg = np.array(catg, dtype=np.uint64).T

        # limit stored catg data to maxfrag len
        stacked[sidx, :catg.shape[0], :] = catg[:maxfrag, :]
        sidx += 1

    # drop the empty rows in case there are fewer loci than the size of array
    newstack = stacked[stacked.sum(axis=2) > 0]
    assert not np.any(newstack.sum(axis=1) == 0), "no zero rows"
    return newstack


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")#, log_file="/tmp/test.log")

    TEST = ip.load_json("../../pedtest/NEW.json")
    TEST.run("4", force=True, quiet=False)
   
    # TEST = ip.load_json("/tmp/TEST3.json")
    # TEST.run("4", force=True, quiet=False)
