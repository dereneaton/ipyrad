#!/usr/bin/env python

"""
Extra utilities for ipa.bpp analyses.
"""

from loguru import logger
import numpy as np
import scipy.stats as stats
from scipy.optimize import minimize


class PriorHelper:
    """

    """
    def __init__(self, gentime_mean, gentime_std, mutrate_mean, mutrate_std):
        self.gentime_mean = gentime_mean
        self.gentime_std = gentime_std
        self.mutrate_mean = mutrate_mean
        self.mutrate_std = mutrate_std

    def get_invgamma_tau_params(self, age_mean, age_std):
        """

        """
        return get_invgamma_params_for_tau(
            age_mean, age_std,
            self.gentime_mean, self.gentime_std,
            self.mutrate_mean, self.mutrate_std,
        )

    def get_invgamma_theta_params(self, ne_mean, ne_std):
        """

        """
        return get_invgamma_params_for_theta(
            ne_mean, ne_std, 
            self.mutrate_mean, self.mutrate_std,
        )

    def draw(self, age_mean, age_std, ne_mean, ne_std):
        """

        """
        pass



def fit_invgamma_params_to_mean_std(mean, std):
    """
    Returns alpha and beta of invgamma distribution fit to a mean
    and std, while constraining alpha >= 2.
    """
    # to prevent overlow divide both numbers by a factor when large.
    factor = len(str(int(mean)))
    if factor > 5:
        factor = 10000
    else:
        factor = 1
    mean /= factor
    std /= factor

    def optim_func(params, mean, std):
        "minimize this function to fit mean, std of invgamma dist."
        # returns mean and variance
        hmean, hvar = stats.invgamma.stats(a=params[0], scale=params[1])
        diff1 = (hmean - mean)**2
        diff2 = (np.sqrt(hvar) - std)**2
        return diff1 + diff2

    # optimize
    model = minimize(
        optim_func,
        x0=(3., 1.),
        args=(mean, std),
        method="L-BFGS-B",
        bounds=[
            (2 + 1e-9, np.inf),
            (0 + 1e-9, np.inf),
        ],
    )

    # multiply scale parameter back by the multiplier
    model.x[1] *= factor
    logger.debug(
        f"shape={model.x[0]}, scale={model.x[1]}, mean={mean}, std={std}")
    return tuple(model.x)



def get_invgamma_params_for_tau(
    age_mean, 
    age_std, 
    gentime_mean=1,
    gentime_std=0.001, 
    mutrate_mean=5e-8,
    mutrate_std=5e-9,
    ):
    """
    Returns the shape and scale parameters for the average percent 
    sequence divergence between the root and tips of a phylogeny as a
    invgamma distributed random variable, given a mean and std for 
    the root node age in years, generation times in years, and 
    mutation rate in per-site-per-generation units.

    Parameters
    ----------
    age_mean (float):
    age_std (float): 
    gentime_mean (float; def=1):
    gentime_std (float; def=0.001): 
    mutrate_mean (float; def=5e-8):
    mutrate_std (float; def=5e-9):
    """
    # get invgamma params for node age 
    a_age, b_age = fit_invgamma_params_to_mean_std(age_mean, age_std)

    # get invgamma params for generation times
    a_gen, b_gen = fit_invgamma_params_to_mean_std(gentime_mean, gentime_std)

    # get invgamma params for mutation rates
    a_mut, b_mut = fit_invgamma_params_to_mean_std(mutrate_mean, mutrate_std)

    # get random variables sampled from each invgamma distribution
    dists = [
        stats.invgamma(a=a_age, scale=b_age).rvs(100000),                
        stats.invgamma(a=a_gen, scale=b_gen).rvs(100000),
        stats.invgamma(a=a_mut, scale=b_mut).rvs(100000),
    ]

    # calculate percent sequence divergence dist from rvs dists
    muts_per_site = (dists[0] * dists[2]) / dists[1]
    
    # get mean and std of this composite invgamma distribution
    s_mean = muts_per_site.mean()
    s_std = muts_per_site.std()

    # fit parameters of the invgamma dist for sequence divergence
    a_seq, b_seq = fit_invgamma_params_to_mean_std(s_mean, s_std)
    return a_seq, b_seq



def get_invgamma_params_for_theta(
    ne_mean, 
    ne_std, 
    mutrate_mean=5e-8,
    mutrate_std=5e-9,
    ):
    """
    Returns the shape and scale parameters for the population mutation 
    rate (Theta=4Neu) across a phylogeny as an invgamma distributed 
    random variable given a mean and std for the effective population
    size, and the per-site-per-generation mutation rate.

    Parameters
    ----------
    ne_mean (float):
    ne_std (float): 
    mutrate_mean (float; def=5e-8):
    mutrate_std (float; def=5e-9):
    """
    # get invgamma params for node age 
    a_ne, b_ne = fit_invgamma_params_to_mean_std(ne_mean, ne_std)

    # get invgamma params for mutation rates
    a_mut, b_mut = fit_invgamma_params_to_mean_std(mutrate_mean, mutrate_std)

    # get random variables sampled from each invgamma distribution
    dists = [
        stats.invgamma(a=a_ne, scale=b_ne).rvs(100000),                
        stats.invgamma(a=a_mut, scale=b_mut).rvs(100000),
    ]

    # calculate percent sequence divergence dist from rvs dists
    theta = dists[0] * dists[1] * 4
    
    # get mean and std of this composite invgamma distribution
    t_mean = theta.mean()
    t_std = theta.std()

    # fit parameters of the invgamma dist for sequence divergence
    a_theta, b_theta = fit_invgamma_params_to_mean_std(t_mean, t_std)
    return a_theta, b_theta



if __name__ == "__main__":

    # from an age estimate (e.g., fossil calibration) this gives you
    # the appropriate invgamma params to set an informative prior on
    # the node in superbpp.
    PARAMS = get_invgamma_params_for_tau(
        age_mean=10000, age_std=5000,
        gentime_mean=1, gentime_std=1,
        )
    print(PARAMS)

    # from a posterior distribution of tau values from the previous 
    # step in a trio analysis, this can convert a tau mean and std to
    # invgamma params to set an informative prior on the node in 
    # superbpp.
    PARAMS = fit_invgamma_params_to_mean_std(0.2, 0.02)
    print(PARAMS)
