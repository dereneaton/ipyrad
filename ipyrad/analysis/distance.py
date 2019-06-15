#!/usr/bin/env python

# py2/3 compat
from __future__ import print_function
from builtins import range

import numpy as np
import pandas as pd
from .snps_extracter import SNPsExtracter
from .snps_imputer import SNPsImputer
from .utils import jsubsample_snps


class Distance(object):
    """
    Calculate pairwise distance between samples with options for imputing 
    missing data based on population information, and for calculating 
    distances from SNPs or sequences and using a simple distance or using 
    a substitution model (JC, HKY). Also can exclude sites based on missing
    data.
    """
    def __init__(
        self, 
        data, 
        imap, 
        minmap, 
        mincov=0.0, 
        impute_method="sample",
        subsample_snps=False,
        # require_nonnegative_symmetric=False,
        random_seed=None,
        quiet=False,
        ):

        # store attributes
        np.random.seed(random_seed)
        self.data = data
        self.imap = imap
        self.minmap = minmap
        self.mincov = mincov
        self.impute_method = impute_method
        self.quiet = quiet
        self.subsample_snps = subsample_snps
        # self.require_nonnegative_symmetric = require_nonnegative_symmetric
        
        # get tools 
        self._se = SNPsExtracter(
            self.data, self.imap, self.minmap, self.mincov, self.quiet)
        self.names = self._se.names


    def run(self):

        # parse the data file to snps and map
        self._se.parse_genos_from_hdf5()

        # subsample one SNP per locus
        if self.subsample_snps:
            mask = jsubsample_snps(self._se.snpsmap, np.random.randint(0, 1e9))
            self._se.snps = self._se.snps[:, mask]

        # impute missing
        newsnps = SNPsImputer(
            self._se.snps, 
            self._se.names, 
            self._se.imap, 
            self.impute_method, 
            self.quiet).run()

        # get distances
        self.dists = pd.DataFrame(
            self._simple_dist(newsnps),
            index=self._se.names,
            columns=self._se.names,
        )


    def _simple_dist(self, snps):
        """
        Pairwise p-distance based on n diffs after missing data imputed.
        """
        ninds = snps.shape[0]
        diffs = np.zeros((ninds, ninds))

        for i in range(ninds):
            for j in range(ninds):
                diffs[i, j] = np.mean(np.abs(snps[i] - snps[j]) ** 2)
                diffs[j, i] = diffs[i, j]
        return diffs

