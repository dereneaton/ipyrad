#!/usr/bin/env python

"popgen tools"

from __future__ import print_function, division
from itertools import chain

import os
import itertools
import math
import numpy as np
import pandas as pd
from ipyrad.analysis.utils import Params
from ipyrad.assemble.utils import IPyradError


class Popgen(object):
    "Analysis functions for calculating theta, Fst, Fis, thetaW, etc."

    def __init__(self, name, data, workdir, mapfile=None):
        
        # i/o paths
        self.workdir = workdir
        self._datafile = data
        self._mapfile = mapfile
        self._check_files()
        self.data = np.zeros()
        self.maparr = np.zeros()

        # init default param settings
        self.params = Params()
        self.popdict = {}
        self.mindict = {}
        self.npops = len(self.popdict)
        self.nboots = 100

        # results dataframes
        self.results = Params()

        # pairwise Fst between all populations
        npops = len(self.popdict)
        arrfst = np.zeros((npops, npops), dtype=np.uint64)
        self.results.fst = pd.DataFrame(
            arrfst
            )

        # individual pi 
        nsamples = len(list(chain(*self.popdict.values())))
        arrpi = np.zeros(nsamples, dtype=np.uint64)
        self.results.pi = pd.DataFrame(
            arrpi
            )

        # population thetas
        npops = len(self.popdict)
        arrtheta = np.zeros(npops, dtype=np.uint64)
        self.results.theta = pd.DataFrame(
            arrtheta
            )

        # parse samples from the data file
        self._check_files()



    def _check_files(self):
        "check input files and file paths"
        # check data file
        if os.path.exists(self.datafile):
            self.datafile = os.path.realpath(self.datafile)
        else:
            raise IPyradError(
                "data file does not exist. Check path: {}"
                .format(self.datafile))

        # check map file
        if self.mapfile:
            self.mapfile = os.path.realpath(self.mapfile)

        # check workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        


    def run(self, force=False, ipyclient=False):
        "calculate the given statistic"
        pass


    def _fst(self):
        """
        Calculate population fixation index Fst using Hudson's estimator.
        Hudson et al. (1992) "Estimation of Levels of Gene Flow From DNA 
        Sequence Data", also returns Fstd with correction for the number of 
        subpopulations (using the number sampled, since the true number is 
        unknown) and number of migrants (Nm) derived from Li (1976b) as 
        described in the Hudson et al. (1992) paper. 
        """       
        # init fst matrix df with named rows and cols
        farr = pd.DataFrame(
            data=np.zeros((self.npops, self.npops)),
            index=self.popdict.keys(),
            columns=self.popdict.keys(),
        ) 
        darr = pd.DataFrame(
            data=np.zeros((self.npops, self.npops)),
            index=self.popdict.keys(),
            columns=self.popdict.keys(),
        ) 
        narr = pd.DataFrame(
            data=np.zeros((self.npops, self.npops)),
            index=self.popdict.keys(),
            columns=self.popdict.keys(),
        ) 
        d = self.npops

        # iterate over pairs of pops and fill Fst values
        pairs = itertools.combinations(self.popdict.keys(), 2)
        for (pop1, pop2) in pairs:
            pop1idx = self.popdict[pop1]
            pop2idx = self.popdict[pop2]
            popaidx = pop1idx + pop2idx

            within1 = list(itertools.combinations(pop1idx, 2))
            within2 = list(itertools.combinations(pop2idx, 2))
            withins = within1 + within2
            allpairs = itertools.combinations(popaidx, 2)
            betweens = itertools.filterfalse(
                lambda x: bool(x in withins),
                allpairs
            )

            diff = [self.data[i] != self.data[j] for (i, j) in withins]
            sums = np.sum(diff, axis=0)
            a = sums / sums.shape[0]
            
            diff = [self.data[i] != self.data[j] for (i, j) in betweens]
            sums = np.sum(diff, axis=0)
            b = sums / sums.shape[0]
            
            farr.loc[pop1, pop2] = abs(1 - (a / b))
            narr.loc[pop1, pop2] = (((d - 1) / d) * (1 / 2) * (a / b - a))
            darr.loc[pop1, pop2] = (
                abs(1 - (a / ((1 / d) * a) + (((d - 1) / d) * b))))
        farr.columns = range(len(self.popdict))
        darr.columns = range(len(self.popdict))
        narr.columns = range(len(self.popdict))
        return farr, darr, narr



    def _fis(self):
        "calculate population inbreeding Fis after filtering"
        pass


    def _pi(self):
        "calculate per-sample heterozygosity after filtering"
        pass


    def _filter_data(self):
        "take input data as phylip seq array and subsample loci from mapfile"
        pass


    def _TajimaD_denom(n, S):
        """
        Tajima's D denominator. I toiled over this to get it right and it is
        known to be working.

        This page has a nice worked example with values for each
        subfunction so you can check your equations:
        https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-
        quantitative-genomics-fall-2005/study-materials/tajimad1.pdf

        :param int N: The number of samples
        :param int S: The number of segregating sites.
        """
        b1 = (n+1)/float(3*(n-1))
        a1 = sum([1./x for x in xrange(1, n)])
        c1 = b1 - (1./a1)
        e1 = c1/a1
        a2 = sum([1./(x**2) for x in xrange(1, n)])
        b2 = (2.*(n**2 + n + 3))/(9*n*(n-1))
        c2 = b2 - (n+2)/(a1*n) + (a2/(a1**2))
        e2 = c2/(a1**2+a2)
        ddenom = math.sqrt(e1*S + e2*S*(S-1))
        return ddenom

