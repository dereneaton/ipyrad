#!/usr/bin/env python

"popgen tools"

from __future__ import print_function, division
from itertools import chain

import os
import h5py
import itertools
import math
import numpy as np
import pandas as pd
from ipyrad.analysis.utils import Params
from ipyrad.assemble.utils import IPyradError
from ipyrad import Assembly

_BAD_IMAP_ERROR = """
Samples in imap not in the hdf5 file: {}"
"""

_SKIP_SAMPLES_WARN = """
Skipping samples in hdf5 not present in imap: {}"
"""

class Popgen(object):
    """
    Analysis functions for calculating theta, Fst, Fis, thetaW, etc.

    Some functions follow Ferretti et al 2012 for calculating stats while
    accounting for missing data:

    Ferretti, L., Raineri, E., & Ramos-Onsins, S. (2012). Neutrality tests for
    sequences with missing data. Genetics, 191(4), 1397-1401.

    Another useful resource for calculating sumstats with missing data:

    Korunes, K., & Samuk, K. (2021). pixy: Unbiased estimation of nucleotide
    diversity and divergence in the presence of missing data. Molecular Ecology
    Resources.
    """

    def __init__(
        self,
        data,
        imap=None,
        minmap=None,
        workdir="analysis-popgen",
        quiet=False,
        ):
        
        # set attributes
        self.imap = (imap if imap else {})
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})
        self.npops = (len(self.imap) if imap else 1)
        self.quiet = quiet
        self.params = Params()
        self.nboots = 100

        # i/o paths
        self.workdir=workdir
        self.mapfile = ""
        # Sets the self.datafile attribute
        self._check_files(data)
        # Sets the self.samples attribute
        self._check_samples()
        #self.data = np.zeros()
        #self.maparr = np.zeros()

        # results dataframes
        self.results = Params()

        # pairwise Fst between all populations
        arrfst = np.zeros((self.npops, self.npops), dtype=np.uint64)
        self.results.fst = pd.DataFrame(
            arrfst
            )

        # individual pi 
        nsamples = len(list(chain(*self.imap.values())))
        arrpi = np.zeros(nsamples, dtype=np.uint64)
        self.results.pi = pd.DataFrame(
            arrpi
            )

        # population thetas
        npops = len(self.imap)
        arrtheta = np.zeros(npops, dtype=np.uint64)
        self.results.theta = pd.DataFrame(
            arrtheta
            )


    def _check_files(self, data):
        "check input files and file paths"

        if isinstance(data, Assembly):
            try:
                # Since v0.9.63-ish the snps_database is stored as an
                # assembly parameter. If it's not there try to construct
                # it from data.outfiles
                self.datafile = data.snps_database
            except AttributeError:
                self.datafile = os.path.join(data.dirs.outfiles + 
                                            f"{data.name}.snps.hdf5")
        else:
            self.datafile = data

        # check data file
        if os.path.exists(self.datafile):
            self.datafile = os.path.realpath(self.datafile)
        else:
            raise IPyradError(
                "data file does not exist. Check path: {}"
                .format(self.datafile))

        # check map file
        #if self.mapfile:
        #    self.mapfile = os.path.realpath(self.mapfile)

        # check workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)


    def _check_samples(self):
        "Read in list of sample names from the datafile"
        
        with h5py.File(self.datafile, 'r') as io5:
            self.samples = [x.decode('utf-8') for x in io5["snps"].attrs["names"]]
        if self.imap:
            # Check agreement between samples in imap and hdf5 file
            imap_samps = list(chain(*self.imap.values()))
            in_imap_not_hdf5 = set(imap_samps).difference(self.samples)
            in_hdf5_not_imap = set(self.samples).difference(imap_samps)

            if in_imap_not_hdf5:
                # Error if you pass in a sample in imap not in hdf5
                raise IPyradError(_BAD_IMAP_ERROR.format(in_imap_not_hdf5))
            if in_hdf5_not_imap:
                if not self.quiet:
                    # Warn here because this is a valid way to remove samples
                    print(_SKIP_SAMPLES_WARN.format(in_hdf5_not_imap))


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
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        darr = pd.DataFrame(
            data=np.zeros((self.npops, self.npops)),
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        narr = pd.DataFrame(
            data=np.zeros((self.npops, self.npops)),
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        d = self.npops

        # iterate over pairs of pops and fill Fst values
        pairs = itertools.combinations(self.imap.keys(), 2)
        for (pop1, pop2) in pairs:
            pop1idx = self.imap[pop1]
            pop2idx = self.imap[pop2]
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
        farr.columns = range(len(self.imap))
        darr.columns = range(len(self.imap))
        narr.columns = range(len(self.imap))
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
        a1 = sum([1./x for x in range(1, n)])
        c1 = b1 - (1./a1)
        e1 = c1/a1
        a2 = sum([1./(x**2) for x in range(1, n)])
        b2 = (2.*(n**2 + n + 3))/(9*n*(n-1))
        c2 = b2 - (n+2)/(a1*n) + (a2/(a1**2))
        e2 = c2/(a1**2+a2)
        ddenom = math.sqrt(e1*S + e2*S*(S-1))
        return ddenom

