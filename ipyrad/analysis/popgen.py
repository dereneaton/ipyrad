#!/usr/bin/env python

"popgen tools"

from __future__ import print_function, division
from itertools import chain

import os
import numpy as np
import pandas as pd
from ipyrad.analysis.utils import Params
from ipyrad.assemble.util import IPyradError


class Popgen(object):
    "Analysis functions for calculating theta, Fst, Fis, etc."

    def __init__(self, name, data, workdir, mapfile=None):
        
        # i/o paths
        self.workdir = workdir
        self.data = data
        self.mapfile = mapfile
        self._check_files()

        # init default param settings
        self.params = Params()
        self.params.popdict = {}
        self.params.mindict = {}
        self.params.nboots = 100

        # results dataframes
        self.results = Params()

        # pairwise Fst between all populations
        npops = len(self.params.popdict)
        arrfst = np.zeros((npops, npops), dtype=np.uint64)
        self.results.fst = pd.DataFrame(
            arrfst
            )

        # individual pi 
        nsamples = len(list(chain(*self.params.popdict.values())))
        arrpi = np.zeros(nsamples, dtype=np.uint64)
        self.results.pi = pd.DataFrame(
            arrpi
            )

        # population thetas
        npops = len(self.params.popdict)
        arrtheta = np.zeros(npops, dtype=np.uint64)
        self.results.theta = pd.DataFrame(
            arrtheta
            )

        # parse samples from the data file
        self._check_files()


        

    def run(self, force, ipyclient):
        "calculate the given statistic"
        pass


    def _fst(self):
        """
        Calculate population fixation index Fst after filtering.
        The variance of allele frequencies between populations.
        Wier and Cockerham (1984) estimator based on average number of 
        pairwise differences: 
           (pi_between - pi_within) / pi_between
        """       


    def _fis(self):
        "calculate population inbreeding Fis after filtering"
        pass


    def _pi(self):
        "calculate per-sample heterozygosity after filtering"
        pass


    def _filter_data(self):
        "take input data as phylip seq array and subsample loci from mapfile"
        pass


    def _check_files(self):
        "check input files and file paths"
        # check data file
        if os.path.exists(self.data):
            self.data = os.path.realpath(self.data)
        else:
            raise IPyradError(
                "data file does not exist. Check path: {}".format(self.data))

        # check map file
        if self.mapfile:
            self.mapfile = os.path.realpath(self.mapfile)

        # check workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
