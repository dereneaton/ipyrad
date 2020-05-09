#!/usr/bin/env python

""" Impute SNPs based on population allele frequencies"""

from __future__ import print_function, division

from copy import deepcopy
import numpy as np

_MISSING_SKLEARN = """
This ipyrad tool requires the library scikit-learn.
You can install it with the following command in a terminal.

conda install scikit-learn -c conda-forge 
"""


class SNPsImputer(object):
    """
    Imputation of missing data based on population allele frequencies. This
    tool is generally meant to be used internally by ipyrad within other 
    analysis tool methods. 

    Parameters:
    -----------
    data: (ndarray)
        A uint8 ndarray of genotype calls that have already been filtered 
        to remove any sites that are not wanted. Only 0,1,2,9 should be 
        in matrix.
    names: list
        Ordered names as extracted from the HDF5 database.
    imap: (dict; default=None)
        Dictionary to assign samples to groups to be used with minmap.
    impute_method: (str; default='sample')
        None, "sample", "simple", "kmeans"

    Functions:
    ----------
    ...
    """
    def __init__(
        self, 
        data, 
        names,
        imap=None,
        impute_method="sample",
        quiet=False,
        ):

        # init attributes
        self.quiet = quiet
        self.snps = deepcopy(data)
        self._mvals = np.sum(self.snps == 9)        

        # data attributes
        self.impute_method = impute_method
        self.imap = imap

        # get names from imap else we'll fill this later with all
        self.names = names


    def _print(self, msg):
        if not self.quiet:
            print(msg)

    def run(self):
        """
        Impute data in-place updating self.snps by filling missing (9) values.
        """
        if self.impute_method == "sample":
            self.snps = self._impute_sample()

        else:
            self.snps[self.snps == 9] = 0
            self._print(
                "Imputation: 'None'; (0, 1, 2) = {:.1f}%, {:.1f}%, {:.1f}%"
                .format(100, 0, 0)            
            )
        return self.snps


    def _impute_sample(self, imap=None):
        """
        Sample derived alleles by their frequency for each population and 
        assign to fill 9 in each column for each pop.
        """
        # override imap
        if not imap:
            imap = self.imap

        # impute data by mean value in each population
        newdata = self.snps.copy()
        for pop, samps in imap.items():

            # sample pop data
            sidxs = sorted(self.names.index(i) for i in samps)
            data = newdata[sidxs, :].copy()

            # number of alleles at each site that are not 9
            nallels = np.sum(data != 9, axis=0) * 2

            # get prob derived at each site using tmp array w/ missing to zero 
            tmp = data.copy()
            tmp[tmp == 9] = 0
            fderived = tmp.sum(axis=0) / nallels

            # sampler
            sampled = np.random.binomial(n=2, p=fderived, size=data.shape)
            data[data == 9] = sampled[data == 9]
            newdata[sidxs, :] = data

        # get all imputed values
        imputed = newdata[np.where(self.snps == 9)]
        self._print(
            "Imputation: 'sampled'; (0, 1, 2) = {:.1f}%, {:.1f}%, {:.1f}%"
            .format(
                100 * np.sum(imputed == 0) / imputed.size,
                100 * np.sum(imputed == 1) / imputed.size,
                100 * np.sum(imputed == 2) / imputed.size,
            )
        )
        return newdata



    def _impute_sample_hier(self, imap=None):
        """
        Sample derived alleles by their frequency for each population and 
        assign to fill 9 in each column for each pop. IF a population has 
        no samples meeting the minmap requirement in the first round of 
        imputation, then a second round is applied in which they sample 
        a genotype based on the overall (non-IMAP) genotype frequencies.
        """
        if 1:
            raise NotImplementedError()

        # override imap
        if not imap:
            imap = self.imap

        # impute data by mean value in each population
        newdata = self.snps.copy()
        for pop, samps in imap.items():

            # sample pop data
            sidxs = sorted(self.names.index(i) for i in samps)
            data = newdata[sidxs, :].copy()

            # number of alleles at each site that are not 9
            nallels = np.sum(data != 9, axis=0) * 2

            # get prob derived at each site using tmp array w/ missing to zero 
            tmp = data.copy()
            tmp[tmp == 9] = 0
            fderived = tmp.sum(axis=0) / nallels

            # sampler
            sampled = np.random.binomial(n=2, p=fderived, size=data.shape)
            data[data == 9] = sampled[data == 9]
            newdata[sidxs, :] = data

        # get all imputed values
        imputed = newdata[np.where(self.snps == 9)]
        self._print(
            "Imputation: 'sampled'; (0, 1, 2) = {:.1f}%, {:.1f}%, {:.1f}%"
            .format(
                100 * np.sum(imputed == 0) / imputed.size,
                100 * np.sum(imputed == 1) / imputed.size,
                100 * np.sum(imputed == 2) / imputed.size,
            )
        )
        return newdata
