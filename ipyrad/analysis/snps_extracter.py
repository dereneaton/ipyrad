#!/usr/bin/env python

"""
Load SNPs array and SNPsMAP from ipyrad SNPs HDF5 database for use in 
created subsampled SNP data sets for ipyrad.analysis tools. 

    Examples:
        pca, structure, treemix, popgen
"""

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard lib
import h5py
import numpy as np
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import jsubsample_snps

"""
.snps should be returned as dtype=... int8, int64?
.check imap for the user, e.g., no dups.
"""


class SNPsExtracter(object):
    """
    Extract .snps and .snpsmap from snps.hdf5 after filtering for indels, 
    bi-allelic, mincov, and minmap.    
    """
    def __init__(
        self, 
        data, 
        imap=None,
        minmap=None,
        mincov=0.0,
        quiet=False,
        ):

        # store params
        self.data = data
        self.imap = imap
        self.imap = (imap if imap else {})
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})
        self.mincov = mincov
        self.quiet = quiet

        # get names from imap, else will be filled/checked against db
        self.names = []
        for key, val in self.imap.items():
            if isinstance(val, (list, tuple)):
                self.names.extend(val)
            elif isinstance(val, str):
                self.names.append(val)        

        # array of (ntaxa, nsnps) as int8, only samples in imap
        self.snps = None
        # array of (nsnps, 2) for subsampling SNPs with subsample_snps()
        self.snpsmap = None

        # check input data
        if self.data.endswith(".seqs.hdf5"):
            raise IPyradError("data should be .snps.hdf5 file not .seqs.hdf5.")
        if self.data.endswith(".vcf"):
            raise NotImplementedError("VCF parsing coming soon.")
        self.parse_names_from_hdf5()



    def _print(self, message):
        if not self.quiet:
            print(message)


    def parse_names_from_hdf5(self):
        """
        Parse ordered sample names from HDF5 and check agains imap names.
        """
        # load arrays from hdf5; could be more efficient...
        io5 = h5py.File(self.data, 'r')

        # all names in database, maybe fewer in this analysis
        self.dbnames = [i.decode() for i in io5["snps"].attrs["names"]]
        # check for sample names not in the database file
        badnames = set(self.names).difference(self.dbnames)
        if badnames:
            raise IPyradError(
                "Samples [{}] are not in data file: {}"
                .format(badnames, self.data))

        # if no names then use database names, else order imap names in dborder
        if not self.names:
            self.names = self.dbnames 
        else:
            self.names = [i for i in self.dbnames if i in self.names]

        # mask snp rows to only include selected samples
        self.sidxs = np.array(sorted([self.dbnames.index(i) for i in self.names]))


    def parse_genos_from_hdf5(self):
        """
        Parse genotype calls from hdf5 snps file and store snpsmap
        for subsampling.
        """
        # load arrays from hdf5; could be more efficient...
        io5 = h5py.File(self.data, 'r')

        # snps is used to filter multi-allel and indel containing 
        snps = io5["snps"][:]

        # snpsmap is used to subsample per locus
        snpsmap = io5["snpsmap"][:, :2]

        # genos are the actual calls we want, after filtering
        genos = io5["genos"][:]  # .sum(axis=2).T

        # report pre-filter
        self._print("Samples: {}".format(len(self.names)))
        self._print("Sites before filtering: {}".format(snps.shape[1]))

        # filter all sites containing an indel in selected samples
        mask0 = np.any(snps[self.sidxs, :] == 45, axis=0)
        self._print("Filtered (indels): {}".format(mask0.sum()))

        # filter all sites w/ multi-allelic
        mask1 = np.sum(genos[:, self.sidxs] == 2, axis=2).sum(axis=1).astype(bool)
        mask1 += np.sum(genos[:, self.sidxs] == 3, axis=2).sum(axis=1).astype(bool)
        self._print("Filtered (bi-allel): {}".format(mask1.sum()))

        # convert genos (e.g., 0/1) to genos sums (e.g., 2)
        genos = genos.sum(axis=2).T
        # convert any summed missing (18) to 9
        genos[genos == 18] = 9

        # filter based on total sample min coverage
        cov = np.sum(genos[self.sidxs, :] != 9, axis=0)
        if isinstance(self.mincov, int):
            mask2 = cov < self.mincov
        elif isinstance(self.mincov, float):
            mask2 = cov < self.mincov * len(self.sidxs)
        else:
            raise IPyradError("mincov should be an int or float.")
        self._print("Filtered (mincov): {}".format(mask2.sum()))

        # filter based on per-population min coverages
        if not self.imap:
            mask3 = np.zeros(snps.shape[1], dtype=np.bool_)
        else:
            mask3 = np.zeros(snps.shape[1], dtype=np.bool_)
            for key, val in self.imap.items():
                mincov = self.minmap[key]
                pidxs = np.array(sorted(self.dbnames.index(i) for i in val))
                try:
                    subarr = genos[pidxs, :]
                except IndexError:
                    raise IPyradError("imap is empty: {} - {}".format(key, val))
                counts = np.sum(subarr != 9, axis=0)
                if isinstance(mincov, float):
                    mask3 += (counts / subarr.shape[0]) < mincov
                elif isinstance(mincov, int):
                    mask3 += counts < mincov
                else:
                    raise IPyradError("minmap dictionary malformed.")
        self._print("Filtered (minmap): {}".format(mask3.sum()))

        # apply mask to snps
        summask = mask0 + mask1 + mask2 + mask3
        allmask = np.invert(summask)
        self._print("Filtered (combined): {}".format(summask.sum()))
        self.snps = genos[self.sidxs, :][:, allmask]

        # bail out if ALL snps were filtered
        if self.snps.size == 0:
            raise IPyradError("No SNPs passed filtering.")

        # report state
        self._print(
            "Sites after filtering: {}".format(self.snps.shape[1])
        )
        self._mvals = np.sum(self.snps == 9)
        self._msites = np.any(self.snps == 9, axis=0).sum()
        self._print(
            "Sites containing missing values: {} ({:.2f}%)"
            .format(
                self._msites, 
                100 * self._msites / self.snps.shape[1],
            )
        )
        self._print(
            "Missing values in SNP matrix: {} ({:.2f}%)"
            .format(
                self._mvals, 
                100 * self._mvals / self.snps.size,                
            )
        )

        # apply mask to snpsmap to map snps
        snpsmap = snpsmap[allmask, :]
        snpsmap[:, 1] = range(snpsmap.shape[0])
        self.snpsmap = snpsmap
        del snpsmap

        # close it up
        io5.close()

        # .snps is an array of 0,1,2 or 9.
        # .snpsmap is ready to subsample .snps to 1-per-locus 


    def subsample_snps(self, random_seed=None, quiet=False):
        """
        Calls jitted functions to subsample 1 SNP per locus/linkage-block
        using snpsmap.
        """
        if not random_seed:
            random_seed = np.random.randint(0, 1e9)
        subarr = self.snps[:, jsubsample_snps(self.snpsmap, random_seed)]
        if not quiet:
            self._print("subsampled {} unlinked SNPs".format(subarr.shape[1]))
        return subarr


    def subsample_loci(self, random_seed=None, quiet=False):
        """
        Calls jitted functions to subsample loci/linkage-blocks with
        replacement to the same number as the original assembly. This does
        not subsample unlinked SNPs per locus.
        """
        raise NotImplementedError("TODO")
        if not random_seed:
            random_seed = np.random.randint(0, 1e9)
        subarr = self.snps[:, jsubsample_snps(self.snpsmap, random_seed)]
        if not quiet:
            self._print(
                "subsampled {} loci w/ replacement"
                .format(subarr.shape[1]))
        return subarr
