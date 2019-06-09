#!/usr/bin/env python

""" Scikit-learn principal componenents analysis for missing data """

from __future__ import print_function, division

import os
import sys
import h5py
import numpy as np
import pandas as pd
from numba import njit

# ipyrad tools
from .snps_extracter import SNPsExtracter
from ipyrad.assemble.utils import IPyradError

# missing imports to be raised on class init
try:
    import toyplot
except ImportError:
    pass

_MISSING_TOYPLOT = ImportError("""
This ipyrad tool requires the plotting library toyplot. 
You can install it with the following command in a terminal.

conda install toyplot -c eaton-lab 
""")

try:
    from sklearn import decomposition 
    from sklearn.impute import SimpleImputer
    from sklearn.cluster import KMeans
except ImportError:
    pass

_MISSING_SKLEARN = """
This ipyrad tool requires the library scikit-learn.
You can install it with the following command in a terminal.

conda install scikit-learn -c conda-forge 
"""


# TODO: could allow LDA as alternative to PCA for supervised (labels) dsets.
# TODO: remove biallel singletons... (option, not sure it's a great idea...)

class PCA(object):
    """
    Principal components analysis of RAD-seq SNPs with iterative
    imputation of missing data.
    
    Parameters:
    -----------
    data: (str, several options)
        A general .vcf file or a .snps.hdf5 file produced by ipyrad.
    workdir: (str; default="./analysis-pca")
        A directory for output files. Will be created if absent.
    ncomponents: (int; default=2)
        Number of PCA components
    kmeans: (int; default=None)
        Number of expected clusters, optionally used for data imputation.
    imap: (dict; default=None)
        Dictionary to assign samples to groups to be used with minmap.
    minmap: (dict; default={})
        Dictionary to assign minimum coverage for groups to filter SNPs
        for inclusion in the analysis. This filter applies before 
        (optional) data imputation step.
    mincov: (float; default=0.5)
        Proportion of total samples that are not N at any site to include
        in data set. 
    impute_method: (str; default='sample')
        None, "sample", "simple", "kmeans"

    Functions:
    ----------
    ...
    """
    def __init__(
        self, 
        data, 
        impute_method=None,
        imap=None,
        minmap=None,
        mincov=0.1,
        ncomponents=None,
        quiet=False,
        ):

        # only check import at init
        if not sys.modules.get("sklearn"):
            raise IPyradError(_MISSING_SKLEARN)
        if not sys.modules.get("toyplot"):
            raise IPyradError(_MISSING_TOYPLOT)

        # init attributes
        self.quiet = quiet
        self.data = os.path.realpath(os.path.expanduser(data))

        # data attributes
        self.ncomponents = ncomponents
        self.impute_method = impute_method
        self.mincov = mincov        
        self.imap = (imap if imap else {})
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})

        # get names from imap else we'll fill this later with all
        self.names = []
        for key, val in self.imap.items():
            if isinstance(val, (list, tuple)):
                self.names.extend(val)
            elif isinstance(val, str):
                self.names.append(val)

        # to be filled
        self.snps = np.array([])
        self.snpsmap = np.array([])
        self.nmissing = 0

        # load .snps and .snpsmap from HDF5
        ext = SNPsExtracter(
            self.data, self.imap, self.minmap, self.mincov, quiet=quiet,
        )
        ext.parse_genos_from_hdf5()
        self.snps = ext.snps
        self.snpsmap = ext.snpsmap
        self.names = ext.names
        self._mvals = ext._mvals

        if self.data.endswith(".vcf"):
            raise NotImplementedError(
                "Sorry, not yet supported. Use .snps.hdf5.")

        # make imap for imputing if not used in filtering.
        if not self.imap:
            self.imap = {'1': self.names}
            self.minmap = {'1': 2}

        # impute missing data
        if (self.impute_method is not False) and self._mvals:
            self._impute_data()


    def _seed(self):   
        return np.random.randint(0, 1e9)        


    def _print(self, msg):
        if not self.quiet:
            print(msg)


    def _impute_data(self):
        """
        Impute data in-place updating self.snps by filling missing (9) values.
        """
        # simple imputer method
        if self.impute_method == "simple":
            self.snps = self._impute_simple()

        elif self.impute_method == "sample":
            self.snps = self._impute_sample()

        elif isinstance(self.impute_method, int):
            self.snps = self._impute_kmeans()

        else:
            self.snps[self.snps == 9] = 0
            self._print(
                "Imputation (null; sets to 0): {:.1f}%, {:.1f}%, {:.1f}%"
                .format(self._mvals, 0, 0)            
            )


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
            "Imputation (sampled by freq. within pops): {:.1f}%, {:.1f}%, {:.1f}%"
            .format(
                100 * np.sum(imputed == 0) / imputed.size,
                100 * np.sum(imputed == 1) / imputed.size,
                100 * np.sum(imputed == 2) / imputed.size,
            )
        )
        return newdata


    def _impute_simple(self, quiet=False, imap=None):
        """
        Assign most frequent value to fill 9 in each column for each pop.
        """
        # override imap
        if not imap:
            imap = self.imap

        # impute data by mean value in each population
        newdata = self.snps.copy()
        model = SimpleImputer(missing_values=9, strategy="most_frequent")
        for pop, samps in imap.items():
            sidxs = sorted(self.names.index(i) for i in samps)
            data = self.snps[sidxs, :]
            model.fit(data)
            newdata[sidxs] = model.transform(data).astype(int)

        # get all imputed values
        imputed = newdata[np.where(self.snps == 9)]
        self._print(
            "Imputation (simple; most freq. within pops): {:.1f}%, {:.1f}%, {:.1f}%"
            .format(
                100 * np.sum(imputed == 0) / imputed.size,
                100 * np.sum(imputed == 1) / imputed.size,
                100 * np.sum(imputed == 2) / imputed.size,
            )
        )
        return newdata


    def _impute_kmeans(self, quiet=False):
        """
        kmeans imputer method takes the following steps:
        
        1. A temporary data set is created by imputing missing data with 
           the most frequent base across all data. 
        2. Clusters are formed using kmeans clustering from unlinked SNPs 
           randomly sampled from this tmp data set. 
        3. Missing data in the original data set are imputed for samples in 
           each clustered group by taking the most frequent base in that 
           cluster. 
        4. If all data were missing for a site in a given cluster then the
           imputed value is the most frequent across all samples. 
        """
        # bail out if no imap dictionary
        if not self.impute_method:
            raise IPyradError(
                "Kmeans imputation method requires 'kmeans' value.")

        # ML models
        fill_model = SimpleImputer(missing_values=9, strategy="most_frequent")
        kmeans_model = KMeans(n_clusters=self.impute_method)
        pca_model = decomposition.PCA(n_components=self.ncomponents)

        # first impute into a copy with most-frequent value across all 
        fill_model.fit(self.snps)
        tmp_data = fill_model.transform(self.snps).astype(int)
        gfreq = fill_model.statistics_

        # cluster data based on subsample of unlinked SNPs from tmp_data
        data = tmp_data[:, subsample_snps(self.snpsmap, self._seed())]
        newdata = pca_model.fit_transform(data)
        kmeans_model.fit(newdata)

        # create a new tmp_imap from kmeans clusters
        labels = np.unique(kmeans_model.labels_)
        tmp_imap = {
            i: np.where(kmeans_model.labels_ == i)[0] for i in labels
        }

        # impute data by mean value in each population. 
        newdata = self.snps.copy()
        for pop, sidxs in tmp_imap.items():
            data = self.snps[sidxs, :]

            # find and fill values that are ALL missing with global most freqx
            allmiss = np.all(data == 9, axis=0)
            data[:, allmiss] = gfreq[allmiss]

            # then use Simple tranformer on rest
            fill_model.fit(data)
            newdata[sidxs] = fill_model.transform(data).astype(int)

        # report
        if not quiet:
            # get all imputed values
            imputed = newdata[np.where(self.snps == 9)]
            self._print(
                "Imputation (kmeans; most freq in clusters): {:.1f}%, {:.1f}%, {:.1f}%"
                .format(
                    100 * np.sum(imputed == 0) / imputed.size,
                    100 * np.sum(imputed == 1) / imputed.size,
                    100 * np.sum(imputed == 2) / imputed.size,
                )
            )
        return newdata


    def run(self, seed=None, subsample=True):
        """
        Decompose genotype array (.snps) into n_components axes. 

        Parameters:
        -----------
        seed: (int)
            Random number seed used if/when subsampling SNPs.
        subsample: (bool)
            Subsample one SNP per RAD locus to reduce effect of linkage.

        Returns:
        --------
        A tuple with two numpy arrays. The first is the new data decomposed
        into principal coordinate space; the second is an array with the 
        variance explained by each PC axis. 
        """
        # update seed. Numba seed cannot be None, so get random int if None
        if seed:
            self.seed = seed
        else:
            self.seed = self._seed()

        # sample one SNP per locus
        if subsample:
            data = self.snps[:, subsample_snps(self.snpsmap, self.seed)]
        else:
            data = self.snps
        self._print(
            "Subsampling SNPs: {}/{}"
            .format(data.shape[1], self.snps.shape[1])
        )

        # decompose pca call
        model = decomposition.PCA(self.ncomponents)
        model.fit(data)
        newdata = model.transform(data)

        # return tuple with new coordinates and variance explained
        return newdata, model.explained_variance_ratio_                


    def run_and_plot_2D(
        self, 
        ax0=0, 
        ax1=1, 
        seed=None, 
        subsample=True, 
        ):
        """
        A convenience function for plotting 2D scatterplot of PCA results.
        """
        # get transformed coordinates and variances
        data, vexp = self.run(seed, subsample)  # , impute_method, kmeans)

        # color map
        colors = toyplot.color.broadcast(
            toyplot.color.brewer.map("Spectral"), shape=len(self.imap),
        )

        # make reverse imap dictionary
        irev = {}
        for pop, vals in self.imap.items():
            for val in vals:
                irev[val] = pop

        # assign styles to populations
        pstyles = {}
        for idx, pop in enumerate(self.imap):
            pstyles[pop] = toyplot.marker.create(
                size=10, 
                shape="o",
                mstyle={
                    "fill": toyplot.color.to_css(colors[idx]),
                    "stroke": "#262626",
                    "fill-opacity": 0.6,
                },
            )

        # assign styled markers to data points
        marks = []
        for name in self.names:
            pop = irev[name]
            mark = pstyles[pop]
            marks.append(mark)

        # plot points with colors x population
        canvas = toyplot.Canvas(400, 300)
        axes = canvas.cartesian(
            grid=(1, 5, 0, 1, 0, 4),
            xlabel="PC{} ({:.1f}%) explained".format(ax0, vexp[ax0] * 100),
            ylabel="PC{} ({:.1f}%) explained".format(ax1, vexp[ax1] * 100),
        )
        mark = axes.scatterplot(
            data[:, ax0],
            data[:, ax1],
            marker=marks,
            title=self.names,
        )

        # add a legend
        if len(self.imap) > 1:
            marks = [(pop, marker) for pop, marker in pstyles.items()]
            canvas.legend(
                marks, 
                corner=("right", 35, 100, min(250, len(pstyles) * 25))
            )

        return canvas, axes, mark



@njit
def subsample_snps(snpsmap, seed):
    "Subsample snps, one per locus, using snpsmap"
    np.random.seed(seed)
    sidxs = np.unique(snpsmap[:, 0])
    subs = np.zeros(sidxs.size, dtype=np.int64)
    idx = 0
    for sidx in sidxs:
        sites = snpsmap[snpsmap[:, 0] == sidx, 1]
        site = np.random.choice(sites)
        subs[idx] = site
        idx += 1
    return subs



if __name__ == "__main__":
    print("Nothing implemented here.")
