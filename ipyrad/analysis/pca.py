#!/usr/bin/env python

""" Scikit-learn principal componenents analysis for missing data """

from __future__ import print_function, division

import os
import sys
import numpy as np
import pandas as pd

# ipyrad tools
from .snps_extracter import SNPsExtracter
from .snps_imputer import SNPsImputer
from ipyrad.analysis.utils import jsubsample_snps
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
    imap: (dict; default=None)
        Dictionary mapping population names to a list of sample names.
    minmap: (dict; default={})
        Dictionary mapping population names to float values (X).
        If a site does not have data across X proportion of samples for
        each population, respectively, the site is filtered from the data set.
    mincov: (float; default=0.5)
        If a site does not have data across this proportion of total samples
        in the data then it is filtered from the data set.
    impute_method: (str; default='sample')
        None, "sample", or an integer for the number of kmeans clusters.
    topcov: (float; default=0.9)
        Affects kmeans method only.    
        The most stringent mincov used as the first iteration in kmeans 
        clustering. Subsequent iterations (niters) are equally spaced between
        topcov and mincov. 
    niters: (int; default=5)
        Affects kmeans method only.        
        kmeans method only.
        Number of iterations of kmeans clustering with decreasing mincov 
        thresholds used to refine population clustering, and therefore to 
        refine the imap groupings used to filter and impute sites.

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
        quiet=False,
        topcov=0.9,
        niters=5,
        #ncomponents=None,
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
        #self.ncomponents = ncomponents
        self.impute_method = impute_method
        self.mincov = mincov        
        self.imap = (imap if imap else {})
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})
        self.topcov = topcov
        self.niters = niters

        # get names from imap else we'll fill this later with all
        # self.names = []
        # for key, val in self.imap.items():
        #     if isinstance(val, (list, tuple)):
        #         self.names.extend(val)
        #     elif isinstance(val, str):
        #         self.names.append(val)

        # to be filled
        self.snps = np.array([])
        self.snpsmap = np.array([])
        self.nmissing = 0

        # coming soon...
        if self.data.endswith(".vcf"):
            raise NotImplementedError(
                "Sorry, not yet supported. Use .snps.hdf5.")

        # load .snps and .snpsmap from HDF5
        first = (True if isinstance(self.impute_method, int) else quiet)
        ext = SNPsExtracter(
            self.data, self.imap, self.minmap, self.mincov, quiet=first,
        )

        # run snp extracter to parse data files
        ext.parse_genos_from_hdf5()       
        self.snps = ext.snps
        self.snpsmap = ext.snpsmap
        self.names = ext.names
        self._mvals = ext._mvals

        # make imap for imputing if not used in filtering.
        if not self.imap:
            self.imap = {'1': self.names}
            self.minmap = {'1': 0.5}
            
        # record missing data per sample
        self.missing = pd.DataFrame({
            "missing": [0.],
            },
            index=self.names,
        )
        miss = np.sum(self.snps == 9, axis=1) / self.snps.shape[1]
        for name in self.names:
            self.missing.missing[name] = round(miss[self.names.index(name)], 2)

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
        # if self.impute_method == "simple":
        # self.snps = SNPsImputer(
        # self.snps, self.names, self.imap, None).run()

        if self.impute_method == "sample":
            self.snps = SNPsImputer(
                self.snps, self.names, self.imap, "sample", self.quiet).run()

        elif isinstance(self.impute_method, int):
            self.snps = self._impute_kmeans(
                self.topcov, self.niters, self.quiet)

        else:
            self.snps[self.snps == 9] = 0
            self._print(
                "Imputation (null; sets to 0): {:.1f}%, {:.1f}%, {:.1f}%"
                .format(100, 0, 0)            
            )


    def _impute_kmeans(self, topcov=0.9, niters=5, quiet=False):

        # the ML models to fit
        pca_model = decomposition.PCA(n_components=None)  # self.ncomponents)
        kmeans_model = KMeans(n_clusters=self.impute_method)

        # start kmeans with a global imap
        kmeans_imap = {'global': self.names}

        # iterate over step values
        iters = np.linspace(topcov, self.mincov, niters)
        for it, kmeans_mincov in enumerate(iters):

            # start message
            kmeans_minmap = {i: self.mincov for i in kmeans_imap}
            self._print(
                "Kmeans clustering: iter={}, K={}, mincov={}, minmap={}"
                .format(it, self.impute_method, kmeans_mincov, kmeans_minmap))

            # 1. Load orig data and filter with imap, minmap, mincov=step
            se = SNPsExtracter(
                self.data, 
                imap=kmeans_imap, 
                minmap=kmeans_minmap, 
                mincov=kmeans_mincov,
                quiet=self.quiet,
            )
            se.parse_genos_from_hdf5()

            # update snpsmap to new filtered data to use for subsampling            
            self.snpsmap = se.snpsmap

            # 2. Impute missing data using current kmeans clusters
            impdata = SNPsImputer(
                se.snps, se.names, kmeans_imap, "sample", self.quiet).run()

            # x. On final iteration return this imputed array as the result
            if it == 4:
                return impdata

            # 3. subsample unlinked SNPs
            subdata = impdata[:, jsubsample_snps(se.snpsmap, self._seed())]

            # 4. PCA on new imputed data values
            pcadata = pca_model.fit_transform(subdata)

            # 5. Kmeans clustering to find new imap grouping
            kmeans_model.fit(pcadata)
            labels = np.unique(kmeans_model.labels_)           
            kmeans_imap = {
                i: [se.names[j] for j in 
                    np.where(kmeans_model.labels_ == i)[0]] for i in labels
            }
            self._print(kmeans_imap)
            self._print("")


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
            data = self.snps[:, jsubsample_snps(self.snpsmap, self.seed)]
            self._print(
                "Subsampling SNPs: {}/{}"
                .format(data.shape[1], self.snps.shape[1])
            )
        else:
            data = self.snps

        # decompose pca call
        model = decomposition.PCA(None)  # self.ncomponents)
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
