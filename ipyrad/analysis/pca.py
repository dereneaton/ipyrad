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


    def _parse_genos_from_hdf5(self):
        """
        Parse genotype calls from hdf5 snps file and store snpsmap
        for subsampling.
        """
        # load arrays from hdf5; could be more efficient...
        io5 = h5py.File(self.data, 'r')

        # snps is used to filter multi-allel and indel containing 
        snps = io5["snps"][:]

        # all names in database, maybe fewer in this analysis
        dbnames = [i.decode() for i in io5["snps"].attrs["names"]]

        # snpsmap is used to subsample per locus
        snpsmap = io5["snpsmap"][:, :2]

        # genos are the actual calls we want, after filtering
        genos = io5["genos"][:]  # .sum(axis=2).T

        # check for sample names not in the database file
        badnames = set(self.names).difference(dbnames)
        if badnames:
            raise IPyradError(
                "Samples [{}] are not in data file: {}"
                .format(badnames, self.data))
        
        # if no names then use database names, else order imap names in dborder
        if not self.names:
            self.names = dbnames 
        else:
            self.names = [i for i in dbnames if i in self.names]

        # report pre-filter
        self._print("Samples: {}".format(len(self.names)))
        self._print("Sites before filtering: {}".format(snps.shape[1]))

        # mask snp rows to only include selected samples
        sidxs = np.array(sorted([dbnames.index(i) for i in self.names]))

        # filter all sites containing an indel in selected samples
        mask0 = np.any(snps[sidxs, :] == 45, axis=0)
        self._print("Filtered (indels): {}".format(mask0.sum()))

        # filter all sites w/ multi-allelic
        mask1 = np.sum(genos[:, sidxs] == 2, axis=2).sum(axis=1).astype(bool)
        mask1 += np.sum(genos[:, sidxs] == 3, axis=2).sum(axis=1).astype(bool)
        self._print("Filtered (bi-allel): {}".format(mask1.sum()))

        # convert genos (e.g., 0/1) to genos sums (e.g., 2)
        genos = genos.sum(axis=2).T
        # convert any summed missing (18) to 9
        genos[genos == 18] = 9

        # filter based on total sample min coverage
        cov = np.sum(genos[sidxs, :] != 9, axis=0)
        if isinstance(self.mincov, int):
            mask2 = cov < self.mincov
        elif isinstance(self.mincov, float):
            mask2 = cov < self.mincov * len(sidxs)
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
                pidxs = np.array(sorted(dbnames.index(i) for i in val))
                subarr = genos[pidxs, :]
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
        self.snps = genos[sidxs, :][:, allmask]

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


    # to match hdf5 function this needs to parse a genotype matrix (.snps)
    # and build a .snpsmap from CHROM, POS.
    def _parse_genos_from_vcf(self):
        """
        Parse genotype calls from input VCF data file. This will be 
        0/1/2 or 9 for missing. In other words, it's the geno file 
        format. 
        """

        # for now, trying to avoid scikit-allel usage.

        # get column names from header line
        with open(self.data) as indata:
            while 1:
                data = indata.readline().strip().split()
                if data[0].upper() == "#CHROM":
                    columns = [data[0][1:]] + data[1:]
                    break

        # parse as a dataframe
        df = pd.read_csv(
            self.data, 
            comment="#",
            sep="\t",
            names=columns,
        )

        # replace INFO column with depth
        df["DEPTH"] 

        # drop unneeded columns
        df = df.drop(columns=["CHROM", "REF", "ALT", "QUAL"])


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



# class PCA(object):
#     "new pca class object"
#     def __init__(self, 
#         data=None, 
#         pops=None,
#         ncomps=10,
#         quiet=True):
#         """ 
#         ipyrad.analysis Baba Class object.

#         Parameters
#         ----------
#         data : Assembly object or path to file
#             Either an ipyrad assembly or a  string path to a .vcf file. If
#             it's a string path then you'll probably want to specify pops as
#             well or else all your dots will be the same color.
            
#         pops : dict or path to file
#             A dictionary specifying the population assignment of each
#             sample. This is optional, since by default if you used a pops
#             file during your assembly the assembly object will include
#             the pops info internally.
#         ncomps : int
#             The number of PCs to calculate. Probably most people won't care
#             to mess with this, but it's simple enough to make it flexible. 

#         Functions
#         ---------
#         run()
#             ...
#         plot()
#             ...

#         """
#         self.quiet = quiet
#         self.ncomponents = ncomps

#         ## parse data as (1) path to data file, or (2) ndarray
#         if isinstance(data, Assembly):
#             self.assembly = data
#             self.pops = data.populations
#             try:
#                 self.data = data.outfiles.vcf
#             except AttributeError as inst:
#                 raise IPyradError(MISSING_VCF_ERROR)  
#         else:
#             ## You need a dummy assembly because we use the machinery
#             ## of _link_populations below to read in the pops data
#             self.assembly = Assembly("ipyrad-pca-tmp", quiet=True)
#             self.data = os.path.realpath(data)
#             self.pops = {}

#         if pops:
#             if isinstance(pops, dict):
#                 ## This is kind of stupid since we're just going to undo this
#                 ## in like 5 lines, but it gets the passed in pops into the
#                 ## same format as an assembly.populations dict, just easier to
#                 ## treat everything the same way.
#                 self.pops = {x: (0, y) for x, y in pops.items()}
#             else:
#                 if not os.path.isfile(pops):
#                     raise IPyradError(
#                         "popfile does not exist - {}".format(pops))

#                 ## If the file you pass in doesn't have the ipyrad minsamp
#                 mindat = [
#                     i.lstrip("#").lstrip().rstrip() for i in
#                     open(pops, 'r').readlines() if i.startswith("#")
#                 ]
#                 if not mindat:
#                     lines = open(pops, 'r').readlines()
#                     p = set([x.split()[1].strip() for x in lines])
#                     with open(pops, 'a') as outfile:
#                         outfile.write(
#                             "# " + " ".join(["{}:1".format(x) for x in p]))

#                 self.assembly.paramsdict["pop_assign_file"] = os.path.realpath(pops)
#                 self.assembly._link_populations()
#                 self.pops = self.assembly.populations
            
#         ## Here the populations continues to maintain info about minsamps,
#         ## which we just get rid of for clarity. Sorry this is dumb, I couldn't
#         ## figure out a clean way to extract from a tuple inside the dict values.
#         tmpdict = {}
#         for samp in self.pops:
#             tmpdict[samp] = self.pops[samp][1]
#         self.pops = tmpdict

#         ## Read in the vcf and extract the samples and the data
#         ## This will set self.samples_vcforder which is a list of sample names
#         ## in the order they appear in the vcf file
#         self._load_calldata()

#         ## If no pops linked yet (either none in the assembly or none passed in)
#         ## then everybody goes into one giant default population.
#         if not self.pops:
#             self.pops = {"All_samples": self.samples_vcforder}

#         if not self.quiet:
#             print("  Using populations:\n{}".format(self.pops))
#         if not self.pops:
#             print("  No populations assigned, so PCA will be monochrome.")

#     ## Load in the vcf and automatically remove multi-allelic snps and biallelic singletons.
#     def _load_calldata(self):
#         callset = allel.read_vcf(self.data, fields=["samples", "GT"])
#         self.samples_vcforder = callset["samples"]

#         gt = allel.GenotypeArray(callset['calldata/GT'])

#         ## All this is for removing multi-allelic snps, and biallelic singletons
#         ac = gt.count_alleles()
#         flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)

#         self.genotypes = gt.compress(flt, axis=0)


#     def get_missing_per_sample(self):
#         return pd.Series(self.genotypes.count_missing(axis=0),\
#                         index=self.samples_vcforder)


#     def get_nsnp_no_missing(self):
#         return np.count_nonzero(pca.genotypes.count_missing(axis=1) == 0)


#     ## Count number of snps that would be retained if n missing sites
#     ## are allowed from 0 to maxm. Returns a series of these values.
#     def _count_miss(self, name=None, maxm=10):
#         tot_snp = len(self.genotypes)
#         missdict = OrderedDict({})
#         for n in range(0, maxm):
#             misscounts = self.genotypes.count_missing(axis=1)
#             misses = misscounts > n
#             nmiss = self.genotypes[:][misses]
#             missdict[n] = tot_snp - nmiss.shape[0]
#         ret = pd.Series(
#             missdict.values(), 
#             index=list(missdict.keys()), 
#             name=name,
#         )
#         return ret


#     ## This returns a dataframe that shows the number of snps you'll retain
#     ## after trimming for certain levels of missing data and/or removing
#     ## the samples with the most missingness. It helps figure out which
#     ## samples to drop and what level of missingness to accept.
#     def missingness(self, maxm=10, nsamps=5, keep_dupes=False):
#         trash_pca = self.copy()

#         ## How much missingness without removing any samples
#         miss_df = pd.DataFrame(trash_pca._count_miss(name="Full", maxm=maxm))

#         miss_samps = trash_pca.get_missing_per_sample()
#         ## Get the first n samples with the most missing data
#         miss_samps = miss_samps.sort_values()[:nsamps]
#         for k, v in miss_samps.iteritems():
#             #print("removing {}".format(k))
#             trash_pca.remove_samples(k)
#             #print(self.genotypes.shape)
#             ret = trash_pca._count_miss(name=k, maxm=maxm)
#             miss_df = pd.concat([miss_df, ret], axis=1)

#         if not keep_dupes:
#             miss_df = miss_df.drop_duplicates()

#         return miss_df


#     ## remove all snps with more than max_missing samples with missing data
#     ## This is a destructive operation.
#     def trim_missing(self, max_missing=100):
#         misscounts = self.genotypes.count_missing(axis=1)
#         misses = misscounts > max_missing
#         self.genotypes = self.genotypes[:][~misses]


#     ## Attempt to fill missing values with the most frequent genotype per
#     ## population. In this case if ./. is the most frequent, then we just
#     ## keep it.
#     ## TODO: Could have other modes of filling missing here, for example
#     ## could make it so ./. doesn't 'count' as a genotype.
#     def fill_missing(self, value='mode'):
#         ## Local function to get the most frequent genotype. If there's
#         ## a tie it will just take the first one
#         def mode(ndarray):
#             modals, counts = np.unique(ndarray, axis=0, return_counts=True)
#             index = np.argmax(counts)
#             return modals[index], counts[index]

#         for pop, samps in self.pops.items():
#             mask = np.isin(self.samples_vcforder, samps)
#             for idx, row in enumerate(self.genotypes):
#                 if self.genotypes[idx][mask].count_missing():
#                     val = mode(row[mask])
#                     #print(gt[idx].to_str(threshold=12))
#                     midx = self.genotypes[idx][mask].is_missing()
#                     tmp = self.genotypes[idx][mask]
#                     tmp[midx] = val[0]
#                     self.genotypes[idx][mask] = tmp


#     def remove_samples(self, samps):
#         ## Allow to just pass in one sample as a string
#         if isinstance(samps, str):
#             samps = [samps]

#         if set(samps) > set(self.samples_vcforder):
#             raise IPyradError("  Trying to remove samples not present in the vcf file: {}".format(samps))

#         ## Remove the samples from the sample list
#         mask = np.isin(self.samples_vcforder, samps)
#         self.samples_vcforder = self.samples_vcforder[~mask]

#         self.genotypes = self.genotypes[:, ~mask]
#         ## Remove biallelic singletons. If you don't do this you get
#         ## a nasty error during svd, like this:
#         ## https://stackoverflow.com/questions/33447808/sklearns-plsregression-valueerror-array-must-not-contain-infs-or-nans
#         ac = self.genotypes.count_alleles()
#         flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
#         self.genotypes = self.genotypes.compress(flt, axis=0)

#         if len(self.samples_vcforder) < self.ncomponents:
#             self.ncomponents = len(self.samples_vcforder)
#             print("  INFO: Number of PCs may not exceed the number of samples.\n  Setting number of PCs = {}".format(self.ncomponents))

        
#     def plot(self, pcs=[1, 2], ax=None, cmap=None, cdict=None, legend=True, title=None, outfile=None):
#         """
#         Do the PCA and plot it.

#         Parameters
#         ---------
#         pcs: list of ints
#         ...
#         ax: matplotlib axis
#         ...
#         cmap: matplotlib colormap
#         ...
#         cdict: dictionary mapping pop names to colors
#         ...
#         legend: boolean, whether or not to show the legend

#         """
#         ## Specify which 2 pcs to plot, default is pc1 and pc2
#         pc1 = pcs[0] - 1
#         pc2 = pcs[1] - 1
#         if pc1 < 0 or pc2 > self.ncomponents - 1:
#             raise IPyradError("PCs are 1-indexed. 1 is min & {} is max".format(self.ncomponents))

#         ## Convert genotype data to allele count data
#         ## We do this here because we might want to try different ways
#         ## of accounting for missing data and "alt" allele counts treat
#         ## missing data as "ref"
#         allele_counts = self.genotypes.to_n_alt()

#         ## Actually do the pca
#         if self.ncomponents > len(self.samples_vcforder):
#             self.ncomponents = len(self.samples_vcforder)
#             print("  INFO: # PCs < # samples. Forcing # PCs = {}".format(self.ncomponents))
#         coords, model = allel.stats.pca(allele_counts, n_components=self.ncomponents, scaler='patterson')

#         self.pcs = pd.DataFrame(coords,
#                                 index=self.samples_vcforder,
#                                 columns=["PC{}".format(x) for x in range(1,self.ncomponents+1)])

#         ## Just allow folks to pass in the name of the cmap they want to use
#         if isinstance(cmap, str):
#             try:
#                 cmap = cm.get_cmap(cmap)
#             except:
#                 raise IPyradError("  Bad cmap value: {}".format(cmap))


#         if not cmap and not cdict:
#             if not self.quiet:
#                 print("  Using default cmap: Spectral")
#             cmap = cm.get_cmap('Spectral')

#         if cmap:
#             if cdict:
#                 print("  Passing in both cmap and cdict defaults to using the cmap value.")
#             popcolors = cmap(np.arange(len(self.pops))/len(self.pops))
#             cdict = {i:j for i, j in zip(self.pops.keys(), popcolors)}

#         fig = ""
#         if not ax:
#             fig = plt.figure(figsize=(6, 5))
#             ax = fig.add_subplot(1, 1, 1)

#         x = coords[:, pc1]
#         y = coords[:, pc2]
#         for pop in self.pops:
#             ## Don't include pops with no samples, it makes the legend look stupid
#             ## TODO: This doesn't prevent empty pops from showing up in the legend for some reason.
#             if len(self.pops[pop]) > 0:
#                 mask = np.isin(self.samples_vcforder, self.pops[pop])
#                 ax.plot(x[mask], y[mask], marker='o', linestyle=' ', color=cdict[pop], label=pop, markersize=6, mec='k', mew=.5)

#         ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
#         ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

#         if legend:
#             ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

#         if fig:
#             fig.tight_layout()

#         if title:
#             ax.set_title(title)

#         if outfile:
#             try:
#                 plt.savefig(outfile, format="png", bbox_inches="tight")
#             except:
#                 print("  Saving pca.plot() failed to save figure to {}".format(outfile))

#         return ax


#     def plot_pairwise_dist(self, labels=None, ax=None, cmap=None, cdict=None, metric="euclidean"):
#         """
#         Plot pairwise distances between all samples

#         labels: bool or list
#                 by default labels aren't included. If labels == True, then labels are read in
#                 from the vcf file. Alternatively, labels can be passed in as a list, should
#                 be same length as the number of samples.
#         """
#         allele_counts = self.genotypes.to_n_alt()
#         dist = allel.stats.pairwise_distance(allele_counts, metric=metric)
#         if not ax:
#             fig = plt.figure(figsize=(5, 5))
#             ax = fig.add_subplot(1, 1, 1)

#         if isinstance(labels, bool):
#             if labels:
#                 labels = list(self.samples_vcforder)
#         elif isinstance(labels, type(None)):
#             pass
#         else:
#             ## If not bool or None (default), then check to make sure the list passed in
#             ## is the right length
#             if not len(labels) == len(self.samples_vcforder):
#                 raise IPyradError(LABELS_LENGTH_ERROR.format(len(labels), len(self.samples_vcforder)))

#         allel.plot.pairwise_distance(dist, labels=labels, ax=ax, colorbar=False)


#     def copy(self):
#         """ returns a copy of the pca analysis object """
#         cp = copy.deepcopy(self)
#         cp.genotypes = allel.GenotypeArray(self.genotypes, copy=True)
#         return cp


MISSING_VCF_ERROR = "  Assembly does not contain a vcf file. Rerun step 7 with `v` included in the `output_formats` parameter."
LABELS_LENGTH_ERROR ="Labels must be same length as number of samples. You have {} labels and {} samples."


#######################################################################
if __name__ == "__main__":
    print("Nothing implemented here.")
