#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142
# pylint: disable=R0915
# pylint: disable=R0914
# pylint: disable=R0912

from __future__ import print_function, division

## ipyrad tools
from ipyrad.assemble.util import IPyradWarningExit, IPyradError, progressbar
from ipyrad import Assembly
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import itertools
import copy
import os

try:
    ## when you have time go back and set attrubutes on toytrees
    import allel
except ImportError:
    raise IPyradWarningExit("""
    Error: pca requires the dependency 'scikit-allel', which we haven't yet
    included in the ipyrad installation. For now, you can install scikit-allel
    using conda with the following command: 

    conda install scikit-allel -c conda-forge
    """)


## set floating point precision in data frames to 3 for prettier printing
pd.set_option('precision', 3)


class PCA(object):
    "new pca class object"
    def __init__(self, 
        data=None, 
        pops=None,
        ncomps=10,
        quiet=True):
        """ 
        ipyrad.analysis Baba Class object.

        Parameters
        ----------
        data : Assembly object or path to file
            Either an ipyrad assembly or a  string path to a .vcf file. If
            it's a string path then you'll probably want to specify pops as
            well or else all your dots will be the same color.
            
        pops : dict or path to file
            A dictionary specifying the population assignment of each
            sample. This is optional, since by default if you used a pops
            file during your assembly the assembly object will include
            the pops info internally.
        ncomps : int
            The number of PCs to calculate. Probably most people won't care
            to mess with this, but it's simple enough to make it flexible. 

        Functions
        ---------
        run()
            ...
        plot()
            ...

        """
        self.quiet = quiet
        self.ncomponents = ncomps

        ## parse data as (1) path to data file, or (2) ndarray
        if isinstance(data, Assembly):
            self.assembly = data
            self.pops = data.populations
            try:
                self.data = data.outfiles.vcf
            except AttributeError as inst:
                raise IPyradError(MISSING_VCF_ERROR)  
        else:
            ## You need a dummy assembly because we use the machinery
            ## of _link_populations below to read in the pops data
            self.assembly = Assembly("ipyrad-pca-tmp", quiet=True)
            self.data = os.path.realpath(data)
            self.pops = {}

        if pops:
            if isinstance(pops, dict):
                ## This is kind of stupid since we're just going to undo this
                ## in like 5 lines, but it gets the passed in pops into the
                ## same format as an assembly.populations dict, just easier to
                ## treat everything the same way.
                self.pops = {x:(0, y) for x, y in pops.items()}
            else:
                if not os.path.isfile(pops):
                    raise IPyradError("popfile does not exist - {}".format(pops))

                ## If the file you pass in doesn't have the stupid ipyrad minsamp
                mindat = [i.lstrip("#").lstrip().rstrip() for i in \
                          open(pops, 'r').readlines() if i.startswith("#")]
                if not mindat:
                    lines = open(pops, 'r').readlines()
                    p = set([x.split()[1].strip() for x in lines])
                    with open(pops, 'a') as outfile:
                        outfile.write("# " + " ".join(["{}:1".format(x) for x in p]))

                self.assembly.paramsdict["pop_assign_file"] = os.path.realpath(pops)
                self.assembly._link_populations()
                self.pops = self.assembly.populations
            
        ## Here the populations continues to maintain info about minsamps,
        ## which we just get rid of for clarity. Sorry this is dumb, I couldn't
        ## figure out a clean way to extract from a tuple inside the dict values.
        tmpdict = {}
        for samp in self.pops:
            tmpdict[samp] = self.pops[samp][1]
        self.pops = tmpdict

        ## Read in the vcf and extract the samples and the data
        ## This will set self.samples_vcforder which is a list of sample names
        ## in the order they appear in the vcf file
        self._load_calldata()

        ## If no pops linked yet (either none in the assembly or none passed in)
        ## then everybody goes into one giant default population.
        if not self.pops:
            self.pops = {"All_samples":self.samples_vcforder}

        if not self.quiet:
            print("  Using populations:\n{}".format(self.pops))
        if not self.pops:
            print("  No populations assigned, so PCA will be monochrome.")


    ## Load in the vcf and automatically remove multi-allelic snps
    ## and biallelic singletons.
    def _load_calldata(self):
        callset = allel.read_vcf(self.data, fields=["samples", "GT"])
        self.samples_vcforder = callset["samples"]

        gt = allel.GenotypeArray(callset['calldata/GT'])

        ## All this is for removing multi-allelic snps, and biallelic singletons
        ac = gt.count_alleles()
        flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)

        self.genotypes = gt.compress(flt, axis=0)


    def get_missing_per_sample(self):
        return pd.Series(self.genotypes.count_missing(axis=0),\
                        index=self.samples_vcforder)


    def get_nsnp_no_missing(self):
        return np.count_nonzero(pca.genotypes.count_missing(axis=1) == 0)


    ## Count number of snps that would be retained if n missing sites
    ## are allowed from 0 to maxm. Returns a series of these values.
    def _count_miss(self, name=None, maxm=10):
        tot_snp = len(self.genotypes)
        missdict = OrderedDict({})
        for n in range(0, maxm):
            misscounts = self.genotypes.count_missing(axis=1)
            misses = misscounts > n
            nmiss = self.genotypes[:][misses]
            missdict[n] = tot_snp - nmiss.shape[0]
        return pd.Series(missdict.values(), index=missdict.keys(), name=name)


    ## This returns a dataframe that shows the number of snps you'll retain
    ## after trimming for certain levels of missing data and/or removing
    ## the samples with the most missingness. It helps figure out which
    ## samples to drop and what level of missingness to accept.
    def missingness(self, maxm=10, nsamps=5, keep_dupes=False):
        trash_pca = self.copy()

        ## How much missingness without removing any samples
        miss_df = pd.DataFrame(trash_pca._count_miss(name="Full", maxm=maxm))

        miss_samps = trash_pca.get_missing_per_sample()
        ## Get the first n samples with the most missing data
        miss_samps = miss_samps.sort_values()[:nsamps]
        for k, v in miss_samps.iteritems():
            #print("removing {}".format(k))
            trash_pca.remove_samples(k)
            #print(self.genotypes.shape)
            ret = trash_pca._count_miss(name=k, maxm=maxm)
            miss_df = pd.concat([miss_df, ret], axis=1)

        if not keep_dupes:
            miss_df = miss_df.drop_duplicates()

        return miss_df


    ## remove all snps with more than max_missing samples with missing data
    ## This is a destructive operation.
    def trim_missing(self, max_missing=100):
        misscounts = self.genotypes.count_missing(axis=1)
        misses = misscounts > max_missing
        self.genotypes = self.genotypes[:][~misses]


    ## Attempt to fill missing values with the most frequent genotype per
    ## population. In this case if ./. is the most frequent, then we just
    ## keep it.
    ## TODO: Could have other modes of filling missing here, for example
    ## could make it so ./. doesn't 'count' as a genotype.
    def fill_missing(self, value='mode'):
        ## Local function to get the most frequent genotype. If there's
        ## a tie it will just take the first one
        def mode(ndarray):
            modals, counts = np.unique(ndarray, axis=0, return_counts=True)
            index = np.argmax(counts)
            return modals[index], counts[index]

        for pop, samps in self.pops.items():
            mask = np.isin(self.samples_vcforder, samps)
            for idx, row in enumerate(self.genotypes):
                if self.genotypes[idx][mask].count_missing():
                    val = mode(row[mask])
                    #print(gt[idx].to_str(threshold=12))
                    midx = self.genotypes[idx][mask].is_missing()
                    tmp = self.genotypes[idx][mask]
                    tmp[midx] = val[0]
                    self.genotypes[idx][mask] = tmp


    def remove_samples(self, samps):
        ## Allow to just pass in one sample as a string
        if isinstance(samps, str):
            samps = [samps]

        if set(samps) > set(self.samples_vcforder):
            raise IPyradError("  Trying to remove samples not present in the vcf file: {}".format(samps))

        ## Remove the samples from the sample list
        mask = np.isin(self.samples_vcforder, samps)
        self.samples_vcforder = self.samples_vcforder[~mask]

        self.genotypes = self.genotypes[:, ~mask]
        ## Remove biallelic singletons. If you don't do this you get
        ## a nasty error during svd, like this:
        ## https://stackoverflow.com/questions/33447808/sklearns-plsregression-valueerror-array-must-not-contain-infs-or-nans
        ac = self.genotypes.count_alleles()
        flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
        self.genotypes = self.genotypes.compress(flt, axis=0)

        if len(self.samples_vcforder) < self.ncomponents:
            self.ncomponents = len(self.samples_vcforder)
            print("  INFO: Number of PCs may not exceed the number of samples.\n  Setting number of PCs = {}".format(self.ncomponents))

        
    def plot(self, pcs=[1, 2], ax=None, cmap=None, cdict=None, legend=True, title=None, outfile=None):
        """
        Do the PCA and plot it.

        Parameters
        ---------
        pcs: list of ints
        ...
        ax: matplotlib axis
        ...
        cmap: matplotlib colormap
        ...
        cdict: dictionary mapping pop names to colors
        ...
        legend: boolean, whether or not to show the legend

        """
        ## Specify which 2 pcs to plot, default is pc1 and pc2
        pc1 = pcs[0] - 1
        pc2 = pcs[1] - 1
        if pc1 < 0 or pc2 > self.ncomponents - 1:
            raise IPyradError("PCs are 1-indexed. 1 is min & {} is max".format(self.ncomponents))

        ## Convert genotype data to allele count data
        ## We do this here because we might want to try different ways
        ## of accounting for missing data and "alt" allele counts treat
        ## missing data as "ref"
        allele_counts = self.genotypes.to_n_alt()

        ## Actually do the pca
        if self.ncomponents > len(self.samples_vcforder):
            self.ncomponents = len(self.samples_vcforder)
            print("  INFO: # PCs < # samples. Forcing # PCs = {}".format(self.ncomponents))
        coords, model = allel.pca(allele_counts, n_components=self.ncomponents, scaler='patterson')

        self.pcs = pd.DataFrame(coords,
                                index=self.samples_vcforder,
                                columns=["PC{}".format(x) for x in range(1,self.ncomponents+1)])

        ## Just allow folks to pass in the name of the cmap they want to use
        if isinstance(cmap, str):
            try:
                cmap = cm.get_cmap(cmap)
            except:
                raise IPyradError("  Bad cmap value: {}".format(cmap))


        if not cmap and not cdict:
            if not self.quiet:
                print("  Using default cmap: Spectral")
            cmap = cm.get_cmap('Spectral')

        if cmap:
            if cdict:
                print("  Passing in both cmap and cdict defaults to using the cmap value.")
            popcolors = cmap(np.arange(len(self.pops))/len(self.pops))
            cdict = {i:j for i, j in zip(self.pops.keys(), popcolors)}

        fig = ""
        if not ax:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(1, 1, 1)

        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in self.pops:
            ## Don't include pops with no samples, it makes the legend look stupid
            ## TODO: This doesn't prevent empty pops from showing up in the legend for some reason.
            if len(self.pops[pop]) > 0:
                mask = np.isin(self.samples_vcforder, self.pops[pop])
                ax.plot(x[mask], y[mask], marker='o', linestyle=' ', color=cdict[pop], label=pop, markersize=6, mec='k', mew=.5)

        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

        if legend:
            ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

        if fig:
            fig.tight_layout()

        if title:
            ax.set_title(title)

        if outfile:
            try:
                plt.savefig(outfile, format="png", bbox_inches="tight")
            except:
                print("  Saving pca.plot() failed to save figure to {}".format(outfile))

        return ax


    def plot_pairwise_dist(self, labels=None, ax=None, cmap=None, cdict=None, metric="euclidean"):
        """
        Plot pairwise distances between all samples

        labels: bool or list
                by default labels aren't included. If labels == True, then labels are read in
                from the vcf file. Alternatively, labels can be passed in as a list, should
                be same length as the number of samples.
        """
        allele_counts = self.genotypes.to_n_alt()
        dist = allel.pairwise_distance(allele_counts, metric=metric)
        if not ax:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_subplot(1, 1, 1)

        if isinstance(labels, bool):
            if labels:
                labels = list(self.samples_vcforder)
        elif isinstance(labels, type(None)):
            pass
        else:
            ## If not bool or None (default), then check to make sure the list passed in
            ## is the right length
            if not len(labels) == len(self.samples_vcforder):
                raise IPyradError(LABELS_LENGTH_ERROR.format(len(labels), len(self.samples_vcforder)))

        allel.plot.pairwise_distance(dist, labels=labels, ax=ax, colorbar=False)


    def copy(self):
        """ returns a copy of the pca analysis object """
        cp = copy.deepcopy(self)
        cp.genotypes = allel.GenotypeArray(self.genotypes, copy=True)
        return cp


MISSING_VCF_ERROR = "  Assembly does not contain a vcf file. Rerun step 7 with `v` included in the `output_formats` parameter."
LABELS_LENGTH_ERROR ="Labels must be same length as number of samples. You have {} labels and {} samples."


#######################################################################
if __name__ == "__main__":
    print("Nothing implemented here.")
