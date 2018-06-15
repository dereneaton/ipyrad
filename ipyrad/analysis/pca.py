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

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import itertools
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
    
        Functions
        ---------
        run()
            ...
        plot()
            ...

        """
        self.quiet = quiet

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
                self.pops = pops
            else:
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


    def remove_samples(self, samps):

        ## Allow to just pass in one sample as a string
        if isinstance(samps, str):
            samps = [samps]

        if set(samps) > set(self.samples_vcforder):
            raise IPyradError("  Trying to remove samples not present in the vcf file: {}".format(samps))

        mask = np.isin(self.samples_vcforder, samps)
        self.samples_vcforder = self.samples_vcforder[~mask]
        self.genotypes = self.genotypes.T[~mask]
        self.genotypes = self.genotypes.T
        

    def plot(self, pcs=[1, 2], ax=None, cmap=None, cdict=None):
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

        """
        ## Specify which 2 pcs to plot, default is pc1 and pc2
        pc1 = pcs[0] - 1
        pc2 = pcs[1] - 1
        if pc1 < 0 or pc2 > 9:
            raise IPyradError("PCs are 1-indexed. 1 is min & 10 is max")

        ## Convert genotype data to allele count data
        ## We do this here because we might want to try different ways
        ## of accounting for missing data and "alt" allele counts treat
        ## missing data as "ref"
        allele_counts = self.genotypes.to_n_alt()

        ## Actually do the pca
        coords, model = allel.stats.pca(allele_counts, n_components=10, scaler='patterson')

        ## Just allow folks to pass in the name of the cmap they want to use
        if isinstance(cmap, str):
            try:
                cmap = cm.get_cmap(cmap)
            except:
                raise IPyradError("  Bad cmap value: {}".format(cmap))


        if not cmap and not cdict:
            print("  Using default cmap: Spectral")
            cmap = cm.get_cmap('Spectral')

        if cmap:
            if cdict:
                print("  Passing in both cmap and cdict defaults to using the cmap value.")
            popcolors = cmap(np.arange(len(self.pops))/len(self.pops))
            cdict = {i:j for i, j in zip(self.pops.keys(), popcolors)}

        if not ax:
            fig = plt.figure(figsize=(6, 5))
            ax = fig.add_subplot(1, 1, 1)

        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in self.pops:
            mask = np.isin(self.samples_vcforder, self.pops[pop])
            ax.plot(x[mask], y[mask], marker='o', linestyle=' ', color=cdict[pop], label=pop, markersize=6, mec='k', mew=.5)

        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

        ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

        if fig:
            fig.tight_layout()


    def copy(self):
        """ returns a copy of the baba analysis object """
        return copy.deepcopy(self)


    MISSING_VCF_ERROR = "  Assembly does not contain a vcf file. Rerun step 7 with `v` included in the `output_formats` parameter."


#######################################################################
if __name__ == "__main__":
    print("Nothing implemented here.")
