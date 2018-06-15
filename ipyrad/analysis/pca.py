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


    def _load_calldata(self):
        callset = allel.read_vcf(self.data, fields=["samples", "GT"])
        self.samples_vcforder = callset["samples"]

        gt = allel.GenotypeArray(callset['calldata/GT'])

        ## All this is for removing multi-allelic snps, and biallelic singletons
        ac = gt.count_alleles()
        flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
        gt_filtered = gt.compress(flt, axis=0)
        ## Get an array of counts of the number of alt alleles
        ## Only keep this around bcz it's all we need for the pca
        self.allele_counts = gt_filtered.to_n_alt()


    def remove_samples(self, samps):

        ## Allow to just pass in one sample as a string
        if isinstance(samps, string):
            samps = [samps]

        if set(samps) > set(self.samples_vcforder):
            raise IPyradError("  Trying to remove samples not present in the vcf file: {}".format(samps))

        mask = np.isin(self.samples_vcforder, samps)
        self.samples_vcforder = self.samples_vcforder[~mask]
        self.allele_counts = self.allele_counts[~mask]
        

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

        ## Actually do the pca
        coords, model = allel.stats.pca(self.allele_counts, n_components=10, scaler='patterson')

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


    def plotwwat(self, 
        show_test_labels=True, 
        use_edge_lengths=True,         
        collapse_outgroup=False, 
        pct_tree_x=0.5, 
        pct_tree_y=0.2,
        subset_tests=None,
        #toytree_kwargs=None,
        *args, 
        **kwargs):

        """ 
        Draw a multi-panel figure with tree, tests, and results 
        
        Parameters:
        -----------
        height: int
        ...

        width: int
        ...

        show_test_labels: bool
        ...

        use_edge_lengths: bool
        ...

        collapse_outgroups: bool
        ...

        pct_tree_x: float
        ...

        pct_tree_y: float
        ...

        subset_tests: list
        ...

        ...

        """

        ## check for attributes
        if not self.newick:
            raise IPyradError("baba plot requires a newick treefile")
        if not self.tests:
            raise IPyradError("baba plot must have a .tests attribute")

        ## ensure tests is a list
        if isinstance(self.tests, dict):
            self.tests = [self.tests]

        ## re-decompose the tree
        ttree = toytree.tree(
            self.newick, 
            orient='down', 
            use_edge_lengths=use_edge_lengths,
            )

        ## subset test to show fewer
        if subset_tests != None:
            #tests = self.tests[subset_tests]
            tests = [self.tests[i] for i in subset_tests]
            boots = self.results_boots[subset_tests]
        else:
            tests = self.tests
            boots = self.results_boots

        ## make the plot
        canvas, axes, panel = baba_panel_plot(
            ttree=ttree,
            tests=tests, 
            boots=boots, 
            show_test_labels=show_test_labels, 
            use_edge_lengths=use_edge_lengths, 
            collapse_outgroup=collapse_outgroup, 
            pct_tree_x=pct_tree_x,
            pct_tree_y=pct_tree_y,
            *args, 
            **kwargs)
        return canvas, axes, panel


    def copy(self):
        """ returns a copy of the baba analysis object """
        return copy.deepcopy(self)


    MISSING_VCF_ERROR = "  Assembly does not contain a vcf file. Rerun step 7 with `v` included in the `output_formats` parameter."


def batch(
    baba,
    ipyclient=None,
    ):
    """
    distributes jobs to the parallel client
    """

    ## parse args
    handle = baba.data
    taxdicts = baba.tests
    mindicts = baba.params.mincov
    nboots = baba.params.nboots

    ## if ms generator make into reusable list
    sims = 0
    if isinstance(handle, types.GeneratorType):
        handle = list(handle)
        sims = 1
    else:
        ## expand locifile path to full path
        handle = os.path.realpath(handle)

    ## parse taxdicts into names and lists if it a dictionary
    #if isinstance(taxdicts, dict):
    #    names, taxdicts = taxdicts.keys(), taxdicts.values()
    #else:
    #    names = []
    names = []
    if isinstance(taxdicts, dict):
        taxdicts = [taxdicts]

    ## an array to hold results (len(taxdicts), nboots)
    tot = len(taxdicts)
    resarr = np.zeros((tot, 7), dtype=np.float64)
    bootsarr = np.zeros((tot, nboots), dtype=np.float64)
    paneldict = {}

    ## TODO: Setup a wrapper to find and cleanup ipyclient
    ## define the function and parallelization to use, 
    ## if no ipyclient then drops back to using multiprocessing.
    if not ipyclient:
        # ipyclient = ip.core.parallel.get_client(**self._ipcluster)
        raise IPyradError("you must enter an ipyparallel.Client() object")
    else:
        lbview = ipyclient.load_balanced_view()

    ## submit jobs to run on the cluster queue
    start = time.time()
    asyncs = {}
    idx = 0

    ## prepare data before sending to engines
    ## if it's a str (locifile) then parse it here just once.
    if isinstance(handle, str):
        with open(handle, 'r') as infile:
            loci = infile.read().strip().split("|\n")
    if isinstance(handle, list):
        pass #sims()

    ## iterate over tests (repeats mindicts if fewer than taxdicts)
    itests = iter(taxdicts)
    imdict = itertools.cycle([mindicts])

    #for test, mindict in zip(taxdicts, itertools.cycle([mindicts])):
    for i in xrange(len(ipyclient)):

        ## next entries unless fewer than len ipyclient, skip
        try:
            test = next(itests)
            mindict = next(imdict)
        except StopIteration:
            continue

        ## if it's sim data then convert to an array
        if sims:
            loci = _msp_to_arr(handle, test)
            args = (loci, test, mindict, nboots)
            print("not yet implemented")
            #asyncs[idx] = lbview.apply_async(dstat, *args)
        else:
            args = [loci, test, mindict, nboots]
            asyncs[idx] = lbview.apply(dstat, *args)
        idx += 1

    ## block until finished, print progress if requested.
    finished = 0
    try:
        while 1:
            keys = [i for (i, j) in asyncs.items() if j.ready()]
            ## check for failures
            for job in keys:
                if not asyncs[job].successful():
                    raise IPyradWarningExit(\
                        " error: {}: {}".format(job, asyncs[job].exception()))
                ## enter results for successful jobs
                else:
                    _res, _bot = asyncs[job].result()
                    
                    ## store D4 results
                    if _res.shape[0] == 1:
                        resarr[job] = _res.T.as_matrix()[:, 0]
                        bootsarr[job] = _bot
                    
                    ## or store D5 results                        
                    else:   
                        paneldict[job] = _res.T

                    ## remove old job
                    del asyncs[job]
                    finished += 1

                    ## submit next job if there is one.
                    try:
                        test = next(itests)
                        mindict = next(imdict)
                        if sims:
                            loci = _msp_to_arr(handle, test)
                            args = (loci, test, mindict, nboots)
                            print("not yet implemented")
                            #asyncs[idx] = lbview.apply_async(dstat, *args)
                        else:
                            args = [loci, test, mindict, nboots]
                            asyncs[idx] = lbview.apply(dstat, *args)
                        idx += 1
                    except StopIteration:
                        pass

            ## count finished and break if all are done.
            #fin = idx - len(asyncs)
            elap = datetime.timedelta(seconds=int(time.time()-start))
            printstr = " calculating D-stats  | {} | "
            progressbar(tot, finished, printstr.format(elap), spacer="")
            time.sleep(0.1)
            if not asyncs:
                print("")
                break

    except KeyboardInterrupt as inst:
        ## cancel all jobs (ipy & multiproc modes) and then raise error
        try:
            ipyclient.abort()
        except Exception:
            pass
        raise inst

    ## dress up resarr as a Pandas DataFrame if 4-part test
    if len(test) == 4:
        if not names:
            names = range(len(taxdicts))
        #print("resarr")
        #print(resarr)
        resarr = pd.DataFrame(resarr, 
            index=names,
            columns=["dstat", "bootmean", "bootstd", "Z", "ABBA", "BABA", "nloci"])

        ## sort results and bootsarr to match if test names were supplied
        resarr = resarr.sort_index()
        order = [list(resarr.index).index(i) for i in names]
        bootsarr = bootsarr[order]
        return resarr, bootsarr
    else:
        ## order results dfs
        listres = []
        for key in range(len(paneldict)):
            listres.append(paneldict[key])
            
        ## make into a multi-index dataframe
        ntests = len(paneldict)
        multi_index = [
            np.array([[i] * 3 for i in range(ntests)]).flatten(),
            np.array(['p3', 'p4', 'shared'] * ntests),
        ]
        resarr = pd.DataFrame(
            data=pd.concat(listres).as_matrix(), 
            index=multi_index,
            columns=listres[0].columns,
            )
        return resarr, None
        #return listres, None  #_res.T, _bot



def dstat(inarr, taxdict, mindict=1, nboots=1000, name=0):
    """ private function to perform a single D-stat test"""

    #if isinstance(inarr, str):
    #    with open(inarr, 'r') as infile:
    #        inarr = infile.read().strip().split("|\n")

    # ## get data as an array from loci file
    # ## if loci-list then parse arr from loci
    if isinstance(inarr, list):
        arr, _ = _loci_to_arr(inarr, taxdict, mindict)
    
    # ## if it's an array already then go ahead
    # elif isinstance(inarr, np.ndarray):
    #     arr = inarr
    # ## if it's a simulation object get freqs from array
    # elif isinstance(inarr, Sim):
    #     arr = _msp_to_arr(inarr, taxdict)

    #elif isinstance(inarr, types.GeneratorType):
    #    arr = _msp_to_arr(inarr, taxdict)
    #elif isinstance(inarr, list):
    #    arr = _msp_to_arr(inarr, taxdict)
    ## get data from Sim object, do not digest the ms generator
    #else:
    #    raise Exception("Must enter either a 'locifile' or 'arr'")

    ## run tests
    #if len(taxdict) == 4:
    if arr.shape[1] == 4:

        ## get results
        res, boots = _get_signif_4(arr, nboots)
    
        ## make res into a nice DataFrame
        res = pd.DataFrame(res, 
            columns=[name],
            index=["Dstat", "bootmean", "bootstd", "Z", "ABBA", "BABA", "nloci"])

    else:
        ## get results
        res, boots = _get_signif_5(arr, nboots)
         ## make int a DataFrame
        res = pd.DataFrame(res,
            index=["p3", "p4", "shared"], 
            columns=["Dstat", "bootmean", "bootstd", "Z", "ABxxA", "BAxxA", "nloci"]
            )

    return res.T, boots



def _loci_to_arr(loci, taxdict, mindict):
    """
    return a frequency array from a loci file for all loci with taxa from 
    taxdict and min coverage from mindict. 
    """

    ## make the array (4 or 5) and a mask array to remove loci without cov
    nloci = len(loci)
    keep = np.zeros(nloci, dtype=np.bool_)
    arr = np.zeros((nloci, 4, 300), dtype=np.float64)

    ## six rows b/c one for each p3, and for the fused p3 ancestor
    if len(taxdict) == 5:
        arr = np.zeros((nloci, 6, 300), dtype=np.float64)

    ## if not mindict, make one that requires 1 in each taxon
    if isinstance(mindict, int):
        mindict = {i: mindict for i in taxdict}
    elif isinstance(mindict, dict):
        mindict = {i: mindict[i] for i in taxdict}
    else:
        mindict = {i: 1 for i in taxdict}

    ## raise error if names are not 'p[int]' 
    allowed_names = ['p1', 'p2', 'p3', 'p4', 'p5']
    if any([i not in allowed_names for i in taxdict]):
        raise IPyradError(\
            "keys in taxdict must be named 'p1' through 'p4' or 'p5'")

    ## parse key names
    keys = sorted([i for i in taxdict.keys() if i[0] == 'p'])
    outg = keys[-1]

    ## grab seqs just for the good guys
    for loc in xrange(nloci):

        ## parse the locus
        lines = loci[loc].split("\n")[:-1]
        names = [i.split()[0] for i in lines]
        seqs = np.array([list(i.split()[1]) for i in lines])

        ## check that names cover the taxdict (still need to check by site)
        covs = [sum([j in names for j in taxdict[tax]]) >= mindict[tax] \
                for tax in taxdict]

        ## keep locus
        if all(covs):
            keep[loc] = True

            ## get the refseq
            refidx = np.where([i in taxdict[outg] for i in names])[0]
            refseq = seqs[refidx].view(np.uint8)
            ancestral = np.array([reftrick(refseq, GETCONS2)[:, 0]])

            ## freq of ref in outgroup
            iseq = _reffreq2(ancestral, refseq, GETCONS2)
            arr[loc, -1, :iseq.shape[1]] = iseq 

            ## enter 4-taxon freqs
            if len(taxdict) == 4:
                for tidx, key in enumerate(keys[:-1]):

                    ## get idx of names in test tax
                    nidx = np.where([i in taxdict[key] for i in names])[0]
                    sidx = seqs[nidx].view(np.uint8)
                   
                    ## get freq of sidx
                    iseq = _reffreq2(ancestral, sidx, GETCONS2)
                   
                    ## fill it in 
                    arr[loc, tidx, :iseq.shape[1]] = iseq

            else:

                ## entere p5; and fill it in
                iseq = _reffreq2(ancestral, refseq, GETCONS2) 
                arr[loc, -1, :iseq.shape[1]] = iseq 
                
                ## enter p1
                nidx = np.where([i in taxdict['p1'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 0, :iseq.shape[1]] = iseq
                
                ## enter p2
                nidx = np.where([i in taxdict['p2'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 1, :iseq.shape[1]] = iseq
                
                ## enter p3 with p4 masked, and p4 with p3 masked
                nidx = np.where([i in taxdict['p3'] for i in names])[0]
                nidy = np.where([i in taxdict['p4'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                sidy = seqs[nidy].view(np.uint8)
                xseq = _reffreq2(ancestral, sidx, GETCONS2)
                yseq = _reffreq2(ancestral, sidy, GETCONS2)
                mask3 = xseq != 0
                mask4 = yseq != 0
                xseq[mask4] = 0
                yseq[mask3] = 0
                arr[loc, 2, :xseq.shape[1]] = xseq
                arr[loc, 3, :yseq.shape[1]] = yseq
                
                ## enter p34 
                nidx = nidx.tolist() + nidy.tolist()
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 4, :iseq.shape[1]] = iseq


    ## size-down array to the number of loci that have taxa for the test
    arr = arr[keep, :, :]

    ## size-down sites to 
    arr = masknulls(arr)

    return arr, keep



## This should be re-written as a dynamic func
def tree2tests(newick, constraint_dict, constraint_exact):
    """
    Returns dict of all possible four-taxon splits in a tree. Assumes
    the user has entered a rooted tree. Skips polytomies.
    """
    ## make tree
    tree = toytree.tree(newick)
    testset = set()

    ## expand constraint_exact if list
    if isinstance(constraint_exact, bool):
        constraint_exact = [constraint_exact] * 4
    elif isinstance(constraint_exact, list):
        if len(constraint_exact) != len(constraint_dict):
            raise Exception("constraint_exact must be bool or [bool, bool, bool, bool]")
    
    ## constraints
    cdict = {"p1":[], "p2":[], "p3":[], "p4":[]}
    if constraint_dict:
        cdict.update(constraint_dict)

    ## traverse root to tips. Treat the left as outgroup, then the right.
    tests = []
    
    ## topnode must have children
    for topnode in tree.tree.traverse("levelorder"):
        for oparent in topnode.children:
            for onode in oparent.traverse("levelorder"):
                if test_constraint(onode, cdict, "p4", constraint_exact[3]):
                    #print(topnode.name, onode.name)
                    
                    ## p123 parent is sister to oparent
                    p123parent = oparent.get_sisters()[0]
                    for p123node in p123parent.traverse("levelorder"):
                        for p3parent in p123node.children:
                            for p3node in p3parent.traverse("levelorder"):
                                if test_constraint(p3node, cdict, "p3", constraint_exact[2]):
                                    #print(topnode.name, onode.name, p3node.name)
                                    
                                    ## p12 parent is sister to p3parent
                                    p12parent = p3parent.get_sisters()[0]
                                    for p12node in p12parent.traverse("levelorder"):
                                        for p2parent in p12node.children:
                                            for p2node in p2parent.traverse("levelorder"):
                                                if test_constraint(p2node, cdict, "p2", constraint_exact[1]):

                                                    ## p12 parent is sister to p3parent
                                                    p1parent = p2parent.get_sisters()[0]
                                                    for p1node in p1parent.traverse("levelorder"):
                                                        #for p1parent in p1node.children:
                                                        #    for p1node in p1parent.traverse("levelorder"):
                                                        if test_constraint(p1node, cdict, "p1", constraint_exact[0]):
                                                            x = (onode.name, p3node.name, p2node.name, p1node.name)
                                                            test = {}
                                                            test['p4'] = onode.get_leaf_names()
                                                            test['p3'] = p3node.get_leaf_names()
                                                            test['p2'] = p2node.get_leaf_names()
                                                            test['p1'] = p1node.get_leaf_names()
                                                            if x not in testset:
                                                                tests.append(test)
                                                                testset.add(x)
    return tests                                            



def test_constraint(node, cdict, tip, exact):
    names = set(node.get_leaf_names())
    const = set(cdict[tip])
    if const:
        if exact:
            #if len(names.intersection(const)) == len(const):
            if names == const:
                return 1
            else:
                return 0
        else:
            if len(names.intersection(const)) == len(names):
                return 1
            else:
                return 0        
    return 1
    
    


#######################################################################
if __name__ == "__main__":
    print("Nothing implemented here.")
