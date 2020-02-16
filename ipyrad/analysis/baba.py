#!/usr/bin/env python

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142
# pylint: disable=R0915
# pylint: disable=R0914
# pylint: disable=R0912

from __future__ import print_function, division

## ipyrad tools
from ipyrad.assemble.write_outputs import reftrick
from ipyrad.assemble.utils import IPyradError, GETCONS, Params
from ipyrad.analysis.utils import progressbar
#from ipyrad.plotting.baba_panel_plot import baba_panel_plot

#import scipy.stats as st  ## used for dfoil
import pandas as pd
import numpy as np
import numba
import itertools
import datetime
import types
import copy
import time
import os

# import tested at init
try: 
    import msprime as ms
except ImportError:
    pass
_MSPRIME_IMPORT= """
This ipyrad analysis tool requires 
You can install it with the following command:

   conda install toytree -c eaton-lab
"""

try:
    import toytree
except ImportError:
    pass
_TOYTREE_IMPORT = """
This ipyrad analysis tool requires 
You can install it with the following command:

   conda install toytree -c eaton-lab
"""
try:
    import toyplot
except ImportError:
    pass
_TOYPLOT_IMPORT = """
This ipyrad analysis tool requires the toyplot package.
You can install it with the following command:

   conda install toyplot -c eaton-lab
"""

## set floating point precision in data frames to 3 for prettier printing
pd.set_option('precision', 3)


class Baba(object):
    "new baba class object"
    def __init__(self, 
        data=None, 
        tests=None, 
        newick=None, 
        nboots=1000, 
        mincov=1):
        """ 
        ipyrad.analysis Baba Class object.

        Parameters
        ----------
        data : string or ndarray
            A string path to a .loci file produced by ipyrad. Alternatively, 
            data can be entered as a Numpy array of float allele frequencies 
            with dimension (nloci, 4 or 5, maxlen). See simulation example 
            in the docs.
            
        tests : dict or list of dicts
            A dictionary mapping Sample names to test taxon names, e.g., 
            test = {'p1': ['a', 'b'], 'p2': ['c'], 'p3': ['e'], 'p4': ['f']}.
            Four taxon tests should have p1-p4 whereas five taxon tests will 
            used if dict keys are p1-p5. Other key names will raise an error. 
            The highest value name (e.g., p5) is the outgroup. 
        
        newick: str
            ...
    
        Functions
        ---------
        run()
            ...
        generate_tests_from_tree()
            ...
        plot()
            ...

        """
        # check external imports
        if not sys.modules.get("toytree"):
            raise ImportError(_TOYTREE_IMPORT)
        if not sys.modules.get("toyplot"):
            raise ImportError(_TOYPLOT_IMPORT)
        if not sys.modules.get("msprime"):
            raise ImportError(_MSPRIME_IMPORT)

        ## parse data as (1) path to data file, or (2) ndarray
        if isinstance(data, str):
            self.data = os.path.realpath(data)
        else:
            self.data = data

        if isinstance(newick, toytree.tree):
            self.newick = newick.tree.write()
        else:
            self.newick = newick

        ## store tests, check for errors
        self.tests = tests

        ## parameters
        self.params = Params()
        self.params.mincov = mincov
        self.params.nboots = nboots
        self.params.quiet = False
        self.params.database = None

        ## results storage
        self.results_table = None
        self.results_boots = None
        


    @property
    def taxon_table(self):
        """
        Returns the .tests list of taxa as a pandas dataframe. 
        By auto-generating this table from tests it means that 
        the table itself cannot be modified unless it is returned 
        and saved. 
        """
        if self.tests:
            keys = sorted(self.tests[0].keys())
            if isinstance(self.tests, list):
                ld = [[(key, i[key]) for key in keys] for i in self.tests]
                dd = [dict(i) for i in ld]
                df = pd.DataFrame(dd)
                return df
            else:
                return pd.DataFrame(pd.Series(self.tests)).T
        else:
            return None



    def run(self, 
        ipyclient=None,
        ):
        """
        Run a batch of dstat tests on a list of tests, where each test is 
        a dictionary mapping sample names to {p1 - p4} (and sometimes p5). 
        Parameters modifying the behavior of the run, such as the number
        of bootstrap replicates (nboots) or the minimum coverage for 
        loci (mincov) can be set in {object}.params.

        Parameters:
        -----------
        ipyclient (ipyparallel.Client object):
            An ipyparallel client object to distribute jobs to a cluster. 
        """
        self.results_table, self.results_boots = batch(self, ipyclient)

        ## skip this for 5-part test results
        if not isinstance(self.results_table, list):
            self.results_table.nloci = np.nan_to_num(self.results_table.nloci)\
                                                 .astype(int)


    def generate_tests_from_tree(self, 
        constraint_dict=None, 
        constraint_exact=False, 
        verbose=True):
        """ 
        Returns a list of all possible 4-taxon tests on a tree (newick file). 
        The number of possible tests can be greatly reduced by setting 
        constraints on the taxon sampling using the constraint_dict arg. 

        Parameters:
        -----------
        constraint_dict: dict
            The constraint dict will limit the tests generated to only include
            the taxa listed in the dict. 

        constraint_exact: bool or list
            If constraint_exact is True then only samples meeting the exact 
            entries in the constraint_dict will be returned, as opposed to all
            subsets of those entries. If a list then different values can be 
            applied to [p1, p2, p3, p4]. For example, if the constraint_dict is
            {"p1": sample1, "p2": sample2, "p3": sample3, "p4": [sample4, sample5]},
            then with constraint_exact==False you get:
            
            sample1, sample2, sample3, sample4
            sample1, sample2, sample3, sample5
            sample1, sample2, sample3, [sample4, sample5]

            and with constraint_exact==True you get only:

            sample1, sample2, sample3, [sample4, sample5]
        """
        if not self.newick:
            raise AttributeError("no newick tree information in {self}.newick")
        tests = tree2tests(self.newick, constraint_dict, constraint_exact)
        if verbose:
            print("{} tests generated from tree".format(len(tests)))
        self.tests = tests


    def plot(self, 
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
        print("Plotting baba results is not implemented in v.0.9.")
        return

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
        #canvas, axes, panel = baba_panel_plot(
        #    ttree=ttree,
        #    tests=tests, 
        #    boots=boots, 
        #    show_test_labels=show_test_labels, 
        #    use_edge_lengths=use_edge_lengths, 
        #    collapse_outgroup=collapse_outgroup, 
        #    pct_tree_x=pct_tree_x,
        #    pct_tree_y=pct_tree_y,
        #    *args, 
        #    **kwargs)
        #return canvas, axes, panel


    def copy(self):
        """ returns a copy of the baba analysis object """
        return copy.deepcopy(self)



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
    for i in range(len(ipyclient)):

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
                    raise IPyradError(\
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
            printstr = " calculating D-stats"
            progressbar(finished, tot, start, message=printstr)
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
    maxlen = np.max(np.array([len(locus.split("\n")[0]) for locus in loci]))
    keep = np.zeros(nloci, dtype=np.bool_)
    arr = np.zeros((nloci, 4, maxlen), dtype=np.float64)

    ## six rows b/c one for each p3, and for the fused p3 ancestor
    if len(taxdict) == 5:
        arr = np.zeros((nloci, 6, maxlen), dtype=np.float64)

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
    for loc in range(nloci):

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
            refseq = seqs[refidx].view(np.uint32)
            ancestral = np.array([reftrick(refseq, GETCONS)[:, 0]])

            ## freq of ref in outgroup
            iseq = _reffreq2(ancestral, refseq, GETCONS)
            try:
                arr[loc, -1, :iseq.shape[1]] = iseq 
            except:
                raise IPyradError("Error in getting outgroup sequences.")
    
            ## enter 4-taxon freqs
            if len(taxdict) == 4:
                for tidx, key in enumerate(keys[:-1]):

                    ## get idx of names in test tax
                    nidx = np.where([i in taxdict[key] for i in names])[0]
                    sidx = seqs[nidx].view(np.uint32)
                   
                    ## get freq of sidx
                    iseq = _reffreq2(ancestral, sidx, GETCONS)
                   
                    ## fill it in 
                    arr[loc, tidx, :iseq.shape[1]] = iseq

            else:

                ## entere p5; and fill it in
                iseq = _reffreq2(ancestral, refseq, GETCONS) 
                arr[loc, -1, :iseq.shape[1]] = iseq 
                
                ## enter p1
                nidx = np.where([i in taxdict['p1'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint32)
                iseq = _reffreq2(ancestral, sidx, GETCONS)
                arr[loc, 0, :iseq.shape[1]] = iseq
                
                ## enter p2
                nidx = np.where([i in taxdict['p2'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint32)
                iseq = _reffreq2(ancestral, sidx, GETCONS)
                arr[loc, 1, :iseq.shape[1]] = iseq
                
                ## enter p3 with p4 masked, and p4 with p3 masked
                nidx = np.where([i in taxdict['p3'] for i in names])[0]
                nidy = np.where([i in taxdict['p4'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint32)
                sidy = seqs[nidy].view(np.uint32)
                xseq = _reffreq2(ancestral, sidx, GETCONS)
                yseq = _reffreq2(ancestral, sidy, GETCONS)
                mask3 = xseq != 0
                mask4 = yseq != 0
                xseq[mask4] = 0
                yseq[mask3] = 0
                arr[loc, 2, :xseq.shape[1]] = xseq
                arr[loc, 3, :yseq.shape[1]] = yseq
                
                ## enter p34 
                nidx = nidx.tolist() + nidy.tolist()
                sidx = seqs[nidx].view(np.uint32)
                iseq = _reffreq2(ancestral, sidx, GETCONS)
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
    for topnode in tree.treenode.traverse("levelorder"):
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
    
    


@numba.jit(nopython=True)
def masknulls(arr):
    nvarr = np.zeros(arr.shape[0], dtype=np.int8)
    trimarr = np.zeros(arr.shape, dtype=np.float64)
    for loc in range(arr.shape[0]):
        nvars = 0
        for site in range(arr.shape[2]):
            col = arr[loc, :, site]
            ## mask cols with 9s
            if not np.any(col == 9):
                ## any non-outgroup shows variation?
                ## todo: check whether BBBBA is ever info?
                if np.any(col[:-1] != col[0]):
                    trimarr[loc, :, nvars] = col
                    nvars += 1
        nvarr[loc] = nvars        
    return trimarr[:, :, :nvarr.max()]



@numba.jit(nopython=True)
def _reffreq2(ancestral, iseq, consdict):
    ## empty arrays
    freq = np.zeros((1, iseq.shape[1]), dtype=np.float64)
    amseq = np.zeros((iseq.shape[0]*2, iseq.shape[1]), dtype=np.uint8)
    
    ## fill in both copies
    for seq in range(iseq.shape[0]):
        for col in range(iseq.shape[1]):  

            ## get this base and check if it is hetero
            base = iseq[seq][col]
            who = consdict[:, 0] == base
            
            ## if not hetero then enter it
            if not np.any(who):
                amseq[seq*2][col] = base
                amseq[seq*2+1][col] = base        
            ## if hetero then enter the 2 resolutions
            else:
                amseq[seq*2][col] = consdict[who, 1][0]
                amseq[seq*2+1][col] = consdict[who, 2][0]

    ## amseq may have N or -, these need to be masked
    for i in range(amseq.shape[1]):
        ## without N or -
        reduced = amseq[:, i][amseq[:, i] != 9]
        counts = reduced != ancestral[0][i]
        if reduced.shape[0]:
            freq[:, i] = counts.sum() / reduced.shape[0]
        else:
            freq[:, i] = 9
    return freq



@numba.jit(nopython=True)
def _prop_dstat(arr):
    
    ## numerator
    abba = ((1.-arr[:, 0]) * (arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))  
    baba = ((arr[:, 0]) * (1.-arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))
    top = abba - baba
    bot = abba + baba

    ## get statistic and avoid zero div  
    sbot = bot.sum()
    if  sbot != 0:
        dst = top.sum() / float(sbot)
    else:
        dst = 0
    
    return abba.sum(), baba.sum(), dst



@numba.jit(nopython=True)
def _get_boots(arr, nboots):
    """
    return array of bootstrap D-stats
    """
    ## hold results (nboots, [dstat, ])
    boots = np.zeros((nboots,))
    
    ## iterate to fill boots
    for bidx in range(nboots):
        ## sample with replacement
        lidx = np.random.randint(0, arr.shape[0], arr.shape[0])
        tarr = arr[lidx]
        _, _, dst = _prop_dstat(tarr)
        boots[bidx] = dst
    
    ## return bootarr
    return boots



@numba.jit(nopython=True)
def _get_signif_4(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """
    abba, baba, dst = _prop_dstat(arr)
    boots = _get_boots(arr, nboots)
    estimate, stddev = (boots.mean(), boots.std())
    zscore = 0.
    if stddev:
        zscore = np.abs(dst) / stddev
    stats = [dst, estimate, stddev, zscore, abba, baba, arr.shape[0]]
    return np.array(stats), boots



@numba.jit(nopython=True)
def _get_signif_5(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """

    statsarr = np.zeros((3, 7), dtype=np.float64)
    bootsarr = np.zeros((3, nboots))

    idx = 0
    for acol in [2, 3, 4]:
        rows = np.array([0, 1, acol, 5])
        tarr = arr[:, rows, :]

        abxa, baxa, dst = _prop_dstat(tarr)
        boots = _get_boots(tarr, nboots)
        estimate, stddev = (boots.mean(), boots.std())
        if stddev:
            zscore = np.abs(dst) / stddev
        else:
            zscore = np.NaN
        stats = [dst, estimate, stddev, zscore, abxa, baxa, arr.shape[0]]

        statsarr[idx] = stats
        bootsarr[idx] = boots
        idx += 1

    return statsarr, bootsarr




######################################################################
## Simulation functions (requires msprime)
######################################################################


class Sim(object):
    def __init__(self, names, sims, nreps, debug):
        self.names = names
        self.sims = sims
        self.nreps = nreps
        self.debug = debug



def _simulate(self, nreps, admix=None, Ns=500000, gen=20):
    """
    Enter a baba.Tree object in which the 'tree' attribute (newick 
    derived tree) has edge lengths in units of generations. You can 
    use the 'gen' parameter to multiply branch lengths by a constant. 

    Parameters:
    -----------

    nreps: (int)
        Number of reps (loci) to simulate under the demographic scenario
    tree: (baba.Tree object)
        A baba.Tree object initialized by calling baba.Tree(*args). 
    admix: (list)
        A list of admixture events to occur on the tree. Nodes must be 
        reference by their index number, and events must occur in time
        intervals when edges exist. Use the .draw() function of the 
        baba.Tree object to see node index numbers and coalescent times.
    Ns: (float)
        Fixed effective population size for all lineages (may allow to vary
        in the future). 
    gen: (int)
        A multiplier applied to branch lengths to scale into units of 
        generations. Example, if all edges on a tree were 1 then you might
        enter 50000 to multiply so that edges are 50K generations long.

    """

    ## node ages
    Taus = np.array(list(set(self.verts[:, 1]))) * 1e4 * gen

    ## The tips samples, ordered alphanumerically
    ## Population IDs correspond to their indexes in pop config
    ntips = len(self.tree)
    #names = {name: idx for idx, name in enumerate(sorted(self.tree.get_leaf_names()))}
    ## rev ladderized leaf name order (left to right on downward facing tree)
    names = {name: idx for idx, name in enumerate(self.tree.get_leaf_names()[::-1])}
    pop_config = [
        ms.PopulationConfiguration(sample_size=2, initial_size=Ns)
        for i in range(ntips)
    ]

    ## migration matrix all zeros init
    migmat = np.zeros((ntips, ntips)).tolist()

    ## a list for storing demographic events
    demog = []

    ## coalescent times
    coals = sorted(list(set(self.verts[:, 1])))[1:]
    for ct in range(len(coals)):
        ## check for admix event before next coalescence
        ## ...
        
        ## print coals[ct], nidxs, time
        nidxs = np.where(self.verts[:, 1] == coals[ct])[0]
        time = Taus[ct+1]

        ## add coalescence at each node
        for nidx in nidxs:
            node = self.tree.search_nodes(name=str(nidx))[0]

            ## get destionation (lowest child idx number), and other
            dest = sorted(node.get_leaves(), key=lambda x: x.idx)[0]
            otherchild = [i for i in node.children if not 
                          i.get_leaves_by_name(dest.name)][0]

            ## get source
            if otherchild.is_leaf():
                source = otherchild
            else:
                source = sorted(otherchild.get_leaves(), key=lambda x: x.idx)[0]
            
            ## add coal events
            event = ms.MassMigration(
                        time=int(time),
                        source=names[source.name], 
                        destination=names[dest.name], 
                        proportion=1.0)
            #print(int(time), names[source.name], names[dest.name])
        
            ## ...
            demog.append(event)
            
            
    ## sim the data
    replicates = ms.simulate(
        population_configurations=pop_config,
        migration_matrix=migmat,
        demographic_events=demog,
        num_replicates=nreps,
        length=100, 
        mutation_rate=1e-8)
    return replicates



## simulates data on 12 taxon tree with two admixture events
def _sim_admix_12(nreps, Ns=500000, gen=20):
    
    # Set the ML values of various parameters
    Taus = np.array([0, 1, 2, 3, 4, 5]) * 1e4 * gen

    # Migration rates C -> B and from IJ -> EF
    m_C_B = 2e-6
    m_IJ_EF = 2e-6
    
    # Population IDs correspond to their indexes in pop_config.
    ntips = len(tree.tree)
    pop_config = [
        ms.PopulationConfiguration(sample_size=2, initial_size=Ns)
        for i in range(ntips)]
    
    ## migration matrix all zeros time=0
    migmat = np.zeros((ntips, ntips)).tolist()
    
    ## set up demography
    demog = [
        ## initial migration from C -> B
        ms.MigrationRateChange(time=0, rate=m_C_B, matrix_index=(1, 2)),
        ms.MigrationRateChange(time=Taus[1], rate=0),

        # merge events at time 1 (b,a), (f,e), (j,i)
        ms.MassMigration(time=Taus[1], source=1, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[1], source=5, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[1], source=9, destination=8, proportion=1.0), 
        
        ## migration from IJ -> EF (backward in time)
        ms.MigrationRateChange(time=Taus[1], rate=m_IJ_EF, matrix_index=(4, 8)), 

        ## merge events at time 2 (c,a), (g,e), (k,i)
        ms.MassMigration(time=Taus[2], source=2, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[2], source=6, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[2], source=10, destination=8, proportion=1.0), 

        ## end migration at ABC and merge
        ms.MigrationRateChange(time=Taus[2], rate=0),
        ms.MassMigration(time=Taus[3], source=3, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[3], source=7, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[3], source=11, destination=8, proportion=1.0),   
        
        ## merge EFJH -> IJKL
        ms.MassMigration(time=Taus[4], source=8, destination=4, proportion=1.0),   
        
        ## merge ABCD -> EFJHIJKL
        ms.MassMigration(time=Taus[5], source=4, destination=0, proportion=1.0),   
    ]

    ## sim the data
    replicates = ms.simulate(
        population_configurations=pop_config,
        migration_matrix=migmat,
        demographic_events=demog,
        num_replicates=nreps,
        length=100, 
        mutation_rate=1e-9)
    
    return replicates



#def _msp_to_arr(simreps, test):
def _msp_to_arr(Sim, test):
    
    ## the fixed tree dictionary
    #fix = {j: [i, i+1] for j, i in zip(list("abcdefghijkl"), range(0, 24, 2))}
    fix = {j: [i, i+1] for j, i in zip(Sim.names, range(0, len(Sim.names)*2, 2))}
    
    ## fill taxdict by test
    keys = ['p1', 'p2', 'p3', 'p4']
    arr = np.zeros((Sim.nreps, 4, 100))
    
    ## unless it's a 5-taxon test
    if len(test) == 5:
        arr = np.zeros((100000, 6, 100))
        keys += ['p5']
    
    ## create array sampler for taxa
    taxs = [test[key] for key in keys]
    idxs = [list(itertools.chain(*[fix[j] for j in i])) for i in taxs]

    ## iterate over reps filling arr
    idx = 0
    for trees in simreps:
        
        ## build genotype array
        shape = trees.get_num_mutations(), trees.get_sample_size()
        garr = np.empty(shape, dtype="u1")
    
        ## fill the garr
        for variant in trees.variants():
            garr[variant.index] = variant.genotypes
        
        if len(test) == 4:
            if garr.shape[0]:
                ## fill my arr with freqs
                for pdx, tax in enumerate(idxs):
                    freq = garr[:, tax]
                    freq = freq.sum(axis=1) / float(freq.shape[1])
                    maxsz = min(freq.shape[0], 100)
                    arr[idx, pdx, :maxsz] = freq[:maxsz]
        else:
            if garr.shape[0]:
                ## get the easy ones
                p1 = garr[:, idxs[0]]
                p2 = garr[:, idxs[1]]
                p5 = garr[:, idxs[4]]
                p34 = garr[:, idxs[2]+idxs[3]]

                ## identity of SNPs is important
                p3 = garr[:, idxs[2]]
                p4 = garr[:, idxs[3]]
                
                ## any rows with data in b are masked in a
                mask3 = np.where(p3.sum(axis=1) == 0)[0]
                mask4 = np.where(p4.sum(axis=1) == 0)[0]
                masked_p3 = p3[mask4]
                masked_p4 = p4[mask3]
                
                ## enter frequencies
                freq = p1
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 0, :maxsz] = freq[:maxsz]
                
                freq = p2
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 1, :maxsz] = freq[:maxsz]
               
                freq = masked_p3
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 2, :maxsz] = freq[:maxsz]               
               
                freq = masked_p4
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 3, :maxsz] = freq[:maxsz]
               
                freq = p34
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 4, :maxsz] = freq[:maxsz]
                 
                freq = p5
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 5, :maxsz] = freq[:maxsz]
        idx += 1

    ## reduce the size of arr to min loci        
    arr = arr[:idx+1]
    
    ## reduce the size of arr to min len
    minl = np.where(np.all(np.all(arr==0, axis=1) == True, axis=0))[0]
    if np.any(minl):
        minl = minl.min()
    else:
        minl = None
    arr = arr[:, :, :minl]
    
    return arr



## combines sim_admix12 + msp_to_arr + baba to return single (stats, boots)
def sim_admix_12_baba(nreps, test, mindict, nboots):
    sims = _sim_admix_12(nreps)
    arr = _msp_to_arr(sims, test)
    stats, boots = baba(arr, test, mindict, nboots)
    return stats, boots



class TreeParser:
    def __init__():
        print("Stub class to allow baba v0.7 to load")


#######################################################################
if __name__ == "__main__":

    ## test input files
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/" \
              +"vib_half_64tip_c85d6m4p99.loci"

    # ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006',
                 'sulcatum_D9_MEX_003',
                 'jamesonii_D12_PWS_1636',
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    ## calculate dstats
