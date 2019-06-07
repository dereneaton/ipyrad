#!/usr/bin/env python

"D-statistic calculations"

# py2/3 compat
from __future__ import print_function, division
from builtins import range

import os
import sys
import time
import copy
import types
import itertools
from collections import OrderedDict

#import scipy.stats as st  ## used for dfoil
import pandas as pd
import numpy as np
import numba
import h5py

## ipyrad tools
from ipyrad.analysis.utils import Params, progressbar, IPyradError
from ipyrad.assemble.write_outputs import reftrick

# import tested at init
try:
    import toytree
except ImportError:
    pass
_TOYTREE_IMPORT = """
This ipyrad analysis tool requires 
You can install it with the following command:

   conda install toytree -c eaton-lab
"""

# from ipyrad.plotting.baba_panel_plot import baba_panel_plot
# set floating point precision in data frames to 3 for prettier printing
# pd.set_option('precision', 3)

# Notes: treegenerateor working, test others, including toytree integratin.


class Baba:
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
    def __init__(self, 
        data=None,
        imap=None,
        minmap=1,
        newick=None,
        nboots=1000,
        ):

        # check imports
        if not sys.modules.get("toytree"):
            raise ImportError(_TOYTREE_IMPORT)

        # parse data as (1) path to data file, or (2) ndarray
        if isinstance(data, str):
            self.data = os.path.realpath(os.path.expanduser(data))
        else:
            self.data = data

        # check dtype of newick/tree entry
        self.newick = newick
        if isinstance(newick, toytree.Toytree.ToyTree):
            self.newick = newick.newick          

        # store tests
        self.imap = imap
        self.minmap = minmap

        # parameters
        self.params = Params()
        self.params.nboots = nboots
        self.params.quiet = False
        self.params.database = None

        # results storage
        self.results_table = None
        self.results_boots = None
       
        # cluster attributes
        self.ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }

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
            print("no tests generated.")
            return None



    def run(self, force=False, ipyclient=None, show_cluster=False, auto=False):
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

        # distribute jobs in a wrapped cleaner function
        pool = Parallel()

        batch(self, ipyclient)

        ## skip this for 5-part test results
        if not isinstance(self.results_table, list):
            self.results_table.nloci = (
                np.nan_to_num(self.results_table.nloci).astype(int))



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
            raise AttributeError("no newick tree")

        tp = TreeParser(self.newick, constraint_dict, constraint_exact)
        tests = tp.testset
        if verbose:
            print("{} tests generated from tree".format(len(tests)))
        self.tests = tests



    def plot(self, 
        show_test_labels=True, 
        use_edge_lengths=False,         
        collapse_outgroup=False, 
        pct_tree_x=0.5, 
        pct_tree_y=0.2,
        subset_tests=None,
        prune_tree_to_tests=False,
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

        """

        ## check for attributes
        if not self.newick:
            raise IPyradError("baba plot requires a newick treefile")
        if not self.tests:
            raise IPyradError("baba plot must have a .tests attribute")

        ## ensure tests is a list
        if isinstance(self.tests, dict):
            self.tests = [self.tests]

        # re-decompose the tree
        ttree = toytree.tree(self.newick)

        # subset test to show fewer
        if subset_tests is not None:
            #tests = self.tests[subset_tests]
            tests = [self.tests[i] for i in subset_tests]
            boots = self.results_boots[subset_tests]
        else:
            tests = self.tests
            boots = self.results_boots

        ## if prune tree
        if prune_tree_to_tests:
            alltesttaxa = set(itertools.chain(*self.taxon_table.values[0]))
            ttree = ttree.drop_tips([i for i in ttree.get_tip_labels()
                                     if i not in alltesttaxa])
            ttree.tree.ladderize()

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


    def _run(self, force=False, ipyclient=None):
        "Function to distribute jobs to ipyclient"

        # load balancer
        lbview = ipyclient.load_balanced_view()

        # check that tests are OK
        if not self.tests:
            raise IPyradError("no tests found")
        if isinstance(self.tests, dict):
            self.tests = [self.tests]

        # check that mindict is OK
        if isinstance(self.minmap, int):
            self.minmap = {i: self.minmap for i in self.imap}
        if not self.minmap:
            self.minmap = {i: 1 for i in self.imap}



        # send jobs to the client (but not all at once b/c njobs can be huge)
        rasyncs = {}
        idx = 0
        for i in range(len(ipyclient)):

            # next entries unless fewer than len ipyclient, skip
            try:
                test = next(itests)
                mindict = next(imdict)
            except StopIteration:
                continue

            rasyncs[idx] = lbview.apply(dstat, *[loci, test, mindict, self.params.nboots])
            idx += 1


def write_tmp_h5(baba):
    "Reduce VCF to temp h5 that jobs will slice from"

    # load in the VCF: if this gets huge we could hdf5 it...
    with open(baba.data) as indata:
        for i in indata:
            if i[:6] == "#CHROM":
                colnames = i[1:].strip().split()
                break

    # below here could be done in chunks...
    df = pd.read_csv(baba.data, comment="#", sep="\t", names=colnames)

    # drop superfluous columns
    df = df.drop(columns=[
        "QUAL", "FILTER", "INFO", "FORMAT", "REF", "ALT"])

    # reduce geno calls to only genos
    for column in df.columns[3:]:
        df[column] = df[column].apply(lambda x: x.split(":")[0])
    
    # set missing data to NaN
    df.iloc[:, 3:] = df.iloc[:, 3:].applymap(sumto)

    # save as hdf5
    r = 54321
    with h5py.File("baba-{}.h5".format(r), "w") as io5:
        # save chrom as an int index instead of string
        io5["CHROM"] = pd.factorize(df.CHROM)[0]

        # save loc as int 
        io5["LOC"] = [4, 4, 10]
    
    



def sumto(value):
    "used in pd.DataFrame applymap to convert genos to derived sums"
    if value == "./.":
        return np.nan
    else:
        return sum((int(i) for i in value.split("/") if i in ("0", "1"))) / 2.
    

def batch(baba, ipyclient=None):
    """
    distributes jobs to the parallel client
    """
    # parse args
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
        pass  #sims()

    ## iterate over tests (repeats mindicts if fewer than taxdicts)
    if not taxdicts:
        print("no tests found")
        return
    else:
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

    # store instead of return...
    self.results_table, self.results_boots




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




def loci_to_arr(self, loci):
    "Converts sequence data in loci file to binary where outgroup is 0"

    # array dimensions
    nloci = len(loci)
    maxsnps = 10
    testshape = (6 if len(self.imap[0]) > 4 else 4) 

    # make an array
    arr = np.zeros((nloci, testshape, maxsnps), dtype=np.float64)
    
    # get the outgroup sample
    keys = sorted([i for i in self.imap if i[0] == 'p'])
    outg = keys[-1]

    # iterate over loci and store 
    pass




def vcf_to_arr(self, loci):
    "Converts VCF SNP data to binary array where outgroup is 0"









def _loci_to_arr(loci, taxdict, mindict):
    """
    return a frequency array from a loci file for all loci with taxa from 
    taxdict and min coverage from mindict. 
    """

    ## get max length of loci
    maxlen = np.max(np.array([len(locus.split("\n")[0]) for locus in loci]))

    ## make the array (4 or 5) and a mask array to remove loci without cov
    nloci = len(loci)
    keep = np.zeros(nloci, dtype=np.bool_)
    # arr = np.zeros((nloci, 4, 300), dtype=np.float64)
    arr = np.zeros((nloci, 4, maxlen), dtype=np.float64)

    ## six rows b/c one for each p3, and for the fused p3 ancestor
    if len(taxdict) == 5:
        # arr = np.zeros((nloci, 6, 300), dtype=np.float64)
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



class TreeParser:
    def __init__(self, newick, constraint_dict, constraint_exact):
        "Traverses tree to build test sets given constraint options."

        # store sets of four-taxon splits
        self.testset = set()
        self.hold = [0, 0, 0, 0]

        # tree to traverse
        self.tree = toytree.tree(newick)
        if not self.tree.is_rooted(): 
            raise IPyradError(
                "generate_tests_from_tree(): tree must be rooted and resolved")

        # constraints
        self.cdict = OrderedDict((i, []) for i in ["p1", "p2", "p3", "p4"])
        if constraint_dict:
            self.cdict.update(constraint_dict)

        # constraint setting
        self.xdict = constraint_exact
        if isinstance(self.xdict, bool):
            self.xdict = [self.xdict] * 4
        if isinstance(self.xdict, list):
            if len(self.xdict) != len(self.cdict):
                raise Exception(
                    "constraint_exact must be bool or list of bools length N")

        # get tests
        self.loop()


    def loop(self, node, idx):
        "getting closer...."
        for topnode in node.traverse():
            for oparent in topnode.children:
                for onode in oparent.traverse():
                    if self.test_constraint(onode, 3):                       
                        self.hold[3] = onode.idx

                        node2 = oparent.get_sisters()[0]
                        for topnode2 in node2.traverse():
                            for oparent2 in topnode2.children:
                                for onode2 in oparent2.traverse():
                                    if self.test_constraint(onode2, 2):                       
                                        self.hold[2] = onode2.idx

                                        node3 = oparent2.get_sisters()[0]
                                        for topnode3 in node3.traverse():
                                            for oparent3 in topnode3.children:
                                                for onode3 in oparent3.traverse():
                                                    if self.test_constraint(onode3, 1):                       
                                                        self.hold[1] = onode3.idx

                                                        node4 = oparent3.get_sisters()[0]
                                                        for topnode4 in node4.traverse():
                                                            for oparent4 in topnode4.children:
                                                                for onode4 in oparent4.traverse():
                                                                    if self.test_constraint(onode4, 0):
                                                                        self.hold[0] = onode4.idx
                                                                        self.testset.add(tuple(self.hold))


    def test_constraint(self, node, idx):
        names = set(node.get_leaf_names())
        const = set(list(self.cdict.values())[idx])
        if const:
            if self.xdict[idx]:
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


# This should be re-written as a dynamic func
def tree2tests(newick, constraint_dict, constraint_exact):
    """
    Returns dict of all possible four-taxon splits in a tree. Assumes
    the user has entered a rooted tree. Skips polytomies.
    """
    # make tree
    tree = toytree.tree(newick)
    if not tree.is_rooted(): 
        raise IPyradError(
            "Input tree must be rooted to use generate_tests_from_tree()")

    # store results 
    testset = set()

    # constraints fill in empty 
    cdict = OrderedDict((i, []) for i in ["p1", "p2", "p3", "p4"])
    if constraint_dict:
        cdict.update(constraint_dict)

    # expand constraint_exact if list
    if isinstance(constraint_exact, bool):
        constraint_exact = [constraint_exact] * 4

    if isinstance(constraint_exact, list):
        if len(constraint_exact) != len(cdict):
            raise Exception(
                "constraint_exact must be bool or list of bools of length N")
    
    # traverse root to tips. Treat the left as outgroup, then the right.
    tests = []
    
    # topnode must have children. All traversals use default "levelorder"
    for topnode in tree.treenode.traverse():
        
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
    amseq = np.zeros((iseq.shape[0] * 2, iseq.shape[1]), dtype=np.uint8)
    
    ## fill in both copies
    for seq in range(iseq.shape[0]):
        for col in range(iseq.shape[1]):

            ## get this base and check if it is hetero
            base = iseq[seq][col]
            who = consdict[:, 0] == base
            
            ## if not hetero then enter it
            if not np.any(who):
                amseq[seq * 2][col] = base
                amseq[seq * 2 + 1][col] = base        
            ## if hetero then enter the 2 resolutions
            else:
                amseq[seq * 2][col] = consdict[who, 1][0]
                amseq[seq * 2 + 1][col] = consdict[who, 2][0]

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
    abba = ((1. - arr[:, 0]) * (arr[:, 1]) * (arr[:, 2]) * (1. - arr[:, 3]))
    baba = ((arr[:, 0]) * (1. - arr[:, 1]) * (arr[:, 2]) * (1. - arr[:, 3]))
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


#######################################################################

if __name__ == "__main__":

    ## test input files
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/"\
              + "vib_half_64tip_c85d6m4p99.loci"

    # ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006',
                 'sulcatum_D9_MEX_003',
                 'jamesonii_D12_PWS_1636',
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    ## calculate dstats
