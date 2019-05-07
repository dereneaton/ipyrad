#!/usr/bin/env python

"Implementation of a twiist like method for reference or anonymous RAD loci"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard
import os
import time
import itertools
# import subprocess as sps

# third party
import pandas as pd
# import toyplot
# import numpy as np
from .utils import progressbar
from ..core.Parallel import Parallel

try:
    import toytree
except ImportError:
    MISSING_IMPORTS = """
To use the ipa.structure module you must install two additional 
libraries which can be done with the following conda command. 

conda install toytree -c eaton-lab
"""



class CladeWeights(object):
    """
    Perform tree weighting on a tree_table.

    Parameters:
    -----------
    ...

    method (str):
        "individual" or "union".

    """
    def __init__(
            self, 
            data, 
            clades,
            name="test", 
            workdir="analysis-clade_weights",
            minsupport=0,
            ):

        # store attrs
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.minsupport = minsupport
        if os.path.exists(data):
            self.tree_table = pd.read_csv(data, index_col=0)
        else:
            self.tree_table = None

        # make sure clades is a dictionary
        if not isinstance(clades, dict):
            raise TypeError("clades argument should be a dictionary.")

        # fill clades dictionary with tests
        self.clades = {}
        for key, value in clades.items():
            
            # make sure all clades values are lists or tuples of strings
            if isinstance(value, str):
                raise TypeError("clades values should be lists or tuples")

            elif isinstance(value, list):
                self.clades[key] = value

            elif isinstance(value, tuple):
                self.clades[key] = value
            
            else:
                raise TypeError("clades argument is malformed.")

        # clades must have >1 tip for individual
        if any([len(i) < 2 for i in self.clades.values()]):
            raise Exception(
                "All clades must have >1 samples in a list or tuple of lists. " + \
                "e.g., ['a', 'b'] or (['a', 'b'], ['c'])"
            )

        # results data frame
        self.clade_weights = None

        # parallelization
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


    def run(self, ipyclient=None, force=False, show_cluster=False, auto=False):
        "Run command for clade weight sampling."
        pool = Parallel(
            tool=self, 
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force},
        )
        pool.wrap_run()


    def _run(self, ipyclient, force=False):
        "distribute jobs in parallel client"

        # make outdir if it doesn't exist and write table to it
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # results file
        clade_weights_path = os.path.join(
            self.workdir, 
            self.name + ".clade_weights.csv")

        # do not overwrite clade weights table
        if os.path.exists(clade_weights_path):
            if not force:
                print((
        "clade_weights table loaded from {}; Use force to instead overwrite."
        .format(clade_weights_path)))
                self.clade_weights = pd.read_csv(clade_weights_path, index_col=0)
                return

        # make empty clade weights table
        self.clade_weights = pd.DataFrame(
            {i: [0.] * self.tree_table.shape[0] for i in self.clades.keys()}
        )

        # setup cluster
        lbview = ipyclient.load_balanced_view()
        time0 = time.time()
        rasyncs = {}

        # distribute jobs on client in chunks
        # TODO: test adaptive chunk sizes.
        chunksize = 20
        for idx in self.tree_table.index[::chunksize]:
            chunk = self.tree_table.tree[idx:idx + chunksize].tolist()

            rasyncs[idx] = lbview.apply(
                clade_weights, *(chunk, self.clades, idx))

        # caculateing clade weights
        done = 0
        nwindows = len(rasyncs)
        message = "calculating tree weights | {}".format(self.name)
        while 1:
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                if rasyncs[idx].successful():
                    # fill with results
                    res = rasyncs[idx].get()
                    # print(res)
                    self.clade_weights.iloc[idx: idx + chunksize] = res.values
                    del rasyncs[idx]
                    done += 1                    
                else:
                    raise Exception(rasyncs[idx].get())
            # progress
            progressbar(done, nwindows, time0, message)
            time.sleep(0.5)
            if not rasyncs:
                print("")
                break

        # progress bar
        self.clade_weights.to_csv(clade_weights_path)        
        return self.clade_weights




def clade_weights(treelist, clades, idx=0):
    """
    Calculate clade weight for [(clade1),(clade1),other,other] by iteratively
    subsampling quartet trees from these clades.
    """
    # make empty clade weights table
    # TODO: should non filled rows be NANs?
    clade_weights = pd.DataFrame(
        {i: [0.] * len(treelist) for i in clades.keys()}
    )

    # iterate over trees
    for tidx, tree in enumerate(treelist):

        # get tips for this subtree
        tree = toytree.tree(tree)
        tips = set(tree.get_tip_labels())        

        # iterate over clades to test
        for name, clade in clades.items():

            # stat counters
            idx = 0
            tsum = 0          

            # union test support for grouping of two clades
            if isinstance(clade, tuple):
                # reduce clades to sampled tips
                tree = toytree.tree(tree)
                tips = set(tree.get_tip_labels())
                iclade1 = tips.intersection(set(clade[0]))
                iclade2 = tips.intersection(set(clade[1]))
                iclade = iclade1.union(iclade2)
                oclade = tips.difference(iclade)

                # get pair samplers
                osamp = itertools.combinations(oclade, 2)
               
                # iterate over quartets
                for c1 in clade[0]:
                    for c2 in clade[1]:
                        for opair in osamp:
                            quartet = set([c1, c2] + list(opair))
                            todrop = set(tree.get_tip_labels()) - quartet
                            dt = tree.drop_tips(todrop)
                            tsum += clade_true(dt.unroot(), iclade)
                            idx += 1
            
            # individual test support for individual clades
            else:
                # get pair samplers
                iclade = tips.intersection(set(clade))
                isamp = itertools.combinations(iclade, 2)
                oclade = tips.difference(iclade)
                osamp = itertools.combinations(oclade, 2)
            
                # iterate over quartets
                for ipair in isamp:
                    for opair in osamp:
                        quartet = set(list(ipair) + list(opair))
                        todrop = set(tree.get_tip_labels()) - quartet
                        dt = tree.drop_tips(todrop)
                        tsum += clade_true(dt.unroot(), iclade)
                        idx += 1
                
            # store result
            if idx:
                clade_weights.loc[tidx, name] = tsum / idx
    return clade_weights



def clade_true(tre, clade):
    """
    Calculate if quartet tre has a split supporting ['a', 'b'] in clade.
    """
    # check for a resolved node
    nodes = [not i.is_leaf() for i in tre.treenode.children]
    
    if True in nodes:
        node = tre.treenode.children[nodes.index(True)]
        tips = [i.name for i in node.children]
        if sum([i in clade for i in tips]) == 1:
            return False
        return True
    else:
        return False
    # root is in a clade
    # tips = [j.name for (i, j) in zip(leaves, tre.treenode.children) if i]
    # if sum([i in clade for i in tips]) == 1:
    #     return False
    # else:
    #     return True
