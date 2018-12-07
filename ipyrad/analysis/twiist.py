#!/usr/bin/env python

"Implementation of a twiist like method for reference or anonymous RAD loci"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard
import os
import sys
import time
import datetime
import itertools
import subprocess as sps

# third party
import pandas as pd
import numpy as np
import toytree



class Twisst:
    """
    Perform tree weighting on a tree_table.
    """
    def __init__(
        self, 
        data, 
        imap, 
        minmap,
        name=None, 
        workdir="analysis-twisst",
        minsnps=1,
        ipyclient=None,
        ):

        # store attrs
        self.imap = imap
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})
        self.name = (name if name else "-".join(sorted(imap.keys())))
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.minsnps = minsnps

        # max clades currently is 4
        assert len(imap) <= 4, "support limited to 4 imap clades currently."

        # parse tree table csv file
        self.tree_table = pd.read_csv(data, index_col=0)

        # attrs to fill
        self.rmap = {}
        self.tree_weights = pd.DataFrame({
            "minmap": np.repeat(False, self.tree_table.shape[0]),
            "subtree": np.nan,
            "nnnodes": np.nan,
            "abcd": np.nan,
            "acbd": np.nan,
            "adbc": np.nan,
            "unk": np.nan
        }, 
        columns=["minmap", "nnodes", "abcd", "acbd", "adbc", "unk", "subtree"],
        )

        # reverse of imap
        for k, v in self.imap.items():
            for i in v:
                self.rmap[i] = k

        # store order of imaps in results (abcd; acbd, adbc)
        abcd = sorted(self.imap.keys())
        self.wmap = {
            "abcd": abcd, 
            "acbd": [abcd[0], abcd[2], abcd[1], abcd[3]],
            "adbd": [abcd[0], abcd[3], abcd[1], abcd[2]],
        }


    def run(self, ipyclient, force=False):
        "distribute jobs in parallel client"

        # do not overwrite tree table
        tree_weights_path = os.path.join(
            self.workdir, self.name + ".tree_weights.csv")
        if os.path.exists(tree_weights_path):
            if not force:
                print((
        "tree_weights table loaded from {}; Use force to instead overwrite."
        .format(tree_weights_path)))
                return

        # setup
        lbview = ipyclient.load_balanced_view()
        time0 = time.time()
        rasyncs = {}

        # initial progress ticker to run during job submission
        print("\rbuilding database...", end="")

        # distribute jobs on client
        done = 0        
        for idx in self.tree_table.index:
            newick = self.tree_table.tree[idx]
            if not pd.isna(newick):
                args = (newick, self.imap, self.rmap, self.minmap)
                rasyncs[idx] = lbview.apply(engine_process, *args)
            else:
                done += 1

        # progress bar
        sys.stdout.flush()
        self._track_progress_and_store_results(rasyncs, time0, done)

        # make outdir if it doesn't exist and write table to it
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.tree_weights.to_csv(tree_weights_path)



    def _track_progress_and_store_results(self, rasyncs, time0, done):
        # track progress and collect results.
        nwindows = self.tree_weights.shape[0]
        while 1:
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                if rasyncs[idx].successful():
                    self.tree_weights.loc[idx, :] = rasyncs[idx].get()
                    del rasyncs[idx]
                    done += 1
                else:
                    raise Exception(rasyncs[idx].get())
            # progress
            progressbar(done, nwindows, time0, "calculating tree weights")
            time.sleep(0.5)
            if not rasyncs:
                break

        # frequencies excluding unresolved
        self.tree_weights["fabcd"] = (self.tree_weights["abcd"] / 
            self.tree_weights[["abcd", "acbd", "adbc"]].sum(axis=1))
        self.tree_weights["facbd"] = (self.tree_weights["acbd"] / 
            self.tree_weights[["abcd", "acbd", "adbc"]].sum(axis=1))
        self.tree_weights["fadbc"] = (self.tree_weights["adbc"] / 
            self.tree_weights[["abcd", "acbd", "adbc"]].sum(axis=1))

        # frequencies with unresolved
        self.tree_weights["uabcd"] = (self.tree_weights["abcd"] / 
            self.tree_weights[["abcd", "acbd", "adbc", "unk"]].sum(axis=1))
        self.tree_weights["uacbd"] = (self.tree_weights["acbd"] / 
            self.tree_weights[["abcd", "acbd", "adbc", "unk"]].sum(axis=1))
        self.tree_weights["uadbc"] = (self.tree_weights["adbc"] / 
            self.tree_weights[["abcd", "acbd", "adbc", "unk"]].sum(axis=1))
        self.tree_weights["uunk"] = (self.tree_weights["unk"] / 
            self.tree_weights[["abcd", "acbd", "adbc", "unk"]].sum(axis=1))



def calculate_weights(imap, subtree):

    # get a generator of unique nkey samples: e.g., quartets 
    tips = sorted(subtree.get_tip_labels(), key=lambda x: x.rsplit("-", 1)[0])
    groups = itertools.groupby(tips, lambda x: x.rsplit("-", 1)[0])
    quarts = itertools.product(*(set(j) for (i, j) in groups))
    
    # record results
    arr = np.zeros(4, dtype=int)
    
    # prune tree and measure topo dists
    for idx, quart in enumerate(quarts):
        
        # get pruned tree
        droptips = set(subtree.get_tip_labels()) - set(quart)
        dtree = subtree.drop_tips(droptips)
                
        # get treenodes of pruned tips always in alphanumeric order of imaps
        tips = sorted(dtree.get_tip_labels())
        nodes = [dtree.treenode.search_nodes(name=i)[0] for i in tips]
        
        # get distances between tips in nnodes
        dists = np.array([
            nodes[0].get_distance(i, topology_only=True) for i in nodes
        ])
        
        # get tip that is closest to the first tip in the list
        one = dists[dists != 0].min()
        nidx = np.argsort(dists)[1]
        
        # if multiple closest then polytomy is unresolved
        if dists[dists == one].size != 1:
            arr[3] += 1
        # store split info
        else:
            arr[nidx - 1] += 1
    return arr
    



def engine_process(newick, imap, rmap, minmap):
    """
    Load toytree for an interval, rename tips to clades with imap, prune
    samples from tree that are not in the test, test for minmap sampling, 
    calculate weights for tree hypothesis, and return tuple of 
    (minmap, weight, subtree).
    """
    # load toytree
    ttree = toytree.tree(newick)
    
    # incremental counter for subtree names from imap
    intmap = {i:0 for i in minmap}

    # keep track of names to drop from tree (not in imap)
    drops = []

    # map names to imap
    for node in ttree.treenode.traverse():
        if node.is_leaf():
            if node.name in rmap:
                imapname = rmap[node.name]
                node.name = "{}-{}".format(imapname, intmap[imapname])
                intmap[imapname] += 1
            else:
                drops.append(node.name)

    # test whether intmap meets minmap filter
    minfilter = [intmap[i] >= minmap[i] for i in minmap]

    # default values
    minmap = False
    weight = np.nan
    nnodes = np.nan
    subtree = np.nan
    abcd = acbd = adbc = unknown = np.nan

    # apply filter
    if not all(minfilter):
        return minmap, nnodes, abcd, acbd, adbc, unknown, subtree
    else:

        # order is defined by alphanumeric order of imap
        abcd = tuple(sorted(imap.keys()))
        acbd = (abcd[0], abcd[2], abcd[1], abcd[3])
        adbc = (abcd[0], abcd[3], abcd[1], abcd[2])

        # calculate weights
        subtree = ttree.drop_tips(drops).collapse_polytomies().unroot()
        abcd, acbd, adbc, null = calculate_weights(imap, subtree)

        # return values
        return True, subtree.nnodes, abcd, acbd, adbc, null, subtree.write()        


def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
        .format(hashes + nohash, int(progress), elapsed, message),
        end="")
    sys.stdout.flush()    

