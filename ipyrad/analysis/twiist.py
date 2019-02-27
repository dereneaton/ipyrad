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
# import subprocess as sps

# third party
import pandas as pd
import numpy as np
import toytree
import toyplot


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
        self.rmap = {}        
        for k, v in self.imap.items():
            for i in v:
                self.rmap[i] = k

        # store order of imaps in results (abcd; acbd, adbc)
        abcd = sorted(self.imap.keys())
        self.wmap = {
            "abcd": ",".join(abcd[:2]) + "|" + ",".join(abcd[2:]), 
            "acbd": ",".join([abcd[0], abcd[2]]) + "|" + ",".join([abcd[1], abcd[3]]),
            "adbc": ",".join([abcd[0], abcd[3]]) + "|" + ",".join([abcd[1], abcd[2]]),
            "unknown": "unknown",
        }


    # TODO: collapse engines if there is a quit or failure...
    def run(self, ipyclient, force=False):
        "distribute jobs in parallel client"

        # do not overwrite tree table
        tree_weights_path = os.path.join(
            self.workdir, 
            self.name + ".tree_weights.csv")
        if os.path.exists(tree_weights_path):
            if not force:
                print((
        "tree_weights table loaded from {}; Use force to instead overwrite."
        .format(tree_weights_path)))
                self.tree_weights = pd.read_csv(tree_weights_path, index_col=0)
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
        message = "calculating tree weights | {}".format(self.name)
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
            progressbar(done, nwindows, time0, message)
            time.sleep(0.5)
            if not rasyncs:
                print("")
                break

        # frequencies excluding unresolved
        sums = self.tree_weights[["abcd", "acbd", "adbc"]].sum(axis=1)       
        self.tree_weights["fabcd"] = np.nan
        self.tree_weights["facbd"] = np.nan
        self.tree_weights["fadbc"] = np.nan
        self.tree_weights.loc[sums > 0, "fabcd"] = (
            self.tree_weights.loc[sums > 0, "abcd"] / sums[sums > 0])
        self.tree_weights.loc[sums > 0, "facbd"] = (
            self.tree_weights.loc[sums > 0, "acbd"] / sums[sums > 0])
        self.tree_weights.loc[sums > 0, "fadbc"] = (
            self.tree_weights.loc[sums > 0, "adbc"] / sums[sums > 0])

        # frequencies with unresolved
        sums = self.tree_weights[["abcd", "acbd", "adbc", "unk"]].sum(axis=1)       
        self.tree_weights["uabcd"] = np.nan
        self.tree_weights["uacbd"] = np.nan
        self.tree_weights["uadbc"] = np.nan
        self.tree_weights["uunk"] = np.nan        
        self.tree_weights.loc[sums > 0, "uabcd"] = (
            self.tree_weights.loc[sums > 0, "abcd"] / sums[sums > 0])
        self.tree_weights.loc[sums > 0, "uacbd"] = (
            self.tree_weights.loc[sums > 0, "acbd"] / sums[sums > 0])
        self.tree_weights.loc[sums > 0, "uadbc"] = (
            self.tree_weights.loc[sums > 0, "adbc"] / sums[sums > 0])
        self.tree_weights.loc[sums > 0, "uunk"] = (
            self.tree_weights.loc[sums > 0, "unk"] / sums[sums > 0])



    def draw_tree_weights(self, 
        minsnps=4, 
        window=25, 
        min_periods=2, 
        label=None, 
        signif=3,
        include_unresolved=False,
        ):
        """
        Draw a sliding window barplot of subtree weights. 
        Parameters:
        -----------
        minsnps: exclude windows with less than minsnps.
        window: N rows to group when computing rolling means of tree weights.
        min_periods: minimum number of rows with data in a window. 
        label: plot label.
        signif: bars highlight regions with greater propotion of one subtree.
        """

        # include unresolved or not
        if include_unresolved:
            tests = ["uabcd", "uacbd", "uadbc", "uunk"]
            tlabels = ["abcd", "acbd", "adbc", "unknown"]
            nlines = 3
        else:
            tests = ["fabcd", "facbd", "fadbc"]
            tlabels = ["abcd", "acbd", "adbc"]
            nlines = 2

        # make dataframe copy and censor low snps windows
        df = self.tree_weights.loc[:, tests].copy()
        df[self.tree_table.nsnps < minsnps] = np.nan

        # get rolling window means 
        fills = df.rolling(
            window=window, 
            min_periods=min_periods, 
            win_type="boxcar", 
            center=True).mean()

        # toyplot drawing
        canvas = toyplot.Canvas(width=900, height=250)
        axes = canvas.cartesian(
            label=label,
            xlabel="Position (Mb)",
            ylabel="Subtree weighting",
        )
        m = axes.fill(fills, 
            baseline="stacked", 
            opacity=0.4, 
            title=[self.wmap[i] for i in tlabels],
        )

        # add line boundaries between fill regions
        for line in range(nlines):
            axes.plot(
                np.sum([fills[i] for i in tests[:line + 1]], axis=0),
                stroke_width=1.25,
            )

        # mark significantly deviating regions
        for col in tests:
            deviant = fills[col].mean() + (signif * fills[col].std())
            tops = np.where(fills[col] > deviant)[0]
            m = axes.scatterplot(
                tops, 
                np.repeat(1.05, tops.size), 
                marker="r1x0.25", size=5)

        # style axes
        indexer = np.arange(0, self.tree_table.index.max(), 250)
        axes.x.ticks.show = True
        axes.x.ticks.locator = toyplot.locator.Explicit(
            locations=np.arange(0, self.tree_table.index.max(), 250),
            labels=(
                self.tree_table.start[indexer] / int(1e6)
                ).astype(int),
        )
        return canvas, axes


def calculate_weights(imap, subtree):
    "calculate tree weights by subsampling quartet trees for imap clades"

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

        # get splits in the tree that are length (2, 2)
        cache = dtree.treenode.get_cached_content("name")
        edges = [i for i in dtree.treenode.get_edges(cache) if len(i[0]) == 2]
        
        # unresolved subtree
        if not edges:
            arr[3] += 1
            
        # store resolved tree
        else:
            split = [i for i in edges[0] if tips[0] in i][0]
            split.remove(tips[0])
            nidx = tips.index(split.pop())
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
    intmap = {i: 0 for i in minmap}

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

