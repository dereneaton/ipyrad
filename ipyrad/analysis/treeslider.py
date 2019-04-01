#!/usr/bin/env python

"Sliding window (or sampling window) for phylo inference"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard
import os
import sys
import time
import datetime
import tempfile

# third party
from numba import njit
import pandas as pd
import numpy as np
import toytree
from .raxml import Raxml as raxml
from ..core.Parallel import Parallel

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


"""
Infer whole scaffold if windowsize = 0, None
scaffold_idx = 0 is default.
extract PHY and trim for a given window entered...
"""


class TreeSlider():
    """
    Performs phylo inference across RAD data sampled in windows. Uses the
    hdf5 database output from ipyrad as input (".seqs.hdf5"). If no window 
    size is entered then entire scaffolds are used as windows. 

    Parameters:
    ------------
    name: name prefix for output files.
    workdir: directory for output files.
    data: .loci file
    imap: optional dictionary mapping of sample names to new names. 
    minmap: optional dictionary of minimum sampling per imap group.
    minsnps: minimum number of SNPs to include window in analysis.
    """
    def __init__(
        self,
        data,
        name=None,
        workdir="./analysis-treeslider",
        window_size=None, 
        slide_size=None,
        scaffold_idxs=None,
        minsnps=1,
        bootstraps=0,
        # imap=None,
        # minmap=None,
        ):

        # store attributes
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))

        # work
        self.scaffold_idxs = scaffold_idxs
        self.window_size = window_size
        self.slide_size = slide_size
        self.minsnps = minsnps
        self.boots = bootstraps

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

        # if not window then slide is set to window
        if (not self.window_size) or (not self.slide_size):
            self.slide_size = self.window_size

        # # parse attributes
        # self.imap = imap
        # self.minmap = minmap
        # self.rmap = {}
        # if self.imap:
        #     for k, v in self.imap.items():
        #         for i in v:
        #             self.rmap[i] = k

        # ensure workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # # fill mindict
        # if not minmap:
        #     if imap:
        #         self.minmap = {i: 1 for i in self.imap}

        # parsed attributes
        self.scaffold_table = None
        self.tree_table = None
        self.phymap = None
        self._pnames = None
        self._parameter_check()

        # default to all scaffolds if none entered.
        self._parse_scaffolds()
        if not self.scaffold_idxs:
            self.scaffold_idxs = self.scaffold_table.index.to_list()
        if isinstance(self.scaffold_idxs, (list, tuple, set)):
            self.scaffold_idxs = sorted(self.scaffold_idxs)

        # build the tree table from the scaffolds, windows, and slides.
        self._parse_tree_table()


    def _parameter_check(self):
        assert os.path.exists(self.data), "database file not found"
        assert self.data.endswith(".seqs.hdf5"), (
            "data must be '.seqs.hdf5' file.")


    def _parse_scaffolds(self):
        "get chromosome lengths from the database"
        with h5py.File(self.data) as io5:

            # parse formatting from db
            self._pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
            self._longname = 1 + max([len(i) for i in self._pnames])

            # parse names and lengths from db
            scafnames = [i.decode() for i in io5["scaffold_names"][:]]
            scaflens = io5["scaffold_lengths"][:]
            self.scaffold_table = pd.DataFrame(
                data={
                    "scaffold_name": scafnames,
                    "scaffold_length": scaflens,
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )



    def _parse_tree_table(self):
        "Build DataFrame for storing results"
        dfs = []
        for scaffold in self.scaffold_idxs:

            # get the length of this scaffold
            chromlen = (
                self.scaffold_table.loc[scaffold, "scaffold_length"])
            
            # get start positions
            if self.window_size:
                starts = np.arange(0, chromlen - self.window_size, self.slide_size)
            else:
                starts = [0]

            # get end positions
            if self.window_size:
                ends = np.arange(self.window_size, chromlen, self.slide_size)
            else:
                ends = [chromlen]

            # build to table
            df = pd.DataFrame(data={
                "scaffold": scaffold,
                "start": starts,
                "end": ends,
                "nsnps": np.nan,
                "nsamplecov": np.nan,
                "tree": np.nan,
                }, 
                columns=["scaffold", "start", "end", "nsnps", "nsamplecov", "tree"], 
            )
            dfs.append(df)

        # concat data from one or more scaffolds
        self.tree_table = pd.concat(dfs)
        self.tree_table.reset_index(drop=True, inplace=True)



    def _parse_scaffold_phymap(self, scaffold_idx):
        "scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table"
        with h5py.File(self.data) as io5:
            colnames = io5["phymap"].attrs["columns"]

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == self.scaffold_idx + 1

            # load dataframe of this scaffold
            self.phymap = pd.DataFrame(
                data=io5["phymap"][mask, :],
                columns=[i.decode() for i in colnames],
            )
     


    def run(
        self, 
        ipyclient=None, 
        force=False, 
        show_cluster=False, 
        auto=False):
        """
        Run command for tree slider.
        """
        pool = Parallel(
            tool=self, 
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force},
            )
        pool.wrap_run()


    def _run(
        self, 
        force=False,
        ipyclient=None,
        ):
        """
        Hidden func to run parallel job.
        """

        # use user name else create one
        if not self.name:
            self.name = "all-scaffolds"

        # get outfile name
        tree_table_path = os.path.join(
            self.workdir,
            "{}.tree_table.csv".format(self.name))

        # do not overwrite tree table
        if os.path.exists(tree_table_path):
            if not force:
                print((
        "Existing tree table loaded from {}; Use force to instead overwrite."
        .format(tree_table_path)))
                return

        # load balance parallel jobs 2-threaded
        lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::2])

        # distribute jobs on client
        time0 = time.time()
        rasyncs = {}

        # initial progress ticker to run during job submission
        print(
            "building database: nwindows={}; minsnps={}"
            .format(
                self.tree_table.shape[0],
                self.minsnps,
            ))

        # subset phymap for this scaffold
        io5 = h5py.File(self.data, 'r')
        scaffs = io5["phymap"][:, 0]

        # submit jobs: (fname, scafidx, minpos, maxpos, minsnps, boots)
        jidx = 0
        for scaff in self.tree_table.scaffold.unique():

            # load phymap for each scaff at a time
            phymap = io5["phymap"][scaffs == scaff + 1, :]
            subtree = self.tree_table[self.tree_table.scaffold == scaff]

            # submit jobs each window at a time
            for idx in subtree.index:

                # get window margins
                start, stop = (
                    subtree.loc[idx, ["start", "end"]].astype(int))

                # subset phymap for this window range
                mask = (phymap[:, 3] > start) & (phymap[:, 4] < stop)
                cmap = phymap[mask, :]
                wmin = cmap[:, 1].min()
                wmax = cmap[:, 2].max()

                # store async result
                args = (self.data, wmin, wmax, self.minsnps, self.boots)
                rasyncs[jidx] = lbview.apply(remote_tree_inference, *args)
                jidx += 1

        # close database
        io5.close()

        # track progress and save result table
        message = "inferring raxml trees"
        self._track_progress_and_store_results(rasyncs, time0, message)
        self.tree_table.to_csv(tree_table_path)


    def _track_progress_and_store_results(self, rasyncs, time0, message):
        # track progress and collect results.
        nwindows = self.tree_table.shape[0]
        done = 0
        while 1:
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                if rasyncs[idx].successful():
                    self.tree_table.iloc[idx, 3:] = rasyncs[idx].get()
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



def remote_tree_inference(db, wmin, wmax, minsnps, boots, threads=2):
    "remote job to run raxml inference."

    # get seqarr for phy indices
    phystring = None
    tree = np.nan
    nsamplecov = np.nan
    nsnps = np.nan

    # parse phylip string if there is sequence in this window
    if wmax > wmin:
        nsnps, nsamplecov, phystring = parse_phystring(db, wmin, wmax)

        # infer a tree if there is variation in this window
        if nsamplecov >= 4 and nsnps >= minsnps:
            fname = os.path.join(
                tempfile.gettempdir(), str(os.getpid()) + ".tmp")
            with open(fname, 'w') as temp:
                temp.write(phystring)

            # init raxml object and run with blocking
            rax = raxml(
                data=fname,
                name="temp_" + str(os.getpid()),
                workdir=tempfile.gettempdir(),
                T=max(1, threads),
                N=max(1, boots),
            )
            rax.run(force=True, quiet=True, block=True)

            if os.path.exists(rax.trees.bipartitions):
                tree = toytree.tree(rax.trees.bipartitions).newick
            else:
                tree = toytree.tree(rax.trees.bestTree).newick

    return nsnps, nsamplecov, tree



def parse_phystring(database, wmin, wmax):
    "Returns phystring else 0 if filtered."

    with h5py.File(database, 'r') as io5:

        # select all sequence data in window range
        seqarr = io5["phy"][:, wmin:wmax]  # cmap[:, 1].min():cmap[:, 2].max()]
        pnames = np.array([
            i.decode() for i in io5["phymap"].attrs["phynames"]
        ])
        longname = 1 + max([len(i) for i in pnames])

    # calculate stats on seqarr
    all_ns = np.all(seqarr == 78, axis=1)
    nsamplecov = seqarr.shape[0] - all_ns.sum()
    nsnps = count_snps(seqarr)

    # dress up as phylip format, drop all N samples
    pnames = pnames[~all_ns]
    pseqs = seqarr[~all_ns, :]
    phylip = ["{} {}".format(len(pnames), pseqs.shape[1])]
    for name, seq in zip(pnames, pseqs):
        phylip.append("{} {}{}".format(
            name,
            " " * (longname - len(name)),
            b"".join(seq.view("S1")).decode(),
        ))
    phystring = ("\n".join(phylip))
    return nsnps, nsamplecov, phystring



def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
          .format(hashes + nohash, int(progress), elapsed, message),
          end="")
    sys.stdout.flush()



@njit
def count_snps(seqarr):
    nsnps = 0
    for site in range(seqarr.shape[1]):
        # make new array
        catg = np.zeros(4, dtype=np.int16)

        ncol = seqarr[:, site]
        for idx in range(ncol.shape[0]):
            if ncol[idx] == 67:    # C
                catg[0] += 1
            elif ncol[idx] == 65:  # A
                catg[1] += 1
            elif ncol[idx] == 84:  # T
                catg[2] += 1
            elif ncol[idx] == 71:  # G
                catg[3] += 1
            elif ncol[idx] == 82:  # R
                catg[1] += 1       # A
                catg[3] += 1       # G
            elif ncol[idx] == 75:  # K
                catg[2] += 1       # T
                catg[3] += 1       # G
            elif ncol[idx] == 83:  # S
                catg[0] += 1       # C
                catg[3] += 1       # G
            elif ncol[idx] == 89:  # Y
                catg[0] += 1       # C
                catg[2] += 1       # T
            elif ncol[idx] == 87:  # W
                catg[1] += 1       # A
                catg[2] += 1       # T
            elif ncol[idx] == 77:  # M
                catg[0] += 1       # C
                catg[1] += 1       # A
        # get second most common site
        catg.sort()

        # if invariant e.g., [0, 0, 0, 9], then nothing (" ")
        if catg[2] > 1:
            nsnps += 1
    return nsnps



# proc = subprocess.Popen([
#     self.raxml_binary, 
#     "--msa", fname, 
#     "--model", "JC", 
#     "--threads", "1", 
#     "--redo",
#     ], 
#     stderr=subprocess.PIPE, 
#     stdout=subprocess.PIPE,
# )
# out, _ = proc.communicate()
# if proc.returncode:
#     raise Exception("raxml error: {}".format(out.decode()))
# tre = toytree.tree(fname + ".raxml.bestTree")