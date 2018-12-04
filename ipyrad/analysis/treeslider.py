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
import itertools

# third party
from numba import njit
import pandas as pd
import numpy as np
import toytree
from .raxml import Raxml as raxml

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class TreeSlider():
    """
    Performs phylo inference across RAD data sampled in windows. Uses the 
    hdf5 database output from ipyrad as input (".seqs.hdf5"). 

    Parameters:
    ------------
    name: name prefix for output files.
    workdir: directory for output files.
    data: .loci file
    minmap: dictionary of minimum sampling to include window in analysis.
    minsnps: minimum number of SNPs to include window in analysis.
    """
    def __init__(
        self, 
        data,
        name="test",
        workdir="./analysis-treeslider",
        window_size=50000, 
        slide_size=10000,
        scaffold_idx=None,
        imap=None,
        minmap=None,
        minsnps=2,
        ):
       
        # store attributes
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))
        self.window_size = window_size
        self.slide_size = slide_size
        self.minsnps = minsnps
               
        # parse attributes
        self.imap = imap
        self.minmap = minmap
        self.rmap = None
        if self.imap:
            for k, v in self.imap.items():
                for i in v:
                    self.rmap[i] = k
        if not scaffold_idx:
            scaffold_idx = 0
        self.scaffold_idx = scaffold_idx
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # fill mindict
        if not minmap:
            if imap:
                self.minmap = {i: 1 for i in self.imap}               

        # parsed attributes
        self.scaffold_table = None
        self.tree_table = None
        self.phymap = None
        self._pnames = None

        self._parameter_check()
        self._parse_scaffolds()


    def _parameter_check(self):
        assert os.path.exists(self.data), "database file not found"
        assert self.data.endswith(".seqs.hdf5"), (
            "data must be '.seqs.hdf5' file.")


    def _parse_scaffolds(self):
        # get chromosome lengths from the database:
        with h5py.File(self.data) as io5:
            self._pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
            self._longname = 1 + max([len(i) for i in self._pnames])
            self.scaffold_table = pd.DataFrame(
                data={
                "scaffold_name": [i.decode() for i in io5["scaffold_names"][:]],
                "scaffold_length": io5["scaffold_lengths"][:],
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )

    def _parse_scaffold_phymap(self):

        # scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table
        with h5py.File(self.data) as io5:
            colnames = io5["phymap"].attrs["columns"]
            mask = io5["phymap"][:, 0] == self.scaffold_idx + 1
            self.phymap = pd.DataFrame(
                data=io5["phymap"][mask, :],
                columns=[i.decode() for i in colnames],
            )
            

    def _parse_tree_table(self):
        chromlen = self.scaffold_table.loc[self.scaffold_idx, "scaffold_length"]
        self.tree_table = pd.DataFrame(
            data = {
            "start": np.arange(0, chromlen - self.window_size, self.slide_size), 
            "end": np.arange(self.window_size, chromlen, self.slide_size),
            "nsnps": np.nan,
            "nsamplecov": np.nan,
            "tree": np.nan,
            }, 
            columns=["start", "end", "nsnps", "nsamplecov", "tree"], 
        )


    def run(self, scaffold_idx=0, window_size=None, slide_size=None, ipyclient=None):
        
        # re-parse the tree-table in case window and slide were updated: 
        self.scaffold_idx = (scaffold_idx if scaffold_idx else self.scaffold_idx)
        self.window_size = (window_size if window_size else self.window_size)
        self.slide_size = (slide_size if slide_size else self.slide_size)
        self._parse_tree_table()
        self._parse_scaffold_phymap()

        # load balance parallel jbos
        lbview = ipyclient.load_balanced_view()

        # distribute jobs on client
        time0 = time.time()
        rasyncs = {}
        for idx in range(20): #self.tree_table.index
            
            # get window margins
            start, stop = self.tree_table.loc[idx, ["start", "end"]].astype(int)
            
            # store async result
            args = (self, start, stop)
            rasyncs[idx] = lbview.apply(engine_process, *args)

        # track progress and save result table
        self._track_progress_and_store_results(rasyncs, time0)
        self.tree_table.to_csv(
            os.path.join(self.workdir, self.name + ".tree_table.csv"),
            )


    def _track_progress_and_store_results(self, rasyncs, time0):
        # track progress and collect results.
        nwindows = self.tree_table.shape[1]
        done = 0
        while 1:
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                if rasyncs[idx].successful():
                    self.tree_table.loc[idx, :] = rasyncs[idx].get()
                    del rasyncs[idx]
                    done += 1
                else:
                    raise Exception(rasyncs[idx].get())
            # progress
            progressbar(done, nwindows, time0, "inferring raxml trees")
            time.sleep(0.5)
            if not rasyncs:
                break


def engine_process(self, start, stop):
    
    # pull up local phymap
    with h5py.File(self.data) as io5:
        mask = io5["phymap"][:, 0] == self.scaffold_idx + 1
        phymap = io5["phymap"][mask, :]

    # get phy index for window margins
    mask = (phymap[:, 3] > start) & (phymap[:, 4] < stop)
    cmap = phymap[mask, :]
                
    # get seqarr for phy indices
    phystring = None
    tree = np.nan
    nsamplecov = np.nan
    nsnps = np.nan

    # parse phylip string if there is sequence in this window
    if cmap.size:
        nsamplecov, nsnps, phystring = parse_phystring(self, cmap)
    
        # infer a tree if there is variation in this window
        if nsamplecov >= 4 and nsnps >= self.minsnps:
            fname = os.path.join(
                tempfile.gettempdir(), str(os.getpid()) + ".tmp")
            with open(fname, 'w') as temp:
                temp.write(phystring)

            # init raxml object and run with blocking
            rax = raxml(
                data=fname, 
                name="temp_" + str(os.getpid()), 
                workdir=tempfile.gettempdir(),
                T=1,
                )
            rax.run(force=True, quiet=True, block=True)
            tree = toytree.tree(rax.trees.bestTree).newick

    return nsnps, nsamplecov, tree




def parse_phystring(self, cmap):
    "Returns phystring else 0 if filtered."

    with h5py.File(self.data) as io5:
        seqarr = io5["phy"][:, cmap[:, 1].min():cmap[:, 2].max()]

    # calculate stats on seqarr
    allNs = np.all(seqarr == 78, axis=1)
    nsamplecov = seqarr.shape[0] - allNs.sum()
    nsnps = count_snps(seqarr)

    # dress up as phylip format
    pnames = self._pnames[~allNs]
    pseqs = seqarr[~allNs, :]
    phylip = ["{} {}".format(len(pnames), pseqs.shape[1])]
    for name, seq in zip(pnames, pseqs):
        phylip.append("{} {}{}".format(
            name, 
            " " * (self._longname - len(name)), 
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