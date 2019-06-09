#!/usr/bin/env python

"""
Extract, filter and format a locus from a scaff, start, end tuple for 
refmapped assembly data.
"""

# py2/3 compat
from __future__ import print_function

# standard lib
import os
import h5py
import numpy as np
import pandas as pd

from .utils import count_snps
from ipyrad.assemble.utils import IPyradError



class WindowExtracter(object):
    """
    Extract and filter a sequence alignment from a genomic window. You can 
    subselect or filter to include certain samples, and to require that data
    are not missing across some number of them, or to relabel them...

    Parameters:
    -----------
    ...
    """
    def __init__(
        self, 
        data, 
        workdir="analysis-window_extracter",
        scaffold_idx=None, 
        start=None, 
        end=None, 
        exclude=None,
        quiet=False,
        ):
        
        # store params
        self.data = data
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.scaffold_idx = scaffold_idx
        self.start = start
        self.end = end
        self.exclude = (exclude if exclude else [])

        # file to write to
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.outfile = os.path.join(
            self.workdir, 
            "scaf{}-{}-{}.phy".format(
                self.scaffold_idx,
                self.start,
                self.end,
            )
        )

        # parse scaffold info
        self.scaffold_table = None
        self.names = []

        # todo: get sidxs here...
        self._parse_scaffolds()

        if self.scaffold_idx is not None:
            # check requested window is in scaff
            self._check_window()

            # pull the sequence from the database
            self.phymap = None
            self._parse_scaffold_phymap(self.scaffold_idx)         

            # report stats on window (ntaxa all missing; nsnps, ntrimmed sites, ..)
            self._get_stats()

            try:
                self.stats = pd.DataFrame({
                    "Scaffold": [self.scaffold_table.scaffold_name[self.scaffold_idx]],
                    "Start": [self.start],
                    "End": [self.end],
                    "nSites": [self._sites],
                    "nSNPs": [self._nsnps],
                    "Missing": [round(self._propn, 2)],
                    "Samples": [self._nsamplecov],
                    "Empty/drop samples": [len(self._dropped)],
                })

            except AttributeError:
                self.stats = None
                print("No data in selected window.")


    def run(self, force=False):
        """
        Write sequence alignment to a file 
        """
        if force or (not os.path.exists(self.outfile)):
            write_phydict_to_phy(self._phydict, self.outfile)
        print("Wrote data to {}".format(self.outfile))


    def _get_stats(self):

        # get mask to select window array region
        mask = (self.phymap.pos0 > self.start) & (self.phymap.pos1 < self.end)
        cmap = self.phymap.values[mask, :]

        # is there any data at all?
        if not cmap.size:
            return
        wmin = cmap[:, 1].min()
        wmax = cmap[:, 2].max()

        # extract array from window        
        with h5py.File(self.data) as io5:
            seqarr = io5["phy"][self.sidxs, wmin:wmax]
            
        # is there any data at all?
        if not seqarr.size:
            return

        # calculate stats on seqarr
        all_ns = np.all(seqarr == 78, axis=1)
        self._nsamplecov = seqarr.shape[0] - all_ns.sum()
        self._nsnps = count_snps(seqarr)
        self._sites = seqarr.shape[1]

        # drop all N samples(?)
        pnames = self._pnames[~all_ns]
        pseqs = seqarr[~all_ns, :]
        self._propn = np.sum(seqarr == 78) / seqarr.size
        self._dropped = self._pnames[all_ns]

        # return dictionary
        self._phydict = {}
        for name, seq in zip(pnames, pseqs):
            name = name + " " * (self._longname - len(name))
            seq = b"".join(seq.view("S1")).decode()
            self._phydict[name] = seq


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


    def _check_window(self):
        "check that scaf start end is in scaffold_table"
        assert self.end > self.start, "end must be > start"
        assert self.end <= self.scaffold_table.scaffold_length[self.scaffold_idx], \
            "end is beyond scaffold length"


    def _parse_scaffolds(self):
        "get chromosome lengths from the database"
        with h5py.File(self.data) as io5:

            # parse formatting from db
            self._pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
            self.names = [i.strip() for i in self._pnames]
            self.sidxs = [
                i for (i, j) in enumerate(self.names) if j not in self.exclude]
            self.names = [
                j for (i, j) in enumerate(self.names) if i in self.sidxs]
            self._pnames = self._pnames[self.sidxs]
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



def write_phydict_to_phy(phydict, outhandle):

    # dress up as phylip format
    phylip = ["{} {}".format(len(phydict), len(next(iter(phydict.values()))))]
    for name, seq in phydict.items():
        phylip.append("{} {}".format(name, seq))
    phystring = ("\n".join(phylip))

    # write to temp file
    with open(outhandle, 'w') as temp:
        temp.write(phystring)
