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
from ipyrad.assemble.utils import GETCONS
from ipyrad.assemble.write_outputs import reftrick


class WindowExtracter(object):
    """
    Extract and filter a sequence alignment from a genomic window. You can 
    subselect or filter to include certain samples, and to require that data
    are not missing across some number of them, or to relabel them...

    Parameters:
    -----------
    data (str):
        A seqs.hdf5 file from ipyrad.
    name (str):
        Prefix name used for outfiles. If None then it is created from the 
        scaffold, start and end positions.
    workdir (str):
        Location for output files to be written. Created if it doesn't exist.
    scaffold_idx (int):


    """
    def __init__(
        self, 
        data, 
        name=None,
        workdir="analysis-window_extracter",
        scaffold_idx=None, 
        start=None, 
        end=None, 
        mincov=4,
        exclude=None,
        imap=None,
        quiet=False,
        ):

        # store params
        self.data = data
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.scaffold_idx = scaffold_idx
        self.start = (start if start else 0)
        self.end = end
        self.exclude = (exclude if exclude else [])
        self.mincov = mincov
        self.imap = imap
        self.quiet = quiet
        # self.minmap = minmap

        # file to write to
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # use provided name else auto gen a name
        if not name:
            self.name = "scaf{}-{}-{}".format(
                self.scaffold_idx,
                self.start,
                self.end
            )
        else:
            self.name = name

        # parse scaffold info
        self.scaffold_table = None
        self.names = []

        # gets names, pnames, scaffold_table, ...
        self._parse_scaffolds()

        # update end to scaff len if not entered
        if not self.end:
            self.end = int(self.scaffold_table.iloc[self.scaffold_idx, -1])

        # allow user to init with None to get scaffold names.
        if self.scaffold_idx is not None:

            # check requested window is in scaff
            self._check_window()

            # pull the sequence from the database
            self.phymap = None
            self._extract_phymap_df()
            self._init_stats()

            # report stats on window (ntaxa missing; nsnps, ntrimmed sites)
            self._extract_seqarr()
            self._calc_initial_stats()
            self._imap_consensus_reduce()
            self._filter_seqarr()
            self._calc_filtered_stats()

            # no data
            if not self.seqarr.size:
                print("No data in selected window.")


    def _print(self, message):
        if not self.quiet:
            print(message)


    def run(self, force=False):
        """
        Write sequence alignment to a file 
        """
        # bail if user never selected a window.
        if self.scaffold_idx is None:
            return "No scaffold selected"

        # make outfile path name
        self.outfile = os.path.join(
            self.workdir, self.name + ".phy")

        # check for force/overwrite
        if force or (not os.path.exists(self.outfile)):
            # convert to file format
            self.write_to_phy()
            self._print("Wrote data to {}".format(self.outfile))
        else:
            raise IOError(
                "Output file already exists. Use force to overwrite.")


    def _extract_seqarr(self):
        """
        Extracts seqarr of the full selected window and computes stats on
        the window.
        """
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
            self.seqarr = io5["phy"][self.sidxs, wmin:wmax]
            
        # is there any data at all?
        if not self.seqarr.size:
            return


    def _calc_initial_stats(self):
        # enter initial stats
        self.stats.loc["prefilter", "sites"] = self.seqarr.shape[1]
        self.stats.loc["prefilter", "snps"] = count_snps(self.seqarr)
        self.stats.loc["prefilter", "missing"] = round(
            np.sum(self.seqarr == 78) / self.seqarr.size, 2)
        self.stats.loc["prefilter", "samples"] = self.seqarr.shape[0]


    def _imap_consensus_reduce(self):
        """
        Get consensus sequence for all samples in clade. 
        """
        # skip if no imap
        if not self.imap:
            return 

        # empty array of shape imap groups
        iarr = np.zeros((len(self.imap), self.seqarr.shape[1]), dtype=np.uint8)
        inames = np.array(sorted(self.imap.keys()))

        # iterate over imap groups
        for ikey, ivals in self.imap.items():
            
            # get subarray for this group
            sidxs = [np.where(self.names == i)[0][0] for i in ivals]
            subarr = self.seqarr[sidxs, :]

            # get consensus sequence
            cons = reftrick(subarr, GETCONS)[:, 0]

            # insert to iarr
            iidx = np.where(inames == ikey)[0][0]
            iarr[iidx] = cons

        # save as new data
        iarr[iarr == 0] = 78
        self.seqarr = iarr
        self.names = inames
        self._longname = 1 + max([len(i) for i in self.names])
        self._pnames = np.array([
            "{}{}".format(name, " " * (self._longname - len(name)))
            for name in self.names
        ])               


    def _filter_seqarr(self):
        """
        Apply filters to remove taxa.
        """
        # drop sites that are too many Ns given (global) mincov
        keep = np.sum(self.seqarr != 78, axis=0) >= self.mincov
        self.seqarr = self.seqarr[:, keep]

        # drop samples that are all Ns
        keep = np.invert(np.all(self.seqarr == 78, axis=1))
        self.seqarr = self.seqarr[keep, :]
        self.names = self.names[keep]
        self._pnames = self._pnames[keep]


    def _calc_filtered_stats(self):
        # update stats
        self.stats.loc["postfilter", "sites"] = self.seqarr.shape[1]
        self.stats.loc["postfilter", "snps"] = count_snps(self.seqarr)
        self.stats.loc["postfilter", "missing"] = round(
            np.sum(self.seqarr == 78) / self.seqarr.size, 2)
        self.stats.loc["postfilter", "samples"] = self.seqarr.shape[0]


    def _extract_phymap_df(self):
        """
        This extracts the phymap DataFrame.
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table. 
        """
        with h5py.File(self.data) as io5:
            colnames = io5["phymap"].attrs["columns"]

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == self.scaffold_idx + 1

            # load dataframe of this scaffold
            self.phymap = pd.DataFrame(
                data=io5["phymap"][mask, :],
                columns=[i.decode() for i in colnames],
            )


    def _init_stats(self):
        # stats table
        scaf = [self.scaffold_table.scaffold_name[self.scaffold_idx]]
        self.stats = pd.DataFrame({
            "scaffold": scaf,
            "start": [self.start],
            "end": [self.end],
            "sites": [0],
            "snps": [0],
            "missing": [0],
            "samples": [0],
        }, index=["prefilter", "postfilter"],
        )


    def _check_window(self):
        "check that scaf start end is in scaffold_table"
        assert self.end > self.start, "end must be > start"
        assert self.end <= self.scaffold_table.scaffold_length[self.scaffold_idx], \
            "end is beyond scaffold length"


    def _parse_scaffolds(self):
        "get chromosome lengths from the database"
        with h5py.File(self.data) as io5:

            # get sample names
            self._pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
            self.names = [i.strip() for i in self._pnames]

            # filter to only the included samples
            self.sidxs = [
                i for (i, j) in enumerate(self.names) if j not in self.exclude]
            self.names = np.array([
                j for (i, j) in enumerate(self.names) if i in self.sidxs])
            self._pnames = self._pnames[self.sidxs]

            # format names to include spacer for phylip, etc.
            self._longname = 1 + max([len(i) for i in self._pnames])
            self._pnames = np.array([
                "{}{}".format(name, " " * (self._longname - len(name)))
                for name in self._pnames
            ])               

            # parse scaf names and lengths from db
            scafnames = [i.decode() for i in io5["scaffold_names"][:]]
            scaflens = io5["scaffold_lengths"][:]
            self.scaffold_table = pd.DataFrame(
                data={
                    "scaffold_name": scafnames,
                    "scaffold_length": scaflens,
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )


    # def write_to_fasta(self):
    #     """
    #     Split into separate loci, only include samples in a locus if they are
    #     not empty, apply filters for locus inclusion. 
    #     """       
        
    #     # make outfile path name
    #     self.outfile = os.path.join(
    #         self.workdir, self.name + ".fasta")

    #     # build loci
    #     locs = []
    #     for idx in self.phymap.index:
    #         start = self.phymap.phy0[idx]
    #         end = self.phymap.phy1[idx]

    #         loc = []
    #         for sidx in range(self.seqarr.shape[0]):
    #             aseq = self.seqarr[sidx, start:end]
    #             if not np.all(aseq == 78):
    #                 name = self.names[sidx]
    #                 seq = bytes(aseq).decode()
    #                 scaf = self.scaffold_table.scaffold_name[self.scaffold_idx]
    #                 header = ">{} | {}:{}-{}".format(name, scaf, start, end)
    #                 loc.append("{}\n{}".format(header, seq))
    #         locs.append("\n".join(loc))

    #     # write to file
    #     with open(self.outfile, 'w') as out:
    #         out.write("\n\n".join(locs))


    # def write_to_loci(self):
    #     """
    #     Split into separate loci, only include samples in a locus if they are
    #     not empty, apply filters for locus inclusion. 
    #     """           
    #     raise NotImplementedError("coming soon")
    #     # build loci
    #     locs = []
    #     for idx in self.phymap.index:
    #         start = self.phymap.phy0[idx]
    #         end = self.phymap.phy1[idx]

    #         for sidx in range(self.seqarr.shape[0]):
    #             aseq = self.seqarr[sidx, start:end]
    #             if not np.all(aseq == 78):
    #                 name = self.names[sidx]
    #                 seq = bytes(aseq).decode()
    #                 scaf = self.scaffold_table.scaffold_name[self.scaffold_idx]
    #                 header = ">{} | {}:{}-{}".format(name, scaf, start, end)
    #                 locs.append("{}\n{}".format(header, seq))

    #     # write to file
    #     with open(self.outfile, 'w') as out:
    #         out.write("\n\n".join(locs))        


    def write_to_phy(self):

        # build phy
        phy = []
        for idx, name in enumerate(self._pnames):
            seq = bytes(self.seqarr[idx]).decode()
            phy.append("{} {}".format(name, seq))
       
        # write to temp file
        ntaxa, nsites = self.seqarr.shape
        with open(self.outfile, 'w') as out:
            out.write("{} {}\n".format(ntaxa, nsites))
            out.write("\n".join(phy))
