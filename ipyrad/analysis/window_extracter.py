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
import itertools
import numpy as np
import pandas as pd

from .utils import count_snps
from ipyrad.assemble.utils import GETCONS, IPyradError
from ipyrad.assemble.write_outputs import reftrick, NEXHEADER


"""
    TODO:

    'format' argument in .run() to write nexus format.

    Select scaffold(s) by index number. If unsure, leave this
    empty when loading a file and then check the .scaffold_table to view
    the indices of scaffolds. 

    To concatenate data from multiple scaffolds
    you can enter a list or slice, e.g., [0, 1, 2] or [0:10]. Default=all.

    To concatenate data from multiple scaffolds
    you can enter a list or slice, e.g., [0, 1, 2] or [0:10]. Default=all.
"""


class WindowExtracter(object):
    """
    Extract and filter a sequence alignment from a genomic window. You can 
    subselect or filter to include certain samples, and to require that data
    are not missing across some number of them, or to relabel them.

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
        Select scaffold by index number. If unsure, leave this
        empty when loading a file and then check the .scaffold_table to view
        the indices of scaffolds. 
    start (int):
        The starting position on a scaffold to obtain data from (default=0)
    end (int):
        The final position on a scaffold to obtain data from (default=end).
    mincov (int):
        The minimum number of individuals that must have data at a site for it
        to be included in the concatenated alignment (default=4). 
    exclude (list):
        A list of sample names to exclude from the data set. Samples can also
        be excluded by using an imap dictionary and not including them.
    imap (dict):
        A dictionary mapping group names (keys) to lists of samples names 
        (values) to be included in the analysis. 
    consensus_reduce (bool):
        The samples in imap groups will be reduced to a consensus sequence 
        that randomly samples two alleles based on the frequency of alleles 
        in the group. This can reduce overall missing data in alignments. 
    quiet (bool):
        Do not print progress information.

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
        minmap=None,
        consensus_reduce=False,
        quiet=False,
        **kwargs
        ):

        # store params
        self.data = data
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.scaffold_idx = scaffold_idx  # (int(scaffold_idx) if scaffold_idx else None)
        self.start = (start if start else 0)
        self.end = end
        self.exclude = (exclude if exclude else [])
        self.mincov = mincov
        self.minmap = minmap
        self.imap = imap
        self.consensus_reduce = consensus_reduce
        self.quiet = quiet

        # file to write to
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # parse scaffold info
        self.scaffold_table = None
        self.names = []

        # gets names, pnames, scaffold_table, ...
        self._parse_scaffolds()

        # update end to scaff len if not entered
        if not self.end:
            # get length of this scaffold
            if isinstance(scaffold_idx, int):
                self.end = int(self.scaffold_table.iloc[self.scaffold_idx, -1])
            # multiple scaffolds
            else:
                self.end = None

        # output prefix name
        self._get_name(name)

        # set parameters as ints or floats
        self._set_filters_type()

        # stats is overwritten if fillseqar runs
        self.stats = "No stats because no scaffolds selected."

        # this will not do anything if no scaffolds selected.
        self._fill_seqarr()


    def _set_filters_type(self):
        """
        Set mincov and minmap to ints. This requires converting floats to ints
        by multiplying by the number of samples, and/or by the number that will
        remain after consensus reduction.
        """
        # consensus sampling changes the number of samples
        nsamples = len(self.sidxs)
        if self.consensus_reduce:
            nsamples = len(self.imap)

        # global filter can be int or float
        self.mincov = (
            self.mincov if isinstance(self.mincov, int) 
            else int(self.mincov * nsamples)
        )

        # set population filter
        if self.minmap:
            for ikey, ivals in self.imap.items():
                # get int value entered by user
                mincov = self.minmap[ikey]

                # get minmap as an int 
                if self.consensus_reduce:
                    self.minmap[ikey] = (
                        mincov if isinstance(mincov, int) 
                        else int(mincov * len(ivals))
                    )
                else:
                    self.minmap[ikey] = (
                        mincov if isinstance(mincov, int) 
                        else int(mincov * 1.0)
                    )                    

        # warning checks for user
        if self.consensus_reduce:
            if self.mincov > nsamples:
                raise IPyradError((
                    "mincov ({}) cannot be greater than number of samples\n "
                    "after consensus_reduce ({}). This will filter all data."
                    ).format(self.mincov, nsamples)
                )
            if self.minmap:
                if max(self.minmap.values()) > 1:
                    raise IPyradError(
                        "minmap int value cannot be >1 when using "
                        "consensus_reduce or all data will be filtered."
                    )


    def _fill_seqarr(self):
        # allow user to init with None to get scaffold names.
        if self.scaffold_idx is not None:

            # check requested window is in scaff
            # self._check_window()

            # pull the sequence from the database
            self.phymap = None
            self._extract_phymap_df()
            self._init_stats()

            # get seqs
            self._extract_seqarr()
            if not self.seqarr.size:
                self._print("No data in selected window.")
                return 

            # report stats on window (ntaxa missing; nsnps, ntrimmed sites)
            self._calc_initial_stats()
            self._imap_consensus_reduce()
            self._filter_seqarr()
            if not self.seqarr.size:
                self._print("No data in filtered window.")
            else:
                self._calc_filtered_stats()


    def _print(self, message):
        if not self.quiet:
            print(message)


    def _get_name(self, name):
        """
        sets output prefix name. If 
        """
        # use provided name else auto gen a name (scaff-start-end)
        if not name:
            if isinstance(self.scaffold_idx, int):
                self.name = "scaf{}-{}-{}".format(
                    self.scaffold_idx,
                    self.start,
                    self.end
                )
            else:
                self.name = "r{}".format(np.random.randint(0, 1e9))
        else:
            self.name = name


    def _extract_phymap_df(self):
        """
        This extracts the phymap DataFrame.
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table. 
        """
        with h5py.File(self.data) as io5:
            colnames = io5["phymap"].attrs["columns"]

            # mask to select this scaff
            #mask = np.zeros(io5["phymap"].shape[0], dtype=np.bool_)
            #for scaff in self.scaffold_idx:
            #    mask += io5["phymap"][:, 0] == scaff + 1
            mask = io5["phymap"][:, 0] == self.scaffold_idx + 1

            # load dataframe of this scaffold
            self.phymap = pd.DataFrame(
                data=io5["phymap"][:][mask],
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


    def run(self, force=False, nexus=False):
        """
        Write sequence alignment to a file 
        """
        # bail if user never selected a window.
        if self.scaffold_idx is None:
            return "No scaffold selected"

        # make outfile path name
        if nexus:           
            self.outfile = os.path.join(
                self.workdir, self.name + ".nex")
        else:
            self.outfile = os.path.join(
                self.workdir, self.name + ".phy")

        # check for force/overwrite
        if force or (not os.path.exists(self.outfile)):

            # convert to file format
            if nexus:
                self._write_to_nex()
            else:
                self._write_to_phy()
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
        self.seqarr = np.zeros((len(self.names), 0))
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
        if (not self.imap) or (not self.consensus_reduce):
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
        Apply filters to remove sites from alignment and to drop taxa if 
        they end up having all Ns.
        """
        # drop sites that are too many Ns given (global) mincov
        drop = np.sum(self.seqarr != 78, axis=0) < self.mincov

        # drop sites that are too many Ns in minmap pops
        if self.minmap:
            for ikey, ivals in self.imap.items():

                # imap drops sites if mincov is below nsamples in group
                if not self.consensus_reduce:
                    sidxs = [np.where(self.names == i)[0][0] for i in ivals]
                    subarr = self.seqarr[sidxs, :]
                    drop += np.sum(subarr != 78, axis=0) < self.minmap[ikey]

                # imap could drop sites in consens if minmap is (1,0,1,1,0)
                else:
                    sidxs = np.where(self.names == ikey)[0][0]
                    subarr = self.seqarr[sidxs, :]
                    drop += np.sum(subarr != 78, axis=0) < self.minmap[ikey]

        # apply site filter
        keep = np.invert(drop)        
        self.seqarr = self.seqarr[:, keep]

        # drop samples that are only Ns after removing lowcov sites
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

            # auto-generate exclude from imap difference
            if self.imap:
                imapset = set(itertools.chain(*self.imap.values()))
                self.exclude = set(self.names).difference(imapset)

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


    def _write_to_phy(self):

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


    def _write_to_nex(self):

        # write the header
        lines = []
        lines.append(
            NEXHEADER.format(self.seqarr.shape[0], self.seqarr.shape[1])
        )

        # grab a big block of data
        sidx = 0
        for block in range(0, self.seqarr.shape[1], 100):           
            # store interleaved seqs 100 chars with longname+2 before
            stop = min(block + 100, self.seqarr.shape[1])
            for idx, name in enumerate(self._pnames):  

                # py2/3 compat --> b"TGCGGG..."
                seqdat = self.seqarr[idx, block:stop]
                lines.append(
                    "  {}{}\n".format(
                        name,
                        bytes(seqdat).decode()
                        )
                    )
            lines.append("\n")
        lines.append("  ;\nend;")
        # print intermediate result and clear
        with open(self.outfile, 'w') as out:
            out.write("".join(lines))


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
