#!/usr/bin/env python

"""
Extract, filter and format a locus from a scaffold, start, and end 
tuple. Intended for extracting concatenated sequence windows.
"""

# py2/3 compat
from __future__ import print_function

# standard lib
import os
import h5py
import itertools
import numpy as np
import pandas as pd
from copy import copy
from numba import njit

from .utils import count_snps
from ipyrad.assemble.utils import GETCONS, IPyradError
from ipyrad.assemble.write_outputs import NEXHEADER


"""
    TODO:
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
    scaffold_idxs (int, list, or range):
        Select scaffold(s) by index number. If unsure, leave this
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
    minmap (dict):
        A dictionary mapping group names (keys) to integers or floats to act
        as a filter requiring that at least N (or N%) of samples in this 
        group have data for a locus to be retained in the dataset. When using
        consensus_reduce=True the imap applies to the reduced data set, i.e.,
        it applies to the groups (keys) so that all values must be <= 1.
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
        scaffold_idxs=None, 
        start=None, 
        end=None, 
        mincov=4,
        rmincov=0,
        exclude=None,
        imap=None,
        minmap=None,
        consensus_reduce=False,
        quiet=False,
        **kwargs
        ):
        # report bad arguments
        if kwargs:
            print(
                "Warning: Some arg names are not recognized and may have "
                "changed. Please check the documentation:\n"
                "{}".format(kwargs))

        # store params
        self.data = data
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.scaffold_idxs = scaffold_idxs
        self.start = (int(start) if start is not None else 0)
        self.end = (int(end) if end is not None else None)
        self.exclude = (exclude if exclude else [])
        self.mincov = mincov
        self.rmincov = float(rmincov if rmincov else 0.0)
        self.imap = imap
        self.minmap = minmap
        self.consensus_reduce = consensus_reduce
        self.quiet = quiet

        # rmincov must be float
        assert rmincov < 1.0, "rmincov must be a float between 0 and 1."

        # file to write to
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # require imap for consensus
        if self.consensus_reduce:

            # must tell it how to reduce samples
            if self.imap is None:
                raise IPyradError(
                    "consensus_reduce option requires an imap dictionary")

        # imap requires a minmap so default set to 0
        if self.imap:
            if self.minmap is None:
                self.minmap = {i: 0 for i in self.imap}

        # global values
        self._name = name
        self._filters_checked = False
        self.scaffold_table = None
        self.names = []
        self.pnames = []

        # per scaff values that are updated each time.
        self._scaffold_idx = None
        self._names = []
        self._pnames = []
        self.phymap = []

        # gets .names, .pnames, .scaffold_table from the database file.
        self._parse_scaffolds()

        # if iterable of size 1 then extract single value. [1] -> 1
        if isinstance(self.scaffold_idxs, (list, tuple, np.ndarray, range)):
            if len(self.scaffold_idxs) == 1:
                self.scaffold_idxs = int(self.scaffold_idxs[0])

        # None for debugging and examining .names and .scaffold_table.
        if self.scaffold_idxs is None:
            self._scaffold_idx = None

        # run a single locus: single prep
        elif isinstance(self.scaffold_idxs, (int, str)):
            self._scaffold_idx = self.scaffold_idxs
            self._single_prep()
            self.phymap = self._phymap
            self.pnames = self._pnames

        # concatenate many loci: runs single_prep for each scaffold in list
        elif isinstance(self.scaffold_idxs, (list, tuple, np.ndarray, range)):

            # suppress messages
            self.quiet = True

            # set generic name
            self.name = 'concat'  # (name if name else "concatenated")

            # limit range input to scaffold_table
            if isinstance(scaffold_idxs, range):
                self.scaffold_idxs = range(
                    self.scaffold_table.index.start,
                    min([
                        self.scaffold_table.index.stop,
                        self.scaffold_idxs.stop,
                    ]),
                    self.scaffold_idxs.step,
                )

            # store chunks
            blocks = []
            names = []
            pnames = []
            stats = []
            phymaps = []

            # iterate over scaffs 
            for scaff in self.scaffold_idxs:
                self._scaffold_idx = scaff
                self.start = self.end = None
                self._single_prep()

                # collect data from this scaff/loc
                if self.seqarr.size:
                    blocks.append(self.seqarr)
                    names.append(self._names)
                    pnames.append(self._pnames)
                    stats.append(self.stats)
                    phymaps.append(self._phymap)

            # if no data passed filtering for any loci then bail out
            if not stats:
                print("no data passed filtering")
                return 

            # concat chunks and stats
            self.phymap = pd.concat(phymaps, ignore_index=True)
            self._stats = pd.concat(stats)
            self.names = sorted(set(itertools.chain(*names)))
            self.pnames = sorted(set(itertools.chain(*pnames)))            

            # special order: if 'reference' it goes first
            if "reference" in self.names:
                ridx = self.names.index("reference")
                self.names.pop(ridx)
                self.pnames.pop(ridx)
                self.names = ["reference"] + self.names
                self.pnames = ["reference"] + self.pnames

            # TODO: print warning if some samples were dropped from imap.
            # ...

            # fill concat seqarr allowing for missing taxa in each block
            nsites = int(self._stats.sites.postfilter.sum())
            self.seqarr = np.zeros((len(self.names), nsites), dtype=np.uint8)
            idx = 0           
            for (inames, iblock) in zip(names, blocks):
                # get index of each name in this block
                for pdx in range(len(inames)):
                    # get string name for index in this block
                    iname = inames[pdx]
                    # get index of this string in concat name list order
                    nidx = self.names.index(iname)
                    # fill concat seqarr at this index
                    self.seqarr[nidx, idx:idx + iblock.shape[1]] = iblock[pdx]
                idx += iblock.shape[1]

            # fill remaining missing with Ns and calc missing on concat
            self.seqarr[self.seqarr == 0] = 78
            missing = np.sum(self.seqarr == 78) / self.seqarr.size

            # recalc stats as a concat block
            totals = pd.DataFrame({
                "scaffold": ["concatenated"], 
                "start": [0], 
                "end": [nsites],
                "sites": [nsites],
                "snps": [self._stats.snps.postfilter.sum()],
                "missing": [missing],
                "samples": [len(self.names)],
            })
            self.stats = totals
            self.quiet = 0  # quiet




    def _single_prep(self):
        """
        This load the full data for all samples from the database and stores
        to .names, .pnames, .scaffold_table. 

        Then it subselects ._scaffold_idx and 
        Applies to a single scaffold/locus.
        """
        # output prefix name, e.g., scaff0-0-1000
        self._get_name(self._name)

        # convert minmap to ints based on current imap.
        # only needs to be done once unless consensus reduce, then always.
        if self.consensus_reduce or (not self._filters_checked):
            self._set_filters_type()           

        # stats is overwritten if fillseqar runs
        self.stats = "No stats because no scaffolds selected."

        # pull the region meta-data from the database
        self._extract_phymap_df()
        self._init_stats()

        # get seqs in the selected region and for all taxa in imap
        self._extract_seqarr()
        if not self.seqarr.size:
            self._print("No data in selected window.")
            return 

        # get stats on window (ntaxa missing; nsnps, ntrimmed sites)
        self._calc_initial_stats()

        # apply optional consensus reduction
        self._new_imap_consensus_reduce()

        # apply filters to seqarr, updates names, pnames, etc.
        self._new_filter_seqarr()

        # recalculate stats
        if not self.seqarr.size:
            self._print("No data in filtered window.")
        else:
            self._calc_filtered_stats()

        # reset ending
        self.end = None


    def _parse_scaffolds(self):
        """
        Called at the beginning of ._single_prep()
        Get chromosome lengths for All from the database and get names.
        """
        with h5py.File(self.data, 'r') as io5:

            # get sample names
            try:
                self.pnames = np.array([
                    i.decode() for i in io5["phymap"].attrs["phynames"]
                ])
            except AttributeError:
                # If "names" aren't encoded as bytes then it's an older version
                # of the snps.hdf5 file, so allow for this.
                self.pnames = [i for i in io5["phymap"].attrs["phynames"]]

            self.allnames = [i.strip() for i in self.pnames]

            # auto-generate exclude from imap difference
            if self.imap:
                imapset = set(itertools.chain(*self.imap.values()))
                self.exclude = set(self.allnames).difference(imapset)

            # filter to only the included samples
            self.sidxs = [
                i for (i, j) in enumerate(self.allnames) if j not in self.exclude]
            self.names = np.array([
                j for (i, j) in enumerate(self.allnames) if i in self.sidxs])
            self.pnames = self.pnames[self.sidxs]

            # format names to include spacer for phylip, etc.
            self._longname = 1 + max([len(i) for i in self.pnames])
            self.pnames = np.array([
                "{}{}".format(name, " " * (self._longname - len(name)))
                for name in self.pnames
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


    def _get_name(self, name):
        """
        sets output prefix name. If 
        """
        # use provided name else auto gen a name (scaff-start-end)
        if not name:
            start = (int(self.start) if self.start else "")
            end = (int(self.end) if self.end else "")
            if isinstance(self._scaffold_idx, int):
                self.name = "scaf{}-{}-{}".format(
                    self._scaffold_idx,
                    start,
                    end
                )
            else:
                self.name = "r{}".format(np.random.randint(0, 1e9))
        else:
            self.name = name


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
        self._minmap = copy(self.minmap)
        self._mincov = (
            self.mincov if isinstance(self.mincov, int) 
            else int(self.mincov * nsamples)
        )

        # warning checks for user
        if self.consensus_reduce:

            # global mincov cannot be greater than imap
            if self._mincov > nsamples:
                raise IPyradError((
                    "mincov ({}) cannot be greater than the number of samples "
                    "after consensus_reduce (the number of imap keys: {}) "
                    "or it will filter all data."
                    ).format(self._mincov, nsamples)
                )

        # re-set population filters as integers
        if self.minmap and self.imap:
            for ikey, ivals in self.imap.items():

                # get int value entered by user
                imincov = self.minmap[ikey]
                self._minmap[ikey] = (
                    imincov if isinstance(imincov, int) 
                    else int(imincov * len(self.imap[ikey]))
                )

        # record that this has been done
        self._filters_checked = True


    def _extract_phymap_df(self):
        """
        This extracts the phymap DataFrame.
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table. 
        """
        with h5py.File(self.data, 'r') as io5:
            colnames = io5["phymap"].attrs["columns"]

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == self._scaffold_idx + 1

            # load dataframe of this scaffold
            self._phymap = pd.DataFrame(
                data=io5["phymap"][:][mask],
                columns=[i.decode() for i in colnames],
            )


    def _init_stats(self):
        # stats table
        scaf = self.scaffold_table.loc[self._scaffold_idx, "scaffold_name"]
        self.stats = pd.DataFrame({
            "scaffold": [scaf],
            "start": [self.start],
            "end": [self.end],
            "sites": [0],
            "snps": [0],
            "missing": [0],
            "samples": [0],
        }, index=["prefilter", "postfilter"],
        )


    def _extract_seqarr(self):
        """
        Extracts seqarr of the full selected window and computes stats on
        the window. This initially includes all taxa in the imap. 
        The phymap contains all selected scaffolds and is used to project 
        the selected positions e.g., 100-200 into their column indices in the
        phy alignment e.g., 0-100. 

        If we selected that we wanted scaffolds 0-20 then this will extract
        sequence data to concatenate from all of these scaffolds.

        If we selected a 200bp window on scaffold 0 then it will extract only
        that subset window from one scaffold.

        # Example phymap where we want to select scaf=0, start=5500, end=20000
                chroms  phy0    phy1    pos0    pos1
            0   1         0      66      3942    4008
            1   1         66     152     5318    5404        <-- start
            2   1         152    216     5634    5698
            3   1         216    266     6236    6286
            4   1         266    338     20667   20739       <-- end
            5   2         400    500     20800   20900
        """

        # get mask to select window array region. No mask for denovo
        self.seqarr = np.zeros((len(self.names), 0))

        # select a single scaffold to pull start-stop alignment from.        
        if self.end:

            # get the chrom/loc that includes pos=start]
            mask1 = self._phymap.pos1 >= self.start
            mask2 = self._phymap.pos0 <= self.end
            mask = mask1 & mask2
            block = self._phymap.loc[mask, :]

            # bail out if block is empty
            if not block.size:
                self.seqarr = np.array([])
                return 

            # get start pos as phy position (how far past pos0 is start)
            # wmin_offset = self.start - block.iloc[0, 3]
            # wmin = int(block.iloc[0, 1] + wmin_offset)
            wmin_offset = max(0, self.start - block.iloc[0, 3])
            wmin = int(block.iloc[0, 1] + wmin_offset)

            wmax_offset = int(max(0, block.iloc[-1, -1] - self.end))
            wmax = int(block.iloc[-1, 2] - wmax_offset)

            # extract sequences
            with h5py.File(self.data, 'r') as io5:
                self.seqarr = io5["phy"][self.sidxs, wmin:wmax]

        # if no end argument then select the entire scaffold(s)/locus
        else:

            # if no hits to this scaffold then skip it
            if self._phymap.size:

                # phy start and end of selected chroms
                wmin = self._phymap.iloc[:, 1].min()
                wmax = self._phymap.iloc[:, 2].max()

                # extract array from window
                with h5py.File(self.data, 'r') as io5:
                    self.seqarr = io5["phy"][self.sidxs, wmin:wmax]


    def _calc_initial_stats(self):
        # enter initial stats
        self.stats.loc["prefilter", "sites"] = self.seqarr.shape[1]
        self.stats.loc["prefilter", "snps"] = count_snps(self.seqarr)
        self.stats.loc["prefilter", "missing"] = round(
            np.sum(self.seqarr == 78) / self.seqarr.size, 2)
        self.stats.loc["prefilter", "samples"] = self.seqarr.shape[0]


    def _new_imap_consensus_reduce(self):

        if not self.consensus_reduce:
            return

        iarr = np.zeros((len(self.imap), self.seqarr.shape[1]), dtype=np.uint8)
        inames = np.array(sorted(self.imap.keys()))

        # iterate over imap groups
        for ikey, ivals in self.imap.items():

            # get subarray for this group
            match = [np.where(self.names == i)[0] for i in ivals]
            sidxs = [i[0] for i in match if i.size]        
            subarr = self.seqarr[sidxs, :]

            # get consensus sequence
            cons = consens_sample(subarr, GETCONS)

            # insert to iarr
            iidx = np.where(inames == ikey)[0][0]
            iarr[iidx] = cons

        # save as new data
        iarr[iarr == 0] = 78
        self.con_seqarr = iarr
        self.con_snames = inames
        _longname = 1 + max([len(i) for i in self.con_snames])
        self.con_pnames = np.array([
            "{}{}".format(name, " " * (_longname - len(name)))
            for name in self.con_snames
        ])


    def _new_filter_seqarr(self):
        """
        Apply row and column filters to the seqarr based on imap/minmap, 
        mincov and rmincov (_mincov, _minmap, ...)
        """
        # convert dash (-) to Ns
        self.seqarr[self.seqarr == 45] = 78

        # GET SITES that are too many Ns in minmap pops
        if self.imap:
            imapdrop = np.zeros(self.seqarr.shape[1], dtype=bool)
            for ikey, ivals in self.imap.items():

                # imap drops sites if mincov is below nsamples in group
                match = [np.where(self.names == i)[0] for i in ivals]
                sidxs = [i[0] for i in match if i.size]
                subarr = self.seqarr[sidxs, :]
                imapdrop += np.sum(subarr != 78, axis=0) < self._minmap[ikey]

        # replace data with consensus reduced to apply filters
        if self.consensus_reduce:
            self.seqarr = self.con_seqarr
            self.wnames = self.con_snames
            self.wpnames = self.con_pnames
        else:
            self.wnames = self.names.copy()
            self.wpnames = self.pnames.copy()

        # drop SITES that don't meet imap/minmap filter
        if self.imap:
            self.seqarr = self.seqarr[:, np.invert(imapdrop)]

        # drop SITES that don't meet mincov filter (applied after cons_reduc)
        mincovdrop = np.sum(self.seqarr != 78, axis=0) < self._mincov
        self.seqarr = self.seqarr[:, np.invert(mincovdrop)]

        # drop SAMPLES that are only Ns after removing lowcov sites
        if self.seqarr.size:
            rcovp = np.sum(self.seqarr != 78, axis=1) / self.seqarr.shape[1]
            keep = rcovp >= self.rmincov
            self.seqarr = self.seqarr[keep, :]
            self._names = self.wnames[keep]
            self._pnames = self.wpnames[keep]


    def run(self, force=False, nexus=False):
        """
        Write sequence alignment to a file 
        """
        # bail if user never selected a window.
        if self._scaffold_idx is None:
            return "No scaffold selected"

        # bail if seqarr is empty
        if not self.seqarr.size:
            raise IPyradError("No data in selected region.")

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


    def _print(self, message):
        if not self.quiet:
            print(message)


    def _calc_filtered_stats(self):
        """
        Fill stats df with seqarr info.
        """
        self.stats.loc["postfilter", "sites"] = self.seqarr.shape[1]
        self.stats.loc["postfilter", "snps"] = count_snps(self.seqarr)
        self.stats.loc["postfilter", "missing"] = round(
            np.sum(self.seqarr == 78) / self.seqarr.size, 2)
        self.stats.loc["postfilter", "samples"] = self.seqarr.shape[0]


    def _check_window(self):
        """
        Check that scaf start and end is in scaffold_table.
        """
        assert self.end > self.start, "end must be > start"
        scaffend = self.scaffold_table.scaffold_length[self._scaffold_idx]
        assert self.end <= scaffend, ("end is beyond scaffold length")


    def _write_to_phy(self):
        """
        Writes the .seqarr matrix as a string to .outfile.
        """
        # build phy
        phy = []
        # print('TEST', len(self._pnames), len(self.pnames))
        for idx, name in enumerate(self.pnames):
            seq = bytes(self.seqarr[idx]).decode()
            phy.append("{} {}".format(name, seq))

        # write to temp file
        ntaxa, nsites = self.seqarr.shape
        with open(self.outfile, 'w') as out:
            out.write("{} {}\n".format(ntaxa, nsites))
            out.write("\n".join(phy))


    def _write_to_nex(self):
        """
        Writes concatenated alignment to nex format...
        """
        # write the header
        lines = []
        lines.append(
            NEXHEADER.format(self.seqarr.shape[0], self.seqarr.shape[1])
        )

        # grab a big block of data
        # sidx = 0
        for block in range(0, self.seqarr.shape[1], 100):           
            # store interleaved seqs 100 chars with longname+2 before
            stop = min(block + 100, self.seqarr.shape[1])
            for idx, name in enumerate(self.pnames):  

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



@njit()
def consens_sample(iseq, consdict):
    "Returns the most common base at each site in order."

    altrefs = np.zeros(iseq.shape[1], dtype=np.uint8)

    for col in range(iseq.shape[1]):
        # expand colums with ambigs and remove N-
        fcounts = np.zeros(111, dtype=np.int64)
        counts = np.bincount(iseq[:, col])  #, minlength=90)
        fcounts[:counts.shape[0]] = counts

        # set N and - to zero, wish numba supported minlen arg
        fcounts[78] = 0
        fcounts[45] = 0

        # add ambig counts to true bases
        for aidx in range(consdict.shape[0]):
            nbases = fcounts[consdict[aidx, 0]]
            for _ in range(nbases):
                fcounts[consdict[aidx, 1]] += 1
                fcounts[consdict[aidx, 2]] += 1
            fcounts[consdict[aidx, 0]] = 0

        # work-around for numba bug floating point error
        if np.any(fcounts):
            frac = fcounts / fcounts.sum()
            frac[frac > 0] -= 0.00000001

            # now get counts from the modified counts arr
            alt = np.argmax(np.random.multinomial(1, pvals=frac))
            altrefs[col] = alt
    return altrefs














    # def _imap_consensus_reduce(self):
    #     """
    #     Get consensus sequence for all samples in clade. 
    #     And resets .names, ._pnames as alphanumeric imap keys for writing.
    #     """
    #     # skip if no imap
    #     if (not self.imap) or (not self.consensus_reduce):
    #         self.wnames = self.names
    #         self.wpnames = self.pnames
    #         return

    #     # empty array of shape imap groups
    #     iarr = np.zeros((len(self.imap), self.seqarr.shape[1]), dtype=np.uint8)
    #     inames = np.array(sorted(self.imap.keys()))

    #     # iterate over imap groups
    #     for ikey, ivals in self.imap.items():

    #         # TODO: SHOULD MINMAP FILTER APPLY HERE??

    #         # get subarray for this group
    #         match = [np.where(self.names == i)[0] for i in ivals]
    #         sidxs = [i[0] for i in match if i.size]        
    #         subarr = self.seqarr[sidxs, :]

    #         # get consensus sequence
    #         cons = consens_sample(subarr, GETCONS)

    #         # insert to iarr
    #         iidx = np.where(inames == ikey)[0][0]
    #         iarr[iidx] = cons

    #     # save as new data --------- (TODO HERE) -------------
    #     iarr[iarr == 0] = 78
    #     self.seqarr = iarr
    #     self.wnames = inames
    #     self._longname = 1 + max([len(i) for i in self.wnames])
    #     self.wpnames = np.array([
    #         "{}{}".format(name, " " * (self._longname - len(name)))
    #         for name in self.wnames
    #     ])

    # def _filter_seqarr(self):
    #     """
    #     Apply filters to remove sites from alignment and to drop taxa if 
    #     they end up having all Ns.
    #     """
    #     # drop SITES that are too many Ns given (global) mincov
    #     drop = np.sum(self.seqarr != 78, axis=0) < self._mincov

    #     # drop SITES that are too many Ns in minmap pops
    #     if self._minmap and self.imap:
    #         for ikey, ivals in self.imap.items():

    #             # imap drops sites if mincov is below nsamples in group
    #             if not self.consensus_reduce:
    #                 match = [np.where(self.names == i)[0] for i in ivals]
    #                 sidxs = [i[0] for i in match if i.size]
    #                 subarr = self.seqarr[sidxs, :]
    #                 drop += np.sum(subarr != 78, axis=0) < self._minmap[ikey]

    #             # imap could drop sites in consens if minmap is (1,0,1,1,0)
    #             else:
    #                 sidxs = np.where(self.wnames == ikey)[0][0]
    #                 subarr = self.seqarr[sidxs, :]
    #                 drop += np.sum(subarr != 78, axis=0) < self._minmap[ikey]

    #     # apply site filter
    #     keep = np.invert(drop)        
    #     self.seqarr = self.seqarr[:, keep]

    #     # convert dash (-) to Ns
    #     self.seqarr[self.seqarr == 45] = 78

    #     # drop SAMPLES that are only Ns after removing lowcov sites
    #     rcovp = np.sum(self.seqarr != 78, axis=1) / self.seqarr.shape[1]
    #     keep = rcovp >= self.rmincov
    #     # keep = np.invert(np.all(self.seqarr == 78, axis=1))
    #     self.seqarr = self.seqarr[keep, :]
    #     self._names = self.wnames[keep]
    #     self._pnames = self.wpnames[keep]







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
