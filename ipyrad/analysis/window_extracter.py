#!/usr/bin/env python

"""
Extract, filter and format a locus from a scaffold, start, and end 
tuple. Intended for extracting concatenated sequence windows.
"""

# standard lib
import os
import itertools
from copy import copy
from typing import List, Dict, Union, Optional

import h5py
import numpy as np
import pandas as pd
from numba import njit
from loguru import logger

from ipyrad.analysis.utils import count_snps
from ipyrad.assemble.utils import GETCONS, IPyradError
from ipyrad.assemble.write_outputs_converter import NEXHEADER


# pylint: disable=too-many-branches, too-many-statements, no-self-use

class NoDataInWindowError(IPyradError):
    pass


class BaseWindowExtracter:
    """
    Extract and filter a sequence alignment from a genomic window. You 
    can subselect samples using an imap dictionary, or filter based on
    sample coverage using global (mincov) or by-population (minmap) 
    minimums. 

    Parameters:
    -----------
    data (str):
        A '.seqs.hdf5' database file from ipyrad.
    name (str):
        Prefix name used for outfiles. If None then it is created from
        the scaffold, start and end positions.
    workdir (str):
        Dir for output files. Created if it doesn't exist.
    scaffold_idxs (int, list, or range):
        Subsample scaffold(s) by index number. If unsure, leave this
        empty when loading a file and then check the .scaffold_table 
        to view the indices of scaffolds. Scaffolds are ordered by 
        their order in the reference genome file.
    start (int):
        Scaffold start position to extract data from (default=0)
    end (int):
        Scaffold end position to extract data from (None = find end)
    mincov (int):
        Minimum number of individuals that must have data at a site 
        for it to be included in the concatenated alignment (def=4). 
    exclude (list):
        A list of sample names to exclude from the data set. Samples 
        can also be excluded by using an imap dictionary and not 
        including them.
    imap (dict):
        A dictionary mapping group names (keys) to lists of sample 
        names (values) to be included in the analysis. This can be 
        used for 3 things: (1) to select samples to extract data for;
        (2) to filter based on sample coverage in groups (minmap); 
        or (3) to use consensus_reduce=True to reduce the dataset to a 
        consensus sequence for each group.
    minmap (dict):
        A dictionary mapping group names (keys) to integers or floats 
        to act as a filter requiring that at least N (or N%) of samples
        in this group have data for a locus to be retained in the 
        dataset. When using consensus_reduce=True the minmap applies to 
        the reduced data set, i.e., it applies to the groups (keys) so 
        that all values must be <= 1.
    consensus_reduce (bool):
        The samples in imap groups will be reduced to a consensus 
        sequence that randomly samples two alleles based on the 
        frequency of alleles in the group. This can reduce overall 
        missing data in alignments.
    """
    def __init__(
        self, 
        data: str, 
        name: Optional[str]="extracted",
        workdir: Optional[str]="analysis-window_extracter",
        scaffold_idxs: Optional[List[int]]=None, 
        start: Optional[int]=None, 
        end: Optional[int]=None, 
        mincov: Union[float, int]=4,
        rmincov: Union[float, int]=0,
        exclude: Optional[List[str]]=None,
        imap: Optional[Dict[str, List[str]]]=None,
        minmap: Optional[Dict[str, Union[float,int]]]=None,
        consensus_reduce: bool=False,
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
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.scaffold_idxs = scaffold_idxs
        self.start = (int(start) if start is not None else 0)
        self.end = (int(end) if end is not None else None)
        self.exclude = set(exclude if exclude else [])
        self.mincov = mincov
        self.rmincov = float(rmincov if rmincov else 0.0)
        self.imap = imap
        self.minmap = minmap
        self.consensus_reduce = consensus_reduce

        # rmincov must be float
        assert rmincov < 1.0, "rmincov must be a float between 0 and 1."

        # sorted names (order in the hdf5) mapped to padded names
        # filled by parse_scaffolds(). This will exclude names of 
        # excluded samples, and will use imap names if consensus.
        self.snames: np.ndarray = None
        self.sidxs: np.ndarray = None
        self.pnames: Dict[str, str] = None
        self.scaffold_table: pd.DataFrame = None
        self.phymap: pd.DataFrame = None
        self.stats: pd.DataFrame = None
        self.outfile: str = None

        # fills: snames, sidxs, scaffold_table
        self._get_scaffold_table()


    def _get_scaffold_table(self):
        """
        Called at the beginning of ._single_prep()
        Get chromosome lengths for All from the database and get names.
        """
        with h5py.File(self.data, 'r') as io5:

            # get sample names and get them as padded names
            snames = io5["phymap"].attrs["phynames"][:].astype(str)

            # auto-update exclude from imap difference
            if self.imap is not None:
                imapset = set(itertools.chain(*self.imap.values()))
                self.exclude.update(set(snames).difference(imapset))
                logger.debug(
                    "dropping samples that are either not in the imap dict, "
                    f"or are in the exclude list: {self.exclude}")

            # filter to only the included samples
            self.sidxs = np.array([
                i for (i, j) in enumerate(snames) if j not in self.exclude
            ])
            self.snames = np.array([
                j for (i, j) in enumerate(snames) if i in self.sidxs])

            # parse scaf names and lengths from db
            self.scaffold_table = pd.DataFrame(
                columns=["scaffold_name", "scaffold_length"],
                data={
                    "scaffold_name": io5["scaffold_names"][:].astype(str),
                    "scaffold_length": io5["scaffold_lengths"][:],
                }, 
            )


    def _get_scaffold_idxs(self):
        """
        Scaffold_idxs can be entered by the user as None (all), int, or 
        List[int] and here we get the list of ints for each selected
        scaffold.
        """
        if self.scaffold_idxs is None:
            self.scaffold_idxs = self.scaffold_table.index.tolist()
        elif isinstance(self.scaffold_idxs, (int, str)):
            self.scaffold_idxs = [self.scaffold_idxs]
        elif isinstance(self.scaffold_idxs, (list, tuple, np.ndarray, range)):
            self.scaffold_idxs = list(self.scaffold_idxs)
        else:
            raise IPyradError('scaffold_idxs argument format not recognized.')

        # only allow start and end if working on a single scaffold
        if len(self.scaffold_idxs) > 1:
            if self.start or self.end:
                logger.warning(
                    "setting start and end to None (use whole scaffolds), "
                    "position selectors can only be used when selecting a "
                    "single scaffold idx at a time."
                )
                self.start = None
                self.end = None


    def _write_to_phy(self, seqarr, names):
        """
        Writes the .seqarr matrix as a string to .outfile.
        """
        # get padded names
        longname = max(len(i) for i in names)
        pnames = [i.ljust(longname + 5) for i in names]

        # build phy
        phy = []
        for idx, name in enumerate(pnames):
            seq = b"".join(seqarr[idx].view("S1")).decode()
            phy.append(f"{name} {seq}")

        # write to temp file
        ntaxa, nsites = seqarr.shape
        self.outfile = os.path.join(self.workdir, self.name + ".phy")
        with open(self.outfile, 'w') as out:
            out.write(f"{ntaxa} {nsites}\n")
            out.write("\n".join(phy))
        logger.info(f"wrote seq array ({ntaxa, nsites}) to {self.outfile}")


    def _write_to_nex(self, seqarr, names):
        """
        Writes concatenated alignment to nex format...
        """
        # get padded names
        longname = max(len(i) for i in names)
        pnames = [i.ljust(longname + 5) for i in names]

        # write the header
        lines = []
        lines.append(NEXHEADER.format(seqarr.shape[0], seqarr.shape[1]))

        # grab a big block of data
        # sidx = 0
        for block in range(0, seqarr.shape[1], 100):
            # store interleaved seqs 100 chars with longname+2 before
            stop = min(block + 100, seqarr.shape[1])
            for idx, name in enumerate(pnames):
                seqdata = b"".join(seqarr[idx, block:stop].view("S1")).decode()
                lines.append(f"  {name}{seqdata}\n")
            lines.append("\n")
        lines.append("  ;\nend;")
        # print intermediate result and clear
        ntaxa, nsites = seqarr.shape        
        self.outfile = os.path.join(self.workdir, self.name + ".nex")        
        with open(self.outfile, 'w') as out:
            out.write("".join(lines))
        logger.info(f"wrote seq array ({ntaxa, nsites}) to {self.outfile}")



class WindowExtracter(BaseWindowExtracter):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # imap requires a minmap so default set to 0 if absent
        if self.imap is None:
            self.imap = {i: i for i in self.snames}
        if isinstance(self.minmap, (int, float)):
            self.minmap = {i: self.minmap for i in self.imap}
        if self.minmap is None:
            self.minmap = {i: 0 for i in self.imap}
        assert list(self.imap.keys()) == list(self.minmap.keys()), (
            "Keys in the imap dict must match keys in the minmap dict."
        )
        self._get_scaffold_idxs()
        self._convert_min_filters_to_floats()


    def run(self, nexus=False):
        """
        Apply filtering to selected genomic region, extract sequence 
        data and write to an output file in phylip format (default)
        or nexus (option).
        """
        # ensure outdir exists
        os.makedirs(self.workdir, exist_ok=True)

        # store chunks
        blocks = []
        names = set([])
        fstats = []
        phymaps = []

        # this COULD be parallelized, but we also want to support a serial
        # version since wex itself is sometimes distributed in parallel,
        # such as in treeslider. 
        for scaffidx in self.scaffold_idxs:
            phymap = self._extract_phymap_for_scaffidx(scaffidx)
            seqarr = self._extract_seqarr(phymap)
            if not seqarr.size:
                logger.debug(f"No data selected in scaffold idx={scaffidx}")
                continue
            stats = self._get_initial_stats(scaffidx, seqarr)
            fseqarr, sample_keep = self._filter_seqarr(seqarr)

            # collect data from this scaff/loc
            if fseqarr.size:
                stats = self._get_filtered_stats(stats, fseqarr)
                fstats.append(stats)
                names.update(set(sample_keep))
                blocks.append(fseqarr)
                phymaps.append(phymap)

        # bail if seqarr is empty
        if not blocks:
            raise NoDataInWindowError(
                "No data passed filtering in the selected region.")

        # concat chunks
        self.stats = pd.concat(fstats)
        self.phymap = pd.concat(phymaps, ignore_index=True)
        seqarr = np.concatenate(blocks, axis=1)

        # remove concat sample row if below rmincov
        rcovprops = np.sum(seqarr != 78, axis=1) / seqarr.shape[1]
        row_keep = rcovprops >= self.rmincov
        seqarr = seqarr[row_keep, :]
        names = [i for (i, j) in zip(self.snames, row_keep) if j]

        # make outfile path name
        if nexus:           
            self._write_to_nex(seqarr, names)
        else:
            self._write_to_phy(seqarr, names)


    def _convert_min_filters_to_floats(self):
        """
        Set mincov and minmap to ints. This requires converting floats
        to ints by multiplying by the number of samples.
        """
        # number of samples
        nsamples = len(self.sidxs)

        # global filter can be int or float
        self._minmap = copy(self.minmap)
        self._mincov = (
            self.mincov if isinstance(self.mincov, int) 
            else int(self.mincov * nsamples)
        )

        # re-set population filters as integers
        for ikey in self.imap:
            imincov = self.minmap[ikey]
            self._minmap[ikey] = (
                imincov if isinstance(imincov, int) 
                else int(imincov * len(self.imap[ikey]))
            )


    def _extract_phymap_for_scaffidx(self, scaffidx):
        """
        This extracts the phymap DataFrame.
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table. 
        """
        with h5py.File(self.data, 'r') as io5:
            colnames = io5["phymap"].attrs["columns"].astype(str)

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == scaffidx + 1

            # load dataframe of this scaffold
            phymap = pd.DataFrame(
                data=io5["phymap"][:][mask],
                columns=colnames,
            )
        return phymap


    def _extract_seqarr(self, phymap):
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
        seqarr = np.zeros(shape=(len(self.snames), 0), dtype=np.uint8)

        # select a single scaffold to pull start-stop alignment from.        
        if self.end:

            # get the chrom/loc that includes pos=start]
            mask1 = phymap.pos1 >= self.start
            mask2 = phymap.pos0 <= self.end
            mask = mask1 & mask2
            block = phymap.loc[mask, :]

            # bail out if block is empty
            if not block.size:
                return np.array([])

            # get start pos as phy position (how far past pos0 is start)
            # wmin_offset = self.start - block.iloc[0, 3]
            # wmin = int(block.iloc[0, 1] + wmin_offset)
            wmin_offset = max(0, self.start - block.iloc[0, 3])
            wmin = int(block.iloc[0, 1] + wmin_offset)

            wmax_offset = int(max(0, block.iloc[-1, -1] - self.end))
            wmax = int(block.iloc[-1, 2] - wmax_offset)

            # extract sequences
            with h5py.File(self.data, 'r') as io5:
                seqarr = io5["phy"][self.sidxs, wmin:wmax]

        # if no end argument then select the entire scaffold(s)/locus
        else:
            # if no hits to this scaffold then skip it
            if phymap.size:
                # phy start and end of selected chroms
                wmin = phymap.iloc[:, 1].min()
                wmax = phymap.iloc[:, 2].max()

                # extract array from window
                with h5py.File(self.data, 'r') as io5:
                    seqarr = io5["phy"][self.sidxs, wmin:wmax]
        return seqarr


    def _get_initial_stats(self, scaffidx, seqarr):
        """
        Returns a dataframe with initial stats for scaffidx seqarr  
        """
        scaffname = self.scaffold_table.loc[scaffidx, "scaffold_name"]    
        scafflen = self.scaffold_table.loc[scaffidx, "scaffold_length"]    
        stats = pd.DataFrame(
            index=[scaffidx],
            data={
                "scaff_name": [scaffname],
                "scaff_len": [scafflen],
                "start": [self.start],
                "end": [self.end],
                "sites_pre": seqarr.shape[1],
                "sites_post": 0,
                "snps_pre": count_snps(seqarr),
                "snps_post": 0,
                "missing_pre": round(np.sum(seqarr == 78) / seqarr.size, 2),
                "missing_post": 0.0,
                "samples_pre": seqarr.shape[0],
                "samples_post": 0,
            }, 
        )
        return stats


    def _filter_seqarr(self, seqarr):
        """
        Apply row and column filters to the seqarr based on imap/minmap, 
        mincov and rmincov (_mincov, _minmap, ...)
        """
        # convert dash (-) to Ns
        fseqarr = seqarr.copy()
        fseqarr[fseqarr == 45] = 78

        # GET SITES that are too many Ns in minmap pops
        imapdrop = np.zeros(fseqarr.shape[1], dtype=bool)
        for pop in self.imap:
            # imap drops sites if mincov is below nsamples in group
            samples = self.imap[pop]
            match = [np.where(self.snames == i)[0] for i in samples]
            sidxs = [i[0] for i in match if i.size]
            subarr = fseqarr[sidxs, :]
            imapdrop += np.sum(subarr != 78, axis=0) < self._minmap[pop]

        # drop SITES that don't meet imap/minmap filter
        fseqarr = fseqarr[:, np.invert(imapdrop)]

        # drop SITES that don't meet mincov filter
        mincovdrop = np.sum(fseqarr != 78, axis=0) < self._mincov
        fseqarr = fseqarr[:, np.invert(mincovdrop)]

        # MASK SAMPLES that are only Ns after removing lowcov sites.
        # dropping of samples will/may occur later on the concat array.
        if fseqarr.size:
            rcovprops = np.sum(fseqarr != 78, axis=1) / fseqarr.shape[1]
            row_mask = rcovprops < self.rmincov
            fseqarr[row_mask, :] = 78
            sample_keep = set(self.snames[np.invert(row_mask)])
            return fseqarr, sample_keep
        return np.array([]), set([])


    def _get_filtered_stats(self, stats, fseqarr):
        """
        Fill stats df with post-filtered seqarr info.
        """       
        stats.sites_post = fseqarr.shape[1]
        stats.snps_post = count_snps(fseqarr)
        stats.missing_post = round(np.sum(fseqarr == 78) / fseqarr.size, 2)
        stats.samples_post = fseqarr.shape[0] - np.all(fseqarr == 78, axis=0).sum()
        return stats



class WindowExtracterConsensus(BaseWindowExtracter):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # consensus reduce requires imap
        if self.imap is None:
            raise IPyradError(
                "consensus_reduce option requires an imap dictionary")
        if isinstance(self.minmap, (int, float)):
            self.minmap = {i: self.minmap for i in self.imap}
        if self.minmap is None:
            self.minmap = {i: 0 for i in self.imap}
        assert list(self.imap.keys()) == list(self.minmap.keys()), (
            "Keys in the imap dict must match keys in the minmap dict."
        )

        # imap consensus reduced ordered names (row order of seqarr)
        self.inames = np.array(sorted(self.imap.keys()))

        # gets .names, .pnames, .scaffold_table from the database file.
        self._get_scaffold_idxs()
        self._convert_min_filters_to_floats()


    def _convert_min_filters_to_floats(self):
        """
        Set mincov and minmap to ints. This requires converting floats
        to ints by multiplying by the number of samples.
        """
        # consensus sampling changes the number of samples
        nsamples = len(self.imap)

        # global filter can be int or float
        self._minmap = copy(self.minmap)
        self._mincov = (
            self.mincov if isinstance(self.mincov, int) 
            else int(self.mincov * nsamples)
        )
        # global mincov cannot be greater than imap
        if self._mincov > nsamples:
            raise IPyradError((
                "mincov ({}) cannot be greater than the number of samples "
                "after consensus_reduce (the number of imap keys: {}) "
                "or it will filter all data."
                ).format(self._mincov, nsamples)
            )
        # re-set population filters as integers
        for ikey in self.imap:
            imincov = self.minmap[ikey]
            self._minmap[ikey] = (
                imincov if isinstance(imincov, int) 
                else int(imincov * len(self.imap[ikey]))
            )


    def run(self, nexus=False):
        """
        Iterates over scaffolds and concats results.
        """
        # ensure outdir exists
        os.makedirs(self.workdir, exist_ok=True)

        # store chunks
        blocks = []
        names = set([])
        fstats = []
        phymaps = []

        for scaffidx in self.scaffold_idxs:
            phymap = self._extract_phymap_for_scaffidx(scaffidx)
            seqarr = self._extract_seqarr(phymap)
            if not seqarr.size:
                # logger.debug(f"No data selected in scaffold idx={scaffidx}")
                continue
            stats = self._get_initial_stats(scaffidx, seqarr)

            # rseqarr is nrows=len(imap), 
            rseqarr, isample_keep = self._imap_consensus_reduce(seqarr)

            # collect data from this scaff/loc
            if rseqarr.size:
                stats = self._get_filtered_stats(stats, rseqarr)
                fstats.append(stats)
                names.update(set(isample_keep))
                blocks.append(rseqarr)
                phymaps.append(phymap)

        # bail if seqarr is empty
        if not blocks:
            raise NoDataInWindowError(
                "No data passed filtering in the selected region.")

        # concat chunks
        self.stats = pd.concat(fstats)
        self.phymap = pd.concat(phymaps, ignore_index=True)
        seqarr = np.concatenate(blocks, axis=1)

        # remove concat sample row if below rmincov
        rcovprops = np.sum(seqarr != 78, axis=1) / seqarr.shape[1]
        row_keep = rcovprops >= self.rmincov
        seqarr = seqarr[row_keep, :]
        names = [i for i, j in zip(sorted(self.imap), row_keep) if j]

        # make outfile path name
        if nexus:           
            self._write_to_nex(seqarr, names)
        else:
            self._write_to_phy(seqarr, names)


    def _imap_consensus_reduce(self, seqarr):
        """
        Reduces seqarr to nrows=len(imap) by sampling consensus sites.
        If no data is present for a pop at a site it returns N. The 
        filter to remove sites or rows below mincov is in the next 
        step, not here.
        """
        # pop seqarray to fill
        iarr = np.zeros((len(self.imap), seqarr.shape[1]), dtype=np.uint8)

        # iterate over imap groups
        imapdrop = np.zeros(iarr.shape[1], dtype=bool)
        for pop in self.inames:
            samples = self.imap[pop]
            match = [np.where(self.snames == i)[0] for i in samples]
            sidxs = [i[0] for i in match if i.size]        
            subarr = seqarr[sidxs, :]

            # get consensus sequence
            cons = consens_sample(subarr, GETCONS)

            # store to drop sites with sample cov below min setting
            imapdrop += np.sum(subarr != 78, axis=0) < self._minmap[pop]

            # store to arr
            iidx = np.where(self.inames == pop)[0][0]
            iarr[iidx] = cons

        # set any remaining empty sites to N.
        iarr[iarr == 0] = 78

        # drop SITES that don't meet imap/minmap filter
        iarr = iarr[:, np.invert(imapdrop)]

        # drop SITES that don't meet mincov filter
        mincovdrop = np.sum(iarr != 78, axis=0) < self._mincov
        iarr = iarr[:, np.invert(mincovdrop)]

        # MASK SAMPLES that are only Ns after removing lowcov sites.
        # dropping of samples will/may occur later on the concat array.
        if iarr.size:
            rcovprops = np.sum(iarr != 78, axis=1) / iarr.shape[1]
            row_mask = rcovprops < self.rmincov
            iarr[row_mask, :] = 78
            isample_keep = set(self.inames[np.invert(row_mask)])
            return iarr, isample_keep
        return np.array([]), set([])


    def _extract_phymap_for_scaffidx(self, scaffidx):
        """
        This extracts the phymap DataFrame.
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table. 
        """
        with h5py.File(self.data, 'r') as io5:
            colnames = io5["phymap"].attrs["columns"].astype(str)

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == scaffidx + 1

            # load dataframe of this scaffold
            phymap = pd.DataFrame(
                data=io5["phymap"][:][mask],
                columns=colnames,
            )
        return phymap


    def _extract_seqarr(self, phymap):
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
        seqarr = np.zeros(shape=(len(self.snames), 0), dtype=np.uint8)

        # select a single scaffold to pull start-stop alignment from.        
        if self.end:

            # get the chrom/loc that includes pos=start]
            mask1 = phymap.pos1 >= self.start
            mask2 = phymap.pos0 <= self.end
            mask = mask1 & mask2
            block = phymap.loc[mask, :]

            # bail out if block is empty
            if not block.size:
                return np.array([])

            # get start pos as phy position (how far past pos0 is start)
            # wmin_offset = self.start - block.iloc[0, 3]
            # wmin = int(block.iloc[0, 1] + wmin_offset)
            wmin_offset = max(0, self.start - block.iloc[0, 3])
            wmin = int(block.iloc[0, 1] + wmin_offset)

            wmax_offset = int(max(0, block.iloc[-1, -1] - self.end))
            wmax = int(block.iloc[-1, 2] - wmax_offset)

            # extract sequences
            with h5py.File(self.data, 'r') as io5:
                seqarr = io5["phy"][self.sidxs, wmin:wmax]

        # if no end argument then select the entire scaffold(s)/locus
        else:
            # if no hits to this scaffold then skip it
            if phymap.size:
                # phy start and end of selected chroms
                wmin = phymap.iloc[:, 1].min()
                wmax = phymap.iloc[:, 2].max()

                # extract array from window
                with h5py.File(self.data, 'r') as io5:
                    seqarr = io5["phy"][self.sidxs, wmin:wmax]
        return seqarr


    def _get_initial_stats(self, scaffidx, seqarr):
        """
        Returns a dataframe with initial stats for scaffidx seqarr  
        """
        scaffname = self.scaffold_table.loc[scaffidx, "scaffold_name"]    
        scafflen = self.scaffold_table.loc[scaffidx, "scaffold_length"]    
        stats = pd.DataFrame(
            index=[scaffidx],
            data={
                "scaff_name": [scaffname],
                "scaff_len": [scafflen],
                "start": [self.start],
                "end": [self.end],
                "sites_pre": seqarr.shape[1],
                "sites_post": 0,
                "snps_pre": count_snps(seqarr),
                "snps_post": 0,
                "missing_pre": round(np.sum(seqarr == 78) / seqarr.size, 2),
                "missing_post": 0.0,
                "samples_pre": seqarr.shape[0],
                "samples_post": 0,
            }, 
        )
        return stats


    def _get_filtered_stats(self, stats, fseqarr):
        """
        Fill stats df with post-filtered seqarr info.
        """       
        stats.sites_post = fseqarr.shape[1]
        stats.snps_post = count_snps(fseqarr)
        stats.missing_post = round(np.sum(fseqarr == 78) / fseqarr.size, 2)
        stats.samples_post = fseqarr.shape[0] - np.all(fseqarr == 78, axis=0).sum()
        return stats



@njit()
def consens_sample(iseq, consdict):
    """
    Returns a sampled base at each site, randomly sampled by the 
    frequency of called bases among samples in the imap population
    at that base.
    """
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



def window_extracter(
    data: str, 
    name: Optional[str]="extracted",
    workdir: Optional[str]="analysis-window_extracter",
    scaffold_idxs: Optional[List[int]]=None, 
    start: Optional[int]=None, 
    end: Optional[int]=None, 
    mincov: Union[float, int]=4,
    rmincov: Union[float, int]=0,
    exclude: Optional[List[str]]=None,
    imap: Optional[Dict[str, List[str]]]=None,
    minmap: Optional[Dict[str, Union[float,int]]]=None,
    consensus_reduce: bool=False,
    **kwargs,
    ) -> 'BaseWindowExtracter':
    """
    Extract and filter a sequence alignment from a genomic window. You 
    can subselect samples using an imap dictionary, or filter based on
    sample coverage using global (mincov) or by-population (minmap) 
    minimums. 

    Parameters:
    -----------
    data (str):
        A '.seqs.hdf5' database file from ipyrad.
    name (str):
        Prefix name used for outfiles. If None then it is created from
        the scaffold, start and end positions.
    workdir (str):
        Dir for output files. Created if it doesn't exist.
    scaffold_idxs (int, list, or range):
        Subsample scaffold(s) by index number. If unsure, leave this
        empty when loading a file and then check the .scaffold_table 
        to view the indices of scaffolds. Scaffolds are ordered by 
        their order in the reference genome file.
    start (int):
        Scaffold start position to extract data from (default=0)
    end (int):
        Scaffold end position to extract data from (None = find end)
    mincov (int):
        Minimum number of individuals that must have data at a site 
        for it to be included in the concatenated alignment (def=4). 
    exclude (list):
        A list of sample names to exclude from the data set. Samples 
        can also be excluded by using an imap dictionary and not 
        including them.
    imap (dict):
        A dictionary mapping group names (keys) to lists of sample 
        names (values) to be included in the analysis. This can be 
        used for 3 things: (1) to select samples to extract data for;
        (2) to filter based on sample coverage in groups (minmap); 
        or (3) to use consensus_reduce=True to reduce the dataset to a 
        consensus sequence for each group.
    minmap (dict):
        A dictionary mapping group names (keys) to integers or floats 
        to act as a filter requiring that at least N (or N%) of samples
        in this group have data for a locus to be retained in the 
        dataset. When using consensus_reduce=True the minmap applies to 
        the reduced data set, i.e., it applies to the groups (keys) so 
        that all values must be <= 1.
    consensus_reduce (bool):
        The samples in imap groups will be reduced to a consensus 
        sequence that randomly samples two alleles based on the 
        frequency of alleles in the group. This can reduce overall 
        missing data in alignments.
    """  
    if consensus_reduce:
        tool = WindowExtracterConsensus(
            data=data,
            name=name,
            workdir=workdir,
            scaffold_idxs=scaffold_idxs,
            start=start,
            end=end,
            mincov=mincov,
            rmincov=rmincov,
            exclude=exclude,
            imap=imap,
            minmap=minmap,
            **kwargs,
        )
        return tool
    tool = WindowExtracter(
        data=data,
        name=name,
        workdir=workdir,
        scaffold_idxs=scaffold_idxs,
        start=start,
        end=end,
        mincov=mincov,
        rmincov=rmincov,
        exclude=exclude,
        imap=imap,
        minmap=minmap,
        **kwargs,        
    )
    return tool



if __name__ == "__main__":
    
    import ipyrad.analysis as ipa
    ipa.set_loglevel("INFO")

    DATA = "/home/deren/Documents/ipyrad/sandbox/refdata_outfiles/refdata.seqs.hdf5"

    wex = window_extracter(
        data=DATA,
        imap={
            "A": ['reference'],
            "B1": ['SLH_AL_1016'],
            "B": ['SLH_AL_1000', 'SLH_AL_1015'],
            "C": ['SLH_AL_1013', 'SLH_AL_1014'],
            "D": ['digested'],
        },
        scaffold_idxs=[0, 1],
        start=0,
        end=1000000,
        mincov=4,
        minmap=0,
        rmincov=0.05,
        consensus_reduce=True,
    )
    wex.run()
    print(wex.stats.T)
    print(wex.phymap)
