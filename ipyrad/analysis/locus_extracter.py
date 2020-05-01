#!/usr/bin/env python

"""
Extract, filter, count SNPs and reorganize filtered loci into a 3-d array
(with maxlen limit) for simple indexing. 
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
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError, GETCONS
from ..assemble.write_outputs import NEXHEADER



class LocusExtracter(object):
    """
    Extract and filter a sequence alignment from a genomic window. You can 
    subselect or filter to include certain samples, and to require that data
    are not missing across some number of them, or to relabel them.

    This tool is primarily for INTERNAL use by ipyrad. Users are more likely
    to want to access the window_extracter tool, which provides similar 
    functionality in addition to more options.  

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
        workdir="analysis-locus_extracter",
        mincov=4,
        exclude=None,
        imap=None,
        minmap=None,
        minsnps=0,
        maxmissing=1.0,
        minlen=50,
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
        self.exclude = (exclude if exclude else [])
        self.mincov = mincov
        self.imap = imap
        self.minmap = minmap
        self.minsnps = minsnps
        self.consensus_reduce = consensus_reduce
        self.maxmissing = maxmissing
        self.minlen = minlen
        self.quiet = quiet

        # hardcoded in locus extracter
        self.rmincov = 0.1

        # minmap defaults to 0 if empty and imap
        if self.imap:
            if self.minmap is None:
                self.minmap = {i: 0 for i in self.imap}

        # file to write to
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # require imap for consensus
        if consensus_reduce and (imap is None):
            raise IPyradError(
                "consensus_reduce option requires an imap dictionary")

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

        # global values
        self.filters = {"mincov": 0, "minmap": 0, "minsnps": 0, "maxmissing": 0}
        self._name = name
        self.loci = []
        self.scaffold_table = None
        self.names = []
        self.pnames = []

        # per scaff values that are updated each time.
        self._scaffold_idx = None
        self._names = []
        self._pnames = []

        # get sample names, sidxs, and scaf names        
        self._parse_scaffolds_meta()

        # get start and stop positions in phymap
        self._load_phymap()

        # set minmap and mincov to ints based on nsamples/imap
        self._set_filters_type()

        # records idx, len, nsnps, nsites for each kept locus.
        self._init_stats()



    def run(self, ipyclient=None, force=False, show_cluster=False, auto=False):
        """
        Distribute tree slider jobs in parallel. 

        Parameters:
        -----------
        ipyclient: (type=ipyparallel.Client); Default=None. 
            If you started an ipyclient manually then you can 
            connect to it and use it to distribute jobs here.

        force: (type=bool); Default=False.
            Force overwrite of existing output with the same name.

        show_cluster: (type=bool); Default=False.
            Print information about parallel connection.

        auto: (type=bool); Default=False.
            Let ipyrad automatically manage ipcluster start and shutdown. 
            This will connect to all avaiable cores by default, but can 
            be modified by changing the parameters of the .ipcluster dict 
            associated with this tool.
        """

        # wrap analysis in parallel client
        pool = Parallel(
            tool=self, 
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force},
            )
        pool.wrap_run()



    def _run(self, force=False, ipyclient=None):
        """
        Parallelize locus filtering
        """
        # load balance parallel jobs 2-threaded
        lbview = ipyclient.load_balanced_view()

        # print Nloci to start.
        print("[locus filter] full data: {}".format(self.phymap.shape[0]))

        # submit jobs in blocks of 2K loci
        chunksize = 2000
        rasyncs = {}
        for lidx in range(0, self.phymap.shape[0], chunksize):

            # select block
            block = self.phymap.iloc[lidx:lidx + chunksize, :]

            # submit job
            args = (self, block)
            rasyncs[lidx] = lbview.apply(remote_filter_loci, *args)

        # wait for jobs to finish
        ipyclient.wait()

        # collect results as an ordered list
        chunk = 0
        self.loci = []
        self.smask = []
        for lidx in sorted(rasyncs.keys()):
            loci, filters, farr, smask = rasyncs[lidx].get()

            # store loci
            self.loci.extend(loci)

            # store filtermap
            self.phymap.iloc[chunk: chunk + chunksize, -1] = farr

            # store filters sums
            for key in self.filters:
                self.filters[key] += filters[key]

            # store retrieval mask
            self.smask.append(smask)

            # advance counter
            chunk += chunksize

        # update stats from smask
        self.smask = np.concatenate(self.smask)
        self.sample_stats.loci = self.smask.sum(axis=0)

        # print results
        print("[locus filter] post filter: {}".format(len(self.loci)))



    def _load_phymap(self):
        # load locidx, seqstart, seqend, genomestart, genomeend
        with h5py.File(self.data, 'r') as io5:
            colnames = io5["phymap"].attrs["columns"]
            self.phymap = pd.DataFrame(
                data=io5["phymap"][:],
                columns=[i.decode() for i in colnames],
            )
        self.phymap.loc[:, "filtered"] = False



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
                    "mincov ({}) cannot be greater than number of samples\n "
                    "after consensus_reduce ({}). This will filter all data."
                    ).format(self._mincov, nsamples)
                )

            # local minmap applies before 
            if self.minmap:
                if max(self.minmap.values()) > 1:
                    raise IPyradError(
                        "minmap int value cannot be >1 when using "
                        "consensus_reduce or all data will be filtered."
                    )

        # re-set population filters as integers
        if self.minmap and self.imap:
            for ikey, ivals in self.imap.items():

                # get int value entered by user
                imincov = self.minmap[ikey]

                # get minmap as an int 
                if self.consensus_reduce:
                    self._minmap[ikey] = (
                        imincov if isinstance(imincov, int) 
                        else int(imincov * len(ivals))
                    )
                else:
                    self._minmap[ikey] = (
                        imincov if isinstance(imincov, int) 
                        else int(imincov * 1.0)
                    )



    def _init_stats(self):

        self.sample_stats = pd.DataFrame(
            data=[0] * len(self.names),
            index=self.names, 
            columns=["loci"],
        )
        # stats table
        # scaf = self.scaffold_table.loc[self._scaffold_idx, "scaffold_name"]
        # self.stats = pd.DataFrame({
        #     "scaffold": [scaf],
        #     "start": [self.start],
        #     "end": [self.end],
        #     "sites": [0],
        #     "snps": [0],
        #     "missing": [0],
        #     "samples": [0],
        # }, index=["prefilter", "postfilter"],
        # )

    # def _calc_initial_stats(self):
    #     # enter initial stats
    #     self.stats.loc["prefilter", "sites"] = self.seqarr.shape[1]
    #     self.stats.loc["prefilter", "snps"] = count_snps(self.seqarr)
    #     self.stats.loc["prefilter", "missing"] = round(
    #         np.sum(self.seqarr == 78) / self.seqarr.size, 2)
    #     self.stats.loc["prefilter", "samples"] = self.seqarr.shape[0]




    def _imap_consensus_reduce(self):
        """
        Called on REMOTE.
        Get consensus sequence for all samples in clade. 
        And resets .names, ._pnames as alphanumeric imap keys for writing.
        """
        # skip if no imap
        if not self.consensus_reduce:
            return 

        # empty array of shape imap groups
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

        # save as new data --------- (TODO HERE) -------------
        iarr[iarr == 0] = 78
        self.con_seqarr = iarr



    def _filter_seqarr(self):
        """
        Called on REMOTE.
        Apply filters to remove sites from alignment and to drop taxa if 
        they end up having all Ns.
        """
        # track filters
        filters = {
            "mincov": 0, "minmap": 0, "minsnps": 0, 
            "maxmissing": 0, "minlen": 0,
        }

        # convert dash (-) to Ns
        self.seqarr[self.seqarr == 45] = 78

        # drop SITES that are too many Ns in minmap pops
        if self.imap:
            imapdrop = np.zeros(self.seqarr.shape[1], dtype=bool)
            for ikey, ivals in self.imap.items():

                # imap drops sites if minmap is below nsamples in group
                match = [np.where(self.names == i)[0] for i in ivals]
                sidxs = [i[0] for i in match if i.size]
                subarr = self.seqarr[sidxs, :]
                imapdrop += np.sum(subarr != 78, axis=0) < self._minmap[ikey]

        # replace data with consensus reduced to apply filters
        if self.consensus_reduce:
            self.seqarr = self.con_seqarr

        # drop SITES that don't meet imap/minmap filter
        if self.imap:        
            self.seqarr = self.seqarr[:, np.invert(imapdrop)]

        # drop SITES that don't meet mincov filter (applied after cons_reduc)
        mincovdrop = np.sum(self.seqarr != 78, axis=0) < self._mincov
        self.seqarr = self.seqarr[:, np.invert(mincovdrop)]

        ### TODO: apply ordered filters

        # keep sample rows: only relevant if locus passes filters.
        keep = np.ones(self.seqarr.shape[0]).astype(bool)

        # all sites were filtered?
        if not self.seqarr.size:
            filters["mincov"] = 1
            return filters, keep

        # filtered to be too short
        if self.seqarr.shape[1] < self.minlen:
            filters["minlen"] = 1
            return filters, keep

        # drop LOCUS if too many Ns globally AFTER site filtering.
        missprop = np.sum(self.seqarr == 78) / self.seqarr.size
        if missprop > self.maxmissing:
            filters["maxmissing"] = 1
            return filters, keep

        # drop LOCUS if too few SNPs
        nsnps = count_snps(self.seqarr)
        if nsnps < self.minsnps:
            filters["minsnps"] = 1
            return filters, keep

        # mark SAMPLES that are only Ns after removing lowcov sites
        # for masking by filling the sample_filter array.
        rcovp = np.sum(self.seqarr != 78, axis=1) / self.seqarr.shape[1]
        keep = rcovp >= self.rmincov
        return filters, keep


    # def _calc_filtered_stats(self):
    #     # update stats
    #     self.stats.loc["postfilter", "sites"] = self.seqarr.shape[1]
    #     self.stats.loc["postfilter", "snps"] = count_snps(self.seqarr)
    #     self.stats.loc["postfilter", "missing"] = round(
    #         np.sum(self.seqarr == 78) / self.seqarr.size, 2)
    #     self.stats.loc["postfilter", "samples"] = self.seqarr.shape[0]


    def _parse_scaffolds_meta(self):
        """
        Called at the beginning of ._single_prep()
        Get chromosome lengths for All from the database and get names.
        """
        with h5py.File(self.data, 'r') as io5:

            # get sample names
            self.pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
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

            # get wpnames for cons
            if self.consensus_reduce:
                self.wnames = sorted(self.imap.keys())
                self._longname = 1 + max([len(i) for i in self.wnames])
                self.wpnames = np.array([
                    "{}{}".format(name, " " * (self._longname - len(name)))
                    for name in self.wnames
                ])
            else:
                self.wnames = self.names
                self.wpnames = self.pnames


    def get_shape(self, lidx, include_empty_rows=False):
        # build phy
        seqarr = self.loci[lidx]
        if not include_empty_rows:
            keep = self.smask[lidx]
            seqarr = seqarr[keep, :]
        return seqarr.shape


    def get_locus(self, lidx, include_empty_rows=False, include_names=True):
        """
        Writes the .seqarr matrix as a string to .outfile.
        """
        # build phy
        seqarr = self.loci[lidx]
        pnames = self.wpnames

        # drop SAMPLES that are only Ns after removing lowcov sites
        if not include_empty_rows:
            keep = self.smask[lidx]
            seqarr = seqarr[keep, :]
            pnames = self.wpnames[keep]

        # return array without names
        if not include_names:
            return seqarr

        # convert to bytes and write spacer names
        phy = []
        for idx, name in enumerate(pnames):
            seq = bytes(seqarr[idx]).decode()
            phy.append("{} {}".format(name, seq))
        return phy


    def get_locus_phy(self, lidx):
        # write to string       
        phy = self.get_locus(lidx)
        ntaxa, nsites = self.seqarr.shape
        stringout = "{} {}\n{}".format(ntaxa, nsites, "\n".join(phy))
        return stringout


    def _write_to_nex(self):

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



def remote_filter_loci(self, block):

    # new fresh list
    self.loci = []
    self.filters = {"mincov": 0, "minmap": 0, "minsnps": 0, "maxmissing": 0}
    self.farr = np.zeros(block.shape[0], dtype=np.bool)
    self.smask = []
    # self.smask = np.zeros((block.shape[0], len(self.names)), dtype=np.bool)

    # open h5 for reading seqarray
    with h5py.File(self.data, 'r') as io5:

        # iterate over the locus index in each block
        for idx, lidx in enumerate(block.index):

            # get sampling positions
            self.start = self.phymap.loc[lidx, "phy0"]
            self.end = self.phymap.loc[lidx, "phy1"]

            # extract sequence
            self.seqarr = io5["phy"][self.sidxs, self.start:self.end]

            # applies SITE filters to seqarr and returns loc and row filters
            self._imap_consensus_reduce()
            filters, smask = self._filter_seqarr()

            # only keep if it passed filtering
            if not any(list(filters.values())):
                self.loci.append(self.seqarr)
                self.smask.append(smask)
            else:
                self.farr[idx] = True

            # update stats
            for key in self.filters:
                self.filters[key] += filters[key]

    return self.loci, self.filters, self.farr, self.smask
