#!/usr/bin/env python

"""Edge trimming/masking of partially filtered loci.

This class is called during filtering of loci, after the duplicate
and minsamples filters have been applied, but before maxsnps filters,
because we first want to clean up (mask) edges of loci where there
may be poor information. These are masked based on minimum sample
coverages.

Order of operations:
1. Trim to `hackers.trim_loci_min_sites`
2. Trim to remove RE overhang len if detected, or `params.trim_loci`, 
whichever is greater.
"""

from typing import TypeVar, Tuple, List
from dataclasses import dataclass
import numpy as np
from numba import njit
from ipyrad.assemble.utils import comp

Assembly = TypeVar("Assembly")

@dataclass
class EdgeTrimmer:
    """class for trimming loci during filtering.

    Although we used to store trim positions for (R1>, <R1, R2>, <R2)
    we now only do so for (R1>, <R2) for PE data, or (R1>, <R1) for
    SE data. Finding and trimming edges at insert positions proved
    too messy in reference-based assemblies where partially overlapping
    reads can tile to create many possible insert positions.

    Note
    ----
    This runs on a remote engine so any `print` statements will be
    reported to `logger.DEBUG`. The print statements are only 
    executed if `debug=True` in this class.
    """
    data: Assembly
    names: List[str]
    seqs: np.ndarray
    nidxs: List[str]
    debug: bool = False

    exclude_ref: bool = False
    """: If True then the ref seq does not contribute to finding edges."""
    trim1: List[int] = None
    """: 0-indexed pos where (start, end) should be trimmed of lowcov edges."""
    trim2: List[int] = None
    """: 0-indexed pos where (start, end) should be trimmed of RE sites."""
    sumlens: List[int] = None
    """: Sum len of trimmed bases from each end."""
    min_sample_coverage_by_site: int = 0
    """: Minimum sample coverage required at locus edge (used for site cov
    trimming, not just sample cov trimming.)"""

    def __post_init__(self):
        # reference is present and should be excluded from trimming info.
        self.exclude_ref = self.data.is_ref and self.data.hackers.exclude_reference

        # edges (len, -len) to trim off low coverage edges (e.g., (5, -5))
        self.trim1 = [0, 0]
        # trims (len, len) to trim off RE overhangs or user trims.
        self.trim2 = [0, 0]
        # record of the summed len of trimmed bases from each side used
        # later to keep sites aligned with positional metadata (depths).
        self.sumlens = [0, 0]

        # e.g., trim_loci_min_sites=10, but nsamples=4, then this is
        # set to 4. BUT, if the reference is in seqs but should not be
        # included then we further subtract one.
        self.min_sample_coverage_by_site = min(
            self.data.hackers.trim_loci_min_sites,
            self.seqs.shape[0] - int(self.exclude_ref)
        )

    def run(self) -> Tuple[List, List, List, bool]:
        """Runs the steps involved in trimming and masking seqs.

        Returns a boolean where True means locus should be filtered.
        """
        # sets .edges array based on sample coverage, and returns bool filter
        if self._set_edges_from_coverage():
            if self.debug:
                print(f"edge trim filtered by coverage: {self.trim1}")
            return None, None, None, True

        # sets .trimseq by slicing .edges from .seqs
        self._trim_seqs_edges_1()

        # sets .trim from RE overhangs
        self._set_trim_from_re_overhangs()

        # updates .trim from user-entered .trim_loci param.
        self._set_trim_from_params()

        # update .edges for .trim lengths.
        self._trim_seqs_edges_2()

        # check that all data was not filtered.
        if self._check_edges():
            if self.debug:
                print(f"edge trim filtered by overhang: {self.trim2}")            
            return None, None, None, True

        # get subsampled names, seqs, nidxs by removing any samples
        # that no longer have data after trimming, which can arise
        # when data are highly tiled.
        # print(self.trim1, self.trim2)
        over_trimmed = self._remove_empty_samples()
        if self.debug and over_trimmed:
            print(f"edge trim filtered by subsampling. Nsamples={len(self.names)}")
        return self.names, self.seqs, self.nidxs, over_trimmed

    def _set_edges_from_coverage(self) -> bool:
        """Fills self.trim1 with trim lens from each edge based on coverage.

        Finds the farthest left and right position where data is
        present (not N or -) and sample coverage >= minsamp, i.e.,
        it calculates per-site sample coverage taking into account
        missing base calls.
        """
        # get nsamples coverage at every site
        site_sample_covs = np.sum((self.seqs != 78) & (self.seqs != 45), axis=0)

        # get left and right most positions from jit'd functions.
        left = locus_left_trim(site_sample_covs, self.min_sample_coverage_by_site)
        right = locus_right_trim(site_sample_covs, self.min_sample_coverage_by_site)

        # return True if no positions had sufficient coverage and
        # this locus should be filtered, else return False
        if right <= left:
            return True

        # else, convert edges 1 from the pos to a -len from the end.
        # if 0 we replace with None for slicing.
        self.trim1[0] = left
        self.trim1[1] = self.seqs.shape[1] - right
        return False

    def _trim_seqs_edges_1(self):
        """Updates .seqs by trimming slice from .trim1.

        end is a length, so we make negative, or None if it is zero.
        """
        start = self.trim1[0]
        end = -self.trim1[1] if self.trim1[1] else None
        self.seqs = self.seqs[:, start:end]

    def _set_trim_from_re_overhangs(self):
        """Fills self.trim with extra edge sizes that should be trimmed.

        Fuzzy search for restriction overhang sequence to be trimmed.
        This applies to the .trimseq array that has already been
        trimmed based on sample coverage.

        TODO: this func could be jit'd for speed improvements.
        """
        for cutter in self.data.params.restriction_overhang:
            if not cutter:
                continue

            # convert to uint8 array
            forward = np.array(list(cutter.encode()), dtype=np.uint8)
            revcomp = np.array(list(comp(cutter).encode()), dtype=np.uint8)[::-1]

            # compare match over cut size skipping Ns and allow 0.25 diffs
            slx = slice(0, forward.shape[0])
            matching = self.seqs[:, slx] == forward
            mask = np.where(
                (self.seqs[:, slx] != 78) & (self.seqs[:, slx] != 45))
            matchset = matching[mask] # non-masked sites that match the RE
            prop = float(matchset.sum()) / matchset.size
            if prop >= 0.75:
                self.trim2[0] = len(forward)
                # print(f"forward={forward}\nmatching={matching}\nmask={mask}\nmatchset={matchset}\nprop={prop:.2f}")
                continue

            # compare match over cut size skipping Ns and allow 0.25 diffs
            slx = slice(self.seqs.shape[1] - revcomp.shape[0], self.seqs.shape[1])
            matching = self.seqs[:, slx] == revcomp
            mask = np.where(
                (self.seqs[:, slx] != 78) & (self.seqs[:, slx] != 45))
            matchset = matching[mask]
            prop = float(matchset.sum()) / matchset.size
            if prop >= 0.75:
                self.trim2[1] = len(revcomp)
            # print(f"revcomp={revcomp}, matching={matching}, mask={mask}, matchset={matchset}, prop={prop:.2f}")

    def _set_trim_from_params(self):
        """Sets .trim2 values based on self.params.trimloci.

        The number at each position of .trimloci is the number of
        bases to trim from the current edge.
        """
        # FIXME (maybe): params.trim_loci still takes a 4-part argument even
        # though we now only use two args (the first and last). If
        # that gets changed then it needs to be updated here. I decided
        # not to change it yet since we may later try to implement middle
        # trimming/masking again...
        self.trim2[0] = max(self.trim2[0], self.data.params.trim_loci[0])
        self.trim2[1] = max(self.trim2[1], abs(self.data.params.trim_loci[3]))

    def _trim_seqs_edges_2(self):
        """Updates .seqs by trimming slice from .trim2.

        This is the final seqs that will be returned.
        """
        start = self.trim2[0]
        end = -self.trim2[1] if self.trim2[1] else None
        self.seqs = self.seqs[:, start:end]
        # self.edges[1] -= self.trims[1]
        # self.edges[2] += self.trims[2]

    def _check_edges(self) -> bool:
        """Return True if no locus length is left untrimmed"""
        if self.seqs.shape[1] < self.data.params.filter_min_trim_len:
            return True
        return False

    def _remove_empty_samples(self) -> Tuple[List[str], np.ndarray, List[str]]:
        """Subsample names, seqs, sidxs to remove samples w/ no data after trimming."""

        # create mask with True is sample has data, False otherwise
        keepmask = np.zeros(len(self.names), dtype=bool)
        for sidx in range(len(self.names)):
            if np.all((self.seqs[sidx] == 78) | (self.seqs[sidx] == 45)):
                keepmask[sidx] = False
            else:
                keepmask[sidx] = True

        # apply mask to each
        self.seqs = self.seqs[keepmask, :]
        self.names = [i for i, j in zip(self.names, keepmask) if j]
        self.nidxs = [i for i, j in zip(self.nidxs, keepmask) if j]

        # save the sum of trimmed lens from each side.
        self.sumlens[0] = self.trim1[0] + self.trim2[0]
        self.sumlens[1] = self.trim1[1] + self.trim2[1]

        # return True if the locus should be filtered.
        return len(self.names) < self.data.params.min_samples_locus

@njit
def locus_left_trim(site_sample_covs: np.ndarray, min_cov_left: int) -> int:
    """Return the leftmost position where sample coverage >= the min
    allowed sample coverage provided in the hackers dict."""
    leftmost = np.where(site_sample_covs >= min_cov_left)[0]
    if leftmost.size:
        return leftmost.min()
    return 0

@njit
def locus_right_trim(site_sample_covs: np.ndarray, min_cov_right: int) -> int:
    """Return the rightmost position where sample coverage >= the min
    allowed sample coverage provided in the hackers dict."""
    rightmost = np.where(site_sample_covs >= min_cov_right)[0]
    if rightmost.size:
        return rightmost.max() + 1
    return 0


if __name__ == "__main__":

    import ipyrad as ip

    # --
    DATA = ip.Assembly(name="test")
    DATA.params.assembly_method = "reference"
    DATA.params.datatype = "pair3rad"
    DATA.params.min_samples_locus = 4
    DATA.params.filter_min_trim_len = 10
    DATA.params.restriction_overhang = ('ATCGG', 'TAGCTT')
    # DATA.params.trim_loci = (2, 2, 2, 2)
    DATA.hackers.trim_loci_min_sites = 5
    DATA.hackers.exclude_reference = False

    # SEQS array loaded from chunked .loci already passed a few filters
    SEQS = np.array([
        list(b"AAAATCGGACCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"NNNATCGGACCTTNNNNTTCCAANNNNNNNNNTACAAGCTA"),
        list(b"NNNATCGGTCCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"NNNATCGGNCCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"AAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTT"),
    ]).astype(np.uint8)

    NAMES = [f"sample-{i}" for i in range(len(SEQS))]
    NIDXS = ["..." for i in SEQS]

    # load Edges class
    test = EdgeTrimmer(DATA, NAMES, SEQS, NIDXS, debug=True)
    names, seqs, nidxs, over = test.run()

    if over:
        print('filtered')
    else:
        for ridx in range(seqs.shape[0]):
            print(bytes(seqs[ridx]).decode())
    # print(test.minsites_left, test.minsites_right)
    # print(test.data.hackers.trim_loci_min_sites)
    # print(test.exclude_ref)

    # # 1: trim for site coverage
    # test.get_edges()
    # print(test.edges)
    # arr = test.seqs[:, test.edges[0]:test.edges[3]]
    # print(arr.view("S1"))
