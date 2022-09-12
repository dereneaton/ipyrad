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

New idea for denovo catgs
--------------------------
Build a table for SNP position and offsets per sample. This info
needs to be recorded here before we can mask all of the indel (dash)
characters.

>>>     pre r1 aligns  |       internal aligns     |  post r2 aligns
>>> ---------NNNNNNNNNNTTTTTTTTTT...  --- ...  TTTTNNNNNNNNNN----
>>> ---------NNNNNNNNNNTTTTTTTTTT... ---- ...  TATTNNNNNNNNNN----
>>> ---------NNNNNNNNNNTTTTTTTTTT...  --- ...  TATTNNNNNNNNNN----
>>> ---------NNNNNNNNNNTTTTTTTTTT...  --- ...  T-TTNNNNNNNNNN----
>>> -----NNNNNNNNNNTTTTTTTTTATTTT...    - ...  TTTTNNNNNNNNNN----
>>> TTTTTTTTTTTTTTTTTTTTTTTTATTTT...    - ...  TTTTTTTTNNNNNNNNNN
>>>                         *                   *

To align the catgs with the reads requires keeping track of the start
position (offset) based on the number of bases trimmed from start.
Then, we must also keep all dash inserts in the locus until after catgs
are enumerated to advance for each one. For RAD and ddRAD this is easy,
but for GBS and other dtypes where reads can reverse comp match this
adds further complexity. The orientation (+/-) tells us whether to
advance from the left versus right to get the positions.

>>> TTTTTTTTTT... [trimmed left: 20] [offset: 0]
>>> TTTTTTTTTT... [trimmed left: 20] [offset: 0]
>>> TTTTTTTTTT... [trimmed left: 20] [offset: 0]
>>> TTTTTTTTTT... [trimmed left: 20] [offset: 0]
>>> TTTTTATTTT... [trimmed left: 16] [offset: 4]
>>> TTTTTATTTT... [trimmed left: 0]  [offset: 20]

But we only need to store for the SNP. The offset from start adds
to the position, whereas any internal indels subtract from it. An
indel at the SNP position equates to NaN for its catg position.

>>> Locus      Pos       S0   S1  S2   S3   S4   S5   S6   ...
>>>   0          5        5    5   5    5    9   25   NaN  ...
>>>   0        200      197  196 197  NaN  203  219   NaN  ...
"""

from typing import TypeVar, Tuple, List
from dataclasses import dataclass, field
import numpy as np
from numba import njit
from loguru import logger
from ipyrad.assemble.utils import comp, IPyradError

logger = logger.bind(name="ipyrad")
Assembly = TypeVar("Assembly")


@dataclass
class Locus:
    """Locus raw and filtered/trimmed data.

    """
    #########################################
    ## INIT Attributes
    #########################################
    lidx: int
    """: Index of this locus in the raw unfiltered loci file."""
    data: Assembly
    """: Assembly object with params settings."""
    names: List[str]
    """: List of names in alphanumeric order. Mutable. Samples can be dropped."""
    nidxs: List[str]
    """: string with information about the chrom/locus positions."""
    seqs: np.ndarray
    """: 2-D array of the samples x sites as an np.uint8 array."""

    #########################################
    ## Trimming attributes
    #########################################
    tseqs: np.ndarray = None
    """: 2-D array of the samples x sites as an np.uint8 array after trimming."""
    trimmed: Tuple[int,int] = (0, 0)
    """: Number of sites in trimmed from each end of .seqs. Index as X[t[0],-t[1]]"""
    pe_insert: Tuple[int,int] = (0, 0)
    """: Number of sites on either side of PE insert in .tseqs. Index as Y[t[0],-t[1]]"""
    site_sample_covs: np.ndarray = None
    """: Counts of sample coverage (non N-n) at each site (pre-trimming)"""
    _trim1: List[int] = field(default_factory=list)
    """: 0-indexed pos where (start, end) should be trimmed of lowcov edges, e.g., (5, -5)."""
    _trim2: List[int] = field(default_factory=list)
    """: 0-indexed pos where (start, end) should be trimmed of RE sites, """
    _min_sample_coverage_by_site: int = 0
    """: Minimum sample coverage required at locus edge (used for site cov
    trimming, not just sample cov trimming.)"""

    #########################################
    ## Filtering attributes
    #########################################
    is_dup: bool = False
    """: filter this locus b/c a sample name is present >1 time."""
    over_trimmed: bool = False
    """: filter this locus if trimming made it below minsample cov."""

    def __post_init__(self):
        # e.g., trim_loci_min_sites=10, but nsamples=4, then this is
        # set to 4. BUT, if the reference is in seqs but should not be
        # included then we further subtract one.
        self.min_sample_coverage_by_site = min(
            self.data.hackers.trim_loci_min_sites,
            self.seqs.shape[0] - int(self.data.drop_ref)
        )

    def get_sample_covs(self):
        """Fills site_sample_covs with data."""
        self.site_sample_covs = np.sum((self.seqs != 78) & (self.seqs != 45), axis=0)

    def run(self) -> Tuple[List, List, List, bool]:
        """Runs the steps involved in trimming and masking seqs."""
        if self.site_sample_covs is None:
            raise IPyradError("first run Locus.get_sample_covs()")

        # sets .edges array based on sample coverage, and returns bool filter
        bad_trim = self._set_edges_from_coverage()
        if bad_trim:
            print(f"edge trim filtered by coverage: {self._trim1}")
            self.over_trimmed = True

        # sets .trimseq by slicing .edges from .seqs
        self._trim_seqs_edges_1()

        # sets .trim from RE overhangs
        self._set_trim_from_re_overhangs()

        # updates .trim from user-entered .trim_loci param.
        self._set_trim_from_params()

        # update .edges for .trim lengths.
        self._trim_seqs_edges_2()

        # check that all data was not filtered.
        bad_trim = self._check_edges()
        if bad_trim:
            print(f"edge trim filtered by overhang: {self._trim2}")
            self.over_trimmed = True

        # get subsampled names, seqs, nidxs by removing any samples
        # that no longer have data after trimming, which can arise
        # when data are highly tiled.
        bad_trim = self._remove_empty_samples()
        if bad_trim:
            print(f"edge trim filtered by subsampling. Nsamples={len(self.names)}")
            self.over_trimmed = True

    def _set_edges_from_coverage(self) -> bool:
        """Fills self.trim1 with trim lens from each edge based on coverage.

        Finds the farthest left and right position where data is
        present (not N or -) and sample coverage >= min_sample_coverage_by_site
        WHICH USES self.data.hackers.trim_loci_min_sites.
        i.e., it calculates per-site sample coverage taking into account
        missing base calls.
        """
        self._trim1 = [0, 0]

        # get left and right most positions from jit'd functions.
        left = locus_left_trim(self.site_sample_covs, self.min_sample_coverage_by_site)
        right = locus_right_trim(self.site_sample_covs, self.min_sample_coverage_by_site)

        # return True if no positions had sufficient coverage and
        # this locus should be filtered, else return False
        if right <= left:
            return True

        # else, convert edges 1 from the pos to a -len from the end.
        # if 0 we replace with None for slicing.
        self._trim1 = (left, self.seqs.shape[1] - right)
        return False

    def _trim_seqs_edges_1(self):
        """Updates .seqs by trimming slice from .trim1.

        end is a length, so we make negative, or None if it is zero.
        """
        start = self._trim1[0]
        end = -self._trim1[1] if self._trim1[1] else None
        self.tseqs = self.seqs[:, start:end]

    def _set_trim_from_re_overhangs(self):
        """Fills self.trim with extra edge sizes that should be trimmed.

        Fuzzy search for restriction overhang sequence to be trimmed.
        This applies to the .trimseq array that has already been
        trimmed based on sample coverage.

        TODO: this func could be jit'd for speed improvements.
        """
        self._trim2 = [0, 0]
        for cutter in self.data.params.restriction_overhang:
            if not cutter:
                continue

            # convert to uint8 array
            forward = np.array(list(cutter.encode()), dtype=np.uint8)
            revcomp = np.array(list(comp(cutter).encode()), dtype=np.uint8)[::-1]

            # compare match over cut size skipping Ns and allow 0.25 diffs
            slx = slice(0, forward.shape[0])
            matching = self.tseqs[:, slx] == forward
            mask = np.where(
                (self.tseqs[:, slx] != 78) & (self.tseqs[:, slx] != 45))
            matchset = matching[mask] # non-masked sites that match the RE
            prop = float(matchset.sum()) / matchset.size
            if prop >= 0.75:
                self._trim2[0] = len(forward)
                # print(f"forward={forward}\nmatching={matching}\nmask={mask}\nmatchset={matchset}\nprop={prop:.2f}")
                continue

            # compare match over cut size skipping Ns and allow 0.25 diffs
            slx = slice(self.tseqs.shape[1] - revcomp.shape[0], self.tseqs.shape[1])
            matching = self.tseqs[:, slx] == revcomp
            mask = np.where(
                (self.tseqs[:, slx] != 78) & (self.tseqs[:, slx] != 45))
            matchset = matching[mask]
            prop = float(matchset.sum()) / matchset.size
            if prop >= 0.75:
                self._trim2[1] = len(revcomp)
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
        self._trim2[0] = max(self._trim2[0], self.data.params.trim_loci[0])
        self._trim2[1] = max(self._trim2[1], abs(self.data.params.trim_loci[3]))

    def _trim_seqs_edges_2(self):
        """Updates .seqs by trimming slice from .trim2.

        This is the final seqs that will be returned.
        """
        start = self._trim2[0]
        end = -self._trim2[1] if self._trim2[1] else None
        self.tseqs = self.tseqs[:, start:end]
        # self.edges[1] -= self.trims[1]
        # self.edges[2] += self.trims[2]

    def _check_edges(self) -> bool:
        """Return True if no locus length is left untrimmed"""
        if self.tseqs.shape[1] < self.data.params.filter_min_trim_len:
            return True
        return False

    def _remove_empty_samples(self) -> Tuple[List[str], np.ndarray, List[str]]:
        """Subsample names, seqs, sidxs to remove samples w/ no data after trimming."""

        # create mask with True if sample has data, False otherwise
        keepmask = np.zeros(len(self.names), dtype=bool)
        for sidx in range(len(self.names)):
            if np.all((self.tseqs[sidx] == 78) | (self.tseqs[sidx] == 45)):
                keepmask[sidx] = False
            else:
                keepmask[sidx] = True

        # apply mask to each
        self.tseqs = self.tseqs[keepmask, :]
        self.seqs = self.seqs[keepmask, :]
        self.names = [i for i, j in zip(self.names, keepmask) if j]
        self.nidxs = [i for i, j in zip(self.nidxs, keepmask) if j]

        # save the sum of trimmed lens from each side. For use w/ catgs later.
        trim0 = self._trim1[0] + self._trim2[0]
        trim1 = self._trim1[1] + self._trim2[1]
        self.trimmed = (trim0, trim1)
        # self.sumlens[0] = self.trim1[0] + self.trim2[0]
        # self.sumlens[1] = self.trim1[1] + self.trim2[1]

        # return True if the locus should be filtered.
        return len(self.names) < self.data.params.min_samples_locus

    def _set_pe_inserts(self) -> None:
        """Fills .pe_insert Tuple with .tseqs index to mask."""
        if not 110 in self.tseqs[0]:
            return
        # find position of the 110 spacer
        pos = np.where(self.tseqs[0] == 110)[0][0]

        # get rightmost position in tseqs with site cov left of pos
        left_ = locus_right_trim(
            self.site_sample_covs[self.trimmed[0]:pos + 4],
            self.min_sample_coverage_by_site)

        # get leftmost pos with site cov right of pos, and add pos len to it
        # right_ = locus_left_trim(
            # self.site_sample_covs[(pos + 4):-self.trimmed[1]],
            # self.min_sample_coverage_by_site)
        # right_ += (pos + 4) - 1
        right_ = locus_left_trim(
            self.site_sample_covs[self.trimmed[0] + (pos+4):],
            self.min_sample_coverage_by_site)

        # index positions on .tseqs where pe_insert is located
        print(right_)
        self.pe_insert = (left_, pos + 4 + right_)

    def mask_inserts(self, mask_int=78) -> None:
        """Mask paired-end insert.

        If the 'nnnn' separator is present then there are likely messy
        aligned edges near to it that should be cleaned. This function
        will identify the positions to the left and right where it is
        min_sites depth. Between these regions is masked. This applies
        after edge trimming, and so there SHOULD ALWAYS be sites to the
        left and right that meet the min sites criterion. Indels outside
        of the pair insert region are not affected by this masking.

        Note
        ----
        The index of the masked region is stored to .pe_insert, which
        indexes relative to .tseqs (post edge-trimmed array).

        Example
        --------
        >>> # BEFORE   8866620000000000000000000000001011111127888
        >>>    get region   **********************************
        >>> GGTGGAAGTCATCAGTNNNN-NNNNNNnnnn--------NNN-NNNNNNNNCGATGTGAG
        >>> GGTGGAAGTCATCAGTTNNNNNNNNNNnnnnNNNNNNNNNNC-CCTAAACCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNN-NNNNNNnnnn--------NNN-NNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNN----NNNNNnnnn--------NNN-NNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNN----NNNNNnnnn--------NNN-NNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNN----NNNNNnnnn--------NNNNNNNNNNCCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNN-NNNNNNnnnn--------NNN--NNNNNNNCGATGTGAG
        >>> GGTGGAAGTCATCAGTTNNNNNNNNNNnnnn--------NNN-NNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNN-NNNNNNnnnn--------NNN-NNNNNNNCCGATGTGAG
        >>>
        >>> # AFTER
        >>>    mask it      **********************************
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        >>> GGTGGAAGTCATCAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGATGTGAG
        """
        self._set_pe_inserts()
        # before = []
        # for x in range(self.seqs.shape[0]):
        #     before.append(bytes(self.seqs[x, (left_-5):(right_ + 5)]).decode())
        # logger.warning(f"{left_},{right_}\n" + "\n".join(before))

        # apply mask to set insert to Ns
        self.tseqs[:, self.pe_insert[0]:self.pe_insert[1]] = mask_int

        # after = []
        # for x in range(self.seqs.shape[0]):
        #     after.append(bytes(self.seqs[x, (left_-5):(right_ + 5)]).decode())
        # logger.warning("\n" + "\n".join(after) + "\n\n")


# These functions are used in different funcs on either .seqs or .tseqs
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
    ip.set_log_level("DEBUG")

    # --
    DATA = ip.Assembly(name="test")
    DATA.params.assembly_method = "reference"
    DATA.params.datatype = "pair3rad"
    DATA.params.min_samples_locus = 4
    DATA.params.filter_min_trim_len = 10
    DATA.params.restriction_overhang = ('','')  #('ATCGG', 'TAGCTT')
    # DATA.params.trim_loci = (2, 2, 2, 2)
    #DATA.hackers.trim_loci_min_sites = 4
    DATA.hackers.exclude_reference = False

    # SEQS array loaded from chunked .loci already passed a few filters
    SEQS1 = np.array([
        list(b'NNNNNNNNNNATCGGTGACGTTTTGTCCTTATTYTCGTAATTKATAATWATATGCCCTCTTTCCACTTCACGACGGTTCCAAACCATCGTGATCTTAGGYGTCACCTGGTRTTATTAACAAAAWATATATTRAATTATATATTATAGTGTAAATTCAAAAAAATAAATTTATGTGTCACATAAAAAGTGTCCAGARCCATTAAAAATATTTGCTTTTTTTAGTAACTGTCAATAATCACAAACATCTTATTGAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTGACGTTTTGTCCTTATTTTCGTARTYGATAATWATATGCCCTCTTTCCACTTCACGACGGTTCCAAACCATCGTGATCTTAGGCGTCACCTGGTRTTWTTAACAAAAAATATATTGAATTAWATATTATAGTGTAAATTCAAAAAAATAAATTTATGTGTCACATAAAAAGTGTCCAGAGCCATTAAAAATATTTGCTTTTTTTAGTAACTGTCAATAATCACAAACATCTTATTGAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTGACGTTTTGTCCTTATTYTCGTAATTGATAATWATATGCCCTCTTTCCACTTCACGACGGTTCCAAACCATCGTGATCTTAGGCGTCACCTGGTRTTATTAACAAAAAATATATTGAATTATAT---ATAGTGTAAATTCAAAAAAATAAATTTATGTGTCACATAAAAAGTGTCCAGAGCCATTAAAAATATTTGCTTTTTTTAGTAACTGTCAATAATCACAAACATYTTATTGAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTGACGTTTTGTCCTTATTCTCGTAATTGATAATAATATGCCCTCTTTCCACTTCACGACGGTTCCAAACCATCGTGATCTTAGGCGTCACCTGGTATTATTAACAAAAAATATATTGAATTAAATATTATAGTGTAAATTCAAAAAAATAAATTTATGTGTCACATAAAAAGTGTCCAGAGCCATTAAAAATATTTGCTTTTTTTAGTAACTGTCAATAATCACAAACATCTTATTGAAGCTANNNNNNNNNN'),
    ]).astype(np.uint8)

    SEQS2 = np.array([
        list(b'NNNNNNNNNNATCGGTAGGCTTCCTGAGGGCCTCCTACCTGAAAGGAATTTCTCAACTAGTCAACATAATGACAAAAGCACCTATCAGCCACGAGCCCAAAGGAAGAAGCACCATAACAAGCCAGAGATAAGAGATCAGATATTGACCACGCNNNNNNNNNNNnnnnNNNNNNNNNNCAGACCAAATTAAATATGTCTCACAAATCCATTATGGTATTTAACCTATTGCTGACAGTACCAAAATAAGTTGCAATCAAACAGAAAAATAATTCAAATATAGAAGGATCAGTGAAATCCTCTGGCATTTGAGATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTAGGCTTCCTGAGGGCCTCCTACCTGAAAGGAATTTCTSAACTAGTCAACATAATGACAAAAGCACCTATCAGCCACGAGCCCAAAGGAAGAAGCACCATAACAAGTCAGAGATAAGAGATCAGATATTGACCACGCNNNN-NNNNNNnnnnNNNNNNNNNNCAGACCAAATTAAATATGTCTCACAAATCCATTATGGTATTTAACCTATTKCTGACAGTACCAAAATAAGTTGCAATCAAACAGAAAAATAATTCAAATATAGAAGGATCAGTGAAATCCTCTGGCATTTGAGATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTAGGCTTCCTGAGGGCCTCCTACCTGAAAGGAATTTCTSAACTAGTCAACATAATGACAAAAGCACCTATCAGCCACGAGCCCAAAGGAAGAAGCACCATAACAAGTCAGAGATAAGAGATCAGATATTGACCACGNNNN--NNNNNNnnnnNNNNNNNNNNCAGACCAAATTAAATATGTCTCACAAATCCATTATGGTATTTAACCTATTGCTGACAGTACCAAAATAAGTTGCAATCAAACAGAAAAATAATTCAAATATAGAAGGATCAGTGAAATCCTCTGGCATTTGAGATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTAGGCTTCCTGAGGGCCTCCTACCTGAAAGGAATTTCTCAAYTAGTCAACWTAATGACAAAAGCACCTATCAGCCACGAGCCCAAAGGAAGAAGCACCATAACAAGYCAGAGATAAGAGATCAGATATTGACCACGCNNNN-NNNNNNnnnnNNNNNNNNNNCAGACCAAATTAAATATGTCTCACAAATCCATTATGGTATTTAACCTATTGCTGACAGTACCAAAATAAGTTGCAATCAAACAGAAAAATAATTCAAATATAGAAGGATCAGTGAAATCCTCTGGCATTTGAGATAAAGCTANNNNNNNNNN'),
    ]).astype(np.uint8)

    SEQS3 = np.array([
        list(b'NNNNNNNNNNATCGGTCCACCTCAGTATCGTAGTAGTACAATTTTTGAAGATGCAACACCCGAGATGGTTAGAGACTTCTTTTGGGATGATAAATTTCGACCAACGTTTGACCCCATGCTCATAAATTCTGAAACACTTGAAGAGTCTCGTACNNNNNNNNNNNnnnnNNNNNNNNNNN--NNAGTGGCATGTAACTGACTTTTGCTGGTTGTCTATGCAGTTTCCGTTCTTCTGTAGCGACCGAGAGTACATAATCGGCCGTAGGATATGGGAGTGTGAAAGAACATTTTACTGCGTGACAAAGGTATAGACTATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTCCACCTCAGTATCGTAGTAGTACAATTTTTGAAGATGCAACACCCGAGATGGTTAGAGACTTCTTTTGGGATGATAAATTTCGACCAACGTTTGACCCCATGCTCATAAATTCTGAAACACTTGAAGAGTCTCGTNNNN---NNNNNNnnnnNNN---NNNNN--NNAGTGGCATGTAACTGACTTTTGCTGGTTGTCTATGCAGTTTCCGTTCTTCTGTAGCGACCGAGAGTACATAATCGGCCGTAGGATATGGGAGTGTGAAAGAACATTTTACTGCGTGACAAAGGTATAGACTATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTCCACCTCAGTATCGTAGTAGTACAATTTTTGAAGATGCAACACCCGAGATGGTTAGAGACTTCTTTTGGGATGATAAATTTCGACCAACGTTTGACCCCATGCTCATAAATTCTGAAACACTTGAAGAGTCTCGTNNNN---NNNNNNnnnnNNN---NNNNN--NNAGTGGCATGTAACTGACTTTTGCTGGTTGTCTATGCAGTTTCCGTTCTTCTGTAGCGACCGAGAGTACATAATCGGCCGTAGGATATGGGAGTGTGAAAGAACATTTTACTGCGTGACAAAGGTATAGACTATAAAGCTANNNNNNNNNN'),
        list(b'NNNNNNNNNNATCGGTCCACCTCAGTATCGTAGTAGTACAATTTTTGAAGATGCAACACCCGAGATGGTWAGAGACTTCTTYTGGGATGATAAATTTCGACCAACGTTTGACCCYATGCTCATAAATTCTGAAACRCTTGAAGAGTCTCGTACTNNNNNNNNNNnnnnNNN---NNNNNNNRGTAKRGCATGTAACTGACTTTTGYTGGYTGTYTATGCAGTTTCCSTTCTTCTGYAGYGACCGAGAGTACATAATCGGCCGTAGGATATGGGAGTGTGAAAGAACATTTTACTGCGTGACAAAGGTATAGACTATAAAGCTANNNNNNNNNN'),
    ]).astype(np.uint8)

    SEQS = np.array([
        list(b"AAAATCGGACCTTNNNNTTCCAACCGNnnnnNTACAAGCTA"),
        list(b"NNNATCGGACCTTNNNNTTCCAANNNNnnnnNTACAAGCTA"),
        list(b"NNNATCGGTCCTTNNNNTTCCAACCGNnnnnNTACAAGCTA"),
        list(b"NNNATCGGNCCTTNNNNTTCCAACCGNnnnnNTACAAGCTA"),
        list(b"AAAANNNNNNNNNNNNNNNNNNNNNNNnnnnNNNNNNNTTT"),
    ]).astype(np.uint8)

    SEQS = SEQS3

    NAMES = [f"sample-{i}" for i in range(len(SEQS))]
    NIDXS = ["..." for i in SEQS]

    # load Edges class
    DATA.drop_ref = False
    test = Locus(0, DATA, NAMES, NIDXS, SEQS) #debug=True)
    test.get_sample_covs()
    test.run()
    test.mask_inserts()
    print(test.trimmed)
    print(test.pe_insert)

    pos = test.pe_insert[0]
    end = test.pe_insert[1]
    ttt = test.trimmed[0]
    eee = test.trimmed[1]
    for nidx, name in enumerate(NAMES):
        seq = bytes(test.seqs[nidx]).decode()
        if pos:
            print(f"{seq[:ttt+5]}...{seq[pos-10+ttt:end+10+ttt]}...{seq[-(eee + 5):]}")
        else:
            print(f"{seq[:ttt+5]}...{seq[-(eee + 5):]}")
    # print("----")
    # for nidx, name in enumerate(NAMES):
    #     seq = bytes(test.seqs[nidx]).decode()
    #     if pos:
    #         print(f"{seq[:ttt+5]}...{seq[pos+ttt:end+ttt+4+1]}...{seq[-(eee + 5):]}")
    #     else:
    #         print(f"{seq[:ttt+5]}...{seq[-(eee + 5):]}")
    print("----")
    for nidx, name in enumerate(NAMES):
        seq = bytes(test.tseqs[nidx]).decode()
        if pos:
            print(f"{seq[:5]}...{seq[(pos-10):(end+10)]}...{seq[-5:]}")
        else:
            print(f"{seq[:5]}...{seq[-5:]}")
        # print(bytes(test.tseqs[nidx][pos-10:end+10]).decode())

    print(test.pe_insert)
    print(test.site_sample_covs[(ttt+pos-1):(ttt+end+6)])

    # if over:
    #     print('filtered')
    # else:
    #     for ridx in range(seqs.shape[0]):
    #         print(bytes(seqs[ridx]).decode())
    # print(test.minsites_left, test.minsites_right)
    # print(test.data.hackers.trim_loci_min_sites)
    # print(test.exclude_ref)

    # # 1: trim for site coverage
    # test.get_edges()
    # print(test.edges)
    # arr = test.seqs[:, test.edges[0]:test.edges[3]]
    # print(arr.view("S1"))
