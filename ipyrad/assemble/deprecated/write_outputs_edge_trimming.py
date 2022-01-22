#!/usr/bin/env python

"""
Locus edge trimming based on presence of cutsites, 
locus trimming params, and site trimming hacker param.
"""

import numpy as np
from numba import njit
from ipyrad.assemble.utils import comp


class Edges:
    """
    Trims edges of overhanging sequences, cutsites, and pair inserts.
    """
    def __init__(self, data, seqs):
        self.data = data
        self.seqs = seqs

        # params
        self.bad = False
        self.exclude_ref = self.data.hackers.exclude_reference
        self.edges = np.array([0, 0, 0, self.seqs.shape[1]])
        self.trims = np.array([0, 0, 0, 0])  # self.seqs.shape[1]])
        self.minlen = self.data.params.filter_min_trim_len

        # -1 to site coverage if ref is excluded from the count
        self.minsites_left = self.data.hackers.trim_loci_min_sites
        self.minsites_right = self.data.hackers.trim_loci_min_sites
        if self.data.is_ref and self.exclude_ref:
            self.minsites_left -= 1
            self.minsites_right -= 1
        self.minsites_left = min(self.minsites_left, self.seqs.shape[0])
        self.minsites_right = min(self.minsites_right, self.seqs.shape[0])        

        # to be filled
        self.trimseq = None


    def get_edges(self):
        """
        Find the proper edges of the locus based on sample coverage,
        and trimming parameters.        
        """
        # get .edges of good locus or .bad
        self.trim_for_coverage()

        # fill trimseq with the trimmed sequence array
        self.trimseq = self.seqs[:, self.edges[0]:self.edges[3]]

        # apply edge filtering to locus
        if not self.bad:
            try:
                self.trim_overhangs()
                self.trim_param_trim_loci()
            except Exception as inst:
                self.bad = True
                # print(f"bad trim edges:\n{self.seqs[:, :15]}")

        # check that locus has anything left
        self.trim_check()


    def trim_for_coverage(self):
        """
        trim edges to where data is not N or - and site coverage 
        meets the argument minimums.
        """
        # how much cov is there at each site?
        mincovs = np.sum((self.seqs != 78) & (self.seqs != 45), axis=0)

        # locus left trim
        self.edges[0] = locus_left_trim(self.minsites_left, mincovs)
        self.edges[3] = locus_right_trim(self.minsites_right, mincovs)
        if self.edges[3] <= self.edges[0]:
            self.bad = True

        # find insert region for paired data to mask it...
        self.edges[1] = 0
        self.edges[2] = 0


    def trim_overhangs(self):
        """
        fuzzy match to trim the restriction_overhangs from r1 and r2
        """
        # trim left side for overhang
        for rcut in self.data.params.restriction_overhang:

            # skip if None
            if not rcut:
                continue

            # convert to integers
            cutter = np.array(list(rcut.encode()))
            rcutter = np.array(list(comp(rcut).encode()))[::-1]

            # compare match over cut size skipping Ns and allow .25 diffs
            slx = slice(0, cutter.shape[0])
            matching = self.trimseq[:, slx] == cutter
            mask = np.where(
                (self.trimseq[:, slx] != 78) & (self.trimseq[:, slx] != 45))
            matchset = matching[mask]
            if float(matchset.sum()) / matchset.size >= 0.75:
                self.trims[0] = len(cutter)

            # trim right side for overhang
            slx = slice(
                self.trimseq.shape[1] - rcutter.shape[0], self.trimseq.shape[1])
            matching = self.trimseq[:, slx] == rcutter
            mask = np.where(
                (self.trimseq[:, slx] != 78) & (self.trimseq[:, slx] != 45))
            matchset = matching[mask]
            if float(matchset.sum()) / matchset.size >= 0.75:
                self.trims[3] = len(rcutter)


    def trim_param_trim_loci(self):
        """
        user entered hard trims
        """
        self.trims[0] = max(
            [self.trims[0], self.data.params.trim_loci[0]])
        self.trims[1] = (
            self.trims[1] - self.data.params.trim_loci[1]
            if self.trims[1] else 0)
        self.trims[2] = (
            self.trims[2] + self.data.params.trim_loci[2]
            if self.trims[2] else 0)
        self.trims[3] = max(
            [self.trims[3], self.data.params.trim_loci[3]])


    def trim_check(self):
        """
        Stores self.bad boolean based on trim limitations.
        """
        self.edges[0] += self.trims[0]
        self.edges[1] -= self.trims[1]
        self.edges[2] += self.trims[2]
        self.edges[3] -= self.trims[3]

        # checks
        if any(self.edges < 0):
            self.bad = True
        if self.edges[3] <= self.edges[0]:
            self.bad = True
        if self.edges[1] > self.edges[2]:
            self.bad = True
        # check total length including insert
        if (self.edges[3] - self.edges[0]) < self.minlen:
            self.bad = True


# -------------------------------------------------------------
# jitted Edges functions (njit = nopython mode)
# -------------------------------------------------------------

@njit
def locus_left_trim(minsamp, mincovs):
    leftmost = np.where(mincovs >= minsamp)[0]
    if leftmost.size:
        return leftmost.min()
    return 0

@njit
def locus_right_trim(minsamp, mincovs):
    rightmost = np.where(mincovs >= minsamp)[0]
    if rightmost.size:
        return rightmost.max() + 1
    return 0



if __name__ == "__main__":

    import ipyrad as ip

    # --
    DATA = ip.Assembly(name="test")
    DATA.params.assembly_method = "reference"
    DATA.params.datatype = "pair3rad"
    DATA.params.min_samples_locus = 6
    DATA.params.restriction_overhang = ('ATCGG', 'TAGCTT')
    # DATA.params.trim_loci = (2, 2, 2, 2)
    DATA.hackers.trim_loci_min_sites = 4
    DATA.hackers.exclude_reference = False

    SEQS = np.array([
        list(b"AAAATCGGACCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"AAAATCGGACCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"AAAATCGGACCTTNNNNTTCCAANNNNNNNNNTACAAGCTA"),
        list(b"NNNATCGGACCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
        list(b"NNNATCGGNCCTTNNNNTTCCAACCGNNNNNNTACAAGCTA"),
    ]).astype(np.uint8)
    test = Edges(DATA, SEQS)
    print(test.minsites_left, test.minsites_right)
    print(test.data.hackers.trim_loci_min_sites)
    print(test.exclude_ref)

    # 1: trim for site coverage
    test.get_edges()
    print(test.edges)
    arr = test.seqs[:, test.edges[0]:test.edges[3]]
    print(arr.view("S1"))
