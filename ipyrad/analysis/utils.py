#!/usr/bin/env python

"""utility functions for the analysis tools
"""

from typing import Dict, List
import os
import numpy as np
from numba import njit, prange


def popfile_to_imap(path: str) -> Dict[str, List[str]]:
    """Parse popfile to an imap dictionary.

    The popfile should be formatted with whitespace separated
    samplename, popname lines.

    Parameters
    ----------
    path: str
        The path to a popfile.

    Example
    -------
    >>> imap = ipa.popfile_to_imap('popfile.txt')

    popfile format example
    ----------------------
    sample_A1   pop_A
    sample_A2   pop_A
    sample_B1   pop_B
    sample_B2   pop_B

    imap format example
    -------------------
    >>> imap = {
    >>>     'pop_A': ['sample_A1', 'sample_A2'],
    >>>     'pop_B': ['sample_B1', 'sample_B2'],
    >>> }
    """
    # TODO: support loading from a URL also.
    popfile = os.path.realpath(os.path.expanduser(path))
    imap = {}
    with open(popfile, 'r') as indata:
        data = [i.strip().split() for i in indata.readlines()]
        for i in data:
            if i[0] not in imap:
                imap[i[0]] = [i[1]]
            else:
                imap[i[0]].append(i[1])
    return imap


@njit
def jsubsample_snps(snpsmap: np.ndarray, seed: int) -> np.ndarray:
    """Subsample snps, one per locus, using snpsmap."""
    np.random.seed(seed)
    lidxs = np.unique(snpsmap[:, 0])
    keep = np.zeros(lidxs.size, dtype=np.int64)
    for sidx, lidx in enumerate(lidxs):
        sites = snpsmap[snpsmap[:, 0] == lidx, 1]
        site = np.random.choice(sites)
        keep[sidx] = site
    return keep


@njit
def jsubsample_loci(snpsmap, seed):
    """
    Return SNPs from re-sampled loci (shape = (nsample, ...can change)
    """
    np.random.seed(seed)

    # the number of unique loci with SNPs in this subset
    lidxs = np.unique(snpsmap[:, 0])

    # resample w/ replacement N loci
    lsample = np.random.choice(lidxs, len(lidxs))

    # the size of array to fill
    size = 0
    for lidx in lsample:
        size += snpsmap[snpsmap[:, 0] == lidx].shape[0]

    # fill with data
    subs = np.zeros(size, dtype=np.int64)
    idx = 0
    for lidx in lsample:
        block = snpsmap[snpsmap[:, 0] == lidx, 1]
        subs[idx: idx + block.size] = block
        idx += block.size
    return len(lidxs), subs


@njit(parallel=True)
def get_spans(maparr, spans):
    """
    Get span distance for each locus in original seqarray. This
    is used to create re-sampled arrays in each bootstrap to sample
    unlinked SNPs from. Used on snpsphy or str or ...
    """
    start = 0
    end = 0
    for idx in prange(1, spans.shape[0] + 1):
        lines = maparr[maparr[:, 0] == idx]
        if lines.size:
            end = lines[:, 3].max()
            spans[idx - 1] = [start, end]
        else:
            spans[idx - 1] = [end, end]
        start = spans[idx - 1, 1]

    # drop rows with no span (invariant loci)
    spans = spans[spans[:, 0] != spans[:, 1]]
    return spans

@njit
def count_snps(seqarr):
    """
    Count the number of SNPs in a np.uint8 seq array. This is used
    in window_extracter.
    """
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
