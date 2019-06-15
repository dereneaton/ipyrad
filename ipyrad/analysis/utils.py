#!/usr/bin/env python

"utility functions for the analysis tools"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard lib
import datetime
import time
import sys
import os

# third party
import numpy as np
from numba import njit, prange



class ProgressBar(object):
    """
    Print pretty progress bar
    """       
    def __init__(self, njobs, start=None, message=""):
        self.njobs = njobs
        self.start = (start if start else time.time())
        self.message = message
        self.finished = 0
        
    @property
    def progress(self):
        return 100 * (self.finished / float(self.njobs))

    @property
    def elapsed(self):
        return datetime.timedelta(seconds=int(time.time() - self.start))
 
    def update(self):
        # build the bar
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))

        # print to stderr
        print("\r[{}] {:>3}% {} | {:<12} ".format(*[
            hashes + nohash,
            int(self.progress),
            self.elapsed,
            self.message,
        ]), end="")
        sys.stdout.flush()



@njit
def jsubsample_snps(snpsmap, seed):
    "Subsample snps, one per locus, using snpsmap"
    np.random.seed(seed)
    sidxs = np.unique(snpsmap[:, 0])
    subs = np.zeros(sidxs.size, dtype=np.int64)
    idx = 0
    for sidx in sidxs:
        sites = snpsmap[snpsmap[:, 0] == sidx, 1]
        site = np.random.choice(sites)
        subs[idx] = site
        idx += 1
    return subs



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



def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
        .format(hashes + nohash, int(progress), elapsed, message),
        end="")
    sys.stdout.flush()    


class Params(object):
    "A dict-like object for storing params values with a custom repr"
    def __init__(self):
        self._i = 0

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __iter__(self):
        return self
    
    def __next__(self):
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if self._i > len(keys) - 1:
            self._i = 0
            raise StopIteration
        else:
            self._i += 1
            return keys[self._i - 1]
        

    def update(self, dict):
        self.__dict__.update(dict)


    def __repr__(self):
        "return simple representation of dict with ~ shortened for paths"
        _repr = ""
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if keys:
            _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
            for key in keys:
                _val = str(self[key]).replace(os.path.expanduser("~"), "~")
                _repr += _printstr.format(key, _val)
        return _repr
