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
import numba
import numpy as np


@numba.jit(nopython=True)
def get_spans(maparr): #, spans):
    """ 
    Get span distance for each locus in original seqarray. This
    is used to create re-sampled arrays in each bootstrap to sample
    unlinked SNPs. 
    """

    # test this; slower with empty spans inside?
    spans = np.zeros((maparr[-1, 0], 2), np.uint64)

    # start at 0, finds change at 1-index of map file
    bidx = 1
    spans = np.zeros((maparr[-1, 0], 2), np.uint64)

    # read through marr and record when locus id changes
    for idx in range(1, maparr.shape[0]):
        cur = maparr[idx, 0]
        if cur != bidx:
            # idy = idx + 1
            spans[cur - 2, 1] = idx
            spans[cur - 1, 0] = idx
            bidx = cur
    spans[-1, 1] = maparr[-1, -1]
    return spans



def progressbar(njobs, finished, start, msg):
    "simple progress bar for ipyrad analysis tools"

    # measure progress
    if njobs:
        progress = 100 * (finished / float(njobs))
    else:
        progress = 100

    # build the bar
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))

    # timestamp
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} | {} |".format(*[
        hashes + nohash,
        int(progress),
        elapsed,
        msg[0],
        msg[1],
    ]), end="")
    sys.stdout.flush()


class Params(object):
    "A dict-like object for storing params values with a custom repr"
    
    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())      
        _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
        for key in keys:
            _val = str(self[key]).replace(os.path.expanduser("~"), "~")
            _repr += _printstr.format(key, _val)
        return _repr