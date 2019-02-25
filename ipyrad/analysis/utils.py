#!/usr/bin/env python

"utility functions for the analysis tools"

# py2/3 compat
from __future__ import print_function
#from builtins import range

# standard lib
import datetime
import time
import sys
import os

# third party
from numba import njit, prange


# parallel and prange give a xncores speedup?
njit(parallel=True)
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



def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
        .format(hashes + nohash, int(progress), elapsed, message),
        end="")
    sys.stdout.flush()    


# New params class is iterable returning keys
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


# def progressbar(njobs, finished, start, msg):
#     "simple progress bar for ipyrad analysis tools"

#     # measure progress
#     if njobs:
#         progress = 100 * (finished / float(njobs))
#     else:
#         progress = 100

#     # build the bar
#     hashes = '#' * int(progress / 5.)
#     nohash = ' ' * int(20 - len(hashes))

#     # timestamp
#     elapsed = datetime.timedelta(seconds=int(time.time() - start))
#     print("\r[{}] {:>3}% {} | {:<12} | {} |".format(*[
#         hashes + nohash,
#         int(progress),
#         elapsed,
#         msg[0],
#         msg[1],
#     ]), end="")
#     sys.stdout.flush()