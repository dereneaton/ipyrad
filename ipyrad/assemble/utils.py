#!/usr/bin/env python

""" 
Various sequence manipulation util functions used by different
parts of the pipeline
"""

from __future__ import print_function
try:
    from itertools import izip, takewhile
except ImportError:
    from itertools import takewhile
    izip = zip

import os
import sys
import time
import socket
import string
import tempfile
import datetime

from loguru import logger
import pandas as pd
import numpy as np


BADCHARS = (
    string.punctuation
    .replace("_", "")
    .replace("-", "")
    .replace(".", "") + " "
)




class AssemblyProgressBar(object):
    """
    Print pretty progress bar with printings specific to Assembly object

    Parameters:
    -----------
    jobs: dict
    start: time.time()
    message: str
    data: ipyrad.Assembly
    """       
    def __init__(self, jobs, start, message, data):
        self.jobs = jobs
        self.start = (start if start else time.time())
        self.message = message
        self.cli = data._cli
        self.spacer = data._spacer
        self.finished = 0

        # filled by check_success
        self.results = {}


    @property
    def progress(self):
        "returns the percent progress as a float"
        if not self.jobs:
            return 0
        return 100 * (self.finished / float(len(self.jobs)))

    @property
    def elapsed(self):
        "returns the elapsed time in nice format"
        return datetime.timedelta(seconds=int(time.time() - self.start))


    def update(self):
        "flushes progress bar at current state to STDOUT"
        # build the bar
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))

        # print to stderr
        if self.cli:
            print("\r{}[{}] {:>3}% {} | {:<12} ".format(
                self.spacer,
                hashes + nohash,
                int(self.progress),
                self.elapsed,
                self.message[0],
            ), end="")
        else:
            print("\r{}[{}] {:>3}% {} | {:<12} | {} |".format(*[
                self.spacer,
                hashes + nohash,
                int(self.progress),
                self.elapsed,
                self.message[0],
                self.message[1],
            ]), end="")
        sys.stdout.flush()


    def block(self):
        """
        Tracks completion of asynchronous result objects in a while 
        loop until they are finished, checking every 0.5 seconds.
        Prints progress either continuously or occasionally, depending
        on TTY output type.
        """
        # get TTY output type of STDOUT
        isatty = os.isatty(1)

        # store the current value
        cur_finished = 0

        # loop until all jobs are finished
        while 1:
            self.finished = sum([i.ready() for i in self.jobs.values()])
            time.sleep(1)
            
            # flush progress and end
            if self.finished == len(self.jobs):
                self.update()
                print("")
                break

            # -- decide whether to flush based on tty ---
            # if value changed then print
            if cur_finished != self.finished:
                self.update()

            # else, only print every 10 minutes if not in tty
            elif not isatty:
                if not int(self.elapsed.seconds) % 60:
                    self.update()
            
            # normal tty print every second
            else:
                self.update()

            # update cur_finished to match finished
            cur_finished = self.finished



    def check(self):
        """
        Will log and raise an error with traceback.
        """
        # check for failures:
        for job in self.jobs:
            if not self.jobs[job].successful():
                # raise the exception from the job and catch it
                logger.debug(job, self.jobs)
                try:
                    self.results[job] = self.jobs[job].get()
                except Exception as inst:
                    # report it to the logger file.
                    # logger.exception("exception on remote engine:")
                    # raise
                    # report it to stdout
                    raise IPyradError("Exception on remote engine:") from inst



class VsearchProgressBar(AssemblyProgressBar):
    """
    Subset Progress bar for step 6 clustering progress.
    """
    @property
    def progress(self):
        "get progress from vsearch stdout"
        logger.warning(self.jobs)
        logger.warning(self.jobs[0])
        logger.warning(self.jobs[0].stdout)        
        if self.jobs[0].stdout:
            return int(self.jobs[0].stdout.split()[-1])
        return 0



LOGFORMAT = (
    "{time:hh:mm:ss} <level>{level: <7}</level> <white>|</white> "
    "<magenta>PID:{process.id}</magenta> <white>|</white> "
    "<cyan>{file}:{line}</cyan> <white>|</white> "
    "<level>{message}</level>"
)


def set_loglevel(loglevel="DEBUG", logfile=None):
    """
    Config and start the logger
    """
    if logfile is None:
        logfile = os.path.join(tempfile.gettempdir(), "ipyrad-log.txt")

    config = {
        "handlers": [
            dict(
                sink=logfile,
                rotation="50 MB",
                format=LOGFORMAT,
                level=loglevel,
                enqueue=True,
                colorize=False,
            ),
        ]
    }
    logger.configure(**config)
    logger.enable("ipyrad")



class IPyradError(Exception):
    def __init__(self, *args, **kwargs):
        # raise the exception with this string message and a traceback
        Exception.__init__(self, *args, **kwargs)


# class IPyradError(Exception):
#     """
#     Exception handler that does clean exit for CLI, but also prints
#     the traceback and cleaner message for API.
#     """
#     def __init__(self, *args, **kwargs):
#         # raise the exception with this string message and a traceback
#         Exception.__init__(self, *args, **kwargs)

#         # but suppress traceback and exit on CLI
#         if not ipyrad.__interactive__:
#             # clean exit for CLI that still exits as an Error (e.g. for HPC)
#             sys.tracebacklimit = 0
#             SystemExit(1)


## utility functions/classes
class Params(object):
    """ 
    A dict-like object for storing params values with a custom repr
    that shortens file paths, and which makes attributes easily viewable
    through tab completion in a notebook while hiding other funcs, attrs, that
    are in normal dicts. 
    """
    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())
        if keys:
            _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
            for key in keys:
                _val = str(self[key]).replace(os.path.expanduser("~"), "~")
                _repr += _printstr.format(key, _val)
        return _repr




class ObjDict(dict):
    """
    Object dictionary allows calling dictionaries in a more
    pretty and Python fashion for storing Assembly data
    """
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __repr__(self):
        result = ""
        if "outfiles" in self.keys():
            dirs_order = ["fastqs", "edits", "clusts", "consens", "outfiles"]
            for key in dirs_order:
                result += key + " : " + self[key] + "\n"
        else:
            for key in sorted(self):
                result += key + " : " + str(self[key]) + "\n"
        return result



CDICT = {i: j for i, j in zip("CATG", "0123")}


# ## used for geno output
# VIEW = {
#     "R": ("G", "A"),
#     "K": ("G", "T"),
#     "S": ("G", "C"),
#     "Y": ("T", "C"),
#     "W": ("T", "A"),
#     "M": ("C", "A"),
#     "A": ("X", "X"),
#     "T": ("X", "X"),
#     "G": ("X", "X"),
#     "C": ("X", "X"),
#     "N": ("X", "X"),
#     "-": ("X", "X"),
#     }

# ## used in hetero() func of consens_se.py
# TRANS = {
#     ('G', 'A'): "R",
#     ('G', 'T'): "K",
#     ('G', 'C'): "S",
#     ('T', 'C'): "Y",
#     ('T', 'A'): "W",
#     ('C', 'A'): "M",
#     }

# # used in write_outfiles.write_geno
TRANSFULL = {
    ('G', 'A'): "R",
    ('G', 'T'): "K",
    ('G', 'C'): "S",
    ('T', 'C'): "Y",
    ('T', 'A'): "W",
    ('C', 'A'): "M",
    ('A', 'C'): "M",
    ('A', 'T'): "W",
    ('C', 'T'): "Y",
    ('C', 'G'): "S",
    ('T', 'G'): "K",
    ('A', 'G'): "R",
    }

# TRANSINT = {
#     (71, 65): 82,
#     (71, 84): 75,
#     (71, 67): 83,
#     (84, 67): 89,
#     (84, 65): 87,
#     (67, 65): 77,
#     (65, 67): 77,
#     (65, 84): 87,
#     (67, 84): 89,
#     (67, 71): 83,
#     (84, 71): 75,
#     (65, 71): 82,
#     }


# GENOOPTS = {
#     b"C": (b"S", b"Y", b"M"),
#     b"A": (b"R", b"W", b"M"),
#     b"T": (b"K", b"Y", b"W"), 
#     b"G": (b"R", b"K", b"S"),
# }


## used for resolving ambiguities
AMBIGS = {
    "R": ("G", "A"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("T", "C"),
    "W": ("T", "A"),
    "M": ("C", "A"),
    }




def chroms2ints(data, intkeys):
    """
    Parse .fai to get a dict with {chroms/scaffolds: ints}, or reversed.
    """
    fai = pd.read_csv(
        data.params.reference_sequence + ".fai",
        names=['scaffold', 'length', 'start', 'a', 'b'],
        sep="\t",
    )
    # Allow CHROM to take integer values, here cast them to str
    fai["scaffold"] = fai["scaffold"].astype(str)

    faidict = {j: i for i, j in enumerate(fai.scaffold)}
    if intkeys:
        revdict = {j: i for i, j in faidict.items()}
        return revdict
    return faidict




def ambigcutters(seq):
    """
    Returns both resolutions of a cut site that has an ambiguous base in
    it, else the single cut site
    """
    resos = []
    if any([i in list("RKSYWM") for i in seq]):
        for base in list("RKSYWM"):
            if base in seq:
                resos.append(seq.replace(base, AMBIGS[base][0]))
                resos.append(seq.replace(base, AMBIGS[base][1]))
        return resos
    else:
        return [seq, ""]



def splitalleles(consensus):
    """ 
    Takes diploid consensus alleles with phase data stored as a mixture
    of upper and lower case characters and splits it into 2 alleles 
    """

    ## store two alleles, allele1 will start with bigbase
    allele1 = list(consensus)
    allele2 = list(consensus)

    hidx = [i for (i, j) in enumerate(consensus) if j in "RKSWYMrkswym"]

    ## do remaining h sites
    for idx in hidx:
        hsite = consensus[idx].encode()
        if hsite.isupper():
            allele1[idx] = PRIORITY[hsite].decode()
            allele2[idx] = MINOR[hsite].decode()
        else:
            allele1[idx] = MINOR[hsite.upper()].decode()
            allele2[idx] = PRIORITY[hsite.upper()].decode()

    ## convert back to strings
    allele1 = "".join(allele1)
    allele2 = "".join(allele2)

    return allele1, allele2


# used by clustmap
def comp(seq):
    """ returns a seq with complement. Preserves little n's for splitters."""
    ## makes base to its small complement then makes upper
    return seq.replace("A", 't')\
              .replace('T', 'a')\
              .replace('C', 'g')\
              .replace('G', 'c')\
              .replace('n', 'Z')\
              .upper()\
              .replace("Z", "n")

# used by clustmap
def bcomp(seq):
    """ returns a seq with complement. Preserves little n's for splitters."""
    ## makes base to its small complement then makes upper
    return seq.replace(b"A", b't')\
              .replace(b'T', b'a')\
              .replace(b'C', b'g')\
              .replace(b'G', b'c')\
              .replace(b'n', b'Z')\
              .upper()\
              .replace(b"Z", b"n")


# used by rawedit
def fullcomp(seq):
    """ returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences"""
    ## this is surely not the most efficient...
    seq = seq.replace("A", 'u')\
             .replace('T', 'v')\
             .replace('C', 'p')\
             .replace('G', 'z')\
             .replace('u', 'T')\
             .replace('v', 'A')\
             .replace('p', 'G')\
             .replace('z', 'C')

    ## No complement for S & W b/c complements are S & W, respectively
    seq = seq.replace('R', 'u')\
             .replace('K', 'v')\
             .replace('Y', 'b')\
             .replace('M', 'o')\
             .replace('u', 'Y')\
             .replace('v', 'M')\
             .replace('b', 'R')\
             .replace('o', 'K')

    seq = seq.replace('r', 'u')\
             .replace('k', 'v')\
             .replace('y', 'b')\
             .replace('m', 'o')\
             .replace('u', 'y')\
             .replace('v', 'm')\
             .replace('b', 'r')\
             .replace('o', 'k')
    return seq


# used by consens
## Alleles priority dict. The key:vals are the same as the AMBIGS dict
## except it returns just one base, w/ the order/priority being (C>A>T>G)
## This dict is used to impute lower case into consens to retain allele
## order for phase in diploids
PRIORITY = {
    b"M": b"C",
    b"Y": b"C",
    b"S": b"C",
    b"W": b"A",
    b"R": b"A",
    b"K": b"T",
}

# The inverse of priority
MINOR = {
    b"M": b"A",
    b"Y": b"T",
    b"S": b"G",
    b"W": b"T",
    b"R": b"G",
    b"K": b"G",
}

# used by write_outputs
# convert byte to list of alleles as ASCII strings
BTS = {
    b"R": ["G", "A"],
    b"K": ["G", "T"],
    b"S": ["G", "C"],
    b"Y": ["T", "C"],
    b"W": ["T", "A"],
    b"M": ["C", "A"],
    b"A": ["A", "A"],
    b"T": ["T", "T"],
    b"G": ["G", "G"],
    b"C": ["C", "C"],
    b"N": ["N", "N"],
    b"-": ["-", "-"]
    }

DUCT = {
    "R": ["G", "A"],
    "K": ["G", "T"],
    "S": ["G", "C"],
    "Y": ["T", "C"],
    "W": ["T", "A"],
    "M": ["C", "A"],
    "A": ["A", "A"],
    "T": ["T", "T"],
    "G": ["G", "G"],
    "C": ["C", "C"],
    "N": ["N", "N"],
    "-": ["-", "-"]
}    

# GETGENO = np.array([
#     list(b"RGA"),
#     list(b"KGT"),
#     list(b"SGC"),
#     list(b"YTC"),
#     list(b"WTA"),
#     list(b"MCA"),
#     list(b"AAA"),
#     list(b"TTT"),
#     list(b"CCC"),
#     list(b"GGG"),
# ], dtype=np.uint8)

# used in baba.py / write_outfiles..py
## with N and - masked to 9
GETCONS = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    [78, 9, 9],
    [45, 9, 9],
    ], dtype=np.uint8)


DCONS = {
    82: (71, 65),
    75: (71, 84),
    83: (71, 67),
    89: (84, 67),
    87: (84, 65),
    77: (67, 65),
    78: (9, 9),
    45: (9, 9),
    67: (67, 67),
    65: (65, 65),
    84: (84, 84),
    71: (71, 71),
}


def clustdealer(pairdealer, optim):
    """ 
    Return 'optim' clusters from 'pairdealer' iterators, and returns 
    1 if this includes the end of the iterator, else returns 0.
    """
    ccnt = 0
    chunk = []
    while ccnt < optim:
        # try refreshing taker, else quit
        try:
            taker = takewhile(lambda x: x[0] != b"//\n", pairdealer)
            oneclust = [b"".join(next(taker))]
        except StopIteration:
            return 1, chunk

        # load one cluster
        while 1:
            try:
                oneclust.append(b"".join(next(taker)))
            except StopIteration:
                break
        chunk.append(b"".join(oneclust))
        ccnt += 1
    return 0, chunk



def get_threaded_view(ipyclient, split=True):
    """ gets optimum threaded view of ids given the host setup """
    ## engine ids
    ## e.g., [0, 1, 2, 3, 4, 5, 6, 7, 8]
    eids = ipyclient.ids

    ## get host names
    ## e.g., ['a', 'a', 'b', 'b', 'a', 'c', 'c', 'c', 'c']
    dview = ipyclient.direct_view()
    hosts = dview.apply_sync(socket.gethostname)

    ## group ids into a dict by their hostnames
    ## e.g., {a: [0, 1, 4], b: [2, 3], c: [5, 6, 7, 8]}
    hostdict = {i: [] for i in hosts}  # defaultdict(list)
    for host, eid in zip(hosts, eids):
        hostdict[host].append(eid)

    ## Now split threads on the same host into separate proc if there are many
    hostdictkeys = list(hostdict.keys())
    for key in hostdictkeys:
        gids = hostdict[key]
        maxt = len(gids)
        if len(gids) >= 4:
            maxt = 2
        ## if 4 nodes and 4 ppn, put one sample per host
        if (len(gids) == 4) and (len(hosts) >= 4):
            maxt = 4
        if len(gids) >= 6:
            maxt = 3
        if len(gids) >= 8:
            maxt = 4
        if len(gids) >= 16:
            maxt = 4
        ## split ids into groups of maxt
        threaded = [gids[i:i + maxt] for i in range(0, len(gids), maxt)]
        lth = len(threaded)
        ## if anything was split (lth>1) update hostdict with new proc
        if lth > 1:
            hostdict.pop(key)
            for hostid in range(lth):
                hostdict[str(key) + "_" + str(hostid)] = threaded[hostid]

    ## make sure split numbering is correct
    #threaded = hostdict.values()
    #assert len(ipyclient.ids) <= len(list(itertools.chain(*threaded)))
    return hostdict


##############################################################
def detect_cpus():
    """
    Detects the number of CPUs on a system. This is better than asking
    ipyparallel since ipp has to wait for Engines to spin up.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.get("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else:  # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.get("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    return 1  # Default
