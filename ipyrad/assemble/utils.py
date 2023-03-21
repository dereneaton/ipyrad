#!/usr/bin/env python

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

import os
import socket
import string
from itertools import takewhile
import pandas as pd
import numpy as np


# used in demux_raw.py and demux_sorted.py to fix sample names.
BADCHARS = (
    string.punctuation
    .replace("_", "")
    .replace("-", "")
    .replace(".", "") + " "
)

# used in demux_raw.py to resolve ambiguous cutters
AMBIGS = {
    "R": ("G", "A"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("T", "C"),
    "W": ("T", "A"),
    "M": ("C", "A"),
}

# used in consens_utils
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

# used in consens_utils to encode IUPAC ints 
TRANS = {
    (71, 65): 82,
    (71, 84): 75,
    (71, 67): 83,
    (84, 67): 89,
    (84, 65): 87,
    (67, 65): 77,
    (65, 67): 77,
    (65, 84): 87,
    (67, 84): 89,
    (67, 71): 83,
    (84, 71): 75,
    (65, 71): 82,
}


CIGARDICT = {
    '-': "I",
    '.': "S",
}


CDICT = dict(zip("CATG", "0123"))


class IPyradError(Exception):
    """Raise a custom exception that will report with traceback.

    This is used to catch and report internal errors in the code, 
    and the traceback will include the source error and error type
    for debugging.
    """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class NoHighDepthClustersError(IPyradError):
    """Raise a custom exception that will report with traceback.
    """


class IPyradExit(SystemExit):
    """Return code 1 to exit with an error message but NO TRACEBACK.

    This is used to catch and report common user mistakes to return
    an error message but not burden them with an ugly traceback.
    """

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

# used by step4.stackarray, step7.split_clusters, ...
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
            taker = takewhile(lambda x: x[0] != "//\n", pairdealer)
            oneclust = ["".join(next(taker))]
        except StopIteration:
            return 1, chunk

        # load one cluster
        while 1:
            try:
                oneclust.append("".join(next(taker)))
            except StopIteration:
                break
        chunk.append("".join(oneclust))
        ccnt += 1
    return 0, chunk

# used by cluster_across.build_clusters
def fullcomp(seq):
    """ 
    returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences
    """
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


# ================================================================
# above here has been checked.
# ================================================================


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


def chroms2ints(data, keys_as_ints: bool = False):
    """Parse .fai to get a dict with {chroms/scaffolds: ints}, or reversed.

    The chrom indices are 1-indexed.

    Parameters
    ----------
    data: Assembly
    intkeys: if True then return {int: name} else {name: int}
    """
    # load reference genome info as a dataframe.
    fai = pd.read_csv(
        str(data.params.reference_sequence) + ".fai",
        names=['scaffold', 'length', 'start', 'a', 'b'],
        sep="\t",
    )

    # get CHROM (scaffold names) as strings.
    fai["scaffold"] = fai["scaffold"].astype(str)

    # get dict mapping {str: int} using 1-indexed chrom indices.
    faidict = {j: i + 1 for i, j in enumerate(fai.scaffold)}

    # return dict as is, or reversed.
    if keys_as_ints:
        return {j: i for i, j in faidict.items()}
    return faidict


def get_fai_values(data, value):
    """Returns the fai table from the reference as an array."""
    reference_file = data.params.reference_sequence
    fai = pd.read_csv(
        reference_file + ".fai",
        names=['scaffold', 'length', 'sumsize', 'a', 'b'],
        sep="\t",
    )
    return fai[value].values


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



NEXHEADER = """#nexus
begin data;
  dimensions ntax={ntax} nchar={nchar};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""

NEXCLOSER = """  ;
end;
"""

STRDICT = {
    'A': '0', 
    'T': '1', 
    'G': '2', 
    'C': '3', 
    'N': '-9', 
    '-': '-9',
}

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
