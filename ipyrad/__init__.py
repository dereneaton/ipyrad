#!/usr/bin/env python

# bring nested functions to top for API access
from .core.assembly import Assembly, merge
from .core.parallel import get_client as _get_client, cluster_info
from .core.startup import Bins as _Bins
from .core.load import load_json

# Dunders
__version__ = "0.8.0-dev"
__author__ = "Deren Eaton & Isaac Overcast"

# CLI __main__ changes to 0
__interactive__ = 1

# log file
__debugfile__ = "./ipyrad_log.txt"
__debugflag__ = "./.debug"

# check binaries
bins = _Bins()

# check hard installs
try:
    import pysam
except ImportError:
    print("""
You must first install 'pysam' with either conda or pip, e.g.,: 

    conda install bioconda::pysam

    or 

    pip install pysam
""")
