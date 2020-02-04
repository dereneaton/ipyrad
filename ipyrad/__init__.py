#!/usr/bin/env python

# bring nested functions to top for API access
from .core.assembly import Assembly, merge
from .core.load import load_json
from .core.Parallel import cluster_info
import os as _os
import sys as _sys
import subprocess as _sps

# Dunders
__version__ = "0.9.32"
__author__ = "Deren Eaton & Isaac Overcast"

# CLI __main__ changes to 0
__interactive__ = 1

# get binaries from conda/bin or conda/env/bin
class _Bins:
    pass


_IMPORT_ERROR = """
Missing requirement: {}

Please run 'conda install {} -c bioconda' or to install
all requirements run 'conda upgrade ipyrad -c bioconda'.
"""

# check binaries
bins = _Bins()
bins.muscle = _os.path.join(_sys.prefix, "bin", "muscle")
bins.samtools = _os.path.join(_sys.prefix, "bin", "samtools")
bins.bedtools = _os.path.join(_sys.prefix, "bin", "bedtools")
bins.vsearch = _os.path.join(_sys.prefix, "bin", "vsearch")
bins.bwa = _os.path.join(_sys.prefix, "bin", "bwa")

for binary, path in bins.__dict__.items():

    # check for conda version
    if not _os.path.exists(path):
        setattr(bins, binary, binary)

        # if not then check for binary in PATH (less reliable versioned...)
        cmd = ['which', binary]
        proc = _sps.Popen(cmd, stderr=_sps.STDOUT, stdout=_sps.PIPE)
        errmsg = proc.communicate()[0]
        if proc.returncode:
            print(errmsg.decode())
            raise ImportError(_IMPORT_ERROR.format(binary, binary))


# if user installed with pip then the following may be missing:
try:
    import pysam
except ImportError:
    print("""
You must first install 'pysam' with either conda or pip, e.g.,: 

    conda install pysam -c bioconda

    or 

    pip install pysam
""")
