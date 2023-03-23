#!/usr/bin/env python

"""API level classes for ipyrad Assembly.

Examples
--------
>>> import ipyrad as ip
>>> data = ip.Assembly("test")
>>> data.params.sorted_fastq_path = "..."
>>> data.run('1', cores=4)

>>> data1 = data.branch("test2", subsamples=['a', 'b'])
>>> data1.run('2', cores=4)
>>> data2 = data.branch("test2", subsamples=['a', 'b'])
>>> data2.run('2', cores=4)

>>> fulldata = data1.merge(data2)
>>> fulldata.run('3', cores=4)
"""

# bring nested functions to top for API access
from ipyrad.core.assembly import Assembly
from ipyrad.core.logger_setup import set_log_level
from ipyrad.core.load_json import load_json
from ipyrad.core.merge import merge
from ipyrad.core.cluster import Cluster

__version__ = "1.0.0-alpha-2"
__author__ = "Deren Eaton & Isaac Overcast"

# configure the logger
set_log_level("INFO")
