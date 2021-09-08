#!/usr/bin/env python

"""
API level Classes are made available, binary paths are checked,
and the logger is activated in WARNING mode, but can be modified
by set_loglevel.
"""

# bring nested functions to top for API access
from ipyrad.core.assembly import Assembly
from ipyrad.core.logger_setup import set_log_level
from ipyrad.core.load_json import load_json
from ipyrad.core.merge import merge

__version__ = "1.0.0-alpha"
__author__ = "Deren Eaton & Isaac Overcast"

# configure the logger
set_log_level("WARNING")
