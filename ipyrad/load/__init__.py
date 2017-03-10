#!/usr/bin/env python2

""" The 'load' module contains functions for loading saved assembly objects
in JSON format. 

Functions 
----------
ip.load.load_assembly() 
ip.load_json()
ip.save_json()
ip.test_assembly()
"""

from .load import test_assembly
from .load import save_json
from .load import load_json
