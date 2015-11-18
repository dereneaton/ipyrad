#!/usr/bin/env ipython2

from . import assemble
from . import plotting
from . import dstats
#from . import core

## create an ipyparallel controller with an exit call to close at sys
## exit, and a unique cluster-id name
from ipyrad.core.parallel import ipcontroller_init
from ipyrad.core.parallel import ipcontroller_set

global __IPNAME__
ipcontroller_init()

## load the functional modules
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import get_params_info
from ipyrad.core.load_dataobj import load_assembly


__version__ = "0.0.66"


