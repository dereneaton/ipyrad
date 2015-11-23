#!/usr/bin/env ipython2

from . import assemble
from . import plotting
from . import dstats
#from . import core


## load the functional modules
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import get_params_info
from ipyrad.core.load_dataobj import load_assembly


__version__ = "0.0.66"


