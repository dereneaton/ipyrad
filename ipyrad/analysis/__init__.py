#!/usr/bin/env python

## version is the same as ipyrad
from ipyrad import __version__


## analysis tools will have a class object that is upper case, which is 
## called by a convenience function which is lower case, and has the 
## same name as the module (file), such that when we import the function
## it clobbers the file name as a import. 

from .baba import Baba as baba
from .tree import Tree as tree
from .bpp import bpp
from .tetrad import tetrad
from .structure import structure
#from .treemix import treemix

#from ..plotting.baba_panel_plot import baba_panel_plot
#from ..plotting.tree_panel_plot import tree_panel_plot
#from ..plotting.share_panel_plot import share_panel_plot
#from .structure import structure



#import tetrad.tetrad as tetrad #
#import baba.baba as baba
#import baba
#from .tetrad.tetrad import tetrad

