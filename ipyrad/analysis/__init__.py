#!/usr/bin/env python

## version is the same as ipyrad
from ipyrad import __version__

## set interactive global. Turned off by CLI programs.
__interactive__ = 1

## analysis tools will have a class object that is upper case, which is 
## called by a convenience function which is lower case, and has the 
## same name as the module (file), such that when we import the function
## it clobbers the file name as a import. 


#from .tree import Tree as tree
#from toytree import tree
from .tetrad import Tetrad as tetrad
from .structure import Structure as structure
from .baba import Baba as baba
from .bpp import Bpp as bpp
from .raxml import Raxml as raxml

#from .treemix import treemix

#from ..plotting.baba_panel_plot import baba_panel_plot
#from ..plotting.tree_panel_plot import tree_panel_plot
#from ..plotting.share_panel_plot import share_panel_plot
#from .structure import structure



#import tetrad.tetrad as tetrad #
#import baba.baba as baba
#import baba
#from .tetrad.tetrad import tetrad

