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
#from .tetrad import Tetrad as tetrad
from .raxml import Raxml as raxml
from .structure import Structure as structure
from .tetrad import Tetrad as tetrad
from .bucky import Bucky as bucky

#from .baba import Baba as baba
from .bpp import Bpp as bpp
from .sratools import SRA as sratools
from .twiist import Twiist as twiist
from .treeslider import TreeSlider as treeslider
#from .popgen import Popgen as popgen
# from .treemix import Treemix as treemix
