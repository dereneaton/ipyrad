#!/usr/bin/env python

## version is the same as ipyrad
from ipyrad import __version__


## analysis tools will have a class object that is upper case, which is 
## called by a convenience function which is lower case, and has the 
## same name as the module (file), such that when we import the function
## it clobbers the file name as a import. 

from .bpp import bpp
from .tetrad import tetrad
from .baba import Tree as tree
from .baba import Baba as baba
#from .structure import structure
#from .treemix import treemix



#import tetrad.tetrad as tetrad #
#import baba.baba as baba
#import baba
#from .tetrad.tetrad import tetrad



