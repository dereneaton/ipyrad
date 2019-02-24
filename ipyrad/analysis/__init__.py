#!/usr/bin/env python

## version is the same as ipyrad
from ipyrad import __version__

## set interactive global. Turned off by CLI programs.
__interactive__ = 1

## analysis tools will have a class object that is upper case, which is 
## called by a convenience function which is lower case, and has the 
## same name as the module (file), such that when we import the function
## it clobbers the file name as the import. 

# conda install raxml -c bioconda
from .raxml import Raxml as raxml
from .twiist import Twisst as twisst
from .treeslider import TreeSlider as treeslider

# conda install structure clumpp -c ipyrad
from .structure import Structure as structure

# conda install tetrad -c eaton-lab
from .tetrad import Tetrad as tetrad

# conda install bucky -c ipyrad
# conda install mrbayes -c BioBuilds
from .bucky import Bucky as bucky

# conda install bpp -c ipyrad
from .bpp import Bpp as bpp

# conda install sratools -c bioconda
from .sratools import SRA as sratools

#from .baba import Baba as baba

#from .popgen import Popgen as popgen
#from .treemix import Treemix as treemix
