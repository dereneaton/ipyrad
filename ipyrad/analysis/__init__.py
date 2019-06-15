#!/usr/bin/env python

# version is the same as ipyrad
from ipyrad import __version__

# analysis tools will have a class object that is upper case, which is 
# called by a convenience function which is lower case, and has the 
# same name as the module (file), such that when we import the function
# it clobbers the file name as the import. 

# conda install mrbayes -c bioconda
# conda install raxml -c bioconda
# conda install scikit-learn -c bioconda
# conda install sratools -c bioconda
# conda install treemix -c bioconda
# conda install toytree -c eaton-lab
# conda install tetrad -c eaton-lab
# conda install structure -c ipyrad
# conda install clumpp -c ipyrad
# conda install bucky -c ipyrad
# conda install bpp -c ipyrad

# tested in 0.9.10
from .raxml import Raxml as raxml
from .mrbayes import MrBayes as mrbayes
from .window_extracter import WindowExtracter as window_extracter
from .treemix import Treemix as treemix
from .pca import PCA as pca
from .tetrad import Tetrad as tetrad
from .treeslider import TreeSlider as treeslider
from .snps_extracter import SNPsExtracter as snps_extracter
from .snps_imputer import SNPsImputer as snps_imputer
from .distance import Distance as distance

# testing
from .structure import Structure as structure
from .clade_weights import CladeWeights as clade_weights
from .digest_genome import DigestGenome as digest_genome
# from .baba import Baba as baba
# from .sratools import SRA as sratools
# from .momi import Momi as momi
# from .eems import Eems as eems

# from .twisst import Twisst as twisst
from .bucky import Bucky as bucky
from .bpp import Bpp as bpp
#from .popgen import Popgen as popgen
