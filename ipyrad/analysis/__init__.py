#!/usr/bin/env python

"""
ipyrad analysis tools -- a set of tools and tool-wrappers for 
evolutionary analysis of RAD-seq data that is simple, reproducible,
and aware of missing data.
"""

# version is the same as ipyrad
from ipyrad import __version__
from ipyrad.analysis.logger_setup import set_loglevel

# analysis tools will have a class object that is upper case, which is 
# called by a convenience function which is lower case, and has the 
# same name as the module (file), such that when we import the function
# it clobbers the file name as the import. 

# conda install mrbayes -c bioconda       # conda install py3.7 is messy...
# conda install raxml -c bioconda
# conda install fasttree -c bioconda      # needs OMP version support...
# conda install scikit-learn -c bioconda
# conda install sratools -c bioconda
# conda install treemix -c bioconda
# conda install toytree -c eaton-lab
# conda install tetrad -c eaton-lab
# conda install structure -c ipyrad
# conda install clumpp -c ipyrad
# conda install bucky -c ipyrad
# conda install bpp -c ipyrad

# tested in v.1.0 
from .raxml import Raxml as raxml
from .window_extracter import window_extracter
from .tree_slider import TreeSlider as tree_slider
from .digest_genome import DigestGenome as digest_genome
from .sratools import SRA as sratools
from .snps_imputer import SNPsImputer as snps_imputer
from .snps_extracter import SNPsExtracter as snps_extracter
from .download import Download as download
from .utils import popfile_to_imap
from .astral import Astral as astral
from .snaq import Snaq as snaq

# TESTING
from .pca import PCA as pca
# from .bucky import Bucky as bucky
# from .bpp import Bpp as bpp
# from .fasttree import Fasttree as fasttree


# analysis tools uses WARNING logger by default.
set_loglevel("INFO")





# # tested in 0.9.10

# from .mrbayes import MrBayes as mrbayes
# from .treemix import Treemix as treemix



# from .treeslider import TreeSlider as treeslider
# from .distance import Distance as distance
# from .structure import Structure as structure
# from .vcf_to_hdf5 import VCFtoHDF5 as vcf_to_hdf5  # tetrad version ahead.
#
# # testing
# from .locus_extracter import LocusExtracter as locus_extracter
# from .tetrad import Tetrad as tetrad
# from .window_extracter import WindowExtracter as window_extracter
# from .clade_weights import CladeWeights as clade_weights

# from .bucky import Bucky as bucky
# from .bpp import Bpp as bpp
# from .fasttree import Fasttree as fasttree
# from .baba import Baba as baba
# from .baba2 import Baba as baba2
# from .coverage import Coverage as coverage
# from .astral import Astral as astral
# from .snaq import Snaq as snaq

# from .momi import Momi as momi
# from .eems import Eems as eems
# from .popgen import Popgen as popgen
