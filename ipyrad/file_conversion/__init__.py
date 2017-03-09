#!/usr/bin/env ipython2

""" mostly a collection of loci2x conversion scripts """

#from os.path import dirname, basename, isfile
#import glob
#MODULES = glob.glob(dirname(__file__)+"/*.py")
#__all__ = [basename(f)[:-3] for f in MODULES if isfile(f)]

from .loci2bpp import loci2bpp
from .loci2cf import loci2cf
#from .loci2migrate import loci2migrate

#from .loci2gphocs import loci2gphocs

#from . import loci2gphocs
#from . import loci2treemix
