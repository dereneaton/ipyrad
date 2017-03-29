#!/usr/bin/env python2

#from . import coverageplots
#from ipyrad.plotting.coverageplots import depthplot
#from ipyrad.plotting.shareplot import shareplot

## we want to have a module with prep functions
## import ipyrad.plotting as iplot


## iplot.share_matrix_plot(share)
## iplot.share_panel_plot(newick, share)
## iplot.share_panel_plot(newick, share, data)


## iplot.baba_panel_plot(newick, tests, boots, show_tips=True)
## iplot.abbba_panel_plot(newick, tests, boots, show_tips=True)

#import shareplot.shareplot as shareplot
from .baba_panel_plot import baba_panel_plot


