#!/usr/bin/env python2

""" plot coverage depths for samples """

from __future__ import print_function
import toyplot
import toyplot.html
import toyplot.svg
#import toyplot.pdf
import numpy as np
from collections import OrderedDict


# pylint: disable=E1101


def depthplot(data, samples=None, dims=(0, 0), xmax=50, outprefix=None):
    """ plots histogram of coverages across clusters"""

    ## select samples to be plotted
    if not samples:
        samples = data.samples.keys()
        samples.sort()
    subsamples = OrderedDict([(i, data.samples[i]) for i in samples])

    ## get canvas dimensions based on n-samples
    if any(dims):
        ## user-supplied dimensions (...)
        print("userdims")
    else:
        if len(subsamples) <= 4:
            ## set dimension to N samples 
            dims = (1, len(subsamples))
        else:
            dims = (len(subsamples)/4, 4)

    ## create canvas
    canvas = toyplot.Canvas(width=250*dims[1], height=200*dims[0])

    ## fill in each panel of canvas with a sample
    for panel, sample in enumerate(subsamples):
        axes = canvas.axes(grid=(dims[0], dims[1], panel), gutter=20)
        axes.x.domain.xmax = xmax
        axes.label.text = sample

        ## statistical called bins
        statdat = subsamples[sample].depths["statmin"]
        sdat = np.histogram(statdat, range(50))

        ## majrule called bins
        maj = data.samples[sample].depths["mjmin"]
        majdat = maj[maj < data.paramsdict["mindepth_statistical"]]
        mdat = np.histogram(majdat, range(50))

        ## excluded bins
        tots = data.samples[sample].depths["total"]
        totsdat = tots[tots < data.paramsdict["mindepth_majrule"]]
        edat = np.histogram(totsdat, range(50))

        # ## set ymax using highest bin...
        # #ymax = ...
        # heights = np.column_stack((sdat,mdat,edat))
        axes.bars(sdat)
        axes.bars(mdat)
        axes.bars(edat)

    ## return objects to be saved...
    if outprefix:
        toyplot.html.render(canvas, fobj=outprefix+".html")
        toyplot.svg.render(canvas, fobj=outprefix+".svg")
        #toyplot.pdf.render(canvas, fobj=outprefix+".pdf")


if __name__ == "__main__":
    pass

#   import ipyrad as ip
#   import numpy as np
#   TEST = ip.Sample()
#   TEST.depths["total"] = np.random.normal(20, 3, 2000)

#   #DATA = ip.Assemble()
#   #DATA.link_samples(TEST)

#   coverageplot(TEST)


