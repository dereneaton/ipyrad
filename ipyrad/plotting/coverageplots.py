#!/usr/bin/env python2

""" plot coverage depths for samples """

from __future__ import print_function
import toyplot
#import toyplot.html
#import toyplot.svg
#import toyplot.pdf
import numpy as np
from collections import OrderedDict


# pylint: disable=E1101


def depthplot(data, samples=None, dims=(None,None), canvas=(None,None), 
              xmax=50, log=False, outprefix=None, use_maxdepth=False):
    """ plots histogram of coverages across clusters"""

    ## select samples to be plotted, requires depths info
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
    if any(canvas):
        print("usercanvas")
        canvas = toyplot.Canvas(width=canvas[0], height=canvas[1])
    else:
        canvas = toyplot.Canvas(width=200*dims[1], height=150*dims[0])

    ## get all of the data arrays
    for panel, sample in enumerate(subsamples):
        ## statistical called bins
        statdat = subsamples[sample].depths
        statdat = statdat[statdat >= data.paramsdict["mindepth_statistical"]]
        if use_maxdepth:
            statdat = {i:j for (i, j) in statdat if \
                                        i < data.paramsdict["maxdepth"]}
        sdat = np.histogram(statdat, range(50))

        ## majrule called bins
        statdat = subsamples[sample].depths
        statdat = statdat[statdat < data.paramsdict["mindepth_statistical"]]
        statdat = statdat[statdat >= data.paramsdict["mindepth_majrule"]]
        if use_maxdepth:
            statdat = statdat[statdat < data.paramsdict["maxdepth"]]
        mdat = np.histogram(statdat, range(50))

        ## excluded bins
        tots = data.samples[sample].depths
        tots = tots[tots < data.paramsdict["mindepth_majrule"]]
        if use_maxdepth:
            tots = tots[tots < data.paramsdict["maxdepth"]]
        edat = np.histogram(tots, range(50))

        ## fill in each panel of canvas with a sample
        axes = canvas.cartesian(grid=(dims[0], dims[1], panel), gutter=25)
        axes.x.domain.xmax = xmax
        axes.label.text = sample
        if log:
            axes.y.scale = "log"

        # heights = np.column_stack((sdat,mdat,edat))
        axes.bars(sdat)
        axes.bars(edat)
        axes.bars(mdat)


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


