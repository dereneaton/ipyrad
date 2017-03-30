#!/usr/bin/env python2

""" plot coverage depths for samples """

from __future__ import print_function
import toyplot
import toyplot.html
import toyplot.svg
import numpy as np
import itertools

try:
    import toyplot.pdf
except ImportError:
    pass

try:
    import ete3 as ete
except ImportError:
    try:
        import ete2 as ete
    except ImportError:
        raise ImportError(ETE_IMPORT_ERROR)


# pylint: disable=E1101

#samples=None, dims=(None,None), canvas=(None,None), 
#              xmax=50, log=False, outprefix=None, use_maxdepth=False)

def shareplot(locifile, tree, **kwargs):

    ## load in the loci data
    with open(locifile, 'r') as locidata:
        loci = locidata.read().split("|\n")[:-1]

    ## load in the tree from a string
    tree = ete.Tree(tree)

    ## get array
    lxs, names = _getarray(loci, tree)

    ## get share matrix
    share = _countmatrix(lxs)

    ## plot the figure
    canvas, axes = _plotshare(share, names, **kwargs)

    ## returning will auto display on a notebook if not saved to var
    return canvas, axes



def _plotshare(share, names, **kwargs):
    """ make toyplot matrix fig"""

    ## set the colormap
    colormap = toyplot.color.LinearMap(toyplot.color.brewer.palette("Spectral"), 
                                 domain_min=share.min(), domain_max=share.max())

    ## set up canvas
    if not kwargs.get('width'):
        width=900
    else:
        width = kwargs['width']
    canvas = toyplot.Canvas(width=width, height=width*0.77778)

    ## order the dta
    table = canvas.matrix((share, colormap), 
                          bounds=(50, canvas.height-100,
                                  50, canvas.height-100), 
                          step=5, tshow=False, lshow=False)

    ## put a box around the table
    table.body.grid.vlines[..., [0, -1]] = 'single'
    table.body.grid.hlines[[0, -1], ...] = 'single'

    ## make hover info on grid
    for i, j in itertools.product(range(len(share)), repeat=2):
        table.body.cell(i,j).title = "%s, %s : %s" % (names[i], names[j], int(share[i,j]))

    ## create barplot
    axes = canvas.cartesian(bounds=(665, 800, 90, 560))

    ## make a hover for barplot
    zf = zip(names[::-1], share.diagonal()[::-1])
    barfloater = ["%s: %s" % (i, int(j)) for i, j in zf]

    ## plot bars
    axes.bars(share.diagonal()[::-1], along='y', title=barfloater)

    ## hide spine, move labels to the left, 
    ## use taxon names, rotate angle, align
    axes.y.spine.show = False
    axes.y.ticks.labels.offset = 0
    axes.y.ticks.locator = toyplot.locator.Explicit(range(len(names)), 
                                                    labels=names[::-1])
    axes.y.ticks.labels.angle = -90
    axes.y.ticks.labels.style = {"baseline-shift":0,
                                 "text-anchor":"end", 
                                 "font-size":"8px"}

    ## rotate xlabels, align with ticks, change to thousands, move up on canvas
    ## show ticks, and hide popup coordinates
    axes.x.ticks.labels.angle = 90
    axes.x.ticks.labels.offset = 20
    axes.x.ticks.locator = toyplot.locator.Explicit(
        range(0, int(share.max()), 
                 int(share.max() / 10)), 
        ["{}".format(i) for i in range(0, int(share.max()),
                                           int(share.max() / 10))])
    axes.x.ticks.labels.style = {"baseline-shift":0,
                                 "text-anchor":"end", 
                                 "-toyplot-anchor-shift":"15px"}
    axes.x.ticks.show = True

    ## add labels
    label_style = {"font-size": "16px", "font-weight": "bold"}
    canvas.text(300, 60, "Matrix of shared RAD loci", style=label_style)
    canvas.text(700, 60, "N RAD loci per sample", style=label_style)

    return canvas, axes



def _getarray(loci, tree):
    """ 
    parse the loci file list and return presence/absence matrix
    ordered by the tips on the tree
    """

    ## order tips
    tree.ladderize()

    ## get tip names
    snames = tree.get_leaf_names()

    ## make an empty matrix
    lxs = np.zeros((len(snames), len(loci)), dtype=np.int)

    ## fill the matrix
    for loc in xrange(len(loci)):
        for seq in loci[loc].split("\n")[:-1]:
            lxs[snames.index(seq.split()[0]), loc] += 1

    return lxs, snames



def _countmatrix(lxs):
    """ fill a matrix with pairwise data sharing """
    
    ## an empty matrix
    share = np.zeros((lxs.shape[0], lxs.shape[0]))

    ## fill above
    names = range(lxs.shape[0])
    for row in lxs:
        for samp1, samp2 in itertools.combinations(names, 2):
            shared = lxs[samp1, lxs[samp2] > 0].sum()
            share[samp1, samp2] = shared

    ## mirror below
    ##share[]

    ## fill diagonal with total sample coverage
    for row in xrange(len(names)):
        share[row, row] = lxs[row].sum()

    return share


ETE_IMPORT_ERROR = """\
  Missing requirement = ete3
  Please run `conda install -c etetoolkit ete3` to install.
"""