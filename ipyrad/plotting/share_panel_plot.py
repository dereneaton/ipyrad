#!/usr/bin/env python2

""" plot coverage depths for samples """

from __future__ import print_function


import toyplot
import toyplot.html
import toyplot.svg
import numpy as np
import itertools

import ipyrad.analysis as ipa
#from ipyrad.analysis import tree

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



## decided to use toyplot styling instead of ape-like args.
dwargs = {
    #"height": 600,
    #"width": 1000,
    "vsize": 0, 
    "vlshow": False, 
    "ewidth": 3, 
    "vlstyle": {"font-size": "18px"}, 
    "cex": "14px", 
    "pct_tree_y": 0.3, 
    "pct_tree_x": 0.7, 
    "lwd_lines": 1,
    }   



## setup panels
class Panel(object):
    def __init__(self, matrix, newick=None, counts=None):
        ## starting tree position will be changed if adding panel plots
        self.matrix = matrix
        if newick:
            self.tree = tree(newick)
            self.tree.verts = self.tree.verts.astype(np.float)
        else:
            self.tree = tree
        self.counts = counts

        ## spacing
        self.xmin_tree = 0.
        self.xmax_tree = 0.        

        self.xmin_mat = 0.
        self.xmax_mat = 0.
        self.ymin_mat = 0.

        self.xmin_count = 0.
        self.xmax_count = 0.


    @property
    def matrix_bounds(self):
        return (self.xmin_mat, self.xmax_mat, self.ymin_mat, self.ymax_mat)

    ## fixed properties
    @property
    def ymax_mat(self):
        return self.xmax_mat - self.xmin_mat



def share_matrix(locifile, tree=None, nameorder=None):
    """ 
    returns a matrix of shared RAD-seq data 
    
    Parameters:
    -----------
    locifile (str):
        Path to a ipyrad .loci file. 
    tree (str):
        Path to Newick file or a Newick string representation of
        a tree. If used, names will be ordered by the ladderized
        tip order. 
    nameorder (list):
        If a tree is not provided you can alternatively enter 
        the sample order as a list here. The tree argument will
        override this argument.

    Returns:
    --------
    matrix (numpy.array):
        A uint64 numpy array of the number of shared loci between
        all pairs of samples.
    """

    ## load in the loci data
    with open(locifile, 'r') as locidata:
        loci = locidata.read().split("|\n")[:-1]

    ## load in the tree from a string
    if tree:
        tree = ete.Tree(tree)
        tree.ladderize()
        snames = tree.get_leaf_names()
        lxs, names = _getarray(loci, snames)
    elif nameorder:
        lxs, names = _getarray(loci, nameorder)
    else:
        raise IOError("must provide either tree or nameorder argument")

    ## get share matrix
    share = _countmatrix(lxs)

    return share



def share_panel_plot(
    matrix,
    newick,
    counts=None,
    *args, 
    **kwargs):

    ## update kwargs
    dwargs.update(kwargs)

    ## parse the tree
    panel = Panel(newick, matrix, counts)

    ## setup the canvas
    canvas = toyplot.Canvas(height=dwargs['height'], width=dwargs['width'])
    axes = canvas.cartesian(bounds=("5%", "95%", "5%", "95%"))
    axes.show = False

    panel = _panel_tree(axes, panel, dwargs)
    panel = _panel_matrix(axes, panel, dwargs)

    ## returning will auto display on a notebook if not stored
    return canvas, axes



def _panel_tree(axes, panel, dwargs):
   
    ## adjust X-axis: boots X is 75% of tree X
    panel.xmax_tree += panel.tree.verts[:, 0].min() * dwargs["pct_tree_x"]
    panel.tree.verts[:, 0] += panel.xmax_tree

    ## get spacer between panels
    #panel.ztot = 0.15 * panel.xmin_tree

    ## add spacer for tip names
    #pctt = 0.2 * (panel.ymin_tree + panel.ymin_test)
    #panel.ymin_tree += pctt / 2.
    #panel.ymin_test += pctt / 2.
    #panel.tree.verts[:, 1] += pctt / 2.
            
    ## add the tree/graph ------------------------------------------------
    _ = axes.graph(panel.tree.edges, 
                   vcoordinates=panel.tree.verts, 
                   ewidth=dwargs["ewidth"], 
                   ecolor=toyplot.color.near_black, 
                   vlshow=dwargs["vlshow"],
                   vsize=dwargs["vsize"],
                   vlstyle=dwargs["vlstyle"],
                   vlabel=panel.tree.names) #.values())   
    return panel





def _plotshare(share, names, **kwargs):
    """ make toyplot matrix fig"""

    ## set the colormap
    colormap = toyplot.color.LinearMap(
        toyplot.color.brewer.palette("Spectral"), 
        domain_min=share.min(), 
        domain_max=share.max())

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



def _getarray(loci, snames):
    """ 
    parse loci list and return presence/absence matrix
    ordered by the tips on the tree or list of names.
    """
    ## make an empty matrix
    lxs = np.zeros((len(snames), len(loci)), dtype=np.uint64)

    ## fill the matrix
    for loc in xrange(len(loci)):
        for seq in loci[loc].split("\n"):
            if "//" not in seq:
                lxs[snames.index(seq.split()[0][:]), loc] += 1

    return lxs, snames



def _countmatrix(lxs):
    """ fill a matrix with pairwise data sharing """
    
    ## an empty matrix
    share = np.zeros((lxs.shape[0], lxs.shape[0]), dtype=np.uint64)

    ## fill above
    names = range(lxs.shape[0])
    for row in lxs:
        for samp1, samp2 in itertools.combinations(names, 2):
            shared = lxs[samp1, lxs[samp2] > 0].sum()
            share[samp1, samp2] = shared

    ## mirror below
    share += share.T

    ## fill diagonal with total sample coverage
    for row in xrange(len(names)):
        share[row, row] = lxs[row].sum()

    return share




ETE_IMPORT_ERROR = """\
  Missing requirement = ete3
  Please run `conda install -c etetoolkit ete3` to install.
"""