#!/usr/bin/env python2

"""
a panel plot function for baba results 
"""

#from ipyrad.analysis import tree # import Tree as tree
from __future__ import print_function

import numpy as np
#import ete3 as ete
import toyplot
import itertools


class Panel(object):
    def __init__(self, tree, edges, verts, names, tests, boots):

        ## starting tree position will be changed if adding panel plots
        self.tree = tree
        self.edges = edges
        self.verts = verts
        self.names = names
        self.tests = tests
        self.boots = np.array(boots)
        
        self.xmin_tree = 0.
        self.xmin_hist = 0.
        
        self.ymin_tree = 0.
        self.ymin_test = 0.
        self.ymin_text = 0.        
        

        self.kwargs = {
            "vsize": 0, 
            "vlshow": False, 
            "ewidth": 3,             
            "pct_tree_y": 0.3, 
            "pct_tree_x": 0.7, 
            "lwd_lines": 1,
            #"cex": "14px",             
            "vlstyle": {"font-size": "18px"}, 
            "style_tip_labels": {
                    "font-size": "14px",
                    "text-anchor" : "start",
                    "-toyplot-anchor-shift":"0",
                    },  
        }


    ## auto-generated properties
    @property
    def ztot(self):
        return 0.15 * self.xmin_tree

    @property    
    def xmin_zs(self):
        return self.xmin_tree * 0.5 
    
    @property
    def yhs(self):
        return np.linspace(self.ymin_tree, self.ymin_test, len(self.tests)+1)
    @property
    def ysp(self):
        return (self.yhs[0] - self.yhs[1])  ## the distance between tests
    @property
    def ytrim(self):
        return 0.70 * self.ysp  ## calculate coords, top&bot 30% of bars cutoff 
    
    @property
    def zright(self):
        return self.xmin_tree - self.ztot
    @property
    def zleft(self):
        return self.xmin_zs + self.ztot / 2.



    ## plotting functions
    def panel_tree(self, axes):
        _panel_tree(self, axes)

    def panel_rect(self, axes):
        _rectangle(self, axes)

    def panel_tip_labels(self, axes):
        _panel_tip_labels(self, axes)

    def panel_test_labels(self, axes):
        pass

    def panel_test_rects(self, axes):
        pass



    
def _panel_rect(panel, axes):

    ## hard-coded color palette for tips
    cdict = {"p1": toyplot.color.Palette()[0], 
             "p2": toyplot.color.Palette()[1],
             "p3": toyplot.color.near_black, 
             "p4": toyplot.color.Palette()[-1],}
    
    ## iterate over tests putting in rects
    for idx, test in enumerate(panel.tests):
        ## check taxdict names
        dictnames = list(itertools.chain(*test.values()))
        badnames = [i for i in dictnames if i not in panel.tree.names]
        if badnames:
            #raise IPyradError(
            print("Warning: found names not in tree:\n -{}"\
                  .format("\n -".join(list(badnames))))

        ## add dashed grid line for tests
        gedges = np.array([[0, 1], [2, 3]])
        gverts = np.array([
                    [panel.xmin_tree, panel.yhs[idx] - panel.ysp/2.],
                    [panel.tree.verts[:, 0].max(), panel.yhs[idx] - panel.ysp/2.],
                    [panel.xmin_zs + panel.ztot/2., panel.yhs[idx] - panel.ysp/2.],
                    [panel.xmin_tree - panel.ztot, panel.yhs[idx] - panel.ysp/2.],
                    ])
        axes.graph(gedges, 
                   vcoordinates=gverts,
                   ewidth=dwargs["lwd_lines"], 
                   ecolor=toyplot.color.Palette()[-1],
                   vsize=0, 
                   vlshow=False,
                   estyle={"stroke-dasharray": "5, 5"},
                   )                              

        ## add test rectangles
        for tax in ["p1", "p2", "p3", "p4"]:
            spx = [panel.tree.tree.search_nodes(name=i)[0] for i in test[tax]]
            spx = np.array([panel.tree.verts[i.idx, 0] for i in spx])
            #spx = np.array([i.x for i in spx], dtype=np.float)
            #spx += panel.xmin_tree
            spx.sort()
            
            ## fill rectangles while leaving holes for missing taxa
            for i in xrange(spx.shape[0]):
                if i == 0:
                    xleft = spx[i] - 0.25
                    xright = spx[i] + 0.25
                if i != spx.shape[0]-1:
                    if spx[i+1] - spx[i] < 1.5:
                        xright += 1
                    else:
                        axes.rects(xleft, xright, 
                            panel.yhs[idx] - panel.ytrim, 
                            panel.yhs[idx+1] + panel.ytrim, color=cdict[tax]) 
                        xleft = spx[i+1] - 0.25
                        xright = spx[i+1] + 0.25
                else:
                    axes.rects(xleft, xright,
                            panel.yhs[idx] - panel.ytrim, 
                            panel.yhs[idx+1] + panel.ytrim, color=cdict[tax]) 
    return panel



def _rectangle(panel, axes):

    ## hard-coded color palette for tips
    cdict = {"p1": toyplot.color.Palette()[0], 
             "p2": toyplot.color.Palette()[1],
             "p3": toyplot.color.near_black, 
             "p4": toyplot.color.Palette()[-1],}

    ## space for tests (pct_tree_y will scale tests space to tree space)
    panel.xmin_test = panel.xmin_tree
    panel.xmax_test = panel.xmax_tree
    panel.ymin_test = panel.ymax_tree - (panel.ymax_tree / float(panel.kwargs["pct_tree_y"]))
    panel.ymax_test = panel.ymin_tree 

    ## debug the panel
    debug = 1
    if debug:
        print('x_test', panel.xmin_test, panel.xmax_test)
        print('y_test', panel.ymin_test, panel.ymax_test)
        axes.rects(
            panel.xmin_test,
            panel.xmax_test,
            panel.ymin_test,
            panel.ymax_test,
            color=toyplot.color.Palette()[3], 
            opacity=0.5,
            )


def null():    
    ## iterate over tests putting in rects
    for idx, test in enumerate(panel.tests):
        ## check taxdict names
        dictnames = list(itertools.chain(*panel.test.values()))
        badnames = [i for i in dictnames if i not in panel.names.values()]
        if badnames:
            print("Warning: found names not in tree:\n -{}"\
                  .format("\n -".join(list(badnames))))

        ## add dashed grid line for tests
        gedges = np.array([[0, 1], [2, 3]])
        gverts = np.array([
                    [panel.xmin_tree, panel.yhs[idx] - panel.ysp/2.],
                    [panel.tree.verts[:, 0].max(), panel.yhs[idx] - panel.ysp/2.],
                    [panel.xmin_zs + panel.ztot/2., panel.yhs[idx] - panel.ysp/2.],
                    [panel.xmin_tree - panel.ztot, panel.yhs[idx] - panel.ysp/2.],
                    ])
        axes.graph(gedges, 
                   vcoordinates=gverts,
                   ewidth=dwargs["lwd_lines"], 
                   ecolor=toyplot.color.Palette()[-1],
                   vsize=0, 
                   vlshow=False,
                   estyle={"stroke-dasharray": "5, 5"},
                   )               




def _panel_tip_labels(panel, axes):
    ## calculate coords
    names = [i for i in panel.names.values() if not isinstance(i, int)]
    #spx = [panel.tree.search_nodes(name=i)[0] for i in names]
    #spx = np.array([panel.verts[i.idx, 0] for i in spx])

    ## add to axes
    xpos = panel.verts[:, 0][panel.verts[:, 1] == panel.verts[:, 1].min()]
    ypos = [panel.ymax_tree] * len(names)
    print(panel.verts)
    print('x', xpos)
    print('y', ypos)
    print('n', names)
    _ = axes.text(
                xpos, 
                ypos,
                names,
                angle=-90, 
                color=toyplot.color.near_black,
                style=panel.kwargs["style_tip_labels"],
                ) 



def _panel_tree(panel, axes):
    ## add the tree/graph ---------------------------
    panel.xmin_tree = panel.verts[:, 0].min()    
    panel.xmax_tree = panel.verts[:, 0].max()
    panel.ymin_tree = panel.verts[:, 1].min()
    panel.ymax_tree = panel.verts[:, 1].max()

    ## debug the panel
    debug = 1
    if debug:
        print('x_tree', panel.xmin_tree, panel.xmax_tree)
        print('y_tree', panel.ymin_tree, panel.ymax_tree)
        axes.rects(
            panel.xmin_tree,
            panel.xmax_tree,
            panel.ymin_tree,
            panel.ymax_tree,
            color=toyplot.color.Palette()[2], 
            opacity=0.5,
            )

    ## add the graph
    _ = axes.graph(panel.edges, 
                   vcoordinates=panel.verts, 
                   ewidth=panel.kwargs["ewidth"], 
                   ecolor=toyplot.color.near_black, 
                   vlshow=panel.kwargs["vlshow"],
                   vsize=panel.kwargs["vsize"],
                   vlstyle=panel.kwargs["vlstyle"],
                   vlabel=panel.names.keys(),
                   )



## the main function.
def baba_panel_plot(
    tree,
    edges, 
    verts, 
    names,
    tests, 
    boots,
    show_tip_labels=True, 
    show_test_labels=True, 
    use_edge_lengths=False, 
    collapse_outgroup=False, 
    pct_tree_x=0.4, 
    pct_tree_y=0.2,
    *args, 
    **kwargs):
    """
    signature...
    """

    ## create Panel plot object and set height & width
    panel = Panel(tree, edges, verts, names, tests, boots)
    if not kwargs.get("width"):
        panel.kwargs["width"] = min(1000, 50*len(panel.tree))
    if not kwargs.get("height"):
        panel.kwargs["height"] = min(1000, 50*len(panel.tests))   

    ## update defaults with kwargs & update size based on ntips & ntests
    kwargs.update(dict(pct_tree_x=pct_tree_x, pct_tree_y=pct_tree_y))
    panel.kwargs.update(kwargs)

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=panel.kwargs['height'], width=panel.kwargs['width'])
    axes = canvas.cartesian(bounds=("10%", "90%", "10%", "90%"))    
    axes.show = False

    ## partition panels for pct_x pct_y
    ## set size of tree
    panel.tests_height = 0.
    panel.xmin_tree = panel.kwargs["width"] - (panel.kwargs["width"] * pct_tree_x)
    panel.ymin_tree = panel.kwargs["height"] * pct_tree_y
    ## set remaining y-height with tests + names
    panel.ymin_ = 0.


    ## add panels to axes
    panel.panel_tree(axes)
    panel.panel_rect(axes)
    panel.panel_tip_labels(axes)
    #panel = _panel_tree(axes, panel, dwargs)
    #panel = _panel_rect(axes, panel, dwargs)
    #panel = _panel_hist(axes, panel, dwargs)
    #panel = _panel_test_labels(axes, panel, show_test_labels)
    #panel = _panel_tree_labels(axes, panel, show_tip_labels)

    return canvas, axes
    


    
    
def _panel_hist(axes, panel, dwargs):
    ## get bounds on hist
    rmax = np.max(np.abs(panel.boots))
    rmax = round(min(1.0, max(0.2, rmax + 0.05*rmax)), 1)
    rmin = round(max(-1.0, min(-0.2, -1 * rmax)), 1)
    allzs = np.abs(np.array([i.mean() / i.std() for i in panel.boots]))
    zmax = max(3., float(np.math.ceil(allzs.max())))

    ## add histograms, and hist axes
    for idx in xrange(panel.boots.shape[0]):
        bins = 30
        mags, xpos = np.histogram(panel.boots[idx], 
                                  bins=bins, 
                                  range=(rmin, rmax), 
                                  density=True)
        ## get hist colors
        thisz = allzs[idx]
        if thisz > 3:
            if panel.boots[idx].mean() < 0:
                color = toyplot.color.Palette()[0]
            else:
                color = toyplot.color.Palette()[1]
        else:
            color = toyplot.color.Palette()[-1]

        ## plot z's within range 
        #zright = panel.xmin_tree - panel.ztot
        #zleft = panel.xmin_zs + panel.ztot/2.
        zprop = thisz / zmax
        zmaxlen = panel.zright - panel.zleft
        zlen = zmaxlen * zprop
        axes.rects(
            panel.zright - zlen, panel.zright,
            panel.yhs[idx] - panel.ytrim, 
            panel.yhs[idx+1] + panel.ytrim,
            color=toyplot.color.Palette()[2])

        ## get hist xspans, heights in range
        xran = np.linspace(panel.xmin_hist, 
                           panel.xmin_zs - panel.ztot/2.,
                           bins+1)
        mags = mags / mags.max()
        mags = (mags * panel.ysp) * 0.8
        yline = panel.yhs[idx+1]
        heights = np.column_stack([[yline for i in mags], mags])
        axes.bars(
            xran[:-1], 
            heights, 
            baseline="stacked",
            color=["white", color],
            )

        ## add x-line to histograms
        gverts = np.array([
            [xran[0], yline], 
            [xran[-1], yline],
            ])
        gedges = np.array([[0, 1]])
        axes.graph(
            gedges, 
            vcoordinates=gverts, 
            ewidth=1, 
            ecolor=toyplot.color.near_black,
            vsize=0, 
            vlshow=False,
            )

    ## add xline to z-plots
    gverts = np.array([
        [panel.xmin_zs + panel.ztot/2., yline], 
        [panel.xmin_tree - panel.ztot, yline], 
        ])
    gedges = np.array([[0, 1]])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        vsize=0, 
        vlshow=False,
        )

    ## add dashed y-line at 0 to histograms
    here = np.where(xpos > 0)[0].min()
    zero = xran[here-1]
    gverts = np.array([
        [zero, panel.ymin_test],
        [zero, panel.ymin_tree - panel.ytrim / 4.],
        ])
    gedges = np.array([[0, 1]])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        eopacity=1.,
        vsize=0, 
        vlshow=False,
        estyle={"stroke-dasharray": "5, 5"},
        )

    ## add solid y-line to z-plots
    gverts = np.array([
        [panel.xmin_tree - panel.ztot, panel.ymin_test + panel.ytrim/4.],
        [panel.xmin_tree - panel.ztot, panel.ymin_tree - panel.ytrim/4.],
        ])
    gedges = np.array([[0, 1]])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        eopacity=1.,
        vsize=0, 
        vlshow=False,
        )

    ## add tick-marks to x-lines (hists and z-plots)
    ticklen = panel.ytrim / 4.
    gedges = np.array([[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]])
    gverts = np.array([
        [panel.xmin_hist, panel.ymin_test - ticklen],
        [panel.xmin_hist, panel.ymin_test],
        [zero, panel.ymin_test - ticklen],                              
        [zero, panel.ymin_test],
        [panel.xmin_zs - panel.ztot / 2., panel.ymin_test - ticklen],
        [panel.xmin_zs - panel.ztot / 2., panel.ymin_test],
        [panel.xmin_zs + panel.ztot / 2., panel.ymin_test - ticklen],
        [panel.xmin_zs + panel.ztot / 2., panel.ymin_test],
        [panel.xmin_tree - panel.ztot, panel.ymin_test - ticklen],
        [panel.xmin_tree - panel.ztot, panel.ymin_test],
        ])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        eopacity=1.,
        vsize=0, 
        vlshow=False,
        )

    ## add tick labels
    labels = [rmin, 0, rmax, zmax, 0]
    axes.text(
        gverts[:, 0][::2], 
        [panel.ymin_test - panel.ysp] * len(labels),
        labels, 
        color=toyplot.color.near_black,
        style={
            "font-size": dwargs["cex"],
            "text-anchor" : "middle",
            "-toyplot-anchor-shift":"0",
            },
        )

    ## add baba abba labels
    axes.text(
        [gverts[:, 0][0], gverts[:, 0][4]],
        [panel.ymin_test - panel.ysp * 2] * 2,
        ["BABA", "ABBA"], 
        color=toyplot.color.near_black,
        style={
            "font-size": dwargs["cex"], #"12px",
            "text-anchor" : "middle",
            "-toyplot-anchor-shift":"0",
            },
        )     

    ## add bootstrap and z-score titles
    axes.text(
        [zero, panel.zright - 0.5 * zmaxlen],
        [panel.ymin_tree + panel.ysp / 2.] * 2,
        ["Bootstrap D-statistics", "Z-scores"], 
        color=toyplot.color.near_black,
        style={
            "font-size": dwargs["cex"], #"12px",
            "text-anchor" : "middle",
            "-toyplot-anchor-shift":"0",
            },
        )
    return panel
    

   
def _panel_test_labels(axes, panel, show_test_labels):
    ## add test-label
    if show_test_labels:
        if isinstance(show_test_labels, bool) or show_test_labels == 1: 
            labels = range(1, len(panel.tests) + 1)
        elif isinstance(show_test_labels, list):
            labels = test_labels
        else:
            raise IPyradError("  label_tests must be a list or boolean")
        axes.text(
            [panel.tree.verts[:, 0].max() + 1] * len(panel.tests),
            panel.yhs[:-1] - panel.ysp / 2., 
            labels,
                color=toyplot.color.near_black,
                style={
                    "font-size": dwargs["cex"],
                    "text-anchor" : "start",
                    "-toyplot-anchor-shift":"0",
                    },
                )   
    return panel
    
    
    



