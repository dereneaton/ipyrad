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


## color palette
COLORS = {"p1": toyplot.color.Palette()[0], 
          "p2": toyplot.color.Palette()[1],
          "p3": toyplot.color.near_black, 
          "p4": toyplot.color.Palette()[-1],}


class Panel(object):
    def __init__(self, ttree, tests, boots, alpha):
        #tree, edges, verts, names, tests, boots, alpha):

        ## starting tree position will be changed if adding panel plots
        #self.tree = tree
        #self.edges = edges
        #self.verts = verts
        #self.names = names
        #self.tests = tests
        ## get attributes from tree object
        for attr in ['tree', 'edges', 'verts', 'names', '_coords', '_lines', '_orient']:
            self.__dict__[attr] = ttree.__dict__[attr]

        self.tests = tests
        self.boots = boots
        try:
            self.allzs = np.abs([i.mean() / i.std() for i in self.boots])
        except Exception:
            self.allzs = []
        self.alpha = alpha

        self.xmin_tree = 0.
        self.xmin_results = 0.
        self.xmin_zscores = 0.        
        self.xmin_hist = 0.

        self.xmax_tree = 0.
        self.xmax_results = 0.
        self.xmax_zscores = 0.
        self.xmax_hist = 0.

        self.ymin_tree = 0.
        self.ymin_test = 0.

        self.ymax_tree = 0.
        self.ymax_test = 0.

        ## default kwargs
        self.kwargs = {
            "vsize": 0, 
            "vlshow": False, 
            "vstyle": {"stroke": "#292724"},            
            "ewidth": 3,             
            "pct_tree_y": 0.3, 
            "pct_tree_x": 0.7, 
            "lwd_lines": 1,
            "vlstyle": {"font-size": "14px"}, 
            #"style_tip_labels": {"font-size": "12px"},
            "style_tip_labels": {"font-size": "12px",
                                 "text-anchor":"start", 
                                 "-toyplot-anchor-shift":"10px", 
                                 "fill": "#292724"},
            "style_test_labels": {"font-size": "12px"},
            "style_results_labels": {"font-size": "12px"},            
            "debug": False,
            "tree_style": "c",            
            "show_tip_labels": True,
        }
        self.check_test_names()


    def check_test_names(self):
        ## check taxdict names
        badnames_set = set()
        for test in self.tests:
            dictnames = list(itertools.chain(*test.values()))
            badnames = [i for i in dictnames if i not in self.names.values()]
            for name in badnames:
                badnames_set.add(name)
        if badnames_set:
            print("Warning: found names not in tree:\n -{}"\
                  .format("\n -".join(list(badnames_set))))        

    @property
    def xpos(self):
        return self.verts[:, 0][self.verts[:, 1] == self.verts[:, 1].min()]
    @property
    def xspacer(self):
        return abs(self.xpos[1] - self.xpos[0])
    @property
    def ypos(self):
        return np.linspace(self.ymin_tree, self.ymin_test, len(self.tests)+2)
    @property
    def yspacer(self):
        return abs(self.ypos[0] - self.ypos[1])  ## the distance between tests
    @property
    def ytrim(self):
        return 0.25 * self.yspacer  ## calculate coords, top&bot 30% of bars cutoff 

    ## plotting functions
    def panel_tree(self, axes):
        _panel_tree(self, axes)

    def panel_test(self, axes):
        _panel_test(self, axes)

    ## these just squeeze in on their own
    def panel_tip_labels(self, axes):
        _panel_tip_labels(self, axes)

    ## if results then plot them
    def panel_results(self, axes):
        _panel_zscores(self, axes)




def _panel_zscores(panel, axes):
    ## space for the results
    tree_width = (panel.xmax_tree - panel.xmin_tree) 
    total_width = (tree_width / float(panel.kwargs["pct_tree_x"]))
    panel.xmin_results = panel.xmin_tree - (total_width - tree_width)
    panel.xmax_results = panel.xmin_tree
    panel.xmin_zscores = (panel.xmin_results + panel.xmax_results) / 2.
    panel.xmax_zscores = panel.xmax_results - (2. * panel.xspacer)
    panel.xmin_hists = panel.xmin_results
    panel.xmax_hists = panel.xmin_zscores - (2. * panel.xspacer)

    if panel.kwargs["debug"]:
        axes.rects(
            panel.xmin_zscores,
            panel.xmax_zscores,
            panel.ymin_test,
            panel.ymax_test,
            color=toyplot.color.Palette()[4], 
            opacity=0.5,
            )
        axes.rects(
            panel.xmin_hists,
            panel.xmax_hists,
            panel.ymin_test,
            panel.ymax_test,
            color=toyplot.color.Palette()[5], 
            opacity=0.5,
            )

    for idx, test in enumerate(panel.tests):
        ## add dashed grid line for tests
        gedges = np.array([[0, 1]])
        gverts = np.array([
                    [panel.xmin_zscores, panel.ypos[idx+1]],
                    [panel.xmax_zscores, panel.ypos[idx+1]],
                    ])
        axes.graph(gedges, 
                   vcoordinates=gverts,
                   ewidth=panel.kwargs["lwd_lines"], 
                   ecolor=toyplot.color.Palette()[-1],
                   vsize=0,
                   vlshow=False,
                   estyle={"stroke-dasharray": "5, 5"},
                   )    

    ## add z-scores bars to lines
    rmax = np.max(np.abs(panel.boots))
    rmax = round(min(1.0, max(0.2, rmax + 0.05*rmax)), 1)
    rmin = round(max(-1.0, min(-0.2, -1 * rmax)), 1)
    zmax = max(3., 1 + float(np.math.ceil(panel.allzs.max())))
    for idx, test in enumerate(panel.tests):
        thisz = panel.allzs[idx]
        zprop = thisz / zmax
        zmaxlen = panel.xmax_zscores - panel.xmin_zscores
        zlen = zmaxlen * zprop        
        axes.rects(
            panel.xmax_zscores, 
            panel.xmax_zscores - zlen,
            panel.ypos[idx+1] - panel.ytrim, 
            panel.ypos[idx+1] + panel.ytrim, 
            color=toyplot.color.Palette()[3],
            )
    ## add axes spines & labels
    ## add xline to z-plots
    gverts = np.array([
        [panel.xmin_zscores, panel.ymin_test],   ## left xbar 
        [panel.xmax_zscores, panel.ymin_test],   ## right xbar
        [panel.xmax_zscores, panel.ymin_tree - panel.ytrim],   ## top ybar
        [panel.xmax_zscores, panel.ymin_test + panel.ytrim],   ## bottom ybar
        [panel.xmin_zscores, panel.ymin_test - panel.ytrim],   ## left tick
        [panel.xmax_zscores, panel.ymin_test - panel.ytrim],   ## right tick
        [panel.xmin_hists, panel.ymin_test],     ## left xbar
        [panel.xmax_hists, panel.ymin_test],     ## right xbar
        [panel.xmin_hists, panel.ymin_test - panel.ytrim],   ## left tick
        [panel.xmax_hists, panel.ymin_test - panel.ytrim],   ## right tick
        ])
    gedges = np.array([
        [0, 1], 
        [2, 3], 
        [4, 0],
        [5, 1],
        [6, 7],
        [8, 6],
        [9, 7],
        ])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        vsize=0, 
        vlshow=False,
        ) 
    ## add text to axes
    tipstyle = {"text-anchor":"middle", "-toyplot-anchor-shift":"0"}
    tipstyle.update(panel.kwargs["style_results_labels"])
    axes.text(
        [panel.xmin_zscores, panel.xmax_zscores, panel.xmin_hists, panel.xmax_hists],
        [panel.ymin_test - panel.yspacer] * 4, 
        [str(zmax), "0.0", str(rmin), str(rmax)],
        color=toyplot.color.near_black,
        style=tipstyle
        )
    ## add text labels
    midpoint_z = panel.xmax_zscores - (panel.xmax_zscores - panel.xmin_zscores) / 2.
    midpoint_h = panel.xmax_hists - (panel.xmax_hists - panel.xmin_hists) / 2.
    axes.text(
        [midpoint_z, midpoint_h],
        #[panel.ymin_test - 2 * panel.yspacer] * 2, 
        [panel.ymin_tree + panel.yspacer / 2.] * 2,
        ["Z-scores", "D-statistics"],
        color=toyplot.color.near_black,
        style=tipstyle
        )

    ## add histograms, and hist axes
    addhists = range(panel.boots.shape[0])[::-1]
    for idx in addhists:
        bins = 30
        mags, xpos = np.histogram(panel.boots[idx], 
                                  bins=bins, 
                                  range=(rmin, rmax), 
                                  density=True)
        ## get hist colors
        thisz = panel.allzs[idx]
        if thisz >= panel.alpha:
            if panel.boots[idx].mean() < 0:
                color = toyplot.color.Palette()[0]
            else:
                color = toyplot.color.Palette()[1]
        else:
            color = toyplot.color.Palette()[-1]

        ## get hist xspans, heights in range
        xran = np.linspace(panel.xmin_hists, panel.xmax_hists, bins+1)
        mags = mags / mags.max()
        mags = (mags * panel.yspacer) * 0.8
        yline = panel.ypos[idx+1] - 1.5*panel.ytrim
        heights = np.column_stack([[yline for i in mags], mags])
        axes.fill(
            xran[:-1], 
            heights, 
            baseline="stacked",
            color=["white", color],
            #style={"stroke": [toyplot.color.to_css(color)]},
            )

    ## add dashed midline to histograms
    gverts = np.array([
        [midpoint_h, panel.ymin_test + panel.ytrim], 
        [midpoint_h, panel.ymin_tree - panel.ytrim],         
        ])
    gedges = np.array([[0, 1]])
    axes.graph(
        gedges, 
        vcoordinates=gverts, 
        ewidth=1, 
        ecolor=toyplot.color.near_black,
        vsize=0, 
        vlshow=False,
        estyle={"stroke-dasharray": "5, 5"},
        )     




def _panel_test(panel, axes):
    ## space for tests (pct_tree_y will scale tests space to tree space)
    panel.xmin_test = panel.xmin_tree
    panel.xmax_test = panel.xmax_tree
    panel.ymax_test = panel.ymin_tree
    tree_height = (panel.ymax_tree - panel.ymin_tree) 
    total_height = (tree_height / float(panel.kwargs["pct_tree_y"]))
    panel.ymin_test = panel.ymin_tree - (total_height - tree_height)

    ## debug the panel
    if panel.kwargs["debug"]:
        print("---rect---")
        print('x_tree', panel.xmin_tree, panel.xmax_tree)
        print('y_tree', panel.ymin_tree, panel.ymax_tree)
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

    ## add features
    add_test_lines(panel, axes)
    add_test_numbers(panel, axes)
    add_rectangles(panel, axes)



def add_rectangles(panel, axes):
    for idx, test in enumerate(panel.tests):    
        ## add test rectangles
        for tax in ["p1", "p2", "p3", "p4"]:
            spx = [panel.tree.search_nodes(name=i)[0] for i in test[tax]]
            spx = np.array([panel.verts[i.idx, 0] for i in spx])
            spx.sort()
        
            ## fill rectangles while leaving holes for missing taxa
            for i in xrange(spx.shape[0]):
                ## if first
                if i == 0:
                    xleft = spx[i] - 0.25
                    xright = spx[i] + 0.25
                ## if not last
                if i != spx.shape[0]-1:
                    if spx[i+1] - spx[i] < 1.5:
                        xright += 1
                    else:
                        axes.rects(xleft, xright, 
                            panel.ypos[idx+1] - panel.ytrim, 
                            panel.ypos[idx+1] + panel.ytrim, 
                            color=COLORS[tax]) 
                        xleft = spx[i+1] - 0.25
                        xright = spx[i+1] + 0.25
                else:
                    axes.rects(xleft, xright,
                            panel.ypos[idx+1] - panel.ytrim, 
                            panel.ypos[idx+1] + panel.ytrim, 
                            color=COLORS[tax]) 


def add_test_lines(panel, axes):
    ## iterate over tests putting in rects
    for idx, test in enumerate(panel.tests):
        ## add dashed grid line for tests
        gedges = np.array([[0, 1]])#, [2, 3]])
        gverts = np.array([
                    [panel.xmin_tree, panel.ypos[idx+1]],
                    [panel.xmax_tree, panel.ypos[idx+1]],
                    ])
                    #[panel.xmin_zscores, panel.ypos[idx+1]], 
                    #[panel.xmax_zscores, panel.ypos[idx+1]], 
                    #])
        axes.graph(gedges, 
                   vcoordinates=gverts,
                   ewidth=panel.kwargs["lwd_lines"], 
                   ecolor=toyplot.color.Palette()[-1],
                   vsize=0,
                   vlshow=False,
                   estyle={"stroke-dasharray": "5, 5"},
                   )                              


def add_test_numbers(panel, axes):
    ## allow for test names
    #if isinstance(show_test_labels, bool) or show_test_labels == 1: 
    #    labels = range(1, len(panel.tests) + 1)
    #elif isinstance(show_test_labels, list):
    #    labels = test_labels
    #else:
    #    raise IPyradError("  label_tests must be a list or boolean")
    labels = range(1, len(panel.tests) + 1)

    ## update styling, text-anchor and anchor-shift are a bit hidden from users
    tipstyle = {"text-anchor":"start", "-toyplot-anchor-shift":"0"}
    tipstyle.update(panel.kwargs["style_test_labels"])

    ## add text
    axes.text(
        [panel.xmax_tree + panel.xspacer] * len(panel.tests), 
        panel.ypos[1:-1], 
        labels,
        color=toyplot.color.near_black,
        style=tipstyle
        )   
    


def _panel_tip_labels(panel, axes):
    ## get coordinates of text
    # names = [i for i in panel.names.values() if not isinstance(i, int)]
    # ypos = [panel.ymin_test] * len(names)

    # ## update styling, text-anchor and anchor-shift are a bit hidden from users
    # tipstyle = {"text-anchor":"start", "-toyplot-anchor-shift":"0"}
    # tipstyle.update(panel.kwargs["style_tip_labels"])

    # ## plot on axes
    # _ = axes.text(panel.xpos, ypos, names,
    #         angle=-90, 
    #         color=toyplot.color.near_black,
    #         style=tipstyle,
    #         ) 

    ## get coordinates of text
    names = panel.tree.get_leaf_names()
    if panel._orient in ["right"]:
        xpos = [panel.verts[:, 0].max()] * len(names)
        ypos = range(len(panel.tree))[::-1] 
        angle = 0.
    elif panel._orient in ['down']:
        xpos = range(len(panel.tree))[::-1]
        ypos = [panel.ymin_test] * len(names) 
        #ypos = [panel.verts[:, 1].min()] * len(names)
        angle = -90.

    tipstyle = {"font-size": "12px",
                "text-anchor":"start", 
                "-toyplot-anchor-shift":"10px", 
                "fill": "#292724"}
    tipstyle.update(panel.kwargs["style_tip_labels"])

    ## plot on axes
    _ = axes.text(xpos, ypos, names,
            angle=angle,
            style=tipstyle,
            ) 


def _panel_tree(panel, axes):
    ## space for the tree
    panel.xmin_tree = panel.verts[:, 0].min()    
    panel.xmax_tree = panel.verts[:, 0].max()
    panel.ymin_tree = panel.verts[:, 1].min()
    panel.ymax_tree = panel.verts[:, 1].max()

    ## debug the panel
    if panel.kwargs["debug"]:
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
    if panel.kwargs["tree_style"] in ["c", "cladogram"]:
        _ = axes.graph(panel.edges, 
                       vcoordinates=panel.verts, 
                       ewidth=panel.kwargs["ewidth"], 
                       ecolor=toyplot.color.near_black, 
                       vlshow=panel.kwargs["vlshow"],
                       vsize=panel.kwargs["vsize"],
                       vlstyle=panel.kwargs["vlstyle"],
                       vlabel=panel.names.keys(),
                       )
    else:
        _ = axes.graph(panel._lines, 
                       vcoordinates=panel._coords, 
                       ewidth=panel.kwargs["ewidth"], 
                       ecolor=toyplot.color.near_black, 
                       #vlshow=panel.kwargs["vlshow"],
                       vsize=0.,
                       vlshow=False,
                       #vlstyle=panel.kwargs["vlstyle"],
                       #vlabel=panel.names.keys(),
                       )
        _ = axes.graph(panel.edges, 
                       vcoordinates=panel.verts, 
                       ewidth=0.,
                       ecolor=toyplot.color.near_black, 
                       vlshow=panel.kwargs["vlshow"],
                       vsize=panel.kwargs["vsize"],
                       vlstyle=panel.kwargs["vlstyle"],
                       vlabel=panel.names.keys(),
                       )




## the main function.
def baba_panel_plot(
    ttree,
    #tree,
    #edges, 
    #verts, 
    #names,
    tests, 
    boots,
    show_tip_labels=True, 
    show_test_labels=True, 
    use_edge_lengths=False, 
    collapse_outgroup=False, 
    pct_tree_x=0.4, 
    pct_tree_y=0.2,
    alpha=3.0,
    *args, 
    **kwargs):
    """
    signature...
    """

    ## create Panel plot object and set height & width
    bootsarr = np.array(boots)
    #panel = Panel(tree, edges, verts, names, tests, bootsarr, alpha)
    panel = Panel(ttree, tests, bootsarr, alpha)
    if not kwargs.get("width"):
        panel.kwargs["width"] = min(1000, 50*len(panel.tree))
    if not kwargs.get("height"):
        panel.kwargs["height"] = min(1000, 50*len(panel.tests))   

    ## update defaults with kwargs & update size based on ntips & ntests
    kwargs.update(dict(pct_tree_x=pct_tree_x, pct_tree_y=pct_tree_y))
    panel.kwargs.update(kwargs)

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=panel.kwargs['height'], width=panel.kwargs['width'])
    axes = canvas.cartesian(bounds=("10%", "90%", "5%", "95%"))    
    axes.show = False

    ## add panels to axes
    panel.panel_tree(axes)
    panel.panel_test(axes)
    panel.panel_tip_labels(axes)
    if isinstance(boots, np.ndarray):
        panel.panel_results(axes)
    return canvas, axes, panel
    


    
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
            panel.ypos[idx] - panel.ytrim, 
            panel.ypos[idx+1] + panel.ytrim,
            color=toyplot.color.Palette()[2])


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
    


    
    



