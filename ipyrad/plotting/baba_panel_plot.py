#!/usr/bin/env python2

"""
a panel plot function for baba results 
"""

import ipyrad.analysis as ipa
import numpy as np
import ete3 as ete
import toyplot
import itertools


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



class Panel(object):
    def __init__(self, newick, tests, bootsarr):
        ## starting tree position will be changed if adding panel plots
        self.tree = ipa.tree(newick)
        self.tree.verts = self.tree.verts.astype(np.float)
        self.tests = tests
        self.boots = np.array(bootsarr)
        
        self.xmin_tree = 0.
        self.xmin_hist = 0.
        
        self.ymin_tree = 0.
        self.ymin_test = 0.
        self.ymin_text = 0.        

        self.ztot = 0.

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



## the main function.
def baba_panel_plot(
    newick,
    taxdicts, 
    bootsarr, 
    show_tip_labels=True, 
    show_test_labels=True, 
    use_edge_lengths=False, 
    collapse_outgroup=False, 
    pct_tree_x=0.3, 
    pct_tree_y=0.7,
    *args, 
    **kwargs):
    """
    signature...
    """

    ## update defaults with kwargs & update size based on ntips & ntests
    dwargs.update(kwargs)
    
    ## create Panel plot object and set height & width
    panel = Panel(newick, taxdicts, bootsarr)
    if not dwargs.get("width"):
        dwargs["width"] = min(1000, 50*len(panel.tree.tree))
    if not dwargs.get("height"):
        dwargs["height"] = min(1000, 50*len(panel.tests))   

    #print dwargs
    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=dwargs['height'], width=dwargs['width'])
    axes = canvas.cartesian(bounds=("5%", "95%", "5%", "95%"))    
    axes.show = False
    panel = _panel_tree(axes, panel, dwargs)
    panel = _panel_rect(axes, panel, dwargs)
    panel = _panel_hist(axes, panel, dwargs)
    panel = _panel_test_labels(axes, panel, show_test_labels)
    panel = _panel_tree_labels(axes, panel, show_tip_labels)

    return canvas, axes
    

def _panel_tree_labels(axes, panel, show_tip_labels):
    if show_tip_labels:
        ## calculate coords
        nams = [i for i in panel.tree.names if not isinstance(i, int)]
        spx = [panel.tree.tree.search_nodes(name=i)[0] for i in nams]
        spx = np.array([i.x for i in spx], dtype=np.float)
        spx += panel.xmin_tree
        if panel.tests:
            spy = [panel.yhs[-1] - panel.ysp / 2.] * len(nams)
        else:
            spy = [panel.ymin_test - 0.5] * len(nams)
        _ = axes.text(spx, spy, nams,
                        angle=-90, 
                        color=toyplot.color.near_black,
                        style={
                            "font-size": dwargs["cex"],
                            "text-anchor" : "start",
                            "-toyplot-anchor-shift":"0",
                            },
                        ) 
    return panel

    
    
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
    

    
def _panel_rect(axes, panel, dwargs):

    ## colors for tips
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
            spx = np.array([i.x for i in spx], dtype=np.float)
            spx += panel.xmin_tree
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
    
    
    

def _panel_tree(axes, panel, dwargs):
    ## panel: adjust y-axis: boots Y is N% of tree Y
    pcy = 1 - dwargs["pct_tree_y"]
    newmn = panel.tree.verts[:, 1].max() * pcy
    newhs = np.linspace(newmn, 
                        panel.tree.verts[:, 1].max(), 
                        len(set(panel.tree.verts[:, 1])))
    panel.ymin_tree += panel.tree.verts[:, 1].max() * pcy
    panel.tree.verts[:, 1] = [newhs[int(i)] for i in panel.tree.verts[:, 1]]

    ## adjust X-axis: boots X is 75% of tree X
    ## how much do I need to add to make the tree be 60%?
    panel.xmin_tree += panel.tree.verts[:, 0].max() * dwargs["pct_tree_x"]
    panel.tree.verts[:, 0] += panel.tree.verts[:, 0].max() * dwargs["pct_tree_x"]

    ## get spacer between panels
    panel.ztot = 0.15 * panel.xmin_tree

    ## add spacer for tip names
    pctt = 0.2 * (panel.ymin_tree + panel.ymin_test)
    panel.ymin_tree += pctt / 2.
    panel.ymin_test += pctt / 2.
    panel.tree.verts[:, 1] += pctt / 2.
            
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


