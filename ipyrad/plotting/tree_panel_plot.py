#!/usr/bin/env python2

"""
a panel plot function for baba results 
"""

#import numpy as np
#import ete3 as ete
from __future__ import print_function
import toyplot
#import itertools
## avoid circular import by NOT import tree
#from ipyrad.analysis.tree import Tree as tree


class Panel(object):
    def __init__(self, ttree):
        ## get attributes from tree object

        ## should I just use a super class of ttree?
        self.tree = ttree.tree
        self.edges = ttree.edges
        self.verts = ttree.verts
        self._coords = ttree._coords
        self._lines = ttree._lines
        self._orient = ttree._orient
        self.tip_labels = ttree.tip_labels
        self.node_labels = ttree.node_labels
        #self.names = ttree.names

        ## defaults panel edges
        self.xmin_tree = 0.
        self.xmax_tree = 0.        
        self.ymin_tree = 0.
        self.ymax_tree = 0.        

        ## default kwargs
        self.kwargs = {
            ## edge defaults
            "ewidth": 3, 

            ## node marker defaults
            "vsize": 0, 
            "vmarker": "o",
            "vstyle": {"stroke": "#292724"},

            ## node label defaults
            "vlshow": False, 
            "vlstyle": {"font-size": "14px"}, 

            ## tip label defaults
            "show_tip_labels": True,
            "tip_labels": None,
            "color_tip_labels": None,  
            "style_tip_labels": {"font-size": "12px",
                                 "text-anchor":"start", 
                                 "-toyplot-anchor-shift":"10px", 
                                 "fill": "#292724"},
            ## tree style and axes
            "tree_style": "p",
            "show_axes": False,
            "debug": False,
            }



    def _panel_tip_labels(panel, axes):

        ## get tip labels. Order is from top to bottom on right-facing
        if panel.kwargs["tip_labels"]:
            names = panel.kwargs["tip_labels"]
        elif panel.kwargs["tip_labels"] == False:
            names = ["" for _ in panel.tip_labels] 
        else:
            names = panel.tip_labels

        ## get coordinates of text from top to bottom (right-facing)
        if panel._orient in ["right"]:
            xpos = [panel.verts[:, 0].max()] * len(names)
            ypos = range(len(panel.tree))[::-1]
            angle = 0.
        elif panel._orient in ['down']:
            xpos = range(len(panel.tree))[::-1]
            ypos = [panel.verts[:, 1].min()] * len(names)
            angle = -90.
        tipstyle = {"font-size": "12px",
                    "text-anchor":"start", 
                    "-toyplot-anchor-shift":"10px", 
                    "fill": "#292724"}
        tipstyle.update(panel.kwargs["style_tip_labels"])

        ## tip color overrides tipstyle[fill]
        if panel.kwargs.get("color_tip_labels"):
            tipstyle.pop("fill")

        ## plot on axes. color is added from top to bottom (right-facing)
        _ = axes.text(xpos, ypos, names,
                angle=angle,
                style=tipstyle,
                color=panel.kwargs["color_tip_labels"] #[::-1],
                ) 



    def _panel_tree(self, axes):
        if self._orient in ["right"]:
            self.xmax_tree = self.verts[:, 0].max()

        ## add the tree/graph ------------------------------------------------
        if self.kwargs["tree_style"] in ["c", "cladogram"]:
            _ = axes.graph(self.edges, 
                           vcoordinates=self.verts, 
                           ewidth=self.kwargs["ewidth"], 
                           ecolor=toyplot.color.near_black, 
                           #estyle=... add round edges 
                           vlabel=self.node_labels.keys(),#names.keys(),
                           #vlabel=self.kwargs["vlabels"],                           
                           vlshow=self.kwargs["vlshow"],
                           vlstyle=self.kwargs["vlstyle"],
                           vsize=self.kwargs["vsize"],
                           )
        else:
            ## add lines for phylogram
            _ = axes.graph(self._lines, 
                           vcoordinates=self._coords, 
                           ewidth=self.kwargs["ewidth"], 
                           ecolor=toyplot.color.near_black, 
                           vlshow=False,
                           vsize=0.,
                           )
            ## add vertices for phylogram
            nodestyle = {
                'fill': toyplot.color.to_css(toyplot.color.Palette()[0]), 
                'stroke': toyplot.color.to_css(toyplot.color.Palette()[0])
                }
            nodestyle.update(self.kwargs["vstyle"])
            _ = axes.graph(self.edges, 
                           vcoordinates=self.verts, 
                           ewidth=0.,
                           vmarker=self.kwargs["vmarker"],
                           vlabel=self.kwargs["vlabel"],
                           #vlabel=self.names.keys(),
                           vlshow=self.kwargs["vlshow"],
                           vlstyle=self.kwargs["vlstyle"], 
                           vsize=self.kwargs["vsize"],
                           vstyle=nodestyle,
                           )
        


## the main function.
def tree_panel_plot(ttree,
    print_args=False,
    *args, 
    **kwargs):
    """
    signature...
    """

    ## create Panel plot object and set height & width
    panel = Panel(ttree)          #tree, edges, verts, names)
    if not kwargs.get("width"):
        panel.kwargs["width"] = min(1000, 25*len(panel.tree))
    if not kwargs.get("height"):
        panel.kwargs["height"] = panel.kwargs["width"]

    ## update defaults with kwargs & update size based on ntips & ntests
    panel.kwargs.update(kwargs)

    ## magic node label arguments overrides others
    if panel.kwargs["show_node_support"]:
        nnodes = sum(1 for i in panel.tree.traverse()) - len(panel.tree)
        ## set node values
        supps = [int(panel.tree.search_nodes(idx=j)[0].support) \
                 for j in range(nnodes)]
        if not panel.kwargs["vsize"]:
            panel.kwargs["vsize"] = 20
        sizes = [panel.kwargs["vsize"] for j in range(nnodes)]
        ## add leaf values
        supps += [""] * len(panel.tree)
        sizes += [0] * len(panel.tree)
        ## override args
        panel.kwargs["vlabel"] = supps
        panel.kwargs["vsize"] = sizes
        panel.kwargs["vlshow"] = True
        #panel.kwargs["vmarker"] = 's'  ## square
        ## if unrooted then hide root node scores
        if len(panel.tree.children) > 2:
            supps[0] = ""
            sizes[0] = 0
        #print(panel.kwargs["vlabels"])
        #print(panel.kwargs["vsize"])
    elif panel.kwargs.get("vlabel"):
        panel.kwargs["vlabel"] = panel.kwargs["vlabel"]
        panel.kwargs["vlshow"] = True
    else:
        panel.kwargs["vlabel"] = panel.node_labels.keys() #names.keys()

    ## debugger / see all options
    if print_args:
        print(panel.kwargs)

    ## maybe add panels for plotting tip traits in the future
    ## ...

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=panel.kwargs['height'], width=panel.kwargs['width'])
    axes = canvas.cartesian(bounds=("10%", "90%", "10%", "90%"))    
    axes.show = panel.kwargs["show_axes"]
    
    ## add panel plots to the axis
    panel._panel_tree(axes)
    if panel.kwargs["show_tip_labels"]:
        panel._panel_tip_labels(axes)

    return canvas, axes, panel
    


