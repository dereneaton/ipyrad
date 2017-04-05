#!/usr/bin/env python2

"""
a panel plot function for baba results 
"""

import numpy as np
import ete3 as ete
import toyplot
import itertools
## avoid circular import by NOT import tree
#from ipyrad.analysis.tree import Tree as tree


class Panel(object):
    def __init__(self, ttree):#tree, edges, verts, coords, lines, names):
        ## get attributes from tree object
        for attr in ['tree', 'edges', 'verts', 'names', '_coords', '_lines', '_orient']:
            self.__dict__[attr] = ttree.__dict__[attr]

        ## defaults panel edges
        self.xmin_tree = 0.
        self.xmax_tree = 0.        
        self.ymin_tree = 0.
        self.ymax_tree = 0.        

        ## default kwargs
        self.kwargs = {
            "vsize": 0, 
            "vmarker": "o",
            "vlshow": False, 
            "ewidth": 3, 
            "vlstyle": {"font-size": "14px"}, 
            "vstyle": {"stroke": "#292724"},
            "style_tip_labels": {"font-size": "12px",
                                 "text-anchor":"start", 
                                 "-toyplot-anchor-shift":"10px", 
                                 "fill": "#292724"},
            "show_tip_labels": True,
            "tree_style": "p",
            "show_axes": False,
            "debug": False,
            }


    def _panel_tip_labels(panel, axes):
        ## get coordinates of text
        names = panel.tree.get_leaf_names()
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

        ## plot on axes
        _ = axes.text(xpos, ypos, names,
                angle=angle,
                style=tipstyle,
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
                           vlabel=self.names.keys(),
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
            nodestyle = {'fill': toyplot.color.to_css(toyplot.color.Palette()[0]), 
                         'stroke': toyplot.color.to_css(toyplot.color.Palette()[0])}
            nodestyle.update(self.kwargs["vstyle"])
            _ = axes.graph(self.edges, 
                           vcoordinates=self.verts, 
                           ewidth=0.,
                           vmarker=self.kwargs["vmarker"],
                           vlabel=self.kwargs["vlabels"],
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
        supps = [int(panel.tree.search_nodes(idx=j)[0].support) for j in range(nnodes)]
        if not panel.kwargs["vsize"]:
            panel.kwargs["vsize"] = 20
        sizes = [panel.kwargs["vsize"] for j in range(nnodes)]
        ## add leaf values
        supps += [""] * len(panel.tree)
        sizes += [0] * len(panel.tree)
        ## override args
        panel.kwargs["vlabels"] = supps
        panel.kwargs["vsize"] = sizes
        panel.kwargs["vlshow"] = True
        #panel.kwargs["vmarker"] = 's'  ## square
        ## if unrooted then hide root node scores
        if len(panel.tree.children) > 2:
            supps[0] = ""
            sizes[0] = 0
        #print(panel.kwargs["vlabels"])
        #print(panel.kwargs["vsize"])
    else:
        panel.kwargs["vlabels"] = panel.names.keys()

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
    


