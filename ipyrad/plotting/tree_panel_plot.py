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


xxx = {
    "font-size": "14px",
    "text-anchor" : "start",
    "-toyplot-anchor-shift":"0",
    }


## todo: maybe can flip stuff if orient is different
class Panel(object):
    def __init__(self, tree, edges, verts, names):
        ## starting tree position will be changed if adding panel plots
        self.tree = tree
        self.edges = edges
        self.verts = verts
        self.names = names

        ## defaults starting
        self.xmin_tree = 0.
        self.xmin_tips = 0.        
        self.ymin_tree = 0.
        self.ymin_tips = 0.        
        self.ztot = 0.

        ## global kwargs
        self.kwargs = {
            "vsize": 0, 
            "vlshow": False, 
            "ewidth": 3, 
            "vlstyle": {"font-size": "14px"}, 
            "tip_label_style": {"font-size": "12px", 
                                "text-anchor" : "start"},
            "tip_label_color": "#292724",
            "tip_label_offset": None,
            "pct_tree_y": 1.0, 
            "pct_tree_x": 0.5, 
            "lwd_lines": 1,
            }   


    def _panel_tip_labels(self, axes):
        ## if panel.orient...
        orient = 1  ## 'right'

        ## get x,y coords of tip labels ------------------------
        names = [i for i in self.names.values() if not isinstance(i, int)]
        spy = [self.tree.search_nodes(name=i)[0] for i in names]
        spy = np.array([self.verts[i.idx, orient] for i in spy])
        spx = [self.xmin_tips] * len(names)

        ## add to plot. Semantics assume orientation='right'
        _ = axes.text(spx, spy, names,
                    angle=0, 
                    color=self.kwargs["tip_label_color"],
                    style=self.kwargs["tip_label_style"],
                    )


    def _panel_tree(self, axes):
        ## add the tree/graph ------------------------------------------------
        _ = axes.graph(self.edges, 
                       vcoordinates=self.verts, 

                       ## edges 
                       ewidth=self.kwargs["ewidth"], 
                       ecolor=toyplot.color.near_black, 
                       #estyle=...
                       
                       ## vert labels 
                       vlabel=self.names.keys(),
                       vlshow=self.kwargs["vlshow"],
                       vlstyle=self.kwargs["vlstyle"],

                       ## vert icons
                       vsize=self.kwargs["vsize"],
                       ## ...
                       )




## the main function.
def tree_panel_plot(
    tree, 
    edges, 
    verts, 
    names, 
    show_tip_labels,
    use_edge_lengths,
    orient,
    print_args=False,
    *args, 
    **kwargs):
    """
    signature...
    """

    ## create Panel plot object and set height & width
    panel = Panel(tree, edges, verts, names)
    if not kwargs.get("width"):
        panel.kwargs["width"] = min(1000, 50*panel.verts[:, 0].max())
    if not kwargs.get("height"):
        panel.kwargs["height"] = min(1000, 25*(panel.verts[:, 1].max()-1))
    ## update defaults with kwargs & update size based on ntips & ntests
    panel.kwargs.update(kwargs)

    ## partition panels by pct_x pct_y
    if show_tip_labels:
        ## the largest tree-x value is what percent of what we want it to be?
        tree_max_x_should_be = panel.kwargs["width"] * panel.kwargs["pct_tree_x"]
        tree_max_x_is = panel.verts[:, 0].max()
        tree_max_x_muliply_by = float(tree_max_x_should_be) / tree_max_x_is
        panel.verts[:, 0] *= tree_max_x_muliply_by
        panel.xmax_tree = panel.verts[:, 0].max()
        if panel.kwargs["tip_label_offset"]:
            panel.xmin_tips = panel.xmax_tree + panel.kwargs["tip_label_offset"]
        else:
            panel.xmin_tips = panel.xmax_tree + (0.05 * panel.xmax_tree)
        #print(panel.xmax_tree, panel.xmin_tips, panel.verts)

    ## debugger / see all options
    if print_args:
        print(panel.kwargs)

    ## maybe add panels for plotting tip traits in the future
    ## ...

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=panel.kwargs['height'], width=panel.kwargs['width'])
    axes = canvas.cartesian(bounds=("10%", "90%", "10%", "90%"))    
    axes.show = False
    
    ## add panel plots to the axis
    panel._panel_tree(axes)
    if show_tip_labels:
        panel._panel_tip_labels(axes)

    return canvas, axes
    


