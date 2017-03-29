#!/usr/bin/env python2

from __future__ import print_function, division


import numpy as np
import ete3 as ete
import copy
import toyplot
from ..plotting.tree_panel_plot import tree_panel_plot


class Tree(object):
    def __init__(self, newick=None, admix=None, orient='right'):

        ## use default newick string if not given
        if newick:
            ## parse and resolve polytomies
            tree = ete.Tree(newick)
            tree.resolve_polytomy()
            self.newick = tree.write()
            self.admix = admix
        else:
            self.newick = "((((a,b),c),d), ((((e,f),g),h) , (((i,j),k),l)));"
        ## parse newick, assigns idx to nodes, returns tre, edges, verts, names
        tree, edges, verts, names = decompose_tree(self.newick, orient=orient)

        ## parse admixture events
        self.admix = admix
        self._check_admix()

        ## store values
        self.tree = tree
        self.edges = edges
        self.verts = verts
        self.names = names.values()  ## in tree plot vlshow order

        ## hidden attr
        self._orient = orient


    def simulate(self, nreps=1000, admix=None, Ns=int(5e5), gen=10):
        sims = _simulate(self, nreps, admix, Ns, gen)
        debug = 0  ## msprime debug
        names = self.tree.get_leaf_names()[::-1]
        Sims = Sim(names, sims, nreps, debug)
        return Sims


    def root(self, outgroup=None, wildcard=None):
        ## set names or wildcard as the outgroup
        if outgroup:
            outs = [i for i in self.tree.get_leaf_names() if i in outgroup]
        elif wildcard:
            outs = [i for i in self.tree.get_leaves() if wildcard in i.name]
        else:
            raise IPyradError(\
                "must enter either a list of outgroup names or a wildcard selector")
        if len(outs) > 1:
            out = self.tree.get_common_ancestor(outs)
        else:
            out = outs[0]
        self.tree.set_outgroup(out)
        ## we split a branch to root it, so let's double each edge so that the 
        ## distance remains the same (i.e., it would be 0.5 but we make it 1).
        self.tree.children[0].dist = 1.
        self.tree.children[1].dist = 1.
        newick = self.tree.write()
        self.__init__(newick=newick, admix=self.admix, orient=self._orient)


    def _check_admix(self):
        ## raise an error if admixture event is not possible in time period
        if self.admix:
            for event in self.admix:
                pass #print(event)

    def draw(
        self, 
        show_tip_labels=True, 
        use_edge_lengths=False, 
        orient="right",
        *args,
        **kwargs):
        """
        plot the tree using toyplot.graph. 

        Parameters:
        -----------
            taxdicts: dict
                Show tests as colored rectangles.
            bootsarr: ndarray
                Show bootstrap distributions (requires taxdicts as well)
            yaxis: bool
                Show the y-axis.
            use_edge_lengths: bool
                Use edge lengths from newick tree.
            show_tips: bool
                Show tip names from tree.
            pct_tree_y: float
                proportion of canvas y-axis showing tree
            ...
        """
        ## re-decompose tree for new orient and edges args
        tree, edges, verts, names = decompose_tree(
            self.newick, 
            orient=orient, 
            use_edge_lengths=use_edge_lengths)
        verts = verts.astype(np.float)

        ## pass to panel plotter
        canvas, axes = tree_panel_plot(tree, edges, verts, names, 
                                       show_tip_labels, 
                                       use_edge_lengths, 
                                       orient, 
                                       *args, 
                                       **kwargs)
        return canvas, axes


## DEPRECATED
def _draw(self, 
    show_tips, 
    use_edge_lengths, 
    orient, 
    yaxis, 
    *args, 
    **kwargs):

    ## update kwargs from defaults
    args = {"height": min(1000, 15*len(self.tree)),
            "width": min(1000, 15*len(self.tree)),
            "vsize": 0, 
            "vlshow": False, 
            "ewidth": 3, 
            "vlstyle": {"font-size": "18px"}, 
            "cex": "14px", 
            "pct_tree_y": 0.3, 
            "pct_tree_x": 0.7, 
            "lwd_lines": 1,
            }

    ## starting tree position will be changed if adding panel plots
    xmin_tree = 0.
    xmax_tree = 0.
    ymin_tree = 0.
    xmin_tips = 0.
    ymin_tips = 0.

    ## if orient is up or down then swap min/max x/y
    ## ...

    ## re-decompoase the tree in case orientation changed
    tree, edges, verts, names = decompose_tree(
        self.newick, 
        orient=orient, 
        use_edge_lengths=False)
    verts = verts.astype(np.float)

    ## space y-axis to fill canvas height



    ## add spacer for tip names (root of tree is at zero for right and down trees)
    pctt = 0.2 * (xmin_tree + xmin_tips)
    xmin_tree += pctt / 2.
    xmin_tips += pctt / 2.
    verts[:, 1] += pctt / 2.

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=args['height'], width=args['width'])
    axes = canvas.cartesian(bounds=("5%", "95%", "5%", "95%"))
            
    ## add the tree/graph ------------------------------------------------
    _ = axes.graph(edges, 
                    vcoordinates=verts, 
                    ewidth=args["ewidth"], 
                    ecolor=toyplot.color.near_black, 
                    vlshow=args["vlshow"],
                    vsize=args["vsize"],
                    vlstyle=args["vlstyle"],
                    vlabel=names.values())   

    ## add names to tips --------------------------------------------------
    if show_tips:
        ## calculate coords
        nams = [i for i in names.values() if not isinstance(i, int)]
        spx = [tree.search_nodes(name=i)[0] for i in nams]
        spx = np.array([i.x for i in spx], dtype=np.float)
        spx += xmin_tree
        if taxdicts:
            spy = [yhs[-1] - ysp / 2.] * len(nams)
        else:
            spy = [ymin_test - 0.5] * len(nams)
        _ = axes.text(spx, spy, nams,
                        angle=-90, 
                        color=toyplot.color.near_black,
                        style={
                            "font-size": args["cex"],
                            "text-anchor" : "start",
                            "-toyplot-anchor-shift":"0",
                            },
                        ) 

    ## plot admix lines ---------------------------------
    if self.admix:
        for event in self.admix:
            ## get event
            source, sink, _, _, _ = event

            ## get nodes from tree
            source = self.tree.search_nodes(name=source)[0]
            sink = self.tree.search_nodes(name=sink)[0]

            ## get coordinates
            fromx = np.max([source.up.x, source.x]) - np.abs(source.up.x - source.x) / 2.
            fromy = source.y + (source.up.y - source.y) / 2.
            tox = np.max([sink.up.x, sink.x]) - np.abs(sink.up.x - sink.x) / 2.
            toy = sink.y + (sink.up.y - sink.y) / 2.
                
            ## if show_tips:
            if show_tips:
                fromy += spacer
                toy += spacer

            ## plot
            mark = axes.plot([fromx, tox], [fromy, toy], 
                            color=toyplot.color.Palette()[1], 
                            style={"stroke-width": 3, 
                                   "stroke-dasharray": "2, 2"},
                            )
                
    ## hide x and hide/show y axies
    axes.x.show = False
    if yaxis:
        axes.y.show = True
    else:
        axes.y.show = False    

    ## return plotting 
    return canvas, axes


################################################################################
## TREE FUNCTIONS ##############################################################
################################################################################

def _collapse_outgroup(tree, taxdicts):
    """ collapse outgroup in ete Tree for easier viewing """
    ## check that all tests have the same outgroup
    outg = taxdicts[0]["p4"]
    if not all([i["p4"] == outg for i in taxdicts]):
        raise Exception("no good")
   
    ## prune tree, keep only one sample from outgroup
    tre = ete.Tree(tree.write(format=1)) #tree.copy(method="deepcopy")
    alltax = [i for i in tre.get_leaf_names() if i not in outg]
    alltax += [outg[0]]
    tre.prune(alltax)
    tre.search_nodes(name=outg[0])[0].name = "outgroup"
    tre.ladderize()

    ## remove other ougroups from taxdicts
    taxd = copy.deepcopy(taxdicts)
    newtaxdicts = []
    for test in taxd:
        #test["p4"] = [outg[0]]
        test["p4"] = ["outgroup"]
        newtaxdicts.append(test)

    return tre, newtaxdicts


## converts newick to (edges, vertices)
def decompose_tree(newick, use_edge_lengths=True, orient=0):
    """ decomposes tree into component parts for plotting """

    ## invert and short name to arg so it is similar to ape
    ig = use_edge_lengths == False

    ## get tree
    tre = ete.Tree(newick=newick)
    tre.ladderize()

    ## map numeric values to internal nodes from root to tips
    names = {}
    idx = 0
    for node in tre.traverse("preorder"):
        if not node.is_leaf():
            if node.name:
                names[idx] = node.name
            else:
                names[idx] = idx
                node.name = str(idx)
            node.idx = idx
            idx += 1
            
    ## map number to the tips, these will be the highest numbers
    for node in tre.get_leaves(): 
        names[idx] = node.name
        node.idx = idx
        #node.name = str(idx)
        idx += 1

    ## create empty edges and coords arrays
    edges = np.zeros((idx-1, 2), dtype=int)
    verts = np.zeros((idx, 2), dtype=float)

    ## postorder: first children and then parents. This moves up the list .
    nidx = 0
    tip_num = len(tre.get_leaves()) - 1
    for node in tre.traverse("postorder"):
        if node.is_leaf():
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = tip_num 
            node.y = 0
            tip_num -= 1
            verts[node.idx] = [node.x, node.y]                      ##
        
        elif node.is_root():
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y]                      ##
        
        else:
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y]                      ##
        nidx += 1

    ## invert for sideways trees (needs to update node.x and node.y)
    if orient in ['down', 0]:
        verts[:, 0] = verts[:, 0] * -1
    if orient in ['left', 1]:
        verts = verts[:, [1, 0]]
    if orient in ['up', 2]:
        verts[:, 0] = verts[:, 0] * -1
        verts[:, 1] = verts[:, 1] * -1
    if orient in ['right', 3]:
        verts = verts[:, [1, 0]]
        verts[:, 0] = verts[:, 0] * -1

    ## if inverted then make left-most tree be at zero
    verts[:, 0] += abs(verts[:, 0].min())

    return tre, edges, verts, names


## DEPRECATED
def cladogram(newick, use_edge_lengths=True, invert=False):
    
    ## invert and short name to arg so it is similar to ape
    ig = use_edge_lengths == False

    ## get tree
    tre = ete.Tree(newick=newick)
    tre.ladderize()

    ## map numeric values to internal nodes from root to tips
    ## preorder: first parent and then children. These indices will
    ## be used in the int edge array to map edges.
    names = {}
    idx = 0
    for node in tre.traverse("preorder"):
        if not node.is_leaf():
            if node.name:
                names[idx] = node.name
            else:
                names[idx] = idx
                node.name = str(idx)
            node.idx = idx
            idx += 1

    ## map number to the tips, these will be the highest numbers
    for node in sorted(tre.get_leaves(), key=lambda x: x.name):
        names[idx] = node.name
        node.idx = idx
        idx += 1

    ## create empty edges and coords arrays
    edges = np.zeros((idx-1, 2), dtype=int)
    verts = np.zeros((idx, 2), dtype=float)
    
    ## postorder: first children and then parents. This moves up the list .
    nidx = 0
    tip_num = len(tre.get_leaves()) - 1
    for node in tre.traverse("postorder"):
        #for nidx in range(idx)[::-1]:
        #node = tre.search_nodes(idx=nidx)[0]
        if node.is_leaf():
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = tip_num 
            node.y = 0
            tip_num -= 1
            verts[node.idx] = [node.x, node.y]
        
        elif node.is_root():
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y]
        
        else:
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y] 
        nidx += 1

    ## invert for sideways trees
    if invert:
        verts[:, 1] = np.abs(verts[:, 1] - tlen)

    return tre, edges, verts, names


