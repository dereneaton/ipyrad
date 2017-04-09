#!/usr/bin/env python2

from __future__ import print_function, division


import numpy as np
import ete3 as ete
import copy
import toyplot
from ..plotting.tree_panel_plot import tree_panel_plot
from ..assemble.util import IPyradError


class Tree(object):
    def __init__(self, newick=None, admix=None, **kwargs):
        #, orient='right', use_edge_lengths=True):

        ## use default newick string if not given
        if newick:
            ## check that tree can be parsed
            self.tree = ete.Tree(newick)
            self.tree.ladderize()
            self.newick = self.tree.write()
            self.admix = admix
        else:
            self.newick = "((((a,b),c),d), ((((e,f),g),h) , (((i,j),k),l)));"
            self.tree = ete.Tree(self.newick)
            self.tree.ladderize()

        ## parse newick, assigns idx to nodes, returns tre, edges, verts, names
        self._decompose_tree(**kwargs)

        ## parse admixture events
        self.admix = admix
        self._check_admix()


    def __str__(self):
        return self.tree.__str__()


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

        ## we split a branch to root it, so let's double each edge so that the 
        ## distance remains the same (i.e., it would be 0.5 but we make it 1).
        self.tree.set_outgroup(out)
        self.tree.children[0].dist = 2. * self.tree.children[0].dist
        self.tree.children[1].dist = 2. * self.tree.children[1].dist 
        self.tree.support = 100
        newick = self.tree.write()
        self.__init__(newick=newick, 
                      admix=self.admix, 
                      #orient=self._orient,
                      #use_edge_lengths=self._use_edge_lengths
                      )


    def _check_admix(self):
        ## raise an error if admixture event is not possible in time period
        if self.admix:
            for event in self.admix:
                pass #print(event)


    def _decompose_tree(self, orient='right', use_edge_lengths=True):
        _decompose_tree(self, orient, use_edge_lengths)


    ## this is the user-interface where all options should be visible 
    ## and documented
    def draw(
        self, 
        show_tip_labels=True, 
        show_node_support=False,
        use_edge_lengths=False, 
        orient="right",
        print_args=False,
        *args,
        **kwargs):
        """
        plot the tree using toyplot.graph. 

        Parameters:
        -----------
            show_tip_labels: bool
                Show tip names from tree.
            use_edge_lengths: bool
                Use edge lengths from newick tree.
            show_node_support: bool
                Show support values at nodes using a set of default 
                options. 

            ...
        """
        ## re-decompose tree for new orient and edges args
        self._decompose_tree(orient=orient, use_edge_lengths=use_edge_lengths)

        ## update kwargs with entered args and all other kwargs
        dwargs = {}
        dwargs["show_tip_labels"] = show_tip_labels
        dwargs["show_node_support"] = show_node_support
        dwargs.update(kwargs)

        ## pass to panel plotter
        canvas, axes, panel = tree_panel_plot(self, print_args, **dwargs)
        return canvas, axes, panel


## DEPRECATED
def __admix():

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




def _decompose_tree(ttree, orient='right', use_edge_lengths=True): 
    """ decomposes tree into component parts for plotting """

    ## set attributes
    ttree._orient = orient
    ttree._use_edge_lengths = use_edge_lengths
    ult = use_edge_lengths == False

    ## map numeric values to internal nodes from root to tips
    names = {}
    idx = 0
    for node in ttree.tree.traverse("preorder"):
        if not node.is_leaf():
            if node.name:
                names[idx] = node.name
            else:
                names[idx] = idx
                node.name = str(idx)
            node.idx = idx
            idx += 1
            
    ## map number to the tips, these will be the highest numbers
    for node in ttree.tree.get_leaves(): 
        names[idx] = node.name
        node.idx = idx
        idx += 1

    ## create empty edges and coords arrays
    ttree.node_labels = names
    ttree.tip_labels = ttree.tree.get_leaf_names()
    #self.tip_labels = self.tree.get_leaf_names()[::-1]
    #self.node_labels = self.names
    ttree.edges = np.zeros((idx - 1, 2), dtype=int)
    ttree.verts = np.zeros((idx, 2), dtype=float)
    ttree._lines = []        # np.zeros((ntips-1), dtype=int)
    ttree._coords = []       # np.zeros((idx * 2 - ntips), dtype=float)

    ## postorder: first children and then parents. This moves up the list .
    nidx = 0
    tip_num = len(ttree.tree.get_leaves()) - 1
    
    ## tips to root to fill in the verts and edges
    for node in ttree.tree.traverse("postorder"):
        if node.is_leaf():
            ## set the xy-axis positions of the tips
            node.y = ttree.tree.get_distance(node)
            if ult:
                node.y = 0. 
            node.x = tip_num
            tip_num -= 1
            
            ## edges connect this vert to
            ttree.verts[node.idx] = [node.x, node.y]
            ttree.edges[nidx] = [node.up.idx, node.idx]

        elif node.is_root():
            node.y = ttree.tree.get_distance(node)
            if ult:
                node.y = -1 * node.get_farthest_leaf(True)[1] - 1
            node.x = sum(i.x for i in node.children) / float(len(node.children))
            ttree.verts[node.idx] = [node.x, node.y]
        
        else:
            ## create new nodes left and right
            node.y = ttree.tree.get_distance(node)
            if ult:
                node.y = -1 * node.get_farthest_leaf(True)[1] - 1
            node.x = sum(i.x for i in node.children) / float(len(node.children))
            ttree.edges[nidx, :] = [node.up.idx, node.idx]
            ttree.verts[node.idx] = [node.x, node.y]
        nidx += 1
        
    ## root to tips to fill in the coords and lines
    cidx = 0
    for node in ttree.tree.traverse():
        ## add yourself
        if not node.is_leaf():
            ttree._coords += [[node.x, node.y]]
            pidx = cidx
            cidx += 1
            for child in node.children:
                ## add children
                ttree._coords += [[child.x, node.y], [child.x, child.y]]
                ttree._lines += [[pidx, cidx]]    ## connect yourself to newx
                ttree._lines += [[cidx, cidx+1]]  ## connect newx to child
                cidx += 2
    ttree._coords = np.array(ttree._coords, dtype=float)
    ttree._lines = np.array(ttree._lines, dtype=int)

    ## invert for sideways trees
    if ttree._orient in ['up', 0]:
        pass
    if ttree._orient in ['left', 1]:
        ttree.verts[:, 1] = ttree.verts[:, 1] * -1
        ttree.verts = ttree.verts[:, [1, 0]]
        ttree._coords[:, 1] = ttree._coords[:, 1] * -1
        ttree._coords = ttree._coords[:, [1, 0]]
    if ttree._orient in ['down', 0]:
        ttree.verts[:, 1] = ttree.verts[:, 1] * -1
        ttree._coords[:, 1] = ttree._coords[:, 1] * -1
    if ttree._orient in ['right', 3]:
        ttree.verts = ttree.verts[:, [1, 0]]
        ttree._coords = ttree._coords[:, [1, 0]]


