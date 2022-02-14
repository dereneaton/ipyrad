#!/usr/bin/env python

"""
Neighbor-joining tree distance based tree inference method.
"""

import numpy as np
import pandas as pd
from numba import njit
try:
    import toytree
except ImportError:
    pass


# TODO: Implement BioNJ algorithm.
#  Gascuel, O (1997) "BIONJ: an improved version of the NJ 
#  algorithm based on a simple model of sequence data."


class NeighborJoining(object):
    """Construct a neighbor-joining tree from a symmetric distance matrix.

    Parameters
    ----------
    matrix: numpy.ndarray or pandas.DataFrame
        A symmetric matrix with float or int values representing distances 
        between samples. See the ipa.distance module for distance matrix 
        calculations.
    steps: int or None
        Number of iterations to run. This allows returning intermediate dist 
        matrices and trees, and is useful for teaching demonstrations.
    """
    def __init__(
        self, 
        matrix,
        steps=None,
        ):
        self.matrix = matrix


    def run(self):
        """Run inference method to infer distance tree.

        """
        # store tip labels, e.g.. ['a', 'b', 'c', ...]
        idxs = np.arange(self.matrix.shape[0])
        
        # require matrix as float array
        mat = np.array(self.matrix).astype(np.float64)
        
        # iteratively reduce matrix dims by joining nodes until shape=3
        nodes = []
        while mat.shape[0] > 3:
            
            # loop
            fidx, gidx, dist_fu, dist_gu, mat, idxs = jiter_matrix(mat, idxs)

            # return f, g, u, fdist, gdist, nidx
            nodes.append([fidx, gidx, dist_fu, dist_gu, max(idxs)])
        
        # final dist: get dist from internal to newnode
        v0 = (1 / 2.) * mat[0, -1]   
        v1 = (1 / (2 * (3 - 2))) * (mat[-1, :].sum() - mat[0, :].sum())
        dist_vw = v0 + v1
        dist_dw = mat[0, -1] - dist_vw
        dist_ew = mat[1, -1] - dist_vw
        nodes.append([idxs[0], idxs[1], dist_dw, dist_ew, dist_vw])
        
        # build tree from nodes info
        tree = build_tree(nodes)

        # relabel tips
        if isinstance(self.matrix, pd.DataFrame):
            tree = tree.set_node_data(
                feature="name", 
                mapping=dict(enumerate(self.matrix.index)),
            )
        
        # return nodes
        self.tree = tree


@njit
def jcalc_newmat(mat, mask, f, g):
    """
    Get distance of all nodes to the new node and return
    new matrix with dim - 1.
    """
    # get new matrix with dim-1 to fill
    ntips = mat.shape[0] - 1
    newmat = np.zeros((ntips, ntips))
    
    # fill newmat with dist between remaining tips
    newmat[:ntips - 1, :ntips - 1] = mat[mask, :][:, mask]  
    
    # fill final row/col with dists to newnode
    val = (1 / 2.) * (mat[f][mask] + mat[g][mask] - mat[f, g])
    newmat[-1, :-1] = val
    newmat[:-1, -1] = val 
    return newmat
        
        

@njit
def jcalc_q_matrix(mat):
    """
    Get the q-matrix distance of each tip to all other nodes.
    Q(i,j) = (n-2)d(i,j) - Sum(d(i,k) - Sum(d(j,k)    
    """
    qmat = np.zeros(mat.shape)
    n = mat.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if i != j:
                v0 = (n - 2) * mat[i, j]
                v1 = mat[i, :].sum() 
                v2 = mat[j, :].sum()
                val = v0 - v1 - v2
                qmat[i, j] = val
                qmat[j, i] = val
    return qmat



@njit
def jcalc_newdist(mat, f, g):
    """
    Get distance from tips (f, g) to new node
    """
    n = mat.shape[0]
    v0 = (1 / 2.) * mat[f, g]
    v1 = 1 / (2. * (n - 2))
    v2 = mat[f, :].sum()
    v3 = mat[g, :].sum()
    val = v0 + v1 * (v2 - v3)
    return val



@njit
def jiter_matrix(mat, idxs):
    """
    Returns one agglomeration
    """
    # the q-matrix compares all-by-all
    qmat = jcalc_q_matrix(mat)

    # join the nodes with lowest q-score to create new node u
    min_idxs = np.where(qmat == qmat.min())
    f = min_idxs[0][0]
    g = min_idxs[1][0]

    # get distance of new node from children
    dist_fu = jcalc_newdist(mat, f, g)
    dist_gu = mat[f, g] - dist_fu

    # get mask to hide rows and columns of agglomerated indices
    mask0 = np.arange(mat.shape[0]) != f
    mask1 = np.arange(mat.shape[0]) != g
    mask = mask0 & mask1

    # create new matrix removing (f, g) and adding (u)
    mat = jcalc_newmat(mat, mask, f, g)

    # get index labels
    fidx = idxs[f]
    gidx = idxs[g]

    # reset labels
    nidxs = np.zeros(mat.shape[0], dtype=np.int64)
    nidxs[:-1] = idxs[mask]
    nidxs[-1] = max(idxs) + 1
    idxs = nidxs

    # return f, g, u, fdist, gdist, nidx
    return fidx, gidx, dist_fu, dist_gu, mat, idxs



def build_tree(nodes):
    """
    Build a newick string and parse into a tree. Nodes are nested
    in order in the nodes list. Each element contains:
    [node_index0, node_index1, dist0, dist1, parent_index]
    except the final element in which the last element is the 
    internal edge length.
    """
    # store internal node indices to newick strings
    relate = {}
    
    # iterate over internal nodes in order
    for n in nodes[:-1]:
        
        # expand internal nodes to newick strings
        if n[0] in relate:
            n[0] = relate[n[0]]
        if n[1] in relate:
            n[1] = relate[n[1]]
            
        # make newick of the new node pair
        newick = "({}:{},{}:{})".format(n[0], n[2], n[1], n[3])
        
        # store internal node 
        relate[n[4]] = newick
    
    # create final unrooted polytomy
    n0, n1, dist0, dist1, iedge = nodes[-1]    
    newick = "({}:{}, {}:{}, {}:{});".format(newick, iedge, n0, dist0, n1, dist1)
    return toytree.tree(newick)



if __name__ == "__main__":

    import pandas as pd
    import ipyrad.analysis as ipa
    import toyplot.browser

    mat = pd.DataFrame({
        "a": [0, 5, 9, 9, 8],
        "b": [5, 0, 10, 10, 9],
        "c": [9, 10, 0, 8, 7],
        "d": [9, 10, 8, 0, 3],
        "e": [8, 9, 7, 3, 0],
        }, index=list("abcde"),
    )

    # get distance matrix from seq data
    # dist = ipa.distance(mat)
    # dist.run()

    # use this distance matrix to infer a tree
    tool = NeighborJoining(mat)
    tool.run()

    # render the tree
    canvas, axes, mark = tool.tree.draw()

    # display tree
    print(tool.tree)
    toyplot.browser.show(canvas)

