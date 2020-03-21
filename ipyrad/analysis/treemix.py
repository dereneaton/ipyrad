#!/usr/bin/env python

# py2/3 compat
from __future__ import print_function

# standar lib
import os
import sys
import gzip
import subprocess

# third party
import numpy as np
import pandas as pd
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import Params
from .snps_extracter import SNPsExtracter
# from ipyrad.assemble.write_outfiles import reftrick, GETCONS


# missing imports to be raised on class init
try:
    import toytree
except ImportError:
    pass
# missing imports to be raised on class init
try:
    import toyplot
except ImportError:
    pass

_MISSING_TOYTREE = ImportError("""
This ipyrad tool requires the plotting library toytree. 
You can install it with the following command in a terminal.

conda install toytree -c eaton-lab 
""")

_MISSING_TREEMIX = ImportError("""
This ipyrad tool requires the progam TREEMIX. See recommended installation 
instructions here:

https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-treemix.html
""")



class Treemix(object):
    """
    Treemix analysis utility function.

    Parameters:
    -----------
    data: str or tuple
        There are *two* options for your input file. 
        1. '.treemix.in.gz': This is a file that is already formatted for
            running in treemix. If you use this no futher formatting will 
            be performed on the data (e.g., imap and minmap will be ignored).
        2. '.snps.hdf5': This file is produced by ipyrad and contains the 
            snps information as well as information about their location and 
            linkage to be used to subsample unlinked SNPs. This input is 
            preferable since we can easily re-sample unlinked SNPs over 
            multiple bootstrap replicates, and can easily include or exclude
            samples using imap or minmap options. 

    name: str
        The name for this run. An alias for '-n'.

    workdir: str
        The output directory for results. An alias for '-w'. 

    imap: dict
        A dictionary mapping population names to sample names. This defines
        how individuals will be grouped to calculate allele frequencies. Only
        samples present in the imap will be inluded in the analysis.

    minmap: dict
        A dictionary mapping population names to minimum sample numbers. You
        can filter SNPs from your data set to include only those that are 
        sampled across at least N individuals in each pop. 

    Attributes:
    -----------
    params: dict
        parameters for this treemix run
    command: str
        returns the command string to run treemix 

    Functions:
    ----------
    write_treemix_file()
        writes a gzipped input file to be used in the treemix program.
    run()
        submits a treemix job to run locally or on a ipyparallel cluster. 
    """    

    # init object for params
    def __init__(
        self,
        data,
        name="test",
        workdir="analysis-treemix", 
        imap=None,
        minmap=None,
        seed=None,
        quiet=False,
        raise_root_error=False,
        binary=None,
        *args, 
        **kwargs):

        # path attributes
        self.name = name
        self.data = data

        # if not imap then it will be set to 1
        self.minmap = minmap
        self.imap = imap

        # others
        self.binary = os.path.join(sys.prefix, "bin", "treemix")
        self.binary = (binary if binary else self.binary)
        self.raise_root_error = raise_root_error
        self._find_binary()

        # params dict
        self.params = Params()
        self.params.k = 0
        self.params.m = 0
        self.params.g = (None, None)
        self.params.bootstrap = 0
        self.params.cormig = 0
        self.params.climb = 0
        self.params.noss = 0
        self.params.seed = (seed if seed else np.random.randint(0, int(1e9)))
        self.params.root = None
        self.params.se = 0
        self.params.global_ = 0

        # get snps and snpmap
        ext = SNPsExtracter(self.data, self.imap, self.minmap, quiet=quiet)
        ext.parse_genos_from_hdf5()
        self.snps = ext.subsample_snps(seed)
        self.names = ext.names
        self.nsites = self.snps.shape[1]
        self.pops = self.imap.keys()
        self.sidxs = {
            pop: [self.names.index(i) for i in self.imap[pop]]
            for pop in self.pops
        }

        # make workdir if it does not exist
        if workdir:
            self.workdir = os.path.abspath(os.path.expanduser(workdir))
        else:
            self.workdir = os.path.join(
                os.path.abspath(os.path.curdir),
                "analysis-treemix",
            )
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        ## set params
        notparams = set(
            ["workdir", "name", "data", "minmap", "imap", "seed", "quiet"]
        )
        for key in set(kwargs.keys()) - notparams:
            self.params[key] = kwargs[key]

        # results files
        self.files = Params()
        self.files.tree = os.path.join(self.workdir, self.name + ".treeout.gz")
        self.files.cov = os.path.join(self.workdir, self.name + ".cov.gz")
        self.files.llik = os.path.join(self.workdir, self.name + ".llik")

        # results
        self.results = Params()
        self.results.tree = ""
        self.results.admixture = []
        self.results.cov = []
        self.results.llik = None


    @property
    def _command_list(self):
        """ build the command list """

        # base args
        cmd = [
            self.binary, 
            "-i", os.path.join(self.workdir, self.name + ".treemix.in.gz"),
            "-o", os.path.join(self.workdir, self.name),
        ]

        # addon params
        args = []
        for key, val in self.params.__dict__.items():
            if key == "g":
                if val[0]:
                    args += ["-" + key, str(val[0]), str(val[1])]
            elif key == "global_":
                if val:
                    args += ["-" + key[:-1]]
            elif key in ["se", "global", "noss"]: 
                if val:
                    args += ["-" + key]
            else:
                if val:
                    args += ["-" + key, str(val)]
        return cmd + args


    @property
    def command(self):
        """ returns command as a string """
        return " ".join(self._command_list)


    def write_treemix_file(self, quiet=False):
        """
        Write genos to treemix gzipped format:
        A   B   C   D
        0,2 2,0 2,0 0,2
        0,2 1,1 3,0 0,3
        0,2 2,0 3,0 0,2
        ...
        """
        outfile = os.path.join(self.workdir, self.name + ".treemix.in.gz")
        outf = open(outfile, 'w')

        # write the headers        
        popnames = sorted(self.imap)
        outf.write(" ".join(popnames) + "\n")

        # create 0,5 pairs for ancestral derived counts
        poptuples = {}
        for pop in popnames:
            ances = np.sum(self.snps[self.sidxs[pop], :] == 0, axis=0) * 2
            deriv = np.sum(self.snps[self.sidxs[pop], :] == 2, axis=0) * 2
            heter = np.sum(self.snps[self.sidxs[pop], :] == 1, axis=0)
            ances += heter
            deriv += heter

            arr = ["{},{}".format(i, j) for i, j in zip(ances, deriv)]
            poptuples[pop] = arr

        # save to file
        np.savetxt(
            outf, 
            np.vstack([poptuples[i] for i in popnames]).T,
            delimiter=" ", 
            fmt="%s",
        )

        # close file handle
        if not quiet:
            print("wrote treemix input file to {}".format(outfile))
        outf.close()


    def draw_tree(self, axes=None):
        """
        Returns a treemix plot on a toyplot.axes object. 
        """
        # create a toytree object from the treemix tree result
        tre = toytree.tree(newick=self.results.tree)

        # draw on axes or create new ones
        if axes:
            canvas = None
            tre.draw(
                axes=axes,
                use_edge_lengths=True,
                tip_labels_align=True,
                edge_align_style={"stroke-width": 1},
                scalebar=True,
            )
        else:
            canvas, axes = tre.draw(
                axes=axes,
                use_edge_lengths=True,
                tip_labels_align=True,
                edge_align_style={"stroke-width": 1},
                scalebar=True,
            )

        # get coords 
        for admix in self.results.admixture:
            # parse admix event
            pidx, pdist, cidx, cdist, weight = admix
            a = _get_admix_point(tre, pidx, pdist)
            b = _get_admix_point(tre, cidx, cdist)

            axes.graph(
                np.array([[0, 1]]),
                vcoordinates=np.array([[a[0], a[1]], [b[0], b[1]]]),
                tmarker=">",
                ewidth=8 * weight,
                eopacity=0.8,
                vlshow=False,
                vsize=0,
            )   
        return canvas, axes


    def draw_cov(self, axes=None):

        # get results
        cov = self.results.cov
        tre = toytree.tree(self.results.tree)

        # names spaced in order
        lnames = toyplot.locator.Explicit(
            locations=range(len(tre.get_tip_labels())),
            labels=tre.get_tip_labels()[::-1],
        )

        # get a colormap and plot the matrix
        cmap = toyplot.color.diverging.map(
            "BlueRed", 
            cov.min(),
            cov.max(),
        )

        canvas, table = toyplot.matrix(
            (cov, cmap),
            width=400, 
            height=400, 
            bshow=True,
            tshow=False,
            lshow=False,
            rlocator=lnames,
            blocator=lnames,      
        )
        return canvas, table


    def run(self, quiet=True):

        # call command
        self.write_treemix_file(quiet=True)

        # call treemix and catch root errors
        try:
            self._call_treemix()
        except KeyError:
            oldroot = self.params.root
            self.params.root = None
            self._call_treemix()
            self.params.root = oldroot

        # parse treemix results files.
        self._parse_results()

 
    def _call_treemix(self):
        """
        Calls command on subprocess.
        """
        # build command line arg
        proc = subprocess.Popen(
            self._command_list,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
        )

        # run command and capture output
        self.stdout, _ = proc.communicate()
        self.stdout = self.stdout.decode()

        # error was raised, check for root error
        if proc.returncode:

            # trim to error message
            msg = self.stdout.strip().split("\n")[-1]

            # skip root error here... gonna treat it as a keyerror
            if ("ERROR in placing root" in msg) and (not self.raise_root_error):
                raise KeyError("ROOT ERROR")
            else:
                raise IPyradError(msg)


    def _parse_results(self):
        """
        Parse results. 

        files.tree --> results.tree, results.admixture
        files.cov --> results.cov
        files.llik --> results.llik
        """
        # get tree and admix from output files 
        with gzip.open(self.files.tree) as tmp:
            data = tmp.readlines()

            # store the tree
            k = 0
            if self.params.noss:
                k += 1
            self.results.tree = data[k].strip().decode()
            self.results.admixture = []

            # get admix events
            for adx in data[k + 1:]:
                dat = [i.decode() for i in adx.strip().split()]
                weight, jweight, jse, pval, clade1, clade2 = dat
                self.results.admixture.append(
                    (clade1, clade2, weight, jweight, jse, pval)
                )

        # get a toytree
        tre = toytree.tree(self.results.tree)
        names = tre.get_tip_labels()[::-1]

        # order admixture 
        for aidx in range(len(self.results.admixture)):
            admix = self.results.admixture[aidx]

            source = toytree.tree(admix[0] + ";")
            if len(source) == 1:
                name = admix[0].split(":")[0]
                sodx = tre.treenode.search_nodes(name=name)[0]
                sodx = sodx.idx
            else:
                lvs = source.get_tip_labels()
                sodx = tre.treenode.get_common_ancestor(lvs).idx

            sink = toytree.tree(admix[1] + ";")
            if len(sink) == 1:
                name = admix[1].split(":")[0]
                sidx = tre.treenode.search_nodes(name=name)[0]
                sidx = sidx.idx
            else:
                lvs = sink.get_tip_labels()
                sidx = tre.treenode.get_common_ancestor(lvs).idx

            self.results.admixture[aidx] = (
                int(sodx),
                float(admix[0].rsplit(":", 1)[1]), 
                int(sidx),
                float(admix[1].rsplit(":", 1)[1]), 
                float(admix[2]), 
                )

        # parse the cov -------------------------
        dat = pd.read_csv(
            self.files.cov,
            sep=" ", 
            header=None, 
            index_col=0, 
            skiprows=1,
        )
        # sort into tree label order
        dat.columns = dat.index
        dat = dat.loc[names]
        dat = dat.T.loc[names]
        self.results.cov = dat.values
        
        # parse the llik -------------------------
        with open(self.files.llik) as indat:
            self.results.llik = float(indat.readlines()[-1].split()[-1])



    def _find_binary(self):
        # check for binary
        list_binaries = [self.binary]

        # check user binary first, then backups
        for binary in list_binaries:
            proc = subprocess.Popen(["which", binary],
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT).communicate()
            # if a binary was found then stop
            if proc[0]:
                return binary

        # if not binaries found
        raise Exception(_MISSING_TREEMIX)



# plotting functions
def _get_admix_point(tre, idx, dist):
    ## parent coordinates
    px, py = tre._coords.verts[idx]
    ## child coordinates
    cx, cy = tre._coords.verts[tre.treenode.search_nodes(idx=idx)[0].up.idx]
    ## angle of hypotenuse
    theta = np.arctan((px - cx) / (py - cy))
    ## new coords along the hypot angle
    horz = np.sin(theta) * dist
    vert = np.cos(theta) * dist
    
    ## change x
    a = tre._coords.verts[idx, 0]
    b = tre._coords.verts[idx, 1] 
    a -= abs(horz)
    if py < cy:
        b += abs(vert)
    else:
        b -= abs(vert)
    return a, b
