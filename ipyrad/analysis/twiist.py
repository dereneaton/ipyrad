#!/usr/bin/env python

"Implementation of a twiist like method for reference or anonymous RAD loci"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard
import os
import glob
import uuid
import tempfile
import itertools
import subprocess as sps

# third party
import pandas as pd
import numpy as np
import toytree

from .raxml import Raxml as raxml
from ipyrad.assemble.utils import IPyradError
import ipyrad as ip


class Twiist():
    """
    Performs phylo inference across loci sampled in windows. If denovo then
    span is the number of randomly sampled loci, and window is the number of 
    random samples. If reference then window is the size of the window selected
    and span is the 
    """
    def __init__(
        self,
        name, 
        data,
        imap, 
        minmap=None,
        reference=False,
        window=None, 
        slide=None,
        minsnps=1,
        randomseed=None,
        ):

        ## init random seed
        self.randomseed = randomseed
        np.random.seed(self.randomseed)
        
        ## store attributes
        self.name = name
        self.data = data
        self.imap = imap
        self.rmap = {}
        self.results_table = None
        self.span = span
        self.window = window
        self.minsnps = minsnps

        # reference data
        self.reference = reference
        self.chromlendict = {}

        # reverse of imap dictionary
        for k, v in self.imap.items():
            for i in v:
                self.rmap[i] = k
        self.ntests = ntests
        self.minsnps = minsnps
        
        if not minmap:
            minmap = {i: 1 for i in self.imap}
        self.minmap = minmap
        
        # store all samples for this test
        self.samples = list(itertools.chain(*[i for i in self.imap.values()]))
        
        ## get idxs of loci for this test
        if self.reference:
            self.check_reference_is_indexed()
            self.idxs = self.get_ref_locus_idxs()

        else:
            self.idxs = self.get_denovo_locus_idxs()
                

    def check_reference_is_indexed(self):
        if not os.path.exists(self.reference):
            raise IOError(
                "Reference file not found: {}".format(self.reference))

        # If reference index exists then bail out unless force
        if not os.path.exists(self.reference + ".fai"):

            # simple samtools index for grabbing ref seqs
            cmd = [ip.bins.samtools, "faidx", self.reference]
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            error = proc.communicate()[0].decode()

            # error handling
            if error:
                if "please use bgzip" in error:
                    raise IPyradError("BGZIP reference file not supported")
                else:
                    raise IPyradError(error)

        # store  chrom lengths
        table = pd.read_csv(
            self.reference + ".fai", 
            names=['scaffold', 'length', 'start', 'a', 'b'],
            sep="\t",
        )
        self.chromlendict = {
            table.scaffold[i]: table.length[i] for i in range(table.shape[0])
        }

    

    def get_ref_locus_idxs(self):
        idxs = []

        with open(self.data) as indata:
            liter = (indata.read().strip().split("|\n"))

        for idx, loc in enumerate(liter):
            lines = loc.split("\n")

            snpline = lines[-1]
            locidx, chidx, pos = snpline.split("|")[1].split(":")            
            names = [i.split()[0] for i in lines[:-1]]

            ## check coverage
            coverage = 0
            for node in self.imap:
                mincov = self.minmap[node]
                if sum([i in names for i in self.imap[node]]) >= mincov:
                    coverage += 1
            
            ## check that the coverage meets threshold
            ## refinfo is not being added correctly....
            if coverage == len(self.imap.keys()):
                pos1, pos2 = pos.split('-')
                refinfo = (idx, chidx, pos1, pos2)
                idxs.append(refinfo) 
                print(refinfo)
        return idxs


    
    def get_denovo_locus_idxs(self):
        """ finds loci with sufficient sampling for this test"""

        ## store idx of passing loci
        idxs = []

        ## open handle
        with open(self.data) as indata:
            liter = (indata.read().strip().split("|\n"))

        ## put chunks into a list
        for idx, loc in enumerate(liter):

            ## parse chunk
            lines = loc.split("\n")[:-1]
            names = [i.split()[0] for i in lines]

            ## check coverage
            coverage = 0
            for node in self.imap:
                mincov = self.minmap[node]
                if sum([i in names for i in self.imap[node]]) >= mincov:
                    coverage += 1
            if coverage == 4:
                idxs.append(idx)      

        ## concatenate into a phylip file
        return idxs
    


    def sample_loci(self, window):
        """ finds loci with sufficient sampling for this test"""

        ## store idx of passing loci
        if not self.reference:
            idxs = np.random.choice(self.idxs, self.ntests)
        else:
            idxs = self.get_window_idxs(window)

        ## open handle, make a proper generator to reduce mem
        with open(self.data) as indata:
            liter = (indata.read().strip().split("|\n"))

        ## store data as dict
        seqdata = {i: "" for i in self.samples}

        ## put chunks into a list
        for idx, loc in enumerate(liter):
            if idx in idxs:
                ## parse chunk
                lines = loc.split("\n")[:-1]
                names = [i.split()[0] for i in lines]
                seqs = [i.split()[1] for i in lines]
                dd = {i: j for (i, j) in zip(names, seqs)}

                ## add data to concatenated seqdict
                for name in seqdata:
                    if name in names:
                        seqdata[name] += dd[name]
                    else:
                        seqdata[name] += "N" * len(seqs[0])
                        
        ## concatenate into a phylip file
        return seqdata
  


    def get_window_idxs(self, window):
        "Returns locus idxs that are on same chrom and within window"

        # get start position
        locidx, chridx, pos1, pos2 = self.idxs[window]

        # get end position
        endwindow = pos2 + self.window

        idxs = []
        for tup in self.idxs:
            if tup[1] == chridx:
                if int(tup[3]) < int(endwindow):
                    idxs.append(tup[0])
        return idxs



    def run_tree_inference(self, nexus, idx):
        """
        Write nexus to tmpfile, runs phyml tree inference, and parses
        and returns the resulting tree. 
        """
        ## create a tmpdir for this test
        tmpdir = tempfile.tempdir
        tmpfile = os.path.join(tempfile.NamedTemporaryFile(
            delete=False,
            prefix=str(idx),
            dir=tmpdir,
        ))

        ## write nexus to tmpfile
        tmpfile.write(str.encode(nexus))
        tmpfile.flush()

        ## infer the tree
        rax = raxml(name=str(idx), data=tmpfile.name, workdir=tmpdir, N=1, T=2)
        rax.run(force=True, block=True, quiet=True)

        ## clean up
        tmpfile.close()

        ## return tree order
        order = get_order(toytree.tree(rax.trees.bestTree))
        return "".join(order)


    
    def run_denovo(self, ipyclient):
        """
        parallelize calls to worker function.
        """
        ## connect to parallel client
        lbview = ipyclient.load_balanced_view()
        
        ## iterate over tests
        asyncs = []
        for window in range(self.ntests): 
            
            ## submit jobs to run
            args = (window, self)
            rasync = lbview.apply(denovo_worker, *args)
            asyncs.append(rasync)
            
        ## wait for jobs to finish
        ipyclient.wait()

        ## check for errors
        for rasync in asyncs:
            if not rasync.successful():
                raise Exception("Error: {}".format(rasync.result()))

        ## return results as df
        results = [i.result() for i in asyncs]
        self.results_table = pd.DataFrame(results)
    
    

    def run_reference(self, ipyclient):
        "parallelize worker function for reference data"

        # connect to client
        lbview = ipyclient.load_balanced_view()
        rasyncs = []
        for chrom, chromlen in self.chromlendict.items():
            # skip chroms that are smaller than window size
            if chromlen > self.window:
                for window in range(0, chromlen, self.window):
                    # check if there is data in this window
                    pass


    def plot(self):
        """
        return a toyplot barplot of the results table.
        """
        if self.results_table is None:
            return "no results found"
        else:
            bb = self.results_table.sort_values(
                by=["ABCD", "ACBD"], 
                ascending=[False, True],
                )

            ## make a barplot
            import toyplot
            c = toyplot.Canvas(width=600, height=200)
            a = c.cartesian()
            m = a.bars(bb)
            return c, a, m



def denovo_worker(self, window):

    # sample random N (window) loci and return as a dict 
    fullseqs = self.sample_loci(window)

    # find all iterations of samples for this quartet
    liters = itertools.product(*self.imap.values())

    ## run tree inference for each iteration of sampledict
    hashval = uuid.uuid4().hex
    weights = []




def worker(self, window):
    """ 
    Calculates the quartet weights for the test at a random
    subsampled chunk of loci.
    """

    # sample random N (window) loci and return as a dict 
    fullseqs = self.sample_loci(window)

    # find all iterations of samples for this quartet
    liters = itertools.product(*self.imap.values())

    ## run tree inference for each iteration of sampledict
    hashval = uuid.uuid4().hex
    weights = []
    for ridx, lidx in enumerate(liters):
        
        ## get subalignment for this iteration and make to nex
        a, b, c, d = lidx
        sub = {}
        for i in lidx:
            if self.rmap[i] == "p1":
                sub["A"] = fullseqs[i]
            elif self.rmap[i] == "p2":
                sub["B"] = fullseqs[i]
            elif self.rmap[i] == "p3":
                sub["C"] = fullseqs[i]
            else:
                sub["D"] = fullseqs[i]
                
        ## write as nexus file
        nex = []
        for tax in list("ABCD"):
            nex.append(">{}         {}".format(tax, sub[tax]))
            
        ## check for too much missing or lack of variants
        nsites, nvar = count_var(nex)

        ## only run test if there's variation present
        if nvar > self.minsnps:
               
            ## format as nexus file
            nexus = "{} {}\n".format(4, len(fullseqs[a])) + "\n".join(nex)    

            ## infer ML tree
            treeorder = self.run_tree_inference(
                nexus, "{}.{}".format(hashval, ridx))

            ## add to list
            weights.append(treeorder)

    ## cleanup - remove all files with the hash val
    rfiles = glob.glob(os.path.join(tempfile.tempdir, "*{}*".format(hashval)))
    for rfile in rfiles:
        if os.path.exists(rfile):
            os.remove(rfile)

    ## return result as weights for the set topologies.
    trees = ["ABCD", "ACBD", "ADBC"]
    wdict = {i: float(weights.count(i)) / len(weights) for i in trees}
    return wdict



def get_order(tre):
    """
    return tree order
    """
    anode = tre.treenode.search_nodes(">A")[0]
    sister = anode.get_sisters()[0]
    sisters = (anode.name[1:], sister.name[1:])
    others = [i for i in list("ABCD") if i not in sisters]
    return sorted(sisters) + sorted(others)


def count_var(nex):
    """
    count number of sites with cov=4, and number of variable sites.
    """
    arr = np.array([list(i.split()[-1]) for i in nex])
    miss = np.any(arr == "N", axis=0)
    nomiss = arr[:, ~miss]
    nsnps = np.invert(np.all(nomiss == nomiss[0, :], axis=0)).sum()
    return nomiss.shape[1], nsnps


