#!/usr/bin/env python

from __future__ import print_function
import os
import gzip
import copy
import subprocess
import itertools
import numpy as np
import pandas as pd
from ipyrad.assemble.util import Params
from ipyrad.analysis.tetrad import get_spans
from ipyrad.assemble.write_outfiles import reftrick, GETCONS

try:
    import toytree
except ImportError:
    print("""
    Required dependency 'toytree' is missing. You can install it
    with the following command:

        conda install toytree -c eaton-lab

    """)
    raise

## alias
OPJ = os.path.join


class Treemix(object):
    """
    Treemix analysis utility function.

    Parameters:
    -----------
    data: str
        The .phy formated sequence file. An alias for '-s'. 
    name: str
        The name for this run. An alias for '-n'.
    workdir: str
        The output directory for results. An alias for '-w'. 
    imap: dict
        A dictionary mapping 


    Attributes:
    -----------
    params: dict
        parameters for this raxml run
    cmd: 
        returns the command string to run raxml

    Functions:
    ----------
    write_output_file()
        writes a gzipped output file in treemix format
    run()
        submits a raxml job to locally or on an ipyparallel client cluster. 

    """    

    ## init object for params
    def __init__(self,
        data,
        name="test",
        workdir="analysis-treemix", 
        imap=None,
        mapfile=None,
        minmap=None,
        *args, 
        **kwargs):

        ## path attributes
        self.name = name
        self.data = data
        self.mapfile = mapfile

        ## if not imap then it will be set to 1
        if imap:
            self.imap = imap

        ## if not minmap then none
        if minmap:
            self.minmap = minmap

        ## params dict
        self.params = Params()
        self.params.k = 0
        self.params.m = 0
        self.params.g = (None, None)
        self.params.bootstrap = 0
        self.params.cormig = 0
        self.params.climb = 0
        self.params.noss = 0
        self.params.seed = np.random.randint(0, int(1e9))
        self.params.root = None
        self.params.se = 0
        self.params.global_ = 0

        ## if mapfile then parse it to an array
        if self.mapfile:
            with open(self.mapfile) as inmap:
                maparr = np.genfromtxt(inmap)[:, [0, 3]].astype(np.uint64)
                spans = np.zeros((maparr[-1, 0], 2), np.uint64)
                spans = get_spans(maparr, spans)
                self.maparr = spans
                self.nsites = spans.shape[0]
        else:
            self.maparr = None
            
        ## make workdir if it does not exist
        if workdir:
            self.workdir = os.path.abspath(os.path.expanduser(workdir))
        else:
            self.workdir = OPJ(os.path.curdir, "analysis-structure")
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        ## set params
        notparams = set(["workdir", "name", "data", "minmap", "imap"])
        for key in set(kwargs.keys()) - notparams:
            self.params[key] = kwargs[key]

        ## check binary
        self._get_binary()

        ## results files
        self.files = Params()
        self.files.cov = OPJ(self.workdir, self.name+".cov.gz")
        self.files.covse = OPJ(self.workdir, self.name+".covse.gz")
        self.files.edges = OPJ(self.workdir, self.name+".edges.gz")
        self.files.llik = OPJ(self.workdir, self.name+".llik")
        self.files.modelcov = OPJ(self.workdir, self.name+".modelcov.gz")
        self.files.treeout = OPJ(self.workdir, self.name+".treeout.gz")        
        self.files.vertices = OPJ(self.workdir, self.name+".vertices.gz")

        ## results
        self.results = Params()
        self.results.tree = ""
        self.results.admixture = []
        self.results.cov = []
        self.results.covse = []
        self.results.modelcov = []   


    @property
    def _command_list(self):
        """ build the command list """

        ## base args
        cmd = [self.params.binary, 
                "-i", OPJ(self.workdir, self.name+".treemix.in.gz"),
                "-o", OPJ(self.workdir, self.name),
                ]

        ## addon params
        args = []
        for key, val in self.params:
            if key not in ["minmap", "binary"]:
                if key == "g":
                    if val[0]:
                        args += ["-"+key, str(val[0]), str(val[1])]
                elif key == "global_":
                    if val:
                        args += ["-"+key[:-1]]
                elif key in ["se", "global", "noss"]: 
                    if val:
                        args += ["-"+key]
                else:
                    if val:
                        args += ["-"+key, str(val)]

        return cmd+args



    @property
    def command(self):
        """ returns command as a string """
        return " ".join(self._command_list)
    

    def _get_binary(self):
        self.params.binary = "treemix"


    def _subsample(self):
        """ returns a subsample of unlinked snp sites """
        spans = self.maparr
        samp = np.zeros(spans.shape[0], dtype=np.uint64)
        for i in xrange(spans.shape[0]):
            samp[i] = np.random.randint(spans[i, 0], spans[i, 1], 1)
        return samp


 
    def copy(self, name):
        """ 
        Returns a copy of the treemix object with the same parameter settings
        but with the files attributes cleared, and with a new 'name' attribute. 
        
        Parameters
        ----------
        name (str):
            A name for the new copied treemix bject that will be used for the 
            output files created by the object. 

        """

        ## make deepcopy of self.__dict__ but do not copy async objects
        subdict = {i:j for i, j in self.__dict__.iteritems() if i != "asyncs"}
        newdict = copy.deepcopy(subdict)

        ## make back into a Treemix object
        #if name == self.name:
        #    raise Exception("new object name must be different from its parent")

        newobj = Treemix(
            data=newdict["data"],
            name=name,
            workdir=newdict["workdir"],
            imap={i:j for i, j in newdict["imap"].items()},
            mapfile=newdict['mapfile'],
            minmap={i:j for i, j in newdict["minmap"].items()},
            seed=np.random.randint(0, int(1e9)),
            )

        ## update special dict attributes but not files
        for key, val in newobj.params.__dict__.iteritems():
            newobj.params.__setattr__(key, self.params.__getattribute__(key))
        #for key, val in newobj.filters.__dict__.iteritems():
        #    newobj.filters.__setattr__(key, self.filters.__getattribute__(key))

        ## new object must have a different name than it's parent
        return newobj



    def write_output_file(self, quiet=False):

        ## read in snpsfile and subsample
        with open(self.data) as ifile:
            
            ## parse file for header and data
            _data = ifile.readlines()
            ntaxa, nsites = [int(i) for i in _data[0].strip().split()]

            ## parse names and seqs
            names = []
            seqs = []
            for line in _data[1:]:
                name, seq = line.strip().split()
                names.append(name)
                seqs.append(list(seq))
            seqarr = np.array(seqs)
            
            ## index names for paired rows
            idxmap = {key:[names.index(i) for i in self.imap[key]] for key in self.imap}
            idxmap = {k: list(itertools.chain(*[[i*2, i*2+1] for i in v])) \
                      for k, v in idxmap.iteritems()}
            
            ## get allele counts
            pcounts, scounts = _get_counts(seqarr)
            
            ## [optional] subsample by mapfile
            if self.mapfile:
                sub = self._subsample()
                pcounts = pcounts[:, sub]
                scounts = scounts[:, sub]
            
            ## mask by minmap 
            sub = _minsample(pcounts, scounts, self.minmap, idxmap)
            pcounts = pcounts[:, sub]
            scounts = scounts[:, sub]
            
            ## print to file
            fdict = {}
            for key, idxs in idxmap.iteritems():
                counts = [
                    "{},{}".format(i,j) for i,j in \
                         zip(np.sum(scounts[idxs, :], axis=0),
                             np.sum(pcounts[idxs, :], axis=0))
                    ]
                fdict[key] = counts
            ## order fdict names
            fnames = sorted(fdict.keys())

        nwritten = sub.sum()
        if not nwritten:
            print("Error: No SNPs passed filtering")
            return 

        if not quiet:
            print('ntaxa {}; nSNPs total {}; nSNPs written {}'\
                  .format(ntaxa, nsites, nwritten))

        with open(OPJ(self.workdir, self.name+".treemix.in.gz"), 'w') as outf:
            outf.write(" ".join(fnames)+"\n")
            farr = np.array([fdict[i] for i in fnames]).T
            np.savetxt(outf, farr, delimiter=" ", fmt="%s")



    def run(self, ipyclient=None, quiet=True):

        ## call command
        self.write_output_file(quiet=True)
        self._call_treemix()
        self._parse_results()

        #if ipyclient:
        #    lbview = ipyclient.load_balanced_view()
        #    command = self._command_list, 
        #    self.async = lbview.apply()
        
        ## check for errors
        ## ...
        ## parse result files
        #with lbview(temp_flags=self.async):
        #    lbview.apply(self._parse_results())



    def _call_treemix(self):
        proc = subprocess.Popen(
            self._command_list,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE)
        stdout, err = proc.communicate()
        return stdout



    def _parse_results(self):
        ## get tree and admix from output files 
        with gzip.open(self.files.treeout) as tmp:
            data = tmp.readlines()

            ## store the tree
            k = 0
            if self.params.noss:
                k += 1
            self.results.tree = data[k].strip()
            self.results.admixture = []

            ## get admix events
            for adx in data[k+1:]:
                weight, jweight, jse, pval, clade1, clade2 = adx.strip().split()
                self.results.admixture.append((clade1, clade2, weight, jweight, jse, pval))

        ## get a toytree
        tre = toytree.tree(self.results.tree)

        ## order admixture 
        for aidx in range(len(self.results.admixture)):
            admix = self.results.admixture[aidx]

            source = toytree.tree(admix[0]+";")
            if len(source.tree) == 1:
                sodx = tre.tree.search_nodes(name=source.tree.name)[0].idx
            else:
                sodx = tre.tree.get_common_ancestor(source.get_tip_labels()).idx

            sink = toytree.tree(admix[1]+";")
            if len(sink.tree) == 1:
                sidx = tre.tree.search_nodes(name=sink.tree.name)[0].idx
            else:
                sidx = tre.tree.get_common_ancestor(sink.get_tip_labels()).idx

            self.results.admixture[aidx] = (
                int(sodx),
                float(admix[0].rsplit(":", 1)[1]), 
                int(sidx),
                float(admix[1].rsplit(":", 1)[1]), 
                float(admix[2]), 
                #float(admix[3]), 
                #float(admix[4]),
                #float(admix[5]),
                )

        ## parse the cov
        names = tre.get_tip_labels()
        self.results.cov = _parse_matrix(self.files.cov, names)
        self.results.covse = _parse_matrix(self.files.covse, names)
        self.results.modelcov = _parse_matrix(self.files.modelcov, names)



def _parse_matrix(mfile, names):
    dat = pd.read_csv(mfile,
                  sep=" ", 
                  header=None, 
                  index_col=0, 
                  skiprows=1)

    ## sort into tree label order
    dat.columns = dat.index
    dat = dat.ix[names]
    dat = dat.T.ix[names]
    return dat.as_matrix()



def _minsample(pcounts, scounts, minmap, idxmap):
    ## mask for minsamp
    fcounts = pcounts + scounts
    minmapped = np.array(
        [fcounts[idxs, :].sum(axis=0) >= minmap[key]*2 for 
         key, idxs in idxmap.iteritems()]
    )
    mask = np.invert(np.any(minmapped == False, axis=0))

    ## now subsample
    return mask



def _get_counts(seqarr):

    outarr = reftrick(seqarr[:].view(np.uint8), GETCONS).view("S1")

    ## allele arrays
    seq0 = seqarr.copy()
    seq1 = seqarr.copy()

    seq0[seq0=="R"] = "G"
    seq0[seq0=="K"] = "G"
    seq0[seq0=="S"] = "G"
    seq0[seq0=="Y"] = "T"
    seq0[seq0=="W"] = "T"
    seq0[seq0=="M"] = "C"

    seq1[seq1=="R"] = "A"
    seq1[seq1=="K"] = "T"
    seq1[seq1=="S"] = "C"
    seq1[seq1=="Y"] = "C"
    seq1[seq1=="W"] = "A"
    seq1[seq1=="M"] = "A"

    ## fill array 
    ntaxa, nsites = seqarr.shape
    priarr = np.zeros((ntaxa*2, nsites), dtype=np.uint8)
    secarr = np.zeros((ntaxa*2, nsites), dtype=np.uint8)

    ## add 2 if allele is the outgroup
    priarr[::2] += seq0 == outarr[:, 0]
    priarr[1::2] += seq1 == outarr[:, 0]
    secarr[::2] += seq0 == outarr[:, 1]
    secarr[1::2] += seq1 == outarr[:, 1]

    ## invert so zero is ancestral and >1 is derived
    return priarr, secarr

    

def _call_treemix(command_list):
    """ call the command as sps """
    proc = subprocess.Popen(
        command_list,
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE
        )
    comm = proc.communicate()
    return comm


