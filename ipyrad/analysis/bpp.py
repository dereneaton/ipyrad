#!/usr/bin/python 

""" convert loci file to bpp format input files """

import os
import sys
import glob
import copy
import itertools
import subprocess
import numpy as np
from toytree import ete3mini as ete
from collections import Counter
from ipyrad.assemble.util import DUCT
#import pandas as pd



class Bpp(object):
    """
    BPP analysis utility function for creating input files, setting parameters, 
    and submitting bpp jobs to run on a parallel cluster. Converts loci 
    file format data to bpp file format, i.e., concatenated phylip-like
    format, and produces imap and ctl input files for bpp. The main 
    function to run jobs is submit_bpp_jobs(). 

    Parameters:
    -----------
    locifile:
        A .loci file produced by ipyrad.
    imap:
        A Python dictionary with 'species' names as keys, and lists of sample
        names for the values. Any sample that is not included in the imap
        dictionary will be filtered out of the data when converting the .loci
        file into the bpp formatted sequence file. Each species in the imap
        dictionary must also be present in the input 'guidetree'.
    guidetree:
        A newick string species tree hypothesis [e.g., (((a,b),(c,d)),e);]
        All taxa in the imap dictionary must also be present in the guidetree.

    Attributes:
    -----------
    params (dict):
        parameters for bpp analses used to create .ctl file
    filters (dict):
        parameters used to filter loci that will be included in the analysis.

    Functions:
    ----------
    submit_bpp_jobs():
        See documentation string for this function.
    write_bpp_files():
        See documentation string for this function.

    Optional parameters (object.params):
    --------------------
    infer_sptree:
        Default=0, only infer parameters on a fixed species tree. If 1, then the
        input tree is treated as a guidetree and tree search is employed to find
        the best tree. The results will include support values for the inferred
        topology.
    infer_delimit:
        Default=0, no delimitation. If 1 then splits in the tree that separate
        'species' will be collapsed to test whether fewer species are a better
        fit to the data than the number in the input guidetree.
    delimit_alg:
        Species delimitation algorithm is a two-part tuple. The first value
        is the algorithm (0 or 1) and the second value is a tuple of arguments
        for the given algorithm. See other ctl files for examples of what the
        delimitation line looks like. This is where you can enter the params
        (e.g., alpha, migration) for the two different algorithms.
        For example, the following args would produce the following ctl lines:
         alg=0, epsilon=5
         > delimit_alg = (0, 5)
         speciesdelimitation = 1 0 5

         alg=1, alpha=2, migration=1
         > delimit_alg = (1, 2, 1)
         speciesdelimitation = 1 1 2 1

         alg=1, alpha=2, migration=1, diagnosis=0, ?=1
         > delimit_alg = (1, 2, 1, 0, 1)
         speciesdelimitation = 1 1 2 1 0 1
    seed:
        A random number seed at start of analysis.
    burnin:
        Number of burnin generations in mcmc
    nsample:
        Number of mcmc generations to run.
    sampfreq:
        How often to sample from the mcmc chain.
    thetaprior:
        Prior on theta (4Neu), gamma distributed. mean = a/b. e.g., (5, 5)
    tauprior
        Prior on root tau, gamma distributed mean = a/b. Last number is dirichlet
        prior for other taus. e.g., (4, 2, 1)
    usedata:
        If false inference proceeds without sequence data (can be used to test
        the effect of priors on the tree distributions).
    cleandata:
        If 1 then sites with missing or hetero characters are removed.
    finetune:
        See bpp documentation.

    """    

    ## init object for params
    def __init__(self,
        locifile, 
        guidetree, 
        imap, 
        workdir, 
        *args, 
        **kwargs):

        ## path attributes
        self.workdir = workdir
        self.asyncs = []
        self._kwargs = {
                "maxloci": None,
                "minmap": None,
                "minsnps": 0,
                "infer_sptree": 0,
                "infer_delimit": 0,
                "delimit_alg": (0, 5),
                "seed": 12345,
                "burnin": 1000,
                "nsample": 10000,
                "sampfreq": 2,
                "thetaprior": (2, 2000),
                "tauprior": (2, 2000, 1),
                "usedata": 1,
                "cleandata": 0,
                "finetune": (0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                "copied": False,
            }
        self._kwargs.update(kwargs)

        ## check workdir
        if self.workdir:
            self.workdir = os.path.abspath(self.workdir)
        else:
            self.workdir = os.path.curdir        
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        ## set the guidetree
        self.tree = ete.Tree(guidetree)

        ## parsing attributes
        self.imap = {i:j for i, j in imap.items()}
        if not self.imap:
            self.imap = {i:i for i in tree.get_leaf_names()}

        ## update stats if alleles instead of loci 
        if ('.alleles.loci' in locifile) and (not self._kwargs['copied']):
            ## add 0/1 to names
            keys = self.imap.keys()
            for key in keys:
                oldvals = self.imap[key]
                newvals = []
                for val in oldvals:
                    newvals += [val+"_0", val+"_1"]
                self.imap[key] = newvals

            ## double the minmap
            if self._kwargs['minmap']:
                self._kwargs["minmap"] = \
                {key: val*2 for key, val in self._kwargs['minmap'].items()}

        ## checks
        assert isinstance(self.imap, dict), "you must enter an IMAP dictionary"
        assert set(self.imap.keys()) == set(self.tree.get_leaf_names()), \
               "IMAP keys must match guidetree names"

        ## filters
        self.filters = Params()
        self.filters.minmap = self._kwargs["minmap"]
        self.filters.maxloci = self._kwargs["maxloci"]
        self.filters.minsnps = self._kwargs["minsnps"]

        ## set bpp parameters with defaults
        self.params = Params()
        notparams = set(["workdir", "maxloci", "minmap", "minsnps"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## results files
        self.files = Params()
        self.files.locifile = locifile
        self.files.mcmcfiles = []
        self.files.outfiles = []
        


    def submit_bpp_jobs(self, 
        prefix,
        nreps, 
        ipyclient, 
        seed=12345,
        quiet=False,
        randomize_order=False,
        ):
        """
        Submits bpp jobs to run on a cluster (ipyparallel Client). 
        The seed for the random number generator if not set is randomly 
        drawn, and if multiple reps are submitted (nreps>1) then each will 
        draw a subsequent random seeds after that. An ipyclient connection 
        is required. Asynchronous result objects are stored in the bpp 
        object submitting the jobs. 
        
        Parameters:
        -----------
        prefix (str):
            a prefix name for output files produced by these bpp runs.
        nreps (int):
            submits nreps replicate jobs to the cluster each with a different 
            random seed drawn starting from the starting seed. 
        ipyclient (ipyparallel.Client)
            an ipyparallel.Client object connected to a running cluster. 
        seed (int):
            an integer used to initiate the random number generator
        quiet (bool):
            whether to print that the jobs have been submitted
        randomize_order (bool):
            if True then when maxloci is set this will randomly sample a 
            different set of N loci in each replicate, rather than sampling
            just the first N loci < maxloci. 
        """
        ## initiate random seed
        np.random.seed(seed)

        ## prepare files
        self._name = prefix
        self._write_seqfile()
        self._write_mapfile()

        ## load-balancer
        lbview = ipyclient.load_balanced_view()

        ## send jobs
        asyncs = []
        for job in xrange(nreps):

            ## change seed and ctl for each rep
            self.params.seed = np.random.randint(0, 1e9, 1)[0]
            ctlfile = self._write_ctlfile(rep=job)
            if randomize_order:
                self._write_seqfile(randomize_order)

            ## submit to engines
            async = lbview.apply(_call_bpp, ctlfile)
            asyncs.append(async)
            self.asyncs.append(async)

        if not quiet:
            sys.stderr.write("submitted {} bpp jobs [{}] ({} loci)\n"\
                             .format(nreps, prefix, self._nloci))



    def write_bpp_files(self, prefix, randomize_order=False, quiet=False):
        """ 
        Writes bpp files (.ctl, .seq, .imap) to the working directory. 

        Parameters:
        ------------
        prefix (str):
            a prefix name for the output files
        randomize_order (bool):
            whether to randomize the locus order, this will allow you to 
            sample different subsets of loci in different replicates when
            using the filters.maxloci option.
        quiet (bool):
            whether to print info to stderr when finished.
        """

        ## remove any old jobs with this same job name
        oldjobs = glob.glob(os.path.join(self.workdir, prefix+"*.ctl.txt"))
        for job in oldjobs:
            os.remove(job)

        ## check params types
        ## ...

        ## write tmp files for the job
        self._name = prefix
        self._write_seqfile(True)
        self._write_mapfile(True)
        self._write_ctlfile()

        if not quiet:
            sys.stderr.write("input files created for job {} ({} loci)\n"\
                             .format(prefix, self._nloci))



    def _write_mapfile(self, name=False):
        ## write the imap file:
        if name:
            self.mapfile = os.path.join(self.workdir, self._name+".imapfile.txt")
        else:
            self.mapfile = os.path.join(self.workdir, "tmp.imapfile.txt")

        with open(self.mapfile, 'w') as mapfile:
            data = ["{:<30} {}".format(val, key) for key \
                    in sorted(self.imap) for val in self.imap[key]]
            mapfile.write("\n".join(data))



    def _write_seqfile(self, name=False, randomize_order=False):

        ## handles
        if name: 
            self.seqfile = os.path.join(self.workdir, self._name+".seqfile.txt")
        else:
            self.seqfile = os.path.join(self.workdir, "tmp.seqfile.txt")
        seqfile = open(self.seqfile, 'w')
        with open(self.files.locifile) as infile:
            loci = infile.read().strip().split("|\n")
            nloci = len(loci)
            if randomize_order:
                np.random.shuffle(loci)

        ## all samples
        samples = list(itertools.chain(*self.imap.values()))

        ## iterate over loci, printing to outfile
        nkept = 0
        for iloc in xrange(nloci):
            lines = loci[iloc].split("//")[0].split()
            names = lines[::2]
            names = ["^"+i for i in names]
            seqs = [list(i) for i in lines[1::2]]
            seqlen = len(seqs[0])
            ## whether to skip this locus based on filters below
            skip = 0

            ## if minmap filter for sample coverage
            if self.filters.minmap:
                covd = {}
                for group, vals in self.imap.items():
                    covd[group] = sum(["^"+i in names for i in vals])
                ## check that coverage is good enough
                if not all([covd[group] >= self.filters.minmap[group] for group \
                            in self.filters.minmap]):
                    skip = 1

            ## too many loci?
            if (not skip) and (self.filters.maxloci):
                if nkept >= self.filters.maxloci:
                    skip = 1

            ## if minsnps filter for snps
            if (not skip) and (self.filters.minsnps):
                npis = 0
                arr = np.array(seqs)
                arr = np.zeros((arr.shape[0]*2, arr.shape[1]), dtype="S1")
                for row in xrange(len(seqs)):
                    fillrow = 2 * row
                    arr[fillrow] = [DUCT[i][0] for i in seqs[row]]
                    arr[fillrow+1] = [DUCT[i][1] for i in seqs[row]]
                    
                for col in xrange(arr.shape[1]):
                    bases = arr[:, col]
                    bases = bases[bases != "N"]
                    bases = bases[bases != "-"]
                    counts = Counter(bases)
                    counts = counts.most_common(2)
                    if len(counts) == 2:
                        if counts[1][1] > 1:
                            npis += 1
                if npis < self.filters.minsnps:
                    skip = 1

            ## build locus as a string
            if not skip:
                ## convert to phylip with caret starter and replace - with N.
                data = ["{:<30} {}".format(i, "".join(k).replace("-", "N")) for \
                    (i, k) in zip(names, seqs) if i[1:] in samples]

                ## if not empty, write to the file
                if data:
                    seqfile.write("{} {}\n\n{}\n\n"\
                               .format(len(data), seqlen, "\n".join(data)))
                    nkept += 1

        ## close up shop
        self._nloci = nkept
        del loci
        seqfile.close()    



    def _write_ctlfile(self, rep=None):
        """ write outfile with any args in argdict """

        ## A string to store ctl info
        ctl = []

        ## write the top header info
        ctl.append("seed = {}".format(self.params.seed))
        ctl.append("seqfile = {}".format(self.seqfile))
        ctl.append("Imapfile = {}".format(self.mapfile))

        path = os.path.join(self.workdir, self._name)
        if isinstance(rep, int):
            mcmcfile = "{}-r{}.mcmc.txt".format(path, rep)
            outfile = "{}-r{}.out.txt".format(path, rep)
        else:
            mcmcfile = "{}.mcmc.txt".format(path, rep)
            outfile = "{}.out.txt".format(path, rep)

        if mcmcfile not in self.files.mcmcfiles:
            self.files.mcmcfiles.append(mcmcfile)
        if outfile not in self.files.outfiles:
            self.files.outfiles.append(outfile)

        ctl.append("mcmcfile = {}".format(mcmcfile))
        ctl.append("outfile = {}".format(outfile))

        ## number of loci (checks that seq file exists and parses from there)
        ctl.append("nloci = {}".format(self._nloci))
        ctl.append("usedata = {}".format(self.params.usedata))
        ctl.append("cleandata = {}".format(self.params.cleandata))

        ## infer species tree
        if self.params.infer_sptree:
            ctl.append("speciestree = 1 0.4 0.2 0.1")
        else:
            ctl.append("speciestree = 0")

        ## infer delimitation (with algorithm 1 by default)
        ctl.append("speciesdelimitation = {} {} {}"\
                   .format(self.params.infer_delimit, 
                           self.params.delimit_alg[0],
                           " ".join([str(i) for i in self.params.delimit_alg[1:]]) 
                           )
                   )
        ## get tree values
        nspecies = str(len(self.imap))
        species = " ".join(sorted(self.imap))
        ninds = " ".join([str(len(self.imap[i])) for i in sorted(self.imap)])
        ctl.append(SPECIESTREE.format(nspecies, species, ninds, self.tree.write(format=9)))

        ## priors
        ctl.append("thetaprior = {} {}".format(*self.params.thetaprior))
        ctl.append("tauprior = {} {} {}".format(*self.params.tauprior))

        ## other values, fixed for now
        ctl.append("finetune = 1: {}".format(" ".join([str(i) for i in self.params.finetune])))
        #CTL.append("finetune = 1: 1 0.002 0.01 0.01 0.02 0.005 1.0")
        ctl.append("print = 1 0 0 0")
        ctl.append("burnin = {}".format(self.params.burnin))
        ctl.append("sampfreq = {}".format(self.params.sampfreq))
        ctl.append("nsample = {}".format(self.params.nsample))

        ## write out the ctl file
        if isinstance(rep, int):
            ctlhandle = "{}-r{}.ctl.txt".format(os.path.join(self.workdir, self._name), rep)
        else:
            ctlhandle = "{}.ctl.txt".format(os.path.join(self.workdir, self._name))
        with open(ctlhandle, 'w') as out:
            out.write("\n".join(ctl))

        return ctlhandle



    def copy(self):
        """ 
        returns a copy of the bpp object with the same parameter settings
        but with the files.mcmcfiles and files.outfiles attributes cleared. 
        """
        ## make deepcopy of self.__dict__ but do not copy async objects
        subdict = {i:j for i,j in self.__dict__.iteritems() if i != "asyncs"}
        newdict = copy.deepcopy(subdict)

        ## make back into a bpp object
        newobj = Bpp(
            locifile=newdict["files"].locifile,
            workdir=newdict["workdir"],
            guidetree=newdict["tree"].write(),
            imap={i:j for i, j in newdict["imap"].items()},
            copied=True,
            )

        ## update special dict attributes but not files
        for key, val in newobj.params.__dict__.iteritems():
            newobj.params.__setattr__(key, self.params.__getattribute__(key))
        for key, val in newobj.filters.__dict__.iteritems():
            newobj.filters.__setattr__(key, self.filters.__getattribute__(key))

        return newobj



def _call_bpp(ctlfile):
    ## call the command
    proc = subprocess.Popen(
        ["bpp", ctlfile], 
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE
        )
    comm = proc.communicate()
    return comm



class Params(object):
    """ 
    A dict-like object for storing params values with a custom repr
    """
    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())
        _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
        for key in keys:
            _repr += _printstr.format(key, str(self[key]))
        return _repr



class Result_00(object):
    """ parse bpp results object for 00 scenario """
    pass



def _get_bpp_results(ofiles):
    
    ## final all outfiles associated with this job name
    ofiles = glob.glob()

    cols = []
    for ofile in ofiles:
        with open(ofile) as infile:
            dat  = infile.read()
        lastbits = dat.split("bpp.mcmc.txt\n\n")[1:]
        results = lastbits[0].split("\n\n")[0].split()
        dat = np.array(results[3:]).reshape(8, 4)
        cols.append(dat[:, 3].astype(float))
    cols = np.array(cols)
    cols = cols.sum(axis=0) / 10.
    dat[:, 3] = cols.astype(str)
    dd = pd.DataFrame(dat[:, 1:])
    dd.columns = ["delim", "prior", "posterior"]
    nspecies = 1 + np.array([list(i) for i in dat[:, 1]], dtype=int).sum(axis=1)
    dd["nspecies"] = nspecies
    return dd
    
    


## GLOBALS
IMAP_REQUIRED = """\
  An IMAP dictionary is required as input to group samples into 'species'.
  Example: {"A":[tax1, tax2, tax3], "B":[tax4, tax5], "C":[tax6, tax7]}
  """

KEYS_DIFFER = """
  MINMAP and IMAP keys must be identical.
  """

PDREAD_ERROR = """\
  Error in trait file: all data must be quantitative (int or floats)
  and missing data must be coded as NA. See the example notebook. If
  sample names are in a data column this will cause an error, try
  passing the argument 'index_col=0' to pandas.read_csv().
  """

SPECIESTREE = """\
species&tree = {} {}
                 {}
                 {}"""  

