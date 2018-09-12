#!/usr/bin/env python

"convert loci file to bpp format input files"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard lib
import os
import sys
import glob
import copy
import subprocess as sps
from collections import Counter

import numpy as np
import pandas as pd
from ipyrad.analysis.utils import Params, progressbar
from ipyrad.assemble.utils import IPyradError, IPyradWarningExit

# DUCT from where?



# try:
#     ## when you have time go back and set attrubutes on toytrees
#     from toytree import ete3mini as ete
# except ImportError:
#     raise IPyradWarningExit("""
#     Error: bpp requires the dependency 'toytree', which we haven't yet
#     included in the ipyrad installation. For now, you can install toytree
#     using conda with the following command: 

#     conda install toytree -c eaton-lab
#     """)

MISSING_IMPORTS = """
You are missing required packages to use ipa.bpp.
First run the following conda install command:

  conda install -c ipyrad bpp
"""


class Bpp(object):
    """
    BPP analysis utility function for creating input files, setting parameters, 
    and submitting bpp jobs to run on a parallel cluster. Converts loci 
    file format data to bpp file format, i.e., concatenated phylip-like
    format, and produces imap and ctl input files for bpp. The main 
    functions are 'write_bpp_files()' and 'run()'.

    Parameters:
    -----------
    name: str
        A name for this analysis object.

    data: str
        The path to a .loci or .alleles.loci file produced by ipyrad.

    imap: dict
        A Python dictionary with 'species' names as keys, and lists of sample
        names for the values. Any sample that is not included in the imap
        dictionary will be filtered out of the data when converting the .loci
        file into the bpp formatted sequence file. Each species in the imap
        dictionary must also be present in the input 'guidetree'.

    guidetree: str
        A newick string species tree hypothesis [e.g., (((a,b),(c,d)),e);]
        All taxa in the imap dictionary must also be present in the guidetree.
        Tree can also be a filename of a newick string.

    load_existing_results: bool
        If True then any existing results files saved in the working
        directory and with the entered name will be loaded and attached 
        to this object. This is useful if returning to a notebook later and 
        you want to summarize results. 


    Attributes:
    -----------
    params (dict):
        parameters for bpp analses used to create .ctl file
    filters (dict):
        parameters used to filter loci that will be included in the analysis.
    workdir (str)
        Default output path is "./analysis-bpp" which will be created if it 
        does not exist. If you enter a different output directly it will 
        similarly be created if it does not exist.

    Functions:
    ----------
    run():
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
        name,
        data=None,
        workdir="analysis-bpp", 
        guidetree=None, 
        imap=None, 
        load_existing_results=False,
        *args, 
        **kwargs):

        # store args
        self.name = name
        self.data = data
        self.workdir = workdir
        self.guidetree = guidetree
        self.imap = imap

        # update kwargs 
        self.asyncs = []
        self._kwargs = {
            "binary": "bpp",
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

        # run checks
        self.check_binary()
        self.check_files()


    def check_binary(self):
        ## check that bpp is installed and in path
        for binary in [self._kwargs["binary"]]:
            cmd = ['which', binary]
            proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
            comm = proc.communicate()[0]
            if not comm:
                raise IPyradError(MISSING_IMPORTS)

    # needs closer checking since the py3 reboot
    def check_files(self):
        # support for legacy args
        if self._kwargs.get("locifile"):
            self.data = self._kwargs.get("locifile")
        if not self.data:
            raise IPyradWarningExit(
                "must enter a 'data' argument (an ipyrad .loci file).")

        # set the guidetree
        if not self.guidetree:
            raise IPyradWarningExit(
                "must enter a 'guidetree' argument (a newick file or string).")
        self.tree = ete.Tree(self.guidetree)

        # check workdir
        if self.workdir:
            self.workdir = os.path.abspath(os.path.expanduser(self.workdir))
        else:
            self.workdir = os.path.join(os.path.curdir, "analysis-bpp")
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # parsing imap dictionary, or create simple 1-1 mapping
        if not self.imap:
            self.imap = {i: [i] for i in self.tree.get_leaf_names()}
        else:
            self.imap = {}
            for key, val in self.imap.items():
                if isinstance(val, (int, str)):
                    self.imap[key] = [str(val)]
                elif isinstance(val, list):
                    self.imap[key] = val
                else:
                    raise IPyradWarningExit(
                        "imap dictionary is not properly formatted")

        # update stats if alleles instead of loci 
        if not self._kwargs["minmap"]:
            self._kwargs["minmap"] = {i: 1 for i in self.tree.get_leaf_names()}

        if ('.alleles.loci' in self.data) and (not self._kwargs['copied']):
            # add 0/1 to names
            keys = self.imap.keys()
            for key in keys:
                oldvals = self.imap[key]
                newvals = []
                for val in oldvals:
                    newvals += [val + "_0", val + "_1"]
                self.imap[key] = newvals

            # double the minmap (copied attribute protects from double 2X)
            self._kwargs["minmap"] = \
                {key: val * 2 for key, val in self._kwargs['minmap'].items()}

        ## checks
        assert isinstance(self.imap, dict), "you must enter an IMAP dictionary"
        assert set(self.imap.keys()) == set(self.tree.get_leaf_names()), \
               "IMAP keys must match guidetree names: \n{}\n{}"\
               .format(self.imap.keys(), self.tree.get_leaf_names())

        ## filters
        self.filters = Params()
        self.filters.minmap = self._kwargs["minmap"]
        self.filters.maxloci = self._kwargs["maxloci"]
        self.filters.minsnps = self._kwargs["minsnps"]

        ## set bpp parameters with defaults
        self.params = Params()
        notparams = set(["workdir", "maxloci", "minmap", "minsnps", "copied"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## results files
        self.files = Params()
        self.files.data = self.data
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []
        
        ## load existing results files for this named bpp object if they exist
        if self.load_existing_results:
            self._load_existing_results(self.name, workdir=self.workdir)



    def _load_existing_results(self, name, workdir):
        """
        Load existing results files for an object with this workdir and name. 
        This does NOT reload the parameter settings for the object...
        """
        ## get mcmcs
        path = os.path.realpath(os.path.join(self.workdir, self.name))
        mcmcs = glob.glob(path+"_r*.mcmc.txt")
        outs = glob.glob(path+"_r*.out.txt")
        trees = glob.glob(path+"_r*.tre")

        for mcmcfile in mcmcs:
            if mcmcfile not in self.files.mcmcfiles:
                self.files.mcmcfiles.append(mcmcfile)
        for outfile in outs:
            if outfile not in self.files.outfiles:
                self.files.outfiles.append(outfile)
        for tree in trees:
            if tree not in self.files.treefiles:
                self.files.treefiles.append(tree)        


    def run(
        self,
        ipyclient, 
        nreps=1, 
        quiet=False,
        randomize_order=False,
        force=False,
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
        nreps (int):
            submits nreps replicate jobs to the cluster each with a different 
            random seed drawn starting from the starting seed. 
        ipyclient (ipyparallel.Client)
            an ipyparallel.Client object connected to a running cluster. 
        quiet (bool):
            whether to print that the jobs have been submitted
        randomize_order (bool):
            if True then when maxloci is set this will randomly sample a 
            different set of N loci in each replicate, rather than sampling
            just the first N loci < maxloci. 
        force (bool):
            Overwrite existing files with the same name. Default=False, skip
            over existing files.
        """

        ## is this running algorithm 00?
        is_alg00 = (not self.params.infer_sptree) and (not self.params.infer_delimit)

        ## clear out pre-existing files for this object
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []
        self.asyncs = []

        ## initiate random seed
        np.random.seed(self.params.seed)

        ## load-balancer
        lbview = ipyclient.load_balanced_view()

        ## send jobs
        for job in range(nreps):

            ## make repname and make ctl filename
            self._name = "{}_r{}".format(self.name, job)
            ctlhandle = os.path.realpath(
                os.path.join(self.workdir, "{}.ctl.txt".format(self._name)))

            ## skip if ctlfile exists
            if (not force) and (os.path.exists(ctlhandle)):
                print("Named ctl file already exists. Use force=True to" \
                    +" overwrite\nFilename:{}".format(ctlhandle))
            else:
                ## change seed and ctl for each rep, this writes into the ctl
                ## file the correct name for the other files which share the 
                ## same rep number in their names.
                #self.params._seed = np.random.randint(0, 1e9, 1)[0]
                self._write_mapfile()
                #if randomize_order:
                self._write_seqfile(randomize_order=randomize_order)
                ctlfile = self._write_ctlfile()
                
                ## submit to engines
                rasync = lbview.apply(_call_bpp, *(self._kwargs["binary"], ctlfile, is_alg00))
                self.asyncs.append(rasync)

                ## save tree file if alg 00
                if is_alg00:
                    self.files.treefiles.append(
                        ctlfile.rsplit(".ctl.txt", 1)[0] + ".tre")

        if self.asyncs and (not quiet):
            sys.stderr.write("submitted {} bpp jobs [{}] ({} loci)\n"\
                             .format(nreps, self.name, self._nloci))



    def write_bpp_files(self, randomize_order=False, quiet=False):
        """ 
        Writes bpp files (.ctl, .seq, .imap) to the working directory. 

        Parameters:
        ------------
        randomize_order (bool):
            whether to randomize the locus order, this will allow you to 
            sample different subsets of loci in different replicates when
            using the filters.maxloci option.
        quiet (bool):
            whether to print info to stderr when finished.
        """

        ## remove any old jobs with this same job name
        self._name = self.name
        oldjobs = glob.glob(os.path.join(self.workdir, self._name+"*.ctl.txt"))
        for job in oldjobs:
            os.remove(job)

        ## check params types
        ## ...

        ## write tmp files for the job
        self._write_seqfile(randomize_order=randomize_order)
        self._write_mapfile()#name=True)
        self._write_ctlfile()

        if not quiet:
            sys.stderr.write("input files created for job {} ({} loci)\n"\
                             .format(self._name, self._nloci))



    def _write_mapfile(self):#, name=False):
        ## write the imap file:
        self.mapfile = os.path.realpath(
            os.path.join(self.workdir, self._name+".imapfile.txt"))
        with open(self.mapfile, 'w') as mapfile:
            data = ["{:<30} {}".format(val, key) for key \
                    in sorted(self.imap) for val in self.imap[key]]
            mapfile.write("\n".join(data))



    def _write_seqfile(self, randomize_order=False):

        ## handles
        self.seqfile = os.path.realpath(
            os.path.join(self.workdir, self._name+".seqfile.txt"))
        seqfile = open(self.seqfile, 'w')
        with open(self.files.data) as infile:
            loci = infile.read().strip().split("|\n")
            nloci = len(loci)
            if randomize_order:
                np.random.shuffle(loci)

        ## all samples
        samples = []
        for mapl in self.imap.values():
            if isinstance(mapl, list):
                samples += mapl
            else:
                samples.append(mapl)

        ## iterate over loci, printing to outfile
        nkept = 0
        for iloc in range(nloci):
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
                for row in range(len(seqs)):
                    fillrow = 2 * row
                    arr[fillrow] = [DUCT[i][0] for i in seqs[row]]
                    arr[fillrow+1] = [DUCT[i][1] for i in seqs[row]]
                    
                for col in range(arr.shape[1]):
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



    def _write_ctlfile(self):#, rep=None):
        """ write outfile with any args in argdict """

        ## A string to store ctl info
        ctl = []

        ## write the top header info
        ctl.append("seed = {}".format(self.params.seed))
        ctl.append("seqfile = {}".format(self.seqfile))
        ctl.append("Imapfile = {}".format(self.mapfile))

        path = os.path.realpath(os.path.join(self.workdir, self._name))
        mcmcfile = "{}.mcmc.txt".format(path)
        outfile = "{}.out.txt".format(path)
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
        ctlhandle = os.path.realpath(
                "{}.ctl.txt".format(os.path.join(self.workdir, self._name)))
        # if isinstance(rep, int):
        #     ctlhandle = os.path.realpath(
        #         "{}-r{}.ctl.txt".format(os.path.join(self.workdir, self._name), rep))
        # else:
        #     ctlhandle = os.path.realpath(
        #         "{}.ctl.txt".format(os.path.join(self.workdir, self._name)))
        with open(ctlhandle, 'w') as out:
            out.write("\n".join(ctl))

        return ctlhandle



    def copy(self, name, load_existing_results=False):
        """ 
        Returns a copy of the bpp object with the same parameter settings
        but with the files.mcmcfiles and files.outfiles attributes cleared, 
        and with a new 'name' attribute. 
        
        Parameters
        ----------
        name (str):
            A name for the new copied bpp object that will be used for the 
            output files created by the object. 

        """

        ## make deepcopy of self.__dict__ but do not copy async objects
        subdict = {i:j for i,j in self.__dict__.iteritems() if i != "asyncs"}
        newdict = copy.deepcopy(subdict)

        ## make back into a bpp object
        if name == self.name:
            raise Exception("new object must have a different 'name' than its parent")
        newobj = Bpp(
            name=name,
            data=newdict["files"].data,
            workdir=newdict["workdir"],
            guidetree=newdict["tree"].write(),
            imap={i:j for i, j in newdict["imap"].items()},
            copied=True,
            load_existing_results=load_existing_results,
            )

        ## update special dict attributes but not files
        for key, val in newobj.params.__dict__.iteritems():
            newobj.params.__setattr__(key, self.params.__getattribute__(key))
        for key, val in newobj.filters.__dict__.iteritems():
            newobj.filters.__setattr__(key, self.filters.__getattribute__(key))

        ## new object must have a different name than it's parent
        return newobj



    def summarize_results(self, individual_results=False):
        """ 
        Prints a summarized table of results from replicate runs, or,
        if individual_result=True, then returns a list of separate
        dataframes for each replicate run. 
        """

        ## return results depending on algorithm

        ## algorithm 00
        if (not self.params.infer_delimit) & (not self.params.infer_sptree):
            if individual_results:
                ## return a list of parsed CSV results
                return [_parse_00(i) for i in self.files.outfiles]
            else:
                ## concatenate each CSV and then get stats w/ describe
                return pd.concat(
                    [pd.read_csv(i, sep='\t', index_col=0) \
                    for i in self.files.mcmcfiles]).describe().T

        ## algorithm 01
        if self.params.infer_delimit & (not self.params.infer_sptree):
            return _parse_01(self.files.outfiles, individual=individual_results)

        ## others
        else:
            return "summary function not yet ready for this type of result"



def _call_bpp(binary, ctlfile, is_alg00):

    ## call the command and block until job finishes
    proc = sps.Popen(["bpp", ctlfile], stderr=sps.STDOUT, stdout=sps.PIPE)
    comm = proc.communicate()

    ## Look for the ~/Figtree.tre file that bpp creates.
    ## This has to be done here to make sure it is instantly run
    ## when the job finishes so other reps won't write over it.
    ## Kludge due to bpp writing forcing the file to $HOME.
    if is_alg00:
        default_figtree_path = os.path.join(os.path.expanduser("~"), "FigTree.tre")
        new_figtree_path = ctlfile.rsplit(".ctl.txt", 1)[0] + ".tre"
        try:
            if os.path.exists(default_figtree_path):
                os.rename(default_figtree_path, new_figtree_path)
        except Exception:
            pass

    return comm


class Result_00(object):
    """ parse bpp results object for 00 scenario """
    pass


def _parse_00(ofile):
    """
    return 00 outfile as a pandas DataFrame
    """
    with open(ofile) as infile:
        ## read in the results summary from the end of the outfile
        arr = np.array(
            [" "] + infile.read().split("Summary of MCMC results\n\n\n")[1:][0]\
            .strip().split())

        ## reshape array 
        rows = 12
        cols = (arr.shape[0] + 1) / rows
        arr = arr.reshape(rows, cols)
        
        ## make into labeled data frame
        df = pd.DataFrame(
            data=arr[1:, 1:], 
            columns=arr[0, 1:], 
            index=arr[1:, 0],
            ).T
        return df
   

def _parse_01(ofiles, individual=False):
    """ 
    a subfunction for summarizing results
    """

    ## parse results from outfiles
    cols = []
    dats = []
    for ofile in ofiles:

        ## parse file
        with open(ofile) as infile:
            dat  = infile.read()
        lastbits = dat.split(".mcmc.txt\n\n")[1:]
        results = lastbits[0].split("\n\n")[0].split()

        ## get shape from ...
        shape = (((len(results) - 3) / 4), 4)
        dat = np.array(results[3:]).reshape(shape)
        cols.append(dat[:, 3].astype(float))
        
    if not individual:
        ## get mean results across reps
        cols = np.array(cols)
        cols = cols.sum(axis=0) / len(ofiles) #10.
        dat[:, 3] = cols.astype(str)

        ## format as a DF
        df = pd.DataFrame(dat[:, 1:])
        df.columns = ["delim", "prior", "posterior"]
        nspecies = 1 + np.array([list(i) for i in dat[:, 1]], dtype=int).sum(axis=1)
        df["nspecies"] = nspecies
        return df
    
    else:
        ## get mean results across reps
        #return cols
        res = []
        for i in range(len(cols)):
            x = dat
            x[:, 3] = cols[i].astype(str)
            x = pd.DataFrame(x[:, 1:])
            x.columns = ['delim', 'prior', 'posterior']
            nspecies = 1 + np.array([list(i) for i in dat[:, 1]], dtype=int).sum(axis=1)
            x["nspecies"] = nspecies
            res.append(x)
        return res


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
