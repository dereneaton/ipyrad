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
import tempfile
import requests
import subprocess as sps

import numpy as np
import pandas as pd

from .utils import Params
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError
from ..analysis.locus_extracter import LocusExtracter

try:
    import toytree
except ImportError:
    pass
_MISSING_TOYTREE = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install toytree -c eaton-lab
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

    randomize_order (bool):
        if True then when maxloci is set this will randomly sample a 
        different set of N loci in each replicate, rather than sampling
        just the first N loci < maxloci. 

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

    # init object for params
    def __init__(
        self,
        name,
        data,
        workdir="analysis-bpp", 
        guidetree=None, 
        imap=None, 
        minmap=None,
        maxloci=100,
        minsnps=0,
        maxmissing=1.0,
        minlen=50,
        reps_resample_loci=False,
        load_existing_results=False,
        *args, 
        **kwargs):

        # results files
        self.files = Params()
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []

        # store args
        self.name = name
        self.data = os.path.realpath(os.path.expanduser(data))
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.guidetree = guidetree
        self.imap = imap
        self.minmap = minmap
        self.maxloci = maxloci
        self.minsnps = minsnps
        self.maxmissing = maxmissing
        self.minlen = minlen
        self.reps_resample_loci = reps_resample_loci
        self.load_existing_results = load_existing_results

        # parallelization
        self.ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }        

        # update kwargs 
        self.asyncs = []
        self._kwargs = {
            "binary": os.path.join(sys.prefix, "bin", "bpp"),
            "infer_sptree": 0,
            "infer_delimit": 0,
            "infer_delimit_args": (0, 2),
            "speciesmodelprior": 1,
            "seed": 12345,
            "burnin": 10000,
            "nsample": 100000,
            "sampfreq": 2,
            "thetaprior": (3, 0.002, "E"),
            "tauprior": (3, 0.002),
            "phiprior": (1, 1),
            "usedata": 1,
            "cleandata": 0,
            "finetune": (0.01, 0.02, 0.03, 0.04, 0.05, 0.01, 0.01),
            "copied": False,
        }
        self._check_kwargs(kwargs)
        self._kwargs.update(kwargs)

        # check that tmp binary is in /tmpdir
        self._check_binary()

        # run checks on args 
        self._check_args()

        # load existing results files for this named bpp object if they exist
        if self.load_existing_results:
            self._load_existing_results(self.name, workdir=self.workdir)



    def _check_binary(self):
        """
        Check for required software. If BPP is not present then a precompiled
        binary is downloaded into the tmpdir.
        """
        # check that toytree is installed
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # platform specific bpp binaries
        platform = "linux"
        if sys.platform != "linux":
            platform = "macos"
        dirname = "bpp-4.1.4-{}-x86_64".format(platform)

        # look for existing binary in tmpdir
        self._kwargs["binary"] = os.path.join(
            tempfile.gettempdir(), dirname, "bin", "bpp"
        )

        # check that bpp is installed and in path            
        cmd = ['which', self._kwargs["binary"]]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if comm:
            return 

        # bpp not found in /tmp, download it.        
        tarname = dirname + ".tar.gz"
        url = "https://github.com/bpp/bpp/releases/download/v4.1.4/" + tarname
        res = requests.get(url, allow_redirects=True)
        tmptar = os.path.join(tempfile.gettempdir(), tarname)
        with open(tmptar, 'wb') as tz:
            tz.write(res.content)

        # decompress tar file 
        cmd = ["tar", "zxvf", tmptar, "-C", tempfile.gettempdir()]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()
        if proc.returncode:
            print(comm[0], comm[1])

        # check that binary now can be found
        cmd = ['which', self._kwargs["binary"]]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if comm:
            return 
        raise IPyradError("bpp binary not found.")


    def _check_kwargs(self, kwargs):
        # support for legacy args
        for kwarg in kwargs:
            if kwarg not in self._kwargs:
                print(
                    "argument {} is either incorrect or no longer supported "
                    "please check the latest documentation".format(kwarg))


    def _check_args(self):
        """
        Check that data is a SEQS HDF5.
        """
        # check for data input
        if 'seqs.hdf5' not in self.data:
            raise IPyradError(
                "'data' argument must be an ipyrad .seqs.hdf5 file.")

        # set the guidetree
        if not self.guidetree:
            raise IPyradError(
                "must enter a 'guidetree' argument (a newick file or string).")
        self.tree = toytree.tree(self.guidetree)

        # check workdir
        if not self.workdir:
            self.workdir = os.path.join(os.path.curdir, "analysis-bpp")
        self.workdir = os.path.realpath(os.path.expanduser(self.workdir))           
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # parsing imap dictionary, or create simple 1-1 mapping
        if not self.imap:
            self.imap = {i: [i] for i in self.tree.get_tip_labels()}

        # check that all tree tip labels are represented in imap
        itips = set(self.imap.keys())
        ttips = set(self.tree.get_tip_labels())
        if itips.difference(ttips):
            raise IPyradError(
                "guidetree tips not in IMAP dictionary: {}"
                .format(itips.difference(ttips)))
        if ttips.difference(itips):
            raise IPyradError(
                "IMAP keys not in guidtree: {}"
                .format(ttips.difference(itips)))

        # check that minmap is OK or set it.
        # ...

        # checks
        assert isinstance(self.imap, dict), "you must enter an IMAP dictionary"
        assert set(self.imap.keys()) == set(self.tree.get_tip_labels()), (
            "IMAP keys must match guidetree names: \n{}\n{}"
            .format(self.imap.keys(), self.tree.get_tip_labels()))


    def _load_existing_results(self, name, workdir):
        """
        Load existing results files for an object with this workdir and name. 
        This does NOT reload the parameter settings for the object...
        """
        ## get mcmcs
        path = os.path.realpath(os.path.join(self.workdir, self.name))
        mcmcs = glob.glob("{}_r*.mcmc.txt".format(path))
        outs = glob.glob("{}_r*.out.txt".format(path))
        trees = glob.glob("{}_r*.tre".format(path))

        for mcmcfile in mcmcs:
            if mcmcfile not in self.files.mcmcfiles:
                self.files.mcmcfiles.append(mcmcfile)
        for outfile in outs:
            if outfile not in self.files.outfiles:
                self.files.outfiles.append(outfile)
        for tree in trees:
            if tree not in self.files.treefiles:
                self.files.treefiles.append(tree)        


    @property
    def _algorithm(self):
        if self._kwargs["infer_sptree"]:
            if self._kwargs["infer_delimit"]:
                return "11"
            return "10"
        else:
            if self._kwargs["infer_delimit"]:
                return "01"
            return "00"


    def run(self, ipyclient=None, force=False, show_cluster=False, auto=False, nreps=1, dry_run=False):
        """
        Submits bpp jobs to run on a cluster (ipyparallel Client). 
        The seed for the random number generator if not set is randomly 
        drawn, and if multiple reps are submitted (nreps>1) then each will 
        draw a subsequent random seeds after that. An ipyclient connection 
        is required. Asynchronous result objects are stored in the bpp 
        object submitting the jobs. 

        Parameters:
        -----------
        ipyclient: (type=ipyparallel.Client); Default=None. 
            If you started an ipyclient manually then you can 
            connect to it and use it to distribute jobs here.

        force: (type=bool); Default=False.
            Force overwrite of existing output with the same name.

        show_cluster: (type=bool); Default=False.
            Print information about parallel connection.

        auto: (type=bool); Default=False.
            Let ipyrad automatically manage ipcluster start and shutdown. 
            This will connect to all avaiable cores by default, but can 
            be modified by changing the parameters of the .ipcluster dict 
            associated with this tool.
        """
        pool = Parallel(
            tool=self,
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force, "nreps": nreps, "dry_run": dry_run},
            )
        pool.wrap_run()


    def _run(self, force, ipyclient, nreps, dry_run):
        "Distribute bpp jobs in parallel."

        # clear out pre-existing files for this object
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []
        self.asyncs = []

        # load-balancer
        lbview = ipyclient.load_balanced_view()

        # apply locus extracter filtering
        self.lex = LocusExtracter(
            data=self.data, 
            imap=self.imap,
            minmap=self.minmap,
            mincov=len(self.imap),  # ENFORCE at least 1 per spp.
            minsnps=self.minsnps,
            maxmissing=self.maxmissing,
            minlen=self.minlen,
        )
        self.lex.run(ipyclient=ipyclient, force=True, show_cluster=False)
        self.lex.wpnames = np.array(["^" + i for i in self.lex.wpnames])
        self.maxloci = min([self.maxloci, len(self.lex.loci)])

        # print BPP header
        print("[bpp v4.1.4]")  

        # initiate random seed 
        np.random.seed(self._kwargs["seed"])

        # static set of sampled loci if not resampling
        self._lidxs = np.random.choice(
                range(len(self.lex.loci)), 
                size=self.maxloci, 
                replace=False,
        )

        # track jobs
        # njobs = len(jobs)
        # printstr = "running {} structure jobs".format(njobs)
        # prog = ProgressBar(njobs, None, printstr)

        # replicate jobs can: ------
        #    - 1. differ by random seed AND by the subset of loci
        #    - 2. differ only by the random seed and use same subset of loci.
        for job in range(nreps):

            # make repname and make ctl filename
            self._name = "{}_r{}".format(self.name, job)
            ctlhandle = os.path.realpath(
                os.path.join(self.workdir, "{}.ctl.txt".format(self._name)))

            # skip if ctlfile exists
            if (not force) and (os.path.exists(ctlhandle)):
                print("Named ctl file already exists. Use force=True to" 
                      " overwrite\nFilename:{}".format(ctlhandle))

            # submit job to run
            else:
                # write imap groupings to the imapfile
                self._write_mapfile()

                # write the seq aligns to the .txt file
                self._write_seqfile()

                # change seed for each rep. CTL has other file paths.
                self._seed = np.random.randint(0, 1e9)
                ctlfile = self._write_ctlfile()

                # submit to engines
                if not dry_run:
                    rasync = lbview.apply(_call_bpp, *(self._kwargs["binary"], ctlfile))
                    self.asyncs.append(rasync)

        # report on the number of submitted jobs
        if self.asyncs:  # and (not quiet):
            print("distributed {} bpp jobs [{}] ({} loci)\n"
                  .format(nreps, self.name, self.maxloci))
        else:
            print("wrote {} bpp ctl files [{}] ({} loci)\n"
                  .format(nreps, self.name, self.maxloci))

        # block until jobs are done with a progress bar.
        # track progress of running job if fewer jobs than cores.
        # else track completed jobs.
        # ...TODO.
        


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
        # remove any old jobs with this same job name
        self._name = self.name
        oldjobs = glob.glob(
            os.path.join(self.workdir, self._name + "*.ctl.txt"))

        for job in oldjobs:
            os.remove(job)

        # check params types
        # ...
        # self._lidxs = np.random.choice(
            # range(len(self.lex.loci)), 
            # size=self.maxloci, 
            # replace=False,
        # )

        # write tmp files for the job
        self._write_seqfile(randomize_order=randomize_order)
        self._write_mapfile()  # name=True)
        self._write_ctlfile()

        # report to user
        if not quiet:
            sys.stderr.write(
                "input files created for job {} ({} loci)\n"
                .format(self._name, self._nloci)
            )


    def _write_mapfile(self):
        """
        Writes the IMAP formatted file for bpp from the ipa IMAP
        """
        # get outfile path
        self.mapfile = os.path.realpath(
            os.path.join(
                self.workdir, self._name + ".imapfile.txt"
            )
        )

        # get longest name in the file
        longname = 0
        for key in sorted(self.imap.keys()):
            for name in self.imap[key]:
                longname = max(longname, len(name))
        formatstr = "{:<" + str(longname + 2) + "} {}"

        # open handle for writing
        with open(self.mapfile, 'w') as mapfile:
            data = [
                formatstr.format(val, key) for key in 
                sorted(self.imap) for val in self.imap[key]
            ]
            mapfile.write("\n".join(data))


    def _write_seqfile(self):
        """
        Write filtered loci from locus extracter
        """
        # set path to output seq data file
        self.seqfile = os.path.realpath(
            os.path.join(
                self.workdir, self._name + ".seqfile.txt"
            )
        )

        # sample loci that meet the filtering requirements
        if self.reps_resample_loci:
            lidxs = np.random.choice(
                range(len(self.lex.loci)), 
                size=self.maxloci, 
                replace=False,
            )
        else:
            lidxs = self._lidxs

        # iterate over loci, printing to outfile
        with open(self.seqfile, 'w') as seqfile:
            for lidx in lidxs:
                header = "{} {}".format(*self.lex.get_shape(lidx))
                locus = "\n".join(self.lex.get_locus(lidx))
                seqfile.write("{}\n{}\n\n".format(header, locus))


    def _write_ctlfile(self):
        """ write outfile with any args in argdict """

        # get full path to out files for this repname
        path = os.path.realpath(os.path.join(self.workdir, self._name))
        mcmcfile = "{}.mcmc.txt".format(path)
        outfile = "{}.out.txt".format(path)

        # store files for this rep
        if mcmcfile not in self.files.mcmcfiles:
            self.files.mcmcfiles.append(mcmcfile)
        if outfile not in self.files.outfiles:
            self.files.outfiles.append(outfile)

        # expand options to fill ctl file
        ctlstring = CTLFILE.format(**{
            "seqfile": self.seqfile,
            "mapfile": self.mapfile,
            "mcmcfile": mcmcfile,
            "outfile": outfile,

            "nloci": self.maxloci,
            "usedata": self._kwargs["usedata"],
            "cleandata": self._kwargs["cleandata"],

            "infer_sptree": int(self._kwargs["infer_sptree"]),
            "infer_delimit": int(self._kwargs["infer_delimit"]),
            "infer_delimit_args": (
                " ".join([str(i) for i in self._kwargs["infer_delimit_args"]])
                if self._kwargs["infer_delimit"]
                else ""),
            "nsp": len(self.imap),
            "spnames": " ".join(sorted(self.imap)),
            "spcounts": " ".join([str(len(self.imap[i])) for i in sorted(self.imap)]),
            "spnewick": self.tree.write(tree_format=9),
            "speciesmodelprior": self._kwargs["speciesmodelprior"],

            "thetaprior": " ".join([str(i) for i in self._kwargs["thetaprior"]]),
            "tauprior": " ".join([str(i) for i in self._kwargs["tauprior"]]),
            "phiprior": " ".join([str(i) for i in self._kwargs["phiprior"]]),

            "seed": self._seed,
            "finetune": " ".join([str(i) for i in self._kwargs["finetune"]]),
            "burnin": self._kwargs["burnin"],
            "sampfreq": self._kwargs["sampfreq"],
            "nsample": self._kwargs["nsample"],
        })

        # write out the ctl file
        ctlhandle = os.path.realpath(
            "{}.ctl.txt".format(os.path.join(self.workdir, self._name)))
        with open(ctlhandle, 'w') as out:
            out.write(ctlstring)

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

        # make deepcopy of self.__dict__ but do not copy async objects
        subdict = {i: j for i, j in self.__dict__.items() if i != "asyncs"}
        newdict = copy.deepcopy(subdict)

        # make back into a bpp object
        if name == self.name:
            raise Exception(
                "new object must have a different 'name' than its parent")

        newobj = Bpp(
            name=name,
            data=newdict["files"].data,
            workdir=newdict["workdir"],
            guidetree=newdict["tree"].write(),
            imap={i: j for i, j in newdict["imap"].items()},
            copied=True,
            load_existing_results=load_existing_results,
            )

        # update special dict attributes but not files
        for key, val in newobj.params.__dict__.items():
            newobj.params.__setattr__(key, self._kwargs.__getattribute__(key))
        for key, val in newobj.filters.__dict__.items():
            newobj.filters.__setattr__(key, self.filters.__getattribute__(key))

        # new object must have a different name than it's parent
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
                    [pd.read_csv(i, sep='\t', index_col=0) 
                    for i in self.files.mcmcfiles]
                    ).describe().T

        ## algorithm 01
        if self.params.infer_delimit & (not self.params.infer_sptree):
            return _parse_01(self.files.outfiles, individual=individual_results)

        ## others
        else:
            print("summary function not yet ready for this type of result")
            return 0



def _call_bpp(binary, ctlfile):  #, is_alg00):

    # call the command and block until job finishes
    cmd = [binary, "--cfile", ctlfile]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    comm = proc.communicate()
    if proc.returncode:
        raise IPyradError(comm[0])

    # Look for the ~/Figtree.tre file that bpp creates.
    # This has to be done here to make sure it is instantly run
    # when the job finishes so other reps won't write over it.
    # Kludge due to bpp harcoded to write FigTree file to CWD.
    # if is_alg00:
    #     default_figtree_path = os.path.join(
    #         os.path.expanduser("~"), "FigTree.tre")
    #     new_figtree_path = ctlfile.rsplit(".ctl.txt", 1)[0] + ".tre"
    #     try:
    #         if os.path.exists(default_figtree_path):
    #             os.rename(default_figtree_path, new_figtree_path)
    #     except Exception:
    #         pass
    # return comm[0]


class Result_00(object):
    """ parse bpp results object for 00 scenario """
    pass


def _parse_00(ofile):
    """
    return 00 outfile as a pandas DataFrame
    """
    with open(ofile) as infile:
        # read in the results summary from the end of the outfile
        arr = np.array(
            [" "] + infile.read().split("Summary of MCMC results\n\n\n")[1:][0]\
            .strip().split())

        # reshape array 
        rows = 12
        cols = (arr.shape[0] + 1) / rows
        arr = arr.reshape(rows, cols)

        # make into labeled data frame
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

    # parse results from outfiles
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


# GLOBALS
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





CTLFILE = """
* I/O 
seqfile = {seqfile}
Imapfile = {mapfile}
mcmcfile = {mcmcfile}
outfile = {outfile}

* DATA
nloci = {nloci}
usedata = {usedata}
cleandata = {cleandata}

* MODEL
speciestree = {infer_sptree}
speciesdelimitation = {infer_delimit} {infer_delimit_args}
speciesmodelprior = {speciesmodelprior}
species&tree = {nsp} {spnames}
               {spcounts}
               {spnewick}

* PRIORS
thetaprior = {thetaprior}
tauprior = {tauprior}
phiprior = {phiprior}

* MCMC PARAMS
seed = {seed}
finetune = 1: {finetune}
print = 1 0 0 0
burnin = {burnin}
sampfreq = {sampfreq}
nsample = {nsample}
"""







