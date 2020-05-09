#!/usr/bin/env python

"convert loci file to bpp format input files"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard lib
import os
import sys
import glob
import time
import copy
import tempfile
import requests
import itertools
import subprocess as sps

import numpy as np
import pandas as pd
import scipy.stats as ss

from .utils import Params, ProgressBar
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError
from ..analysis.locus_extracter import LocusExtracter

try:
    import toytree
    # from toytree.utils import bpp2newick
except ImportError:
    pass
_MISSING_TOYTREE = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install toytree -c eaton-lab
"""

# TODO: raise error on locus_filter result 0


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

    reps_sample_loci (bool):
        if True then when maxloci is set this will randomly sample a 
        different set of N loci in each replicate, rather than sampling
        just the first N loci < maxloci. 

    infer_sptree:
        Default=0, only infer parameters on a fixed species tree. If 1, then the
        input tree is treated as a guidetree and tree search is employed to find
        the best tree. The results will include support values for the inferred
        topology.

    infer_delimit:
        Default=0, no delimitation. If 1 then splits in the tree that separate
        'species' will be collapsed to test whether fewer species are a better
        fit to the data than the number in the input guidetree.

    infer_delimit_args:
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

    maxloci (int):
        The max number of loci that will be used in an analysis.
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
        data=None,
        workdir="analysis-bpp", 
        guidetree=None, 
        imap=None, 
        minmap=None,
        maxloci=100,
        minsnps=0,
        maxmissing=1.0,
        minlen=50,
        reps_resample_loci=False,
        # load_existing_results=False,
        *args, 
        **kwargs):

        # results files
        self.files = Params()
        self.files.mcmcfiles = []
        self.files.outfiles = []
        self.files.treefiles = []

        # store args
        self.data = data
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.guidetree = guidetree
        self.imap = imap
        self.minmap = minmap
        self.maxloci = maxloci
        self.minsnps = minsnps
        self.maxmissing = maxmissing
        self.minlen = minlen
        self.reps_resample_loci = reps_resample_loci
        # self.load_existing_results = load_existing_results

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
        self.kwargs = {
            "binary": None,
            "infer_sptree": 0,
            "infer_delimit": 0,
            "infer_delimit_args": (0, 2),
            "speciesmodelprior": 1,
            "seed": 12345,
            "burnin": 10000,
            "nsample": 100000,
            "sampfreq": 2,
            "thetaprior": (3, 0.002),
            "tauprior": (3, 0.002),
            "phiprior": (1, 1),
            "usedata": 1,
            "cleandata": 0,
            "finetune": (0.01, 0.02, 0.03, 0.04, 0.05, 0.01, 0.01),
            "copied": False,
        }

        # binary is needed for running or loading and combining results
        self._check_binary()

        # update and check kwargs if data else do nothing which allows 
        # loading dummy objects for summarizing existing results.
        if data:
            self._check_kwargs(kwargs)
            self.kwargs.update(kwargs)
            self._check_args()


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
        self.kwargs["binary"] = os.path.join(
            tempfile.gettempdir(), dirname, "bin", "bpp"
        )

        # check that bpp is installed and in path            
        cmd = ['which', self.kwargs["binary"]]
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
        cmd = ['which', self.kwargs["binary"]]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if comm:
            return 
        raise IPyradError("bpp binary not found.")


    def _check_kwargs(self, kwargs):
        # support for legacy args
        for kwarg in kwargs:
            if kwarg not in self.kwargs:
                print(
                    "argument {} is either incorrect or no longer supported "
                    "please check the latest documentation".format(kwarg))

            # check type
            if kwarg in ['nloci', 'burnin', 'sampfreq', 'nsample']:
                kwargs[kwarg] = int(kwargs[kwarg])


    def _check_args(self):
        """
        Check that data is a SEQS HDF5.
        """
        # expand path to data
        self.data = os.path.realpath(os.path.expanduser(self.data))

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
        # get mcmcs
        path = os.path.realpath(os.path.join(self.workdir, self.name))
        mcmcs = sorted(glob.glob("{}_r*.mcmc.txt".format(path)))
        outs = sorted(glob.glob("{}_r*.out.txt".format(path)))
        trees = sorted(glob.glob("{}_r*.tre".format(path)))

        for mcmcfile in mcmcs:
            if mcmcfile not in self.files.mcmcfiles:
                self.files.mcmcfiles.append(mcmcfile)
        for outfile in outs:
            if outfile not in self.files.outfiles:
                self.files.outfiles.append(outfile)
        for tree in trees:
            if tree not in self.files.treefiles:
                self.files.treefiles.append(tree)        
        print("[ipa.bpp] found {} existing result files".format(len(mcmcs)))


    @property
    def _algorithm(self):
        if self.kwargs["infer_sptree"]:
            if self.kwargs["infer_delimit"]:
                return "11"
            return "10"
        else:
            if self.kwargs["infer_delimit"]:
                return "01"
            return "00"


    def run(self, ipyclient=None, force=False, show_cluster=True, auto=False, nreps=1, dry_run=False):
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

        # raise error if no loci passed filtering
        if not self.maxloci:
            raise IPyradError("No loci passed filtering. Check params.")

        # print BPP header
        print("[ipa bpp] bpp v4.1.4")

        # initiate random seed 
        np.random.seed(self.kwargs["seed"])

        # static set of sampled loci if not resampling
        self._lidxs = np.random.choice(
            range(len(self.lex.loci)), 
            size=self.maxloci, 
            replace=False,
        )

        # replicate jobs
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
                # get random subset if reps_resample_loci=True
                self._write_seqfile()

                # change seed for each rep. CTL has other file paths.
                self._seed = np.random.randint(0, 1e9)
                ctlfile = self._write_ctlfile()

                # submit to engines
                if not dry_run:
                    args = (self.kwargs["binary"], ctlfile, self._algorithm)
                    rasync = lbview.apply(_call_bpp, *args)
                    self.asyncs.append(rasync)

        # report on the files written
        if not self.asyncs:  # and (not quiet):
            print("[ipa.bpp] wrote {} bpp ctl files (name={}, nloci={})"
                  .format(nreps, self.name, self.maxloci))

        # report on the number of submitted jobs
        else:
            print("[ipa.bpp] distributed {} bpp jobs (name={}, nloci={})"
                  .format(nreps, self.name, self.maxloci))

            # setup progress bar
            rep = 0
            nits = int(self.kwargs["nsample"])
            prog = ProgressBar(nits, None, "progress on rep {}".format(rep))
            prog.finished = 0
            prog.update()

            # block until jobs are done with a progress bar.
            spacer = 0
            while 1:

                # get file to check for results
                checkresult = self.files.mcmcfiles[rep]

                # check for job progress periodically
                if spacer == 5:
                    if os.path.exists(checkresult):
                        with open(checkresult, 'rb') as infile:
                            nlines = 0
                            for line in infile:
                                nlines += 1
                        prog.finished = nlines
                    spacer = 0

                # break between checking progress       
                prog.update()
                time.sleep(5)
                spacer += 1

                # finished jobs
                finished = [i.ready() for i in self.asyncs]

                # check if a different rep should be tracked.
                if not all(finished):
                    if finished[rep]:
                        rep = [i for (i, j) in enumerate(finished) if not j][0]
                    prog.message = "progress on rep {}".format(rep)

                # all jobs finished
                else:
                    prog.message = "progress on all reps"
                    prog.finished = nits
                    prog.update()
                    print("")
                    break

            # check for job failures and raise an error:
            for job in self.asyncs:
                if not job.successful():
                    print(job.result())


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
            "usedata": self.kwargs["usedata"],
            "cleandata": self.kwargs["cleandata"],

            "infer_sptree": int(self.kwargs["infer_sptree"]),
            "infer_delimit": int(self.kwargs["infer_delimit"]),
            "infer_delimit_args": (
                " ".join([str(i) for i in self.kwargs["infer_delimit_args"]])
                if self.kwargs["infer_delimit"]
                else ""),
            "nsp": len(self.imap),
            "spnames": " ".join(sorted(self.imap)),
            "spcounts": " ".join([str(len(self.imap[i])) for i in sorted(self.imap)]),
            "spnewick": self.tree.write(tree_format=9),
            "speciesmodelprior": self.kwargs["speciesmodelprior"],

            "thetaprior": " ".join([str(i) for i in self.kwargs["thetaprior"]]),
            "tauprior": " ".join([str(i) for i in self.kwargs["tauprior"]]),
            "phiprior": " ".join([str(i) for i in self.kwargs["phiprior"]]),
            "estimate_theta": ("E" if self._algorithm == "00" else ""),

            "seed": self._seed,
            "finetune": " ".join([str(i) for i in self.kwargs["finetune"]]),
            "burnin": self.kwargs["burnin"],
            "sampfreq": self.kwargs["sampfreq"],
            "nsample": self.kwargs["nsample"],
        })

        # write out the ctl file
        ctlhandle = os.path.realpath(
            "{}.ctl.txt".format(os.path.join(self.workdir, self._name)))
        with open(ctlhandle, 'w') as out:
            out.write(ctlstring)

        return ctlhandle


    def copy(self, name):
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
        if name == self.name:
            raise Exception(
                "new object must have a different 'name' than its parent")
        newself = copy.deepcopy(self)
        newself.name = name
        newself.files.mcmcfiles = []
        newself.files.outfiles = []
        newself.asyncs = []
        return newself


    def summarize_results(self, algorithm, individual_results=False):
        """ 
        Prints a summarized table of results from replicate runs, or,
        if individual_result=True, then returns a list of separate
        dataframes for each replicate run. 
        """
        # reports number of results found
        self._load_existing_results(self.name, self.workdir)
        assert algorithm in ["00", "01", "10", "11"]
        print("[ipa.bpp] summarizing algorithm '{}' results".format(algorithm))

        # algorithms supported
        if algorithm == "00":
            return self._summarize_00(individual_results)
        if algorithm == "10":
            return self._summarize_10(individual_results)
        if algorithm == "01":
            return self._summarize_01(individual_results)



    def _summarize_00(self, individual_results):
        """
        Combines MCMC files together and writes a ctl file then calls
        bpp with the print=-1 option which means "read in" so that it 
        will compute a new posterior table...
        """
        # load out tables of summarized posteriors
        tables = []
        for ofile in self.files.outfiles:
            with open(ofile, 'r') as infile:
                lines = infile.readlines()[-12:]
                data = [i.strip().split() for i in lines]
                index = [i[0] for i in data[1:]]
                df = pd.DataFrame(
                    data=[i[1:] for i in data[1:]],
                    columns=data[0],
                    index=index,
                )
                tables.append(df.astype(float))

        # load mcmc tables of posteriors
        dfs = [
            pd.read_csv(i, sep='\t', index_col=0) 
            for i in self.files.mcmcfiles
        ]

        # return a list of parsed CSV results
        if individual_results:
            return tables, dfs

        # concatenate each CSV and then get stats w/ describe
        else:
            print('[ipa.bpp] combining mcmc files')

            # new file handles
            cf = self.files.mcmcfiles[0].rsplit("_", 1)[0] + "_concat.mcmc.txt"
            of = self.files.mcmcfiles[0].rsplit("_", 1)[0] + "_concat.out.txt"            

            # existing and new ctl files
            ctlfile = os.path.join(self.workdir, self.name + "_r0.ctl.txt")
            newctl = os.path.join(self.workdir, self.name + "_tmp.ctl.txt")

            # write a concatenated mcmc file
            concat = pd.concat(dfs, ignore_index=True)
            concat.to_csv(cf, sep="\t", float_format="%.6f")

            # write a tmp ctl file with print=-1 and mcmcfile=cf
            with open(ctlfile, 'r') as infile:
                with open(newctl, 'w') as outfile:
                    cdat = infile.readlines()
                    for line in cdat:
                        if 'mcmcfile' in line:
                            line = "mcmcfile = {}\n".format(cf)
                        if 'outfile' in line:
                            line = "outfile = {}\n".format(of)
                        if 'print' in line:
                            line = "print = -1\n"
                        outfile.write(line)

            # run bpp on the new ctlfile
            _call_bpp(self.kwargs["binary"], newctl, "00")

            # cleanup tmp file
            os.remove(newctl)

            # load the new table
            with open(ofile, 'r') as infile:
                lines = infile.readlines()[-12:]
                data = [i.strip().split() for i in lines]
                index = [i[0] for i in data[1:]]
                table = pd.DataFrame(
                    data=[i[1:] for i in data[1:]],
                    columns=data[0],
                    index=index,
                )
            return table, concat



    def _summarize_01(self, individual_results):
        dfs = []
        for ofile in self.files.outfiles:
            with open(ofile, 'r') as infile:
                dat = infile.read().split("posterior\n")[1]
                table, dat = dat.split("Order of ancestral nodes:")
                data = [i.strip().split() for i in table.strip().split("\n")]
                df = pd.DataFrame(
                    data=data,
                    columns=["x", "delim", "prior", "posterior"]
                )
                df = df.drop(columns=['x'])
                df["nspecies"] = [i.count("1") + 1 for i in df["delim"]]
                df["posterior"] = df["posterior"].astype(float)
                dfs.append(df)

        if individual_results:
            return dfs
        else:
            df = dfs[0]
            for odf in dfs[1:]:
                df["posterior"] += odf.posterior
            df["posterior"] /= len(dfs)
            return df



    def _summarize_10(self, individual_results):
        """
        Returns a tuple with (tree, treedist) where tree is a Toytree 
        containing the MJrule tree with support values on the edges and 
        treedist is a multitree containing the posterior distribution of trees. 
        """
        # store results 
        trees = []
        treelists = []

        # get best trees
        for treefile in self.files.outfiles:
            with open(treefile, 'r') as infile:

                # jump to end of file to get besttree
                for line in infile:
                    pass
                newick = line.split(";")[0] + ";"

                # get majority-rule tree
                # while 1:
                #     line = next(infile)
                #     if line.startswith("(C) Majority"):
                #         newick = next(infile).strip()
                #         break

                # extract proper newick from bpp last line
                newick = newick.replace(" #", "")
                tree = toytree.tree(newick)

                # convert support values to ints
                for node in tree.treenode.traverse():
                    node.support = int(round(node.support * 100))
                trees.append(tree)

        # get posteriors
        for treefile in self.files.mcmcfiles:
            post = []
            with open(treefile, 'r') as infile:
                btrees = infile.readlines()
                for newick in btrees:
                    # newick = bpp2newick(tre)
                    post.append(newick)
                treelists.append(post)

        # return results
        if individual_results:
            return trees, [toytree.mtree(i) for i in treelists]
        else:
            return trees, toytree.mtree(list(itertools.chain(*treelists)))



    def draw_priors(self, gentime_min, gentime_max, mutrate_min, mutrate_max, invgamma=True, seed=123):
        """

        """
        import toyplot

        # setup canvas
        canvas0 = toyplot.Canvas(width=925, height=300)
        ax0 = canvas0.cartesian(
            bounds=(50, 275, 50, 250), 
            xlabel="prior on mutation rates (x10^-8)",
        )
        ax1 = canvas0.cartesian(
            bounds=(350, 575, 50, 250), 
            xlabel="prior on theta (4Neu)",
        )
        ax2 = canvas0.cartesian(
            bounds=(650, 875, 50, 250),
            xlabel="prior on Ne",
        )
        for ax in (ax0, ax1, ax2):
            ax.y.ticks.labels.show = False

        # distribution of mutation_rates ---------------------------------
        mean = (mutrate_max + mutrate_min) / 2.
        var = ((mutrate_max - mutrate_min) ** 2) / 16
        a = mean ** 2 / var
        b = mean / var
        muts_rvs = ss.gamma.rvs(
            a, **{"scale": 1 / b, 'random_state': 123, "size": 1000})

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax0.fill(x * 1e8, y, opacity=0.33, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax0.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0] * 1e8, 1), round(x[-1] * 1e8, 1)))


        # distribution of prior on theta ---------------------------------
        # invgamma_a = self.kwargs["thetaprior"][0]
        # invgamma_b = self.kwargs["thetaprior"][1]
        # mean = invgamma_b / (invgamma_a - 1.)
        # var = (invgamma_b ** 2) / (((invgamma_a - 1) ** 2) * (invgamma_a - 2))
        # a = mean ** 2 / var
        # b = mean / var
        a = self.kwargs["thetaprior"][0]
        b = self.kwargs["thetaprior"][1]
        if invgamma:
            b = 1 / b

        theta_rvs = ss.gamma.rvs(
            a, **{"scale": 1 / b, 'random_state': 123, "size": 1000})

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax1.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax1.label.text = (
                    "95% CI: {:.3f} - {:.3f}"
                    .format(x[0], x[-1]))


        # distribution of effective population sizes --------------------
        ne_rvs = (theta_rvs / (muts_rvs * 4))
        mean, var, std = ss.bayes_mvs(ne_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax2.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax2.label.text = (
                    "95% CI: {:.0f} - {:.0f}"
                    .format(x[0], x[-1]))


        # plot  ---------------------------------------------------------
        canvas1 = toyplot.Canvas(width=925, height=300)
        ax0 = canvas1.cartesian(
            bounds=(50, 275, 50, 250), 
            xlabel="prior on mutation rates (x10^-8)",
        )
        ax1 = canvas1.cartesian(
            bounds=(350, 575, 50, 250), 
            xlabel="prior on theta (4Neu)",
        )
        ax2 = canvas1.cartesian(
            bounds=(650, 875, 50, 250),
            xlabel="prior on Ne",
        )
        for ax in (ax0, ax1, ax2):
            ax.y.ticks.labels.show = False

        # distribution of generation times -------------------------------
        mean = (gentime_max + gentime_min) / 2.
        var = ((gentime_max - gentime_min) ** 2) / 16
        a = mean ** 2 / var
        b = mean / var
        gens_rvs = ss.gamma.rvs(
            a, **{"scale": 1 / b, 'random_state': 123, "size": 1000})

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax0.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[0])

            if cix == 0.95:
                ax0.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(round(x[0], 1), round(x[-1], 1)))


        # distribution of prior on tau ----------------------------------
        a = self.kwargs["tauprior"][0]
        b = self.kwargs["tauprior"][1]        
        if invgamma:
            b = 1 / b
        tau_rvs = ss.gamma.rvs(
            a, **{"scale": 1 / b, 'random_state': 123, "size": 1000})

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax1.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[1])

            if cix == 0.95:
                ax1.label.text = (
                    "95% CI: {:.4f} - {:.4f}"
                    .format(x[0], x[-1]))


        # distribution of divergence times in years ----------------------
        div_rvs = (gens_rvs * tau_rvs) / muts_rvs
        div_rvs /= 1e6
        mean, var, std = ss.bayes_mvs(div_rvs)
        a = mean.statistic ** 2 / var.statistic
        b = mean.statistic / var.statistic

        # draw dist
        for cix in (0.99, 0.95, 0.5):
            edge = (1 - cix) / 2.
            x = np.linspace(
                ss.gamma.ppf(edge, a, **{"scale": 1 / b}),
                ss.gamma.ppf(1 - edge, a, **{"scale": 1 / b}),
                100)
            y = ss.gamma.pdf(x, a, **{"scale": 1 / b})
            ax2.fill(x, y, opacity=0.25, color=toyplot.color.Palette()[2])

            if cix == 0.95:
                ax2.label.text = (
                    "95% CI: {:.1f} - {:.1f}"
                    .format(x[0], x[-1]))

        return canvas0, canvas1



    def draw_posteriors(self, ):
        pass



class Transformer(object):
    """
    ...
    """
    def __init__(self, df, gentime_min, gentime_max, mutrate_min, mutrate_max, seed=123):

        self.df = df
        self.seed = seed

        self.gentime_min = gentime_min
        self.gentime_max = gentime_max
        self.mutrate_min = mutrate_min
        self.mutrate_max = mutrate_max

        self.gentime_mean = (gentime_max + gentime_min) / 2.
        self.gentime_var = ((gentime_max - gentime_min) ** 2) / 16
        self.gentime_a = self.gentime_mean ** 2 / self.gentime_var
        self.gentime_b = self.gentime_mean / self.gentime_var

        self.mutrate_mean = (mutrate_max + mutrate_min) / 2.
        self.mutrate_var = ((mutrate_max - mutrate_min) ** 2) / 16
        self.mutrate_a = self.mutrate_mean ** 2 / self.mutrate_var
        self.mutrate_b = self.mutrate_mean / self.mutrate_var

        self._sample_gentime_rvs()
        self._sample_mutrate_rvs()


    def _sample_gentime_rvs(self):
        self.gentime_rvs = ss.gamma.rvs(
            self.gentime_a, 
            **{
                "scale": 1 / self.gentime_b, 
                "random_state": self.seed, 
                "size": self.df.shape[0],
            },
        )


    def _sample_mutrate_rvs(self):
        self.mutrate_rvs = ss.gamma.rvs(
            self.mutrate_a, 
            **{
                "scale": 1 / self.mutrate_b, 
                "random_state": self.seed, 
                "size": self.df.shape[0],
            },
        )


    def _get_gentime_x(self, nvalues=100):
        xvals = np.linspace(
            ss.gamma.ppf(0.0001, self.gentime_a, **{"scale": 1 / self.gentime_b}),
            ss.gamma.ppf(0.9999, self.gentime_a, **{"scale": 1 / self.gentime_b}),
            nvalues,
        )
        return xvals


    def _get_mutrate_x(self, nvalues=100):
        xvals = np.linspace(
            ss.gamma.ppf(0.0001, self.mutrate_a, **{"scale": 1 / self.mutrate_b}),
            ss.gamma.ppf(0.9999, self.mutrate_a, **{"scale": 1 / self.mutrate_b}),
            nvalues,
        )
        return xvals


    def _sample_tau(self, colname):
        # check that it is a tau column
        if "tau" not in colname:
            raise IPyradError("not a tau value: {}".format(colname))

        # get mean, var, std
        mvs = ss.bayes_mvs(self.df[colname])
        a = mvs[0].statistic ** 2 / mvs[1].statistic
        b = mvs[0].statistic / mvs[1].statistic

        # sampled taus
        self.tau_rvs = ss.gamma.rvs(
            a, 
            **{
                'scale': 1 / b, 
                "random_state": self.seed, 
                "size": self.df.shape[0],
            },
        )


    def _sample_theta(self, colname):
        # check that it is a theta column
        if "theta" not in colname:
            raise IPyradError("not a theta value: {}".format(colname))

        # get mean, var, std
        mvs = ss.bayes_mvs(self.df[colname])
        a = mvs[0].statistic ** 2 / mvs[1].statistic
        b = mvs[0].statistic / mvs[1].statistic

        # sampled taus
        self.theta_rvs = ss.gamma.rvs(
            a, 
            **{
                'scale': 1 / b, 
                "random_state": self.seed, 
                "size": self.df.shape[0],
            },
        )


    def _transform_tau(self, colname):

        # sampled taus for this column parameter
        self._sample_tau(colname)

        # get div time as (tau * gentime) / mutrate
        # i.e., (gens * (time/gen)) / (muts/site/gen)
        self.div_rvs = (self.tau_rvs * self.gentime_rvs) / self.mutrate_rvs


    def _transform_theta(self, colname):

        # sampled taus for this column parameter
        self._sample_theta(colname)

        # get Ne as (theta / 4*u)
        self.ne_rvs = self.theta_rvs / (self.mutrate_rvs * 4)


    def transform(self, colname):

        if "tau" in colname:
            self._transform_tau(colname)
            mean, var, std = ss.bayes_mvs(self.div_rvs)
            print("mean: {}".format(mean))
            print("95% CI: {}-{}".format(mean.minmax[0], mean.minmax[1]))
            draw_dist(
                mean.statistic, var.statistic, "Divergence time")

        if "theta" in colname:
            self._transform_theta(colname)
            mean, var, std = ss.bayes_mvs(self.ne_rvs)
            print("mean: {}".format(mean))
            print("95% CI: {}-{}".format(mean.minmax[0], mean.minmax[1]))
            canvas, axes = draw_dist(
                mean.statistic, var.statistic, "Effective population size")
        return canvas, axes




def _call_bpp(binary, ctlfile, alg):
    """
    Remote function call of BPP binary
    """
    # call the command and block until job finishes
    cmd = [binary, "--cfile", ctlfile]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    comm = proc.communicate()
    if proc.returncode:
        raise IPyradError(comm[0])

    # bpp writes a Figtree.tre result to CWD of the ipyparallel engine (ugh.)
    # so we need to catch it quickly and move it somewhere relevant. 
    if alg == "00":
        figfile = "Figtree.tre"
        newfigpath = ctlfile.replace(".ctl.txt", ".nex")
        if os.path.exists(figfile):
            os.rename(figfile, newfigpath)

    if os.path.exists("./SeedUsed"):
        os.remove("./SeedUsed")



def draw_dist(mean, var, xlabel):
    """
    Tranformer class subfunction to draw posterior density.
    """
    import toyplot
    a = mean ** 2 / var
    b = mean / var

    # set up plots
    canvas = toyplot.Canvas(width=400, height=300)
    axes = canvas.cartesian(ylabel="density", xlabel=xlabel)

    # confidence intervals shaded
    for civ in [0.99, 0.95, 0.50]:

        # get 100 values evenly spaced across 99% 
        edge = (1 - civ) / 2.
        x = np.linspace(
            ss.gamma.ppf(edge, a, **{'scale': 1 / b}),
            ss.gamma.ppf(1 - edge, a, **{'scale': 1 / b}), 
            100)

        # plot values across range of gamma
        axes.fill(
            x,  # / 1e6,
            ss.gamma.pdf(x, a, **{'scale': 1 / b}),
            opacity=0.25,
            color=toyplot.color.Palette()[0],
        )

    axes.y.ticks.labels.show = False
    axes.x.ticks.show = True
    return canvas, axes



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
thetaprior = {thetaprior} {estimate_theta}
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
