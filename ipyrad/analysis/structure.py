#!/usr/bin/env python

"convenience wrappers for running structure in a jupyter notebook"

# py2/3 compat
from __future__ import print_function
from builtins import range

# standard lib
import os
import re
import sys
import glob
import time
import subprocess as sps

# third party
import numpy as np
import pandas as pd

# ipyrad utils
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError
from .utils import Params, ProgressBar
from .snps_extracter import SNPsExtracter

# tODO: add subsample_snps = False as an option.


MISSING_IMPORTS = """
To use the ipa.structure module you must install two additional 
libraries which can be done with the following conda command. 

conda install structure clumpp -c ipyrad
"""


class Structure(object):
    """ 
    Create and return an ipyrad.analysis Structure Object. This object allows
    you to easily enter parameter setting to submit structure jobs to run in 
    parallel on a cluster. 


    Parameters
    -----------
    name (str):
        A prefix name for all output files. 

    data (str):
        A .snps.hdf5 file from ipyrad. If you do not have this file see the 
        ipyrad.analysis docs to create this file from a VCF output.

    workdir (str):
        Directory for output files; will be created if not present.

    ... (common ipyrad-analysis params supported)

    Attributes:
    ----------
    mainparams (dict):
        A dictionary with the mainparams used by STRUCTURE
    extraparams (dict):
        A dictionary with the extraparams used by STRUCTURE
    clumppparams (dict):
        A ditionary with the parameter settings used by CLUMPP
    header (pandas.DataFrame):
        Returns the header columns of the str file
    result_files (list):
        Returns a list of result files for finished STRUCTURE jobs submitted 
        by this object. 
    asyncs: (list):
        A list of asynchronous result objects for each job that was 
        submitted to the ipyclient. These can be used for debugging if 
        a job fails.


    Functions:
    ----------
    run(*args, **kwargs):
        Submits independent replicate jobs to run on a cluster.
    get_clumpp_table(kpop):
        Returns a table of results for K=kpop permuted across all replicates.

    """
    def __init__(
        self, 
        name, 
        data, 
        workdir="./analysis-structure",
        imap=None,
        minmap=None,
        mincov=0.0,
        minmaf=0.0,
        quiet=False,
        load_only=False,
        subsample_snps=True,
        ):

        # printing strategy
        self.quiet = quiet

        # get path to saved files and load any existing files
        self.name = name
        self.workdir = os.path.abspath(os.path.expanduser(workdir))

        # check attribute for existing results at this name.
        if self.result_files:
            self._print(
                "{} previous results loaded for run [{}]"
                .format(len(self.result_files), self.name))

        # the snps database file contains data and names, etc.
        self.data = os.path.abspath(os.path.expanduser(data))

        # filtering parameters
        self.imap = imap
        self.minmap = minmap
        self.mincov = mincov
        self.minmaf = minmaf
        self.subsample_snps = subsample_snps

        # run checks
        self.STRUCTURE = os.path.join(sys.prefix, "bin", "structure")
        self.CLUMPP = os.path.join(sys.prefix, "bin", "CLUMPP")
        self._check_binaries()
        self._setup_dirs()

        # load the database file for filtering/extracting later
        self._ext = SNPsExtracter(
            data=self.data,
            imap=self.imap, 
            minmap=self.minmap, 
            mincov=self.mincov,
            minmaf=self.minmaf,
        )       

        # can skip parsing the file if load=True
        self._load_only = load_only
        if not self._load_only:
            self._ext.parse_genos_from_hdf5()

        # header columns from imap, user can modify after init.
        self._get_header()

        # params
        self.mainparams = _MainParams()
        self.extraparams = _ExtraParams()
        self.clumppparams = _ClumppParams()
        self.rasyncs = {}

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


    def _check_binaries(self):
        "check for structure and clumpp"
        for binary in [self.STRUCTURE, self.CLUMPP]:
            if not os.path.exists(binary):
                raise IPyradError(MISSING_IMPORTS)


    def _setup_dirs(self):
        # make workdir if it does not exist
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)


    def _print(self, value):
        if not self.quiet:
            print(value)


    def _get_header(self):
        """
        Build header columns of the STRUCTURE input file.
        """
        # names are alphanumeric (x2) except for "reference" which is at top
        if "reference" not in self._ext.names:
            labels = sorted(self._ext.names * 2)
        else:
            labels = sorted(self._ext.names[1:] * 2)
            labels = ["reference", "reference"] + labels

        # make reverse imap dict for extracting popdata
        if self.imap:
            rdict = {}
            for key, val in self.imap.items():
                if isinstance(val, (str, int)):
                    rdict[val] = key
                elif isinstance(val, (list, tuple)):
                    for tax in val:
                        rdict[tax] = key

        popdata = [""] * len(labels)  # [rdict[i] for i in self.labels]
        popflag = [""] * len(labels)  # ["0"] * len(self.labels)
        locdata = [""] * len(labels)
        phenotype = [""] * len(labels)        

        self.header = pd.DataFrame(
            [labels, popdata, popflag, locdata, phenotype],
            index=["labels", "popdata", "popflag", "locdata", "phenotype"]).T


    # def check_files(self):
    #     "check file format and get quick stats."


    @property
    def result_files(self):
        "returns a list of files that have finished structure"
        reps = os.path.join(
            self.workdir, 
            self.name + "-K-*-rep-*_f")
        repfiles = sorted(glob.glob(reps))
        return repfiles


    def run(
        self, 
        kpop, 
        nreps, 
        seed=12345,
        ipyclient=None, 
        force=False, 
        quiet=False,
        show_cluster=False, 
        auto=False):
        """
        Distribute structure jobs in parallel. 

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

        # bail out if object was init with load=True
        if self._load_only:
            sys.stderr.write(
                "To call .run() you must re-init the structure object without load_only=True."
            )
            return

        # print software info
        

        # load the parallel client
        pool = Parallel(
            tool=self, 
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={
                "force": force, "nreps": nreps,
                "kpop": kpop, "seed": seed, "quiet": False},
            )
        pool.wrap_run()


    def _run(self, kpop, nreps=1, seed=12345, force=False, ipyclient=None, quiet=False):
        """ 
        submits a job to run on the cluster and returns an asynchronous result
        object. K is the number of populations, randomseed if not set will be 
        randomly drawn, ipyclient if not entered will raise an error. If nreps
        is set then multiple jobs will be started from new seeds, each labeled
        by its replicate number. If force=True then replicates will be 
        overwritten, otherwise, new replicates will be created starting with 
        the last file N found in the workdir. 

        Parameters:
        -----------
        kpop: (int or list)
            The MAXPOPS parameter in structure, i.e., the number of populations
            assumed by the model (K). You can enter multiple integers as a list.

        nreps: (int):
            Number of independent replicate runs starting from distinct 
            random seeds.

        ipyclient: (ipyparallel.Client Object)
            An ipyparallel client connected to an ipcluster instance. This is 
            used to distribute parallel jobs. If you don't know what this is
            then you should use the option 'auto=True' instead.

        auto: (bool):
            Automatically start an ipcluster instance on this node connected
            to all available cores to use for multiprocessing, and shutdown 
            the ipcluster when jobs are finished. This will enforce blocking
            until all submitted jobs are finished.

        seed: (int):
            Random number seed.

        force: (bool):
            If force is true then old replicates are removed and new reps start
            from rep-0. Otherwise, new reps start at end of existing rep numbers.

        quiet: (bool)
            Whether to print number of jobs submitted to stderr

        Example: 
        ---------
        import ipyparallel as ipp

        # init structure object
        s = ipa.structure(
            name="test",
            data="mydata.str", 
            mapfile="mydata.snps.map",
            workdir="structure-results",
        )

        # modify some basic params
        s.mainparams.numreps = 100000
        s.mainparams.burnin = 10000

        # submit many jobs
        s.run(
            kpop=[2,3,4,5,6], 
            nreps=10, 
            auto=True,
            )
        """
        # initiate starting seed
        np.random.seed(seed)

        # start load balancer
        if ipyclient:
            lbview = ipyclient.load_balanced_view()

        # build new requested jobs
        if isinstance(kpop, int):
            kpop = [kpop]
        jobs = []
        for k in kpop:
            for rep in range(nreps):
                jobs.append((k, rep))

        # get previous results 
        res = {}
        for i in self.result_files:
            d1 = i.split("-")[-3:]
            tup = (int(d1[0]), int(d1[-1][0]))
            res[tup] = i

        # if force then remove old files and leave tups in jobs
        if force:
            for i in res.values():
                os.remove(i)

        # if not force then remove tups from jobs and use existing results
        else:
            for i in res:
                jobs.remove(i)

        # track jobs
        njobs = len(jobs)
        printstr = "running {} structure jobs".format(njobs)
        prog = ProgressBar(njobs, None, printstr)

        # bail out now if all jobs are completed
        if not jobs:
            print("{} finished jobs. No further jobs to run.".format(len(res)))
            return

        # print initial progress bar
        prog.finished = 0
        prog.update()

        # submit jobs to queue
        for job in jobs:
            # sample random seed for this rep
            self.extraparams.seed = np.random.randint(0, 1e9, 1)[0]

            # prepare files (randomly subsamples snps in each rep)
            k, rep = job
            mname, ename, sname = self.write_structure_files(k, rep)

            # build args with new tmp file strings
            args = [
                self.STRUCTURE,
                mname, ename, sname,
                self.name, 
                self.workdir,
                self.extraparams.seed, 
                self.ntaxa, 
                self.nsites,
                k,
                rep,
            ]

            # call structure (submit job to queue)
            rasync = lbview.apply(_call_structure, *(args))
            name = "{}-{}".format(k, rep)
            self.rasyncs[name] = rasync
            prog.update()

        # track progress...
        while 1:
            fins = [i for i in self.rasyncs if self.rasyncs[i].ready()]
            for i in fins:
                prog.finished += 1
                del self.rasyncs[i]
            prog.update()
            time.sleep(0.9)

            if not self.rasyncs:
                print("")
                break


    def write_structure_files(self, kpop, rep=1):
        """ 
        Prepares input files for running structure. Users typically do not need
        to call this function since it is called internally by .run(). But it
        is optionally available here in case users wish to generate files and 
        run structure separately.
        """

        # check params
        self.mainparams.numreps = int(self.mainparams.numreps)
        self.mainparams.burnin = int(self.mainparams.burnin)

        # write tmp files for the job. Rando avoids filename conflict.
        mname = os.path.join(
            self.workdir, 
            "tmp-{}-{}-{}.mainparams.txt".format(self.name, kpop, rep))
        ename = os.path.join(
            self.workdir, 
            "tmp-{}-{}-{}.extraparams.txt".format(self.name, kpop, rep))
        sname = os.path.join(
            self.workdir, 
            "tmp-{}-{}-{}.strfile.txt".format(self.name, kpop, rep))
        tmp_m = open(mname, 'w')
        tmp_e = open(ename, 'w')
        tmp_s = open(sname, 'w')

        # write params files
        tmp_m.write(self.mainparams._asfile())
        tmp_e.write(self.extraparams._asfile())

        # get sequence array
        if self.subsample_snps:
            subs = self._ext.subsample_snps(quiet=True)
        else:
            subs = self.snps.copy()

        arr = np.zeros(
            (self.header.shape[0], subs.shape[1]),
            dtype=np.int8,
        )

        # update object for subsampled array
        self.ntaxa = subs.shape[0]
        self.nsites = subs.shape[1]

        # build split row genotype matrix (ugh, what a terrible format)
        for idx in range(subs.shape[0]):

            sidx = idx * 2

            # on odd numbers subtract both
            arr[sidx] = subs[idx].copy()
            arr[sidx + 1] = subs[idx].copy()

            arr[sidx][(arr[sidx] == 1) | (arr[sidx] == 2)] -= 1
            arr[sidx + 1][arr[sidx + 1] == 2] -= 1
        arr[arr == 9] = -9

        # convert to dataframe for writing
        df = pd.concat([self.header, pd.DataFrame(arr)], axis=1)
        df.to_csv(tmp_s, sep="\t", header=False, index=False)

        # close tmp files
        tmp_m.close()
        tmp_e.close()
        tmp_s.close()
        return mname, ename, sname



    def get_clumpp_table(self, kvalues, max_var_multiple=0, quiet=False):
        """
        Returns a dictionary of results tables for making structure barplots.
        This calls the same functions used in get_evanno_table() to call 
        CLUMPP to permute replicates.

        Parameters:
        -----------
        kvalues : list or int
            A kvalue or list of kvalues to run CLUMPP on and return a 
            results table. 

        max_var_multiple: int
            A multiplier value to use as a filter for convergence of runs. 
            Default=0=no filtering. As an example, if 10 replicates 
            were run then the variance of the run with the minimum variance is
            used as a benchmark. If other runs have a variance that is N times 
            greater then that run will be excluded. Remember, if replicate runs 
            sampled different distributions of SNPs then it is not unexpected that 
            they will have very different variances. However, you may still want 
            to exclude runs with very high variance since they likely have 
            not converged. 

        Returns:
        --------
        table : dict or pd.DataFrame
            A dictionary of dataframes with admixture proportions.
        """
        ## do not allow bad vals
        if max_var_multiple:
            if max_var_multiple < 1:
                raise ValueError('max_var_multiple must be >1')

        if isinstance(kvalues, int):
            return _get_clumpp_table(self, kvalues, max_var_multiple, quiet)
        else:
            tabledict = {}
            for kpop in kvalues:
                table = _get_clumpp_table(self, kpop, max_var_multiple, quiet)
                tabledict[kpop] = table
            return tabledict



    def get_evanno_table(self, kvalues, max_var_multiple=0, quiet=False):
        """
        Calculates the Evanno table from results files for tests with 
        K-values in the input list kvalues. The values lnPK, lnPPK,
        and deltaK are calculated. The max_var_multiplier arg can be used
        to exclude results files based on variance of the likelihood as a 
        proxy for convergence. 

        Parameters:
        -----------
        kvalues : list
            The list of K-values for which structure was run for this object.
            e.g., kvalues = [3, 4, 5]

        max_var_multiple: int
            A multiplier value to use as a filter for convergence of runs. 
            Default=0=no filtering. As an example, if 10 replicates 
            were run then the variance of the run with the minimum variance is
            used as a benchmark. If other runs have a variance that is N times 
            greater then that run will be excluded. Remember, if replicate runs 
            sampled different distributions of SNPs then it is not unexpected that 
            they will have very different variances. However, you may still want 
            to exclude runs with very high variance since they likely have 
            not converged. 

        quiet: bool
            Suppresses printed messages about convergence.

        Returns:
        --------
        table : pandas.DataFrame
            A data frame with LPK, LNPPK, and delta K. The latter is typically
            used to find the best fitting value of K. But be wary of over
            interpreting a single best K value. 
        """
        ## do not allow bad vals
        if max_var_multiple:
            if max_var_multiple < 1:
                raise ValueError('max_variance_multiplier must be >1')

        table = _get_evanno_table(self, kvalues, max_var_multiple, quiet)
        return table



def _call_structure(STRUCTURE, mname, ename, sname, name, workdir, seed, ntaxa, nsites, kpop, rep):
    "make the subprocess call to structure"
    # create call string
    outname = os.path.join(workdir, "{}-K-{}-rep-{}".format(name, kpop, rep))

    cmd = [
        STRUCTURE,
        "-m", mname, 
        "-e", ename, 
        "-K", str(kpop),
        "-D", str(seed), 
        "-N", str(ntaxa), 
        "-L", str(nsites),
        "-i", sname, 
        "-o", outname,
    ]

    # call the shell function
    proc = sps.Popen(cmd, stdout=sps.PIPE, stderr=sps.STDOUT)
    comm = proc.communicate()

    if proc.returncode:
        raise IPyradError(comm[0])

    # cleanup
    oldfiles = [mname, ename, sname]
    for oldfile in oldfiles:
        if os.path.exists(oldfile):
            os.remove(oldfile)
    return comm



class _MainParams(Params):
    """
    A dictionary object of mainparams parameter arguments to STRUCTURE. 
    See STRUCTURE docs for details on their function. Modify by setting as 
    an object or dict, e.g.:

    struct.mainparams.popflag = 1
    struct.mainparams["popflag"] = 1
    """
    def __init__(self):
        self.burnin = int(10000)
        self.numreps = int(50000)
        self.ploidy = 2
        self.missing = -9
        self.onerowperind = 0
        self.label = 1
        self.popdata = 0
        self.popflag = 0
        self.locdata = 0
        self.phenotype = 0
        self.extracols = 0
        self.markernames = 0
        self.recessivealleles = 0
        self.mapdistances = 0
        self.phased = 0
        self.phaseinfo = 0
        self.markovphase = 0
        self.notambiguous = -999

    def _asfile(self):
        return _MAINPARAMS.format(**self.__dict__)



class _ExtraParams(Params):
    """
    A dictionary object of extraparams parameter arguments to STRUCTURE. 
    See STRUCTURE docs for details on their function. Modify by setting as 
    an object or dict, e.g.:

    struct.extraparams.noadmix = 1
    struct.extraparams["noadmix"] = 1
    """
    def __init__(self):
        self.noadmix = 0
        self.linkage = 0
        self.usepopinfo = 0
        self.locprior = 0
        self.freqscorr = 1
        self.onefst = 0
        self.inferalpha = 1
        self.popalphas = 0
        self.alpha = 1.0

        self.inferlambda = 0
        self.popspecificlambda = 0
        self.lambda_ = 1.0

        self.fpriormean = 0.01
        self.fpriorsd = 0.05
        self.unifprioralpha = 1
        self.alphamax = 10.0
        self.alphapriora = 1.0
        self.alphapriorb = 2.0
        self.log10rmin = -4.0
        self.log10rmax = 1.0
        self.log10rpropsd = 0.1
        self.log10rstart = -2.0

        self.gensback = 2
        self.migrprior = 0.01
        self.pfrompopflagonly = 0

        self.locispop = 0
        self.locpriorinit = 1.0
        self.maxlocprior = 20.0

        self.printnet = 1               ## do we want these to print ?
        self.printlambda = 1            ##
        self.printqsum = 1              ##
        self.sitebysite = 0
        self.printqhat = 0
        self.updatefreq = 10000
        self.printlikes = 0
        self.intermedsave = 0
        self.echodata = 0       
        self.ancestdist = 0
        self.numboxes = 1000
        self.ancestpint = 0.90

        self.computeprob = 1
        self.admburnin = 500
        self.alphapropsd = 0.025
        self.startatpopinfo = 0
        self.randomize = 0
        self.seed = 12345
        self.metrofreq = 10
        self.reporthitrate = 0


    def _asfile(self):
        return _EXTRAPARAMS.format(**self.__dict__)



class _ClumppParams(Params):
    """
    A dictionary object of params arguments to CLUMPP.
    See CLUMPP docs for details on their function. Modify by setting as 
    an object or dict, e.g.:

    struct.clumppparams.datatype = 1
    struct.clumpparams["datatype"] = 1
    """

    def __init__(self):
        self.datatype = 0
        self.indfile = 0
        self.outfile = 0
        self.popfile = 0
        self.miscfile = 0
        #self.kpop = 3
        #self.c = 3
        #self.r = 10
        self.m = 3
        self.w = 1
        self.s = 2
        self.greedy_option = 2
        self.repeats = 50000
        self.permutationsfile = 0
        self.print_permuted_data = 0
        self.permuted_datafile = 0
        self.print_every_perm = 0
        self.every_permfile = 0
        self.permfile = 0
        self.print_random_inputorder = 0
        self.random_inputorderfile = 0
        self.override_warnings = 0
        self.order_by_run = 1


    def _asfile(self):
        return _CLUMPPARAMS.format(**self.__dict__)



def _get_clumpp_table(self, kpop, max_var_multiple, quiet):
    "private function to clumpp results"

    ## concat results for k=x
    reps, excluded = _concat_reps(self, kpop, max_var_multiple, quiet)
    if reps:
        ninds = reps[0].inds
        nreps = len(reps)
    else:
        ninds = nreps = 0
    if not reps:
        return "no result files found"

    # return if kpop is 1
    if kpop == 1:
        if not quiet:
            print("Nothing to permute or plot for kpop=1, but these results can be used for Evanno.")
            return

    clumphandle = os.path.join(self.workdir, "tmp.clumppparams.txt")
    self.clumppparams.kpop = kpop
    self.clumppparams.c = ninds
    self.clumppparams.r = nreps
    with open(clumphandle, 'w') as tmp_c:
        tmp_c.write(self.clumppparams._asfile())
    
    ## create CLUMPP args string
    outfile = os.path.join(self.workdir, 
        "{}-K-{}.clumpp.outfile".format(self.name, kpop))
    indfile = os.path.join(self.workdir, 
        "{}-K-{}.clumpp.indfile".format(self.name, kpop))
    miscfile = os.path.join(self.workdir, 
        "{}-K-{}.clumpp.miscfile".format(self.name, kpop))

    # shorten filenames because clumpp can't handle names > 50 chars.
    clumphandle = clumphandle.replace(os.path.realpath('.'), '.', 1)
    clumphandle = clumphandle.replace(os.path.expanduser('~'), '~', 1)
    indfile = indfile.replace(os.path.realpath('.'), '.', 1)
    indfile = indfile.replace(os.path.expanduser('~'), '~', 1)
    outfile = outfile.replace(os.path.realpath("."), '.', 1)    
    outfile = outfile.replace(os.path.expanduser('~'), '~', 1)    
    miscfile = miscfile.replace(os.path.realpath("."), '.', 1)  
    miscfile = miscfile.replace(os.path.expanduser('~'), '~', 1)    

    cmd = [
        self.CLUMPP, clumphandle, 
        "-i", indfile,
        "-o", outfile, 
        "-j", miscfile,
        "-r", str(nreps), 
        "-c", str(ninds), 
        "-k", str(kpop),
    ]

    # call clumpp
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    comm = proc.communicate()   

    # cleanup
    for rfile in [indfile, miscfile]:
        if os.path.exists(rfile):
            os.remove(rfile)

    # parse clumpp results file
    ofile = os.path.join(
        self.workdir, 
        "{}-K-{}.clumpp.outfile".format(self.name, kpop))

    # load clumpp outfile as pandas df
    if os.path.exists(ofile):

        # load table and select columns with ancestry proportions
        table = pd.read_csv(ofile, delim_whitespace=True, header=None)
        table = table.iloc[:, 5:]
        table.columns = range(table.shape[1])

        # set index to names based on header (data file subset by imap)
        try:
            table.index = self.header.labels[::2].values           
        except ValueError as inst:
            print("N samples has changed, be sure to load imap dictionary.")
            raise inst

        # print report to user
        if not quiet:
            print(
                "[K{}] {}/{} results permuted across replicates (max_var={})."
                .format(kpop, nreps, nreps + excluded, max_var_multiple))

        # return the final table
        return table

    else:
        # TODO: shouldn't we raise an error here?
        sys.stderr.write("No files ready for {}-K-{} in {}\n"\
                         .format(self.name, kpop, self.workdir))
        if len(outfile) > 50:
            print("""
    This error may be caused by the length of your output filename. For some 
    reason Clumpp cannot handle filenames longer than 50 characters...
    {}
    """.format(" ".join(cmd)), file=sys.stderr)



def _concat_reps(self, kpop, max_var_multiple, quiet, **kwargs):
    """
    Combine structure replicates into a single indfile, 
    returns nreps, ninds. Excludes reps with too high of 
    variance (set with max_variance_multiplier) to exclude
    runs that did not converge. 
    """
   
    ## make an output handle
    outf = os.path.join(self.workdir, 
        "{}-K-{}.clumpp.indfile".format(self.name, kpop))
    
    ## combine replicates and write to indfile
    excluded = 0
    reps = []
    with open(outf, 'w') as outfile:
        repfiles = glob.glob(
            os.path.join(
                self.workdir, 
                self.name + "-K-{}-rep-*_f".format(kpop)))

        ## get result as a Rep object
        for rep in repfiles:
            result = Rep(rep, kpop=kpop)
            reps.append(result)

        ## exclude results with variance NX above (min) 
        newreps = []
        if len(reps) > 1:
            min_var_across_reps = np.min([i.var_lnlik for i in reps])
        else:
            min_var_across_reps = reps[0].var_lnlik

        ## iterate over reps
        for rep in reps:

            ## store result w/o filtering
            if not max_var_multiple:
                newreps.append(rep)
                outfile.write(rep.stable)

            ## use max-var-multiple as a filter for convergence                
            else:
                #print(
                #    rep.var_lnlik, 
                #    min_var_across_reps, 
                #    rep.var_lnlik / min_var_across_reps, 
                #    max_var_multiple)
                ## e.g., repvar is 1.05X minvar. We keep it if maxvar <= 1.05
                if (rep.var_lnlik / min_var_across_reps) <= max_var_multiple:
                    newreps.append(rep)
                    outfile.write(rep.stable)
                else:
                    excluded += 1

    return newreps, excluded



def _get_evanno_table(self, kpops, max_var_multiple, quiet):
    """
    Calculates Evanno method K value scores for a series
    of permuted clumpp results. 
    """

    ## iterate across k-vals
    kpops = sorted(kpops)
    replnliks = []

    for kpop in kpops:
        ## concat results for k=x
        reps, excluded = _concat_reps(self, kpop, max_var_multiple, quiet)

        ## report if some results were excluded
        if excluded:
            if not quiet:
                sys.stderr.write(
                "[K{}] {} reps excluded (not converged) see 'max_var_multiple'.\n"\
                .format(kpop, excluded))

        if reps:
            ninds = reps[0].inds
            nreps = len(reps)
        else:
            ninds = nreps = 0
        if not reps:
            print("no result files found")

        ## all we really need is the lnlik
        replnliks.append([i.est_lnlik for i in reps])

    ## compare lnlik and var of results
    if len(replnliks) > 1:
        lnmean = [np.mean(i) for i in replnliks]
        lnstds = [np.std(i, ddof=1) for i in replnliks]
    else:
        lnmean = replnliks
        lnstds = np.nan

    tab = pd.DataFrame(
        index=kpops,
        data={
            "Nreps": [len(i) for i in replnliks],
            "lnPK": [0] * len(kpops),
            "lnPPK": [0] * len(kpops),
            "deltaK": [0] * len(kpops),
            "estLnProbMean": lnmean, 
            "estLnProbStdev": lnstds,
        }
        )

    ## calculate Evanno's
    for kpop in kpops[1:]:
        tab.loc[kpop, "lnPK"] = tab.loc[kpop, "estLnProbMean"] \
                              - tab.loc[kpop-1, "estLnProbMean"]

    for kpop in kpops[1:-1]:
        tab.loc[kpop, "lnPPK"] = abs(tab.loc[kpop+1, "lnPK"] 
                                     - tab.loc[kpop, "lnPK"])
        tab.loc[kpop, "deltaK"] = (abs(
                                    tab.loc[kpop+1, "estLnProbMean"] - \
                                    2.0 * tab.loc[kpop, "estLnProbMean"] + \
                                    tab.loc[kpop-1, "estLnProbMean"]) / \
                                   tab.loc[kpop, "estLnProbStdev"])
        
    ## return table
    return tab



class Rep(object):
    """ parsed structure result file object """
    def __init__(self, repfile, kpop):
        self.repfile = repfile
        self.est_lnlik = 0
        self.mean_lnlik = 0
        self.var_lnlik = 0
        self.alpha = 0
        self.inds = 0
        self.kpop = kpop

        ## get table string        
        psearch = re.compile(r"\)   :  ")
        dsearch = re.compile(r"\)\s+\d+ :  ")
        self.stable = self.parse(psearch, dsearch)

        ## record if it is high variance
        self.variance_multiple = abs(self.var_lnlik / self.mean_lnlik)

    def parse(self, psearch, dsearch):
        """ parse an _f structure output file """
        stable = ""
        with open(self.repfile) as orep:
            dat = orep.readlines()
            for line in dat:
                ## stat lines
                if "Estimated Ln Prob of Data" in line:
                    self.est_lnlik = float(line.split()[-1])
                if "Mean value of ln likelihood" in line:
                    self.mean_lnlik = float(line.split()[-1])
                if "Variance of ln likelihood" in line:
                    self.var_lnlik = float(line.split()[-1])
                if "Mean value of alpha" in line:
                    self.alpha = float(line.split()[-1])

                ## matrix lines
                nonline = psearch.search(line)
                popline = dsearch.search(line)

                #if ")   :  " in line:
                if nonline:
                    ## check if sample is supervised...
                    abc = line.strip().split()
                    outstr = "{}{}{}".format(
                        " ".join([abc[0], abc[0], abc[2], 
                                  abc[0].split('.')[0]]),
                        " :  ",
                        " ".join(abc[4:])
                    )
                    self.inds += 1
                    stable += outstr + "\n"

                elif popline:
                    ## check if sample is supervised...
                    abc = line.strip().split()
                    prop = ["0.000"] * self.kpop
                    pidx = int(abc[3]) - 1
                    prop[pidx] = "1.000"
                    outstr = "{}{}{}".format(
                        " ".join([abc[0], abc[0], abc[2], 
                                  abc[0].split('.')[0]]),
                        " :  ",
                        " ".join(prop)
                    )
                    self.inds += 1
                    stable += outstr+"\n"

            stable += "\n"
        return stable




_MAINPARAMS = """
#define BURNIN         {burnin}                   //
#define NUMREPS        {numreps}                  //

#define PLOIDY         {ploidy}                   //
#define MISSING        {missing}                  //
#define ONEROWPERIND   {onerowperind}              //     

#define LABEL          {label}                     //
#define POPDATA        {popdata}                   //
#define POPFLAG        {popflag}                   //
#define LOCDATA        {locdata}                   //
#define PHENOTYPE      {phenotype}                 //
#define EXTRACOLS      {extracols}                 //

#define MARKERNAMES       {markernames}           //
#define RECESSIVEALLELES  {recessivealleles}      //
#define MAPDISTANCES      {mapdistances}          //
#define PHASED            {phased}                //
#define PHASEINFO         {phaseinfo}             //
#define MARKOVPHASE       {markovphase}           //
#define NOTAMBIGUOUS      {notambiguous}          //
"""


_EXTRAPARAMS = """
#define NOADMIX             {noadmix}                    //
#define LINKAGE             {linkage}                    //
#define USEPOPINFO          {usepopinfo}                 //
#define LOCPRIOR            {locprior}                   //
#define FREQSCORR           {freqscorr}                  //
#define ONEFST              {onefst}                     //
#define INFERALPHA          {inferalpha}                 // 
#define POPALPHAS           {popalphas}                  // 
#define ALPHA               {alpha}                      // 

#define INFERLAMBDA         {inferlambda}                  //
#define POPSPECIFICLAMBDA   {popspecificlambda}            //
#define LAMBDA              {lambda_}                       // 
             
#define FPRIORMEAN          {fpriormean}                   //
#define FPRIORSD            {fpriorsd}                     // 
#define UNIFPRIORALPHA      {unifprioralpha}               // 
#define ALPHAMAX            {alphamax}                     // 
#define ALPHAPRIORA         {alphapriora}                  // 
#define ALPHAPRIORB         {alphapriorb}                  //
#define LOG10RMIN           {log10rmin}                    //
#define LOG10RMAX           {log10rmax}                    //
#define LOG10RPROPSD        {log10rpropsd}                 //
#define LOG10RSTART         {log10rstart}                  //

#define GENSBACK            {gensback}                      //
#define MIGRPRIOR           {migrprior}                      //
#define PFROMPOPFLAGONLY    {pfrompopflagonly}              //

#define LOCISPOP            {locispop}                      //
#define LOCPRIORINIT        {locpriorinit}                  //
#define MAXLOCPRIOR         {maxlocprior}                   //

#define PRINTNET            {printnet}                     // 
#define PRINTLAMBDA         {printlambda}                  //
#define PRINTQSUM           {printqsum}                    //
#define SITEBYSITE          {sitebysite}                   // 
#define PRINTQHAT           {printqhat}                    // 
#define UPDATEFREQ          {updatefreq}                   // 
#define PRINTLIKES          {printlikes}                   // 
#define INTERMEDSAVE        {intermedsave}                 //
#define ECHODATA            {echodata}                     // 
#define ANCESTDIST          {ancestdist}                   // 
#define NUMBOXES            {numboxes}                     // 
#define ANCESTPINT          {ancestpint}                   // 

#define COMPUTEPROB         {computeprob}                    //
#define ADMBURNIN           {admburnin}                      // 
#define ALPHAPROPSD         {alphapropsd}                    // 
#define STARTATPOPINFO      {startatpopinfo}                 //
#define RANDOMIZE           {randomize}                      //
#define SEED                {seed}                           //
#define METROFREQ           {metrofreq}                      // 
#define REPORTHITRATE       {reporthitrate}                  //  

"""

_CLUMPPARAMS = """
DATATYPE                      {datatype}               #
INDFILE                       {indfile}                #
POPFILE                       {popfile}                #
OUTFILE                       {outfile}                #
MISCFILE                      {miscfile}               #
K                             {kpop}                   #
C                             {c}                      #
R                             {r}                      # 
M                             {m}                      #
W                             {w}                      #
S                             {s}                      # 
GREEDY_OPTION                 {greedy_option}          #
REPEATS                       {repeats}                #
PERMUTATIONFILE               {permutationsfile}       #
PRINT_PERMUTED_DATA           {print_permuted_data}    #
PERMUTED_DATAFILE             {permuted_datafile}      #
PRINT_EVERY_PERM              {print_every_perm}       #
EVERY_PERMFILE                {every_permfile}             #
PRINT_RANDOM_INPUTORDER       {print_random_inputorder}    # 
RANDOM_INPUTORDERFILE         {random_inputorderfile}      #
OVERRIDE_WARNINGS             {override_warnings}          #
ORDER_BY_RUN                  {order_by_run}               #
"""
