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
import subprocess as sps

# third party
import numpy as np
import pandas as pd
from ipyrad.analysis.utils import get_spans, Params
from ipyrad.assemble.utils import IPyradError


MISSING_IMPORTS = """
To use the ipa.structure module you must install two additional 
libraries which can be done with the following conda command. 

  conda install -c ipyrad structure clumpp
  conda install -c eaton-lab toytree
"""


# These are almost all default values.
class Structure(object):
    """ 
    Create and return an ipyrad.analysis Structure Object. This object allows
    you to easily enter parameter setting to submit structure jobs to run in 
    parallel on a cluster. 


    Parameters
    -----------
    name (str):
        A prefix name for all output files. 
    strfile (str):
        The path to a .str or .ustr file formatted to run in STRUCTURE. The 
        first is expected to include all SNPs for a data set, and is meant to
        be combined with a mapfile (see next), whereas the ustr file contains
        a random sample of unlinked SNPs. 
    mapfile (str):
        The path to a .snps.map file from ipyrad. This has information about 
        which SNPs are linked (from the same RAD locus) which allow for sampling
        unlinked SNPs. 

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
    def __init__(self, name, data, workdir=None, mapfile=None):

        # store args
        self.name = name
        self.data = data
        self.workdir = workdir
        self.mapfile = mapfile

        # run checks
        self.check_binaries()
        self.setupdirs()
        self.check_files()

        # params
        self.mainparams = _MainParams()
        self.extraparams = _ExtraParams()
        self.clumppparams = _ClumppParams()
        self.asyncs = []


    def check_binaries(self):
        "check for structure and clumpp"
        for binary in ['structure', 'clumpp']:
            cmd = ["which", binary]
            proc = sps.Popen(cmd, stdout=sps.PIPE, stderr=sps.PIPE)
            stdout = proc.communicate()[0]
            if not stdout:
                raise IPyradError(MISSING_IMPORTS) 


    def setupdirs(self):
        ## make workdir if it does not exist
        if self.workdir:
            self.workdir = os.path.abspath(os.path.expanduser(self.workdir))
        else:
            self.workdir = os.path.join(
                os.path.abspath('.'), 
                "analysis-structure")
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)


    def check_files(self):
        "check that strfile exists, print and parse some info from it"

        # expand path
        self.data = os.path.abspath(os.path.expanduser(self.data))
        with open(self.data) as ifile:
            lines = ifile.readlines()
            self.ntaxa = len(lines) // 2
            self.nsites = len(lines[0].strip().split()[1:])
            self.labels = [i.split('\t')[0].strip() for i in lines][::2]
            self.popdata = [i.split('\t')[1] for i in lines][::2]
            self.popflag = [i.split('\t')[2] for i in lines][::2]
            self.locdata = [i.split('\t')[3] for i in lines][::2]
            self.phenotype = [i.split('\t')[4] for i in lines][::2]
            del lines

        # if mapfile then parse it to an array
        if self.mapfile:
            with open(self.mapfile) as inmap:
                maparr = np.genfromtxt(inmap)[:, [0, 3]].astype(np.uint64)
                spans = np.zeros((maparr[-1, 0], 2), np.uint64)
                spans = get_spans(maparr, spans)
                self.maparr = spans
                self.nsites = spans.shape[0]
        else:
            self.maparr = None
            

    def _subsample(self):
        "returns a subsample of unlinked snp sites"
        spans = self.maparr
        samp = np.zeros(spans.shape[0], dtype=np.uint64)
        for i in range(spans.shape[0]):
            samp[i] = np.random.randint(spans[i, 0], spans[i, 1], 1)
        return samp


    @property
    def header(self):
        "returns a header dataframe for viewing"
        header = pd.DataFrame(
            [self.labels, self.popdata, self.popflag, self.locdata, self.phenotype],
            index=["labels", "popdata", "popflag", "locdata", "phenotype"]).T
        return header


    @property
    def result_files(self):
        "returns a list of files that have finished structure"
        reps = os.path.join(
            self.workdir, 
            self.name + "-K-*-rep-*_f")
        repfiles = glob.glob(reps)
        return repfiles


    def run(
        self,
        kpop, 
        nreps, 
        ipyclient=None,
        seed=12345, 
        force=False,
        quiet=False, 
        ):

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
        kpop: (int)
            The MAXPOPS parameter in structure, i.e., the number of populations
            assumed by the model (K). 

        nreps: (int):
            Number of independent runs starting from distinct seeds.

        ipyclient: (ipyparallel.Client Object)
            An ipyparallel client connected to an ipcluster instance. This is 
            used to manage parallel jobs. If not present a single job will
            run and block until finished (i.e., code is not parallel).

        seed: (int):
            Random number seed used for subsampling unlinked SNPs if a mapfile
            is linked to the Structure Object. 

        force: (bool):
            If force is true then old replicates are removed and new reps start
            from rep-0. Otherwise, new reps start at end of existing rep numbers.

        quiet: (bool)
            Whether to print number of jobs submitted to stderr

        Example: 
        ---------
        import ipyparallel as ipp
        import ipyrad.analysis as ipa

        ## get parallel client
        ipyclient = ipp.Client()

        ## get structure object
        s = ipa.structure(
                name="test",
                data="mydata.str", 
                mapfile="mydata.snps.map",
                workdir="structure-results",
                )

        ## modify some basic params
        s.mainparams.numreps = 100000
        s.mainparams.burnin = 10000

        ## submit many jobs
        for kpop in [3, 4, 5]:
            s.run(
                kpop=kpop, 
                nreps=10, 
                ipyclient=ipyclient,
                )

        ## block until all jobs finish
        ipyclient.wait()
        """

        # initiate starting seed
        np.random.seed(seed)

        # start load balancer
        if ipyclient:
            lbview = ipyclient.load_balanced_view()

        # remove old jobs with this same name
        handle = os.path.join(
            self.workdir, 
            self.name + "-K-{}-*".format(kpop))
        oldjobs = glob.glob(handle)
        if force or (not oldjobs):
            for job in oldjobs:
                os.remove(job)
            repstart = 0
            repend = nreps
        else:
            repstart = max([int(i.split("-")[-1][:-2]) for i in oldjobs])
            repend = repstart + nreps

        ## check that there is a ipcluster instance running
        for rep in range(repstart, repend):

            ## sample random seed for this rep
            self.extraparams.seed = np.random.randint(0, 1e9, 1)[0]

            ## prepare files (randomly subsamples snps if mapfile)
            mname, ename, sname = self.write_structure_files(kpop, rep)
            args = [
                mname, ename, sname,
                self.name, 
                self.workdir,
                self.extraparams.seed, 
                self.ntaxa, 
                self.nsites,
                kpop, 
                rep]

            if ipyclient:
                ## call structure
                rasync = lbview.apply(_call_structure, *(args))
                self.asyncs.append(rasync)

            else:
                if not quiet:
                    sys.stderr.write("submitted 1 structure job [{}-K-{}]\n"\
                                 .format(self.name, kpop))
                comm = _call_structure(*args)
                return comm

        if ipyclient:
            if not quiet:
                sys.stderr.write("submitted {} structure jobs [{}-K-{}]\n"\
                                .format(nreps, self.name, kpop))



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

        ## subsample SNPs as unlinked if a mapfile is present.
        ## & write pop data to the tmp_s file if present
        assert len(self.popdata) == len(self.labels), \
            "popdata list must be the same length as the number of taxa"

        with open(self.data) as ifile:
            _data = ifile.readlines()
            ## header
            header = np.array([i.strip().split("\t")[:5] for i in _data])
            ## seqdata
            seqdata = np.array([i.strip().split("\t")[5:] for i in _data])

            ## enter popdata into seqfile if present in self
            if any(self.popdata):
                ## set popdata in header
                header[::2, 1] = self.popdata
                header[1::2, 1] = self.popdata
                
                ## set flag to all 1s if user entered popdata but no popflag
                if not any(self.popflag):
                    self.popflag = [1 for i in self.popdata]
                    header[:, 2] = 1
                else:
                    header[::2, 2] = self.popflag
                    header[1::2, 2] = self.popflag

            ## subsample SNPs if mapfile is present
            if isinstance(self.maparr, np.ndarray):
                seqdata = seqdata[:, self._subsample()]
                
            ## write fullstr
            fullstr = np.concatenate([header, seqdata], axis=1)
            np.savetxt(tmp_s, fullstr, delimiter="\t", fmt="%s")

        ## close tmp files
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



def _call_structure(mname, ename, sname, name, workdir, seed, ntaxa, nsites, kpop, rep):
    "make the subprocess call to structure"

    # create call string
    outname = os.path.join(workdir, "{}-K-{}-rep-{}".format(name, kpop, rep))

    cmd = [
        "structure", 
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
        self.burnin = int(250000)
        self.numreps = int(1e6)
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

    clumphandle = os.path.join(self.workdir, "tmp.clumppparams.txt")
    self.clumppparams.kpop = kpop
    self.clumppparams.c = ninds
    self.clumppparams.r = nreps
    with open(clumphandle, 'w') as tmp_c:
        tmp_c.write(self.clumppparams._asfile())
    
    ## create CLUMPP args string
    outfile = os.path.join(self.workdir, 
                "{}-K-{}.outfile".format(self.name, kpop))
    indfile = os.path.join(self.workdir, 
                "{}-K-{}.indfile".format(self.name, kpop))
    miscfile = os.path.join(self.workdir, 
                "{}-K-{}.miscfile".format(self.name, kpop))

    # shorten filenames because clumpp can't handle names > 50 chars.
    for filename in [clumphandle, indfile, outfile]:
        filename = filename.replace(os.path.expanduser('~'), '~', 1)
    cmd = [
        "CLUMPP", clumphandle, 
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
        "{}-K-{}.outfile".format(self.name, kpop))
    if os.path.exists(ofile):
        csvtable = pd.read_csv(ofile, delim_whitespace=True, header=None)
        table = csvtable.loc[:, 5:]
    
        ## apply names to cols and rows
        table.columns = range(table.shape[1])
        table.index = self.labels
        if not quiet:
            sys.stderr.write(
                "[K{}] {}/{} results permuted across replicates (max_var={}).\n"\
                .format(kpop, nreps, nreps + excluded, max_var_multiple))
        return table

    else:
        # TODO: shouldn't we raise an error here?
        sys.stderr.write("No files ready for {}-K-{} in {}\n"\
                         .format(self.name, kpop, self.workdir))
        if len(outfile) > 50:
            print("""
    This error may be caused by the length of your output filename. For some 
    reason Clumpp cannot handle filenames longer than 50 characters...
    """, file=sys.stderr)





def _concat_reps(self, kpop, max_var_multiple, quiet, **kwargs):
    """
    Combine structure replicates into a single indfile, 
    returns nreps, ninds. Excludes reps with too high of 
    variance (set with max_variance_multiplier) to exclude
    runs that did not converge. 
    """
   
    ## make an output handle
    outf = os.path.join(self.workdir, 
        "{}-K-{}.indfile".format(self.name, kpop))
    
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
        dsearch = re.compile(r"\)    \d :  ")
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
                    stable += outstr+"\n"

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
