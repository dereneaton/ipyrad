#!/usr/bin/env python 

"ipyrad.analysis wrapper for parallel BUCKy concordance tree analyses"

from __future__ import print_function
from builtins import range

import os
import glob
import time
import numpy as np
import subprocess as sps
from collections import Counter
from ipyrad.analysis.utils import progressbar, Params
from ipyrad.assemble.utils import IPyradError


MISSING_IMPORTS = """
To use the ipa.bucky module you must install two additional 
libraries which can be done with the following conda command. 

conda install -c ipyrad bucky
conda install -c bioconda mrbayes
"""


class Bucky(object):
    """
    Bucky analysis utility function for creating input files, setting parameters, 
    and submitting jobs to run on a parallel cluster. Converts loci/alleles
    file format data to nexus file format, The main functions are write_nex_files()
    and run()

    Parameters:
    -----------
    name:
        A name for this analysis that will be a prefix for all output files.
    data:
        A .loci or .alleles.loci file produced by ipyrad.
    workdir:
        A directory for all output files. Default is "./analysis-bucky"
    samples:
        A list of samples to include in the analysis. If empty then all samples
        present in the loci file will be used, and loci will only be included
        in the analysis if data is present for all samples. 
    mb_mcmc_ngen:
        ....
    mb_mcmc_burnin:
        ...
    mb_mcmc_sample_freq:
        ...
    bucky_alpha:
        ...
    bucky_nchains:
        ...
    bucky_nreps:
        ...
    bucky_niter:
        ...

    Attributes:
        ...u
    """
    def __init__(self, name, data, workdir=None, **kwargs):

        ## store attributes
        self.name = name
        self.data = data
        self._kwargs = {
            "minsnps": 0,
            "maxloci": None,
            "seed": None,
            "mb_mcmc_ngen": int(1e6),
            "mb_mcmc_burnin": int(1e5),
            "mb_mcmc_sample_freq": int(1e3),
            "bucky_alpha": [0.1, 1.0, 10.0],
            "bucky_nchains": 4,
            "bucky_nreps": 4,
            "bucky_niter": int(1e6),
            "copied": False,
            "samples": None,
        }
        self._kwargs.update(kwargs)

        ## check binaries
        self.check_binaries()

        ## check workdir
        if workdir:
            self.workdir = os.path.realpath(
                os.path.abspath(os.path.expanduser(workdir)))
        else:
            self.workdir = os.path.realpath(
                os.path.join(os.path.curdir, "analysis-bucky"))
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        ## set bucky parameters with defaults
        self.params = Params()
        notparams = set(["workdir", "samples", "copied"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## set defaults
        self.params.seed = np.random.randint(1e9)
        if self._kwargs["seed"]:
            self.params.seed = self._kwargs["seed"]
        self.samples = []
        if self._kwargs["samples"]:
            self.samples = set(self._kwargs["samples"])
        self._alleles = 0
        if ".alleles" in self.data:
            self._alleles = 1

        ## results files
        self.files = Params()
        self.files.data = data
        self.files.nexfiles = []
        self.files.sumfiles = []
        self.files.buckyfiles = []        

        ## accessible results
        ## get from property functions
        self.results = Params()
        self.results.concordance_trees = Params()
        self.results.population_trees = Params()
        self.results.concordance_factors = Params()

        ## during init, if results files exist then make them accessible
        ## otherwise they are only set when they are made. 
        #_buckies = glob.glob(os.path.join(self.workdir, self.name, "*.bucky*"))
        ## ...

    def check_binaries(self):
        "check for structure and clumpp"
        for binary in ['mb', 'mbsum', 'bucky']:
            cmd = ["which", binary]
            proc = sps.Popen(cmd, stdout=sps.PIPE, stderr=sps.PIPE)
            stdout = proc.communicate()[0]
            if not stdout:
                raise IPyradError(MISSING_IMPORTS) 


    @property
    def _results(self):
        """
        parse results from .bucky results files for this {workdir}/{name}.
        """
        pass


    @property
    def _trees(self):
        """ 
        parse trees from .bucky results files for this {workdir}/{name}.
        """
        pass


    def copy(self, name):
        """
        Returns a new bucky object with the same parameters as its parent. A 
        new name must be provided that is different from the parent name. 
        """
        pass


    def write_nexus_files(self, force=False, quiet=False):
        """
        Write nexus files to {workdir}/{name}/[0-N].nex, If the directory already
        exists an exception will be raised unless you use the force flag which 
        will remove all files in the directory. 

        Parameters:
        -----------
        force (bool):
            If True then all files in {workdir}/{name}/*.nex* will be removed. 

        """
        # clear existing files 
        existing = glob.glob(os.path.join(self.workdir, self.name, "*.nex"))
        if any(existing):
            if force:
                for rfile in existing:
                    os.remove(rfile)
            else:
                path = os.path.join(self.workdir, self.name)
                raise IPyradError(EXISTING_NEX_FILES.format(path))

        ## parse the loci or alleles file
        with open(self.files.data) as infile:
            loci = iter(infile.read().strip().split("|\n"))

        ## use entered samples or parse them from the file
        if not self.samples:
            with open(self.files.data) as infile:
                samples = set(
                    (i.split()[0] for i in infile.readlines() if "//" not in i)                    
                )
        else:
            samples = set(self.samples)

        ## keep track of how many loci pass filtering
        totn = len(samples)
        nloci = 0

        ## this set is just used for matching, then we randomly
        ## subsample for real within the locus so it varies 
        if self._alleles:
            msamples = {i + rbin() for i in samples}
        else:
            msamples = samples

        ## write subsampled set of loci
        for loc in loci:
            ## get names and seqs from locus
            dat = loc.split("\n")[:-1]
            try:
                names = [i.split()[0] for i in dat]
                snames = set(names)
                seqs = np.array([list(i.split()[1]) for i in dat])
            except IndexError:
                print(ALLELESBUGFIXED)
                continue

            ## check name matches
            if len(snames.intersection(msamples)) == totn:

                ## prune sample names if alleles. Done here so it is randomly
                ## different in every locus which allele is selected from 
                ## each sample (e.g., 0 or 1)
                if self._alleles:
                    _samples = [i + rbin() for i in samples]
                else:
                    _samples = samples

                ## re-order seqs to be in set order
                seqsamp = seqs[[names.index(tax) for tax in _samples]]

                ## resolve ambiguities randomly if .loci file otherwise
                ## sample one of the alleles if .alleles file.
                if not self._alleles:
                    seqsamp = _resolveambig(seqsamp)

                ## find parsimony informative sites
                if _count_PIS(seqsamp, self.params.minsnps):
                    ## keep the locus
                    nloci += 1

                    ## remove empty columns given this sampling
                    copied = seqsamp.copy()
                    copied[copied == "-"] == "N"
                    rmcol = np.all(copied == "N", axis=0)
                    seqsamp = seqsamp[:, ~rmcol]

                    ## write nexus file
                    if self._alleles:
                        ## trim off the allele number
                        samps = [i.rsplit("_", 1)[0] for i in _samples]
                        mdict = dict(zip(samps, ["".join(i) for i in seqsamp]))
                    else:
                        mdict = dict(zip(_samples, ["".join(i) for i in seqsamp]))
                    self._write_nex(mdict, nloci)

                    ## quit early if using maxloci
                    if nloci == self.params.maxloci: 
                        break


        ## print data size
        if not quiet:
            path = os.path.join(self.workdir, self.name)
            path = path.replace(os.path.expanduser("~"), "~")
            print("wrote {} nexus files to {}".format(nloci, path))



    def run(self, steps=None, ipyclient=None, force=False, quiet=False):
        """
        Submits an ordered list of jobs to a load-balancer to complete 
        the following tasks, and reports a progress bar:
        (1) Write nexus files for each locus
        (2) Run mrBayes on each locus to get a posterior of gene trees
        (3) Run mbsum (a bucky tool) on the posterior set of trees
        (4) Run Bucky on the summarized set of trees for all alpha values.

        Parameters:
        -----------
        ipyclient (ipyparallel.Client())
            A connected ipyparallel Client object used to distribute jobs
        force (bool):
            Whether to overwrite existing files with the same name and workdir
            if they exist. Default is False.
        quiet (bool):
            Whether to suppress progress information. Default is False.
        steps (list):
            A list of integers of steps to perform. This is useful if a 
            job was interrupted, or you created a new bucky object copy, 
            or you wish to run an analysis under a new set of parameters, 
            after having run it once. For example, if you finished running
            steps 1 and 2 (write nexus files and infer mrbayes posteriors), 
            but you want to rerun steps 3 and 4 with new settings, then you
            could enter `steps=[3,4]` and also `force=True` to run steps 3 
            and 4 with a new set of parameters. Default argument is None 
            which means run all steps. 
        """

        ## require ipyclient
        if not ipyclient:
            raise IPyradError("an ipyclient object is required")

        ## check the steps argument
        if not steps:
            steps = [1, 2, 3, 4]
        if isinstance(steps, (int, str)):
            steps = [int(i) for i in [steps]]
        if isinstance(steps, list):
            if not all(isinstance(i, int) for i in steps):
                raise IPyradError("steps must be a list of integers")

        ## run steps ------------------------------------------------------
        ## TODO: wrap this function so it plays nice when interrupted.
        if 1 in steps:
            self.write_nexus_files(force=force, quiet=quiet)
        if 2 in steps:
            self.run_mrbayes(force=force, quiet=quiet, ipyclient=ipyclient)
        if 3 in steps:
            self.run_mbsum(force=force, quiet=quiet, ipyclient=ipyclient)
        if 4 in steps:
            self.run_bucky(force=force, quiet=quiet, ipyclient=ipyclient)

        ## make sure jobs are done if waiting (TODO: maybe make this optional)
        ipyclient.wait()



    def _write_nex(self, mdict, nlocus):
        """ 
        function that takes a dictionary mapping names to sequences, 
        and a locus number, and writes it as a NEXUS file with a mrbayes 
        analysis block given a set of mcmc arguments.
        """

        ## create matrix as a string
        max_name_len = max([len(i) for i in mdict])
        namestring = "{:<" + str(max_name_len + 1) + "} {}\n"
        matrix = ""
        for i in mdict.items():
            matrix += namestring.format(i[0], i[1])

        ## ensure dir
        minidir = os.path.realpath(os.path.join(self.workdir, self.name))
        if not os.path.exists(minidir):
            os.makedirs(minidir)

        ## write nexus block
        handle = os.path.join(minidir, "{}.nex".format(nlocus))
        with open(handle, 'w') as outnex:
            outnex.write(NEXBLOCK.format(**{
                "ntax": len(mdict), 
                "nchar": len(list(mdict.values())[0]), 
                "matrix": matrix,
                "ngen": self.params.mb_mcmc_ngen, 
                "sfreq": self.params.mb_mcmc_sample_freq, 
                "burnin": self.params.mb_mcmc_burnin, 
                })) 



    def run_mbsum(self, ipyclient, force=False, quiet=False):
        """
        Sums two replicate mrbayes runs for each locus
        """
        minidir = os.path.realpath(os.path.join(self.workdir, self.name))
        trees1 = glob.glob(os.path.join(minidir, "*.run1.t"))
        trees2 = glob.glob(os.path.join(minidir, "*.run2.t"))

        ## clear existing files 
        existing = glob.glob(os.path.join(self.workdir, self.name, "*.sumt"))
        if any(existing):
            if force:
                for rfile in existing:
                    os.remove(rfile)
            else:
                path = os.path.join(self.workdir, self.name)
                raise IPyradError(EXISTING_SUMT_FILES.format(path))

        ## load balancer
        lbview = ipyclient.load_balanced_view()

        ## submit each to be processed
        asyncs = []
        for tidx in range(len(trees1)):
            rep1 = trees1[tidx]
            rep2 = trees2[tidx]
            outname = os.path.join(minidir, str(tidx) + ".sumt")
            rasync = lbview.apply(call_mbsum, *(rep1, rep2, outname))
            asyncs.append(rasync)

        ## track progress
        start = time.time()
        printstr = "sum replicate runs"
        while 1:
            ready = [i.ready() for i in asyncs]
            if not quiet:            
                progressbar(sum(ready), len(ready), start, printstr)
            if len(ready) == sum(ready):
                if not quiet:
                    print("")
                break
            else:
                time.sleep(0.1)

        ## check success
        for rasync in asyncs:
            if not rasync.successful():
                raise IPyradError(rasync.result())



    def run_mrbayes(self, ipyclient, force=False, quiet=False):
        """
        calls the mrbayes block in each nexus file.
        """
        ## get all the nexus files for this object
        minidir = os.path.realpath(os.path.join(self.workdir, self.name))
        nexus_files = glob.glob(os.path.join(minidir, "*.nex"))

        ## clear existing files 
        existing = glob.glob(os.path.join(minidir, "*.nex.*"))
        if any(existing):
            if force:
                for rfile in existing:
                    os.remove(rfile)
            else:
                raise IPyradError(EXISTING_NEXdot_FILES.format(minidir))

        ## load balancer
        lbview = ipyclient.load_balanced_view()

        ## submit each to be processed
        asyncs = []
        for nex in nexus_files:
            rasync = lbview.apply(call_mb, nex)
            asyncs.append(rasync)

        ## track progress
        start = time.time()
        printstr = "infer gene-tree posteriors"
        while 1:
            ready = [i.ready() for i in asyncs]
            if not quiet:            
                progressbar(sum(ready), len(ready), start, printstr)
            if len(ready) == sum(ready):
                if not quiet:
                    print("")
                break
            else:
                time.sleep(0.1)

        ## check success
        for rasync in asyncs:
            if not rasync.successful():
                raise IPyradError(rasync.result())
            


    def run_bucky(self, ipyclient, force=False, quiet=False, subname=False):
        """
        Runs bucky for a given set of parameters and stores the result 
        to the ipa.bucky object. The results will be stored by default
        with the name '{name}-{alpha}' unless a argument is passed for
        'subname' to customize the output name. 

        Parameters:
        -----------
        subname (str):
            A custom name prefix for the output files produced by the bucky
            analysis and output into the {workdir}/{name} directory.
        force (bool):
            If True then existing result files with the same name prefix
            will be overwritten. 
        quiet (bool):
            If True the progress bars will be suppressed. 
        ipyclient (ipyparallel.Client)
            An active ipyparallel client to distribute jobs to.

        """
        ## check for existing results files
        minidir = os.path.realpath(os.path.join(self.workdir, self.name))
        infiles = glob.glob(os.path.join(minidir, "*.sumt"))
        outroot = os.path.realpath(os.path.join(self.workdir, self.name))

        ## build alpha list
        if isinstance(self.params.bucky_alpha, list):
            alphas = self.params.bucky_alpha
        else:
            alphas = [self.params.bucky_alpha]

        ## load balancer
        lbview = ipyclient.load_balanced_view()

        ## submit each to be processed
        asyncs = []
        for alpha in alphas:
            pathname = os.path.join(outroot, "CF-a" + str(alpha))
            if (os.path.exists(pathname)) and (force != True):
                print("BUCKy results already exist for this filepath. " + \
                      "Use force to overwrite")
            else:
                args = [
                    alpha, 
                    self.params.bucky_nchains, 
                    self.params.bucky_nreps, 
                    self.params.bucky_niter, 
                    pathname,
                    infiles]
                rasync = lbview.apply(call_bucky, *args)
                asyncs.append(rasync)

        ## track progress
        start = time.time()
        printstr = "infer CF posteriors"
        while 1:
            ready = [i.ready() for i in asyncs]
            if not quiet:            
                progressbar(sum(ready), len(ready), start, printstr)
            if len(ready) == sum(ready):
                if not quiet:
                    print("")
                break
            else:
                time.sleep(0.1)

        ## check success
        for rasync in asyncs:
            if not rasync.successful():
                raise IPyradError(rasync.result())



#################################################

def call_bucky(alpha, nchains, nreps, niter, outname, infiles):
    "call bucky binary command"
    # build command string
    cmd = [
        "bucky", 
        "-a", str(alpha),
        "-c", str(nchains),
        "-k", str(nreps),
        "-n", str(int(niter)), 
        "-o", outname,
    ]
    for ifile in infiles:
        cmd.append(ifile)

    # call bucky
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    stdout = proc.communicate()
    if proc.returncode:
        return stdout


def call_mb(infile):
    "call mrbayes on a nex file"
    # call mrbayes
    cmd = ['mb', infile]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    stdout = proc.communicate()
        
    # check for errors
    if proc.returncode:
        return stdout


def call_mbsum(rep1, rep2, outname):
    "calls mbsum for each pair of replicates"
    cmd = [
        "mbsum", 
        "-n", "0", 
        "-o", outname,
        rep1, 
        rep2,
    ]
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    stdout = proc.communicate()

    ## check for errors
    if proc.returncode:
        return stdout


def rbin():
    "returns a random binomial as a string with an underscore before it"
    return "_{}".format(np.random.binomial(1, 0.5))


def _resolveambig(subseq):
    """ 
    Randomly resolves iupac hetero codes. This is a shortcut
    for now, we could instead use the phased alleles in RAD loci.
    """
    N = []
    for col in subseq:
        rand = np.random.binomial(1, 0.5)
        N.append([_AMBIGS[i][rand] for i in col])
    return np.array(N)



def _count_PIS(seqsamp, N):
    """ filters for loci with >= N PIS """
    counts = [Counter(col) for col in seqsamp.T if not ("-" in col or "N" in col)]
    pis = [i.most_common(2)[1][1] > 1 for i in counts if len(i.most_common(2))>1]
    if sum(pis) >= N:
        return sum(pis)
    else:
        return 0      


## GLOBALS
NEXBLOCK = """\
#NEXUS
begin data;
dimensions ntax={ntax} nchar={nchar};
format datatype=dna interleave=yes gap=- missing=N;
matrix
{matrix}
    ;

begin mrbayes;
set autoclose=yes nowarn=yes;
lset nst=6 rates=gamma;
mcmc ngen={ngen} samplefreq={sfreq} printfreq={ngen};
sump burnin={burnin};
sumt burnin={burnin};
end;
"""

EXISTING_NEX_FILES = """\
Nexus files linked to this object (i.e., same workdir & name)
already exist at {}. 
To remove the existing files and write new files use the
argument force=True. 
"""

EXISTING_NEXdot_FILES = """\
Result files linked to this object (i.e., same workdir & name)
already exist at {}. 
To remove the existing files and write new mb result files use 
the argument force=True. 
"""

EXISTING_SUMT_FILES = """\
Sumtree (sumt) files linked to this object (i.e., same workdir & name)
already exist at {}. 
To remove the existing files and write new sumt files use the
argument force=True. 
"""

ALLELESBUGFIXED = """\
Warning: encountered an error in the alleles file format. This 
is a bug that was fixed in v.0.7.2. Rerun step 7 on this data
set to ensure that the alleles file is properly formatted. 
"""


## a dictionary mapping ambiguous characters
_AMBIGS = {
    "R": ("G", "A"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("T", "C"),
    "W": ("T", "A"),
    "M": ("C", "A"), 
    "A": ("A", "A"), 
    "T": ("T", "T"), 
    "G": ("G", "G"), 
    "C": ("C", "C"), 
    "-": ("-", "-"), 
    "N": ("N", "N"),
    }
