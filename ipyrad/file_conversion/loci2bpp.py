#!/usr/bin/python 

""" convert loci file to bpp format input files """

import os
import sys
import itertools
import pandas as pd
from ipyrad.assemble.util import IPyradWarningExit
#import numpy as np


def loci2bpp(name, locifile, imap, guidetree,
              minmap=None,
              maxloci=None,
              infer_sptree=0,
              infer_delimit=0,
              delimit_alg=(0, 5),
              seed=12345,
              burnin=1000,
              nsample=10000,
              sampfreq=2,
              thetaprior=(5, 5),
              tauprior=(4, 2, 1),
              traits_df=None,
              nu=0,
              kappa=0,
              useseqdata=1,
              usetraitdata=1,
              cleandata=0,
              wdir=None,
              finetune=(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
              verbose=0):
    """
    Converts loci file format to bpp file format, i.e., concatenated phylip-like
    format, and produces imap and ctl input files for bpp.

    Parameters:
    -----------
    name:
      A prefix name for output files that will be produced
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
        All species in the imap dictionary must also be present in the guidetree.

    Optional parameters:
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
        Species delimitation algorithm. This is a two-part tuple. The first value
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
    traits_df:
        A pandas DataFrame with trait data properly formatted. This means only
        quantitative traits are included, and missing values are NaN.
        The first column contains sample names, with "Indiv" as the header.
        The following columns have a header row with trait names. This script
        will write a CSV trait file with trait values mean-standardized, with
        NaN replaced by "NA", and with sample not present in IMAP removed.
    nu:
        A prior on phenotypic trait variance (0) for iBPP analysis.
    kappa:
        A prior on phenotypic trait mean (0) for iBPP analysis.
    useseqdata:
        If false inference proceeds without sequence data (can be used to test
        the effect of priors on the tree distributions).
    usetraitdata:
        If false inference proceeds without trait data (can be used to test
        the effect of priors on the trait distributions).
    cleandata:
        If 1 then sites with missing or hetero characters are removed.
    wdir:
        A working directory to write files to.
    finetune:
        See bpp documentation.
    verbose:
        If verbose=1 the ctl file text will also be written to screen (stderr).

    """
    ## check args
    if not imap:
        raise IPyradWarningExit(IMAP_REQUIRED)
    if minmap:
        if minmap.keys() != imap.keys():
            raise IPyradWarningExit(KEYS_DIFFER)

    ## working directory, make sure it exists
    if wdir:
        wdir = os.path.abspath(wdir)
        if not os.path.exists(wdir):
            raise IPyradWarningExit(" working directory (wdir) does not exist")
    else:
        wdir = os.path.curdir

    ## if traits_df then we make '.ibpp' files
    prog = 'bpp'
    if isinstance(traits_df, pd.DataFrame):
        prog = 'ibpp'
    outfile = os.path.join(wdir, "{}.{}.seq.txt".format(name, prog))
    mapfile = os.path.join(wdir, "{}.{}.imap.txt".format(name, prog))

    ## open outhandles
    fout = open(outfile, 'w')
    fmap = open(mapfile, 'w')

    ## parse the loci file
    with open(locifile, 'r') as infile:
        ## split on "//" for legacy compatibility
        loci = infile.read().strip().split("|\n")
        nloci = len(loci)

    ## all samples
    samples = list(itertools.chain(*imap.values()))

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
        if minmap:
            covd = {}
            for group, vals in imap.items():
                covd[group] = sum(["^"+i in names for i in vals])
            ## check that coverage is good enough
            if not all([covd[group] >= minmap[group] for group in minmap]):
                skip = 1

        ## too many loci?
        if maxloci:
            if nkept >= maxloci:
                skip = 1

        ## build locus as a string
        if not skip:
            ## convert to phylip with caret starter and replace - with N.
            data = ["{:<30} {}".format(i, "".join(k).replace("-", "N")) for \
                (i, k) in zip(names, seqs) if i[1:] in samples]

            ## if not empty, write to the file
            if data:
                fout.write("{} {}\n\n{}\n\n"\
                           .format(len(data), seqlen, "\n".join(data)))
                nkept += 1

    ## close up shop
    fout.close()

    ## write the imap file:
    data = ["{:<30} {}".format(val, key) for key \
            in sorted(imap) for val in imap[key]]
    fmap.write("\n".join(data))
    fmap.close()

    ## write ctl file
    write_ctl(name, imap, guidetree, nkept,
              infer_sptree, infer_delimit, delimit_alg,
              seed, burnin, nsample, sampfreq,
              thetaprior, tauprior, traits_df, nu, kappa,
              cleandata, useseqdata, usetraitdata, wdir,
              finetune, verbose)

    ## print message?
    sys.stderr.write("new files created ({} loci, {} species, {} samples)\n"\
                     .format(nkept, len(imap.keys()),
                             sum([len(i) for i in imap.values()])))
    sys.stderr.write("  {}.{}.seq.txt\n".format(name, prog))
    sys.stderr.write("  {}.{}.imap.txt\n".format(name, prog))
    sys.stderr.write("  {}.{}.ctl.txt\n".format(name, prog))
    if isinstance(traits_df, pd.DataFrame):
        sys.stderr.write("  {}.{}.traits.txt\n".format(name, prog))

    ## return the ctl file string
    return os.path.abspath(
        "{}.{}.ctl.txt".format(os.path.join(wdir, name), prog))



def write_ctl(name, imap, guidetree, nloci,
              infer_sptree, infer_delimit, delimit_alg,
              seed, burnin, nsample, sampfreq,
              thetaprior, tauprior, traits_df, nu0, kappa0,
              cleandata, useseqdata, usetraitdata, wdir,
              finetune, verbose):

    """ write outfile with any args in argdict """

    ## A string to store ctl info
    CTL = []

    ## check the tree (can do this better once we install ete3 w/ ipyrad)
    if not guidetree.endswith(";"):
        guidetree += ";"

    ## if traits_df then we make '.ibpp' files
    prog = 'bpp'
    if isinstance(traits_df, pd.DataFrame):
        prog = 'ibpp'

    ## write the top header info
    CTL.append("seed = {}".format(seed))
    CTL.append("seqfile = {}.{}.seq.txt".format(os.path.join(wdir, name), prog))
    CTL.append("Imapfile = {}.{}.imap.txt".format(os.path.join(wdir, name), prog))
    CTL.append("mcmcfile = {}.{}.mcmc.txt".format(os.path.join(wdir, name), prog))
    CTL.append("outfile = {}.{}.out.txt".format(os.path.join(wdir, name), prog))
    if isinstance(traits_df, pd.DataFrame):
        CTL.append("traitfile = {}.{}.traits.txt".format(os.path.join(wdir, name), prog))

    ## number of loci (checks that seq file exists and parses from there)
    CTL.append("nloci = {}".format(nloci))
    CTL.append("usedata = {}".format(useseqdata))
    CTL.append("cleandata = {}".format(cleandata))

    ## infer species tree
    if infer_sptree:
        CTL.append("speciestree = 1 0.4 0.2 0.1")
    else:
        CTL.append("speciestree = 0")

    ## infer delimitation (with algorithm 1 by default)
    CTL.append("speciesdelimitation = {} {} {}"\
               .format(infer_delimit, delimit_alg[0],
                       " ".join([str(i) for i in delimit_alg[1:]])))

    ## if using iBPP (if not traits_df, we assume you're using normal bpp (v.3.3+)
    if isinstance(traits_df, pd.DataFrame):
        ## check that the data frame is properly formatted
        try:
            traits_df.values.astype(float)
        except Exception:
            raise IPyradWarningExit(PDREAD_ERROR)

        ## subsample to keep only samples that are in IMAP, we do not need to
        ## standarize traits b/c ibpp does that for us.
        samples = sorted(list(itertools.chain(*imap.values())))
        didx = [list(traits_df.index).index(i) for i in traits_df.index if i not in samples]
        dtraits = traits_df.drop(traits_df.index[didx])

        ## mean standardize traits values after excluding samples
        straits = dtraits.apply(lambda x: (x - x.mean()) / (x.std()))

        ## convert NaN to "NA" cuz that's what ibpp likes, and write to file
        ftraits = straits.fillna("NA")
        traitdict = ftraits.T.to_dict("list")

        ## get reverse imap dict
        rev = {val:key for key in sorted(imap) for val in imap[key]}

        ## write trait file
        traitfile = "{}.{}.traits.txt".format(os.path.join(wdir, name), prog)
        with open(traitfile, 'w') as tout:
            tout.write("Indiv\n")
            tout.write("\t".join(
                ['Species'] + list(ftraits.columns))+"\n"
                )
            #for key in sorted(traitdict):
            #    tout.write("\t".join([key, rev[key]] + \
            #        ["^"+str(i) for i in traitdict[key]])+"\n"
            #        )
            nindT = 0
            for ikey in sorted(imap.keys()):
                samps = imap[ikey]
                for samp in sorted(samps):
                    if samp in traitdict:
                        tout.write("\t".join([samp, rev[samp]] + \
                            [str(i) for i in traitdict[samp]])+"\n"
                        )
                        nindT += 1

        #    tout.write("Indiv\n"+"\t".join(["Species"]+\
        #    ["t_{}".format(i) for i in range(len(traitdict.values()[0]))])+"\n")
        #    for key in sorted(traitdict):
        #        print >>tout, "\t".join([key, rev[key]] + \
        #                                [str(i) for i in traitdict[key]])
        #ftraits.to_csv(traitfile)

        ## write ntraits and nindT and traitfilename
        CTL.append("ntraits = {}".format(traits_df.shape[1]))
        CTL.append("nindT = {}".format(nindT))  #traits_df.shape[0]))
        CTL.append("usetraitdata = {}".format(usetraitdata))
        CTL.append("useseqdata = {}".format(useseqdata))

        ## trait priors
        CTL.append("nu0 = {}".format(nu0))
        CTL.append("kappa0 = {}".format(kappa0))

        ## remove ibpp incompatible options
        CTL.remove("usedata = {}".format(useseqdata))
        CTL.remove("speciestree = {}".format(infer_sptree))

    ## get tree values
    nspecies = str(len(imap))
    species = " ".join(sorted(imap))
    ninds = " ".join([str(len(imap[i])) for i in sorted(imap)])

    ## write the tree
    CTL.append("""\
species&tree = {} {}
                 {}
                 {}""".format(nspecies, species, ninds, guidetree))


    ## priors
    CTL.append("thetaprior = {} {}".format(*thetaprior))
    CTL.append("tauprior = {} {} {}".format(*tauprior))

    ## other values, fixed for now
    CTL.append("finetune = 1: {}".format(" ".join([str(i) for i in finetune])))
    #CTL.append("finetune = 1: 1 0.002 0.01 0.01 0.02 0.005 1.0")
    CTL.append("print = 1 0 0 0")
    CTL.append("burnin = {}".format(burnin))
    CTL.append("sampfreq = {}".format(sampfreq))
    CTL.append("nsample = {}".format(nsample))

    ## write out the ctl file
    with open("{}.{}.ctl.txt".format(os.path.join(wdir, name), prog), 'w') as out:
        out.write("\n".join(CTL))

    ## if verbose print ctl
    if verbose:
        sys.stderr.write("ctl file\n--------\n"+"\n".join(CTL)+"\n--------\n\n")




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
