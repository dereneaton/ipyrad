
import os
import sys
import itertools
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
              traitdict=None,
              nu=0,
              kappa=0,
              useseqdata=1,
              usetraitdata=1,
              cleandata=0,
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
        Prior on theta (4Neu), gamma distributed. mean = a/b.
    tauprior
        Prior on tau (branch lengths), gamma distributed. mean = a/b. Last number
        is dirichlet distribution of ...
    traitdict:
        A trait dictionary for iBPP analyses. Keys should be sample names and
        values should be lists of trait values. Only quantitative traits.
        Missing values should be listed as NA.
    nu0:
        A prior on ? for iBPP trait analyses...
    kappa0:
        A prior on ? for iBPP trait analyses...
    useseqdata:
        If false inference proceeds without sequence data (can be used to test
        the effect of priors on the tree distributions).
    usetraitdata:
        If false inference proceeds without trait data (can be used to test
        the effect of priors on the trait distributions).
    cleandata:
        If 1 then sites with missing or hetero characters are removed.


    """
    ## check args
    if not imap:
        raise IPyradWarningExit(IMAP_REQUIRED)
    if minmap:
        if minmap.keys() != imap.keys():
            raise IPyradWarningExit(KEYS_DIFFER)

    ## if traitdict then we make '.ibpp' files
    prog = 'bpp'
    if traitdict:
        prog = 'ibpp'
    outfile = "{}.{}.seq.txt".format(name, prog)
    mapfile = "{}.{}.imap.txt".format(name, prog)

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
            ## convert to phylip with caret starter
            data = ["{:<30} {}".format(i, "".join(k)) for (i, k) in \
                    zip(names, seqs) if i[1:] in samples]

            ## if not empty, write to the file
            if data:
                fout.write("{} {}\n\n{}\n\n".format(len(data), seqlen, "\n".join(data)))
                nkept += 1

    ## close up shop
    fout.close()

    ## write the imap file:
    data = ["{:<30} {}".format(val, key) for key in sorted(imap) for val in imap[key]]
    fmap.write("\n".join(data))
    fmap.close()

    ## write ctl file
    write_ctl(name, imap, guidetree, nkept,
              infer_sptree, infer_delimit, delimit_alg,
              seed, burnin, nsample, sampfreq,
              thetaprior, tauprior, traitdict, nu, kappa,
              cleandata, useseqdata, usetraitdata, verbose)

    ## print message?
    sys.stderr.write("new files created ({} loci, {} species, {} samples)\n"\
                     .format(nkept, len(imap.keys()),
                             sum([len(i) for i in imap.values()])))
    sys.stderr.write("  {}.{}.seq.txt\n".format(name, prog))
    sys.stderr.write("  {}.{}.imap.txt\n".format(name, prog))
    sys.stderr.write("  {}.{}.ctl.txt\n".format(name, prog))
    if traitdict:
        sys.stderr.write("  {}.{}.traits.txt\n".format(name, prog))

    ## return the ctl file string
    return os.path.abspath(os.path.join(
        os.path.curdir, "{}.{}.ctl.txt".format(name, prog))
        )



def write_ctl(name, IMAP, guidetree, nloci,
              infer_sptree, infer_delimit, delimit_alg,
              seed, burnin, nsample, sampfreq,
              thetaprior, tauprior, traitdict, nu, kappa,
              cleandata, useseqdata, usetraitdata, verbose):

    """ write outfile with any args in argdict """

    ## A string to store ctl info
    CTL = []

    ## check the tree (can do this better once we install ete3 w/ ipyrad)
    if not guidetree.endswith(";"):
        guidetree += ";"

    ## if traitdict then we make '.ibpp' files
    prog = 'bpp'
    if traitdict:
        prog = 'ibpp'

    ## write the top header info
    CTL.append("seed = {}".format(seed))
    CTL.append("seqfile = {}.{}.seq.txt".format(name, prog))
    CTL.append("Imapfile = {}.{}.imap.txt".format(name, prog))
    CTL.append("mcmcfile = {}.{}.mcmc.txt".format(name, prog))
    CTL.append("outfile = {}.{}.out.txt".format(name, prog))
    if traitdict:
        CTL.append("traitfile = {}.{}.traits.txt".format(name, prog))

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

    ## if using iBPP (if not traits, we assume you're using normal bpp (v.3.3+)
    if traitdict:
        ## write ntraits and nindT and traitfilename
        CTL.append("ntraits = {}".format(len(traitdict.values()[0])))
        CTL.append("nindT = {}".format(len(traitdict)))
        CTL.append("usetraitdata = {}".format(usetraitdata))
        CTL.append("useseqdata = {}".format(useseqdata))

        ## trait priors
        CTL.append("nu0 = {}".format(nu))
        CTL.append("kappa0 = {}".format(kappa))

        ## remove ibpp incompatible options
        CTL.remove("usedata = {}".format(useseqdata))
        CTL.remove("speciestree = {}".format(infer_sptree))

        ## get reverse imap dict
        rev = {val:key for key in sorted(IMAP) for val in IMAP[key]}

        ## write trait file
        with open("{}.{}.traits.txt".format(name, prog), 'w') as tout:
            tout.write("Indiv\n"+"\t".join(["Species"]+\
            ["t_{}".format(i) for i in range(len(traitdict.values()[0]))])+"\n")
            for key in sorted(traitdict):
                print >>tout, "\t".join([key, rev[key]] + \
                                        [str(i) for i in traitdict[key]])

    ## get tree values
    nspecies = str(len(IMAP))
    species = " ".join(IMAP.keys())
    ninds = " ".join([str(len(i)) for i in IMAP.values()])

    ## write the tree
    CTL.append("""\
species&tree = {} {}
                 {}
                 {}""".format(nspecies, species, ninds, guidetree))


    ## priors
    CTL.append("thetaprior = {} {}".format(*thetaprior))
    CTL.append("tauprior = {} {} {}".format(*tauprior))

    ## other values, fixed for now
    CTL.append("finetune = 1: 1 0.002 0.01 0.01 0.02 0.005 1.0")
    CTL.append("print = 1 0 0 0")
    CTL.append("burnin = {}".format(burnin))
    CTL.append("sampfreq = {}".format(sampfreq))
    CTL.append("nsample = {}".format(nsample))

    ## write out the ctl file
    with open("{}.{}.ctl.txt".format(name, prog), 'w') as out:
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
