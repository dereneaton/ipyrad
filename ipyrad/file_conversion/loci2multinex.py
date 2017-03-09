#!/usr/bin/python 

""" 
convert ipyrad loci file to nexus files for running a batch of 
mrbayes jobs for bucky
"""


import os
import sys
import itertools
import pandas as pd
from ipyrad.assemble.util import IPyradWarningExit
#import numpy as np

import numpy as np
import glob
from collections import Counter


def loci2multinex(name, 
                  locifile, 
                  subsamples=None,
                  outdir=None,
                  maxloci=None,
                  minSNPs=1,
                  seed=12345,
                  mcmc_burnin=int(1e6),
                  mcmc_ngen=int(2e6),
                  mcmc_sample_freq=1000,
                  ):

    """
    Converts loci file format to multiple nexus formatted files, one for 
    each locus, and writes a mrbayes block in the nexus information. The 
    mrbayes block will be set to run 2 replicate chains, for [mcmc_ngen]
    generations, skipping [burnin] steps, and sampling every 
    [mcmc_sample_freq] steps. 


    Parameters:
    -----------
    name: (str)
        A prefix name for output files that will be produced
    locifile: (str)
        A .loci file produced by ipyrad.
    maxloci: (int)
        Limit the number of loci to the first N loci with sufficient sampling
        to be included in the analysis. 
    minSNPs: (int)
        Only include loci that have at least N parsimony informative SNPs.
    seed: (int)
        Random seed used for resolving ambiguities.
    burnin: (int)
        mrbayes nexus block burnin parameter used for 'sump burnin' and 'sumt burnin'.
        The number of generations to skip before starting parameter and tree sampling. 
    mcmc_ngen: (int)
        mrbayes nexus block 'mcmc ngen' and 'mcmc printfreq' parameters. We don't really
        to have any info printed to screen, so these values are set equal. This is the 
        length of the chains that will be run. 
    mcmc_sample_freq: (int)
        mrbayes nexus block 'mcmc samplefreq' parameter. The frequency of sampling from
        the mcmc chain. 

    """

    ## workdir is the top level directory (e.g., analysis-bucky)
    if outdir:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.curdir

    ## enforce full path names
    outdir = os.path.realpath(outdir)

    ## outdir is a named directory within this (e.g., analysis-bucky/subs/)
    outdir = os.path.join(outdir, "bucky-{}".format(name))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        ## remove {number}.nex files in this folder
        ofiles = glob.glob(os.path.join(outdir, "[0-9].nex*"))
        for ofile in ofiles:
            os.remove(ofile)

    ## parse the locifile to a list
    with open(locifile) as infile:
        loci = infile.read().strip().split("|\n")

    ## convert subsamples to a set
    if not subsamples:
        ## get all sample names from loci
        with open(locifile) as infile:
            subs = set((i.split()[0] for i in infile.readlines() if "//" not in i))
    else:   
        subs = set(subsamples)

    ## keep track of how many loci pass
    lens = len(subs)
    nlocus = 0
    
    ## create subsampled data set
    for loc in loci:
        dat = loc.split("\n")[:-1]

        ## get names and seq from locus
        names = [i.split()[0] for i in dat]
        seqs = np.array([list(i.split()[1]) for i in dat])

        ## check that locus has required samples for each subtree
        if len(set(names).intersection(set(subs))) == lens:
            ## order the same way every time
            seqsamp = seqs[[names.index(tax) for tax in subs]]
            seqsamp = _resolveambig(seqsamp)
            pis = _count_PIS(seqsamp, minSNPs)

            if pis:
                nlocus += 1
                ## remove invariable columns given this subsampling
                copied = seqsamp.copy()
                copied[copied == "-"] = "N"
                rmcol = np.all(copied == "N", axis=0)
                seqsamp = seqsamp[:, ~rmcol]

                ## write to a nexus file
                mdict = dict(zip(subs, [i.tostring() for i in seqsamp]))
                nexmake(mdict, nlocus, outdir, mcmc_burnin, mcmc_ngen, mcmc_sample_freq)

    print "wrote {} nexus files to {}".format(nlocus, outdir)



def nexmake(mdict, nlocus, dirs, mcmc_burnin, mcmc_ngen, mcmc_sample_freq):
    """ 
    function that takes a dictionary mapping names to 
    sequences, and a locus number, and writes it as a NEXUS
    file with a mrbayes analysis block.
    """
    ## create matrix as a string
    max_name_len = max([len(i) for i in mdict])
    namestring = "{:<" + str(max_name_len+1) + "} {}\n"
    matrix = ""
    for i in mdict.items():
        matrix += namestring.format(i[0], i[1])
    
    ## write nexus block
    handle = os.path.join(dirs, "{}.nex".format(nlocus))
    with open(handle, 'w') as outnex:
        outnex.write(NEXBLOCK.format(**{
            "ntax": len(mdict), 
            "nchar": len(mdict.values()[0]), 
            "matrix": matrix,
            "ngen": mcmc_ngen, 
            "sfreq": mcmc_sample_freq, 
            "burnin": mcmc_burnin, 
            })) 



## a dictionary mapping ambiguous characters
_AMBIGS = {"R": ("G", "A"),
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
          "N": ("N", "N")}
            


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




