#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401

from __future__ import print_function, division


from ipyrad.assemble.write_outfiles import reftrick, GETCONS
from ipyrad.assemble.util import *

import scipy.stats as st
import ipyparallel as ipp
import ipyrad as ip
import pandas as pd
import numpy as np
import numba
import itertools
import datetime
import time
import os
import toyplot


## prettier printing
# pd.options.display.float_format = '{:.4f}'.format

# @numba.jit('i4(i4[:])')
# def sum1d(array):
#     """ a sum function that is typed for speed in numba"""
#     sumn = 0.0
#     for i in range(array.shape[0]):
#         sumn += array[i]
#     return sumn

# @numba.jit('f4(i4[:], i4[:])')
# def jcalc_d12(abbba, babba):
#     """ D12 calc for fixed differences from pdf (pandas data frame)"""
#     return sum1d(abbba-babba)/sum1d(abbba+babba)

# @numba.jit('f4(i4[:], i4[:])')
# def jcalc_d1(abbaa, babaa):
#     """ D1 calc for fixed differences from pdf (pandas data frame)"""
#     return sum1d(abbaa-babaa)/sum1d(abbaa+babaa)

# @numba.jit('f4(i4[:], i4[:])')
# def jcalc_d2(ababa, baaba):
#     """ D2 calc for fixed differences from pdf (pandas data frame)"""
#     return sum1d(ababa-baaba)/sum1d(ababa+baaba)


# @numba.jit('f4[:,:](i4[:,:], i4)')#, nopython=True)
# def jtestloop(vals, nboots):
#     """ fast numba testloop"""
#     ## create empty results array
#     barr = np.zeros((nboots, 3), dtype=np.float32)
#     ## fill array
#     for iboot in xrange(nboots):
#         samples = np.random.randint(0, vals.shape[0], vals.shape[0])
#         ## create empty boot array
#         bootarr = np.zeros((vals.shape[0], 19), dtype=np.int32)
#         ## fill the boots array
#         for irand in xrange(vals.shape[0]):
#             bootarr[irand] = vals[samples[irand]]
#         ## calculate Dstats from bootarr and insert to barr
#         barr[iboot][0] += jcalc_d12(bootarr[:, 8], bootarr[:, 12])
#         barr[iboot][1] += jcalc_d12(bootarr[:, 7], bootarr[:, 11])
#         barr[iboot][2] += jcalc_d12(bootarr[:, 6], bootarr[:, 10])
#     return barr


# @numba.jit('f4[:,:](i4[:,:], i4[:])', nopython=True)
# def jtestloop2(vals, rands):
#     """ fast numba testloop"""
#     ## create empty results array
#     barr = np.zeros((rands.shape[0], 3), dtype=np.float32)
#     ## fill array
#     for iboot in xrange(rands.shape[0]):
#         #samples = np.random.randint(0, vals.shape[0], vals.shape[0])
#         ## create empty boot array
#         bootarr = np.zeros((vals.shape[0], 19), dtype=np.int32)
#         ## fill the boots array
#         for irand in xrange(vals.shape[0]):
#             bootarr[irand] = vals[rands[irand]]
#         ## calculate Dstats from bootarr and insert to barr
#         barr[iboot][0] += jcalc_d12(bootarr[:, 8], bootarr[:, 12])
#         barr[iboot][1] += jcalc_d12(bootarr[:, 7], bootarr[:, 11])
#         barr[iboot][2] += jcalc_d12(bootarr[:, 6], bootarr[:, 10])
#     return barr


# ## call function to get test statistics
# def jdstat_part(pdf, nboots):
#     """ Function to perform bootstrap resampling to measure
#     significance of partitioned D-statistics. """
#     ## dict to store boot results with column order D12, D1, D2
#     barr = np.zeros((nboots, 3), dtype=np.float32)

#     ## do bootstrap resampling with replacement
#     for iboot in xrange(nboots):
#         samples = np.random.randint(0, pdf.shape[0], pdf.shape[0])
#         bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
#         barr[iboot] = [calc_d12(bootdf), calc_d1(bootdf), calc_d2(bootdf)]

#     ## array for full data results
#     rarr = np.zeros((9,), dtype=np.float16)
#     rarr[0:3] = [calc_d12(pdf), calc_d1(pdf), calc_d2(pdf)]
#     rarr[3:6] = [barr]


#     results["D_12"] = calc_d12(pdf)
#     results["D_1"] = calc_d1(pdf)
#     results["D_2"] = calc_d2(pdf)

#     ## get standard deviation & Z from boots
#     results["D12sd"] = np.std(boots["D12"])
#     results["Z12"] = abs(results["D_12"])/float(results["D12sd"])
#     results["D1sd"] = np.std(boots["D1"])
#     results["Z1"] = abs(results["D_1"])/float(results["D1sd"])
#     results["D2sd"] = np.std(boots["D2"])
#     results["Z2"] = abs(results["D_2"])/float(results["D2sd"])
#     return pd.Series(results)


# ## Functions to calculate partitioned D-statistics
# def calc_d12(pdf):
#     """ D12 calc for fixed differences from pdf (pandas data frame)"""
#     return sum(pdf.ABBBA-pdf.BABBA)/float(sum(pdf.ABBBA+pdf.BABBA))

# def calc_d1(pdf):
#     """ D1 calc for fixed differences from pdf (pandas data frame)"""
#     return sum(pdf.ABBAA-pdf.BABAA)/float(sum(pdf.ABBAA+pdf.BABAA))

# def calc_d2(pdf):
#     """ D2 calc for fixed differences from pdf (pandas data frame)"""
#     return sum(pdf.ABABA-pdf.BAABA)/float(sum(pdf.ABABA+pdf.BAABA))



## Functions to calculate D-foil
def calc_dfo(pdf):
    """ DFO calc for fixed differences from pdf """
    nleft = pdf.BABAA+pdf.BBBAA+pdf.ABABA+pdf.AAABA
    nright = pdf.BAABA+pdf.BBABA+pdf.ABBAA+pdf.AABAA
    return sum(nleft-nright)/float(sum(nleft+nright))

def calc_dil(pdf):
    """ DIL calc for fixed differences from pdf """
    nleft = pdf.ABBAA+pdf.BBBAA+pdf.BAABA+pdf.AAABA
    nright = pdf.ABABA+pdf.BBABA+pdf.BABAA+pdf.AABAA
    return sum(nleft-nright)/float(sum(nleft+nright))

def calc_dfi(pdf):
    """ DFI calc for fixed differences from pdf """
    nleft = pdf.BABAA+pdf.BABBA+pdf.ABABA+pdf.ABAAA
    nright = pdf.ABBAA+pdf.ABBBA+pdf.BAABA+pdf.BAAAA
    return sum(nleft-nright)/float(sum(nleft+nright))

def calc_dol(pdf):
    """ DOL calc for fixed differences from pdf """
    nleft = pdf.BAABA+pdf.BABBA+pdf.ABBAA+pdf.ABAAA
    nright = pdf.ABABA+pdf.ABBBA+pdf.BABAA+pdf.BAAAA
    return sum(nleft-nright)/float(sum(nleft+nright))



# ## call function to get test statistics
# def dstat_part(pdf, nboots):
#     """ Function to perform bootstrap resampling to measure
#     significance of partitioned D-statistics. """
#     ## dict to store boot results with column order D12, D1, D2
#     barr = np.zeros((nboots, 3), dtype=np.float16)

#     ## do bootstrap resampling with replacement
#     for iboot in xrange(nboots):
#         samples = np.random.randint(0, pdf.shape[0], pdf.shape[0])
#         bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
#         barr[iboot] = [calc_d12(bootdf), calc_d1(bootdf), calc_d2(bootdf)]

#     ## array for full data results
#     rarr = np.zeros((9,), dtype=np.float16)
#     rarr[0:3] = [calc_d12(pdf), calc_d1(pdf), calc_d2(pdf)]
#     rarr[3:6] = [barr]


#     results["D_12"] = calc_d12(pdf)
#     results["D_1"] = calc_d1(pdf)
#     results["D_2"] = calc_d2(pdf)

#     ## get standard deviation & Z from boots
#     results["D12sd"] = np.std(boots["D12"])
#     results["Z12"] = abs(results["D_12"])/float(results["D12sd"])
#     results["D1sd"] = np.std(boots["D1"])
#     results["Z1"] = abs(results["D_1"])/float(results["D1sd"])
#     results["D2sd"] = np.std(boots["D2"])
#     results["Z2"] = abs(results["D_2"])/float(results["D2sd"])
#     return pd.Series(results)



def dstat_foil(pdf, nboots):
    """ Function to perform boostrap resampling on Dfoil stats """
    ## dict to store results
    results = {}

    ## dict to store bootstrap reps
    boots = {"DFO": [],
             "DIL": [],
             "DFI": [],
             "DOL": []}

    ## do bootstrap resampling with replacement
    for _ in xrange(nboots):
        samples = np.random.randint(0, len(pdf), len(pdf))
        bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
        boots["DFO"].append(calc_dfo(bootdf))
        boots["DIL"].append(calc_dil(bootdf))
        boots["DFI"].append(calc_dfi(bootdf))
        boots["DOL"].append(calc_dol(bootdf))

    ## calculate on full data
    results["DFO"] = calc_dfo(pdf)
    results["DIL"] = calc_dil(pdf)
    results["DFI"] = calc_dfi(pdf)
    results["DOL"] = calc_dol(pdf)

    ## get standard deviation & Z from boots
    results["DFOsd"] = np.std(boots["DFO"])
    results["Z_DFO"] = abs(results["DFO"])/float(results["DFOsd"])
    results["DILsd"] = np.std(boots["DIL"])
    results["Z_DIL"] = abs(results["DIL"])/float(results["DILsd"])
    results["DFIsd"] = np.std(boots["DFI"])
    results["Z_DFI"] = abs(results["DFI"])/float(results["DFIsd"])
    results["DOLsd"] = np.std(boots["DOL"])
    results["Z_DOL"] = abs(results["DOL"])/float(results["DOLsd"])
    return pd.Series(results)


## Functions to calculate Dfoil with chi-square test """
def x_dfo(pdf):
    """ calculate DFO significance by chi-square test """
    nleft = [pdf.BABAA[i]+pdf.BBBAA[i]+pdf.ABABA[i]+pdf.AAABA[i] \
              for i in range(len(pdf))]
    nright = [pdf.BAABA[i]+pdf.BBABA[i]+pdf.ABBAA[i]+pdf.AABAA[i] \
              for i in range(len(pdf))]
    getd = [(i-j)/float(i+j) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    xstat = [((i-j)**2/float(i+j)) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    sig = [1.-scipy.stats.chi2.cdf(x, 1) for x in xstat]
    return [np.mean(getd), np.std(getd), np.mean(sig)]



def x_dil(pdf):
    """ calculate DIL significance by chi-square test """
    nleft = [pdf.ABBAA[i]+pdf.BBBAA[i]+pdf.BAABA[i]+pdf.AAABA[i] \
              for i in xrange(len(pdf))]
    nright = [pdf.ABABA[i]+pdf.BBABA[i]+pdf.BABAA[i]+pdf.AABAA[i] \
              for i in xrange(len(pdf))]
    getd = [(i-j)/float(i+j) if (i+j) > 0 else 0 for \
            i, j in zip(nleft, nright)]
    xstat = [((i-j)**2/float(i+j)) if (i+j) > 0 else 0 for \
            i, j in zip(nleft, nright)]
    sig = [1.-scipy.stats.chi2.cdf(x, 1) for x in xstat]
    return [np.mean(getd), np.std(getd), np.mean(sig)]



def x_dfi(pdf):
    """ calculate DFI significane by chi-square test """
    nleft = [pdf.BABAA[i]+pdf.BABBA[i]+pdf.ABABA[i]+pdf.ABAAA[i] \
              for i in xrange(len(pdf))]
    nright = [pdf.ABBAA[i]+pdf.ABBBA[i]+pdf.BAABA[i]+pdf.BAAAA[i] \
              for i in xrange(len(pdf))]
    getd = [(i-j)/float(i+j) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    xstat = [((i-j)**2/float(i+j)) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    sig = [1.-scipy.stats.chi2.cdf(x, 1) for x in xstat]
    return [np.mean(getd), np.std(getd), np.mean(sig)]



def x_dol(pdf):
    """ calculate DOL significance by chi-square test """
    nleft = [pdf.BAABA[i]+pdf.BABBA[i]+pdf.ABBAA[i]+pdf.ABAAA[i] \
              for i in xrange(len(pdf))]
    nright = [pdf.ABABA[i]+pdf.ABBBA[i]+pdf.BABAA[i]+pdf.BAAAA[i] \
               for i in xrange(len(pdf))]
    getd = [(i-j)/float(i+j) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    xstat = [((i-j)**2/float(i+j)) if (i+j) > 0 else 0 for \
             i, j in zip(nleft, nright)]
    sig = [1.-scipy.stats.chi2.cdf(x, 1) for x in xstat]
    return [np.mean(getd), np.std(getd), np.mean(sig)]



# def loci2pdf(loci, where=None, ntotal=None):
#     """
#     takes ms output file created using dfoil_sim.py and
#     creates a table of site counts similar to what the dfoil_sim.py
#     script attempts to do, but correctly.

#     Parameters
#     ----------
#     loci : list
#         list of loci
#     ntotal : int
#         total number of sites simulated, since ms does not output
#         invariant sites this is needed to calc AAAAA

#     Returns
#     -------
#     results : pandas.Dataframe
#         A DataFrame with results

#     """
#     ## site patterns
#     sitep = ["total",
#              "AAAAA", "AAABA", "AABAA", "AABBA",
#              "ABAAA", "ABABA", "ABBAA", "ABBBA",
#              "BAAAA", "BAABA", "BABAA", "BABBA",
#              "BBAAA", "BBABA", "BBBAA", "BBBBA", "locID", "pos"]

#     ## Create DataFrame
#     lcounts = pd.DataFrame(0, columns=sitep,
#                               index=xrange(loci.shape[0]),
#                               dtype=np.int32)
#     ## counter for position
#     pos = 0
#     ## iterate over loci
#     for iloc in xrange(loci.shape[0]):
#         ## get real length
#         ntotal = loci[iloc][0].astype("S1").tostring().find('9')
#         ## get site patterns in this locus
#         counts = loci[iloc][:][:, :ntotal].astype("S1")
#         ## for each site in this locus
#         counts[counts == '0'] = 'B'
#         counts[counts == '1'] = 'A'
#         for site in counts.T:
#             #print(site)
#             if site[-1] not in ['9', 'B']:
#                 lcounts[site.tostring()][iloc] += 1
#         ## fill in meta info
#         #lcounts["AAAAx"][iloc] = ntotal-lcounts.iloc[iloc].values.sum()
#         lcounts["total"][iloc] = ntotal
#         lcounts["pos"][iloc] = int(pos)
#         lcounts["locID"][iloc] = where[iloc]+1
#         pos += ntotal
#     i = 0
#     while os.path.exists(
#             os.path.join(
#               os.path.curdir, "dstat_%s.csv") % i):
#         i += 1
#     handle = os.path.join(os.path.curdir, "dstat_%s.csv") % i
#     lcounts.to_csv(handle, sep="\t")
#     return lcounts

# def ms2loci(handle, maxlen=200):
#     """ converts ms output file to loci list """
#     ## read in the input file
#     with open(handle, 'r') as infile:
#         indata = infile.read()

#     ## split by locus, skip first chunk which contains ms code and random seeds
#     loci = indata.strip().split("//")[1:]
#     farr = np.ones((len(loci), 5, maxlen), dtype="int8")

#     ## iterate
#     for iloc in xrange(farr.shape[0]):
#         arr = np.int8([list(j) for j in loci[iloc].strip().split("\n")[2:]])
#         farr[iloc][:, :arr.shape[1]] = arr
#     return farr

# def _loci2loci(handle, taxonlist):
#     """
#     Converts .loci file to a numpy array of loci with sufficient taxon sampling
#     for the test in 'taxonlist'. Returns the float array with frequency of SNP
#     in each row/taxon. 

#     Params:
#       - handle:
#         A .loci file handle
#       - taxonlist:
#         A list/array of lists/arrays: [[p1],[p2],[p3],[p4],[outg]] specifying
#         a four or five taxon test to perform.

#     Converts loci file to a binary array with XXXXA values for up to 5 taxa
#     stored as floats, such that if multiple individuals were listed in a taxon
#     position they are represented as a SNP frequency. If only four taxa were
#     entered then the values are XXX9A, and the fourth values will be ignored
#     in all computation. The loci array dimensions is (nloci, 5, maxlen), however,
#     nloci is the number of loci that have sufficient taxon sampling to be included,
#     and so it is typically less than the full number in the .loci file.
#     """

#     ## read in the input file
#     with open(handle, 'r') as infile:
#         ## split on "//" for legacy compatibility
#         loci = infile.read().strip().split("//")[:-1]

#     ## create emtpy array to fill
#     nloci = len(loci)
#     farr = np.ones((nloci, 5, maxlen), dtype=np.float64)
#     taxc = np.zeros((nloci,))

#     ## iterate over loci to find those which have taxon sampling
#     for iloc in xrange(nloci):
#         lines = loci[iloc].split("\n", 1)[1].split()
#         names = [i[1:] for i in lines[::2]]
#         seqs = np.array([list(i) for i in lines[1::2]])
#         seqlen = seqs.shape[1]

#         taxi = sum([i in names for i in taxonlist])
#         taxc[iloc] += taxi
#         if taxi == len(taxonlist):
#             arr = np.zeros((5, maxlen), dtype=np.float64)
#             ## find most frequent allele among outgroups and call that the
#             ## outgroup allele (A). How do we break ties? For simplicity, we'll
#             ## consistently choose lowest base to break ties (e.g., A over C)
#             arr[-1].fill(1)

#             ## fill fake data columns with 9s
#             arr[:, seqlen-maxlen:].fill(9)

#             ## get outgroup values
#             outvals = seqs[names.index(taxonlist[4])]
#             for itax in xrange(4):
#                 ## make 1s all sites that match to outg
#                 tmparr = np.int8(seqs[names.index(taxonlist[itax])] == outvals)
#                 ## make 9s all sites that have (N-RKSMW)
#                 #tmparr[
#                 arr[itax][:tmparr.shape[0]] = tmparr
#             farr[iloc] = arr

#     ## warn if no SNPs are found
#     ## warn if no loci have sampling of all taxa

#     #print(np.histogram(taxc, range(7)))
#     ## return array that includes np.ones for loci w/o taxa
#     return farr[taxc == len(taxonlist)], taxc

##############################################################


def partd(handle, test, mindict, nboots):
    pass


def dfoil(handle, test, mindict, nboots):
    pass


def batch(handle, taxdicts, mindicts=None, nboots=100, ipyclient=None, quiet=False):
    """
    parallel mode
    """

    ## an array to hold results (len(taxdicts), nboots)
    tot = len(taxdicts)
    resarr = np.zeros((tot, 4), dtype=np.float64)
    bootarr = np.zeros((tot, nboots), dtype=np.float64)

    ## if no ipyclient then assume Default is running, else raise error
    if not ipyclient:
        pass

    ## submit jobs to run on the cluster queue
    else:
        ## get client
        start = time.time()
        lbview = ipyclient.load_balanced_view()
        asyncs = {}
        idx = 0
    
        ## iterate over tests (repeats mindicts if fewer than taxdicts)
        for test, mindict in zip(taxdicts, itertools.cycle([mindicts])):
            asyncs[idx] = lbview.apply(baba, *(handle, test, mindict, nboots))
            idx += 1

        ## block until finished, print progress if requested.
        while 1:
            keys = [i for (i, j) in asyncs.items() if j.ready()]
            ## check for failures
            for job in keys:
                if not asyncs[job].successful():
                    raise IPyradWarningExit(\
                        " error: {}: {}".format(job, asyncs[job].exception()))
                ## enter results for successful jobs
                else:
                    resarr[job], bootarr[job] = asyncs[job].result()
                    del asyncs[job]

            ## count finished
            fin = tot - len(asyncs) 
            elap = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(tot, fin, " calculating D-stats  | {} | ".format(elap))
            time.sleep(0.1)
            if not asyncs:
                print("")
                break

        return resarr, bootarr



def baba(inarr, taxdict, mindict=None, nboots=100):
    """
    Takes input loci file or frequency array and a dictionary describing a 
    four or five taxon test to perform, and return D-statistic and related 
    results, and a bootstrap array. 


    Parameters
    ----------
    inarr : string or ndarray
        A string path to a .loci file produced by ipyrad. Alternatively, data
        can be entered as an ndarray. Numpy array float allele frequencies 
        with dimension (nloci, 4 or 5, maxlen). See docs.
        
    test : dict
        A dictionary mapping Sample names to test taxon names, e.g., 
        test = {'p1': ['a', 'b'], 'p2': ['c', 'd'], 'p3': ['e'], 'p4': ['f']}.
        Four taxon tests should have p1-p4 whereas five taxon tests will 
        be performed if dict keys are p1-p5. Other names will raise an error. 
        The highest value name (e.g., p5) is the outgroup. 

    mindict : dict
        A dictionary mapping minimum sample values to taxon names (e.g., p1 and 
        p2, like above). Loci that do not meet the minimum sample coverage will
        be excluded from calculations and bootstraps. If no mindict is entered
        then a default dict is used with values=1.

    nboots: int
        The number of non-parametric bootstrap replicates to perform in which 
        loci are re-sampled with replacement to the same number as in the 
        original data set (after applying mindict filtering). 

    Returns
    -------
    out : tuple of ndarrays
        The first element is a ndarray with [dstat, bootmean, bootstd, Z-score], 
        the second element is a ndarray with the bootstrap dstats. 

    boots : ndarray
        An array of dstat values over nboots replicates. See plotting functions
        for usage.
    """

    ## get data as an array from loci file
    if isinstance(inarr, str):
        arr = _loci_to_arr(inarr, taxdict, mindict)
    elif isinstance(inarr, np.ndarray):
        arr = inarr
    else:
        raise Exception("Must enter either a 'locifile' or 'arr'")

    ## run tests
    if len(taxdict) == 4:
        res, boots = _get_signif_4(arr, nboots)
    else:
        res, boots = _get_signif_5(arr, nboots)
    return res, boots



## can reduce array size significantly which would speed things up if we reduce
## to minlen being Nsnps and remove all cols without SNPs.
def _loci_to_arr(locifile, taxdict, mindict):
    """
    return a frequency array from a loci file for all loci with taxa from 
    taxdict and min coverage from mindict. 
    """

    ## read in the input file
    with open(locifile, 'r') as infile:
        loci = infile.read().strip().split("|\n")
        nloci = len(loci)

    ## get max loc length
    maxlen = 0
    for iloc in xrange(nloci):
        lines = loci[iloc].split("\n")[:-1]
        _maxl = len(lines[0]) 
        maxlen = max(maxlen, _maxl)

    ## make the array (4 or 5)
    arr = np.zeros((nloci, len(taxdict), maxlen), dtype=np.float64)
    
    ## if not mindict, make one that requires 1 in each taxon
    if not mindict:
        mindict = {i:1 for i in taxdict}

    ## raise error if names are not 'p[int]' 
    allowed_names = ['p1', 'p2', 'p3', 'p4', 'p5']
    if any([i not in allowed_names for i in taxdict]):
        raise NameError("keys in taxdict must be named 'p1' through 'p4' or 'p5'")

    ## parse key names
    keys = sorted([i for i in taxdict.keys() if i[0] == 'p'])
    outg = keys[-1]

    ## grab seqs just for the good guys
    for loc in xrange(nloci):    

        ## parse the locus
        lines = loci[loc].split("\n")[:-1]
        names = [i.split()[0] for i in lines]
        seqs = np.array([list(i.split()[1]) for i in lines])

        ## check that names cover the taxdict
        covs = [sum([j in names for j in taxdict[tax]]) >= mindict[tax] \
                for tax in taxdict]
        if all(covs):
            
            ## get the refseq 
            ref = np.where([i in taxdict[outg] for i in names])[0]
            refseq = seqs[ref].view(np.uint8)
            ancestral = np.array([reftrick(refseq, GETCONS)[:, 0]])

            ## and fill it in
            iseq = _reffreq(ancestral, refseq, GETCONS)
            arr[loc, -1, :iseq.shape[1]] = iseq 

            ## fill each other tax freq in taxdict
            for tidx, key in enumerate(keys[:-1]):

                ## get idx of names in test tax
                nidx = np.where([i in taxdict[key] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                ## get freq of sidx
                iseq = _reffreq(ancestral, sidx, GETCONS)
                ## fill it in 
                arr[loc, tidx, :iseq.shape[1]] = iseq

    ## size-down array to the number of loci that have taxa for the test
    ## this was already done in COV, this would filter for loci that 
    ## contain at least some variants, but that is not required.
    #who = np.all(np.all(arr==0, axis=2), axis=1)
    #arr = arr[~who]

    ## remove empty sites now that we subsampled
    maxlen = np.where(np.all(np.all(arr==0, axis=1), axis=0))[0].min()
    arr = arr[:, :, :maxlen]

    return arr



@numba.jit(nopython=True)
def _reffreq(refseq, iseq, consdict):
    ## empty arrays
    freq = np.zeros((1, iseq.shape[1]), dtype=np.float64)
    amseq = np.zeros((iseq.shape[0]*2, iseq.shape[1]), dtype=np.uint8)
    
    ## fill in both copies
    for seq in xrange(iseq.shape[0]):
        for col in xrange(iseq.shape[1]):  

            ## expand colums with ambigs and remove N-
            base = iseq[seq][col]
            who = consdict[:, 0] == base
            
            ## resolve heteros and enter into 
            if not np.any(who):
                amseq[seq*2][col] = base
                amseq[seq*2+1][col] = base        
            else:
                amseq[seq*2][col] = consdict[who, 1][0]
                amseq[seq*2+1][col] = consdict[who, 2][0]

    ## get as frequencies
    amseq = (refseq == amseq).astype(np.float64)
    for i in xrange(amseq.shape[0]):
        freq += amseq[i]   
    return freq / np.float64(amseq.shape[0])



@numba.jit(nopython=True)
def _prop_dstat(arr):
    
    ## numerator
    abba = ((1.-arr[:, 0]) * (arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))  
    baba = ((arr[:, 0]) * (1.-arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))
    top = abba - baba
    bot = abba + baba

    ## get statistic and avoid zero div   
    if bot.sum() != 0:
        dstat = top.sum() / float(bot.sum())
    else:
        dstat = 0
    
    return abba.sum(), baba.sum(), dstat



@numba.jit(nopython=True)
def _get_boots(arr, nboots):
    """
    return array of bootstrap D-stats
    """
    ## hold results (nboots, [dstat, ])
    boots = np.zeros((nboots,))
    
    ## iterate to fill boots
    for bidx in xrange(nboots):
        ## sample with replacement
        lidx = np.random.randint(0, arr.shape[0], arr.shape[0])
        tarr = arr[lidx]
        _, _, dstat = _prop_dstat(tarr)
        boots[bidx] = dstat
    
    ## return bootarr
    return boots



#@numba.jit(nopython=True)
def _get_signif_4(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """
    ## serial execution
    abba, baba, dstat = _prop_dstat(arr)
    boots = _get_boots(arr, nboots)
    e, s = (boots.mean(), boots.std())
    #z = np.abs(dstat) / s
    z = np.abs(e) / s
    stats = np.array([dstat, e, s, z])
    return stats, boots



@numba.jit(nopython=True)
def _get_signif_5(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """
    abba, baba, dstat = _prop_dstat(arr)
    boots = _get_boots(arr, nboots)
    e, s = (boots.mean(), boots.std())
    z = np.abs(dstat) / s
    stats = np.array([dstat, e, s, z])
    return stats, boots    



## function to convert ms simulated trees to arr 
## assumes a fixed 12 taxon tree
## assumes simlen = 100 
def _msp_to_arr(nreps, simreps, test):
    
    ## the fixed tree dictionary
    fix = {j: [i, i+1] for j, i in zip(list("abcdefghijkl"), range(0, 24, 2))}
    
    ## fill taxdict by test
    keys = ['p1', 'p2', 'p3', 'p4']
    taxs = [test[key] for key in keys]
    idxs = [list(itertools.chain(*[fix[j] for j in i])) for i in taxs]

    ## array to fill, limit to 100 len
    arr = np.zeros((nreps, sum([len(i) for i in idxs]), 100))
    
    ## iterate over reps filling arr
    for idx, trees in enumerate(simreps):
        
        ## build genotype array
        shape = trees.get_num_mutations(), trees.get_sample_size()
        garr = np.empty(shape, dtype="u1")
    
        ## fill the garr
        for variant in trees.variants():
            garr[variant.index] = variant.genotypes
        garr = garr.T
        
        ## fill my arr with freqs
        for pdx, tax in enumerate(idxs):
            freq = garr[tax]
            freq = freq.sum(axis=0) / float(freq.shape[0])
            maxsz = min(freq.shape[0], 100)
            arr[idx, pdx, :maxsz] = freq[:maxsz]
            
    ## reduce the size of arr to min  
    minl = np.where(np.all(np.all(arr==0, axis=1) == True, axis=0))[0].min()
    arr = arr[:, :, :minl]
    
    return arr


#######################################################################
## plotting funcs                                                ######
#######################################################################


def bootplot(resarr, bootarr, alpha=0.05, *args, **kwargs):
    """
    plot the distribution of bootstrap replicates and significance
    """
    ## grab/update settable defaults
    args = {"height": 400, 
            "width": 1000, 
            "label-font-size": "16px", 
            "tick-font-size": "12px"}
    args.update(kwargs)

    ## storage for all the data
    cutoff = st.norm.ppf(1-(alpha)/2) ## 1.96 for 0.05
    ntests = resarr.shape[0]
    xlines = np.zeros(ntests)

    ## get borders canvas
    canvas = toyplot.Canvas(height=args['height'], width=args['width'])
    bmin = canvas.height * 0.05
    bmax = canvas.height * 0.35
    hmin = canvas.height * 0.45
    hmax = canvas.height * 0.95
    wmin = canvas.width * 0.1
    wmax = canvas.width * 0.9

    ## space between plots (min 50)
    xlines = np.linspace(wmin, wmax, ntests+1)
    spacer = 0.75 * (xlines[1] - xlines[0])
    
    ## get line plot scale
    rmax = np.max(np.abs(bootarr))
    rmax = round(min(1.0, max(0.2, rmax+0.1*rmax)), 1)
    rmin = round(max(-1.0, min(-0.2, -1 * rmax)), 1)
    
    ## add the rest
    for idx in xrange(ntests):
        ## new data
        res = resarr[idx]
        boot = bootarr[idx]
        hist = np.histogram(boot, bins=50, range=(rmin, rmax), density=True)
        
        ## get p-value from z-score
        sign = res[3] > cutoff
        
        ## next axes
        dims = (xlines[idx], xlines[idx]+spacer, hmin, hmax)
        axes = canvas.cartesian(bounds=dims)
        axes.bars(hist, along='y', color=toyplot.color.Palette()[sign])

        ## style leftmost edge
        if idx == 0:
            ## add histograms y-label
            axes.y.label.text = "D-stat bootstraps"
            axes.y.label.style = {"font-size": args['label-font-size'],
                                  "fill": toyplot.color.near_black}
            axes.y.label.offset = 30 #wmin / 2. ## 40

            ## axes style
            axes.y.ticks.show = True
            axes.y.ticks.labels.style = {"font-size": args['tick-font-size']}
            axes.y.ticks.labels.offset = 10
        else:        
            ## styling left most
            axes.y.ticks.show = False
            axes.y.ticks.labels.show = False

        ## shared axis style
        axes.x.show = False
        axes.padding = 0.5

    ## add dash through histograms
    dims = (xlines[0], xlines[-1], hmin, hmax)# canvas.height)
    axes = canvas.cartesian(bounds=dims)
    axes.hlines(y = 0, style={"stroke-dasharray": "5, 10"})
    axes.show = False

    ## add bar plots
    dims = (xlines[0], xlines[-1], bmin, bmax)
    axes = canvas.cartesian(bounds=dims)
    axes.bars(xlines[1:], resarr[:, 3], opacity=0.75)

    ## bars axis styling
    axes.padding = 0.5
    axes.y.ticks.labels.offset = 10
    axes.y.ticks.labels.style = {"font-size": args['tick-font-size']}
    axes.x.show = False#True
    axes.hlines(y = cutoff, style={"stroke-dasharray": "5, 10"})
    
    ## bars y-label
    axes.y.label.text = "Z-score"
    axes.y.label.offset = 30 #40
    axes.y.label.style = {"font-size": args['label-font-size'],
                          "fill": toyplot.color.near_black}

    return canvas
    



def panelplot():

    ## setup canvas height in three parts 
    canvas = toyplot.Canvas(width=1000, height=1000)

    ## plot min/max
    pwmin = canvas.width * 0.05
    pwmax = canvas.width * 0.95
    phmin = canvas.height * 0.05
    phmax = canvas.height * 0.95
    pwidth = pwmax - pwmin
    pheight = phmax - phmin

    ## tree plot min/max ----------------------------------------------
    div_tree_xmin = pwmin
    div_tree_xmax = pwmin + pwidth * 0.25
    div_tree_ymin = phmin
    div_tree_ymax = phmin + pheight * 0.50
    dims = (div_tree_xmin, div_tree_xmax, 
            div_tree_ymin, div_tree_ymax)
    div_tree = canvas.cartesian(bounds=dims)
    div_tree.graph(edges, 
                   vcoordinates=coord, 
                   ewidth=3, 
                   ecolor=toyplot.color.near_black, 
                   vlshow=False,
                   vsize=0,
                   along='y')
    div_tree.show=False

    ## separator between tree and blocks for names --------------------
    namespace = pwidth * 0.15
    div_sep_xmax = div_tree_xmax + namespace 
    dims = (div_tree_xmax, div_sep_xmax, 
            div_tree_ymin, div_tree_ymax)
    div_sep = canvas.cartesian(bounds=dims)
    div_sep.fill([0, 100], [100, 100], color=toyplot.color.Palette()[1])
    div_sep.show = False

    ## get blocks for tests
    tests = range(10)
    ntests = len(tests)

    ## block spacers
    blocks = np.linspace(div_sep_xmax, pwidth, ntests+1)
    spacer = 0.75 * (blocks[1] - blocks[0])
    for bidx in xrange(ntests):
        ## create block 
        dims = (blocks[bidx], blocks[bidx]+spacer, 
                div_tree_ymin, div_tree_ymax)
        div_block = canvas.cartesian(bounds=dims)
        print dims
        
        ## functions to fill block based on taxonomy of test
        div_block.fill([0, 100], [100, 100])
        div_block.show = False 
     
    ## add bar plots ---------------------------------------------------
    div_z_ymin = div_tree_ymax + pheight * 0.05
    div_z_ymax = div_tree_ymax + pheight * 0.15
    fudge = 3.
    dims = (blocks[0], blocks[-1]-spacer/fudge, div_z_ymin, div_z_ymax)
    div_z = canvas.cartesian(bounds=dims)
    #print dims
    div_z.bars(blocks[1:], blocks[1:]+spacer, resarr[:, 3], 
               opacity=0.75, color=toyplot.color.Palette()[2])

    ## bars axis styling
    div_z.padding = 0.5
    div_z.y.ticks.labels.offset = 10
    div_z.y.ticks.labels.style = {"font-size": "12px"}# args['tick-font-size']}
    div_z.x.show = False#True
    div_z.y.spine.style = {"stroke-width": 2}# "2.5px"}
    ## bars y-label
    div_z.y.label.text = "Z-score"
    div_z.y.label.offset = "40px"
    div_z.y.label.style = {"font-size": "16px", #args['label-font-size'],
                           "fill": toyplot.color.near_black}
    ## add cutoff line
    cutoff = 3.3
    div_z.hlines(y = cutoff, 
                 style={"stroke-width": 2, #"2.5px", 
                        "stroke-dasharray": "5, 10"})

    ## plot histogram distributions --------------------------------------
    rmax = np.max(np.abs(bootsarr))
    rmax = round(min(1.0, max(0.2, rmax+0.1*rmax)), 1)
    rmin = round(max(-1.0, min(-0.2, -1 * rmax)), 1)

    ## space between plots 
    div_hist_ymin = div_z_ymax + pheight * 0.05
    div_hist_ymax = pheight

    ## add the rest
    for idx in xrange(ntests):
        ## new data
        res = resarr[idx]
        boot = bootsarr[idx]
        hist = np.histogram(boot, bins=50, range=(rmin, rmax), density=True)

        ## get p-value from z-score 
        sign = res[3] > cutoff

        ## make histogram cartesians
        dims = (blocks[idx], blocks[idx]+spacer, div_hist_ymin, div_hist_ymax)
        div_hist = canvas.cartesian(bounds=dims)
        div_hist.bars(hist, along='y', color=toyplot.color.Palette()[sign])

        ## next axes
        #dims = (blocks[idx], blocks[idx]+spacer, div_hist_ymin, div_hist_ymax)
        #axes = canvas.cartesian(bounds=dims)
        #axes.bars(hist, along='y', color=toyplot.color.Palette()[sign])

        ## style leftmost edge
        if idx == 0:
            ## add histograms y-label
            div_hist.y.label.text = "D-stat bootstraps"
            div_hist.y.label.offset = "40px" 
            div_hist.y.label.style = {"font-size": "16px", #args['label-font-size'],
                                      "fill": toyplot.color.near_black}

            ## axes style
            div_hist.y.ticks.show = True
            div_hist.y.ticks.labels.style = {"font-size": "12px", 
                                             "fill": toyplot.color.near_black}
            #args['tick-font-size']}
            div_hist.y.ticks.labels.offset = 10
            div_hist.y.spine.style = {"stroke-width": 2}#"2.5px"}
        else:        
            ## styling left most
            div_hist.y.ticks.show = False
            div_hist.y.ticks.labels.show = False
            
        ## shared axis style
        div_hist.x.show = False
        div_hist.padding = 0.5

    ## add dash through histograms
    dims = (blocks[0], blocks[-1], div_hist_ymin, div_hist_ymax)
    div_hist = canvas.cartesian(bounds=dims)
    div_hist.hlines(y = 0, style={"stroke-width": 2, #"2.5px", 
                                  "stroke-dasharray": "5, 10"})
    div_hist.show = False

    return canvas



#######################################################################
if __name__ == "__main__":

    ## test input files
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/" \
              +"vib_half_64tip_c85d6m4p99.loci"

    # ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006',
                 'sulcatum_D9_MEX_003',
                 'jamesonii_D12_PWS_1636',
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    ## calculate dstats
