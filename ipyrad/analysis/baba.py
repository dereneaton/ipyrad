#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401

from __future__ import print_function, division
import pandas as pd
import numpy as np
import scipy.stats
import numba
import sys
import os

## prettier printing
pd.options.display.float_format = '{:.4f}'.format

@numba.jit('i4(i4[:])')
def sum1d(array):
    """ a sum function that is typed for speed in numba"""
    sumn = 0.0
    for i in range(array.shape[0]):
        sumn += array[i]
    return sumn

@numba.jit('f4(i4[:], i4[:])')
def jcalc_d12(abbba, babba):
    """ D12 calc for fixed differences from pdf (pandas data frame)"""
    return sum1d(abbba-babba)/sum1d(abbba+babba)

@numba.jit('f4(i4[:], i4[:])')
def jcalc_d1(abbaa, babaa):
    """ D1 calc for fixed differences from pdf (pandas data frame)"""
    return sum1d(abbaa-babaa)/sum1d(abbaa+babaa)

@numba.jit('f4(i4[:], i4[:])')
def jcalc_d2(ababa, baaba):
    """ D2 calc for fixed differences from pdf (pandas data frame)"""
    return sum1d(ababa-baaba)/sum1d(ababa+baaba)



@numba.jit('f4[:,:](i4[:,:], i4)')#, nopython=True)
def jtestloop(vals, nboots):
    """ fast numba testloop"""
    ## create empty results array
    barr = np.zeros((nboots, 3), dtype=np.float32)
    ## fill array
    for iboot in xrange(nboots):
        samples = np.random.randint(0, vals.shape[0], vals.shape[0])
        ## create empty boot array
        bootarr = np.zeros((vals.shape[0], 19), dtype=np.int32)
        ## fill the boots array
        for irand in xrange(vals.shape[0]):
            bootarr[irand] = vals[samples[irand]]
        ## calculate Dstats from bootarr and insert to barr
        barr[iboot][0] += jcalc_d12(bootarr[:, 8], bootarr[:, 12])
        barr[iboot][1] += jcalc_d12(bootarr[:, 7], bootarr[:, 11])
        barr[iboot][2] += jcalc_d12(bootarr[:, 6], bootarr[:, 10])
    return barr




@numba.jit('f4[:,:](i4[:,:], i4[:])', nopython=True)
def jtestloop2(vals, rands):
    """ fast numba testloop"""
    ## create empty results array
    barr = np.zeros((rands.shape[0], 3), dtype=np.float32)
    ## fill array
    for iboot in xrange(rands.shape[0]):
        #samples = np.random.randint(0, vals.shape[0], vals.shape[0])
        ## create empty boot array
        bootarr = np.zeros((vals.shape[0], 19), dtype=np.int32)
        ## fill the boots array
        for irand in xrange(vals.shape[0]):
            bootarr[irand] = vals[rands[irand]]
        ## calculate Dstats from bootarr and insert to barr
        barr[iboot][0] += jcalc_d12(bootarr[:, 8], bootarr[:, 12])
        barr[iboot][1] += jcalc_d12(bootarr[:, 7], bootarr[:, 11])
        barr[iboot][2] += jcalc_d12(bootarr[:, 6], bootarr[:, 10])
    return barr



## call function to get test statistics
def jdstat_part(pdf, nboots):
    """ Function to perform bootstrap resampling to measure
    significance of partitioned D-statistics. """
    ## dict to store boot results with column order D12, D1, D2
    barr = np.zeros((nboots, 3), dtype=np.float32)

    ## do bootstrap resampling with replacement
    for iboot in xrange(nboots):
        samples = np.random.randint(0, pdf.shape[0], pdf.shape[0])
        bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
        barr[iboot] = [calc_d12(bootdf), calc_d1(bootdf), calc_d2(bootdf)]

    ## array for full data results
    rarr = np.zeros((9,), dtype=np.float16)
    rarr[0:3] = [calc_d12(pdf), calc_d1(pdf), calc_d2(pdf)]
    rarr[3:6] = [barr]


    results["D_12"] = calc_d12(pdf)
    results["D_1"] = calc_d1(pdf)
    results["D_2"] = calc_d2(pdf)

    ## get standard deviation & Z from boots
    results["D12sd"] = np.std(boots["D12"])
    results["Z12"] = abs(results["D_12"])/float(results["D12sd"])
    results["D1sd"] = np.std(boots["D1"])
    results["Z1"] = abs(results["D_1"])/float(results["D1sd"])
    results["D2sd"] = np.std(boots["D2"])
    results["Z2"] = abs(results["D_2"])/float(results["D2sd"])
    return pd.Series(results)



## Functions to calculate partitioned D-statistics
def calc_d12(pdf):
    """ D12 calc for fixed differences from pdf (pandas data frame)"""
    return sum(pdf.ABBBA-pdf.BABBA)/float(sum(pdf.ABBBA+pdf.BABBA))

def calc_d1(pdf):
    """ D1 calc for fixed differences from pdf (pandas data frame)"""
    return sum(pdf.ABBAA-pdf.BABAA)/float(sum(pdf.ABBAA+pdf.BABAA))

def calc_d2(pdf):
    """ D2 calc for fixed differences from pdf (pandas data frame)"""
    return sum(pdf.ABABA-pdf.BAABA)/float(sum(pdf.ABABA+pdf.BAABA))



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



## call function to get test statistics
def dstat_part(pdf, nboots):
    """ Function to perform bootstrap resampling to measure
    significance of partitioned D-statistics. """
    ## dict to store boot results with column order D12, D1, D2
    barr = np.zeros((nboots, 3), dtype=np.float16)

    ## do bootstrap resampling with replacement
    for iboot in xrange(nboots):
        samples = np.random.randint(0, pdf.shape[0], pdf.shape[0])
        bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
        barr[iboot] = [calc_d12(bootdf), calc_d1(bootdf), calc_d2(bootdf)]

    ## array for full data results
    rarr = np.zeros((9,), dtype=np.float16)
    rarr[0:3] = [calc_d12(pdf), calc_d1(pdf), calc_d2(pdf)]
    rarr[3:6] = [barr]


    results["D_12"] = calc_d12(pdf)
    results["D_1"] = calc_d1(pdf)
    results["D_2"] = calc_d2(pdf)

    ## get standard deviation & Z from boots
    results["D12sd"] = np.std(boots["D12"])
    results["Z12"] = abs(results["D_12"])/float(results["D12sd"])
    results["D1sd"] = np.std(boots["D1"])
    results["Z1"] = abs(results["D_1"])/float(results["D1sd"])
    results["D2sd"] = np.std(boots["D2"])
    results["Z2"] = abs(results["D_2"])/float(results["D2sd"])
    return pd.Series(results)



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



def loci2pdf(loci, where=None, ntotal=None):
    """
    takes ms output file created using dfoil_sim.py and
    creates a table of site counts similar to what the dfoil_sim.py
    script attempts to do, but correctly.

    Parameters
    ----------
    loci : list
        list of loci
    ntotal : int
        total number of sites simulated, since ms does not output
        invariant sites this is needed to calc AAAAA

    Returns
    -------
    results : pandas.Dataframe
        A DataFrame with results

    """
    ## site patterns
    sitep = ["total",
             "AAAAA", "AAABA", "AABAA", "AABBA",
             "ABAAA", "ABABA", "ABBAA", "ABBBA",
             "BAAAA", "BAABA", "BABAA", "BABBA",
             "BBAAA", "BBABA", "BBBAA", "BBBBA", "locID", "pos"]

    ## Create DataFrame
    lcounts = pd.DataFrame(0, columns=sitep,
                              index=xrange(loci.shape[0]),
                              dtype=np.int32)
    ## counter for position
    pos = 0
    ## iterate over loci
    for iloc in xrange(loci.shape[0]):
        ## get real length
        ntotal = loci[iloc][0].astype("S1").tostring().find('9')
        ## get site patterns in this locus
        counts = loci[iloc][:][:, :ntotal].astype("S1")
        ## for each site in this locus
        counts[counts == '0'] = 'B'
        counts[counts == '1'] = 'A'
        for site in counts.T:
            #print(site)
            if site[-1] not in ['9', 'B']:
                lcounts[site.tostring()][iloc] += 1
        ## fill in meta info
        #lcounts["AAAAx"][iloc] = ntotal-lcounts.iloc[iloc].values.sum()
        lcounts["total"][iloc] = ntotal
        lcounts["pos"][iloc] = int(pos)
        lcounts["locID"][iloc] = where[iloc]+1
        pos += ntotal
    i = 0
    while os.path.exists(
            os.path.join(
              os.path.curdir, "dstat_%s.csv") % i):
        i += 1
    handle = os.path.join(os.path.curdir, "dstat_%s.csv") % i
    lcounts.to_csv(handle, sep="\t")
    return lcounts



def ms2loci(handle, maxlen=200):
    """ converts ms output file to loci list """
    ## read in the input file
    with open(handle, 'r') as infile:
        indata = infile.read()

    ## split by locus, skip first chunk which contains ms code and random seeds
    loci = indata.strip().split("//")[1:]
    farr = np.ones((len(loci), 5, maxlen), dtype="int8")

    ## iterate
    for iloc in xrange(farr.shape[0]):
        arr = np.int8([list(j) for j in loci[iloc].strip().split("\n")[2:]])
        farr[iloc][:, :arr.shape[1]] = arr
    return farr



def loci2loci(handle, taxonlist, maxlen=200):
    """
    Converts .loci file to a binary loci array.
    Params:
      - handle:
        A .loci file handle
      - taxonlist:
        A list/array of lists/arrays: [[p1],[p2],[p3],[p4],[outg]] specifying
        a four or five taxon test to perform.

    Converts loci file to a binary array with XXXXA values for up to 5 taxa
    stored as floats, such that if multiple individuals were listed in a taxon
    position they are represented as a SNP frequency. If only four taxa were
    entered then the values are XXX9A, and the fourth values will be ignored
    in all computation. The loci array dimensions is (nloci, 5, maxlen), however,
    nloci is the number of loci that have sufficient taxon sampling to be included,
    and so it is typically less than the full number in the .loci file.
    """

    ## read in the input file
    with open(handle, 'r') as infile:
        ## split on "//" for legacy compatibility
        loci = infile.read().strip().split("//")[:-1]

    ## create emtpy array to fill
    nloci = len(loci)
    farr = np.ones((nloci, 5, maxlen), dtype="np.float64")
    taxc = np.zeros((nloci,))

    ## iterate over loci to find those which have taxon sampling
    for iloc in xrange(nloci):
        lines = loci[iloc].split("\n", 1)[1].split()
        names = [i[1:] for i in lines[::2]]
        seqs = np.array([list(i) for i in lines[1::2]])
        seqlen = seqs.shape[1]

        taxi = sum([i in names for i in taxonlist])
        taxc[iloc] += taxi
        if taxi == len(taxonlist):
            arr = np.zeros((5, maxlen), dtype="np.float64")
            ## find most frequent allele among outgroups and call that the
            ## outgroup allele (A). How do we break ties? For simplicity, we'll
            ## consistently choose lowest base to break ties (e.g., A over C)
            arr[-1].fill(1)

            ## fill fake data columns with 9s
            arr[:, seqlen-maxlen:].fill(9)

            ## get outgroup values
            outvals = seqs[names.index(taxonlist[4])]
            for itax in xrange(4):
                ## make 1s all sites that match to outg
                tmparr = np.int8(seqs[names.index(taxonlist[itax])] == outvals)
                ## make 9s all sites that have (N-RKSMW)
                #tmparr[
                arr[itax][:tmparr.shape[0]] = tmparr
            farr[iloc] = arr

    ## warn if no SNPs are found
    ## warn if no loci have sampling of all taxa

    #print(np.histogram(taxc, range(7)))
    ## return array that includes np.ones for loci w/o taxa
    return farr[taxc == len(taxonlist)], taxc



if __name__ == "__main__":
    ## test input files
    # MSFILE = "/home/deren/Dropbox/dfoiled_copy/" \
    #         +"mysim_3_1_L10000_W100_test.sites.ms.tmp"
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/" \
              +"vib_half_64tip_c85d6m4p99.loci"

    # ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006',
                 'sulcatum_D9_MEX_003',
                 'jamesonii_D12_PWS_1636',
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    # ## get binary loci list from file
    # LOCI = ms2loci(MSFILE)
    # #print(LOCI)
    # PDF = loci2pdf(LOCI, 100)
    # print(PDF)
    # NBOOTS = 100
    # RANDS = np.random.randint(0, NBOOTS, NBOOTS).astype(np.int32)
    # #print(dstat_part(PDF, nboots=100))
    # #print(dstat_foil(PDF, nboots=100))


    # sys.exit()
    LOCI, TAXC = loci2loci(LOCIFILE, TAXONLIST)
    print(LOCI.shape)
    print(LOCI[:2])

    # ## get data frame of site counts by loci
    NTOTAL = 100
    PDF = loci2pdf(LOCI, where=np.where(TAXC==5)[0])
    # print(PDF)

    ## calculate dstats
