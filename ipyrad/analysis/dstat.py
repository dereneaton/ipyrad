#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
from __future__ import print_function
import pandas as pd
import numpy as np
import scipy.stats
#rom collections import OrderedDict



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
    ## dict to store results
    results = {}
    
    ## dict to store bootstrap reps
    boots = {"D12":[], 
             "D1": [],
             "D2": []}
    
    ## do bootstrap resampling with replacement
    for _ in xrange(nboots):
        samples = np.random.randint(0, len(pdf), len(pdf))
        bootdf = pd.DataFrame([pdf.loc[i] for i in samples])
        boots["D12"].append(calc_d12(bootdf))
        boots["D1"].append(calc_d1(bootdf))
        boots["D2"].append(calc_d2(bootdf))
        
    ## calculate on full data
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


### converters 
def convert(binary):
    """ translates binary (00010) to pattern (AAABA) """
    return binary.replace("0", "A").replace("1", "B")



def loci2pdf(loci, ntotal):
    """ takes ms output file created using dfoil_sim.py and 
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
    ## Create locus dict
    lcounts = {}
    pos = 0

    ## site patterns
    sitep = ["total",
             "AAAAA", "AAABA", "AABAA", "AABBA",
             "ABAAA", "ABABA", "ABBAA", "ABBBA",
             "BAAAA", "BAABA", "BABAA", "BABBA",
             "BBAAA", "BBABA", "BBBAA", "BBBBA"]

    ## iterate over loci
    for i in xrange(len(loci)):   #, loc in enumerate(loci):
        ## create pattern counter for this locus
        lcounts[i] = pd.Series({site:0 for site in sitep})
        ## get site patterns in this locus
        counts = np.array([list(x)[:ntotal] for x in loci[i]]).T
        ## for each site in this locus
        for site in counts:
            sitestr = site.tostring()
            pattern = convert(sitestr)
            if pattern[-1] != 'B':
                lcounts[i][pattern] += 1
        ## fill in meta info
        lcounts[i]["AAAAA"] = int(ntotal)-lcounts[i].values.sum()
        lcounts[i]["total"] = ntotal
        lcounts[i]["#chrom"] = "rad_"+str(i)
        lcounts[i]["pos"] = pos
        pos += ntotal
    return pd.DataFrame(lcounts).T



def ms2loci(handle):
    """ converts ms output file to loci list """
    ## read in the input file
    with open(handle, 'r') as infile:
        indata = infile.read()

    ## split by locus, skip first chunk which contains ms
    ## code and random seeds
    loci = indata.strip().split("//")[1:]
    loci = [i.strip().split("\n")[2:] for i in loci]
    return loci


## convert loci file to binary loci list

def loci2loci(handle, taxonlist):
    """ converts loci file to a binary loci list """
    ## read in the input file
    with open(handle, 'r') as infile:
        indata = infile.read()

    ## parse taxon list
    tax1 = taxonlist[0]
    tax2 = taxonlist[1]
    tax3 = taxonlist[2]
    tax4 = taxonlist[3]
    outg = taxonlist[4]

    ## split on "//" for legacy compatibility
    loci = indata.strip().split("//")[:-1]
    loci[0] = " |1|\n" + loci[0]

    ## iterate over loci to find those which have taxon sampling
    for loc in loci[:200]:
        lines = loc.split("\n", 1)[1].split()
        names = [i[1:] for i in lines[::2]]
        seqs = np.array([list(i) for i in lines[1::2]])

        if all([i in names for i in taxonlist]):
            arr = np.array([seqs[names.index(tax1)],
                            seqs[names.index(tax2)],
                            seqs[names.index(tax3)],
                            seqs[names.index(tax4)],
                            seqs[names.index(outg)]
                           ])
    return arr






if __name__ == "__main__":

    ## test input files
    MSFILE = "/home/deren/Desktop/dfoiled/" \
            +"mysim_3_1_L10000_W100_test.sites.ms.tmp"
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/" \
              +"vib_half_64tip_c85d6m40p99.loci"

    ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006', 
                 'sulcatum_D9_MEX_003', 
                 'jamesonii_D12_PWS_1636', 
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    ## get binary loci list from file
    #LOCI = ms2loci(MSFILE)
    LOCI = loci2loci(LOCIFILE, TAXONLIST)
    print(LOCI)

    ## get data frame of site counts by loci
    NTOTAL = 10000
    #PDF = loci2pdf(LOCI, NTOTAL)
    #print(PDF)

    ## calculate dstats
    #print(dstat_foil(PDF, nboots=100))    
    #print(dstat_part(PDF, nboots=100))




