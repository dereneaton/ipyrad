#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401

from __future__ import print_function, division


from ipyrad.assemble.write_outfiles import reftrick, GETCONS2
from ipyrad.assemble.util import *

import scipy.stats as st
import ipyparallel as ipp
import ipyrad as ip
import pandas as pd
import numpy as np
import numba
import itertools
import datetime
import types
import copy
import time
import os

## non-included imports
import toyplot
import ete3 as ete
try: 
    import msprime as ms
except ImportError:
    pass


## prettier printing
# pd.options.display.float_format = '{:.4f}'.format


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



def partd(handle, test, mindict, nboots):
    pass


def dfoil(handle, test, mindict, nboots):
    pass


def batch(handle, taxdicts, mindicts=None, nboots=1000, ipyclient=None, quiet=False):
    """
    parallel mode
    """

    ## if ms generator make into reusable list
    sims = 0
    if isinstance(handle, types.GeneratorType):
        handle = list(handle)
        sims = 1

    ## parse taxdicts into names and lists if it a dictionary
    if isinstance(taxdicts, dict):
        names, taxdicts = taxdicts.keys(), taxdicts.values()
    else:
        names = []

    ## an array to hold results (len(taxdicts), nboots)
    tot = len(taxdicts)
    resarr = np.zeros((tot, 6), dtype=np.float64)
    bootsarr = np.zeros((tot, nboots), dtype=np.float64)

    ## if no ipyclient then assume Default is running, else raise error
    if not ipyclient:
        raise IPyradError("  must provide an ipyclient Object")

    ## submit jobs to run on the cluster queue
    else:
        ## get client
        start = time.time()
        lbview = ipyclient.load_balanced_view()
        asyncs = {}
        idx = 0
    
        ## iterate over tests (repeats mindicts if fewer than taxdicts)
        for test, mindict in zip(taxdicts, itertools.cycle([mindicts])):
            ## if it's sim data then convert to an array
            if sims:
                arr = _msp_to_arr(handle, test)
                args = (arr, test, mindict, nboots)
                asyncs[idx] = lbview.apply(baba, *args)
            else:
                args = (handle, test, mindict, nboots)
                asyncs[idx] = lbview.apply(baba, *args)
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
                    _res, _bot = asyncs[job].result()
                    resarr[job] = _res.as_matrix()[:, 0]
                    bootsarr[job] = _bot
                    del asyncs[job]

            ## count finished
            fin = tot - len(asyncs) 
            elap = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(tot, fin, " calculating D-stats  | {} | ".format(elap))
            time.sleep(0.1)
            if not asyncs:
                print("\n")
                break

        ## dress up resarr as a Pandas DataFrame
        if not names:
            names = range(len(taxdicts))
        resarr = pd.DataFrame(resarr, 
                index=names,
                columns=["dstat", "bootmean", "bootstd", "ABBA", "BABA", "Z"])

        ## sort results dataframe and bootsarr to match
        resarr = resarr.sort_index()
        order = [list(resarr.index).index(i) for i in names]
        bootsarr = bootsarr[order]
        return resarr, bootsarr



def baba(inarr, taxdict, mindict=1, nboots=1000, name=0):
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
        arr, _ = _loci_to_arr(inarr, taxdict, mindict)
    elif isinstance(inarr, np.ndarray):
        arr = inarr
    #elif isinstance(inarr, types.GeneratorType):
    #    arr = _msp_to_arr(inarr, taxdict)
    #elif isinstance(inarr, list):
    #    arr = _msp_to_arr(inarr, taxdict)
    ## get data from Sim object, do not digest the ms generator
    elif isinstance(inarr, Sim):
        arr = _msp_to_arr(inarr, taxdict)
    else:
        raise Exception("Must enter either a 'locifile' or 'arr'")

    ## run tests
    if len(taxdict) == 4:

        ## get results
        res, boots = _get_signif_4(arr, nboots)
    
        ## make res into a nice DataFrame
        res = pd.DataFrame(res, 
                columns=[name],
                index=["dstat", "bootmean", "bootstd", "abba", "baba", "Z"])

    else:
        ## get results
        res, boots = _get_signif_5(arr, nboots)

        ## make int a DataFrame
        res = pd.DataFrame(res,
            index=["p3", "p4", "shared"], 
            columns=["dstat", "bootmean", "bootstd", "abxxa", "baxxa", "Z"]
            )

    return res, boots



def _loci_to_arr(locifile, taxdict, mindict):
    """
    return a frequency array from a loci file for all loci with taxa from 
    taxdict and min coverage from mindict. 
    """

    ## read in the input file
    with open(locifile, 'r') as infile:
        loci = infile.read().strip().split("|\n")
        nloci = len(loci)

    ## make the array (4 or 5) and a mask array to remove loci without cov
    keep = np.zeros(nloci, dtype=np.bool_)
    arr = np.zeros((nloci, 4, 300), dtype=np.float64)
    if len(taxdict) == 5:
        arr = np.zeros((nloci, 6, 300), dtype=np.float64)

    ## if not mindict, make one that requires 1 in each taxon
    if not mindict:
        mindict = {i:1 for i in taxdict}
    if isinstance(mindict, int):
        mindict = {i:mindict for i in taxdict}

    ## raise error if names are not 'p[int]' 
    allowed_names = ['p1', 'p2', 'p3', 'p4', 'p5']
    if any([i not in allowed_names for i in taxdict]):
        raise IPyradError(\
            "keys in taxdict must be named 'p1' through 'p4' or 'p5'")

    ## parse key names
    keys = sorted([i for i in taxdict.keys() if i[0] == 'p'])
    outg = keys[-1]

    ## grab seqs just for the good guys
    for loc in xrange(nloci):    

        ## parse the locus
        lines = loci[loc].split("\n")[:-1]
        names = [i.split()[0] for i in lines]
        seqs = np.array([list(i.split()[1]) for i in lines])

        ## check that names cover the taxdict (still need to check by site)
        covs = [sum([j in names for j in taxdict[tax]]) >= mindict[tax] \
                for tax in taxdict]

        ## keep locus
        if all(covs):
            keep[loc] = True

            ## get the refseq
            refidx = np.where([i in taxdict[outg] for i in names])[0]
            refseq = seqs[refidx].view(np.uint8)
            ancestral = np.array([reftrick(refseq, GETCONS2)[:, 0]])

            ## freq of ref in outgroup
            iseq = _reffreq2(ancestral, refseq, GETCONS2)
            arr[loc, -1, :iseq.shape[1]] = iseq 

            ## enter 4-taxon freqs
            if len(taxdict) == 4:
                for tidx, key in enumerate(keys[:-1]):

                    ## get idx of names in test tax
                    nidx = np.where([i in taxdict[key] for i in names])[0]
                    sidx = seqs[nidx].view(np.uint8)
                   
                    ## get freq of sidx
                    iseq = _reffreq2(ancestral, sidx, GETCONS2)
                   
                    ## fill it in 
                    arr[loc, tidx, :iseq.shape[1]] = iseq

            else:

                ## entere p5; and fill it in
                iseq = _reffreq2(ancestral, refseq, GETCONS2) 
                arr[loc, -1, :iseq.shape[1]] = iseq 
                
                ## enter p1
                nidx = np.where([i in taxdict['p1'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 0, :iseq.shape[1]] = iseq
                
                ## enter p2
                nidx = np.where([i in taxdict['p2'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 1, :iseq.shape[1]] = iseq
                
                ## enter p3 with p4 masked, and p4 with p3 masked
                nidx = np.where([i in taxdict['p3'] for i in names])[0]
                nidy = np.where([i in taxdict['p4'] for i in names])[0]
                sidx = seqs[nidx].view(np.uint8)
                sidy = seqs[nidy].view(np.uint8)
                xseq = _reffreq2(ancestral, sidx, GETCONS2)
                yseq = _reffreq2(ancestral, sidy, GETCONS2)
                mask3 = xseq != 0
                mask4 = yseq != 0
                xseq[mask4] = 0
                yseq[mask3] = 0
                arr[loc, 2, :xseq.shape[1]] = xseq
                arr[loc, 3, :yseq.shape[1]] = yseq
                
                ## enter p34 
                nidx = nidx.tolist() + nidy.tolist()
                sidx = seqs[nidx].view(np.uint8)
                iseq = _reffreq2(ancestral, sidx, GETCONS2)
                arr[loc, 4, :iseq.shape[1]] = iseq


    ## size-down array to the number of loci that have taxa for the test
    arr = arr[keep, :, :]

    ## size-down sites to 
    arr = masknulls(arr)

    return arr, keep



@numba.jit(nopython=True)
def masknulls(arr):
    nvarr = np.zeros(arr.shape[0], dtype=np.int8)
    trimarr = np.zeros(arr.shape, dtype=np.float64)
    for loc in xrange(arr.shape[0]):
        nvars = 0
        for site in xrange(arr.shape[2]):
            col = arr[loc, :, site]
            ## mask cols with 9s
            if not np.any(col == 9):
                ## any non-outgroup shows variation?
                ## todo: check whether BBBBA is ever info?
                if np.any(col[:-1] != col[0]):
                    trimarr[loc, :, nvars] = col
                    nvars += 1
        nvarr[loc] = nvars        
    return trimarr[:, :, :nvarr.max()]



@numba.jit(nopython=True)
def _reffreq2(ancestral, iseq, consdict):
    ## empty arrays
    freq = np.zeros((1, iseq.shape[1]), dtype=np.float64)
    amseq = np.zeros((iseq.shape[0]*2, iseq.shape[1]), dtype=np.uint8)
    
    ## fill in both copies
    for seq in xrange(iseq.shape[0]):
        for col in xrange(iseq.shape[1]):  

            ## get this base and check if it is hetero
            base = iseq[seq][col]
            who = consdict[:, 0] == base
            
            ## if not hetero then enter it
            if not np.any(who):
                amseq[seq*2][col] = base
                amseq[seq*2+1][col] = base        
            ## if hetero then enter the 2 resolutions
            else:
                amseq[seq*2][col] = consdict[who, 1][0]
                amseq[seq*2+1][col] = consdict[who, 2][0]

    ## amseq may have N or -, these need to be masked
    for i in xrange(amseq.shape[1]):
        ## without N or -
        reduced = amseq[:, i][amseq[:, i] != 9]
        counts = reduced != ancestral[0][i]
        if reduced.shape[0]:
            freq[:, i] = counts.sum() / reduced.shape[0]
        else:
            freq[:, i] = 9
    return freq



@numba.jit(nopython=True)
def _prop_dstat(arr):
    
    ## numerator
    abba = ((1.-arr[:, 0]) * (arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))  
    baba = ((arr[:, 0]) * (1.-arr[:, 1]) * (arr[:, 2]) * (1.-arr[:, 3]))
    top = abba - baba
    bot = abba + baba

    ## get statistic and avoid zero div  
    sbot = bot.sum()
    if  sbot != 0:
        dstat = top.sum() / float(sbot)
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



@numba.jit(nopython=True)
def _get_signif_4(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """
    abba, baba, dstat = _prop_dstat(arr)
    boots = _get_boots(arr, nboots)
    e, s = (boots.mean(), boots.std())
    z = 0.
    if s:
        z = np.abs(e) / s
    stats = np.array([dstat, e, s, abba, baba, z])
    return stats, boots



@numba.jit(nopython=True)
def _get_signif_5(arr, nboots):
    """
    returns a list of stats and an array of dstat boots. Stats includes
    z-score and two-sided P-value. 
    """

    statsarr = np.zeros((3, 6), dtype=np.float64)
    bootsarr = np.zeros((3, nboots))

    idx = 0
    for acol in [2, 3, 4]:
        rows = np.array([0, 1, acol, 5])
        tarr = arr[:, rows, :]

        abxa, baxa, dstat = _prop_dstat(tarr)
        boots = _get_boots(tarr, nboots)
        e, s = (boots.mean(), boots.std())
        if s:
            z = np.abs(dstat) / s
        else:
            z = np.NaN
        stats = np.array([dstat, e, s, abxa, baxa, z])

        statsarr[idx] = stats
        bootsarr[idx] = boots
        idx += 1

    return statsarr, bootsarr



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
    wmin = canvas.width * 0.15
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
        ## new data from resarr dataframe
        #res = resarr.iloc[idx]
        boot = bootarr[idx]
        hist = np.histogram(boot, bins=50, range=(rmin, rmax), density=True)
        
        ## get p-value from z-score
        sign = resarr.Z[idx] > cutoff
        if sign:
            if resarr.bootmean[idx] > 0:
                histcolor = toyplot.color.Palette()[0]
            else:
                histcolor = toyplot.color.Palette()[1]
        else:
            histcolor = toyplot.color.Palette()[-1]
        
        ## next axes
        dims = (xlines[idx], xlines[idx]+spacer, hmin, hmax)
        axes = canvas.cartesian(bounds=dims)
        axes.bars(hist, along='y', color=histcolor)

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
    axes.hlines(y=0, style={"stroke-dasharray": "5, 10"})
    axes.show = False

    ## add bar plots
    dims = (xlines[0], xlines[-1], bmin, bmax)
    axes = canvas.cartesian(bounds=dims)
    axes.bars(xlines[1:], resarr.Z, opacity=0.75, color=toyplot.color.Palette()[2])

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
    



def panelplot(tests, resarr, bootsarr, tree):

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
        print(dims)
        
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
    rmax = round(min(1.0, max(0.1, rmax+0.1*rmax)), 1)
    rmin = round(max(-1.0, min(-0.1, -1 * rmax)), 1)

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



class Sim(object):
    def __init__(self, names, sims, nreps, debug):
        self.names = names
        self.sims = sims
        self.nreps = nreps
        self.debug = debug



class Tree(object):
    def __init__(self, newick=None, admix=None):

        ## use default newick string if not given
        if newick:
            self.newick = newick
            self.admix = admix
        else:
            self.newick = "((((a,b),c),d), ((((e,f),g),h) , (((i,j),k),l)));"
        ## parse newick, assigns idx to nodes, returns tre, edges, verts, names
        tree, edges, verts, names = cladogram(self.newick)

        ## parse admixture events
        self.admix = admix
        self._check_admix()

        ## store values
        self.tree = tree
        self.edges = edges
        self.verts = verts
        self.names = names.values()  ## in tree plot vlshow order
        #self.ladderized_tipnames = \
        #    [self.verts[self.tree.search_nodes(name=name)[0].idx, 0]
        #     for name in self.tree.get_leaf_names()]



    def _check_admix(self):
        ## raise an error if admixture event is not possible in time period
        if self.admix:
            for event in self.admix:
                pass #print(event)


    def simulate(self, nreps=1000, admix=None, Ns=int(5e5), gen=10):
        sims = _simulate(self, nreps, admix, Ns, gen)
        debug = 0  ## msprime debug
        names = self.tree.get_leaf_names()[::-1]
        Sims = Sim(names, sims, nreps, debug)
        return Sims


    def draw(
            self, 
            yaxis=False, 
            show_tips=False, 
            use_edge_lengths=True, 
            taxdicts=None, 
            bootsarr=None,
            collapse_outgroup=False,
            test_labels=False,
            **kwargs):
        """
        plot the tree using toyplot.graph. 

        Parameters:
        -----------
            taxdicts: dict
                Show tests as colored rectangles.
            bootsarr: ndarray
                Show bootstrap distributions (requires taxdicts as well)
            yaxis: bool
                Show the y-axis.
            use_edge_lengths: bool
                Use edge lengths from newick tree.
            show_tips: bool
                Show tip names from tree.
            pct_tree_y: float
                proportion of canvas y-axis showing tree
            ...
        """
        ## build toyplot plot
        canvas, axes = _draw(self, yaxis, show_tips, use_edge_lengths, taxdicts, 
                             bootsarr, collapse_outgroup, test_labels, **kwargs)
        return canvas, axes



def _draw(self, yaxis, show_tips, use_edge_lengths, taxdicts, bootsarr, collapse_outgroup, test_labels, **kwargs):

    ## update kwargs from defaults
    args = {"height": min(1000, 15*len(self.tree)),
            "width": min(1000, 15*len(self.tree)),
            "vsize": 0, 
            "vlshow": False, 
            "ewidth": 3, 
            "vlstyle": {"font-size": "18px"}, 
            "cex": "14px", 
            "pct_tree_y": 0.3, 
            "pct_tree_x": 0.7, 
            "lwd_lines": 1,
                }
    ## default if only a tree plot
    if not taxdicts:
        args.update({
            "vsize": 20,
            "vlshow": True,
            })
    args.update(kwargs)

    ## collapse outgroup
    if collapse_outgroup:
        tree, taxdicts = _collapse_outgroup(self.tree, taxdicts)
        newick = tree.write(format=1)
    else:
        newick = self.newick

    ## starting tree position will be changed if adding panel plots
    xmin_tree = 0.
    ymin_tree = 0.
    ymin_test = 0.
    ymin_text = 0.

    ## convert vert to length 1s if not using edges
    verts = copy.deepcopy(self.verts)

    #if not use_edge_lengths:
    tree, edges, verts, names = cladogram(newick, use_edge_lengths=False)
    verts = verts.astype(np.float)

    ## ensure types
    bootsarr = np.array(bootsarr)

    ## relocate tree for panel plots
    if np.any(bootsarr) or taxdicts:
            
        ## adjust Y-axis: boots Y is 125% of tree Y
        pcy = 1 - args["pct_tree_y"]
        newmn = verts[:, 1].max() * pcy
        newhs = np.linspace(newmn, verts[:, 1].max(), len(set(verts[:, 1])))
        ymin_tree += verts[:, 1].max() * pcy
        verts[:, 1] = [newhs[int(i)] for i in verts[:, 1]]

        ## adjust X-axis: boots X is 75% of tree X
        if np.any(bootsarr):
            ## how much do I need to add to make the tree be 60%?
            xmin_tree += verts[:, 0].max() * args["pct_tree_x"]
            verts[:, 0] += verts[:, 0].max() * args["pct_tree_x"]

        ## get spacer between panels
        ztot = 0.15 * xmin_tree

    ## add spacer for tip names
    pctt = 0.2 * (ymin_tree + ymin_test)
    ymin_tree += pctt / 2.
    ymin_test += pctt / 2.
    verts[:, 1] += pctt / 2.

    ## create a canvas and a single cartesian coord system
    canvas = toyplot.Canvas(height=args['height'], width=args['width'])
    axes = canvas.cartesian(bounds=("5%", "95%", "5%", "95%"))
            
    ## add the tree/graph ------------------------------------------------
    _ = axes.graph(edges, 
                    vcoordinates=verts, 
                    ewidth=args["ewidth"], 
                    ecolor=toyplot.color.near_black, 
                    vlshow=args["vlshow"],
                    vsize=args["vsize"],
                    vlstyle=args["vlstyle"],
                    vlabel=names.values())   

    ## add rects for test taxa -------------------------------------------
    if taxdicts:
        yhs = np.linspace(ymin_tree, ymin_test, len(taxdicts)+1)
        ## calculate coords, top and bot 30% of bars is cutoff 
        xmin_hists = 0.
        xmin_zs = (xmin_tree * 0.5) 
        ysp = (yhs[0] - yhs[1])  ## yyy
        ytrim = 0.70 * ysp

        ## colors for tips
        cdict = {"p1": toyplot.color.Palette()[0], 
                 "p2": toyplot.color.Palette()[1],
                 "p3": toyplot.color.near_black, 
                 "p4": toyplot.color.Palette()[-1],}

        ## iterate over tests putting in rects
        for idx, test in enumerate(taxdicts):
            ## check taxdict names
            dictnames = list(itertools.chain(*test.values()))
            badnames = [i for i in dictnames if i not in names.values()]
            if badnames:
                #raise IPyradError(
                print("Warning: found names not in tree:\n -{}"\
                      .format("\n -".join(list(badnames))))

            ## add dashed grid line for tests
            gedges = np.array([[0, 1], [2, 3]])
            gverts = np.array([
                        [xmin_tree, yhs[idx]-ysp/2.],
                        [verts[:, 0].max(), yhs[idx]-ysp/2.],
                        [xmin_zs+ztot/2., yhs[idx]-ysp/2.],
                        [xmin_tree-ztot, yhs[idx]-ysp/2.],
                        ])
            axes.graph(gedges, 
                        vcoordinates=gverts,
                        ewidth=args["lwd_lines"], 
                        ecolor=toyplot.color.Palette()[-1],
                        vsize=0, 
                        vlshow=False,
                        estyle={"stroke-dasharray": "5, 5"},
                        )                              

            ## add test rectangles
            for tax in ["p1", "p2", "p3", "p4"]:
                spx = [tree.search_nodes(name=i)[0] for i in test[tax]]
                spx = np.array([i.x for i in spx], dtype=np.float)
                spx += xmin_tree
                spx.sort()
                ## fill rectangles while leaving holes for missing taxa
                for i in xrange(spx.shape[0]):
                    if i == 0:
                        xleft = spx[i] - 0.25
                        xright = spx[i] + 0.25
                    if i != spx.shape[0]-1:
                        if spx[i+1] - spx[i] < 1.5:
                            xright += 1
                        else:
                            axes.rects(xleft, xright, 
                                yhs[idx] - ytrim, 
                                yhs[idx+1] + ytrim, color=cdict[tax]) 
                            xleft = spx[i+1] - 0.25
                            xright = spx[i+1] + 0.25
                    else:
                        axes.rects(xleft, xright,
                                yhs[idx] - ytrim, 
                                yhs[idx+1] + ytrim, color=cdict[tax]) 
                    
        ## add test-label
        if test_labels:
            if isinstance(test_labels, bool) or test_labels == 1: 
                labels = range(1, len(taxdicts) + 1)
            elif isinstance(test_labels, list):
                labels = test_labels
            else:
                raise IPyradError("  label_tests must be a list or boolean")
            axes.text(
                [verts[:, 0].max() + 1] * len(taxdicts),
                yhs[:-1] - ysp / 2., 
                labels,
                    color=toyplot.color.near_black,
                    style={
                        "font-size": args["cex"],
                        "text-anchor" : "start",
                        "-toyplot-anchor-shift":"0",
                        },
                    )                

        ## add hists
        if np.any(bootsarr):
            ## get bounds on hist
            rmax = np.max(np.abs(bootsarr))
            rmax = round(min(1.0, max(0.2, rmax+0.05*rmax)), 1)
            rmin = round(max(-1.0, min(-0.2, -1 * rmax)), 1)
            allzs = np.abs(np.array([i.mean() / i.std() for i in bootsarr]))
            zmax = max(3., float(np.math.ceil(allzs.max())))

            ## add histograms, and hist axes
            for idx in xrange(bootsarr.shape[0]):
                bins = 30
                mags, xpos = np.histogram(bootsarr[idx], 
                        bins=bins, range=(rmin, rmax), density=True)

                ## get hist colors
                thisz = allzs[idx]
                if thisz > 3:
                    if bootsarr[idx].mean() < 0:
                        color = toyplot.color.Palette()[0]
                    else:
                        color = toyplot.color.Palette()[1]
                else:
                    color = toyplot.color.Palette()[-1]

                ## plot z's within range 
                zright = xmin_tree - ztot
                zleft = xmin_zs + ztot/2.
                zprop = thisz / zmax
                zmaxlen = zright - zleft
                zlen = zmaxlen * zprop
                axes.rects(
                    zright-zlen, zright,
                    yhs[idx] - ytrim, 
                    yhs[idx+1] + ytrim,
                    color=toyplot.color.Palette()[2])

                ## get hist xspans, heights in range
                xran = np.linspace(xmin_hists, xmin_zs - ztot/2., bins+1)
                mags = mags / mags.max()
                mags = (mags * ysp) * 0.8
                yline = yhs[idx+1]
                heights = np.column_stack([[yline for i in mags], mags])
                axes.bars(
                    xran[:-1], 
                    heights, 
                    baseline="stacked",
                    color=["white", color],
                    )

                ## add x-line to histograms
                gverts = np.array([
                    [xran[0], yline], 
                    [xran[-1], yline],
                    ])
                gedges = np.array([[0, 1]])
                axes.graph(
                    gedges, 
                    vcoordinates=gverts, 
                    ewidth=1, 
                    ecolor=toyplot.color.near_black,
                    vsize=0, 
                    vlshow=False,
                    )

            ## add xline to z-plots
            gverts = np.array([
                [xmin_zs + ztot/2., yline], 
                [xmin_tree - ztot, yline], 
                ])
            gedges = np.array([[0, 1]])
            axes.graph(
                gedges, 
                vcoordinates=gverts, 
                ewidth=1, 
                ecolor=toyplot.color.near_black,
                vsize=0, 
                vlshow=False,
                )

            ## add dashed y-line at 0 to histograms
            here = np.where(xpos > 0)[0].min()
            zero = xran[here-1]
            gverts = np.array([
                [zero, ymin_test],
                [zero, ymin_tree - ytrim / 4.],
                ])
            gedges = np.array([[0, 1]])
            axes.graph(
                gedges, 
                vcoordinates=gverts, 
                ewidth=1, 
                ecolor=toyplot.color.near_black,
                eopacity=1.,
                vsize=0, 
                vlshow=False,
                estyle={"stroke-dasharray": "5, 5"},
                )

            ## add solid y-line to z-plots
            gverts = np.array([
                [xmin_tree - ztot, ymin_test + ytrim/4.],
                [xmin_tree - ztot, ymin_tree - ytrim/4.],
                ])
            gedges = np.array([[0, 1]])
            axes.graph(
                gedges, 
                vcoordinates=gverts, 
                ewidth=1, 
                ecolor=toyplot.color.near_black,
                eopacity=1.,
                vsize=0, 
                vlshow=False,
                )

            ## add tick-marks to x-lines (hists and z-plots)
            ticklen = ytrim / 4.
            gedges = np.array([[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]])
            gverts = np.array([
                [xmin_hists, ymin_test - ticklen],
                [xmin_hists, ymin_test],
                [zero, ymin_test - ticklen],                              
                [zero, ymin_test],
                [xmin_zs - ztot / 2., ymin_test - ticklen],
                [xmin_zs - ztot / 2., ymin_test],
                [xmin_zs + ztot / 2., ymin_test - ticklen],
                [xmin_zs + ztot / 2., ymin_test],
                [xmin_tree - ztot, ymin_test - ticklen],
                [xmin_tree - ztot, ymin_test],
                ])
            axes.graph(
                gedges, 
                vcoordinates=gverts, 
                ewidth=1, 
                ecolor=toyplot.color.near_black,
                eopacity=1.,
                vsize=0, 
                vlshow=False,
                )

            ## add tick labels
            labels = [rmin, 0, rmax, zmax, 0]
            axes.text(
                gverts[:, 0][::2], 
                [ymin_test - ysp] * len(labels),
                labels, 
                color=toyplot.color.near_black,
                style={
                    "font-size": args["cex"],
                    "text-anchor" : "middle",
                    "-toyplot-anchor-shift":"0",
                    },
                )

            ## add baba abba labels
            axes.text(
                [gverts[:, 0][0], gverts[:, 0][4]],
                [ymin_test - ysp * 2] * 2,
                ["BABA", "ABBA"], 
                color=toyplot.color.near_black,
                style={
                    "font-size": args["cex"], #"12px",
                    "text-anchor" : "middle",
                    "-toyplot-anchor-shift":"0",
                    },
                )     

            ## add bootstrap and z-score titles
            axes.text(
                [zero, zright - 0.5 * zmaxlen],
                [ymin_tree + ysp / 2.] * 2,
                ["Bootstrap D-statistics", "Z-scores"], 
                color=toyplot.color.near_black,
                style={
                    "font-size": args["cex"], #"12px",
                    "text-anchor" : "middle",
                    "-toyplot-anchor-shift":"0",
                    },
                )

    ## add names to tips --------------------------------------------------
    if show_tips:
        ## calculate coords
        nams = [i for i in names.values() if not isinstance(i, int)]
        spx = [tree.search_nodes(name=i)[0] for i in nams]
        spx = np.array([i.x for i in spx], dtype=np.float)
        spx += xmin_tree
        if taxdicts:
            spy = [yhs[-1] - ysp / 2.] * len(nams)
        else:
            spy = [ymin_test - 0.5] * len(nams)
        _ = axes.text(spx, spy, nams,
                        angle=-90, 
                        color=toyplot.color.near_black,
                        style={
                            "font-size": args["cex"],
                            "text-anchor" : "start",
                            "-toyplot-anchor-shift":"0",
                            },
                        ) 

    ## plot admix lines ---------------------------------
    if self.admix:
        for event in self.admix:
            ## get event
            source, sink, _, _, _ = event

            ## get nodes from tree
            source = self.tree.search_nodes(name=source)[0]
            sink = self.tree.search_nodes(name=sink)[0]

            ## get coordinates
            fromx = np.max([source.up.x, source.x]) - np.abs(source.up.x - source.x) / 2.
            fromy = source.y + (source.up.y - source.y) / 2.
            tox = np.max([sink.up.x, sink.x]) - np.abs(sink.up.x - sink.x) / 2.
            toy = sink.y + (sink.up.y - sink.y) / 2.
                
            ## if show_tips:
            if show_tips:
                fromy += spacer
                toy += spacer

            ## plot
            mark = axes.plot([fromx, tox], [fromy, toy], 
                            color=toyplot.color.Palette()[1], 
                            style={"stroke-width": 3, 
                                   "stroke-dasharray": "2, 2"},
                            )
                
    ## hide x and hide/show y axies
    axes.x.show = False
    if yaxis:
        axes.y.show = True
    else:
        axes.y.show = False    

    ## return plotting 
    return canvas, axes




def _collapse_outgroup(tree, taxdicts):
    """ collapse outgroup in ete Tree for easier viewing """
    ## check that all tests have the same outgroup
    outg = taxdicts[0]["p4"]
    if not all([i["p4"] == outg for i in taxdicts]):
        raise Exception("no good")
   
    ## prune tree, keep only one sample from outgroup
    tre = ete.Tree(tree.write(format=1)) #tree.copy(method="deepcopy")
    alltax = [i for i in tre.get_leaf_names() if i not in outg]
    alltax += [outg[0]]
    tre.prune(alltax)
    tre.search_nodes(name=outg[0])[0].name = "outgroup"
    tre.ladderize()

    ## remove other ougroups from taxdicts
    taxd = copy.deepcopy(taxdicts)
    newtaxdicts = []
    for test in taxd:
        #test["p4"] = [outg[0]]
        test["p4"] = ["outgroup"]
        newtaxdicts.append(test)

    return tre, newtaxdicts




## convertes newick to (edges, vertices)
def cladogram(newick, use_edge_lengths=True, invert=False):
    
    ## invert and short name to arg so it is similar to ape
    ig = use_edge_lengths == False

    ## get tree
    tre = ete.Tree(newick=newick)
    tre.ladderize()

    ## map numeric values to internal nodes from root to tips
    ## preorder: first parent and then children. These indices will
    ## be used in the int edge array to map edges.
    names = {}
    idx = 0
    for node in tre.traverse("preorder"):
        if not node.is_leaf():
            if node.name:
                names[idx] = node.name
            else:
                names[idx] = idx
                node.name = str(idx)
            node.idx = idx
            idx += 1

    ## map number to the tips, these will be the highest numbers
    for node in sorted(tre.get_leaves(), key=lambda x: x.name):
        names[idx] = node.name
        node.idx = idx
        idx += 1

    ## create empty edges and coords arrays
    edges = np.zeros((idx-1, 2), dtype=int)
    verts = np.zeros((idx, 2), dtype=float)
    
    ## postorder: first children and then parents. This moves up the list .
    nidx = 0
    tip_num = len(tre.get_leaves()) - 1
    for node in tre.traverse("postorder"):
        #for nidx in range(idx)[::-1]:
        #node = tre.search_nodes(idx=nidx)[0]
        if node.is_leaf():
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = tip_num 
            node.y = 0
            tip_num -= 1
            verts[node.idx] = [node.x, node.y]
        
        elif node.is_root():
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y]
        
        else:
            edges[nidx-1, :] = node.up.idx, node.idx
            node.x = sum(i.x for i in node.children) / 2.
            if ig:
                node.y = node.get_farthest_leaf(ig)[1] + 1
            else:
                node.y = node.get_farthest_leaf()[1]
            verts[node.idx] = [node.x, node.y] 
        nidx += 1

    ## invert for sideways trees
    if invert:
        verts[:, 1] = np.abs(verts[:, 1] - tlen)

    return tre, edges, verts, names



######################################################################
## Simulation functions (requires msprime)
######################################################################

def _simulate(self, nreps, admix=None, Ns=500000, gen=20):
    """
    Enter a baba.Tree object in which the 'tree' attribute (newick 
    derived tree) has edge lengths in units of generations. You can 
    use the 'gen' parameter to multiply branch lengths by a constant. 

    Parameters:
    -----------

    nreps: (int)
        Number of reps (loci) to simulate under the demographic scenario
    tree: (baba.Tree object)
        A baba.Tree object initialized by calling baba.Tree(*args). 
    admix: (list)
        A list of admixture events to occur on the tree. Nodes must be 
        reference by their index number, and events must occur in time
        intervals when edges exist. Use the .draw() function of the 
        baba.Tree object to see node index numbers and coalescent times.
    Ns: (float)
        Fixed effective population size for all lineages (may allow to vary
        in the future). 
    gen: (int)
        A multiplier applied to branch lengths to scale into units of 
        generations. Example, if all edges on a tree were 1 then you might
        enter 50000 to multiply so that edges are 50K generations long.

    """

    ## node ages
    Taus = np.array(list(set(self.verts[:, 1]))) * 1e4 * gen

    ## The tips samples, ordered alphanumerically
    ## Population IDs correspond to their indexes in pop config
    ntips = len(self.tree)
    #names = {name: idx for idx, name in enumerate(sorted(self.tree.get_leaf_names()))}
    ## rev ladderized leaf name order (left to right on downward facing tree)
    names = {name: idx for idx, name in enumerate(self.tree.get_leaf_names()[::-1])}
    pop_config = [
        ms.PopulationConfiguration(sample_size=2, initial_size=Ns)
        for i in range(ntips)
    ]

    ## migration matrix all zeros init
    migmat = np.zeros((ntips, ntips)).tolist()

    ## a list for storing demographic events
    demog = []

    ## coalescent times
    coals = sorted(list(set(self.verts[:, 1])))[1:]
    for ct in xrange(len(coals)):
        ## check for admix event before next coalescence
        ## ...
        
        ## print coals[ct], nidxs, time
        nidxs = np.where(self.verts[:, 1] == coals[ct])[0]
        time = Taus[ct+1]

        ## add coalescence at each node
        for nidx in nidxs:
            node = self.tree.search_nodes(name=str(nidx))[0]

            ## get destionation (lowest child idx number), and other
            dest = sorted(node.get_leaves(), key=lambda x: x.idx)[0]
            otherchild = [i for i in node.children if not 
                          i.get_leaves_by_name(dest.name)][0]

            ## get source
            if otherchild.is_leaf():
                source = otherchild
            else:
                source = sorted(otherchild.get_leaves(), key=lambda x: x.idx)[0]
            
            ## add coal events
            event = ms.MassMigration(
                        time=int(time),
                        source=names[source.name], 
                        destination=names[dest.name], 
                        proportion=1.0)
            #print(int(time), names[source.name], names[dest.name])
        
            ## ...
            demog.append(event)
            
            
    ## sim the data
    replicates = ms.simulate(
        population_configurations=pop_config,
        migration_matrix=migmat,
        demographic_events=demog,
        num_replicates=nreps,
        length=100, 
        mutation_rate=1e-8)
    return replicates



## simulates data on 12 taxon tree with two admixture events
def _sim_admix_12(nreps, Ns=500000, gen=20):
    
    # Set the ML values of various parameters
    Taus = np.array([0, 1, 2, 3, 4, 5]) * 1e4 * gen

    # Migration rates C -> B and from IJ -> EF
    m_C_B = 2e-6
    m_IJ_EF = 2e-6
    
    # Population IDs correspond to their indexes in pop_config.
    ntips = len(tree.tree)
    pop_config = [
        ms.PopulationConfiguration(sample_size=2, initial_size=Ns)
        for i in range(ntips)]
    
    ## migration matrix all zeros time=0
    migmat = np.zeros((ntips, ntips)).tolist()
    
    ## set up demography
    demog = [
        ## initial migration from C -> B
        ms.MigrationRateChange(time=0, rate=m_C_B, matrix_index=(1, 2)),
        ms.MigrationRateChange(time=Taus[1], rate=0),

        # merge events at time 1 (b,a), (f,e), (j,i)
        ms.MassMigration(time=Taus[1], source=1, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[1], source=5, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[1], source=9, destination=8, proportion=1.0), 
        
        ## migration from IJ -> EF (backward in time)
        ms.MigrationRateChange(time=Taus[1], rate=m_IJ_EF, matrix_index=(4, 8)), 

        ## merge events at time 2 (c,a), (g,e), (k,i)
        ms.MassMigration(time=Taus[2], source=2, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[2], source=6, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[2], source=10, destination=8, proportion=1.0), 

        ## end migration at ABC and merge
        ms.MigrationRateChange(time=Taus[2], rate=0),
        ms.MassMigration(time=Taus[3], source=3, destination=0, proportion=1.0), 
        ms.MassMigration(time=Taus[3], source=7, destination=4, proportion=1.0), 
        ms.MassMigration(time=Taus[3], source=11, destination=8, proportion=1.0),   
        
        ## merge EFJH -> IJKL
        ms.MassMigration(time=Taus[4], source=8, destination=4, proportion=1.0),   
        
        ## merge ABCD -> EFJHIJKL
        ms.MassMigration(time=Taus[5], source=4, destination=0, proportion=1.0),   
    ]

    ## sim the data
    replicates = ms.simulate(
        population_configurations=pop_config,
        migration_matrix=migmat,
        demographic_events=demog,
        num_replicates=nreps,
        length=100, 
        mutation_rate=1e-9)
    
    return replicates


#def _msp_to_arr(simreps, test):
def _msp_to_arr(Sim, test):
    
    ## the fixed tree dictionary
    #fix = {j: [i, i+1] for j, i in zip(list("abcdefghijkl"), range(0, 24, 2))}
    fix = {j: [i, i+1] for j, i in zip(Sim.names, range(0, len(Sim.names)*2, 2))}
    
    ## fill taxdict by test
    keys = ['p1', 'p2', 'p3', 'p4']
    arr = np.zeros((Sim.nreps, 4, 100))
    
    ## unless it's a 5-taxon test
    if len(test) == 5:
        arr = np.zeros((100000, 6, 100))
        keys += ['p5']
    
    ## create array sampler for taxa
    taxs = [test[key] for key in keys]
    idxs = [list(itertools.chain(*[fix[j] for j in i])) for i in taxs]

    ## iterate over reps filling arr
    idx = 0
    for trees in simreps:
        
        ## build genotype array
        shape = trees.get_num_mutations(), trees.get_sample_size()
        garr = np.empty(shape, dtype="u1")
    
        ## fill the garr
        for variant in trees.variants():
            garr[variant.index] = variant.genotypes
        
        if len(test) == 4:
            if garr.shape[0]:
                ## fill my arr with freqs
                for pdx, tax in enumerate(idxs):
                    freq = garr[:, tax]
                    freq = freq.sum(axis=1) / float(freq.shape[1])
                    maxsz = min(freq.shape[0], 100)
                    arr[idx, pdx, :maxsz] = freq[:maxsz]
        else:
            if garr.shape[0]:
                ## get the easy ones
                p1 = garr[:, idxs[0]]
                p2 = garr[:, idxs[1]]
                p5 = garr[:, idxs[4]]
                p34 = garr[:, idxs[2]+idxs[3]]

                ## identity of SNPs is important
                p3 = garr[:, idxs[2]]
                p4 = garr[:, idxs[3]]
                
                ## any rows with data in b are masked in a
                mask3 = np.where(p3.sum(axis=1) == 0)[0]
                mask4 = np.where(p4.sum(axis=1) == 0)[0]
                masked_p3 = p3[mask4]
                masked_p4 = p4[mask3]
                
                ## enter frequencies
                freq = p1
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 0, :maxsz] = freq[:maxsz]
                
                freq = p2
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 1, :maxsz] = freq[:maxsz]
               
                freq = masked_p3
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 2, :maxsz] = freq[:maxsz]               
               
                freq = masked_p4
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 3, :maxsz] = freq[:maxsz]
               
                freq = p34
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 4, :maxsz] = freq[:maxsz]
                 
                freq = p5
                freq = freq.sum(axis=1) / float(freq.shape[1])
                maxsz = min(freq.shape[0], 100)
                arr[idx, 5, :maxsz] = freq[:maxsz]
        idx += 1

    ## reduce the size of arr to min loci        
    arr = arr[:idx+1]
    
    ## reduce the size of arr to min len
    minl = np.where(np.all(np.all(arr==0, axis=1) == True, axis=0))[0]
    if np.any(minl):
        minl = minl.min()
    else:
        minl = None
    arr = arr[:, :, :minl]
    
    return arr



## combines sim_admix12 + msp_to_arr + baba to return single (stats, boots)
def sim_admix_12_baba(nreps, test, mindict, nboots):
    sims = _sim_admix_12(nreps)
    arr = _msp_to_arr(sims, test)
    stats, boots = baba(arr, test, mindict, nboots)
    return stats, boots




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
