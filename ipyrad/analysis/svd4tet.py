#!/usr/bin/env ipython2

""" 
SVD-quartet like tree inference. Modelled on the following papers:

Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under 
the coalescent, Bioinformatics, 30(23): 3317-3324.

Chifman, J. and L. Kubatko. 2015. Identifiability of the unrooted species 
tree topology under the coalescent model with time-reversible substitution 
processes, site-specific rate variation, and invariable sites, Journal of 
Theoretical Biology 374: 35-47

"""

# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0103
# pylint: disable=C0301



from __future__ import print_function, division
import os
import sys
import glob
import h5py
import time
import random
import ipyrad
import numba
import datetime
import itertools
import subprocess
import numpy as np
import pandas as pd
import ipyparallel as ipp
from fractions import Fraction

from ipyrad.assemble.util import ObjDict, IPyradWarningExit, progressbar
#from collections import Counter, OrderedDict

## ete3 is an extra dependency not included with ipyrad
try:
    import ete3
except ImportError:
    try:
        import ete2 as ete3
    except ImportError:
        raise IPyradWarningExit("""
    svd4tet requires the dependency `ete3`. You can install
    it with the command `conda install -c etetoolkit ete3`
    """)

## set the logger
import logging
LOGGER = logging.getLogger(__name__)

## debug numba code
#numba.NUMBA_DISABLE_JIT = 1

## The 16 x 16 matrix of site counts. 
MKEYS = """
    AAAA AAAC AAAG AAAT  AACA AACC AACG AACT  AAGA AAGC AAGG AAGT  AATA AATC AATG AATT
    ACAA ACAC ACAG ACAT  ACCA ACCC ACCG ACCT  ACGA ACGC ACGG ACGT  ACTA ACTC ACTG ACTT
    AGAA AGAC AGAG AGAT  AGCA AGCC AGCG AGCT  AGGA AGGC AGGG AGGT  AGTA AGTC AGTG AGTT
    ATAA ATAC ATAG ATAT  ATCA ATCC ATCG ATCT  ATGA ATGC ATGG ATGT  ATTA ATTC ATTG ATTT

    CAAA CAAC CAAG CAAT  CACA CACC CACG CACT  CAGA CAGC CAGG CAGT  CATA CATC CATG CATT
    CCAA CCAC CCAG CCAT  CCCA CCCC CCCG CCCT  CCGA CCGC CCGG CCGT  CCTA CCTC CCTG CCTT
    CGAA CGAC CGAG CGAT  CGCA CGCC CGCG CGCT  CGGA CGGC CGGG CGGT  CGTA CGTC CGTG CGTT
    CTAA CTAC CTAG CTAT  CTCA CTCC CTCG CTCT  CTGA CTGC CTGG CTGT  CTTA CTTC CTTG CTTT

    GAAA GAAC GAAG GAAT  GACA GACC GACG GACT  GAGA GAGC GAGG GAGT  GATA GATC GATG GATT
    GCAA GCAC GCAG GCAT  GCCA GCCC GCCG GCCT  GCGA GCGC GCGG GCGT  GCTA GCTC GCTG GCTT
    GGAA GGAC GGAG GGAT  GGCA GGCC GGCG GGCT  GGGA GGGC GGGG GGGT  GGTA GGTC GGTG GGTT
    GTAA GTAC GTAG GTAT  GTCA GTCC GTCG GTCT  GTGA GTGC GTGG GTGT  GTTA GTTC GTTG GTTT

    TAAA TAAC TAAG TAAT  TACA TACC TACG TACT  TAGA TAGC TAGG TAGT  TATA TATC TATG TATT
    TCAA TCAC TCAG TCAT  TCCA TCCC TCCG TCCT  TCGA TCGC TCGG TCGT  TCTA TCTC TCTG TCTT
    TGAA TGAC TGAG TGAT  TGCA TGCC TGCG TGCT  TGGA TGGC TGGG TGGT  TGTA TGTC TGTG TGTT
    TTAA TTAC TTAG TTAT  TTCA TTCC TTCG TTCT  TTGA TTGC TTGG TTGT  TTTA TTTC TTTG TTTT
"""



def get_seqarray(data, boot):
    """ 
    Takes an Assembly object and looks for a phylip file unless the path 
    argument is used then it looks in path. If boots then it grabs the 
    resampled seqarray, else it samples the original seqarray. 
    """

    ## read in the file
    spath = open(data.outfiles.svdinput, 'r')
    line = spath.readline().strip().split()
    ntax = int(line[0])
    nbp = int(line[1])
    #LOGGER.info("array shape: (%s, %s)", ntax, nbp)

    ## make a tmp seq array
    tmpseq = np.zeros((ntax, nbp), dtype="S1")

    ## get initial array
    if not boot:
        print("  loading array [{} taxa x {} bp]".format(ntax, nbp))        
        ## use 'w' to make initial array
        with h5py.File(data.svd.h5in, 'w') as io5:
            ## create array storage for both real seq and for later bootstraps
            ## seqdata will be stored as int8
            seqarr = io5.create_dataset("seqarr", (ntax, nbp), dtype=np.uint8)
            io5.create_dataset("bootarr", (ntax, nbp), dtype=np.uint8)
            ## fill the tmp array from the input phy
            for line, seq in enumerate(spath.readlines()):
                tmpseq[line] = np.array(list(seq.split()[-1]))
            ## save array to disk so it can be easily accessed by slicing
            seqarr[:] = tmpseq.view(np.uint8)
            del tmpseq
            
    else:
        ## use 'r+' to read and write to existing array
        with h5py.File(data.svd.h5in, 'r+') as io5:        
            ## load in the seqarr
            tmpseq = pd.DataFrame(io5["seqarr"][:])
            ## fill the boot array with a re-sampled phy w/ replacement
            io5["bootarr"][:] = np.uint8(tmpseq.sample(n=tmpseq.shape[1], 
                                                       replace=True, axis=1))
            del tmpseq
    return data



## FROM THE ITERTOOLS RECIPES COOKCOOK
def random_combination(iterable, nquarts):
    """
    Random selection from itertools.combinations(iterable, r). 
    Use this if not sampling all possible quartets.
    """
    pool = tuple(iterable)
    size = len(pool)
    indices = random.sample(xrange(size), nquarts)
    return tuple(pool[i] for i in indices)



def random_product(iter1, iter2):
    """ random sampler for equa_splits func"""
    pool1 = tuple(iter1)
    pool2 = tuple(iter2)
    ind1 = random.sample(pool1, 2)
    ind2 = random.sample(pool2, 2)
    return tuple(ind1+ind2)



MUL = lambda x, y: x*y

def n_choose_k(n, k):
    """ calculate the number of quartets as n-choose-k. This is used
    in equal splits to decide whether a split should be exhaustively sampled
    or randomly sampled. Edges near tips can be exhaustive while highly nested
    edges probably have too many quartets
    """
    return int(reduce(MUL, (Fraction(n-i, i+1) for i in range(k)), 1))



def equal_splits(data, nquarts, ipyclient):
    """ sample quartets even across splits of the starting tree """
    ## choose chunker for h5 arr
    chunk = (nquarts // (2 * len(ipyclient))) + (nquarts % (2 * len(ipyclient)))
    LOGGER.info("E: nquarts = %s, chunk = %s", nquarts, chunk)

    ## get starting tree
    tre = ete3.Tree(".tmptre")
    tre.unroot()
    print("  starting tree:")
    print(tre)

    ## randomly sample all splits of tree
    splits = [([z.name for z in i], 
               [z.name for z in j]) \
               for (i, j) in tre.get_edges()]
    ## only keep internal splits
    splits = [i for i in splits if all([len(j) > 1 for j in i])]
    N = len(data.samples)
    if len(splits) < ((N * (N-1)) // 2):
        print("  starting tree is unresolved, sample more quartets")

    ## turn each into an iterable split sampler
    ## if the nquartets for that split is small, then sample all of them
    ## if it is UUUUUGE, then make it a random sampler from that split
    qiters = []
    ## how many quartets are we gonna sample from each quartet?
    squarts = nquarts // len(splits)

    for split in splits:
        ## if small then sample all
        if n_choose_k(len(split[0]), 2) * n_choose_k(len(split[1]), 2) < squarts*2:
            qiter = itertools.product(
                        itertools.combinations(split[0], 2), 
                        itertools.combinations(split[1], 2))
        ## or create random sampler
        else:
            qiter = (random_product(split[0], split[1]) for _ in xrange(nquarts))
        qiters.append(qiter)

    ## make qiters infinitely cycling
    qiters = itertools.cycle(qiters)

    ## iterate over qiters sampling from each, if one runs out, keep 
    ## sampling from remaining qiters. Keep going until samples is filled
    with h5py.File(data.svd.h5in, 'a') as io5:
        ## create data sets
        del io5["samples"]
        del io5["quartets"]
        smps = io5.create_dataset("samples", (nquarts, 4), 
                           dtype=np.uint16, chunks=(chunk, 4),
                           compression='gzip')
        io5.create_dataset("quartets", (nquarts, 4), 
                            dtype=np.uint16, chunks=(chunk, 4),
                            compression='gzip')

        ## fill 1000 at a time for efficiency
        i = 0
        while i < nquarts:
            qdat = []
            while len(qdat) < min(1000, smps.shape[0]): 
                qiter = qiters.next()
                try:
                    qdat.append(qiter.next())
                except StopIteration:
                    print(len(qdat))
                    continue
            dat = np.array(qdat, dtype=np.uint16)
            smps[i:i+dat.shape[0]] = dat
            i += 1000

    ## save attrs to data object
    data.svd.nquarts = nquarts
    data.svd.chunk = chunk

    return data




def get_quartets(data, method, nquarts, ipyclient):
    """ 
    Find all quartets of samples and store in a large array 
    """

    ## calculate how many quartets to generate
    if method == 'all':
        nquarts = n_choose_k(len(data.samples), 4)
        #qiter = itertools.combinations(range(len(data.samples)), 4)
        #nquarts = sum(1 for _ in qiter)

    ## create a chunk size for sampling from the array of quartets. This should
    ## be relatively large so that we don't spend a lot of time doing I/O, but
    ## small enough that jobs finish every few hours since that is how the 
    ## checkpointing works.
    breaks = 2
    if nquarts < 5000:
        breaks = 1
    if nquarts > 100000:
        breaks = 4
    if nquarts > 500000:
        breaks = 8

    chunk = (nquarts // (breaks * len(ipyclient))) + \
            (nquarts % (breaks * len(ipyclient)))
    #LOGGER.info("nquarts = %s, chunk = %s", nquarts, chunk)

    ## 'samples' stores the indices of the quartet. 
    ## `quartets` stores the correct quartet in the order (1,2|3,4)
    ## `weights` stores the calculated weight of the quartet in 'quartets'
    ## we gzip this for now, but check later if this has a big speed cost
    print("  populating array with {} quartets".format(nquarts))
    with h5py.File(data.svd.h5in, 'a') as io5:
        ## create data sets
        smps = io5.create_dataset("samples", (nquarts, 4), 
                           dtype=np.uint16, chunks=(chunk, 4),
                           compression='gzip')
        io5.create_dataset("quartets", (nquarts, 4), 
                            dtype=np.uint16, chunks=(chunk, 4),
                            compression='gzip')

        ## populate array with all possible quartets. This allows us to 
        ## sample from the total, and also to continue from a checkpoint
        qiter = itertools.combinations(range(len(data.samples)), 4)
        i = 0

        ## fill 1000 at a time for efficiency
        while i < nquarts:
            if method != "all":
                qiter = []
                while len(qiter) < min(1000, smps.shape[0]):
                    qiter.append(
                        random_combination(range(len(data.samples)), 4))
                dat = np.array(qiter)
            else:
                dat = np.array(list(itertools.islice(qiter, 1000)))
            smps[i:i+dat.shape[0]] = dat
            i += 1000

    ## save attrs to data object
    data.svd.nquarts = nquarts
    data.svd.chunk = chunk

    return data



#############################################################################
#############################################################################
## SVD computation functions
## 
#############################################################################
#############################################################################

@numba.jit('f4(f4[:])', nopython=True)
def get_weights(scores):
    """ 
    calculates quartet weights from ordered svd scores. Following 
    description from Avni et al. 
    """
    scores.sort()
    scores = scores[::-1]
    if scores[2]:
        weight = np.float32((scores[2]-scores[0]) / 
                            (np.exp(scores[2]-scores[1]) * scores[2]))
    else:
        #weight = np.float32(0.00001)
        weight = scores[2]
    
    if np.isnan(weight):
        weight = np.float32(0.01)
        
    return weight



@numba.jit('Tuple((f4,f4,f4))(u4[:,:])', nopython=True)
def abba_baba(mat):
    """ 
    calculate dstats from the count array and return as a float tuple 
    """
    ## get all the abba sites from mat
    baba = 0
    for i in range(16):
        if i % 5:
            baba += mat[i, i]
    baba = np.float32(baba)
    
    ## get all the baba sites from mat
    abba = np.float32(\
            mat[1, 4] + mat[2, 8] + mat[3, 12] +\
            mat[4, 1] + mat[6, 9] + mat[7, 13] +\
            mat[8, 2] + mat[9, 6] + mat[11, 14] +\
            mat[12, 3] + mat[13, 7] + mat[14, 11])
     
    ## calculate D, protect from ZeroDivision
    denom = abba + baba
    if denom:
        dstat = (abba-baba)/denom
    else:
        dstat = 0
    
    return abba, baba, dstat

        

@numba.jit('u1[:,:](u1[:,:],b1[:])', nopython=True)
def removeNs(seqchunk, mask):
    """ removes ncolumns from snparray prior to matrix calculation"""
    ## mask columns that contain Ns
    for idx in xrange(seqchunk.shape[1]):
        if np.sum(seqchunk[:, idx] == 78):
            mask[idx] = False
    ## apply mask
    newarr = seqchunk[:, mask]
    ## return smaller Nmasked array
    return newarr



#@numba.jit()
def reduce_arr(arr, map, sidx):
    """ 
    If a full seqarray and map file are provided then the array is reduced 
    to a single unlinked snp from each locus for this quartet.
    """
    pass



@numba.jit('u4[:,:](u1[:,:])', nopython=True)
def chunk_to_matrix(narr):
    """ 
    numba compiled code to get matrix fast.
    arr is a 4 x N seq matrix converted to np.int8
    I convert the numbers for ATGC into their respective index for the MAT
    matrix, and leave all others as high numbers, i.e., -==45, N==78. 
    """

    ## get seq alignment and create an empty array for filling
    mat = np.zeros((16, 16), dtype=np.uint32)

    ## replace ints with small ints that index their place in the 
    ## 16x16. If not replaced, the existing ints are all very large
    ## and the column will be excluded.
    for x in xrange(narr.shape[1]):
        i = narr[:, x]
        ## convert to index values
        i[i == 65] = 0  ## already zero
        i[i == 67] = 1
        i[i == 71] = 2
        i[i == 84] = 3

        if np.sum(i) < 16:
            mat[i[0]*4:(i[0]+4)*4]\
               [i[1]]\
               [i[2]*4:(i[2]+4)*4]\
               [i[3]] += 1
    return mat



#@numba.jit('(u2[:,:],u1[:,:],u1[:,:])')#, nopython=True)
#@numba.jit(nopython=True)
def nworker(smpchunk, seqview, tests):
    """ The workhorse function. All numba. """

    ## get the input arrays ready
    rquartets = np.zeros((smpchunk.shape[0], 4), dtype=np.uint16)
    rweights = np.zeros(smpchunk.shape[0], dtype=np.float32)
    mask = np.ones(seqview.shape[1], dtype=np.bool_)
    rdstats = np.zeros((smpchunk.shape[0], 3), dtype=np.float32)

    ## fill arrays with results using numba funcs
    for idx in xrange(smpchunk.shape[0]):
        ## get seqchunk for 4 samples (4, ncols) and remove Ns
        seqchunk = seqview[smpchunk[idx]]
        seqnon = removeNs(seqchunk, mask)

        ## get svdscores for each arrangement of seqchunk
        qscores = np.zeros(3, dtype=np.float32)
        mats = np.zeros((3, 16, 16), dtype=np.uint32)
        for test in range(3):
            mats[test] = chunk_to_matrix(seqnon[tests[test]])
            ## get svd scores
            tmpscore = np.linalg.svd(mats[test].astype(np.float32))[1]
            qscores[test] = np.sqrt(tmpscore[11:]).sum()

        ## sort to find the best qorder
        best = np.where(qscores == qscores.min())[0]
        bidx = tests[best][0]
        LOGGER.info("""
            best: %s, 
            bidx: %s, 
            qscores: %s, 
            mats %s
            """, best, bidx, qscores, mats)
        rquartets[idx] = smpchunk[idx][bidx]#[0]

        ## get weights from the three scores sorted
        rweights[idx] = get_weights(qscores)

        ## get dstat from the best (correct) matrix 
        ## (or should we get all three?) [option]
        rdstats[idx] = abba_baba(mats[best][0])

    #return 
    return rquartets, rweights, rdstats 


#############################################################################
#############################################################################
## End svd funcs
#############################################################################
#############################################################################


def run_qmc(data, boot):
    """ runs quartet max-cut """

    ## make call lists
    cmd1 = " ".join(
            [ipyrad.bins.qmc,
            " qrtt="+data.svd.qdump, 
            " weights=off"+
            " otre=.tmptre"])

    cmd2 = " ".join(
            [ipyrad.bins.qmc,
            " qrtt="+data.svd.qdump, 
            " weights=on"+
            " otre=.tmpwtre"])

    ## run them
    for cmd in [cmd1, cmd2]:
        #LOGGER.info(cmd)
        try:
            subprocess.check_call(cmd, shell=True,
                                       stderr=subprocess.STDOUT,
                                       stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as inst:
            LOGGER.error("Error in wQMC: \n({}).".format(inst))
            LOGGER.error(subprocess.STDOUT)
            raise inst

    ## read in the tmp files since qmc does not pipe
    tmptre = open(".tmptre").read().strip()
    tmpwtre = open(".tmpwtre").read().strip()    

    ## convert int names back to str names
    tmptre = ete3.Tree(tmptre)
    tmptre = renamer(data, tmptre)
    tmpwtre = ete3.Tree(tmpwtre)
    tmpwtre = renamer(data, tmpwtre)

    ## save the boot tree
    if boot:
        with open(data.svd.tboots, 'a') as outboot:
            outboot.write(tmptre+"\n")
        with open(data.svd.wboots, 'a') as outboot:
            outboot.write(tmpwtre+"\n")

    ## store full data trees to Assembly
    else:
        with open(data.svd.tre, 'w') as outtree:
            outtree.write(tmptre)
        with open(data.svd.wtre, 'w') as outtree:
            outtree.write(tmpwtre)

    ## save assembly with new tree
    data.save()



def renamer(data, tre):
    """ renames newick from numbers to sample names"""
    ## order the numbered tree tip lables
    names = tre.get_leaves()
    names.sort(key=lambda x: int(x.name))

    ## order the sample labels in the same order they are 
    ## in the seq file (e.g., .snp, .phy)
    snames = data.samples.keys()
    snames.sort()

    ## replace numbered names with snames
    for (tip, sname) in zip(names, snames):
        tip.name = sname

    ## return with only topology and leaf labels
    return tre.write(format=9)



def dump(data):
    """ 
    prints the quartets to a file formatted for wQMC 
    """
    ## open the h5 database
    io5 = h5py.File(data.svd.h5out, 'r')

    ## create an output file for writing
    data.svd.qdump = os.path.join(data.dirs.svd, data.name+"_quartets.txt")
    outfile = open(data.svd.qdump, 'w')

    ## todo: should pull quarts order in randomly
    for idx in range(0, data.svd.nquarts, data.svd.chunk):
        quarts = [list(j) for j in io5["quartets"][idx:idx+data.svd.chunk]]
        weight = io5["weights"][idx:idx+data.svd.chunk]
        chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
                                                in zip(quarts, weight)]
        outfile.write("\n".join(chunk)+"\n")

    ## close output file and h5 database
    outfile.close()
    io5.close()



def insert_to_array(data, start, results):
    """
    inputs results from workers into hdf4 array
    """
    qrts, wgts, dsts = results

    with h5py.File(data.svd.h5out, 'r+') as out:
        chunk = data.svd.chunk
        out['quartets'][start:start+chunk] = qrts
        out['weights'][start:start+chunk] = wgts
        out["dstats"][start:start+chunk] = dsts
    ## save checkpoint
    #data.svd.checkpoint_arr = np.where(ww == 0)[0].min()



#def insert_to_array(data, result):
#     """ 
#     Takes a tmpfile output from finished worker enters it into the 
#     full h5 array, and deletes the tmpfile
#     """
#     out5 = h5py.File(data.svd.h5out, 'r+')

#     with h5py.File(result) as inh5:
#         qrts = inh5['quartets'][:]
#         wgts = inh5['weights'][:]

#         chunk = data.svd.chunk
#         start = int(result.split("_")[-1][:-3])
#         out5['quartets'][start:start+chunk] = qrts
#         out5['weights'][start:start+chunk] = wgts
#     out5.close()

#     ## remove tmp h5path
#     os.remove(result)



def svd_obj_init(data, method):
    """ creates svd attribute to Assembly object """
    ## create ObjDict for storage
    data.svd = ObjDict()

    ## add array path
    data.svd.method = method
    data.svd.h5in = os.path.join(data.dirs.svd, data.name+"_input.h5")
    data.svd.h5out = os.path.join(data.dirs.svd, data.name+"_output.h5")

    ## original tree paths
    data.svd.tre = os.path.join(data.dirs.svd, data.name+"_svd4tet.tre")
    data.svd.wtre = os.path.join(data.dirs.svd, data.name+"_svd4tet.w.tre")

    ## bootstrap tree paths
    data.svd.tboots = os.path.join(data.dirs.svd, data.name+"_svd4tet.boots")
    data.svd.wboots = os.path.join(data.dirs.svd, data.name+"_svd4tet.w.boots")

    ## bootstrap labeled o.g. trees paths
    data.svd.btre = os.path.join(data.dirs.svd, data.name+"_svd4tet.support.tre")
    data.svd.bwtre = os.path.join(data.dirs.svd, data.name+"_svd4tet.w.support.tre")

    ## NHX formatted tre with rich information
    data.svd.nhx = os.path.join(data.dirs.svd, data.name+"_svd4tet.nhx")
    data.svd.wnhx = os.path.join(data.dirs.svd, data.name+"_svd4tet.w.nhx")

    ## checkpoints
    data.svd.checkpoint_boot = 0
    data.svd.checkpoint_arr = 0

    ## save to object and return
    data.save()
    return data



def svd4tet(data, nboots=100, method="all", nquarts=None, force=False):
    """ 
    API wrapper for svd4tet analysis

    data 
        ipyrad Assembly object
    nboots
        number of non-parametric bootstrap replicates to run
    method
        all, random, or equal. Default is all, which samples all possible
        quartets. For very large trees (>50 tips) this may take too long, 
        in which case you should use random or equal. The arg nquarts 
        determines how many quartets will be samples. In random, nquarts
        are sampled and used. In equal, a starting tree is inferred and
        the random quartets are drawn so that they are spread ~equally
        across splits of the tree. 
    nquarts 
        The numer of random quartets sampled in random or equal method. 
        Default is 10000, or all if all < 10000. 
    force
        Overwrite existing
    """

    ## check that method was entered correctly
    assert method in ["all", "random", "equal"], \
        "method type not recognized, must be one of ['all', 'random', 'equal']"

    if method != "all":
        ## require nquarts if method not all
        assert nquarts, "if method != all, must enter a value for nquarts"
        ## don't allow nquarts to be greater than all
        totalquarts = n_choose_k(len(data.samples), 4)
        if nquarts > totalquarts:
            print("  nquarts > total quartets, switching to method='all'")
            method = "all"
        if nquarts < 500:
            print("  few possible quartets, only method='all' available")
            method = "all"

    ## launch ipclient, assumes ipyparallel is running
    ipyclient = ipp.Client(timeout=10)

    ## protects it from KBD
    try:
        run(data, nboots, method, nquarts, force, ipyclient)

    except (KeyboardInterrupt, SystemExit):

        ## protect from KBD while saving
        try:
            ## cancel submitted jobs
            #ipyclient.abort()
            ## kill running jobs
            #ipyclient.close()

            ## remove any abandoned tmp arrays 
            abandon = glob.glob(os.path.join(data.dirs.svd, "*_tmp_*.h5"))
            for afile in abandon:
                os.remove(afile)

        except KeyboardInterrupt:
            pass

    finally:
        ## checkpoint the state and save
        LOGGER.info("\n  saving checkpoints to [Assembly].svd")
        LOGGER.info("  array checkpoint: %s", data.svd.checkpoint_arr)
        LOGGER.info("  boot checkpoint: %s", data.svd.checkpoint_boot)
        data.save()        



def run(data, nboots, method, nquarts, force, ipyclient):
    """ 
    Run svd4tet inference on a sequence or SNP alignment for all samples in
    the Assembly. 

    By default the job starts from 0 or where it last left off, unless 
    force=True, then it starts from 0. 
    """

    ## load svd attributes if they exist
    fresh = 0
    if not force:
        try:
            if data.svd.checkpoint_boot or data.svd.checkpoint_arr:
                print("  loading from svd checkpoint")
                print("  array checkpoint: {}".format(data.svd.checkpoint_arr))
                print("  boots checkpoint: {}".format(data.svd.checkpoint_boot))
                print("  sampling method: {}".format(data.svd.method))
                ## require method to be same as loaded type
                assert method == data.svd.method, \
                    "loaded object method={}, cannot change methods midstream"+\
                    " use force argument to start new run with new method."
            else:
                fresh = 1
        except (AttributeError, IOError):
            fresh = 1

    ## if svd results do not exist or force then restart
    if force or fresh:
        ## make an analysis directory if it doesn't exist
        data.dirs.svd = os.path.realpath(
                            os.path.join(
                                data.dirs.project, data.name+"_analysis_svd"))
        if not os.path.exists(data.dirs.svd):
            try:
                os.mkdir(data.dirs.svd)
            except OSError:
                ## if not there then create new svd directory
                data.dirs.svd = os.path.join(
                                    os.path.curdir, data.name+"_analysis_svd")
                os.mkdir(data.dirs.svd)
                print("  output directory created at: {}".format(data.dirs.svd))

        ## init new svd ObjDict
        data = svd_obj_init(data, method)

        ## get the real seq array into hdf5 h5in
        data = get_seqarray(data, boot=False)

        ## make quartet arrays into hdf5. Allow subsetting samples eventually.
        ## and optimize chunk value given remaining quartets and ipyclient    
        if method == "equal":
            ## print equal header
            print("  loading {} random quartet samples for starting tree inference"\
                  .format(nquarts))
            ## grab test number for starting tree
            data = get_quartets(data, method, nquarts, ipyclient)
            print("  inferring {} x 3 quartet trees for starting tree"\
                  .format(nquarts))
            ## infer starting tree
            inference(data, ipyclient, bidx=0)
            ## sample quartets from starting tree
            print("  loading {} equal-splits quartets from starting tree"\
                  .format(nquarts))            
            data = equal_splits(data, nquarts, ipyclient)
            ## remove starting tree tmp files
            tmps = [data.svd.tre, data.svd.wtre, data.svd.tboots, 
                    data.svd.wboots, data.svd.btre, data.svd.bwtre]
            for tmp in tmps:
                try:
                    os.remove(tmp)
                except OSError:
                    continue

        ## will sample all or random set of quartets    
        else:
            if method == "random":
                print("  loading {} random quartet samples"\
                      .format(nquarts))
            else:   
                nquarts = n_choose_k(len(data.samples), 4)                
                print("  loading all {} possible quartets"\
                      .format(nquarts))
            data = get_quartets(data, method, nquarts, ipyclient)


    ## run the full inference 
    if not data.svd.checkpoint_boot:
        print("  inferring {} x 3 quartet trees".format(nquarts))
        inference(data, ipyclient, bidx=0)
    else:
        print("  full inference finished")
        elapsed = datetime.timedelta(seconds=int(0)) #time.time()-start))
        progressbar(20, 20, " | {}".format(elapsed))


    ## run the bootstrap replicates
    if nboots:
        start = time.time()
        print("  running {} bootstrap replicates".format(nboots))
    
        ## get current boot
        for bidx in range(data.svd.checkpoint_boot, nboots):
        
            if data.svd.checkpoint_arr == 0:
                data = get_seqarray(data, boot=True)
                #LOGGER.info("  new boot array sampled")
                data.svd.checkpoint_boot = bidx
            ## start boot inference
            elapsed = datetime.timedelta(seconds=int(time.time()-start))            
            progressbar(nboots, bidx, " | {}".format(elapsed))
            inference(data, ipyclient, bidx=True)
        elapsed = datetime.timedelta(seconds=int(time.time()-start))            
        progressbar(20, 20, " | {}".format(elapsed))

        ## write outputs with bootstraps
        write_outputs(data, with_boots=1)

    else:
        ## write outputs without bootstraps
        write_outputs(data, with_boots=0)



def write_outputs(data, with_boots=1):
    """ write final tree files """

    ## create tree with support values on edges
    write_supports(data, with_boots)

    ## make tree figure
    data.svd.quicktre = data.svd.nhx.rsplit(".", 1)[0]+"_quartet_sampling.pdf"

    ## hide error message during tree plotting
    #save_stdout = sys.stdout           
    ## file-like obj to catch stdout
    #sys.stdout = cStringIO.StringIO()  
    ## run func with stdout hidden
    #ipyclient = ipp.Client(**args)
    ## resets stdout
    #sys.stdout = save_stdout    
    quickfig(data.svd.nhx, data.svd.quicktre)

    ## print finished
    print("\n  Final quartet-joined and weighted quartet-joined tree files:\
          \n  {}\n  {}".format(
          data.svd.tre, 
          data.svd.wtre))

    if with_boots:
        print("\n  Bootstrap tree files:\
              \n  {}\n  {}".format(
              data.svd.tboots, 
              data.svd.wboots))

        print("\n  Final trees with bootstraps as edge labels: \n  {}\n  {}".format(
              data.svd.btre, 
              data.svd.bwtre))

    print("\n  Final trees with edge info in NHX format: \n  {}\n  {}".format(
          data.svd.nhx,
          data.svd.wnhx))

    qtre = ete3.Tree(data.svd.tre)
    qtre.unroot()
    print("\n  Quick view of unrooted topology from unweighted analysis")
    print(qtre)

    print("\n  Quartet sampling visualized on a tree: \n  {}".format(
          data.svd.quicktre))
    print("  ")

    docslink = "ipyrad.readthedocs.org/cookbook.html"    
    print("  * For tips on plotting trees see: {}".format(docslink))

    citelink = "ipyrad.readthedocs.org/svd4tet.html"
    print("  * For tips on citing this software see: {}".format(citelink))
    print("")

    return data



PIECOLORS = ['#a6cee3',
             '#1f78b4']

def layout(node):
    """ layout for ete3 tree plotting fig """
    if node.is_leaf():
        nameF = ete3.TextFace(node.name, tight_text=False, 
                                         fgcolor="#262626", fsize=8)
        ete3.add_face_to_node(nameF, node, column=0, position="aligned")
        node.img_style["size"] = 0
              
    else:
        if not node.is_root():
            node.img_style["size"] = 0
            node.img_style["shape"] = 'square'
            node.img_style["fgcolor"] = "#262626"   
            if "quartets_total" in node.features:
                ete3.add_face_to_node(ete3.PieChartFace(
                                    percents=[float(node.quartets_sampled_prop), 
                                    100-float(node.quartets_sampled_prop)],
                                    width=15, height=15, colors=PIECOLORS),
                                    node, column=0, position="float-behind")  
            if "bootstrap" in node.features:
                ete3.add_face_to_node(ete3.AttrFace(
                                  "bootstrap", fsize=7, fgcolor='red'),
                                  node, column=0, position="branch-top")  
        else:
            node.img_style["size"] = 0
            if node.is_root():
                node.img_style["size"] = 0
                node.img_style["shape"] = 'square'
                node.img_style["fgcolor"] = "#262626"   



def quickfig(input_tre, outname):
    """ make a quick ete3 fig. Plots total quartets """
    ts = ete3.TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.mode = 'r'
    ts.draw_guiding_lines = True
    ts.show_scale = False
    ts.scale = 25

    tre = ete3.Tree(input_tre)
    tre.ladderize()
    tre.convert_to_ultrametric(tree_length=len(tre)//2)
    tre.render(file_name=outname, h=40*len(tre), tree_style=ts)



def get_total(tre, node):
    """ get total number of quartets possible for a split """
    tots = set(tre.get_leaf_names())
    down = set(node.get_leaf_names())
    up = tots - down
    return n_choose_k(len(down), 2) * n_choose_k(len(up), 2)
    

    
def get_sampled(data, tre, node, names):
    """ get how many quartets were sampled that are informative for a split"""
    ## get leaves up and down
    tots = set(tre.get_leaf_names())
    down = set(node.get_leaf_names())
    up = tots - down

    ## get up and down as index
    idxd = set([names.index(i) for i in down])
    idxu = set([names.index(i) for i in up])

    sampled = 0
    with h5py.File(data.svd.h5out, 'r') as io5:
        qrts = io5["quartets"][:]
        ## iterate over quartets 
        for qrt in qrts:
            sqrt = set(qrt)
            if len(sqrt.intersection(idxd)) > 1:
                if len(sqrt.intersection(idxu)) > 1:
                    sampled += 1

    return sampled



def write_supports(data, with_boots):
    """ writes support values as edge labels on unrooted tree """
    ## get name indices
    names = data.samples.keys()
    names.sort()

    ## get unrooted best trees
    otre = ete3.Tree(data.svd.tre, format=0)
    otre.unroot()
    for node in otre.traverse():
        node.add_feature("bootstrap", 0)
        node.add_feature("quartets_total", get_total(otre, node))
        node.add_feature("quartets_sampled", get_sampled(data, otre, node, names))
        try:
            prop = 100*(float(node.quartets_sampled) / node.quartets_total)
        except ZeroDivisionError:
            prop = 0.0
        node.add_feature("quartets_sampled_prop", prop)
        node.dist = 0
        node.support = 0

    wtre = ete3.Tree(data.svd.wtre, format=0)
    wtre.unroot()
    for node in wtre.traverse():
        node.add_feature("bootstrap", 0)
        node.add_feature("quartets_total", get_total(wtre, node))
        node.add_feature("quartets_sampled", get_sampled(data, wtre, node, names))
        try:
            prop = 100*(float(node.quartets_sampled) / node.quartets_total)
        except ZeroDivisionError:
            prop = 0.0
        node.add_feature("quartets_sampled_prop", prop)
        node.dist = 0
        node.support = 0

    ## get unrooted boot trees
    if with_boots:
        oboots = open(data.svd.tboots, 'r').readlines()
        wboots = open(data.svd.wboots, 'r').readlines()
        oboots = [ete3.Tree(btre.strip()) for btre in oboots]
        wboots = [ete3.Tree(btre.strip()) for btre in wboots]    
        _ = [btre.unroot() for btre in oboots]
        _ = [btre.unroot() for btre in wboots]

        ## get and set support values 
        for tre, boots in zip([otre, wtre], [oboots, wboots]):
            for btre in boots:
                common = tre.compare(btre, unrooted=True)
                for bnode in common["common_edges"]:
                    ## check monophyly of each side of split
                    a = tre.check_monophyly(bnode[0], target_attr='name', unrooted=True)
                    b = tre.check_monophyly(bnode[1], target_attr='name', unrooted=True)
                    ## if both sides are monophyletic
                    if a[0] and b[0]:
                        ## find which is the 'bottom' node, to attach support to
                        node = list(tre.get_monophyletic(bnode[0], target_attr='name'))
                        node.extend(list(tre.get_monophyletic(bnode[1], target_attr='name')))
                        ## add +1 suport to (edge dist) to this edge
                        if not node[0].is_leaf():
                            node[0].dist += 1
                            node[0].support += 1
                            node[0].bootstrap += 1

        ## change support values to percentage
        for tre in [otre, wtre]:
            for node in tre.traverse():
                node.dist = int(100 * (node.dist / len(wboots)))
                node.support = int(100 * (node.support / len(wboots)))
                node.bootstrap = int(100 * (node.bootstrap / len(wboots)))

        ## return as newick string w/ support as edge labels (lengths)
        with open(data.svd.btre, 'w') as outtre:
            outtre.write(otre.write(format=5))

        with open(data.svd.bwtre, 'w') as outtre:
            outtre.write(wtre.write(format=5))

        features = ["bootstrap", "quartets_total", "quartets_sampled", "quartets_sampled_prop"]            
                    
    else:
        features = ["quartets_total", "quartets_sampled", "quartets_sampled_prop"]

    ## return as NHX format with extra info
    with open(data.svd.nhx, 'w') as outtre:
        outtre.write(wtre.write(format=0, features=features))

    with open(data.svd.wnhx, 'w') as outtre:
        outtre.write(wtre.write(format=0, features=features))

    return otre, wtre



def inference(data, ipyclient, bidx):
    """ run inference and store results """

    ## how many chunks are we gonna be doing?
    njobs = sum(1 for _ in iter(xrange(
                data.svd.checkpoint_arr, 
                data.svd.nquarts, 
                data.svd.chunk)))

    ## an iterator to distribute those chunks
    #jobiter = iter(xrange(data.svd.checkpoint_arr, 
    #                      data.svd.nquarts, 
    #                      data.svd.chunk))
    jobiter = iter(xrange(0, data.svd.nquarts, 100))
    njobs = sum(1 for _ in jobiter)
    jobiter = iter(xrange(0, data.svd.nquarts, 100))    
    data.svd.chunk = 100


    LOGGER.info("chunksize: %s, start: %s, total: %s, njobs: %s", \
            data.svd.chunk, data.svd.checkpoint_arr, data.svd.nquarts, njobs)

    ## a distributor for engine jobs
    lbview = ipyclient.load_balanced_view()

    ## create h5 storage data sets for the correct quartets & and their weights
    with h5py.File(data.svd.h5out, 'w') as io5:
        io5.create_dataset("quartets", (data.svd.nquarts, 4), 
                            dtype=np.uint16, chunks=(data.svd.chunk, 4))
        io5.create_dataset("weights", (data.svd.nquarts,), 
                            dtype=np.float32, chunks=(data.svd.chunk,))
        io5.create_dataset("dstats", (data.svd.nquarts, 3), 
                            dtype=np.float32, chunks=(data.svd.chunk, 3))

    ## A view to the inh5 hdf5 array will be left open since each time we 
    ## submit a job we will be slicing a chunk of quartets from disk. 
    ## This is also where the seqarray is, but we assume that can be held
    ## in memory just fine (this could change!). One view of the array is fine
    ## since ipyparallel does NOT copy arrays in memory when sending to engines.
    inh5 = h5py.File(os.path.realpath(data.svd.h5in), 'r')
    if not data.svd.checkpoint_boot:
        seqs = inh5["seqarr"][:]
    else:
        seqs = inh5["bootarr"][:]

    ## the three indexed resolutions of each quartet
    tests = np.array([[0, 1, 2, 3], 
                      [0, 2, 1, 3], 
                      [0, 3, 1, 2]], dtype=np.uint8)

    ## start progress bar timer and submit initial n jobs
    start = time.time()
    res = {}
    for _ in xrange(min(len(ipyclient), njobs)):
        ## get chunk of quartet samples
        qidx = jobiter.next()
        smps = inh5["samples"][qidx:qidx+data.svd.chunk]
        print(qidx, smps.shape, seqs.shape)
        ## send chunk off to be worked on
        res[qidx] = lbview.apply(nworker, *[smps, seqs, tests])

    #ipyclient.wait_interactive()
    done = 0
    while 1:
        ## print progress unless bootstrapping, diff progbar for that.
        if not bidx:
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            progressbar(njobs, done, " | {}".format(elapsed))

        ## check for finished jobs
        curkeys = res.keys()
        finished = [i.ready() for i in res.values()]

        ## remove finished and submit new jobs
        if any(finished):
            for ikey in curkeys:
                if res[ikey].completed:
                    if res[ikey].successful():
                        LOGGER.info("cleanup key %s", ikey)
                        ## track finished
                        done += 1
                        ## insert results into hdf5 data base
                        results = res[ikey].get(0)
                        insert_to_array(data, ikey, results)
                        ## purge memory of the old one
                        del res[ikey]
                    else:
                        ## print error if something went wrong
                        meta = res[ikey].metadata
                        if meta.error:
                            LOGGER.error("""\
                        stdout: %s
                        stderr: %s 
                        error: %s""", meta.stdout, meta.stderr, meta.error)
                        del res[ikey]

                ## submit new jobs
                try:
                    qidx = jobiter.next()
                    smps = inh5["samples"][qidx:qidx+data.svd.chunk]
                    ## send chunk off to be worked on
                    res[qidx] = lbview.apply(nworker, *[smps, seqs, tests])
                    print(qidx, smps.shape, seqs.shape)

                except StopIteration:
                    continue
        else:
            time.sleep(0.1)

        ## know when to finish
        if njobs == done:
            break

    ## final progress bar
    if not bidx:
        progressbar(njobs, done, " | {}".format(elapsed))
        print("")

    ## convert to txt file for wQMC
    dump(data)    
    sys.exit(2)

    ## run quartet joining algorithm
    if not bidx:
        run_qmc(data, boot=0)
    else:
        run_qmc(data, boot=1)

    ## reset the checkpoint_arr
    data.svd.checkpoint_arr = 0



if __name__ == "__main__":

    ## imports
    import ipyrad.analysis as ipa
    #import ipyrad as ip
    #import ipyparallel as ipp

    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/cli/cli.json")
    DATA = ipyrad.load_json("~/Documents/ipyrad/tests/iptutorial/cli.json")
    ## run
    ipa.svd4tet.wrapper(DATA, nboots=10, method='equal', nquarts=50, force=True)

