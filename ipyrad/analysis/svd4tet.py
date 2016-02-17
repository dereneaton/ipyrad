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

from __future__ import print_function, division
import os
import h5py
import random
import itertools
import numpy as np
from ipyrad.assemble.util import *
from collections import Counter, OrderedDict


import logging
LOGGER = logging.getLogger(__name__)



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



def get_seqarray(data, path=None):
    """ Takes a data object connected to an ipyclient and looks for a 
    phylip file. Or should I use another format?
    """

    ## allow path override of object
    if path:
        phylip = open(path, 'r')
    else:
        ## turn the phylip file into an array that can be indexed
        phylip = open(data.outfiles.phylip, 'r')
        line = phylip.readline().strip().split()
        ntax = int(line[0])
        nbp = int(line[1])

    ## make a seq array
    seqarray = np.zeros((ntax, nbp), dtype="S1")

    with h5py.File(data._svd.path, 'w') as io5:
        io5.create_dataset("seqarray", (ntax, nbp), dtype="S1")        

        for line, seq in enumerate(phylip.readlines()):
            seqarray[line] = np.array(list(seq.split()[-1]))

        ## save array to disk for so it can be easily accessed from 
        ## many engines on arbitrary nodes 
        io5["seqarray"][:] = seqarray
        del seqarray

    return data



def get_quartets(data, samples):
    """ load all quartets into an array """

    ## calculate how many quartets to generate
    qiter = itertools.combinations(range(len(samples)), 4)
    nquarts = sum(1 for _ in qiter)

    ## load an array of the right size
    ## first 4 are the sample indices, then the quartet, and then the weight.
    if nquarts < 1000:
        chunk = nquarts
    else:
        chunk = 1000

    print("populating array with {} quartets".format(nquarts))
    with h5py.File(data._svd.path, 'a') as io5:
        io5.create_dataset("samples", (nquarts, 4), 
                           dtype=np.int16, chunks=(chunk, 4))
        io5.create_dataset("quartets", (nquarts, 4), 
                            dtype=np.int16, chunks=(chunk, 4))        
        io5.create_dataset("weights", (nquarts, 1), 
                            dtype=np.float16, chunks=(chunk, 1))

        
        ## populate array with all possible quartets. This allows us to 
        ## sample from the total, and also to continue from a checkpoint
        qiter = itertools.combinations(range(len(samples)), 4)
        i = 0
        while i < nquarts:
            dat = np.array(list(itertools.islice(qiter, 1000)))
            io5["samples"][i:i+dat.shape[0]] = dat
            i += 1000

    ## set to data object
    data._svd.nquarts = nquarts
    data._svd.chunk = chunk

    return data



## FROM THE ITERTOOLS RECIPES COOKCOOK
def random_combination(iterable, nquarts):
    """
    Random selection from itertools.combinations(iterable, r). 
    Use this if not sampling all possible quartets.
    """
    pool = tuple(iterable)
    size = len(pool)
    indices = sorted(random.sample(xrange(size), nquarts))
    return tuple(pool[i] for i in indices)



def seq_to_matrix(arr, sidx):
    """ 
    Takes a 4-taxon alignment and generates a 16x16 count matrix. Ignores 
    ambiguous sites, sites with missing data, and sites with indels. The 
    order of the samples (sidx) is important. 

    MDICT is an orderedict with the following keys that will be reshaped to 
    be 16x16: 


    """
    ## select only the relevant rows
    arr = arr[sidx, :]
    
    ## convert all nulls to Ns
    for null in list("RSKYWM-"):
        arr[arr == null] = "N"

    ## mask any columns with Ns
    arr = arr[:, ~np.any(arr == "N", axis=0)]

    ## get dict of observed patterns
    counts = tablestack(arr)

    ## fill mdict with ordered patterns
    mdict = OrderedDict([(i, 0) for i in MKEYS.split()])
    for patt in counts:
        mdict[patt] = counts[patt]

    ## shape into 16x16 array
    matrix = np.array(mdict.values()).reshape(16, 16)
    return matrix



def tablestack(arr):
    """ makes a count dict of each unique array element """
    ## goes by 10% at a time to minimize memory overhead. Is possible it skips
    ## the last chunk, but this shouldn't matter.
    table = Counter()
    for i in xrange(0, arr.shape[1], arr.shape[1]//10):
        tmp = Counter([j.tostring() for j in arr[:, i:i+arr.shape[1]//10].T])
        table.update(tmp)
    return table



def worker(args):
    """ gets data from seq array and puts results into results array """

    ## parse args
    data, qidx = args

    ## get 1000 quartets from qidx
    io5 = h5py.File(data._svd.path, 'r+')
    quartets = io5["samples"][qidx:qidx+1000]
    seqarray = io5["seqarray"][:]

    ## get data from seqarray
    results = [svd4tet(seqarray, qtet) for qtet in quartets]

    ## put results into h5
    io5["quartets"][qidx:qidx+1000] = results
    io5.close()



def svd4tet(arr, samples):
    """ calc rank. From Chibatko and Kiffman (2014)

    Our proposal for inferring the true species-level relationship within a 
    sample of four taxa is thus the following. For each of the three possible 
    splits, construct the matrix FlatL1|L2(P) and compute SVD (L1|L2). The 
    split with the smallest score is taken to be the true split. 

    """
    ## get the three resolution of four taxa
    split0 = [samples[0], samples[1], samples[2], samples[3]]
    split1 = [samples[0], samples[2], samples[1], samples[3]]
    split2 = [samples[0], samples[3], samples[1], samples[2]]
    splits = [split0, split1, split2]

    ## get the three matrices
    mats = [seq_to_matrix(arr, sidx) for sidx in splits]

    ## calculate which is closest to a rank 10 matrix
    scores = [np.linalg.svd(m, full_matrices=1, compute_uv=0) for m in mats]
    score = [np.sqrt(sco[11:].sum()) for sco in scores]

    ## return the split with the lowest score
    idx = score.index(min(score))
    return splits[idx]



def dump():
    """ prints the quartets to a file formatted for qmaxcut"""

    outfile = open("test.svd.leaves", 'w')

    for i in range(0, 1000, 10):
        chunk = [list(i) for i in io5["quartets"][i:i+10]]
        chunk = ["{},{}|{},{}".format(*i) for i in chunk]
        outfile.write(" ".join(chunk)+" ")
    outfile.close()




def main(data, ipyclient, path=None, checkpoint=0):
    """ 
    main funcs
    """
    ## store some svd stats in data
    data._svd = ObjDict()
    data._svd.path = os.path.join(data.dirs.outfiles, data.name+"_svd.h5")
    data._svd.checkpoint = 0

    ## check for checkpoint restart
    if checkpoint:
        ## scan the hdf5 array to find where it left off writing
        raise IPyradWarningExit("not yet implemented")

    ## get the seq array into hdf5
    print("loading array")
    data = get_seqarray(data, path)

    ## make quartet arrays into hdf5. Allow subsetting samples eventually.
    ## and optimize chunk value given remaining quartets and ipyclient    
    print("loading quartets")
    samples = data.samples
    data = get_quartets(data, samples)

    ## load a job queue
    jobs = []
    for qidx in xrange(data._svd.checkpoint, nquarts, 1000):
        jobs.append([data, qidx])

    ## make a distributor for engines
    lbview = ipyclient.load_balanced_view()
    results = lbview.map_async(worker, jobs)
    results.get()

    print("  done")
    return data



## step 1 create an h5 array with all possible quartets






