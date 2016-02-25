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
# pylint: disable=C0103

from __future__ import print_function, division
import os
import sys
import h5py
import time
import random
import ipyrad
import itertools
import subprocess
import numpy as np
import ipyparallel as ipp

from ipyrad.assemble.util import *
from collections import Counter, OrderedDict

import logging
LOGGER = logging.getLogger(__name__)


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



def get_seqarray(data, path=None, snps=False):
    """ 
    Takes an Assembly object and looks for a phylip file unless the path 
    argument is used then it looks in path. 
    """

    ## allow path override of object
    if path:
        spath = open(path, 'r')
    else:
        ## turn the phylip file into an array that can be indexed
        if snps:
            spath = open(data.outfiles.snp, 'r')
        else:
            spath = open(data.outfiles.phy, 'r')            
    line = spath.readline().strip().split()
    ntax = int(line[0])
    nbp = int(line[1])

    ## make a seq array
    tmpseq = np.zeros((ntax, nbp), dtype="S1")
    with h5py.File(data._svd.h5, 'w') as io5:
        seqarray = io5.create_dataset("seqarray", (ntax, nbp), dtype="S1")        

        for line, seq in enumerate(spath.readlines()):
            tmpseq[line] = np.array(list(seq.split()[-1]))

        ## save array to disk for so it can be easily accessed from 
        ## many engines on arbitrary nodes 
        seqarray[:] = tmpseq
        del tmpseq

    return data



def get_quartets(data, samples):
    """ load all quartets into an array """

    ## calculate how many quartets to generate
    qiter = itertools.combinations(range(len(samples)), 4)
    nquarts = sum(1 for _ in qiter)

    ## create a chunk size for sampling from the array of quartets. This should
    ## be relatively large so that we don't spend a lot of time doing I/O, but
    ## small enough that jobs finish every few hours since that is how the 
    ## checkpointing works.
    if nquarts < 1000:
        chunk = nquarts // 20
    else:
        chunk = 1000

    ## 'samples' stores the indices of the quartet. 
    ## `quartets` stores the correct quartet in the order (1,2|3,4)
    ## `weights` stores the calculated weight of the quartet in 'quartets'
    print("  populating array with {} quartets".format(nquarts))
    with h5py.File(data._svd.h5, 'a') as io5:
        ## create data sets
        smps = io5.create_dataset("samples", (nquarts, 4), 
                           dtype=np.uint16, chunks=(chunk, 4))
        #qrts = io5.create_dataset("quartets", (nquarts, 4), 
        #                    dtype=np.uint16, chunks=(chunk, 4))        
        #wgts = io5.create_dataset("weights", (nquarts,), 
        #                    dtype=np.float16, chunks=(chunk,))
        
        ## populate array with all possible quartets. This allows us to 
        ## sample from the total, and also to continue from a checkpoint
        qiter = itertools.combinations(range(len(samples)), 4)
        i = 0
        ## fill 1000 at a time for efficiency
        while i < nquarts:
            dat = np.array(list(itertools.islice(qiter, 1000)))
            smps[i:i+dat.shape[0]] = dat
            i += 1000

    ## save attrs to data object
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
    """
    ## select only the relevant rows
    arr = arr[sidx, :]
    
    ## convert all nulls to Ns
    for null in list("RSKYWM-"):
        arr[arr == null] = "N"

    ## mask any columns for this 4-taxon alignment that has Ns
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
    chunk = data._svd.chunk

    ## get 1000 quartets from qidx
    lpath = os.path.realpath(data._svd.h5)
    with h5py.File(lpath, 'r') as io5in:
        smps = io5in["samples"][qidx:qidx+chunk]
        seqs = io5in["seqarray"][:]

        ## get data from seqarray
        results = [svd4tet(seqs, qtet) for qtet in smps]
        rquartets = np.array([i[0] for i in results])
        rweights = np.array([i[1] for i in results])
        #print("rquartets", rquartets, rquartets.shape)
        #print("rweights", rweights, rweights.shape)

    tmpchunk = os.path.join(data.dirs.svd, data.name+"_tmp_{}.h5".format(qidx))
    with h5py.File(tmpchunk, 'w') as io5out:
        #io5out["weights"][qidx:qidx+chunk] = np.array(rweights)
        #io5out["quartets"][qidx:qidx+chunk] = np.array(rquartets)
        io5out["weights"] = np.array(rweights)
        io5out["quartets"] = np.array(rquartets)

    return tmpchunk
    #with h5py.File(lpath, 'r') as io5:
    #    print(io5["quartets"][:])



def svd4tet(arr, samples, weights=0):
    """ calc rank. From Chibatko and Kiffman (2014)

    Our proposal for inferring the true species-level relationship within a 
    sample of four taxa is thus the following. For each of the three possible 
    splits, construct the matrix FlatL1|L2(P) and compute SVD (L1|L2). The 
    split with the smallest score is taken to be the true split. 

    """
    ## get the three resolution of four taxa
    splits = [[samples[0], samples[1], samples[2], samples[3]],
              [samples[0], samples[2], samples[1], samples[3]],
              [samples[0], samples[3], samples[1], samples[2]]]

    ## get the three matrices
    mats = [seq_to_matrix(arr, sidx) for sidx in splits]

    ## calculate which is closest to a rank 10 matrix
    scores = [np.linalg.svd(m, full_matrices=1, compute_uv=0) for m in mats]
    score = [np.sqrt(sco[11:].sum()) for sco in scores]

    ## get splits sorted in the order of score
    sortres = [(x, y) for (x, y) in sorted(zip(score, splits))]
    rsplits = [i[1] for i in sortres]
    rweight = [i[0] for i in sortres]

    ## calculate weights for quartets (Avni et al. 2014). Here I use svd-scores
    ## for the weights rather than genetic distances, but could do either.
    weight = 0
    if weights:
        weight = get_weight_svd(rweight)
        #weight = get_weight_snir(arr, rsplits)

    ## return the split with the lowest score
    return rsplits[0], weight


def get_weight_svd(ssplits):
    """ 
    Caulcate weight similar to Avni et al. but using SVD score distance from 
    rank 10 instead of JC distance. Experimental...
    """
    dl, dm, dh = ssplits
    print(dl, dm, dh)

    calc = (dh-dl) / (np.exp(dh-dm) * dh)
    print(calc)
    return calc



def get_distances(arr, splits):
    """ 
    Calculate JC distance from the 4-taxon shared sequence data. This func is
    not used. Finish writing it if you end up using gen. distances.
    """

    ## distances given (0,1|2,3)
    d_ab = np.bool(arr[splits[0], :] != arr[splits[1], :]).sum() / arr.shape[1]
    d_cd = np.bool(arr[splits[2], :] != arr[splits[3], :]).sum() / arr.shape[1]
    dist0 = d_ab + d_cd

    ## distances given (0,2|1,3)
    d_ac = np.bool(arr[splits[0], :] != arr[splits[2], :]).sum() / arr.shape[1]
    d_bd = np.bool(arr[splits[1], :] != arr[splits[3], :]).sum() / arr.shape[1]
    dist1 = d_ac + d_bd

    ## distances given (0,3|1,2)
    d_ad = np.bool(arr[splits[0], :] != arr[splits[3], :]).sum() / arr.shape[1]
    d_bc = np.bool(arr[splits[1], :] != arr[splits[2], :]).sum() / arr.shape[1]
    dist2 = d_ad + d_bc

    return dist0, dist1, dist2



def get_weight_snir(arr, splits):
    """ 
    Calculate the quartet weight as described in Avni et al. (2014) eq.1.
    Essentially our goal is to downweight quartets with very long edges
    """

    ## calculate JC distances
    distances = get_distances(arr, splits)

    ## sort distances
    dl, dm, dh = sorted(distances)

    ## calculate weight
    return (dh-dl) / (np.exp(dh-dm) * dh)



def run_qmc(data, useweights=False):
    """ runs quartet max-cut """

    ## use weight info
    weights = "on"
    if not useweights:
        weights = "off"

    ## output handle
    data._svd.tree = os.path.join(data.dirs.svd, data.name+"_svd4tet.tre")

    ## make call list
    cmd = [ipyrad.bins.qmc,
           " qrtt="+data._svd.qdump, 
           " weights="+weights, 
           " otre="+data._svd.tree]
    cmd = " ".join(cmd)

    ## run it
    LOGGER.info(cmd)
    try:
        subprocess.check_call(cmd, shell=True,
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error("Error in wQMC: \n({}).".format(inst))
        LOGGER.error(subprocess.STDOUT)
        LOGGER.error(cmd)
        raise inst

    ## convert int names back to str names
    








def dump(data):
    """ 
    prints the quartets to a file formatted for QMC in random order.
    """
    ## open the h5 database
    io5 = h5py.File(data._svd.results, 'r')

    ## create an output file for writing
    data._svd.qdump = os.path.join(data.dirs.svd, data.name+"_quartets.txt")
    outfile = open(data._svd.qdump, 'w')

    ## should pull quarts order in randomly
    for idx in range(0, 1000, 100):
        quarts = [list(j) for j in io5["quartets"][idx:idx+100]]
        weight = io5["weights"][idx:idx+100]
        chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
                                                     in zip(quarts, weight)]
        outfile.write("\n".join(chunk)+"\n")

    ## close output file and h5 database
    outfile.close()
    io5.close()



def insert_to_array(data, resqrt, reswgt, path):
    """ 
    Takes a tmpfile output from finished worker, enters it into the 
    full h5 array, and deletes the tmpfile
    """

    with h5py.File(path) as inh5:
        qrts = inh5['quartets'][:]
        wgts = inh5['weights'][:]

        chunk = data._svd.chunk
        start = int(path.split("_")[-1][:-3])
        resqrt[start:start+chunk] = qrts
        reswgt[start:start+chunk] = wgts

    os.remove(path)


def main(data, ipyclient, path=None, checkpoint=0):
    """ 
    main funcs
    """
    ## make an analysis directory if it doesn't exit
    data.dirs.svd = os.path.realpath(
                        os.path.join(
                            data.dirs.project, data.name+"_analysis_svd"))
    if not os.path.exists(data.dirs.svd):
        os.mkdir(data.dirs.svd)

    ## store some svd stats in data
    data._svd = ObjDict()
    data._svd.h5 = os.path.join(data.dirs.svd, data.name+"_input.h5")
    data._svd.results = os.path.join(data.dirs.svd, data.name+"_output.h5")    
    data._svd.checkpoint = 0

    ## checkpoint info will be retrieved from h5 if it exists
    if checkpoint:
        ## scan the hdf5 array to find where it left off writing
        data._svd.checkpoint = 0        
        raise IPyradWarningExit("not yet implemented")

    ## get the seq array into hdf5
    #LOGGER.info("loading array")
    print("loading array")
    data = get_seqarray(data, path)

    ## make quartet arrays into hdf5. Allow subsetting samples eventually.
    ## and optimize chunk value given remaining quartets and ipyclient    
    #LOGGER.info("loading quartets")
    print("loading quartets")
    samples = data.samples
    data = get_quartets(data, samples)

    ## load a job queue
    # jobs = []
    # for qidx in xrange(data._svd.checkpoint, 
    #                    data._svd.nquarts, 
    #                    data._svd.chunk):
    #     jobs.append([data, qidx])
    # ## THIS ONE WORKS GREAT, but loads all into mem
    # lbview = ipyclient.load_balanced_view()
    # res = lbview.map_async(worker, jobs)
    # res.wait_interactive()
    # print('done')
    # sys.exit()

    njobs = len(xrange(data._svd.checkpoint, 
                       data._svd.nquarts, data._svd.chunk))
    jobiter = iter(xrange(data._svd.checkpoint, 
                          data._svd.nquarts, data._svd.chunk))

    ## make a distributor for engines
    LOGGER.info("sending jobs to Engines")
    lbview = ipyclient.load_balanced_view()

    ## open a view to the super h5 array
    resh5 = h5py.File(data._svd.results, 'w')
    resqrt = resh5.create_dataset("quartets", (data._svd.nquarts, 4), 
                                  dtype=np.uint16, chunks=(data._svd.chunk, 4))
    reswgt = resh5.create_dataset("weights", (data._svd.nquarts,), 
                                  dtype=np.float16, chunks=(data._svd.chunk,))

    ## submit initial n jobs
    res = {}
    for i in range(len(ipyclient)):
        res[i] = lbview.apply(worker, [data, jobiter.next()])

    ## iterate over remaining jobs
    keys = res.keys()
    finished = 0
    while res.keys():
        time.sleep(1)
        for key in keys:
            try:
                ## query for finished results
                result = res[key].get(0)
                ## put it into the super array
                insert_to_array(data, resqrt, reswgt, result)
                ## submit a new job to the queue
                del res[key]
                finished += 1
                print(r"  {:.1f}% complete".format(100*(finished / njobs)))
                try:
                    res[key] = lbview.apply(worker, 
                                    [data, jobiter.next()])                
                except StopIteration:
                    continue

            except (ipp.error.TimeoutError, KeyError):
                continue
    resh5.close()

    ## convert to txt file for wQMC
    dump(data)    

    ## run quartet joining algorithm
    run_qmc(data)

    ## 
    print("  done")
    return data

    #dview = ipyclient[:]
    #dview.block = False

    ## a list of jobs to run
    #todo = [dview.apply_async(data, jobiter.next()) for _ in range(14)]

    #results = lbview.map_async(worker, jobs)
    #res = {}
    #for idx, job in enumerate(jobs):
    #    res[idx] = lbview.apply_async(worker, job)
    #results.get()

    # for idx in range(len(ipyclient)):
    #     #res[str(idx)] = lbview.apply_async(worker, [data, jobiter.next()])
    #     res[idx] = lbview.apply(worker, [data, jobiter.next()])

    # while 1:
    #     if not res.keys():
    #         print("breaking")
    #         break
    #     else:
    #         keys = res.keys()
    #         print('elapsed', res[keys[0]].elapsed)
    #         for idx in keys:
    #             try:
    #                 print(idx)
    #                 print("found this:", res[idx].get(0))

    #             except ipp.error.TimeoutError:
    #                 ## keep on looking
    #                 continue
    #             else:
    #                 try:
    #                     del res[idx]
    #                     res[idx] = lbview.apply_async(worker, 
    #                                                  [data, jobiter.next()])
    #                 except StopIteration:
    #                     pass
            
    #     time.sleep(0.25)

    #results.wait_interactive()
    #worker(jobs[0])

    #results = lbview.apply_async(worker, jobs)
    #results.wait_interactive()

    #while 1:
        ## hang out a while
    #    time.sleep(1)
        ## query for a finished job




    ## Dump database results into a file
    #dump(data)

    ## Run quartet max-cut on the quartets
    #run_qmc(data)



if __name__ == "__main__":

    ## imports
    import ipyrad.analysis as ipa
    import ipyrad as ip
    import ipyparallel as ipp
    import time

    IPYCLIENT = ipp.Client()
    time.sleep(2)

    DATA = ip.load_json("~/Documents/ipyrad/tests/cli/cli")

    ipa.svd4tet.main(DATA, IPYCLIENT)

