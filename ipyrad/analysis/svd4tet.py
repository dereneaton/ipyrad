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
import glob
import h5py
import time
import random
import ipyrad
import itertools
import subprocess
import numpy as np
import pandas as pd
import ipyparallel as ipp
from numba import jit

from ipyrad.assemble.util import ObjDict, IPyradWarningExit, progressbar
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

MAT = np.zeros((16, 16), dtype=np.int32)



# def regen_boot_array(data, resh5):
#     """ bootstrap replicates re-using same quartet order as original """

#     ## create the bootstrap sampled seqarray
#     with h5py.File(data.svd.h5, 'r+') as io5:
#         print("  resampling bootstrap array")
#         ## get original seqarray as a DF
#         seqarray = pd.DataFrame(resh5["seqarray"][:])
#         ## get boot array by randomly sample columns w/ replacement
#         bootarr = seqarray.sample(n=seqarray.shape[1], replace=True, axis=1)
#         resh5["bootarr"][:] = bootarr
#         del bootarr
#         ## save the boot array 
#         print("  done resampling bootstrap array")



def get_seqarray(data, dtype, boot):
    """ 
    Takes an Assembly object and looks for a phylip file unless the path 
    argument is used then it looks in path. 
    """

    ## check that input file is correct
    assert dtype in ['snp', 'usnp', 'phy'], \
        "only snp, usnp or phy data types allowed"
    ## read in the file
    spath = open(data.outfiles[dtype], 'r')
    line = spath.readline().strip().split()
    ntax = int(line[0])
    nbp = int(line[1])
    LOGGER.info("array shape: (%s,%s)", ntax, nbp)

    ## make a seq array
    tmpseq = np.zeros((ntax, nbp), dtype="S1")

    ## get initial array
    if not boot:
        ## use 'w' to make initial array
        with h5py.File(data.svd.h5in, 'w') as io5:
            ## create array storage for both real seq and for later bootstraps
            seqarr = io5.create_dataset("seqarr", (ntax, nbp), dtype="S1")
            io5.create_dataset("bootarr", (ntax, nbp), dtype="S1")
            ## fill the tmp array from the input phy
            for line, seq in enumerate(spath.readlines()):
                tmpseq[line] = np.array(list(seq.split()[-1]))
            ## save array to disk so it can be easily accessed from 
            ## many engines on arbitrary nodes 
            seqarr[:] = tmpseq
            del tmpseq
    else:
        ## use 'r+' to read and write to existing array
        with h5py.File(data.svd.h5in, 'r+') as io5:        
            ## load in the seqarr
            tmpseq = pd.DataFrame(io5["seqarr"][:])
            ## fill the boot array with a re-sampled phy w/ replacement
            io5["bootarr"][:] = np.array(tmpseq.sample(n=tmpseq.shape[1], 
                                         replace=True, axis=1), dtype="S1")
            del tmpseq
    return data



def get_quartets(data, samples):
    """ load all quartets into an array """

    ## calculate how many quartets to generate
    qiter = itertools.combinations(range(1, len(samples)+1), 4)
    nquarts = sum(1 for _ in qiter)

    ## create a chunk size for sampling from the array of quartets. This should
    ## be relatively large so that we don't spend a lot of time doing I/O, but
    ## small enough that jobs finish every few hours since that is how the 
    ## checkpointing works.
    if nquarts < 1000:
        chunk = nquarts // 20
    else:
        chunk = 500

    ## 'samples' stores the indices of the quartet. 
    ## `quartets` stores the correct quartet in the order (1,2|3,4)
    ## `weights` stores the calculated weight of the quartet in 'quartets'
    print("  populating array with {} quartets".format(nquarts))
    with h5py.File(data.svd.h5in, 'a') as io5:
        ## create data sets
        smps = io5.create_dataset("samples", (nquarts, 4), 
                            dtype=np.uint16, chunks=(chunk, 4))
        io5.create_dataset("quartets", (nquarts, 4), 
                            dtype=np.uint16, chunks=(chunk, 4))        
        
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
    data.svd.nquarts = nquarts
    data.svd.chunk = chunk

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
    """ wrapper for numba func """
    narr = arr[sidx, :].view(np.int8)
    mat = np.zeros((16, 16), dtype=np.int32)
    mat = jseq_to_matrix(narr, mat)
    return mat



@jit('i4[:,:](i1[:,:],i4[:,:])', nopython=True)
def jseq_to_matrix(narr, mat):
    """ 
    numba compiled code to get matrix fast.
    arr is a 4 x N seq matrix converted to np.int8
    I convert the numbers for ATGC into their respective index for the MAT
    matrix, and leave all others as high numbers, i.e., -==45, N==78. 
    """

    for x in xrange(narr.shape[1]):
        i = narr.T[x]
        ## convert to index values
        i[i == 65] = 0
        i[i == 67] = 1
        i[i == 71] = 2
        i[i == 84] = 3

        if np.sum(i) < 16:
            mat[i[0]*4:(i[0]+4)*4]\
               [i[1]]\
               [i[2]*4:(i[2]+4)*4]\
               [i[3]] += 1

    return mat


# def seq_to_matrix(arr, sidx):
#     """ 
#     Takes a 4-taxon alignment and generates a 16x16 count matrix. Ignores 
#     ambiguous sites, sites with missing data, and sites with indels. The 
#     order of the samples (sidx) is important. 
#     """
#     ## select only the relevant rows
#     arr = arr[sidx, :]
    
#     ## convert all nulls to Ns
#     for null in list("RSKYWM-"):
#         arr[arr == null] = "N"

#     ## mask any columns for this 4-taxon alignment that has Ns
#     arr = arr[:, ~np.any(arr == "N", axis=0)]

#     ## get dict of observed patterns
#     counts = tablestack(arr)

#     ## fill mdict with ordered patterns
#     mdict = OrderedDict([(i, 0) for i in MKEYS.split()])
#     for patt in counts:
#         mdict[patt] = counts[patt]

#     ## shape into 16x16 array
#     matrix = np.array(mdict.values()).reshape(16, 16)
#     LOGGER.info("matrix: %s", matrix)
#     return matrix



# def tablestack(arr):
#     """ makes a count dict of each unique array element """
#     ## goes by 10% at a time to minimize memory overhead. Is possible it skips
#     ## the last chunk, but this shouldn't matter. Check tho. 
#     table = Counter()
#     for i in xrange(0, arr.shape[1], arr.shape[1]//10):
#         tmp = Counter([j.tostring() for j in arr[:, i:i+arr.shape[1]//10].T])
#         table.update(tmp)
#     return table



def worker(args):
    """ gets data from seq array and puts results into results array """

    ## parse args
    LOGGER.info("I'm inside a worker")
    data, qidx = args
    chunk = data.svd.chunk

    ## get n chunk quartets from qidx
    lpath = os.path.realpath(data.svd.h5in)
    with h5py.File(lpath, 'r') as io5in:
        #start = time.time()
        smps = io5in["samples"][qidx:qidx+chunk]

        if not data.svd.checkpoint_boot:
            seqs = io5in["seqarr"][:]
        else:
            seqs = io5in["bootarr"][:]            
        #LOGGER.info("takes %s secs to read smps & seqs", (time.time() - start))

        ## get data from seqarray
        #start = time.time()        
        results = [svd4tet(seqs, qtet) for qtet in smps]
        rquartets = np.array([i[0] for i in results])
        rweights = np.array([i[1] for i in results])

    tmpchunk = os.path.join(data.dirs.svd, data.name+"_tmp_{}.h5".format(qidx))
    with h5py.File(tmpchunk, 'w') as io5out:
        #io5out["weights"][qidx:qidx+chunk] = np.array(rweights)
        #io5out["quartets"][qidx:qidx+chunk] = np.array(rquartets)
        io5out["weights"] = np.array(rweights)
        io5out["quartets"] = np.array(rquartets)

    return tmpchunk



def svd4tet(arr, samples):
    """ calc rank. From Kubatko and Chiffman (2014)

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
    #mats = [seq_to_matrix(arr, sidx) for sidx in splits]
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
    calc = (dh-dl) / (np.exp(dh-dm) * dh)
    #print(calc)
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
        LOGGER.info(cmd)
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
    ## ... dendropy or ete2 code

    ## save the tree
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




def dump(data):
    """ 
    prints the quartets to a file formatted for QMC in random order.
    """
    ## open the h5 database
    io5 = h5py.File(data.svd.h5out, 'r')

    ## create an output file for writing
    data.svd.qdump = os.path.join(data.dirs.svd, data.name+"_quartets.txt")
    outfile = open(data.svd.qdump, 'w')

    ## todo: should pull quarts order in randomly
    for idx in range(0, 1000, 100):
        quarts = [list(j) for j in io5["quartets"][idx:idx+100]]
        weight = io5["weights"][idx:idx+100]
        chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
                                                     in zip(quarts, weight)]
        outfile.write("\n".join(chunk)+"\n")

    ## close output file and h5 database
    outfile.close()
    io5.close()



def insert_to_array(data, result):
    """ 
    Takes a tmpfile output from finished worker, enters it into the 
    full h5 array, and deletes the tmpfile
    """

    out5 = h5py.File(data.svd.h5out, 'r+')

    with h5py.File(result) as inh5:
        qrts = inh5['quartets'][:]
        wgts = inh5['weights'][:]

        chunk = data.svd.chunk
        start = int(result.split("_")[-1][:-3])
        out5['quartets'][start:start+chunk] = qrts
        out5['weights'][start:start+chunk] = wgts
    out5.close()

    ## remove tmp h5path
    os.remove(result)



def svd_obj_init(data):
    """ creates svd attribute to Assembly object """
    ## create ObjDict for storage
    data.svd = ObjDict()

    ## add array path
    data.svd.h5in = os.path.join(data.dirs.svd, data.name+"_input.h5")
    data.svd.h5out = os.path.join(data.dirs.svd, data.name+"_output.h5")

    ## original tree paths
    data.svd.tre = os.path.join(data.dirs.svd, data.name+"_svd4tet.tre")
    data.svd.wtre = os.path.join(data.dirs.svd, data.name+"_svd4tet.wtre")

    ## bootstrap tree paths
    data.svd.tboots = os.path.join(data.dirs.svd, 
                                   data.name+"_svd4tet.tre.boots")
    data.svd.wboots = os.path.join(data.dirs.svd, 
                                   data.name+"_svd4tet.wtre.boots")

    ## bootstrap labeled o.g. trees paths
    data.svd.btre = os.path.join(data.dirs.svd, 
                                 data.name+"_svd4tet.tre.boots.support")
    data.svd.bwtre = os.path.join(data.dirs.svd, 
                                 data.name+"_svd4tet.wtre.boots.support")
    ## checkpoints
    data.svd.checkpoint_boot = 0
    data.svd.checkpoint_arr = 0

    ## save to object and return
    data.save()
    return data




def wrapper(data, samples=None, dtype='snp', nboots=100, force=False):
    """ wraps main in try/except statement """

    ## launch ipclient, assumes ipyparallel is running
    ipyclient = ipp.Client(timeout=10)

    ## protects it from KBD
    try:
        main(data, samples, ipyclient, dtype, nboots, force)

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




def main(data, samples, ipyclient, dtype, nboots, force):
    """ 
    Run svd4tet inference on a sequence or SNP alignment for all samples 
    the Assembly. 

    By default the job starts from 0 or where it last left off, unless 
    force=True, then it starts from 0. 
    """

    ## subset the samples
    if samples:
        samples = ipyrad.core.assembly._get_samples(data, samples)

    ## load svd attributes if they exist
    fresh = 0
    if not force:
        try:
            if data.svd.checkpoint_boot or data.svd.checkpoint_arr:
                print("  loading from svd checkpoint")
                print("  array checkpoint {}".format(data.svd.checkpoint_arr))
                print("  boots checkpoint {}".format(data.svd.checkpoint_boot))
            else:
                fresh = 1
        except (AttributeError, IOError):
            print("  starting new svd analysis")
            fresh = 1

    ## if svd results do not exist or force then restart
    if force or fresh:
        ## make an analysis directory if it doesn't exit
        data.dirs.svd = os.path.realpath(
                            os.path.join(
                                data.dirs.project, data.name+"_analysis_svd"))
        if not os.path.exists(data.dirs.svd):
            try:
                os.mkdir(data.dirs.svd)
            except OSError:
                ## if path doesn't exist (i.e., changed computers for analysis)
                ## then create new svd directory
                data.dirs.svd = os.path.join(
                                    os.path.curdir, data.name+"_analysis_svd")
                os.mkdir(data.dirs.svd)
                print("  new dir/ created: {}".format(data.dirs.svd))

        ## init new svd
        data = svd_obj_init(data)

        ## get the real seq array into hdf5 h5in
        print("  loading array")
        data = get_seqarray(data, dtype=dtype, boot=False)

        ## make quartet arrays into hdf5. Allow subsetting samples eventually.
        ## and optimize chunk value given remaining quartets and ipyclient    
        print("  loading quartets")
        samples = data.samples
        data = get_quartets(data, samples)

    ## run the full inference 
    if not data.svd.checkpoint_boot:
        print("  inferring quartets for the full data")
        inference(data, ipyclient, bidx=0)

    ## run the bootstrap replicates
    print("  running {} bootstrap replicates".format(nboots))
    ## get current boot
    for bidx in range(data.svd.checkpoint_boot, nboots):
        ## get a new bootarray if starting a new boot
        LOGGER.info("  boot = {}".format(bidx))
        if data.svd.checkpoint_arr == 0:
            data = get_seqarray(data, dtype=dtype, boot=True)
            LOGGER.info("  new boot array sampled")
            data.svd.checkpoint_boot = bidx
        ## start boot inference
        progressbar(nboots, bidx)
        inference(data, ipyclient, bidx=True)
    progressbar(nboots, nboots)

    ## create tree with support values on edges
    ## ...

    ## print finished
    print("  Finished. Final tree files: \n{}\n{}\n{}\n{}\n{}\n{}".format(
          data.svd.tre, 
          data.svd.tboots,
          data.svd.wtre, 
          data.svd.wboots,
          data.svd.btre,
          data.svd.bwtre,
          ))

    return data



def inference(data, ipyclient, bidx):
    """ run inference and store results """

    ## a distributor of chunks
    njobs = sum(1 for _ in iter(xrange(data.svd.checkpoint_arr, 
                                       data.svd.nquarts, data.svd.chunk)))
    jobiter = iter(xrange(data.svd.checkpoint_arr, 
                          data.svd.nquarts, data.svd.chunk))
    LOGGER.info("chunksize: %s, start: %s, total: %s, njobs: %s", \
            data.svd.chunk, data.svd.checkpoint_arr, data.svd.nquarts, njobs)

    ## make a distributor for engines
    lbview = ipyclient.load_balanced_view()
    LOGGER.info("sending jobs to %s Engines", len(ipyclient))

    ## open a view to the super h5 array
    with h5py.File(data.svd.h5out, 'w') as out5:
        out5.create_dataset("quartets", (data.svd.nquarts, 4), 
                            dtype=np.uint16, chunks=(data.svd.chunk, 4))
        out5.create_dataset("weights", (data.svd.nquarts,), 
                            dtype=np.float16, chunks=(data.svd.chunk,))

    ## submit initial n jobs
    assert len(ipyclient) > 0, "No ipyparallel Engines found"
    res = {}
    for i in range(len(ipyclient)):
        try:
            res[i] = lbview.apply(worker, [data, jobiter.next()])
        except StopIteration:
            continue

    ## iterate over remaining jobs
    keys = res.keys()
    finished = 0

    while res.keys():
        time.sleep(3)
        progressbar(njobs, finished)
        for key in keys:
            try:
                ## query for finished results
                result = res[key].get(0)
                ## put it into the super array
                insert_to_array(data, result)
                ## delete result, update checkpoint
                del res[key]
                finished += 1
                ## update the minimum quartets finished/filled.
                with h5py.File(data.svd.h5out, 'r') as tmp5:
                    ww = tmp5["weights"][:]
                    try:
                        data.svd.checkpoint_arr = np.where(ww == 0)[0].min()
                        LOGGER.info("arr saved at %s", data.svd.checkpoint_arr)
                    except (ValueError, AttributeError):
                        ## array is full (no zeros)
                        pass
                ## submit new jobs
                try:
                    res[key] = lbview.apply(worker, 
                                    [data, jobiter.next()])
                    LOGGER.info("new job added to Engine %s", key)
                except StopIteration:
                    continue

            except (ipp.error.TimeoutError, KeyError):
                continue
    progressbar(njobs, finished)                
    print("")

    ## convert to txt file for wQMC
    dump(data)    

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
    DATA = ipyrad.load_json("~/Documents/RADmissing/rad1/half_min20.json")
    ## run
    ipa.svd4tet.wrapper(DATA)

