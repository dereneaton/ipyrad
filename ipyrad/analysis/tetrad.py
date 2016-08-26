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
# pylint: disable=E1103
# pylint: disable=F0401
# pylint: disable=W0212
# pylint: disable=W0142
# pylint: disable=C0103
# pylint: disable=C0301
# pylint: disable=R0914
# pylint: disable=R0915


from __future__ import print_function, division
import os
import json
import h5py
import time
import numba
import random
import socket
import datetime
import itertools
import subprocess
import numpy as np
import ipyrad as ip
from bitarray import bitarray
from fractions import Fraction
from collections import defaultdict
from ipyrad.assemble.util import ObjDict, IPyradWarningExit, progressbar

## for our desired form of parallelism we will limit 1 thread per cpu
numba.config.NUMBA_DEFAULT_NUM_THREADS = 1

## debug numba code
#numba.config.NUMBA_DISABLE_JIT = 1

## ete3 is an extra dependency not included with ipyrad
## replace with biopython asap
try:
    import ete3
except ImportError:
    try:
        import ete2 as ete3
    except ImportError:
        raise IPyradWarningExit("""
    tetrad requires the dependency `ete3`. You can install
    it with the command `conda install -c etetoolkit ete3`
    Sorry for the inconvenience, this will be incorporated into the
    ipyrad installation eventually.
    """)

# try:
#     import skbio.tree as sktree
#     from io import StringIO
# except ImportError:
#     raise IPyradWarningExit("""
#     tetrad requires the dependency `biopython`. You can install
#     it with the command `conda install -c anaconda biopython`. 
#     Sorry for the inconvenience, this will be incorporated into the
#     ipyrad installation eventuall.
#     """)        

## set the logger
import logging
LOGGER = logging.getLogger(__name__)

## The 16 x 16 matrix of site counts (The phylogenetic invariants). 
## It's here just to look at. 
PHYLO_INVARIANTS = """
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



#############################################################################
#############################################################################
## Quartet inference Class Object
#############################################################################
#############################################################################


class Quartet(object):
    """
    The main tetrad object for storing data and checkpointing. It is 
    initialized with a name, and args to command line (e.g., sampling method, 
    starting tree, nboots, etc.). 
    """

    def __init__(self, name, wdir=os.path.curdir, method='all'):
        ## version is ipyrad version
        self._version = ip.__version__

        ## name this assembly
        self.name = name
        self.dirs = os.path.realpath(wdir)
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)

        ## store default cluster information 
        self._ipcluster = {
            "cluster_id" : "",
            "profile" : "default", 
            "engines" : "Local", 
            "quiet" : 0, 
            "timeout" : 60, 
            "cores" : ip.assemble.util.detect_cpus()}

        ## Sampling method attributes
        self.method = method
        self.nboots = 0
        self.nquartets = 0
        self.chunksize = 0
        self.resolve = 0

        ## store samples from the seqarray
        self.samples = []

        ## self.populations ## if we allow grouping samples
        ## (haven't done this yet)

        ## hdf5 data bases init and delete existing
        self.h5in = os.path.join(self.dirs, self.name+".input.h5")
        self.h5out = os.path.join(self.dirs, self.name+".output.h5")

        ## input files
        self.files = ObjDict()
        self.files.seqfile = None
        self.files.mapfile = None
        self.files.treefile = None
        self.files.qdump = None

        ## store tree file paths
        self.trees = ObjDict()
        ## the full trees
        self.trees.tre = os.path.join(self.dirs, self.name+".full.tre")
        ## the extended majority rule consensus tree w/ support values
        self.trees.cons = os.path.join(self.dirs, self.name+".consensus.tre")
        ## all bootstrap trees 
        self.trees.boots = os.path.join(self.dirs, self.name+".boots")        
        ## NHX formatted tre with rich information
        self.trees.nhx = os.path.join(self.dirs, self.name+".nhx.tre")     
        ## a file for tree and quartet stats
        self.trees.stats = os.path.join(self.dirs, self.name+".stats.txt")
        ## checkpointing information
        self.checkpoint = ObjDict()
        self.checkpoint.boots = 0
        self.checkpoint.arr = 0



    def refresh(self):
        """ 
        Remove all existing results files and reinit the h5 arrays 
        so that the Quartet object is just like fresh from a CLI start
        """

        ## clear any existing results files
        oldfiles = [self.files.qdump] + self.trees.values()
        for oldfile in oldfiles:
            if oldfile:
                if os.path.exists(oldfile):
                    os.remove(oldfile) 

        ## io5 input has 'samples' and 'seqarr' keys which can be left
        ## io5 output has 'qboots', 'qstats', 'quartets', 'weights'
        with h5py.File(self.h5out, 'w') as io5:
            ## remove old arrs
            for key in io5.keys():
                del io5[key]
            ## create fresh ones
            io5.create_dataset("quartets", (self.nquartets, 4), 
                                dtype=np.uint16, chunks=(self.chunksize, 4))
            io5.create_dataset("weights", (self.nquartets,), 
                                dtype=np.float64, chunks=(self.chunksize, ))
            io5.create_dataset("qstats", (self.nquartets, 4), 
                                dtype=np.uint32, chunks=(self.chunksize, 4))
            io5.create_group("qboots")        

        ## reset metadata
        self.checkpoint.array = 0
        self.checkpoint.boots = 0



    def parse_names(self):
        """ parse names from the seqfile """
        ## parse samples from the sequence file
        self.samples = []
        with iter(open(self.files.seqfile, 'r')) as infile:
            infile.next().strip().split()
            while 1:
                try:
                    self.samples.append(infile.next().split()[0])
                except StopIteration:
                    break
        ## names are in the order of the sequences in seqfile
        #self.samples = sorted(self.samples)



    ## Filling the h5in seqarray
    def init_seqarray(self):
        """ 
        Fills the seqarr with the full data set, and creates a bootsarr copy
        with the following modifications:

        1) converts "-" into "N"s, since they are similarly treated as missing. 
        2) randomly resolve ambiguities (RSKWYM)
        3) convert to uint8 for smaller memory load and faster computation
        """

        ## read in the file
        try:
            spath = open(self.files.seqfile, 'r')
        except IOError:
            raise IPyradWarningExit(NO_SNP_FILE\
                                    .format(self.files.seqfile))
        line = spath.readline().strip().split()
        ntax = int(line[0])
        nbp = int(line[1])

        ## make a tmp seq array
        print("  loading seq array [{} taxa x {} bp]".format(ntax, nbp))        
        tmpseq = np.zeros((ntax, nbp), dtype=np.uint8)
    
        ## create array storage for real seq and the tmp bootstrap seqarray
        with h5py.File(self.h5in, 'w') as io5:
            io5.create_dataset("seqarr", (ntax, nbp), dtype=np.uint8)
            io5.create_dataset("bootsarr", (ntax, nbp), dtype=np.uint8)
            io5.create_dataset("bootsmap", (nbp, 2), dtype=np.uint32)

            ## if there is a map file, load it into the bootsmap
            if self.files.mapfile:
                with open(self.files.mapfile, 'r') as inmap:
                    ## parse the map file from txt and save as dataset
                    maparr = np.genfromtxt(inmap, dtype=np.uint32)
                    io5["bootsmap"][:] = maparr[:, [0, 3]]

                    ## parse the span info from maparr and save to dataset
                    spans = np.zeros((maparr[-1, 0], 2), np.uint64)
                    spans = get_spans(maparr, spans)
                    io5.create_dataset("spans", data=spans)
                    print("  max unlinked SNPs per quartet: {}".format(spans.shape[0]))

            ## fill the tmp array from the input phy
            for line, seq in enumerate(spath.readlines()):
                tmpseq[line] = np.array(list(seq.split()[-1])).view(np.uint8)

            ## convert '-' into 'N'
            tmpseq[tmpseq == 45] = 78

            ## save array to disk so it can be easily accessed by slicing
            ## This unmodified array is used again later for sampling boots
            io5["seqarr"][:] = tmpseq

            ## resolve ambiguous IUPAC codes
            if self.resolve:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## save modified array to disk            
            io5["bootsarr"][:] = tmpseq

            ## memory cleanup
            #del tmpseq

            ## get initial array
            LOGGER.info("original seqarr \n %s", io5["seqarr"][:, :20])
            LOGGER.info("original bootsarr \n %s", io5["bootsarr"][:, :20])
            LOGGER.info("original bootsmap \n %s", io5["bootsmap"][:20, :])



    def sample_bootseq_array(self):
        """ 
        Takes the seqarray and re-samples columns chunks based on linkage 
        from the maps array. Saves the new map information in the bootmap. 

        The seqarray is 'cleaned up' in several ways: 
        1) converts "-" into "N"s, since they are similarly treated as missing. 
        2) randomly resolve ambiguities (RSKWYM)
        3) convert to uint8 for smaller memory load and faster computation
        """

        ## use 'r+' to read and write to existing array
        with h5py.File(self.h5in, 'r+') as io5:        
            ## load in the seqarr and maparr
            seqarr = io5["seqarr"][:]

            ## resample columns with replacement
            newarr = np.zeros(seqarr.shape, dtype=np.uint8)
            cols = np.random.randint(0, seqarr.shape[1], seqarr.shape[1])
            tmpseq = shuffle_cols(seqarr, newarr, cols)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self.resolve:
                tmpseq = resolve_ambigs(tmpseq)
        
            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## fill the boot array with a re-sampled phy w/ replacement
            io5["bootsarr"][:] = tmpseq
            del tmpseq



    def sample_bootseq_array_map(self):
        """ Re-samples loci with replacement to fill the bootarr """

        ## open view to the in database
        with h5py.File(self.h5in, 'r+') as io5:
            ## load the original data (seqarr and spans)
            seqarr = io5["seqarr"][:]
            spans = io5["spans"][:]            

            ## get size of the new locus re-samples array
            nloci = spans.shape[0]
            loci = np.random.choice(nloci, nloci)
            arrlen = get_shape(spans, loci)

            ## create a new bootsarr and maparr to fill
            del io5["bootsarr"]
            del io5["bootsmap"]
            newbarr = np.zeros((seqarr.shape[0], arrlen), dtype=np.uint8)
            newbmap = np.zeros((arrlen, 2), dtype=np.uint32)
            newbmap[:, 1] = np.arange(1, arrlen+1)
            
            ## fill the new arrays            
            tmpseq, tmpmap = fill_boot(seqarr, newbarr, newbmap, spans, loci)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self.resolve:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## store data sets
            io5.create_dataset("bootsmap", data=tmpmap)
            io5.create_dataset("bootsarr", data=tmpseq)

            LOGGER.info("resampled bootsarr \n %s", io5["bootsarr"][:, :10])
            LOGGER.info("resampled bootsmap \n %s", io5["bootsmap"][:10, :])



    ## Functions to fill h5in with samples
    def store_N_samples(self):
        """ Find all quartets of samples and store in a large array """
        ## create a chunk size for sampling from the array of quartets. This should
        ## be relatively large so that we don't spend a lot of time doing I/O, but
        ## small enough that jobs finish every few hours for checkpointing
        breaks = 2
        if self.nquartets < 5000:
            breaks = 1
        if self.nquartets > 100000:
            breaks = 4
        if self.nquartets > 500000:
            breaks = 8

        cpus = self._ipcluster["cores"]
        self.chunksize = (self.nquartets // (breaks * cpus) + \
                         (self.nquartets % (breaks * cpus)))
        LOGGER.info("nquarts = %s, chunk = %s", self.nquartets, self.chunksize)

        ## 'samples' stores the indices of the quartet. 
        ## `quartets` stores the correct quartet in the order (1,2|3,4)
        ## `weights` stores the calculated weight of the quartet in 'quartets'
        ## we gzip this for now, but check later if this has a big speed cost

        ## create h5 OUT empty arrays
        with h5py.File(self.h5out, 'w') as io5:
            io5.create_dataset("quartets", (self.nquartets, 4), 
                                dtype=np.uint16, chunks=(self.chunksize, 4))
            io5.create_dataset("weights", (self.nquartets,), 
                                dtype=np.float64, chunks=(self.chunksize, ))
            io5.create_dataset("qstats", (self.nquartets, 4), 
                                dtype=np.uint32, chunks=(self.chunksize, 4))
            io5.create_group("qboots")


        ## append to h5 IN array (which also has seqarray) and fill it
        with h5py.File(self.h5in, 'a') as io5:
            ## create data sets
            io5.create_dataset("samples", (self.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self.chunksize, 4),
                               compression='gzip')

            ## populate array with all possible quartets. This allows us to 
            ## sample from the total, and also to continue from a checkpoint
            qiter = itertools.combinations(xrange(len(self.samples)), 4)
            i = 0

            ## fill chunksize at a time for efficiency
            while i < self.nquartets:
                if self.method != "all":
                    ## grab the next random 1000
                    qiter = []
                    while len(qiter) < min(self.chunksize, io5["samples"].shape[0]):
                        qiter.append(
                            random_combination(range(len(self.samples)), 4))
                    dat = np.array(qiter)
                else:
                    ## grab the next ordered chunksize
                    dat = np.array(list(itertools.islice(qiter, self.chunksize)))

                ## store to h5 
                io5["samples"][i:i+self.chunksize] = dat[:io5["samples"].shape[0] - i]
                i += self.chunksize



    def store_equal_samples(self):
        """ 
        sample quartets evenly across splits of the starting tree, and fills
        in remaining samples with random quartet samples. Uses a hash dict to 
        not sample the same quartet twice, so for very large trees this can 
        take a few minutes to find millions of possible quartet samples. 
        """
        
        ## choose chunker for h5 arr
        breaks = 2
        if self.nquartets < 5000:
            breaks = 1
        if self.nquartets > 100000:
            breaks = 4
        if self.nquartets > 500000:
            breaks = 8

        cpus = self._ipcluster["cores"]
        self.chunksize = (self.nquartets // (breaks * cpus) + \
                         (self.nquartets % (breaks * cpus)))
        LOGGER.info("nquarts = %s, chunk = %s", self.nquartets, self.chunksize)

        ## create h5 OUT empty arrays
        with h5py.File(self.h5out, 'w') as io5:
            io5.create_dataset("quartets", (self.nquartets, 4), 
                                dtype=np.uint16, chunks=(self.chunksize, 4))
            io5.create_dataset("weights", (self.nquartets,), 
                                dtype=np.float64, chunks=(self.chunksize, ))
            io5.create_dataset("qstats", (self.nquartets, 4), 
                                dtype=np.uint32, chunks=(self.chunksize, 4))
            io5.create_group("qboots")

        ## get starting tree, unroot, randomly resolve, ladderize
        tre = ete3.Tree(self.files.treefile, format=0)
        tre.unroot()
        tre.resolve_polytomy(recursive=True)
        tre.ladderize()

        ## randomly sample all splits of tree and convert tip names to indices
        splits = [([self.samples.index(z.name) for z in i], 
                   [self.samples.index(z.name) for z in j]) \
                   for (i, j) in tre.get_edges()]
    
        ## only keep internal splits (no single tips edges)
        ## this seemed to cause problems with unsampled tips
        splits = [i for i in splits if all([len(j) > 1 for j in i])]

        ## turn each into an iterable split sampler
        ## if the nquartets for that split is small, then sample all of them
        ## if it is big, then make it a random sampler from that split
        qiters = []

        ## how many min quartets are we gonna sample from each split?
        squarts = self.nquartets // len(splits)

        ## how many iterators can be sampled to saturation?
        nsaturation = 0

        for split in splits:
            ## if small number at this split then sample all possible sets
            ## we will exhaust this quickly and then switch to random for 
            ## the larger splits.
            if n_choose_k(len(split[0]), 2) * n_choose_k(len(split[1]), 2) < squarts*2:
                qiter = (i+j for (i, j) in itertools.product(
                            itertools.combinations(split[0], 2), 
                            itertools.combinations(split[1], 2)))
                nsaturation += 1

            ## else create random sampler across that split, this is slower
            ## because it can propose the same split repeatedly and so we 
            ## have to check it against the 'sampled' set.
            else:
                qiter = (random_product(split[0], split[1]) for _ \
                         in xrange(self.nquartets))
                nsaturation += 1

            ## store all iterators into a list
            qiters.append(qiter)

        #for split in splits:
        #    print(split)

        ## make qiters infinitely cycling
        qiters = itertools.cycle(qiters)
        cycler = itertools.cycle(range(len(splits)))

        ## store visiting quartets
        sampled = set()

        ## iterate over qiters sampling from each, if one runs out, keep 
        ## sampling from remaining qiters. Keep going until samples is filled
        with h5py.File(self.h5in, 'a') as io5:
            ## create data sets
            io5.create_dataset("samples", (self.nquartets, 4), 
                                      dtype=np.uint16, 
                                      chunks=(self.chunksize, 4),
                                      compression='gzip')

            ## fill chunksize at a time for efficiency
            i = 0
            empty = set()
            edge_targeted = 0
            random_target = 0

            ## keep filling quartets until nquartets are sampled
            while i < self.nquartets:
                qdat = []
                ## keep filling this chunk until its full
                while len(qdat) < self.chunksize:
                    ## grab the next iterator
                    qiter = qiters.next()
                    cycle = cycler.next()

                    ## sample from iterator
                    try:
                        qrtsamp = qiter.next()
                        if tuple(qrtsamp) not in sampled:
                            qdat.append(qrtsamp)
                            sampled.add(qrtsamp)
                            edge_targeted += 1
                        #else:
                        #    print('repeat')
                        
                    ## unless iterator is empty, then skip it
                    except StopIteration:
                        empty.add(cycle)

                    ## break when all edge samplers are empty
                    if len(empty) == nsaturation:
                        break

                ## if array is not full then add random samples
                while len(qdat) < self.chunksize:
                    qrtsamp = random_combination(range(len(self.samples)), 4)
                    if tuple(qrtsamp) not in sampled:
                        qdat.append(qrtsamp)
                        sampled.add(qrtsamp)
                        random_target += 1

                ## stick chunk into h5 array
                dat = np.array(qdat, dtype=np.uint16)
                io5["samples"][i:i+self.chunksize] = dat[:io5["samples"].shape[0] - i]
                i += self.chunksize

            print("  equal sampling: {} edge quartets, {} random quartets "\
                  .format(edge_targeted, random_target))



    def run_qmc(self, boot):
        """ runs quartet max-cut on a quartets file """

        cmd = " ".join(
                [ip.bins.qmc,
                " qrtt="+self.files.qdump,
                " weights=on"+
                " otre=.tmpwtre"])

        ## run them
        try:
            subprocess.check_call(cmd, shell=True,
                                       stderr=subprocess.STDOUT,
                                       stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as inst:
            LOGGER.error("Error in wQMC: \n({}).".format(inst))
            LOGGER.error(subprocess.STDOUT)
            raise inst

        ## read in the tmp files since qmc does not pipe
        intmpwtre = open(".tmpwtre", 'r')

        ## convert int names back to str names
        tmpwtre = self.renamer(ete3.Tree(intmpwtre.read().strip()))

        ## save the boot tree
        if boot:
            with open(self.trees.boots, 'a') as outboot:
                outboot.write(tmpwtre+"\n")

        ## store full data trees to Assembly
        else:
            with open(self.trees.tre, 'w') as outtree:
                outtree.write(tmpwtre)

        ## save JSON file checkpoint
        intmpwtre.close()
        self.save()



    def dump_qmc(self):
        """ 
        Makes a reduced array that excludes quartets with no information and 
        prints the quartets and weights to a file formatted for wQMC 
        """
        ## open the h5 database
        io5 = h5py.File(self.h5out, 'r')

        ## create an output file for writing
        self.files.qdump = os.path.join(self.dirs, self.name+".quartets.txt")
        LOGGER.info("qdump file %s", self.files.qdump)
        outfile = open(self.files.qdump, 'w')

        ## todo: should pull quarts order in randomly?
        for idx in xrange(0, self.nquartets, self.chunksize):
            ## get mask of zero weight quartets
            mask = io5["weights"][idx:idx+self.chunksize] != 0
            LOGGER.info("exluded = %s, mask shape %s", 
                        self.chunksize - mask.shape[0], mask.shape)

            ## apply mask
            LOGGER.info('q shape %s', io5["quartets"][idx:idx+self.chunksize].shape)
            masked_quartets = io5["quartets"][idx:idx+self.chunksize, :][mask, :]


            quarts = [list(j) for j in masked_quartets]
            weight = io5["weights"][idx:idx+self.chunksize][mask]

            ## format and print
            chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
                                                    in zip(quarts, weight)]
            outfile.write("\n".join(chunk)+"\n")

        ## close output file and h5 database
        outfile.close()
        io5.close()



    def renamer(self, tre):
        """ renames newick from numbers to sample names"""
        ## get the tre with numbered tree tip labels
        names = tre.get_leaves()

        ## replace numbered names with snames
        for name in names:
            name.name = self.samples[int(name.name)]

        ## return with only topology and leaf labels
        return tre.write(format=9)



    def write_output_splash(self):
        """ write final tree files """

        ## print stats file location:
        print(STATSOUT.format(opr(self.trees.stats)))

        ## print finished tree information ---------------------
        print(FINALTREES.format(opr(self.trees.tre)))

        ## print bootstrap information --------------------------
        if self.nboots:
            ## get consensus, map values to tree edges, record stats file
            self._compute_tree_stats()
            ## print bootstrap info
            print(BOOTTREES.format(opr(self.trees.cons),
                                   opr(self.trees.boots))) 

        ## print the ASCII tree only if its small
        if len(self.samples) < 200:
            if self.nboots:
                wctre = ete3.Tree(self.trees.cons, format=0)
                #wctre.unroot()
                #wctre.ladderize()
                print(wctre.get_ascii(show_internal=True, 
                                      attributes=["dist", "name"]))
                print("")
            else:
                qtre = ete3.Tree(self.trees.tre, format=0)
                qtre.unroot()
                print(qtre.get_ascii())
                print("")

        ## print PDF filename & tips -----------------------------
        docslink = "ipyrad.readthedocs.org/cookbook.html"    
        citelink = "ipyrad.readthedocs.org/tetrad.html"
        print(LINKS.format(docslink, citelink))



    def _compute_tree_stats(self):
        """ writes support values as edge labels on unrooted tree """

        ## get name indices
        names = self.samples

        ## get majority rule consensus tree of weighted Q bootstrap trees
        if self.nboots:
            fulltre = ete3.Tree(self.trees.tre, format=0)
            with open(self.trees.boots, 'r') as inboots:
                wboots = [fulltre] + \
                  [ete3.Tree(i.strip(), format=0) for i in inboots.readlines()]
            wctre, wcounts = consensus_tree(wboots, names=names)


        ## build stats file
        with open(self.trees.stats, 'w') as ostats:

            ## print Quartet info
            ostats.write("## Analysis info\n")
            ostats.write("{:<30}  {:<20}\n".format("Name", self.name))
            ostats.write("{:<30}  {:<20}\n".format("Sampling_method", self.method))
            ostats.write("{:<30}  {:<20}\n".format("Sequence_file", self.files.seqfile))
            ostats.write("{:<30}  {:<20}\n".format("Map_file", self.files.mapfile))
            used_treefile = [self.files.treefile if self.method == 'equal' else None][0]
            ostats.write("{:<30}  {:<20}\n".format("Guide_tree", used_treefile))
            ostats.write("\n")

            ## calculate Quartet stats
            ostats.write("## Quartet statistics (coming soon)\n")
            ostats.write("{:<30}  {:<20}\n".format("N_sampled_quartets", self.nquartets))
            proportion = 100*(self.nquartets / float(n_choose_k(len(self.samples), 4)))
            ostats.write("{:<30}  {:<20.1f}\n".format("percent_sampled_of_total", proportion))
            mean_loci = 0
            mean_snps = 0
            mean_weight = 0
            mean_dstat = 0
            ostats.write("{:<30}  {:<20}\n".format("Mean_N_loci_per_split", mean_loci))
            ostats.write("{:<30}  {:<20}\n".format("Mean_SNPs_per_split", mean_snps))
            ostats.write("{:<30}  {:<20}\n".format("Mean_quartet_weight", mean_weight))
            ostats.write("{:<30}  {:<20}\n".format("Mean_abba_baba", mean_dstat))
            ostats.write("\n")

            ## print tree output files info
            ostats.write("## Tree files\n")
            ostats.write("{:<30}  {:<20}\n".format("Initial_tree", self.trees.tre))
            ostats.write("{:<30}  {:<20}\n".format("bootstrap_replicates", self.trees.boots))
            ostats.write("{:<30}  {:<20}\n".format("extended_majrule_consens", self.trees.cons))            
            ostats.write("\n")

            ## print bootstrap splits
            if self.nboots:
                ostats.write("## splits observed in {} trees\n".format(len(wboots)))
                for i, j in enumerate(self.samples):
                    ostats.write("{:<3} {}\n".format(i, j))
                ostats.write("\n")
                for split, freq in wcounts:
                    if split.count('1') > 1:
                        ostats.write("{}   {:.2f}\n".format(split, round(freq, 2)))
                ostats.write("\n")


        ## parallelized this function because it can be slogging
        ipyclient = ip.core.parallel.get_client(**self._ipcluster)
        lbview = ipyclient.load_balanced_view()
        
        ## store results in dicts
        qtots = {}
        qsamp = {}
        tots = set(wctre.get_leaf_names())
        ## iterate over node traversal. 
        for node in wctre.traverse():
            ## this is slow, needs to look at every sampled quartet
            ## so we send it be processed on an engine
            qtots[node] = lbview.apply(_get_total, *(tots, node))
            qsamp[node] = lbview.apply(_get_sampled, *(self, tots, node, names))

        ## wait for jobs to finish
        ipyclient.wait()

        ## put results into tree
        for node in wctre.traverse():
            ## this is fast, just calcs n_choose_k
            total = qtots[node].result()
            sampled = qsamp[node].result()
            ## store the results to the tree            
            node.add_feature("quartets_possible", total)
            node.add_feature("quartets_sampled", sampled)

        features = ["quartets_total", "quartets_sampled"]

        ## return as NHX format with extra info
        with open(self.trees.cons, 'w') as outtre:
            outtre.write(wctre.write(format=0, features=features))



    ## not currently being used and its ugly, just provide a cookbook
    # def quickfig(self):
    #     """ make a quick ete3 fig. Plots total quartets """
    #     ts = ete3.TreeStyle()
    #     ts.layout_fn = layout
    #     ts.show_leaf_name = False
    #     ts.mode = 'r'
    #     ts.draw_guiding_lines = True
    #     ts.show_scale = False
    #     ts.scale = 25

    #     tre = ete3.Tree(self.trees.nhx)
    #     tre.ladderize()
    #     tre.convert_to_ultrametric(tree_length=len(tre)//2)
    #     tre.render(file_name=self.trees.pdf, h=40*len(tre), tree_style=ts)



    def save(self):
        """ save a JSON file representation of Quartet Class for checkpoint"""

        ## save each attribute as dict
        fulldumps = json.dumps(self.__dict__, 
                               sort_keys=False, 
                               indent=4, 
                               separators=(",", ":"),
                               )

        ## save to file, make dir if it wasn't made earlier
        assemblypath = os.path.join(self.dirs, self.name+".tet.json")
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)
    
        ## protect save from interruption
        done = 0
        while not done:
            try:
                with open(assemblypath, 'w') as jout:
                    jout.write(fulldumps)
                done = 1
            except (KeyboardInterrupt, SystemExit): 
                print('.')
                continue        



    def insert_to_array(self, start, results, bidx):
        """ inputs results from workers into hdf4 array """
        qrts, wgts, qsts = results

        with h5py.File(self.h5out, 'r+') as out:
            chunk = self.chunksize
            out['quartets'][start:start+chunk] = qrts
            out['weights'][start:start+chunk] = wgts

            if bidx:
                out["qboots/b{}".format(bidx-1)][start:start+chunk] = qsts
            else:
                out["qstats"][start:start+chunk] = qsts

        ## save checkpoint
        #data.svd.checkpoint_arr = np.where(ww == 0)[0].min()


    ########################################################################
    ## Main functions
    ########################################################################
    def run(self, force=0, quiet=0):
        """ 
        Run quartet inference on a SNP alignment. If checkpoint values exist 
        it will continue from where it left off unless force=True to force a
        a restart using the same parameter values. The analysis launches an
        ipcluster instance by default unless newclient=0, in which case it 
        expects to find an ipcluster instance running under the "default"
        profile. 
        """

        ## wrap everything in a try statement so we can ensure that it will
        ## save if interrupted and we will clean up the 
        ## client instance at the end. If it was created then we kill it. If
        ## using an existing client then we simply clean the memory space.
        try:
            ## find an ipcluster instance
            ipyclient = ip.core.parallel.get_client(**self._ipcluster)

            ## print a message about the cluster status
            ## if MPI setup then we are going to wait until all engines are
            ## ready so that we can print how many cores started on each 
            ## host machine exactly. 
            if not quiet:
                if self._ipcluster["engines"] == "MPI":
                    hosts = ipyclient[:].apply_sync(socket.gethostname)
                    for hostname in set(hosts):
                        print("  host compute node: [{} cores] on {}"\
                              .format(hosts.count(hostname), hostname))
                    print("")
                ## if Local setup then we know that we can get all the cores for 
                ## sure and we won't bother waiting for them to start, since 
                ## they'll start grabbing jobs once they're started. 
                else:
                    _cpus = min(ip.assemble.util.detect_cpus(), 
                                self._ipcluster["cores"])
                    print("  local compute node: [{} cores] on {}\n"\
                          .format(_cpus, socket.gethostname()))

            ## run the full inference or print finished prog bar if it's done
            if not self.checkpoint.boots:
                print("  inferring {} x 3 induced quartet trees".format(self.nquartets))
                self.inference(0, ipyclient)

            ## run the bootstrap replicates -------------------------------
            if self.nboots:
                if self.nboots == self.checkpoint.boots:
                    print("  continuing {} bootstrap replicates from boot {}"\
                          .format(self.nboots, self.checkpoint.boots))  
                else:
                    print("  running {} bootstrap replicates".format(self.nboots))              
    
                ## load from current boot
                for bidx in xrange(self.checkpoint.boots+1, self.nboots):
                    ## get resampled array and set checkpoint
                    if self.checkpoint.arr == 0:
                        if self.files.mapfile:
                            self.sample_bootseq_array_map()
                        else:
                            self.sample_bootseq_array() 

                    ## start boot inference, (1-indexed !!!)
                    self.inference(bidx, ipyclient)
                    self.checkpoint.boots = bidx

                ## write outputs with bootstraps
                self.write_output_splash()

            else:
                ## write outputs without bootstraps
                self.write_output_splash()


        ## handle exceptions so they will be raised after we clean up below
        except KeyboardInterrupt as inst:
            LOGGER.info("assembly interrupted by user.")
            print("\n  Keyboard Interrupt by user. Cleaning up...")

        except IPyradWarningExit as inst:
            LOGGER.info("IPyradWarningExit: %s", inst)
            print("  IPyradWarningExit: {}".format(inst))

        except Exception as inst:
            LOGGER.info("caught an unknown exception %s", inst)
            print("\n  Exception found: {}".format(inst))

        ## close client when done or interrupted
        finally:
            try:
                ## save the Assembly
                self.save()                
                
                ## can't close client if it was never open
                if ipyclient:

                    ## if CLI, stop jobs and shutdown
                    if self._ipcluster["cluster_id"]:
                        ipyclient.abort()
                        ipyclient.close()
                    ## if API, stop jobs and clean queue
                    else:
                        ipyclient.abort()
                        ipyclient.purge_everything()
            
            ## if exception is close and save, print and ignore
            except Exception as inst2:
                LOGGER.error("shutdown warning: %s", inst2)



    def inference(self, bidx, ipyclient):
        """ 
        Inference sends slices of jobs to the parallel engines for computing
        and collects the results into the output hdf5 array as they finish. 
        """

        ## an iterator to distribute sampled quartets in chunks
        njobs = sum(1 for _ in \
                xrange(self.checkpoint.arr, self.nquartets, self.chunksize))
        jobiter = iter(xrange(self.checkpoint.arr, self.nquartets, self.chunksize))
        LOGGER.info("chunksize: %s, start: %s, total: %s, njobs: %s", \
                self.chunksize, self.checkpoint.arr, self.nquartets, njobs)

        ## if bootstrap create an output array for results unless we are 
        ## restarting an existing analysis, then use the one already present
        with h5py.File(self.h5out, 'r+') as out:
            if 'b{}'.format(bidx) not in out["qboots"].keys():
                out["qboots"].create_dataset("b{}".format(bidx), 
                                             (self.nquartets, 4), 
                                             dtype=np.uint32, 
                                             chunks=(self.chunksize, 4))

        ## a distributor for engine jobs
        lbview = ipyclient.load_balanced_view()

        ## the three indexed resolutions of each quartet
        tests = np.array([[0, 1, 2, 3], 
                          [0, 2, 1, 3], 
                          [0, 3, 1, 2]], dtype=np.uint8)

        ## start progress bar timer and submit initial n jobs
        start = time.time()
        res = {}
        for _ in xrange(njobs):
            ## get chunk of quartet samples and send to a worker engine
            qidx = jobiter.next()
            LOGGER.info('qidx: %s', qidx)
            with h5py.File(self.h5in, 'r') as inh5:
                smps = inh5["samples"][qidx:qidx+self.chunksize]
            res[qidx] = lbview.apply(nworker, *[self, smps, tests])

        ## keep adding jobs until the jobiter is empty
        done = 0
        while 1:
            ## print progress unless bootstrapping, diff progbar for that.
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            if not bidx:
                progressbar(njobs, done, " initial tree | {}".format(elapsed))
            else:
                progressbar(njobs, done, " boot {:<7} | {}".format(bidx, elapsed))

            ## check for finished jobs
            curkeys = res.keys()
            finished = [i.ready() for i in res.values()]

            ## remove finished and submit new jobs
            if any(finished):
                for ikey in curkeys:
                    if res[ikey].ready():
                        if res[ikey].successful():
                            LOGGER.info("cleanup key %s", ikey)
                            ## track finished
                            done += 1
                            ## insert results into hdf5 data base
                            results = res[ikey].get(0)
                            self.insert_to_array(ikey, results, bidx)
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
                        ## send chunk off to be worked on
                        qidx = jobiter.next()
                        with h5py.File(self.h5in, 'r') as inh5:
                            smps = inh5["samples"][qidx:qidx+self.chunksize]
                        res[qidx] = lbview.apply(nworker, *[self, smps, tests])

                    ## if no more jobs then just wait until these are done
                    except StopIteration:
                        continue
            else:
                time.sleep(0.01)

            ## done is counted on finish, so this means we're done
            if njobs == done:
                break

        ## final progress bar
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        if not bidx:
            progressbar(njobs, done, " initial tree | {}".format(elapsed))
        else:
            progressbar(njobs, done, " boot {:<7} | {}".format(bidx, elapsed))
        print("")

        ## convert to txt file for wQMC
        self.dump_qmc()

        ## send to qmc
        if not bidx:
            self.run_qmc(0)
            #lbview.apply(self.run_qmc, 0)
        else:
            self.run_qmc(1)            
            #lbview.apply(self.run_qmc, 1)

        ## reset the checkpoint_arr
        self.checkpoint.arr = 0




#########################################################################
## Sampling functions 
#########################################################################
MUL = lambda x, y: x*y

## FROM THE ITERTOOLS RECIPES COOKCOOK
## TODO: replace random.sample with numpy random so that 
## our random seed stays the same
def random_combination(iterable, nquartets):
    """
    Random selection from itertools.combinations(iterable, r). 
    Use this if not sampling all possible quartets.
    """
    pool = tuple(iterable)
    size = len(pool)
    indices = random.sample(xrange(size), nquartets)
    return tuple(pool[i] for i in indices)



def random_product(iter1, iter2):
    """ random sampler for equal_splits func"""
    pool1 = tuple(iter1)
    pool2 = tuple(iter2)
    ind1 = random.sample(pool1, 2)
    ind2 = random.sample(pool2, 2)
    return tuple(ind1+ind2)



def n_choose_k(n, k):
    """ calculate the number of quartets as n-choose-k. This is used
    in equal splits to decide whether a split should be exhaustively sampled
    or randomly sampled. Edges near tips can be exhaustive while highly nested
    edges probably have too many quartets
    """
    return int(reduce(MUL, (Fraction(n-i, i+1) for i in range(k)), 1))


#############################################################################
#############################################################################
## SVD computation functions
## 
#############################################################################
#############################################################################

@numba.jit('f8(f8[:])', nopython=True)
def get_weights(scores):
    """ 
    calculates quartet weights from ordered svd scores. Following 
    description from Avni et al. 
    """
    ## lowest to highest [best, ils1, ils2]
    scores.sort()
    ## calculate weight given the svd scores
    if scores[2]:
        weight = (scores[2]-scores[0]) / \
                 (np.exp(scores[2]-scores[1]) * scores[2])
    else:
        weight = 0
    return weight



@numba.jit('u4[:](u4[:,:])', nopython=True)
def count_snps(mat):
    """ 
    calculate dstats from the count array and return as a float tuple 
    """

    ## get [aabb, baba, abba, aaab] 
    snps = np.zeros(4, dtype=np.uint32)

    ## get concordant (aabb) pis sites
    snps[0] = np.uint32(\
           mat[0, 5] + mat[0, 10] + mat[0, 15] + \
           mat[5, 0] + mat[5, 10] + mat[5, 15] + \
           mat[10, 0] + mat[10, 5] + mat[10, 15] + \
           mat[15, 0] + mat[15, 5] + mat[15, 10])

    ## get discordant (baba) sites
    for i in range(16):
        if i % 5:
            snps[1] += mat[i, i]
    
    ## get discordant (abba) sites
    snps[2] = mat[1, 4] + mat[2, 8] + mat[3, 12] +\
              mat[4, 1] + mat[6, 9] + mat[7, 13] +\
              mat[8, 2] + mat[9, 6] + mat[11, 14] +\
              mat[12, 3] + mat[13, 7] + mat[14, 11]

    ## get autapomorphy sites
    snps[3] = (mat.sum() - np.diag(mat).sum()) - snps[2]

    return snps


        
## TODO: either use pure numpy here or guvectorized func
@numba.jit('u1[:,:](u1[:,:],b1[:],u4[:])', nopython=True)
def subsample_snps(seqchunk, nmask, maparr):
    """ 
    removes ncolumns from snparray prior to matrix calculation, and 
    subsamples 'linked' snps (those from the same RAD locus) such that
    for these four samples only 1 SNP per locus is kept. This information
    comes from the 'map' array (map file). 
    """
    ## mask columns that contain Ns
    rmask = np.ones(seqchunk.shape[1], dtype=np.bool_)
    #LOGGER.info("rmask : %s %s", rmask.shape, rmask.sum())

    for idx in xrange(rmask.shape[0]):
        if nmask[idx]: 
            rmask[idx] = False
    #LOGGER.info("rmasked : %s %s", rmask.shape, rmask.sum())    

    ## apply mask
    newarr = seqchunk[:, rmask]

    ## return smaller Nmasked array
    return newarr



@numba.jit('u1[:,:](u1[:,:],b1[:],u4[:])', nopython=True)
def subsample_snps_map(seqchunk, nmask, maparr):
    """ 
    removes ncolumns from snparray prior to matrix calculation, and 
    subsamples 'linked' snps (those from the same RAD locus) such that
    for these four samples only 1 SNP per locus is kept. This information
    comes from the 'map' array (map file). 
    """
    ## mask columns that contain Ns
    rmask = np.ones(seqchunk.shape[1], dtype=np.bool_)

    ## apply mask to the mapfile
    last_snp = 0
    for idx in xrange(rmask.shape[0]):
        if nmask[idx]:
            ## mask if Ns
            rmask[idx] = False
        else:
            ## also mask if SNP already sampled 
            this_snp = maparr[idx]
            if maparr[idx] == last_snp:
                rmask[idx] = False
            ## record this snp
            last_snp = this_snp  
    
    ## apply mask
    newarr = seqchunk[:, rmask]
    
    ## return smaller Nmasked array
    return newarr



@numba.jit('u4[:,:,:](u1[:,:])', nopython=True)
def chunk_to_matrices(narr):
    """ 
    numba compiled code to get matrix fast.
    arr is a 4 x N seq matrix converted to np.int8
    I convert the numbers for ATGC into their respective index for the MAT
    matrix, and leave all others as high numbers, i.e., -==45, N==78. 
    """

    ## get seq alignment and create an empty array for filling
    mats = np.zeros((3, 16, 16), dtype=np.uint32)        

    ## replace ints with small ints that index their place in the 
    ## 16x16. If not replaced, the existing ints are all very large
    ## and the column will be excluded.
    for x in xrange(narr.shape[1]):
        i = narr[:, x]
        if np.sum(i) < 16:
            mats[0, i[0]*4:(i[0]+4)*4]\
                    [i[1]]\
                    [i[2]*4:(i[2]+4)*4]\
                    [i[3]] += 1
                
    ## get matrix 2
    mats[1, 0:4, 0:4] = mats[0, 0].reshape(4, 4)
    mats[1, 0:4, 4:8] = mats[0, 1].reshape(4, 4)
    mats[1, 0:4, 8:12] = mats[0, 2].reshape(4, 4)
    mats[1, 0:4, 12:16] = mats[0, 3].reshape(4, 4)
    mats[1, 4:8, 0:4] = mats[0, 4].reshape(4, 4)
    mats[1, 4:8, 4:8] = mats[0, 5].reshape(4, 4)
    mats[1, 4:8, 8:12] = mats[0, 6].reshape(4, 4)
    mats[1, 4:8, 12:16] = mats[0, 7].reshape(4, 4)
    mats[1, 8:12, 0:4] = mats[0, 8].reshape(4, 4)
    mats[1, 8:12, 4:8] = mats[0, 9].reshape(4, 4)
    mats[1, 8:12, 8:12] = mats[0, 10].reshape(4, 4)
    mats[1, 8:12, 12:16] = mats[0, 11].reshape(4, 4)
    mats[1, 12:16, 0:4] = mats[0, 12].reshape(4, 4)
    mats[1, 12:16, 4:8] = mats[0, 13].reshape(4, 4)
    mats[1, 12:16, 8:12] = mats[0, 14].reshape(4, 4)
    mats[1, 12:16, 12:16] = mats[0, 15].reshape(4, 4)
    
    ## get matrix 3
    mats[2, 0:4, 0:4] = mats[0, 0].reshape(4, 4).T
    mats[2, 0:4, 4:8] = mats[0, 1].reshape(4, 4).T
    mats[2, 0:4, 8:12] = mats[0, 2].reshape(4, 4).T
    mats[2, 0:4, 12:16] = mats[0, 3].reshape(4, 4).T
    mats[2, 4:8, 0:4] = mats[0, 4].reshape(4, 4).T
    mats[2, 4:8, 4:8] = mats[0, 5].reshape(4, 4).T
    mats[2, 4:8, 8:12] = mats[0, 6].reshape(4, 4).T
    mats[2, 4:8, 12:16] = mats[0, 7].reshape(4, 4).T
    mats[2, 8:12, 0:4] = mats[0, 8].reshape(4, 4).T
    mats[2, 8:12, 4:8] = mats[0, 9].reshape(4, 4).T
    mats[2, 8:12, 8:12] = mats[0, 10].reshape(4, 4).T
    mats[2, 8:12, 12:16] = mats[0, 11].reshape(4, 4).T
    mats[2, 12:16, 0:4] = mats[0, 12].reshape(4, 4).T
    mats[2, 12:16, 4:8] = mats[0, 13].reshape(4, 4).T
    mats[2, 12:16, 8:12] = mats[0, 14].reshape(4, 4).T
    mats[2, 12:16, 12:16] = mats[0, 15].reshape(4, 4).T  
                
    return mats




@numba.jit(nopython=True)
def calculate(seqnon, tests):
    """ groups together several numba compiled funcs """

    ## create empty matrices
    #LOGGER.info("tests[0] %s", tests[0])
    #LOGGER.info('seqnon[[tests[0]]] %s', seqnon[[tests[0]]])
    mats = chunk_to_matrices(seqnon[tests[0]])

    ## epmty svdscores for each arrangement of seqchunk
    qscores = np.zeros(3, dtype=np.float64)

    for test in range(3):
        ## get svd scores
        tmpscore = np.linalg.svd(mats[test].astype(np.float64))[1]
        #qscores[test] = np.sqrt(tmpscore[11:]).sum()
        qscores[test] = np.sqrt(np.sum(tmpscore[11:]**2))

    ## sort to find the best qorder
    best = np.where(qscores == qscores.min())[0]
    bidx = tests[best][0]
    qsnps = count_snps(mats[best][0])

    # LOGGER.info("""
    #     best: %s, 
    #     bidx: %s, 
    #     qscores: %s, 
    #     qsnps: %s
    #     mats \n %s
    #     """, best, bidx, qscores, qsnps, mats)

    return bidx, qscores, qsnps



def nworker(data, smpchunk, tests):
    """ The workhorse function. Not numba. """

    ## open the seqarray view, the modified array is in bootsarr
    inh5 = h5py.File(data.h5in, 'r')
    seqview = inh5["bootsarr"][:]

    ## choose function based on mapfile arg
    if data.files.mapfile:
        subsample = subsample_snps_map
        maparr = inh5["bootsmap"][:]
    else:
        subsample = subsample_snps
        maparr = np.zeros((2, 2), dtype=np.uint32)

    ## create an N-mask array of all seq cols
    nall_mask = seqview[:] == 78

    ## tried numba compiling everythign below here, but was not faster
    ## than making nmask w/ axis arg in numpy

    ## get the input arrays ready
    rquartets = np.zeros((smpchunk.shape[0], 4), dtype=np.uint16)
    rweights = np.zeros(smpchunk.shape[0], dtype=np.float64)
    rdstats = np.zeros((smpchunk.shape[0], 4), dtype=np.uint32)

    ## record how many quartets have no information
    excluded = 0

    ## fill arrays with results using numba funcs
    for idx in xrange(smpchunk.shape[0]):
        ## get seqchunk for 4 samples (4, ncols) 
        sidx = smpchunk[idx]
        seqchunk = seqview[sidx]

        ## get N-containing columns in 4-array
        nmask = nall_mask[sidx].sum(axis=0, dtype=np.bool_)
        #LOGGER.info('not N-masked sites: %s', nmask.sum())

        ## remove Ncols from seqchunk & sub-sample unlinked SNPs
        #LOGGER.info("seqchunk %s", seqchunk.shape)
        seqnon = subsample(seqchunk, nmask, maparr[:, 0])
        #LOGGER.info("seqnon sites %s", seqnon.shape)
        #LOGGER.info("before sub: %s, after %s", seqchunk.shape, seqnon.shape)

        ## get matrices if there are any shared SNPs
        if seqnon.shape[1]:
            ## returns best-tree index, qscores, and qstats
            bidx, qscores, qstats = calculate(seqnon, tests)

            ## get weights from the three scores sorted. 
            ## Only save to file if the quartet has information
            rdstats[idx] = qstats 
            
            iwgt = get_weights(qscores)
            if iwgt:
                rweights[idx] = iwgt
                rquartets[idx] = smpchunk[idx][bidx]
                LOGGER.info("""\n
                    ------------------------------------
                    bidx: %s
                    qstats: %s, 
                    weight: %s, 
                    scores: %s
                    ------------------------------------
                    """,
                    bidx, qstats, rweights[idx], qscores)
            else:
                excluded += 1

    LOGGER.warning("excluded quartets %s", excluded)    
    #return 
    return rquartets, rweights, rdstats 




########################################################################
## GLOBALS
########################################################################

MIDSTREAM_MESSAGE = """
    loaded object method={}
    cannot change sampling methods midstream
    use force argument to start new run with new method
"""
## require method to be same as loaded type
## assert method == data..method, MIDSTREAM_MESSAGE.format(method)

LOADING_MESSAGE = """\
  Continuing checkpointed analysis: {}
    sampling method: {}
    bootstrap checkpoint: {}
    array checkpoint: {}
"""

LOADING_RANDOM = """\
    loading {} random quartet samples to infer a starting tree 
    inferring {} x 3 quartet trees
"""

LOADING_STARTER = """\
    loading {} equal-splits quartets from starting tree
    """

NO_SNP_FILE = """\
    Cannot find SNP file. You entered: '{}'. 
    """


AMBIGS = {"R":("G", "A"),
          "K":("G", "T"),
          "S":("G", "C"),
          "Y":("T", "C"),
          "W":("T", "A"),
          "M":("C", "A")}

## convenience functions
def opr(path):
    """ shorthand for realpath """
    return os.path.realpath(path)


@numba.jit(nopython=True)
def shuffle_cols(seqarr, newarr, cols):
    """ used in bootstrap resampling without a map file """
    for idx in xrange(cols.shape[0]):
        newarr[:, idx] = seqarr[:, cols[idx]]
    return newarr


## TODO: this could be numbified by making AMBIGS into two arrays
## IF SO, pay attention to the numba random seed being different from numpy
def resolve_ambigs(tmpseq):
    """ returns a seq array with 'RSKYWM' randomly replaced with resolved bases"""
    ## iterate over the bases 'RSKWYM': [82, 83, 75, 87, 89, 77]
    for ambig in np.uint8([82, 83, 75, 87, 89, 77]):
        ## get all site in this ambig
        idx, idy = np.where(tmpseq == ambig)
        ## get the two resolutions of the ambig
        res1, res2 = AMBIGS[ambig.view("S1")]
        ## randomly sample half those sites
        halfmask = np.random.choice([True, False], idx.shape[0])
        ## replace ambig bases with their resolutions
        for i in xrange(halfmask.shape[0]):
            if halfmask[i]:
                tmpseq[idx[i], idy[i]] = np.array(res1).view(np.uint8)
            else:
                tmpseq[idx[i], idy[i]] = np.array(res2).view(np.uint8)
    return tmpseq



@numba.jit(nopython=True)
def get_spans(maparr, spans):
    """ get span distance for each locus in original seqarray """
    ## start at 0, finds change at 1-index of map file
    bidx = 0
    
    ## read through marr and record when locus id changes
    for idx in xrange(maparr.shape[0]):
        cur = maparr[idx, 0]
        if cur != bidx:
            spans[cur-1, 1] = idx+1
            spans[cur, 0] = idx+1
    return spans



@numba.jit(nopython=True)
def get_shape(spans, loci):
    """ get shape of new bootstrap resampled locus array """
    width = 0
    for idx in xrange(loci.shape[0]):
        width += spans[loci[idx], 1] - spans[loci[idx], 0]
    return width
    


@numba.jit(nopython=True)
def fill_boot(seqarr, newboot, newmap, spans, loci):
    """ fills the new bootstrap resampled array """
    ## column index
    cidx = 0
  
    ## resample each locus
    for i in xrange(loci.shape[0]):
        
        ## grab a random locus's columns
        x1 = spans[loci[i]][0]
        x2 = spans[loci[i]][1]
        cols = seqarr[:, x1:x2]

        ## randomize columns within colsq
        cord = np.random.choice(cols.shape[1], cols.shape[1], replace=False)
        rcols = cols[:, cord]
        
        ## fill bootarr with n columns from seqarr
        ## the required length was already measured
        newboot[:, cidx:cidx+cols.shape[1]] = rcols

        ## fill bootmap with new map info
        newmap[cidx: cidx+cols.shape[1], 0] = i+1
        
        ## advance column index
        cidx += cols.shape[1]

    ## return the concatenated cols
    return newboot, newmap



## thoughts on this... we need to re-estimate chunksize given whatever
## new parallel setup is passed in. Start from arr checkpoint to end. 
def load_json(path):
    """ Load a json serialized Quartet Class object """

    ## load the JSON string and try with name+.json
    if not path.endswith(".tet.json"):
        path += ".tet.json"

    ## expand user
    path = path.replace("~", os.path.expanduser("~"))

    ## load the json file
    try:
        with open(path, 'r') as infile:
            fullj = _byteify(json.loads(infile.read(),
                            object_hook=_byteify), 
                        ignore_dicts=True)
    except IOError:
        raise IPyradWarningExit("""\
    Cannot find checkpoint (.test.json) file at: {}""".format(path))

    ## create a new Quartet Class
    newobj = Quartet(fullj["name"], fullj["dirs"], fullj["method"])

    ## fill in the same attributes
    for key in fullj:
        newobj.__setattr__(key, fullj[key])

    newobj.files = ObjDict(newobj.files)
    newobj.trees = ObjDict(newobj.trees)
    newobj.checkpoint = ObjDict(newobj.checkpoint)

    return newobj



def _byteify(data, ignore_dicts=False):
    """
    converts unicode to utf-8 when reading in json files
    """
    if isinstance(data, unicode):
        return data.encode("utf-8")

    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]

    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
        }
    return data





########################################################################
## Plotting functions
########################################################################


def consensus_tree(trees, names=None, cutoff=0.0):
    """ 
    An extended majority rule consensus function for ete3. 
    Modelled on the similar function from scikit-bio tree module. If 
    cutoff=0.5 then it is a normal majority rule consensus, while if 
    cutoff=0.0 then subsequent non-conflicting clades are added to the tree.
    """

    ## find which clades occured with freq > cutoff
    namedict, clade_counts = _find_clades(trees, names=names)

    ## filter out the < cutoff clades
    fclade_counts = _filter_clades(clade_counts, cutoff)

    ## build tree
    consens_tree, _ = _build_trees(fclade_counts, namedict)
    ## make sure no singleton nodes were left behind
    return consens_tree, clade_counts



def _filter_clades(clade_counts, cutoff):
    """ 
    A subfunc of consensus_tree(). Removes clades that occur 
    with freq < cutoff.
    """

    ## store clades that pass filter
    passed = []
    clades = np.array([list(i[0]) for i in clade_counts], dtype=np.int8)
    counts = np.array([i[1] for i in clade_counts], dtype=np.float64)
    
    for idx in xrange(clades.shape[0]):
        conflict = False
    
        if counts[idx] < cutoff:
            continue
            
        if np.sum(clades[idx]) > 1:
            # check the current clade against all the accepted clades to see if
            # it conflicts. A conflict is defined as:
            # 1. the clades are not disjoint
            # 2. neither clade is a subset of the other
            # OR:
            # 1. it is inverse of clade (affects only <fake> root state)
            # because at root node it mirror images {0011 : 95}, {1100 : 5}.
            for aidx in passed:
                #intersect = clade.intersection(accepted_clade)
                summed = clades[idx] + clades[aidx]
                intersect = np.max(summed) > 1
                subset_test0 = np.all(clades[idx] - clades[aidx] >= 0)
                subset_test1 = np.all(clades[aidx] - clades[idx] >= 0)
                #invert_test = np.bool_(clades[aidx]) != np.bool_(clades[idx])

                #if np.all(invert_test):
                #    counts[aidx] += counts[idx]
                #    conflict = True
                if intersect:
                    if (not subset_test0) and (not subset_test1):
                        conflict = True

        if conflict == False:
            passed.append(idx)

    ## rebuild the dict
    rclades = []#j for i, j in enumerate(clade_counts) if i in passed]
    ## set the counts to include mirrors
    for idx in passed:
        rclades.append((clades[idx], counts[idx]))
    return rclades



def _find_clades(trees, names):
    """ 
    A subfunc of consensus_tree(). Traverses trees to count clade occurrences.
    Names are ordered by names, else they are in the order of the first
    tree. 
    """
    ## index names from the first tree
    if not names:
        names = trees[0].get_leaf_names()
    ndict = {j:i for i, j in enumerate(names)}
    namedict = {i:j for i, j in enumerate(names)}

    ## store counts
    clade_counts = defaultdict(int)
    ## count as bitarray clades in each tree
    for tree in trees:
        tree.unroot()
        for node in tree.traverse('postorder'):
            bits = bitarray('0'*len(tree))
            for child in node.iter_leaf_names():
                bits[ndict[child]] = 1
            ## if parent is root then mirror flip one child (where bit[0]=0)
            if not node.is_root():
                if node.up.is_root():
                    if bits[0]:
                        bits.invert()
            clade_counts[bits.to01()] += 1

    ## convert to freq
    for key, val in clade_counts.items():
        clade_counts[key] = val / float(len(trees))

    ## return in sorted order
    clade_counts = sorted(clade_counts.items(), 
                          key=lambda x: x[1],
                          reverse=True)
    return namedict, clade_counts



def _build_trees(fclade_counts, namedict):
    """ 
    A subfunc of consensus_tree(). Build an unrooted consensus tree 
    from filtered clade counts. 
    """

    ## storage
    nodes = {}
    idxarr = np.arange(len(fclade_counts[0][0]))
    queue = []

    ## create dict of clade counts and set keys
    countdict = defaultdict(int)
    for clade, count in fclade_counts:
        mask = np.int_(list(clade)).astype(np.bool)
        ccx = idxarr[mask]
        queue.append((len(ccx), frozenset(ccx)))
        countdict[frozenset(ccx)] = count

    while queue:
        queue.sort()
        (clade_size, clade) = queue.pop(0)
        new_queue = []
    
        # search for ancestors of clade
        for (_, ancestor) in queue:
            if clade.issubset(ancestor):
                # update ancestor such that, in the following example:
                # ancestor == {1, 2, 3, 4}
                # clade == {2, 3}
                # new_ancestor == {1, {2, 3}, 4}
                new_ancestor = (ancestor - clade) | frozenset([clade])          
                countdict[new_ancestor] = countdict.pop(ancestor)
                ancestor = new_ancestor
            
            new_queue.append((len(ancestor), ancestor))
   
        # if the clade is a tip, then we have a name
        if clade_size == 1:
            name = list(clade)[0]
            name = namedict[name]
        else:
            name = None 
        
        # the clade will not be in nodes if it is a tip
        children = [nodes.pop(c) for c in clade if c in nodes]
        node = ete3.Tree(name=name)    
        for child in children:
            node.add_child(child)
        if not node.is_leaf():
            node.dist = int(round(100*countdict[clade]))
        else:
            node.dist = int(0) 
        
        nodes[clade] = node
        queue = new_queue
    tre = nodes.values()[0]
    tre.unroot()
    ## return the tree and other trees if present
    return tre, list(nodes.values())

    

def _get_total(tots, node):
    """ get total number of quartets possible for a split """
    down = set(node.get_leaf_names())
    up = tots - down
    return n_choose_k(len(down), 2) * n_choose_k(len(up), 2)
    

    
def _get_sampled(data, tots, node, names):
    """ get how many quartets were sampled that are informative for a split"""
    ## get leaves up and down
    down = set(node.get_leaf_names())
    up = tots - down

    ## get up and down as index
    idxd = set([names.index(i) for i in down])
    idxu = set([names.index(i) for i in up])

    ## find how many sampled quartets span each edge
    sampled = 0

    ## do chunks at a time in case qrts is huge
    idx = 0
    with h5py.File(data.h5out, 'r') as io5:
        qrts = io5["quartets"][idx:idx+data.chunksize]
        for qrt in qrts:
            sqrt = set(qrt)
            if len(sqrt.intersection(idxd)) > 1:
                if len(sqrt.intersection(idxu)) > 1:
                    sampled += 1
            idx += data.chunksize

    return sampled



## GLOBALS #############################################################

STATSOUT = """
  Statistics for sampling, discordance, and tree support:
    > {}
    """

FINALTREES = """\
  Full tree inferred from by weighted quartet-joining of the SNP supermatrix
    > {}
    """

BOOTTREES = """\
  Extended majority-rule consensus with support as edge lengths:
    > {}

  All bootstrap trees:
    > {}
    """

ASCII_TREE = """\
  ASCII view of unrooted topology from the weighted analysis
    {}
    """

LINKS = """\
  * For tips on plotting trees in R: {}     
  * For tips on citing this software: {} 
    """

########################################################################
## JUNK

        # ## get unrooted weighted quartets tree
        # wtre = ete3.Tree(self.trees.wtre, format=0)
        # wtre.unroot()
        # for node in wtre.traverse():
        #     node.add_feature("bootstrap", 0)
        #     node.add_feature("quartets_total", _get_total(wtre, node))
        #     node.add_feature("quartets_sampled", _get_sampled(self, wtre, node, names))
        #     try:
        #         prop = 100*(float(node.quartets_sampled) / node.quartets_total)
        #     except ZeroDivisionError:
        #         prop = 0.0
        #     node.add_feature("quartets_sampled_prop", prop)
        #     node.dist = 0
        #     node.support = 0


        # ## get unrooted boot trees
        # if with_boots:
        #     oboots = open(self.trees.tboots, 'r').readlines()
        #     wboots = open(self.trees.wboots, 'r').readlines()
        #     oboots = [ete3.Tree(btre.strip()) for btre in oboots]
        #     wboots = [ete3.Tree(btre.strip()) for btre in wboots]    
        #     _ = [btre.unroot() for btre in oboots]
        #     _ = [btre.unroot() for btre in wboots]

        #     ## get and set support values 
        #     for tre, boots in zip([otre, wtre], [oboots, wboots]):
        #         for btre in boots:
        #             common = tre.compare(btre, unrooted=True)
        #             for bnode in common["common_edges"]:
        #                 ## check monophyly of each side of split
        #                 a = tre.check_monophyly(bnode[0], target_attr='name', unrooted=True)
        #                 b = tre.check_monophyly(bnode[1], target_attr='name', unrooted=True)
        #                 ## if both sides are monophyletic
        #                 if a[0] and b[0]:
        #                     ## find which is the 'bottom' node, to attach support to
        #                     node = list(tre.get_monophyletic(bnode[0], target_attr='name'))
        #                     node.extend(list(tre.get_monophyletic(bnode[1], target_attr='name')))
        #                     ## add +1 suport to (edge dist) to this edge
        #                     if not node[0].is_leaf():
        #                         node[0].dist += 1
        #                         node[0].support += 1
        #                         node[0].bootstrap += 1

        #     ## change support values to percentage
        #     for tre in [otre, wtre]:
        #         for node in tre.traverse():
        #             node.dist = int(100 * (node.dist / len(wboots)))
        #             node.support = int(100 * (node.support / len(wboots)))
        #             node.bootstrap = int(100 * (node.bootstrap / len(wboots)))

        #     ## return as newick string w/ support as edge labels (lengths)
        #     with open(self.trees.tbtre, 'w') as outtre:
        #         outtre.write(otre.write(format=5))

        #     with open(self.trees.wbtre, 'w') as outtre:
        #         outtre.write(wtre.write(format=5))
        #     features = ["bootstrap", "quartets_total", "quartets_sampled", "quartets_sampled_prop"]            
        # else:
        #     
        



if __name__ == "__main__":

    ## imports
    import ipyrad.analysis as ipa
    #import ipyrad as ip
    #import ipyparallel as ipp

    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/cli/cli.json")
    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/iptutorial/cli.json")
    ## run


