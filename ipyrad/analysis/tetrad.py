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
import copy
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
from ipyrad.assemble.util import IPyradWarningExit, progressbar

## when you have time go back and set attrubutes on toytrees
from toytree import ete3mini as ete3

## numba debugging
#numba.config.NUMBA_DEFAULT_NUM_THREADS = 1
#numba.config.NUMBA_DISABLE_JIT = 1

## set the logger
import logging
LOGGER = logging.getLogger(__name__)

## a reminder of the 16 x 16 matrix of site counts (The phylogenetic invariants). 
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
## Tetrad inference Class Object
#############################################################################
#############################################################################



class Params(object):
    """ 
    A dict-like object for storing params values with a custom repr
    """
    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())
        _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
        for key in keys:
            _repr += _printstr.format(key, str(self[key]))
        return _repr




class Tetrad(object):
    """
    The main tetrad object for storing data and checkpointing. It is 
    initialized with a name, and args to command line (e.g., sampling method, 
    starting tree, nboots, etc.). 
    """

    def __init__(self,
        name, 
        workdir="analysis-tetrad",
        seqfile=None, 
        mapfile=None, 
        method='all', 
        guidetreefile=None, 
        nboots=0, 
        nquartets=0, 
        resolve=True, 
        initarr=True, 
        load=False,
        quiet=False,
        *args, 
        **kwargs):

        ## interactive or CLI
        self.kwargs = kwargs

        ## name this assembly
        self.name = name
        self.dirs = os.path.realpath(workdir)
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)

        ## store default cluster information 
        self._ipcluster = {
            "cluster_id" : "",
            "profile" : "default", 
            "engines" : "Local", 
            "quiet" : 0, 
            "timeout" : 60, 
            "cores" : 0}

        ## Sampling method attributes
        self.method = method
        self.nboots = nboots
        self.nquartets = nquartets
        self._chunksize = 0
        self._resolve = resolve

        ## store samples from the seqarray
        self.samples = []

        ## self.populations ## if we allow grouping samples
        ## (haven't done this yet)

        ## hdf5 data bases init and delete existing
        self.database = Params()
        self.database.input = os.path.join(self.dirs, self.name+".input.h5")
        self.database.output = os.path.join(self.dirs, self.name+".output.h5")        

        ## input files
        self.files = Params()
        self.files.seqfile = seqfile
        self.files.mapfile = mapfile
        self.files.guidetreefile = guidetreefile
        self.files.qdump = None
        self.files.stats = None

        ## store tree file paths (init as None)
        self.trees = Params()
        self._tmp = None          ## the temporary qmc tree         
        self.trees.tree = None    #os.path.join(self.dirs, self.name+".full.tre")
        self.trees.cons = None    #os.path.join(self.dirs, self.name+".consensus.tre")
        self.trees.boots = None   #os.path.join(self.dirs, self.name+".boots")        
        self.trees.nhx = None     #os.path.join(self.dirs, self.name+".nhx.tre")     
        
        ## stats is written to os.path.join(self.dirs, self.name+".stats.txt")
        self.stats = Params()
        self.stats.n_quartets_sampled = self.nquartets
        #self.stats.prop_quartets_sampled = None

        ## checkpointing information
        self.checkpoint = Params()
        self.checkpoint.boots = 0
        self.checkpoint.arr = 0

        ## init the seq data and samples
        if load:
            self._load_json(self.name, self.dirs)

        elif seqfile:
            if initarr:
                self._init_seqarray(quiet=quiet)
                self._parse_names()
        else:
            raise IPyradWarningExit("must enter a seqfile argument.")

        ## if quartets not entered then sample all
        total = n_choose_k(len(self.samples), 4)
        if self.method != "all":
            if int(self.nquartets) >= total:
                self.method = "all"
                self.nquartets = total
                print("nquartets > total: switching to method='all' ")
            if not self.nquartets:
                raise IPyradWarningExit("must enter nquartets value w/ method='random'")
        else:
            self.nquartets = total
            
        

    def refresh(self):
        """ 
        Remove all existing results files and reinit the h5 arrays 
        so that the tetrad object is just like fresh from a CLI start.
        """

        ## clear any existing results files
        oldfiles = [self.files.qdump] + \
                    self.database.__dict__.values() + \
                    self.trees.__dict__.values()
        for oldfile in oldfiles:
            if oldfile:
                if os.path.exists(oldfile):
                    os.remove(oldfile)
        ## reinit
        self.__init__(
            name=self.name, 
            seqfile=self.files.seqfile, 
            mapfile=self.files.mapfile,
            workdir=self.dirs,
            method=self.method,
            guidetreefile=self.files.guidetreefile,
            resolve=self._resolve, 
            nboots=self.nboots, 
            nquartets=self.nquartets, 
            initarr=True, 
            quiet=True,
            )



    def _parse_names(self):
        ## parse samples from the sequence file
        self.samples = []
        with iter(open(self.files.seqfile, 'r')) as infile:
            infile.next().strip().split()
            while 1:
                try:
                    self.samples.append(infile.next().split()[0])
                except StopIteration:
                    break



    def _init_seqarray(self, quiet=False):
        """ 
        Fills the seqarr with the full data set, and creates a bootsarr copy
        with the following modifications:

        1) converts "-" into "N"s, since they are similarly treated as missing. 
        2) randomly resolve ambiguities (RSKWYM)
        3) convert to uint8 for smaller memory load and faster computation
        """

        ## read in the seqfile
        try:
            spath = open(self.files.seqfile, 'r')
        except IOError:
            raise IPyradWarningExit(NO_SNP_FILE.format(self.files.seqfile))
        line = spath.readline().strip().split()
        ntax = int(line[0])
        nbp = int(line[1])

        ## make a tmp seq array
        if not quiet:
            print("loading seq array [{} taxa x {} bp]".format(ntax, nbp))        
        tmpseq = np.zeros((ntax, nbp), dtype=np.uint8)
    
        ## create array storage for real seq and the tmp bootstrap seqarray
        with h5py.File(self.database.input, 'w') as io5:
            io5.create_dataset("seqarr", (ntax, nbp), dtype=np.uint8)
            io5.create_dataset("bootsarr", (ntax, nbp), dtype=np.uint8)
            io5.create_dataset("bootsmap", (nbp, 2), dtype=np.uint32)

            ## if there is a map file, load it into the bootsmap
            if self.files.mapfile:
                with open(self.files.mapfile, 'r') as inmap:
                    ## parse the map file from txt and save as dataset
                    maparr = np.genfromtxt(inmap, dtype=np.uint64)
                    io5["bootsmap"][:] = maparr[:, [0, 3]]

                    ## parse the span info from maparr and save to dataset
                    spans = np.zeros((maparr[-1, 0], 2), np.uint64)
                    spans = get_spans(maparr, spans)
                    io5.create_dataset("spans", data=spans)
                    if not quiet:
                        print("max unlinked SNPs per quartet (nloci): {}"\
                              .format(spans.shape[0]))
            else:
                io5["bootsmap"][:, 0] = np.arange(io5["bootsmap"].shape[0])

            ## fill the tmp array from the input phy
            for line, seq in enumerate(spath.readlines()):
                tmpseq[line] = np.array(list(seq.split()[-1])).view(np.uint8)

            ## convert '-' or '_' into 'N'
            tmpseq[tmpseq == 45] = 78
            tmpseq[tmpseq == 95] = 78            

            ## save array to disk so it can be easily accessed by slicing
            ## This unmodified array is used again later for sampling boots
            io5["seqarr"][:] = tmpseq

            ## resolve ambiguous IUPAC codes
            if self._resolve:
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



    def _sample_bootseq_array(self):
        ## Takes the seqarray and re-samples columns and saves to bootsarr. 
        ## use 'r+' to read and write to existing array
        with h5py.File(self.database.input, 'r+') as io5:        
            ## load in the seqarr and maparr
            seqarr = io5["seqarr"][:]

            ## resample columns with replacement
            newarr = np.zeros(seqarr.shape, dtype=np.uint8)
            cols = np.random.randint(0, seqarr.shape[1], seqarr.shape[1])
            tmpseq = shuffle_cols(seqarr, newarr, cols)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self._resolve:
                tmpseq = resolve_ambigs(tmpseq)
        
            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## fill the boot array with a re-sampled phy w/ replacement
            io5["bootsarr"][:] = tmpseq
            del tmpseq



    def _sample_bootseq_array_map(self):
        ## Re-samples loci with replacement to fill the bootarr
        with h5py.File(self.database.input, 'r+') as io5:
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
            if self._resolve:
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



    ## Functions to fill h5 with samples
    def _store_N_samples(self, ncpus):
        """ 
        Find all quartets of samples and store in a large array
        Create a chunk size for sampling from the array of quartets. 
        This should be relatively large so that we don't spend a lot of time 
        doing I/O, but small enough that jobs finish often for checkpointing.
        """
        breaks = 2
        if self.nquartets < 5000:
            breaks = 1
        if self.nquartets > 100000:
            breaks = 4
        if self.nquartets > 500000:
            breaks = 8

        ## chunk up the data
        self._chunksize = (self.nquartets // (breaks * ncpus) + \
                          (self.nquartets % (breaks * ncpus)))
        LOGGER.info("nquarts = %s, chunk = %s", self.nquartets, self._chunksize)

        ## 'samples' stores the indices of the quartet. 
        ## `quartets` stores the correct quartet in the order (1,2|3,4)
        ## `weights` stores the weight of the quartet in 'quartets'
        ## we gzip this for now, but check later if this has a big speed cost

        ## create h5 OUT empty arrays
        with h5py.File(self.database.output, 'w') as io5:
            io5.create_dataset("quartets", 
                               (self.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self._chunksize, 4))
            io5.create_dataset("qstats", 
                               (self.nquartets, 4), 
                               dtype=np.uint32, 
                               chunks=(self._chunksize, 4))
            io5.create_group("qboots")


        ## append to h5 IN array (which also has seqarray) and fill it
        with h5py.File(self.database.input, 'a') as io5:
            ## create data sets
            io5.create_dataset("samples", 
                               (self.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self._chunksize, 4),
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
                    while len(qiter) < min(self._chunksize, io5["samples"].shape[0]):
                        qiter.append(
                            random_combination(range(len(self.samples)), 4))
                    dat = np.array(qiter)
                else:
                    ## grab the next ordered chunksize
                    dat = np.array(list(itertools.islice(qiter, self._chunksize)))

                ## store to h5 
                io5["samples"][i:i+self._chunksize] = dat[:io5["samples"].shape[0] - i]
                i += self._chunksize



    ## NOT YET READY FOR USE (TESTING)
    def _store_equal_samples(self, ncpus):
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

        self._chunksize = (self.nquartets // (breaks * ncpus) + \
                         (self.nquartets % (breaks * ncpus)))
        LOGGER.info("nquarts = %s, chunk = %s", self.nquartets, self._chunksize)

        ## create h5 OUT empty arrays
        with h5py.File(self.database.output, 'w') as io5:
            io5.create_dataset("quartets", 
                               (self.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self._chunksize, 4))
            io5.create_dataset("qstats", 
                               (self.nquartets, 4), 
                               dtype=np.uint32, 
                               chunks=(self._chunksize, 4))
            io5.create_group("qboots")

        ## get starting tree, unroot, randomly resolve, ladderize
        tre = ete3.Tree(self.files.guidetreefile, format=0)
        #tre = toytree.tree(self.files.guidetreefile, format=0)
        tre.tree.unroot()
        tre.tree.resolve_polytomy(recursive=True)
        tre.tree.ladderize()

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
        with h5py.File(self.database.input, 'a') as io5:
            ## create data sets
            io5.create_dataset("samples", 
                               (self.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self._chunksize, 4),
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
                while len(qdat) < self._chunksize:
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
                while len(qdat) < self._chunksize:
                    qrtsamp = random_combination(range(len(self.samples)), 4)
                    if tuple(qrtsamp) not in sampled:
                        qdat.append(qrtsamp)
                        sampled.add(qrtsamp)
                        random_target += 1

                ## stick chunk into h5 array
                dat = np.array(qdat, dtype=np.uint16)
                io5["samples"][i:i+self._chunksize] = dat[:io5["samples"].shape[0] - i]
                i += self._chunksize

            print("  equal sampling: {} edge quartets, {} random quartets "\
                  .format(edge_targeted, random_target))



    def _run_qmc(self, boot):
        """ runs quartet max-cut on a quartets file """

        ## convert to txt file for wQMC
        self._tmp = os.path.join(self.dirs, ".tmpwtre")
        cmd = [ip.bins.qmc, "qrtt="+self.files.qdump, "otre="+self._tmp] 

        ## run them
        proc = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        res = proc.communicate()
        if proc.returncode:
            #LOGGER.error("Error in QMC: \n({}).".format(res))
            LOGGER.error(res)
            raise IPyradWarningExit(res[1])

        ## read in the tmp files since qmc does not pipe
        with open(self._tmp) as intree:
            ## convert int names back to str names renamer returns a newick str
            #tmp = toytree.tree(intree.read().strip())
            tmp = ete3.Tree(intree.read().strip())
            tmpwtre = self._renamer(tmp)#.tree)

        ## save the tree
        if boot:
            self.trees.boots = os.path.join(self.dirs, self.name+".boots")
            with open(self.trees.boots, 'a') as outboot:
                outboot.write(tmpwtre+"\n")
        else:
            self.trees.tree = os.path.join(self.dirs, self.name+".tree")
            with open(self.trees.tree, 'w') as outtree:
                outtree.write(tmpwtre)

        ## save JSON file checkpoint
        self._save()



    def _dump_qmc(self):
        """ 
        Makes a reduced array that excludes quartets with no information and 
        prints the quartets and weights to a file formatted for wQMC 
        """
        ## open the h5 database
        io5 = h5py.File(self.database.output, 'r')

        ## create an output file for writing
        self.files.qdump = os.path.join(self.dirs, self.name+".quartets.txt")
        LOGGER.info("qdump file %s", self.files.qdump)
        outfile = open(self.files.qdump, 'w')

        ## todo: should pull quarts order in randomly? or doesn't matter?
        for idx in xrange(0, self.nquartets, self._chunksize):
            ## get mask of zero weight quartets
            #mask = io5["weights"][idx:idx+self.chunksize] != 0
            #weight = io5["weights"][idx:idx+self.chunksize][mask]
            #LOGGER.info("exluded = %s, mask shape %s", 
            #            self._chunksize - mask.shape[0], mask.shape)
            #LOGGER.info('q shape %s', io5["quartets"][idx:idx+self._chunksize].shape)
            masked_quartets = io5["quartets"][idx:idx+self._chunksize, :]#[mask, :]
            quarts = [list(j) for j in masked_quartets]

            ## format and print
            #chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
            #                                        in zip(quarts, weight)]
            chunk = ["{},{}|{},{}".format(*i) for i in quarts]
            outfile.write("\n".join(chunk)+"\n")


        ## close output file and h5 database
        outfile.close()
        io5.close()



    def _renamer(self, tre):
        """ renames newick from numbers to sample names"""
        ## get the tre with numbered tree tip labels
        names = tre.get_leaves()

        ## replace numbered names with snames
        for name in names:
            name.name = self.samples[int(name.name)]

        ## return with only topology and leaf labels
        return tre.write(format=9)



    def _finalize_stats(self):
        """ write final tree files """

        ## print stats file location:
        #print(STATSOUT.format(opr(self.files.stats)))

        ## print finished tree information ---------------------
        print(FINALTREES.format(opr(self.trees.tree)))

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
                #wctre = toytree.tree(self.trees.cons, format=0)
                #wctre.unroot()
                wctre.ladderize()
                print(wctre.get_ascii(show_internal=True, 
                                      attributes=["dist", "name"]))
                print("")
            else:
                qtre = ete3.Tree(self.trees.tree, format=0)
                #qtre = toytree.tree(self.trees.tree, format=0)
                qtre.tree.unroot()
                print(qtre.get_ascii())
                print("")

        ## print PDF filename & tips -----------------------------
        docslink = "toytree.readthedocs.io/"    
        citelink = "ipyrad.readthedocs.io/tetrad.html"
        print(LINKS.format(docslink, citelink))



    def _compute_tree_stats(self):
        """ writes support values as edge labels on unrooted tree """
        compute_tree_stats(self)



    def _save(self):
        """ save a JSON file representation of Tetrad Class for checkpoint"""

        ## save each attribute as dict
        fulldict = copy.deepcopy(self.__dict__)
        for i, j in fulldict.items():
            if isinstance(j, Params):
                fulldict[i] = j.__dict__
        fulldumps = json.dumps(fulldict,
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



    def _insert_to_array(self, start, results):
        """ inputs results from workers into hdf4 array """
        qrts, wgts, qsts = results
        #qrts, wgts = results
        #print(qrts)

        with h5py.File(self.database.output, 'r+') as out:
            chunk = self._chunksize
            out['quartets'][start:start+chunk] = qrts
            ##out['weights'][start:start+chunk] = wgts

            ## entered as 0-indexed !
            if self.checkpoint.boots:
                key = "qboots/b{}".format(self.checkpoint.boots-1)
                out[key][start:start+chunk] = qsts
            else:
                out["qstats"][start:start+chunk] = qsts



    ## thoughts on this... we need to re-estimate chunksize given whatever
    ## new parallel setup is passed in. Start from arr checkpoint to end. 
    def _load_json(self, name, workdir, quiet=False):
        """ Load a json serialized Tetrad Class object """

        ## load the JSON string and try with name+.json
        path = os.path.join(workdir, name)
        if not path.endswith(".tet.json"):
            path += ".tet.json"

        ## expand user
        path = path.replace("~", os.path.expanduser("~"))

        ## load the json file as a dictionary
        try:
            with open(path, 'r') as infile:
                fullj = _byteify(json.loads(infile.read(),
                                object_hook=_byteify), 
                            ignore_dicts=True)
        except IOError:
            raise IPyradWarningExit("""\
        Cannot find checkpoint (.tet.json) file at: {}""".format(path))

        ## set old attributes into new tetrad object
        self.name = fullj["name"]
        self.files.seqfile = fullj["files"]["seqfile"]
        self.files.mapfile = fullj["files"]["mapfile"]        
        self.dirs = fullj["dirs"]
        self.method = fullj["method"]
        self._init_seqarray(quiet=quiet)
        self._parse_names()

        ## fill in the same attributes
        for key in fullj:
            ## fill Params a little different
            if key in ["files", "database", "trees", "stats", "checkpoint"]:
                filler = fullj[key]
                for ikey in filler:
                    self.__dict__[key].__setattr__(ikey, fullj[key][ikey])
            else:
                self.__setattr__(key, fullj[key])



    ########################################################################
    ## Main functions
    ########################################################################
    def run(self, force=0, quiet=0):
        """ 
        Run quartet inference on a SNP alignment and distribute work
        across an ipyparallel cluster (ipyclient). Unless passed an 
        ipyclient explicitly, it looks for a running ipcluster instance
        running from the "default" profile, and will raise an exception
        if one is not found within a set time limit. Parameter settings
        influencing the run (e.g., nquartets, sampling method) should
        be set on the tetrad Class object itself. 

        Parameters
        ----------
        force (bool):
            ...
        quiet (bool):
            ...

        """

        ## wrap everything in a try statement so we can ensure that it will
        ## save if interrupted and we will clean up the 
        try:
            ## find an ipcluster instance
            ipyclient = ip.core.parallel.get_client(**self._ipcluster)

            ## print a message about the cluster status
            if not quiet:
                print(ip.cluster_info(ipyclient))

            ## grab 2 engines from each host (2 multi-thread jobs per host)
            rdict = ipyclient[:].apply(socket.gethostname).get_dict()
            hosts = set(rdict.values())
            hostdict = {host: [i for i in rdict if rdict[i] == host] for host in hosts}
            targets = list(itertools.chain(*[hostdict[i][:2] for i in hostdict]))
            lbview = ipyclient.load_balanced_view(targets=targets)

            ## get or init quartet sampling ---------------------------
            if not self._chunksize:
                self.nquartets = n_choose_k(len(self.samples), 4)
                if self.method != 'equal':
                    ## store N sampled quartets into the h5 array
                    self._store_N_samples(ncpus=len(lbview))
                else:
                    self._store_equal_samples(ncpus=len(lbview))

            ## calculate invariants for the full array ----------------
            start = time.time()            
            if not self.trees.tree:
                print("inferring {} induced quartet trees".format(self.nquartets))
                self.inference(start, lbview)
                print("")

            ## calculate for bootstraps -------------------------------            
            start = time.time()
            while self.checkpoint.boots < self.nboots:
                if self.files.mapfile:
                    self._sample_bootseq_array_map()
                else:
                    self._sample_bootseq_array() 

                ## start boot inference, (1-indexed !!!)
                self.checkpoint.boots += 1
                self.inference(start, lbview)

            ## write outputs with bootstraps
            print("")
            self.files.stats = os.path.join(self.dirs, self.name+"_stats.txt")
            if not self.kwargs.get("cli"):
                self._compute_tree_stats()
            else:
                self._finalize_stats()               


        ## handle exceptions so they will be raised after we clean up below
        except KeyboardInterrupt as inst:
            LOGGER.info("assembly interrupted by user.")
            print("\nKeyboard Interrupt by user. Cleaning up...")
            raise

        except IPyradWarningExit as inst:
            LOGGER.info("IPyradWarningExit: %s", inst)
            print(inst)
            raise 

        except Exception as inst:
            LOGGER.info("caught an unknown exception %s", inst)
            print("\n  Exception found: {}".format(inst))
            raise

        ## close client when done or interrupted
        finally:
            try:
                ## save the Assembly
                self._save()                
                
                ## can't close client if it was never open
                if ipyclient:

                    ## if CLI, stop jobs and shutdown
                    if 'ipyrad-cli' in self._ipcluster["cluster_id"]:
                        LOGGER.info("  shutting down engines")
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        LOGGER.info("  finished shutdown")
                    else:
                        if not ipyclient.outstanding:
                            ipyclient.purge_everything()
                        else:
                            ## nanny: kill the engines left running, report kill.
                            ipyclient.shutdown(hub=True, block=False)
                            ipyclient.close()
                            print("\n  warning: ipcluster shutdown and must be restarted")
            
            ## if exception is close and save, print and ignore
            except Exception as inst2:
                LOGGER.error("shutdown warning: %s", inst2)




    def inference(self, start, lbview):
        """ 
        Inference sends slices of jobs to the parallel engines for computing
        and collects the results into the output hdf5 array as they finish. 
        """

        ## an iterator to distribute sampled quartets in chunks
        gen = xrange(self.checkpoint.arr, self.nquartets, self._chunksize)
        njobs = sum(1 for _ in gen)
        jobiter = iter(gen)
        LOGGER.info("chunksize: %s, start: %s, total: %s, njobs: %s", \
            self._chunksize, self.checkpoint.arr, self.nquartets, njobs)

        ## if bootstrap create an output array for results unless we are 
        ## restarting an existing boot, then use the one already present
        key = "b{}".format(self.checkpoint.boots)
        with h5py.File(self.database.output, 'r+') as out:
            if key not in out["qboots"].keys():
                out["qboots"].create_dataset(key, 
                                            (self.nquartets, 4), 
                                            dtype=np.uint32, 
                                            chunks=(self._chunksize, 4))

        ## initial progress bar
        elapsed = datetime.timedelta(seconds=int(time.time()-start))
        if not self.checkpoint.boots:
            progressbar(1, 0, 
                " initial tree | {} | ".format(elapsed), spacer=0)
        else:
            progressbar(self.nboots, self.checkpoint.boots, 
                " boot {:<7} | {} | ".format(self.checkpoint.boots, elapsed), 
                spacer=0)

        ## distribute jobs across nodes
        res = {}
        for _ in xrange(njobs):
            ## get chunk of quartet samples and send to a worker engine
            qidx = jobiter.next()
            LOGGER.info('submitting chunk: %s', qidx)
            with h5py.File(self.database.input, 'r') as inh5:
                smps = inh5["samples"][qidx:qidx+self._chunksize]
                res[qidx] = lbview.apply(nworker, *[self, smps, TESTS])

        ## keep adding jobs until the jobiter is empty
        done = 0
        while 1:
            ## check for finished jobs
            curkeys = res.keys()
            finished = [i.ready() for i in res.values()]

            ## remove finished and submit new jobs
            if any(finished):
                for ikey in curkeys:
                    if res[ikey].ready():
                        if res[ikey].successful():
                            LOGGER.info("collecting results chunk: %s, tool %s ms", ikey, res[ikey].elapsed*1e3)
                            ## track finished
                            done += 1
                            ## insert results into hdf5 data base
                            results = res[ikey].get(0)
                            LOGGER.info("%s", results[1])
                            self._insert_to_array(ikey, results) #, bidx)
                            ## purge memory of the old one
                            del res[ikey]
                        else:
                            ## print error if something went wrong
                            raise IPyradWarningExit(""" error in 'inference'\n{}
                                """.format(res[ikey].exception()))

                    ## submit new jobs
                    try:
                        ## send chunk off to be worked on
                        qidx = jobiter.next()
                        with h5py.File(self.database.input, 'r') as inh5:
                            smps = inh5["samples"][qidx:qidx+self._chunksize]
                        res[qidx] = lbview.apply(nworker, *[self, smps, TESTS])

                    ## if no more jobs then just wait until these are done
                    except StopIteration:
                        continue
            else:
                time.sleep(0.01)

            ## print progress unless bootstrapping, diff progbar for that.
            elapsed = datetime.timedelta(seconds=int(time.time()-start))
            if not self.checkpoint.boots:
                progressbar(njobs, done, 
                    " initial tree | {} | ".format(elapsed), spacer=0)
            else:
                progressbar(self.nboots, self.checkpoint.boots, 
                    " boot {:<7} | {} | ".format(self.checkpoint.boots, elapsed), 
                    spacer=0)

            ## done is counted on finish, so this means we're done
            if njobs == done:
               break

        ## dump quartets to a file
        self._dump_qmc()

        ## send to qmc
        if not self.checkpoint.boots:
            self._run_qmc(0)
        else:
            self._run_qmc(1)            

        ## reset the checkpoint_arr
        self.checkpoint.arr = 0



################################################################
## STATS FUNCTIONS
################################################################
def compute_tree_stats(self):
    """ 
    compute stats for stats file and NHX tree features
    """

    ## get name indices
    names = self.samples

    ## get majority rule consensus tree of weighted Q bootstrap trees
    if self.nboots:
        fulltre = ete3.Tree(self.trees.tree, format=0)
        #fulltre = toytree.tree(self.trees.tree, format=0)
        #[toytree.tree(i.strip(), format=0).tree for i in inboots.readlines()]
        with open(self.trees.boots, 'r') as inboots:
            #wboots = [fulltre.tree] + \
            wboots = [fulltre] + \
            [ete3.Tree(i.strip(), format=0) for i in inboots.readlines()]
        wctre, wcounts = consensus_tree(wboots, names=names)
        self.trees.cons = os.path.join(self.dirs, self.name + ".cons")
        with open(self.trees.cons, 'w') as ocons:
            ocons.write(wctre.write(format=0))

    ## build stats file and write trees
    self.trees.nhx = os.path.join(self.dirs, self.name + ".nhx")
    with open(self.files.stats, 'w') as ostats:

        ## print Tetrad info
        ostats.write("## Analysis info\n")
        ostats.write("{:<30}  {:<20}\n".format("Name", self.name))
        ostats.write("{:<30}  {:<20}\n".format("Sampling_method", self.method))
        ostats.write("{:<30}  {:<20}\n".format("Sequence_file", self.files.seqfile))
        ostats.write("{:<30}  {:<20}\n".format("Map_file", self.files.mapfile))
        used_treefile = [self.files.guidetreefile if self.method == 'equal' else None][0]
        ostats.write("{:<30}  {:<20}\n".format("Guide_tree", used_treefile))
        ostats.write("\n")

        ## get Tetrad stats
        ostats.write("## quartet statistics (coming soon!!)\n")
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
        ostats.write("{:<30}  {:<20}\n".format("Initial_tree", self.trees.tree))
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
    with open(self.trees.nhx, 'w') as outtre:
        outtre.write(wctre.write(format=0, features=features))

    ## close the client view
    ipyclient.close()




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
    """ get the number of quartets as n-choose-k. This is used
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

# @numba.jit('f8(f8[:])', nopython=True)#, cache=True)
# def get_weights(scores):
#     """ 
#     gets quartet weights from ordered svd scores. Following 
#     description from Avni et al. 
#     """
#     ## lowest to highest [best, ils1, ils2]
#     scores.sort()
#     ## get weight given the svd scores
#     if scores[2]:
#         weight = (scores[2]-scores[0]) / \
#                  (np.exp(scores[2]-scores[1]) * scores[2])
#     else:
#         weight = 0
#     return weight



@numba.jit('u4[:](u4[:,:])', nopython=True)#, cache=True)
def count_snps(mat):
    """ 
    get dstats from the count array and return as a float tuple 
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


# @numba.jit('u1[:,:](u1[:,:],b1[:],u4[:])', nopython=True)
# def subsample_snps(seqchunk, nmask, maparr):
#     """ 
#     removes ncolumns from snparray prior to matrix calculation, and 
#     subsamples 'linked' snps (those from the same RAD locus) such that
#     for these four samples only 1 SNP per locus is kept. This information
#     comes from the 'map' array (map file). 
#     """
#     ## mask columns that contain Ns
#     rmask = np.ones(seqchunk.shape[1], dtype=np.bool_)
#     #LOGGER.info("rmask : %s %s", rmask.shape, rmask.sum())

#     for idx in xrange(rmask.shape[0]):
#         if nmask[idx]: 
#             rmask[idx] = False
#     #LOGGER.info("rmasked : %s %s", rmask.shape, rmask.sum())    

#     ## apply mask
#     newarr = seqchunk[:, rmask]

#     ## return smaller Nmasked array
#     return newarr



@numba.jit('b1[:](u1[:,:],b1[:],u4[:])', nopython=True)#, cache=True)
def subsample_snps_map(seqchunk, nmask, maparr):
    """ 
    removes ncolumns from snparray prior to matrix calculation, and 
    subsamples 'linked' snps (those from the same RAD locus) such that
    for these four samples only 1 SNP per locus is kept. This information
    comes from the 'map' array (map file). 
    """
    ## mask columns that contain Ns
    rmask = np.zeros(seqchunk.shape[1], dtype=np.bool_)

    ## apply mask to the mapfile
    last_loc = -1
    for idx in xrange(maparr.shape[0]):
        if maparr[idx] != last_loc:
            if not nmask[idx]:
                rmask[idx] = True
            last_loc = maparr[idx]
    
    ## apply mask
    #newarr = seqchunk[:, rmask]
    
    ## return smaller Nmasked array
    return rmask



# @numba.jit('u1[:,:](u1[:,:],b1[:],u4[:])', nopython=True)
# def subsample_snps_map(seqchunk, nmask, maparr):
#     """ 
#     removes ncolumns from snparray prior to matrix calculation, and 
#     subsamples 'linked' snps (those from the same RAD locus) such that
#     for these four samples only 1 SNP per locus is kept. This information
#     comes from the 'map' array (map file). 
#     """
#     ## mask columns that contain Ns
#     rmask = np.ones(seqchunk.shape[1], dtype=np.bool_)

#     ## apply mask to the mapfile
#     last_snp = 0
#     for idx in xrange(rmask.shape[0]):
#         if nmask[idx]:
#             ## mask if Ns
#             rmask[idx] = False
#         else:
#             ## also mask if SNP already sampled 
#             this_snp = maparr[idx]
#             if maparr[idx] == last_snp:
#                 rmask[idx] = False
#             ## record this snp
#             last_snp = this_snp  
    
#     ## apply mask
#     newarr = seqchunk[:, rmask]
    
    ## return smaller Nmasked array
#    return newarr



@numba.jit('u4[:,:,:](u1[:,:],u4[:],b1[:])', nopython=True)#, cache=True)
def chunk_to_matrices(narr, mapcol, nmask):
    """ 
    numba compiled code to get matrix fast.
    arr is a 4 x N seq matrix converted to np.int8
    I convert the numbers for ATGC into their respective index for the MAT
    matrix, and leave all others as high numbers, i.e., -==45, N==78. 
    """

    ## get seq alignment and create an empty array for filling
    mats = np.zeros((3, 16, 16), dtype=np.uint32)

    ## replace ints with small ints that index their place in the 
    ## 16x16. This no longer checks for big ints to exclude, so resolve=True
    ## is now the default, TODO. 
    last_loc = -1
    for idx in xrange(mapcol.shape[0]):
        if not nmask[idx]:
            if not mapcol[idx] == last_loc:
                i = narr[:, idx]
                mats[0, (4*i[0])+i[1], (4*i[2])+i[3]] += 1      
                last_loc = mapcol[idx]

    ## fill the alternates
    x = np.uint8(0)
    for y in np.array([0, 4, 8, 12], dtype=np.uint8):
        for z in np.array([0, 4, 8, 12], dtype=np.uint8):
            mats[1, y:y+np.uint8(4), z:z+np.uint8(4)] = mats[0, x].reshape(4, 4)
            mats[2, y:y+np.uint8(4), z:z+np.uint8(4)] = mats[0, x].reshape(4, 4).T
            x += np.uint8(1)

    return mats



@numba.jit(nopython=True)#, cache=True)
def calculate(seqnon, mapcol, nmask, tests):
    """ groups together several numba compiled funcs """

    ## create empty matrices
    #LOGGER.info("tests[0] %s", tests[0])
    #LOGGER.info('seqnon[[tests[0]]] %s', seqnon[[tests[0]]])
    mats = chunk_to_matrices(seqnon, mapcol, nmask)

    ## empty svdscores for each arrangement of seqchunk
    svds = np.zeros((3, 16), dtype=np.float64)
    qscores = np.zeros(3, dtype=np.float64)
    ranks = np.zeros(3, dtype=np.float64)

    for test in range(3):
        ## get svd scores
        svds[test] = np.linalg.svd(mats[test].astype(np.float64))[1]
        ranks[test] = np.linalg.matrix_rank(mats[test].astype(np.float64))

    ## get minrank, or 11
    minrank = int(min(11, ranks.min()))
    for test in range(3):
        qscores[test] = np.sqrt(np.sum(svds[test, minrank:]**2))

    ## sort to find the best qorder
    best = np.where(qscores == qscores.min())[0]
    #best = qscores[qscores == qscores.min()][0]
    bidx = tests[best][0]
    qsnps = count_snps(mats[best][0])

    return bidx, qsnps
    #, qscores
    #return bidx, qscores, qsnps



def nworker(data, smpchunk, tests):
    """ The workhorse function. Not numba. """
    
    ## tell engines to limit threads
    #numba.config.NUMBA_DEFAULT_NUM_THREADS = 1
    
    ## open the seqarray view, the modified array is in bootsarr
    with h5py.File(data.database.input, 'r') as io5:
        seqview = io5["bootsarr"][:]
        maparr = io5["bootsmap"][:]

    ## create an N-mask array of all seq cols
    nall_mask = seqview[:] == 78

    ## tried numba compiling everythign below here, but was not faster
    ## than making nmask w/ axis arg in numpy

    ## get the input arrays ready
    rquartets = np.zeros((smpchunk.shape[0], 4), dtype=np.uint16)
    rweights = None
    #rweights = np.ones(smpchunk.shape[0], dtype=np.float64)
    rdstats = np.zeros((smpchunk.shape[0], 4), dtype=np.uint32)

    #times = []
    ## fill arrays with results using numba funcs
    for idx in xrange(smpchunk.shape[0]):
        ## get seqchunk for 4 samples (4, ncols) 
        sidx = smpchunk[idx]
        seqchunk = seqview[sidx]

        ## get N-containing columns in 4-array, and invariant sites.
        nmask = np.any(nall_mask[sidx], axis=0)
        nmask += np.all(seqchunk == seqchunk[0], axis=0)

        ## get matrices if there are any shared SNPs
        ## returns best-tree index, qscores, and qstats
        #bidx, qscores, qstats = calculate(seqchunk, maparr[:, 0], nmask, tests)
        bidx, qstats = calculate(seqchunk, maparr[:, 0], nmask, tests)        
        
        ## get weights from the three scores sorted. 
        ## Only save to file if the quartet has information
        rdstats[idx] = qstats 
        rquartets[idx] = smpchunk[idx][bidx]

    return rquartets, rweights, rdstats 
    #return rquartets, rweights, rdstats 




########################################################################
## GLOBALS
########################################################################

## the three indexed resolutions of each quartet
TESTS = np.array([[0, 1, 2, 3], 
                  [0, 2, 1, 3], 
                  [0, 3, 1, 2]], dtype=np.uint8)


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
    inferring {} quartet trees
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


@numba.jit(nopython=True)#, cache=True)
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



@numba.jit(nopython=True)#, cache=True)
def get_spans(maparr, spans):
    """ get span distance for each locus in original seqarray """
    ## start at 0, finds change at 1-index of map file
    bidx = 1
    spans = np.zeros((maparr[-1, 0], 2), np.uint64)
    ## read through marr and record when locus id changes
    for idx in xrange(1, maparr.shape[0]):
        cur = maparr[idx, 0]
        if cur != bidx:
            idy = idx + 1
            spans[cur-2, 1] = idx
            spans[cur-1, 0] = idx
            bidx = cur
    spans[-1, 1] = maparr[-1, -1]
    return spans



@numba.jit(nopython=True)#, cache=True)
def get_shape(spans, loci):
    """ get shape of new bootstrap resampled locus array """
    width = 0
    for idx in xrange(loci.shape[0]):
        width += spans[loci[idx], 1] - spans[loci[idx], 0]
    return width
    


@numba.jit(nopython=True)#, cache=True)
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
                invert_test = np.bool_(clades[aidx]) != np.bool_(clades[idx])

                if np.all(invert_test):
                    counts[aidx] += counts[idx]
                    conflict = True
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
            # if not node.is_root():
            #     if node.up.is_root():
            #         if bits[0]:
            #             bits.invert()
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
        #node = toytree.tree(name=name).tree
        for child in children:
            node.add_child(child)
        if not node.is_leaf():
            node.dist = int(round(100*countdict[clade]))
            node.support = int(round(100*countdict[clade]))
        else:
            node.dist = int(100) 
            node.support = int(100)
        
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
    with h5py.File(data.database.output, 'r') as io5:
        qrts = io5["quartets"][idx:idx+data._chunksize]
        for qrt in qrts:
            sqrt = set(qrt)
            if len(sqrt.intersection(idxd)) > 1:
                if len(sqrt.intersection(idxu)) > 1:
                    sampled += 1
            idx += data._chunksize

    return sampled



## GLOBALS #############################################################

STATSOUT = """
Statistics for sampling, discordance, and tree support:
  {}
"""

FINALTREES = """\
Best tree inferred from the full SNP array:
  {}
"""

BOOTTREES = """\
Extended majority-rule consensus with support values:
  {}

All bootstrap trees:
  {}
"""

ASCII_TREE = """\
  {}
"""

LINKS = """\
* For tips on plotting trees: {}     
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
    #import ipyrad.analysis as ipa
    #import ipyrad as ip
    #import ipyparallel as ipp
    pass

    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/cli/cli.json")
    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/iptutorial/cli.json")
    ## run


