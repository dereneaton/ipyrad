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
# pylint: disable=R0914



from __future__ import print_function, division
import os
import sys
import json
import h5py
import time
import numba
import random
import ipyrad
import cStringIO
import datetime
import itertools
import subprocess
import numpy as np
import pandas as pd
import ipyparallel as ipp
from fractions import Fraction

from ipyrad.assemble.util import ObjDict, IPyradWarningExit, progressbar
from collections import OrderedDict

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

## The 16 x 16 matrix of site counts. This is just for looking at. 
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



#############################################################################
#############################################################################
## Quartet inference Class Object
#############################################################################
#############################################################################


class Quartet(object):
    """
    The main svd4tet object for storing data and checkpointing. It is 
    initialized with a name, and args to command line (e.g., sampling method, 
    starting tree, nboots, etc.). 
    """

    def __init__(self, name, outdir=os.path.curdir, method='all'):
        ## version is ipyrad version
        self._version = ipyrad.__version__

        ## name this assembly
        self.name = name
        self.dirs = os.path.realpath(outdir)
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)

        ## store cluster information the same as ipyrad Assembly objects
        self._ipcluster = {}
        for ipkey in ["id", "profile", "engines", "cores"]:
            self._ipcluster[ipkey] = None

        ## store cpu information, a more convenient call to ipyrad funcs
        self.cpus = 4

        ## Sampling method attributes
        self.method = method
        self.nboots = 0
        self.nquartets = 0
        self.chunksize = 0
        self.resolve_ambigs = 0

        ## store samples from the seqarray
        self.samples = []

        ## self.populations ## if we allow grouping samples (not yet)

        ## hdf5 data bases init and delete existing
        self.h5in = os.path.join(self.dirs, self.name+".input.h5")
        self.h5out = os.path.join(self.dirs, self.name+".output.h5")

        ## input files
        self.files = ObjDict()
        self.files.seqfile = None
        self.files.mapfile = None
        self.files.treefile = None
        self.files.qdump = None

        ## store tree data
        self.trees = ObjDict()
        ## the full trees
        self.trees.ttre = os.path.join(self.dirs, self.name+".tre")
        self.trees.wtre = os.path.join(self.dirs, self.name+".w.tre")
        ## list of bootstrap trees 
        self.trees.tboots = os.path.join(self.dirs, self.name+".boots")        
        self.trees.wboots = os.path.join(self.dirs, self.name+".w.boots")        
        ## trees with boot support as edge length
        self.trees.tbtre = os.path.join(self.dirs, self.name+".support.tre")        
        self.trees.wbtre = os.path.join(self.dirs, self.name+".w.support.tre")     
        ## NHX formatted tre with rich information
        self.trees.tnhx = os.path.join(self.dirs, self.name+".nhx.tre")     
        self.trees.wnhx = os.path.join(self.dirs, self.name+".w.nhx.tre")     
        ## PDF
        self.trees.pdf = os.path.join(self.dirs, self.name+".pdf")

        ## checkpointing information
        self.checkpoint = ObjDict()
        self.checkpoint.boots = 0
        self.checkpoint.arr = 0



    ## convenience function
    #def opj(self, path):
    #    return os.path.join(self.dirs, path)


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
        ## make sure sorted
        self.samples = sorted(self.samples)



    ## parallel launcher and wrapper to kill remote engines at exit
    def _launch(self, nwait):
        """ 
        Creates a client for a given profile to connect to the running 
        clusters. 
        """
        #save_stdout = sys.stdout           
        try: 
            clusterargs = [self._ipcluster['id'], self._ipcluster["profile"]]
            argnames = ["cluster_id", "profile"]
            args = {key:value for key, value in zip(argnames, clusterargs)}

            ## wait for at least 1 engine to connect
            for _ in range(nwait):
                try:
                    ## using this wrap to avoid ipyparallel's ugly warnings
                    ## save orig stdout
                    save_stdout = sys.stdout 
                    save_stderr = sys.stderr
                    ## file-like obj to catch stdout
                    sys.stdout = cStringIO.StringIO()
                    sys.stderr = cStringIO.StringIO()                    
                    ## run func with stdout hidden
                    ipyclient = ipp.Client(**args)
                    ## resets stdout
                    sys.stdout = save_stdout
                    sys.stderr = save_stderr
                    break

                except IOError as inst:
                    time.sleep(0.1)

            ## check that at least one engine has connected
            for _ in range(300):
                initid = len(ipyclient)
                time.sleep(0.1)
                if initid:
                    break
        except KeyboardInterrupt as inst:
            ## ensure stdout is reset even if Exception was raised            
            sys.stdout = save_stdout
            raise inst

        except IOError as inst:
            ## ensure stdout is reset even if Exception was raised
            sys.stdout = save_stdout
            print(inst)
            raise inst
        return ipyclient



    def _clientwrapper(self, func, args, nwait):
        """ wraps a call with error messages for when ipyparallel fails"""
        ## emtpy error string
        inst = ""

        ## wrapper to ensure closure of ipyparallel
        try:
            ipyclient = ""
            ipyclient = self._launch(nwait)
            args.append(ipyclient)
            func(*args)

        except Exception as inst:
            ## Caught unhandled exception, print and reraise
            LOGGER.error(inst)
            print("\n  Caught unknown exception - {}".format(inst))
            raise  ## uncomment raise to get traceback

        ## close client when done or interrupted
        finally:
            try:
                self.save()                
                ## can't close client if it was never open
                if ipyclient:
                    ## if CLI, stop jobs and shutdown
                    ipyclient.abort()                        
                    ipyclient.close()
            ## if exception is close and save, print and ignore
            except Exception as inst2:
                LOGGER.error("shutdown warning: %s", inst2)
            if inst:
                IPyradWarningExit(inst)



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
        spath = open(self.files.seqfile, 'r')
        line = spath.readline().strip().split()
        ntax = int(line[0])
        nbp = int(line[1])

        ## make a tmp seq array
        print("  loading array [{} taxa x {} bp]".format(ntax, nbp))        
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
            if self.resolve_ambigs:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## save modified array to disk            
            io5["bootsarr"][:] = tmpseq

            ## memory cleanup
            del tmpseq

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
            if self.resolve_ambigs:
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
            if self.resolve_ambigs:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## store data sets
            io5.create_dataset("bootsmap", data=tmpmap)
            io5.create_dataset("bootsarr", data=tmpseq)

            LOGGER.info("resampled bootsarr \n %s", io5["bootsarr"][:, :20])
            LOGGER.info("resampled bootsmap \n %s", io5["bootsmap"][:20, :])




    ## Functions to fill h5in with samples
    def store_N_samples(self):
        """ Find all quartets of samples and store in a large array """
        ## print header
        # if self.method == 'all':
        #     print("    loading all {} possible quartets".format(self.nquartets))        
        # else:
        #     print("    loading {} random quartet samples".format(self.nquartets))        

        ## create a chunk size for sampling from the array of quartets. This should
        ## be relatively large so that we don't spend a lot of time doing I/O, but
        ## small enough that jobs finish every few hours since that is how the 
        ## checkpointing works.
        breaks = 2
        if self.nquartets < 5000:
            breaks = 1
        if self.nquartets > 100000:
            breaks = 4
        if self.nquartets > 500000:
            breaks = 8

        self.chunksize = (self.nquartets // (breaks * self.cpus)) + \
                         (self.nquartets % (breaks * self.cpus))
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
                                dtype=np.float32, chunks=(self.chunksize, ))
            io5.create_dataset("dstats", (self.nquartets, 3), 
                                dtype=np.float32, chunks=(self.chunksize, 3))

        ## append to h5 IN array (which also has seqarray) and fill it
        with h5py.File(self.h5in, 'a') as io5:
            ## create data sets
            io5.create_dataset("samples", (self.nquartets, 4), 
                               dtype=np.uint16, chunks=(self.chunksize, 4),
                               compression='gzip')

            ## populate array with all possible quartets. This allows us to 
            ## sample from the total, and also to continue from a checkpoint
            qiter = itertools.combinations(xrange(len(self.samples)), 4)
            i = 0

            ## fill 1000 at a time for efficiency
            while i < self.nquartets:
                if self.method != "all":
                    qiter = []
                    while len(qiter) < min(1000, io5["samples"].shape[0]):
                        qiter.append(
                            random_combination(range(len(self.samples)), 4))
                    dat = np.array(qiter)
                else:
                    dat = np.array(list(itertools.islice(qiter, 1000)))
                io5["samples"][i:i+dat.shape[0]] = dat
                i += 1000



    def store_equal_samples(self):
        """ 
        sample quartets even across splits of the starting tree 
        """
        
        ## choose chunker for h5 arr
        self.chunksize = (self.nquartets // (2 * len(self.cpus))) +\
                         (self.nquartets % (2 * len(self.cpus)))
        LOGGER.info("E: nquarts = %s, chunk = %s", self.nquartets, self.chunksize)

        ## get starting tree
        tre = ete3.Tree(".tmptre")
        tre.unroot()
        print("  starting tree: \n {}".format(tre))

        ## randomly sample all splits of tree
        splits = [([z.name for z in i], 
                   [z.name for z in j]) \
                   for (i, j) in tre.get_edges()]
    
        ## only keep internal splits
        splits = [i for i in splits if all([len(j) > 1 for j in i])]
        N = len(self.samples)
        if len(splits) < ((N * (N-1)) // 2):
            print("  starting tree is unresolved, sample more quartets")

        ## turn each into an iterable split sampler
        ## if the nquartets for that split is small, then sample all of them
        ## if it is UUUUUGE, then make it a random sampler from that split
        qiters = []
        ## how many quartets are we gonna sample from each quartet?
        squarts = self.nquartets // len(splits)

        for split in splits:
            ## if small number at this split then sample all
            if n_choose_k(len(split[0]), 2) * n_choose_k(len(split[1]), 2) < squarts*2:
                qiter = itertools.product(
                            itertools.combinations(split[0], 2), 
                            itertools.combinations(split[1], 2))
            ## else create random sampler for this split
            else:
                qiter = (random_product(split[0], split[1]) for _ in xrange(self.nquartets))
            qiters.append(qiter)

        ## make qiters infinitely cycling
        qiters = itertools.cycle(qiters)

        ## iterate over qiters sampling from each, if one runs out, keep 
        ## sampling from remaining qiters. Keep going until samples is filled
        with h5py.File(self.h5in, 'a') as io5:
            ## create data sets
            del io5["samples"]
            io5.create_dataset("samples", (self.nquartets, 4), 
                                      dtype=np.uint16, 
                                      chunks=(self.chunksize, 4),
                                      compression='gzip')

            ## fill 1000 at a time for efficiency
            i = 0
            while i < self.nquartets:
                qdat = []
                while len(qdat) < min(1000, io5["samples"].shape[0]): 
                    qiter = qiters.next()
                    try:
                        qdat.append(qiter.next())
                    except StopIteration:
                        print(len(qdat))
                        continue
                dat = np.array(qdat, dtype=np.uint16)
                io5["samples"][i:i+dat.shape[0]] = dat
                i += 1000



    ## functions to get samples for the 'equal' method
    def get_equal_samples(self, ipyclient):
        """ get starting tree for 'equal' sampling method """

        ## infer starting tree if one wasn't provided
        if not self.files.treefile:
            ## grab the minimum needed for a good tree
            self.nquartets = len(self.samples)**2.8
            self.store_N_samples()
            print("""\
    loading {} random quartet samples to infer a starting tree 
    inferring {} x 3 quartet trees
    """.format(self.nquartets**2.8, 3*(self.nquartets**2.8)))

            ## run inference functions on sampled quartets 
            self._clientwrapper(self.inference, [0], 45)

        ## sample quartets from starting tree
        print("""\
    loading {} equal-splits quartets from starting tree
    """.format(self.nquartets))
        self.store_equal_samples()
    
        ## remove starting tree tmp files
        tmps = [self.tre, self.wtre, self.tboots, 
                self.wboots, self.tbtre, self.wbtre]
        for tmp in tmps:
            if os.path.exists(tmp):
                os.remove(tmp)



    def run_qmc(self, boot):
        """ runs quartet max-cut on a quartets file """

        ## make call lists
        cmd1 = " ".join(
                [ipyrad.bins.qmc,
                " qrtt="+self.files.qdump,
                " weights=off"+
                " otre=.tmptre"])

        cmd2 = " ".join(
                [ipyrad.bins.qmc,
                " qrtt="+self.files.qdump,
                " weights=on"+
                " otre=.tmpwtre"])

        ## run them
        for cmd in [cmd1, cmd2]:
            try:
                subprocess.check_call(cmd, shell=True,
                                           stderr=subprocess.STDOUT,
                                           stdout=subprocess.PIPE)
            except subprocess.CalledProcessError as inst:
                LOGGER.error("Error in wQMC: \n({}).".format(inst))
                LOGGER.error(subprocess.STDOUT)
                raise inst

        ## read in the tmp files since qmc does not pipe
        intmptre = open(".tmptre", 'r')
        intmpwtre = open(".tmpwtre", 'r')

        ## convert int names back to str names
        tmptre = self.renamer(ete3.Tree(intmptre.read().strip()))
        tmpwtre = self.renamer(ete3.Tree(intmpwtre.read().strip()))

        ## save the boot tree
        if boot:
            with open(self.trees.tboots, 'a') as outboot:
                outboot.write(tmptre+"\n")
            with open(self.trees.wboots, 'a') as outboot:
                outboot.write(tmpwtre+"\n")

        ## store full data trees to Assembly
        else:
            with open(self.trees.ttre, 'w') as outtree:
                outtree.write(tmptre)
            with open(self.trees.wtre, 'w') as outtree:
                outtree.write(tmpwtre)

        ## save JSON file checkpoint
        intmptre.close()
        intmpwtre.close()
        self.save()



    def dump_qmc(self):
        """ prints the quartets to a file formatted for wQMC """
        ## open the h5 database
        io5 = h5py.File(self.h5out, 'r')

        ## create an output file for writing
        self.files.qdump = os.path.join(self.dirs, self.name+".quartets.txt")
        outfile = open(self.files.qdump, 'w')
        LOGGER.info("qdump file %s", self.files.qdump)

        ## todo: should pull quarts order in randomly
        for idx in xrange(0, self.nquartets, self.chunksize):
            quarts = [list(j) for j in io5["quartets"][idx:idx+self.chunksize]]
            weight = io5["weights"][idx:idx+self.chunksize]
            chunk = ["{},{}|{},{}:{}".format(*i+[j]) for i, j \
                                                    in zip(quarts, weight)]
            outfile.write("\n".join(chunk)+"\n")

        ## close output file and h5 database
        outfile.close()
        io5.close()



    def renamer(self, tre):
        """ renames newick from numbers to sample names"""
        ## order the numbered tree tip lables
        names = tre.get_leaves()
        names.sort(key=lambda x: int(x.name))

        ## order the sample labels in the same order they are 
        ## in the seq file (e.g., .snp, .phy)
        snames = self.samples
        snames.sort()

        ## replace numbered names with snames
        for (tip, sname) in zip(names, snames):
            tip.name = sname

        ## return with only topology and leaf labels
        return tre.write(format=9)



    def write_output_splash(self, with_boots=1):
        """ write final tree files """
        ## create tree with support values on edges
        self.write_supports(with_boots)

        ## hide error message during tree plotting
        #save_stdout = sys.stdout           
        ## file-like obj to catch stdout
        #sys.stdout = cStringIO.StringIO()  
        ## run func with stdout hidden
        #ipyclient = ipp.Client(**args)
        ## resets stdout
        #sys.stdout = save_stdout    
        #self.quickfig()

        ## print finished tree information ---------------------
        print("""
  Final quartet-joined and weighted quartet-joined (.w.) tree files:
    - {}
    - {}
    """.format(opr(self.trees.ttre), opr(self.trees.wtre)))

        ## print bootstrap information --------------------------
        if with_boots:
            print("""\
  Bootstrap trees:
    - {}
    - {}

  Final tree with bootstrap support as edge lengths:
    - {}
    - {}
    """.format(opr(self.trees.tboots), opr(self.trees.wboots), 
               opr(self.trees.tbtre), opr(self.trees.wbtre)))

        ## print rich information--------------------------------
        print("""\
  Final tree with rich information in NHX format:
    - {}
    - {}
    """.format(opr(self.trees.tnhx), opr(self.trees.wnhx)))

        ## print ASCII tree --------------------------------------
        qtre = ete3.Tree(self.trees.tnhx, format=0)
        qtre.unroot()
        print("""\
  ASCII view of unrooted topology from unweighted analysis
    {}
    """.format(qtre.get_ascii(attributes=["name", "support"]), show_internal=True))

        ## print PDF filename & tips -----------------------------
        docslink = "ipyrad.readthedocs.org/cookbook.html"    
        citelink = "ipyrad.readthedocs.org/svd4tet.html"
        print("""\
  * For tips on plotting these trees in R see: 
    - {}     
  * For tips on citing this software see: 
    - {} 

    """.format(docslink, citelink))




    def write_supports(self, with_boots):
        """ writes support values as edge labels on unrooted tree """
        ## get name indices
        names = self.samples
        names.sort()

        ## get unrooted best trees
        otre = ete3.Tree(self.trees.ttre, format=0)
        otre.unroot()
        for node in otre.traverse():
            node.add_feature("bootstrap", 0)
            node.add_feature("quartets_total", get_total(otre, node))
            node.add_feature("quartets_sampled", get_sampled(self, otre, node, names))
            try:
                prop = 100*(float(node.quartets_sampled) / node.quartets_total)
            except ZeroDivisionError:
                prop = 0.0
            node.add_feature("quartets_sampled_prop", prop)
            node.dist = 0
            node.support = 0

        ## get unrooted weighted quartets tree
        wtre = ete3.Tree(self.trees.wtre, format=0)
        wtre.unroot()
        for node in wtre.traverse():
            node.add_feature("bootstrap", 0)
            node.add_feature("quartets_total", get_total(wtre, node))
            node.add_feature("quartets_sampled", get_sampled(self, wtre, node, names))
            try:
                prop = 100*(float(node.quartets_sampled) / node.quartets_total)
            except ZeroDivisionError:
                prop = 0.0
            node.add_feature("quartets_sampled_prop", prop)
            node.dist = 0
            node.support = 0

        ## get unrooted boot trees
        if with_boots:
            oboots = open(self.trees.tboots, 'r').readlines()
            wboots = open(self.trees.wboots, 'r').readlines()
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
            with open(self.trees.tbtre, 'w') as outtre:
                outtre.write(otre.write(format=5))

            with open(self.trees.wbtre, 'w') as outtre:
                outtre.write(wtre.write(format=5))
            features = ["bootstrap", "quartets_total", "quartets_sampled", "quartets_sampled_prop"]            
        else:
            features = ["quartets_total", "quartets_sampled", "quartets_sampled_prop"]

        ## return as NHX format with extra info
        with open(self.trees.tnhx, 'w') as outtre:
            outtre.write(wtre.write(format=0, features=features))

        with open(self.trees.wnhx, 'w') as outtre:
            outtre.write(wtre.write(format=0, features=features))


    ## not currently being used...
    def quickfig(self):
        """ make a quick ete3 fig. Plots total quartets """
        ts = ete3.TreeStyle()
        ts.layout_fn = layout
        ts.show_leaf_name = False
        ts.mode = 'r'
        ts.draw_guiding_lines = True
        ts.show_scale = False
        ts.scale = 25

        tre = ete3.Tree(self.trees.nhx)
        tre.ladderize()
        tre.convert_to_ultrametric(tree_length=len(tre)//2)
        tre.render(file_name=self.trees.pdf, h=40*len(tre), tree_style=ts)



    def save(self):
        """ save a JSON file representation of Quartet Class for checkpoint"""

        ## save each attribute as dict
        fulldumps = json.dumps(self.__dict__, 
                               sort_keys=False, 
                               indent=4, 
                               separators=(",", ":"),
                               )

        ## save to file, make dir if it wasn't made earlier
        assemblypath = os.path.join(self.dirs, self.name+".svd.json")
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



    def insert_to_array(self, start, results):
        """ inputs results from workers into hdf4 array """
        qrts, wgts, dsts = results

        with h5py.File(self.h5out, 'r+') as out:
            chunk = self.chunksize
            out['quartets'][start:start+chunk] = qrts
            out['weights'][start:start+chunk] = wgts
            out["dstats"][start:start+chunk] = dsts
        ## save checkpoint
        #data.svd.checkpoint_arr = np.where(ww == 0)[0].min()


    ########################################################################
    ## Main functions
    ########################################################################
    def run(self, force):
        """ 
        Run svd4tet inference on a sequence or SNP alignment for all samples in
        the Assembly. By default the job starts from 0 or where it last left 
        off, unless force=True, then it starts from 0. This passes args to 
        the function 'inference' which does the real work. 
        """
        ## Check for existing (loaded) Quartet Class Object
        fresh = 0
        if not force:
            if self.checkpoint.boots or self.checkpoint.arr:
                print(LOADING_MESSAGE.format(
                    self.method, self.checkpoint.arr, self.checkpoint.boots))
            else:
                ## if no checkpoints then start fresh
                fresh = 1

        ## if svd results do not exist or force then restart
        if force or fresh:
            ## if nquartets not provided then get all possible
            if not self.nquartets:
                self.nquartets = n_choose_k(len(self.samples), 4)
            ## store N sampled quartets into the h5 array
            self.store_N_samples()

            ## infer a starting tree from the N sampled quartets 
            ## this could take a long time. Calls the parallel client.
            if self.method == "equal":
                self.store_equal_samples()
            
        ## run the full inference or print finished prog bar if it's done
        if not self.checkpoint.boots:
            print("  inferring {} x 3 quartet trees".format(self.nquartets))
            self._clientwrapper(self.inference, [0], 45)

        else:
            print("  full inference finished")
            elapsed = datetime.timedelta(seconds=int(0)) 
            progressbar(20, 20, " | {}".format(elapsed))


        ## run the bootstrap replicates -------------------------------
        if self.nboots:
            print("  running {} bootstrap replicates".format(self.nboots))  
    
            ## load from current boot
            for bidx in xrange(self.checkpoint.boots, self.nboots):
                ## get resampled array and set checkpoint
                if self.checkpoint.arr == 0:
                    if self.files.mapfile:
                        self.sample_bootseq_array_map()
                    else:
                        self.sample_bootseq_array() 
                    self.checkpoint.boots = bidx

                ## start boot inference, (1-indexed !!!)
                self._clientwrapper(self.inference, [bidx+1], 45)            

            ## ensure that all qmc jobs are finished
            #ipyclient = ipp.Client()
            #ipyclient.wait()

            ## write outputs with bootstraps
            self.write_output_splash(with_boots=1)

        else:
            ## write outputs without bootstraps
            self.write_output_splash(with_boots=0)



    def inference(self, bidx, ipyclient):
        """ run inference and store results """
        ## an iterator to distribute sampled quartets in chunks
        jobiter = iter(xrange(self.checkpoint.arr, self.nquartets, self.chunksize))
        njobs = sum(1 for _ in jobiter)
        jobiter = iter(xrange(self.checkpoint.arr, self.nquartets, self.chunksize))
        LOGGER.info("chunksize: %s, start: %s, total: %s, njobs: %s", \
                self.chunksize, self.checkpoint.arr, self.nquartets, njobs)

        ## a distributor for engine jobs
        lbview = ipyclient.load_balanced_view()

        ## the three indexed resolutions of each quartet
        tests = np.array([[0, 1, 2, 3], 
                          [0, 2, 1, 3], 
                          [0, 3, 1, 2]], dtype=np.uint8)

        ## start progress bar timer and submit initial n jobs
        start = time.time()
        res = {}
        for _ in xrange(min(self.cpus, njobs)):
            ## get chunk of quartet samples and send to a worker engine
            qidx = jobiter.next()
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
                            self.insert_to_array(ikey, results)
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
                        with h5py.File(self.h5in, 'r') as inh5:
                            smps = inh5["samples"][qidx:qidx+self.chunksize]
                        ## send chunk off to be worked on
                        res[qidx] = lbview.apply(nworker, *[self, smps, tests])
                        #print(qidx, smps.shape, seqs.shape)

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
    """ random sampler for equa_splits func"""
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

@numba.jit('f4(f4[:])', nopython=True)
def get_weights(scores):
    """ 
    calculates quartet weights from ordered svd scores. Following 
    description from Avni et al. 
    """
    ## lowest to highest [best, ils1, ils2]
    scores.sort()
    ## calculate weight given the svd scores
    if scores[2]:
        weight = np.float32((scores[2]-scores[0]) / 
                            (np.exp(scores[2]-scores[1]) * scores[2]))
    else:
        weight = np.float32(0.00001)
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

        

@numba.jit('u1[:,:](u1[:,:],b1[:],b1[:],u4[:])', nopython=True)
def subsample_snps(seqchunk, rmask, nmask, maparr):
    """ 
    removes ncolumns from snparray prior to matrix calculation, and 
    subsamples 'linked' snps (those from the same RAD locus) such that
    for these four samples only 1 SNP per locus is kept. This information
    comes from the 'map' array (map file). 
    """
    ## mask columns that contain Ns
    for idx in xrange(rmask.shape[0]):
        if nmask[idx]: 
            rmask[idx] = False
    
    ## apply mask
    newarr = seqchunk[:, rmask]
    
    ## return smaller Nmasked array
    return newarr



@numba.jit('u1[:,:](u1[:,:],b1[:],b1[:],u4[:])', nopython=True)
def subsample_snps_map(seqchunk, rmask, nmask, maparr):
    """ 
    removes ncolumns from snparray prior to matrix calculation, and 
    subsamples 'linked' snps (those from the same RAD locus) such that
    for these four samples only 1 SNP per locus is kept. This information
    comes from the 'map' array (map file). 
    """
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
        if np.sum(i) < 16:
            mat[i[0]*4:(i[0]+4)*4]\
               [i[1]]\
               [i[2]*4:(i[2]+4)*4]\
               [i[3]] += 1
    return mat



## TODO: store stats for each quartet (nSNPs, etc.)
#@numba.jit(nopython=True)
def nworker(data, smpchunk, tests):
    """ The workhorse function. All numba. """

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

    ## get the input arrays ready
    rquartets = np.zeros((smpchunk.shape[0], 4), dtype=np.uint16)
    rweights = np.zeros(smpchunk.shape[0], dtype=np.float32)
    rmask = np.ones(seqview.shape[1], dtype=np.bool_)
    rdstats = np.zeros((smpchunk.shape[0], 3), dtype=np.float32)

    ## fill arrays with results using numba funcs
    for idx in xrange(smpchunk.shape[0]):
        ## get seqchunk for 4 samples (4, ncols) 
        sidx = smpchunk[idx]
        seqchunk = seqview[sidx]

        ## get N-containing columns in 4-array
        nmask = nall_mask[sidx].sum(axis=0, dtype=np.bool_)        

        ## remove Ncols from seqchunk & sub-sample unlinked SNPs
        seqnon = subsample(seqchunk, rmask, nmask, maparr[:, 0])
        #LOGGER.info("before sub: %s, after %s", seqchunk.shape, seqnon.shape)

        ## get svdscores for each arrangement of seqchunk
        qscores = np.zeros(3, dtype=np.float32)
        mats = np.zeros((3, 16, 16), dtype=np.uint32)
        for test in range(3):
            #LOGGER.error("%s", seqnon[tests[test]].shape)
            mats[test] = chunk_to_matrix(seqnon[tests[test]])

            ## get svd scores
            tmpscore = np.linalg.svd(mats[test].astype(np.float32))[1]
            qscores[test] = np.sqrt(tmpscore[11:]).sum()

        ## sort to find the best qorder
        best = np.where(qscores == qscores.min())[0]
        #LOGGER.error("\n%s", mats[best][0])
        bidx = tests[best][0]
        #LOGGER.info("""
        #    best: %s, 
        #    bidx: %s, 
        #    qscores: %s, 
        #    mats %s
        #    """, best, bidx, qscores, mats)

        ## get weights from the three scores sorted. 
        ## Only save to file if the quartet has information
        iweight = get_weights(qscores)
        rweights[idx] = iweight
        rquartets[idx] = smpchunk[idx][bidx]            
        #    LOGGER.error("zero weight: %s :\n %s", qscores, mats[best][0])

        ## get dstat from the best (correct) matrix 
        ## (or should we get all three?) [option]
        rdstats[idx] = abba_baba(mats[best][0])

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
    Loading from saved checkpoint:
      sampling method: {}
      array checkpoint: {}
      bootstrap checkpoint: {}

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



def load_json(path):
    """ Load a json serialized Quartet Class object """

    ## load the JSON string and try with name+.json
    fullj = ""
    checkfor = [path+".svd.json", path]
    for inpath in checkfor:
        inpath = inpath.replace("~", os.path.expanduser("~"))
        try:
            ## load in the JSON object
            with open(inpath, 'r') as infile:
                fullj = json.loads(infile.read())
        except IOError:
            pass

    ## raise an error if json not loaded
    if not fullj:
        raise IPyradWarningExit("""
    Cannot find checkpoint (JSON) file at: {}
    """.format(inpath))

    ## create a new Quartet Class
    newobj = Quartet(fullj["name"], fullj["dirs"], fullj["method"])

    ## fill in the same attributes
    for key in fullj:
        newobj.__setattr__(key, fullj[key])

    return newobj




########################################################################
## Plotting functions
########################################################################



PIECOLORS = ['#a6cee3', '#1f78b4']
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
    with h5py.File(data.h5out, 'r') as io5:
        qrts = io5["quartets"][:]
        ## iterate over quartets 
        for qrt in qrts:
            sqrt = set(qrt)
            if len(sqrt.intersection(idxd)) > 1:
                if len(sqrt.intersection(idxu)) > 1:
                    sampled += 1
    return sampled





if __name__ == "__main__":

    ## imports
    import ipyrad.analysis as ipa
    #import ipyrad as ip
    #import ipyparallel as ipp

    #DATA = ipyrad.load_json("~/Documents/ipyrad/tests/cli/cli.json")
    DATA = ipyrad.load_json("~/Documents/ipyrad/tests/iptutorial/cli.json")
    ## run
    ipa.svd4tet.wrapper(DATA, nboots=10, method='equal', nquarts=50, force=True)

