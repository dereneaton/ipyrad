#!/usr/bin/env python

""" 
SVD-quartet like tree inference. Modelled on the following papers:

Chifman, J. and L. Kubatko. 2014. Quartet inference from SNP data under 
the coalescent, Bioinformatics, 30(23): 3317-3324.

Chifman, J. and L. Kubatko. 2015. Identifiability of the unrooted species 
tree topology under the coalescent model with time-reversible substitution 
processes, site-specific rate variation, and invariable sites, Journal of 
Theoretical Biology 374: 35-47
"""

# py2/3 compat
from __future__ import print_function, division
from builtins import range

# standard lib
import os
import sys
import json
import h5py
import time
import copy
import itertools
import subprocess as sps
from fractions import Fraction
from collections import defaultdict

# third party
import numba
import ctypes
import datetime
import numpy as np
import ipyrad as ip
from ipyrad.analysis.utils import progressbar, Params
from ipyrad.assemble.util import IPyradError


# # TOYTREE is required to run TETRAD
# try:
#     from toytree import ete3mini as ete3
# except ImportError:
#     raise IPyradError("""
#     Error: tetrad requires the dependency 'toytree', which we haven't yet
#     included in the ipyrad installation. For now, you can install toytree
#     using conda with the following command: 

#     conda install toytree -c eaton-lab
#     """)


class Tetrad(object):
    """
    The main tetrad object for storing data and checkpointing. It is 
    initialized with a name, and args to command line (e.g., sampling method, 
    starting tree, nboots, etc.). 
    """

    def __init__(self,
        name, 
        data=None,
        workdir="analysis-tetrad",
        mapfile=None, 
        guidetree=None, 
        method='all', 
        nquartets=0, 
        nboots=0, 
        resolve_ambigs=True, 
        load=False,
        quiet=False,
        save_invariants=False,
        *args, 
        **kwargs):

        ## check additional arguments from kwargs.
        self.kwargs = {
            "initarr": True,
            "cli": False,
            }
        self.kwargs.update(kwargs)

        ## name this assembly
        self.samples = []
        self.name = name
        self.dirs = os.path.abspath(os.path.expanduser(workdir))
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)

        ## store ipcluster information 
        self._ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }

        ## Sampling method attributes 
        self.params = Params()
        self.params.method = method
        self.params.nboots = nboots
        self.params.nquartets = nquartets
        self.params.resolve_ambigs = resolve_ambigs
        self.params.save_invariants = save_invariants

        ## private attributes
        self._chunksize = 0
        self._tmp = None

        ## self.populations ## if we allow grouping samples
        ## (haven't done this yet)

        ## hdf5 data bases init and delete existing
        self.database = Params()
        self.database.input = os.path.join(self.dirs, self.name + ".input.h5")
        self.database.output = os.path.join(self.dirs, self.name + ".output.h5")        

        ## input files
        self.files = Params()
        self.files.data = data
        self.files.mapfile = mapfile
        self.files.tree = guidetree
        self.files.qdump = None
        self.files.stats = None 

        ## fill in files
        if mapfile:
            self.files.mapfile = os.path.abspath(os.path.expanduser(mapfile))
        if guidetree:
            self.files.tree = os.path.abspath(os.path.expanduser(guidetree))
        ## check file paths:
        if self.files.mapfile:
            if not os.path.exists(self.files.mapfile):
                raise IOError("file path {} not found".format(self.files.mapfile))
        if self.files.tree:
            if not os.path.exists(self.files.tree):
                raise IOError("file path {} not found".format(self.files.tree))

        ## load tree file paths if they exist, or None if empty
        self.trees = Params()
        self.trees.tree = os.path.join(self.dirs, self.name + ".full.tre")
        self.trees.cons = os.path.join(self.dirs, self.name + ".consensus.tre")
        self.trees.boots = os.path.join(self.dirs, self.name + ".boots")        
        self.trees.nhx = os.path.join(self.dirs, self.name + ".nhx.tre")
        for key, val in self.trees.__dict__.items():
            if not os.path.exists(val):
                self.trees.__dict__[key] = None

        ## stats is written to os.path.join(self.dirs, self.name+".stats.txt")
        #self.stats = Params()
        #self.stats.n_quartets_sampled = self.params.nquartets

        ## checkpointing information
        self.checkpoint = Params()
        self.checkpoint.boots = 0
        self.checkpoint.arr = 0

        ## init the by loading existing or parsing new data --------------------
        if load:
            self._load(self.name, self.dirs)
            self._parse_names()
        elif data and os.path.exists(data):
            if self.kwargs["initarr"]:
                self._init_seqarray(quiet=quiet)
                self._parse_names()
        else:
            raise IPyradError("must enter a data (sequence file) argument.")

        ## check quartets ----------------------------------------------------
        ## depending on the quartet sampling method selected the number of 
        ## quartets that must be sampled will be calculated, or error raised.
        total = n_choose_k(len(self.samples), 4)
        if self.params.method == "all":
            self.params.nquartets = total
        else:
            if not self.params.nquartets:
                self.params.nquartets = int(len(self.samples) ** 2.8)
                print("using default setting for 'random' nquartets = N**2.8 ")

            if self.params.nquartets > total:
                self.params.method = "all"
                print(" Warning: nquartets > total possible quartets "
                + "({})\n Changing to sampling method='all'".format(total))

    ## INIT FUNCTIONS -----------------------------------------------
    def _parse_names(self):
        """ parse sample names from the sequence file"""
        self.samples = []
        with iter(open(self.files.data, 'r')) as infile:
            infile.next().strip().split()
            while 1:
                try:
                    self.samples.append(infile.next().split()[0])
                except StopIteration:
                    break


    def _init_seqarray(self, quiet=False):
        """ 
        Fills the seqarr with the full SNP data set while keeping memory 
        requirements super low, and creates a bootsarr copy with the 
        following modifications:

        1) converts "-" into "N"s, since they are similarly treated as missing. 
        2) randomly resolve ambiguities (RSKWYM)
        3) convert to uint8 for smaller memory load and faster computation. 
        """

        ## read in the data (seqfile)
        spath = open(self.files.data, 'r')
        line = spath.readline().strip().split()
        ntax = int(line[0])
        nbp = int(line[1])
        tmpseq = np.zeros((ntax, nbp), dtype=np.uint8)
        if not quiet:
            print("loading seq array [{} taxa x {} bp]".format(ntax, nbp))        
    
        ## create array storage for original seqarray, the map used for
        ## subsampling unlinked SNPs (bootsmap) and an array that will be
        ## refilled for each bootstrap replicate (bootsarr).
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
            for line, seq in enumerate(spath):
                tmpseq[line] = np.array(list(seq.split()[-1])).view(np.uint8)

            ## convert '-' or '_' into 'N'
            tmpseq[tmpseq == 45] = 78
            tmpseq[tmpseq == 95] = 78            

            ## save array to disk so it can be easily accessed by slicing
            ## This hardly modified array is used again later for sampling boots
            io5["seqarr"][:] = tmpseq

            ## resolve ambiguous IUPAC codes, this is done differently each rep.
            ## everything below here (resolve, index) is done each rep.
            if self.params.resolve_ambigs:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices, nothing else matters.
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## save modified array to disk            
            io5["bootsarr"][:] = tmpseq

        ## cleanup files and mem
        del tmpseq
        spath.close()


    def _refresh(self):
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

        ## store old ipcluster info
        oldcluster = copy.deepcopy(self._ipcluster)

        ## reinit the tetrad object data.
        self.__init__(
            name=self.name, 
            data=self.files.data, 
            mapfile=self.files.mapfile,
            workdir=self.dirs,
            method=self.params.method,
            guidetree=self.files.tree,
            resolve_ambigs=self.params.resolve_ambigs,
            save_invariants=self.params.save_invariants,
            nboots=self.params.nboots, 
            nquartets=self.params.nquartets, 
            initarr=True, 
            quiet=True,
            cli=self.kwargs.get("cli")
            )

        ## retain the same ipcluster info
        self._ipcluster = oldcluster

    ## BOOTSTRAP SEQARRAY SAMPLING FUNCTIONS ------------------------
    def _sample_bootseq_array(self):
        """
        Takes the seqarray and re-samples columns and saves to bootsarr.
        """
        ## use 'r+' to read and write to existing array. This is super 
        ## similar to what is called in __init__. 
        with h5py.File(self.database.input, 'r+') as io5:  
            ## load in the seqarr and maparr
            seqarr = io5["seqarr"][:]

            ## resample columns with replacement
            newarr = np.zeros(seqarr.shape, dtype=np.uint8)
            cols = np.random.randint(0, seqarr.shape[1], seqarr.shape[1])
            tmpseq = shuffle_cols(seqarr, newarr, cols)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self.params.resolve_ambigs:
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
        """
        Re-samples loci with replacement to fill the bootarr sampling
        only a single SNP from each locus according to the maparr. 
        """
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
            if self.params.resolve_ambigs:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## store data sets
            io5.create_dataset("bootsmap", data=tmpmap)
            io5.create_dataset("bootsarr", data=tmpseq)
            #LOGGER.info("resampled bootsarr \n %s", io5["bootsarr"][:, :10])
            #LOGGER.info("resampled bootsmap \n %s", io5["bootsmap"][:10, :])

    ## QUARTET TAXON SAMPLING FUNCTIONS -----------------------------
    def _store_N_samples(self, start, ipyclient, quiet=False):
        """ 
        Find all quartets of samples and store in a large array
        A chunk size is assigned for sampling from the array of quartets
        based on the number of cpus available. This should be relatively 
        large so that we don't spend a lot of time doing I/O, but small 
        enough that jobs finish often for checkpointing.
        """

        breaks = 2
        if self.params.nquartets < 5000:
            breaks = 1
        if self.params.nquartets > 100000:
            breaks = 8
        if self.params.nquartets > 500000:
            breaks = 16
        if self.params.nquartets > 5000000:
            breaks = 32

        ## chunk up the data
        ncpus = len(ipyclient)    
        self._chunksize = (self.params.nquartets // (breaks * ncpus) \
                        + (self.params.nquartets % (breaks * ncpus)))

        ## create h5 OUT empty arrays
        ## 'quartets' stores the inferred quartet relationship (1 x 4)
        ## This can get huge, so we need to choose the dtype wisely. 
        ## the values are simply the index of the taxa, so uint16 is good.
        with h5py.File(self.database.output, 'w') as io5:
            io5.create_dataset("quartets", 
                               (self.params.nquartets, 4), 
                               dtype=np.uint16, 
                               chunks=(self._chunksize, 4))
            ## group for bootstrap invariant matrices ((16, 16), uint32)
            ## these store the actual matrix counts. dtype uint32 can store
            ## up to 4294967295. More than enough. uint16 max is 65535.
            ## the 0 boot is the original seqarray.
            io5.create_group("invariants")

        ## append to h5 IN array (which has the seqarray, bootsarr, maparr)
        ## and fill it with all of the quartet sets we will ever sample.
        ## the samplign method will vary depending on whether this is random, 
        ## all, or equal splits (in a separate but similar function). 
        with h5py.File(self.database.input, 'a') as io5:
            try:
                io5.create_dataset(
                    name="quartets", 
                    shape=(self.params.nquartets, 4), 
                    dtype=np.uint16, 
                    chunks=(self._chunksize, 4),
                    compression='gzip')
            except RuntimeError:
                raise IPyradError(
                    "database file already exists for this analysis, "
                  + "you must run with the force flag to overwrite")
            
        ## submit store job to write into self.database.input
        if self.params.method == "all":
            rasync = ipyclient[0].apply(store_all, self)
        elif self.params.method == "random":
            rasync = ipyclient[0].apply(store_random, self)
        elif self.params.method == "equal":
            rasync = ipyclient[0].apply(store_equal, self) 

        ## progress bar 
        printstr = "generating q-sets | {} | "
        prog = 0        
        while 1:
            elapsed = datetime.timedelta(seconds=int(time.time() - start))
            if not quiet:
                if rasync.stdout:
                    prog = int(rasync.stdout.strip().split()[-1])
                progressbar(self.params.nquartets, prog,
                            printstr.format(elapsed), spacer="")
            if not rasync.ready():
                time.sleep(0.1)
            else:
                break

        if not rasync.successful():
            raise IPyradError(rasync.result())
        if not quiet:
            print("")

    ## QMC QUARTET-JOINING FUNCTIONS ---------------------------------
    def _dump_qmc(self):
        """
        Writes the inferred quartet sets from the database to a text 
        file to be used as input for QMC. Quartets that had no information
        available (i.e., no SNPs) were written to the database as 0,0,0,0
        and are excluded here from the output.
        """

        ## open the h5 database
        with h5py.File(self.database.output, 'r') as io5:

            ## create an output file for writing
            self.files.qdump = os.path.join(self.dirs, self.name+".quartets.txt")
            with open(self.files.qdump, 'w') as qdump:

                ## pull from db
                for idx in range(0, self.params.nquartets, self._chunksize):
                    qchunk = io5["quartets"][idx:idx+self._chunksize, :]
                    quarts = [tuple(j) for j in qchunk if np.any(j)]

                    ## shuffle and format for qmc
                    np.random.shuffle(quarts)
                    chunk = ["{},{}|{},{}".format(*i) for i in quarts]
                    qdump.write("\n".join(chunk)+"\n")


    def _run_qmc(self, boot):
        """
        Runs quartet max-cut QMC on the quartets qdump file.
        """

        ## build command
        self._tmp = os.path.join(self.dirs, ".tmptre")
        cmd = [ip.bins.qmc, "qrtt=" + self.files.qdump, "otre=" + self._tmp]

        ## run it
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        res = proc.communicate()
        if proc.returncode:
            raise IPyradWarningExit(res[1])

        ## parse tmp file written by qmc into a tree and rename it
        with open(self._tmp, 'r') as intree:
            tre = ete3.Tree(intree.read().strip())
            names = tre.get_leaves()
            for name in names:
                name.name = self.samples[int(name.name)]
            tmptre = tre.write(format=9)

        ## save the tree to file
        if boot:
            self.trees.boots = os.path.join(self.dirs, self.name+".boots")
            with open(self.trees.boots, 'a') as outboot:
                outboot.write(tmptre + "\n")
        else:
            self.trees.tree = os.path.join(self.dirs, self.name+".tree")
            with open(self.trees.tree, 'w') as outtree:
                outtree.write(tmptre)

        ## save the file
        self._save()


    def _compute_stats(self, start, ipyclient, quiet=False):
        "Compute sampling stats and consens trees"
        
        ## get name indices
        names = self.samples

        ## make a consensus from bootstrap reps.
        if self.checkpoint.boots:
            tre = ete3.Tree(self.trees.tree, format=0)
            tre.unroot()
            with open(self.trees.boots, 'r') as inboots:
                bb = [ete3.Tree(i.strip(), format=0) for i in inboots.readlines()]
                bb = [tre] + bb

            ## calculate consensus supports
            ctre, counts = consensus_tree(bb, names=names)
            self.trees.cons = os.path.join(self.dirs, self.name + ".cons")
            with open(self.trees.cons, 'w') as ocons:
                ocons.write(ctre.write(format=0))

        else:
            ctre = ete3.Tree(self.trees.tree, format=0)
            ctre.unroot()

        ## build stats file and write trees
        self.trees.nhx = os.path.join(self.dirs, self.name + ".nhx")
        lbview = ipyclient.load_balanced_view()
        qtots = {}
        qsamp = {}
        tots = sum(1 for i in ctre.iter_leaves())
        totn = set(ctre.get_leaf_names())

        ## iterate over node traversal
        for node in ctre.traverse():
            ## this is slow, needs to look at every sampled quartet
            ## so we send it to be processed on engines
            qtots[node] = lbview.apply(get_total, *(tots, node))
            qsamp[node] = lbview.apply(get_sampled, *(self, totn, node))

        ## wait for jobs to finish (+1 to lenjob is for final progress printer)
        alljobs = qtots.values() + qsamp.values()
        lenjobs = len(alljobs) + 1
        printstr = "calculating stats | {} | "
        done = 0
        while 1:
            if not quiet:
                done = sum([i.ready() for i in alljobs])
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(lenjobs, done, 
                    printstr.format(elapsed), spacer="")
            if (lenjobs - 1) == done:
                break
            else:
                time.sleep(0.1)
        ## store results in the tree object
        for node in ctre.traverse():
            total = qtots[node].result()
            sampled = qsamp[node].result()
            node.add_feature("quartets_total", total)
            node.add_feature("quartets_sampled", sampled)
        features = ["quartets_total", "quartets_sampled"]

        ## update final progress
        elapsed = datetime.timedelta(seconds=int(time.time()-start))        
        progressbar(1, 1, printstr.format(elapsed), spacer="")
        if not quiet:
            print("")

        ## write tree in NHX format 
        with open(self.trees.nhx, 'w') as outtre:
            outtre.write(ctre.write(format=0, features=features))


    def _save(self):
        """
        Save a JSON serialized tetrad instance to continue from a checkpoint.
        """

        ## save each attribute as dict
        fulldict = copy.deepcopy(self.__dict__)
        for i, j in fulldict.items():
            if isinstance(j, Params):
                fulldict[i] = j.__dict__
        fulldumps = json.dumps(
            fulldict,
            sort_keys=False, 
            indent=4, 
            separators=(",", ":"),
            )

        ## save to file, make dir if it wasn't made earlier
        assemblypath = os.path.join(self.dirs, self.name + ".tet.json")
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


    def _load(self, name, workdir, quiet=False):
        "Load JSON serialized tetrad instance to continue from a checkpoint."

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
            raise IPyradError("""\
        Cannot find checkpoint (.tet.json) file at: {}""".format(path))

        ## set old attributes into new tetrad object
        self.name = fullj["name"]
        self.files.data = fullj["files"]["data"]
        self.files.mapfile = fullj["files"]["mapfile"]        
        self.dirs = fullj["dirs"]
        self._init_seqarray(quiet=quiet)
        self._parse_names()

        ## fill in the same attributes
        for key in fullj:
            ## fill Params a little different
            if key in ["files", "params", "database", 
                       "trees", "stats", "checkpoint"]:
                filler = fullj[key]
                for ikey in filler:
                    self.__dict__[key].__setattr__(ikey, fullj[key][ikey])
            else:
                self.__setattr__(key, fullj[key])


    def _inference(self, start, ipyclient, quiet):
        """
        Sends slices of quartet sets to parallel engines for computing, 
        enters results into output database, sends finished quartet sets
        to QMC for tree inference, and prints progress bars.
        """

        ## load-balancer for single-threaded execution jobs
        lbview = ipyclient.load_balanced_view()

        ## an iterator that grabs quartet chunk start positions
        jobs = range(self.checkpoint.arr, self.params.nquartets, self._chunksize)

        ## if this is a bootstrap then init a new boot array in the database
        ## max val is 65535 in here if uint16
        bootkey = "boot{}".format(self.checkpoint.boots)
        with h5py.File(self.database.output, 'r+') as io5:
            if bootkey not in io5["invariants"].keys():
                io5["invariants"].create_dataset(
                    bootkey, 
                    (self.params.nquartets, 16, 16),
                    dtype=np.uint16,
                    chunks=(self._chunksize, 16, 16))

        ## start progress bar if new or skip if bootstrapping
        elapsed = datetime.timedelta(seconds=int(time.time() - start))
        if self.checkpoint.boots:
            printstr = "bootstrap trees   | {} | "
        else:
            printstr = "initial tree      | {} | "
            if not quiet:
                progressbar(1, 0, printstr.format(elapsed), spacer="")

        ## submit jobs distriuted across the cluster.
        asyncs = {}
        for job in jobs:
            asyncs[job] = lbview.apply(nworker, *(self, job))

        ## wait for jobs to finish, catch results as they return and
        ## enter them into the HDF5 database to keep memory low.
        done = 0
        while 1:
            ## gather finished jobs
            finished = [i for i, j in asyncs.iteritems() if j.ready()]

            ## iterate over finished list
            for key in finished:
                rasync = asyncs[key]
                if rasync.successful():
                    ## store result
                    done += 1
                    results = rasync.result()
                    self._insert_to_array(key, results)
                    ## purge from memory
                    del asyncs[key]
                else:
                    raise IPyradError(rasync.result())

            # progress bar is different if first vs boot tree
            if not quiet:
                if not self.checkpoint.boots:
                    progressbar(len(jobs), done, start, printstr)
                else:
                    progressbar(
                        self.params.nboots, 
                        self.checkpoint.boots, 
                        printstr, 
                        start)

            ## done is counted on finish, so this means we're done
            if len(asyncs) == 0:
                break
            else:
                time.sleep(0.1)

        ## dump quartets into a text file for QMC
        self._dump_qmc()

        ## send to QMC
        if not self.checkpoint.boots:
            self._run_qmc(0)
        else:
            self._run_qmc(1)

        ## reset the checkpoint arr
        self.checkpoint.arr = 0

        ## print spacer if finished first tree or last boot.
        if (not self.checkpoint.boots) and (not quiet):
            print("")

        elif (self.checkpoint.boots == self.params.nboots) and (not quiet):
            print("")


    def _insert_to_array(self, chunk, results):
        "Enters results arrays into the HDF5 database."

        ## two result arrs
        chunksize = self._chunksize
        qrts, invs = results

        ## enter into db
        with h5py.File(self.database.output, 'r+') as io5:
            io5['quartets'][chunk:chunk + chunksize] = qrts

            ## entered as 0-indexed !
            if self.params.save_invariants:
                if self.checkpoint.boots:
                    key = "invariants/boot{}".format(self.checkpoint.boots)
                    io5[key][chunk:chunk + chunksize] = invs
                else:
                    io5["invariants/boot0"][chunk:chunk + chunksize] = invs


    def run(self, force=False, quiet=False, ipyclient=None):
        """
        Parameters
        ----------
        force (bool):
            Overwrite existing results for object with the same name
            and workdir as this one.
        verbose (int):
            0=primt nothing; 1=print progress bars; 2=print pringress
            bars and cluster information.
        ipyclient (ipyparallel.Client object):
            A connected ipyclient object. If ipcluster instance is 
            not running on the default profile then ...
        """

        ## force overwrite needs to clear out the HDF5 database
        if force:
            self._refresh()

        ## print nquartet statement
        if not quiet:
            print("inferring {} quartet tree sets".format(self.params.nquartets))

        ## wrap the run in a try statement to ensure we properly shutdown
        ## and cleanup on exit or interrupt. 
        inst = None
        try:
            ## find and connect to an ipcluster instance given the information
            ## in the _ipcluster dictionary if a connected client was not given.
            if not ipyclient:
                args = self._ipcluster.items() + [("spacer", "")]
                ipyclient = ip.core.parallel.get_client(**dict(args))

            ## print the cluster connection information
            if not quiet:
                ip.cluster_info(ipyclient)

            ## store ipyclient engine pids to the dict so we can 
            ## hard-interrupt them later if assembly is interrupted. 
            ## Only stores pids of engines that aren't busy at this moment, 
            ## otherwise it would block here while waiting to find their pids.
            self._ipcluster["pids"] = {}
            for eid in ipyclient.ids:
                engine = ipyclient[eid]
                if not engine.outstanding:
                    pid = engine.apply(os.getpid).get()
                    self._ipcluster["pids"][eid] = pid            

            ## fill the input array with quartets to sample --------------------
            start = time.time()
            if (not self.checkpoint.boots) and (not self.trees.tree):
                self._store_N_samples(start, ipyclient, quiet=quiet)

            ## calculate invariants for the full seqarray ----------------------
            start = time.time()            
            if not self.trees.tree:
                self._inference(start, ipyclient, quiet=quiet)
            else:
                if not quiet:
                    print("initial tree already inferred")

            ## calculate invariants for each bootstrap rep ----------------------
            start = time.time()
            if self.params.nboots:
                if self.checkpoint.boots: #<= self.params.nboots:
                    if not quiet:
                        print("{} bootstrap trees already inferred"\
                              .format(self.checkpoint.boots))

                while self.checkpoint.boots < self.params.nboots:
                    ## resample the bootseq array
                    if self.files.mapfile:
                        self._sample_bootseq_array_map()
                    else:
                        self._sample_bootseq_array()

                    ## start boot inference 
                    self.checkpoint.boots += 1
                    self._inference(start, ipyclient, quiet=quiet)

            ## write output stats -----------------------------------------------
            #self.files.stats = os.path.join(self.dirs, self.name+"_stats.txt")
            start = time.time()
            self._compute_stats(start, ipyclient, quiet=quiet)

        ## handle exceptions so they will be raised after we clean up below
        except KeyboardInterrupt as inst:
            print("\nKeyboard Interrupt by user. Cleaning up...")

        except IPyradError as inst:
            print("\nError encountered: {}".format(inst))

        except Exception as inst:
            print("\nUnknown exception encountered: {}".format(inst))

        ## close client when done or interrupted
        finally:
            try:
                ## save the Assembly
                self._save()                
                
                ## can't close client if it was never open
                if ipyclient:

                    ## send SIGINT (2) to all engines
                    ipyclient.abort()
                    time.sleep(1)
                    for engine_id, pid in self._ipcluster["pids"].items():
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                        time.sleep(0.25)
                    
                    ## if CLI, stop jobs and shutdown
                    if 'ipyrad-cli' in self._ipcluster["cluster_id"]:
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                    else:
                        if not ipyclient.outstanding:
                            ipyclient.purge_everything()
                        else:
                            ## nanny: kill everything, something bad happened
                            ipyclient.shutdown(hub=True, block=False)
                            ipyclient.close()
                            print("\nwarning: ipcluster shutdown and must be restarted")
                
                ## reraise the error now that we're cleaned up
                #if inst:
                #    raise inst

            ## if exception during shutdown then we really screwed up
            except Exception as inst2:
                print("warning: error during shutdown:\n{}".format(inst2))



#################################################
## worker function to run jit funcs
#################################################

def nworker(data, chunk):
    """
    Worker to distribute work to jit funcs. Wraps everything on an 
    engine to run single-threaded to maximize efficiency for 
    multi-processing.
    """

    ## set the thread limit on the remote engine
    oldlimit = set_mkl_thread_limit(1)

    ## open seqarray view, the modified arr is in bootstarr
    with h5py.File(data.database.input, 'r') as io5:
        seqview = io5["bootsarr"][:]
        maparr = io5["bootsmap"][:, 0]
        smps = io5["quartets"][chunk:chunk + data._chunksize]

        ## create an N-mask array of all seq cols
        nall_mask = seqview[:] == 78

    ## init arrays to fill with results
    rquartets = np.zeros((smps.shape[0], 4), dtype=np.uint16)
    rinvariants = np.zeros((smps.shape[0], 16, 16), dtype=np.uint16)

    ## fill arrays with results as we compute them. This iterates
    ## over all of the quartet sets in this sample chunk. It would
    ## be nice to have this all numbified.
    for idx in range(smps.shape[0]):
        sidx = smps[idx]
        seqs = seqview[sidx]

        ## these axis calls cannot be numbafied, but I can't 
        ## find a faster way that is JIT compiled, and I've
        ## really, really, really tried. Tried again now that
        ## numba supports axis args for np.sum. Still can't 
        ## get speed improvements by numbifying this loop.
        nmask = np.any(nall_mask[sidx], axis=0)
        nmask += np.all(seqs == seqs[0], axis=0) 

        ## here are the jitted funcs
        bidx, invar = calculate(seqs, maparr, nmask, TESTS)

        ## store results
        rquartets[idx] = smps[idx][bidx]
        rinvariants[idx] = invar

    ## reset thread limit
    set_mkl_thread_limit(oldlimit)

    ## return results...
    return rquartets, rinvariants


#################################################
## jit compiled invariant count functions
#################################################


@numba.jit(nopython=True)
def calculate(seqnon, mapcol, nmask, tests):
    """
    Groups together several other numba funcs.
    """

    ## get the invariants matrix
    mats = chunk_to_matrices(seqnon, mapcol, nmask)

    ## empty arrs to fill
    svds = np.zeros((3, 16), dtype=np.float64)
    scor = np.zeros(3, dtype=np.float64)
    rank = np.zeros(3, dtype=np.float64)

    ## why svd and rank?
    for test in range(3):
        svds[test] = np.linalg.svd(mats[test].astype(np.float64))[1]
        rank[test] = np.linalg.matrix_rank(mats[test].astype(np.float64))

    ## get minrank, or 11
    minrank = int(min(11, rank.min()))
    for test in range(3):
        scor[test] = np.sqrt(np.sum(svds[test, minrank:]**2))

    ## sort to find the best qorder
    best = np.where(scor == scor.min())[0]
    bidx = tests[best][0]

    return bidx, mats[0]



@numba.jit('u4[:,:,:](u1[:,:],u4[:],b1[:])', nopython=True)
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
    for idx in range(mapcol.shape[0]):
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


##################################################
## quartet sampling funcs to fill the database
##################################################


def store_all(self):
    """
    Populate array with all possible quartets. This allows us to 
    sample from the total, and also to continue from a checkpoint
    """
    with h5py.File(self.database.input, 'a') as io5:
        fillsets = io5["quartets"]

        ## generator for all quartet sets
        qiter = itertools.combinations(range(len(self.samples)), 4)
        i = 0
        while i < self.params.nquartets:
            ## sample a chunk of the next ordered N set of quartets
            dat = np.array(list(itertools.islice(qiter, self._chunksize)))
            end = min(self.params.nquartets, dat.shape[0] + i)
            fillsets[i:end] = dat[:end - i]
            i += self._chunksize

            ## send progress update to stdout on engine
            print(min(i, self.params.nquartets))


def store_random(self):
    """
    Populate array with random quartets sampled from a generator.
    Holding all sets in memory might take a lot, but holding a very
    large list of random numbers for which ones to sample will fit 
    into memory for most reasonable sized sets. So we'll load a 
    list of random numbers in the range of the length of total 
    sets that can be generated, then only keep sets from the set 
    generator if they are in the int list. I did several tests to 
    check that random pairs are as likely as 0 & 1 to come up together
    in a random quartet set. 
    """

    with h5py.File(self.database.input, 'a') as io5:
        fillsets = io5["quartets"]

        ## set generators
        qiter = itertools.combinations(range(len(self.samples)), 4)
        rand = np.arange(0, n_choose_k(len(self.samples), 4))
        np.random.shuffle(rand)
        rslice = rand[:self.params.nquartets]
        rss = np.sort(rslice)
        riter = iter(rss)
        del rand, rslice

        ## print progress update 1 to the engine stdout
        print(self._chunksize)

        ## set to store
        rando = riter.next()
        tmpr = np.zeros((self.params.nquartets, 4), dtype=np.uint16)
        tidx = 0
        while 1:
            try:
                for i, j in enumerate(qiter):
                    if i == rando:
                        tmpr[tidx] = j
                        tidx += 1
                        rando = riter.next()

                    ## print progress bar update to engine stdout
                    if not i % self._chunksize:
                        print(min(i, self.params.nquartets))

            except StopIteration:
                break
        ## store into database
        fillsets[:] = tmpr
        del tmpr


def store_equal(self):
    """
    Takes a tetrad class object and populates array with random 
    quartets sampled equally among splits of the tree so that 
    deep splits are not overrepresented relative to rare splits, 
    like those near the tips. 
    """

    with h5py.File(self.database.input, 'a') as io5:
        fillsets = io5["quartets"]

        ## require guidetree
        if not os.path.exists(self.files.tree):
            raise IPyradError(
                "To use sampling method 'equal' requires a guidetree")
        tre = ete3.Tree(self.files.tree)
        tre.unroot()
        tre.resolve_polytomy(recursive=True)

        ## randomly sample internals splits
        splits = [([self.samples.index(z.name) for z in i],
                   [self.samples.index(z.name) for z in j]) \
                   for (i, j) in tre.get_edges()]

        ## only keep internal splits, not single tip edges
        splits = [i for i in splits if all([len(j) > 1 for j in i])]

        ## how many min quartets shoudl be equally sampled from each split
        squarts = self.params.nquartets // len(splits)

        ## keep track of how many iterators are saturable.
        saturable = 0

        ## turn each into an iterable split sampler
        ## if the nquartets for that split is small, then sample all, 
        ## if it is big then make it a random sampler for that split.
        qiters = []

        ## iterate over splits sampling quartets evenly
        for idx, split in enumerate(splits):
            ## if small number at this split then sample all possible sets
            ## we will exhaust this quickly and then switch to random for 
            ## the larger splits.
            total = n_choose_k(len(split[0]), 2) * n_choose_k(len(split[1]), 2)
            if total < squarts*2:
                qiter = (i + j for (i, j) in itertools.product(
                    itertools.combinations(split[0], 2), 
                    itertools.combinations(split[1], 2)))
                saturable += 1

            ## else create random sampler across that split, this is slower
            ## because it can propose the same split repeatedly and so we 
            ## have to check it against the 'sampled' set.
            else:
                qiter = (random_product(split[0], split[1]) for _ \
                         in range(self.params.nquartets))

            ## store all iterators into a list
            qiters.append((idx, qiter))

        ## create infinite cycler of qiters
        qitercycle = itertools.cycle(qiters)

        ## store visited quartets
        sampled = set()

        ## fill chunksize at a time
        i = 0
        empty = set()
        edge_targeted = 0
        random_targeted = 0

        ## keep filling quartets until nquartets are sampled.
        while i < self.params.nquartets:
            ## grab the next iterator
            cycle, qiter = qitercycle.next()

            ## sample from iterators, store sorted set.
            try:
                qrtsamp = tuple(sorted(qiter.next()))
                if qrtsamp not in sampled:
                    sampled.add(qrtsamp)
                    edge_targeted += 1
                    i += 1
                    ## print progress bar update to engine stdout
                    if not i % self._chunksize:
                        print(min(i, self.params.nquartets))                    

            except StopIteration:
                empty.add(cycle)
                if len(empty) == saturable:
                    break


        ## if array is not full then add random samples
        while i <= self.params.nquartets:
            newset = tuple(sorted(np.random.choice(
                range(len(self.samples)), 4, replace=False)))
            if newset not in sampled:
                sampled.add(newset)
                random_targeted += 1
                i += 1
                ## print progress bar update to engine stdout
                if not i % self._chunksize:
                    print(min(i, self.params.nquartets))

        ## store into database
        print(self.params.nquartets)
        fillsets[:] = np.array(tuple(sampled))
        del sampled



########################################
## some itertools cookbook recipes
########################################


def n_choose_k(n, k):
    """ 
    Get the number of quartets as n-choose-k. This is used in equal 
    splits to decide whether a split should be exhaustively sampled
    or randomly sampled. Edges near tips can be exhaustive while highly 
    nested edges can do with less sampling.
    """
    mulfunc = lambda x, y: x * y
    return int(reduce(mulfunc, (Fraction(n  -i, i+1) for i in range(k)), 1))


def random_combination(nsets, n, k):
    """
    Returns nsets unique random quartet sets sampled from
    n-choose-k without replacement combinations.
    """
    sets = set()
    while len(sets) < nsets:
        newset = tuple(sorted(np.random.choice(n, k, replace=False)))
        sets.add(newset)
    return tuple(sets)


def random_product(iter1, iter2):
    """ 
    Random sampler for equal_splits functions
    """
    iter4 = np.concatenate([
        np.random.choice(iter1, 2, replace=False),
        np.random.choice(iter2, 2, replace=False)
        ])
    return iter4


########################################
## JIT compiled parsing functions
########################################

@numba.jit(nopython=True)
def get_spans(maparr, spans):
    """ 
    Get span distance for each locus in original seqarray. This
    is used to create re-sampled arrays in each bootstrap to sample
    unlinked SNPs. 
    """
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


@numba.jit(nopython=True)
def resolve_ambigs(tmpseq):
    """ 
    Randomly resolve ambiguous bases. This is applied to each boot
    replicate so that over reps the random resolutions don't matter.
    Sites are randomly resolved, so best for unlinked SNPs since 
    otherwise linked SNPs are losing their linkage information... 
    though it's not like we're using it anyways.
    """

    ## the order of rows in GETCONS
    for aidx in xrange(6):
        #np.uint([82, 75, 83, 89, 87, 77]):
        ambig, res1, res2 = GETCONS[aidx]

        ## get true wherever tmpseq is ambig
        idx, idy = np.where(tmpseq == ambig)
        halfmask = np.random.choice(np.array([True, False]), idx.shape[0])

        for col in xrange(idx.shape[0]):
            if halfmask[col]:
                tmpseq[idx[col], idy[col]] = res1
            else:
                tmpseq[idx[col], idy[col]] = res2
    return tmpseq


@numba.jit(nopython=True)#, cache=True)
def shuffle_cols(seqarr, newarr, cols):
    """ 
    Used in bootstrap resampling when no map file is present.
    """
    for idx in xrange(cols.shape[0]):
        newarr[:, idx] = seqarr[:, cols[idx]]
    return newarr


@numba.jit(nopython=True)#, cache=True)
def get_shape(spans, loci):
    """ 
    Returns shape of new bootstrap resampled locus array. The 
    shape can change because when we resample loci the number
    of SNPs in the data set can change.
    """
    width = 0
    for idx in xrange(loci.shape[0]):
        width += spans[loci[idx], 1] - spans[loci[idx], 0]
    return width


@numba.jit(nopython=True)#, cache=True)
def fill_boot(seqarr, newboot, newmap, spans, loci):
    """ 
    Fills the new bootstrap SNP array and map array with
    new data based on the resampled loci for this boot.
    """
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


#########################################
## globals
#########################################

## used by resolve_ambigs
GETCONS = np.array([[82, 71, 65],
                    [75, 71, 84],
                    [83, 71, 67],
                    [89, 84, 67],
                    [87, 84, 65],
                    [77, 67, 65]], dtype=np.uint8)


## the three indexed resolutions of each quartet
TESTS = np.array([[0, 1, 2, 3], 
                  [0, 2, 1, 3], 
                  [0, 3, 1, 2]], dtype=np.uint8)


##########################################
## custom exception messages
##########################################

MIDSTREAM_MESSAGE = """
    loaded object method={}
    cannot change sampling methods midstream
    use force argument to start new run with new method
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


############################################
## Function to limit threading to single
############################################

def set_mkl_thread_limit(cores):
    """
    set mkl thread limit and return old value so we can reset
    when finished. 
    """
    if "linux" in sys.platform:
        mkl_rt = ctypes.CDLL('libmkl_rt.so')
    else:
        mkl_rt = ctypes.CDLL('libmkl_rt.dylib')
    oldlimit = mkl_rt.mkl_get_max_threads()
    mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(cores)))
    return oldlimit


#############################################
## Tree and statistics operations
#############################################


def get_total(tots, node):
    """ get total number of quartets possible for a split"""
    if (node.is_leaf() or node.is_root()):
        return 0
    else:
        ## Get counts on down edges. 
        ## How to treat polytomies here?
        if len(node.children) > 2:
            down_r = node.children[0]
            down_l = node.children[1]
            for child in node.children[2:]:
                down_l += child
        else:
            down_r, down_l = node.children
        lendr = sum(1 for i in down_r.iter_leaves())
        lendl = sum(1 for i in down_l.iter_leaves())

        ## get count on up edge sister
        up_r = node.get_sisters()[0]
        lenur = sum(1 for i in up_r.iter_leaves())

        ## everyone else
        lenul = tots - (lendr + lendl + lenur)

        ## return product
        return lendr * lendl * lenur * lenul



def get_sampled(data, totn, node):
    """ get total number of quartets sampled for a split"""
    ## convert tip names to ints
    names = sorted(totn)
    cdict = {name: idx for idx, name in enumerate(names)}
    
    ## skip some nodes
    if (node.is_leaf() or node.is_root()):
        return 0
    else:
        ## get counts on down edges
        if len(node.children) > 2:
            down_r = node.children[0]
            down_l = node.children[1]
            for child in node.children[2:]:
                down_l += child
        else:
            down_r, down_l = node.children

        lendr = set(cdict[i] for i in down_r.get_leaf_names())
        lendl = set(cdict[i] for i in down_l.get_leaf_names())

        ## get count on up edge sister
        up_r = node.get_sisters()[0]
        lenur = set(cdict[i] for i in up_r.get_leaf_names())

        ## everyone else
        lenul = set(cdict[i] for i in totn) - set.union(lendr, lendl, lenur)

    idx = 0
    sampled = 0
    with h5py.File(data.database.output, 'r') as io5:
        end = io5["quartets"].shape[0]
        while 1:
            ## break condition
            if idx >= end:
                break

            ## counts matches
            qrts = io5["quartets"][idx:idx+data._chunksize]
            for qrt in qrts:
                sqrt = set(qrt)
                if all([sqrt.intersection(i) for i in [lendr, lendl, lenur, lenul]]):
                    sampled += 1

            ## increase span
            idx += data._chunksize
    return sampled



def consensus_tree(trees, names=None, cutoff=0.0):
    """ 
    An extended majority rule consensus function for ete3. 
    Modelled on the similar function from scikit-bio tree module. If 
    cutoff=0.5 then it is a normal majority rule consensus, while if 
    cutoff=0.0 then subsequent non-conflicting clades are added to the tree.
    """

    ## find which clades occured with freq > cutoff
    namedict, clade_counts = find_clades(trees, names=names)

    ## filter out the < cutoff clades
    fclade_counts = filter_clades(clade_counts, cutoff)

    ## build tree
    consens_tree, _ = build_trees(fclade_counts, namedict)
    ## make sure no singleton nodes were left behind
    return consens_tree, clade_counts



def filter_clades(clade_counts, cutoff):
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



def find_clades(trees, names):
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
            #bits = bitarray('0'*len(tree))
            bits = np.zeros(len(tree), dtype=np.bool_)
            for child in node.iter_leaf_names():
                bits[ndict[child]] = True
            ## if parent is root then mirror flip one child (where bit[0]=0)
            # if not node.is_root():
            #     if node.up.is_root():
            #         if bits[0]:
            #             bits.invert()
            bitstring = "".join([np.binary_repr(i) for i in bits])
            clade_counts[bitstring] += 1

    ## convert to freq
    for key, val in clade_counts.items():
        clade_counts[key] = val / float(len(trees))

    ## return in sorted order
    clade_counts = sorted(clade_counts.items(), 
                          key=lambda x: x[1],
                          reverse=True)
    return namedict, clade_counts



def build_trees(fclade_counts, namedict):
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

    
###############################################
## Save/Load from JSON operations checkpoint.
###############################################

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

