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
from builtins import range, str

# standard lib
import os
import sys
import json
import h5py
import time
import copy
import itertools
import subprocess as sps

# third party
import toytree
import numba
import ctypes
import datetime
import numpy as np
from scipy.special import comb

import ipyrad as ip
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import Params, get_spans



class Tetrad(object):
    """
    The main tetrad object for saving/loading data and checkpointing with JSON,
    connecting to parallel client, and storing results for eacy access.

    Params: 
        name (str): a string to use as prefix for outputs.
        data (str): a phylip formatted file (.snps or .phy from ipyrad).
        mapfile (str): [opt] a mapfile from ipyrad (.snpsmap).
        method (str): [opt; def=all] 'all', 'random', or 'equal'.
        guidetree (str): [opt] a newick file for guiding 'equal' method.
        nquartets (int): [opt] number of samples for 'random' method.
        nboots (int): [opt] number of non-parametric bootstrap replicates.
        resolve_ambigs (bool): [def=True] Whether to include ambiguous sites.
        load (bool): if True object is loaded from [workdir]/[name].json.
        quiet (bool): if True progress is not printed to stdout.
        save_invariants (bool): if True invariants array is written to hdf5.

    Functions:
        run()

    Attributes:
        params : optional params can be set after object instantiation
        trees : access trees from object after analysis is finished
        samples: names of samples in the data set
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
        save_invariants=False,
        *args, 
        **kwargs):      

        # check additional (hidden) arguments from kwargs.
        self.quiet = False
        self.kwargs = {
            "initarr": True,
            "cli": False,
            }
        self.kwargs.update(kwargs)

        # are we in the CLI?
        self._cli = False
        if self.kwargs.get("cli"):
            self._cli = True

        # name, sample, and workdir
        self.name = name
        self.samples = []
        self.dirs = os.path.abspath(os.path.expanduser(workdir))
        if not os.path.exists(self.dirs):
            os.mkdir(self.dirs)

        # required arguments 
        params = (method, nboots, nquartets, resolve_ambigs, save_invariants)
        self._init_params(*params)
        self._init_database()
        self._init_file_handles(data, mapfile, guidetree)

        # is this a new analysis, or loading an existing object? If loading 
        # then checkpoints and file handles are all pulled in. 
        if load:
            self._load(self.name, self.dirs)
        else:
            #self._init_file_handles(data, mapfile, guidetree)
            self._parse_names()
            if self.kwargs["initarr"]:
                self._init_seqarray()

        # default ipcluster information for finding a running Client
        self.ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }

        # Sampling method attributes (should I put path info in here too?)

        # private attributes
        self._chunksize = 1
        self._tmp = None

        # self.populations ## if we allow grouping samples
        # (haven't done this yet)
        self._count_quartets()

        # stats is written to os.path.join(self.dirs, self.name+".stats.txt")
        self.stats = Params()
        self.stats.n_quartets_sampled = self.params.nquartets
        #self.stats.avg

        # checkpointing information
        self.checkpoint = Params()
        self.checkpoint.boots = 0
        self.checkpoint.arr = 0


    @property
    def _spacer(self):
        """ return print spacer for CLI versus API """
        if self._cli:
            return "  "
        return ""


    def _init_params(self, method, nboots, nquartets, resolve_ambigs, save_invariants):
        self.params = Params()
        self.params.method = method
        self.params.nboots = nboots
        self.params.nquartets = nquartets
        self.params.resolve_ambigs = resolve_ambigs
        self.params.save_invariants = save_invariants


    def _init_database(self):
        self.database = Params()
        self.database.input = os.path.join(
            self.dirs, self.name + ".input.h5")
        self.database.output = os.path.join(
            self.dirs, self.name + ".output.h5")        


    def _init_file_handles(self, data, mapfile, guidetree):
        "set handles and clear old data unless loaded existing data"

        # input files
        self.files = Params()       
        self.files.qdump = None
        self.files.data = os.path.abspath(data)
        self.files.mapfile = None
        self.files.guidetree = None
        if mapfile:
            self.files.mapfile = os.path.abspath(mapfile)
        if guidetree:
            self.files.guidetree = os.path.abspath(guidetree)

        # set tree file paths
        self.trees = Params()
        self.trees.tree = os.path.join(self.dirs, self.name + ".tre")
        self.trees.cons = os.path.join(self.dirs, self.name + ".cons.tre")
        self.trees.boots = os.path.join(self.dirs, self.name + ".boots.tre")        
        self.trees.nhx = os.path.join(self.dirs, self.name + ".nhx.tre")

        # check user supplied paths:
        for sfile in self.files:
            if self.files[sfile]:
                if not os.path.exists(self.files[sfile]):
                    raise IOError(
                        "file path not found: {}".format(self.files[sfile])
                        )

        # if not loading files
        if not (data and os.path.exists(data)):
            raise IPyradError(
                "must enter a data (e.g., snp file) argument or use load=True")

        # remove any existing results files
        for sfile in self.trees:
            if self.trees[sfile]:
                if os.path.exists(self.trees[sfile]):
                    os.remove(self.trees[sfile])


    def _count_quartets(self):
        """
        Depending on the quartet sampling method selected the number of 
        quartets that must be sampled will be calculated, or error raised.
        """
        total = int(comb(len(self.samples), 4))
        if self.params.method == "all":
            self.params.nquartets = total
        else:
            if not self.params.nquartets:
                self.params.nquartets = int(len(self.samples) ** 2.8)
                self._print(
                    "using default setting for 'random' nquartets = N**2.8 ")

            if self.params.nquartets > total:
                self.params.method = "all"
                self._print(
                    " Warning: nquartets > total possible quartets " + 
                    + "({})\n Changing to sampling method='all'".format(total))


    def _load_file_paths(self):
        "load file paths if they exist, or None if empty"
        for key, val in self.trees.__dict__.items():
            if not os.path.exists(val):
                self.trees.__dict__[key] = None


    def _parse_names(self):
        "parse sample names from the sequence file"
        self.samples = []
        with iter(open(self.files.data, 'r')) as infile:
            next(infile)  
            while 1:
                try:
                    self.samples.append(next(infile).split()[0])
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
        # read in the data (seqfile)
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

            # if there is a map file, load it into the bootsmap
            if self.files.mapfile:
                with open(self.files.mapfile, 'r') as inmap:
                    
                    # parse the map file from txt and save as dataset
                    maparr = np.genfromtxt(inmap, dtype=np.uint64)
                    maparr[:, 1] = 0
                    maparr = maparr.astype(int)
                    io5["bootsmap"][:] = maparr[:, [0, 3]]

                    # parse the span info from maparr and save to dataset
                    spans = np.zeros((maparr[-1, 0], 2), dtype=np.uint64)
                    spans = get_spans(maparr, spans)
                    io5.create_dataset("spans", data=spans)
                    if not quiet:
                        print("max unlinked SNPs per quartet (nloci): {}"\
                              .format(spans.shape[0]))
            else:
                io5["bootsmap"][:, 0] = np.arange(io5["bootsmap"].shape[0])

            ## fill the tmp array from the input phy
            for line, seq in enumerate(spath):
                tmpseq[line] = (
                    np.array(list(seq.split()[-1]))
                    .astype(bytes)
                    .view(np.uint8)
                )              

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
        # clear any existing results files
        oldfiles = [self.files.qdump] + \
            list(self.database.__dict__.values()) + \
            list(self.trees.__dict__.values())
        for oldfile in oldfiles:
            if oldfile:
                if os.path.exists(oldfile):
                    os.remove(oldfile)

        # store old ipcluster info
        oldcluster = copy.deepcopy(self.ipcluster)

        # reinit the tetrad object data.
        self.__init__(
            name=self.name, 
            data=self.files.data, 
            mapfile=self.files.mapfile,
            workdir=self.dirs,
            method=self.params.method,
            guidetree=self.files.guidetree,
            resolve_ambigs=self.params.resolve_ambigs,
            save_invariants=self.params.save_invariants,
            nboots=self.params.nboots, 
            nquartets=self.params.nquartets, 
            initarr=True, 
            quiet=True,
            cli=self.kwargs.get("cli")
            )

        # retain the same ipcluster info
        self.ipcluster = oldcluster


    def _store_N_samples(self, ipyclient):
        """ 
        Find all quartets of samples and store in a large array
        A chunk size is assigned for sampling from the array of quartets
        based on the number of cpus available. This should be relatively 
        large so that we don't spend a lot of time doing I/O, but small 
        enough that jobs finish often for checkpointing.
        """
        # chunking
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
        self._chunksize = max([
            1, 
            sum([
                (self.params.nquartets // (breaks * ncpus)),
                (self.params.nquartets % (breaks * ncpus)),
            ])
        ])

        ## create h5 OUT empty arrays
        ## 'quartets' stores the inferred quartet relationship (1 x 4)
        ## This can get huge, so we need to choose the dtype wisely. 
        ## the values are simply the index of the taxa, so uint16 is good.
        with h5py.File(self.database.output, 'w') as io5:
            io5.create_dataset(
                name="quartets", 
                shape=(self.params.nquartets, 4), 
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
                    "database file already exists for this analysis, you must"
                  + "run with the force flag to overwrite")
            
        # submit store job to write into self.database.input
        if self.params.method == "all":
            rasync = ipyclient[0].apply(store_all, self)
        elif self.params.method == "random":
            rasync = ipyclient[0].apply(store_random, self)
        elif self.params.method == "equal":
            rasync = ipyclient[0].apply(store_equal, self) 

        ## progress bar 
        printstr = ("generating q-sets", "")
        prog = 0        
        while 1:
            if rasync.stdout:
                prog = int(rasync.stdout.strip().split()[-1])
            self._progressbar(
                self.params.nquartets, prog, self.start, printstr)

            if not rasync.ready():
                time.sleep(0.1)
            else:
                break

        if not rasync.successful():
            raise IPyradError(rasync.result())
        self._print("")


    def _print(self, message):
        if not self.quiet:
            print(message)


    def _progressbar(self, njobs, finished, start, message):
        # measure progress
        if not self.quiet:
            if njobs:
                progress = 100 * (finished / float(njobs))
            else:
                progress = 100

            # build the bar
            hashes = '#' * int(progress / 5.)
            nohash = ' ' * int(20 - len(hashes))

            # timestamp
            elapsed = datetime.timedelta(seconds=int(time.time() - start))

            # print to stderr
            if self.kwargs["cli"]:
                print("\r[{}] {:>3}% {} | {:<12} ".format(
                    hashes + nohash,
                    int(progress),
                    elapsed,
                    message[0],
                ), end="")
            else:
                print("\r[{}] {:>3}% {} | {:<12} ".format(*[
                    hashes + nohash,
                    int(progress),
                    elapsed,
                    message[0],
                ]), end="")
            sys.stdout.flush()


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


    def _load(self, name, workdir="analysis-tetrad"):
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
        self._init_seqarray()
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


    def _get_parallel(self, ipyclient):
        if not ipyclient:
            #args = list(self.ipcluster.items()) + [("spacer", "")]
            ipyclient = ip.core.parallel.get_client(self) # **dict(args))

            # if THAT fails, launch a custom ipcluster in cli mode
            # self.ipcluster["cluster_id"] = 'ipyrad-cli'
            # ... cli = True
            # ... ipcluster_register()

        # print the cluster connection information
        if not self.quiet:
            ip.cluster_info(ipyclient)

        ## store ipyclient engine pids to the dict so we can 
        ## hard-interrupt them later if assembly is interrupted. 
        ## Only stores pids of engines that aren't busy at this moment, 
        ## otherwise it would block here while waiting to find their pids.
        self.ipcluster["pids"] = {}
        for eid in ipyclient.ids:
            engine = ipyclient[eid]
            if not engine.outstanding:
                pid = engine.apply(os.getpid).get()
                self.ipcluster["pids"][eid] = pid
        return ipyclient


    def _cleanup_parallel(self, ipyclient):
        try:
            ## save the Assembly
            self._save()                
            
            ## can't close client if it was never open
            if ipyclient:

                ## send SIGINT (2) to all engines
                ipyclient.abort()
                time.sleep(1)
                for engine_id, pid in self.ipcluster["pids"].items():
                    if ipyclient.queue_status()[engine_id]["tasks"]:
                        os.kill(pid, 2)
                    time.sleep(0.25)
                
                ## if CLI, stop jobs and shutdown
                if 'ipyrad-cli' in self.ipcluster["cluster_id"]:
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

        ## if exception during shutdown then we really screwed up
        except Exception as inst2:
            print("warning: error during shutdown:\n{}".format(inst2))


    def run(self, force=False, quiet=False, ipyclient=None):
        """
        Parameters
        ----------
        force (bool):
            Overwrite existing results for object with the same name
            and workdir as this one.
        quiet (int):
            0=primt nothing; 1=print progress bars; 2=print pringress
            bars and cluster information.
        ipyclient (ipyparallel.Client object):
            A connected ipyclient object. If ipcluster instance is 
            not running on the default profile then ...
        """
        # force overwrite needs to clear out the HDF5 database
        self.quiet = quiet
        if force:
            self._refresh()

        # print nquartet statement
        self._print(
            "inferring {} quartet tree sets".format(self.params.nquartets))

        # wrap the run in a try statement to ensure we properly shutdown
        try:
            # find and connect to an ipcluster instance given the information
            # in the ipcluster dict. Connect to running one, or launch new. 
            ipyclient = self._get_parallel(ipyclient)

            # fill the input array with quartets to sample
            self.start = time.time()
            self._store_N_samples(ipyclient)

            # calculate invariants for the full seqarray 
            self.start = time.time()
            if os.path.exists(self.trees.tree):
                print("initial tree already inferred")
            else:
                Inference(self, ipyclient, self.start).run()
                    
            # calculate invariants for each bootstrap rep 
            self.start = time.time()
            if self.params.nboots:
                if self.checkpoint.boots == self.params.nboots:
                    print("{} bootstrap trees already inferred"
                          .format(self.checkpoint.boots))
                else:
                    while self.checkpoint.boots < self.params.nboots:
                        self.checkpoint.boots += 1
                        Inference(self, ipyclient, self.start).run()

            # write output stats
            TreeStats(self, ipyclient).run()

        # handle exceptions so they will be raised after we clean up below
        except KeyboardInterrupt as inst:
            print("\nKeyboard Interrupt by user. Cleaning up...")

        except IPyradError as inst:
            print("\nError encountered: {}\n[see trace below]\n".format(inst))
            raise 

        except Exception as inst:
            print("\nException encountered: {}\n[see trace below]\n".format(inst))
            raise

        # close client when done or interrupted
        finally:
            self._cleanup_parallel(ipyclient)


# class with functions run on the remote engines
class Inference:
    def __init__(self, tet, ipyclient, start=None, quiet=False):
        self.tet = tet
        self.quiet = quiet
        self.boot = bool(self.tet.checkpoint.boots)
        self.start = (start if start else time.time())
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.jobs = range(
            self.tet.checkpoint.arr, 
            self.tet.params.nquartets, 
            self.tet._chunksize)
        self.printstr = ("initial tree   ", "")
        if self.boot:
            self.printstr = ("bootstrap trees ", "")

    def run(self):
        self.get_snp_array()
        self.init_new_invariants_array()
        self.fill_database()
        self.write_tree()


    def init_new_invariants_array(self):
        "if this is a bootstrap then init a new boot array in the database"
        ## max val is 65535 in here if uint16
        bootkey = "boot{}".format(self.tet.checkpoint.boots)
        with h5py.File(self.tet.database.output, 'r+') as io5:
            if bootkey not in io5["invariants"].keys():
                io5["invariants"].create_dataset(
                    name=bootkey, 
                    shape=(self.tet.params.nquartets, 16, 16),
                    dtype=np.uint16,
                    chunks=(self.tet._chunksize, 16, 16))


    def get_snp_array(self):
        # resample the bootseq array
        if self.tet.files.mapfile:
            self.sample_bootseq_array_map()
        else:
            self.sample_bootseq_array()


    def insert_to_hdf5(self, chunk, results):
        #two result arrs
        chunksize = self.tet._chunksize
        qrts, invs = results

        ## enter into db
        with h5py.File(self.tet.database.output, 'r+') as io5:
            io5['quartets'][chunk:chunk + chunksize] = qrts

            ## entered as 0-indexed !
            if self.tet.params.save_invariants:
                if self.boot:
                    key = "invariants/boot{}".format(self.tet.checkpoint.boots)
                    io5[key][chunk:chunk + chunksize] = invs
                else:
                    io5["invariants/boot0"][chunk:chunk + chunksize] = invs


    def fill_database(self):
        # submit jobs distriuted across the cluster.
        asyncs = {}
        for job in self.jobs:
            asyncs[job] = self.lbview.apply(nworker, *(self.tet, job))

        # wait for jobs to finish, catch results as they return and
        # enter into HDF5 database and delete to keep memory low.
        done = 0
        while 1:
            ## gather finished jobs
            finished = [i for i, j in asyncs.items() if j.ready()]

            ## iterate over finished list
            for key in finished:
                rasync = asyncs[key]
                if rasync.successful():

                    # store result and purge it
                    done += 1
                    results = rasync.result()
                    self.insert_to_hdf5(key, results)
                    del asyncs[key]
                else:
                    raise IPyradError(rasync.result())

            # progress bar is different if first vs boot tree
            if not self.quiet:
                if not self.boot:
                    self.tet._progressbar(
                        len(self.jobs), 
                        done, 
                        self.start, 
                        self.printstr,
                    )
                else:
                    self.tet._progressbar(
                        self.tet.params.nboots, 
                        self.tet.checkpoint.boots, 
                        self.start,
                        self.printstr)

            ## done is counted on finish, so this means we're done
            if len(asyncs):
                time.sleep(0.1)
            else:
                break


    def write_tree(self):
        # dump quartets into a text file for QMC
        self.dump_qmc()

        # send to QMC
        if not self.boot:
            self.run_qmc(0)
            self.tet._print("")

        else:
            self.run_qmc(1)
            if self.tet.checkpoint.boots == self.tet.params.nboots:
                self.tet._print("")

        ## reset the checkpoint arr
        #self.tet.checkpoint.arr = 0


    def sample_bootseq_array(self):
        "Takes the seqarray and re-samples columns and saves to bootsarr."
        ## use 'r+' to read and write to existing array. This is super 
        ## similar to what is called in __init__. 
        with h5py.File(self.tet.database.input, 'r+') as io5:  
            ## load in the seqarr and maparr
            seqarr = io5["seqarr"][:]

            ## resample columns with replacement
            newarr = np.zeros(seqarr.shape, dtype=np.uint8)
            cols = np.random.randint(0, seqarr.shape[1], seqarr.shape[1])
            tmpseq = shuffle_cols(seqarr, newarr, cols)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self.tet.params.resolve_ambigs:
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
        """
        Re-samples loci with replacement to fill the bootarr sampling
        only a single SNP from each locus according to the maparr. 
        """
        with h5py.File(self.tet.database.input, 'r+') as io5:
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
            newbmap[:, 1] = np.arange(1, arrlen + 1)
            
            ## fill the new arrays            
            tmpseq, tmpmap = fill_boot(seqarr, newbarr, newbmap, spans, loci)

            ## resolve ambiguous bases randomly. We do this each time so that
            ## we get different resolutions.
            if self.tet.params.resolve_ambigs:
                tmpseq = resolve_ambigs(tmpseq)

            ## convert CATG bases to matrix indices
            tmpseq[tmpseq == 65] = 0
            tmpseq[tmpseq == 67] = 1
            tmpseq[tmpseq == 71] = 2
            tmpseq[tmpseq == 84] = 3

            ## store data sets
            io5.create_dataset("bootsmap", data=tmpmap)
            io5.create_dataset("bootsarr", data=tmpseq)


    def dump_qmc(self):
        """
        Writes the inferred quartet sets from the database to a text 
        file to be used as input for QMC. Quartets that had no information
        available (i.e., no SNPs) were written to the database as 0,0,0,0
        and are excluded here from the output.
        """
        ## open the h5 database
        with h5py.File(self.tet.database.output, 'r') as io5:

            ## create an output file for writing
            self.tet.files.qdump = os.path.join(
                self.tet.dirs, 
                self.tet.name + ".quartets.txt")
            with open(self.tet.files.qdump, 'w') as qdump:

                ## pull from db
                for idx in range(
                    0, 
                    self.tet.params.nquartets, 
                    self.tet._chunksize):
                    qchunk = io5["quartets"][idx:idx + self.tet._chunksize, :]
                    quarts = [tuple(j) for j in qchunk if np.any(j)]

                    ## shuffle and format for qmc
                    np.random.shuffle(quarts)
                    chunk = ["{},{}|{},{}".format(*i) for i in quarts]
                    qdump.write("\n".join(chunk) + "\n")


    def run_qmc(self, boot):
        "Runs quartet max-cut QMC on the quartets qdump file."
        # build command
        self._tmp = os.path.join(self.tet.dirs, ".tmptre")
        cmd = [
            ip.bins.qmc, 
            "qrtt=" + self.tet.files.qdump, 
            "otre=" + self._tmp,
        ]

        # run QMC on quartets input
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        res = proc.communicate()
        if proc.returncode:
            raise IPyradError(res[1])

        # parse tmp file written by QMC into a tree and rename tips
        ttre = toytree.tree(self._tmp)
        tips = ttre.treenode.get_leaves()
        for tip in tips:
            tip.name = self.tet.samples[int(tip.name)]
        newick = ttre.treenode.write(format=9)      

        # save the tree to file
        if boot:
            with open(self.tet.trees.boots, 'a') as outboot:
                outboot.write(newick + "\n")
        else:
            with open(self.tet.trees.tree, 'w') as outtree:
                outtree.write(newick)

        # save the new checkpoint and file paths
        self.tet._save()


class TreeStats:
    def __init__(self, tet, ipyclient):
        self.tet = tet
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.start = (self.tet.start if self.tet.start else time.time())
        self.samples = self.tet.samples


    def run(self):
        self.build_bootstrap_consensus()
        #self.build_nhx_stats()


    def build_bootstrap_consensus(self):
        "Compute sampling stats and consens trees"   
        # make a consensus from bootstrap reps.
        if self.tet.checkpoint.boots:
            boottrees = toytree.mtree(self.tet.trees.boots)
            self.ctre = boottrees.get_consensus_tree()
            self.ctre.write(self.tet.trees.cons)


    def build_nhx_stats(self):
        "Compute quartet sampling stats, and maybe others"
        ## build stats file and write trees
        qtots = {}
        qsamp = {}
        tots = len(self.ctre)
        totn = set(self.ctre.get_tip_labels())

        ## iterate over node traversal
        for node in self.ctre.treenode.traverse():
            # this is slow, needs to look at every sampled quartet
            # so we send it to be processed on engines
            qtots[node] = self.lbview.apply(get_total, *(tots, node))

            # TODO: error here on pickling...
            qsamp[node] = self.lbview.apply(get_sampled, *(self, totn, node))

        ## wait for jobs to finish (+1 to lenjob is for final progress printer)
        alljobs = list(qtots.values()) + list(qsamp.values())
        printstr = ("calculating stats", "")
        done = 0
        while 1:
            done = [i.ready() for i in alljobs]
            self.tet._progressbar(
                sum(done), len(done), self.tet.start, printstr)
            if all(done):
                break
            time.sleep(0.1)

        ## store results in the tree object
        for node in self.ctre.treenode.traverse():
            total = qtots[node].result()
            sampled = qsamp[node].result()
            node.add_feature("quartets_total", total)
            node.add_feature("quartets_sampled", sampled)
        features = ["quartets_total", "quartets_sampled"]

        ## update final progress
        self._progressbar(0, 0, 0, 0, True)

        ## write tree in NHX format 
        with open(self.tet.trees.nhx, 'w') as outtre:
            outtre.write(self.ctre.tree.write(format=0, features=features))        



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
    # TODO: is there a way to make this work for non-MKL (e.g., BLAS)?
    # or ideally to work more generally for both? Maybe just try/except, 
    # maybe OPM_NUMTHREADS?
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
        # TODO: new numba funcs supported, maybe this time...
        nmask = np.any(nall_mask[sidx], axis=0)
        nmask += np.all(seqs == seqs[0], axis=0) 

        ## here are the jitted funcs
        bidx, invar = calculate(seqs, maparr, nmask, TESTS)

        ## store results
        rquartets[idx] = smps[idx][bidx]
        rinvariants[idx] = invar

    # reset thread limit
    set_mkl_thread_limit(oldlimit)

    # return results...
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
                mats[0, (4 * i[0]) + i[1], (4 * i[2]) + i[3]] += 1      
                last_loc = mapcol[idx]

    ## fill the alternates
    x = np.uint8(0)
    for y in np.array([0, 4, 8, 12], dtype=np.uint8):
        for z in np.array([0, 4, 8, 12], dtype=np.uint8):
            mats[1, y:y + np.uint8(4), z:z + np.uint8(4)] = mats[0, x].reshape(4, 4)
            mats[2, y:y + np.uint8(4), z:z + np.uint8(4)] = mats[0, x].reshape(4, 4).T
            x += np.uint8(1)

    return mats


##################################################
## quartet sampling funcs to fill the database
##################################################

# replace self here with tet to be more clear
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


# not yet updated
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
        #rand = np.arange(0, n_choose_k(len(self.samples), 4))
        rand = np.arange(0, int(comb(len(self.samples), 4)))
        np.random.shuffle(rand)
        rslice = rand[:self.params.nquartets]
        rss = np.sort(rslice)
        riter = iter(rss)
        del rand, rslice

        ## print progress update 1 to the engine stdout
        print(self._chunksize)

        ## set to store
        rando = next(riter)
        tmpr = np.zeros((self.params.nquartets, 4), dtype=np.uint16)
        tidx = 0
        while 1:
            try:
                for i, j in enumerate(qiter):
                    if i == rando:
                        tmpr[tidx] = j
                        tidx += 1
                        rando = next(riter)

                    ## print progress bar update to engine stdout
                    if not i % self._chunksize:
                        print(min(i, self.params.nquartets))

            except StopIteration:
                break
        ## store into database
        fillsets[:] = tmpr
        del tmpr


# not yet updated for toytree or py3
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
        tre = toytree.etemini.TreeNode(self.files.tree)
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
            #total = n_choose_k(len(split[0]), 2) * n_choose_k(len(split[1]), 2)
            total = int(comb(len(split[0]), 2)) * int(comb(len(split[1]), 2))
            if total < squarts * 2:
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
            cycle, qiter = next(qitercycle)

            ## sample from iterators, store sorted set.
            try:
                qrtsamp = tuple(sorted(next(qiter)))
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

# def n_choose_k(n, k):
#     """ 
#     Get the number of quartets as n-choose-k. This is used in equal 
#     splits to decide whether a split should be exhaustively sampled
#     or randomly sampled. Edges near tips can be exhaustive while highly 
#     nested edges can do with less sampling.
#     """
#     mulfunc = lambda x, y: x * y
#     return int(reduce(mulfunc, (Fraction(n  -i, i+1) for i in range(k)), 1))


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
def resolve_ambigs(tmpseq):
    """ 
    Randomly resolve ambiguous bases. This is applied to each boot
    replicate so that over reps the random resolutions don't matter.
    Sites are randomly resolved, so best for unlinked SNPs since 
    otherwise linked SNPs are losing their linkage information... 
    though it's not like we're using it anyways.
    """

    ## the order of rows in GETCONS
    for aidx in range(6):
        #np.uint([82, 75, 83, 89, 87, 77]):
        ambig, res1, res2 = GETCONS[aidx]

        ## get true wherever tmpseq is ambig
        idx, idy = np.where(tmpseq == ambig)
        halfmask = np.random.choice(np.array([True, False]), idx.shape[0])

        for col in range(idx.shape[0]):
            if halfmask[col]:
                tmpseq[idx[col], idy[col]] = res1
            else:
                tmpseq[idx[col], idy[col]] = res2
    return tmpseq


@numba.jit(nopython=True)
def shuffle_cols(seqarr, newarr, cols):
    "Used in bootstrap resampling when no map file is present."
    for idx in range(cols.shape[0]):
        newarr[:, idx] = seqarr[:, cols[idx]]
    return newarr


@numba.jit(nopython=True)
def get_shape(spans, loci):
    """ 
    Returns shape of new bootstrap resampled locus array. The 
    shape can change because when we resample loci the number
    of SNPs in the data set can change.
    """
    width = 0
    for idx in range(loci.size):
        width += spans[loci[idx], 1] - spans[loci[idx], 0]
    return width


@numba.jit(nopython=True)
def fill_boot(seqarr, newboot, newmap, spans, loci):
    """ 
    Fills the new bootstrap SNP array and map array with
    new data based on the resampled loci for this boot.
    """
    ## column index
    cidx = 0
  
    ## resample each locus
    for i in range(loci.shape[0]):
        
        ## grab a random locus's columns
        x1 = spans[loci[i]][0]
        x2 = spans[loci[i]][1]
        cols = seqarr[:, x1:x2]

        ## randomize columns within colsq
        cord = np.random.choice(cols.shape[1], cols.shape[1], replace=False)
        rcols = cols[:, cord]
        
        ## fill bootarr with n columns from seqarr
        ## the required length was already measured
        newboot[:, cidx:cidx + cols.shape[1]] = rcols

        ## fill bootmap with new map info
        newmap[cidx: cidx + cols.shape[1], 0] = i + 1
        
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
            qrts = io5["quartets"][idx:idx + data._chunksize]
            for qrt in qrts:
                sqrt = set(qrt)
                if all([sqrt.intersection(i) for i in [lendr, lendl, lenur, lenul]]):
                    sampled += 1

            ## increase span
            idx += data._chunksize
    return sampled


    
###############################################
## Save/Load from JSON operations checkpoint.
###############################################

def _byteify(data, ignore_dicts=False):
    """
    converts unicode to utf-8 when reading in json files
    """
    if isinstance(data, str):
        return data.encode("utf-8")

    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]

    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.items()
        }
    return data
