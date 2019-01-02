#!/usr/bin/env python

""" demultiplex raw sequence data given a barcode map."""

# py2/3 compatible imports
from __future__ import print_function
from builtins import range
try:
    from itertools import izip, islice
except ImportError:
    from itertools import islice
    izip = zip

# external imports
import os
import io
import gzip
import glob
import time
import shutil
import pickle
import numpy as np
import subprocess as sps
from collections import Counter

# ipyrad imports
from ipyrad.core.sample import Sample
from ipyrad.assemble.utils import IPyradError, ambigcutters


class Step1:
    def __init__(self, data, force, ipyclient):
        # store attrs
        self.data = data
        self.force = force
        self.ipyclient = ipyclient

        # check input data files
        self.sfiles = self.data.params.sorted_fastq_path
        self.rfiles = self.data.params.raw_fastq_path
        self.select_method()
        self.setup_dirs()


    def run(self):
        if self.method == "link_fastqs":
            FileLinker(self).run()
        else:
            Demultiplexer(self).run()


    def setup_dirs(self):
        "create output directory, tmp directory."
        # assign fastq dir
        self.data.dirs.fastqs = os.path.join(
            self.data.params.project_dir,
            self.data.name + "_fastqs")       
        self.data.dirs.fastqs = os.path.realpath(self.data.dirs.fastqs)

        # remove existing if force flag
        if self.force:
            if os.path.exists(self.data.dirs.fastqs):
                shutil.rmtree(self.data.dirs.fastqs)
        # bail out if overwrite necessary but no force flag
        else:
            if os.path.exists(self.data.dirs.fastqs):
                raise IPyradError(
                    "Fastq dir {} already exists: use force flag to overwrite".
                    format(self.data.dirs.fastqs))

        # ensure project dir exists
        if not os.path.exists(self.data.params.project_dir):
            os.mkdir(self.data.params.project_dir)

        # ensure fastq dir exists
        if not os.path.exists(self.data.dirs.fastqs):
            os.mkdir(self.data.dirs.fastqs)


    def select_method(self):
        "Checks file paths and for existing samples and returns function"
        # do not allow both a sorted_fastq_path and a raw_fastq
        if self.sfiles and self.rfiles:
            raise IPyradError(NOT_TWO_PATHS)

        # but also require that at least one exists
        if not (self.sfiles or self.rfiles):
            raise IPyradError(NO_SEQ_PATH_FOUND)

        # print headers
        if self.data._cli:
            if self.sfiles:
                self.data._print(
                    "\n{}Step 1: Loading sorted fastq data to Samples"
                    .format(self._spacer))
            else:
                self.data._print(
                    "\n{}Step 1: Demultiplexing fastq data to Samples"
                    .format(self._spacer))

        # if Samples already exist then bail out unless force
        if self.data.samples:
            if not self.force:
                raise IPyradError(
                    SAMPLES_EXIST.format(
                        len(self.data.samples), self.data.name))
            else:
                if glob.glob(self.sfiles):
                    self.method = "link_fastqs"
                else:
                    self.method = "demultiplex"                   

        # Creating new Samples
        else:
            if glob.glob(self.sfiles):
                self.method = "link_fastqs"
            else:
                self.method = "demultiplex"


class FileLinker:
    def __init__(self, step):
        self.data = step.data
        self.input = step.sfiles
        self.fastqs = glob.glob(self.input)
        self.ftuples = []

        # parallel distributor
        self.ipyclient = step.ipyclient
        self.lbview = self.ipyclient.load_balanced_view()


    def run(self):
        # checks for bad names and fills self.ftuples with file handles
        self.check_files()
        # distributes jobs to parallel
        self.remote_run_linker()


    def check_files(self):
        # Assert files are not .bz2 format
        if any([i.endswith(".bz2") for i in self.fastqs]):
            raise IPyradError(NO_SUPPORT_FOR_BZ2.format(self.sfiles))

        # filter out any files without proper file endings. Raise if None
        endings = ("gz", "fastq", "fq")
        self.fastqs = [i for i in self.fastqs if i.split(".")[-1] in endings]
        if not self.fastqs:
            raise IPyradError(NO_FILES_FOUND_PAIRS
                .format(self.data.params.sorted_fastq_path))

        # link pairs into tuples
        if 'pair' in self.data.params.datatype:
            # check that names fit the paired naming convention
            # trying to support flexible types (_R2_, _2.fastq)
            r1_try1 = [i for i in self.fastqs if "_R1_" in i]
            r1_try2 = [i for i in self.fastqs if i.endswith("_1.fastq.gz")]
            r1_try3 = [i for i in self.fastqs if i.endswith("_R1.fastq.gz")]

            r2_try1 = [i for i in self.fastqs if "_R2_" in i]
            r2_try2 = [i for i in self.fastqs if i.endswith("_2.fastq.gz")]
            r2_try3 = [i for i in self.fastqs if i.endswith("_R2.fastq.gz")]

            r1s = [r1_try1, r1_try2, r1_try3]
            r2s = [r2_try1, r2_try2, r2_try3]

            # check that something was found
            if not r1_try1 + r1_try2 + r1_try3:
                raise IPyradError(
                    "Paired filenames are improperly formatted. See Docs.")

            # find the one with the right number of R1s
            for idx, tri in enumerate(r1s):
                if len(tri) == len(self.fastqs) / 2:
                    break
            r1_files = r1s[idx]
            r2_files = r2s[idx]

            if len(r1_files) != len(r2_files):
                raise IPyradError(R1_R2_name_error.format(
                    len(r1_files), len(r2_files)))
            self.ftuples = [(i, j) for i, j in zip(r1_files, r2_files)]

        # data are not paired, create empty tuple pair
        else:
            # print warning if _R2_ is in names when not paired
            idx = 0
            if any(["_R2_" in i for i in self.fastqs]):
                print(NAMES_LOOK_PAIRED_WARNING)
            self.ftuples = [(i, "") for i in self.fastqs]


    def remote_run_linker(self):
        "read in fastq files and count nreads for stats and chunking in s2."

        # local counters 
        createdinc = 0

        # iterate over input files
        for ftup in self.ftuples:
            
            # remove file extension from name
            sname = get_name_from_file(ftup[0], None, None)

            # Create new Sample Class objects with names from files
            if sname not in self.data.samples:
                newsamp = Sample(sname)
                newsamp.stats.state = 1
                newsamp.barcode = None
                newsamp.files.fastqs = [ftup]
                self.data.samples[sname] = newsamp
                createdinc += 1

        # send jobs to engines for counting with cat/zcat | wc
        rasyncs = {}
        if createdinc:
            for sample in self.data.samples.values():
                gzipped = bool(sample.files.fastqs[0][0].endswith(".gz"))
                rasyncs[sample.name] = self.lbview.apply(
                    zbufcountlines, 
                    *(sample.files.fastqs[0][0], gzipped)
                )

        # wait for link jobs to finish if parallel
        start = time.time()
        printstr = ("loading reads       ", "s1")
        while 1:
            fin = [i.ready() for i in rasyncs.values()]
            self.data._progressbar(len(fin), sum(fin), start, printstr)
            time.sleep(0.1)
            if len(fin) == sum(fin):
                self.data._print("")
                break

        # collect link job results           
        for sname in rasyncs:
            rasync = rasyncs[sname]
            if rasync.successful():
                res = rasyncs[sname].get() / 4
                self.data.samples[sname].stats.reads_raw = res
                self.data.samples[sname].stats_dfs.s1["reads_raw"] = res
                self.data.samples[sname].state = 1
            else:
                raise IPyradError(rasync.get())

        # print if data were linked
        if createdinc:
            # double for paired data
            if 'pair' in self.data.params.datatype:
                createdinc = createdinc * 2
            if self.data._cli:
                self.data._print("{}{} fastq files loaded to {} Samples."
                    .format(self.data._spacer, createdinc, len(self.samples)))

        # save step-1 stats. We don't want to write this to the fastq dir, b/c
        # it is not necessarily inside our project dir. Instead, we'll write 
        # this file into our project dir in the case of linked_fastqs.
        self.data.stats_dfs.s1 = self.data._build_stat("s1")
        self.data.stats_files.s1 = os.path.join(
            self.data.params.project_dir,
            self.data.name + '_s1_demultiplex_stats.txt')
        with open(self.data.stats_files.s1, 'w') as outfile:
            (self.data.stats_dfs.s1
                .fillna(value=0)
                .astype(np.int)
                .to_string(outfile))


class Demultiplexer:
    def __init__(self, step):
        self.data = step.data
        self.input = step.rfiles
        self.fastqs = glob.glob(self.input)
        self.ipyclient = step.ipyclient
        # single engine jobs
        self.iview = self.ipyclient.load_balanced_view(targets=[0])

        # limited multi-engine jobs
        if len(self.ipyclient.ids) >= 12:
            targets = self.ipyclient.ids[::4]
        else:
            targets = self.ipyclient.ids[:4]
        self.lbview = self.ipyclient.load_balanced_view(targets=targets)

        # we better have barcodes...
        if not self.data.barcodes:
            self.data._link_barcodes()

        # attrs filled by check_files
        self.ftuples = []
        self.chunksdict = {}
        self.longbar = None
        self.check_files()

        # attrs filled by get_barcode_dict
        self.cutters = None
        self.matchdict = {}       
        self.get_barcode_dict()        

        # store stats for each file handle (grouped results of chunks)
        self.stats = Stats()

    def run(self):
        # Estimate size of files to plan parallelization. 
        self.setup_for_splitting()

        # work load; i.e., is there one giant file or many small files?
        self.splitfiles()

        # process the files or chunked file bits        
        self.remote_run_barmatch()

        # concatenate chunks
        self.concatenate_chunks()

        # store stats and create Sample objects in Assembly
        self.store_stats()

    # init functions -----------------------------------
    def check_files(self):
        "Check that data files are present and formatted correctly"

        # check for data using glob for fuzzy matching
        if not self.fastqs:
            raise IPyradError(
                NO_RAWS.format(self.data.params.raw_fastq_path))

        # find longest barcode
        if not os.path.exists(self.data.params.barcodes_path):
            raise IPyradError(
                "Barcodes file not found. You entered: '{}'".format(
                    self.data.params.barcodes_path))

        # Handle 3rad multi-barcodes. Gets len of the first one. 
        blens = [len(i.split("+")[0]) for i in self.data.barcodes.values()]
        if len(set(blens)) == 1:
            self.longbar = (blens[0], 'same')
        else:
            self.longbar = (max(blens), 'diff')

        # For 3rad we need to add the length info for barcodes_R2
        if "3rad" in self.data.params.datatype:
            blens = [
                len(i.split("+")[1]) for i in self.data.barcodes.values()]
            self.longbar = (self.longbar[0], self.longbar[1], max(blens))

        # gather raw sequence filenames (people want this to be flexible ...)
        if 'pair' in self.data.params.datatype:
            firsts = [i for i in self.fastqs if "_R1_" in i]
            if not firsts:
                raise IPyradError(
                    "First read files names must contain '_R1_'.")
            seconds = [i.replace("_R1_", "_R2_") for i in firsts]
            self.ftuples = list(zip(firsts, seconds))
        else:
            self.ftuples = list(zip(self.fastqs, iter(int, 1)))


    def get_barcode_dict(self):
        # returns a list of both resolutions of cut site 1
        # (TGCAG, ) ==> [TGCAG, ]
        # (TWGC, ) ==> [TAGC, TTGC]
        # (TWGC, AATT) ==> [TAGC, TTGC]
        self.cutters = [
            ambigcutters(i) for i in self.data.params.restriction_overhang
        ]
        assert self.cutters, "Must enter a restriction_overhang for demultiplexing."

        # get matchdict
        self.matchdict = inverse_barcodes(self.data)


    def setup_for_splitting(self):
        # create a tmpdir for chunked_files and a chunk optimizer 
        self.tmpdir = os.path.realpath(
            os.path.join(self.data.dirs.fastqs, "tmpdir")
        )
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)

        # chunk into 16 pieces
        self.nreads = estimate_nreads(self.data, self.ftuples[0][0])
        self.optim = int(self.nreads / 16)

        # if more files than cpus or optim<8M: no chunking
        self.do_file_split = 0
        omin = int(8e6)
        if (len(self.ftuples) > len(self.ipyclient)) or (self.optim > omin):
            self.do_file_split = 1


    def splitfiles(self):
        "sends raws to be chunked"

        # send slices N at a time. The dict chunkfiles stores a tuple of 
        # rawpairs dictionary to store asyncresults for sorting jobs
        printstr = ('chunking large files', 's1')
        start = time.time()
        done = 0
        njobs = (
            32 * len(self.ftuples) if "pair" in self.data.params.datatype
            else 16 * len(self.ftuples)
        )
        rasyncs = {}
        chunksdict = {}
        for fidx, ftup in enumerate(self.ftuples):

            # get file handle w/o basename for stats output
            handle = os.path.splitext(os.path.basename(ftup[0]))[0]
            if not self.do_file_split:
                chunksdict[handle] = [ftup]

            # chunk file into 4 bits using zcat_make_temps                
            else:               
                args = (self.data, ftup, fidx, self.tmpdir, self.optim, start)
                rasyncs[handle] = self.iview.apply(zcat_make_temps, *args)

        # track progress until finished
        if rasyncs:
            while 1:
                done = len(glob.glob(os.path.join(self.tmpdir, "chunk*_*_*")))
                self.data._progressbar(
                    njobs, min(njobs, done - 1), start, printstr)
                time.sleep(0.5)
                if all([i.ready() for i in rasyncs.values()]):
                    break

            # store results
            for key, val in rasyncs.items():
                chunksdict[key] = val.get()       

            # clean up                    
            self.ipyclient.purge_everything()                    
            self.data._progressbar(njobs, njobs, start, printstr)
            self.data._print("")

        # return value
        self.chunksdict = chunksdict
        

    def remote_run_barmatch(self):
        "Submit chunks to be sorted by barmatch() and collect stats"
        # progress bar info
        start = time.time()
        printstr = ("sorting reads       ", "s1")

        # chunkfiles is a dict with {handle: chunkslist, ...}. The func barmatch
        # writes results to samplename files with PID number, and also writes a 
        # pickle for chunk specific results with fidx suffix, which it returns.
        rasyncs = {}
        ridx = 0
        for handle, ftuplist in self.chunksdict.items():
            for fidx, ftuple in enumerate(ftuplist):
                args = (
                    self.data, ftuple, self.longbar, 
                    self.cutters, self.matchdict, fidx,
                    )
                rasync = self.lbview.apply(barmatch, args)
                rasyncs[ridx] = (handle, rasync)
                ridx += 1

            # get ready to receive stats: 'total', 'cutfound', 'matched'
            self.stats.perfile[handle] = np.zeros(3, dtype=np.int)

        # collect and store results as jobs finish
        njobs = len(rasyncs)
        done = 0
        while 1:
            # get list of ridx numbers for finished jobs
            finished = [i for (i, j) in rasyncs.items() if j[1].ready()]

            # cleanup finished ridx jobs and grab stats
            for ridx in finished:
                handle, rasync = rasyncs[ridx]
                if rasync.successful():
                    pkl = rasync.get()
                    self.stats.fill_from_pickle(pkl, handle)
                    del rasyncs[ridx]
                    done += 1
                else:
                    raise IPyradError(rasync.get())

            # print progress
            self.data._progressbar(njobs, done, start, printstr)
            time.sleep(0.1)
            if njobs == done:
                self.data._print("")
                break


    def concatenate_chunks(self):
        """ 
        If multiple chunk files match to the same sample name but with 
        different barcodes (i.e., they are technical replicates) then this
        will assign all the files to the same sample name file.
        """
        # collate files progress bar
        start = time.time()
        printstr = ("writing/compressing ", "s1")
        self.data._progressbar(10, 0, start, printstr) 

        # get all the files
        ftmps = glob.glob(os.path.join(
            self.data.dirs.fastqs, 
            "tmpdir", 
            "tmp_*.fastq"))

        # a dict to assign tmp files to names/reads
        r1dict = {}
        r2dict = {}
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            r1dict[sname] = []
            r2dict[sname] = []

        # assign to name keys
        for ftmp in ftmps:
            base, orient, _ = ftmp.rsplit("_", 2)
            sname = base.rsplit("/", 1)[-1].split("tmp_", 1)[1]
            if orient == "R1":
                r1dict[sname].append(ftmp)
            else:
                r2dict[sname].append(ftmp)

        ## concatenate files
        snames = []
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            snames.append(sname)

        writers = []
        for sname in set(snames):
            tmp1s = sorted(r1dict[sname])
            tmp2s = sorted(r2dict[sname])
            writers.append(
                self.iview.apply(
                    collate_files, 
                    *(self.data, sname, tmp1s, tmp2s))
                )

        total = len(writers)
        while 1:
            ready = [i.ready() for i in writers]
            self.data._progressbar(total, sum(ready), start, printstr)
            time.sleep(0.1)
            if all(ready):
                self.data._print("")
                break
        

    def store_stats(self):
        "Write stats and stores to Assembly object."

        # out file
        self.data.stats_files.s1 = os.path.join(
            self.data.dirs.fastqs, 's1_demultiplex_stats.txt')
        outfile = open(self.data.stats_files.s1, 'w')

        # write the header for file stats ------------------------------------
        outfile.write("{:<35}  {:>13}{:>13}{:>13}\n"
            .format("raw_file", "total_reads", "cut_found", "bar_matched"))

        # write the file stats
        r1names = sorted(self.stats.perfile)
        for fname in r1names:
            dat = self.stats.perfile[fname]
            outfile.write(
                "{:<35}  {:>13}{:>13}{:>13}\n"
                .format(fname, dat[0], dat[1], dat[2])
            )
            # repeat for pairfile
            if 'pair' in self.data.params.datatype:
                fname = fname.replace("_R1_", "_R2_")
                outfile.write(
                    "{:<35}  {:>13}{:>13}{:>13}\n"
                    .format(fname, dat[0], dat[1], dat[2])
                )

        # spacer, how many records for each sample --------------------------
        outfile.write("\n{:<35}  {:>13}\n"
            .format("sample_name", "total_reads"))

        # names alphabetical. Write to file. Will save again below to Samples.
        snames = set()
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            snames.add(sname)
            
        for sname in sorted(list(snames)):
            outfile.write("{:<35}  {:>13}\n"
                .format(sname, self.stats.fsamplehits[sname]))

        ## spacer, which barcodes were found -----------------------------------
        outfile.write('\n{:<35}  {:>13} {:>13} {:>13}\n'
            .format("sample_name", "true_bar", "obs_bar", "N_records"))

        ## write sample results
        for sname in sorted(self.data.barcodes):
            if "-technical-replicate-" in sname:
                fname = sname.rsplit("-technical-replicate", 1)[0]  
            else:
                fname = sname
                
            # write perfect hit
            hit = self.data.barcodes[sname]
            offhitstring = ""
        
            # write off-n hits
            # sort list of off-n hits  
            if fname in self.stats.fdbars:
                offkeys = list(self.stats.fdbars.get(fname))
                for offhit in offkeys[::-1]:
                    # exclude perfect hit
                    if offhit not in self.data.barcodes.values():
                        offhitstring += (
                            "{:<35}  {:>13} {:>13} {:>13}\n"
                            .format(sname, hit, offhit, 
                                self.stats.fbarhits[offhit] / 2)
                            )
                        #sumoffhits += fbarhits[offhit]
            
                # write string to file
                outfile.write("{:<35}  {:>13} {:>13} {:>13}\n"
                    .format(sname, hit, hit, 
                        self.stats.fbarhits[hit] / 2))
                outfile.write(offhitstring)
            
        # write misses
        misskeys = list(self.stats.fmisses.keys())
        misskeys.sort(key=self.stats.fmisses.get)
        for key in misskeys[::-1]:
            outfile.write('{:<35}  {:>13}{:>13}{:>13}\n'
                .format("no_match", "_", key, self.stats.fmisses[key]))
        outfile.close()        

        # Link Sample with this data file to the Assembly object
        for sname in snames:

            # make the sample
            sample = Sample(sname)

            # allow multiple barcodes if its a replicate. 
            barcodes = []
            for n in range(500):
                fname = sname + "-technical-replicate-{}".format(n)
                fbar = self.data.barcodes.get(fname)
                if fbar:
                    barcodes.append(fbar)
            if barcodes:
                sample.barcode = barcodes
            else:
                sample.barcode = self.data.barcodes[sname]

            # file names        
            if 'pair' in self.data.params.datatype:
                sample.files.fastqs = [(
                    os.path.join(
                        self.data.dirs.fastqs, sname + "_R1_.fastq.gz"),
                    os.path.join(
                        self.data.dirs.fastqs, sname + "_R2_.fastq.gz"), 
                    )]
            else:
                sample.files.fastqs = [
                    (os.path.join(
                        self.data.dirs.fastqs, 
                        sname + "_R1_.fastq.gz",
                        ),
                    ""),
                ]

            # fill in the summary stats
            sample.stats["reads_raw"] = int(self.stats.fsamplehits[sname])
            # fill in the full df stats value
            sample.stats_dfs.s1["reads_raw"] = int(
                self.stats.fsamplehits[sname])

            # Only link Sample if it has data
            if sample.stats["reads_raw"]:
                sample.stats.state = 1
                self.data.samples[sample.name] = sample
            else:
                print("Excluded sample: no data found for", sname)

        # initiate s1 key for data object
        self.data.stats_dfs.s1 = self.data._build_stat("s1")

        # cleanup
        shutil.rmtree(self.tmpdir)


# this class is created and run inside the barmatch function() that is run
# on remote engines for parallelization.
class BarMatch:
    def __init__(self, data, ftuple, longbar, cutters, matchdict, fidx):
        """
        Sorts reads to samples based on barcodes and writes stats to a pickle.
        """
        # store attrs
        self.data = data
        self.longbar = longbar
        self.cutters = cutters
        self.ftuple = ftuple
        self.matchdict = matchdict
        self.fidx = fidx

        # when to write to disk
        self.chunksize = int(1e6) 
        self.epid = os.getpid()
        self.filestat = np.zeros(3, dtype=int)
        
        # store reads per sample (group technical replicates)
        self.samplehits = {}
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            self.samplehits[sname] = 0

        # store all barcodes observed
        self.barhits = {}
        for barc in self.matchdict:
            self.barhits[barc] = 0

        # store reads and bars matched to samples
        self.read1s = {} 
        self.read2s = {} 
        self.dbars = {} 
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            self.read1s[sname] = []
            self.read2s[sname] = []
            self.dbars[sname] = set()

        # store counts of what didn't match to samples
        self.misses = {}
        self.misses['_'] = 0


    def run(self):
        self.demux = self.get_matching_function()
        self.open_read_generators()
        pkl = self.sort_reads()
        self.close_read_generators()
        return pkl


    def get_matching_function(self):
        if self.longbar[1] == 'same':
            if self.data.params.datatype == '2brad':
                return getbarcode1
            else:
                return getbarcode2
        else:
            return getbarcode3


    def open_read_generators(self):

        # get file type
        if self.ftuple[0].endswith(".gz"):
            self.ofile1 = gzip.open(self.ftuple[0], 'r')
        else:
            self.ofile1 = open(self.ftuple[0], 'r')

        # create iterators 
        fr1 = iter(self.ofile1) 
        quart1 = izip(fr1, fr1, fr1, fr1)
        
        # create second read iterator for paired data
        if self.ftuple[1]:
            if self.ftuple[0].endswith(".gz"):
                self.ofile2 = gzip.open(self.ftuple[1], 'r')
            else:
                self.ofile2 = open(self.ftuple[1], 'r')

            # create iterators
            fr2 = iter(self.ofile2)  
            quart2 = izip(fr2, fr2, fr2, fr2)
            self.quarts = izip(quart1, quart2)
        else:
            self.quarts = izip(quart1, iter(int, 1))


    def close_read_generators(self):
        self.ofile1.close()
        if self.ftuple[1]:
            self.ofile2.close()


    def sort_reads(self):
        while 1:
            # read in four lines of data and increase counter
            try:
                read1, read2 = next(self.quarts)
                read1 = list(read1)
                self.filestat[0] += 1
            except StopIteration:
                break
        
            # use barcode finding function to separate read from barcode
            if '3rad' not in self.data.params.datatype:
                # Parse barcode. Uses the parsing function selected above.
                barcode = self.demux(self.cutters, read1, self.longbar)
            else:
                # Here we're just reusing the findbcode function for R2
                # and reconfiguring the longbar tuple to have the maxlen for
                # the R2 barcode.
                barcode1 = find3radbcode(self.cutters, self.longbar, read1)
                barcode2 = find3radbcode(self.cutters, self.longbar, read2)
                barcode = barcode1 + "+" + barcode2

            # ensure barcode is string
            try:
                barcode = barcode.decode()
            except AttributeError:
                pass          
       
            # find if it matches 
            sname_match = self.matchdict.get(barcode)

            if sname_match:

                # add to observed set of bars
                self.dbars[sname_match].add(barcode)
                self.filestat[1:] += 1

                self.samplehits[sname_match] += 1
                self.barhits[barcode] += 1
                if barcode in self.barhits:
                    self.barhits[barcode] += 1
                else:
                    self.barhits[barcode] = 1
        
                # trim off barcode
                lenbar = len(barcode)
                if '3rad' in self.data.params.datatype:
                    ## Iff 3rad trim the len of the first barcode
                    lenbar = len(barcode1)
        
                # for 2brad we trim the barcode AND the synthetic overhang
                # The `+1` is because it trims the newline
                if self.data.params.datatype == '2brad':
                    overlen = len(self.cutters[0][0]) + lenbar + 1
                    read1[1] = read1[1][:-overlen] + "\n"
                    read1[3] = read1[3][:-overlen] + "\n"
                else:
                    read1[1] = read1[1][lenbar:]
                    read1[3] = read1[3][lenbar:]
        
                # Trim barcode off R2 and append. Only 3rad datatype
                # pays the cpu cost of splitting R2
                if '3rad' in self.data.params.datatype:
                    read2 = list(read2)
                    read2[1] = read2[1][len(barcode2):]
                    read2[3] = read2[3][len(barcode2):]
        
                # append to sorted reads list
                self.read1s[sname_match].append(b"".join(read1).decode())
                if 'pair' in self.data.params.datatype:
                    self.read2s[sname_match].append(b"".join(read2).decode()) 

            else:
                self.misses["_"] += 1
                if barcode:
                    self.filestat[1] += 1

            # Write to each sample file (pid's have different handles)
            if not self.filestat[0] % int(1e6):
                
                # write reads to file
                writetofile(self.data, self.read1s, 1, self.epid)
                if 'pair' in self.data.params.datatype:
                    writetofile(self.data, self.read2s, 2, self.epid)
                
                # clear out lits of sorted reads
                for sname in self.data.barcodes:
                    if "-technical-replicate-" in sname:
                        sname = sname.rsplit("-technical-replicate", 1)[0]
                    self.read1s[sname] = []
                    self.read2s[sname] = []             

        ## write the remaining reads to file
        writetofile(self.data, self.read1s, 1, self.epid)
        if 'pair' in self.data.params.datatype:
            writetofile(self.data, self.read2s, 2, self.epid)

        ## return stats in saved pickle b/c return_queue is too small
        ## and the size of the match dictionary can become quite large
        samplestats = [self.samplehits, self.barhits, self.misses, self.dbars]
        pklname = os.path.join(
            self.data.dirs.fastqs, 
            "tmpdir",
            "tmp_{}_{}.p".format(self.epid, self.fidx))
        with open(pklname, 'wb') as wout:
            pickle.dump([self.filestat, samplestats], wout)
        return pklname


# used inside BarMatch to store stats nicely.
class Stats:
    def __init__(self):
        # stats for each raw input file
        self.perfile = {}

        # stats for each sample
        self.fdbars = {}
        self.fsamplehits = Counter()
        self.fbarhits = Counter()
        self.fmisses = Counter()


    def fill_from_pickle(self, pkl, handle):

        # load in stats pickle
        with open(pkl, 'rb') as infile:
            filestats, samplestats = pickle.load(infile)

        ## pull new stats
        self.perfile[handle] += filestats

        ## update sample stats
        samplehits, barhits, misses, dbars = samplestats
        self.fsamplehits.update(samplehits)
        self.fbarhits.update(barhits)
        self.fmisses.update(misses)
        self.fdbars.update(dbars)


# -------------------------------------
# EXTERNAL FUNCS 
# -------------------------------------

def barmatch(args):
    # run procesor
    bar = BarMatch(*args)
    # writes reads t ofile and writes stats to pickle
    pkl = bar.run()
    return pkl


def get_name_from_file(fname, splitnames, fields):
    "Grab Sample names from demultiplexed input fastq file names"

    # allowed extensions
    file_extensions = [".gz", ".fastq", ".fq", ".fasta", ".clustS", ".consens"]
    base, _ = os.path.splitext(os.path.basename(fname))

    # remove read number from name
    base = base.replace("_R1_.", ".")\
               .replace("_R1_", "")\
               .replace("_R1.", ".")

    # remove extensions, retains '.' in file names.
    while 1:
        tmpb, tmpext = os.path.splitext(base)
        if tmpext in file_extensions:        
            base = tmpb
        else:
            break

    # split on 'splitnames' and select field 'fields'
    if fields:
        namebits = base.split(splitnames)
        base = []
        for field in fields:
            try:
                base.append(namebits[field])
            except IndexError:
                pass
        base = splitnames.join(base)

    # don't allow empty names
    if not base:
        raise IPyradError("""
    Found invalid/empty filename in link_fastqs. Check splitnames argument.
    """)
    return base


def zbufcountlines(filename, gzipped):
    "Fast line counter using unix utils"
    if gzipped:
        cmd1 = ["gunzip", "-c", filename]
    else:
        cmd1 = ["cat", filename]
    cmd2 = ["wc"]
    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, stderr=sps.PIPE)
    proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, stderr=sps.PIPE)
    res = proc2.communicate()[0]
    if proc2.returncode:
        raise IPyradError("error zbufcountlines {}:".format(res))
    nlines = int(res.split()[0])
    return nlines


def find3radbcode(cutters, longbar, read):
    "find barcode sequence in the beginning of read"
    ## default barcode string
    for ambigcuts in cutters:
        for cutter in ambigcuts:
            ## If the cutter is unambiguous there will only be one.
            if not cutter:
                continue
            search = read[1][:int(longbar[0] + len(cutter) + 1)]
            splitsearch = search.decode().rsplit(cutter, 1)
            if len(splitsearch) > 1:
                return splitsearch[0]
    ## No cutter found
    return splitsearch[0] 


def getbarcode1(cutters, read1, longbar):
    "find barcode for 2bRAD data"
    #+1 is for the \n at the end of the sequence line
    lencut = len(cutters[0][0]) + 1
    return read1[1][:-lencut][-longbar[0]:]

def getbarcode2(_, read1, longbar):
    "finds barcode for invariable length barcode data"
    return read1[1][:longbar[0]]

def getbarcode3(cutters, read1, longbar):
    "find barcode sequence in the beginning of read"
    ## default barcode string
    for cutter in cutters[0]:
        ## If the cutter is unambiguous there will only be one.
        if not cutter:
            continue
        search = read1[1][:int(longbar[0] + len(cutter) + 1)]
        barcode = search.rsplit(cutter, 1)
        if len(barcode) > 1:
            return barcode[0]
    ## No cutter found
    return barcode[0] 



def writetofile(data, dsort, read, pid):
    "Writes sorted data to tmp files"
    if read == 1:
        rrr = "R1"
    else:
        rrr = "R2"

    # appends to file for each sample, avoids parallel fighting by using 
    # pid assigned file handle.
    for sname in dsort:
        handle = os.path.join(
            data.dirs.fastqs, 
            "tmpdir",
            "tmp_{}_{}_{}.fastq".format(sname, rrr, pid))
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))
            # b"".join(dsort[sname]).decode())



def collate_files(data, sname, tmp1s, tmp2s):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    ## out handle
    out1 = os.path.join(data.dirs.fastqs, "{}_R1_.fastq.gz".format(sname))
    out = io.BufferedWriter(gzip.open(out1, 'w'))

    ## build cmd
    cmd1 = ['cat']
    for tmpfile in tmp1s:
        cmd1 += [tmpfile]

    ## compression function
    proc = sps.Popen(['which', 'pigz'], 
        stderr=sps.PIPE, 
        stdout=sps.PIPE).communicate()
    if proc[0].strip():
        compress = ["pigz"]
    else:
        compress = ["gzip"]

    ## call cmd
    proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
    proc2 = sps.Popen(compress, 
        stdin=proc1.stdout, 
        stderr=sps.PIPE, 
        stdout=out)
    err = proc2.communicate()
    if proc2.returncode:
        raise IPyradError("error in collate_files R1 %s", err)
    proc1.stdout.close()
    out.close()

    ## then cleanup
    for tmpfile in tmp1s:
        os.remove(tmpfile)

    if 'pair' in data.params.datatype:
        ## out handle
        out2 = os.path.join(data.dirs.fastqs, "{}_R2_.fastq.gz".format(sname))
        out = io.BufferedWriter(gzip.open(out2, 'w'))

        ## build cmd
        cmd1 = ['cat']
        for tmpfile in tmp2s:
            cmd1 += [tmpfile]

        ## call cmd
        proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
        proc2 = sps.Popen(compress, 
            stdin=proc1.stdout,
            stderr=sps.PIPE, 
            stdout=out)
        err = proc2.communicate()
        if proc2.returncode:
            raise IPyradError("error in collate_files R2 %s", err)
        proc1.stdout.close()
        out.close()

        ## then cleanup
        for tmpfile in tmp2s:
            os.remove(tmpfile)




def inverse_barcodes(data):
    """ Build full inverse barcodes dictionary """
    matchdict = {}
    bases = set("CATGN")
    poss = set()

    ## do perfect matches
    for sname, barc in data.barcodes.items():
        ## remove -technical-replicate-N if present
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        matchdict[barc] = sname
        poss.add(barc)

        if data.params.max_barcode_mismatch:
            ## get 1-base diffs
            for idx1, base in enumerate(barc):
                diffs = bases.difference(base)
                for diff in diffs:
                    lbar = list(barc)
                    lbar[idx1] = diff
                    tbar1 = "".join(lbar)
                    if tbar1 not in poss:
                        matchdict[tbar1] = sname                    
                        poss.add(tbar1)
                    else:
                        if matchdict.get(tbar1) != sname:
                            print("""\
        Note: barcodes {}:{} and {}:{} are within {} base change of each other
            Ambiguous barcodes that match to both samples will arbitrarily
            be assigned to the first sample. If you do not like this idea 
            then lower the value of max_barcode_mismatch and rerun step 1\n"""
            .format(sname, barc, 
                matchdict[tbar1], data.barcodes[matchdict[tbar1]],
                data.params.max_barcode_mismatch))

                ## if allowing two base difference things get big
                ## for each modified bar, allow one modification to other bases
                if data.params.max_barcode_mismatch > 1:
                    for idx2, _ in enumerate(tbar1):
                        ## skip the base that is already modified
                        if idx2 != idx1:
                            for diff in bases.difference(tbar1[idx2]):
                                ltbar = list(tbar1)
                                ltbar[idx2] = diff
                                tbar2 = "".join(ltbar)
                                if tbar2 not in poss:
                                    matchdict[tbar2] = sname                    
                                    poss.add(tbar2)
                                else:
                                    if matchdict.get(tbar2) != sname:
                                        print("""\
        Note: barcodes {}:{} and {}:{} are within {} base change of each other\
             Ambiguous barcodes that match to both samples will arbitrarily
             be assigned to the first sample. If you do not like this idea 
             then lower the value of max_barcode_mismatch and rerun step 1\n"""
             .format(sname, barc, 
                     matchdict[tbar2], data.barcodes[matchdict[tbar2]],
                     data.params.max_barcode_mismatch))
    return matchdict



def estimate_nreads(data, testfile):
    """ 
    Estimate a reasonable optim value by grabbing a chunk of sequences, 
    decompressing and counting them, to estimate the full file size.
    """
    ## count the len of one file and assume all others are similar len
    insize = os.path.getsize(testfile)
    tmp_file_name = os.path.join(
        data.params.project_dir, 
        "tmp-step1-count.fq")

    if testfile.endswith(".gz"):
        infile = gzip.open(testfile)
        outfile = gzip.open(tmp_file_name, 'w', compresslevel=5)
    else:
        infile = open(testfile)
        outfile = open(tmp_file_name, 'w')
        
    ## We'll take the average of the size of a file based on the
    ## first 10000 reads to approximate number of reads in the main file
    dat = b"".join(islice(infile, 40000))
    outfile.write(dat)
    outfile.close()
    infile.close()

    ## Get the size of the tmp file
    tmp_size = os.path.getsize(tmp_file_name)

    # divide by the tmp file size and multiply by 10000 to approximate
    # the size of the input .fq files
    inputreads = int(insize / tmp_size) * 10000
    os.remove(tmp_file_name)

    return inputreads


# def _cleanup_and_die(data):
#     """ cleanup func for step 1 """
#     tmpfiles = glob.glob(os.path.join(data.dirs.fastqs, "tmp_*_R*.fastq"))
#     tmpfiles += glob.glob(os.path.join(data.dirs.fastqs, "tmp_*.p"))
#     for tmpf in tmpfiles:            
#         os.remove(tmpf)




# def demux2(data, chunkfiles, cutters, longbar, matchdict, ipyclient):
#     """ 
#     Submit chunks to be sorted by the barmatch() function then 
#     calls putstats().
#     """

#     ## parallel stuff, limit to 1/4 of available cores for RAM limits.
#     start = time.time()
#     printstr = ("sorting reads       ", "s1")
#     lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::4])

#     ## store statcounters and async results in dicts
#     perfile = {}
#     filesort = {}
#     total = 0
#     done = 0 

#     ## chunkfiles is a dict with {handle: chunkslist, ...}. The func barmatch
#     ## writes results to samplename files with PID number, and also writes a 
#     ## pickle for chunk specific results with fidx suffix, which it returns.
#     for handle, rawtuplist in chunkfiles.items():
#         for fidx, rawtuple in enumerate(rawtuplist):
#             # get args for job
#             args = (data, rawtuple, cutters, longbar, matchdict, fidx)

#             # submit the job --------------------
#             rasync = lbview.apply(barmatch, *args)
#             filesort[total] = (handle, rasync)
#             total += 1

#             # get ready to receive stats: 'total', 'cutfound', 'matched'
#             perfile[handle] = np.zeros(3, dtype=np.int)

#     ## stats for each sample
#     fdbars = {}
#     fsamplehits = Counter()
#     fbarhits = Counter()
#     fmisses = Counter()

#     ## a tuple to hold my dictionaries
#     statdicts = perfile, fsamplehits, fbarhits, fmisses, fdbars

#     ## wait for jobs to finish
#     while 1:
#         fin = [i for i, j in filesort.items() if j[1].ready()]
#         data._progressbar(total, done, start, printstr)
#         time.sleep(0.1)
#         if total == done:
#             break

#         ## cleanup
#         for key in fin:
#             tup = filesort[key]
#             if tup[1].successful():
#                 pfile = tup[1].result()
#                 handle = tup[0]
#                 if pfile:
#                     ## check if this needs to return data
#                     putstats(pfile, handle, statdicts)
#                     ## purge to conserve memory
#                     del filesort[key]
#                     done += 1
#             else:
#                 ip.logger.debug(tup[1].exception())
#                 raise IPyradError(tup[1].exception())
#     print("")
#     return statdicts         




# used by splitfiles()
def zcat_make_temps(data, ftup, num, tmpdir, optim, start):
    """ 
    Call bash command 'cat' and 'split' to split large files into 4 bits.
    """
    # read it, is it gzipped?
    catcmd = ["cat"]
    if ftup[0].endswith(".gz"):
        catcmd = ["gunzip", "-c"]

    # get reading commands for r1s, r2s
    cmd1 = catcmd + [ftup[0]]
    cmd2 = catcmd + [ftup[1]]

    # make name prefix
    chunk1 = os.path.join(tmpdir, "chunk1_{}_".format(num))
    chunk2 = os.path.join(tmpdir, "chunk2_{}_".format(str(num)))
 
    # command to split and write to prefix
    cmd3 = ["split", "-a", "4", "-l", str(int(optim) * 4), "-", chunk1]
    cmd4 = ["split", "-a", "4", "-l", str(int(optim) * 4), "-", chunk2]

    # start 'split ... | gunzip -c rawfile'
    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc3 = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc1.stdout)
    res = proc3.communicate()[0]
    if proc3.returncode:
        raise IPyradError("error in zcat_make_temps:\n{}".format(res))
    
    # grab output handle results from read1s
    chunks1 = sorted(glob.glob(chunk1 + "*"))

    # repeat for paired reads
    if "pair" in data.params.datatype:
        proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc4 = sps.Popen(cmd4, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc2.stdout)
        res = proc4.communicate()[0]
        if proc4.returncode:
            raise IPyradError("error in zcat_make_temps:\n{}".format(res))            
        chunks2 = sorted(glob.glob(chunk2 + "*"))   
    else:
        chunks2 = [0] * len(chunks1)

    # ensure r1==r2
    assert len(chunks1) == len(chunks2), "Different number of R1 and R2 files"

    # ensure full progress bar b/c estimates njobs could be off
    return list(zip(chunks1, chunks2))



## EXPERIMENTAL; not yet implemented
def barmatch2(data, tups, cutters, longbar, matchdict, fnum):
    """
    cleaner barmatch func...
    """

    ## how many reads to store before writing to disk
    waitchunk = int(1e6)
    ## pid name for this engine
    epid = os.getpid()

    ## counters for total reads, those with cutsite, and those that matched
    filestat = np.zeros(3, dtype=np.int)
    ## store observed sample matches
    samplehits = {}
    ## dictionaries to store first and second reads until writing to file
    dsort1 = {} 
    dsort2 = {} 
    ## dictionary for all bars matched in sample
    dbars = {} 

    ## fill for sample names
    for sname in data.barcodes:
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]
        samplehits[sname] = 0
        dsort1[sname] = []
        dsort2[sname] = []
        dbars[sname] = set()

    ## store observed bars
    barhits = {}
    for barc in matchdict:
        barhits[barc] = 0

    ## store others
    misses = {}
    misses['_'] = 0

    ## build func for finding barcode
    getbarcode = get_barcode_func(data, longbar)

    ## get quart iterator of reads
    if tups[0].endswith(".gz"):
        ofunc = gzip.open
    else:
        ofunc = open

    ## create iterators 
    ofile1 = ofunc(tups[0], 'r')
    fr1 = iter(ofile1) 
    quart1 = izip(fr1, fr1, fr1, fr1)
    if tups[1]:
        ofile2 = ofunc(tups[1], 'r')
        fr2 = iter(ofile2)  
        quart2 = izip(fr2, fr2, fr2, fr2)
        quarts = izip(quart1, quart2)
    else:
        quarts = izip(quart1, iter(int, 1))

    ## go until end of the file
    while 1:
        try:
            read1, read2 = next(quarts)
            read1 = list(read1)
            filestat[0] += 1
        except StopIteration:
            break
    
        barcode = ""
        ## Get barcode_R2 and check for matching sample name
        if '3rad' in data.params.datatype:
            ## Here we're just reusing the findbcode function
            ## for R2, and reconfiguring the longbar tuple to have the
            ## maxlen for the R2 barcode
            ## Parse barcode. Use the parsing function selected above.
            barcode1 = find3radbcode(cutters=cutters, 
                                longbar=longbar, read1=read1)
            barcode2 = find3radbcode(cutters=cutters, 
                                longbar=(longbar[2], longbar[1]), read1=read2)
            barcode = barcode1 + "+" + barcode2
        else:
            ## Parse barcode. Uses the parsing function selected above.
            barcode = getbarcode(cutters, read1, longbar)
   
        ## find if it matches 
        sname_match = matchdict.get(barcode)

        if sname_match:
            #sample_index[filestat[0]-1] = snames.index(sname_match) + 1
            ## record who matched
            dbars[sname_match].add(barcode)
            filestat[1] += 1
            filestat[2] += 1
            samplehits[sname_match] += 1
            barhits[barcode] += 1
            if barcode in barhits:
                barhits[barcode] += 1
            else:
                barhits[barcode] = 1
    
            ## trim off barcode
            lenbar = len(barcode)
            if '3rad' in data.params.datatype:
                ## Iff 3rad trim the len of the first barcode
                lenbar = len(barcode1)
    
            if data.params.datatype == '2brad':
                overlen = len(cutters[0][0]) + lenbar + 1
                read1[1] = read1[1][:-overlen] + "\n"
                read1[3] = read1[3][:-overlen] + "\n"
            else:
                read1[1] = read1[1][lenbar:]
                read1[3] = read1[3][lenbar:]
    
            ## Trim barcode off R2 and append. Only 3rad datatype
            ## pays the cpu cost of splitting R2
            if '3rad' in data.params.datatype:
                read2 = list(read2)
                read2[1] = read2[1][len(barcode2):]
                read2[3] = read2[3][len(barcode2):]
    
            ## append to dsort
            dsort1[sname_match].append("".join(read1))
            if 'pair' in data.params.datatype:
                dsort2[sname_match].append("".join(read2))

        else:
            misses["_"] += 1
            if barcode:
                filestat[1] += 1

        ## how can we make it so all of the engines aren't trying to write to
        ## ~100-200 files all at the same time? This is the I/O limit we hit..
        ## write out at 100K to keep memory low. It is fine on HPC which can 
        ## write parallel, but regular systems might crash
        if not filestat[0] % waitchunk:
            ## write the remaining reads to file"
            writetofile(data, dsort1, 1, epid)
            if 'pair' in data.params.datatype:
                writetofile(data, dsort2, 2, epid)
            ## clear out dsorts
            for sample in data.barcodes:
                if "-technical-replicate-" in sname:
                    sname = sname.rsplit("-technical-replicate", 1)[0]
                dsort1[sname] = []
                dsort2[sname] = []
            ## reset longlist
            #longlist = np.zeros(waitchunk, dtype=np.uint32)                

    ## close open files
    ofile1.close()
    if tups[1]:
        ofile2.close()

    ## write the remaining reads to file
    writetofile(data, dsort1, 1, epid)
    if 'pair' in data.params.datatype:
        writetofile(data, dsort2, 2, epid)

    ## return stats in saved pickle b/c return_queue is too small
    ## and the size of the match dictionary can become quite large
    samplestats = [samplehits, barhits, misses, dbars]
    outname = os.path.join(data.dirs.fastqs, "tmp_{}_{}.p".format(epid, fnum))
    with open(outname, 'wb') as wout:
        pickle.dump([filestat, samplestats], wout)

    return outname


## GLOBALS
NO_RAWS = """\
    No data found in {}. Fix path to data files.
    """

OVERWRITING_FASTQS = """\
{spacer}[force] overwriting fastq files previously created by ipyrad.
{spacer}This _does not_ affect your original/raw data files."""



NOT_TWO_PATHS = """\
    Error: Must enter a raw_fastq_path or sorted_fastq_path, but not both.
    """
NO_SEQ_PATH_FOUND = """\
    Error: Step 1 requires that you enter one of the following:
        (1) a sorted_fastq_path
        (2) a raw_fastq_path + barcodes_path
    """
SAMPLES_EXIST = """\
    Error: {} Samples already found in Assembly {}.
    (Use force argument to overwrite)
    """
NO_SUPPORT_FOR_BZ2 = """\
    Found bz2 formatted files in 'sorted_fastq_path': {}
    ipyrad does not support bz2 files. The only supported formats for samples
    are .gz, .fastq, and .fq. The easiest thing to do is probably go into
    your sorted_fastq_path directory and issue this command `bunzip2 *`. You
    will probably also need to update your params file to reflect the fact
    that sample raw files now probably end with .fq or .fastq.
    """
NO_FILES_FOUND_PAIRS = """\
    No files found in 'sorted_fastq_path': {}
    Check that file names match the required convention for paired datatype
    i.e., paired file names should be identical save for _R1_ and _R2_
    (note the underscores before AND after R*).
    """
R1_R2_name_error = """\
    Paired file names must be identical except for _R1_ and _R2_.
    We detect {} R1 files and {} R2 files.
    """
NAMES_LOOK_PAIRED_WARNING = """\
    Warning: '_R2_' was detected in a file name, which suggests the data may
    be paired-end. If so, you should set the parameter 'datatype' to a paired
    option (e.g., pairddrad or pairgbs) and re-run step 1, which will require
    using the force flag (-f) to overwrite existing data.
    """
