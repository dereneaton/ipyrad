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
from ipyrad.assemble.utils import IPyradError, ambigcutters, BADCHARS
from ipyrad.assemble.pair_fastqs import (
    get_fastq_tuples_dict_from_paths_list,
    get_paths_list_from_fastq_str,
)


class Step1:
    def __init__(self, data, force, ipyclient):
        # store attrs
        self.data = data
        self.force = force
        self.ipyclient = ipyclient
        self.skip = False

        # check input data files
        self.sfiles = self.data.params.sorted_fastq_path
        self.rfiles = self.data.params.raw_fastq_path
        self.print_headers()
        self.select_method()
        self.setup_dirs()


    def run(self):
        if not self.skip:
            if self.method == "link_fastqs":
                FileLinker(self).run()
            else:
                Demultiplexer(self).run()
            self.data.save()


    def setup_dirs(self):
        "create output directory, tmp directory."
        # assign fastq dir
        self.data.dirs.fastqs = os.path.join(
            self.data.params.project_dir,
            self.data.name + "_fastqs")       
        self.data.dirs.fastqs = os.path.realpath(self.data.dirs.fastqs)

        # Do NOT delete any directory if you're just linking sorted fastqs
        # This allows you to reuse _fastqs from previous assemblies.
        if self.method == "link_fastqs":
            pass
        # remove existing if force flag.
        elif self.force:
            if os.path.exists(self.data.dirs.fastqs):
                shutil.rmtree(self.data.dirs.fastqs)

        # bail out if overwrite necessary but no force flag.
        else:
            if os.path.exists(self.data.dirs.fastqs):
                raise IPyradError(
                    "Fastq dir {} already exists: use force to overwrite"
                    .format(self.data.dirs.fastqs))
                self.skip = True

        # ensure project dir exists
        if not os.path.exists(self.data.params.project_dir):
            os.mkdir(self.data.params.project_dir)

        # ensure fastq dir exists, but don't make the directory if linking
        # because it's empty and unused
        if not os.path.exists(self.data.dirs.fastqs) and\
            not self.method == "link_fastqs":
            os.mkdir(self.data.dirs.fastqs)


    def print_headers(self):
        # print headers
        if self.data._cli:
            if self.sfiles:
                self.data._print(
                    "\n  Step 1: Loading sorted fastq data to Samples")
            else:
                self.data._print(
                    "\n  Step 1: Demultiplexing fastq data to Samples")


    def select_method(self):
        "Checks file paths and for existing samples and returns function"
        # do not allow both a sorted_fastq_path and a raw_fastq
        if self.sfiles and self.rfiles:
            raise IPyradError(NOT_TWO_PATHS)

        # but also require that at least one exists
        if not (self.sfiles or self.rfiles):
            raise IPyradError(NO_SEQ_PATH_FOUND)

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
    """
    Loads Samples from file names and check sample names for bad chars.
    """
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
            raise IPyradError(NO_SUPPORT_FOR_BZ2.format(self.input))

        # filter out any files without proper file endings. Raise if None
        endings = ("gz", "fastq", "fq")
        self.fastqs = [i for i in self.fastqs if i.split(".")[-1] in endings]
        if not self.fastqs:
            raise IPyradError(
                NO_FILES_FOUND_PAIRS
                .format(self.data.params.sorted_fastq_path))

        # simple reality check to verify PE data has an even number of files
        if 'pair' in self.data.params.datatype:
            if len(self.fastqs) % 2:
                raise IPyradError(PE_ODD_NUMBER_OF_FILES)

        # get list of expanded paths
        paths = get_paths_list_from_fastq_str(self.fastqs)

        # get dict of {basename: (Path-R1, Path-R2)} from filenames
        self.ftuples = get_fastq_tuples_dict_from_paths_list(paths)

        # # link pairs into tuples
        # if 'pair' in self.data.params.datatype:
        #     # check that names fit the paired naming convention
        #     r1s = [i for i in self.fastqs if "_R1_" in i]
        #     r2s = [i for i in self.fastqs if "_R2_" in i]

        #     # file checks
        #     if not r1s:
        #         raise IPyradError(
        #             "No fastqs files found. File names must contain '_R1_' "
        #             "(and '_R2_' for paired data). See Docs.")
        #     if len(r1s) != len(r2s):
        #         raise IPyradError(
        #             R1_R2_name_error.format(len(r1s), len(r2s)))

        #     # store tuples                    
        #     self.ftuples = []
        #     for r1file in r1s:
        #         r2file = r1file.replace("_R1_", "_R2_")
        #         if not os.path.exists(r2file):
        #             raise IPyradError(
        #                 "Expected R2 file {} to match R1 file {}"
        #                 .format(r1file, r2file)
        #                 )
        #         self.ftuples.append((r1file, r2file))

        # # data are not paired, create empty tuple pair
        # else:
        #     # print warning if _R2_ is in names when not paired
        #     if any(["_R2_" in i for i in self.fastqs]):
        #         print(NAMES_LOOK_PAIRED_WARNING)
        #     self.ftuples = [(i, "") for i in self.fastqs]


    def remote_run_linker(self):
        "read in fastq files and count nreads for stats and chunking in s2."

        # local counters
        createdinc = 0

        # iterate over input files
        for sname, ftup in self.ftuples.items():
            # remove file extension from name
            # sname = get_name_from_file(ftup[0], None, None)
            # print(sname, ftup)

            # Create new Sample Class objects with names from files
            if sname not in self.data.samples:
                newsamp = Sample(sname)
                newsamp.stats.state = 1
                newsamp.barcode = None
                newsamp.files.fastqs = [(str(ftup[0]), str(ftup[1]))]
                self.data.samples[sname] = newsamp
                createdinc += 1

        # send jobs to engines for counting with cat/zcat | wc
        rasyncs = {}
        if createdinc:
            for sample in self.data.samples.values():

                # get zip var
                gzipped = bool(sample.files.fastqs[0][0].endswith(".gz"))

                # submit job to count lines and store async
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
            res = rasyncs[sname].get() / 4
            self.data.samples[sname].stats.reads_raw = res
            self.data.samples[sname].stats_dfs.s1["reads_raw"] = res
            self.data.samples[sname].state = 1

        # print if data were linked
        if createdinc:
            # double for paired data
            if 'pair' in self.data.params.datatype:
                createdinc = createdinc * 2
            if self.data._cli:
                self.data._print(
                    "{} fastq files loaded to {} Samples."
                    .format(
                        createdinc, 
                        len(self.data.samples),
                    )
                )

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
                .astype(np.int64)
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

        # re-parse the barcodes file in case hackers options changed
        # e.g., (i7 demuxing or merge technical replicates)
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

        # i7 tags there will be only one barcode, so this overrides "datatype"
        # so that if you are using pair3rad if doesn't cause problems.
        # For pair3rad we need to add the length info for barcodes_R2
        if not self.data.hackersonly.demultiplex_on_i7_tags:
            if "3rad" in self.data.params.datatype:
                blens = [
                    len(i.split("+")[1]) for i in self.data.barcodes.values()
                ]
                # Record if bar1 and bar2 are different lengths
                if self.longbar[0] != max(blens):
                    self.longbar = (self.longbar[0], 'diff', max(blens))
                else:
                    self.longbar = (self.longbar[0], 'same', max(blens))

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
        """
        Checks sample names and replaces bad chars in dict with _
        And returns a list of both resolutions of cut site 1 for ambigs.
        # (TGCAG, ) ==> [TGCAG, ]
        # (TWGC, ) ==> [TAGC, TTGC]
        # (TWGC, AATT) ==> [TAGC, TTGC]
        """
        # expand ambigs
        self.cutters = [
            ambigcutters(i) for i in self.data.params.restriction_overhang
        ]
        assert self.cutters, "Must enter a restriction_overhang for demultiplexing."

        # get matchdict
        self.matchdict = inverse_barcodes(self.data)


    def setup_for_splitting(self, omin=int(8e6)):
        """
        Decide to split or not based on whether 1/16th of file size is 
        bigger than omin, which is default to 8M reads.
        """
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
        # for each file submitted we expect it to create 16 or 32 files.
        if rasyncs:
            while 1:
                # break when all jobs are finished
                if all([i.ready() for i in rasyncs.values()]):
                    break

                # ntemp files written or being written
                done = len(glob.glob(os.path.join(self.tmpdir, "chunk*_*_*")))
                self.data._progressbar(njobs, done, start, printstr)
                time.sleep(0.5)

            # store results
            for key, val in rasyncs.items():
                chunksdict[key] = val.get()       

            # clean up                    
            self.ipyclient.purge_everything()                    
            self.data._progressbar(10, 10, start, printstr)
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
                    self.data,
                    ftuple,
                    self.longbar,
                    self.cutters,
                    self.matchdict,
                    fidx,
                    )
                rasync = self.lbview.apply(barmatch, args)
                rasyncs[ridx] = (handle, rasync)
                ridx += 1

            # get ready to receive stats: 'total', 'cutfound', 'matched'
            self.stats.perfile[handle] = np.zeros(3, dtype=np.int64)

        # collect and store results as jobs finish
        njobs = len(rasyncs)
        done = 0
        while 1:
            # get list of ridx numbers for finished jobs
            finished = [i for (i, j) in rasyncs.items() if j[1].ready()]

            # cleanup finished ridx jobs and grab stats
            for ridx in finished:
                handle, rasync = rasyncs[ridx]
                pkl = rasync.get()
                self.stats.fill_from_pickle(pkl, handle)
                del rasyncs[ridx]
                done += 1

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
        outfile.write(
            "{:<35}  {:>13}{:>13}{:>13}\n"
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
        outfile.write(
            "\n{:<35}  {:>13}\n"
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
                                int(self.stats.fbarhits[offhit] / 2))
                            )
                        #sumoffhits += fbarhits[offhit]
            
                # write string to file
                outfile.write("{:<35}  {:>13} {:>13} {:>13}\n"
                    .format(sname, hit, hit, 
                        int(self.stats.fbarhits[hit] / 2)))
                outfile.write(offhitstring)
            
        # write misses
        misskeys = list(self.stats.fmisses.keys())
        misskeys.sort(key=self.stats.fmisses.get)
        for key in misskeys[::-1]:
            outfile.write('{:<35}  {:>13} {:>13} {:>13}\n'
                .format("no_match", "_", key, self.stats.fmisses[key]))
        outfile.close()        

        # Link Sample with this data file to the Assembly object
        for sname in snames:

            # make the sample
            sample = Sample(sname)

            # allow multiple barcodes if its a replicate. 
            barcodes = []
            for n in range(500):
                fname = "{}-technical-replicate-{}".format(sname, n)
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
        """
        Gzips are always bytes so let's use rb to make unzipped also bytes.
        """
        # get file type
        if self.ftuple[0].endswith(".gz"):
            self.ofile1 = gzip.open(self.ftuple[0], 'rb')
        else:
            self.ofile1 = open(self.ftuple[0], 'rb')

        # create iterators 
        fr1 = iter(self.ofile1) 
        quart1 = izip(fr1, fr1, fr1, fr1)

        # create second read iterator for paired data
        if self.ftuple[1]:
            if self.ftuple[0].endswith(".gz"):
                self.ofile2 = gzip.open(self.ftuple[1], 'rb')
            else:
                self.ofile2 = open(self.ftuple[1], 'rb')

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

            # i7 barcodes (get from name string instead of sequence)
            if self.data.hackersonly.demultiplex_on_i7_tags:
                barcode = read1[0].decode().rsplit(":", 1)[-1].split("+")[0].strip()

            else:
                # COMBINATORIAL BARCODES (BCODE1+BCODE2)
                if '3rad' in self.data.params.datatype:
                    barcode1 = find3radbcode(self.cutters, self.longbar[0], read1)
                    barcode2 = find3radbcode(self.cutters, self.longbar[2], read2)
                    barcode = barcode1 + "+" + barcode2

                # USE BARCODE PARSER: length or splitting
                else:
                    # Parse barcode. Uses the parsing function selected above.
                    barcode = self.demux(self.cutters, read1, self.longbar)

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
                lenbar1 = len(barcode)
                if '3rad' in self.data.params.datatype:
                    ## Iff 3rad trim the len of the first barcode
                    lenbar1 = len(barcode1)
                    lenbar2 = len(barcode2)

                # no trim on i7 demux
                if self.data.hackersonly.demultiplex_on_i7_tags:
                    lenbar1 = lenbar2 = 0

                # for 2brad we trim the barcode AND the synthetic overhang
                # The `+1` is because it trims the newline
                if self.data.params.datatype == '2brad':
                    overlen = len(self.cutters[0][0]) + lenbar1 + 1
                    read1[1] = read1[1][:-overlen] + b"\n"
                    read1[3] = read1[3][:-overlen] + b"\n"
                else:
                    read1[1] = read1[1][lenbar1:]
                    read1[3] = read1[3][lenbar1:]

                # Trim barcode off R2 and append. Only 3rad datatype
                # pays the cpu cost of splitting R2
                if '3rad' in self.data.params.datatype:
                    read2 = list(read2)
                    read2[1] = read2[1][lenbar2:]
                    read2[3] = read2[3][lenbar2:]

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
                write_to_file(self.data, self.read1s, 1, self.epid)
                if 'pair' in self.data.params.datatype:
                    write_to_file(self.data, self.read2s, 2, self.epid)

                # clear out lits of sorted reads
                for sname in self.data.barcodes:
                    if "-technical-replicate-" in sname:
                        sname = sname.rsplit("-technical-replicate", 1)[0]
                    self.read1s[sname] = []
                    self.read2s[sname] = []             

        ## write the remaining reads to file
        write_to_file(self.data, self.read1s, 1, self.epid)
        if 'pair' in self.data.params.datatype:
            write_to_file(self.data, self.read2s, 2, self.epid)

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


# CALLED BY FILELINKER
def get_name_from_file(fname, splitnames, fields):
    "Grab Sample names from demultiplexed input fastq file names"

    # allowed extensions
    file_extensions = [".gz", ".fastq", ".fq", ".fasta", ".clustS", ".consens"]
    base, _ = os.path.splitext(os.path.basename(fname))

    # remove read number from name
    base = base.replace("_R1_.", ".")\
               .replace("_R1_", "")\
               .replace("_R1.", ".")

    # To test running pe data as concatenated SE
    # 3/8/20 iao
    #base = base.replace("_R2_.", ".")\
    #           .replace("_R2_", "")\
    #           .replace("_R2.", ".")

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

    # replace any bad characters from name with _
    base = "".join([
        i.replace(i, "_") if i in BADCHARS else i for i in base
    ])        

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
    # default barcode string
    for ambigcuts in cutters:
        for cutter in ambigcuts:
            # If the cutter is unambiguous there will only be one.
            if not cutter:
                continue
            search = read[1][:int(longbar + len(cutter) + 1)]
            splitsearch = search.decode().rsplit(cutter, 1)
            if len(splitsearch) > 1:
                return splitsearch[0]
    # No cutter found
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

        # bytes-strings!
        search = read1[1][:int(longbar[0] + len(cutter) + 1)]

        try:
            search = search.decode()
        except (AttributeError, TypeError):
            pass

        try:
            cutter = cutter.decode()
        except (AttributeError, TypeError):
            pass

        try:
            barcode = search.rsplit(cutter, 1)
        except (AttributeError, TypeError):
            barcode = search.decode().rsplit(cutter, 1)

        if len(barcode) > 1:
            return barcode[0]
    ## No cutter found
    return barcode[0] 



def write_to_file(data, dsort, read, pid):
    "Writes sorted data to tmp files"
    if read == 1:
        rrr = "R1"
    else:
        rrr = "R2"

    # appends to file for each sample, avoids parallel fighting by using 
    # pid assigned file handle.
    for sname in dsort:

        # file out handle
        handle = os.path.join(
            data.dirs.fastqs, 
            "tmpdir",
            "tmp_{}_{}_{}.fastq".format(sname, rrr, pid))

        # append to this sample name
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))



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

    # do perfect matches
    for sname, barc in data.barcodes.items():
        
        # remove -technical-replicate-N if present
        if "-technical-replicate-" in sname:
            sname = sname.rsplit("-technical-replicate", 1)[0]

        # store {barcode: name} mapping
        matchdict[barc] = sname

        # record that this barcodes has been seen
        poss.add(barc)

        # get x-off barcodes 
        if data.params.max_barcode_mismatch:
            
            # iterate over bases in the barcode
            for idx1, base in enumerate(barc):
                
                # get bases that are not this one
                diffs = bases.difference(base)
                
                # iter over the ways this barcode has other bases as this pos.
                for diff in diffs:
                    lbar = list(barc)
                    lbar[idx1] = diff
                    tbar1 = "".join(lbar)

                    # if this new barcode has not been observed store it.
                    if tbar1 not in poss:
                        matchdict[tbar1] = sname
                        poss.add(tbar1)

                    # if it has been seen in another taxon, problem.
                    else:
                        print("""\n
        Warning: 
        Sample: {} ({})
        is within {} base changes of sample ({})
        Ambiguous barcodes that match to both samples will arbitrarily 
        be assigned to the first sample. If you do not like this then 
        lower the value of max_barcode_mismatch and rerun (recommended).
        """.format(
            sname,
            barc, 
            data.params.max_barcode_mismatch + 1,
            matchdict.get(tbar1), 
            )
        )

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
    try:
        dat = b"".join(islice(infile, 40000))
    except TypeError:
        dat = "".join(islice(infile, 40000))
    outfile.write(dat)
    outfile.close()
    infile.close()

    ## Get the size of the tmp file
    tmp_size = os.path.getsize(tmp_file_name)

    # divide by the tmp file size and multiply by 10000 to approximate
    # the size of the input .fq files
    # If the input file is less than 40000 lines long the islice returns
    # 0, so here just return the full file size because this is toy data.
    try:
        inputreads = int(insize / tmp_size) * 10000
    except ZeroDivisionError:
        inputreads = insize
    os.remove(tmp_file_name)

    return inputreads


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
    proc1 = sps.Popen(
        cmd1, 
        stderr=sps.STDOUT, 
        stdout=sps.PIPE, 
        universal_newlines=True)
    proc3 = sps.Popen(
        cmd3, 
        stderr=sps.STDOUT, 
        stdout=sps.PIPE, 
        stdin=proc1.stdout, 
        universal_newlines=True)
    res = proc3.communicate()[0]
    if proc3.returncode:
        raise IPyradError("error in zcat_make_temps:\n{}".format(res))
    
    # grab output handle results from read1s
    chunks1 = sorted(glob.glob(chunk1 + "*"))

    # repeat for paired reads
    if "pair" in data.params.datatype:
        proc2 = sps.Popen(
            cmd2, 
            stderr=sps.STDOUT, 
            stdout=sps.PIPE, 
            universal_newlines=True)
        proc4 = sps.Popen(
            cmd4, 
            stderr=sps.STDOUT, 
            stdout=sps.PIPE, 
            stdin=proc2.stdout, 
            universal_newlines=True)
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
PE_ODD_NUMBER_OF_FILES = """\
    Paired-end datatype indicated by `datatype` parameter, but
    `sorted_fastq_path` contains an odd number of files. Please check files
    in this path to ensure it includes _only_ R1/R2 paired-end .gz files.
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


if __name__ == "__main__":

    import ipyrad as ip

    # test to load samples w/ names "*_R1_*", "*_R2_*"
    # data = ip.Assembly("TEST")
    # data.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    # data.params.project_dir = "/tmp/9test"
    # data.run("1", force=True, auto=True)
    # print(data.stats)

    # test to load samples w/ names "*_R1.*", "*_R2.*"
    data = ip.Assembly("TEST")
    data.params.sorted_fastq_path = "../../pedtest/DEMUX_fastqs/integ*.gz"
    data.params.project_dir = "/tmp/9test"
    data.run("1", force=True, auto=True)
    print(data.stats)

    # test to load samples w/ names "*_1.*", "*_2.*"
    # ...
