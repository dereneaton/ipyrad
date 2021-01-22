#!/usr/bin/env python

""" 
A simpler form of the previous demultiplex.py script. 
Raw sequence is assigned to separate fastq files for each sample 
in a barcode map.
"""

# py2/3 compatible imports
from __future__ import print_function
from builtins import range

# external imports
import os
import io
import gzip
import glob
import shutil
import subprocess as sps

from loguru import logger
import numpy as np
from ipyrad.core.sample import Sample
from ipyrad.assemble.demux_utils import BarMatch, Stats
from ipyrad.assemble.utils import IPyradError, ambigcutters, BADCHARS
from ipyrad.assemble.utils import AssemblyProgressBar



class Step1:
    """
    Select step1 function to get fastq reads for each sample.
    """
    def __init__(self, data, force, ipyclient):
        # store attrs
        self.data = data
        self.force = force
        self.ipyclient = ipyclient

        # check input data files
        self.sfiles = self.data.params.sorted_fastq_path
        self.rfiles = self.data.params.raw_fastq_path
        self.print_headers()
        self.select_method()
        self.setup_dirs()


    def run(self):
        """ 
        Select and run the step1 mode: loading or demux'ing data
        """
        if self.method == "link_fastqs":
            FileLinker(self).run()
        else:
            SimpleDemux(self).run()
        self.data.save()


    def print_headers(self):
        """
        print the header for the CLI
        """
        # print headers
        if self.data._cli:
            if self.sfiles:
                self.data._print(
                    "\n  Step 1: Loading sorted fastq data to Samples")
            else:
                self.data._print(
                    "\n  Step 1: Demultiplexing fastq data to Samples")


    def select_method(self):
        """
        Checks file paths and for existing samples and returns function
        """
        # do not allow both a sorted_fastq_path and a raw_fastq
        if self.sfiles and self.rfiles:
            raise IPyradError(
                "Must enter a raw_fastq_path or sorted_fastq_path, "
                "but not both.")

        # but also require that at least one exists
        if not (self.sfiles or self.rfiles):
            raise IPyradError(
                "Step 1 requires that you enter one of the following:\n"
                "    (1) a sorted_fastq_path\n"
                "    (2) a raw_fastq_path + barcodes_path\n")

        # if Samples already exist then bail out unless force
        if self.data.samples:
            if not self.force:
                raise IPyradError((
                    "{} Samples already found in Assembly {}.\n"
                    "(Use force argument to overwrite)")
                    .format(len(self.data.samples), self.data.name)
                )
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
        logger.debug("Step 1 method: {}".format(self.method))


    def setup_dirs(self):
        """
        create output directory, tmp directory.
        """
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

        # ensure project dir exists
        if not os.path.exists(self.data.params.project_dir):
            os.mkdir(self.data.params.project_dir)

        # ensure fastq dir exists, but don't make the directory if linking
        # because it's empty and unused
        if not os.path.exists(self.data.dirs.fastqs):
            if not self.method == "link_fastqs":
                os.mkdir(self.data.dirs.fastqs)





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
        self.lbview = step.ipyclient.load_balanced_view()


    def run(self):
        """
        checks files and then counts reads in each one.
        """
        # checks for bad names and fills self.ftuples with file handles
        self.check_files()
        # distributes jobs to parallel
        self.remote_run_linker()


    def check_files(self):
        """
        Check file endings, check for PE matching
        """
        # Assert files are not .bz2 format
        if any([i.endswith(".bz2") for i in self.fastqs]):
            raise IPyradError((
                "Found bz2 formatted files in 'sorted_fastq_path': {} "
                "ipyrad does not support bz2 files. The only supported "
                "formats for samples are .gz, .fastq, and .fq. The easiest "
                "thing to do is probably go into your sorted_fastq_path "
                "directory and issue this command `bunzip2 *`. You "
                "will probably also need to update your params file to "
                "reflect the fact that sample raw files now probably end "
                "with .fq or .fastq.")
                .format(self.fastqs))

        # filter out any files without proper file endings. Raise if None
        endings = ("gz", "fastq", "fq")
        self.fastqs = [i for i in self.fastqs if i.split(".")[-1] in endings]
        if not self.fastqs:
            raise IPyradError((
                "No files found in 'sorted_fastq_path': {}\n"
                "Check that file names match the required convention for "
                "paired datatype, i.e., paired file names should be "
                "identical save for _R1_ and _R2_ (note the underscores "
                "before AND after R*).")
                .format(self.data.params.sorted_fastq_path))

        # link pairs into tuples
        if 'pair' in self.data.params.datatype:
            # check that names fit the paired naming convention
            r1s = [i for i in self.fastqs if "_R1_" in i]
            r2s = [i for i in self.fastqs if "_R2_" in i]

            # file checks
            if not r1s:
                raise IPyradError(
                    "No fastqs files found. File names must contain '_R1_' "
                    "(and '_R2_' for paired data). See Docs.")
            if len(r1s) != len(r2s):
                raise IPyradError((
                    "Paired file names must be identical except for "
                    "_R1_ and _R2_. We detect {} R1 files and {} R2 files."
                    ).format(len(r1s), len(r2s)))

            # store tuples                    
            self.ftuples = []
            for r1file in r1s:
                r2file = r1file.replace("_R1_", "_R2_")
                if not os.path.exists(r2file):
                    raise IPyradError(
                        "Expected R2 file {} to match R1 file {}"
                        .format(r1file, r2file)
                        )
                self.ftuples.append((r1file, r2file))

        # data are not paired, create empty tuple pair
        else:
            # print warning if _R2_ is in names when not paired
            if any(["_R2_" in i for i in self.fastqs]):
                message = (
                    "'_R2_' was detected in a file name, which suggests the "
                    "data may be paired-end. If so, you should set the "
                    "parameter 'datatype' to a paired option (e.g., "
                    "pairddrad or pairgbs) and re-run step 1, which will "
                    "require using the force flag (-f) to overwrite "
                    "existing data.")
                logger.warning(message)
            self.ftuples = [(i, "") for i in self.fastqs]
        logger.debug("ftuples: {}".format(self.ftuples))


    def remote_run_linker(self):
        """
        Read in fastq files and count nreads for stats and chunking in s2.
        """
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
                # submit job to count lines and store async
                rasyncs[sample.name] = self.lbview.apply(
                    zbufcountlines, sample.files.fastqs[0][0])

        # wait for all to finish
        printstr = ("loading reads       ", "s1")
        prog = AssemblyProgressBar(rasyncs, None, printstr, self.data)
        prog.block()
        prog.check()

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
                    .format(createdinc, len(self.data.samples))
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
                .astype(np.int)
                .to_string(outfile))





class SimpleDemux:
    """
    Demultiplexer that uses only a single engine and processes one
    file at a time, since I/O limits are generally limiting.
    """
    def __init__(self, step):

        # store the step object inputs
        self.data = step.data
        self.input = step.rfiles
        self.fastqs = glob.glob(self.input)

        # get parallel info
        self.iview = step.ipyclient.load_balanced_view(targets=[0])

        # re-parse the barcodes file in case hackers options changed
        # e.g., (i7 demuxing or merge technical replicates)
        self.data._link_barcodes()

        # attrs filled by check_files
        self.ftuples = []
        self.longbar = None
        self.check_files()

        # attrs filled by get_barcode_dict
        self.cutters = None
        self.matchdict = {}       
        self.get_barcode_dict()        

        # store stats for each file handle (grouped results of chunks)
        self.stats = Stats()    


    def check_files(self):
        """
        Check .fastqs files are present and formatted correctly as .ftuples
        """
        # check for data using glob for fuzzy matching
        if not self.fastqs:
            raise IPyradError(
                "No data found in {}. Fix path to data files."
                .format(self.data.params.raw_fastq_path))

        # find longest barcode
        if not os.path.exists(self.data.params.barcodes_path):
            raise IPyradError(
                "Barcodes file not found. You entered: '{}'"
                .format(self.data.params.barcodes_path))

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
        logger.debug('longbar: {}'.format(self.longbar))            
        logger.debug('ftuples: {}'.format(self.ftuples))



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
        assert self.cutters, (
            "Must enter a restriction_overhang for demultiplexing.")

        # get matchdict
        self.matchdict = inverse_barcodes(self.data)
        logger.debug("cutters: {}".format(self.cutters))
        logger.debug("matchdict: {}...".format(str(self.matchdict)[:50]))



    def remote_run_barmatch(self):
        """
        Submit chunks to be sorted by barmatch() and collect stats
        """
        # chunkfiles is a dict with {handle: chunkslist, ...}. The func barmatch
        # writes results to samplename files with PID number, and also writes a 
        # pickle for chunk specific results with fidx suffix, which it returns.
        rasyncs = {}
        for fidx, ftuple in enumerate(self.ftuples):

            # get filename 
            handle = ftuple[0].rsplit("_R1", 1)[0]

            # args to barmatch
            args = (
                self.data, ftuple, self.longbar,
                self.cutters, self.matchdict, fidx
            )

            # submit job to engine
            rasyncs[(handle, fidx)] = self.iview.apply(barmatch, args)

            # get ready to receive stats: 'total', 'cutfound', 'matched'
            self.stats.perfile[handle] = np.zeros(3, dtype=np.int)

        # progress bar info
        printstr = ("sorting reads       ", "s1")
        prog = AssemblyProgressBar(rasyncs, None, printstr, self.data)
        prog.block()
        prog.check()

        # collect and store results as jobs finish
        for tup in rasyncs:

            # get the handle
            handle, fidx = tup
            pkl = rasyncs[tup].get()
            
            # fill stats for this handle
            self.stats.fill_from_pickle(pkl, handle)
        del rasyncs



    def concatenate_chunks(self):
        """ 
        If multiple chunk files match to the same sample name but with 
        different barcodes (i.e., they are technical replicates) then this
        will assign all the files to the same sample name file and gzip it.
        """
        # get all the files
        ftmps = glob.glob(os.path.join(self.data.dirs.fastqs, "tmp_*.fastq"))                        

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

        # concatenate files
        snames = []
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            snames.append(sname)

        # submit jobs to remote engines
        rasyncs = {}
        for sname in set(snames):
            tmp1s = sorted(r1dict[sname])
            tmp2s = sorted(r2dict[sname])

            args = (self.data, sname, tmp1s, tmp2s)
            rasyncs[sname] = self.iview.apply(collate_files, *args)

        # collate files progress bar
        printstr = ("writing/compressing ", "s1")
        prog = AssemblyProgressBar(rasyncs, None, printstr, self.data)
        prog.block()
        prog.check()



    def store_stats(self):
        """
        Write stats and stores to Assembly object.
        """
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
        pkls = os.path.join(self.data.dirs.fastqs, "tmpstats_*.pkl")
        for pkl in glob.glob(pkls):
            os.remove(pkl)



    def run(self):
        """
        Run complete class
        """
        self.remote_run_barmatch()
        self.concatenate_chunks()
        self.store_stats()





# CALLED BY FILELINKER
def get_name_from_file(fname, splitnames, fields):
    """
    Grab Sample names from demultiplexed input fastq file names
    """
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
        raise IPyradError(
            "Found invalid/empty filename in link_fastqs. Check "
            "splitnames argument.")
    return base



def zbufcountlines(filename):
    """
    Fast line counter using unix utils
    """
    if filename.endswith(".gz"):
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



# CALLED BY GET_BARCODE_DICT
def inverse_barcodes(data):
    """ 
    Build full inverse barcodes dictionary 
    """
    matchdict = {}
    bases = set("CATGN")
    poss = set()

    # do perfect matches
    for sname in data.barcodes:
        barc = data.barcodes[sname]
        
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
                        logger.warning((
                            "\nSample: {} ({}) is within {} base changes of "
                            "sample ({}). Ambiguous barcodes that match to "
                            "both samples will arbitrarily be assigned to the "
                            "first sample. If you do not like this then lower "
                            "the value of max_barcode_mismatch and rerun "
                            "(recommended).")
                            .format(sname, barc, 
                                data.params.max_barcode_mismatch + 1,
                                matchdict.get(tbar1)
                            ))

                # if allowing two base difference things get big
                # for each modified bar, allow one modification to other bases
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



def barmatch(args):
    """
    function to call BarMatch.run() and return 
    """
    # run procesor
    tool = BarMatch(*args)
    # writes reads t ofile and writes stats to pickle
    pkl = tool.run()
    return pkl



# CALLED BY CONCATENATE_CHUNKS
def collate_files(data, sname, tmp1s, tmp2s):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    # out handle
    out1 = os.path.join(data.dirs.fastqs, "{}_R1_.fastq.gz".format(sname))
    out2 = os.path.join(data.dirs.fastqs, "{}_R2_.fastq.gz".format(sname))

    # build cmd
    cmd1 = ['cat']
    for tmpfile in tmp1s:
        cmd1 += [tmpfile]

    # get compression function
    proc = sps.Popen(['which', 'pigz'], stderr=sps.PIPE, stdout=sps.PIPE)
    pigz_bin = proc.communicate()
    if pigz_bin[0].strip():
        compress = ["pigz"]
    else:
        compress = ["gzip"]

    # call cmd
    out = io.BufferedWriter(gzip.open(out1, 'w'))    
    proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
    proc2 = sps.Popen(compress, stdin=proc1.stdout, stderr=sps.PIPE, stdout=out)
    eout, _ = proc2.communicate()
    if proc2.returncode:
        logger.exception(eout)
        raise IPyradError("error in collate_files R1 {}".format(eout))
    proc1.stdout.close()
    out.close()

    # then cleanup
    for tmpfile in tmp1s:
        os.remove(tmpfile)

    # do R2 files
    if 'pair' in data.params.datatype:     
        # build cmd
        cmd1 = ['cat']
        for tmpfile in tmp2s:
            cmd1 += [tmpfile]

        # call cmd
        out = io.BufferedWriter(gzip.open(out2, 'w'))
        proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE)
        proc2 = sps.Popen(compress, stdin=proc1.stdout, stderr=sps.PIPE, stdout=out)
        err = proc2.communicate()
        if proc2.returncode:
            logger.exception(err[0])
            raise IPyradError("error in collate_files R2 {}".format(err[0]))
        proc1.stdout.close()
        out.close()

        # then cleanup
        for tmpfile in tmp2s:
            os.remove(tmpfile)




if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("INFO")

    # test FileLinker by loading SE fastq files
    # tdata = ip.Assembly("test-emprad")
    # tdata.params.project_dir = "/tmp"
    # tdata.params.sorted_fastq_path = "/home/deren/Documents/ipyrad/sandbox/oak-fastqs/*.gz"
    # tdata.params.datatype = "rad"
    # tdata.run("1", auto=True, force=True)
    # print(tdata.stats.head())

    CURDIR = os.path.dirname(__file__)
    SIM_PREFIX = os.path.join(CURDIR, "../../tests/ipsimdata/")

    # test DEMUX by loading SE fastq files
    # tdata = ip.Assembly("test-simrad")
    # tdata.params.project_dir = "/tmp"
    # tdata.params.raw_fastq_path = SIM_PREFIX + "rad_example_R*.gz"
    # tdata.params.barcodes_path = SIM_PREFIX + "rad_example_barcodes.txt"
    # tdata.params.datatype = "rad"
    # tdata.run("1", auto=True, force=True)
    # print(tdata.stats.head())


    # # test DEMUX by loading PE-ddRAD fastq files
    # tdata = ip.Assembly("test-simpairddrad")
    # tdata.params.project_dir = "/tmp"
    # tdata.params.raw_fastq_path = SIM_PREFIX + "pairddrad_example_R*.gz"
    # tdata.params.barcodes_path = SIM_PREFIX + "pairddrad_example_barcodes.txt"
    # tdata.params.datatype = "pairddrad"
    # tdata.params.restriction_overhang = ("TGCAG", "CGG")
    # tdata.run("1", auto=True, force=True)
    # print(tdata.stats.head())


    # test by loading messy PE-GBS fastq files
    tdata = ip.Assembly("test-amaranth")
    tdata.params.project_dir = "/tmp"
    tdata.params.sorted_fastq_path = "/home/deren/Documents/kmerkit/data/amaranths/hybrid*.gz"
    tdata.params.assembly_method = "denovo"
    tdata.params.datatype = "pair3rad"
    tdata.params.restriction_overhang = ("ATCGG", "TAGCTT")
    tdata.params.filter_adapters = 2
    tdata.run("1", auto=True, force=True)
    print(tdata.stats.head())
   

    # # test by loading messy PE-GBS fastq files
    # tdata = ip.Assembly("test-emppairgbs")
    # tdata.params.project_dir = "/tmp"
    # tdata.params.sorted_fastq_path = "/home/deren/Desktop/messy/Cyathophora/*.gz"
    # tdata.params.barcodes_path = "/home/deren/Desktop/messy/Pedicularis-PEGBS.barcodes.txt"    
    # tdata.params.datatype = "pairgbs"
    # tdata.params.restriction_overhang = ("TGCAG", "TGCAG")
    # tdata.params.filter_adapters = 2
    # tdata.run("1", auto=True, force=True)
    # print(tdata.stats.head())
