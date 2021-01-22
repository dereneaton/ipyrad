#!/usr/bin/env python

"""
Some utilities used in demux.py for demultiplexing.
"""

# py2/3 compat
try:
    from itertools import izip, islice
except ImportError:
    from itertools import islice
    izip = zip

import os
import gzip
import pickle
from collections import Counter
from loguru import logger
import numpy as np



class BarMatch:
    """
    Assign reads to samples based on barcodes and writes stats to 
    a pickled dict object.
    This class is created and run inside the barmatch function() that 
    is run on remote engines for parallelization. It processes a single
    large fastq file. It assigns reads from multiple barcodes to the 
    same sample name if -technical-replicates.
    """
    def __init__(self, data, ftuple, longbar, cutters, matchdict, fidx):
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
        
        # store all barcodes observed
        self.barhits = {}
        for barc in self.matchdict:
            self.barhits[barc] = 0

        # store reads and bars matched to samples
        # store reads per sample (group technical replicates)        
        self.read1s = {} 
        self.read2s = {} 
        self.dbars = {} 
        self.samplehits = {}
        for sname in self.data.barcodes:
            if "-technical-replicate-" in sname:
                sname = sname.rsplit("-technical-replicate", 1)[0]
            self.samplehits[sname] = 0
            self.read1s[sname] = []
            self.read2s[sname] = []
            self.dbars[sname] = set()

        # store counts of what didn't match to samples
        self.misses = {}
        self.misses['_'] = 0

        # get the barcode matching function
        if self.longbar[1] == 'same':
            if self.data.params.datatype == '2brad':
                self.demux = getbarcode1
            else:
                self.demux = getbarcode2
        else:
            self.demux = getbarcode3

        # file i/o handles
        self.ofile1 = None
        self.ofile2 = None
        self.quarts = None


    def run(self):
        """
        1. Open generators that pull in 4 lines at a time (fastq)
        2. 
        """
        self.open_read_generators()
        pkl = self.assign_reads()
        self.close_read_generators()
        return pkl


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
        """
        Close the file handles that were opened as generators.
        """
        self.ofile1.close()
        if self.ftuple[1]:
            self.ofile2.close()


    def assign_reads(self):
        """
        Assign reads to samples using barcode matching and the reads
        to fastq files. The statistics are stored in a dictionary which
        is pickled, and the tmp filename is returned.
        """
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
                barcode = read1[0].decode().rsplit(":", 1)[-1].split("+")[0]

            else:
                # COMBINATORIAL BARCODES (BCODE1+BCODE2)
                if '3rad' in self.data.params.datatype:
                    barcode1 = find3radbcode(self.cutters, self.longbar, read1)
                    barcode2 = find3radbcode(self.cutters, self.longbar, read2)
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
            "tmpstats_{}_{}.pkl".format(self.epid, self.fidx))
        with open(pklname, 'wb') as wout:
            pickle.dump([self.filestat, samplestats], wout)
            logger.debug("dumped stats to {}".format(pklname))

        return pklname




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


def find3radbcode(cutters, longbar, read):
    "find barcode sequence in the beginning of read"
    # default barcode string
    for ambigcuts in cutters:
        for cutter in ambigcuts:
            # If the cutter is unambiguous there will only be one.
            if not cutter:
                continue
            search = read[1][:int(longbar[0] + len(cutter) + 1)]
            splitsearch = search.decode().rsplit(cutter, 1)
            if len(splitsearch) > 1:
                return splitsearch[0]
    # No cutter found
    return splitsearch[0] 


def write_to_file(data, dsort, read, pid):
    """
    Writes sorted data to tmp files
    """
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
            "tmp_{}_{}_{}.fastq".format(sname, rrr, pid))

        # append to this sample name
        with open(handle, 'a') as out:
            out.write("".join(dsort[sname]))
        logger.debug("appending to {}".format(handle))





class Stats:
    """
    Used inside BarMatch to store stats nicely.
    """
    def __init__(self):
        # stats for each raw input file
        self.perfile = {}

        # stats for each sample
        self.fdbars = {}
        self.fsamplehits = Counter()
        self.fbarhits = Counter()
        self.fmisses = Counter()


    def fill_from_pickle(self, pkl, handle):
        """
        ...
        """
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


