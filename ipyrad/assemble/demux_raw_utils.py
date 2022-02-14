#!/usr/bin/env python

"""Some utilities used in demux.py for demultiplexing.
"""

from typing import Dict, Tuple, TypeVar
import io
from pathlib import Path
import gzip
import subprocess as sps
from collections import Counter
from dataclasses import dataclass
from loguru import logger
from ipyrad.assemble.utils import IPyradError

Assembly = TypeVar("Assembly")


class BarMatch:
    """
    Assign reads to samples based on barcode matching.

    This class is created and run inside the barmatch function() that 
    is run on remote engines for parallelization. It processes a single
    large fastq file. It assigns reads from multiple barcodes to the 
    same sample name if -technical-replicates in name.
    """
    def __init__(self, data, barcodes, ftuple, longbar, cutters, matchdict, fidx):
        # store attrs
        self.data = data
        self.longbar = longbar
        self.cutters = cutters
        self.ftuple = ftuple
        self.matchdict = matchdict
        self.barcodes = barcodes
        self.fidx = fidx

        # when to write to disk and unique suffix for tmpfile
        self.chunksize = int(1e6) 
        
        # store reads until chunk is ready to write to disk
        self.read1s = {} 
        self.read2s = {} 

        # file stats
        self.nreads = 0
        self.cut_detected = 0
        self.nmatched = 0

        # store all observed barcodes counts
        self.barcode_counts = Counter()

        # store nreads per sample
        self.sample_to_counts = Counter()

        # store misses
        self.misses = Counter()

        # set defaults for collections while grouping tech reps
        for sname in self.barcodes:
            # if "-technical-replicate-" in sname:
                # sname = sname.rsplit("-technical-replicate", 1)[0]
            self.read1s[sname] = []
            self.read2s[sname] = []

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
        self.assign_reads()
        self.close_read_generators()
        return (
            self.nreads, 
            self.nmatched, 
            self.sample_to_counts, 
            self.barcode_counts,
            self.misses,
        )


    def open_read_generators(self):
        """
        Gzips are always bytes so let's use rb to make unzipped also bytes.
        """
        # get file type
        if self.ftuple[0].suffix == ".gz":
            self.ofile1 = gzip.open(self.ftuple[0], 'r', encoding="utf-8")
        else:
            self.ofile1 = open(self.ftuple[0], 'r', encoding="utf-8")

        # create iterators 
        quart1 = zip(self.ofile1, self.ofile1, self.ofile1, self.ofile1)

        # create second read iterator for paired data
        if self.ftuple[1]:
            if self.ftuple[0].endswith(".gz"):
                self.ofile2 = gzip.open(self.ftuple[1], 'r', encoding="utf-8")
            else:
                self.ofile2 = open(self.ftuple[1], 'r', encoding="utf-8")

            # create iterators
            quart2 = zip(self.ofile2, self.ofile2, self.ofile2, self.ofile2)
            self.quarts = zip(quart1, quart2)
        else:
            self.quarts = zip(quart1, iter(int, 1))

    def close_read_generators(self):
        """Close the file handles that were opened as generators."""
        self.ofile1.close()
        if self.ftuple[1]:
            self.ofile2.close()

    def assign_reads(self):
        """Assign reads to samples using barcode matching.

        Writes reads to reads to fastq files. The statistics are 
        stored in a dictionary which is pickled, and the tmp filename
        is returned.
        """
        while 1:
            # read in four lines of data and increase counter
            try:
                read1, read2 = next(self.quarts)
                read1 = list(read1)
                self.nreads += 1
            except StopIteration:
                break

            # i7 barcodes (get from name string instead of sequence)
            if self.data.hackers.demultiplex_on_i7_tags:
                barcode = read1[0].rsplit(":", 1)[-1].split("+")[0].strip()

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

            # find if it matches 
            sname_match = self.matchdict.get(barcode)
            if sname_match:

                # increment counters
                self.sample_to_counts.update([sname_match])
                self.barcode_counts.update([barcode])
                self.nmatched += 1

                # trim off barcode
                lenbar1 = len(barcode)
                if '3rad' in self.data.params.datatype:
                    # Iff 3rad trim the len of the first barcode
                    lenbar1 = len(barcode1)
                    lenbar2 = len(barcode2)

                # no trim on i7 demux
                if self.data.hackers.demultiplex_on_i7_tags:
                    lenbar1 = lenbar2 = 0

                # for 2brad we trim the barcode AND the synthetic overhang
                # The `+1` is because it trims the newline
                if self.data.params.datatype == '2brad':
                    overlen = len(self.cutters[0][0]) + lenbar1 + 1
                    read1[1] = read1[1][:-overlen] + "\n"
                    read1[3] = read1[3][:-overlen] + "\n"
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
                self.read1s[sname_match].append("".join(read1))
                if self.data.is_pair:
                    self.read2s[sname_match].append("".join(read2))

            else:
                if barcode:
                    self.misses.update([barcode])
                else:
                    self.misses.update("_")

            # Write to each sample file (pid's have different handles)
            # print(self.nmatched)
            if not self.nmatched % int(2e5):

                # tell logger stats
                print(f"progress:\n{self.sample_to_counts}")

                # write reads to file
                write_to_file(self.data, self.read1s, 1, self.fidx)
                if self.data.is_pair:
                    write_to_file(self.data, self.read2s, 2, self.fidx)

                # clear out lits of sorted reads
                for sname in self.barcodes:
                    # if "-technical-replicate-" in sname:
                        # sname = sname.rsplit("-technical-replicate", 1)[0]
                    self.read1s[sname] = []
                    self.read2s[sname] = []             

        print(f"progress: {self.sample_to_counts}")
        ## write the remaining reads to file
        write_to_file(self.data, self.read1s, 1, self.fidx)
        if self.data.is_pair:
            write_to_file(self.data, self.read2s, 2, self.fidx)


def getbarcode1(cutters, read1, longbar):
    "find barcode for 2bRAD data"
    #+1 is for the \n at the end of the sequence line
    lencut = len(cutters[0][0]) + 1
    return read1[1][:-lencut][-longbar[0]:]

def getbarcode2(_, read1, longbar):
    "finds barcode for invariable length barcode data"
    return read1[1][:longbar[0]]

def getbarcode3(cutters, read1, longbar):
    """find barcode sequence in the beginning of read."""
    # default barcode string
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
    """
    find barcode sequence in the beginning of read
    """
    # default barcode string
    for ambigcuts in cutters:
        # print("AMBIGCUTS", ambigcuts)            
        for cutter in ambigcuts:
            # print("cutter", cutter)            
            # If the cutter is unambiguous there will only be one.
            if not cutter:
                continue
            search = read[1][:int(longbar + len(cutter) + 1)]
            try:
                splitsearch = search.decode().rsplit(cutter, 1)
            except (AttributeError, TypeError):
                splitsearch = search.rsplit(cutter, 1)
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
        handle = data.tmpdir / f"tmp_{sname}_{rrr}_{pid}.fastq"

        # append to this sample name
        with open(handle, 'a', encoding="utf-8") as out:
            out.write("".join(dsort[sname]))
        # logger.debug("appending to {}".format(handle))
        # logger.complete()


def collate_files(data, sname, tmp1s, tmp2s):
    """ 
    Collate temp fastq files in tmp-dir into 1 gzipped sample.
    """
    # out handle
    out1 = data.stepdir / f"{sname}_R1.fastq.gz"
    out2 = data.stepdir / f"{sname}_R2.fastq.gz"

    # build cmd
    cmd1 = ['cat']
    for tmpfile in tmp1s:
        cmd1 += [tmpfile]

    # get compression function
    with sps.Popen(['which', 'pigz'], stderr=sps.PIPE, stdout=sps.PIPE) as proc:
        pigz_bin = proc.communicate()
        if pigz_bin[0].strip():
            compress = ["pigz"]
        else:
            compress = ["gzip"]

    # call cmd
    out = io.BufferedWriter(gzip.open(out1, 'w'))    
    with sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.PIPE) as proc1:
        with sps.Popen(compress, stdin=proc1.stdout, stderr=sps.PIPE, stdout=out) as proc2:
            eout, _ = proc2.communicate()
    if proc2.returncode:
        logger.exception(eout)
        raise IPyradError("error in collate_files R1 {}".format(eout))
    # proc1.stdout.close()
    out.close()

    # then cleanup
    for tmpfile in tmp1s:
        tmpfile.unlink()
        # os.remove(tmpfile)

    # do R2 files
    if data.is_pair:
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
            tmpfile.unlink()
            # os.remove(tmpfile)
