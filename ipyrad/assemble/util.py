#!/usr/bin/env ipython2

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

# pylint: disable=E1101
# pylint: disable=W0212
# pylint: disable=W0142

from __future__ import print_function
import os
import sys
import socket
import tempfile
import itertools
import ipyrad
import gzip
from collections import defaultdict

try:
    import subprocess32 as sps
except ImportError:
    import subprocess as sps

import logging
LOGGER = logging.getLogger(__name__)

## a subset of functions to import when importing as *
#__all__ = ["IPyradError", "IPyradParamsError", "IPyradWarningExit",
#           "ObjDict", "comp"]


### custom Exception classes
class IPyradParamsError(Exception):
    """ Exception handler indicating error in parameter entry """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class IPyradError(Exception):
    """ Exception handler indicating error in during assembly """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class IPyradWarningExit(SystemExit):
    """
    Exception handler that does clean exit for CLI, but also prints
    the traceback and cleaner message for API.
    """
    def __init__(self, *args, **kwargs):
        if ipyrad.__interactive__:
            raise IPyradError(*args, **kwargs)
        else:
            SystemExit.__init__(self, *args, **kwargs)


## utility functions/classes
class Params(object):
    """ 
    A dict-like object for storing params values with a custom repr
    that shortens file paths, and which makes attributes easily viewable
    through tab completion in a notebook while hiding other funcs, attrs, that
    are in normal dicts. 
    """
    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())
        if keys:
            _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
            for key in keys:
                _val = str(self[key]).replace(os.path.expanduser("~"), "~")
                _repr += _printstr.format(key, _val)
        return _repr




class ObjDict(dict):
    """
    Object dictionary allows calling dictionaries in a more
    pretty and Python fashion for storing Assembly data
    """
    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __repr__(self):
        result = ""
        if "outfiles" in self.keys():
            dirs_order = ["fastqs", "edits", "clusts", "consens", "outfiles"]
            for key in dirs_order:
                result += key + " : " + self[key] + "\n"
        else:
            for key in sorted(self):
                result += key + " : " + str(self[key]) + "\n"
        return result


### A decorator for use in step5 base calling
def memoize(func):
    """ Memoization decorator for a function taking one or more arguments. """
    class Memodict(dict):
        """ just a dict"""
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            """ this makes it faster """
            ret = self[key] = func(*key)
            return ret

    return Memodict().__getitem__




CDICT = {i:j for i, j in zip("CATG", "0123")}


## used for geno output
VIEW = {"R":("G", "A"),
        "K":("G", "T"),
        "S":("G", "C"),
        "Y":("T", "C"),
        "W":("T", "A"),
        "M":("C", "A"),
        "A":("X", "X"),
        "T":("X", "X"),
        "G":("X", "X"),
        "C":("X", "X"),
        "N":("X", "X"),
        "-":("X", "X"),
        }

## used in hetero() func of consens_se.py
TRANS = {('G', 'A'):"R",
         ('G', 'T'):"K",
         ('G', 'C'):"S",
         ('T', 'C'):"Y",
         ('T', 'A'):"W",
         ('C', 'A'):"M"}

TRANSFULL = {
         ('G', 'A'):"R",
         ('G', 'T'):"K",
         ('G', 'C'):"S",
         ('T', 'C'):"Y",
         ('T', 'A'):"W",
         ('C', 'A'):"M",
         ('A', 'C'):"M",
         ('A', 'T'):"W",
         ('C', 'T'):"Y",
         ('C', 'G'):"S",
         ('T', 'G'):"K",
         ('A', 'G'):"R"}


## used for resolving ambiguities
AMBIGS = {"R":("G", "A"),
          "K":("G", "T"),
          "S":("G", "C"),
          "Y":("T", "C"),
          "W":("T", "A"),
          "M":("C", "A")}



def ambigcutters(seq):
    """
    Returns both resolutions of a cut site that has an ambiguous base in
    it, else the single cut site
    """
    resos = []
    if any([i in list("RKSYWM") for i in seq]):
        for base in list("RKSYWM"):
            if base in seq:
                resos.append(seq.replace(base, AMBIGS[base][0]))
                resos.append(seq.replace(base, AMBIGS[base][1]))
        return resos
    else:
        return [seq, ""]



def splitalleles(consensus):
    """ takes diploid consensus alleles with phase data stored as a mixture
    of upper and lower case characters and splits it into 2 alleles """

    ## store two alleles, allele1 will start with bigbase
    allele1 = list(consensus)
    allele2 = list(consensus)
    hidx = [i for (i, j) in enumerate(consensus) if j in "RKSWYMrkswym"]

    ## do remaining h sites
    for idx in hidx:
        hsite = consensus[idx]
        if hsite.isupper():
            allele1[idx] = PRIORITY[hsite]
            allele2[idx] = MINOR[hsite]
        else:
            allele1[idx] = MINOR[hsite.upper()]
            allele2[idx] = PRIORITY[hsite.upper()]

    ## convert back to strings
    allele1 = "".join(allele1)
    allele2 = "".join(allele2)

    return allele1, allele2



def comp(seq):
    """ returns a seq with complement. Preserves little n's for splitters."""
    ## makes base to its small complement then makes upper
    return seq.replace("A", 't')\
              .replace('T', 'a')\
              .replace('C', 'g')\
              .replace('G', 'c')\
              .replace('n', 'Z')\
              .upper()\
              .replace("Z", "n")



def fullcomp(seq):
    """ returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences"""
    ## this is surely not the most efficient...
    seq = seq.replace("A", 'u')\
             .replace('T', 'v')\
             .replace('C', 'p')\
             .replace('G', 'z')\
             .replace('u', 'T')\
             .replace('v', 'A')\
             .replace('p', 'G')\
             .replace('z', 'C')

    ## No complement for S & W b/c complements are S & W, respectively
    seq = seq.replace('R', 'u')\
             .replace('K', 'v')\
             .replace('Y', 'b')\
             .replace('M', 'o')\
             .replace('u', 'Y')\
             .replace('v', 'M')\
             .replace('b', 'R')\
             .replace('o', 'K')

    seq = seq.replace('r', 'u')\
             .replace('k', 'v')\
             .replace('y', 'b')\
             .replace('m', 'o')\
             .replace('u', 'y')\
             .replace('v', 'm')\
             .replace('b', 'r')\
             .replace('o', 'k')
    return seq



def fastq_touchup_for_vsearch_merge(read, outfile, reverse=False):
    """ option to change orientation of reads and sets Qscore to B """
    
    counts = 0
    with open(outfile, 'w') as out:
        ## read in paired end read files 4 lines at a time
        if read.endswith(".gz"):
            fr1 = gzip.open(read, 'rb')
        else:
            fr1 = open(read, 'rb')
        quarts = itertools.izip(*[iter(fr1)]*4)

        ## a list to store until writing
        writing = []

        while 1:
            try:
                lines = quarts.next()
            except StopIteration:
                break
            if reverse:
                seq = lines[1].strip()[::-1]
            else:
                seq = lines[1].strip()
            writing.append("".join([
                lines[0],
                seq+"\n",
                lines[2],
                "B"*len(seq)
            ]))

            ## write to disk
            counts += 1
            if not counts % 1000:
                out.write("\n".join(writing)+"\n")
                writing = []
        if writing:
            out.write("\n".join(writing))
            
    out.close()
    fr1.close()
                               


def merge_pairs_after_refmapping(data, two_files, merged_out):
    """ 
    A function to merge fastq files produced by bam2fq.
    """

    ## create temp files 
    nonmerged1 = tempfile.NamedTemporaryFile(
        mode='wb',
        dir=data.dirs.edits,
        suffix="_nonmerged_R1_.fastq").name
    nonmerged2 = tempfile.NamedTemporaryFile(
        mode='wb',
        dir=data.dirs.edits,
        suffix="_nonmerged_R2_.fastq").name 

    ## get the maxn and minlen values
    minlen = str(max(32, data.paramsdict["filter_min_trim_len"]))
    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']
 
    ## set the quality scores abritrarily high and orient R2 correctly
    tmp1 = two_files[0][0]
    tmp2 = two_files[0][1]
    fastq_touchup_for_vsearch_merge(tmp1, tmp1+".tu", False)
    fastq_touchup_for_vsearch_merge(tmp2, tmp2+".tu", True)

    ## command string to call vsearch
    cmd = [ipyrad.bins.vsearch,
           "--fastq_mergepairs", tmp1+".tu",
           "--reverse", tmp2+".tu",
           "--fastqout", merged_out,
           "--fastqout_notmerged_fwd", nonmerged1,
           "--fastqout_notmerged_rev", nonmerged2,
           "--fasta_width", "0",
           "--fastq_minmergelen", minlen,
           "--fastq_maxns", str(maxn),
           "--fastq_minovlen", "10",
           "--fastq_maxdiffs", "4",
           "--label_suffix", "_m1",
           "--fastq_qmax", "1000",
           "--threads", "2",
           "--fastq_allowmergestagger"]
    
    ## run vsearch but allow kbd
    proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    try:
        res = proc.communicate()[0]
    except KeyboardInterrupt:
        proc.kill()

    ## cleanup tmp files if job failed or stopped    
    if proc.returncode:
        LOGGER.error("Error: %s %s", cmd, res)
        raise IPyradWarningExit("Error merge pairs:\n %s\n%s", cmd, res)

    ## record how many read pairs were merged
    with open(merged_out, 'r') as tmpf:
        nmerged = sum(1 for i in tmpf.readlines()) // 4

    ## Concat unmerged pairs with a 'nnnn' separator
    with open(merged_out, 'ab') as combout:
        ## read in paired end read files 4 lines at a time
        fr1 = open(nonmerged1, 'rb')
        quart1 = itertools.izip(*[iter(fr1)]*4)
        fr2 = open(nonmerged2, 'rb')
        quart2 = itertools.izip(*[iter(fr2)]*4)
        quarts = itertools.izip(quart1, quart2)

        ## a list to store until writing
        writing = []
        counts = 0

        ## iterate until done
        while 1:
            try:
                read1s, read2s = quarts.next()
            except StopIteration:
                break

            ## store the read
            writing.append("".join([
                read1s[0],
                read1s[1].strip() + "nnnn" + \
                read2s[1], 
                read1s[2],
                read1s[3].strip() + "nnnn" + \
                read2s[3], 
                ]))

            ## count up until time to write
            counts += 1
            if not counts % 10:
                combout.write("".join(writing))
                writing = []
        ## write the remaining
        if writing:
            combout.write("".join(writing))

        ## close handles
        fr1.close()
        fr2.close()
        combout.close()

    ## remove temp files (or do this later)
    rmfiles = [nonmerged1, nonmerged2, tmp1, tmp2, tmp1+".tu", tmp2+".tu"]
    for rmfile in rmfiles:
        if os.path.exists(rmfile):
            os.remove(rmfile)

    return nmerged



def merge_after_pysam(data, clust):
    """
    This is for pysam post-flight merging. The input is a cluster
    for an individual locus. We have to split the clusters, write
    R1 and R2 to files then call merge_pairs(). This is not ideal,
    it's slow, but it works. This is the absolute worst way to do this,
    it bounces all the files for each locus off the disk. I/O _hog_.
    """
    try:
        r1file = tempfile.NamedTemporaryFile(mode='wb', delete=False,
                                            dir=data.dirs.edits,
                                            suffix="_R1_.fastq")
        r2file = tempfile.NamedTemporaryFile(mode='wb', delete=False,
                                            dir=data.dirs.edits,
                                            suffix="_R2_.fastq")

        r1dat = []
        r2dat = []
        for locus in clust:
            sname, seq = locus.split("\n")
            ## Have to munge the sname to make it look like fastq format
            sname = "@" + sname[1:]
            r1, r2 = seq.split("nnnn")
            r1dat.append("{}\n{}\n{}\n{}".format(sname, r1, "+", "B"*(len(r1))))
            r2dat.append("{}\n{}\n{}\n{}".format(sname, r2, "+", "B"*(len(r2))))

        r1file.write("\n".join(r1dat))
        r2file.write("\n".join(r2dat))
        r1file.close()
        r2file.close()

        ## Read in the merged data and format to return as a clust
        merged_file = tempfile.NamedTemporaryFile(mode='wb',
                                            dir=data.dirs.edits,
                                            suffix="_merged.fastq").name

        clust = []
        merge_pairs(data, [(r1file.name, r2file.name)], merged_file, 0, 1)

        with open(merged_file) as infile:
            quarts = itertools.izip(*[iter(infile)]*4)

            while 1:
                try:
                    sname, seq, _, _ = quarts.next()
                    ## Vsearch expects R2 oriented how it would be in a raw data file
                    ## i.e. revcomp, and that's also how it returns it
                    ## but we want to maintain the genomic orientation so R1 and
                    ## R2 are both on the + strand and both in ascending positional order
                    ## so here if the reads don't merge we have to revcomp R2
                    if not "_m1" in sname.rsplit(";", 1)[1]:
                        try:
                            R1, R2 = seq.split("nnnn")
                            seq = R1 + "nnnn" + revcomp(R2)
                        except ValueError as inst:
                            LOGGER.error("Failed merge_after_pysam: {} {}".format(sname, seq))
                            raise

                except StopIteration:
                    break
                ## put sname back
                sname = ">" + sname[1:]
                clust.extend([sname.strip(), seq.strip()])
    except:
        LOGGER.info("Error in merge_pairs post-refmap.")
        raise
    finally:
        for i in [r1file.name, r2file.name, merged_file]:
            if os.path.exists(i):
                log_level = logging.getLevelName(LOGGER.getEffectiveLevel())
                ## uncomment this to shit ALL over the filesystem. 
                ## if not log_level == "DEBUG":
                os.remove(i)

    return clust



def merge_pairs(data, two_files, merged_out, revcomp, merge):
    """
    Merge PE reads. Takes in a list of unmerged files [r1, r2] and the
    filehandle to write merged data to, and it returns the number of reads
    that were merged (overlapping). If merge==0 then only concat pairs (nnnn),
    no merging in vsearch.

    Parameters
    -----------
    two_files (tuple):
        A list or tuple of the [r1, r2] files to be merged. 

    merged_out (str):
        A string file handle for the merged data to be written to. 

    revcomp (bool):
        Whether or not to revcomp the R2s. 

    merge (bool):
        Whether or not to perform vsearch merging. If not then reads are simply
        concatenated with a 'nnnn' separator. 

    Returns
    --------
    If merge is on then the func will return the number of pairs 
    successfully merged, else it returns -1. 
    """

    LOGGER.debug("Entering merge_pairs()")

    ## Return the number of merged pairs
    nmerged = -1

    ## Check input files from inside list-tuple [(r1, r2)]
    for fhandle in two_files[0]:
        if not os.path.exists(fhandle):
            raise IPyradWarningExit("""
    Attempting to merge a file that doesn't exist - {}""".format(fhandle))

    ## If it already exists, clean up the old merged file
    if os.path.exists(merged_out):
        os.remove(merged_out)

    ## if merge then catch nonmerged in a separate file
    if merge:
        nonmerged1 = tempfile.NamedTemporaryFile(mode='wb',
                                            dir=data.dirs.edits,
                                            suffix="_nonmerged_R1_.fastq").name
        nonmerged2 = tempfile.NamedTemporaryFile(mode='wb',
                                            dir=data.dirs.edits,
                                            suffix="_nonmerged_R2_.fastq").name

    ## if not merging then the nonmerged reads will come from the normal edits
    else:
        nonmerged1 = two_files[0][0]
        nonmerged2 = two_files[0][1]

    ## get the maxn and minlen values
    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']
    minlen = str(max(32, data.paramsdict["filter_min_trim_len"]))

    ## we need to gunzip the files if they are zipped (at least for now)
    if merge and two_files[0][0].endswith(".gz"):
        LOGGER.info("gunzipping pairs")
        tmp1 = os.path.splitext(two_files[0][0])[0]+".tmp1"
        tmp2 = os.path.splitext(two_files[0][1])[0]+".tmp2"

        out1 = open(tmp1, 'w')
        out2 = open(tmp2, 'w')
        gun1 = sps.Popen(["gunzip", "-c", two_files[0][0]],
                          stderr=sps.STDOUT, stdout=out1, close_fds=True)
        gun2 = sps.Popen(["gunzip", "-c", two_files[0][1]],
                          stderr=sps.STDOUT, stdout=out2, close_fds=True)
        _ = gun1.communicate()
        _ = gun2.communicate()
        out1.close()
        out2.close()
    else:
        tmp1 = two_files[0][0]
        tmp2 = two_files[0][1]
    try:
        ## If we are actually mergeing and not just joining then do vsearch
        if merge:
            ## create tmp files with high quality scores and with R2 oriented
            cmd = [ipyrad.bins.vsearch,
                   "--fastq_mergepairs", tmp1,
                   "--reverse", tmp2,
                   "--fastqout", merged_out,
                   "--fastqout_notmerged_fwd", nonmerged1,
                   "--fastqout_notmerged_rev", nonmerged2,
                   "--fasta_width", "0",
                   "--fastq_minmergelen", minlen,
                   "--fastq_maxns", str(maxn),
                   "--fastq_minovlen", "20",
                   "--fastq_maxdiffs", "4",
                   "--label_suffix", "_m1",
                   "--fastq_qmax", "1000",
                   "--threads", "2",
                   "--fastq_allowmergestagger"]

            LOGGER.debug("merge cmd: %s", " ".join(cmd))
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            try:
                res = proc.communicate()[0]
            except KeyboardInterrupt:
                proc.kill()

            if proc.returncode:
                LOGGER.error("Error: %s %s", cmd, res)
                ## remove temp files
                rmfiles = [os.path.splitext(two_files[0][0])[0]+".tmp1",
                           os.path.splitext(two_files[0][1])[0]+".tmp2",
                           nonmerged1, nonmerged2]
                for rmfile in rmfiles:
                    if os.path.exists(rmfile):
                        os.remove(rmfile)
                raise IPyradWarningExit("Error merge pairs:\n %s\n%s", cmd, res)

            ## record how many read pairs were merged
            with open(merged_out, 'r') as tmpf:
                #nmerged = len(tmpf.readlines()) // 4
                nmerged = sum(1 for i in tmpf.readlines()) // 4

        ## Combine the unmerged pairs and append to the merge file
        with open(merged_out, 'ab') as combout:
            ## read in paired end read files 4 lines at a time
            if nonmerged1.endswith(".gz"):
                fr1 = gzip.open(nonmerged1, 'rb')
            else:
                fr1 = open(nonmerged1, 'rb')
            quart1 = itertools.izip(*[iter(fr1)]*4)
            if nonmerged2.endswith(".gz"):
                fr2 = gzip.open(nonmerged2, 'rb')
            else:
                fr2 = open(nonmerged2, 'rb')
            quart2 = itertools.izip(*[iter(fr2)]*4)
            quarts = itertools.izip(quart1, quart2)

            ## a list to store until writing
            writing = []
            counts = 0

            ## iterate until done
            while 1:
                try:
                    read1s, read2s = quarts.next()
                except StopIteration:
                    break
                if revcomp:
                    writing.append("".join([
                        read1s[0],
                        read1s[1].strip() + "nnnn" + \
                        comp(read2s[1].strip()[::-1]) + "\n",
                        read1s[2],
                        read1s[3].strip() + "nnnn" + \
                        read2s[3].strip()[::-1] + "\n",
                        ]))
                else:
                    writing.append("".join([
                        read1s[0],
                        read1s[1].strip() + "nnnn" + \
                        read2s[1],
                        read1s[2],
                        read1s[3].strip() + "nnnn" + \
                        read2s[3],
                        ]))

                counts += 1
                if not counts % 10:
                    combout.write("".join(writing)) #+"\n")
                    writing = []

            if writing:
                combout.write("".join(writing))

        ## close handles
        fr1.close()
        fr2.close()
        combout.close()

    except Exception as inst:
        LOGGER.error("Exception in merge_pairs - {}".format(inst))
        raise
    ## No matter what happens please clean up the temp files.
    finally:
        ## if merged then delete the nonmerge tmp files
        if merge:
            ## remove temp files
            rmfiles = [nonmerged1, nonmerged2,
                       os.path.splitext(two_files[0][0])[0]+".tmp1",
                       os.path.splitext(two_files[0][1])[0]+".tmp2"]
            for rmfile in rmfiles:
                if os.path.exists(rmfile):
                    os.remove(rmfile)

    return nmerged


## Doesn't work right now.
def merge_pair_pipes(data, sample, clust):
    r1_inpipe_name = os.path.join(data.tmpdir, 'r1_in')
    r2_inpipe_name = os.path.join(data.tmpdir, 'r2_in')
    unmapped1 = tempfile.NamedTemporaryFile()
    unmapped2 = tempfile.NamedTemporaryFile()
    merged = tempfile.NamedTemporaryFile()
    pipelist = [r1_inpipe_name, r2_inpipe_name]
    try:
        r1reads, r2reads = map(lambda x: x.split("nnnn"), clust)
    except:
        LOGGER.debug("bad clust")
        raise

    def childR1( ):
        LOGGER.debug("Entering c1 {} {}".format(data, sample))
        pipeout = os.open(r1_inpipe_name, os.O_WRONLY)
        for line in r1reads:
            os.write(pipeout, line)
        os.close(pipeout)

    def childR2( ):
        LOGGER.debug("Entering c2 {} {}".format(data, sample))
        pipeout = os.open(r2_inpipe_name, os.O_WRONLY)
        for line in r2reads:
            os.write(pipeout, line)
        os.close(pipeout)

    def vsearch( ):
        LOGGER.debug("entering vsearch")
        cmd = ["/Users/iovercast/miniconda2/bin/vsearch",
               "--fastq_mergepairs", r1_inpipe_name,
               "--reverse", r2_inpipe_name,
               "--fastqout", merged.name,
               "--fastqout_notmerged_fwd", unmapped1,
               "--fastqout_notmerged_rev", unmapped2,
               "--fasta_width", "0",
               "--fastq_minmergelen", "32",
               "--fastq_maxns", "6",
               "--fastq_minovlen", "10",
               "--fastq_maxdiffs", "4",
               "--label_suffix", "_m1",
               "--fastq_qmax", "1000",
               "--threads", "2",
               "--fastq_allowmergestagger"]
        proc = subprocess.Popen(cmd)
        ret = proc.communicate()

    joblist = [childR1, childR2]
    for pipe in pipelist:
        if os.path.exists(pipe):
            os.unlink(pipe)
        if not os.path.exists(pipe):
            os.mkfifo(pipe)

    children = []
    isparent = True
    for job in joblist:
        child = os.fork()
        if child:
            children.append(child)
        else:
            isparent = False
            job()
            os._exit(0)

    if isparent:
        for child in children:
            os.wait3(child)

    with open(merged, 'r'):
        clusts = merged.read()
        LOGGER.error("lenmerged = {}".format(len(clusts)))
        return clusts


def revcomp(sequence):
    "returns reverse complement of a string"
    sequence = sequence[::-1].strip()\
                             .replace("A", "t")\
                             .replace("T", "a")\
                             .replace("C", "g")\
                             .replace("G", "c").upper()
    return sequence



def unhetero(amb):
    " returns bases from ambiguity code"
    amb = amb.upper()
    return AMBIGS.get(amb)




## Alleles priority dict. The key:vals are the same as the AMBIGS dict
## except it returns just one base, w/ the order/priority being (C>A>T>G)
## This dict is used to impute lower case into consens to retain allele
## order for phase in diploids
PRIORITY = {"M": "C",
            "Y": "C",
            "S": "C",
            "W": "A",
            "R": "A",
            "K": "T"}

## The inverse of priority
MINOR = {"M": "A",
         "Y": "T",
         "S": "G",
         "W": "T",
         "R": "G",
         "K": "G"}


DUCT = {"R":["G", "A"],
        "K":["G", "T"],
        "S":["G", "C"],
        "Y":["T", "C"],
        "W":["T", "A"],
        "M":["C", "A"],
        "A":["A", "A"],
        "T":["T", "T"],
        "G":["G", "G"],
        "C":["C", "C"],
        "N":["N", "N"],
        "-":["-", "-"]
        }


def unstruct(amb):
    """ This is copied from pyrad.alignable, and is referenced in
    several of the loci2*.py conversion modules. It duplicates some
    of the effort of unhetero(), but i guess it's fine for now. Probably
    could merge these two functions if you wanted to.
    """
    amb = amb.upper()
    ## returns bases from ambiguity code"
    return DUCT.get(amb)




def clustdealer(pairdealer, optim):
    """ return optim clusters given iterators, and whether it got all or not"""
    ccnt = 0
    chunk = []
    while ccnt < optim:
        ## try refreshing taker, else quit
        try:
            taker = itertools.takewhile(lambda x: x[0] != "//\n", pairdealer)
            oneclust = ["".join(taker.next())]
        except StopIteration:
            #LOGGER.debug('last chunk %s', chunk)
            return 1, chunk

        ## load one cluster
        while 1:
            try:
                oneclust.append("".join(taker.next()))
            except StopIteration:
                break
        chunk.append("".join(oneclust))
        ccnt += 1
    return 0, chunk




def progressbar(njobs, finished, msg="", spacer="  "):
    """ prints a progress bar """
    if njobs:
        progress = 100*(finished / float(njobs))
    else:
        progress = 100
        
    hashes = '#'*int(progress/5.)
    nohash = ' '*int(20-len(hashes))
    if not ipyrad.__interactive__:
        msg = msg.rsplit("|", 2)[0]

    args = [spacer, hashes+nohash, int(progress), msg]
    print("\r{}[{}] {:>3}% {} ".format(*args), end="")
    sys.stdout.flush()




def get_threaded_view(ipyclient, split=True):
    """ gets optimum threaded view of ids given the host setup """
    ## engine ids
    ## e.g., [0, 1, 2, 3, 4, 5, 6, 7, 8]
    eids = ipyclient.ids

    ## get host names
    ## e.g., ['a', 'a', 'b', 'b', 'a', 'c', 'c', 'c', 'c']
    dview = ipyclient.direct_view()
    hosts = dview.apply_sync(socket.gethostname)

    ## group ids into a dict by their hostnames
    ## e.g., {a: [0, 1, 4], b: [2, 3], c: [5, 6, 7, 8]}
    hostdict = defaultdict(list)
    for host, eid in zip(hosts, eids):
        hostdict[host].append(eid)

    ## Now split threads on the same host into separate proc if there are many
    hostdictkeys = hostdict.keys()
    for key in hostdictkeys:
        gids = hostdict[key]
        maxt = len(gids)
        if len(gids) >= 4:
            maxt = 2
        ## if 4 nodes and 4 ppn, put one sample per host
        if (len(gids) == 4) and (len(hosts) >= 4):
            maxt = 4
        if len(gids) >= 6:
            maxt = 3
        if len(gids) >= 8:
            maxt = 4
        if len(gids) >= 16:
            maxt = 4
        ## split ids into groups of maxt
        threaded = [gids[i:i+maxt] for i in xrange(0, len(gids), maxt)]
        lth = len(threaded)
        ## if anything was split (lth>1) update hostdict with new proc
        if lth > 1:
            hostdict.pop(key)
            for hostid in range(lth):
                hostdict[str(key)+"_"+str(hostid)] = threaded[hostid]

    ## make sure split numbering is correct
    #threaded = hostdict.values()
    #assert len(ipyclient.ids) <= len(list(itertools.chain(*threaded)))
    LOGGER.info("threaded_view: %s", dict(hostdict))
    return hostdict




##############################################################
def detect_cpus():
    """
    Detects the number of CPUs on a system. This is better than asking
    ipyparallel since ipp has to wait for Engines to spin up.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    return 1 # Default
