#!/usr/bin/env ipython2

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

from __future__ import print_function
import os
import sys
import glob
import gzip
import tempfile
import itertools
import subprocess

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

class IPyradWarningExit(Exception):
    """ Exception handler indicating error in during assembly """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class ObjDict(dict):
    """ object dictionary allows calling dictionaries in a more 
    pretty and Python fashion for storing Assembly data """
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


AMBIGS = {"R":("G", "A"),
          "K":("G", "T"),
          "S":("G", "C"),
          "Y":("T", "C"),
          "W":("T", "A"),
          "M":("C", "A")}
def ambigcutters(seq):
    """ Returns both resolutions of a cut site that has an ambiguous base in 
    it, else the single cut site """
    resos = []
    if any([i in list("RKSYWM") for i in seq]):
        for base in list("RKSYWM"):
            if base in seq:
                resos.append(seq.replace(base, AMBIGS[base][0]))
                resos.append(seq.replace(base, AMBIGS[base][1]))
        return resos
    else:
        return [seq, ""]


def breakalleles(consensus):
    """ break ambiguity code consensus seqs into two alleles """
    allele1 = ""
    allele2 = ""
    bigbase = ""
    for base in consensus:
        if base in tuple("RKSYWM"):
            unhet1, unhet2 = unhetero(base)
            hetset = set([unhet1, unhet2])
            allele1 += uplow((unhet1, unhet2))
            allele2 += hetset.difference(uplow((unhet1, unhet2))).pop()
            if not bigbase:  
                bigbase = uplow((allele1, allele2))
        elif base in tuple("rksywm"):
            unhet1, unhet2 = unhetero(base)
            hetset = set([unhet1, unhet2])
            allele2 += uplow((unhet1, unhet2))
            allele1 += hetset.difference(uplow((unhet1, unhet2))).pop()
        else:
            allele1 += base
            allele2 += base
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



def getsplits(filename):
    """ Calculate optimum splitting based on file size. Does not unzip files, 
    assumes average rate of compression. This is a fast alternative to counting 
    lines which takes too long on huge files.
    """
    filesize = os.stat(filename).st_size
    if filesize < 10000000:
        optim = 200000
    elif filesize < 4000000000:
        optim = 500000
    elif filesize < 8000000000:
        optim = 4000000
    else:
        optim = 8000000
    return optim




def merge_pairs(data, sample): #, unmerged_files):
    """ 
    Merge PE reads. Takes in a tuple of unmerged files and returns the file 
    name of the merged/combined PE reads and the number of reads that were 
    merged (overlapping)
    """
    #LOGGER.debug("Entering merge_pairs - %s", unmerged_files)

    ## tempnames for merge files
    sample.files.merged = os.path.join(data.dirs.edits,
                                       sample.name+"_merged_.fastq")
    sample.files.nonmerged1 = os.path.join(data.dirs.edits,
                                           sample.name+"_nonmerged_R1_.fastq")
    sample.files.nonmerged2 = os.path.join(data.dirs.edits,
                                           sample.name+"_nonmerged_R2_.fastq")
    sample.files.revcomp = os.path.join(data.dirs.edits,
                                        sample.name+"_revcomp_R2_.fastq")

    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']
    minlen = str(max(32, data.paramsdict["filter_min_trim_len"]))

    # unmerged_files[1])
    assert os.path.exists(sample.files.edits[0][1]), \
           "No paired read file (_R2_ file) found." 

    ## make revcomp file
    cmd = data.bins.vsearch \
      + " --fastx_revcomp "+sample.files.edits[0][1] \
      + " --fastqout "+sample.files.revcomp
    LOGGER.warning(cmd)
    LOGGER.debug(cmd)    
    try:
        subprocess.check_call(cmd, shell=True, 
                                   stderr=subprocess.STDOUT, 
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(subprocess.STDOUT)
        LOGGER.error(cmd)
        raise SystemExit("Error in revcomping: \n ({})".format(inst))

    ## vsearch merging
    cmd = data.bins.vsearch \
      +" --fastq_mergepairs "+sample.files.edits[0][0] \
      +" --reverse "+sample.files.revcomp \
      +" --fastqout "+sample.files.merged \
      +" --fastqout_notmerged_fwd "+sample.files.nonmerged1 \
      +" --fastqout_notmerged_rev "+sample.files.nonmerged2 \
      +" --fasta_width 0 " \
      +" --fastq_allowmergestagger " \
      +" --fastq_minmergelen "+minlen \
      +" --fastq_maxns "+str(maxn) \
      +" --fastq_minovlen 12 " \
      +" --fastq_maxdiffs 4 "

    LOGGER.warning(cmd)
    try:
        subprocess.check_call(cmd, shell=True,
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error("Error in merging pairs: \n({}).".format(inst))
        LOGGER.error(subprocess.STDOUT)
        LOGGER.error(cmd)
        sys.exit("Error in merging pairs: \n({}).".format(inst))

    ## record how many read pairs were merged
    with open(sample.files.merged, 'r') as tmpf:
        nmerged = len(tmpf.readlines())

    LOGGER.debug("Merged pairs - %d", nmerged)

    ## Combine the unmerged pairs and append to the merge file
    with open(sample.files.merged, 'ab') as combout:
        ## read in paired end read files"
        ## create iterators to sample 4 lines at a time
        fr1 = open(sample.files.nonmerged1, 'rb')
        quart1 = itertools.izip(*[iter(fr1)]*4)
        fr2 = open(sample.files.nonmerged2, 'rb')
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
            writing.append("\n".join([
                            read1s[0].strip(),
                            read1s[1].strip()+\
                                "ssss"+comp(read2s[1].strip())[::-1],
                            read1s[2].strip(),
                            read1s[3].strip()+\
                                "ssss"+read2s[3].strip()[::-1]]
                            ))
            counts += 1
            if not counts % 1000:
                combout.write("\n".join(writing)+"\n")
                writing = []
        if writing:
            combout.write("\n".join(writing))
            combout.close()

    os.remove(sample.files.nonmerged1)
    os.remove(sample.files.nonmerged2)

    return sample.files.merged, nmerged


## This is hold-over code from pyrad V3 alignable, it's only used
## by loci2vcf so you could move it there if you like
def most_common(L):
    return max(itertools.groupby(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]


def preview_truncate_fq( data, sample_fastq ):
    """ If we are running in preview mode, truncate the input fq.gz file
    so it'll run quicker, just so we can see if it works. Input is tuple of the file
    names of the sample fq. Function returns a list of one tuple of 1 or 2 elements
    depending on whether you're doing paired or single end. The elements are
    paths to a temp files of the sample fq truncated to some much smaller size.
    """

    ## Return a list of filenames
    truncated_fq = []

    for read in sample_fastq[0]:

        ## If the R2 is empty then exit the loop
        if not read:
            continue
        try:
            if read.endswith(".gz"):
                f = gzip.open(os.path.realpath(read), 'rb')
            else:
                f = open(os.path.realpath(read), 'rb')

            ## create iterators to sample 4 lines at a time 
            quart = itertools.izip(*[iter(f)]*4)

            with tempfile.NamedTemporaryFile( 'w+b', delete=False,
                    dir=os.path.realpath(data.dirs.working),
                    prefix=read+".preview_tmp", suffix=".fq") as tmp_fq:
                try:
                    ## Sample the first x number of reads. On real data 2e6 is enough.  
                    for i in range(data._hackersonly["preview_truncate_length"]):
                        tmp_fq.write( "".join( quart.next() ) )
                except StopIteration:
                    LOGGER.info("preview_truncate_length > size of sample, means "+\
                                "your sample is small, nbd")
                except Exception as e:
                    LOGGER.info("preview truncate length, caught exception {}".format(e))
                    raise
            truncated_fq.append( tmp_fq.name )
            f.close()
        except AttributeError as e:
            ## R2 during SE is passed out as 0
            #truncated_fq.append( 0 )
            pass

    return [tuple(truncated_fq)]


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



def uplow(hsite):
    """ allele precedence used in assigning upper and lower case letters to 
    a consensus sequence to store the the phased allele pattern for diploids. 
    G > T > C > A """
    prec = {('G', 'A'):"G",
            ('A', 'G'):"G",
            ('G', 'T'):"G",
            ('T', 'G'):"G",
            ('G', 'C'):"G",
            ('C', 'G'):"G",
            ('T', 'C'):"T",
            ('C', 'T'):"T",
            ('T', 'A'):"T",
            ('A', 'T'):"T",
            ('C', 'A'):"C",
            ('A', 'C'):"C"}
    bigbase = prec.get(hsite)
    if not bigbase:
        bigbase = hsite[0]
    return bigbase


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
        "-":["-", "-"]}
def unstruct(amb):
    """ This is copied from pyrad.alignable, and is referenced in
    several of the loci2*.py conversion modules. It duplicates some
    of the effort of unhetero(), but i guess it's fine for now. Probably
    could merge these two functions if you wanted to. 
    """
    amb = amb.upper()
    ## returns bases from ambiguity code"
    return DUCT.get(amb)



def zcat_make_temps(args):
    """ 
    Call bash command 'zcat' and 'split' to split large files. The goal
    is to create N splitfiles where N is a multiple of the number of processors
    so that each processor can work on a file in parallel.
    """

    ## split args
    data, raws, num, tmpdir, optim = args
    LOGGER.debug("zcat splittin' %s", os.path.split(raws[0])[-1])

    ## get optimum lines per file
    if not optim:
        optim = getsplits(raws[0])
    LOGGER.info("optim = %s", optim)

    ## is it gzipped
    cat = "cat"
    if raws[0].endswith(".gz"):
        cat = "gunzip -c"

    ### run splitter
    cmd = " ".join([cat, raws[0], "|", "split", "-l", str(optim),
                   "-", os.path.join(tmpdir, "chunk1_"+str(num)+"_")])
    _ = subprocess.check_call(cmd, shell=True)

    chunks1 = glob.glob(os.path.join(tmpdir, "chunk1_"+str(num)+"_*"))
    chunks1.sort()

    if "pair" in data.paramsdict["datatype"]:
        cmd = " ".join([cat, raws[1], "|", "split", "-l", str(optim),
                  "-", os.path.join(tmpdir, "chunk2_"+str(num)+"_")])
        _ = subprocess.check_call(cmd, shell=True)

        chunks2 = glob.glob(os.path.join(tmpdir, "chunk2_"+str(num)+"_*"))
        chunks2.sort()
    
    else:
        chunks2 = [0]*len(chunks1)

    assert len(chunks1) == len(chunks2), \
        "R1 and R2 files are not the same length."

    return [raws[0], zip(chunks1, chunks2)]


