#!/usr/bin/env ipython2

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

from __future__ import print_function
import os
import glob
import itertools
import subprocess

import logging
LOGGER = logging.getLogger(__name__)

def ambigcutters(seq):
    """ returns both resolutions of a cut site that has an ambiguous base in 
    it, else the single cut site """
    resos = []
    ambigs = {"R":("G", "A"),
              "K":("G", "T"),
              "S":("G", "C"),
              "Y":("T", "C"),
              "W":("T", "A"),
              "M":("C", "A")}
    if any([i in list("RKSYWM") for i in seq]):
        for base in list("RKSYWM"):
            if base in seq:
                resos.append(seq.replace(base, ambigs[base][0]))
                resos.append(seq.replace(base, ambigs[base][1]))
        return resos
    else:
        return [seq, ""]

def comp(seq):
    """ returns a seq with small complement"""
    return seq.replace("A", 't')\
           .replace('T', 'a')\
           .replace('C', 'g')\
           .replace('G', 'c')\
           .replace('n', 'Z')\
           .upper()\
           .replace("Z", "n")\
           .replace("S", "s")

## TODO: Delete this? This function isn't used in the codebase
## maybe delete it.
def getoptim(filename):
    """ Calculate optimum splitting based on file size. 
    Does not unzip files, assumes average rate of compression. 
    This is a fast alternative to counting lines which takes 
    too long on huge files.
    """
    filesize = os.stat(filename).st_size
    if filesize < 160000000:
        optim = 40000
    elif filesize < 4000000000:
        optim = 800000
    elif filesize < 8000000000:
        optim = 12000000
    else:
        optim = 24000000
    return optim

def getsplits(filename):
    """ Calculate optimum splitting based on file size. 
    Does not unzip files, assumes average rate of compression. 
    This is a fast alternative to counting lines which takes 
    too long on huge files.
    """
    filesize = os.stat(filename).st_size
    if filesize < 10000000:
        optim = 400000
    elif filesize < 4000000000:
        optim = 1000000
    elif filesize < 8000000000:
        optim = 4000000
    else:
        optim = 8000000
    return optim


def merge_pairs( data, sample, unmerged_files ):
    """ Merge PE reads. Takes in a tuple of unmerged files
    and returns the file name of the merged/combined PE
    reads and the number of reads that were merged (overlapping)
    """
    LOGGER.debug("Entering merge_pairs - %s", unmerged_files)

    ## tempnames for merge files
    merged = os.path.join(data.dirs.edits,
                          sample.name+"_merged_.fastq")
    nonmerged1 = os.path.join(data.dirs.edits,
                          sample.name+"_nonmerged_R1_.fastq")
    nonmerged2 = os.path.join(data.dirs.edits,
                          sample.name+"_nonmerged_R2_.fastq")

    try:
        maxn = sum(data.paramsdict['max_low_qual_bases'])
    except TypeError:
        maxn = data.paramsdict['max_low_qual_bases']

    assert os.path.exists(unmerged_files[1]), \
           "No paired read file (_R2_ file) found." 

    ## vsearch merging
    cmd = data.bins.vsearch \
      +" --fastq_mergepairs "+unmerged_files[0] \
      +" --reverse "+unmerged_files[1] \
      +" --fastqout "+merged \
      +" --fastqout_notmerged_fwd "+nonmerged1 \
      +" --fastqout_notmerged_rev "+nonmerged2 \
      +" --fasta_width 0 " \
      +" --fastq_allowmergestagger " \
      +" --fastq_minmergelen 32 " \
      +" --fastq_maxns "+str(maxn) \
      +" --fastq_minovlen 12 "

    LOGGER.debug( cmd )
    try:
        subprocess.check_call(cmd, shell=True,
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(subprocess.STDOUT)
        LOGGER.error(cmd)
        sys.exit("Error in merging pairs: \n({}).".format(inst))
    ## record how many read pairs were merged
    with open(merged, 'r') as tmpf:
        nmerged = len(tmpf.readlines())

    LOGGER.debug( "Merged pairs - %d", nmerged )
    ## Combine the unmerged pairs and append to the merge file
    with open(merged, 'ab') as combout:
        ## read in paired end read files"
        ## create iterators to sample 4 lines at a time
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

        combout.write("\n".join(writing))

    os.remove( nonmerged1 )
    os.remove( nonmerged2 )

    return merged, nmerged

def most_common(L):
    return max(groupby(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]

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
    trans = {"R":("G", "A"),
             "K":("G", "T"),
             "S":("G", "C"),
             "Y":("T", "C"),
             "W":("T", "A"),
             "M":("C", "A")}
    return trans.get(amb)

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

def unstruct(amb):
    """ This is copied from pyrad.alignable, and is referenced in
    several of the loci2*.py conversion modules. It duplicates some
    of the effort of unhetero(), but i guess it's fine for now. Probably
    could merge these two functions if you wanted to. 
    TODO: Also could make the D dict{} a global so you wouldn't have to 
    recreate it every time this function is called. Could save some cycles.
    """
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"],
         "A":["A","A"],
         "T":["T","T"],
         "G":["G","G"],
         "C":["C","C"],
         "N":["N","N"],
         "-":["-","-"]}
    return D.get(amb)

def zcat_make_temps(args):
    """ call bash command zcat and split to split large files """
    ## split args
    data, raws, num, optim = args

    ## get optimum lines per file
    if not optim:
        optim = getsplits(raws[0])
    #LOGGER.info("optim = %s", optim)

    ## is it gzipped
    cat = "cat"
    if raws[0].endswith(".gz"):
        cat = "gunzip -c"

    ### run splitter
    cmd = " ".join([cat, raws[0], "|", "split", "-l", str(optim),
                   "-", os.path.join(data.dirs.fastqs, "chunk1_"+str(num)+"_")])
    _ = subprocess.call(cmd, shell=True,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE,
                             close_fds=True)
    chunks1 = glob.glob(os.path.join(
                        data.dirs.fastqs, "chunk1_"+str(num)+"_*"))
    chunks1.sort()

    if "pair" in data.paramsdict["datatype"]:
        cmd = " ".join([cat, raws[1], "|", "split", "-l", str(optim),
                  "-", os.path.join(data.dirs.fastqs, "chunk2_"+str(num)+"_")])
        _ = subprocess.call(cmd, shell=True,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE,
                             close_fds=True)
        chunks2 = glob.glob(os.path.join(
                        data.dirs.fastqs, "chunk2_"+str(num)+"_*"))
        chunks2.sort()
        #LOGGER.debug("chunksfiles: %s %s", chunks1, chunks2)
        assert len(chunks1) == len(chunks2), \
            "R1 and R2 files are not the same length."
    else:
        chunks2 = [0]*len(chunks1)

    return [raws[0], zip(chunks1, chunks2)]


