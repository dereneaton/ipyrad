#!/usr/bin/env ipython2

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

from __future__ import print_function
import os
import glob
import subprocess


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


