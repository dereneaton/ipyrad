#!/usr/bin/env ipython2

""" Various sequence manipulation util functions used by different
parts of the pipeline
"""

# pylint: disable=E1101
# pylint: disable=W0212

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

class IPyradWarningExit(SystemExit):
    """ Exception handler indicating error in during assembly """
    def __init__(self, *args, **kwargs):
        SystemExit.__init__(self, *args, **kwargs)


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



def splitalleles(consensus):
    """ takes diploid consensus alleles with phase data stored as a mixture 
    of upper and lower case characters and splits it into 2 alleles """

    ## store two alleles, allele1 will start with bigbase
    allele1 = list(consensus)
    allele2 = list(consensus)
    hidx = [i for (i, j) in enumerate(consensus) if j in list("RKSWYM")]

    ## do remaining h sites
    for idx in hidx:
        hsite = consensus[idx]
        if hsite.isupper:
            allele1[idx] = PRIORITY.get(hsite)
            allele2[idx] = MINOR.get(hsite)
        else:
            allele1[idx] = MINOR.get(hsite)
            allele2[idx] = PRIORITY.get(hsite)

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
      +" --reverse "+sample.files.edits[0][1] \
      +" --fastqout "+sample.files.merged \
      +" --fastqout_notmerged_fwd "+sample.files.nonmerged1 \
      +" --fastqout_notmerged_rev "+sample.files.nonmerged2 \
      +" --fasta_width 0 " \
      +" --fastq_allowmergestagger " \
      +" --fastq_minmergelen "+minlen \
      +" --fastq_maxns "+str(maxn) \
      +" --fastq_minovlen 20 " \
      +" --fastq_maxdiffs 4 " \
      +" --label_suffix _m1" \
      +" --threads 0"

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
        sample.stats.reads_merged = len(tmpf.readlines()) // 4

    LOGGER.info("Merged pairs - %d", sample.stats.reads_merged)

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
                                "nnnn"+read2s[1].strip(),
                                #comp(read2s[1].strip())[::-1],
                            read1s[2].strip(),
                            read1s[3].strip()+\
                                "nnnn"+read2s[3].strip()]
                                #read2s[3].strip()[::-1]]
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

    return sample


## This is hold-over code from pyrad V3 alignable, it's only used
## by loci2vcf so you could move it there if you like
def most_common(L):
    return max(itertools.groupby(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]






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
    LOGGER.info("zcat is using optim = %s", optim)

    ## is it gzipped
    cat = "cat"
    if raws[0].endswith(".gz"):
        cat = "gunzip -c"

    ### run splitter
    ### The -a flag tells split how long the suffix for each split file
    ### should be. It uses lowercase letters of the alphabet, so `-a 4`
    ### will have 26^4 possible tmp file names.
    cmd = " ".join([cat, raws[0], "|", "split", "-a", "4", "-l", str(optim),
                   "-", os.path.join(tmpdir, "chunk1_"+str(num)+"_")])
    _ = subprocess.check_call(cmd, shell=True)

    chunks1 = glob.glob(os.path.join(tmpdir, "chunk1_"+str(num)+"_*"))
    chunks1.sort()

    if "pair" in data.paramsdict["datatype"]:
        cmd = " ".join([cat, raws[1], "|", "split", "-a", "4", "-l", str(optim),
                  "-", os.path.join(tmpdir, "chunk2_"+str(num)+"_")])
        _ = subprocess.check_call(cmd, shell=True)

        chunks2 = glob.glob(os.path.join(tmpdir, "chunk2_"+str(num)+"_*"))
        chunks2.sort()
    
    else:
        chunks2 = [0]*len(chunks1)

    assert len(chunks1) == len(chunks2), \
        "R1 and R2 files are not the same length."

    return [raws[0], zip(chunks1, chunks2)]



def preview_truncate_fq(data, sample_fastq, nlines=None):
    """ 
    If we are running in preview mode, truncate the input fq.gz file so it'll 
    run quicker, just so we can see if it works. Input is tuple of the file
    names of the sample fq, and the # of lines to truncate to. Function 
    returns a list of one tuple of 1 or 2 elements depending on whether 
    you're doing paired or single end. The elements are paths to a temp files 
    of the sample fq truncated to some much smaller size.
    """

    ## Return a list of filenames
    truncated_fq = []

    ## If nlines is not set do something sensible, just truncate to 10%
    ## of the input read file size
    if not nlines:
       nlines = os.path.getsize(sample_fastq[0][0])/10 

    ## grab rawdata tuple pair from fastqs list [(x_R1_*, x_R2_*),]
    ## do not need to worry about multiple appended fastq files b/c preview
    ## mode will only want to sample from one file pair.
    for read in sample_fastq[0]:

        ## If the R2 is empty then exit the loop
        if not read:
            continue
        try:
            if read.endswith(".gz"):
                infile = gzip.open(os.path.realpath(read), 'rb')
            else:
                infile = open(os.path.realpath(read), 'rb')

            ## slice from data some multiple of 4 lines, no need to worry
            ## about truncate length being longer than the file this way.
            optim = nlines*4
            quarts = itertools.islice(infile, optim)

            ## write to a tmp file in the same place zcat_make_tmps would write
            with tempfile.NamedTemporaryFile('w+b', delete=False,
                          dir=os.path.realpath(data.dirs.project),
                          prefix=read+".preview_tmp_", 
                          suffix=".fq") as tmp_fq:
                tmp_fq.write("".join(quarts))
            ## save file name and close input
            truncated_fq.append(tmp_fq.name)
            infile.close()

        except AttributeError as _:
            ## R2 during SE is passed out as 0 
            #truncated_fq.append( 0 )
            raise

        except KeyboardInterrupt as _:
            LOGGER.info("""
    Caught keyboard interrupt during preview mode. Cleaning up preview files.
            """)
            ## clean up preview files
            try:
                truncated_fq.append(tmp_fq.name)
                for truncfile in truncated_fq:
                    if os.path.exists(truncfile):
                        os.remove(truncfile)
            except Exception as inst:
                LOGGER.debug("Error cleaning up truncated fq files: {}"\
                             .format(inst))
                raise
        except Exception as inst:
            LOGGER.debug("Some other stupid error - {}".format(inst))

    return [tuple(truncated_fq)]



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


#### Worker class to hold threaded or non-threaded views
class Worker():
    """ independent threaded CPUs"""
    def __init__(self, name, ipyclient, targets, func, data):
        ## Engine views
        self.client = ipyclient
        self.lview = self.client.load_balanced_view(targets=targets)
        self.dview = self.client.direct_view(targets=targets)
        self.ids = targets
        self.name = name
        
        ## how long to wait between checking for available proc
        self.delay = 0.5
        
        ## store the data and func
        self.data = data
        self.func = func
        
        ## is job finished
        self.result = {}
        
        ## finished results
        self.running = 0
        self.finished = 0
        
    
    def synchronize(self):
        """ checks whether it's finished """
        try:
            self.result[self.name] = self.running.get(0)
        except ipp.error.TimeoutError:
            pass
        except Exception as e:
            self.finished = 1
            #print(data, ': fail', e.traceback)
            self.result[self.name] = e
        else:
            self.finished = 1
            self.running = None
            

    def job(self, data, func):
        """ while job is running store async result in running."""
        self.running = self.lview.apply_async(func, data)
        
                
    def run(self):
        ## starts self.running
        self.job(self.data, self.func)
        while self.running:
            print(self.name, 'pending...')
            ## check for jobs finished, stops self.running
            self.synchronize()
            time.sleep(self.delay)
