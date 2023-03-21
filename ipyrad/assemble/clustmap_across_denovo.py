#!/usr/bin/env python

"""Denovo clustering across samples using vsearch and muscle.

"""

from typing import TypeVar, Iterator
import sys
import gzip
import time
import random
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT

from loguru import logger
import numpy as np
from ipyrad.assemble.clustmap_within_denovo_utils import iter_muscle_alignments
from ipyrad.assemble.utils import IPyradError, fullcomp
from ipyrad.core.progress_bar import AssemblyProgressBar

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")
logger = logger.bind(name="ipyrad")
BIN = Path(sys.prefix) / "bin"
BIN_MUSCLE = str(BIN / "muscle")
BIN_VSEARCH = str(BIN / "vsearch")


class ClustMapAcrossDenovo:
    def __init__(self, step):
        self.lbview = step.lbview
        self.data = step.data
        self.quiet = step.quiet
        self.clust_database = self.data.stepdir / "alignment_database.fa.gz"
        self.samples = step.samples

    def run(self):
        """Runs the set of methods for denovo or reference method."""
        # prepare clustering inputs for hierarchical clustering
        self.concat_consens_files()

        # big clustering
        self.cluster_all()

        # build clusters
        self.build_denovo_clusters()

        # align denovo clusters
        self.align_denovo_clusters()

        # concat aligned files
        self.concat_alignments()

        # in addition to advancing sample state, store its clust_database
        for sample in self.samples.values():
            sample.files.database = self.clust_database
            sample.state = 6
        self.data.save_json()

    def concat_consens_files(self):
        """Prepares ONE large randomized concatenated consens input file.
        """
        jobs = {0: self.lbview.apply(build_concat_files, *(self.data, self.samples))}
        msg = "concatenating inputs"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()

    def cluster_all(self):
        """Cluster. Uses vsearch output to track progress.
        """
        # nthreads=0 defaults to using all cores
        jobs = {0: self.lbview.apply(cluster, self.data)}
        msg = "clustering across"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)        
        prog.block()
        prog.check()

    def build_denovo_clusters(self):
        """build denovo clusters from vsearch clustered seeds
        """
        # filehandles; if not multiple tiers then 'x' is jobid 0
        jobs = {}

        # sort utemp files, count seeds.
        utemp = self.data.tmpdir / "across.utemp"
        jobs[0] = self.lbview.apply(sort_seeds, utemp)
        jobs[1] = self.lbview.apply(count_seeds, utemp)
        msg = "counting clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()
        nseeds = int(jobs[1].get())
        args = (self.data, nseeds)
        jobs[2] = self.lbview.apply(write_denovo_cluster_chunks, *args)
        msg = "building clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()

    def align_denovo_clusters(self):
        """Distributes parallel jobs to align_to_array() function. 
        """
        path = self.data.tmpdir / "across_[0-9]*.unaligned_chunk"
        clustbits = list(path.parent.glob(path.name))
        jobs = {}
        start = time.time()
        for idx, clustbit in enumerate(clustbits):
            jobs[idx] = self.lbview.apply(write_alignments_across, clustbit)            
        msg = "aligning clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()
        logger.debug(f"aligning took: {time.time() - start}")

    def concat_alignments(self):
        """Concatenate aligned chunks.

        This step is not necessary... we just chunk it up again in step 7...
        it's nice having a file as a product, but why bother...
        It creates a header with names of all samples that were present when
        step 6 was completed. 
        """
        # get files
        globpath = self.data.tmpdir.glob("across_*.alignment")
        alignbits = sorted(
            globpath,
            key=lambda x: int(x.name.split("_", 1)[1].split(".")[0]))

        # store path to clust database
        with gzip.open(self.clust_database, 'w') as out:
            snames = sorted(self.samples)
            namestr = "#" + ",".join(["@" + i for i in snames]) + "\n"
            out.write(namestr.encode())
            for alignbit in alignbits:
                with open(alignbit, 'r', encoding="utf-8") as indat:
                    dat = indat.read().strip()
                    if dat:
                        out.write(dat.encode() + b"\n")


def build_concat_files(data, samples):
    """Write a concatenated file of consensus alleles (no RSWYMK/rswymk).

    Orders reads by length and shuffles randomly within length classes.
    """
    # tmp and output file paths
    catfile = data.tmpdir / "across-catcons.gz"
    allhaps = data.tmpdir / "across-cathaps.fa"
    allsort = data.tmpdir / "across-catsort.fa"
    allshuf = data.tmpdir / "across-catshuf.fa"

    # concatenate all of the gzipped consens files
    conshandles = [j.files.consens for i, j in samples.items()]
    conshandles.sort()
    cmd = ['cat'] + conshandles
    with open(catfile, 'w', encoding="utf-8") as out:
        with Popen(cmd, stdout=out, close_fds=True) as proc:
            proc.communicate()

    # a string of sed substitutions for temporarily replacing hetero sites
    # skips lines with '>', so it doesn't affect taxon names
    subs = ["/>/!s/W/A/g", "/>/!s/w/A/g", "/>/!s/R/A/g", "/>/!s/r/A/g",
            "/>/!s/M/A/g", "/>/!s/m/A/g", "/>/!s/K/T/g", "/>/!s/k/T/g",
            "/>/!s/S/C/g", "/>/!s/s/C/g", "/>/!s/Y/C/g", "/>/!s/y/C/g"]
    subs = ";".join(subs)

    # impute pseudo-haplo information to avoid mismatch at hetero sites
    # the read data with hetero sites is put back into clustered data later.
    # pipe passed data from gunzip to sed.
    cmd1 = ["gunzip", "-c", catfile]
    cmd2 = ["sed", subs]
    with Popen(cmd1, stdout=PIPE, close_fds=True) as proc1:
        with open(allhaps, 'w', encoding="utf-8") as out:
            with Popen(cmd2, stdin=proc1.stdout, stdout=out, close_fds=True) as proc2:
                proc2.communicate()
        proc1.stdout.close()

    # now sort the file using vsearch
    cmd1 = [
        BIN_VSEARCH,
        "--sortbylength", allhaps,
        "--fasta_width", "0",
        "--output", allsort,
    ]
    with Popen(cmd1, close_fds=True) as proc1:
        proc1.communicate()

    # shuffle sequences within size classes. Tested seed (8/31/2016)
    # shuffling works repeatably with seed.
    random.seed(data.hackers.random_seed)

    # open an iterator to lengthsorted file and grab two lines at at time
    with open(allshuf, 'w', encoding="utf-8") as out:
        with open(allsort, 'r', encoding="utf-8") as indat:

            # generator yields two lines at a time
            pairgen = zip(iter(indat), iter(indat))
            pair = next(pairgen)
            chunk = ["".join(pair)]
            oldlen = len(pair[1])

            # get all consens reads of the same length and write
            for pair in pairgen:
                if len(pair[1]) != oldlen:
                    random.shuffle(chunk)
                    out.write("".join(chunk))
                    chunk = ["".join(pair)]
                else:
                    chunk.append("".join(pair))

            # do the last chunk
            if chunk:
                random.shuffle(chunk)
                out.write("".join(chunk))

def cluster(data):
    """Runs vsearch cluster_smallmem on the concatenated data. 

    This can be very time consuming for super large and 
    non-overlapping data sets.
    """
    # get files for this jobid
    catshuf = data.tmpdir / "across-catshuf.fa"
    uhaplos = data.tmpdir / "across.utemp"
    hhaplos = data.tmpdir / "across.htemp"

    # parameters that vary by datatype
    # (too low of cov values yield too many poor alignments)
    strand = "plus"
    cov = 0.5         # 0.90
    if data.params.datatype in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.60
    elif data.params.datatype == "pairgbs":
        strand = "both"
        cov = 0.75     # 0.90

    # cluster using our custom sorted file (sorted by len)
    cmd = [
        BIN_VSEARCH,
        "--cluster_smallmem", str(catshuf),
        "--usersort",
        "--strand", strand,
        "--query_cov", str(cov),
        "--minsl", str(0.5),
        "--id", str(data.params.clust_threshold),
        "--userout", str(uhaplos),
        "--notmatched", str(hhaplos),
        "--userfields", "query+target+qstrand",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        "--fasta_width", "0",
        "--minseqlength", str(data.params.filter_min_trim_len),
        "--threads", str(max(1, data.ipcluster['threads'])),
        # "-fulldp", # no longer necessary w/ new vsearch
    ]

    # send command to logger and wait a second so it clears.
    print(" ".join(cmd))
    time.sleep(1.01)

    # get progress from the stdout of the subprocess
    with Popen(cmd, stdout=PIPE, stderr=STDOUT, close_fds=True) as proc:
        proc.communicate()
    print(100)


def count_seeds(uhandle: Path) -> int:
    """Use bash commands to quickly count N seeds from utemp file"""
    with open(uhandle, 'r', encoding="utf-8") as insort:
        cmd1 = ["cut", "-f", "2"]
        cmd2 = ["uniq"]
        cmd3 = ["wc"]
        with Popen(cmd1, stdin=insort, stdout=PIPE, close_fds=True
            ) as proc1:
            with Popen(cmd2, stdin=proc1.stdout, stdout=PIPE, close_fds=True
                ) as proc2:
                with Popen(cmd3, stdin=proc2.stdout, stdout=PIPE, close_fds=True
                    ) as proc3:
                    res = proc3.communicate()
                    nseeds = int(res[0].split()[0])
                    proc1.stdout.close()
                    proc2.stdout.close()
                    proc3.stdout.close()
    return nseeds

def sort_seeds(uhandle: Path) -> None:
    """sort seeds from cluster results"""
    cmd = ["sort", "-k", "2", str(uhandle), "-o", str(uhandle) + ".sort"]
    with Popen(cmd, close_fds=True) as proc:
        proc.communicate()

def iter_build_denovo_clusters(data: Assembly) -> Iterator:
    """Build unaligned clusters from vsearch output. 

    Splits paired clusters into before and after the nnnn separator.
    """
    # load all concat fasta files into a dictionary (memory concerns here...)
    conshandle = data.tmpdir / "across-catcons.gz"
    allcons = {}
    with gzip.open(conshandle, 'rt') as iocons:
        cons = zip(*[iocons] * 2)
        for namestr, seq in cons:
            nnn, sss = [i.strip() for i in (namestr, seq)]
            allcons[nnn[1:]] = sss

    # iterate through usort grabbing seeds and matches
    usort_file = data.tmpdir / "across.utemp.sort"
    with open(usort_file, 'r', encoding="utf-8") as insort:
        lastseed = None
        fseqs = []

        for line in insort:
            hit, seed, orient = line.strip().split()
        
            # this is a new seed.
            if seed != lastseed:
                # should we bother to return singleton?
                if fseqs:
                    yield "\n".join(fseqs)
                fseqs = [f">{seed};*\n{allcons[seed]}"]
                lastseed = seed

            # this is a hit to the seed.
            if orient == "-":
                seq = fullcomp(allcons[hit])[::-1]
            else:
                seq = allcons[hit]
            fseqs.append(f">{hit};{orient}\n{seq}")
    if fseqs:
        yield "\n".join(fseqs)
    del allcons

def write_denovo_cluster_chunks(data: Assembly, nseeds: int):
    """Write built clusters to chunk files."""
    # chunksize to break into parallel aligning jobs
    optim = int(np.ceil(nseeds / (data.ncpus * 4)))
    print(f"chunking align files; optim={optim}")

    # iterate over built clusters
    idx = 0
    clusters = []
    for idx, clust in enumerate(iter_build_denovo_clusters(data)):
        clusters.append(clust)

        if not idx % optim:
            path = data.tmpdir / f"across_{idx}.unaligned_chunk"
            with open(path, 'w', encoding="utf-8") as out:
                out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
                clusters = []

    # write whatever is left over to the clusts file
    if clusters:
        path = data.tmpdir / f"across_{idx}.unaligned_chunk"        
        with open(path, 'w', encoding="utf-8") as out:
            out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")

def iter_alignment_across_format(handle: Path) -> Iterator:
    """Generator to yield filtered, aligned, paired, and sorted clusters.

    This differs from the within-sample function by allowing any
    size of cluster (no upper limit) and not using the derep numbers
    to sort names after alignment.

    Example for testing
    -------------------
    >>> tmpfile = "/tmp/test_tmp_clustmap/1A_0_chunk_0.ali"
    >>> align_gen = iter_alignment_format(tmpfile)
    >>> print(next(align_gen))
    """
    for ali1, ali2 in iter_muscle_alignments(handle, int(1e9)):
        # get dict mapping headers to sequences
        head_to_seq1 = {}
        for line in ali1:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq1:
                    head_to_seq1[key] += line.strip()
                else:
                    head_to_seq1[key] = line.strip()
        head_to_seq2 = {}
        for line in ali2:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq2:
                    head_to_seq2[key] += line.strip()
                else:
                    head_to_seq2[key] = line.strip()

        # sort the first reads
        try:
            keys = sorted(head_to_seq1)
            seed = [i for i in keys if i[-1] == "*"][0]
        except IndexError as inst:
            print(keys)
            raise IndexError(f"{keys}\n\n{ali1}") from inst

        seed = keys.pop(keys.index(seed))
        order = [seed] + keys

        alignment = []
        if head_to_seq2:
            for key in order:
                alignment.append(
                    f"{key}\n{head_to_seq1[key]}nnnn{head_to_seq2[key]}")
        else:
            for key in order:
                alignment.append(f"{key}\n{head_to_seq1[key]}")
        yield "\n".join(alignment)        

def write_alignments_across(handle: Path) -> None:
    """Writes alignments to tmpfiles.
    
    Iterates over chunks of aligned loci from generator and
    writes to a concatenated output file.
    """
    print(f"aligning: {handle}") # engine sends to logger.info
    tmpout = handle.with_suffix(".alignment")
    chunk = []
    with open(tmpout, 'w', encoding="utf-8") as out:
        for idx, alignment in enumerate(iter_alignment_across_format(handle)):
            chunk.append(alignment)
            if not idx % 1_000:
                out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")
                chunk = []
        if chunk:
            out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")

def reconcat_across(data: Assembly, sample: Sample) -> None:
    """Concatenate .aligned chunks into a single .clusters.gz file.
    """
    chunks = list(data.tmpdir.glob(f"{sample.name}_chunk_[0-9].alignment"))
    chunks.sort(key=lambda x: int(x.name.rsplit("_", 1)[-1][:-10]))

    # concatenate finished clusters
    clustfile = data.stepdir / f"{sample.name}.clusters.gz"
    with gzip.open(clustfile, 'w') as out:
        for fname in chunks:
            with open(fname, 'r', encoding="utf-8") as infile:
                dat = infile.read().strip()
                if dat:
                    out.write(dat.encode() + b"\n")



# def threaded_muscle(chunk):
#     """
#     Runs the i/o limited muscle calls on non-blocking threads.
#     """
#     # get all clusters
#     with open(chunk, 'rt') as infile:
#         clusts = infile.read().split("//\n//\n")

#     # submit jobs on THIS ENGINE to run on 5 threads.
#     futures = []
#     with concurrent.futures.ThreadPoolExecutor(5) as pool:
#         for clust in clusts:
#             futures.append(pool.submit(muscle_align, clust))
#         for fut in concurrent.futures.as_completed(futures):
#             assert fut.done() and not fut.cancelled()
#     stacks = [i.result() for i in futures]

#     # write to file when chunk is finished
#     odx = chunk.rsplit("_")[-1]
#     odir = os.path.dirname(chunk)
#     alignfile = os.path.join(odir, "aligned_{}.fa".format(odx))
#     with open(alignfile, 'wt') as outfile:
#         outfile.write("\n//\n//\n".join(stacks))



# def store_alleles(seqs):
#     """
#     Returns a mask selecting columns with lower case calls, and 
#     a boolean of whether or not any exist. This is used to put them 
#     back into alignments after muscle destroys all this info during
#     alignment.
#     """
#     # get shape of the array and new empty array
#     shape = (len(seqs), max([len(i) for i in seqs]))
#     arrseqs = np.zeros(shape, dtype=np.bytes_)

#     # iterate over rows 
#     for row in range(arrseqs.shape[0]):
#         seqsrow = seqs[row]
#         arrseqs[row, :len(seqsrow)] = list(seqsrow)

#     # mask the lo...
#     amask = np.char.islower(arrseqs)
#     if np.any(amask):
#         return amask, True
#     else:
#         return amask, False


# def retrieve_alleles_after_aligning(intarr, amask):
#     """
#     Imputes lower case allele calls back into alignments 
#     while taking account for spacing caused by insertions.
#     """
#     newmask = np.zeros(intarr.shape, dtype=np.bool_)
    
#     for ridx in range(intarr.shape[0]):
#         iarr = intarr[ridx]
#         indidx = np.where(iarr == 45)[0]
        
#         # if no indels then simply use the existing mask
#         if not indidx.size:
#             newmask[ridx] = amask[ridx]
            
#         # if indels that impute 
#         else:
#             allrows = np.arange(amask.shape[1])
#             mask = np.ones(allrows.shape[0], dtype=np.bool_)
#             for idx in indidx:
#                 if idx < mask.shape[0]:
#                     mask[idx] = False
#             not_idx = allrows[mask == 1]
            
#             # fill in new data into all other spots
#             newmask[ridx, not_idx] = amask[ridx, :not_idx.shape[0]]
#     return newmask

