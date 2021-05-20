#!/usr/bin/env python

"""
Denovo clustering across samples using vsearch and muscle.
"""

import os
import sys
import gzip
import glob
import time
import random
import itertools
import subprocess as sps
import concurrent.futures
from loguru import logger
import numpy as np
from ipyrad.assemble.utils import IPyradError, fullcomp
from ipyrad.core.progress_bar import AssemblyProgressBar


BIN_MUSCLE = os.path.join(sys.prefix, "bin", "muscle")
BIN_VSEARCH = os.path.join(sys.prefix, "bin", "vsearch")


class ClustMapAcrossDenovo:
    def __init__(self, step):
        self.lbview = step.lbview
        self.data = step.data
        self.data.samples = step.samples
        self.tmpdir = step.tmpdir
        self.stepdir = step.stepdir
        self.quiet = step.quiet
        self.clust_database = os.path.join(
            self.stepdir, "alignment_database.fa")


    def run(self):
        """
        Runs the set of methods for denovo or reference method
        """
        # DENOVO
        if self.data.params.assembly_method == "denovo":

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

        for sname in self.data.samples:
            self.data.samples[sname].files.database = self.clust_database
            self.data.samples[sname].state = 6
        self.data.save_json()


    def concat_consens_files(self):
        """
        Prepares ONE large randomized concatenated consens input file.
        """
        jobs = {0: self.lbview.apply(build_concat_files, self.data)}
        msg = "concatenating inputs"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()


    def cluster_all(self):
        """
        Hierarchical cluster. Uses vsearch output to track progress.
        """
        # nthreads=0 defaults to using all cores
        jobs = {0: self.lbview.apply(cluster, self.data)}
        msg = "clustering across"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)        
        prog.block()
        prog.check()


    def build_denovo_clusters(self):
        """
        build denovo clusters from vsearch clustered seeds
        """
        # filehandles; if not multiple tiers then 'x' is jobid 0
        jobs = {}

        # sort utemp files, count seeds.
        utemp = os.path.join(self.tmpdir, "across.utemp")
        jobs[0] = self.lbview.apply(sort_seeds, utemp)
        jobs[1] = self.lbview.apply(count_seeds, utemp)
        msg = "counting clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()
        nseeds = int(jobs[1].get())
        jobs[2] = self.lbview.apply(
            build_denovo_clusters, *(self.data, nseeds)
        )
        msg = "building clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()


    def align_denovo_clusters(self):
        """
        Distributes parallel jobs to align_to_array() function. 
        """
        globpath = os.path.join(self.tmpdir, "unaligned.chunk_*")
        clustbits = glob.glob(globpath)
        jobs = {}
        start = time.time()
        for idx, clustbit in enumerate(clustbits):
            jobs[idx] = self.lbview.apply(threaded_muscle, clustbit)
            # args = (self.data, list(self.data.samples.values()), clustbits[idx])
            # jobs[idx] = self.lbview.apply(align_to_array, *args)
        msg = "aligning clusters"
        prog = AssemblyProgressBar(jobs, msg, 6, self.quiet)
        prog.block()
        prog.check()
        logger.debug(f"aligning took: {time.time() - start}")


    def concat_alignments(self):
        """
        This step is not necessary... we just chunk it up again in step 7...
        it's nice having a file as a product, but why bother...
        It creates a header with names of all samples that were present when
        step 6 was completed. 
        """
        # get files
        globpath = os.path.join(self.tmpdir, "aligned_*.fa")
        alignbits = glob.glob(globpath)
        alignbits = sorted(
            alignbits, 
            key=lambda x: int(os.path.basename(x).rsplit("_")[1].split(".")[0])
        )

        # store path to clust database 
        with open(self.clust_database, 'wt') as out:
            snames = sorted(self.data.samples)
            out.write("#" + ",".join(["@" + i for i in snames]) + "\n")
            for alignbit in alignbits:
                with open(alignbit, 'rt') as indat:
                    dat = indat.read()
                    if dat:
                        out.write(dat)


def build_concat_files(data):
    """
    Make a concatenated consens file with sampled alleles 
    (no RSWYMK/rswymk). Orders reads by length and shuffles randomly 
    within length classes
    """
    conshandles = [data.samples[i].files.consens for i in data.samples]
    conshandles.sort()

    # concatenate all of the gzipped consens files
    cmd = ['cat'] + conshandles
    catfile = os.path.join(data.tmpdir, "across-catcons.gz")
    with open(catfile, 'w') as output:
        call = sps.Popen(cmd, stdout=output, close_fds=True)
        call.communicate()

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

    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, close_fds=True)
    allhaps = catfile.replace("-catcons.gz", "-cathaps.fa")
    with open(allhaps, 'w') as output:
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=output, close_fds=True)
        proc2.communicate()
    proc1.stdout.close()

    # now sort the file using vsearch
    allsort = catfile.replace("-catcons.gz", "-catsort.fa")
    cmd1 = [
        BIN_VSEARCH,
        "--sortbylength", allhaps,
        "--fasta_width", "0",
        "--output", allsort,
    ]
    proc1 = sps.Popen(cmd1, close_fds=True)
    proc1.communicate()

    # shuffle sequences within size classes. Tested seed (8/31/2016)
    # shuffling works repeatably with seed.
    random.seed(data.hackers.random_seed)

    # open an iterator to lengthsorted file and grab two lines at at time
    allshuf = catfile.replace("-catcons.gz", "-catshuf.fa")
    outdat = open(allshuf, 'wt')
    indat = open(allsort, 'rt')
    idat = zip(iter(indat), iter(indat))
    done = 0

    chunk = [next(idat)]
    while not done:
        # grab 2-lines until they become shorter (unless there's only one)
        oldlen = len(chunk[-1][-1])
        while 1:
            try:
                dat = next(idat)
            except StopIteration:
                done = 1
                break
            if len(dat[-1]) == oldlen:
                chunk.append(dat)
            else:
                # send the last chunk off to be processed
                random.shuffle(chunk)
                outdat.write("".join(itertools.chain(*chunk)))
                # start new chunk
                chunk = [dat]
                break

    # do the last chunk
    random.shuffle(chunk)
    outdat.write("".join(itertools.chain(*chunk)))
    indat.close()
    outdat.close()


def cluster(data, nthreads=0):
    """
    Runs vsearch cluster_smallmem on the concatenated data. This can 
    be very time consuming for super large and non-overlapping data
    sets.
    """
    # get files for this jobid
    catshuf = os.path.join(data.tmpdir, "across-catshuf.fa")
    uhaplos = os.path.join(data.tmpdir, "across.utemp")
    hhaplos = os.path.join(data.tmpdir, "across.htemp")

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

    cmd = [
        BIN_VSEARCH,
        "-cluster_smallmem", catshuf,
        "-strand", strand,
        "-query_cov", str(cov),
        "-minsl", str(0.5),
        "-id", str(data.params.clust_threshold),
        "-userout", uhaplos,
        "-notmatched", hhaplos,
        "-userfields", "query+target+qstrand",
        "-maxaccepts", "1",
        "-maxrejects", "0",
        "-fasta_width", "0",
        "--minseqlength", str(data.params.filter_min_trim_len),
        "-threads", str(nthreads),
        "-fulldp",
        "-usersort",
    ]

    # send command to logger and wait a second so it clears.
    print(" ".join(cmd))
    time.sleep(1.01)

    # get progress from the stdout of the subprocess
    proc = sps.Popen(
        cmd, stdout=sps.PIPE, stderr=sps.STDOUT, close_fds=True,
    )
    proc.communicate()
    print(100)


def count_seeds(uhandle):
    """
    uses bash commands to quickly count N seeds from utemp file
    """
    with open(uhandle, 'r') as insort:
        cmd1 = ["cut", "-f", "2"]
        cmd2 = ["uniq"]
        cmd3 = ["wc"]
        proc1 = sps.Popen(cmd1, stdin=insort, stdout=sps.PIPE, close_fds=True)
        proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, close_fds=True)
        proc3 = sps.Popen(cmd3, stdin=proc2.stdout, stdout=sps.PIPE, close_fds=True)
        res = proc3.communicate()
        nseeds = int(res[0].split()[0])
        proc1.stdout.close()
        proc2.stdout.close()
        proc3.stdout.close()
    return nseeds


def sort_seeds(uhandle):
    """
    sort seeds from cluster results
    """
    cmd = ["sort", "-k", "2", uhandle, "-o", uhandle + ".sort"]
    proc = sps.Popen(cmd, close_fds=True)
    proc.communicate()


def build_denovo_clusters(data, nseeds):
    """
    Builds unaligned clusters from vsearch output. Splits paired 
    read clusters into before and after the nnnn separator.
    """
    # load all concat fasta files into a dictionary (memory concerns here...)
    conshandle = os.path.join(data.tmpdir, "across-catcons.gz")
    allcons = {}
    with gzip.open(conshandle, 'rt') as iocons:
        cons = zip(*[iocons] * 2)
        for namestr, seq in cons:
            nnn, sss = [i.strip() for i in (namestr, seq)]
            allcons[nnn[1:]] = sss

    # chunksize to break into parallel aligning jobs
    optim = np.ceil(nseeds / (data.ncpus * 4))

    # iterate through usort grabbing seeds and matches
    usort_file = os.path.join(data.tmpdir, "across.utemp.sort")
    with open(usort_file, 'rt') as insort:
        # iterator, seed null, and seqlist null
        isort = iter(insort)
        loci = 0
        lastseed = 0
        fseqs = []
        seqlist = []
        seqsize = 0

        while 1:
            try:
                hit, seed, ori = next(isort).strip().split()
            except StopIteration:
                break
        
            # store hit if still matching to same seed
            if seed == lastseed:
                if ori == "-":
                    seq = fullcomp(allcons[hit])[::-1]
                else:
                    seq = allcons[hit]
                fseqs.append(">{}\n{}".format(hit, seq))

            # store seed and hit (to a new cluster) if new seed.
            else:  
                # store the last fseq, count it, and clear it
                if fseqs:
                    seqlist.append("\n".join(fseqs))
                    seqsize += 1
                    fseqs = []

                # occasionally write to file
                if seqsize >= optim:
                    if seqlist:
                        loci += seqsize
                        pathname = os.path.join(
                            data.tmpdir, f"unaligned.chunk_{loci}")
                        with open(pathname, 'wt') as clustout:
                            clustout.write(
                                "\n//\n//\n".join(seqlist) + "\n//\n//\n")
                        # reset counter and list
                        seqlist = []
                        seqsize = 0

                # store the new seed on top of fseqs
                fseqs.append(">{}\n{}".format(seed, allcons[seed]))
                lastseed = seed

                # store the first hit to the seed
                seq = allcons[hit]
                if ori == "-":
                    seq = fullcomp(seq)[::-1]
                fseqs.append(">{}\n{}".format(hit, seq))

    # write whatever is left over to the clusts file
    if fseqs:
        seqlist.append("\n".join(fseqs))
        seqsize += 1
        loci += seqsize
    if seqlist:
        pathname = os.path.join(data.tmpdir, f"unaligned.chunk_{loci}")
        with open(pathname, 'wt') as clustsout:
            clustsout.write("\n//\n//\n".join(seqlist) + "\n//\n//\n")
    # final progress and cleanup
    del allcons


def muscle_align(clust):
    """
    Run muscle alignment on a subprocess
    """
    # easy peezy for unpaired data
    if "nnnn" not in clust:
        proc = sps.Popen(
            [BIN_MUSCLE, '-in', '-'], 
            stdin=sps.PIPE, 
            stdout=sps.PIPE,
            universal_newlines=True,
        )
        out, _ = proc.communicate(clust)

        # remove the newline breaks in the sequences
        out = "\n>".join([
            i[::-1].replace("\n", "", 1)[::-1] 
            for i in out.strip().split("\n>")
        ])
        return out

    # for paired data, split and align each separate, then rejoin
    # split the read1s and read2s
    inseqs = {}
    for line in clust.split():
        if line.startswith(">"):
            name = line.strip()
            inseqs[name] = ""
        else:
            if "nnnn" in line:
                inseqs[name] = line.strip().split("nnnn", 1)
            else:
                inseqs[name] = (line.strip(), "NNNNNNNNNN")
    
    # make back into strings
    clust1 = "\n".join([f"{key}\n{inseqs[key][0]}" for key in inseqs])
    clust2 = "\n".join([f"{key}\n{inseqs[key][1]}" for key in inseqs])
    
    # run each individually on subprocess
    proc = sps.Popen(
        [BIN_MUSCLE, '-in', '-'], 
        stdin=sps.PIPE, 
        stdout=sps.PIPE,
        universal_newlines=True,
    )
    out1, _ = proc.communicate(clust1)
    proc = sps.Popen(
        [BIN_MUSCLE, '-in', '-'], 
        stdin=sps.PIPE, 
        stdout=sps.PIPE,
        universal_newlines=True,
    )
    out2, _ = proc.communicate(clust2)
    
    # split the clust1s and clust2s to align reads w/ names b/c 
    # muscle doesn't retain the input order.
    inseqs = {}
    for line in out1.split():
        if line.startswith(">"):
            name = line.strip()
            inseqs[name] = ""
        else:
            inseqs[name] += line.strip()

    for line in out2.split():
        if line.startswith(">"):
            name = line.strip()
        else:
            if "nnnn" not in inseqs[name]:
                inseqs[name] += "nnnn" + line.strip()
            else:
                inseqs[name] += line.strip()
    
    return "\n".join([f"{key}\n{inseqs[key]}" for key in sorted(inseqs)])


def threaded_muscle(chunk):
    """
    Runs the i/o limited muscle calls on non-blocking threads.
    """
    # get all clusters
    with open(chunk, 'rt') as infile:
        clusts = infile.read().split("//\n//\n")

    # submit jobs on THIS ENGINE to run on 5 threads.
    futures = []
    with concurrent.futures.ThreadPoolExecutor(5) as pool:
        for clust in clusts:
            futures.append(pool.submit(muscle_align, clust))
        for fut in concurrent.futures.as_completed(futures):
            assert fut.done() and not fut.cancelled()
    stacks = [i.result() for i in futures]

    # write to file when chunk is finished
    odx = chunk.rsplit("_")[-1]
    odir = os.path.dirname(chunk)
    alignfile = os.path.join(odir, "aligned_{}.fa".format(odx))
    with open(alignfile, 'wt') as outfile:
        outfile.write("\n//\n//\n".join(stacks) + "\n//\n//\n")


def align_to_array(data, samples, chunk):
    """
    Opens a tmp clust chunk and iterates over align jobs.
    """
    # data are already chunked, read in the whole thing
    with open(chunk, 'rt') as infile:
        clusts = infile.read().split("//\n//\n")[:-1]    

    # snames to ensure sorted order
    samples.sort(key=lambda x: x.name)

    # create a persistent shell for running muscle in. 
    proc = sps.Popen(["bash"], stdin=sps.PIPE, stdout=sps.PIPE, bufsize=0)

    # iterate over clusters until finished
    allstack = []
    for ldx in range(len(clusts)):
        istack = []
        lines = clusts[ldx].strip().split("\n")
        names = lines[::2]
        seqs = lines[1::2]

        # skip aligning and continue if duplicates present (locus too big)
        # but reshape locs to be same lengths by adding --- to end, this 
        # simplifies handling them in step7 (they're still always filtered)
        unames = set([i.rsplit("_", 1)[0] for i in names])
        if len(unames) < len(names):
            longname = max([len(i) for i in seqs])
            seqs = [i.ljust(longname, "-") for i in seqs]
            istack = [">{}\n{}".format(i[1:], j) for i, j in zip(names, seqs)]
            allstack.append("\n".join(istack))
            continue

        # else locus looks good, align it.
        # is there a paired-insert in any samples in the locus?
        try:

            # try to split cluster list at nnnn separator for each read
            left = [i.split("nnnn")[0] for i in seqs]
            right = [i.split("nnnn")[1] for i in seqs]

            # align separately
            istack1 = muscle_it(proc, names, left)
            istack2 = muscle_it(proc, names, right)

            # combine in order
            for sdx in range(len(istack1)):
                n1, s1 = istack1[sdx].split("\n")
                s2 = istack2[sdx].split("\n")[-1]
                istack.append(n1 + "\n" + s1 + "nnnn" + s2)

        # no insert just align a single locus
        except IndexError:
            istack = muscle_it(proc, names, seqs)

        # store the locus
        if istack:
            allstack.append("\n".join(istack))

    # cleanup
    proc.stdout.close()
    if proc.stderr:
        proc.stderr.close()
    proc.stdin.close()
    proc.wait()

    # write to file when chunk is finished
    odx = chunk.rsplit("_")[-1]
    alignfile = os.path.join(data.tmpdir, "aligned_{}.fa".format(odx))
    with open(alignfile, 'wt') as outfile:
        outfile.write("\n//\n//\n".join(allstack) + "\n//\n//\n")



def muscle_it(proc, names, seqs):
    """
    Align with muscle, ensure name order, and return as string
    """  
    istack = []

    # append counter to names because muscle doesn't retain order
    nnames = [">{};*{}".format(j[1:], i) for i, j in enumerate(names)]
    
    # make back into strings
    cl1 = "\n".join(["\n".join(i) for i in zip(nnames, seqs)])

    # store allele (lowercase) info, returns mask with lowercases
    # amask, abool = store_alleles(seqs)

    # send align1 to the bash shell (TODO: check for pipe-overflow)
    cmd1 = ("echo -e '{}' | {} -quiet -in - ; echo {}"
            .format(cl1, BIN_MUSCLE, "//\n"))
    proc.stdin.write(cmd1.encode())

    # read the stdout by line until splitter is reached
    align1 = []
    for line in iter(proc.stdout.readline, b'//\n'):
        align1.append(line.decode())

    # reorder b/c muscle doesn't keep order
    lines = "".join(align1)[1:].split("\n>")
    dalign1 = dict([i.split("\n", 1) for i in lines])
    keys = sorted(
        dalign1.keys(), 
        key=lambda x: int(x.rsplit("*")[-1])
    )
    seqarr = np.zeros(
        (len(names), len(dalign1[keys[0]].replace("\n", ""))),
        dtype='S1',
    )
    for kidx, key in enumerate(keys):
        concatseq = dalign1[key].replace("\n", "")
        seqarr[kidx] = list(concatseq)

    # get alleles back using fast jit'd function.
    if np.sum(amask):
        intarr = seqarr.view(np.uint8)
        iamask = retrieve_alleles_after_aligning(intarr, amask)
        seqarr[iamask] = np.char.lower(seqarr[iamask])

    # sort in sname (alphanumeric) order. 
    istack = []    
    wkeys = np.argsort([i.rsplit("_", 1)[0] for i in keys])
    for widx in wkeys:
        wname = names[widx]
        istack.append(
            "{}\n{}".format(wname, b"".join(seqarr[widx]).decode()))
    return istack



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

