#!/usr/bin/env python

"""
Reference based clustering across samples based on mapped positions
to a reference genome.
"""

import os
import sys
import glob
import subprocess as sps
import numpy as np
from loguru import logger
from pysam import AlignmentFile, FastaFile
from ipyrad.assemble.utils import IPyradError, chroms2ints
from ipyrad.core.progress_bar import AssemblyProgressBar


BIN_SAMTOOLS = os.path.join(sys.prefix, "bin", "samtools")
BIN_BEDTOOLS = os.path.join(sys.prefix, "bin", "bedtools")


class ClustMapAcrossReference:
    def __init__(self, step):
        self.step = step
        self.samples = step.samples
        self.quiet = step.quiet
        self.data = step.data

        # store path to clust database 
        for sname in step.samples:
            self.data.samples[sname].files.database = os.path.join(
                self.step.stepdir, 
                f"{self.data.name}_clust_database.fa")


    def run(self):
        """
        Runs the core step functions
        """
        # concat
        self.remote_concat_bams()

        # get extents of regions using bedtools merge
        self.remote_build_ref_regions()

        # build clusters from regions
        self.remote_build_ref_clusters()

        # concat aligned files (This is not necessary, chunk again in s7)
        self.remote_concat_alignments()


    def remote_concat_bams(self):
        """
        Merge bam files into a single large sorted indexed bam
        """
        args = (self.step.data, self.step.samples)
        rasyncs = {0: self.step.lbview.apply(concat_bams, *args)}

        # progress bar
        message = "concatenating bams"
        prog = AssemblyProgressBar(rasyncs, message, step=6, quiet=self.quiet)
        prog.block()
        prog.check()


    def remote_build_ref_regions(self):
        """
        call bedtools remotely and track progress
        """
        msg = "fetching regions"
        jobs = {0: self.step.lbview.apply(build_ref_regions, self.data)}
        prog = AssemblyProgressBar(jobs, msg, 6, self.step.quiet)
        prog.block()
        prog.check()
        self.step.regions = jobs[0].get()
        logger.debug('regions: {}...'.format(self.step.regions[:10]))


    def remote_build_ref_clusters(self):
        """
        build clusters and find variants/indels to store
        """       
        # send N jobs each taking chunk of regions
        ncpus = self.data.ncpus
        nloci = len(self.step.regions)
        optim = int((nloci // ncpus) + (nloci % ncpus))
        optim = int(np.ceil(optim / 2))
        logger.debug(f"using optim chunk size: {optim}")

        jobs = {}
        for idx, chunk in enumerate(range(0, nloci, optim)):
            region = self.step.regions[chunk: chunk + optim]
            if region:
                args = (self.data, idx, region)
                jobs[idx] = self.step.lbview.apply(build_ref_clusters, *args)

        # send jobs to func
        msg = "building database"
        prog = AssemblyProgressBar(jobs, msg, 6, self.step.quiet)
        prog.block()
        prog.check()


    def remote_concat_alignments(self):
        """
        concatenate fa chunks.
        """
        msg = "concat alignments"
        args = (self.data, self.samples)
        jobs = {0: self.step.lbview.apply(concat_alignments, *args)}
        prog = AssemblyProgressBar(jobs, msg, 6, self.step.quiet)
        prog.block()
        prog.check()


def concat_bams(data, samples):
    """
    Merge bam files into a single large sorted indexed bam
    """
    # concatenate consens bamfiles for all samples in this assembly
    catbam = os.path.join(data.tmpdir, f"{data.name}.cat.bam")
    cmd1 = [
        BIN_SAMTOOLS, 
        "merge", 
        "-f", 
        catbam,
    ]

    # Use the sample.files.consens to track branching
    for sname in samples:
        cmd1.append(data.samples[sname].files.consens)
    proc = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    err = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError(f"error in: {' '.join(cmd1)}: {err}")

    # sort the bam file
    cmd2 = [
        BIN_SAMTOOLS,
        "sort",
        "-T",
        catbam + '.tmp',
        "-o", 
        os.path.join(data.tmpdir, f"{data.name}.cat.sorted.bam"),
        catbam,
    ]
    proc = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)
    err = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError(f"error in: {' '.join(cmd2)}: {err}")        
    os.remove(catbam)

    # index the bam file
    cmd3 = [
        BIN_SAMTOOLS,
        "index", "-c",
        os.path.join(data.tmpdir, f"{data.name}.cat.sorted.bam",
        ),           
    ]
    proc = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE)
    err = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError(f"error in: {' '.join(cmd3)}: {err}")


def resolve_duplicates(keys, arr):
    """
    Tries to join together duplicate consens reads that were not previously
    collapsed, likely because there was no overlap of the sequences for one 
    or more samples, but there was for others. Joins two consens reads if the
    """
    newkeys = []
    snames = np.array([i.rsplit(":", 2)[0].rsplit("_", 1)[0] for i in keys])
    newarr = np.zeros((len(set(snames)) + 1, arr.shape[1]), dtype="S1")
    
    # put reference into arr
    newarr[0] = arr[0]
    
    # fill rest while merging dups
    nidx = 1
    seen = set()
    for sidx, _ in enumerate(keys):
        sname = snames[sidx]
        if sname not in seen:
            # add to list of seen names
            seen.add(sname)

            # get all rows of data for this sname (+1 b/c ref)
            didx = np.where(snames == sname)[0] + 1
            if didx.size > 1:
                iarr = arr[didx, :].view(np.uint8)
                iarr[iarr == 78] = 0
                iarr[iarr == 45] = 0
                if np.all(np.any(iarr == 0, axis=0)):
                    newarr[nidx] = iarr.max(axis=0).view("S1")
                else:
                    raise IPyradError("duplicate could not be resolved")
                # store key with reference to all dups
                ikeys = [keys[i - 1] for i in didx]
                fidxs = ";".join([i.rsplit("_", 1)[-1] for i in ikeys])
                newkeys.append("{}_{}".format(sname, fidxs))
                
            else:
                # store array data and orig key
                newarr[nidx] = arr[didx]
                newkeys.append(keys[sidx])
                
            nidx += 1
    
    # fill terminal edges with N again since array can increase
    newarr[newarr == b""] = b"N"
    return newkeys, newarr


def build_ref_regions(data):
    """
    Use bedtools to pull in consens reads overlapping some region of ref
    """
    cmd1 = [
        BIN_BEDTOOLS,
        "bamtobed",
        "-i", 
        os.path.join(data.tmpdir, f"{data.name}.cat.sorted.bam")
    ]

    cmd2 = [
        BIN_BEDTOOLS,    
        "merge", 
        "-d", "0",
        "-i", "-",
    ]

    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc2 = sps.Popen(
        cmd2, 
        stdin=proc1.stdout,
        stderr=sps.STDOUT,
        stdout=sps.PIPE,
    )
    result = proc2.communicate()[0].decode()
    if proc2.returncode:
        raise IPyradError(f"error in {' '.join(cmd2)}: {result}")
    regs = [i.split("\t") for i in result.strip().split("\n")]
    return [(i, int(j), int(k)) for i, j, k in regs]


def build_ref_clusters(data, idx, iregion):
    """
    Given a chunk of regions this will pull in the reference for each 
    region and then pull in all consens reads matching to that region. 
    It uses cigar info to align the consens reads with the ref. This 
    also merges consens from the same sample that were not merged 
    earlier, which is why we expect no duplicate samples in the output 
    of reference assemblies.
    """
    # prepare i/o for bamfile with mapped reads
    bamfile = os.path.join(data.tmpdir, f"{data.name}.cat.sorted.bam")
    alignments = AlignmentFile(bamfile, 'rb')

    # dict to map chromosome names to integers
    faidict = chroms2ints(data, False)

    # prepare i/o for pysam reference indexed
    reffai = FastaFile(data.params.reference_sequence)

    # store path to cluster bit
    outbit = os.path.join(data.tmpdir, f"aligned_{idx}.fa")

    # get clusters
    iregions = iter(iregion)
    clusts = []

    while 1:
        # pull in all consens reads mapping to a bed region
        try:
            region = next(iregions)
            reads = alignments.fetch(*region)
        except StopIteration:
            break

        # build a dict to reference seqs and cigars by name
        mstart = 9e12
        mend = 0
        rdict = {}
        for read in reads:
            rstart = read.reference_start
            rend = rstart + read.qlen
            mstart = min(mstart, rstart)
            mend = max(mend, rend)
            rdict[read.qname] = (read.seq, read.cigar, rstart, rend)
        keys = sorted(rdict.keys(), key=lambda x: x.rsplit(":", 2)[0])

        # pull in the reference for this region (1-indexed)
        refs = reffai.fetch(region[0], mstart, mend)

        # make empty array
        rlen = mend - mstart
        arr = np.zeros((len(keys) + 1, rlen), dtype=bytes)
        arr[0] = list(refs.upper())

        # fill arr with remaining samples
        for kidx, key in enumerate(keys):
            seq, cigar, start, end = rdict[key]

            # how far ahead of ref start and short of ref end is this read
            fidx = start - mstart
            eidx = arr.shape[1] - (mend - end)

            # enter into the array, trim end if longer than pulled ref
            arr[kidx + 1, fidx:eidx] = list(seq)[:eidx - fidx]

            # mod sequence according to cigar for indels and ambigs
            # csums is the location of impute on the seq, so it must be 
            # incremented by fidx and not extend past eidx
            for cidx, cig in enumerate(cigar):
                if cig[0] == 4:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[kidx + 1, csums] = arr[kidx + 1, csums].lower()
                if cig[0] == 1:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[kidx + 1, csums] = b"-"

        # fill terminal edges with N
        arr[arr == b""] = b"N"

        # duplicates merge here (only perfect merge on all Ns) and reshape
        # the array to match. This will need to be resolved in catgs...
        # if it does not merge then
        try:
            keys, arr = resolve_duplicates(keys, arr)
        except IPyradError:
            pass

        # get consens seq and variant site index 
        clust = [">reference_{}:{}:{}-{}\n{}".format(
            0, 
            faidict[region[0]] + 1, mstart + 1, mend + 1,   # 1-indexed
            b"".join(arr[0]).decode()
        )]
        for kidx, key in enumerate(keys):    
            clust.append(
                ">{}\n{}".format(key, b"".join(arr[kidx + 1]).decode())
            )
        clusts.append("\n".join(clust))

    # dump to temp file until concat in next step.
    with open(outbit, 'wt') as outfile:
        if clusts:
            outfile.write("\n//\n//\n".join(clusts) + "\n//\n//\n")
    alignments.close()


def concat_alignments(data, samples):
    """
    This step is not necessary... we just chunk it up again in step 7...
    it's nice having a file as a product, but why bother...
    It creates a header with names of all samples that were present when
    step 6 was completed. 
    """
    # get files
    globlist = glob.glob(os.path.join(data.tmpdir, "aligned_*.fa"))
    clustbits = sorted(
        globlist,
        key=lambda x: int(x.rsplit("_", 1)[1].split(".")[0]),
    )

    # write clusters to file with a header that has all samples in db        
    snames = sorted(samples)
    clustdb = data.samples[snames[0]].files.database
    with open(clustdb, 'wt') as out:
        out.write("#{}\n".format(",@".join(snames)))
        for clustfile in clustbits:
            with open(clustfile, 'r') as indata:
                dat = indata.read()
                if dat:
                    out.write(dat)  # + "//\n//\n")


if __name__ == "__main__":

    pass
