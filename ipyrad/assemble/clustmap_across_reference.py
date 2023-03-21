#!/usr/bin/env python

"""Step 6 reference-based.

Reference based clustering across samples based on mapped positions
to a reference genome.
"""

import sys
import gzip
import subprocess as sps
from pathlib import Path
import numpy as np
from loguru import logger
from pysam import AlignmentFile, FastaFile  # pylint: disable=no-name-in-module
from ipyrad.assemble.utils import IPyradError, chroms2ints
from ipyrad.core.progress_bar import AssemblyProgressBar

BIN = Path(sys.prefix) / "bin"
BIN_SAMTOOLS = BIN / "samtools"
BIN_BEDTOOLS = BIN / "bedtools"

logger = logger.bind(name="ipyrad")


class ClustMapAcrossReference:
    def __init__(self, step):
        self.step = step
        self.samples = step.samples
        self.quiet = step.quiet
        self.data = step.data
        self.clust_database = self.data.stepdir / "alignment_database.fa.gz"

    def run(self):
        """Runs the core step functions"""
        # concat
        self.remote_concat_bams()

        # get extents of regions using bedtools merge
        self.remote_build_ref_regions()

        # build clusters from regions
        self.remote_build_ref_clusters()

        # concat aligned files (This is not necessary, chunk again in s7)
        self.concat_alignments()

        # store path to clust database
        for sample in self.samples.values():
            sample.files.database = self.clust_database
        self.data.save_json()

    def remote_concat_bams(self):
        """Merge bam files into a single large sorted indexed bam"""
        args = (self.step.data, self.step.samples)
        rasyncs = {0: self.step.lbview.apply(concat_bams, *args)}
        message = "concatenating bams"
        prog = AssemblyProgressBar(rasyncs, message, step=6, quiet=self.quiet)
        prog.block()
        prog.check()

    def remote_build_ref_regions(self):
        """Call bedtools remotely and track progress."""
        msg = "fetching regions"
        jobs = {0: self.step.lbview.apply(build_ref_regions, self.data)}
        prog = AssemblyProgressBar(jobs, msg, 6, self.step.quiet)
        prog.block()
        prog.check()
        self.step.regions = jobs[0].get()
        # logger.debug('regions: {}...'.format(self.step.regions[:10]))

    def remote_build_ref_clusters(self):
        """Build clusters and find variants/indels to store."""
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

    def concat_alignments(self):
        """Concat denovo alignments.

        This step is not necessary... we just chunk it up again in step 7...
        it's nice having a file as a product, but why bother...
        It creates a header with names of all samples that were present when
        step 6 was completed.
        """
        # get files
        globlist = self.data.tmpdir.glob("aligned_*.fa")
        clustbits = sorted(
            globlist,
            key=lambda x: int(x.stem.rsplit("_", 1)[1]),
        )

        # write clusters to file with a header that has all samples in db
        # store path to clust database
        snames = sorted(self.samples)
        with gzip.open(self.clust_database, 'w') as out:
            namestr = "#" + ",".join(["@" + i for i in snames]) + "\n"
            out.write(namestr.encode())
            for alignbit in clustbits:
                with open(alignbit, 'r', encoding="utf-8") as indat:
                    dat = indat.read().strip()
                    if dat:
                        out.write(dat.encode() + b"\n")


def concat_bams(data, samples):
    """Merge bam files into a single large sorted indexed bam
    """
    catbam = data.tmpdir / f"{data.name}.cat.bam"
    cmd1 = [
        str(BIN_SAMTOOLS),
        "merge",
        "-f",
        str(catbam),
    ]

    # Use the sample.files.consens to track branching
    for sname in samples:
        cmd1.append(str(samples[sname].files.consens))
    proc = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    err = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError(
            f"An error occurred in samtools merge. Check for memory limits.\n"
            f": {' '.join(cmd1)}: {err}"
        )

    # sort the bam file
    cmd2 = [
        str(BIN_SAMTOOLS),
        "sort",
        "-T",
        str(catbam) + '.tmp',
        "-o",
        str(data.tmpdir / f"{data.name}.cat.sorted.bam"),
        str(catbam),
    ]
    proc = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)
    err = proc.communicate()[0].decode()
    if proc.returncode:
        raise IPyradError(
            f"An error occurred in samtools sort:\n{' '.join(cmd2)}: {err}")
    catbam.unlink()

    # index the bam file
    cmd3 = [
        str(BIN_SAMTOOLS),
        "index", "-c",
        str(data.tmpdir / f"{data.name}.cat.sorted.bam")
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
    """Use bedtools to pull in consens reads overlapping some region of ref
    """
    cmd1 = [
        str(BIN_BEDTOOLS),
        "bamtobed",
        "-i",
        str(data.tmpdir / f"{data.name}.cat.sorted.bam")
    ]

    cmd2 = [
        str(BIN_BEDTOOLS),
        "merge",
        "-d", "0",
        "-i", "-",
    ]
    # print(cmd1)
    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc2 = sps.Popen(
        cmd2,
        stdin=proc1.stdout,
        stderr=sps.STDOUT,
        stdout=sps.PIPE,
    )
    # print(cmd2)
    result = proc2.communicate()[0].decode()
    if proc2.returncode:
        raise IPyradError(f"error in {' '.join(cmd2)}: {result}")
    print(result)
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
    bamfile = data.tmpdir / f"{data.name}.cat.sorted.bam"
    alignments = AlignmentFile(str(bamfile), 'rb')

    # dict to map chromosome names to integers
    faidict = chroms2ints(data, False)
    print("CHROMS2INTS CHANGED, CHECK INTS INDEXING")

    # prepare i/o for pysam reference indexed
    reffai = FastaFile(str(data.params.reference_sequence))

    # store path to cluster bit
    outbit = data.tmpdir / f"aligned_{idx}.fa"

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
        print(rdict)
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
    with open(outbit, 'w', encoding="utf-8") as outfile:
        if clusts:
            outfile.write("\n//\n//\n".join(clusts) + "\n//\n//\n")
    alignments.close()


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s6_clustmap_across import Step6

    ip.set_log_level("DEBUG", "/tmp/test.log")

    data = ip.load_json("/tmp/RICHIE.json")
    data = data.branch("RICH2", subsample=list(data.samples)[:10])
    # print(data.stats)
    with ip.Cluster(cores=2) as ipyclient:
        step = Step6(data, True, False, ipyclient)
        step.run()
        # concat_bams(step.data, step.samples)
