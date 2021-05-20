#!/usr/bin/env python

"""
Reference based clustering across samples based on mapped positions
to a reference genome.
"""

import os
import sys
from loguru import logger
import pysam
from ipyrad.assemble.utils import IPyradError, fullcomp, chroms2ints
from ipyrad.core.progress_bar import AssemblyProgressBar


BIN_SAMTOOLS = os.path.join(sys.prefix, "bin", "samtools")
BIN_BEDTOOLS = os.path.join(sys.prefix, "bin", "bedtools")


class ClustMapAcrossReference:
    def __init__(self, step):
        pass


    def run(self):

        # concat
        self.remote_concat_bams()

        # get extents of regions using bedtools merge
        self.remote_build_ref_regions()

        # build clusters from regions
        self.remote_build_ref_clusters()

        # concat aligned files (This is not necessary, chunk again in s7)
        self.concat_alignments()


    def remote_concat_bams(self):
        """
        Merge bam files into a single large sorted indexed bam
        """
        # concatenate consens bamfiles for all samples in this assembly
        catbam = os.path.join(self.data.stepdir, f"{self.data.name}.cat.bam")
        cmd1 = [
            BIN_SAMTOOLS, 
            "merge", 
            "-f", 
            catbam,
        ]

        # Use the sample.files.consens info, rather than data.dirs to allow
        # for merging assemblies after step 5 where data.dirs is invalid/empty.
        for sample in self.samples:
            cmd1.append(sample.files.consens)

        proc = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)

        # progress bar
        while proc.poll() != 0:
            self.data._progressbar(3, 0, start, printstr)
            time.sleep(0.1)

        # parse result
        err = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(
                "error in: {}: {}".format(" ".join(cmd1), err))

        # sort the bam file
        cmd2 = [
            BIN_SAMTOOLS,
            "sort",
            "-T",
            catbam + '.tmp',
            "-o", 
            os.path.join(self.data.stepdir, f"{self.data.name}.cat.sorted.bam"),
            catbam,
        ]
        proc = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)

        # progress bar
        while not proc.poll() == 0:
            self.data._progressbar(3, 1, start, printstr)
            time.sleep(0.1)

        # parse result
        err = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(
                "error in: {}: {}".format(" ".join(cmd2), err))
        os.remove(catbam)

        # index the bam file
        cmd3 = [
            BIN_SAMTOOLS,
            "index", 
            os.path.join(self.data.stepdir, f"{self.data.name}.cat.sorted.bam",
            ),           
        ]
        proc = sps.Popen(cmd3, stderr=sps.STDOUT, stdout=sps.PIPE)

        # progress bar
        while not proc.poll() == 0:
            self.data._progressbar(3, 2, start, printstr)
            time.sleep(0.1)

        # parse result
        err = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(
                "error in: {}: {}".format(" ".join(cmd3), err))
        self.data._progressbar(3, 3, start, printstr)
        self.data._print("")


    def remote_build_ref_regions(self):
        """
        call bedtools remotely and track progress
        """
        msg = "fetching regions"
        jobs = {0: self.ipyclient[0].apply(build_ref_regions, self.data)}
        prog = AssemblyProgressBar(jobs, msg, 6, self.step.quiet)
        prog.block()
        prog.check()
        self.regions = rasync.get()
        logger.debug('regions: {}...'.format(self.regions[:10]))


    def remote_build_ref_clusters(self):
        """
        build clusters and find variants/indels to store
        """       
        # send N jobs each taking chunk of regions
        ncpus = self.data.ncpus
        nloci = len(self.regions)
        optim = int((nloci // ncpus) + (nloci % ncpus))
        optim = int(np.ceil(optim / 2))

        # send jobs to func
        printstr = ("building database   ", "s6")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()
        prog.jobs = {}
        for idx, chunk in enumerate(range(0, nloci, optim)):
            region = self.regions[chunk: chunk + optim]
            if region:
                args = (self.data, idx, region)
                prog.jobs[idx] = self.lbview.apply(build_ref_clusters, *args)

        # print progress while bits are aligning
        prog.block()
        prog.check()



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
    for sidx, key in enumerate(keys):
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
        os.path.join(
            data.dirs.across,
            "{}.cat.sorted.bam".format(data.name)
        )
    ]

    cmd2 = [
        ipyrad.bins.bedtools, 
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
        raise IPyradError(
            "error in {}: {}".format(" ".join(cmd2), result))
    regs = [i.split("\t") for i in result.strip().split("\n")]
    return [(i, int(j), int(k)) for i, j, k in regs]



def build_ref_clusters(data, idx, iregion):
    """
    Given a chunk of regions this will pull in the reference for each region
    and then pull in all consens reads matching to that region. It uses cigar
    info to align the consens reads with the ref. This also merges consens
    from the same sample that were not merged earlier, which is why we expect
    no duplicate samples in the output of reference assemblies.
    """

    # prepare i/o for bamfile with mapped reads
    bamfile = AlignmentFile(
        os.path.join(
            data.dirs.across,
            "{}.cat.sorted.bam".format(data.name)),
        'rb')

    # dict to map chromosome names to integers
    faidict = chroms2ints(data, False)

    # prepare i/o for pysam reference indexed
    reffai = FastaFile(data.params.reference_sequence)

    # store path to cluster bit
    outbit = os.path.join(data.tmpdir, "aligned_{}.fa".format(idx))

    # get clusters
    iregions = iter(iregion)
    clusts = []

    while 1:
        # pull in all consens reads mapping to a bed region
        try:
            region = next(iregions)
            reads = bamfile.fetch(*region)
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
        for idx, key in enumerate(keys):
            seq, cigar, start, end = rdict[key]

            # how far ahead of ref start and short of ref end is this read
            fidx = start - mstart
            eidx = arr.shape[1] - (mend - end)

            # enter into the array, trim end if longer than pulled ref
            arr[idx + 1, fidx:eidx] = list(seq)[:eidx - fidx]

            # mod sequence according to cigar for indels and ambigs
            # csums is the location of impute on the seq, so it must be 
            # incremented by fidx and not extend past eidx
            for cidx, cig in enumerate(cigar):
                if cig[0] == 4:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[idx + 1, csums] = arr[idx + 1, csums].lower()
                if cig[0] == 1:
                    csums = sum(i[1] for i in cigar[:cidx])
                    csums += eidx
                    if csums < fidx:
                        arr[idx + 1, csums] = b"-"

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
        for idx, key in enumerate(keys):    
            clust.append(
                ">{}\n{}".format(key, b"".join(arr[idx + 1]).decode())
            )
        clusts.append("\n".join(clust))

    # dump to temp file until concat in next step.
    with open(outbit, 'w') as outfile:
        if clusts:
            outfile.write("\n//\n//\n".join(clusts) + "\n//\n//\n")

