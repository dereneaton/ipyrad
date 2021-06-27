#!/usr/bin/env python

"""
Reference assembly related functions. Some of these are very similar
to the functions in the denovo_utils but they sometimes differ subtly
and so it was easier to write seperate versions for denovo and ref.
"""

# pylint: disable=no-member, 

from typing import List
import os
import sys
import gzip
import subprocess as sps
import numpy as np
import pysam
from ipyrad.assemble.utils import IPyradError


BIN_BWA = os.path.join(sys.prefix, "bin", "bwa")
BIN_SAMTOOLS = os.path.join(sys.prefix, "bin", "samtools")
BIN_VSEARCH = os.path.join(sys.prefix, "bin", "vsearch")
BIN_BEDTOOLS = os.path.join(sys.prefix, "bin", "bedtools")


def index_ref_with_bwa(data, alt=False):
    """
    Index the reference sequence, unless it already exists
    """
    # get ref file from params, alt ref is for subtraction
    if not alt:
        refseq_file = data.params.reference_sequence
    else:
        refseq_file = data.params.reference_as_filter

    # check that ref file exists
    if not os.path.exists(refseq_file):
        if not alt:
            raise IPyradError(
                f"Assembly method {data.params.assembly_method} "
                "requires that a valid reference_sequence_path. "
                f"You entered {data.params.reference_sequence}")
        raise IPyradError(
            "reference_as_filter requires that you enter a reference "
            "fasta file. The path you entered was not found: "
            f"{data.params.reference_as_filter}")

    # If reference sequence already exists then bail out of this func
    index_files = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if all([os.path.isfile(refseq_file + i) for i in index_files]):
        print("reference is bwa indexed: {}".format(refseq_file))
        return

    # bwa index <reference_file>
    cmd = [BIN_BWA, "index", refseq_file]
    print(" ".join(cmd))
    proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=None)
    error = proc.communicate()[1].decode()

    # error handling for one type of error on stderr
    if proc.returncode:
        if "please use bgzip" in error:
            raise IPyradError((
                "Reference sequence must be de-compressed fasta or bgzip "
                "compressed, your file is probably gzip compressed. The "
                "simplest fix is to gunzip your reference sequence by "
                "running this command: \n"
                "    gunzip {}\n"
                "Then edit your params file to remove the `.gz` from the "
                "end of the path to your reference sequence file and rerun "
                "step 3 with the `-f` flag."
                .format(refseq_file)))
        raise IPyradError(error)


def index_ref_with_sam(data, alt=False):
    """
    Index ref for building scaffolds w/ index numbers in steps 5-6
    """
    # get ref file from params, alt ref is for subtraction
    if not alt:
        refseq_file = data.params.reference_sequence
    else:
        refseq_file = data.params.reference_as_filter

    # check whether it is already indexed
    if not os.path.exists(refseq_file):
        if not alt:
            raise IPyradError((
                "Assembly method {} requires that you enter a "
                "reference_sequence_path. The path you entered was not "
                "found: \n{}")
                .format(data.params.assembly_method, data.params.reference_sequence))
        raise IPyradError((
            "reference_as_filter requires that you enter a reference "
            "fasta file. The path you entered was not found: \n{}")
            .format(data.params.reference_as_filter))

    # If reference index exists then bail out unless force
    if os.path.exists(refseq_file + ".fai"):
        print(f"reference is sam indexed: {refseq_file}")
        return

    # complain if file is bzipped
    if refseq_file.endswith(".gz"):
        raise IPyradError("You must decompress your genome file.") 

    # index the file
    print(f"indexing {refseq_file} with pysam/samtools")
    pysam.faidx(refseq_file)


def mapping_reads_minus(data, sample, nthreads):
    """
    Map reads to the reference-filter fasta to get unmapped fastq
    files to use for downstream analyses.
    """
    if not data.params.reference_as_filter:
        return
    reference = data.params.reference_as_filter
    if not os.path.exists(reference):
        raise IPyradError(f"reference_filter sequence not found: {reference}")

    # input reads are concat if present else trims
    read1 = os.path.join(data.tmpdir, f"{sample.name}_concat_R1.fastq.gz")
    read2 = os.path.join(data.tmpdir, f"{sample.name}_concat_R2.fastq.gz")
    read1 = read1 if os.path.exists(read1) else sample.files.edits[0][0]
    read2 = read2 if os.path.exists(read2) else sample.files.edits[0][1]

    # setup cmd1 (mapping w/ bwa)
    cmd1 = [BIN_BWA, "mem", "-t", str(max(1, nthreads)), "-M", reference]
    cmd1 += [read1, read2 if read2 else ""]
    if data.hackers.bwa_args:
        for arg in data.hackers.bwa_args.split()[::-1]:
            cmd1.insert(2, arg)
    cmd1 += ['-o', os.path.join(data.tmpdir, f"{sample.name}.sam")]
    print(" ".join(cmd1))

    # run cmd1
    proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.DEVNULL)
    error1 = proc1.communicate()[1]
    if proc1.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd1)}\nerror: {error1}")

    # setup cmd2 (sam to bam)
    cmd2 = [BIN_SAMTOOLS, 'view', '-b', '-F', '0x904']
    cmd2 += ['-U', os.path.join(data.tmpdir, f"{sample.name}.unmapped.bam")]
    cmd2 += [os.path.join(data.tmpdir, f"{sample.name}.sam")]
    print(' '.join(cmd2))

    # run cmd2
    proc2 = sps.Popen(cmd2, stderr=sps.PIPE, stdout=sps.DEVNULL)
    error2 = proc2.communicate()[1]
    if proc2.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd2)}\nerror: {error2}")

    # setup cmd3 (bam to fastq unmapped)
    cmd3 = [BIN_SAMTOOLS, 'fastq', '-v', '45']
    cmd3 += ['-1', os.path.join(data.tmpdir, f"{sample.name}.unmapped_R1.fastq")]
    cmd3 += ['-2', os.path.join(data.tmpdir, f"{sample.name}.unmapped_R2.fastq")]
    cmd3 += [os.path.join(data.tmpdir, f"{sample.name}.unmapped.bam")]
    print(' '.join(cmd3))

    # run cmd3
    proc3 = sps.Popen(cmd3, stderr=sps.PIPE, stdout=sps.DEVNULL)
    error3 = proc3.communicate()[1]
    if proc3.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd3)}\nerror: {error3}")


def join_pairs_for_derep_ref(data, sample):
    """ 
    Combines read1 and read2 with a 'nnnn' separator. Because reads
    will be refmapped we DO NOT REVCOMP read2 in 'reference' datatype.
    The joined read pairs are used for dereplication and decloning.
    """
    # input file options
    unmapped1 = os.path.join(data.tmpdir, f"{sample.name}.unmapped_R1.fastq")
    unmapped2 = os.path.join(data.tmpdir, f"{sample.name}.unmapped_R2.fastq")
    concat1 = os.path.join(data.tmpdir, f"{sample.name}.concat_R1.fastq.gz")
    concat2 = os.path.join(data.tmpdir, f"{sample.name}.concat_R2.fastq.gz")
    trim1 = sample.files.edits[0][0]
    trim2 = sample.files.edits[0][1]

    # file precedence
    order1 = (trim1, concat1, unmapped1)
    order2 = (trim2, concat2, unmapped2)
    read1 = [i for i in order1 if os.path.exists(i)][-1]
    read2 = [i for i in order2 if os.path.exists(i)][-1]

    # read in paired end read files 4 lines at a time
    xopen = gzip.open if read1.endswith(".gz") else open
    readio1 = xopen(read1, 'rt')
    readio2 = xopen(read2, 'rt')
    quart1 = zip(*[readio1] * 4)
    quart2 = zip(*[readio2] * 4)
    quarts = zip(quart1, quart2)

    # output file
    out = open(os.path.join(data.tmpdir, f"{sample.name}_joined.fastq"), 'wt')

    # a list to store until writing
    writing = []
    counts = 0

    # iterate until done
    while 1:
        try:
            read1s, read2s = next(quarts)
        except StopIteration:
            break

        writing.append(
            "".join([
                read1s[0],
                read1s[1].strip() + "nnnn" + read2s[1],
                read1s[2],
                read1s[3].strip() + "nnnn" + read2s[3],
            ])
        )

        # keep count
        counts += 1
        if not counts % 5000:
            out.write("".join(writing))
            writing = []
    if writing:
        out.write("".join(writing))

    # close handles
    readio1.close()
    readio2.close()
    out.close()
    print(f"joined read pairs of {sample.name}")


def tag_for_decloning(data, sample):
    """
    Pull i5 tag from Illumina index and insert into the fastq name 
    header so we can use it later to declone even after the fastq 
    indices are gone. THIS IS DONE BEFORE DEREPLICATION, so that 
    identical reads are not collapsed if they have different i5s.

    # e.g.,                                            i7       i5
    @NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    +
    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

    to

    @NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA    
    CAAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    +
    FFFFFFFFBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

    3rad uses random adapters to identify pcr duplicates. For removing
    pcr dups later we need to insert the tag into the sequence here
    so they will not be dereplicated together.
    """
    # paired reads are merged or joined in the merged file
    tmpin = os.path.join(data.tmpdir, f"{sample.name}_joined.fastq")
    tmpout = os.path.join(data.tmpdir, f"{sample.name}_declone.fastq")

    # Remove adapters from head of sequence and write out
    # tmp_outfile is now the input file for the next step
    # first vsearch derep discards the qscore so we iterate pairs
    outfile = open(tmpout, 'wt') 
    with open(tmpin, 'rt') as infile:

        # iterate over 2 lines a time
        duo = zip(*[infile] * 4)
        writing = []
        counts2 = 0

        # a list to store until writing
        while 1:
            try:
                read = next(duo)
            except StopIteration:
                break

            # extract i5 if it exists else use empty string           
            try:
                indices = read[0].split(":")[-1]
                index = indices.split("+")[1].strip()
                assert len(index) == 8
            except AssertionError:
                index = ""
            fake_quality = "B" * len(index)

            # add i5 to the 5' end of the sequence
            newread = [
                read[0],
                index + read[1],
                read[2],
                fake_quality + read[3]
            ]
            writing.append("".join(newread))

            # Write the data in chunks
            counts2 += 1
            if not counts2 % 5000:
                outfile.write("".join(writing))
                writing = []

        if writing:
            outfile.write("".join(writing))
    outfile.close()
    print(f"tagged for decloning: {sample.name}")


def dereplicate_func(data, sample):
    """
    Dereplicates reads and sorts so reads that were highly replicated 
    are at the top, and singletons at bottom, writes output to derep 
    file. Paired data are already joined into a single concat file.
    """
    # file precedence 
    infiles = [
        sample.files.edits[0][0],
        os.path.join(data.tmpdir, f"{sample.name}_concat_R1.fastq.gz"),
        os.path.join(data.tmpdir, f"{sample.name}_joined.fastq"),
        os.path.join(data.tmpdir, f"{sample.name}_declone.fastq"),
    ]
    infile = [i for i in infiles if os.path.exists(i)][-1]

    # datatypes options
    strand = "plus"
    if data.params.datatype in ['gbs', '2brad']:
        strand = "both"

    # do dereplication with vsearch
    cmd = [
        BIN_VSEARCH,
        "--derep_fulllength", infile,
        "--strand", strand,
        "--output", os.path.join(data.tmpdir, f"{sample.name}_derep.fa"),
        "--fasta_width", str(0),
        "--minseqlength",  str(data.params.filter_min_trim_len),
        "--sizeout", 
        "--relabel_md5",
        "--quiet",
        # "--threads", str(nthreads),
        #"--fastq_qmax", "1000",        
    ]

    # decompress argument (IF ZLIB is missing this will not work!!) 
    # zlib is part of the conda installation.
    if infile.endswith(".gz"):
        cmd.append("--gzip_decompress")

    # build PIPEd job
    print(" ".join(cmd))
    proc = sps.Popen(cmd, stderr=sps.PIPE, close_fds=True)
    errmsg = proc.communicate()[1]
    if proc.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd)}\nerror: {errmsg}")


def tag_to_header_for_decloning(data, sample):
    """
    Pull i5 tag from dereplicated sequence and append to (PE joined) 
    dereped seq header.

    # e.g.,      
    0004a51ebbcd442afb6b6d02f5daf553;size=6;*
    TATCGGTCATCGGCTAAGAATAAGAGAAAAAAACAAGTGAATGATAATGAATATGGATATGACTAAAA

    to 

    0004a51ebbcd442afb6b6d02f5daf553;tag=TATCGGTC;size=6;*
    ATCGGCTAAGAATAAGAGAAAAAAACAAGTGAATGATAATGAATATGGATATGACTAAAA
    """
    tmpin = os.path.join(data.tmpdir, f"{sample.name}_derep.fa")
    tmpout = os.path.join(data.tmpdir, f"{sample.name}_tagged.fa")

    # Remove adapters from head of sequence and write out
    # tmp_outfile is now the input file for the next step
    # first vsearch derep discards the qscore so we iterate pairs
    outfile = open(tmpout, 'wt') 
    with open(tmpin, 'rt') as infile:

        # iterate over 2 lines a time
        duo = zip(*[infile] * 2)
        writing = []
        counts2 = 0

        # a list to store until writing
        while 1:
            try:
                read = next(duo)
            except StopIteration:
                break

            # extract i5 if it exists else use empty string           
            tag = read[1][:8]
            name, size = read[0].split(";")

            # add i5 to the 5' end of the sequence
            newread = [
                "{};tag={};{}".format(name, tag, size),
                read[1][8:],
            ]
            writing.append("".join(newread))

            # Write the data in chunks
            counts2 += 1
            if not counts2 % 1000:
                outfile.write("".join(writing))
                writing = []

        if writing:
            outfile.write("".join(writing))
            outfile.close()


def split_derep_pairs_ref(data, sample):
    """
    Takes R1nnnnR2 derep reads from paired data and splits it back into
    separate R1 and R2 parts for read mapping.
    """
    # select tagged if it exists else choose derep
    infiles = [
        os.path.join(data.tmpdir, f"{sample.name}_tagged.fa"),
        os.path.join(data.tmpdir, f"{sample.name}_derep.fa"),
    ]
    infile = [i for i in infiles if os.path.exists(i)][0]

    # output fasta names
    out1 = os.path.join(data.tmpdir, f"{sample.name}_derep_split_R1.fa")
    out2 = os.path.join(data.tmpdir, f"{sample.name}_derep_split_R2.fa")
    print(f"resplitting paired reads to {sample.name}")

    # open files for writing
    splitderep1 = open(out1, 'wt')
    splitderep2 = open(out2, 'wt')

    with open(infile, 'rt') as infile:
        # Read in the infile two lines at a time: (seqname, sequence)
        duo = zip(*[infile] * 2)

        # lists for storing results until ready to write
        split1s = []
        split2s = []

        # iterate over input splitting, saving, and writing.
        idx = 0
        while 1:
            try:
                itera = next(duo)
            except StopIteration:
                break

            # split the duo into separate parts and inc counter
            try:
                part1, part2 = itera[1].split("nnnn")
            except ValueError:
                print("SPLIT THIS {}".format(itera))
            idx += 1

            # R1 needs a newline, but R2 inherits it from the original file
            # store parts in lists until ready to write
            split1s.append("{}{}\n".format(itera[0], part1))
            split2s.append("{}{}".format(itera[0], part2))

            # if large enough then write to file
            if not idx % 10000:
                splitderep1.write("".join(split1s))
                splitderep2.write("".join(split2s))
                split1s = []
                split2s = []

    # write final chunk if there is any
    if any(split1s):
        splitderep1.write("".join(split1s))
        splitderep2.write("".join(split2s))

    # close handles
    splitderep1.close()
    splitderep2.close()


def mapping_reads(data, sample, nthreads):
    """
    Map reads to the reference to get a sorted bam.
    """
    reference = data.params.reference_sequence
    if not os.path.exists(reference):
        raise IPyradError(f"reference_sequence not found: {reference}")

    # input reads are concat if present else trims
    read1 = os.path.join(data.tmpdir, f"{sample.name}_derep_split_R1.fa")
    read2 = os.path.join(data.tmpdir, f"{sample.name}_derep_split_R2.fa")
    if not os.path.exists(read1):
        read1 = os.path.join(data.tmpdir, f"{sample.name}_derep.fa")
        read2 = None
    assert os.path.exists(read1), f"derep files missing: {read1}"

    # setup cmd1 (mapping w/ bwa)
    cmd1 = [BIN_BWA, "mem", "-t", str(max(1, nthreads)), "-M", reference]
    cmd1.append(read1)
    if read2:
        cmd1.append(read2)
    if data.hackers.bwa_args:
        for arg in data.hackers.bwa_args.split()[::-1]:
            cmd1.insert(2, arg)
    cmd1 += ['-o', os.path.join(data.tmpdir, f"{sample.name}.sam")]
    print(" ".join(cmd1))

    # run cmd1
    proc1 = sps.Popen(cmd1, stderr=sps.PIPE, stdout=sps.DEVNULL)
    error1 = proc1.communicate()[1]
    if proc1.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd1)}\nerror: {error1}")

    # setup cmd2 (sam to bam)
    cmd2 = [BIN_SAMTOOLS, 'view', '-b', '-F', '0x904']
    cmd2 += [os.path.join(data.tmpdir, f"{sample.name}.sam")]
    print(' '.join(cmd2))

    # run cmd2
    proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE)

    # setup cmd3 (sort sam)
    cmd3 = [
        BIN_SAMTOOLS, "sort", 
        "-T", os.path.join(data.tmpdir, f"{sample.name}.sam.tmp"),
        "-O", "bam", 
        "-o", os.path.join(data.stepdir, f"{sample.name}.bam"),
    ]
    print(' '.join(cmd3))

    # run cmd3
    proc3 = sps.Popen(
        cmd3, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc2.stdout)
    error3 = proc3.communicate()[1]
    if proc3.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd3)}\nerror: {error3}")
    proc2.stdout.close()

    # setup cmd4 (index bam file)
    cmd4 = [BIN_SAMTOOLS, "index", "-c", os.path.join(data.stepdir, f"{sample.name}.bam")]

    # run cmd4
    proc4 = sps.Popen(cmd4, stderr=sps.PIPE, stdout=sps.DEVNULL)
    error4 = proc4.communicate()[1]
    if proc4.returncode:
        raise IPyradError(f"cmd: {' '.join(cmd4)}\nerror: {error4}")


def build_clusters_from_cigars(data, sample):
    """
    Directly building clusters relative to reference. Uses the function 
    cigared() to impute indels relative to reference. This means add - 
    for insertion and skip* deletions. Skipping is not a good final 
    solution.

    If a read pair is merged completely then we mask the restriction
    overhangs from each end. This is because the read-through can
    recover the uncut sequence in SOME seqs, leading to false SNPS.

                ------ATCGGAGGTA ...  not merged
                ATTCACCG--------      fully merged
                ATTCACCG--------      fully merged
                ATTCACCG--------      fully merged
                ------ATCGGAGGTA ...  not merged
                ------ATCGGAGGTA ...  not merged

                ATTCACNNNNNAGGTA  ==  masked tips for merged pairs
 
    """
    # get all regions with reads as a List of 'scaff\tstart\tend' strs
    regions = bedtools_merge(data, sample)

    # if no matched regions (no clusters) then bail out.
    if not regions:
        return 

    # Make into a generator of tuples
    # bedtools is 0-based. sam is 1-based, so increment the positions here
    regions_split = (i.split("\t") for i in regions)
    regions = ((i, int(j) + 1, int(k) + 1) for (i, j, k) in regions_split)

    # access reads from bam file using pysam
    bamfile = os.path.join(data.stepdir, f"{sample.name}.bam")
    bamin = pysam.AlignmentFile(bamfile, 'rb')

    # output path 
    opath = os.path.join(data.stepdir, f"{sample.name}.clusters.gz")
    out = gzip.open(opath, 'wt')
    idx = 0

    # iterate over all regions to build clusters
    clusters = []
    for reg in regions:
        # uncomment and compare against ref sequence when testing
        # ref = get_ref_region(data.paramsdict["reference_sequence"], *reg)
        reads = bamin.fetch(*reg)

        # store reads in a dict
        rdict = {}

        # paired-end data cluster building
        if data.is_pair:

            # match paired reads together in a dictionary (same qname)
            # for each mapped read in this region.
            for read in reads:
                if read.qname not in rdict:
                    rdict[read.qname] = [read, None]
                else:
                    rdict[read.qname][1] = read

            # sort keys by derep number
            keys = sorted(
                rdict.keys(),
                key=lambda x: int(x.split("=")[-1]), reverse=True)

            # get max RE length for tip-masking
            try:
                tipmask = max([len(i) for i in data.params.restriction_overhang])
            except Exception:
                tipmask = 6

            # build the cluster based on map positions, orientation, cigar
            clust = []
            for key in keys:

                # both pairs must match or we exclude the reads
                align1, align2 = rdict[key]
                if align1 and align2:

                    # create two empty arrays length of the region
                    lref = reg[2] - reg[1]
                    arr1 = np.zeros(lref, dtype="U1")
                    arr2 = np.zeros(lref, dtype="U1")
                    arr1.fill("-")
                    arr2.fill("-")

                    # how far ahead of the start does this read begin
                    seq = cigared(align1.seq, align1.cigar)
                    start1 = align1.reference_start - (reg[1] - 1) 
                    len1 = len(seq)
                    arr1[start1:start1 + len1] = list(seq)
                    
                    seq = cigared(align2.seq, align2.cigar)
                    start2 = align2.reference_start - (reg[1] - 1)
                    len2 = len(seq)
                    arr2[start2:start2 + len2] = list(seq)

                    # tip-masking for merged pairs.
                    if start2 <= (start1 + tipmask):
                        arr2[start1:start1 + tipmask] = "-"
                    if (start1 + len1) >= (start2 + len2 - tipmask):
                        arr1[start2+len2-tipmask:start2+len2] = "-"
                    if start1 != 0:
                        arr1[start1:start1 + tipmask] = "-"
                    if start2 + len2 != lref:
                        arr2[start2+len2-tipmask:start2+len2] = "-"

                    arr3 = join_arrays(arr1, arr2)
                    pairseq = "".join(arr3)

                    ori = "+"
                    if align1.is_reverse:
                        ori = "-"
                    derep = align1.qname.split("=")[-1]

                    # insert decloning tag into name, or not.
                    if data.hackers.declone_PCR_duplicates:
                        tag = align1.qname.split(";")[-2]
                        rname = "{}:{}-{};{};size={};{}".format(
                            reg[0], reg[1], reg[2], tag, derep, ori,
                        )
                    else:
                        rname = "{}:{}-{};size={};{}".format(
                            reg[0], reg[1], reg[2], derep, ori,
                        )
                    clust.append("{}\n{}".format(rname, pairseq))


        # single-end data cluster building
        else:   
            mstart = int(9e12)
            mend = 0

            for read in reads:
                rdict[read.qname] = read
                mstart = min(mstart, read.aend - read.alen)
                mend = max(mend, read.aend)

            # sort keys by derep number
            keys = sorted(
                rdict.keys(),
                key=lambda x: int(x.split("=")[-1]), reverse=True)

            # build the cluster based on map positions, orientation, cigar
            clust = []
            for key in keys:
                align1 = rdict[key]

                #aref = np.array(list(ref[1]))
                lref = mend - mstart
                arr1 = np.zeros(lref, dtype="U1")
                arr1.fill("-")

                # how far ahead of the start does this read begin
                seq = cigared(align1.seq, align1.cigar)
                rstart = (align1.aend - align1.alen) - mstart
                arr1[rstart:rstart + len(seq)] = list(seq)
                aseq = "".join(arr1)

                ori = "+"
                if align1.is_reverse:
                    ori = "-"
                derep = align1.qname.split("=")[-1]
                # Pysam coords are 0 based, but sam files are 1 based, and
                # since we we build sam in step 5, we need to account for the
                # diffrent indexing strategies here by incrementing mstart
                # and mend
                rname = "{}:{}-{};size={};{}".format(
                    reg[0], mstart + 1, mend + 1, derep, ori)
                clust.append("{}\n{}".format(rname, aseq))

        # store this cluster
        if clust:
            clusters.append("\n".join(clust))
            idx += 1

        # if 1000 clusters stored then write to disk
        if not idx % 1000:
            if clusters:
                out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
                clusters = []

    # write final remaining clusters to disk
    if clusters:
        out.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
    out.close()


def join_arrays(arr1, arr2):
    """
    Fastp already performed paired-read correction for overlapping 
    reads. But if diffs persist at overlap somehow then mask them here.
    Also mask bases at cutsites if no insert is present.
    """
    arr3 = np.zeros(arr1.size, dtype="U1")
    for i in range(arr1.size):
        if arr1[i] == arr2[i]:
            arr3[i] = arr1[i]

        elif arr1[i] == "N":
            if arr2[i] == "-":
                arr3[i] = "N"
            else:
                arr3[i] = arr2[i]

        elif arr2[i] == "N":
            if arr1[i] == "-":
                arr3[i] = "N"
            else:
                arr3[i] = arr1[i]

        elif arr1[i] == "-":
            if arr2[i] == "N":
                arr3[i] = "N"
            else:
                arr3[i] = arr2[i]

        elif arr2[i] == "-":
            if arr1[i] == "N":
                arr3[i] = "N"
            else:
                arr3[i] = arr1[i]
        else:
            arr3[i] = "N"
    return arr3


def cigared(sequence, cigartups) -> str:
    """
    modify sequence based on its cigar string
    """
    start = 0
    seq = ""
    for tup in cigartups:
        flag, add = tup
        if flag == 0:
            seq += sequence[start:start + add]
        if flag == 1:
            pass
        if flag == 2:
            seq += "-" * add
            start -= add
        if flag == 4:
            pass
        start += add
    return seq


def bedtools_merge(data, sample) -> List[str]:
    """
    Get all contiguous genomic regions with one or more overlapping
    reads. This is the shell command we'll eventually run

    bedtools bamtobed -i 1A_0.sorted.bam | bedtools merge [-d 100]
        -i <input_bam>  :   specifies the input file to bed'ize
        -d <int>        :   For PE set max distance between reads
    """
    mappedreads = os.path.join(data.stepdir, sample.name + ".bam")

    # command to call `bedtools bamtobed`, and pipe output to stdout
    # Usage:   bedtools bamtobed [OPTIONS] -i <bam>
    # Usage:   bedtools merge [OPTIONS] -i <bam>
    cmd1 = [BIN_BEDTOOLS, "bamtobed", "-i", mappedreads]
    cmd2 = [BIN_BEDTOOLS, "merge", "-i", "-"]

    # If PE the -d flag to tell bedtools how far apart to allow mate pairs.
    if data.is_pair:

        # uses the hackers dict value if user set one, else 500 as a 
        # reasonable max insert size.
        max_insert = data.hackers.max_inner_mate_distance
        if max_insert is None:
            max_insert = 500
        cmd2.insert(2, str(max_insert))
        cmd2.insert(2, "-d")

    # pipe output from bamtobed into merge
    proc1 = sps.Popen(cmd1, stderr=sps.STDOUT, stdout=sps.PIPE)
    proc2 = sps.Popen(cmd2, stderr=sps.STDOUT, stdout=sps.PIPE, stdin=proc1.stdout)
    result = proc2.communicate()[0].decode()
    proc1.stdout.close()

    # check for errors and do cleanup
    if proc2.returncode:
        raise IPyradError(f"error in {cmd2} {result}")

    # Report the number of regions we're returning
    regions = result.strip().split("\n")
    print(f"{sample.name}: max_insert={max_insert}, nregions={len(regions)}")
    return regions
