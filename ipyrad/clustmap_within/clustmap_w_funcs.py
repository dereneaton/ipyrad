#!/usr/bin/env python

"""Functions called by ClustMapBase.py

"""

from typing import Tuple
import sys
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT, DEVNULL
from loguru import logger
import pysam
from ipyrad.schema import Sample
from ipyrad.core import Assembly, IPyradError

BIN = Path(sys.prefix) / "bin"
BIN_BWA = str(BIN / "bwa")
BIN_SAMTOOLS = str(BIN / "samtools")  # indexing
BIN_VSEARCH = str(BIN / "vsearch")    # dereplicating


def dereplicate(data: Assembly, sample: Sample) -> None:
    """Dereplicate reads and sort so that most replicated are on top.

    Paired data are dereplicated as joined (dataNNNNdata) reads.
    """
    infiles = [
        Path(sample.files.trimmed[0][0]),
        data.tmpdir / f"{sample.name}_concat.fastq.gz",
        data.tmpdir / f"{sample.name}_merged.fa",
        data.tmpdir / f"{sample.name}_joined.fastq",
        data.tmpdir / f"{sample.name}_decloned.fa",
    ]
    infile = [i for i in infiles if i.exists()][-1]

    # only check both strands for identical sequences if a single cutter
    # is cutting at both ends of a sequence. Otherwise, RAD sequences
    # are not expected to match on both strands. This means that for
    # GBS and 2BRAD data the user must enter both REs even though the
    # data was prepared with only one.
    res = data.params.restriction_overhang
    if res[0] == res[1]:
        strand = "both"
    else:
        strand = "plus"

    # do dereplication with vsearch
    cmd = [
        BIN_VSEARCH,
        "--fastx_uniques", str(infile),
        "--strand", strand,
        "--fastaout", str(data.tmpdir / f"{sample.name}_derep.fa"),
        "--fasta_width", str(0),
        "--minseqlength", str(data.params.filter_min_trim_len),
        "--sizeout",
        "--relabel_md5",
        "--quiet",
        # "--threads", str(4), multithreading not supported
        # "--fastq_qmax", "1000",
    ]

    # decompress argument (IF ZLIB is missing this will not work!!)
    # zlib is part of the conda installation.
    if infile.suffix == ".gz":
        cmd.append("--gzip_decompress")

    # build PIPEd job
    print(f"@@DEBUG: {' '.join(cmd)}", flush=True)
    with Popen(cmd, stderr=STDOUT, stdout=PIPE, close_fds=True) as proc:
        errmsg = proc.communicate()[0].decode()
        if proc.returncode:
            print(f"@@ERROR: {' '.join(cmd)}\n{errmsg}")
            raise IPyradError(errmsg)

    # remove tmpfile IF IT IS NOT A TRIMMED FILE
    if infile.parent == data.tmpdir:
        infile.unlink()


def tag_inline_for_decloning(data: Assembly, sample: Sample) -> None:
    """Tag reads w/ i5 inline prior to dereplicating.

    Pull i5 tag from Illumina index and insert into the fastq name
    header so we can use it later to declone even after the fastq
    indices are gone. THIS IS DONE BEFORE DEREPLICATION, so that
    identical reads are not collapsed if they have different i5s.

    This obviously increases the time it takes to map reads by quite
    a bit, since identical reads are not collapsed. Too bad.

    3rad uses random adapters to identify pcr duplicates. For removing
    pcr dups later we need to insert the tag into the sequence here
    so they will not be dereplicated together.

    >>>                                                    i7       i5
    >>> >NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA
    >>> AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    the i5 is inserted to the read with "F" quality score

    >>> ********
    >>> >NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA
    >>> CAAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    """
    # The input files are different depending on ref or not.
    if data.is_ref:
        tmpin = data.tmpdir / f"{sample.name}_joined.fastq"
    else:
        tmpin = data.tmpdir / f"{sample.name}_merged.fa"
    tmpout = data.tmpdir / f"{sample.name}_decloned.fa"

    # Remove adapters from head of sequence and write out
    # tmp_outfile is now the input file for the next step
    # first vsearch derep discards the qscore so we iterate pairs
    with open(tmpout, 'w', encoding="utf-8") as out:
        with open(tmpin, 'r', encoding="utf-8") as infile:

            # iterate over 2 lines a time if fasta else 4 for fastq
            if data.is_ref:
                duo = zip(*[infile] * 4)
            else:
                duo = zip(*[infile] * 2)

            # a list to store until writing
            tmp = []
            for idx, read in enumerate(duo):

                # extract i5 if it exists else use empty string
                try:
                    indices = read[0].split(":")[-1]
                    index = indices.split("+")[1].strip()
                    # assert int(index), "failed to parse i5 index"
                except (IndexError, AssertionError):
                    index = ""

                # add i5 to the 5' end of the sequence and insert the
                # length of the i5 index (usually 8) into the header.
                newread = f">{read[0][1:]}{index}{read[1]}"
                tmp.append(newread)

                # Write the data in chunks
                if not idx % 50_000:
                    out.write("".join(tmp))
                    tmp = []
            if tmp:
                out.write("".join(tmp))
    print(f"@@DEBUG: {sample.name} - tagged with inline i5s for decloning")
    tmpin.unlink()


def map_to_reference_as_filter(data: Assembly, sample: Sample) -> Tuple[int, float]:
    """Map reads to the reference-filter fasta.

    Mapped reads are discarded and unmapped reads are kept as data.
    """
    reference = Path(data.params.reference_as_filter)

    # priority: select concats first if they exist, else trim files.
    read1 = data.tmpdir / f"{sample.name}_concat_R1.fastq.gz"
    read2 = data.tmpdir / f"{sample.name}_concat_R2.fastq.gz"
    read1 = read1 if read1.exists() else Path(sample.files.trimmed[0][0])
    read2 = read2 if read2.exists() else Path(sample.files.trimmed[0][1])

    # get threading arg
    nthreads = max(1, data.ipcluster["threads"])
    # if data.hackers.bwa_args:
    #     for arg in data.hackers.bwa_args.split()[::-1]:
    #         cmd1.insert(2, arg)

    ################################################################
    # setup cmd1 (mapping w/ bwa)
    cmd1 = [BIN_BWA, "mem", "-t", str(nthreads), "-M", str(reference)]
    cmd1 += [str(read1), str(read2) if read2 else ""]
    cmd1 += ['-o', str(data.tmpdir / f"{sample.name}_ref_filter.sam")]
    print(f"@@DEBUG: {' '.join(cmd1)}")

    # run cmd1
    with Popen(cmd1, stderr=PIPE, stdout=DEVNULL) as proc1:
        error1 = proc1.communicate()[1].decode()
        if proc1.returncode:
            print(f"@@ERROR: {' '.join(cmd1)}\n{error1}.")
            raise IPyradError(f"cmd: {' '.join(cmd1)}\nerror: {error1}")

    ################################################################
    # setup cmd2 (sam to bam)
    cmd2 = [BIN_SAMTOOLS, 'view', '-b', '-F', '0x904']
    cmd2 += ['-U', str(data.tmpdir / f"{sample.name}_unmapped.bam")]
    cmd2 += [str(data.tmpdir / f"{sample.name}_ref_filter.sam")]
    print(f"@@DEBUG: {' '.join(cmd2)}")

    # run cmd2
    with Popen(cmd2, stderr=PIPE, stdout=DEVNULL) as proc2:
        error2 = proc2.communicate()[1].decode()
        if proc2.returncode:
            print(f"@@ERROR: {' '.join(cmd2)}\n{error2}.")
            raise IPyradError(f"cmd: {' '.join(cmd2)}\nerror: {error2}")

    #################################################################
    # setup cmd3 (bam to fastq unmapped)
    cmd3 = [BIN_SAMTOOLS, 'fastq', '-v', '45']
    cmd3 += ['-1', str(data.tmpdir / f"{sample.name}_unmapped_R1.fastq")]
    cmd3 += ['-2', str(data.tmpdir / f"{sample.name}_unmapped_R2.fastq")]
    cmd3 += [str(data.tmpdir / f"{sample.name}_unmapped.bam")]
    print(f"@@DEBUG: {' '.join(cmd3)}")

    # run cmd3
    with Popen(cmd3, stderr=PIPE, stdout=DEVNULL) as proc3:
        error3 = proc3.communicate()[1].decode()
        if proc3.returncode:
            print(f"@@ERROR: {' '.join(cmd3)}\n{error3}.")
            raise IPyradError(f"cmd: {' '.join(cmd3)}\nerror: {error3}")

    # return the number of filtered reads
    unmapped = data.tmpdir / f"{sample.name}_unmapped_R1.fastq"
    with open(unmapped, 'r', encoding="utf-8") as ion:
        n_unmapped = int(sum(1 for i in ion) / 4)
    n_filtered = sample.stats_s1.reads_passed_filter - n_unmapped
    n_filtered_prop = n_filtered / sample.stats_s1.reads_passed_filter
    print(
        f"@@DEBUG: {sample.name} - proportion reads filtered by mapping "
        f"to reference filter: {n_filtered_prop:.3f}")
    return n_filtered, n_filtered_prop


def concat_multiple_fastqs_from_merged_sample(data: Assembly, sample: Sample) -> None:
    """Concat files if multiple Assemblies were merged between steps 1-2.

    Create a temporary concatenated file for multiple trim input
    files, which arises when Assemblies were merged between steps
    1 and 2.
    """
    # define output files
    concat1 = data.tmpdir / f"{sample.name}_concat_R1.fastq.gz"
    concat2 = data.tmpdir / f"{sample.name}_concat_R2.fastq.gz"

    read1s = [i[0] for i in sample.files.trimmed]
    if len(read1s) > 1:
        cmd = ['cat'] + read1s
        with open(concat1, 'w', encoding="utf-8") as cout:
            with Popen(cmd, stderr=PIPE, stdout=cout, close_fds=True) as proc:
                res = proc.communicate()[1].decode()
                if proc.returncode:
                    logger.error(res)
                    raise IPyradError(f"cmd: {' '.join(cmd)}\nerror: {res}")

    read2s = [i[1] for i in sample.files.trimmed if i[1]]
    if len(read2s) > 1:
        cmd = ['cat'] + read2s
        with open(concat2, 'w', encoding="utf-8") as cout:
            with Popen(cmd, stderr=PIPE, stdout=cout, close_fds=True) as proc:
                res = proc.communicate()[1].decode()
                if proc.returncode:
                    logger.error(res)
                    raise IPyradError(f"cmd: {' '.join(cmd)}\nerror: {res}")


def index_ref_with_bwa(data: Assembly, as_filter: bool = False) -> None:
    """Index the reference sequence, unless it already exists

    """
    # get ref file from params, alt ref is for subtraction
    if not as_filter:
        refseq_file = data.params.reference_sequence
    else:
        refseq_file = data.params.reference_as_filter

    # check that ref file exists
    if not refseq_file.exists():
        if not as_filter:
            raise IPyradError(
                f"You must enter a valid reference sequence in params. "
                f"You entered {data.params.reference_sequence}")
        raise IPyradError(
            "reference_as_filter requires that you enter a reference "
            "fasta file. The path you entered was not found: "
            f"{data.params.reference_as_filter}")

    # If reference sequence already exists then bail out of this func
    suffs = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    # don't use Path.with_suffix here b/c '.fa.ann' double suffix is messy.
    if all(Path(str(refseq_file) + i).exists() for i in suffs):
        print(f"reference is bwa indexed: {refseq_file}")
        return

    # bwa index <reference_file>
    cmd = [str(BIN_BWA), "index", str(refseq_file)]
    print(f"@@DEBUG: {' '.join(cmd)}")
    with Popen(cmd, stderr=PIPE, stdout=None) as proc:
        error = proc.communicate()[1].decode()

    # error handling for one type of error on stderr
    if proc.returncode:
        if "please use bgzip" in error:
            raise IPyradError(
                "Reference sequence must be de-compressed fasta or bgzip "
                "compressed, your file is probably gzip compressed. The "
                "simplest fix is to gunzip your reference sequence by "
                "running this command: \n"
                f"    gunzip {refseq_file}\n"
                "Then edit your params file to remove the `.gz` from the "
                "end of the path to your reference sequence file and rerun "
                "step 3 with the `-f` flag.")
        raise IPyradError(error)


def index_ref_with_sam(data: Assembly, as_filter: bool = False) -> None:
    """Index ref for building scaffolds w/ index numbers in steps 5-6.

    """
    # get ref file from params, alt ref is for subtraction
    if not as_filter:
        refseq_file = data.params.reference_sequence
    else:
        refseq_file = data.params.reference_as_filter

    # check whether it is already indexed
    if not refseq_file.exists():
        if not as_filter:
            raise IPyradError(
                "You must enter a valid reference_sequence_path. You "
                f"entered: {data.params.reference_sequence}")
        raise IPyradError(
            "reference_as_filter requires that you enter a reference "
            "fasta file. The path you entered was not found:\n"
            f"{data.params.reference_as_filter}")

    # If reference index exists then bail out unless force
    if Path(str(refseq_file) + ".fai").exists():
        print(f"@@DEBUG: reference_as_filter is sam indexed: {refseq_file}")
        return

    # complain if file is bzipped
    if refseq_file.suffix == ".gz":
        msg = f"You must decompress the reference_as_filter file: {refseq_file}."
        print(f"@@ERROR: {msg}")
        raise IPyradError(msg)

    # index the file
    print(f"@@DEBUG: indexing {refseq_file} with pysam/samtools")
    pysam.faidx(str(refseq_file))


if __name__ == "__main__":
    pass
