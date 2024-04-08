#!/usr/bin/env python

"""Functions called on remote engines by clustmap_w_denovo.py

"""

from typing import TypeVar, Tuple
import sys
from pathlib import Path
from subprocess import Popen, PIPE, STDOUT

from loguru import logger
from ipyrad.core import IPyradError
from ipyrad.core.utils import comp

# pylint: disable=too-many-branches, too-many-statements, consider-using-with, too-many-lines

logger = logger.bind(name="ipyrad")

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")

BIN = Path(sys.prefix) / "bin"
BIN_VSEARCH = str(BIN / "vsearch")
# BIN_BWA = str(BIN / "bwa")
# BIN_SAMTOOLS = str(BIN / "samtools")


def merge_pairs_with_vsearch(data: Assembly, sample: Sample) -> int:
    """Merge PE reads using vsearch to find overlap."""
    # input files (select only the top one)
    in1 = [
        data.tmpdir / f"{sample.name}_unmapped_R1.fastq",
        data.tmpdir / f"{sample.name}_concat_R1.fq.gz",
        sample.files.trimmed[0],
    ]
    in2 = [
        data.tmpdir / f"{sample.name}_unmapped_R2.fastq",
        data.tmpdir / f"{sample.name}_concat_R2.fq.gz",
        sample.files.trimmed[1],
    ]
    index = min([i for i, j in enumerate(in1) if Path(j).exists()])
    infile1 = in1[index]
    infile2 = in2[index]

    # define output files
    mergedfile = data.tmpdir / f"{sample.name}_merged.fa"
    nonmerged1 = data.tmpdir / f"{sample.name}_nonmerged_R1.fa"
    nonmerged2 = data.tmpdir / f"{sample.name}_nonmerged_R2.fa"

    # get the maxn and minlen values
    try:
        maxn = sum(data.params.max_low_qual_bases)
    except TypeError:
        maxn = data.params.max_low_qual_bases
    minlen = str(max(32, data.params.filter_min_trim_len))

    # vsearch merge can now take gzipped files (v.2.8)
    cmd = [
        BIN_VSEARCH,
        "--fastq_mergepairs", str(infile1),
        "--reverse", str(infile2),
        "--fastaout", str(mergedfile),
        "--fastaout_notmerged_fwd", str(nonmerged1),
        "--fastaout_notmerged_rev", str(nonmerged2),
        "--fasta_width", "0",
        "--fastq_minmergelen", str(minlen),
        "--fastq_maxns", str(maxn),
        "--fastq_minovlen", "10",  # default value,
        "--fastq_maxdiffs", "4",  # <- fastp has done pe overlap correction
        # "--label_suffix", "_m1",
        "--fastq_qmax", "93",     # <- Set high to allow FASTQ+64
        "--threads", "2",
        # "--fastq_allowmergestagger",
    ]
    print(f"@@DEBUG: {' '.join(cmd)}", flush=True)
    with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
        res = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(f"Error merge pairs:\n {cmd}\n{res}")

    # count number of lines / 2 from output
    with open(mergedfile, 'r', encoding="utf-8") as merged:
        nlines = sum(1 for line in merged)
    return int(nlines / 2)


def join_end_to_end(data: Assembly, sample: Sample) -> None:
    """Combine read pairs end to end for dereplication and/or clustering.
    """
    # output path (may already contain overlap merged pairs)
    joined = data.tmpdir / f"{sample.name}_merged.fa"

    # input paths (read pairs not overlap merged)
    nonmerged1 = data.tmpdir / f"{sample.name}_nonmerged_R1.fa"
    nonmerged2 = data.tmpdir / f"{sample.name}_nonmerged_R2.fa"

    # iterate over pairs and write append to out
    tmp = []
    with open(joined, 'a', encoding="utf-8") as out:
        # open both indata io
        with open(
            nonmerged1, 'r', encoding="utf-8") as in1, open(
                nonmerged2, 'r', encoding="utf-8") as in2:

            # generator to sample 2 lines at a time.
            quart1 = zip(*[in1] * 2)
            quart2 = zip(*[in2] * 2)

            # write out and fasta: >header\n{seq1}nnnn{seq2}
            for idx, (read1s, read2s) in enumerate(zip(quart1, quart2)):
                header = read1s[0][1:]
                read1 = read1s[1].strip()
                read2 = comp(read2s[1].strip())[::-1]
                newread = f">{header}{read1}nnnn{read2}"
                tmp.append(newread)

                # write on intervals
                if not idx % 50_000:
                    out.write("\n".join(tmp) + "\n")
                    tmp = []

        # write remaining to outfile
        if tmp:
            out.write("\n".join(tmp))

        # remove temporary files
        nonmerged1.unlink()
        nonmerged2.unlink()


def cluster(data: Assembly, sample: Sample) -> Tuple[Path, Path]:
    """Calls vsearch for clustering with cluster_smallmem.

    cov varies by data type, values were chosen
    based on experience, but could be edited by users.
    """
    # get dereplicated reads for denovo+reference or denovo-reference
    handles = [
        data.tmpdir / f"{sample.name}_derep.fa",
        data.tmpdir / f"{sample.name}_derep_tag.fa",
    ]
    derephandle = [i for i in handles if i.exists()][-1]
    assert derephandle, "bad derep handle"

    # create handles for the outfiles
    uhandle = data.stepdir / f"{sample.name}_matches.tsv"

    # minsl: the percentage of the seed that must be matched
    #    smaller values for RAD/ddRAD where we might want to combine, say 50bp
    #    reads and 100bp reads in the same analysis.
    minsl = 0.5
    # query_cov: the percentage of the query sequence that must match seed
    #    smaller values are needed for gbs where only the tips might overlap
    #    larger values for pairgbs where they should overlap near completely
    #    small minsl and high query cov allows trimmed reads to match to untrim
    #    seed for rad/ddrad/pairddrad.
    cov = 0.75

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

    # get call string
    # we used to do --cluster_smallmem w/ --usersort to sort by depth
    # but now there is an option --cluster_size to accomplish this,
    # however, --cluster_fast sorts by length, and this MAY be better
    # when most sequences are similar depths (e.g., 2-3).
    cmd = [
        BIN_VSEARCH,
        "--cluster_size", str(derephandle),
        "--strand", strand,
        "--query_cov", str(cov),
        "--id", str(data.params.clust_threshold),
        "--minsl", str(minsl),
        "--userout", str(uhandle),
        "--userfields", "query+target+id+gaps+qstrand+qcov",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        "--threads", str(max(1, data.ipcluster['threads'])),
        # "--notmatched", str(temphandle),
        "--fasta_width", "0",
        "--minseqlength", str(data.params.filter_min_trim_len),
        # "-fastq_qmax", "100",
        # "-fulldp",  # no longer relevant, it is now always on.
    ]

    # log cmd and run vsearch
    print(f"@@DEBUG: {' '.join(cmd)}", flush=True)

    with Popen(cmd, stderr=STDOUT, stdout=PIPE, close_fds=True) as proc:
        res = proc.communicate()[0].decode()
    if proc.returncode:
        print(f"@@ERROR: {' '.join(cmd)}\n{res}")
        raise IPyradError(f"cmd {cmd}: {res}")
    # log finished.
    print(f"@@INFO: {sample.name} finished clustering", flush=True)


if __name__ == "__main__":
    pass
