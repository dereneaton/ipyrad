#!/usr/bin/env python

"""Load fastqs per sample and trim with fastp.

"""

from typing import Tuple
import sys
import gzip
import json
from collections import Counter
from pathlib import Path
from subprocess import Popen, STDOUT, PIPE
from loguru import logger

logger = logger.bind(name="ipyrad")
FASTP_BINARY = Path(sys.prefix) / "bin" / "fastp"


class TrimFastqs:
    """Simple read trimming with fastp.

    https://github.com/OpenGene/fastp
    """
    def __init__(
        self,
        sname: str,
        reads: Tuple[Path, Path],
        outdir: Path,
        is_pair: bool,
        restriction_overhang: Tuple[str, str],
        trim_reads: Tuple[int, int],
        filter_min_trim_len: int,
        phred_qscore_offset: int,
        max_low_qual_bases: int,
        filter_adapters: int,
    ):

        # input args are passed as full paths
        self.sname = sname
        self.read1, self.read2 = reads
        self.is_pair = is_pair
        self.restriction_overhang = restriction_overhang
        self.trim_reads = trim_reads
        self.filter_min_trim_len = filter_min_trim_len
        self.phred_qscore_offset = phred_qscore_offset
        self.max_low_qual_bases = max_low_qual_bases
        self.filter_adapters = filter_adapters

        # output file paths (do not bother gzipping since these are tmp files)
        self.out1 = outdir / f"{sname}.trimmed_R1.fastq.gz"
        self.out2 = None
        if self.is_pair:
            self.out2 = outdir / f"{sname}.trimmed_R2.fastq.gz"

        # paths to stats files
        self.json = outdir / f'{sname}.fastp.json'
        self.html = outdir / f'{sname}.fastp.html'

        # get the command
        self.command = []
        self._check_trim_edges()
        self._build_command()

    def _check_trim_edges(self) -> None:
        # by default we trim the adapter lengths from the start of R1
        # and R2 reads based on auto-detection of the end position of
        # the restriction overhang. This is most reliable since users
        # may have demux'd their own data which may or may not include
        # the barcode still.
        trim_front1 = estimate_trim_position(
            self.read1, self.restriction_overhang[0], 200_000, end=True,
        )
        print(
            f"@@DEBUG: {self.sname} inferred length of barcode to trim "
            f"from front of read1: {trim_front1}", flush=True)

        if self.is_pair:
            trim_front2 = estimate_trim_position(
                self.read2, self.restriction_overhang[1], 200_000, end=True,
            )
            print(
                f"@@DEBUG: {self.sname} inferred length of barcode to trim "
                f"from front of read2: {trim_front1}", flush=True)
        else:
            trim_front2 = 0

        # HOwever, the auto trim can be overriden by setting values for
        # the `trim_reads` param to either a positive integer length,
        # or, if set to a negative value then NO TRIMMING will be done.
        # *If it was set to 0 then we use the estimated values above.*
        if self.trim_reads[0]:
            estimated = trim_front1
            if self.trim_reads[0] < 0:
                trim_front1 = 0
            else:
                trim_front1 = abs(self.trim_reads[0])
            print(
                "@@WARNING:"
                f"overriding estimated trim position ({estimated}) and "
                f"instead using user value for .trim_reads[0] = {trim_front1}",
                flush=True,
            )

        if self.trim_reads[1]:
            estimated = trim_front2
            if self.trim_reads[1] < 0:
                trim_front2 = 0
            else:
                trim_front2 = abs(self.trim_reads[1])
            print(
                "@@WARNING:"
                f"overriding estimated trim position ({estimated}) and "
                f"instead using user value for .trim_reads[1] = {trim_front2}",
                flush=True,
            )
        # store values
        self.trim_reads = (trim_front1, trim_front2)

    def _build_command(self):
        """Calls fastp in subprocess and writes tmpfiles to workdir.

        0 = no trim or quality check.
        1 = only trimming
        2 = trimming and quality checks.

        (-A turns off adapter trimming)
        (-Q turns off quality filtering)
        """
        if self.is_pair:
            cmd = [
                str(FASTP_BINARY),
                "-i", str(self.read1),
                "-I", str(self.read2),
                "-o", str(self.out1),
                "-O", str(self.out2),
                "--detect_adapter_for_pe",
                "-w", "2",         # 2 workers + 2 i/o threads.
            ]
        else:
            cmd = [
                str(FASTP_BINARY),
                "-i", str(self.read1),
                "-o", str(self.out1),
                "-w", "3",         # 3 workers + 1 i/o threads.
            ]

        # submit 4-threaded jobs 
        # the preferred fastp threading is 3 + io1 + io2 if PE

        cmd.extend([
            "-5",  # sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-3",  # sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-q", str(15 + self.phred_qscore_offset - 33),  # minqual
            "-l", str(self.filter_min_trim_len),  # minlen
            "-y", "-Y", "50",  # turns on and sets complexity filter to 50
            "-c",              # paired-end base correction
            "--n_base_limit", str(self.max_low_qual_bases),
            "-j", str(self.json),
            "-h", str(self.html),
            "--trim_front1", str(self.trim_reads[0]),
            "--trim_front2", str(self.trim_reads[1]),
        ])

        # TODO: UMIs here...
        # -U, --umi       enable unique molecular identifier (UMI) preprocessing
        # --umi_loc       specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
        # --umi_len       if the UMI is in read1/read2, its length should be provided (int [=0])
        # --umi_prefix    if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])        

        # turn off filters if settings are lower than 2
        # hard coded fasta adapters file with Truseq and AAAAAAA
        # in theory, not using -a turns on auto-detection of adapters, but
        # its not working well currently, so we supply a file for PE data,
        # but THAT's not working for SE data, so we supply a string.
        if self.filter_adapters == 2:
            # from fastp: "For PE data, the adapters can be detected by per-read
            # overlap analysis, which seeks for the overlap of each pair of reads.
            # This method is robust and fast, so normally you don't have to input
            # the adapter sequence even you know it."
            #
            # so, we only worry about SE data, especially GBS or ddRAD where
            # a short fragment read may be: ---------------------------->
            #                     [adapter1][bar][cut][data][rc-cut][rc-bar][adapter2]
            # and so we want to trim not only the adapter but also the
            # revcomp of the end cut and bar if present. SE ddrad will likely
            # not have a bar at other end, but GBS will... The best fix for
            # this is to not use -a so that fastp auto-detect the adapter as
            # the sequence [rc-cut][rc-bar][adapter2].

            # Thus, we do not enter -a for any setting.
            pass
            # if not self.is_pair:
            #     illumina = "AGATCGGAAGAGC"
            #     cmd.extend(["-a", illumina])

        # turns off by user request
        if self.filter_adapters == 1:
            cmd.extend("-A")
        if self.filter_adapters == 0:
            cmd.extend("-Q")
        self.command = cmd

    def run(self):
        """Run the fastp command on a subprocess"""
        with Popen(self.command, stderr=STDOUT, stdout=PIPE) as proc:
            out = proc.communicate()
        if proc.returncode:
            out = out[0].decode()
            print(
                "@@ERROR: "
                f"\ncmd: {' '.join(self.command)}"
                f"\nerr: {out}",
                flush=True
            )
            raise ValueError(out)

        # callback sent to logger.info on completion
        print(f"@@INFO: cmd: {' '.join(self.command)}", flush=True)

    def parse_stats_from_json(self):
        """Get stats from the fastp JSON file."""
        with open(self.json, 'r', encoding="utf-8") as indata:
            jdata = json.loads(indata.read())
        return (jdata, [(str(self.out1), str(self.out2))])


def estimate_trim_position(
    fastq: Path,
    restriction_overhang: str,
    max_len: int = 25,
    max_reads: int = 200_000,
    end: bool = True,
) -> int:
    """Return estimation position of end of cut1 on R1 files.

    At this point read1 files should have one of two formats:
    >>> ACACACTGCAGXXXX....
    >>> bbbbbb^^^^^dddd
    or
    >>> TGCAGXXXX....
    >>> ^^^^^dddd
    where
    b = barcode
    ^ = cut overhang
    d = data

    >>>           ------------------x---------------------->
    >>> [adapter1][barcode1][cutter1][data ..................]

    This function will read the first max_reads reads to find the avg
    position of the last cutter overhang position (i.e., end of barcoe).
    If the cutter overhang is not in >20% of reads then 0 is returned.

    Parameters
    ----------
    end: bool
        If True then the end of the position of the end of the overhang
        is returned, else the start position is returned.

    Note
    ----
    We use rfind from 25 bases in to accommodate cases in which the re
    sequence occurs within the barcode.
    """
    if str(fastq).endswith(".gz"):
        fp = gzip.open(fastq, 'rt', encoding="utf-8")
    else:
        fp = open(fastq, 'rt', encoding="utf-8")
    quart = zip(fp, fp, fp, fp)
    count = range(max_reads)

    # barcodes = Counter()
    positions = Counter()
    for idx, line in zip(count, quart):
        try:
            cut_start = line[1][:max_len].index(restriction_overhang)
            positions[cut_start] += 1
            # barcodes[line[1][:cut_start]] += 1
        except ValueError:
            pass

    # must occur in >20%
    commons = positions.most_common(10)
    if len(commons) > 1:
        best, second = commons[:2]
        if float(best[1] / second[1]) < 2.:
            print(
                "WARNING: a single best cut position could not be estimated "
                f"for {fastq}.", flush=True)
            return 0

    # trim up to Nth position, fastp takes trim N bases...
    # , barcodes.most_common(1)[0][0]
    if end:
        return positions.most_common(1)[0][0] + len(restriction_overhang)
    return positions.most_common(1)[0][0]


if __name__ == "__main__":

    # mixed sample, should not find one.
    print(estimate_trim_position("../../pedtest/small_tmp_R1.fastq.gz", "ATCGG", max_len=35, max_reads=500000, end=True))

    # demux sample w/ barcodes removed, should find 0.
    print(estimate_trim_position("../../pedtest/NEW_fastqs/lachnoglossa-DE123_R1.fastq.gz", "ATCGG", end=True))

    # print(estimate_trim_position("../../pedtest/NEW_fastqs/lachnoglossa-DE123_R2.fastq.gz", "ATCGG"))