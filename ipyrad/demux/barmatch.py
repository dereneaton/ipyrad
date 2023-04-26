#!/usr/bin/env python

"""Classes for matching barcodes for different datatypes.

"""

from typing import Dict, Tuple, List, TypeVar, Iterator
import io
import itertools
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import gzip
from dataclasses import dataclass, field
from loguru import logger
from ipyrad.core.exceptions import IPyradError

logger = logger.bind(name="ipyrad")
Assembly = TypeVar("Assembly")
# CHUNKSIZE = 10_000_000


@dataclass
class BarMatching:
    """Base class for barcode matching.

    See subclasses which have different versions of the function
    `_iter_matched_barcode` to find barcode matches based on i7,
    combinatorial, or single inline barcodes. The subclasses all
    share the functions of this class, which includes iterating
    over the fastq(s), storing stats, and writing to tmp files.
    """
    fastqs: Tuple[str, str]
    """: A tuple with paired R1 and R2 fastq files."""
    barcodes_to_names: Dict[str, str]
    """: Dict matching barcodes to sample names."""
    cuts1: List[str]
    """: List of RE overhangs to match on R1."""
    cuts2: List[str]
    """: List of RE overhangs to match on R2."""
    merge_technical_replicates: bool
    """: ..."""
    outpath: Path
    """: ..."""
    cores: int
    """: ..."""
    chunksize: int

    # stats counters
    barcode_misses: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that don't match."""
    barcode_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that match."""
    sample_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record number of hits per sample."""

    def __post_init__(self):
        self._format_check()

    def _format_check(self):
        """Check that data is appropriate for selected format."""
        pass

    def _iter_fastq_reads(self) -> Iterator[Tuple[List[str], List[str]]]:
        """Yields fastq quartets of lines from fastqs (gzip OK)."""
        # create first read iterator for paired data
        opener = gzip.open if self.fastqs[0].suffix == ".gz" else io.open
        # ofile1 = opener(self.fastqs[0], 'rt', encoding="utf-8")
        ofile1 = opener(self.fastqs[0], 'rb')
        quart1 = zip(ofile1, ofile1, ofile1, ofile1)

        # create second read iterator for paired data
        if self.fastqs[1]:
            # ofile2 = opener(self.fastqs[1], 'rt', encoding="utf-8")
            ofile2 = opener(self.fastqs[1], 'rb')
            quart2 = zip(ofile2, ofile2, ofile2, ofile2)
        else:
            quart2 = iter(int, 1)

        # yield from iterators as 4 items as a time (fastq)
        for read1, read2 in zip(quart1, quart2):
            yield read1, read2

    def _iter_matched_barcode(self):
        """SUBCLASSES REPLACE THIS FUNCTION."""
        raise NotImplementedError("See subclasses.")

    def _iter_matched_chunks(self) -> Iterator[Tuple[List[str], List[str]]]:
        """Stores matched reads until N then writes to file."""
        read1s = {}
        read2s = {}
        nstored = 0

        # iterate over matched reads
        for read1, read2, match in self._iter_matched_barcode():
            # store r1 as 4-line string
            fastq1 = b"".join(read1)
            if match in read1s:
                read1s[match].append(fastq1)
            else:
                read1s[match] = [fastq1]

            # store r2 as 4-line string
            if read2:
                fastq2 = b"".join(read2)
                if match in read2s:
                    read2s[match].append(fastq2)
                else:
                    read2s[match] = [fastq2]

            # write to file when size is big enough and reset.
            nstored += 1
            if nstored > self.chunksize:
                yield read1s, read2s
                read1s = {}
                read2s = {}
                nstored = 0

        # write final chunk if data
        yield read1s, read2s

    def run(self) -> None:
        """Multiprocessed writing is much faster, especially on HPC.

        Some overhead from i/o limitations, but most time here is spent
        on the string concatenation and gzip compression, which can
        happen in parallel on different engines.

        TODO
        ----
        Allow writers to use read1s or read2s dict as shared memory
        to prevent duplication of the Memory usage here on every CPU.
        Example here: https://stackoverflow.com/questions/65980183/processpoolexecutor-on-shared-dataset-and-multiple-arguments
        """
        with ProcessPoolExecutor(max_workers=self.cores) as pool:
            total = 0
            for read1s, read2s in self._iter_matched_chunks():
                nprocessed = min(self.chunksize, sum(len(i) for i in read1s.values()))
                total += nprocessed
                logger.info(
                    f"writing/compressing {nprocessed:.0f} matched reads "
                    f"(total={total:.0f})")

                rasyncs = {}

                # parallel workers cannot write to the same file so lists
                # assigned to technical replicates need to be grouped
                if self.merge_technical_replicates:
                    groups = itertools.groupby(
                        list(read1s), key=lambda x: x.rsplit("-technical-replicate-")[0]
                    )
                    for key, names in groups:
                        names = list(names)
                        if "-technical-replicate-" in names[0]:
                            read1s[key] = list(itertools.chain(*[read1s.pop(i) for i in sorted(names)]))
                            read2s[key] = list(itertools.chain(*[read2s.pop(i) for i in sorted(names)]))
                    logger.warning(list(read1s))

                # both dicts share the same names
                for name in read1s:
                    # # if merging tech reps then remove suffix
                    # if self.merge_technical_replicates:
                    #     fname = name.split("-technical-replicate-")[0]
                    # else:
                    #     fname = name

                    # write to R1 chunk file.
                    path1 = self.outpath / f"{name}_R1.fastq.gz"
                    data = read1s[name]
                    rasyncs[f"{name}_R1"] = pool.submit(write, *(path1, data))

                    # write to R2 chunk file.
                    if read2s:
                        path2 = self.outpath / f"{name}_R2.fastq.gz"
                        data = read2s[name]
                        rasyncs[f"{name}_R2"] = pool.submit(write, *(path2, data))

                # raise exception for any writing errors
                for res in rasyncs:
                    rasyncs[res].result()


@dataclass
class BarMatchingI7(BarMatching):
    """Subclass of Barmatching that matches barcode in i7 header.

    Example 3RAD R1 file with i7 tag in header
    ------------------------------------------
    >>> # asterisk part is the i7 --->                  ********
    >>> @NB551405:60:H7T2GAFXY:4:21612:8472:20380 1:N:0:TATCGGTC+ACCAGGGA
    >>> ATCGGTATGCTGGAGGTGGTGGTGGTGGAGGTGGACGTTACAAGGGTTCTGGTGGTAGCCGATCAG...
    >>> +
    >>> EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE...
    """
    # check that i7's exist for this data
    def _format_check(self):
        read1, _ = next(self._iter_fastq_reads())
        barcode = read1[0].strip().rsplit(b":", 1)[-1].split(b"+")[0]
        if not barcode:
            raise IPyradError(
                "No i7 index exists in this data. Example read1:\n{read1}")

    def _iter_matched_barcode(self) -> Iterator[Tuple[str, str, str]]:
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():
            # pull barcode from header
            barcode = read1[0].strip().rsplit(b":", 1)[-1].split(b"+")[0]
            # look for match
            match = self.barcodes_to_names.get(barcode)

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1
                yield read1, read2, match
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1


@dataclass
class BarMatchingSingleInline(BarMatching):
    """Subclass of Barmatching SE or PE data w/ inline barcodes only on R1.

    Example R1 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= restriction overhang.
    >>>
    >>> ********-----
    >>> @E00526:227:H53YNCCX2:8:1202:7710:23354 1:N:0:
    >>> CTGCAACTATCGGAGCGAATGAAAC........GACTCAACATAACGGGTCTGATCATTGAG
    >>> +
    >>> AA<FFJJJJJJJJJJJJJJJJJJJJ........JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    """
    maxlen1: int = 0
    """: Max len of the read1 inline barcodes."""

    def __post_init__(self):
        self.maxlen1 = max([len(i) for i in self.barcodes_to_names])
        self.maxlen1 += len(self.cuts1[0])

    def _iter_matched_barcode(self) -> Iterator[Tuple[str, str, str]]:
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():

            # find barcode from start of R1 (barcode1 + RE1 overhang)
            barcode = cut_matcher(read1[1][:self.maxlen1], self.cuts1)

            # look for matches
            match = self.barcodes_to_names.get(barcode)
            # print(f"Match={match}; barcode={barcode}; read1={read1[1][:20]}")

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1

                # TRIM inline barcode
                read1 = [read1[0], read1[1][len(barcode):], read1[2], read1[3][len(barcode):]]
                yield read1, read2, match
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1


@dataclass
class BarMatchingCombinatorialInline(BarMatching):
    """Subclass of Barmatching for combinatorial inline barcodes.

    Example R1 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= restriction overhang.
    >>>
    >>> ********-----
    >>> @E00526:227:H53YNCCX2:8:1202:7710:23354 1:N:0:
    >>> CTGCAACTATCGGAGCGAATGAAAC........GACTCAACATAACGGGTCTGATCATTGAG
    >>> +
    >>> AA<FFJJJJJJJJJJJJJJJJJJJJ........JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

    Example R2 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= restriction overhang.
    >>>
    >>> ********----
    >>> @E00526:227:H53YNCCX2:8:1202:7446:23354 2:N:0:CGAACTGT+ACAACAGT
    >>> ATGCTGTCGATCCCAACCACCACGC........TTTTTTTCTATCTCAACTATTTACAACAA
    >>> +
    >>> AAFFFJJJJJJJJJFJFJJJJJJ-F........AFJ<JFJJJJAJFFAA-F<A-AAF-AFFJ
    """
    maxlen1: int = 0
    """: Max len of the read1 inline barcodes + re."""
    maxlen2: int = 0
    """: Max len of the read2 inline barcodes + re."""

    def __post_init__(self):
        self.maxlen1 = max([len(i.split(b"_")[0]) for i in self.barcodes_to_names])
        self.maxlen1 += len(self.cuts1[0])
        self.maxlen2 = max([len(i.split(b"_")[1]) for i in self.barcodes_to_names])
        self.maxlen2 += len(self.cuts2[0])

    def _iter_matched_barcode(self):
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        # get a list of cutters and off-by-one's
        for read1, read2 in self._iter_fastq_reads():

            # find barcode from start of R1 (barcode1 + RE1 overhang)
            match_r1 = cut_matcher(read1[1][:self.maxlen1], self.cuts1)

            # pull barcode from start of R2 (barcode2 + RE2 overhang)
            match_r2 = cut_matcher(read2[1][:self.maxlen2], self.cuts2)

            # look for matches
            # barcode = f"{match_r1}_{match_r2}"
            barcode = match_r1 + b"_" + match_r2
            match = self.barcodes_to_names.get(barcode)

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1

                # TRIM INLINE BARCODE(S)
                read1 = [read1[0], read1[1][len(match_r1):], read1[2], read1[3][len(match_r1):]]
                read2 = [read2[0], read2[1][len(match_r2):], read2[2], read2[3][len(match_r2):]]
                yield read1, read2, match
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1


@dataclass
class BarMatch2BRADInline(BarMatching):
    """TODO: Need some test data for this, copied from older ipyrad code.

    # for 2brad we trim the barcode AND the synthetic overhang
    # The `+1` is because it trims the newline
    if self.data.params.datatype == '2brad':
        overlen = len(self.cutters[0][0]) + lenbar1 + 1
        read1[1] = read1[1][:-overlen] + "\n"
        read1[3] = read1[3][:-overlen] + "\n"
    else:
        read1[1] = read1[1][lenbar1:]
        read1[3] = read1[3][lenbar1:]
    """


def cut_matcher(read: str, cutters: List[str]) -> str:
    """Returns the barcode sequence before the detected re overhang.

    This will test each of the input cutters, which are multiple if it
    contains an ambiguity code, and then also tests cutters that are
    differ from the cutter sequence by `offby` (default=1). This
    generator stops when the first valid hit occurs, for speed.
    """
    for cut in cutters:
        pos = read.rfind(cut)
        if pos > 0:
            return read[:pos]
    return b"XXX"


def write(path: Path, data: List[str]) -> None:
    with gzip.open(path, 'a') as out:
        out.write(b"".join(data))


if __name__ == "__main__":
    pass
