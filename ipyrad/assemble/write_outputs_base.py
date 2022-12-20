#!/usr/bin/env python

"""Loads processed .loci chunks and writes to HDF5 snps database.

"""

from typing import Iterator, Dict, List, Tuple, TypeVar
from abc import ABC, abstractmethod
from pathlib import Path

Assembly = TypeVar("Assembly")

class DatabaseWriter(ABC):
    def __init__(self, data: Assembly, samples: Dict[str,"SampleSchema"]):
        self.data: Assembly = data
        """: Assembly object."""
        self.samples: Dict[str,"SampleSchema"] = samples
        """: Dict of SampleSchema objects in Step7."""

        # attributes to be filled.
        self.loci_chunks: List[Path] = None
        """: A sorted list of Paths to loci chunks."""
        self.snames: List[str] = None
        """: A list of names in alphanumeric order, optionally including ref as sample."""
        self.sidxs: Dict[str, str] = {}
        """: A dict mapping names to the seqs database row index."""
        self.nmissing: int = 0
        """: Counter to record number of cells of missing data (N or -)."""
        self.snppad: int = 0
        """: Length of name region in .loci file before sequences start."""

    def _get_sorted_loci_chunks(self) -> None:
        """Gather all loci bits and sort by filename. And get snppad len."""
        self.loci_chunks = sorted(
            Path(self.data.tmpdir).glob("chunk-*.loci"),
            key=lambda x: int(x.name.split('-')[1]))

    def _get_snppad(self) -> None:
        """Set .snppad with length of spacer before SNP info starts."""
        with open(self.loci_chunks[0], 'r', encoding="utf-8") as indata:
            self.snppad = indata.readline().rfind(" ") + 1

    def _get_snames(self) -> None:
        """Get sorted names with 'reference' on top if present."""
        self.snames = sorted(self.samples)
        if self.data.drop_ref:
            self.snames.remove("reference")
        elif self.data.is_ref:
            self.snames.remove("reference")
            self.snames = ['reference'] + self.snames
        self.sidxs = {i: self.snames.index(i) for i in self.snames}

    @abstractmethod
    def _init_datasets(self) -> None:
        """Create a database output file."""

    def _iter_loci(self) -> Iterator[Tuple[Dict[str,str], Tuple[int,str,int,int]]]:
        """Yields loci from each ordered .loci file until empty.

        The 'reference' sample is always included here when is_ref.
        """
        lidx = 0
        for locfile in self.loci_chunks:
            with open(locfile, 'r', encoding="utf-8") as indata:

                # iterate over loci and yield one at a time.
                names_to_seqs = {}
                for line in indata:

                    # within a locus, fill the data.
                    if line[0] != "/":
                        name, seq = line.split()
                        names_to_seqs[name] = bytes(seq, "utf-8")

                    else:
                        # store the snpstring too
                        names_to_seqs['snpstring'] = line[self.snppad:].split("|")[0]

                        # parse reference position from snpstring
                        if self.data.is_ref:
                            line_chunks = line.split("|")
                            chrom_int, chrom_name, pos = line_chunks[1].split(":")
                            chrom_int = int(chrom_int)
                            pos0, pos1 = [int(i) for i in pos.split("-")]
                        else:
                            chrom_int, chrom_name, pos0, pos1 = lidx, "RAD", 0, 0

                        # end of locus, yield the dict.
                        yield names_to_seqs, (chrom_int, chrom_name, pos0, pos1)
                        lidx += 1
                        names_to_seqs = {}                        

    @abstractmethod
    def run(self):
        """Initialize HDF5 database and then fill it in chunks."""


if __name__ == "__main__":
    pass
