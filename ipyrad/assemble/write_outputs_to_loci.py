#!/usr/bin/env python

"""Concatenate .loci chunks and format meta-data in snpstring.

The data in the .loci chunks is already sorted alphanumerically, but
for reference-based assemblies it contains the 'reference' sequence
as the first/top sample. This will be removed if the exclude_reference
is True. If so, the SNP calling and filtering has already taken into
account the expectation that it will be removed. It has only been 
kept thus far so that we can store the REF allele in the HDF5 files.

Denovo
------
write a final (filtered) locus number: |2|

Reference
---------
Clean the snpstring to be: |{chromint}:{chrom}:{start}-{end}|
"""

from typing import Iterator, List
from ipyrad.assemble.write_outputs_base import DatabaseWriter

class LociWriter(DatabaseWriter):
    """Join locus bits into a single file."""

    def _iter_chunk_of_loci(self) -> Iterator[List[str]]:
        """Yields loci from each ordered .loci file until empty.

        Drops the reference sequence if .drop_ref.
        """
        lidx = 0
        for locfile in self.loci_chunks:
            loci = []
            with open(locfile, 'r', encoding="utf-8") as indata:
                loc = []
                for line in indata:
                    if line[0] != "/":
                        # EXCLUDE REF FROM LOCI? OR LEAVE IT?
                        if line.startswith("reference ") & self.data.drop_ref:
                            continue
                        loc.append(line.strip())
                    else:
                        # store the snpstring and meta-info line.
                        if self.data.is_ref:
                            print("DO REF loCI string here..")
                        else:
                            line = line.strip().rsplit("|", 3)[0]
                            loc.append(f"line|{lidx}|")
                        lidx += 1
                loci.append("\n".join(loc))
            # yield the locus for the entire chunk.
            yield loci

    def _init_datasets(self):
        """Not used here, but overwriting from base class."""

    def run(self):
        """Iterate over loci, clean snpstring, end write to final loci file."""
        # fill .loci_chunks with sorted tmp files
        self._get_sorted_loci_chunks()
        # format and write to file.
        with open(self.data.outfiles["loci"], 'w', encoding="utf-8") as out:
            for loci_chunk in self._iter_chunk_of_loci():
                out.write("\n".join(loci_chunk) + "\n")


# def get_nondash_from_left(seq: str) -> int:
#     """Get index of first non dash character from left."""
#     for i, j in enumerate(seq):
#         if j != "-":
#             return i
#     return 0

# def get_nondash_from_right(seq: str) -> int:
#     """Get index of first non dash character from right."""
#     lseq = len(seq)
#     for i, j in enumerate(seq[::-1]):
#         if j != "-":
#             return lseq - i
#     return len(seq)


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7

    # data = ip.load_json("/tmp/TEST5.json")
    data = ip.load_json("/home/deren/Documents/ipyrad/sra-fastqs/cyatho.json")    
    data.hackers.exclude_reference = True

    with ip.Cluster(cores=4) as ipyclient:
        step = Step7(data, quiet=False, force=True, ipyclient=ipyclient)
        step.run()
        LociWriter(step.data, step.samples).run()
