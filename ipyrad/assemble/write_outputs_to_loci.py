#!/usr/bin/env python

"""Concatenate .loci chunk to .loci file.

Cleans the snpstring to be: {chromint}:{chrom}:{start}-{end}
"""

from typing import Iterator, Dict, List
from pathlib import Path

class LociWriter:
    def __init__(self, data):
        self.data = data

        # attributes to be filled.
        self.snames: List[str] = []
        """: A list of names in alphanumeric order, optionally including ref as sample."""
        self.sidxs: Dict[str, str] = {}
        """: A dict mapping names to the seqs database row index."""
        self.drop_ref: bool = (
            self.data.is_ref & self.data.hackers.exclude_reference)
        """: Boolean for whether reference samples exists and should be dropped."""
        self._set_sorted_names()

    def _set_sorted_names(self):
        """Get sorted names.

        Names are in alphanumeric order except 'reference' which is
        put first, if present, unless `hackers.exclude_reference=True`.
        """
        self.snames = sorted(self.data.samples)

        # drop the reference sample, and add it back as first sample
        # unless we're really dropping it for real. As first sample it
        # is easy to use for filtering later.
        self.snames.remove("reference")
        if not self.drop_ref:
            self.snames = ["reference"] + self.snames
        self.sidxs = {sample: i for (i, sample) in enumerate(self.snames)}        

    def _get_sorted_loci_chunks(self) -> None:
        """Gather all loci bits and sort by filename."""
        return sorted(
            Path(self.data.tmpdir).glob("*.loci"),
            key=lambda x: int(str(x).rsplit("-", 1)[-1][:-5]),
        )

    def _iter_loci(self) -> Iterator[List[str]]:
        """Yields loci from each ordered .loci file until empty.

        Drops the reference sequence if .drop_ref.
        """
        for locfile in self._get_sorted_loci_chunks():
            loci = []
            with open(locfile, 'r', encoding="utf-8") as indata:
                loc = []
                for line in indata:
                    if line[0] != "/":
                        # EXCLUDE REF FROM LOCI? OR LEAVE IT?
                        # if line.startswith("reference") & self.drop_ref:
                            # continue
                        loc.append(line.strip())
                    else:
                        line = line.strip().rsplit("|", 2)[0]
                        loc.append(line + "|")
                loci.append("\n".join(loc))
            # yield the locus for the entire chunk.
            yield loci

    def run(self):
        """Iterate over loci, clean snpstring, end write to final loci file."""
        with open(self.data.outfiles["loci"], 'w', encoding="utf-8") as out:
            for loci_chunk in self._iter_loci():
                out.write("\n".join(loci_chunk) + "\n")

if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7
    with ip.Cluster(cores=4) as ipyclient:
        data = ip.load_json("/tmp/TEST5.json")
        step = Step7(data, quiet=False, force=True, ipyclient=ipyclient)
        step._split_clusters()
        step._apply_filters_and_trimming()
        LociWriter(step.data).run()
