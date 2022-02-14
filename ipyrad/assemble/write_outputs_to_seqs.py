#!/usr/bin/env python

"""Loads processed .loci chunks and writes to HDF5 seqs database.

Datasets in seqs.hdf5
---------------------
- phy
- phymap
- scaffold_names
- scaffold_lengths

Format of seqs.hdf5 'phymap'
----------------------------
- 1-indexed index of scaffold/locus names (in 'scaffold_names' str dataset)
- 0-indexed 'phy' start position
- 0-indexed 'phy' end position
- 1-indexed reference start position (just zeros if denovo)
- 1-indexed reference end position (just zeros if denovo)

Example phymap
--------------
>>> [1,   0, 300, 20000, 20300],     # early scaffolds are usually larger
>>> [1, 300, 600, 40000, 40300],     # and will have high ref start-end pos.
>>> [1, 600, 900, 55000, 55300],     # but low indices in the phy array.
>>> ...
>>> [512, 50000, 50300, 3000, 3300], # late scaffolds are usually smaller,
>>> [512, 50300, 50600, 5000, 5300]  # and so have low ref start-end, but
>>>                                  # high indices in the phy array.

Note
----
As much as I don't like it, the scaffolds in phymap are 1-indexed
for compatibility reasons with other formats like VCF, which do not
allow a 0-indexed scaffold or positions.
"""

from typing import Iterator, Dict, List, Tuple, TypeVar
from pathlib import Path
import h5py
import numpy as np
from ipyrad.assemble.utils import get_fai_values

Assembly = TypeVar("Assembly")

class SeqsDatabase:
    def __init__(self, data: Assembly):
        self.data = data
        """: Assembly object."""
        self.name = Path(data.stepdir) / f"{data.name}.seqs.hdf5"
        """: Output HDF5 file name."""

        # attributes to be filled.
        self.snames: List[str] = []
        """: A list of names in alphanumeric order, optionally including ref as sample."""
        self.sidxs: Dict[str, str] = {}
        """: A dict mapping names to the seqs database row index."""
        self.drop_ref: bool = (
            self.data.is_ref & self.data.hackers.exclude_reference)
        """: Boolean for whether reference samples exists and should be dropped."""
        self.nmissing: int = 0
        """: Counter to record number of cells of missing data (N or -)."""

    def _get_sorted_loci_chunks(self) -> None:
        """Gather all loci bits and sort by filename."""
        return sorted(
            Path(self.data.tmpdir).glob("*.loci"),
            key=lambda x: int(str(x).rsplit("-", 1)[-1][:-5]),
        )

    def _init_datasets(self) -> None:
        """Create the datasets in the seqs database.

        The order of samples in the seqs database is alphanumeric
        except for the 'reference' sample, if included as a sample
        using the hackers option, in which case it is listed first.
        """
        with h5py.File(self.name, 'w') as io5:

            # create the datasets.
            io5.create_dataset(
                name="phy",
                shape=(self.data.assembly_stats.nsamples, self.data.assembly_stats.nbases),
                dtype=np.uint8,
            )
            phymap = io5.create_dataset(
                name="phymap",
                shape=(self.data.assembly_stats.nloci, 5),
                dtype=np.uint64,
            )

            # alphanumeric name order except 'reference' on top.
            self.snames = sorted(self.data.samples)

            # drop the reference sample, and add it back as first sample
            # unless we're really dropping it for real. As first sample it
            # is easy to use for filtering later.
            self.snames.remove("reference")
            if not self.drop_ref:
                self.snames = ["reference"] + self.snames
            self.sidxs = {sample: i for (i, sample) in enumerate(self.snames)}

            # store the scaffold names and lengths. The name strings are
            # now always stored as strings in HDF5 attrs.
            if self.data.is_ref:
                io5["scaffold_lengths"] = get_fai_values(self.data, "length")
                io5["scaffold_names"] = get_fai_values(self.data, "scaffold").astype("S")

            # store the name of reference genome file.
            if self.data.is_ref:
                phymap.attrs["reference"] = Path(self.data.params.reference_sequence).name
            else:
                phymap.attrs["reference"] = "pseudoref"

            # store ordered names and column labels
            phymap.attrs["names"] = self.snames
            phymap.attrs["columns"] = ["chroms", "phy0", "phy1", "pos0", "pos1"]

    def _iter_loci(self) -> Iterator[Tuple[Dict[str,str], Tuple[int,str,int,int]]]:
        """Yields loci from each ordered .loci file until empty.
        """
        for locfile in self._get_sorted_loci_chunks():
            with open(locfile, 'r', encoding="utf-8") as indata:

                names_to_seqs = {}
                for line in indata:

                    # within a locus, fill the data.
                    if line[0] != "/":
                        name, seq = line.split()
                        if (name == "reference") & self.drop_ref:
                            continue
                        names_to_seqs[name] = bytes(seq, "utf-8")

                    else:
                        # parse reference position from snpstring
                        chrom_int, chrom_name, pos0, pos1 = 0, "", 0, 0
                        if self.data.is_ref:
                            line_chunks = line.split("|")
                            chrom_int, chrom_name, pos = line_chunks[1].split(":")
                            chrom_int = int(chrom_int)
                            pos0, pos1 = [int(i) for i in pos.split("-")]

                        # end of locus, yield the dict.
                        yield names_to_seqs, (chrom_int, chrom_name, pos0, pos1)

    def _iter_supermatrix_chunks(self) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
        """Yields chunks of supermatrix and phymap.

        yields (phy, phymap, nmissing) for chunks of ~50_000 columns
        of the supermatrix at a time.
        """
        # size of chunks that will be yielded, with some overflow of columns
        chunk_size = 50_000
        nsamps = len(self.snames)

        # make array to fill w/ overflow on column size.
        tmparr = np.zeros((nsamps, chunk_size + 50_000), dtype=np.uint8)
        tmpmap = []

        # iterate over loci from .loci files
        pos = 0
        for names_to_seqs, positions in self._iter_loci():

            # get locus dict and ref positions
            chrom_int, _, pos0, pos1 = positions

            # order the dict by alphanumeric names except with reference first.
            ordered_seqs = (names_to_seqs[i] for i in self.snames if i in names_to_seqs)

            # convert to a uint8 array so we can filter inserts and fill Ns
            loc = np.array([list(i) for i in ordered_seqs], dtype=np.uint8)

            # filter PE inserts (all N in is_ref, or all N or -)
            if self.data.is_ref:
                # if ref was not excluded already then we need to pass over
                # it when looking for inserts since it's inserts are not empty.
                if not self.drop_ref:
                    mask = np.all(loc[1:] == 78, axis=0)
                else:
                    mask = np.all(loc == 78, axis=0)
            else:
                mask = np.all((loc == 45) | (loc == 78), axis=0)
            loc = loc[:, np.invert(mask)]

            # enter data into the tmparr
            for tmp_sidx, name in enumerate(names_to_seqs):
                global_sidx = self.sidxs[name]
                tmparr[global_sidx, pos:pos + loc.shape[1]] = loc[tmp_sidx]

            # enter position data into the tmpmap
            tmpmap.append((chrom_int, pos, pos + loc.shape[1], pos0, pos1))

            # if the tmparr is full, yield it and refresh.
            pos += loc.shape[1]
            if pos >= chunk_size:

                # trim extra, fill empty cells to N (78), and convert map to array.
                tmparr = tmparr[:, :pos]
                tmparr[tmparr == 0] = 78
                tmpmap = np.array(tmpmap, dtype=np.uint64)
                self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
                yield tmparr, tmpmap

                # reset
                tmparr = np.zeros((nsamps, chunk_size + 50_000), dtype=np.uint8)
                tmpmap = []
                pos = 0

        # yield the final leftover chunk that did not reach chunk_size.
        if pos:
            tmparr = tmparr[:, :pos]
            tmparr[tmparr == 0] = 78
            tmpmap = np.array(tmpmap, dtype=np.uint64)
            self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
            yield tmparr, tmpmap

    def _fill_database(self) -> None:
        """Iterates over chunks of data entering into the H5 database.
        """
        start_arr = 0
        start_map = 0
        with h5py.File(self.name, 'a') as io5:
            for arr_chunk, map_chunk in self._iter_supermatrix_chunks():

                # get end positions
                end_arr = start_arr + arr_chunk.shape[1]
                end_map = start_map + map_chunk.shape[0]

                # insert to dataset array
                io5["phy"][:, start_arr:end_arr] = arr_chunk
                io5["phymap"][start_map:end_map] = map_chunk

                # increment start positions
                start_arr += arr_chunk.shape[1]
                start_map += map_chunk.shape[0]

    def _write_stats(self):
        """Write missing values in matrix to stats file, and logger.

        Because this is run on an engine the logger will report what
        is written with print.
        """
        with h5py.File(self.name, 'r') as io5:
            shape = io5["phy"].shape               # pylint: disable=no-member
            size = io5["phy"].size                 # pylint: disable=no-member            
            perc_missing = self.nmissing / size
        print(f"wrote seqs.hdf5, shape={shape}, "
            f"missing={perc_missing:.2f}, path={self.name}.")

        # write stats to the output file
        stats_file = Path(self.data.stepdir) / "s7_assembly_stats.txt"
        with open(stats_file, 'a', encoding="utf-8") as out:
            outstr = (
                f"sequence matrix size: {shape}, "
                f"{perc_missing:.2f}% missing sites.\n"
            )
            out.write(outstr)

    def run(self):
        """Initialize HDF5 database and then fill it in chunks."""
        self._init_datasets()
        self._fill_database()
        self._write_stats()


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    DATA = ip.load_json("/tmp/TEST5.json")
    # uncomment to include the ref
    # DATA.hackers.exclude_reference = False
    print("EXCLUDE REF=", DATA.hackers.exclude_reference)

    with ip.Cluster(4) as ipyclient:
        step = Step7(DATA, ipyclient=ipyclient, force=True, quiet=False)
        step._split_clusters()
        step._apply_filters_and_trimming()
        step._collect_stats()

        # run it.
        db = SeqsDatabase(DATA)
        db.run()

        # print supermatrix parts from start and end of loci file.
        with h5py.File(db.name, 'r') as io5:
            # print first ten sites in first chunk file
            for i, n in enumerate(db.snames):
                print(n, io5["phy"][i, :16].tobytes().decode())
            print()
            for i, n in enumerate(db.snames):
                print(n, io5["phy"][i, -16:].tobytes().decode())
