#!/usr/bin/env python

"""Loads processed .loci chunks and writes to HDF5 seqs database.

Attrs in seqs_hdf5
------------------
- version: int : 1 is current. Older will not have this attr.
- names: List[str] : alphanumerically sorted sample names but with 
    'reference ' as the first sample if present.
- reference: str : filename of reference fasta file.

Datasets in seqs.hdf5
---------------------
- phy              : shape=(nsamples, nsites)     : arr of uint8 seqs
- phymap           : shape=(nsites, 5)            : arr mapping seqs windows to genomic positions or loci.
- scaffold_names   : shape=(nscaffolds or nloci,) : names of scaffolds or locus numbers
- scaffold_lengths : shape=(nscaffolds or nloci,) : lengths of scaffolds or loci

Format of seqs.hdf5 'phymap'
----------------------------
- 0-indexed index of scaffold (or locus idx if denovo)
- 0-indexed 'phy' start position 
- 0-indexed 'phy' end position
- 0-indexed reference start position (always zero if denovo)
- 0-indexed reference end position (always loclen if denovo)

Example phymap
--------------
>>> [0,   0, 300, 20000, 20300],     # early scaffolds are usually larger
>>> [0, 300, 600, 40000, 40300],     # and will have high ref start-end pos.
>>> [0, 600, 900, 55000, 55300],     # but low indices in the phy array.
>>> ...
>>> [512, 50000, 50300, 3000, 3300], # late scaffolds are usually smaller,
>>> [512, 50300, 50600, 5000, 5300]  # and so have low ref start-end, but
>>>                                  # high indices in the phy array.

# Scaff0, idx0, idx300, scaffpos20000, scaffpos20300       # RAD locus 0
# Scaff0, idx300, idx600, scaffpos40000, scaffpos40300     # RAD locus 1
# ...                                                      # RAD loci...
# Scaff512, idx50300, idx50600, scaffpos5000, scaffpos5300 # RAD locus -1

Note
----
The phy array includes extra columns at the end containing 0s. This is
by design, albeit a little wasteful. The data should always be accessed
from this file with respect to the phymap array which specifies which
windows correspond to which loci.
"""

from typing import Dict, TypeVar, Iterator, Tuple
from pathlib import Path
import h5py
import numpy as np
from ipyrad.assemble.utils import get_fai_values
from ipyrad.assemble.write_outputs_base import DatabaseLoader

Assembly = TypeVar("Assembly")
CHUNKSIZE = 50_000

class SeqsDatabase(DatabaseLoader):

    def __init__(self, data: Assembly, samples: Dict[str,"SampleSchema"]):
        self.name = Path(data.stepdir) / f"{data.name}.seqs_hdf5"
        """: Output HDF5 file name."""
        super().__init__(data, samples)

    def _init_datasets(self) -> None:
        """Create the datasets in the seqs database.

        The order of samples in the seqs database is alphanumeric
        except for the 'reference' sample, if included as a sample
        using the hackers option, in which case it is listed first.
        """
        with h5py.File(self.name, 'w') as io5:

            # store meta-data
            io5.attrs["names"] = self.snames
            io5.attrs["version"] = 1
            if self.data.is_ref:
                io5.attrs["reference"] = Path(self.data.params.reference_sequence).name
            else:
                io5.attrs["reference"] = "pseudoref"

            # create the seq and seq positions arr datasets.
            _ = io5.create_dataset(
                name="phy",
                shape=(self.data.assembly_stats.nsamples, self.data.assembly_stats.nbases),
                dtype=np.uint8,
            )
            phymap = io5.create_dataset(
                name="phymap",
                shape=(self.data.assembly_stats.nloci, 5),
                dtype=np.uint64,
            )

            # Store phymap meta-data
            # store the scaffold names and lengths. The name strings are
            # now always stored as strings in HDF5 attrs.
            if self.data.is_ref:
                phymap["scaffold_lengths"] = get_fai_values(self.data, "length")
                phymap["scaffold_names"] = get_fai_values(self.data, "scaffold").astype("S")

            # store ordered names and column labels
            phymap.attrs["indexing"] = [0, 0, 0, 0, 0]
            phymap.attrs["columns"] = ["chroms", "phy0", "phy1", "pos0", "pos1"]

    def _iter_supermatrix_chunks(self) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Yields chunks of supermatrix, snpmap, and reference seq.

        """
        # concat size of chunks that will be yielded, with some overflow of columns
        nsamps = self.data.assembly_stats.nsamples

        # make array to fill w/ overflow on column size.
        tmparr = np.zeros((nsamps, CHUNKSIZE + 5_000), dtype=np.uint8)
        # tmpref = np.zeros((CHUNKSIZE + 5_000), dtype=np.uint8)
        tmpmap = []

        # iterate over loci from .loci files incrementing the locus every
        # time so that the locus IDs in snps.hdf5 match those in seqs.hdf5.
        # keep track of position in array to send a chunk when full.
        pos = 0
        for names_to_seqs, positions in self._iter_loci():

            # get locus dict and ref positions
            chrom_int, _, pos0, pos1 = positions

            # remove SNP string
            _ = names_to_seqs.pop("snpstring")

            # order the dict by alphanumeric names except with reference first.
            # ordered_seqs = (names_to_seqs[i] for i in self.snames if i in names_to_seqs)
            ordered_seqs = [names_to_seqs[i] for i in self.snames if i in names_to_seqs]

            # convert locus to uint8 array and subselect only SNP sites.
            loc = np.array([list(i) for i in ordered_seqs], dtype=np.uint8)

            # create mask to filter PE insert (all N in is_ref, or all N or -)
            if self.data.is_ref:
                # if ref was not excluded already then we need to pass over
                # it when looking for inserts since it's inserts are not empty.
                if not self.data.drop_ref:
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
            if not self.data.is_ref:
                tmpmap.append((chrom_int, pos, pos + loc.shape[1], pos0, loc.shape[1]))
            else:
                tmpmap.append((chrom_int, pos, pos + loc.shape[1], pos0, pos1))

            # if the tmparr is full, yield it and refresh.
            pos += loc.shape[1]
            if pos >= CHUNKSIZE:

                # trim extra, fill empty cells to N (78), and convert map to array.
                tmparr = tmparr[:, :pos]
                tmparr[tmparr == 0] = 78
                tmpmap = np.array(tmpmap, dtype=np.uint64)
                self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
                yield tmparr, tmpmap

                # reset
                tmparr = np.zeros((nsamps, CHUNKSIZE + 50_000), dtype=np.uint8)
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
        self._get_sorted_loci_chunks()
        self._get_snppad()
        self._get_snames()
        self._init_datasets()
        self._fill_database()
        self._write_stats()


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7
    ip.set_log_level("DEBUG")#, log_file="/tmp/test.log")

    DATA = ip.load_json("/home/deren/Documents/ipyrad/sra-fastqs/cyatho.json")
    # DATA = ip.load_json("/tmp/TEST5.json")
    # uncomment to include the ref
    # DATA.hackers.exclude_reference = False
    print("EXCLUDE REF=", DATA.hackers.exclude_reference)

    with ip.Cluster(4) as ipyclient:
        step = Step7(DATA, ipyclient=ipyclient, force=True, quiet=False)
        step._write_cluster_chunks()
        step._apply_filters_and_trimming()
        step._collect_stats()
        step._write_stats_files()

        # run it.
        db = SeqsDatabase(step.data, step.samples)
        db.run()

        # print supermatrix parts from start and end of loci file.
        with h5py.File(db.name, 'r') as _io5:
            # print first ten sites in first chunk file
            for i, n in enumerate(db.snames):
                print(n, _io5["phy"][i, :16].tobytes().decode())
            print()
            for i, n in enumerate(db.snames):
                print(n, _io5["phy"][i, -16:].tobytes().decode())
