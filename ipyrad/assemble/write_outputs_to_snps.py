#!/usr/bin/env python

"""Loads processed .loci chunks and writes to HDF5 snps database.

Datasets in snps.hdf5
---------------------
- snps
- snpsmap
- genos
- pseudoref

Format of snps.hdf5 'snpsmap'
----------------------------
- 1-indexed total locus counter.
- 0-indexed index of nSNPs on this locus.
- 0-indexed position of SNP on this locus.
- 1-indexed index of scaffold (same as 1-index locus-index if denovo).
- 1-indexed index of SNP on reference (just zeros if denovo).

Example snpsmap
---------------
>>> [1, 0, 10, 1, 10000],
>>> [1, 1, 22, 1, 10012],
>>> [1, 2, 65, 1, 10065],
>>> [2, 0, 80, 1, 20000],  # locus 2 has only 1 SNP, at position 20K on scaff 1.
>>> [3, 0, 25, 1, 30000],
>>> [3, 1, 75, 1, 30050],
>>> ...
>>> [10000, 3, 88, 100, 200]  # scaff 100 is small and so SNP pos is low.

Note
----
As much as I don't like it, the loci and scaffolds in snpsmap are
1-indexed for compatibility reasons with other formats like VCF,
which do not allow a 0-indexed scaffold or positions.
"""

from typing import Iterator, Dict, List, Tuple, TypeVar
from pathlib import Path
import h5py
import numpy as np
from numba import njit
from ipyrad.core.schema import Project
from ipyrad.assemble.utils import get_fai_values

Assembly = TypeVar("Assembly")

# used to resolve ambiguous alleles in jit'd funcs: e.g., 82 -> (71, 65) 
AMBIGS = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    ], dtype=np.uint8)

# used when making genos calls for each sample from ALTs.
AMBIGS_FULL = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    [78,  9,  9],
    [45,  9,  9],
    [67, 67, 67],
    [65, 65, 65],
    [84, 84, 84],
    [71, 71, 71],
    ], dtype=np.uint8)


class SnpsDatabase:
    def __init__(self, data: Assembly):
        self.data = data
        """: Assembly object."""
        self.name = Path(data.stepdir) / f"{data.name}.snps.hdf5"
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
        self.snppad: int = 0
        """: Length of name region in .loci file before sequences start."""

    def _get_sorted_loci_chunks(self) -> None:
        """Gather all loci bits and sort by filename. And get snppad len."""
        paths = sorted(
            Path(self.data.tmpdir).glob("*.loci"),
            key=lambda x: int(str(x).rsplit("-", 1)[-1][:-5]),
        )
        with open(paths[0], 'r', encoding="utf-8") as indata:
            self.snppad = indata.readline().rfind(" ") + 1
        return paths

    def _init_datasets(self) -> None:
        """Create the datasets in the snps database.

        The order of samples in the snps database is alphanumeric
        except for the 'reference' sample, if included as a sample
        using the hackers option, in which case it is listed first.
        """
        with h5py.File(self.name, 'w') as io5:

            # stores the snps sequence array.
            io5.create_dataset(
                name="snps",
                shape=(self.proj.assembly_stats.nsamples, self.proj.assembly_stats.nsnps),
                dtype=np.uint8,
            )
            # stores the positions of SNPs on loci and scaffolds.
            snpsmap = io5.create_dataset(
                name="snpsmap",
                shape=(self.proj.assembly_stats.nsnps, 5),
                dtype=np.uint64,
            )
            # stores the REF and alleles for writing ALT/REF in VCF
            io5.create_dataset(
                name="alts",
                shape=(self.proj.assembly_stats.nsnps, 4),
                dtype=np.uint8,
            )
            # stores genotype calls (ALTS) as (0,0, 0,1, 0,2, etc) where 
            # for reference data 0 is the REF allele, else it is the most
            # common allele in denovo.
            io5.create_dataset(
                name="genos",
                shape=(self.proj.assembly_stats.nsamples, self.proj.assembly_stats.nsnps, 2),
                dtype=np.uint8,
            )

            # alphanumeric name order except 'reference' on top.
            self.snames = sorted(self.data.samples)

            # drop the reference sample, and add it back as first sample.
            # EVEN if we plan to exclude it from the sample data later,
            # we will still keep it in the array for now so that we can
            # write it to the HDF5 as the reference seq. NOTE: this is
            # different from how the seqs.hdf5 writing code works.
            self.snames.remove("reference")
            self.snames = ["reference"] + self.snames
            self.sidxs = {sample: i for (i, sample) in enumerate(self.snames)}

            # store the scaffold names and lengths. The name strings are
            # now always stored as strings in HDF5 attrs.
            if self.data.is_ref:
                io5["scaffold_lengths"] = get_fai_values(self.data, "length")
                io5["scaffold_names"] = get_fai_values(self.data, "scaffold").astype("S")

            # store the name of reference genome file.
            if self.data.is_ref:
                snpsmap.attrs["reference"] = Path(self.data.params.reference_sequence).name
            else:
                snpsmap.attrs["reference"] = "pseudoref"

            # store ordered names and column labels
            if self.drop_ref:
                snpsmap.attrs["names"] = [i for i in self.snames if i != "reference"]
            else:
                snpsmap.attrs["names"] = self.snames
            snpsmap.attrs["columns"] = ["locus", "locidx", "locpos", "scaf", "scafpos"]
            snpsmap.attrs["indexing"] = [1, 0, 0, 1, 1]

    def _iter_loci(self) -> Iterator[Tuple[Dict[str,str], Tuple[int,int,int]]]:
        """Yields loci from each ordered .loci file until empty.
        """
        for locfile in self._get_sorted_loci_chunks():
            with open(locfile, 'r', encoding="utf-8") as indata:

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
                        chrom, pos0, pos1 = 0, 0, 0
                        if self.data.is_ref:
                            line_end = line.split("|")[1]
                            chrom = int(line_end.split(":")[0])
                            pos0, pos1 = (
                                int(i) for i in line_end.split(":")[1].split(",")[0].split("-")
                            )

                        # end of locus, yield the dict.
                        yield names_to_seqs, (chrom, pos0, pos1)

    def _iter_supermatrix_chunks(self) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Yields chunks of supermatrix and snpsmap.
        """
        # size of chunks that will be yielded, with some overflow of columns
        chunk_size = 10_000
        nsamps = self.proj.assembly_stats.nsamples

        # make array to fill w/ overflow on column size.
        tmparr = np.zeros((nsamps, chunk_size + 5_000), dtype=np.uint8)
        tmpref = np.zeros((chunk_size + 5_000), dtype=np.uint8)
        tmpmap = []

        # iterate over loci from .loci files incrementing the locus every
        # time so that the locus IDs in snps.hdf5 match those in seqs.hdf5.
        # keep track of position in array to send a chunk when full.
        locidx = 1
        pos = 0
        for names_to_seqs, positions in self._iter_loci():

            # get locus dict and ref positions
            chrom, pos0, _ = positions

            # get snpstring
            snps_mask = np.array(list(names_to_seqs.pop("snpstring"))) != " "

            # if no snps are present the skip this locus
            nsnps = np.sum(snps_mask)
            if not nsnps:
                locidx += 1
                continue

            # order the dict by alphanumeric names except with reference first.
            ordered_seqs = (names_to_seqs[i] for i in self.snames if i in names_to_seqs)

            # convert locus to uint8 array and subselect only SNP sites.
            loc = np.array([list(i) for i in ordered_seqs], dtype=np.uint8)
            snps_sites = loc[:, snps_mask]
            snps_pos = np.where(snps_sites)[0]

            # store the reference sequence if present
            if self.data.is_ref:
                tmpref[pos: pos + snps_sites.shape[1]] = snps_sites[0, :]

            # enter snps data into the tmparr: (0, 'ref'), (1, 'A'), (2, 'B')...
            for tmp_sidx, name in enumerate(names_to_seqs):
                pos_slice = slice(pos, pos + snps_sites.shape[1])
                if self.drop_ref:
                    if name == "reference":
                        continue
                    global_sidx = self.sidxs[name] - 1
                    tmparr[global_sidx, pos_slice] = snps_sites[tmp_sidx - 1]
                else:
                    global_sidx = self.sidxs[name]
                    tmparr[global_sidx, pos_slice] = snps_sites[tmp_sidx]

            # enter position data into the tmpmap
            for idx in range(nsnps):
                spos = snps_pos[idx]
                tmpmap.append((locidx, idx, spos, chrom, pos0 + spos))

            # advance locus counter and supermtarix position index
            locidx += 1
            pos += snps_sites.shape[1]

            # if the tmparr is full, yield it and refresh.
            if pos >= chunk_size:

                # trim extra, fill empty cells to N (78), and convert map to array.
                tmpref = tmpref[:pos]
                tmparr = tmparr[:, :pos]
                tmparr[tmparr == 0] = 78
                tmpmap = np.array(tmpmap, dtype=np.uint64)
                self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
                yield tmparr, tmpmap, tmpref

                # reset
                tmparr = np.zeros((nsamps, chunk_size + 5_000), dtype=np.uint8)
                tmpmap = []
                pos = 0

        # yield the final leftover chunk that did not reach chunk_size.
        if pos:
            tmpref = tmpref[:pos]            
            tmparr = tmparr[:, :pos]
            tmparr[tmparr == 0] = 78
            tmpmap = np.array(tmpmap, dtype=np.uint64)
            self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
            yield tmparr, tmpmap, tmpref

    def _fill_database(self) -> None:
        """Iterates over chunks of data entering into the H5 database.
        """
        start_arr = 0
        start_map = 0
        with h5py.File(self.name, 'a') as io5:
            for arr_chunk, map_chunk, ref_chunk in self._iter_supermatrix_chunks():
                # get end positions
                end_arr = start_arr + arr_chunk.shape[1]
                end_map = start_map + map_chunk.shape[0]

                # insert to dataset array
                io5["snps"][:, start_arr:end_arr] = arr_chunk
                io5["snpsmap"][start_map:end_map] = map_chunk
                if self.data.is_ref:
                    io5["alts"][start_arr:end_arr, 0] = ref_chunk

                # increment start positions
                start_arr += arr_chunk.shape[1]
                start_map += map_chunk.shape[0]

    def _fill_alts(self):
        """get REF, ALT genotypes for use in VCF and GENOS calls.

        """
        with h5py.File(self.name, 'a') as io5:

            # pass a chunk of phy seq array to jit'd function
            phycols = io5["snps"].shape[1]  # pylint: disable=no-member
            for chunk in range(0, phycols, 10_000):

                # get alts from observations in chunk of snps
                cslice = slice(chunk, chunk + 10_000)
                alts = get_ordered_alts(io5["snps"][:, cslice], AMBIGS)

                # if ref, shift to make ref allele first, in chunks.
                if self.data.is_ref:
                    ref = io5["alts"][cslice, 0]
                    alts = get_reordered_alts_by_reference(alts, ref)

                # enter the altschunk into 'alts'
                io5["alts"][cslice, :] = alts

    def _fill_genos(self):
        """Fill genos as indices of REF and ALT alleles.
        """
        with h5py.File(self.name, 'a') as io5:

            # pass a chunk of phy seq array to jit'd function
            phycols = io5["snps"].shape[1]  # pylint: disable=no-member
            for chunk in range(0, phycols, 10_000):

                # slice to get chunk
                cslice = slice(chunk, chunk + 10_000)
                chunk = io5["snps"][:, cslice]
                alts = io5["alts"][cslice, :]

                # get genos and enter into hdf5
                genos = get_genos(chunk, alts, AMBIGS_FULL)
                io5["genos"][cslice, :, :] = genos

    def _write_stats(self):
        """Write missing values in matrix to stats file, and logger.

        Because this is run on an engine the logger will report what
        is written with print.
        """
        with h5py.File(self.name, 'r') as io5:
            shape = io5["snps"].shape               # pylint: disable=no-member
            size = io5["snps"].size                 # pylint: disable=no-member
            perc_missing = self.nmissing / size
        print(f"wrote snps.hdf5, shape={shape}, "
            f"missing={perc_missing:.2f}, path={self.name}.")

        # write stats to the output file
        stats_file = Path(self.data.stepdir) / "s7_assembly_stats.txt"
        with open(stats_file, 'a', encoding="utf-8") as out:
            outstr = (
                f"snps matrix size: {shape}, "
                f"{perc_missing:.2f}% missing sites.\n"
            )
            out.write(outstr)

    def run(self):
        """Initialize HDF5 database and then fill it in chunks."""
        self._init_datasets()
        self._fill_database()
        self._write_stats()
        self._fill_alts()
        self._fill_genos()

@njit
def get_ordered_alts(
    phy: np.ndarray, cons: np.ndarray) -> np.ndarray:
    """Find all observed bases at each site in order of frequency.

    This is used when filling VCF to fill the ALTS index.
    """
    # array to store observed genos for writing to VCF.
    alt_refs = np.zeros((phy.shape[1], 4), dtype=np.uint8)

    # column 0 will have observation, column 1 defaults to '.'
    alt_refs[:, 1] = 46

    # iterate over columns of the matrix
    for cidx in range(phy.shape[1]):
        # count occurrences of each uint8 character.
        fcounts = np.zeros(111, dtype=np.uint64)
        counts = np.bincount(phy[:, cidx])
        fcounts[:counts.shape[0]] = counts

        # remove Ns and -s
        fcounts[78] = 0
        fcounts[45] = 0

        # expand columns ambiguous bases to A,C,T,G
        for aidx in range(cons.shape[0]):
            nbases = fcounts[cons[aidx, 0]]
            for _ in range(nbases):
                fcounts[cons[aidx, 1]] += 1
                fcounts[cons[aidx, 2]] += 1                
            fcounts[cons[aidx, 0]] = 0

        # get next most frequent, save to alt_refs, and clear from fcounts.
        idx = 0
        while 1:
            most_common = np.argmax(fcounts)
            if most_common:
                alt_refs[cidx, idx] = most_common
                fcounts[most_common] = 0
                idx += 1
            else:
                break
    return alt_refs

@njit
def get_reordered_alts_by_reference(
    alts: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Returns ALTs so that the reference is first, and other 
    observations come after.
    """
    # create new empty array
    new_alts = np.zeros(alts.shape, dtype=np.uint8)

    # at all sites where pseudo 0 matches reference, leave it
    matched = alts[:, 0] == ref[:]      # [True, True, False, ...]
    midxs = np.where(matched)[0]        # [0, 1, ...]
    new_alts[midxs, :] = alts[midxs]

    # at other sites, shift order so ref is first
    nidxs = np.where(~matched)[0]          # [20, 33, 50, 88, ...]
    for idx in nidxs:
        pre = alts[idx, :]                 # [84, 67, 0, 0]
        new = ref[idx]                     # 67
        new_alts[idx, 0] = new             # [67, ...]

        # get index that was already filled
        midx = np.where(pre == new)[0][0]  # 1

        # fill other cells of new alts with remaining
        pre_idx = 0
        for newidx in range(1, pre.size):
            if pre_idx != midx:
                new_alts[idx, newidx] = pre[pre_idx]
                pre_idx += 1
    return new_alts

@njit
def get_genos(
    snps: np.ndarray, alts: np.ndarray, cons: np.ndarray) -> np.ndarray:
    """Return genotype as 0/1/2/3 or 9 for missing as, e.g., (0, 1).

    """
    # create genos array for the chunk.
    genos = np.zeros((snps.shape[0], snps.shape[1], 2), dtype=np.uint8)

    # fill an array for the two resolutions of each base.
    snps1 = np.zeros(snps.shape, dtype=np.uint8)
    snps2 = np.zeros(snps.shape, dtype=np.uint8)    

    # iterate over samples b/c can't index 2-d in jit.
    for bidx in range(cons.shape[0]):
        for sidx in range(snps.shape[0]):
            base = cons[bidx, 0]                   # 83
            mask = snps[sidx] == base
            snps1[sidx][mask] = cons[bidx, 1]
            snps2[sidx][mask] = cons[bidx, 2]
            # snps1[sidx, mask] = cons[bidx, 1]      # 71
            # snps2[sidx, mask] = cons[bidx, 2]      # 67
            # [67, 67, 83, 83, 71, 67, 9, 9, ...]

    # for each site convert base to index of alt for that site.
    for aidx, arr in enumerate((snps1, snps2)):
        for cidx in range(arr.shape[1]):
            alt = alts[cidx]                     # [71, 67, 0, 0]
            for sidx in range(arr.shape[0]):
                xpos = np.where(alt == arr[sidx, cidx])[0]
                if xpos.size:
                    genos[sidx, cidx, aidx] = xpos[0]
                else:
                    genos[sidx, cidx, aidx] = 9
    return genos


AMBIGS_FULL = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    [78,  9,  9],
    [45,  9,  9],
    [67, 67, 67],
    [65, 65, 65],
    [84, 84, 84],
    [71, 71, 71],
    ], dtype=np.uint8)



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    TEST = ip.load_json("/tmp/TEST5.json")
    PROJ = Project.parse_file(TEST.json_file)

    # pretend this has been made into a Step7 class object
    TEST.tmpdir = "/tmp/TEST5_tmp_outfiles"
    TEST.stepdir = "/tmp/TEST5_outfiles"
    TEST.samples["reference"] = ip.core.schema.SampleSchema(name="reference")

    # uncomment both of these to include the ref
    TEST.hackers.exclude_reference = False
    PROJ.assembly_stats.nsamples += 1
    print("EXCLUDE REF=", TEST.hackers.exclude_reference)

    # run it.
    db = SnpsDatabase(TEST, PROJ)
    db.run()

    # print supermatrix parts from start and end of loci file.
    with h5py.File(db.name, 'r') as io5:
        names = io5["snpsmap"].attrs["names"][:]
        snps = io5["snps"][:]
        genos = io5["genos"][:]

        print(names)
        # print first ten sites in first chunk file
        for i, name in enumerate(names):            
            print(f"{names[i]:<10}", snps[i, :16].tobytes().decode())            
        print(genos[-1, :16])

        # print()
        # # print last ten sites in last chunk file.
        # for i, name in enumerate(names):
        #     print(f"{names[i]:<10}", snps[i, -16:].tobytes().decode())
