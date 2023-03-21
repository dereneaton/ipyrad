#!/usr/bin/env python

"""Loads processed .loci chunks and writes to HDF5 snps database.

Attrs in snps.hdf5 
------------------
- version : int : 1 is current. Older will not have this attr.
- names: List[str] : alphanumerically sorted sample names but with 
    'reference ' as the first sample if present.
- reference: str : filename of reference fasta file.

Datasets in snps.hdf5
---------------------
- snps    : shape=(nsamples, nsnps)    : arr of diploid geno calls 
- snpsmap : shape=(nsnps, 5)           : arr of SNP positions
- genos   : shape=(nsamples, nsnps, 2) : arr of geno calls (0,1,2)
- alts    : shape=(nsnps, 4)           : arr of alt calls (...)

Attrs of the snpsmap dataset
-----------------------------
- scaffold_names: List[str]
- scaffold_lengths: List[int]
- columns: List[str] : ["loc_idx", "loc_snp_idx", "loc_snp_pos", "scaf_idx", "scaf_pos"]
- indexing: List[int] : shows that each columns is indexed from 0.

Format of data snps.hdf5 'snpsmap'
----------------------------------
- 0-indexed total locus counter.
- 0-indexed index of nSNPs on this locus.
- 0-indexed position of SNP on this locus.
- 0-indexed index of scaffold name (.attrs['scaffold_names'][idx])
- 0-indexed position of SNP on scaffold (same as pos of SNP on locus if denovo)

Example snpsmap
---------------
The indexing of rows in .snpsmap corresponds to SNPs in .snps (i.e.,
the column indexing of .snps). The first three columns of .snpsmap 
describe the location of the SNP in each discrete RAD locus (i.e., 
it matches up with the .loci file). The fourth and fifth columns
describe the location of SNPs on scaffolds.

>>> [0, 0, 10, 0, 10000],
>>> [0, 1, 22, 0, 10012],
>>> [0, 2, 65, 0, 10065],
>>> [1, 0, 80, 0, 20000],
>>> [2, 0, 25, 0, 30000],
>>> [2, 1, 75, 0, 30050],
>>> ...
>>> [10000, 3, 88, 100, 200]  # scaff 100 is small and so SNP pos is low.
# Loc0, SNP0, Pos10, Scaff0, ScaffPos10000
# Loc0, SNP1, Pos22, Scaff0, ScaffPos10000
# Loc1, SNP0, Pos80, Scaff0, ScaffPos20000

Note
----
The total locus counter index and the scaffold number index are both
0-indexed as of v1. This is different from older ipyrad versions. 
Also, for denovo data, the last two columns will match the first 
and third columns, i.e., map[:, [3, 4]] == map[:, [0, 2]].
"""

from typing import TypeVar, Dict, Iterator, Tuple
from pathlib import Path
import h5py
import numpy as np
from numba import njit
from ipyrad.assemble.utils import get_fai_values
from ipyrad.assemble.write_outputs_base import DatabaseWriter

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

# size of concatenated data processed at one time in memory.
CHUNKSIZE = 10_000


class SnpsDatabaseWriter(DatabaseWriter):

    def __init__(self, data: Assembly, samples: Dict[str, "SampleSchema"]):
        # self.name = Path(data.stepdir) / f"{data.name}.snps.hdf5"
        self.name = data.outfiles['snps_database']
        """: Output HDF5 file name."""
        super().__init__(data, samples)

    def _init_datasets(self) -> None:
        """Create the datasets in the snps database.

        The order of samples in the snps database is alphanumeric
        except for the 'reference' sample which if present is listed
        first. The hackers.exclude reference does not yet apply here,
        except in that it already affect the .nsamples attribute so
        that the datasets are the right size.
        """
        with h5py.File(self.name, 'w') as io5:

            # store meta-data
            io5.attrs["names"] = self.snames
            io5.attrs["version"] = 1
            io5.attrs["nsnps"] = self.data.assembly_stats.nsnps
            if self.data.is_ref:
                io5.attrs["reference"] = Path(self.data.params.reference_sequence).name
            else:
                io5.attrs["reference"] = "pseudoref"

            # stores the snps sequence array.
            io5.create_dataset(
                name="snps",
                shape=(self.data.assembly_stats.nsamples, self.data.assembly_stats.nsnps),
                dtype=np.uint8,
            )
            # stores the positions of SNPs on loci and scaffolds.
            snpsmap = io5.create_dataset(
                name="snpsmap",
                shape=(self.data.assembly_stats.nsnps, 5),
                dtype=np.uint64,
            )
            # stores the REF and alleles for writing ALT/REF in VCF
            io5.create_dataset(
                name="alts",
                shape=(self.data.assembly_stats.nsnps, 4),
                dtype=np.uint8,
            )
            # stores genotype calls (ALTS) as (0,0, 0,1, 0,2, etc) where 
            # for reference data 0 is the REF allele, else it is the most
            # common allele in denovo.
            io5.create_dataset(
                name="genos",
                shape=(self.data.assembly_stats.nsnps, self.data.assembly_stats.nsamples, 2),
                dtype=np.uint8,
            )

            # indices of sorted sample names 
            self.sidxs = {sample: i for (i, sample) in enumerate(self.snames)}

            # store the scaffold names and lengths. The name strings are
            # now always stored as strings in HDF5 attrs.
            if self.data.is_ref:
                snpsmap.attrs["scaffold_lengths"] = get_fai_values(self.data, "length")
                snpsmap.attrs["scaffold_names"] = get_fai_values(self.data, "scaffold").astype("S")

            # store the name of reference genome file.
            snpsmap.attrs["columns"] = ["loc_idx", "loc_snp_idx", "loc_snp_pos", "scaf_idx", "scaf_pos"]
            snpsmap.attrs["indexing"] = [0, 0, 0, 0, 0]

            # store ordered names and column labels
            # if self.drop_ref:
                # snpsmap.attrs["names"] = [i for i in self.snames if i != "reference"]
            # else:

    def _iter_supermatrix_chunks(self) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Yields chunks of supermatrix, snpmap, and reference seq.

        """
        # concat size of chunks that will be yielded, with some overflow of columns
        nsamps = self.data.assembly_stats.nsamples

        # make array to fill w/ overflow on column size.
        tmparr = np.zeros((nsamps, CHUNKSIZE + 5_000), dtype=np.uint8)
        tmpref = np.zeros((CHUNKSIZE + 5_000), dtype=np.uint8)
        tmpmap = []

        # iterate over loci from .loci files incrementing the locus every
        # time so that the locus IDs in snps.hdf5 match those in seqs.hdf5.
        # keep track of position in array to send a chunk when full.
        locidx = 0
        pos = 0
        for names_to_seqs, positions in self._iter_loci():

            # get locus dict and ref positions: 
            # denovo = (0, 'RAD', 0, 0)         # the last three are fixed.
            # refere = (0, 'scaff13', 500, 800) # all are real data.
            chrom_int, _, pos0, _ = positions

            # get snpstring
            snps_mask = np.array(list(names_to_seqs.pop("snpstring"))) != " "

            # if no snps are present the skip this locus
            nsnps = np.sum(snps_mask)
            if not nsnps:
                locidx += 1
                continue

            # order the dict by alphanumeric names except with reference first.
            # ordered_seqs = (names_to_seqs[i] for i in self.snames if i in names_to_seqs)
            ordered_seqs = [names_to_seqs[i] for i in self.snames if i in names_to_seqs]

            # convert locus to uint8 array and subselect only SNP sites.
            loc = np.array([list(i) for i in ordered_seqs], dtype=np.uint8)
            snps_sites = loc[:, snps_mask]
            snps_pos = np.where(snps_mask)[0]

            # store the reference sequence if present. This is separate from
            # whether or not it will being stored as a sample.
            if self.data.is_ref:
                tmpref[pos: pos + snps_sites.shape[1]] = snps_sites[0, :]

            # enter snps data into the tmparr: (0, 'ref'), (1, 'A'), (2, 'B')...
            for tmp_sidx, name in enumerate(names_to_seqs):
                pos_slice = slice(pos, pos + snps_sites.shape[1])
                if self.data.drop_ref:
                    if name == "reference":
                        continue
                    global_sidx = self.sidxs[name] - 1
                else:
                    global_sidx = self.sidxs[name]
                tmparr[global_sidx, pos_slice] = snps_sites[tmp_sidx]

            # enter position data into the tmpmap
            for idx in range(nsnps):
                spos = snps_pos[idx]
                tmpmap.append((locidx, idx, spos, chrom_int, pos0 + spos))

            # advance locus counter and supermtarix position index
            locidx += 1
            pos += snps_sites.shape[1]

            # if the tmparr is full, yield it and refresh.
            if pos >= CHUNKSIZE:

                # trim extra, fill empty cells to N (78), and convert map to array.
                tmpref = tmpref[:pos]
                tmparr = tmparr[:, :pos]
                tmparr[tmparr == 0] = 78
                tmpmap = np.array(tmpmap, dtype=np.uint64)
                self.nmissing += np.sum((tmparr == 78) | (tmparr == 45))
                yield tmparr, tmpmap, tmpref

                # reset
                tmparr = np.zeros((nsamps, CHUNKSIZE + 5_000), dtype=np.uint8)
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

        This fills 'snps' and 'snpsmap' dsets, and 'alts' if is_ref.
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

                # if reference then use REF as the alt allele.
                if self.data.is_ref:
                    io5["alts"][start_arr:end_arr, 0] = ref_chunk

                # increment start positions
                start_arr += arr_chunk.shape[1]
                start_map += map_chunk.shape[0]

    def _fill_alts(self):
        """get REF, ALT genotypes for use in VCF and GENOS calls.

        ALT is used to fill the VCF database later, and to make the
        GENO calls. It doesn't really need to be stored in the SNPs 
        database long term though, but it is for now.
        """
        with h5py.File(self.name, 'a') as io5:

            # pass a chunk of phy seq array to jit'd function
            snpcols = io5["snps"].shape[1]  # pylint: disable=no-member
            for chunk in range(0, snpcols, CHUNKSIZE):

                # get alts from observations in chunk of snps
                cslice = slice(chunk, chunk + CHUNKSIZE)
                alts = jit_get_ordered_alts(io5["snps"][:, cslice], AMBIGS)

                # if ref, shift to make ref allele first, in chunks.
                if self.data.is_ref:
                    ref = io5["alts"][cslice, 0]
                    alts = jit_get_reordered_alts_by_reference(alts, ref)

                # enter the altschunk into 'alts'
                io5["alts"][cslice, :] = alts

    def _fill_genos(self):
        """Fill genos as indices of REF and ALT alleles.
        """
        with h5py.File(self.name, 'a') as io5:

            # pass a chunk of phy seq array to jit'd function
            snpcols = io5["snps"].shape[1]  # pylint: disable=no-member
            for chunk in range(0, snpcols, CHUNKSIZE):

                # slice to get chunk
                cslice = slice(chunk, chunk + CHUNKSIZE)
                chunk = io5["snps"][:, cslice]
                alts = io5["alts"][cslice, :]

                # get genos and enter into hdf5
                genos = jit_get_genos(chunk, alts, AMBIGS_FULL)
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
                f"{100 * perc_missing:.2f}% missing sites.\n"
            )
            out.write(outstr)

    def run(self):
        """Initialize HDF5 database and then fill it in chunks."""

        # parent class methods
        self._get_sorted_loci_chunks()
        self._get_snppad()
        self._get_snames()

        # child class methods
        self._init_datasets()
        self._fill_database()
        self._write_stats()
        self._fill_alts()
        self._fill_genos()

@njit
def jit_get_ordered_alts(arr: np.ndarray, cons: np.ndarray) -> np.ndarray:
    """Find all observed bases at each site in order of frequency.

    This is used when filling VCF to fill the ALTS index.
    """
    # array to store observed genos for writing to VCF.
    alt_refs = np.zeros((arr.shape[1], 4), dtype=np.uint8)

    # column 0 will have observation, column 1 defaults to '.'
    alt_refs[:, 1] = 46

    # iterate over columns of the matrix
    for cidx in range(arr.shape[1]):
        # count occurrences of each uint8 character.
        fcounts = np.zeros(111, dtype=np.uint64)
        counts = np.bincount(arr[:, cidx])
        fcounts[:counts.shape[0]] = counts

        # remove Ns and -s
        fcounts[78] = 0
        fcounts[45] = 0

        # expand columns ambiguous bases to A,C,T,G
        # iterate over alleles in AMBIGS
        for aidx in range(cons.shape[0]):
            # get the counts of this allele
            nbases = fcounts[cons[aidx, 0]]
            # for each count, add one to its two resolved bases.
            # print(cidx, arr[:, cidx], cons[aidx, 0].tobytes(), nbases)
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
def jit_get_reordered_alts_by_reference(alts: np.ndarray, ref: np.ndarray) -> np.ndarray:
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
def jit_get_genos(snps: np.ndarray, alts: np.ndarray, cons: np.ndarray) -> np.ndarray:
    """Return genotype as 0/1/2/3 or 9 for missing as, e.g., (0, 1).

    why do 
    """
    # create genos array for the chunk shape=(nsnps, nsamples, 2)
    genos = np.zeros((snps.shape[1], snps.shape[0], 2), dtype=np.uint8)

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
            # print(alt)
            for sidx in range(arr.shape[0]):
                xpos = np.where(alt == arr[sidx, cidx])[0]
                if xpos.size:
                    genos[cidx, sidx, aidx] = xpos[0]
                else:
                    genos[cidx, sidx, aidx] = 9
    return genos



if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7    
    ip.set_log_level("INFO")#, log_file="/tmp/test.log")

    DATA = ip.load_json("/home/deren/Documents/ipyrad/sra-fastqs/cyatho.json")
    # DATA = ip.load_json("/tmp/small.json")

    # uncomment to include the ref
    DATA.hackers.exclude_reference = True
    print("EXCLUDE REF=", DATA.hackers.exclude_reference)

    # run it.
    with ip.Cluster(4) as ipyclient:
        step = Step7(DATA, ipyclient=ipyclient, force=True, quiet=False)
        step._write_cluster_chunks()
        step._apply_filters_and_trimming()
        step._collect_stats()
        step._write_stats_files()

        db = SnpsDatabase(step.data, step.samples)
        db.run()

        # # print supermatrix parts from start and end of loci file.
        # with h5py.File(db.name, 'r') as IO5:
        #     NAMES = IO5["snpsmap"].attrs["names"][:]
        #     SNPS = IO5["snps"][:]
        #     GENOS = IO5["genos"][:]
        #     ALTS = IO5["alts"][:]

        #     print(NAMES)
        #     # print first ten sites in first chunk file
        #     for i, n in enumerate(NAMES):            
        #         print(f"{NAMES[i]:<10}", SNPS[i, :16].tobytes().decode())            
        #     print(ALTS[:16])
        #     # print(GENOS[-1, :16])
            
