#!/usr/bin/env python

"""Class run on remote engine to filter/trim loci in step7.

If denovo:
    - SNPs are called among the samples in data.samples.
    - Trimming removes low coverage edges and RE overhangs.
    - Positional information is shifted for trimming.

If reference w/ hackers.exclude_reference=False:
    - SNPs are called among the samples including ref data.
    - Trimming removes low coverage edges and RE overhangs.
    - Positional information is shifted for trimming.    
    - stats on coverage and polymorphism include the ref.    
    - Loci outputs includes all samples + reference.

If reference w/ hackers.exclude_reference=True:
    - SNPs are called among the samples but NOT including ref data.
    - Trimming removes low coverage edges and RE overhangs.
    - Positional information is shifted for trimming.
    - stats on coverage and polymorphism do not include ref.    
    - Loci outputs includes all samples + reference (because we still
    need the ref information for filling VCF later. It is still good
    to have it in .loci for visual validation. But ref seq will be 
    excluded from the other output files as a sample.)
"""

from typing import Iterable, List, Tuple, Dict, TypeVar
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from numba import njit
from ipyrad.assemble.write_outputs_new_edge_trim import EdgeTrimmer
from ipyrad.assemble.utils import chroms2ints

Assembly = TypeVar("Assembly")
AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)

# pylint: disable=too-many-branches, too-many-statements


@dataclass
class ChunkProcess:
    data: Assembly
    chunksize: int
    chunkfile: str

    nsamples: int = 0
    """: Number of samples that is updated given drop_ref."""
    drop_ref: bool = False
    """: boolean for whether to include ref sample when finding SNPs"""
    filters: pd.DataFrame = None
    """: DataFrame with filters applied to each locus."""
    stats: Dict[str, Dict] = None
    """: Dict of Dicts storing polymorphism statistics."""
    loci_passed_filters: List[str] = field(default_factory=list)
    """: List of loci that will be written to disk."""
    pnames: Dict[str,str] = None
    """: Dict of padded names for loci file."""

    def __post_init__(self):
        # get boolean for whether 'reference' will be written to files.
        self.drop_ref = self.data.is_ref & self.data.hackers.exclude_reference

        # get dict with padding for the samples that will be written
        self.pnames = self._get_padded_names()

        # filters dataframe is size of nloci (chunksize)
        self.filters = pd.DataFrame(
            index=range(self.chunksize),
            columns=["dups", "minsamp", "maxind", "maxvar", "maxshared"],
            data=False,
            dtype=bool,
        )

        # get tmp list of sample names that will be written to file.
        # This is may or may not include reference, since when drop_ref
        # is on the ref should not be included in stats.
        samples = list(self.data.samples)
        if self.drop_ref:
            samples.remove("reference")

        # sets stats defaults, using nsamples that will be written
        self.stats = {
            'nbases': 0,
            'sample_cov': {i: 0 for i in samples},
            'locus_cov': {i: 0 for i in range(1, len(samples) + 1)},
            'var_sites': {i: 0 for i in range(2000)},
            'pis_sites': {i: 0 for i in range(2000)},
            'var_props': {},
            'pis_props': {},
        }

    def _get_padded_names(self) -> Dict[str,str]:
        """Return a dict mapping names to padded left-aligned names
        and padding for the SNP line.

        Name padding is created relative to ALL samples that could 
        be present in any locus, including the reference, since if 
        present if is always written to the .loci outputs.

        Note
        ----
        A similar function to this is present in other format writers
        which may exclude the 'reference' sample and thus pad the names
        slightly differently (e.g., phy and snps outputs).
        """
        names = list(self.data.samples)
        longname = max(len(i) for i in names)
        pnames = {i: i.ljust(longname + 5) for i in names}
        pnames["snpstring"] = "//".ljust(longname + 5)
        return pnames

    def run(self):
        """Iterates over all loci in the chunkfile applying filters.

        - Writes the filtered loci to {chunkfile}.loci
        - Writes the filters array to {chunkfile}.csv
        """
        # keep lists of variation among loci that passed filtering for
        # making histograms at the end.
        pis_list = []
        var_list = []

        # yields only the data of samples in data.samples from each locus
        # INCLUDING the reference if is_ref even if excluding later.
        for lidx, names, seqs, nidxs, is_dup in self._iter_loci():

            # filter if duplicate samples are present in locus
            if is_dup:
                self.filters.loc[lidx, "dups"] = True
                continue

            # filter if sample coverage is too low in locus
            if self._filter_minsamp_pops(names):
                self.filters.loc[lidx, "minsamp"] = True
                continue

            # get seqs with edges trimmed, and names and nidxs potentially
            # subsampled if a sample no longer had data after trimming,
            # which can occur sometimes with highly tiled data.
            # The `trimmer.sumlens` is used later to update positions
            # to make sure new locus sites align with depths info positions.
            trimmer = EdgeTrimmer(self.data, names, seqs, nidxs, debug=True)
            names, seqs, nidxs, over_trimmed = trimmer.run()

            # filter if sample coverage is too low after trimming edges.
            if over_trimmed:
                self.filters.loc[lidx, "minsamp"] = True
                continue

            # [denovo-only] FIXME: are the inserts Ns or -s here?...
            # ...

            # [reference-only] mask the insert region made up on 'N's.
            # this is used for counting nsites to not include the insert
            # region when computing % variable sites, etc.
            if self.data.is_pair and self.data.params.min_samples_locus > 1:
                if self.drop_ref:
                    insert_mask = np.all(seqs[1: :] == 78, axis=0)
                else:
                    insert_mask = np.all(seqs[:] == 78, axis=0)
                # insert_mask = seqs[:, np.invert(insert)]

            # filter is number of indels is too high.
            if self._filter_maxindels(seqs):
                self.filters.loc[lidx, "maxindels"] = True
                continue

            # find all variable sites in seqs.
            snpsarr = self._get_snps_array(seqs)

            # calculate stats on the proportion of sites polymorphic.
            # we'll store these in a bit if the locus passes remaining filters.
            nsnps = np.sum(snpsarr > 0)
            npis = np.sum(snpsarr == 2)
            nsites = seqs.shape[1] - insert_mask.sum()
            prop_var = nsnps / nsites
            prop_pis = npis / nsites

            # filter for max polymorphic sites per locus.
            if prop_var > self.data.params.max_snps_locus:
                self.filters.loc[lidx, "maxsnps"] = True
                continue

            # filter for max shared polymorphisms per site.
            if self._filter_maxshared(seqs):
                self.filters.loc[lidx, "maxshared"] = True
                continue

            # -----------------------
            # LOCUS PASSED FILTERING
            # -----------------------

            # store statistics. If drop_ref this will not 
            pis_list.append(prop_pis)
            var_list.append(prop_var)
            for name in names:
                if name == "reference":
                    if self.drop_ref:
                        continue
                self.stats['sample_cov'][name] += 1
            self.stats['locus_cov'][len(names) - int(self.drop_ref)] += 1
            self.stats['var_sites'][nsnps] += 1
            self.stats['pis_sites'][npis] += 1
            self.stats['nbases'] += nsites

            # convert to string and store it.
            locus = self._to_locus(names, seqs, nidxs, snpsarr, trimmer)
            self.loci_passed_filters.append(locus)

        # clean up large stats dicts by dropping where values = 0
        self.stats['var_sites'] = {
            i: j for (i, j) in self.stats['var_sites'].items() if j}
        self.stats['pis_sites'] = {
            i: j for (i, j) in self.stats['pis_sites'].items() if j}

        # if no loci passed filtering then don't write to file.
        if self.loci_passed_filters:

            # write to output chunk file.
            with open(self.chunkfile + ".loci", 'w', encoding="utf-8") as out:
                out.write("\n".join(self.loci_passed_filters) + "\n")

            # save dataframe of filters for loci that were removed for stats
            mask = self.filters.sum(axis=1).astype(bool).values
            self.filters.loc[mask, :].to_csv(self.chunkfile + ".csv")

            # calculate histograms for polymorphism stats, using evenly spaced
            # bins but with the inclusion of one extra bin >0 but very low, 
            # so that invariant can be distinguished from low number variants.
            nice_bins = [0, 0.001] + [round(i, 2) for i in np.linspace(0.01, 0.25, 25)]
            mags, bins = np.histogram(var_list, bins=nice_bins)
            self.stats['var_props'] = dict(zip(bins, mags))
            mags, bins = np.histogram(pis_list, bins=nice_bins)
            self.stats['pis_props'] = dict(zip(bins, mags))

    def _filter_maxindels(self, seqs: np.ndarray) -> bool:
        """Max size of internal indels. Denovo vs. Ref, single versus paired."""
        # get max indels for read1, read2
        inds = maxind_numba(seqs)
        if inds > self.data.params.max_indels_locus:
            return True
        return False

    def _filter_minsamp_pops(self, names: List[str]) -> bool:
        """Return True if locus is filtered by min_samples_locus.

        The format of the .populations dict is {'popname': minsamp, ...}
        and can be set in CLI by entering a popfile that is parsed,
        or in the API by setting it directly on Assembly.populations.
        """
        # default: no population information
        if not self.data.populations:
            nsamples = len(names) - int(self.drop_ref)
            if nsamples < self.data.params.min_samples_locus:
                return True
            return False

        # use populations
        minfilters = []
        for pop in self.data.populations:
            samps = self.data.populations[pop][1]
            minsamp = self.data.populations[pop][0]
            if len(set(samps).intersection(set(names))) < minsamp:
                minfilters.append(pop)
        if any(minfilters):
            return True
        return False

    def _filter_maxshared(self, seqs: np.ndarray) -> bool:
        """Return True if too many samples share a polymorphism at a site.

        This is a paralog filter. If a polymorphism is shared across
        N% of samples it is more likely to be a fixed difference at
        paralogous sites that were erroneously treated a orthologs,
        than it is likely to be a true widespread polymorphism. Or
        at least that's the idea.'
        """
        n_heteros = count_maxheteros_numba(seqs)
        n_samples = seqs.shape[0] - int(self.drop_ref)
        max_heteros = int(self.data.params.max_shared_h_locus * n_samples)
        if n_heteros > max_heteros:
            return True
        return False

    def _iter_loci(self) -> Iterable[Tuple[int, List[str], np.ndarray, List[str], bool]]:
        """Generator to yield loci names, seqs, and nidxs (meta strings).

        Iterates over the chunk file yielding each time it reaches the
        end of a locus (which are separated by //). Names and sequences
        are extracted from the fasta data, as well as 'nidx' strings
        which record the locus name and position (reference info).
        Also checks whether a name exists >1 time in locus, indicating
        a poor clustering that led to duplicate.
        """
        with open(self.chunkfile, 'r', encoding="utf-8") as io_chunk:
            names = []
            nidxs = [] # reference mapping info
            seqs = []
            lidx = 0
            store_sequence = False

            for line in io_chunk:

                # yield the current locus info and reset
                if line[0] == "/":
                    if names:
                        is_dup = len(names) != len(set(names))
                        seqarr = np.array(seqs).astype(np.uint8)
                        yield (lidx, names, seqarr, nidxs, is_dup)

                        # reset.
                        names.clear()
                        nidxs.clear()
                        seqs.clear()
                        lidx += 1

                # get sample name and mapping info
                elif line[0] == ">":
                    name, nidx = line[1:].strip().rsplit("_", 1)
                    if name not in self.data.samples:
                        continue
                    names.append(name)
                    nidxs.append(nidx)
                    store_sequence = True

                # store the sequence if the name was stored.
                elif store_sequence:
                    seqs.append(list(line.strip().encode("utf-8")))
                    store_sequence = False

            # final locus has no ending //, so if there is data then yield last
            if names:
                is_dup = len(names) != len(set(names))
                seqarr = np.array(seqs).astype(np.uint8)
                yield (lidx, names, seqarr, nidxs, is_dup)

    def _get_snps_array(self, seqs: np.ndarray) -> np.ndarray:
        """Return an array with two types of SNP calls.

        The returned array is shape=(nsites, 2) where the first
        column indicates whether or not the site is variable, while
        the second column indicates whether or not it is phylog.
        informative (minor allele is present in more than 1 sample).
        The two different records are used to make the snpstring in
        the .loci format.

        The core function is jit'd for speed.
        """
        return snpcount_numba(seqs, int(self.drop_ref))

    def _to_locus(
        self,
        names: List[str],
        seqs: np.ndarray,
        nidxs: List[str],
        snpsarr: np.ndarray,
        trimmer: EdgeTrimmer,
        ) -> str:
        """Return a locus as a .loci formatted str.

        The most important thing here is the updating of the 'nidx'
        metadata, which records the scaffold/locus ID and POSITION
        where this locus starts on that scaffold. Because seqs has
        been trimmed, we need to use the information in trimmer to
        update these position coordinates.

        snpstring
        ---------
        denovo: |{chromint}|
        ref:    |{chromint}:{chromname}:{start}-{end}|
        """
        # convert uint8: 0->32 " ", 1->46 '.', 2->42 '*'
        snpsarr[snpsarr == 0] = 32
        snpsarr[snpsarr == 1] = 45
        snpsarr[snpsarr == 2] = 42
        snpstring = bytes(snpsarr).decode()

        # [denovo]: ...
        # TODO...
        if not self.data.is_ref:
            pass

        # [reference]: all samples have the same nidx = 'scaff:start-end'
        # TODO: is the above True? check with empirical data...
        else:
            # get dict mapping {int: str} for chroms names.
            faidict = chroms2ints(self.data, keys_as_ints=True)

            # get original position of locus on reference.
            refpos = ":".join(nidxs[0].rsplit(":", 2)[-2:])
            chrom, pos = refpos.split(":")
            ostart, end = pos.split("-")
            chrom_int = int(chrom)
            chrom_name = str(faidict[chrom_int])

            # shift locus position to reflect trimmed sites from each edge.
            start = int(ostart) + trimmer.sumlens[0]
            end = int(end) - trimmer.sumlens[1]

            # validate that positions still align after trimming.
            assert end - start == seqs.shape[1], (
               f"trim positions alignment error: {end-start} != {seqs.shape[1]}")

            # shift reference position string from edge trims
            nidbits = []
            for bit in nidxs[1:]:
                bkey = []
                for cbit in bit.split(";"):
                    cidx, _, pos = cbit.split(":")

                    # start pos of sample is its consens start pos + ostart
                    # where ostart is the ref locus start pos after trim. So
                    # how far ahead of ref start does the consens read start.
                    posplus = int(pos.split("-")[0]) - int(ostart)
                    bkey.append(f"{cidx}:{posplus}")
                nidbits.append("-".join(bkey))
            # print(nidbits, "is this relevant to empirical data?")
            refpos = f"{chrom}:{chrom_name}:{start}-{end}"
            nidxstring = ",".join(nidbits)

        # write sample data for the locus
        locus = []
        for idx, name in enumerate(names):
            name_space = self.pnames[name]
            seq = seqs[idx].tobytes().decode()
            locus.append(f"{name_space}{seq}")

        # write SNP string for the locus
        name_space = self.pnames["snpstring"]
        locus.append(f"{name_space}{snpstring}|{refpos}|{nidxstring}|")
        return "\n".join(locus)

@njit
def maxind_numba(seqs: np.ndarray) -> int:
    """Return the size of the largest indel.

    Indels are only counted here as '-' sites (45), not 'N' (78).
    """
    inds = 0
    for row in range(seqs.shape[0]):
        where = np.where(seqs[row] != 45)[0]
        left = np.min(where)
        right = np.max(where)
        obs = np.sum(seqs[row, left:right] == 45)
        if obs > inds:
            inds = obs
    return inds

@njit
def snpcount_numba(seqs: np.ndarray, rowstart: int = 0) -> np.ndarray:
    """Return the SNP array (see get_snps_array docstring).

    Parameters
    ----------
    seqs: ndarray
        A locus sequence array shape (ntaxa, nsites) in np.uint8.
    rowstart: int
        Taxon row to start on. Default if 0 (iter over all taxa),
        but when excluding the reference as counting towards
        identifying variants then the first row is skipped (the
        reference sample is always first row).
    """

    # record for every site as 0, 1, or 2, where 0 indicates the site
    # is invariant, 1=autapomorphy, and 2=synapomorphy.
    snpsarr = np.zeros(seqs.shape[1], dtype=np.uint8)

    # record the frequency of minor allele to filter by maf.
    # Nah, recommend users can do this after in ipa.window_extracter...?
    # snpfreqs = np.zeros(seqs.shape[1], dtype=np.float)

    # iterate over all loci
    for site in range(seqs.shape[1]):

        # count Cs As Ts and Gs at each site (up to 65535 sample depth)
        catg = np.zeros(4, dtype=np.uint16)

        # select the site column (potentially skipping first sample if ref.)
        ncol = seqs[rowstart:, site]

        # iterate over bases in the site column recording
        for idx in range(ncol.shape[0]):
            if ncol[idx] == 67:    # C
                catg[0] += 1
            elif ncol[idx] == 65:  # A
                catg[1] += 1
            elif ncol[idx] == 84:  # T
                catg[2] += 1
            elif ncol[idx] == 71:  # G
                catg[3] += 1
            elif ncol[idx] == 82:  # R
                catg[1] += 1       # A
                catg[3] += 1       # G
            elif ncol[idx] == 75:  # K
                catg[2] += 1       # T
                catg[3] += 1       # G
            elif ncol[idx] == 83:  # S
                catg[0] += 1       # C
                catg[3] += 1       # G
            elif ncol[idx] == 89:  # Y
                catg[0] += 1       # C
                catg[2] += 1       # T
            elif ncol[idx] == 87:  # W
                catg[1] += 1       # A
                catg[2] += 1       # T
            elif ncol[idx] == 77:  # M
                catg[0] += 1       # C
                catg[1] += 1       # A

        # sort counts so we can find second most common site.
        catg.sort()

        # if invariant      [0, 0, 0, 9] -> 0
        # if autapomorphy   [0, 0, 1, 8] -> 1
        # if synapomorphy   [0, 0, 2, 7] -> 2
        if catg[2] == 0:
            pass
        elif catg[2] == 1:
            snpsarr[site] = 1
        else:
            snpsarr[site] = 2
    return snpsarr

@njit
def count_maxheteros_numba(seqs: np.ndarray) -> int:
    """Return max number of samples with a shared polymorphism.
    """
    counts = np.zeros(seqs.shape[1], dtype=np.uint16)
    for fidx in range(seqs.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(seqs[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts.max()


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")
    # from dataclasses import asdict

    data = ip.load_json("/tmp/TEST5.json")
    data.params.restriction_overhang = ("TGCAG", "CGG")
    data.hackers.exclude_reference = True

    with ip.Cluster(cores=2) as ipyclient:
        step = ip.assemble.s7_assemble.Step7(data, True, True, ipyclient)
        step._split_clusters()
        step._apply_filters_and_trimming()
        step._collect_stats()

        # proc = ChunkProcess(step.data, 100, '/tmp/TEST5_tmp_outfiles/chunk-0')
        # proc.run()

        # lidx, names, seqs, nidxs, is_dup  = next(proc._iter_loci())
        # trimmer = EdgeTrimmer(step.data, names, seqs, nidxs, debug=True)
        # names, seqs, nidxs, over_trimmed = trimmer.run()

        # print(trimmer)

        # print(list(proc.data.samples))
        # proc.run()
        # print(proc.filters.head())
    # print(asdict(proc))

    # confirm that SNPs are being called correct.
