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

from typing import Iterator, List, Tuple, Dict, TypeVar
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd
import numpy as np
from numba import njit
from ipyrad.assemble.write_outputs_new_edge_trim import Locus
from ipyrad.assemble.utils import chroms2ints

Assembly = TypeVar("Assembly")
AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)

# pylint: disable=too-many-branches, too-many-statements


@dataclass
class ChunkProcess:
    data: Assembly
    """: Assembly object with paths and params."""
    samples: Dict[str, 'SampleSchema']
    """: Dict of samples in Step7 object."""
    chunksize: int
    """: int with max number of alignment in chunk file."""
    chunkfile: Path
    """: Chunk of alignments file."""
    nsamples: int = 0
    """: Number of samples that is updated given drop_ref."""
    filters: pd.DataFrame = None
    """: DataFrame with filters applied to each locus."""
    stats: Dict[str, Dict] = None
    """: Dict of Dicts storing polymorphism statistics."""
    loci_passed_filters: List[str] = field(default_factory=list)
    """: List of loci that will be written to disk."""
    pnames: Dict[str,str] = None
    """: Dict of padded names for loci file."""
    snames: List[str] = None
    """: Alphanumerically sorted names but with 'reference' on top if present."""
    shift_tables: List = field(default_factory=list)
    """: ...."""
    start_lidx: int = 0
    """: Locus index at the start of this chunk."""

    def __post_init__(self):
        # fill some attrs
        self._get_lidx()
        self._get_snames()
        self._get_pnames()

        # filters dataframe is size of nloci (chunksize)
        self.filters = pd.DataFrame(
            index=range(self.chunksize),
            columns=["dups", "minsamp", "maxind", "maxvar", "maxshared"],
            data=False,
            dtype=bool,
        )

        # sets stats defaults, using nsamples that will be written
        self.stats = {
            'nbases': 0,
            'sample_cov': {i: 0 for i in self.samples},
            'locus_cov': {i: 0 for i in range(1, len(self.samples) + 1)},
            'var_sites': {i: 0 for i in range(2000)},
            'pis_sites': {i: 0 for i in range(2000)},
            'var_props': {},
            'pis_props': {},
        }

    # ======== INIT FUNCTIONS =================================== ##
    def _get_lidx(self) -> None:
        """Get raw locus index at start of this chunkfile."""
        self.start_lidx = int(self.chunkfile.name.split("-")[-1][:-4])

    def _get_snames(self) -> None:
        """Get sorted names with 'reference' on top if present."""
        self.snames = sorted(self.samples)
        if self.data.drop_ref:
            self.snames.remove("reference")
        elif self.data.is_ref:
            self.snames.remove("reference")
            self.snames = ['reference'] + self.snames

    def _get_pnames(self) -> None:
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
        # get length of names for padding.
        longname = max(len(i) for i in self.snames)
        self.pnames = {i: i.ljust(longname + 5) for i in self.snames}
        self.pnames["snpstring"] = "//".ljust(longname + 5)

    # ==== MAIN FUNCTIONS ======================================== ##
    def _iter_loci(self) -> Iterator[Locus]:
        """Generator to yield loci names, seqs, and nidxs (meta strings).

        Iterates over the chunk file yielding each time it reaches the
        end of a locus (which are separated by //). Names and sequences
        are extracted from the fasta data, as well as 'nidx' strings
        which record the locus name and position (reference info).
        Also checks whether a name exists >1 time in locus, indicating
        a poor clustering that led to duplicate.
        """
        lidx = self.start_lidx
        with open(self.chunkfile, 'r', encoding="utf-8") as io_chunk:
            names = []
            nidxs = [] # reference mapping info
            seqs = []
            store_sequence = False

            for line in io_chunk:

                # yield the current locus info and reset
                if line[0] == "/":
                    if names:
                        is_dup = len(names) != len(set(names))
                        seqarr = np.array(seqs).astype(np.uint8)
                        yield Locus(lidx, self.data, names, nidxs, seqarr, is_dup=is_dup)

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
                yield Locus(lidx, self.data, names, nidxs, seqarr, is_dup=is_dup)

    def run(self):
        """Fills .loci_passed_filters and .stats by applying trimming
        and filters to all loci in _iter_loci.
        """
        # keep lists of variation among loci that passed filtering for
        # making histograms at the end.
        pis_list = []
        var_list = []

        # yields only the data of samples in data.samples from each locus
        # INCLUDING the reference if is_ref even if excluding later.
        for locus in self._iter_loci():

            # filter if duplicate samples are present in locus
            if locus.is_dup:
                self.filters.loc[locus.lidx, "dups"] = True
                continue

            # filter if sample coverage is too low in locus
            if self._filter_minsamp_pops(locus.names):
                self.filters.loc[locus.lidx, "minsamp"] = True
                continue

            # get sample coverages that will be used for trimming.
            locus.get_sample_covs()

            # apply edge trimming to locus to update:
            # creates .tseqs as trimmed version of .seqs
            # possibly subsamples .names to remove empty samples.
            # fills .nidxs with refpos info
            # fills .trimmed with n bases trimmed from each end.
            # fills .over_trimmed boolean if sample coverage becomes too low.
            locus.run()
            if locus.over_trimmed:
                self.filters.loc[locus.lidx, "minsamp"] = True
                continue

            # [reference-only] mask the insert region made up on 'N's.
            # this is used for counting nsites to not include the insert
            # region when computing % variable sites, etc.
            # if self.data.is_pair and self.data.params.min_samples_locus > 1:
                # if self.data.drop_ref:
                    # insert_mask = np.all(seqs[1: :] == 78, axis=0)
                # else:
                    # insert_mask = np.all(seqs[:] == 78, axis=0)
                # insert_mask = seqs[:, np.invert(insert)]

            # mask PE insert region to be Ns in .tseqs
            locus.mask_inserts()

            # find all variable sites in .tseqs. At this point this can
            # call SNPs within the PE insert region (on edges usually)
            # which we will want to mask after masking inserts.
            snpsarr = self._get_snps_array(locus.tseqs)

            # [denovo-only] Store shifted positions of each SNP
            # (see docstring). This uses indels (dashes) in the orig
            # .seqs array (before masking the insert).
            shift_table = self._get_shift_table(locus, snpsarr)

            # filter is number of indels is too high.
            if self._filter_maxindels(locus.tseqs):
                self.filters.loc[locus.lidx, "maxindels"] = True
                continue

            # calculate stats on the proportion of sites polymorphic.
            # we'll store these in a bit if the locus passes remaining filters.
            nsnps = np.sum(snpsarr > 0)
            npis = np.sum(snpsarr == 2)

            # mask PE insert region when counting nsites
            if locus.pe_insert[1]:
                insert_mask_size = locus.pe_insert[1] - locus.pe_insert[0]
            else:
                insert_mask_size = 0
            nsites = locus.tseqs.shape[1] - insert_mask_size
            prop_var = nsnps / nsites
            prop_pis = npis / nsites

            # filter for max polymorphic sites per locus.
            if prop_var > self.data.params.max_snps_locus:
                self.filters.loc[locus.lidx, "maxsnps"] = True
                continue

            # filter for max shared polymorphisms per site.
            if self._filter_maxshared(locus.tseqs):
                self.filters.loc[locus.lidx, "maxshared"] = True
                continue

            # -----------------------
            # LOCUS PASSED FILTERING
            # -----------------------

            # store statistics. If drop_ref this will not
            pis_list.append(prop_pis)
            var_list.append(prop_var)
            for name in locus.names:
                if name == "reference":
                    if self.data.drop_ref:
                        continue
                self.stats['sample_cov'][name] += 1
            self.stats['locus_cov'][len(locus.names) - int(self.data.drop_ref)] += 1
            self.stats['var_sites'][nsnps] += 1
            self.stats['pis_sites'][npis] += 1
            self.stats['nbases'] += nsites

            # convert to string and store it.
            locstr = self._to_locus(locus, snpsarr)
            self.loci_passed_filters.append(locstr)
            self.shift_tables.extend(shift_table)

        # clean up large stats dicts by dropping where values = 0
        self.stats['var_sites'] = {
            i: j for (i, j) in self.stats['var_sites'].items() if j}
        self.stats['pis_sites'] = {
            i: j for (i, j) in self.stats['pis_sites'].items() if j}

        # - Writes the filtered loci to {chunkfile}.loci
        # - Writes the filters array to {chunkfile}.csv
        # if no loci passed filtering then don't write to file.
        if self.loci_passed_filters:

            # write to output chunk file.
            with open(self.chunkfile.with_suffix(".loci"), 'w', encoding="utf-8") as out:
                out.write("\n".join(self.loci_passed_filters) + "\n")

            # save dataframe of filters for loci that were removed for stats
            mask = self.filters.sum(axis=1).astype(bool).values
            self.filters.loc[mask, :].to_csv(self.chunkfile.with_suffix(".csv"))

            # calculate histograms for polymorphism stats, using evenly spaced
            # bins but with the inclusion of one extra bin >0 but very low,
            # so that invariant can be distinguished from low number variants.
            nice_bins = [0, 0.001] + [round(i, 2) for i in np.linspace(0.01, 0.25, 25)]
            mags, bins = np.histogram(var_list, bins=nice_bins)
            self.stats['var_props'] = dict(zip(bins, mags))
            mags, bins = np.histogram(pis_list, bins=nice_bins)
            self.stats['pis_props'] = dict(zip(bins, mags))

            # ...
            with open(self.chunkfile.with_suffix('.npy'), 'wb') as out:
                np.save(out, np.array(self.shift_tables, dtype=np.uint64))

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
            nsamples = len(names) - int(self.data.drop_ref)
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
        n_samples = seqs.shape[0] - int(self.data.drop_ref)
        max_heteros = int(self.data.params.max_shared_h_locus * n_samples)
        if n_heteros > max_heteros:
            return True
        return False

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
        return snpcount_numba(seqs, int(self.data.drop_ref))

    def _get_shift_table(self, locus: Locus, snpsarr: np.ndarray) -> List[List[int]]:
        """Return table with positions of SNPs in .seqs.

        Record the shift for every sample that has data for a locus,
        even if it is N or dash at the SNP site. This is stored by
        the snames sample index (Sidx), the consens locus index for
        that sample (Cidx) from the nidx string, and then the locus
        and position of the trimmed alignment (output .loci fragment).
        Finally, the shift (see below) is recorded.

        Format
        ------
        >>>  Sidx    Cidx    Locus      Pos      Shift
        >>>    0      10       0        10        0
        >>>    0     200       0       100       10
        >>>    ...
        >>>   100    150K     500K     300       100
        >>> # should store as dtype=uint64 to be safe.

        Shift
        -----
        >>>    ...............                    .=trimmed
        >>>  0 -----NNNNNNNNNNCCCCCCCCCCCCCCCCC
        >>>  1 -----NNNNNNNNNNCCCCCCCCCCCCCCCCC
        >>>  2 NNNNNNNNNNTTTTTCCCCCTCCCCCCCCCCC
        >>>  3 NNNNNNNNNNTTTTTCC--CTCCCCCCCCCCC
        >>>                        *              *=SNP
        >>>  sample 0 shift at SNP = 0 [consens has not been trimmed]
        >>>  sample 2 shift at SNP = 5 [5 bases trimmed from consens]
        >>>  sample 3 shift at SNP = 3 [5 trimmed - 2 skipped indels]

        The "shift" records how many sites a position in the alignment
        is offset from its position in the individual consens sequence.
        It is used to fetch the depths of base calls for VCF output
        from the step5 HDF5 files. Step7 trimming and alignment indels
        cause the shifts.

        To find the shift we first check the orientation and start
        from 5' or 3' end. Then, we count how many non N or dash
        bases were trimmed from that side. This is the base shift.
        Then, for each SNP we subtract from the base shift the number
        of indels between the start trim position and the SNP.
        """
        # for each SNP in locus.
        shift_table = []
        pxs = np.where(snpsarr != 0)[0]
        for pos in pxs:
            for nidx, name in enumerate(locus.names):
                sidx = self.snames.index(name)
                nidxstring = locus.nidxs[nidx]
                cidx = int(nidxstring[:-2])
                ori = nidxstring[-1]

                # FIXME
                if ori == "-":
                    raise NotImplementedError('todo')
                else:
                    # get baseshift - 10 (N padding) from SEQS
                    base_shift = np.argmax((locus.seqs[nidx] != 78) & (locus.seqs[nidx] != 45))
                    base_shift -= 10

                    # get indel shift from TSEQS
                    # subtract for each indels between start and SNP
                    indels = np.sum(locus.tseqs[nidx, :pos] == 45)

                    # get shift and store full record
                    shift = base_shift - indels
                    shift_table.append([sidx, cidx, locus.lidx, pos, shift])
        return shift_table

    def _to_locus(self, locus: Locus, snpsarr: np.ndarray) -> str:
        """Return a locus as a .loci formatted str.

        This is a tmpfile, not the final .loci output. It includes
        additional metadata that is stored in the snpstring.

        The most important thing here is the updating of the 'nidx'
        metadata, which records the scaffold/locus ID and POSITION
        where this locus starts on that scaffold. Because seqs has
        been trimmed, we need to use the information in trimmer to
        update these position coordinates.

        snpstring
        ---------
        denovo: |{...}|{lefttrim,righttrim;}
        ref:    |{chromint}:{chromname}:{start}-{end}|
        """
        # convert uint8: 0->32 " ", 1->46 '.', 2->42 '*'
        snpsarr[snpsarr == 0] = 32
        snpsarr[snpsarr == 1] = 45
        snpsarr[snpsarr == 2] = 42
        snpstring = bytes(snpsarr).decode()

        # [denovo]: ...
        if not self.data.is_ref:
            refpos = str(locus.lidx)
            nidxstring = ";".join(locus.nidxs)

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
            # TODO: IF REFERENCE?
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

        # [applies to both denovo and reference]
        # write sample data for the locus in sorted name order
        locdict = {}
        for idx, name in enumerate(locus.names):
            sidx = self.snames.index(name)
            name_space = self.pnames[name]

            # convert terminal dashes into Ns
            # ----AAATTTCCCCNNNNN--AAATTTCCCC---- (before)
            # NNNNAAATTTCCCCNNNNN--AAATTTCCCCNNNN (after)
            seq = locus.tseqs[idx]
            # print("-" in seq, seq[:10])
            if 45 in seq:
                nondash = np.where(seq != 45)[0]
                left = nondash.min()
                right = nondash.max()
                seq[:left] = 78
                seq[right + 1:] = 78
            seq = seq.tobytes().decode()
            locdict[sidx] = f"{name_space}{seq}"

        # convert dict to an ordered list
        loclist = [locdict[i] for i in sorted(locdict)]

        # add SNP string and NIDX info to the end of the locus list
        name_space = self.pnames["snpstring"]
        loclist.append(f"{name_space}{snpstring}|{refpos}|{nidxstring}|")
        return "\n".join(loclist)

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
def snpcount_numba(seqs: np.ndarray, rowstart: int=0) -> np.ndarray:
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
