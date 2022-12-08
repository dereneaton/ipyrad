#!/usr/bin/env python

"""Load SNPs array and SNPsMAP from ipyrad SNPs HDF5 database.

These data can be used to create subsampled SNP data sets, or for
imputation, in other ipyrad.analysis tools.

This tool is used in: pca, structure, treemix, popgen, baba

Note
----
snpsmap columns:
    0: 1-indexed scaff id
    1: 0-indexed snpidx
    2: 1-indexed scaffpos
    3: 0-indexed orig. scaffidx
    4: snpcounter]

Example
-------
>>> import ipyrad.analysis as ipa
>>> tool = ipa.snps_extracter(DATA)
>>> tool.run()
>>> ...
"""

from typing import Optional, Dict, List, Union, Tuple
from pathlib import Path
import numpy as np
import pandas as pd
from loguru import logger
from ipyrad.core.cluster import Cluster
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import jsubsample_snps, jsubsample_loci
from ipyrad.analysis.progress import ProgressBar
from ipyrad.analysis.database import SnpsDatabase

logger = logger.bind(name="ipa")


# how many cols of SNPs to load in at once from snps, genos, snpsmap
CHUNKSIZE = 10_000
STATS_HEADER = [
    "samples",
    "pre_filter_snps",
    "pre_filter_percent_missing",
    "filter_by_indels_present",
    "filter_by_non_biallelic",
    "filter_by_mincov",
    "filter_by_minmap",
    "filter_by_invariant_after_subsampling",
    "filter_by_minor_allele_frequency",
    "post_filter_snps",
    "post_filter_snp_containing_linkage_blocks",
    "post_filter_percent_missing",
]

class SNPsExtracter:
    """Extract SNP data from HDF5 after filtering.

    This tool is used to subsample and filter SNPs for linkage and
    missingness. The filtered SNPs are accessible as a numpy array.
    This tool is primarily for internal use by ipyrad-analysis, and
    is used in other tools such as ipa.pca, ipa.popgen, etc.

    Parameters
    ----------
    data: str
        path to the .snps.hdf5 file produced by ipyrad (or created by
        converting a vcf file to .snps.hdf5 using ipyrad.analysis).
    imap: dict[str, List[str]]
        dict mapping population names to a list of sample names. This
        is used to select which samples from the data will be
        included in the filtered dataset.
    minmap: dict
        dict mapping population names to a float or integer
        representing the minimum samples required to be present per
        population for a SNP to be retained in the filtered dataset.
    mincov: int or float
        the minimum samples required globally for a SNP to be retained
        in the dataset. Int is treated as a min number of samples,
        float is treated as a min proportion of samples.
    minmaf: int or float
        minimum minor allele frequency. The minor allele at a variant
        must be present across at least n% of samples (haploid base
        calls with data) at a SNP to retain the SNP in the dataset.
        The default is zero, meaning that any SNP will be retained,
        greater values require minor variants to be shared across more
        samples. An int value of 1 will remove sites heterozygous in
        a single diploid, whereas 2 would would remove a site hetero
        in two samples, or homozygous in only one sample. Float values
        take into account the proportion of samples with missing data
        in each site so that it is a proportion of the alleles present.
    ncores: int
        If ncores > 1 the work is parallelized on this many processors.
    """
    def __init__(
        self,
        data: Union[str, Path],
        imap: Optional[Dict] = None,
        minmap: Optional[Dict] = None,
        mincov: Union[float,int] = 0.0,
        minmaf: Union[float,int] = 0.0,
        # rmincov: float=0.0,
        ):

        # store params
        self.data = Path(data)
        """String file path to a snps HDF5 file."""
        self.imap = imap if imap else {}
        """Dict mapping population names to a list of existing sample names."""
        self.minmap = minmap if minmap else {i: 1 for i in self.imap}
        """Dict mapping population names to mincov numbers or proportions per population."""
        self.mincov = mincov
        """The minimum sample coverage across all samples or subsamples (if imap)."""
        self.maf = minmaf
        """The minimum frequency of the minor allele else SNP is filtered."""

        # attributes to be filled
        self.nsnps: int = None
        """Number of SNPs before filters are applied."""
        self.dbnames: np.ndarray=None
        """Array of ordered names in the HDF5 snps array."""
        self.names: List[str]=None
        """List of ordered names in the IMAP dict (may be subset of dbnames)."""
        self.sidxs: np.ndarray=None
        """Array of int indices of names in .dbnames that are in .names."""
        self.mask: np.ndarray=None
        """Array of (nsnps,) where True=filtered."""
        self.snps: np.ndarray=None
        """Array of (nsamples, nsnps) w/ diploid genotype calls as uint8."""
        self.genos: np.ndarray=None
        """Array of (nsamples, nsnps) w/ diploid genotype calls as integers."""
        self.snpsmap: np.ndarray=None
        """Array of shape=(nsnps, 2) with locus and site indices."""
        self.stats: pd.Series=None
        """DataFrame with filtering statistics."""

        self._set_names_from_imap()
        self._set_nsnps_and_check_h5_file_format()
        self._set_names_sidxs_and_dbnames()

    def _set_names_from_imap(self) -> None:
        """Fill names, nsnps, and check input H5 file."""
        self.names = []
        for _, val in self.imap.items():
            if isinstance(val, (list, tuple, np.ndarray)):
                self.names.extend(val)
            elif isinstance(val, str):
                self.names.append(val)
        # report repeated names
        for name in self.names:
            if self.names.count(name) > 1:
                raise IPyradError(f"Name {name} is repeated in imap.")

    def _set_nsnps_and_check_h5_file_format(self) -> None:
        """Check input data is proper format, and get nsnps."""
        if self.data.suffix == ".vcf":
            raise IPyradError("input should be hdf5, see the vcf_to_hdf5 tool.")
        # this will raise an error msg if using an outdated version
        with SnpsDatabase(self.data, 'r') as io5:
            self.nsnps = int(io5.attrs['nsnps'])

    def _set_names_sidxs_and_dbnames(self) -> None:
        """Check names in h5 file match those in imap dict. Also sets names."""
        # this will raise an error msg if using an outdated version
        with SnpsDatabase(self.data, 'r') as io5:

            # get all names in database as List[str]
            self.dbnames = list(io5.attrs["names"])

            # raise error if any imap sample names not in database names
            badnames = set(self.names).difference(self.dbnames)
            if badnames:
                raise IPyradError(
                    f"Samples [{badnames}] are not in data file: {self.data}")

            # if no imap names then use db names, else subset to imap
            # names but order into db names order
            if not self.names:
                self.names = self.dbnames
            else:
                self.names = [i for i in self.dbnames if i in self.names]

            # update sidxs to mask rows not in selected samples
            sidxs = sorted([self.dbnames.index(i) for i in self.names])
            self.sidxs = np.array(sidxs)

    ################################################################
    ## the main user function
    ################################################################
    def run(self, cores: int=1, log_level: str="INFO", ipyclient: "Client"=None):
        """Loads SNP data from HDF5, applies filters, and logs stats.

        The filtered statistics are saved in the .stats attribute, and
        are printed if the log_level if above the `ipa.set_log_level()`
        (default='INFO'). This operation can be parallelized for speed
        improvements on large datasets by entering a value > 1 for
        the `cores` arg.

        Parameters
        ----------
        cores: int
            Number of cores to parallelize filtering. Default=1.
        log_level: str
            verbosity of outputs to the logger (mostly used internally).
        """
        if cores == 1:
            self._run(log_level)
            return

        # run with an existing ipyclient or start and wrap a new one.
        if ipyclient is not None:
            self._run(log_level=log_level, ipyclient=ipyclient)
        else:
            with Cluster(cores=cores, logger_name="ipa") as client:
                self._run(log_level=log_level, ipyclient=client)

    def _run(self, log_level: str="INFO", ipyclient=None):
        """Parse genotype calls from HDF5 snps file.

        This runs the filtering steps to extract and filter SNP data
        and stores the snps and snpsmap as numpy arrays.

        This performs in a chunked way to avoid memory limits on super
        large datasets.
        """
        # df for reporting filter stats
        stats = pd.Series(index=STATS_HEADER, dtype=int)

        # store masks to be concatenated at end
        mask_arrs = []
        genos_arrs = []
        snpsmap_arrs = []
        snps_arrs = []
        nmissing = 0
        ntotal = 0

        # run `get_masks_chunk()` on a single core or in parallel
        if ipyclient is None:
            jobs = {}
            for start in range(0, self.nsnps, CHUNKSIZE):
                jobs[start] = self._get_masks_chunk(start)
        else:
            lbview = ipyclient.load_balanced_view()
            prog = ProgressBar({}, "SNP filtering")
            for start in range(0, self.nsnps, CHUNKSIZE):
                prog.jobs[start] = lbview.apply(self._get_masks_chunk, start)
            prog.block()
            prog.check()
            jobs = prog.jobs

        # collect results from chunked jobs
        for job in jobs:
            try:
                snpsmap, snps, genos, masks, nmiss, ntot = jobs[job].get()
            except (NameError, AttributeError):
                snpsmap, snps, genos, masks, nmiss, ntot = jobs[job]
            snpsmap_arrs.append(snpsmap)
            snps_arrs.append(snps)
            genos_arrs.append(genos)
            mask_arrs.append(masks)
            nmissing += nmiss
            ntotal += ntot

        if not mask_arrs:
            raise IPyradError("No SNPs found.")

        # concatenate masks into supermask array, and del list
        mask = np.concatenate(mask_arrs)
        self.mask = mask.sum(axis=1).astype(bool)
        self.snpsmap = np.concatenate(snpsmap_arrs)
        self.genos = np.concatenate(genos_arrs, axis=1).astype(np.uint8)
        self.snps = np.concatenate(snps_arrs, axis=1)

        # assign unique site index to every snpsmap site
        self.snpsmap[:, 1] = range(self.snpsmap.shape[0])

        # record missing pre-impute (TODO: move to ?)
        if self.genos.size:
            missing_cells = np.sum(self.genos == 9)
            missing_percent = missing_cells / self.genos.size
        else:
            missing_percent = 1.

        # report stats
        stats.samples = len(self.names)
        stats.pre_filter_snps = self.nsnps
        stats.pre_filter_percent_missing = 100 * (nmissing / ntotal)
        stats.filter_by_indels_present = mask[:, 0].sum()
        stats.filter_by_non_biallelic = mask[:, 1].sum()
        stats.filter_by_mincov = mask[:, 2].sum()
        stats.filter_by_minmap = mask[:, 3].sum()
        stats.filter_by_invariant_after_subsampling = mask[:, 4].sum()
        stats.filter_by_minor_allele_frequency = mask[:, 5].sum()
        stats.post_filter_snps = self.snpsmap.shape[0]
        stats.post_filter_snp_containing_linkage_blocks = np.unique(self.snpsmap[:, 0]).size
        stats.post_filter_percent_missing = 100 * missing_percent
        logger.log(log_level, f"filter statistics:\n{stats}")
        self.stats = stats

    ################################################################
    ## processing functions
    ################################################################
    def _get_masks_chunk(self, start: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
        """A single chunk iteration to load data and calc filters from h5.

        This could be parallelized. It only reads from h5.
        """
        chunkslice = slice(start, start + CHUNKSIZE)
        with SnpsDatabase(self.data, 'r') as io5:

            # snps is used to filter multi-allel and indel containing.
            snps = io5["snps"]
            # genos are the actual calls we want, after filtering.
            genos = io5["genos"]
            # snpsmap is the position information of SNPs on loci.
            snpsmap = io5["snpsmap"]

            # get size of the mask to create
            start = chunkslice.start
            end = min(chunkslice.stop, snps.shape[1])
            nsnps = end - start

            # select slice for this chunk and sample set
            snpsmap = snpsmap[start:end, :2]
            snps = snps[self.sidxs, start:end]
            genos = genos[start:end, self.sidxs, :].astype(np.uint8)

            # measure number of missing cells
            nmissing = np.sum(genos == 9)
            ntotal = genos.size

            # get filter masks and diploid genotypes
            masks, diplos = self._masks_filter(nsnps, snps, genos)

            # flatten the mask to a boolean for each SNP
            flat_mask = np.invert(masks.sum(axis=1).astype(bool))

            # apply mask to arrays
            snpsmap = snpsmap[flat_mask]
            snps = snps[:, flat_mask]
            diplos = diplos[:, flat_mask]
        return snpsmap, snps, diplos, masks, nmissing, ntotal

    def _masks_filter(self, nsnps: int, snps: np.ndarray, genos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return arrays with filter masks and diploid genotypes.

        This is called wtihin `_get_masks_chunk()`.
        The filter masks are not yet applied to the genotype array. The
        snps and genos arrays have already been subsampled to include
        only this chunk of snps, and only the selected rows of
        subsamples.
        """
        # masks are boolean shape=(6, chunksize)
        masks = np.zeros((nsnps, 6), dtype=bool)

        # mask0 is True if an indel is present
        masks[:, 0] = np.any(snps == 45, axis=0)

        # mask1 is True if a third or fourth allele is present.
        masks[:, 1] = np.sum(genos == 2, axis=2).sum(axis=1).astype(bool)
        masks[:, 1] += np.sum(genos == 3, axis=2).sum(axis=1).astype(bool)

        # mask2 is True if sample coverage is below mincov.
        # mask missing calls from genotype array
        genomask = np.ma.array(data=genos, mask=(genos == 9))

        # count number of non-masked haplotypes in each site [2, 2, 4, 0, ...]
        # here a zero indicates the the site is fully masked (i.e., missing)
        # for the selected set of samples.
        nhaplos = (~genomask.mask).sum(axis=2).sum(axis=1)

        # This accomodates missing haploid calls (0/9) since it counts alleles
        if isinstance(self.mincov, int):
            masks[:, 2] = nhaplos < (2 * self.mincov)
        elif isinstance(self.mincov, float):
            masks[:, 2] = nhaplos < (2 * self.mincov * len(self.sidxs))
        else:
            raise IPyradError("mincov should be an int or float.")

        # mask3 is True if any pop sample cov is below imap/minmap
        for pop, samps in self.imap.items():
            # get min cov for this pop
            mincov = self.minmap[pop]

            # get the indices for samples in this pop
            imap_sidxs = [self.names.index(i) for i in samps]

            # select the samples from the geno array
            subarr = genomask[:, imap_sidxs, :]

            # get number of haplotypes skipping masked missing
            nhaplos = (~subarr.mask).sum(axis=2).sum(axis=1)
            if isinstance(mincov, int):
                masks[:, 3] += nhaplos < (2 * mincov)
            elif isinstance(mincov, float):
                masks[:, 3] += nhaplos < (2 * mincov * len(imap_sidxs))
            else:
                raise IPyradError("minmap dictionary malformed.")

        # mask4 is True if site is not a variant w/ current sampling---
        # get the common diploid call at each site (0, 1, 2); anything
        # above this (e.g., 3, 4) necessarily includes a 2 or 3 geno
        # call, which is already filtered by the bi-allele filter.
        diplo_common = (genomask
            .sum(axis=2)
            .mean(axis=1)
            .round()
            .astype(int)
            .data
        )
        # get diploid calls at each site
        diplos = genomask.sum(axis=2).data.T

        # mask4 is sites where all samples match common geno.
        masks[:, 4] = np.all(diplo_common == diplos, axis=0)

        # mask5 is True if maf is below minmaf setting -----------
        called_0 = (genomask == 0).sum(axis=2).sum(axis=1).data
        called_1 = (genomask == 1).sum(axis=2).sum(axis=1).data

        # this suppresses a divide by zero error which represents sites
        # with no observations of alleles 0 or 1. This is OK since such
        # sites will already be filtered as either non-biallelic or as
        # below mincov (all missing for these subsamples), or as
        # invariant (again due to all missing for these subsamples).
        with np.errstate(divide='ignore', invalid='ignore'):
            if isinstance(self.maf, int):
                freqs = called_1
            else:
                freqs = called_1 / (called_0 + called_1)
                freqs[freqs > 0.5] = 1 - (freqs[freqs > 0.5])
            masks[:, 5] = freqs < self.maf

        # set 9 for missing values in diploid genotype array
        diplos[snps == 78] = 9
        return masks, diplos

    ################################################################
    ## Subsampling/linkage functions
    ################################################################
    def subsample_snps(self, random_seed:Optional[int]=None, log_level: str="INFO") -> np.ndarray:
        """Return an array with 1 SNP sampled per locus/linkage-block.

        Uses numba jit function for speed, and snpsmap array to find
        linkage information.
        """
        rng = np.random.default_rng(random_seed)
        subarr = self.snps[:, jsubsample_snps(self.snpsmap, rng.integers(2**31))]
        logger.log(log_level, f"subsampled {subarr.shape[1]} unlinked SNPs.")
        return subarr

    def subsample_genos(self, random_seed:Optional[int]=None, log_level: str="INFO") -> np.ndarray:
        """Return an array with 1 SNP geno sampled per locus/linkage-block.

        Uses numba jit function for speed, and snpsmap array to find
        linkage information.
        """
        rng = np.random.default_rng(random_seed)
        subarr = self.genos[:, jsubsample_snps(self.snpsmap, rng.integers(2**31))]
        logger.log(log_level, f"subsampled {subarr.shape[1]} unlinked SNPs.")
        return subarr

    def subsample_loci(
        self,
        random_seed: Optional[int]=None,
        return_sites: bool=False,
        # invariant_loci: Optional[int]=None,
        log_level: str="INFO",
        ) -> np.ndarray:
        """Return an array of snps by re-sampling loci w/ replacement.

        Calls jitted functions to subsample loci/linkage-blocks with
        replacement to the same number as the original assembly.
        This does not subsample unlinked SNPs per locus, but instead
        re-samples linked SNPs for use in bootstrapping. Read below
        to be sure this is doing what you want it to do.

        Parameters
        ----------
        random_seed: int
            Random number generator seed.
        return_sites: bool
            If True sites np.uint8 values are returned, else diploid
            genotypes (0,1,2,9) are returned.
        invariant_loci: Optional[int]
            An optional additional number of invariant loci to resample
            from. When one of these is randomly sampled no SNPs are
            added to the array. This emulates sampling loci from a
            genome more closely. The number of invariant loci can be
            estimated from running locus_extracter instead of
            snps_extracter...?
        log_level: str
            The logging level of this function can be modified to
            effects its logged message. This is useful if it will be
            run many times repeatedly, set to "WARNING" to suppress.

        Note
        ----
        THIS RESAMPLES ONLY LOCI IN THE SNPS HDF5 FILE, MEANING THAT IT
        DOES NOT RESAMPLE INVARIANT LOCI. AND, IN FACT IT RESAMPLES
        LOCI AFTER FILTERING HAS REDUCED THE LOCI TO ONLY THOSE THAT
        INCLUDE THE SELECTED SAMPLES.
            So, even though a full data set (loci file) may include
        3000 loci, perhaps only 2000 of those will include a particular
        set of four samples, and only 1000 of those might be variable.
        This func will return a random resampling with replacment of
        those 1000 loci.
        """
        rng = np.random.default_rng(random_seed)
        nloci, lidxs = jsubsample_loci(self.snpsmap, rng.integers(2**31))
        if return_sites:
            subarr = self.snps[:, lidxs]
        else:
            subarr = self.genos[:, lidxs]
        logger.log(
            log_level,
            f"subsampled {subarr.shape[1]} SNPs from {nloci} variable "
            "loci w/ replacement.",
        )
        return subarr

    def get_population_geno_counts(
        self,
        subsample:bool=False,
        random_seed: Optional[int]=None,
        log_level="INFO",
        ):
        """Return dataframe with genos in treemix format.

        A   B   C   D
        0,2 2,0 2,0 0,2
        0,2 1,1 3,0 0,3
        0,2 2,0 3,0 0,2
        ...

        Parameters
        ----------
        subsample: bool
            If True then 1 SNP is randomly sampled per linkage block.
        random_seed: Optional[int]
            Seed used for random subsampling of unlinked SNPs.
        """
        # check for required imap groupings
        if not self.imap:
            imap = {i: [i] for i in self.names}
        else:
            imap = self.imap.copy()

        # write the headers
        data = pd.DataFrame(columns=sorted(imap))

        # get data
        if subsample:
            genos = self.subsample_genos(random_seed=random_seed, log_level=log_level)
        else:
            genos = self.genos

        # create 0,5 pairs for ancestral derived counts
        for pop in data.columns:
            samp = [self.names.index(i) for i in imap[pop]]
            ances = np.sum(genos[samp, :] == 0, axis=0) * 2
            deriv = np.sum(genos[samp, :] == 2, axis=0) * 2
            heter = np.sum(genos[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter
            data.loc[:, pop] = [
                "{},{}".format(i, j) for i, j in zip(ances, deriv)
            ]
        return data

    def get_population_geno_frequency(
        self,
        subsample: bool=False,
        random_seed=None,
        log_level="INFO",
        # imap: Optional[Dict]=None
        ):
        """Return a dataframe with genotype frequencies as in construct format.

        The population names will the names in the .imap dictionary
        if one is present, otherwise every individual will be treated
        as a separate population.

             0    1    2   ...
        A    0.0  0.5  1.0
        B    1.0  0.5  0.0
        C    0.0  0.0  0.0

        You can optionally impute missing data before running this
        command, or not -- missing values will be filled with NaN,
        which is fine since construct will try to infer their values.

        Parameters
        ----------
        subsample: bool
            If True then max of 1 SNP per linkage group is sampled.
        random_seed: int
            Random number generator seed.
        log_level: str
            A logging level name ("DEBUG", "INFO", "WARNING", ...).
        """
        # imap: Dict[str, List[str]]
        #     A dictionary mapping a new population name to a list of
        #     existing sample names. The allele frequencies will be of
        #     all samples that have data at each SNP in each population.

        # check for required imap groupings
        if not self.imap:
            imap = {i: [i] for i in self.names}
        else:
            imap = self.imap.copy()

        # write the headers
        data = pd.DataFrame(columns=sorted(imap))

        # subsample dataset to use.
        if subsample:
            genos = self.subsample_genos(random_seed=random_seed, log_level=log_level)
        else:
            genos = self.genos

        # create 0,5 pairs for ancestral derived counts
        for pop in data.columns:
            samp = [self.names.index(i) for i in imap[pop]]
            ances = np.sum(genos[samp, :] == 0, axis=0) * 2
            deriv = np.sum(genos[samp, :] == 2, axis=0) * 2
            heter = np.sum(genos[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter
            # allow divide by zero to be set as NaN
            with np.errstate(divide='ignore', invalid='ignore'):
                data.loc[:, pop] = deriv / (deriv + ances)
        return data.T



if __name__ == "__main__":

    import toytree
    import ipcoal
    import ipyrad.analysis as ipa
    import ipyrad as ip
    ipa.set_log_level("DEBUG")
    ip.set_log_level("DEBUG")

    tree = toytree.rtree.unittree(10, 1e6)
    model = ipcoal.Model(tree, Ne=1e5, nsamples=6)
    model.sim_loci(50, 100)
    model.apply_missing_mask(0.5)
    model.write_snps_to_hdf5(name="test", outdir="/tmp")

    # write popfile and load back as an imap
    model.write_popfile(name='test', outdir="/tmp", diploid=True)
    IMAP = ipa.popfile_to_imap("/tmp/test.popfile.tsv")

    # model.write_snps_to_hdf5(name="test", outdir="/tmp", diploid=True)
    tool = ipa.snps_extracter("/tmp/test.snps.hdf5",
        imap=IMAP,
        minmap={i:1 for i in IMAP},
        minmaf=0.1,
    )
    tool.run(cores=2)
    # print(tool.subsample_snps().view())

    # ipa.snps_imputer(tool.genos, tool.names, tool.imap, inplace=True).run()
    # print(tool.subsample_genos())
