#!/usr/bin/env python

"""
Load SNPs array and SNPsMAP from ipyrad SNPs HDF5 database for use in 
created subsampled SNP data sets for ipyrad.analysis tools. 

    Examples:
        pca, structure, treemix, popgen, baba
"""

# standard lib
import h5py
import numpy as np
import pandas as pd
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import jsubsample_snps, jsubsample_loci

"""
.snps should be returned as dtype=... int8, int64?
.check imap for the user, e.g., no dups.
"""


class SNPsExtracter:
    """
    Extract .snps and .snpsmap from snps.hdf5 after filtering for 
    indels, bi-allelic, mincov, and minmap.

    Parameters
    ------------
    data (str)
        path to the .snps.hdf5 file produced by ipyrad (or created by 
        converting a vcf file to .snps.hdf5 using ipyrad.analysis).
    imap: dict[str, List[str]]
        dictionary mapping population names to a list of sample names.
    minmap: (dict)
        a dictionary mapping population names to a float or integer 
        representing the minimum samples required to be present for 
        a SNP to be retained in the filtered dataset.
    mincov: int or float
        the minimum samples required globally for a SNP to be retained 
        in the dataset.
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
    """
    def __init__(
        self, 
        data, 
        imap=None,
        minmap=None,
        mincov=0.0,
        minmaf=0.0,
        quiet=False,
        ):

        # store params
        self.data = data
        self.imap = imap
        self.imap = (imap if imap else {})
        self.minmap = (minmap if minmap else {i: 1 for i in self.imap})
        self.mincov = mincov
        self.quiet = quiet
        self.maf = minmaf

        # get names from imap, else will be filled/checked against db
        self.names = []
        for key in self.imap:
            if isinstance(self.imap[key], (list, tuple, np.ndarray)):
                self.names.extend(self.imap[key])
            elif isinstance(self.imap[key], str):
                self.names.append(self.imap[key])

        # boolean mask of sites that are filtered
        self.supermask = None
        # array of (ntaxa, nsnps) as int8, only samples in imap
        self.snps = None
        # array of (nsnps, 2) for subsampling SNPs with subsample_snps()
        self.snpsmap = None

        # check input data
        if self.data.endswith(".seqs.hdf5"):
            raise IPyradError("data should be .snps.hdf5 file not .seqs.hdf5.")
        if self.data.endswith(".vcf"):
            raise IPyradError("input should be hdf5, see the vcf_to_hdf5 tool.")
        self.parse_names_from_hdf5()

        # get number of snps
        with h5py.File(self.data, 'r') as io5:
            self.nsnps = io5['snpsmap'].shape[0]


 
    def _print(self, message):
        if not self.quiet:
            print(message)


    def parse_names_from_hdf5(self):
        """
        Parse ordered sample names from HDF5 and check agains imap names.
        """
        # load arrays from hdf5; could be more efficient.
        with h5py.File(self.data, 'r') as io5:

            # all names in database, maybe fewer in this analysis
            try:
                self.dbnames = [i.decode() for i in io5["snps"].attrs["names"]]
            except AttributeError:
                self.dbnames = [i for i in io5["snps"].attrs["names"]]

            # check for sample names not in the database file
            badnames = set(self.names).difference(self.dbnames)
            if badnames:
                raise IPyradError(
                    "Samples [{}] are not in data file: {}"
                    .format(badnames, self.data))

            # if no names then use database names, else order imap in dborder
            if not self.names:
                self.names = self.dbnames 
            else:
                self.names = [i for i in self.dbnames if i in self.names]

            # mask snp rows to only include selected samples
            self.sidxs = np.array(
                sorted([self.dbnames.index(i) for i in self.names]))


    def parse_genos_from_hdf5(self, return_as_characters=False):
        """
        Parse genotype calls from hdf5 snps file and store snpsmap
        for subsampling. This performs in a chunked way to avoid memory
        limits on super large datasets.

        snpsmap columns:
            0: 1-indexed scaff id
            1: 0-indexed snpidx
            2: 1-indexed scaffpos
            3: 0-indexed orig. scaffidx
            4: snpcounter]
        """

        # how many cols of SNPs to load in at once from snps, genos, snpsmap
        chunksize = 10000

        # store masks to be concatenated at end
        maskarrs = []

        # build filter for this chunksize
        for chunkidx in range(0, self.nsnps, chunksize):

            # send slice/chunk to get masks
            chunkslice = slice(chunkidx, chunkidx + chunksize)
            masks = self.get_masks_chunk(chunkslice)
            maskarrs.append(masks)

        # get stats from individual masks
        print(maskarrs)


        # concatenate masks into supermask array, and del list


        # report mask statistics


        # apply mask to create snpsmap


    def defunct(self):
        # snpsmap is used to subsample per locus 
        # [:, 1] is overwrit as a range over non-filtered sites below.
        snpsmap = io5["snpsmap"][chunkslice, :2]

        # .snps is an array of 0,1,2 or 9.
        # .snpsmap is ready to subsample .snps to 1-per-locus 


        # report pre-filter
        self._print("Samples: {}".format(len(self.names)))
        self._print("Sites before filtering: {}".format(snps.shape[1]))
        self._print("Filtered (indels): {}".format(mask0.sum()))
        self._print("Filtered (bi-allel): {}".format(mask1.sum()))
        self._print("Filtered (mincov): {}".format(mask2.sum()))
        self._print("Filtered (minmap): {}".format(mask3.sum()))
        self._print("Filtered (subsample invariant): {}".format(mask4.sum()))
        # report mask5 filters even if they overlap with other filters.
        self._print("Filtered (minor allele frequency): {}".format(mask5.sum()))


        # Calculate combined filters -------------------------------------
        # apply the filters
        summask = mask0 + mask1 + mask2 + mask3 + mask4 + mask5
        allmask = np.invert(summask)

        # store filtered int genotype calls (the output used by most tools)
        self.snps = diplo[self.sidxs, :][:, allmask]

        totalsnps = summask.sum()
        self._print("Filtered (combined): {}".format(totalsnps))

        # bail out if ALL snps were filtered
        if self.snps.size == 0:
            raise IPyradError("No SNPs passed filtering.")

        # report state
        self._print(
            "Sites after filtering: {}".format(self.snps.shape[1])
        )
        self._mvals = np.sum(self.snps == 9)
        self._msites = np.any(self.snps == 9, axis=0).sum()
        self._print(
            "Sites containing missing values: {} ({:.2f}%)"
            .format(
                self._msites, 
                100 * self._msites / self.snps.shape[1],
            )
        )
        self._print(
            "Missing values in SNP matrix: {} ({:.2f}%)"
            .format(
                self._mvals, 
                100 * self._mvals / self.snps.size,                
            )
        )

        # apply mask to snpsmap to map snps
        snpsmap = snpsmap[allmask, :]
        snpsmap[:, 1] = range(snpsmap.shape[0])
        self.snpsmap = snpsmap
        del snpsmap

        # final report
        self._print(
            "SNPs (total): {}\nSNPs (unlinked): {}".format(
                self.snps.shape[1], 
                np.unique(self.snpsmap[:, 0]).size
            )
        )

        # overwrite geno calls with str data.
        # store the filtered SNP calls (actual A,C,T or G) 
        if return_as_characters:
            self.snps = snps[self.sidxs, :][:, allmask].view("S1")
            # if the b'' bothers you then you can call .astype(str)




    def get_masks_chunk(self, chunkslice):
        """
        A single chunk iteration of parse_genos_from_hdf5()       
        """
        # masks are boolean shape=(5, chunksize)
        masks = np.zeros((5, chunkslice.stop - chunkslice.start), dtype=np.bool_)  

        # load arrays from hdf5; could be more efficient...
        with h5py.File(self.data, 'r') as io5:

            # snps is used to filter multi-allel and indel containing 
            snps = io5["snps"][:, chunkslice]

            # genos are the actual calls we want, after filtering
            genos = io5["genos"][chunkslice, :, :]

            # mask0 is True if an indel is present
            masks[0] = np.any(snps[self.sidxs, :] == 45, axis=0)


            # mask1 is True if a third or fourth allele is present.
            masks[1] = np.sum(genos[:, self.sidxs] == 2, axis=2).sum(axis=1).astype(bool)
            masks[1] += np.sum(genos[:, self.sidxs] == 3, axis=2).sum(axis=1).astype(bool)


            # mask missing calls from genotype array
            genomask = np.ma.array(data=genos, mask=(genos == 9))

            # mask2 is True if sample coverage is below mincov. This method
            # accomodates missing haploid calls (0/9) by counting alleles.
            nhaplos = (~genomask).sum(axis=2).sum(axis=1)
            if isinstance(self.mincov, int):
                masks[2] = nhaplos < (2 * self.mincov)
            elif isinstance(self.mincov, float):
                masks[2] = nhaplos < (2 * self.mincov * len(self.sidxs))
            else:
                raise IPyradError("mincov should be an int or float.")


            # mask3 is True if any pop sample cov in below imap/minmap -------
            if self.imap:
                for pop in self.imap:
                    
                    # get min cov for this pop
                    mincov = self.minmap[pop]

                    # get the indices for samples in this pop
                    pidxs = np.array(
                        sorted(self.dbnames.index(i) for i in self.imap[pop])
                    )

                    # select the samples from the geno array
                    subarr = genomask[:, pidxs, :]

                    # get number of haplotypes skipping masked missing
                    nhaplos = (~subarr).sum(axis=2).sum(axis=1)                   
                    if isinstance(mincov, int):
                        masks[3] += nhaplos < (2 * mincov)
                    elif isinstance(mincov, float):
                        masks[3] += nhaplos < (2 * mincov * len(pidxs))
                    else:
                        raise IPyradError("minmap dictionary malformed.")


            # mask4 is True if site is not a variant w/ current sampling------
            # get the common diploid call at each site (0, 1, 2) 
            diplo_common = (
                genomask[:, self.sidxs, :]
                .sum(axis=2)
                .mean(axis=1)
                .round()
                .astype(int)
                .data
            )
            # get diploid calls at each site
            diplos = genomask[:, self.sidxs, :].sum(axis=2).data.T
            # mask sites where all samples match common geno.
            masks[4] = np.all(diplo_common == diplos, axis=0)


            # mask5 is True if maf is below minmaf setting -------------------
            called_0 = (genomask == 0).sum(axis=2).sum(axis=1).data
            called_1 = (genomask == 0).sum(axis=2).sum(axis=1).data            
            if isinstance(self.maf, int):
                freqs = called_1
            else:
                freqs = called_1 / (called_0 + called_1)
                freqs[freqs > 0.5] = 1 - (freqs[freqs > 0.5])
            masks[5] = freqs < self.maf

        return masks




    def old_parse_genos_from_hdf5(self, return_as_characters):
        """
        A single chunk iteration of parse_genos_from_hdf5()       
        """
        # load arrays from hdf5; could be more efficient...
        with h5py.File(self.data, 'r') as io5:

            # snps is used to filter multi-allel and indel containing 
            snps = io5["snps"][:]

            # snpsmap is used to subsample per locus 
            # [:, 1] is overwrit as a range over non-filtered sites below.
            snpsmap = io5["snpsmap"][:, :2]

            # genos are the actual calls we want, after filtering
            genos = io5["genos"][:]  # .sum(axis=2).T

            # report pre-filter
            self._print("Samples: {}".format(len(self.names)))
            self._print("Sites before filtering: {}".format(snps.shape[1]))


            # filter all sites containing an indel in selected samples--------
            # mask0 is True if an indel is present
            mask0 = np.any(snps[self.sidxs, :] == 45, axis=0)
            self._print("Filtered (indels): {}".format(mask0.sum()))


            # filter all sites w/ multi-allelic in the selected samples-------
            # mask1 is True if a third allele is present.
            mask1 = np.sum(genos[:, self.sidxs] == 2, axis=2).sum(axis=1).astype(bool)
            mask1 += np.sum(genos[:, self.sidxs] == 3, axis=2).sum(axis=1).astype(bool)
            self._print("Filtered (bi-allel): {}".format(mask1.sum()))


            # filter based on mincov of subsamples (default = 0.0)------------
            # mask2 is True if sample coverage is below mincov.
            # convert genos (e.g., 1/1) to genos sums (e.g., 2)
            # TODO: if haploid missing calls are allowed (0|9) this is problem
            diplo = genos.sum(axis=2).T

            # convert any summed missing (18) to 9
            diplo[diplo == 18] = 9
            cov = np.sum(diplo[self.sidxs, :] != 9, axis=0)

            if isinstance(self.mincov, int):
                mask2 = cov < self.mincov

            elif isinstance(self.mincov, float):
                mask2 = cov < self.mincov * len(self.sidxs)

            else:
                raise IPyradError("mincov should be an int or float.")
            self._print("Filtered (mincov): {}".format(mask2.sum()))


            # mask3 is True if any pop sample cov in below imap/minmap -------
            # filter based on per-population min coverages
            mask3 = np.zeros(snps.shape[1], dtype=np.bool_)
            if self.imap:
                for key, val in self.imap.items():
                    mincov = self.minmap[key]
                    pidxs = np.array(sorted(self.dbnames.index(i) for i in val))
                    try:
                        subarr = diplo[pidxs, :]
                    except IndexError as err:
                        raise IPyradError(
                            "imap is empty: {} - {}".format(key, val)
                        ) from err
                    counts = np.sum(subarr != 9, axis=0)
                    if isinstance(mincov, float):
                        mask3 += (counts / subarr.shape[0]) < mincov
                    elif isinstance(mincov, int):
                        mask3 += counts < mincov
                    else:
                        raise IPyradError("minmap dictionary malformed.")
            self._print("Filtered (minmap): {}".format(mask3.sum()))


            # mask4 is True if site is not a variant w/ current sampling------
            # filter based on whether subsample still is variable at site
            # after masking missing data values, while also collapsing data
            # back into diploid genotype calls where haplo-missing is masked.
            marr = np.ma.array(
                data=diplo[self.sidxs, :], 
                mask=diplo[self.sidxs, :] == 9,
            )

            # round to the most common call (0, 1, 2)
            common = marr.mean(axis=0).round().astype(int)

            # mask if any non-mask data is not the common genotype and invert
            mask4 = np.invert(np.any(marr != common, axis=0).data)
            self._print("Filtered (subsample invariant): {}".format(mask4.sum()))

            
            # mask5 is True if maf is below minmaf setting -------------------
            # maf filter: applies to all gene copies (haploids)
            if isinstance(self.maf, int):
                freqs = marr.sum(axis=0)
            else:
                freqs = marr.sum(axis=0) / (2 * np.sum(marr.mask == False, axis=0))
                freqs[freqs > 0.5] = 1 - (freqs[freqs > 0.5])
            mask5 = freqs <= self.maf

            # only report filters unique to mask5
            # masks = mask0 + mask1 + mask2 + mask3 + mask4
            # drop = mask5[~masks].sum()
            # self._print("Filtered (minor allele frequency): {}".format(drop))

            # report mask5 filters even if they overlap with other filters.
            self._print("Filtered (minor allele frequency): {}".format(mask5.sum()))


            # Calculate combined filters -------------------------------------
            # apply the filters
            summask = mask0 + mask1 + mask2 + mask3 + mask4 + mask5
            allmask = np.invert(summask)

            # store filtered int genotype calls (the output used by most tools)
            self.snps = diplo[self.sidxs, :][:, allmask]

            totalsnps = summask.sum()
            self._print("Filtered (combined): {}".format(totalsnps))

            # bail out if ALL snps were filtered
            if self.snps.size == 0:
                raise IPyradError("No SNPs passed filtering.")

            # report state
            self._print(
                "Sites after filtering: {}".format(self.snps.shape[1])
            )
            self._mvals = np.sum(self.snps == 9)
            self._msites = np.any(self.snps == 9, axis=0).sum()
            self._print(
                "Sites containing missing values: {} ({:.2f}%)"
                .format(
                    self._msites, 
                    100 * self._msites / self.snps.shape[1],
                )
            )
            self._print(
                "Missing values in SNP matrix: {} ({:.2f}%)"
                .format(
                    self._mvals, 
                    100 * self._mvals / self.snps.size,                
                )
            )

            # apply mask to snpsmap to map snps
            snpsmap = snpsmap[allmask, :]
            snpsmap[:, 1] = range(snpsmap.shape[0])
            self.snpsmap = snpsmap
            del snpsmap

            # final report
            self._print(
                "SNPs (total): {}\nSNPs (unlinked): {}".format(
                    self.snps.shape[1], 
                    np.unique(self.snpsmap[:, 0]).size
                )
            )

            # overwrite geno calls with str data.
            # store the filtered SNP calls (actual A,C,T or G) 
            if return_as_characters:
                self.snps = snps[self.sidxs, :][:, allmask].view("S1")
                # if the b'' bothers you then you can call .astype(str)

        # .snps is an array of 0,1,2 or 9.
        # .snpsmap is ready to subsample .snps to 1-per-locus 




    def subsample_snps(self, random_seed=None, quiet=False):
        """
        Calls jitted functions to subsample 1 SNP per locus/linkage-block
        using snpsmap.
        """
        if not random_seed:
            random_seed = np.random.randint(0, 1e9)
        subarr = self.snps[:, jsubsample_snps(self.snpsmap, random_seed)]
        if not quiet:
            self._print("subsampled {} unlinked SNPs".format(subarr.shape[1]))
        return subarr


    def subsample_loci(self, random_seed=None, quiet=False):
        """
        Calls jitted functions to subsample loci/linkage-blocks with
        replacement to the same number as the original assembly. This does
        not subsample unlinked SNPs per locus, but instead re-samples linked
        SNPs for use in bootstrapping (e.g., BABA analysis).

        THIS RESAMPLES ONLY LOCI IN THE SNPS HDF5 FILE, MEANING THAT IT DOES
        NOT RESAMPLE INVARIANT LOCI. AND IN FACT IT RESAMPLES LOCI AFTER
        FILTERING HAS REDUCED THE LOCI TO ONLY THOSE THAT INCLUDE THE 
        SELECTED SAMPLES. 
        So a full data set (e.g., loci file) may include 
        3000 loci of which only 2000 include data for the a particular set
        of four samples, and only 1000 of those are variable. This will 
        return a random resampling with replacment of those 1000 loci.
        """
        if not random_seed:
            random_seed = np.random.randint(0, 1e9)
        nloci, lidxs = jsubsample_loci(self.snpsmap, random_seed)
        subarr = self.snps[:, lidxs]
        if not quiet:
            self._print(
                "subsampled {} SNPs from {} variable loci w/ replacement"
                .format(subarr.shape[1], nloci))
        return subarr


    def get_population_geno_counts(self, quiet=False):
        """
        Returns DF with genos in treemix format.
        A   B   C   D
        0,2 2,0 2,0 0,2
        0,2 1,1 3,0 0,3
        0,2 2,0 3,0 0,2
        ...
        """
        # check for required imap groupings
        assert self.imap, "imap dictionary is required to get counts."

        # write the headers
        df = pd.DataFrame(columns=sorted(self.imap))

        # create 0,5 pairs for ancestral derived counts
        for pop in df.columns:
            samp = [self.names.index(i) for i in self.imap[pop]]
            ances = np.sum(self.snps[samp, :] == 0, axis=0) * 2
            deriv = np.sum(self.snps[samp, :] == 2, axis=0) * 2
            heter = np.sum(self.snps[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter

            df.loc[:, pop] = [
                "{},{}".format(i, j) for i, j in zip(ances, deriv)
            ]
        return df



    def get_population_geno_frequency(self, subsample=False, individuals=False, seed=None, quiet=False):
        """
        Returns DF with geno freqs as in construct format.
             0    1    2   ...
        A    0.0  0.5  1.0
        B    1.0  0.5  0.0
        C    0.0  0.0  0.0

        You can optionally impute missing data before running this command
        or not and missing values will be filled with NaN, which is fine since
        construct will try to infer their values.

        Parameters:
        -----------
        subsample (bool):
            If True then max of 1 SNP per linkage group is sampled.
        quiet (bool):
            Print verbose output.
        """
        # check for required imap groupings
        assert self.imap, "imap dictionary is required to get counts."

        # write the headers
        df = pd.DataFrame(columns=sorted(self.imap))

        # create 0,5 pairs for ancestral derived counts
        for pop in df.columns:
            samp = [self.names.index(i) for i in self.imap[pop]]
            ances = np.sum(self.snps[samp, :] == 0, axis=0) * 2
            deriv = np.sum(self.snps[samp, :] == 2, axis=0) * 2
            heter = np.sum(self.snps[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter
            df.loc[:, pop] = deriv / (deriv + ances)

        # subsample unlinked SNPs
        if subsample:
            if not seed:
                seed = np.random.randint(int(1e9))
            df = df.iloc[jsubsample_snps(self.snpsmap, seed), :]

        return df.T






    # def subsample_loci_full(self, nloci, random_seed=None, quiet=False):
    #     """
    #     Calls jitted functions to subsample loci/linkage-blocks with
    #     replacement to the same number as the original assembly. This does
    #     not subsample unlinked SNPs per locus, but instead re-samples linked
    #     SNPs for use in bootstrapping (e.g., BABA analysis).
    #     THIS RESAMPLES FROM ALL LOCI IN THE DATA SET, WHETHER IT IS 
    #     VARIABLE OR NOT.
    #     So a full data set (e.g., loci file) may include 
    #     3000 loci of which only 2000 include data for the a particular set
    #     of four samples. This will return 2000 a random resampling with 
    #     replacment of those 2000 loci.
    #     """

    #     # set random seed 
    #     if not random_seed:
    #         random_seed = np.random.randint(0, 1e9)

    #     # expand snpsmap to include invariant loci that can be resampled.


    #     # subsample nloci with replacment
    #     nloci, lidxs = jsubsample_loci(self.snpsmap, random_seed)
    #     subarr = self.snps[:, lidxs]


    #     if not quiet:
    #         self._print(
    #             "subsampled {} SNPs from {} loci w/ replacement"
    #             .format(subarr.shape[1], nloci))           
    #     return subarr
