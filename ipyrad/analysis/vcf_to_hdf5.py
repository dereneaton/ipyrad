#!/usr/bin/env python

"""Convert VCF to HDF5 database format for SNP analyses.

This stores 'snps', 'snpsmap', and 'genos' arrays in the HDF5.

snpsmap columns:
    0: linkage block (ints as factors, so 0 vs 1 indexing doesn't matter)
    1: index on linkage block (0-indexed).
    2: position on assembled scaff alignment (?-indexed) (NOT USED CURRENTLY?)
    3: original chrom int (1-indexed)
    4. genomic position on original chrom (1-indexed); OR global snp counter.
"""

import os
import gzip
import tempfile

from loguru import logger
import h5py
import numpy as np
import pandas as pd
from numba import njit

from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import GETCONS, IPyradError


class VCFtoHDF5:
    """Create a snps.hdf5 file conversion from a VCF file.

    For ipyrad assembled RAD seq data this will use RAD loci as the
    grouping of SNPs into linkage blocks for subsampling. If VCF is 
    from some other assembly (e.g., WGS) then you can use the 
    ld_block_size arg to group SNPs into linkage blocks for subsampling
    analyses.

    Parameters
    ----------
    data: str
        Path to input VCF file.
    name: str
        Prefix name for output files.
    workdir: str
        Directory name for output files. Will be created if absent.
    ld_block_size: int
        By default each scaffold will be treated as an unlinked locus.
        To assign SNPs on the same locus to be on different "loci" you can
        set a block size. For example, 50000 will assign SNPs to 50Kb windows
        and then analyses will use this information to subsample SNPs.
    scaffolds: array-like
        A list of the VCF scaffold (CHROM) names as strings that you want to 
        include in the HDF5 database. Default is None which means use all 
        scaffolds. If you are not sure of the CHROM names then take a look
        at the VCF file you are trying to convert.
    chunksize: int
        The size of VCF chunks to read in at one time. This only affects 
        memory usage and speed. Default value is 5000. Use larger values 
        if the number of SNPs is super large to increase speed greatly.
    log_level: str
        DEBUG is most verbose, INFO is less, and WARNING least.
    """
    def __init__(
        self,
        data,
        name="test",
        workdir="./analysis-vcf2hdf5",
        ld_block_size=None,
        scaffolds=None,
        chunksize=2500,
        ):

        # attrs
        self.data = data
        self.name = (name if name else "test")
        self.workdir = (workdir if workdir else tempfile.gettempdir())
        self.names = []
        self.nsamples = 0
        self.nsnps = 0
        self.hlines = 0
        self.ld_block_size = ld_block_size
        self.chunksize = chunksize
        self.database = ""
        self.scaffolds = ([] if not scaffolds else [str(i) for i in scaffolds])

        # check for data file
        self.database = os.path.join(self.workdir, self.name + ".snps.hdf5")
        assert os.path.exists(self.data), "file {} not found".format(self.data)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # vcf format info
        self.source = ""
        self.reference = ""
        self.nscaffolds = 0


    def run(self, force=False):
        """Parse and convert data to HDF5 file format."""
        logger.info("Indexing VCF to HDF5 database file")
        if os.path.exists(self.database) and not force:
            logger.error("hdf5 file exists. Use `force=True` to overwrite.")
            return

        # get sample names, count header lines and nsnps
        self.get_meta()

        # init the database to fill
        self.init_database()

        # fill snps matrix
        self.build_chunked_matrix()

        # report on new database
        with h5py.File(self.database, 'r') as io5:
            self.nscaffolds = io5["snpsmap"][-1, 0]
        logger.info(f"HDF5: {self.nsnps} SNPs; {self.nscaffolds} linkage groups")
        logger.info(f"SNP database written to {self.database}")


    def get_meta(self):
        """Skip and count ## lines, and get names from first # line in VCF."""
        # store a list of chrom names
        chroms = set()

        if self.data.endswith(".gz"):
            infile = gzip.open(self.data, 'r')
        else:
            infile = open(self.data, 'r')

        # get data header line
        for dat in infile:

            # split on space
            try:
                data = dat.decode().split()
            except AttributeError:
                data = dat.split()

            # parse meta data lines
            if data[0][0] == "#":
                if data[0][1] != "#":
                    # parse names from header line in their order
                    self.names = data[9:]
                    self.nsamples = len(self.names)
                else:
                    # store header line count and look for format str
                    self.hlines += 1
                    if "source=" in data[0].lower():
                        self.source = data[0].split("source=")[-1].strip()
                    if "reference=" in data[0].lower():
                        self.reference = data[0].split("reference=")[-1].strip()

            # meta snps data
            else:
                # only store SNP and CHROM if CHROM is in self.scaffolds
                if not self.scaffolds:
                    self.nsnps += 1
                    chroms.add(data[0])
                else:
                    if str(data[0]) in self.scaffolds:
                        self.nsnps += 1
                        chroms.add(data[0])                       

        # close file handle
        infile.close()

        # convert chroms into a factorized list
        logger.info("VCF: {} SNPs; {} scaffolds".format(self.nsnps, len(chroms)))        


    def init_database(self):
        """
        # load vcf file as a pandas dataframe in chunks.
        """
        # init the database file
        with h5py.File(self.database, 'w') as io5:

            # core data sets (should SNPs be S1?)
            io5.create_dataset("genos", (self.nsnps, self.nsamples, 2), np.uint8)
            io5.create_dataset("snps", (self.nsamples, self.nsnps), np.uint8)
            io5.create_dataset("snpsmap", (self.nsnps, 5), np.uint32)
            io5["snps"].attrs["names"] = [i.encode() for i in self.names]
            io5["genos"].attrs["names"] = [i.encode() for i in self.names]
            io5["snpsmap"].attrs["columns"] = [
                b"locus", b"locidx", b"locpos", b"scaf", b"scafpos",
            ]


    def build_chunked_matrix(self):
        """Fill HDF5 database with VCF data in chunks at a time."""
        # chunk retriever
        vcfchunker = pd.read_csv(
            self.data, 
            sep="\t", 
            header=None,
            skiprows=self.hlines + 1, 
            chunksize=int(self.chunksize),
            index_col=False,  # prevent from interpreting int CHROM as index
            # dtype="string",
        )

        # counters and/or notes on indexing
        cdict = {
            'chromint': 0,    # re-code and store original chroms as ints
            'curscaff': -1,   # keep track of continuation over chunks
            'loc': 0,         # chrom name/idx
            'snpidx': 0,      # 0-indexed index of snp on chrom
            'pos': 0,         # 1-indexed position of snp on chrom
            'snps': 0,        # 0-indexed global SNP counter
        }

        # progress bar 
        prog = AssemblyProgressBar(
            jobs={}, message="converting VCF to HDF5", quiet=False)
        prog.update()

        # open h5
        with h5py.File(self.database, 'a') as io5:

            # iterate over chunks of the file
            ichunk = 0
            for cdf in vcfchunker:

                # filter to selected chroms
                if self.scaffolds:
                    mask = cdf.loc[:, 0].isin(self.scaffolds)
                    cdf = cdf.loc[mask]
                    if not cdf.size:
                        continue

                # get sub arrays as ints
                genos, snps = chunk_to_arrs(cdf, self.nsamples)

                # get sub snpsmap and advance counters (cols)
                snpsmap, cdict = self.get_snpsmap_external(cdf, cdict)

                # write to hdf5
                io5['snps'][:, ichunk:ichunk + cdf.shape[0]] = snps.T
                io5['genos'][ichunk:ichunk + cdf.shape[0], :] = genos
                io5['snpsmap'][ichunk:ichunk + cdf.shape[0], :] = snpsmap

                # print progress
                ichunk += cdf.shape[0]
                prog.finished = ichunk
                prog.update()

            # close h5 handle
            self._print("")



    def get_snpsmap_external(self, cdf, cdict):
        """
        Construct snpsmap from vcf chunk assuming it was assembled
        from an external software tool.
        """
        # print("\nNEWCHUNK-----------------------------------")
        # convert snps back to "S1" view to enter data...
        nsnps = cdf.shape[0]
        snpsmap = np.zeros((nsnps, 5), dtype=np.uint32)      

        # start BLOCK and POS as empty
        cdf["BLOCK"] = 0
        cdf["SNPIDX"] = 0
        cdf["CHROMINT"] = 0

        # iterate over scaffolds, which could load in one big chrom that
        # extends into the next chunk, or many small chroms.
        for group, scaff in cdf.groupby(0):

            # reset positional counter if this is a new scaff on same chunk
            if group != cdict['curscaff']:
                cdict['loc'] += 1
                cdict['pos'] = scaff[1].min()
                cdict['snpidx'] = 0
                cdict['curscaff'] = group
                cdict['chromint'] += 1

                # store original chrom as an int
                cdf.loc[scaff.index, "CHROMINT"] = cdict['chromint']

            # print('\n\ngroup', group)
            # loop and break scaff into linkage blocks if it is too large.
            # a scaff can finish within a chunk, but most chunks will end
            # with an unfinished scaff.
            while 1:

                # does the next block finish before end of scaff?
                ends = (cdict['pos'] + self.ld_block_size) <= scaff[1].max()

                # grab subset of scaff containing only the next ld block
                if ends:
                    mask0 = scaff[1] >= cdict["pos"]
                    mask1 = scaff[1] < (cdict["pos"] + self.ld_block_size)
                    bmask = mask0 & mask1                   
                    block = scaff.loc[bmask]

                    # if empty block then slide window further to find data
                    if not block.size:
                        cdict['pos'] += self.ld_block_size

                    # if block then store loc and snpidx and advance counts
                    else:
                        # store new ld block and snp index
                        cdf.loc[block.index, "BLOCK"] = cdict['loc']
                        cdf.loc[block.index, "SNPIDX"] = (
                            range(cdict['snpidx'], cdict['snpidx'] + block.shape[0])
                        )

                        # advance current position
                        cdict['loc'] += 1
                        cdict['snpidx'] = 0
                        cdict['pos'] += self.ld_block_size
                        cdict['snps'] += block.shape[0]
                        cdict['curscaff'] = block.iloc[-1, 0]

                # grab all remainder of this scaff
                else:
                    block = scaff[scaff[1] >= cdict['pos']]
                    cdf.loc[block.index, "BLOCK"] = cdict['loc']
                    cdf.loc[block.index, "SNPIDX"] = (
                        range(cdict['snpidx'], cdict['snpidx'] + block.shape[0])
                    )

                    # reset counters
                    cdict['pos'] = block[1].min()
                    cdict['snpidx'] = block.shape[0]
                    cdict['snps'] += block.shape[0]
                    cdict['curscaff'] = block.iloc[-1, 0]
                    break
                    
            # Note: upon leaving this loop we are likely in the middle
            # of a linkage block, so we will start the next one using 
            # the current cdict['loc'] and cdict['pos'].

        # store it (CHROMS/BLOCKS are stored 1-indexed !!!!!!)
        snpsmap[:, 0] = cdf.BLOCK     # recoded blocks (as ints)
        snpsmap[:, 1] = cdf.SNPIDX    # index of SNPs on each block
        snpsmap[:, 2] = cdf[1]        # position of SNPs on original chroms
        snpsmap[:, 3] = cdf.CHROMINT  # original chroms (as ints)
        snpsmap[:, 4] = range(cdict['snps'], cdict['snps'] + nsnps)
        # print(snpsmap)

        return snpsmap, cdict



    def get_snpsmap_ipyrad(self, chunkdf, col0, col1, col2, col4):
        """
        DEPRECATED: just always use the 'external' tool.

        A potential use for this tool would be to 're-code' your RAD loci
        by 'concatenating' RAD loci that are within N bp of each other 
        to occur on the same linkage block. But not that useful generally, 
        since 1 snp per RAD locus is probably usually good enough.
        """
        # convert snps back to "S1" view to enter data...
        nsnps = chunkdf.shape[0]
        snpsmap = np.zeros((nsnps, 5), dtype=np.uint32)

        # snpsmap: if ipyrad denovo it's easy, and they should just use hdf5.
        if "ipyrad" in self.source:

            if "pseudo-ref" in self.reference:

                # print warning that we're not ussng ld_block_size
                if self.ld_block_size:
                    self._print(
                        "\nThis appears to be a denovo assembly, "
                        "ld_block_size arg is being ignored.")

                # get locus index
                snpsmap[:, 0] = (
                    chunkdf[0].factorize()[0].astype(np.uint32) + col0)

                # get snp index counter possibly continuing from last chunk
                snpsmap[:, 1] = np.concatenate(
                    [range(i[1].shape[0]) for i in chunkdf.groupby(0)])
                snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += col1

                # get snp pos counter possibly continuing from last chunk
                snpsmap[:, 2] = chunkdf[1] + col2
                snpsmap[:, 3] = snpsmap[:, 0]

                # get total snp counter, always continuing from any previous chunk
                snpsmap[:, 4] = range(col4, snpsmap.shape[0] + col4)

            # snpsmap: if ipyrad ref VCF the per-RAD loc info is available too
            else:
                # skip to generic vcf method if ld_block_size is set:
                if not self.ld_block_size:
                    snpsmap[:, 0] = (
                        chunkdf[0].factorize()[0].astype(np.uint32) + col0)

                    # add ldx counter from last chunk
                    snpsmap[:, 1] = np.concatenate(
                        [range(i[1].shape[0]) for i in chunkdf.groupby(0)])
                    snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += col1
                    snpsmap[:, 2] = chunkdf[1] + col2
                    snpsmap[:, 3] = snpsmap[:, 0]
                    snpsmap[:, 4] = range(col4, snpsmap.shape[0] + col4)

        # snpsmap: for other program's VCF's we need ldsize arg to chunk.
        else:
            if not self.ld_block_size:
                raise IPyradError(
                    "You must enter an ld_block_size estimate for this VCF.")
        return snpsmap, col0, col1, col2, col4




def get_genos(gstr):
    """
    Extract geno from vcf string: e.g., (0/1:...). 
    First pulls the indices then uses refalt to pull genos.
    """
    gen = gstr.split(":")[0]
    try:
        return gen[0], gen[2]
    except IndexError:
        return 9, 9
    # return gstr[0], gstr[2]


def return_g(gstr, i):
    "returns the genotype str from vcf at one position (0/1) -> 0"
    gen = gstr.split(":")[0]
    try:
        return int(gen[i])
    except:
        return 9


# vectorized version of return g
v_return_g = np.vectorize(return_g, otypes=[np.int8])


def chunk_to_arrs(chunkdf, nsamples):
    """
    Read in chunk of VCF and convert to numpy arrays
    """
    # nsnps in this chunk
    nsnps = chunkdf.shape[0]

    # load chrom and pos columns as a (nsnps, 2) shape array
    # convert to ints: pos is already int, chrom converted to int factors
    arrpos = chunkdf.iloc[:, [0, 1]].values
    arrpos[:, 0] = pd.factorize(arrpos[:, 0])[0]
    arrpos = arrpos.astype(np.int64)

    # get geno calls as int8 (0/1/2/3/9)
    ref = chunkdf.iloc[:, 3].astype(bytes).view(np.int8).values
    alts = chunkdf.iloc[:, 4].astype(bytes)
    sas = np.char.replace(alts, b",", b"")
    alts1 = np.zeros(alts.size, dtype=np.int8)
    alts2 = np.zeros(alts.size, dtype=np.int8)
    alts3 = np.zeros(alts.size, dtype=np.int8)
    lens = np.array([len(i) for i in sas])
    alts1[lens == 1] = [i[0] for i in sas[lens == 1]]
    alts2[lens == 2] = [i[1] for i in sas[lens == 2]]
    alts3[lens == 3] = [i[2] for i in sas[lens == 3]]

    # genotypes as int8 
    g0 = v_return_g(chunkdf.iloc[:, 9:], 0)
    g1 = v_return_g(chunkdf.iloc[:, 9:], 2)
    genos = np.zeros((nsnps, nsamples, 2), dtype=np.int8)
    genos[:, :, 0] = g0
    genos[:, :, 1] = g1

    # numba func to fill
    snps = jfill_snps(nsnps, nsamples, ref, g0, g1, alts1, alts2, alts3)
    return genos, snps


@njit()
def jfill_snps(nsnps, nsamples, ref, g0, g1, alts1, alts2, alts3):
    """
    jit-compiled function to fill the snps array based on input
    genotype calls and the ref and alts columns.
    """

    # fill snps
    snps = np.zeros((nsnps, nsamples), dtype=np.int8)

    # fill snps in rows by indexing genos from ref,alt with g0,g1
    for ridx in range(snps.shape[0]):

        # get it
        tmpr = ref[ridx]
        tmp0 = g0[ridx]
        tmp1 = g1[ridx]

        # missing set to 78
        tmpsnps = snps[ridx]
        tmpsnps[tmp0 == 9] = 78
        snps[ridx] = tmpsnps

        # 0/0 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 + tmp1) == 0] = tmpr
        snps[ridx] = tmpsnps

        # 1/1 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 1) & (tmp1 == 1)] = alts1[ridx]
        snps[ridx] = tmpsnps   

        # 2/2 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 2) & (tmp1 == 2)] = alts2[ridx]
        snps[ridx] = tmpsnps   

        # 3/3 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 3) & (tmp1 == 3)] = alts3[ridx]
        snps[ridx] = tmpsnps 

    # fill ambiguity sites 
    ambs = np.where(g0 != g1)
    for idx in range(ambs[0].size):

        # row, col indices of the ambiguous site in the snps mat
        row = ambs[0][idx]
        col = ambs[1][idx]

        # get genos (0123) from the geno matrices
        a0 = g0[row, col]
        a1 = g1[row, col]        
        alls = sorted([a0, a1])

        # get the alleles (CATG) from the ref/alt matrices
        if alls[0] == 0:
            b0 = ref[row]
            if alls[1] == 1:
                b1 = alts1[row]
            elif alls[1] == 2:
                b1 = alts2[row]
            else:
                b1 = alts3[row]
        elif alls[0] == 1:
            b0 = alts1[row]
            if alls[1] == 2:
                b1 = alts2[row]
            else:
                b1 = alts3[row]
        elif alls[0] == 2:
            b0 = alts2[row]
            b1 = alts3[row]

        # convert allele tuples into an ambiguity byte
        fill = np.argmax((GETCONS[:, 2] == b0) & (GETCONS[:, 1] == b1))
        snps[row, col] = GETCONS[fill, 0]

    # return the three arrays
    return snps




if __name__ == "__main__":

    # quick test
    import h5py
    import ipcoal
    import toytree
    import ipyrad.analysis as ipa

    # get a tree topology
    tree = toytree.rtree.unittree(ntips=5, treeheight=1e5, seed=123)

    # setup simulation of loci
    mod = ipcoal.Model(tree=tree, Ne=1e5, nsamples=2)
    mod.sim_loci(nloci=100, nsites=1000)

    # write vcf
    mod.write_vcf(name="test-short", outdir="/tmp", diploid=True)

    # setup converter tool
    tool = ipa.vcf_to_hdf5(
        data="/tmp/test-short.vcf", 
        name="test-short", 
        workdir="/tmp",
        ld_block_size=500,
        chunksize=500,
    )

    # write converted file
    tool.run(force=True)

    # check h5 database
    with h5py.File("/tmp/test-short.snps.hdf5", 'r') as io5:
        print(io5['snpsmap'][:50])
