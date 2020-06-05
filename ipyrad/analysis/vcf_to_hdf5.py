#!/usr/bin/env python

"convert VCF to database format for SNP analyses"

# py2/3 compat
from __future__ import print_function, division
from builtins import range

import os
import tempfile

import h5py
import numpy as np
import pandas as pd
from numba import njit

from .utils import ProgressBar
from ..assemble.utils import TRANSFULL, GETCONS, IPyradError


class VCFtoHDF5(object):
    """
    Creates a temporary snps.hdf5 file conversion of the VCF file.
    For ipyrad assembled RAD seq data this will use RAD loci as the
    grouping of SNPs into linkage blocks for subsampling. If VCF is from
    some other assembly (e.g., WGS) then you can use the ld_block_size arg to
    group SNPs into linkage blocks for subsampling analyses.
    """
    def __init__(
        self,
        data,
        name="test",
        workdir="./analysis-vcf2hdf5",
        ld_block_size=None,
        quiet=False,
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
        self.database = ""
        self.quiet = quiet

        # check for data file
        self.database = os.path.join(self.workdir, self.name + ".snps.hdf5")
        assert os.path.exists(self.data), "file {} not found".format(self.data)
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # vcf format info
        self.source = ""
        self.reference = ""


    def run(self):
        """
        Parse and convert data to HDF5 file format.
        """
        # print message
        self._print("Indexing VCF to HDF5 database file")

        # get sample names, count header lines and nsnps
        self.get_meta()

        # init the database to fill
        self.init_database()

        # fill snps matrix
        self.build_chunked_matrix()

        # report on new database
        with h5py.File(self.database, 'r') as io5:
            self.nscaffolds = io5["snpsmap"][-1, 0]
            # self.nlinkagegroups = io5["snpsmap"][-1, 3]

        self._print(
            "HDF5: {} SNPs; {} linkage group"
            .format(
                self.nsnps,
                self.nscaffolds,
            )
        )
        self._print(
            "SNP database written to {}"
            .format(self.database)
        )


    def _print(self, msg):
        if not self.quiet:
            print(msg)


    def get_meta(self):
        """
        Skip and count ## lines, and get names from first # line in VCF.
        """
        # store a list of chrom names
        chroms = set()

        if self.data.endswith(".gz"):
            import gzip
            infile = gzip.open(self.data)
        else:
            infile = open(self.data)

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
                self.nsnps += 1
                chroms.add(data[0])

        # close file handle
        infile.close()

        # convert chroms into a factorized list
        self._print("VCF: {} SNPs; {} scaffolds".format(self.nsnps, len(chroms)))        


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
        """
        Fill HDF5 database with VCF data in chunks at a time.
        """

        # chunk retriever
        self.df = pd.read_csv(
            self.data, 
            sep="\t", 
            skiprows=self.hlines, 
            chunksize=int(1e5),
            index_col=False,  # prevent from interpreting int CHROM as index
        )

        # open h5
        with h5py.File(self.database, 'a') as io5:
            prog = ProgressBar(self.nsnps, 0, "converting VCF to HDF5")
            prog.finished = 0
            prog.update()

            # iterate over chunks of the file
            xx = 0
            lastchrom = "NULL"
            e0 = 0  # 1-indexed new-locus index, will advance in get_snps/lastchrom
            e1 = 0  # 0-indexed snps-per-loc index
            e2 = 0  # 0-indexed snps-per-loc position
            e3 = 0  # 1-indexed original-locus index, TODO, advancer
            e4 = 0  # 0-indexed global snps counter
            for chunkdf in self.df:

                # get sub arrays
                genos, snps = chunk_to_arrs(chunkdf, self.nsamples)

                # get sub snpsmap
                snpsmap, lastchrom = self.get_snpsmap(
                    chunkdf, lastchrom=lastchrom, e0=e0, e1=e1, e2=e2, e4=e4)

                # store sub arrays
                e0 = snpsmap[-1, 0].astype(int)
                e1 = snpsmap[-1, 1].astype(int) + 1
                e2 = snpsmap[-1, 2].astype(int) + 1
                e4 = snpsmap[-1, 4].astype(int) + 1

                # write to HDF5
                io5['snps'][:, xx:xx + chunkdf.shape[0]] = snps.T
                io5['genos'][xx:xx + chunkdf.shape[0], :] = genos
                io5['snpsmap'][xx:xx + chunkdf.shape[0], :] = snpsmap
                xx += chunkdf.shape[0]

                # print progress
                prog.finished = xx
                prog.update()

            # return with last chunk
            self.df = chunkdf

            # close h5 handle
            self._print("")


    def build_matrix(self):
        """
        Fill database with VCF data
        """

        # load vcf as dataframe (chunking makes this an iterator)
        self.df = pd.read_csv(self.data, sep="\t", skiprows=self.hlines)

        # get ref and alt alleles as a string series
        refalt = (self.df.REF + self.df.ALT).apply(str.replace, args=(",", ""))

        # genos array to fill from geno indices and refalt
        genos = np.zeros((self.nsnps, self.nsamples, 2), dtype="u1")
        snps = np.zeros((self.nsnps, self.nsamples), dtype="S1")
        snpsmap = np.zeros((self.nsnps, 5), dtype="u4")

        # iterate over samples indices
        for sidx in range(self.nsamples):

            # store geno calls for this sample (+9 goes to geno columns in vcf)
            glist = self.df.iloc[:, sidx + 9].apply(get_genos).apply(sorted)
            genos[:, sidx, :] = pd.concat([
                glist.apply(lambda x: x[0]),
                glist.apply(lambda x: x[1]),
            ], axis=1)

            # iterate over geno indices to get alleles
            for gidx in range(genos.shape[0]):

                # this sample's geno
                sgeno = genos[gidx, sidx]
                genos[gidx, sidx] = sgeno

                # bail if geno is missing (TODO: or indel).
                if sgeno[0] == 9:
                    call = b"N"
                    snps[gidx, sidx] = call
                    continue

                # convert to geno alleles
                call0 = refalt[gidx][sgeno[0]]
                call1 = refalt[gidx][sgeno[1]]

                # get ambiguity code
                if call0 != call1:
                    call = TRANSFULL[(call0, call1)].encode()
                else:
                    call = call0

                # store call
                snps[gidx, sidx] = call

        # convert snps to uint8
        snps = snps.view("u1").T

        # snpsmap: if ipyrad denovo it's easy, and they should just use hdf5.
        if ("ipyrad" in self.source) and ("pseudo-ref" in self.reference):
            snpsmap[:, 0] = self.df["#CHROM"].factorize()[0].astype("u4") + 1
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in self.df.groupby("#CHROM")])
            snpsmap[:, 2] = self.df.POS
            snpsmap[:, 3] = 0
            snpsmap[:, 4] = range(snpsmap.shape[0])
            self.df["BLOCK"] = self.df["#CHROM"]

        # snpsmap: if ipyrad ref VCF the per-RAD loc info is available too
        elif "ipyrad" in self.source:

            # skip to generic vcf method if ld_block_size is set:
            if not self.ld_block_size:
                snpsmap[:, 0] = self.df["#CHROM"].factorize()[0]
                snpsmap[:, 1] = np.concatenate(
                    [range(i[1].shape[0]) for i in self.df.groupby("#CHROM")])
                snpsmap[:, 2] = self.df.POS
                snpsmap[:, 3] = 0
                snpsmap[:, 4] = range(snpsmap.shape[0])
                self.df["BLOCK"] = self.df["#CHROM"]

        # snpsmap: for other program's VCF's we need ldsize arg to chunk.
        else:
            if not self.ld_block_size:
                raise IPyradError(
                    "You must enter an ld_block_size estimate for this VCF.")

        # (TODO: numpy/numba this loop instead of pand)
        # cut it up by block size (unless it's denovo, then skip.)
        if (self.ld_block_size) and ("pseudo-ref" not in self.reference):

            # block and df index counters
            bidx = 0
            dfidx = 0

            # iterate over existing scaffolds
            for _, scaff in self.df.groupby("#CHROM"):

                # start and end of this scaffold
                gpos = scaff.iloc[0, 1]
                end = scaff.POS.max()

                # iterate to break scaffold into linkage blocks
                while 1:

                    # grab a block
                    mask = (scaff.POS >= gpos) & (scaff.POS < gpos + self.ld_block_size)
                    block = scaff[mask]

                    # check for data and sample a SNP
                    if block.size:

                        # enter new block into dataframe
                        self.df.loc[dfidx:dfidx + block.shape[0], "BLOCK"] = bidx
                        dfidx += block.shape[0]
                        bidx += 1

                    # advance counter
                    gpos += self.ld_block_size

                    # break on end of scaff
                    if gpos > end:
                        break

            # store it (CHROMS/BLOCKS are stored 1-indexed !!!!!!)
            snpsmap[:, 0] = self.df.BLOCK + 1
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in self.df.groupby("BLOCK")])
            snpsmap[:, 2] = self.df.POS
            snpsmap[:, 3] = self.df["#CHROM"].factorize()[0] + 1
            snpsmap[:, 4] = range(self.df.shape[0])

        # store to database
        with h5py.File(self.database, 'a') as io5:
            io5['snps'][:] = snps
            io5['genos'][:] = genos
            io5['snpsmap'][:] = snpsmap
        del snps


    def get_snpsmap(self, chunkdf, lastchrom, e0, e1, e2, e4):

        # convert snps back to "S1" view to enter data...
        nsnps = chunkdf.shape[0]
        snpsmap = np.zeros((nsnps, 5), dtype=np.uint32)

        # check whether locus is same as end of last chunk
        currchrom = chunkdf.iloc[0, 0]
        if currchrom != lastchrom:
            e0 += 1
            e1 = 0
            e2 = 0

        # snpsmap: if ipyrad denovo it's easy, and they should just use hdf5.
        if ("ipyrad" in self.source) and ("pseudo-ref" in self.reference):

            # print warning that we're not ussng ld_block_size
            if self.ld_block_size:
                self._print(
                    "\nThis appears to be a denovo assembly, "
                    "ld_block_size arg is being ignored.")

            # get locus index
            snpsmap[:, 0] = (
                chunkdf["#CHROM"].factorize()[0].astype(np.uint32) + e0)

            # get snp index counter possibly continuing from last chunk
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in chunkdf.groupby("#CHROM")])
            snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1

            # get snp pos counter possibly continuing from last chunk
            snpsmap[:, 2] = chunkdf.POS + e2
            snpsmap[:, 3] = snpsmap[:, 0]

            # get total snp counter, always continuing from any previous chunk
            snpsmap[:, 4] = range(e4, snpsmap.shape[0] + e4)


        # snpsmap: if ipyrad ref VCF the per-RAD loc info is available too
        elif "ipyrad" in self.source:

            # skip to generic vcf method if ld_block_size is set:
            if not self.ld_block_size:
                snpsmap[:, 0] = (
                    chunkdf["#CHROM"].factorize()[0].astype(np.uint32) + e0)

                # add ldx counter from last chunk
                snpsmap[:, 1] = np.concatenate(
                    [range(i[1].shape[0]) for i in chunkdf.groupby("#CHROM")])
                snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1
                snpsmap[:, 2] = chunkdf.POS + e2
                snpsmap[:, 3] = snpsmap[:, 0]
                snpsmap[:, 4] = range(e4, snpsmap.shape[0] + e4)

        # snpsmap: for other program's VCF's we need ldsize arg to chunk.
        else:
            if not self.ld_block_size:
                raise IPyradError(
                    "You must enter an ld_block_size estimate for this VCF.")

        # cut it up by block size (unless it's denovo, then skip.)
        if (self.ld_block_size) and ("pseudo-ref" not in self.reference):

            # create a BLOCK column to keep track of original chroms
            chunkdf["BLOCK"] = 0

            # block and df index counters
            original_e0 = e0
            dfidx = e4

            # iterate over existing scaffolds (e.g., could be one big chrom)
            for _, scaff in chunkdf.groupby("#CHROM"):

                # current start and end POS of this scaffold before breaking
                gpos = e2
                end = scaff.POS.max()

                # iterate to break scaffold into linkage blocks
                while 1:

                    # grab a block 
                    mask = (scaff.POS >= gpos) & (scaff.POS < gpos + self.ld_block_size)
                    block = scaff[mask]

                    # check for data and sample a SNP
                    if block.size:

                        # enter new block into dataframe
                        chunkdf.loc[dfidx:dfidx + block.shape[0], "BLOCK"] = e0
                        dfidx += block.shape[0]
                        e0 += 1

                    # advance counter
                    gpos += self.ld_block_size

                    # break on end of scaff
                    if gpos > end:
                        break

            # store it (CHROMS/BLOCKS are stored 1-indexed !!!!!!)
            snpsmap[:, 0] = chunkdf.BLOCK
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in chunkdf.groupby("BLOCK")])
            # add ldx counter from last chunk
            snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1
            snpsmap[:, 2] = chunkdf.POS
            snpsmap[:, 3] = chunkdf["#CHROM"].factorize()[0] + original_e0
            snpsmap[:, 4] = range(e4, chunkdf.shape[0] + e4)
        return snpsmap, currchrom


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
    In development...
    Read in chunk of VCF and convert to numpy arrays
    """
    # nsnps in this chunk
    nsnps = chunkdf.shape[0]

    # chrom, pos as factor, integers
    arrpos = chunkdf.iloc[:, [0, 1]].values
    arrpos[:, 0] = pd.factorize(arrpos[:, 0])[0]
    arrpos = arrpos.astype(np.int64)

    # base calls as int8 (0/1/2/3/9)
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
