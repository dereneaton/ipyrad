#!/usr/bin/env python

"""
Utils for writing VCF with or without depths info.
"""

import os
import glob
import time

from loguru import logger
import h5py
import numpy as np
import pandas as pd

import ipyrad
from ipyrad.assemble.utils import BTS, GETCONS, DCONS, chroms2ints



class FillVCF:
    """
    Incorporate indels and trim amounts when grabbing depths from CATG arrays
    (depth arrays from step 5). Indels are only releveant to denovo data.
    """
    def __init__(self, data, nsnps, sample):

        # input locus bits
        self.locbits = glob.glob(os.path.join(data.tmpdir, "chunk*.loci"))
        self.locbits = sorted(
            self.locbits, key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))
        self.loclines = None

        # input arrays of indels arrays
        self.indbits = glob.glob(os.path.join(data.tmpdir, "chunk*.indels*"))
        if not self.indbits:
            self.indbits = [None] * len(self.locbits)

        # input trim arrays
        self.trimbits = glob.glob(os.path.join(data.tmpdir, "chunk*.npy"))
        self.trimbits = sorted(
            self.trimbits, key=lambda x: int(x.rsplit("-", 1)[-1][:-4]))

        # array to store vcfdepths for this taxon
        self.vcfd = np.zeros((nsnps, 4), dtype=np.uint32)

        # the sample for this comp
        self.sname = sample.name
        self.isref = bool(data.isref)

        # snpsmap has locations of SNPs on trimmed loci, e.g., 
        # no SNPs are on loc 1 and 2, first is on 3 at post-trim pos 11
        # [    3     0    11     1 41935]
        # [    4     0    57     1 56150]
        with h5py.File(data.snps_database, 'r') as io5:
            self.snpsmap = io5['snpsmap'][:, [0, 2]]   

        # TODO: scaffs should be ordered (right?) so no need to load it all!
        # All catgs for this sample (this could be done more mem efficient...)
        with h5py.File(sample.files.database, 'r') as io5:
            self.catgs = io5['catg'][:]
            self.maxlen = self.catgs.shape[1]

        # Sample-level counters
        self.locidx = 0
        self.snpidx = 0


    def run(self):
        "loops over chunked files streaming through all loci for this sample"
        for idx in range(len(self.locbits)):
            self.localidx = 0
            self.locfill(idx)


    def locfill(self, idx):
        "iterates over loci in chunkfile to get and enter catgs for snps"
        # load the arrays for this bit
        edges = np.load(self.trimbits[idx])
        inds = self.indbits[idx]
        if inds:
            inds = np.load(inds)

        # iterate over the chunk of trimmed loci
        self.loclines = iter(open(self.locbits[idx], 'r'))
        while 1:

            # yield increments locidx by 1
            try:
                self.yield_loc()
            except StopIteration:
                break

            # get snps for this locus (1-indexed locus idxs)
            self.locsnps = self.snpsmap[self.snpsmap[:, 0] == self.locidx]

            # get global trim for this locus (0-indexed edge arr)
            self.gtrim = edges[self.localidx - 1]

            # if SNPs and data for this sample enter catgs
            if (self.locsnps.size) and (self.sname in self.names):
                if self.isref:
                    self.ref_enter_catgs()
                else:
                    self.denovo_enter_catgs()
            else:
                # advance SNP counter even though this sample wasn't in SNP
                self.snpidx += self.locsnps.shape[0]


    def ref_enter_catgs(self):

        # map SNP position to pre-trim locus position
        nidx = self.names.index(self.sname)
        sidx = self.sidxs[nidx]
        tups = [[int(j) for j in i.split(":")] for i in sidx.split("-")]

        # SNP is in samples, so get and store catg data for locidx
        # [0] post-trim chrom:start-end of locus
        # [1:] how far ahead of start does this sample start
        # FOR DEBUGGING 
        # seq = seqs[nidx]
        # seqarr = np.array(list(seq))

        # enter each SNP 
        for snp in self.locsnps[:, 1]:
            # in case multiple consens were merged in step 6 of this sample
            for tup in tups:
                cidx, coffset = tup
                pos = snp + (self.gtrim - coffset)
                if (pos >= 0) & (pos < self.maxlen):
                    self.vcfd[self.snpidx] += self.catgs[cidx, pos]
            self.snpidx += 1


    def denovo_enter_catgs(self):
        """
        Grab catg depths for each SNP position -- needs to take into account
        trim from left end, and impution of indels.
        """
        nidx = self.names.index(self.sname)
        sidx = self.sidxs[nidx]
        tups = [[int(j) for j in i.split("-")] for i in sidx.split(":")]

        # SNP is in samples, so get and store catg data for locidx
        # [0] post-trim chrom:start-end of locus
        # [1:] how far ahead of start does this sample start
        # FOR DEBUGGING 
        seq = self.seqs[nidx]

        # enter each SNP 
        for snp in self.locsnps[:, 1]:
            # indels before this SNP
            ishift = seq[:snp].count("-")

            # in case multiple consens were merged in step 6 of this sample
            for tup in tups:
                cidx, coffset = tup
                # pos = snp + (self.gtrim - coffset) - ishift
                pos = snp + coffset - ishift                
                if (pos >= 0) & (pos < self.maxlen):
                    self.vcfd[self.snpidx] += self.catgs[cidx, pos]
            self.snpidx += 1


    def yield_loc(self):
        self.names = []
        self.seqs = []
        while 1:
            line = next(self.loclines)
            if "|\n" not in line:
                # skip if .loci chunk is empty
                try:
                    name, seq = line.split()
                    self.names.append(name)
                    self.seqs.append(seq)
                except ValueError:
                    continue
            else:
                self.locidx += 1
                self.localidx += 1                
                self.sidxs = [i for i in line.rsplit("|", 2)[1].split(',')]
                break




def build_vcf(data, chunksize=1000):
    """
    
    """
    # removed at init of Step function anyway.
    if os.path.exists(data.outfiles.vcf):
        os.remove(data.outfiles.vcf)

    # dictionary to translate locus numbers to chroms
    if data.isref:
        revdict = chroms2ints(data, True)

    # pull locus numbers and positions from snps database
    with h5py.File(data.snps_database, 'r') as io5:

        # iterate over chunks
        for chunk in range(0, io5['genos'].shape[0], chunksize):

            # if reference then psuedo ref is already ordered with REF first.
            pref = io5['pseudoref'][chunk:chunk + chunksize]
            snpmap = io5['snpsmap'][chunk:chunk + chunksize]

            # load array chunks
            if data.isref:
                genos = io5['genos'][chunk:chunk + chunksize, 1:, :]
                snames = data.snames[1:]

                # 1-indexed to 0-indexed (1/9/2019)
                chroms = [revdict[i - 1] for i in snpmap[:, 3]]
                ids = [
                    "loc{}_pos{}".format(i - 1, j) for (i, j) 
                    in snpmap[:, [0, 2]]
                ]

                # reference based positions: pos on scaffold: 4, yes. tested.
                pos = snpmap[:, 4]
                #offset = 1

            else:
                genos = io5['genos'][chunk:chunk + chunksize, :, :]
                snames = data.snames
                chroms = ["RAD_{}".format(i - 1) for i in snpmap[:, 0]]
                ids = [
                    "loc{}_pos{}".format(i - 1, j) for (i, j) 
                    in snpmap[:, [0, 2]]
                ]
                # denovo based positions: pos on locus. tested. works. right.
                # almost. POS is 1 indexed.
                pos = snpmap[:, 2] + 1
                # offset = 0

            # get alt genotype calls
            alts = [
                b",".join(i).decode().strip(",")
                for i in pref[:, 1:].view("S1") 
            ]

            # build df label cols
            df_pos = pd.DataFrame({
                '#CHROM': chroms,
                'POS': pos,            # 1-indexed
                'ID': ids,             # 0-indexed
                'REF': [i.decode() for i in pref[:, 0].view("S1")],
                'ALT': alts,
                'QUAL': [13] * genos.shape[0],
                'FILTER': ['PASS'] * genos.shape[0],
            })

            # get sample coverage at each site
            nsamps = (
                genos.shape[1] - np.any(genos == 9, axis=2).sum(axis=1)
            )

            # store sum of coverage at each site
            asums = []

            # build depth columns for each sample
            df_depth = pd.DataFrame({})
            for sname in snames:

                # build geno strings
                genostrs = [
                    "{}/{}".format(*k) for k in [
                        i for i in [
                            list(j) for j in genos[:, snames.index(sname)]
                        ]
                    ]
                ]

                # change 9's into missing
                genostrs = ["./." if i == "9/9" else i for i in genostrs]

                # genostrs = [
                # b"/".join(i).replace(b"9", b".").decode()
                # for i in genos[:, snames.index(sname)]
                # .astype(bytes)
                # ]

                # build depth and depthsum strings
                dpth = os.path.join(data.tmpdir, sname + ".depths.hdf5")
                with h5py.File(dpth, 'r') as s5:
                    dpt = s5['depths'][chunk:chunk + chunksize]
                    sums = [sum(i) for i in dpt]
                    strs = [
                        ",".join([str(k) for k in i.tolist()])
                        for i in dpt
                    ]

                    # save concat string to name
                    df_depth[sname] = [
                        "{}:{}:{}".format(i, j, k) for (i, j, k) in 
                        zip(genostrs, sums, strs)]

                    # add sums to global list
                    asums.append(np.array(sums))

            # make final columns
            nsums = sum(asums)
            colinfo = pd.Series(
                name="INFO",
                data=[
                    "NS={};DP={}".format(i, j) for (i, j) in zip(nsamps, nsums)
                ])
            colform = pd.Series(
                name="FORMAT",
                data=["GT:DP:CATG"] * genos.shape[0],
                )

            # concat and order columns correctly
            infocols = pd.concat([df_pos, colinfo, colform], axis=1)
            infocols = infocols[
                ["#CHROM", "POS", "ID", "REF", "ALT",
                 "QUAL", "FILTER", "INFO", "FORMAT"]]
            arr = pd.concat([infocols, df_depth], axis=1)

            # debugging                       
            #print(arr.head())
            ## PRINTING VCF TO FILE
            ## choose reference string
            if data.isref:
                reference = data.params.reference_sequence
            else:
                reference = "pseudo-reference (most common base at site)"

            header = VCFHEADER.format(
                date=time.strftime("%Y/%m/%d"),
                version=ipyrad.__version__,
                reference=os.path.basename(reference)
                ) 

            with open(data.outfiles.vcf, 'a') as out:
                if chunk == 0:
                    out.write(header)
                    arr.to_csv(out, sep='\t', index=False)
                else:
                    arr.to_csv(out, sep='\t', index=False, header=False)



VCFHEADER = """\
##fileformat=VCFv4.0
##fileDate={date}
##source=ipyrad_v.{version}
##reference={reference}
##phasing=unphased
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
"""
