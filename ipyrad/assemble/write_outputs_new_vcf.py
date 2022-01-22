#!/usr/bin/env python

"""Writes data from SNPs HDF5 to VCF format.

"""

from typing import TypeVar, Iterator
from pathlib import Path
import time

from loguru import logger
import h5py
import numpy as np
import pandas as pd

import ipyrad
from ipyrad.assemble.utils import BTS, GETCONS, DCONS, chroms2ints


Assembly = TypeVar("Assembly")


class BuildVcf:
    def __init__(self, data: Assembly):
        self.data = data
        self.revdict = chroms2ints(self.data, True)

    def _iter_snps_data(self) -> Iterator:
        """Yield chunks of data from SNPS HDF5."""
        with h5py.File(self.data.outfiles.)

    def _write_to_vcf(self) -> None:
        """Writes chunks to VCF file."""





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
