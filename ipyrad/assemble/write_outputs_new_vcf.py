#!/usr/bin/env python

"""Writes data from SNPs HDF5 to VCF format.

"""

from typing import TypeVar, Iterator, Tuple, List
from pathlib import Path
import time

import h5py
import numpy as np
import pandas as pd

import ipyrad
from ipyrad.assemble.utils import chroms2ints


Assembly = TypeVar("Assembly")
CHUNKSIZE = 100 # 10_000


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


class BuildVcf:
    def __init__(self, data: Assembly):
        self.data = data

        # attributes to be filled.
        self.revdict = chroms2ints(self.data, keys_as_ints=True)
        """: A dict to convert chrom ints to chrom str names."""
        self.snames: List[str] = []
        """: A list of names in alphanumeric order, optionally including ref as sample."""

    def _iter_snps_data_chunk(self) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Yield chunks of data from SNPS HDF5."""
        with h5py.File(self.data.outfiles["snps_database"], 'r') as io5:
            self.snames = list(io5['snpsmap'].attrs["names"])
            for start in range(0, io5['snpsmap'].shape[0], CHUNKSIZE): # pylint: disable=no-member
                snpsmap = io5['snpsmap'][start:start + CHUNKSIZE, :]
                alts = io5['alts'][start:start + CHUNKSIZE, :]
                genos = io5['genos'][:, start:start + CHUNKSIZE, :]
                yield snpsmap, alts, genos

    def _iter_vcf_data_chunk(self) -> Iterator[pd.DataFrame]:
        """Writes chunks to VCF file."""
        for snpsmap, alts, genos in self._iter_snps_data_chunk():

            # get scaffold names, positions, and ID string.
            chrom_names = [self.revdict[i] for i in snpsmap[:, 3]]
            rad_ids = [
                f"loc_{i}_pos_{j}" for (i, j) in zip(snpsmap[:, 0], snpsmap[:, 2])
            ]

            # get ref and alternate calls.
            ref = list(alts[:, 0].tobytes().decode())
            alleles = []
            for idx in range(alts.shape[0]):
                calls = alts[idx, 1:]
                calls = calls[calls > 0]
                calls = ",".join(calls.tobytes().decode())
                alleles.append(calls)

            # build vcf dataframe for first 7 columns.
            vcfdf = pd.DataFrame({
                '#CHROM': chrom_names,                  # chrom str name
                'POS': snpsmap[:, 4],                   # 1-indexed position on scaff
                'ID': rad_ids,                          # RAD locus (x-positioned) and position 0-indexed.
                'REF': ref,                             # reference allele.
                'ALT': alleles,                         # other alleles in order.
                'QUAL': [13] * snpsmap.shape[0],        # arbitrarily high quality
                'FILTER': ['PASS'] * snpsmap.shape[0],  # no filter applied.
            })

            # get sample coverage at each SNP site
            nsamples = genos.shape[0] - np.any(genos == 9, axis=2).sum(axis=0)
            colinfo = pd.Series(
                name="INFO", data=[f"NS={i}" for i in nsamples]
            )

            # get format column: what type of metadata to expect
            colform = pd.Series(
                name="FORMAT", data=["GT"] * snpsmap.shape[0],
            )

            # get genotypes relative to REF/ALTS as strings (0/0, 0/1, ...)
            colgenos = pd.DataFrame({})
            for sname in self.snames:
                genostrs = [
                    "{}/{}".format(*sorted(k)) for k in [
                        list(j) for j in genos[self.snames.index(sname)]
                    ]
                ]
                genostrs = ["./." if i == "9/9" else i for i in genostrs]
                colgenos[sname] = genostrs

            # concat and order columns correctly
            infocols = pd.concat([vcfdf, colinfo, colform], axis=1)
            infocols = infocols[
                ["#CHROM", "POS", "ID", "REF", "ALT",
                 "QUAL", "FILTER", "INFO", "FORMAT"]]
            vcfdf = pd.concat([infocols, colgenos], axis=1)
            yield vcfdf

    def run(self):
        """Write chunks of VCF table to file."""
        with open(self.data.outfiles["vcf"], 'w', encoding="utf-8") as out:

            # print the VCF header.
            if self.data.is_ref:
                reference = Path(self.data.params.reference_sequence).name
            else:
                reference = "pseudo-reference (most common base at site)"
            header = VCFHEADER.format(
                date=time.strftime("%Y/%m/%d"),
                version=ipyrad.__version__,
                reference=reference
            )
            out.write(header + "\n")

            # write data table
            head = True
            for vcfdf in self._iter_vcf_data_chunk():
                vcfdf.to_csv(out, sep='\t', index=False, header=head)
                head = False




# def build_vcf(data, chunksize=1000):
#     """

#     """
#     # removed at init of Step function anyway.
#     if os.path.exists(data.outfiles.vcf):
#         os.remove(data.outfiles.vcf)

#     # dictionary to translate locus numbers to chroms
#     if data.isref:
#         revdict = chroms2ints(data, True)

#     # pull locus numbers and positions from snps database
#     with h5py.File(data.snps_database, 'r') as io5:

#         # iterate over chunks
#         for chunk in range(0, io5['genos'].shape[0], chunksize):

#             # if reference then psuedo ref is already ordered with REF first.
#             pref = io5['pseudoref'][chunk:chunk + chunksize]
#             snpmap = io5['snpsmap'][chunk:chunk + chunksize]

#             # load array chunks
#             if data.isref:
#                 genos = io5['genos'][chunk:chunk + chunksize, 1:, :]
#                 snames = data.snames[1:]

#                 # 1-indexed to 0-indexed (1/9/2019)
#                 chroms = [revdict[i - 1] for i in snpmap[:, 3]]
#                 ids = [
#                     "loc{}_pos{}".format(i - 1, j) for (i, j)
#                     in snpmap[:, [0, 2]]
#                 ]

#                 # reference based positions: pos on scaffold: 4, yes. tested.
#                 pos = snpmap[:, 4]
#                 #offset = 1

#             else:
#                 genos = io5['genos'][chunk:chunk + chunksize, :, :]
#                 snames = data.snames
#                 chroms = ["RAD_{}".format(i - 1) for i in snpmap[:, 0]]
#                 ids = [
#                     "loc{}_pos{}".format(i - 1, j) for (i, j)
#                     in snpmap[:, [0, 2]]
#                 ]
#                 # denovo based positions: pos on locus. tested. works. right.
#                 # almost. POS is 1 indexed.
#                 pos = snpmap[:, 2] + 1
#                 # offset = 0

#             # get alt genotype calls
#             alts = [
#                 b",".join(i).decode().strip(",")
#                 for i in pref[:, 1:].view("S1")
#             ]

#             # build df label cols
#             df_pos = pd.DataFrame({
#                 '#CHROM': chroms,            # scaff str names
#                 'POS': pos,                  # 1-indexed
#                 'ID': ids,                   # 0-indexed
#                 'REF': [i.decode() for i in pref[:, 0].view("S1")],
#                 'ALT': alts,
#                 'QUAL': [13] * genos.shape[0],
#                 'FILTER': ['PASS'] * genos.shape[0],
#             })

#             # get sample coverage at each site
#             nsamps = (
#                 genos.shape[1] - np.any(genos == 9, axis=2).sum(axis=1)
#             )

#             # store sum of coverage at each site
#             asums = []

#             # build depth columns for each sample
#             df_depth = pd.DataFrame({})
#             for sname in snames:

#                 # build geno strings
#                 genostrs = [
#                     "{}/{}".format(*k) for k in [
#                         i for i in [
#                             list(j) for j in genos[:, snames.index(sname)]
#                         ]
#                     ]
#                 ]

#                 # change 9's into missing
#                 genostrs = ["./." if i == "9/9" else i for i in genostrs]

#                 # genostrs = [
#                 # b"/".join(i).replace(b"9", b".").decode()
#                 # for i in genos[:, snames.index(sname)]
#                 # .astype(bytes)
#                 # ]

#                 # build depth and depthsum strings
#                 dpth = os.path.join(data.tmpdir, sname + ".depths.hdf5")
#                 with h5py.File(dpth, 'r') as s5:
#                     dpt = s5['depths'][chunk:chunk + chunksize]
#                     sums = [sum(i) for i in dpt]
#                     strs = [
#                         ",".join([str(k) for k in i.tolist()])
#                         for i in dpt
#                     ]

#                     # save concat string to name
#                     df_depth[sname] = [
#                         "{}:{}:{}".format(i, j, k) for (i, j, k) in
#                         zip(genostrs, sums, strs)]

#                     # add sums to global list
#                     asums.append(np.array(sums))

#             # make final columns
#             nsums = sum(asums)
#             colinfo = pd.Series(
#                 name="INFO",
#                 data=[
#                     "NS={};DP={}".format(i, j) for (i, j) in zip(nsamps, nsums)
#                 ])
#             colform = pd.Series(
#                 name="FORMAT",
#                 data=["GT:DP:CATG"] * genos.shape[0],
#                 )

#             # concat and order columns correctly
#             infocols = pd.concat([df_pos, colinfo, colform], axis=1)
#             infocols = infocols[
#                 ["#CHROM", "POS", "ID", "REF", "ALT",
#                  "QUAL", "FILTER", "INFO", "FORMAT"]]
#             arr = pd.concat([infocols, df_depth], axis=1)

#             # debugging
#             #print(arr.head())
#             ## PRINTING VCF TO FILE
#             ## choose reference string
#             if self.data.is_ref:
#                 reference = data.params.reference_sequence
#             else:
#                 reference = "pseudo-reference (most common base at site)"

#             header = VCFHEADER.format(
#                 date=time.strftime("%Y/%m/%d"),
#                 version=ipyrad.__version__,
#                 reference=os.path.basename(reference)
#                 )

#             with open(data.outfiles.vcf, 'a') as out:
#                 if chunk == 0:
#                     out.write(header)
#                     arr.to_csv(out, sep='\t', index=False)
#                 else:
#                     arr.to_csv(out, sep='\t', index=False, header=False)


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s7_assemble import Step7, SnpsDatabase
    ip.set_log_level("INFO", log_file="/tmp/test.log")

    DATA = ip.load_json("/tmp/TEST5.json")
    DATA.params.output_formats.append("v")

    # uncomment to include the ref
    DATA.hackers.exclude_reference = False
    print("EXCLUDE REF=", DATA.hackers.exclude_reference)

    # run it.
    with ip.Cluster(4) as ipyclient:
        step = Step7(DATA, ipyclient=ipyclient, force=True, quiet=False)
        step._split_clusters()
        step._apply_filters_and_trimming()
        step._collect_stats()

        db = SnpsDatabase(DATA)
        db.run()

        BuildVcf(DATA).run()
