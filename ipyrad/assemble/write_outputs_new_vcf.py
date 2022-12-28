#!/usr/bin/env python

"""Writes data from SNPs HDF5 to VCF format.

"""

from typing import TypeVar, Iterator, Tuple, List
from pathlib import Path
import time

import h5py
import numpy as np
import pandas as pd
from loguru import logger
from pkg_resources import get_distribution
# from ipyrad.assemble.utils import chroms2ints
from ipyrad.assemble.write_outputs_vcf_depths import CombinedDepths


logger = logger.bind(name="ipyrad")
Assembly = TypeVar("Assembly")
# CHUNKSIZE = 100
CHUNKSIZE = 2_000
VERSION = str(get_distribution('ipyrad')).split()[1]

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


class BuildVcfBase:
    def __init__(self, data: Assembly):
        self.data = data

        # attributes to be filled.
        # self.revdict = chroms2ints(self.data, keys_as_ints=True)
        """: A dict to convert chrom ints to chrom str names."""
        self.snames: List[str] = []
        """: A list of names in alphanumeric order, optionally including ref as sample."""

    def _iter_snps_data_chunk(self) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Yield chunks of data from SNPS HDF5."""
        with h5py.File(self.data.outfiles["snps_database"], 'r') as io5:
            self.snames = list(io5.attrs["names"])
            for start in range(0, io5['snpsmap'].shape[0], CHUNKSIZE): # pylint: disable=no-member
                snpsmap = io5['snpsmap'][start:start + CHUNKSIZE, :]
                alts = io5['alts'][start:start + CHUNKSIZE, :]
                genos = io5['genos'][start:start + CHUNKSIZE, :, :]
                yield snpsmap, alts, genos


class BuildVcfReference(BuildVcfBase):
    """Reference only object for building VCF. 
    
    This differs from the BuildVcfDenovo object in ...
    """
    def _iter_vcf_data_chunk(self) -> Iterator[pd.DataFrame]:
        """Writes chunks to VCF file."""
        for snpsmap, alts, genos in self._iter_snps_data_chunk():

            # get scaffold names, positions, and ID string.
            chrom_names = [self.revdict[i] for i in snpsmap[:, 3]]
            rad_ids = [
                f"loc{i}_pos{j}_scaff{x}_pos{y}" 
                for (i, j) in zip(snpsmap[:, 0], snpsmap[:, 2])
            ]


class BuildVcfDenovo(BuildVcfBase):
    """Denovo only object for building VCF. 
    
    This differs from the BuildVcfReference object in ...
    """
    def __init__(self, data: Assembly, snames: List[str]):
        self.snames = snames
        super().__init__(data)
        self.idepths = self._iter_vcf_depths()        

    def _iter_vcf_depths(self):
        """Yield depths for each SNP in order."""
        comb = CombinedDepths(self.data, self.snames)
        comb.open_handles()
        for snp in comb.iter_snp_depths():
            yield snp
        comb.close_handles()

    def _iter_vcf_data_chunk(self) -> Iterator[pd.DataFrame]:
        """Writes chunks to VCF file."""
        for snpsmap, alts, genos in self._iter_snps_data_chunk():

            # get CHROM strings of RAD locus at each SNP site.
            # >>> ['RADtag0', RADtag0', 'RADtag1', ...]
            chrom_names = [f"RADtag{i}" for i in snpsmap[:, 0]]

            # get ID strings for description loc+pos 0-indexed
            # >>> [loc0_snp0, loc0_snp1, loc1_snp0, ...]
            loc_pos_gen = zip(snpsmap[:, 0], snpsmap[:, 1])
            rad_ids = [f"loc{i}_snp{j}" for (i, j) in loc_pos_gen]

            # get ref alleles as a list of strings: 
            # >>> ["A", "T", "C", ...]
            ref = list(alts[:, 0].tobytes().decode())

            # get alt alleles as a list of single or comma-joined strings: 
            # >>> ["T", "A,C", "G,T", ...]
            alleles = []
            for idx in range(alts.shape[0]):
                calls = alts[idx, 1:]
                calls = calls[calls > 0]
                calls = ",".join(calls.tobytes().decode())
                alleles.append(calls)

            # build vcf dataframe for first 7 columns.
            # >>> #CHROM  POS     ID            REF     ALT     QUAL    FILTER
            # >>> RADtag0 8       loc0_snp0       G       A       13      PASS
            # >>> ...
            vcfdf = pd.DataFrame({
                '#CHROM': chrom_names,                  # chrom str name
                'POS': snpsmap[:, 4] + 1,  # *******    # 1-indexed position on scaff
                'ID': rad_ids,                          # RAD locus (x-positioned) and position 0-indexed.
                'REF': ref,                             # reference allele.
                'ALT': alleles,                         # other alleles in order.
                'QUAL': [13] * snpsmap.shape[0],        # arbitrarily high quality
                'FILTER': ['PASS'] * snpsmap.shape[0],  # no filter applied.
            })

            # fill depth information from ...
            # >>> 0/0:10:10,0,0,0  1/1:20:10,10,0,0  0/0:5:5,5,5,5 ...
            sumdepths = []
            indepths = []
            catgs = []
            for _ in range(min(CHUNKSIZE, snpsmap.shape[0])):
                site_depths = next(self.idepths)
                isum = 0
                idepth = []
                icatg = []
                for sample_depth in site_depths:
                    depth = int(sample_depth[0])
                    isum += depth
                    idepth.append(depth)
                    icatg.append(sample_depth[1])
                sumdepths.append(isum)
                indepths.append(idepth)
                catgs.append(icatg)

            # get sample coverage at each SNP site (column 8)
            # >>> INFO
            # >>> NS=10;DP=300
            nsamples = genos.shape[1] - np.any(genos == 9, axis=2).sum(axis=1)
            colinfo = pd.Series(
                name="INFO", data=[f"NS={i};DP={j}" for i, j in zip(nsamples, sumdepths)]
            )

            # get format column: what type of metadata to expect (column 9)
            # >>> FORMAT
            # >>> GT:DP:CATG
            colform = pd.Series(
                name="FORMAT", data=["GT:DP:CATG"] * snpsmap.shape[0],
            )

            # get genotype calls (index of REF+ALTS) as strings
            # this makes up part of columns 9 -> 9+nsamples
            # >>> 0/0, 0/1, 1/1, ...
            colgenos = pd.DataFrame({})
            for sidx, sname in enumerate(self.snames):

                genlist = (sorted(j) for j in genos[:, sidx])
                genostrs = (f"{i[0]}/{i[1]}" for i in genlist)
                genostrs = ["./." if i == "9/9" else i for i in genostrs]
                colgenos[sname] = genostrs

            # get depths [[10, 20, 12, 15, ...], [50, 55, 60, 33, ...]...]
            coldepths = pd.DataFrame(np.array(indepths), columns=self.snames)

            # get catgs [['0,0,0,10', '0,0,5,5', ...], ['50,0,0,0', '0,0,30,0'...]]
            colcatgs = pd.DataFrame(np.array(catgs), columns=self.snames)

            # join colgenos, coldepths and colcatgs
            for sname in colgenos.columns:
                colgenos[sname] = [f"{i}:{j}:{z}" for (i, j, z) in 
                    zip(colgenos[sname], coldepths[sname], colcatgs[sname])
                ]

            # concat and order columns correctly
            infocols = pd.concat([vcfdf, colinfo, colform], axis=1)
            infocols = infocols[
                ["#CHROM", "POS", "ID", "REF", "ALT",
                 "QUAL", "FILTER", "INFO", "FORMAT"]]
            vcfdf = pd.concat([infocols, colgenos], axis=1)
            yield vcfdf

    def run(self):
        """Write chunks of VCF table to file."""
        # TODO: Gzip it? BGZip it?
        with open(self.data.outfiles["vcf"], 'w', encoding="utf-8") as out:

            # print the VCF header.
            if self.data.is_ref:
                reference = Path(self.data.params.reference_sequence).name
            else:
                reference = "pseudo-reference (most common base at site)"
            header = VCFHEADER.format(
                date=time.strftime("%Y/%m/%d"),
                version=VERSION,
                reference=reference
            )
            out.write(header)

            # write data table in chunks at a time.
            head = True
            for vcfdf in self._iter_vcf_data_chunk():
                vcfdf.to_csv(out, sep='\t', index=False, header=head)
                head = False



if __name__ == "__main__":

    # from ipyrad.assemble.s7_assemble import Step7, SnpsDatabase
    import ipyrad as ip
    ip.set_log_level("INFO", log_file="/tmp/test.log")

    # DATA = ip.load_json("../../sra-fastqs/cyatho.json")
    DATA = ip.load_json("../../tests/ipsim.json")
    DATA.tmpdir = DATA.params.project_dir / "ipsim_tmp_outfiles"
    DATA.outfiles['vcf'] = "/tmp/test.vcf"
    DATA.drop_ref = False

    v = BuildVcfDenovo(DATA, DATA.samples)    
    v.run()

    # for i in v._iter_snps_data_chunk():
        # print([x.shape for x in i])

    # for i in v.idepths:
        # print(i)

    # for i in v._iter_vcf_data_chunk():
        # print(i.iloc[:, [2,3,4,7,8,9,10]])

    raise SystemExit()

    # pd.set_option('max_colwidth', 1000)
    # print(y.iloc[:5, :])
    # print(y.iloc[-5:, :])
    # v.run()


    # DATA = ip.load_json("/tmp/TEST5.json")
    # DATA.params.output_formats.append("v")

    # # uncomment to include the ref
    # DATA.hackers.exclude_reference = False
    # print("EXCLUDE REF=", DATA.hackers.exclude_reference)

    # # run it.
    # with ip.Cluster(4) as ipyclient:
    #     step = Step7(DATA, ipyclient=ipyclient, force=True, quiet=False)
    #     step._split_clusters()
    #     step._apply_filters_and_trimming()
    #     step._collect_stats()

    #     db = SnpsDatabase(DATA)
    #     db.run()

    #     BuildVcf(DATA).run()
