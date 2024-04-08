#!/usr/bin/env python

"""Tests of the demultiplexing subcommand.

- test_demux_se
- test_demux_pe_single_barcodes
- test_demux_pe_combinatorial_barcodes
- test_demux_pe_i7
- test_demux_se_i7
- test_demux_pe_mismatch_values
- test_demux_pe_chunksize_too_big
- test_demux_merge_replicates
- test_demux_autodetect_re
- test_demux_autodetect_re_override
- test_demux_bad_data_param
- test_demux_bad_barcode_param
- test_demux_bad_outdir_param

"""

import sys
import shutil
from pathlib import Path
import unittest
import subprocess
import tarfile
import ipyrad as ip
from ipyrad.core.logger_setup import capture_logs


CLI = Path(sys.prefix) / "bin" / "ipyrad"


class TestLoadFastqs(unittest.TestCase):

    def setUp(self):
        # create tmpdir for ipyrad files
        self.testdir = Path("/tmp/ipyrad-tests")
        self.tmpdir = self.testdir / "assembly"
        self.datadir = self.testdir / "ipsimdata"

        # ensure testdir exists
        self.testdir.mkdir(exist_ok=True)

        # ensure tmpdir doesn't exist
        if self.tmpdir.exists():
            shutil.rmtree(self.testdir)

        # TODO: download the ipsimdata.tar.gz archive from URL
        # ...

        # ensure data dir exists
        if not self.datadir.exists():
            with tarfile.open("../../../tests/ipsimdata-1.tar.gz") as itmp:
                itmp.extractall(self.testdir)

    # def test_demux_se(self):
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "rad_example_R1*.gz",
    #         barcodes_path=self.datadir / "rad_example_barcodes.txt",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         cores=7,
    #         chunksize=1e6,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     areads = tool._sample_stats["1A_0"]
    #     self.assertEqual(nreads, 239866)
    #     self.assertEqual(areads, 19862)
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_se_cli(self):
    #     """Demultiplex SE data using the CLI."""
    #     cmd = [
    #         CLI, "demux",
    #         "-d", self.datadir / "rad_example_R1*.gz",
    #         "-b", self.datadir / "rad_example_barcodes.txt",
    #         "-o", self.tmpdir,
    #         "-c", "4",
    #         "-m", "1",
    #     ]
    #     subprocess.run(cmd, check=True)
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_se_small_chunks(self):
    #     """Demultiplex SE data using small chunk size and verify result."""
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "rad_example_R1*.gz",
    #         barcodes_path=self.datadir / "rad_example_barcodes.txt",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         cores=4,
    #         chunksize=701,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     areads = tool._sample_stats["1A_0"]
    #     self.assertEqual(nreads, 239866)
    #     self.assertEqual(areads, 19862)
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_se_tech_reps_not_merged(self):
    #     """Demultiplex SE data with technical replicates not merged and verify result."""
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "rad_example_R1*.gz",
    #         barcodes_path=self.datadir / "rad_example_barcodes_techreps.txt",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         merge_technical_replicates=False,
    #         cores=4,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     areads = tool._sample_stats["1A_0"]
    #     a0reads = tool._sample_stats["1A_0-technical-replicate-0"]
    #     self.assertEqual(nreads, 239866)  # total reads
    #     self.assertEqual(areads, 0)       # rep reads are not summed
    #     self.assertEqual(a0reads, 19862)  # techrep0 reads
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_se_tech_reps_merged(self):
    #     """Demultiplex SE data w/ technical replicates merged and verify result."""
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "rad_example_R1*.gz",
    #         barcodes_path=self.datadir / "rad_example_barcodes_techreps.txt",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         merge_technical_replicates=True,
    #         cores=4,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     areads = tool._sample_stats["1A_0"]
    #     a0reads = tool._sample_stats["1A_0-technical-replicate-0"]
    #     self.assertEqual(nreads, 239866)  # total reads
    #     self.assertEqual(areads, 60041)   # rep reads are summed
    #     self.assertEqual(a0reads, 0)      # techrep0 was cleared
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_pe_single_barcodes(self):
    #     """Demultiplex paired-end data single barcodes and verify result."""
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "pairddrad_example_R*.gz",
    #         barcodes_path=self.datadir / "pairddrad_example_barcodes.txt",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         cores=4,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     self.assertEqual(nreads, 239812)
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_pe_combinatorial_barcodes(self):
    #     """Demultiplex paired-end data comb barcodes and verify result."""
    #     tool = ip.Demux(
    #         fastq_paths=self.datadir / "small_tmp_R*.gz",
    #         barcodes_path=self.datadir / "barcodes-fewer-plate1.csv",
    #         outpath=self.tmpdir,
    #         max_barcode_mismatch=1,
    #         cores=4,
    #     )
    #     tool.run()
    #     nreads = sum(tool._sample_stats.values())
    #     self.assertEqual(nreads, 102528)
    #     shutil.rmtree(self.tmpdir)

    # def test_demux_override_re(self):
    #     """Demultiplex with re overrides by user and verify warnings."""
    #     with capture_logs("INFO") as cap:
    #         tool = ip.Demux(
    #             fastq_paths=self.datadir / "small_tmp_R*.gz",
    #             barcodes_path=self.datadir / "barcodes-fewer-plate1.csv",
    #             outpath=self.tmpdir,
    #             max_barcode_mismatch=1,
    #             re1="ATCG",
    #             re2="CGATC",
    #             cores=4,
    #         )
    #         tool.run()
    #         shutil.rmtree(self.tmpdir)

    #     # check that warning were raised.
    #     self.assertTrue([
    #         "WARNING:ipyrad.demux.demux:user entered ATCG"
    #         in logmsg for logmsg in cap
    #     ])
    #     self.assertTrue([
    #         "WARNING:ipyrad.demux.demux:user entered CGATC"
    #         in logmsg for logmsg in cap
    #     ])

    def test_demux_pe_i7(self):
        """..."""

    def test_demux_se_i7(self):
        """This is impossible, should raise an error."""


if __name__ == "__main__":
    unittest.main()
