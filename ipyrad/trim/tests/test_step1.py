#!/usr/bin/env python

"""Unittests for step 1 scenarios.

The object of step 1 is to:
    - load fastq data from one or more folders
    - identify if the data is single or paired-end
    - identify if technical replicates are present and optionally merge
    - identify restriction overhang is present and should be trimmed
    - run fastp to trim reads with the correct arguments.
    - store results to Samples at:
        - .fastqs : List[Path]
        - .trimmed: List[Tuple[Path, Path | None]]

Tests
-----
1. Load fastqs from a single location
2. Load fastqs from multiple locations
3. Load fastqs from a single location with technical replicates
4. Load fastqs from multiple locations with technical replicates
5. ...

# raise Warning, drop one or more samples, but proceed.
x. One sample fails because no reads pass filtering

# raise IPyradError() exception.
x. The only entered fastq file path does not exist
x. One out of several entered fastq file paths does not exist.
y. All samples fail because no reads pass filtering

"""

import unittest
from pathlib import Path
import ipyrad as ip


class TestLoadFastqs(unittest.TestCase):

    def setUp(self):
        # setup dirs
        self.testdir = Path("/tmp/ipyrad-tests")
        self.tmpdir = self.testdir / "assembly"
        self.datadir = self.testdir / "ipsimdata"

        self.data = ip.Assembly(name="TEST")
        self.data.params.project_dir = self.tmpdir

    # def test_load_single_end_fastqs(self):
    #     """Run basic step1 and check Types of stored results."""
    #     self.data.params.fastq_paths = "../../../sra-fastqs/4*.fastq"
    #     self.data.run('1', force=True, cores=6)

    #     # select a sample
    #     sample = list(self.data.samples.values())[0]

    #     # SE fastqs should be a list with tuples of (Path, '')
    #     self.assertIsInstance(sample.files.fastqs, list)
    #     [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
    #     [self.assertEqual(i[1], "") for i in sample.files.fastqs]

    #     # SE trimmed should a tuple of (Path, '')
    #     self.assertIsInstance(sample.files.trimmed, tuple)
    #     self.assertIsInstance(sample.files.trimmed[0], Path)
    #     self.assertEqual(sample.files.trimmed[1], "")

    # def test_load_single_end_fastqs_from_multiple_paths(self):
    #     """Load samples from more than one location."""
    #     self.data.params.fastq_paths = [
    #         "../../../sra-fastqs/40*.fastq",
    #         "../../../sra-fastqs/2*.fastq",
    #     ]
    #     self.data.run('1', force=True, cores=6)
    #     sample = list(self.data.samples.values())[0]
    #     self.assertIsInstance(sample.files.fastqs, list)
    #     self.assertEqual(len(sample.files.fastqs), 1)
    #     self.assertEqual(len(self.data.samples), 2)
    #     self.assertEqual(self.data.stats.loc["29154_superba_SRR1754715", "raw_reads"], ...)

    # def test_load_paired_end_fastqs(self):
    #     """Load PE data samples one path entered as the only item in a list."""
    #     self.data.params.fastq_paths = ["../../../pedtest/demux_2023-3-28/linea*.gz"]
    #     self.data.run('1', force=True, cores=2)
    #     sample = list(self.data.samples.values())[0]

    #     # PE fastqs should be a list with tuples of (Path, '')
    #     self.assertIsInstance(sample.files.fastqs, list)
    #     [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
    #     [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]

    #     # PE trimmed should a tuple of (Path, '')
    #     self.assertIsInstance(sample.files.trimmed, tuple)
    #     self.assertIsInstance(sample.files.trimmed[0], Path)
    #     self.assertIsInstance(sample.files.trimmed[1], Path)

    # def test_load_single_end_fastqs_with_techs_merged(self):
    #     # select two sample's data as inputs
    #     self.data.params.fastq_paths = [
    #         "../../../sra-fastqs/2*.fastq",
    #         "../../../sra-fastqs/40*.fastq",
    #     ]
    #     # designate that these two samples are replicates to be merged
    #     self.data.params.technical_replicates = {
    #         "TECH": ["29154_superba_SRR1754715", "40578_rex_SRR1754724"],
    #     }
    #     # load and trim the reads
    #     self.data.run('1', force=True, cores=6)
    #     # fastqs should be [(R1, ''), (R1, '')]
    #     sample = list(self.data.samples.values())[0]
    #     self.assertIsInstance(sample.files.fastqs, list)
    #     [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
    #     [self.assertEqual(i[1], "") for i in sample.files.fastqs]
    #     self.assertIsInstance(sample.files.trimmed, tuple)
    #     self.assertIsInstance(sample.files.trimmed[0], Path)
    #     self.assertEqual(sample.files.trimmed[1], "")
    #     self.assertEqual(self.data.stats.loc["TECH", "reads_raw"], 2404936)

    # def test_load_single_end_fastqs_with_techs_merged_raise_bad_name(self):
    #     with self.assertRaises(ValueError):
    #         self.data.params.fastq_paths = [
    #             "../../../sra-fastqs/2*.fastq",
    #             "../../../sra-fastqs/40*.fastq",
    #         ]
    #         self.data.params.technical_replicates = {
    #             "TECH": ["29154_superba_SRR1754715", "XXX"],
    #         }
    #         self.data.run('1', force=True, cores=2)

    # def test_5_load_paired_end_fastqs_with_techs_not_merged(self):
    #     pass

    # def test_user_can_load_dataset_after_fastq_inputs_are_moved(self):
    #     """
    #     Ideally the user can load their Assembly file even if it cannot
    #     find data files anymore as a way to access the stats, or to
    #     re-run only a later step for which data paths are still intact.

    #     Maybe only do this in terms of the orig fastq data paths.
    #     """
    #     pass


if __name__ == "__main__":

    # this prevents logger.exceptions from raising themselves
    # since caught exceptions will have explicit Exception types.
    # ip.set_log_level("CRITICAL")
    unittest.main()
