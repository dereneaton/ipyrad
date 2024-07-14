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
        self.testdir.mkdir(exist_ok=True)
        self.tmpdir = self.testdir / "assembly"
        self.datadir = self.testdir / "ipsimdata"

        self.data1 = ip.Assembly(name="TEST-denovo-se")
        self.data1.params.project_dir = self.tmpdir
        self.data2 = ip.Assembly(name="TEST-denovo-pe")
        self.data2.params.project_dir = self.tmpdir

    def test_load_single_end_fastqs(self):
        """Run basic step1 and check Types of stored results."""
        self.data1.params.fastq_paths = "../../../sra-fastqs/4*.fastq"
        self.data1.run('1', force=True, cores=6)

        # select a sample
        self.assertEqual(len(self.data1.samples), 3)
        sample = list(self.data1.samples.values())[0]

        # SE fastqs should be a list with tuples of (Path, '')
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
        [self.assertEqual(i[1].name, "null") for i in sample.files.fastqs]

        # SE trimmed should a tuple of (Path, '')
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)
        self.assertEqual(sample.files.trimmed[1].name, "null")

        # reload from JSON and test types again
        self.data1 = ip.load_json(self.data1.json_file)
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
        [self.assertEqual(i[1].name, "null") for i in sample.files.fastqs]
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)
        self.assertEqual(sample.files.trimmed[1].name, "null")

    def test_load_single_end_fastqs_from_multiple_paths(self):
        """Run basic step1 and check Types of stored results."""
        self.data = self.data1.branch("tmp")
        self.data.params.fastq_paths = [
            "../../../sra-fastqs/40*.fastq",
            "../../../sra-fastqs/2*.fastq",
        ]

        # run
        self.data.run('1', force=True, cores=6)

        # write to JSON and reload
        self.data = ip.load_json(self.data.json_file)

        # select a sample
        self.assertEqual(len(self.data.samples), 2)
        sample = list(self.data.samples.values())[0]

        # SE fastqs should be a list with tuples of (Path, '')
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
        [self.assertEqual(i[1].name, "null") for i in sample.files.fastqs]

        # SE trimmed should a tuple of (Path, '')
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)
        self.assertEqual(sample.files.trimmed[1].name, "null")

    def test_load_single_end_fastqs_with_techs_merged(self):
        # select two sample's data as inputs
        self.data = self.data1.branch("tmp")
        self.data.params.fastq_paths = [
            "../../../sra-fastqs/2*.fastq",
            "../../../sra-fastqs/40*.fastq",
        ]
        # designate that these two samples are replicates to be merged
        self.data.params.technical_replicates = {
            "TECH": ["29154_superba_SRR1754715", "40578_rex_SRR1754724"],
        }
        # load and trim the reads
        self.data.run('1', force=True, cores=6)
        sample = list(self.data.samples.values())[0]
        # check fastqs types
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
        self.assertEqual(sample.files.fastqs[0][1], Path("null"))
        # check trimmed types
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)
        self.assertEqual(sample.files.trimmed[1], Path("null"))
        # expected merged count for TECH sample
        self.assertEqual(self.data.stats.loc["TECH", "reads_raw"], 2404936)

    def test_load_paired_end_fastqs(self):
        """Load PE data samples one path entered as the only item in a list."""
        self.data2.params.fastq_paths = ["../../../pedtest/demux_2023-3-28/linea*.gz"]
        self.data2.run('1', force=True, cores=2)

        # select a sample
        self.assertEqual(len(self.data2.samples), 1)
        sample = list(self.data2.samples.values())[0]

        # PE fastqs should be a list with tuples of (Path, '')
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]

        # PE trimmed should a tuple of (Path, '')
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)

        # reload from JSON and test types again
        self.data2 = ip.load_json(self.data2.json_file)
        self.assertIsInstance(sample.files.fastqs, list)
        [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
        [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
        self.assertIsInstance(sample.files.trimmed, tuple)
        self.assertIsInstance(sample.files.trimmed[0], Path)
        self.assertIsInstance(sample.files.trimmed[1], Path)

    # def test_load_...

    # def test_load_paired_end_fastqs_with_techs_merged(self):
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
