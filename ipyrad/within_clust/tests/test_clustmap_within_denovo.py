#!/usr/bin/env python


"""Unittests for step 2 denovo scenarios.

The object of step 2 denovo is to:
    - [optionally] map to reference_as_filter to discard matched reads.
    - [optionally] tag sequences with i5 tag from headers to declone PCR dups. (WARNING if i5 tags present but not used?)
    - [optionally] merge (if overlapped) else join PE reads to de-replicate as a pair.
    - dereplicate reads using vsearch.
    - cluster dereplicated reads using vsearch.
    - store files to Samples at:
        - .clustmap : Tuple[Path, Path] = (_derep, _matches)
    - store stats to Samples at:
        - merged_read_pairs
        - merged_read_pairs_prop

Tests
-----
1. Load fastqs from a single location
2. Load fastqs from multiple locations
3. Load fastqs from a single location with technical replicates
4. Load fastqs from multiple locations with technical replicates
5. ...
"""


import unittest
from pathlib import Path
import ipyrad as ip


class TestClustMapWithin(unittest.TestCase):

    def setUp(self):
        # SE dataset
        self.data1 = ip.load_json("/tmp/ipyrad-tests/TEST-denovo-se.json")
        # PE dataset
        self.data2 = ip.load_json("/tmp/ipyrad-tests/TEST-denovo-pe.json.")

    def denovo_cluster_se_data(self):
        self.data1.run("2", cores=6, threads=2, force=True)

        # select a sample
        self.assertEqual(len(self.data.samples), 3)
        sample = list(self.data.samples.values())[0]


    # def test_load_single_end_fastqs(self):
    #     """Run basic step1 and check Types of stored results."""
    #     self.data.params.fastq_paths = "../../../sra-fastqs/4*.fastq"
    #     self.data.run('1', force=True, cores=6)


    #     # SE fastqs should be a list with tuples of (Path, '')
    #     self.assertIsInstance(sample.files.fastqs, list)
    #     [self.assertIsInstance(i[0], Path) for i in sample.files.fastqs]
    #     [self.assertIsInstance(i[1], Path) for i in sample.files.fastqs]
    #     [self.assertEqual(i[1].name, "null") for i in sample.files.fastqs]


if __name__ == "__main__":
    pass
