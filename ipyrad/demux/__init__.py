#!/usr/bin/env python

"""Example usage of the Demux class for demultiplexing reads.

API
---
>>> tool = Demux(
>>>     fastq_paths="../../pedtest/Pedicularis_plate1_R*.fastq.gz",
>>>     barcodes_path="../../pedtest/barcodes-fewer-plate1.csv",
>>>     outpath="../../pedtest/demux_2024-3-16",
>>>     max_barcode_mismatch=1,
>>>     cores=7,
>>>     chunksize=1e6,
>>> )
>>> tool.run()

CLI
---
$ ipyrad demux -d DATA -b BARCODES -o /tmp -m 1 -c 7 --chunksize 1e6

"""
