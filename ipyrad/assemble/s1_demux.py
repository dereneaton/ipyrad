#!/usr/bin/env python

"""Step 1 runs either Demultiplexing or Loading of input data files.

Step1 uses one of the following two classes:
- Demultiplexing: SimpleDemux()
- Loading: FileLinker()

Supported scenarios (se and pe)
-------------------------------
- load demuxed data
- demux on inline barcodes
- demux on i7 barcodes
- demux on combinatorial inline barcodes (pair3rad)
- demux on combinatorial i7 barcodes?

"""

from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.utils import IPyradError
from ipyrad.assemble.demux_sorted import FileLinker
from ipyrad.assemble.demux_raw import SimpleDemux


class Step1(BaseStep):
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, step=1, quiet=quiet, force=force)
        self.ipyclient = ipyclient
        self.pre_check()

    def pre_check(self):
        """Check that either sorted_fastq_path exists, or raw_fastq_path
        and barcodes_path both exist.
        """
        if not self.data.params.sorted_fastq_path:
            if not self.data.params.raw_fastq_path:
                raise IPyradError(
                    "You must enter either a sorted_fastq_path or "
                    "raw_fastq_path in param settings.")
        if self.data.params.raw_fastq_path:
            if not self.data.params.barcodes_path:
                raise IPyradError(
                    "You must enter a barcodes_path in param settings "
                    "when demultiplexing from raw_fastq_path.")

    def run(self):
        """Runs a different Step1 class depending on input data method"""
        if self.data.params.sorted_fastq_path is not None:
            FileLinker(self).run()
        else:
            SimpleDemux(self.data, quiet=quiet, ipyclient=ipyclient).run()


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    # LOADING PRE-DEMUX'D DATA.
    # TEST = ip.Assembly("PEDIC")
    # TEST.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    # TEST.params.project_dir = "/tmp"
    # TEST.run('1', force=True, quiet=True)

    # ONE FASTQ FILE AND UNIQUE BARCODES TEST (SE RAD)
    TEST = ip.Assembly("TEST1")
    TEST.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_R1*.gz"    
    TEST.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes.txt"
    TEST.params.project_dir = "/tmp"
    TEST.params.max_barcode_mismatch = 1
    TEST.run('1', force=True, quiet=True)

    # # MULTIPLE FASTQ FILES AND ONE UNIQUE BARCODES TEST
    # TESTX = TEST.branch("TEST2")
    # TESTX.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_*.gz"
    # TESTX.run('1', force=True, quiet=True)

    # # ONE FASTQ FILE AND TECHNICAL REPLICATES IN BARCODES (hackers on)
    # TESTX = TEST.branch("TEST3")
    # TESTX.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes_techreps.txt"
    # TESTX.run('1', force=True, quiet=True)

    # # ONE FASTQ FILE AND TECHNICAL REPLICATES IN BARCODES (hackers off)    
    # TESTX = TEST.branch("TEST4")
    # TESTX.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes_techreps.txt"
    # TESTX.hackers.merge_technical_replicates = False
    # TESTX.run('1', force=True, quiet=True)

    # DEMUX PAIRED_END TEST
    TESTX = TEST.branch("TEST5")
    TESTX.params.raw_fastq_path = "../../tests/ipsimdata/pairddrad_example_*.gz"    
    TESTX.params.barcodes_path = "../../tests/ipsimdata/pairddrad_example_barcodes.txt"
    TESTX.params.datatype = "pairddrad"
    TESTX.run('1', force=True, quiet=True)

    # DEMUX PAIRED_3RAD TEST (must have combinatorial barcodes)
    # TESTX = TEST.branch("TEST5")
    # TESTX.params.raw_fastq_path = "../../tests/ipsimdata/pairddrad_example_*.gz"    
    # TESTX.params.barcodes_path = "../../tests/ipsimdata/pairddrad_example_barcodes.txt"
    # TESTX.params.datatype = "pair3rad"
    # TESTX.run('1', force=True, quiet=True)    

    # EMPIRICAL DEMUX ON i7 outer tags...
    # TODO

    # EMPIRICAL ...

