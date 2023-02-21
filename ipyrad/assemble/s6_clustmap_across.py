#!/usr/bin/env python

"""Detect orthology across samples and write to fasta.

Denovo: Cluster across samples w/ vsearch and build clusters.
Reference: Extract mapping from bam files and bedtools to build clusters
  and apply filters to identify bad mapping.
"""

from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.clustmap_across_reference import ClustMapAcrossReference
from ipyrad.assemble.clustmap_across_denovo import ClustMapAcrossDenovo


class Step6(BaseStep):
    """Group loci across samples by orthology using clustering in denovo
    or by pulling from map positions in reference.
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 6, quiet, force)
        self.is_ref = self.data.is_ref
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.data.ncpus = len(self.ipyclient.ids)
        self.regions = []

    def run(self):
        """Runs the set of methods for denovo or reference method."""
        if self.data.params.assembly_method == "denovo":
            ClustMapAcrossDenovo(self).run()
        else:
            ClustMapAcrossReference(self).run()
        for _, sample in self.samples.items():
            sample.state = 6
            sample._clear_old_results()
        self.data.save_json()


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")# log_file="/tmp/test.log")
   
    TEST = ip.load_json("../../pedtest/NEW.json")
    TEST.params.min_depth_majrule = 1
    TEST.run("6", force=True, quiet=True)
    print(TEST.stats)

    # TEST = ip.load_json("/tmp/TEST5.json")
    # TEST.run("6", force=True, quiet=False)

    # tdata = ip.load_json("/tmp/test-amaranth.json")
    # tdata.run("6", auto=True, force=True)
    # print(tdata.stats)
