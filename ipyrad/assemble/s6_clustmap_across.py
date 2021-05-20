#!/usr/bin/env python

"""
Denovo: Cluster across samples w/ vsearch and build clusters.
Reference: Extract mapping from bam files and bedtools to build clusters
  and apply filters to identify bad mapping.
"""


from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.clustmap_across_reference import ClustMapAcrossReference
from ipyrad.assemble.clustmap_across_denovo import ClustMapAcrossDenovo


class Step6(BaseStep):
    """
    Group loci across samples by orthology using clustering in denovo
    or by pulling from map positions in reference.
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 6, quiet, force)
        self.is_ref = bool('ref' in self.data.params.assembly_method)
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.data.ncpus = len(self.ipyclient.ids)
        self.data.tmpdir = self.tmpdir
        self.regions = []


    def run(self):
        """
        Runs the set of methods for denovo or reference method
        """
        # DENOVO
        if self.data.params.assembly_method == "denovo":
            ClustMapAcrossDenovo(self).run()
        else:
            ClustMapAcrossReference(self).run()



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG", stderr=False, logfile="/tmp/test.log")
   
    TEST = ip.load_json("/tmp/TEST1.json")
    TEST.run("6", force=True, quiet=False)

    # tdata = ip.load_json("/tmp/test-amaranth.json")
    # tdata.run("6", auto=True, force=True)
    # print(tdata.stats)
