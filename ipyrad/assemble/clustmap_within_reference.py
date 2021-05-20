#!/usr/bin/env python

"""
Reference based assembly methods
"""

from loguru import logger


class ClustMapReference:
    """
    de novo assembly pipeline.
    """
    def __init__(self, step):
        # inherit attrs from step
        self.data = step.data
        self.samples = step.samples
        self.params = step.data.params
        self.tmpdir = step.tmpdir
        self.stepdir = step.stepdir
        self.quiet = step.quiet
        self.lbview = step.lbview
        self.thview = step.thview
        self.is_paired = "pair" in self.params.datatype


    def index_references(self):
        """
        Index reference and/or reference_filter files.
        """
        jobs = {}
        if self.data.params.reference_sequence:
            logger.debug("indexing reference with bwa and samtools")
            rasync1 = self.lbview.apply(index_ref_with_bwa, self.data)
            rasync2 = self.lbview.apply(index_ref_with_sam, self.data)
            jobs['bwa_index_ref'] = rasync1
            jobs['sam_index_ref'] = rasync2

        if self.data.params.reference_as_filter:
            logger.debug("indexing reference_filter with bwa and samtools")
            rasync1 = self.lbview.apply(index_ref_with_bwa, self.data, alt=1)
            rasync2 = self.lbview.apply(index_ref_with_sam, self.data, alt=1)
            jobs['bwa_index_alt'] = rasync1
            jobs['sam_index_alt'] = rasync2

        # track job
        msg = "indexing reference"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()        

