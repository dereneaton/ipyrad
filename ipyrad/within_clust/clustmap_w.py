# !/usr/bin/env python

"""Class to organize within-sample clust/map operations.

"""

from loguru import logger
from ipyrad.core import BaseStep, track_remote_jobs
from ipyrad.within_clust.clustmap_w_funcs import (
    dereplicate,
    tag_inline_for_decloning,
    map_to_reference_as_filter,
)
logger = logger.bind(name="ipyrad")


class ClustMapBase(BaseStep):
    """Load and filter/trim fastq input files to create Samples."""
    def __init__(self, data, force, ipyclient):
        super().__init__(data, step=2, force=force)
        # self.samples: Dict[str, Sample]
        # self.stepdir: Path
        # self.tmpdir: Path
        self.ipyclient = ipyclient

    def run(self):
        """Different methods apply to reference vs denovo."""
        # reference:
        # [unique] .index_references
        # [shared] .mapping_to_reference_filter
        # [unique] .join_pairs_for_derep
        # [shared] .tag_for_decloning
        # [shared] .dereplicate
        # [shared] .tag_back_to_header
        # [unique] .mapping

        # denovo:
        # [unique] .index_reference_as_filter
        # [shared] .remote_map_to_reference_filter
        # [unique] .remote_pair_merge_overlaps_with_vsearch
        # [unique] .remote_pair_join_unmerged_end_to_end
        # [shared] .remote_tag_inline_for_decloning
        # [unique] .check_strand_type
        # [shared] .remote_dereplicate
        # [shared] .tag_back_to_header
        # [unique] .cluster

    def remote_map_to_reference_as_filter(self):
        """Map reads to filter reference to get unmapped fastq.

        # i1: trimmed/{}_trimmed_R[1,2].fastq.gz
        # i0: tmpdir/{}_concat_R[1,2].fastq.gz
        # o: tmpdir/{}_unmapped_R[1,2].fastq
        """
        if not self.data.params.reference_as_filter:
            return

        # run this 4-threaded to be nice to I/O
        threads = self.data.ipcluster['threads']
        thview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::threads])
        rasyncs = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            func = map_to_reference_as_filter
            rasyncs[sname] = thview.apply(func, *args)
        results = track_remote_jobs(rasyncs, self.ipyclient)

        # store stats
        for sname, val in results.items():
            sample = self.samples[sname]
            sample.stats_s2.reads_mapped_to_ref_filter = val[0]
            sample.stats_s2.reads_mapped_to_ref_filter_prop = val[1]

    def remote_tag_inline_for_decloning(self):
        """Moves the i5 tags from the index to start of reads for decloning.

        # i: tmpdir/{}_merged.fa
        # o: tmpdir/{}_decloned.fa
        """
        if (not self.data.is_pair) or (not self.data.params.declone_PCR_duplicates):
            return

        lbview = self.ipyclient.load_balanced_view()
        rasyncs = {}
        for sname, sample in self.samples.items():
            args = (self.data, self.samples[sname])
            func = tag_inline_for_decloning
            rasyncs[sname] = lbview.apply(func, *args)
        track_remote_jobs(rasyncs, self.ipyclient)

    def remote_dereplicate(self):
        """Moves the i5 tags from the index to start of reads for decloning.

        # i3: trimmed/{}_trimmed.fastq.gz         # se data
        # i2: tmpdir/{}_concat.fastq.gz           # se assembly merge
        # i1: tmpdir/{}_merged.fa                 # pe data
        # i0alt: tmpdir/{}_joined.fa              # pe data
        # i0: tmpdir/{}_decloned.fa               # pe w/ declone
        # o: tmpdir/{}_derep.fa
        """
        logger.info("dereplicating identical reads")
        threads = self.data.ipcluster['threads']
        thview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::threads])
        rasyncs = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            func = dereplicate
            rasyncs[sname] = thview.apply(func, *args)
        track_remote_jobs(rasyncs, self.ipyclient)


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    # RUN EXAMPLES IN DENOVO OR REFERENCE MODULES
