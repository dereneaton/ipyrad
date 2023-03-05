#!/usr/bin/env python

"""Reference based assembly method.

"""

from loguru import logger
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.clustmap_within_both import ClustMapBase
from ipyrad.assemble.clustmap_within_reference_utils import (
    index_ref_with_bwa,
    index_ref_with_sam,
    join_pairs_for_derep_ref,
    split_derep_pairs_ref,
    mapping_reads,
    write_clusters,
    # build_clusters_from_cigars,
)

logger = logger.bind(name="ipyrad")


class ClustMapReference(ClustMapBase):
    """Reference mapping assembly pipeline."""
    def __init__(self, step):
        super().__init__(step)

    def index_references(self):
        """Index reference and/or reference_filter files."""
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

    def join_pairs_for_derep(self):
        """Joins pairs end to end for decloning and dereplicating.

        i: tmpdir/[trim, concat, or unmapped].fastq[.gz]
        o: tmpdir/joined.fastq
        """
        if not self.data.is_pair:
            return
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                join_pairs_for_derep_ref,
                *(self.data, self.samples[sname])
            )
        if jobs:
            msg = "joining pairs for derep/declone"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def split_derep_pairs_for_mapping(self):
        """Dereplicate read(pairs) to fasta format using vsearch.

        i: tmpdir/{}_[derep,tagged].fa
        o: tmpdir/{}_derep_split_R[1,2].fa
        """
        if not self.data.is_pair:
            return
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                split_derep_pairs_ref,
                *(self.data, self.samples[sname])
            )
        if jobs:
            msg = "splitting derep pairs for mapping"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def mapping_to_reference(self):
        """Map reads to reference to get mapped bams.

        i: tmpdir/{}_derep_[split]_R[1,2].fa
        o: tmpdir/{}.bam
        """
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.thview.apply(
                mapping_reads,
                *(self.data, self.samples[sname], self.data.ipcluster['threads'])
            )
        if jobs:
            msg = "mapping reads to reference"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def build_clusters_from_cigars(self):
        """Builds clusters.gz files from .bams."""
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.thview.apply(write_clusters, *args)
        if jobs:
            msg = "building clusters"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def run(self):
        """Run steps of the reference-based step3.

        Series of functions to run single or paired reference assembly
        with or without PCR decloning.
        """
        # creates several suffix files in same dir as reference(s)
        self.index_references()
        # creates _concat_R* files if multiple present.
        self.concat_trimmed_files_from_assembly_merge()
        # creates _unmapped_R* files with kept reads
        self.mapping_to_reference_filter()
        # creates _joined files with all pairs merged by 'nnnn'
        self.join_pairs_for_derep()
        # creates _decloned files IF decloning is turned on.
        self.decloning_transfer_tags_inline()
        # creates _derep from merged or decloned.
        self.dereplicate()
        # creates _derep_tag.fa from _derep IF decloning is turned on.
        self.decloning_transfer_tags_to_header()
        # creates _derep_split_R[1,2].fa ...
        self.split_derep_pairs_for_mapping()
        # creates .bam by mapping derep paired files to reference
        self.mapping_to_reference()
        # creates .clusters.gz files in stepdir
        self.build_clusters_from_cigars()
        # ...
        # self.declone_clusters()
        # ...
        self.calculate_sample_stats()


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    TEST = ip.load_json("/tmp/RICHIE.json")
    TEST.ipcluster['threads'] = 6
    TEST.run("34", force=True, quiet=False, cores=1)
    print(TEST.stats)

    # TEST = ip.load_json("/tmp/ama-denovo.json")
    # TEST = TEST.branch("ama-ref")
    # TEST.params.assembly_method = "reference"
    # TEST.params.reference_sequence = "/home/deren/Documents/ipyrad/sandbox/Ahypochondriacus_459_v2.0.fa"

    # with ip.Cluster(cores=2) as ipyclient:
    #     step = ip.assemble.s3_clustmap_within.Step3(TEST, 1, 0, ipyclient)
    #     c = ClustMapReference(step)
    #     c.index_references()
    #     c.concat_trimmed_files_from_assembly_merge()
    #     c.mapping_to_reference_filter()
    #     c.join_pairs_for_derep()


    # TEST.params.reference_sequence = "../../tests/ipsimdata/pairddrad_example_genome.fa"

    # practice with concat of some samples
    # DATA = ip.load_json("/home/deren/Documents/ipyrad/sandbox/ama3rad/assembly/amaranth.json")
    # DATA1 = DATA.branch("test-branch1", subsample=["SLH_AL_1000", "SLH_AL_1009"])
    # DATA2 = DATA.branch("test-branch2", subsample=["SLH_AL_1000", "SLH_AL_1009"])
    # DATA = ip.merge("merge", [DATA1, DATA2])

    # print(DATA.params)
    # DATA.params.reference_as_filter = "../../tests/ipsimdata/rad_example_genome.fa"
    # DATA.ipcluster['threads'] = 4
    # DATA.run("3", force=True, cores=8, quiet=True)
