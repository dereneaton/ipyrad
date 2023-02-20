#!/usr/bin/env python

"""Denovo clustering of reads within samples.

Identify reads from the same loci using clustering in vsearch.
Support for SE and PE denovo and denovo-reference assemblies,
and tagging PCR duplicates.
"""

import time
import itertools
from loguru import logger
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.clustmap_within_reference import ClustMapBase
from ipyrad.assemble.clustmap_within_reference_utils import index_ref_with_bwa
from ipyrad.assemble.clustmap_within_denovo_utils import (
    merge_pairs_with_vsearch,
    join_end_to_end,
    cluster,
    write_clusters,
    muscle_chunker,
    write_alignments,
    reconcat,
)

logger = logger.bind(name="ipyrad")


class ClustMapDenovo(ClustMapBase):
    """de novo within-sample assembly pipeline."""
    def __init__(self, step):
        super().__init__(step)
        self.data.max_indels = 5
        """: max number of gap-openings in a within-sample stack."""

        # sort samples to start largest files clustering first
        self.sorted_samples = sorted(
            self.samples, key=lambda x: self.samples[x].stats_s2.reads_raw)

    def index_reference_as_filter(self):
        """Index reference_filter with BWA for mapping reads.

        This can be used in denovo assemblies for filtering reads from
        the 'reference_as_filter' param.
        """
        jobs = {}
        if not self.data.params.reference_as_filter:
            return
        logger.debug("indexing reference_as_filter with bwa")
        rasync1 = self.lbview.apply(index_ref_with_bwa, self.data, alt=1)
        jobs['bwa_index_alt'] = rasync1
        msg = "indexing reference"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def pair_merge_overlaps_with_vsearch(self):
        """Merge reads based on overlapping using vsearch.

        This returns the number of merged pairs stored to stats.

        # i2: edits/{}_edits_R[1,2].fastq
        # i1: tmpdir/{}_concat_edits_R[1,2].fastq
        # i0: tmpdir/{}_unmapped_R[1,2].fastq
        # o: tmpdir/{}_merged.fa
        # o: tmpdir/{}_nonmerged_R[1,2].fa
        """
        # skip if not pairs
        if not self.data.is_pair:
            return
        logger.info("merging overlapping pairs with vsearch")
        jobs = {}
        for sname in self.sorted_samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(merge_pairs_with_vsearch, *args)
        msg = "merging overlapping paired reads"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

        # store n merged reads
        for sname, sample in self.samples.items():
            sample.stats_s3.merged_pairs = prog.results[sname]
            sample.stats_s3.merged_pairs_prop = (
                prog.results[sname] / sample.stats_s2.reads_raw)

    def pair_join_unmerged_end_to_end(self):
        """Joins end-to-end the unmerged paired reads.

        Concats to the end of this file any vsearch merged pe reads.

        # i: tmpdir/{}_merged.fa
        # i: tmpdir/{}_nonmerged_R[1,2].fa
        # o: tmpdir/{}_merged.fa
        """
        if not self.data.is_pair:
            return
        logger.info("joining unmerged paired reads for derep")
        jobs = {}
        for sname in self.sorted_samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(join_end_to_end, *args)
        msg = "joining non-overlapping pairs"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def cluster_build(self):
        """Cluster reads and build cluster output files.

        Submit clustering/mapping job. All jobs will start in order
        and the tracking progress bar will progress as each group
        finishes. These functions are grouped together so that each
        sample can be passed through without waiting on all samples.
        """
        # sort samples by size to cluster/align largest first
        start = time.time()
        casyncs = {}
        for sname in self.sorted_samples:
            args = (self.data, self.samples[sname])
            casyncs[sname] = self.thview.apply(cluster, *args)

        # submit cluster building job
        basyncs = {}
        for sname in self.sorted_samples:
            args = (self.data, self.samples[sname])
            with self.lbview.temp_flags(after=casyncs[sname]):
                basyncs[sname] = self.lbview.apply(write_clusters, *args)

        # get cluster data: --> /tmp/_utemp.tsv /tmp/htemp.fa
        msg = "clustering"
        prog1 = AssemblyProgressBar(casyncs, msg, 3, self.quiet, start)
        prog1.block()
        prog1.check()

        # build unaligned clusters: --> /tmp/_clust.txt
        start = time.time()
        msg = "building clusters"
        prog2 = AssemblyProgressBar(basyncs, msg, 3, self.quiet, start)
        prog2.block()
        prog2.check()

    def muscle_chunk(self):
        """Break unaligned cluster into chunks for aligning (/tmp/_[0-9].ali)

        This will find either .clusters.txt or .clusters_decloned.txt
        to use as input. Priority for second if present.
        """
        start = time.time()
        hasyncs = {}
        for sname in self.sorted_samples:
            args = (self.data, self.samples[sname])
            hasyncs[sname] = self.thview.apply(muscle_chunker, *args)
        msg = "chunking clusters"
        prog3 = AssemblyProgressBar(hasyncs, msg, 3, self.quiet, start)
        prog3.block()
        prog3.check()

    def muscle_align_chunks(self):
        """Aligns all chunked loci using muscle"""
        # submit largest samples first
        sorted_samples = sorted(
            self.samples,
            key=lambda x: self.samples[x].stats_s2.reads_raw,
            reverse=True,
        )

        # submit 10 aligning jobs per sample.
        jobs = {}
        for sname in sorted_samples:
            jobs[sname] = []
            for idx in range(10):
                handle = self.data.tmpdir / f"{sname}_chunk_{idx}.ali"
                # submit to be aligned if any data in this file.
                if handle.stat().st_size:
                    jobs[sname].append(
                        self.lbview.apply(write_alignments, handle))

        # a list with all aasyncs concatenated for use in progbar
        aasyncs = itertools.chain(*jobs.values())
        aasyncs = dict(enumerate(aasyncs))

        # submit cluster building job for each sample *after* all align jobs
        basyncs = {}
        for sname in sorted_samples:
            args = (self.data, self.samples[sname])
            with self.lbview.temp_flags(after=jobs[sname]):
                rasync = self.lbview.apply(reconcat, *args)
                basyncs[sname] = rasync

        # track job 1 progress
        msg = "aligning clusters"
        prog = AssemblyProgressBar(aasyncs, msg, 3, self.quiet)
        prog.block()
        prog.check()

        # track job 2 progress
        msg = "concat clusters"
        prog = AssemblyProgressBar(basyncs, msg, 3, self.quiet)
        prog.block()
        prog.check()

    def run(self):
        """Run the core functions."""
        self.index_reference_as_filter()
        self.mapping_to_reference_filter()
        self.concat_trimmed_files_from_assembly_merge()
        self.pair_merge_overlaps_with_vsearch()
        self.pair_join_unmerged_end_to_end()      # -> ...
        self.decloning_transfer_tags_inline()     # -> tmp/_decloned.fa
        self.dereplicate()                        # -> tmp/_derep.fa
        self.decloning_transfer_tags_to_header()  # -> tmp/_derep_tag.fa
        self.cluster_build()                      # -> tmp/.clusters.txt
        self.declone_clusters()                   # -> tmp/.clusters_decloned.txt
        self.muscle_chunk()                       # -> tmp/.ali
        self.muscle_align_chunks()                # -> tmp/.alignment
        self.calculate_sample_stats()


if __name__ == "__main__":


    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    # TEST = ip.load_json("/tmp/TEST5.json")
    # TEST = TEST.branch("TEST5-denovo")
    # TEST.params.assembly_method = "denovo"
    # TEST.params.reference_sequence = "../../tests/ipsimdata/pairddrad_example_genome.fa"

    # ClustMapDenovo(TEST)

    # TEST.run("3", force=True, quiet=True)

    # DATA = ip.load_json("/tmp/TEST1.json")
    # DATA.run('3', force=True)

    # STEP = ip.assemble.s3_clustmap_within.Step3(DATA, 1, 0, CLIENT)
    # STEP.samples['1A_0'].concat = "A"
    # TOOL = ClustMapDenovo(STEP)
    # TOOL.run()
#