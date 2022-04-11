#!/usr/bin/env python

"""Denovo clustering of reads within samples.

Identify reads from the same loci using clustering in vsearch. 
Support for SE and PE denovo and denovo-reference assemblies, 
and tagging PCR duplicates.
"""

import time
import itertools
from loguru import logger
import pandas as pd
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.schema import Stats3
from ipyrad.assemble.clustmap_within_reference_utils import index_ref_with_bwa 
from ipyrad.assemble.clustmap_within_denovo_utils import (
    concat_multiple_edits,
    mapping_reads_minus,
    merge_pairs_with_vsearch,
    join_end_to_end,
    tag_seq_for_decloning,
    dereplicate,
    retag_header_after_derep_for_decloning,
    cluster,
    write_clusters,
    muscle_chunker,
    write_alignments,
    reconcat,
    set_sample_stats
)


class ClustMapDenovo:
    """de novo within-sample assembly pipeline."""
    def __init__(self, step):
        # attach all relevent attributes to data (Assembly)
        self.data = step.data
        self.data.max_indels = 8

        # job submitting, parallel, or progress bar relevant
        self.samples = step.samples
        self.quiet = step.quiet
        self.lbview = step.lbview
        self.thview = step.thview

        # set empty stats_s3 on each Sample
        for sample in self.samples.values():
            sample.stats_s3 = Stats3(
                min_depth_maj_during_step3=self.data.params.min_depth_majrule,
                min_depth_stat_during_step3=self.data.params.min_depth_statistical,
            )

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

    def concat_trimmed_files_from_assembly_merge(self):
        """Concat when assemblies were merged before step3.

        # i: edits/{sname}_edits_[R1,R2].fastq
        # o: tmpdir/{sname}_concat_edits_[R1,R2].fastq
        """
        # skip if no merge happened.
        if not any(len(self.samples[i].files.edits) > 1 for i in self.samples):
            return
        # else run the job
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            if len(self.samples[sname].files.edits) > 1:
                jobs[sname] = self.lbview.apply(concat_multiple_edits, *args)
        if jobs:
            msg = "concatenating merged assembly inputs"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def mapping_to_reference_filter(self):
        """Map reads to filter reference to get unmapped fastq.

        # i1: edits/{}_edits_R[1,2].fastq
        # i0: tmpdir/{}_concat_edits_R[1,2].fastq
        # o: tmpdir/{}_unmapped_R[1,2].fastq
        """
        if not self.data.params.reference_as_filter:
            return        
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.thview.apply(mapping_reads_minus, *args)
        if jobs:
            msg = "mapping to reference as filter"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

        # store stats
        for sname, val in prog.results.items():
            sample = self.data.samples[sname]
            sample.stats_s3.reads_mapped_to_ref_filter = val[0]
            sample.stats_s3.reads_mapped_to_ref_filter_prop = val[1]
            logger.debug(
                f"{sname} proportion refmap filtered reads: {val[1]:.4f}")            

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
        for sname in self.samples:
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
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(join_end_to_end, *args)
        msg = "joining non-overlapping pairs"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def decloning_transfer_tags_inline(self):
        """Moves the i5 tags from the index to start of reads for decloning.

        # i: tmpdir/{}_merged.fa
        # o: tmpdir/{}_decloned.fa
        """
        if (not self.data.is_pair) or (not self.data.hackers.declone_PCR_duplicates):
            return
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(tag_seq_for_decloning, *args)
        if jobs:
            msg = "tagging reads for decloning    "
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()            

    def dereplicate(self):
        """Dereplicate sequences (read pairs are merged).

        # i3: edits/{}_edits.fastq                 # se data
        # i2: tmpdir/{}_concat_edits.fastq         # se assembly merge
        # i1: tmpdir/{}_merged.fa                  # pe data
        # i0: tmpdir/{}_decloned.fa                # pe w/ declone
        # o: tmpdir/{}_derep.fa
        """
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.thview.apply(dereplicate, *args)
        msg = "dereplicating"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def decloning_retag_to_header(self):
        """Moves the i5 tags from the index to start of reads for decloning.

        # i: tmpdir/{}_derep.fa
        # o: tmpdir/{}_derep_tag.fa
        """
        if (not self.data.is_pair) or (not self.data.hackers.declone_PCR_duplicates):
            return
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(
                retag_header_after_derep_for_decloning, *args)
        if jobs:
            msg = "moving tags to header"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()            

    def cluster_build_and_chunk(self):
        """Cluster reads and build cluster output files.

        Submit clustering/mapping job. All jobs will start in order
        and the tracking progress bar will progress as each group 
        finishes. These functions are grouped together so that each
        sample can be passed through without waiting on all samples.
        """
        # sort samples by size to cluster/align largest first
        sorted_samples = sorted(
            self.samples, key=lambda x: self.samples[x].stats_s2.reads_raw)
        start = time.time()
        casyncs = {}
        for sname in sorted_samples:
            args = (self.data, self.samples[sname])
            casyncs[sname] = self.thview.apply(cluster, *args)

        # submit cluster building job
        basyncs = {}
        for sname in sorted_samples:
            args = (self.data, self.samples[sname])
            with self.lbview.temp_flags(after=casyncs[sname]):
                basyncs[sname] = self.lbview.apply(write_clusters, *args)

        # submit cluster chunking job
        hasyncs = {}
        for sname in sorted_samples:
            args = (self.data, self.samples[sname])
            with self.lbview.temp_flags(after=basyncs[sname]):
                hasyncs[sname] = self.lbview.apply(muscle_chunker, *args)

        # track job progress
        msg = "clustering"
        prog1 = AssemblyProgressBar(casyncs, msg, 3, self.quiet, start)
        prog1.block()
        prog1.check()

        # track job progress
        start = time.time()
        msg = "building clusters"
        prog2 = AssemblyProgressBar(basyncs, msg, 3, self.quiet, start)
        prog2.block()
        prog2.check()

        # track job progress
        start = time.time()
        msg = "chunking clusters"
        prog3 = AssemblyProgressBar(hasyncs, msg, 3, self.quiet, start)
        prog3.block()
        prog3.check()

        # store declone pcr clusters stats to samples
        for sname, val in prog2.results.items():
            pcr_dups, pcr_dups_prop = val
            sample = self.data.samples[sname]
            sample.stats_s3.pcr_duplicates = pcr_dups
            sample.stats_s3.pcr_duplicates_prop = pcr_dups_prop
            if self.data.hackers.declone_PCR_duplicates:
                logger.info(f"{sname} proportion pcr duplicates: {pcr_dups_prop}")

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
                jobs[sname].append(self.lbview.apply(write_alignments, handle))

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

    def calculate_sample_stats(self):
        """Send samples to calc depths on remote, and then enter stats
        to sample objects non-parallel.
        """
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(set_sample_stats, *args)

        # progress bar
        msg = "calculating stats"
        prog = AssemblyProgressBar(jobs, msg, 3, self.quiet)
        prog.block()
        prog.check()

        # store/save
        self.data.samples = prog.results
        self.data.save_json()

        # write stats text file
        handle = self.data.stepdir / "s3_cluster_stats.txt"
        with open(handle, 'w', encoding="utf-8") as out:
            table = pd.DataFrame({
                i: self.data.samples[i].stats_s3.dict() 
                for i in self.data.samples
            }).T
            table.sort_index(inplace=True)
            table.to_string(out)
            cols=["clusters_total", "clusters_hidepth", "mean_depth_stat"]
            logger.info("\n" + table.loc[:, cols].to_string())

    def run(self):
        """Run the core functions."""
        self.index_reference_as_filter()
        self.mapping_to_reference_filter()
        self.concat_trimmed_files_from_assembly_merge()
        self.pair_merge_overlaps_with_vsearch()
        self.pair_join_unmerged_end_to_end()
        self.decloning_transfer_tags_inline()
        self.dereplicate()
        self.decloning_retag_to_header()
        self.cluster_build_and_chunk()
        self.muscle_align_chunks()
        self.calculate_sample_stats()


if __name__ == "__main__":


    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    TEST = ip.load_json("/tmp/TEST5.json")
    TEST = TEST.branch("TEST5-denovo")
    TEST.params.assembly_method = "denovo"
    TEST.params.reference_sequence = "../../tests/ipsimdata/pairddrad_example_genome.fa"

    ClustMapDenovo(TEST)

    # TEST.run("3", force=True, quiet=True)

    # DATA = ip.load_json("/tmp/TEST1.json")
    # DATA.run('3', force=True)

    # STEP = ip.assemble.s3_clustmap_within.Step3(DATA, 1, 0, CLIENT)
    # STEP.samples['1A_0'].concat = "A"
    # TOOL = ClustMapDenovo(STEP)
    # TOOL.run()
# 