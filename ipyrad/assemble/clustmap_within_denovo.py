#!/usr/bin/env python

"""Denovo clustering of reads within samples.

Identify reads from the same loci using clustering in vsearch. 
Support for SE and PE denovo and denovo-reference assemblies, 
and tagging PCR duplicates.
"""

import os
import time
import itertools
from loguru import logger
import pandas as pd
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.clustmap_within_reference_utils import (
    index_ref_with_bwa, 
    index_ref_with_sam,
)
from ipyrad.assemble.clustmap_within_denovo_utils import (
    concat_multiple_edits,
    merge_pairs_with_vsearch,
    merge_end_to_end,
    dereplicate,
    cluster,
    build_clusters,
    muscle_chunker,
    align_and_parse,
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

    def index_references(self):
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

        i: edits/{sname}_edits_[R1,R2].fastq
        o: tmpdir/{sname}_concat_edits_[R1,R2].fastq
        """
        # skip if no merge happened.
        if not any(len(self.samples[i].files.edits) > 1 for i in self.samples):
            return
        # else run the job
        jobs = {}
        for sname in self.samples:
            if len(self.samples[sname].files.edits) > 1:
                jobs[sname] = self.lbview.apply(
                    concat_multiple_edits,
                    *(self.data, self.samples[sname])
                )
        if jobs:
            msg = "concatenating merged assembly inputs"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def pair_merge_with_vsearch(self):
        """Merge reads based on overlapping using vsearch.

        i: tmpdir/{sname}_concat_edits_[R1,R2].fastq OR edits
        o: tmpdir/{sname}_merged.fastq, tmpdir/{sname}_nonmerged_R[1,2].fastq, 
        """
        # skip if not pairs
        if not self.data.is_pair:
            return
        logger.info("merging paired reads with vsearch")
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                merge_pairs_with_vsearch,
                *(self.data, self.samples[sname])
            )
        msg = "merging overlapping paired reads"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def pair_join_unmerged_end_to_end(self):
        """Joins end-to-end the unmerged paired reads, and concats to the
        end of this file any merged reads from step2.

        i: tmpdir/{sname}_nonmerged_[R1,R2].fastq
        o: tmpdir/{sname}_merged.fastq
        """
        if not self.data.is_pair:
            return
        logger.info("joining unmerged paired reads")            
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                merge_end_to_end, 
                *(self.data, self.samples[sname], True, True)
            )
        msg = "joining non-overlapping pairs"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def decloning_transfer_tags_inline(self):
        """Moves the tags from the index to the reads for decloning.

        TODO: run on test data.
        """
        if (not self.data.is_pair) or (not self.data.hackers.declone_PCR_duplicates):
            return
        raise NotImplementedError("TODO.")

    def dereplicate(self):
        """Dereplicate sequences (read pairs are merged).

        i: tmpdir/{}_[merged,declone].fastq
        o: tmpdir/{}_derep.fa
        """
        jobs = {}
        for sname in self.samples:
            # debugging
            # dereplicate(self.step, self.samples[sname])
            jobs[sname] = self.thview.apply(
                dereplicate, 
                *(self.data, self.samples[sname]))
        msg = "dereplicating"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()


    def decloning_transfer_tags_to_header(self):
        """

        """


    def optional_pair_split_joined_pairs(self):
        """

        """


    def optional_map_to_reference_filter(self):
        """

        """


    def optional_merge_back_after_refmapping(self):
        """

        """


    def cluster_build_and_chunk(self):
        """Cluster reads and build cluster output files.

        Submit clustering/mapping job. All jobs will start in order
        and the tracking progress bar will progress as each group 
        finishes.
        """
        # sort samples by size to cluster/align largest first
        sorted_samples = sorted(
            self.samples, key=lambda x: self.samples[x].stats_s2.reads_raw)
        start = time.time()
        casyncs = {}
        for sname in sorted_samples:
            casyncs[sname] = self.thview.apply(
                cluster, *(self.data, self.samples[sname])
            )

        # submit cluster building job
        basyncs = {}
        for sname in sorted_samples:
            with self.lbview.temp_flags(after=casyncs[sname]):
                basyncs[sname] = self.lbview.apply(
                    build_clusters,
                    *(self.data, self.samples[sname])
                )

        # submit cluster chunking job
        hasyncs = {}
        for sname in sorted_samples:
            with self.lbview.temp_flags(after=basyncs[sname]):
                hasyncs[sname] = self.lbview.apply(
                    muscle_chunker,
                    *(self.data, self.samples[sname])
                )

        # track job progress
        msg = "clustering"
        prog = AssemblyProgressBar(casyncs, msg, 3, self.quiet, start)
        prog.block()
        prog.check()

        # track job progress
        start = time.time()
        msg = "building clusters"
        prog = AssemblyProgressBar(basyncs, msg, 3, self.quiet, start)
        prog.block()
        prog.check()

        # track job progress
        start = time.time()
        msg = "chunking clusters"
        prog = AssemblyProgressBar(hasyncs, msg, 3, self.quiet, start)
        prog.block()
        prog.check()        


    def muscle_align_chunks(self):
        """Aligns all chunked loci using muscle"""
        # submit largest samples first
        sorted_samples = sorted(
            self.samples, key=lambda x: self.samples[x].stats_s2.reads_raw)

        # submit 10 aligning jobs per sample.
        jobs = {}
        for sname in sorted_samples:
            jobs[sname] = []
            for idx in range(10):            
                handle = os.path.join(
                    self.data.tmpdir, f"{sname}_chunk_{idx}.ali")
                args = (
                    handle, 
                    self.data.max_indels, 
                    "gbs" in self.data.params.datatype,
                    self.data.hackers.declone_PCR_duplicates
                )
                jobs[sname].append(self.lbview.apply(align_and_parse, *args))

        # a list with all aasyncs concatenated for use in progbar
        aasyncs = itertools.chain(*jobs.values())
        aasyncs = dict(enumerate(aasyncs))

        # submit cluster building job for each sample *after* all align jobs
        basyncs = {}
        for sname in sorted_samples:
            with self.lbview.temp_flags(after=jobs[sname]):
                rasync = self.lbview.apply(
                    reconcat, 
                    *(self.data, self.samples[sname])
                )
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
            jobs[sname] = self.lbview.apply(
                set_sample_stats, 
                *(self.data, self.samples[sname])
            )

        # progress bar
        msg = "calculating stats"
        prog = AssemblyProgressBar(jobs, msg, 3, self.quiet)
        prog.block()
        prog.check()

        # store/save
        self.data.samples = {i: prog.results[i] for i in prog.results}
        self.data.save_json()

        # write stats text file
        statsdf = pd.DataFrame(
            index=sorted(self.data.samples),
            columns=["clusters_total", "clusters_hidepth", "mean_depth_mj"],
        )
        for sname in self.data.samples:
            statsdict = self.data.samples[sname].stats_s3.dict()
            for i in statsdf.columns:
                statsdf.loc[sname, i] = statsdict[i]                    
        logger.info("\n" + statsdf.to_string())
 

    def run(self):
        """Run the core functions."""
        self.index_references()
        # self.map_to_reference_filter()
        self.concat_trimmed_files_from_assembly_merge()
        self.pair_merge_with_vsearch()
        self.pair_join_unmerged_end_to_end()
        self.decloning_transfer_tags_inline()
        self.dereplicate()
        self.decloning_transfer_tags_to_header()
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