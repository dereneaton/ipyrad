#!/usr/bin/env python

"""
Reference based assembly methods
"""

import os
from loguru import logger
import pandas as pd
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.clustmap_within_reference_utils import (
    index_ref_with_bwa, 
    index_ref_with_sam,
    mapping_reads_minus,
    join_pairs_for_derep_ref,
    tag_for_decloning,
    dereplicate_func,
    tag_to_header_for_decloning,
    split_derep_pairs_ref,
    mapping_reads,
    build_clusters_from_cigars,
)
from ipyrad.assemble.clustmap_within_denovo_utils import (
    concat_multiple_edits,
    # reconcat,
    set_sample_stats
)


class ClustMapReference:
    """
    de novo assembly pipeline.
    """
    def __init__(self, step):
        # inherit attrs from step
        self.data = step.data
        self.samples = step.samples
        self.data.tmpdir = step.tmpdir
        self.data.stepdir = step.stepdir

        self.quiet = step.quiet
        self.lbview = step.lbview
        self.thview = step.thview


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


    def concat_trimmed_files_from_assembly_merge(self):
        """
        If assembly merging was performed after step2 and before step3
        then multiple 'trim' fastq files must be merged within samples
        to be used as inputs for step3. 

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
                logger.debug(f"concatenating merged inputs: {sname}")
                jobs[sname] = self.lbview.apply(
                    concat_multiple_edits,
                    *(self.data, self.samples[sname])
                )
        if jobs:
            msg = "concatenating merged assembly inputs"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()


    def mapping_to_refminus(self):
        """
        map reads to altref to get unmapped fastq files.

        # i: tmpdir/{}_R[1,2]-tmp.fa
        # o: tmpdir/{}-tmp-umap[1,2].FASTQ
        """
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.thview.apply(
                mapping_reads_minus,
                *(self.data, self.samples[sname], self.data.ipcluster['threads'])
            )
        if jobs:
            msg = "joining pairs for dereplication"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()        


    def join_pairs_for_derep(self):
        """
        Joins pairs end to end for decloning and dereplicating.
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


    def tag_for_decloning(self):
        """
        Joins pairs end to end for decloning and dereplicating.
        i: tmpdir/[trim, concat, or unmapped].fastq[.gz]
        o: tmpdir/joined.fastq
        """
        if not self.data.hackers.declone_PCR_duplicates:
            return
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                tag_for_decloning,
                *(self.data, self.samples[sname])
            )
        if jobs:
            msg = "tagging reads for decloning"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()            


    def dereplicate(self):
        """
        Dereplicate read(pairs) to fasta format using vsearch.
        i: tmpdir/{}_[trimmed_R1, joined, declone].fastq
        o: tmpdir/{}_derep.fa
        """
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                dereplicate_func,
                *(self.data, self.samples[sname])
            )
        if jobs:
            msg = "dereplicating reads"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()            
    

    def tag_to_header(self):
        """
        Dereplicate read(pairs) to fasta format using vsearch.
        i: tmpdir/{}_[trimmed_R1, joined, declone].fastq
        o: tmpdir/{}_derep.fa
        """
        if not self.data.hackers.declone_PCR_duplicates:
            return        
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.lbview.apply(
                tag_to_header_for_decloning,
                *(self.data, self.samples[sname])
            )
        if jobs:
            msg = "tagging headers for decloning"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()            
    

    def split_derep_pairs_for_mapping(self):
        """
        Dereplicate read(pairs) to fasta format using vsearch.
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
        """
        map reads to reference to get mapped bams.

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
        """
        Builds clusters.gz files from .bams 
        """
        jobs = {}
        for sname in self.samples:
            jobs[sname] = self.thview.apply(
                build_clusters_from_cigars,
                *(self.data, self.samples[sname]),
            )
        if jobs:
            msg = "building clusters"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()          


    def calculate_sample_stats(self):
        """
        Send samples to calc depths on remote, and then enter stats
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
            columns=[
                "clusters_total", 
                "clusters_hidepth", 
                "mean_depth_mj",
                "mean_depth_stat",
                "reads_mapped_to_ref_prop",
            ],
        )
        for sname in self.data.samples:
            statsdict = self.data.samples[sname].stats_s3.dict()
            for i in statsdf.columns:
                statsdf.loc[sname, i] = statsdict[i]

        handle = os.path.join(self.data.stepdir, 's3_cluster_stats.txt')
        with open(handle, 'w') as outfile:
            statsdf.to_string(outfile)
        logger.info("\n" + statsdf.to_string())


    def run(self):
        """
        Series of functions to run single or paired reference assembly
        with or without PCR decloning.
        """
        self.index_references()
        self.concat_trimmed_files_from_assembly_merge()
        self.mapping_to_refminus()
        self.join_pairs_for_derep()
        self.tag_for_decloning()
        self.dereplicate()
        self.tag_to_header()
        self.split_derep_pairs_for_mapping()
        self.mapping_to_reference()
        self.build_clusters_from_cigars()
        self.calculate_sample_stats()


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG", logfile="/tmp/test.log")
    DATA = ip.load_json("/home/deren/Documents/ipyrad/sandbox/ama3rad/assembly/amaranth.json")

    # practice with concat of some samples
    DATA1 = DATA.branch("test-branch1", subsample=["SLH_AL_1000", "SLH_AL_1009"])
    DATA2 = DATA.branch("test-branch2", subsample=["SLH_AL_1000", "SLH_AL_1009"])    
    DATA = ip.merge("merge", [DATA1, DATA2])
    
    # print(DATA.params)
    DATA.params.reference_as_filter = "../../tests/ipsimdata/rad_example_genome.fa"
    DATA.ipcluster['threads'] = 4
    DATA.run("3", force=True, cores=8, quiet=True)
