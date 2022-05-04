#!/usr/bin/env python

"""Call consensus alleles from clustered stacks within samples.

"""

from typing import Dict, TypeVar
from collections import Counter
import numpy as np
import pandas as pd
from loguru import logger

from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.schema import Stats5
from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.clustmap_within_denovo_utils import iter_clusters
from ipyrad.assemble.s4_joint_estimate import recal_hidepth_cluster_stats
from ipyrad.assemble.consens_utils import (
    Processor,
    concat_catgs,
    concat_denovo_consens,
    concat_reference_consens,
)

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")
logger = logger.bind(name="ipyrad")


class Step5(BaseStep):
    """Run Step3 clustering/mapping using vsearch or bwa"""
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 5, quiet, force)
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.data.max_frag = self.data.hackers.max_fragment_length
        self.keep_masks: Dict[str, np.ndarray] = {}

        # sample paths to be used
        for sname, sample in self.samples.items():
            sample.files.depths = str(self.data.stepdir / f"{sname}.catg.hdf5")
            if self.data.is_ref:
                sample.files.consens = str(self.data.stepdir / f"{sname}.bam")
            else:
                sample.files.consens = str(self.data.stepdir / f"{sname}.consens.gz")

    def run(self):
        """Submit jobs to run either denovo, reference, or complex."""
        self.calculate_depths_and_max_frag()
        self.make_chunks()
        self.process_chunks()
        self.concatenate_chunks()
        self.data.save_json()

    def calculate_depths_and_max_frag(self):
        """Check whether mindepth has changed and calc nclusters and maxlen
        """
        # send jobs to be processed on engines
        jobs = {}
        for sname, sample in self.samples.items():
            kwargs = dict(data=self.data, sample=sample, majrule=True)
            rasync = self.lbview.apply(recal_hidepth_cluster_stats, **kwargs)
            jobs[sname] = rasync
        msg = "calculating depths"
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()

        # check for failures and collect results
        for sname, res in prog.results.items():
            keep_mask, max_frag = res
            self.data.max_frag = int(max(self.data.max_frag, max_frag))
            self.keep_masks[sname] = keep_mask
            logger.debug(f"high depth clusters in {sname}: {keep_mask.sum()}")
        logger.debug(
            f"max_fragment_length will be constrained to {self.data.max_frag}")

    def make_chunks(self):
        """Split clusters into chunks for parallel processing
        """
        msg = "chunking clusters"
        jobs = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample, self.keep_masks[sname])
            rasync = self.lbview.apply(make_chunk_files, *args)
            jobs[sname] = rasync
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()

    def process_chunks(self):
        """Process the cluster chunks into arrays and consens or bam files
        """
        def processor_wrap(data, sample, chunkfile):
            """Send job to Processor class on remote engine."""
            proc = Processor(data, sample, chunkfile)
            proc.run()
            return proc.counters, proc.filters

        # get mean values across all samples
        self.data.error_est = np.mean([
            self.data.samples[i].stats_s4.error_est
            for i in self.data.samples
        ])
        self.data.hetero_est = np.mean([
            self.data.samples[i].stats_s4.error_est
            for i in self.data.samples
        ])

        # submit jobs (can be many per sample)
        msg = "consens calling"
        prog = AssemblyProgressBar({}, msg, 5, self.quiet)
        for sname, sample in self.samples.items():
            chunks = list(self.data.tmpdir.glob(f"{sname}_chunk_[0-9]*"))
            chunks.sort(key=lambda x: int(x.name.split('_')[-1]))
            for cidx, chunk in enumerate(chunks):
                args = (self.data, sample, chunk)
                jid = (sname, cidx)
                prog.jobs[jid] = self.lbview.apply(processor_wrap, *args)

        # send chunks to be processed
        prog.block()
        prog.check()

        # collect results from jobs since where they're grouped by sname
        stats = {i: [Counter(), Counter()] for i in self.data.samples}
        for key, result in prog.results.items():
            sname, _ = key
            stats[sname][0].update(result[0])
            stats[sname][1].update(result[1])

        # ONLY ADVANCE STATE IF CONSENS READS EXIST
        for sname, (counters, filters) in stats.items():
            sample = self.samples[sname]
            if counters['nconsens']:
                sample.state = 5

            prefiltered_by_depth = np.invert(self.keep_masks[sname]).sum()

            # STORE STATS
            sample.stats_s5 = Stats5(
                clusters_total=self.keep_masks[sname].size,
                consensus_total=counters["nconsens"],
                nsites=counters["nsites"],
                nhetero=counters["nheteros"],
                heterozygosity = (0. if counters['nsites'] == 0
                    else counters['nheteros'] / counters['nsites']),
                filtered_by_depth=filters["depth"] + prefiltered_by_depth,
                filtered_by_max_h=filters["maxh"],
                filtered_by_max_n=filters["maxn"],
                filtered_by_max_alleles=filters["maxalleles"],
                min_depth_maj_during_step5=self.data.params.min_depth_majrule,
                min_depth_stat_during_step5=self.data.params.min_depth_statistical,
            )

        # WRITE TO STATS FILE and LOGGER
        stats_file = self.data.stepdir / "s5_stats_consensus.txt"
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=list(self.samples[sname].stats_s5.dict()),
        )
        for sname, sample in self.samples.items():
            istats = sample.stats_s5.dict()
            for column in statsdf.columns:
                statsdf.loc[sname, column] = istats[column]
        with open(stats_file, 'w', encoding="utf-8") as out:
            out.write(statsdf.to_string())
            logger.info("\n" + statsdf.iloc[:, :7].to_string())

    def concatenate_chunks(self):
        """Concatenate processed chunks into final stepdir.

        Relabels chunks by new consensus IDs. This spends most of its
        time storing CATG data that will probably not be used,
        but is important for saving SNP depths info.
        """
        # concat catgs for each sample
        jobs1 = {}
        jobs2 = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample)

            # submit h5 counts concat
            jobs1[sname] = self.lbview.apply(concat_catgs, *args)

            # submit seq file concat
            if not self.data.is_ref:
                jobs2[sname] = self.lbview.apply(concat_denovo_consens, *args)
            else:
                jobs2[sname] = self.lbview.apply(concat_reference_consens, *args)

        jobs = list(jobs1.values()) + list(jobs2.values())
        jobs = dict(enumerate(jobs))
        msg = "indexing alleles"
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()


def make_chunk_files(data: Assembly, sample: Sample, keep_mask: np.ndarray) -> None:
    """Split job into 5_000 hidepth clusters each.
    """
    # open to cluster generator
    clusters = iter_clusters(sample.files.clusters, gzipped=True)

    # load in high depth clusters and then write to chunk
    chunk = []
    sidx = 0
    for idx, clust in enumerate(clusters):
        if not keep_mask[idx]:
            continue
        chunk.append("".join(clust))

        # write to chunk file and reset        
        if len(chunk) == 5_000:
            end = sidx + len(chunk)
            handle = data.tmpdir / f"{sample.name}_chunk_{sidx}_{end}"
            with open(handle, 'w', encoding="utf-8") as out:
                out.write("//\n//\n".join(chunk) + "//\n//\n")
            chunk = []
            sidx += 5000

    # write any remaining
    if chunk:
        end = sidx + len(chunk)        
        handle = data.tmpdir / f"{sample.name}_chunk_{sidx}_{end}"
        with open(handle, 'w', encoding="utf-8") as out:
            out.write("//\n//\n".join(chunk) + "//\n//\n")


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    TEST = ip.load_json("/tmp/TEST5.json")
    TEST.run("5", force=True, quiet=False)
