#!/usr/bin/env python

"""Call consensus alleles from clustered stacks within samples.

Example
-------
>>> # load data as `consens_utils.Processor` object 
>>> data = ip.load_json(JSON_FILE)
>>> with ip.Cluster(cores=2) as ipyclient:
>>>     step = Step5(data, 1, 0, ipyclient)
>>>     proc = step.debug(SAMPLE_NAME)
>>>
>>> # interrogate within generator funcs
>>> for loc in proc._iter_filter_heteros():
>>>     print(loc.cidx, loc.filters)
>>>
>>> # or, run call filter generators 
>>> proc.collect_data()
"""

from typing import Dict, TypeVar
from collections import Counter
import numpy as np
import pandas as pd
from loguru import logger

from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.schema import Stats5
from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.clustmap_within_both import iter_clusters
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
        self.set_s4_params()
        self.process_chunks()
        self.concatenate_chunks()
        self.store_stats()
        self.data.save_json()

    def debug(self, sname: str) -> Processor:
        """Function for developers only.

        Used to debug and check consensus calling algorithms. This
        returns a Processor object for an unchunked file for a sample.
        """
        self.calculate_depths_and_max_frag()
        self.set_s4_params()
        self.make_chunks(chunksize=int(1e9))
        tmpfile = list(self.data.tmpdir.glob(f"{sname}_chunk_*"))[0]
        sample = self.data.samples[sname]
        return Processor(self.data, sample, tmpfile)

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

    def make_chunks(self, chunksize: int=5000):
        """Split clusters into chunks for parallel processing
        """
        msg = "chunking clusters"
        jobs = {}
        for sname, sample in self.samples.items():
            args = (self.data, sample, self.keep_masks[sname], chunksize)
            rasync = self.lbview.apply(make_chunk_files, *args)
            jobs[sname] = rasync
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()

    def set_s4_params(self):
        """Set values from step4 results to the data object."""
        # get mean values across all samples
        self.data.error_est = np.mean([
            self.data.samples[i].stats_s4.error_est
            for i in self.data.samples
        ])
        self.data.hetero_est = np.mean([
            self.data.samples[i].stats_s4.error_est
            for i in self.data.samples
        ])

    def process_chunks(self):
        """Process the cluster chunks into arrays and consens or bam files
        """
        def processor_wrap(data, sample, chunkfile):
            """Send job to Processor class on remote engine."""
            proc = Processor(data, sample, chunkfile)
            proc.run()
            return proc.counters, proc.filters

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

        # return stats in a dict. This will be stored to the objects
        # later after we complete the concatenation steps.
        stats = {i: [Counter(), Counter()] for i in self.data.samples}
        for key, result in prog.results.items():
            sname, _ = key
            stats[sname][0].update(result[0])
            stats[sname][1].update(result[1])

        # store stats to the samples, but don't advance the state yet.
        for sname, (counters, filters) in stats.items():
            sample = self.samples[sname]

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

    def store_stats(self):
        """Collect stats and write to Samples and statsfile."""
        # collect results from jobs since where they're grouped by sname

        # WRITE TO STATS FILE and LOGGER
        stats_file = self.data.stepdir / "s5_stats_consensus.txt"
        snames = sorted(self.samples)
        statsdf = pd.DataFrame(
            index=snames,
            columns=list(self.samples[snames[0]].stats_s5.dict()),
        )
        for sname, sample in self.samples.items():
            istats = sample.stats_s5.dict()
            for column in statsdf.columns:
                statsdf.loc[sname, column] = istats[column]
        with open(stats_file, 'w', encoding="utf-8") as out:
            out.write(statsdf.to_string())
            logger.info("\n" + statsdf.iloc[:, :7].to_string())
        # advance sample states
        for sample in self.samples.values():
            if sample.stats_s5.consensus_total:
                sample.state = 5

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


def make_chunk_files(
    data: Assembly, sample: Sample, keep_mask: np.ndarray, chunksize: int=5000) -> None:
    """Split job into 5_000 hidepth clusters each.

    Parameters
    ----------
    data: Assembly
        params and samples
    sample: Sample
        filepaths and stats
    keep_mask: np.ndarray
        Used to filter which consens reads will be included in chunk
        files used for parallel processing.
    chunksize: int
        Chunksize used for breaking the problem into parallel chunks.
        The chunks are unzipped. Default is 5K, but is set to very
        large when using debug() function to not split files.
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
        if len(chunk) == int(chunksize):
            end = sidx + len(chunk)
            handle = data.tmpdir / f"{sample.name}_chunk_{sidx}_{end}"
            with open(handle, 'w', encoding="utf-8") as out:
                out.write("//\n//\n".join(chunk) + "//\n//\n")
            chunk = []
            sidx += int(chunksize)

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
