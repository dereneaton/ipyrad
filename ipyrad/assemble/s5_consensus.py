#!/usr/bin/env python

"""Call consensus alleles from clustered stacks within samples.

Example
-------
>>> # load data as `consens_utils.ConsensusProcessor` object
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

import ipyparallel as ipp
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.schema import Stats5
from ipyrad.assemble.utils import NoHighDepthClustersError
from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.s4_joint_estimate import recal_hidepth_cluster_stats, make_chunk_files
from ipyrad.assemble.consens_utils import (
    ConsensusProcessor,
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

    def debug(self, sname: str) -> ConsensusProcessor:
        """Function for developers only.

        Used to debug and check consensus calling algorithms. This
        returns a ConsensusProcessor object for an unchunked file for a sample.
        """
        self.calculate_depths_and_max_frag()
        self.make_chunks(chunksize=int(1e9))
        self.set_s4_params()
        tmpfile = list(self.data.tmpdir.glob(f"{sname}_chunk_*"))[0]
        sample = self.samples[sname]
        return ConsensusProcessor(self.data, sample, tmpfile)

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

        # check results. If no high depth clusts then catch case
        for sname in prog.jobs:
            try:
                kmask, max_frag = prog.jobs[sname].get()
                self.keep_masks[sname] = kmask
                self.data.max_frag = int(max(self.data.max_frag, max_frag))
                logger.debug(f"hidepth clusts in {sname}: {kmask.sum()}")
            except ipp.error.RemoteError:
                exc = prog.jobs[sname].exception()
                if exc.ename == "NoHighDepthClustersError":
                    self.keep_masks[sname] = np.zeros(
                        self.data.samples[sname].stats_s3.clusters_total,
                        dtype=np.bool_,
                    )
                    logger.warning(f"No high depth clusters in {sname}.")
                else:
                    raise
        logger.debug(
            f"max_fragment_length will be constrained to {self.data.max_frag}")

    def make_chunks(self, chunksize: int = 5000):
        """Split clusters into chunks for parallel processing
        """
        msg = "chunking clusters"
        jobs = {}
        for sname, sample in self.samples.items():
            args = dict(data=self.data, sample=sample, keep_mask=self.keep_masks[sname], chunksize=chunksize)
            rasync = self.lbview.apply(make_chunk_files, **args)
            jobs[sname] = rasync
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()

    def set_s4_params(self):
        """Set values from step4 results to the data object."""
        # get mean values across all samples
        self.data.error_est = np.nanmean([
            self.samples[i].stats_s4.error_est for i in self.samples
        ])
        self.data.hetero_est = np.nanmean([
            self.samples[i].stats_s4.error_est for i in self.samples
        ])

    def process_chunks(self):
        """Process the cluster chunks into arrays and consens or bam files
        """
        def processor_wrap(data, sample, chunkfile):
            """Send job to Processor class on remote engine."""
            proc = ConsensusProcessor(data, sample, chunkfile)
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
        stats = {i: [Counter(), Counter()] for i in self.samples}
        for key, result in prog.results.items():
            sname, _ = key
            stats[sname][0].update(result[0])
            stats[sname][1].update(result[1])

        # store stats to the samples, but don't advance the state yet.
        for sname, (counters, filters) in stats.items():
            sample = self.samples[sname]

            # successful consens reads
            prefiltered_by_depth = np.invert(self.keep_masks[sname]).sum()

            # STORE STATS
            sample.stats_s5 = Stats5(
                clusters_total=self.keep_masks[sname].size,
                consensus_total=counters["nconsens"],
                nsites=counters["nsites"],
                nhetero=counters["nheteros"],
                heterozygosity=(
                    0. if counters['nsites'] == 0
                    else counters['nheteros'] / counters['nsites']
                ),
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

        # iterate over samples
        for sname, sample in self.samples.items():
            istats = sample.stats_s5.dict()
            for column in statsdf.columns:
                statsdf.loc[sname, column] = istats[column]

        # print stats to logger
        with open(stats_file, 'w', encoding="utf-8") as out:
            out.write(statsdf.to_string())
            logger.info("\n" + statsdf.iloc[:, :7].to_string())

        # advance sample states
        for sample in self.samples.values():
            if sample.stats_s5.consensus_total:
                sample.state = 5
                sample._clear_old_results()
            else:
                sample.state = 4

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

            # skip this sample if no consensus data was found.
            if not sample.stats_s5.consensus_total:
                continue

            # submit h5 counts concat
            args = (self.data, sample)
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


if __name__ == "__main__":

    import ipyrad as ip
    #ip.set_log_level("DEBUG")#, log_file="/tmp/test.log")

    # testing in this region fails b/c some function are in this module
    # instead of another. Fix ithis..

    TEST = ip.load_json("/tmp/RICHIE.json")
    # TEST.ipcluster['threads'] = 2
    # TEST.run("5", force=True, quiet=True, cores=4)
    with ip.Cluster(cores=4) as ipyclient:
        step = Step5(TEST, force=True, quiet=False, ipyclient=ipyclient)
        step.run()
        # step.calculate_depths_and_max_frag()
        # step.set_s4_params()
        # make_chunk_files(step.data, step.samples['RA-225'], step.keep_masks["RA-225"], 1000)
        # step.make_chunks()  # chunksize=int(1e9))
        # cons = step.debug("RA-225")
        # print(cons)
    # print(TEST.stats)

    # for JSON in ["/tmp/TEST1.json", "/tmp/TEST5.json"]:
        # TEST = ip.load_json(JSON)
        # TEST.run("5", force=True, quiet=True)

    # TEST = ip.load_json("../../pedtest/NEW.json")
    # TEST.params.min_depth_majrule = 1
    # TEST.run("5", force=True, quiet=True)
    # print(TEST.stats)

    # TEST = ip.load_json("/tmp/TEST5.json")
    # TEST.run("5", force=True, quiet=False)
