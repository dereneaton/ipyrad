#!/usr/bin/env python

"""

"""

from loguru import logger
import pandas as pd
from ipyrad.core import BaseStep, track_remote_jobs
from ipyrad.within_clust_build.clustbuild_w_funcs import (
    write_tmp_clusters,
    write_alignments,
    reconcat,
)

logger = logger.bind(name="ipyrad")


class Step3(BaseStep):
    def __init__(self, data, force, ipyclient):
        super().__init__(data, step=3, force=force)
        self.ipyclient = ipyclient
        self._stats = {}

    def run(self):
        self._remote_write_muscle_chunk_files()
        self._remote_write_alignments()
        self._remote_concat()
        self._write_json_file()
        self._write_stats_file()

    def _remote_write_muscle_chunk_files(self):
        """Write chunks of X unaligned clusters to /tmp/sname_unaligned.fa"""
        logger.info("building clusters")
        rasyncs = {}
        lbview = self.ipyclient.load_balanced_view()
        func = write_tmp_clusters
        for sname, sample in self.samples.items():
            args = (self.data.tmpdir, sample)
            rasyncs[sname] = lbview.apply(func, *args)
        track_remote_jobs(rasyncs, self.ipyclient)

    def _remote_write_alignments(self):
        """Write chunks of 5K aligned clusters to /tmp/sname_aligned.fa"""
        logger.info("aligning clusters")
        rasyncs = {}
        # paired data will use 2 concurrent subprocesses per engine.
        threads = 2 if self.data.is_pair else 1
        lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::threads])
        func = write_alignments
        for sname, sample in self.samples.items():
            unaligned = sorted(self.data.tmpdir.glob(f"{sname}_unaligned_*.fa"))
            for ufile in unaligned:
                rasyncs[ufile.name] = lbview.apply(func, ufile)
        logger.debug(f"queued {len(rasyncs)} tmpfiles for aligning")
        track_remote_jobs(rasyncs, self.ipyclient)

    def _remote_concat(self):
        """Write chunks of 5K aligned clusters to /tmp/sname_aligned.fa"""
        logger.info("concatenating tmp files and calculating statistics")
        rasyncs = {}
        lbview = self.ipyclient.load_balanced_view()
        func = reconcat
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            rasyncs[sname] = lbview.apply(func, *args)
        track_remote_jobs(rasyncs, self.ipyclient)
        self._stats = {i: j.get() for i, j in rasyncs.items()}

    def _write_json_file(self):
        # store results to samples
        for sname, sample in self.samples.items():
            nclusters, mean, median, std = self._stats[sname]
            sample.stats_s3.clusters = nclusters
            sample.stats_s3.cluster_depth_mean = mean
            sample.stats_s3.cluster_depth_median = median
            sample.stats_s3.cluster_depth_std = std
            sample.stats_s3.cluster_depth_histogram = []
            sample.stats_s3.filtered_bad_alignment = 0
            if nclusters:
                sample.state = 3

        # update the Assembly object
        for sname, sample in self.samples.items():
            self.data.samples[sname] = sample
        self.data.save_json()

    def _write_stats_file(self):
        # write stats files
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=[
                "clusters",
                "cluster_depth_mean",
                "cluster_depth_median",
                "cluster_depth_std",
                # "cluster_depth_histogram",
                "filtered_bad_alignment",
            ],
        )
        for sname, sample in self.samples.items():
            statsdict = sample.stats_s3.model_dump()
            for i in statsdf.columns:
                if statsdict[i]:
                    statsdf.loc[sname, i] = statsdict[i]
        handle = self.data.stepdir / "s3_cluster_stats"
        # self.data.stats_files.s3 = ...
        with open(handle, 'w', encoding="utf-8") as outfile:
            statsdf.fillna(value=0).to_string(outfile)
        logger.info(f"step 3 results written to {handle}")


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    # data = load_json("/tmp/pairgbs_merge.json")
    # sample = data.samples["1A_0"]
    data = ip.load_json("/tmp/pedtest/half-demuxed.json")
    with ip.Cluster(cores=6) as ipyclient:
        tool = Step3(data, True, ipyclient)
        tool.run()
