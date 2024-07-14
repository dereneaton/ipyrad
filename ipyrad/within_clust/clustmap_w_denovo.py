#!/usr/bin/env python

"""ClustMapDenovo child class of ClustMapBase for within-sample denovo.

"""

from typing import TypeVar
from loguru import logger
import pandas as pd
from ipyrad.within_clust.clustmap_w import ClustMapBase
from ipyrad.core import track_remote_jobs
from ipyrad.schema.sample_schema import Stats2
from ipyrad.within_clust.clustmap_w_funcs import index_ref_with_bwa
from ipyrad.within_clust.clustmap_w_denovo_funcs import (
    merge_pairs_with_vsearch,
    join_end_to_end,
    cluster,
)

logger = logger.bind(name="ipyrad")
Client = TypeVar("Client")
Assembly = TypeVar("Assembly")


class ClustMapDenovo(ClustMapBase):
    def __init__(self, data: Assembly, force: bool, ipyclient: Client):
        super().__init__(data, force, ipyclient)
        self._merged_pairs = {}
        self._merged_pairs_prop = {}

    def run(self):
        """Run the core functions."""
        self.remote_index_reference_as_filter()
        self.remote_map_to_reference_as_filter()
        self.remote_pairs_merge_overlaps_with_vsearch()
        self.remote_pairs_join_unmerged_end_to_end()       # -> ...
        # self.remote_decloning_transfer_tags_inline()     # -> tmp/_decloned.fa
        self.remote_dereplicate()                          # -> tmp/_derep.fa
        # self.remote_decloning_transfer_tags_to_header()  # -> tmp/_derep_tag.fa
        self.remote_cluster()                              # -> out/...
        self.write_json_file()
        self.write_stats_file()

    def remote_index_reference_as_filter(self):
        """Index reference_filter for bwa mapping."""
        if not self.data.params.reference_as_filter:
            return

        logger.info("indexing reference_as_filter with bwa")
        func = index_ref_with_bwa
        kwargs = dict(data=self.data, as_filter=True)
        rasyncs = {0: self.ipyclient[0].apply(func, **kwargs)}
        track_remote_jobs(rasyncs, self.ipyclient)

    def remote_pairs_merge_overlaps_with_vsearch(self):
        """Merge reads based on overlapping using vsearch.

        This returns the number of merged pairs stored to stats.

        # i2: trimmed/{}_trimmed_R[1,2].fastq.gz
        # i1: tmpdir/{}_concat_edits_R[1,2].fastq.gz
        # i0: tmpdir/{}_unmapped_R[1,2].fastq
        # o: tmpdir/{}_merged.fa
        # o: tmpdir/{}_nonmerged_R[1,2].fa
        """
        if not self.data.is_pair:
            for sname in self.samples:
                self._merged_pairs[sname] = 0
                self._merged_pairs_prop[sname] = 0
            return

        logger.info("merging overlapping paired reads")
        lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::2])
        rasyncs = {}
        func = merge_pairs_with_vsearch
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            rasyncs[sname] = lbview.apply(func, *args)
        results = track_remote_jobs(rasyncs, self.ipyclient)

        # store results to samples
        for sname, sample in self.samples.items():
            self._merged_pairs[sname] = results[sname]
            self._merged_pairs_prop[sname] = (
                results[sname] / sample.stats_s1.reads_raw)

    def remote_pairs_join_unmerged_end_to_end(self):
        """Joins end-to-end the unmerged paired reads.

        Concats to the end of this file any vsearch merged pe reads.

        # i: tmpdir/{}_merged.fa
        # i: tmpdir/{}_nonmerged_R[1,2].fa
        # o: tmpdir/{}_merged.fa
        """
        if not self.data.is_pair:
            return

        logger.info("joining unmerged paired reads for dereplication")
        lbview = self.ipyclient.load_balanced_view()
        rasyncs = {}
        func = join_end_to_end
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            rasyncs[sname] = lbview.apply(func, *args)
        track_remote_jobs(rasyncs, self.ipyclient)

    def remote_cluster(self):
        """Cluster reads to write cluster info to files."""
        logger.info("clustering within samples")
        threads = (
            self.data.ipcluster['threads'] if self.data.ipcluster['threads']
            else 4)
        lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::threads])
        rasyncs = {}
        func = cluster
        for sname, sample in self.samples.items():
            args = (self.data, sample)
            rasyncs[sname] = lbview.apply(func, *args)

        # debugging kbd
        try:
            track_remote_jobs(rasyncs, self.ipyclient)
        except KeyboardInterrupt:
            logger.warning("job interrupted by user...")

        # store handles for the outfiles
        for sname, sample in self.samples.items():
            derep = self.data.tmpdir / f"{sname}_derep.fa"
            derep = derep.rename(self.data.stepdir / f"{sname}_derep.fa")
            matches = self.data.stepdir / f"{sname}_matches.tsv"
            sample.files.clustmap = (derep, matches)

    def write_json_file(self):
        """Write a stats summary file."""
        for sname, sample in self.samples.items():
            # if derep file size is 0 then no clusters exist
            if not sample.files.clustmap[0].stat().st_size:
                logger.warning(f"sample {sname} contains no clusters.")
                continue

            # advance sample state, stats, and store to Assembly
            sample.state = 2
            sample.stats_s2 = Stats2()

            # -- if pair: updated in remote_pairs_merge_overlap_with_vsearch
            sample.stats_s2.merged_pairs = self._merged_pairs[sname]
            sample.stats_s2.merged_pairs_prop = float(self._merged_pairs_prop[sname])
            # -- not relevant here
            # reads_mapped_to_ref: int = None
            # reads_mapped_to_ref_prop: float = None
            # -- if ref_filter: updated in remote_map_to_reference_as_filter
            # reads_mapped_to_ref_filter: int = None
            # reads_mapped_to_ref_filter_prop: float = None
            self.data.samples[sname] = sample
        self.data.save_json()

    def write_stats_file(self):
        """Write a easily readable tabular stats output file."""
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=[
                "merged_pairs",
                "merged_pairs_prop",
                "reads_mapped_to_ref",
                "reads_mapped_to_ref_prop",
                "reads_mapped_to_ref_filter",
                "reads_mapped_to_ref_filter_prop",
            ],
        )
        for sname, sample in self.samples.items():
            statsdict = sample.stats_s2.model_dump()
            for i in statsdf.columns:
                if statsdict[i]:
                    statsdf.loc[sname, i] = statsdict[i]
        handle = self.data.stepdir / 's2_cluster_prep_stats.txt'
        # self.data.stats_files.s2 = ...
        with open(handle, 'w', encoding="utf-8") as outfile:
            statsdf.fillna(value=0).to_string(outfile)
        logger.info(f"step 2 results written to {handle}")


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    # data = load_json("/tmp/pairgbs_merge.json")
    # data = load_json("/tmp/pedtest/half-demuxed.json")
    # data.ipcluster['threads'] = 4
    # with ip.Cluster(cores=8) as ipyclient:
    #     tool = ClustMapDenovo(data, True, ipyclient)
    #     tool.run()
    data = ip.load_json("/tmp/ipyrad-tests/assembly/TEST-denovo-se.json")
    data.run("2", force=True, cores=6, threads=2)
    print(data.stats)
    print(data.samples["40578_rex_SRR1754724"].files)
    print(data.samples["40578_rex_SRR1754724"].stats_s2)

    # # branch to subsample
    # subs = [
    #     "kansuensis-DE284",
    #     "kansuensis-DE704",
    #     "kansuensis-DE739",
    #     "kansuensis-DE366",
    #     "kansuensis-DE662",
    #     "lineata-DE757",
    # ]

    # # for clust in [0.95, 0.90, 0.85]:
    # for clust in [0.90]:
    #     ndata = data.branch(f"test-c{int(clust * 100)}", subsample=subs)
    #     ndata.params.clust_threshold = clust

    #     with ip.Cluster(cores=8) as ipyclient:
    #         tool = ClustMapDenovo(ndata, True, ipyclient)
    #         tool.run()

        # tool.remote_index_reference_as_filter()
    #     tool.remote_concat_multiple_fastqs_from_merged_sample()
    #     tool.remote_map_to_reference_as_filter()
    #     tool.remote_pairs_merge_overlaps_with_vsearch()
    #     tool.remote_pairs_join_unmerged_end_to_end()
    #     tool.remote_dereplicate()
    #     tool.remote_cluster()
    #     tool.write_json_file()
