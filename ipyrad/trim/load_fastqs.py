#!/usr/bin/env python

"""Load fastqs from one or more folders and join pairs together into
samples if PE data.

"""

from typing import Dict
import gzip
from pathlib import Path
from loguru import logger
import numpy as np
import pandas as pd
from ipyrad.demux.infer_overhang import infer_overhang
from ipyrad.core import BaseStep, track_remote_jobs
from ipyrad.schema import Sample
from ipyrad.schema.sample_schema import Stats1
from ipyrad.demux.pair_fastqs import (
    get_paths_list_from_fastq_str,
    get_fastq_tuples_dict_from_paths_list,
)

logger = logger.bind(name="ipyrad")


class Step1(BaseStep):
    """Load and filter/trim fastq input files to create Samples."""
    def __init__(self, data, force, ipyclient):
        # creates output directory
        super().__init__(data, step=1, force=force)
        self._filenames_to_fastq_tuples = {}
        self._samples: Dict[str, Sample] = {}
        self._results: Dict[str, Dict[str, float]] = {}
        self.ipyclient = ipyclient

    def run(self) -> None:
        self._load_fastqs()
        self._check_re_overhangs()
        self._run_fastp()
        self._create_samples()
        self._write_json_file()
        self._write_stats_file()

    def _load_fastqs(self) -> None:
        """Expand fastq path (e.g., ./data/*.gz) to {sample: (R1, R2)}."""
        fqs = get_paths_list_from_fastq_str(self.data.params.fastq_paths)
        snames_to_fq_tuples = get_fastq_tuples_dict_from_paths_list(fqs)
        self._filenames_to_fastq_tuples = snames_to_fq_tuples

    def _check_re_overhangs(self) -> None:
        """Infer re overhang from kmer analysis and compare w/ params setting."""
        read_r1s = [i[0] for i in self._filenames_to_fastq_tuples.values()]
        re1 = infer_overhang(read_r1s, max_reads=10_000)
        logger.debug(f"inferred R1 restriction overhang as: {re1}")

        read_r2s = [i[1] for i in self._filenames_to_fastq_tuples.values()]
        re2 = infer_overhang(read_r2s, max_reads=10_000)
        logger.debug(f"inferred R2 restriction overhang as: {re2}")

        # set values if no values were previously set
        if not self.data.params.restriction_overhang[0]:
            self.data.params.restriction_overhang = (re1, re2)
            logger.info(f"saving 'restriction overhang' as ({re1}, {re2}) to params.")
        else:
            ure1 = self.data.params.restriction_overhang[0]
            if ure1 != re1:
                logger.warning(
                    "kmer analysis identified the read1 restriction overhang "
                    f"as '{re1}', however, you entered {ure1}. Your entered "
                    "value is being used for now, but we recommend comparing "
                    f"your results with a run using {re1}."
                )
            ure2 = self.data.params.restriction_overhang[1]
            if ure2 != re2:
                logger.warning(
                    "kmer analysis identified the read2 restriction overhang "
                    f"as '{re2}', however, you entered {ure2}. Your entered "
                    "value is being used for now, but we recommend comparing "
                    f"your results with a run using {re2}."
                )

    def _run_fastp(self) -> None:
        """The main function to call fastp in parallel."""

        # function to call ReadTrimming on a remote engine.
        def trim_reads(**kwargs):
            from ipyrad.trim.trim_fastqs import TrimFastqs
            tool = TrimFastqs(**kwargs)
            tool.run()
            return tool.parse_stats_from_json()

        # arguments to trim reads, some updated below.
        kwargs = dict(
            # sname="TODO",
            # reads="TODO",
            # is_pair="TODO",
            outdir=self.data.stepdir,
            restriction_overhang=self.data.params.restriction_overhang,
            trim_reads=self.data.params.trim_reads,
            filter_min_trim_len=self.data.params.filter_min_trim_len,
            phred_qscore_offset=self.data.params.phred_qscore_offset,
            max_low_qual_bases=self.data.params.max_low_qual_bases,
            filter_adapters=self.data.params.filter_adapters,
        )

        # serial execution option
        if self.ipyclient is None:
            logger.warning("No ipyclient, executing serially.")
            for sname, fastq_tuple in self._filenames_to_fastq_tuples.items():
                kwargs['sname'] = sname
                kwargs['reads'] = fastq_tuple
                kwargs['is_pair'] = bool(fastq_tuple[1])
                self._results[sname] = trim_reads(**kwargs)
            return

        # Parallel execution
        lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::4])

        # submit jobs to parallel client
        rasyncs = {}
        for sname, fastq_tuple in self._filenames_to_fastq_tuples.items():
            kwargs['sname'] = sname
            kwargs['reads'] = fastq_tuple
            kwargs['is_pair'] = bool(fastq_tuple[1])
            rasyncs[sname] = lbview.apply(trim_reads, **kwargs)
        self._results = track_remote_jobs(rasyncs, self.ipyclient)

    def _create_samples(self) -> None:
        for sname, result in self._results.items():
            j, filepaths = result
            sample = Sample(name=sname)
            sample.files.fastqs = self._filenames_to_fastq_tuples[sname]
            sample.files.trimmed = filepaths
            sample.stats_s1 = Stats1(
                reads_raw=j['summary']['before_filtering']['total_reads'],
                reads_filtered_by_Ns=j['filtering_result']['too_many_N_reads'],
                reads_filtered_by_low_quality=j['filtering_result']['low_quality_reads'],
                reads_filtered_by_low_complexity=j['filtering_result']['low_complexity_reads'],
                reads_filtered_by_minlen=j['filtering_result']['too_short_reads'],
                adapter_trimmed_reads=j['adapter_cutting']['adapter_trimmed_reads'],
                adapter_trimmed_bases=j['adapter_cutting']['adapter_trimmed_bases'],
                mean_len_R1_before_trimming=j['summary']['before_filtering']['read1_mean_length'],
                mean_len_R1_after_trimming=j['summary']['after_filtering']['read1_mean_length'],
                mean_len_R2_before_trimming=j['summary']['before_filtering'].get("read2_mean_length", None),
                mean_len_R2_after_trimming=j['summary']['after_filtering'].get("read2_mean_length", None),
                reads_passed_filter=j['summary']['after_filtering']['total_reads'],
            )

            # for paired data these values make more sense when divided
            # in half to represent "read pairs" instead of reads.
            pair_keys = [
                "reads_raw",
                "reads_filtered_by_Ns",
                "reads_filtered_by_low_quality",
                "reads_filtered_by_low_complexity",
                "reads_filtered_by_minlen",
                "reads_passed_filter"
            ]
            # if paired
            if sample.files.fastqs[1]:
                for key in pair_keys:
                    setattr(sample.stats_s1, key, int(getattr(sample.stats_s1, key) / 2))
            # save sample to sample dict
            self._samples[sname] = sample

    def _write_json_file(self):
        """Writes samples to the JSON file."""
        # only advance the STATE for samples that were successful
        for sname, sample in self._samples.items():
            if not sample.stats_s1.reads_passed_filter:
                logger.warning(
                    f"sample {sname} has 0 reads after filtering "
                    "and will be excluded.")
            else:
                sample.state = 1
                self.data.samples[sname] = sample

        # put samples into Assembly and save updated JSON
        self.data.save_json()

    def _write_stats_file(self):
        """Write a easily readable tabular stats output file."""
        statsdf = pd.DataFrame(
            index=sorted(self._samples),
            columns=[
                "reads_raw",
                "reads_passed_filter",
                "reads_filtered_by_Ns",
                "reads_filtered_by_low_quality",
                "reads_filtered_by_low_complexity",
                "reads_filtered_by_minlen",
                "adapter_trimmed_reads",
                "adapter_trimmed_bases",
                "mean_len_R1_before_trimming",
                "mean_len_R2_before_trimming",
                "mean_len_R1_after_trimming",
                "mean_len_R2_after_trimming",
            ],
        )
        for sname in self._samples:
            sample = self.data.samples[sname]
            statsdict = sample.stats_s1.dict()
            for i in statsdf.columns:
                if statsdict[i]:
                    statsdf.loc[sname, i] = statsdict[i]
        handle = self.data.stepdir / 's1_read_stats.txt'
        # self.data.stats_files.s1
        with open(handle, 'w', encoding="utf-8") as outfile:
            statsdf.fillna(value=0).to_string(outfile)
        logger.info(f"step 1 results written to {handle}")


def estimate_trim_position(
    fastq: Path,
    restriction_overhang: str,
    max_len: int = 25,
    max_reads: int = 200_000,
) -> int:
    """Return estimation position of end of cut1 on R1 files.

    At this point read1 files should have one of two formats:
    >>> ACACACTGCAGXXXX....
    >>> bbbbbb^^^^^dddd
    or
    >>> TGCAGXXXX....
    >>> ^^^^^dddd
    where
    b = barcode
    ^ = cut overhang
    d = data

    This function will read the first max_reads reads to find the avg
    position of the last cutter overhang position (i.e., end of barcoe).
    If the cutter overhang is not in >20% of reads then 0 is returned.

    Note
    ----
    We use rfind from 25 bases in to accommodate cases in which the re
    sequence occurs within the barcode.
    """
    if str(fastq).endswith(".gz"):
        fp = gzip.open(fastq, 'rt', encoding="utf-8")
    else:
        fp = open(fastq, 'rt', encoding="utf-8")
    quart = zip(fp, fp, fp, fp)
    count = range(max_reads)

    observed = []
    for idx, line in zip(count, quart):
        try:
            match = line[1][:max_len].rindex(restriction_overhang)
            observed.append(match)
        except ValueError:
            pass

    # must occur in >20%
    if len(observed) < max_reads * 0.20:
        return 0
    # trim up to Nth position, fastp takes trim N bases...
    return int(np.round(np.mean(observed)))


if __name__ == "__main__":

    import ipyrad as ip
    # from ipyrad.core2.assembly import Assembly
    ip.set_log_level("DEBUG")

    # data.params.filter_adapters = 2

    # with ip.Cluster(cores=8) as ipyclient:
    #     step = Step1(data, True, ipyclient)
    #     step._load_fastqs()
    #     step._check_re_overhangs()
    #     step._run_fastp()
    #     step._create_samples()
    #     step._write_json_file()
    #     step._write_stats_file()

    # data = Assembly("pairgbs_merge")
    # data.params.fastq_path = "/tmp/demux_pairgbs/*.gz"

    data = ip.Assembly(name="half-demuxed")
    data.params.fastq_paths = [
        "../../pedtest/demux_2023-3-28/linea*.gz",
        "../../pedtest/demux_2023-3-28/kansu*.gz",
        # "../../pedtest/demux_2023-3-28/lachno*.gz",
    ]
    data.params.project_dir = "../../pedtest/"

    with ip.Cluster(cores=5) as ipyclient:
        step = Step1(data, True, ipyclient)
        step.run()

    # FASTQS = list(Path("../../pedtest/NEW_fastqs").glob("*fastq.gz"))
    # pairs = get_filenames_to_paired_fastqs(FASTQS)

    # R1s = [i[0] for i in pairs.values()]
    # print("inferred R1 re overhang is:", infer_overhang(R1s, max_reads=10_000))

    # R2s = [i[1] for i in pairs.values()]
    # print("inferred R2 re overhang is:", infer_overhang(R2s, max_reads=10_000))