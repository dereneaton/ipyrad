#!/usr/bin/env python

"""Trim read (pairs) for adapters and low quality bases using 'fastp'.

This processes the demultiplexed data files, and expects that
the files have not been preprocessed by other tools already.

Further read processing to merge pairs or derep takes place in Step3.

Notes
------
The behavior in v.1.0 changed to automatically trim the restriction
overhangs from sequences during step2. They had previously been
retained because it improved subsqeuent alignments, and made is easier
to tell that reads were oriented properly. However, we now include
NNNNN padding which serves fine for improving edge alignements, and
so removing the overhangs altogether makes for cleaner results. To
suppress this automatic trimming behavior you can set
Assembly.params.trim_reads = ()
"""

from typing import TypeVar
import sys
import json
from pathlib import Path
import subprocess as sps

import pandas as pd
from loguru import logger
from ipyrad.assemble.utils import IPyradError
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.base_step import BaseStep
from ipyrad.core.schema import Stats2

Assembly = TypeVar("Assembly")
logger = logger.bind(name="ipyrad")


class Step2(BaseStep):
    """Run Step2 read trimming and filtering using fastp.

    The program fastp auto-detects the adapter sequences by
    analyzing the first 1M reads. In addition to those, we tell it
    to look for the adapters in the `./adapter.fa` file.
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, step=2, quiet=quiet, force=force)
        self.ctuple_dict = {}

        # threading: default 4-threads per job.
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view(self.ipyclient.ids[::4])

    def run(self):
        """Run the component subfuctions to concatenate input files,
        run the fastp binary, and cleanup tmp files."""
        self._concat_multiple_raws()
        self._distribute_fastp()
        self._write_json_file()
        self._write_stats_file()

    def _concat_multiple_raws(self):
        """Fills the `ctuples` dict with fastq paths {sname: (R1,R2)}.

        Don't bother to parallelize this since its rarely run and
        when it is it's pretty I/O limited.

        Concatenate multiple raw files into a single 'concat' file.
        This is needed when data from multiple step1 (e.g., demux)
        runs are 'merged' into a single assembly object, or when a
        demux run included merged technical replicates.
        """
        for sname, sample in self.data.samples.items():

            # check if user has moved any files since step 1
            for ftuple in sample.files.fastqs:
                for fastq in ftuple:
                    if fastq:
                        if not Path(fastq).exists():
                            raise IOError(
                                f"fastq files from step 1 are missing: {fastq}")

            # if only one set of files then no concatenation is necessary
            if len(sample.files.fastqs) == 1:
                self.ctuple_dict[sname] = sample.files.fastqs[0]

            # otherwise, we need to concatenate multiple fastq files.
            else:
                logger.info(f"concatenating inputs: {sample.files.fastqs}")

                # cat command works for gzip or not. Index 0 of tuples is R1s.
                cmd1 = ["cat"] + [i[0] for i in sample.files.fastqs]

                # write to new concat handle
                conc1 = self.data.tmpdir / f"{sample.name}_R1_concat.fastq"
                if sample.files.fastqs[0][0].endswith(".gz"):
                    conc1 = conc1.with_suffix(".gz")

                # call and write to open outfile
                with open(conc1, 'w', encoding="utf-8") as cout1:
                    with sps.Popen(
                        cmd1, stderr=sps.STDOUT, stdout=cout1, close_fds=True
                        ) as proc1:
                        res1 = proc1.communicate()[0]
                if proc1.returncode:
                    raise IPyradError(f"error in: {cmd1}, {res1}")

                # Only set conc2 if R2 actually exists
                conc2 = 0
                if self.data.is_pair:

                    # out _R2 filehandle
                    conc2 = self.data.tmpdir / f"{sample.name}_R2_concat.fastq"
                    if sample.files.fastqs[0][0].endswith(".gz"):
                        conc2 = conc2.with_suffix(".gz")

                    # write concat results directly to _concat outfile.
                    cmd2 = ["cat"] + [i[1] for i in sample.files.fastqs]
                    with open(conc2, 'w', encoding="utf-8") as cout2:
                        with sps.Popen(
                            cmd2, stderr=sps.STDOUT, stdout=cout2, close_fds=True,
                            ) as proc2:
                            res2 = proc2.communicate()[0]
                    if proc2.returncode:
                        raise IPyradError(
                            "Error concatenating fastq files. Make sure all "
                            f"these files exist: {cmd2}\nError message: {res2}")

                # store new file handles
                self.ctuple_dict[sname] = (conc1, conc2)
                logger.info(f'concat: {self.ctuple_dict[sname]}')

    def _distribute_fastp(self):
        """Distribute N / 2 fastp jobs on N processors.

        """
        logger.info("submitting remote fastp jobs")

        def trim_reads(**kwargs):
            """Remote function for applying fastp read trimming."""
            tool = ReadTrimming(**kwargs)
            tool.run()
            return tool.parse_stats_from_json()

        # submit jobs to run
        jobs = {}
        for sname in self.data.samples:
            kwargs = dict(
                data=self.data,
                sname=sname,
                read1=self.ctuple_dict[sname][0],
                read2=self.ctuple_dict[sname][1],
            )
            # logger.debug(" ".join(ReadTrimming(*args).command))
            jobs[sname] = self.lbview.apply(trim_reads, **kwargs)

        # wait for all to finish
        message = "processing reads"
        prog = AssemblyProgressBar(jobs, message, step=2, quiet=self.quiet)
        prog.block()
        prog.check()

        # store file paths and stats
        for sname, result in prog.results.items():
            jdata, filepaths = result
            # store a list of tuples of strings
            self.data.samples[sname].files.edits = filepaths
            stats = Stats2()
            stats.reads_raw = jdata['summary']['before_filtering']['total_reads']
            stats.reads_filtered_by_Ns = jdata['filtering_result']['too_many_N_reads']
            stats.reads_filtered_by_low_quality = jdata['filtering_result']['low_quality_reads']
            stats.reads_filtered_by_low_complexity = jdata['filtering_result']['low_complexity_reads']
            stats.reads_filtered_by_minlen = jdata['filtering_result']['too_short_reads']
            stats.mean_len_R1_before_trimming = jdata['summary']['before_filtering']['read1_mean_length']
            stats.mean_len_R1_after_trimming = jdata['summary']['after_filtering']['read1_mean_length']
            stats.mean_len_R2_before_trimming = (
                jdata['summary']['before_filtering'].get("read2_mean_length", None))
            stats.mean_len_R2_after_trimming = (
                jdata['summary']['after_filtering'].get("read2_mean_length", None))
            stats.reads_passed_filter = jdata['summary']['after_filtering']['total_reads']

            # store reads numbers as unpaired count unlike in fastp
            pair_keys = [
                "reads_raw",
                "reads_filtered_by_Ns",
                "reads_filtered_by_low_quality",
                "reads_filtered_by_low_complexity",
                "reads_filtered_by_minlen",
                "reads_passed_filter"
            ]
            if self.data.is_pair:
                for key in pair_keys:
                    setattr(stats, key, int(getattr(stats, key) / 2))

            # store to Sample object
            self.data.samples[sname].stats_s2 = stats

    def _write_json_file(self):
        """Writes samples to the JSON file."""
        # only advance the STATE for samples that were successful
        for sname, sample in self.samples.items():
            if not sample.stats_s2.reads_passed_filter:
                logger.warning(
                    f"sample {sname} has 0 reads after filtering "
                    "and will be excluded.")
            else:
                sample.state = 2
                sample._clear_old_results()

        # put samples into Assembly and save updated JSON
        self.data.save_json()

    def _write_stats_file(self):
        """Write a easily readable tabular stats output file."""
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=[
                "reads_raw",
                "reads_passed_filter",
                "reads_filtered_by_Ns",
                "reads_filtered_by_low_quality",
                "reads_filtered_by_low_complexity",
                "reads_filtered_by_minlen",
                "mean_len_R1_before_trimming",
                "mean_len_R2_before_trimming",
                "mean_len_R1_after_trimming",
                "mean_len_R2_after_trimming",
            ],
        )
        for sname in self.samples:
            sample = self.data.samples[sname]
            statsdict = sample.stats_s2.dict()
            for i in statsdf.columns:
                if statsdict[i]:
                    statsdf.loc[sname, i] = statsdict[i]

        handle = self.data.stepdir / 's2_read_trim_stats.txt'
        with open(handle, 'w', encoding="utf-8") as outfile:
            statsdf.fillna(value=0).to_string(outfile)
        logger.info(
            "\n" +
            statsdf.loc[:, ["reads_raw", "reads_passed_filter"]].to_string()
        )


class ReadTrimming:
    """Simple read trimming with fastp.

    https://github.com/OpenGene/fastp
    """
    def __init__(self, data: Assembly, sname: str, read1: Path, read2: Path):

        # input args are passed as full paths
        self.data = data
        self.read1 = read1
        self.read2 = read2 if read2 else None
        self.fastp_binary = Path(sys.prefix) / "bin" / "fastp"

        # output file paths (do not bother gzipping since these are tmp files)
        self.out1 = self.data.stepdir / f"{sname}.trimmed_R1.fastq.gz"
        self.out2 = (
            self.data.stepdir / f"{sname}.trimmed_R2.fastq.gz"
            if self.data.is_pair else "")

        # paths to stats files
        self.json = self.data.stepdir / f'{sname}.fastp.json'
        self.html = self.data.stepdir / f'{sname}.fastp.html'

        # get the command
        self.command = []
        self.build_command()

    def build_command(self):
        """Calls fastp in subprocess and writes tmpfiles to workdir.

        0 = no trim or quality check.
        1 = only trimming
        2 = trimming and quality checks.

        (-A turns off adapter trimming)
        (-Q turns off quality filtering)
        """
        if self.data.is_pair:
            cmd = [
                str(self.fastp_binary),
                "-i", str(self.read1),
                "-I", str(self.read2),
                "-o", str(self.out1),
                "-O", str(self.out2),
                "--detect_adapter_for_pe"
            ]
            # if we merged in step2 this would prevent branching at
            # step 3 to try both denovo and reference, so not doing it.
            # if self.data.params.assembly_method:
                # cmd += ['-m', '--merged_out', self.out1 + ".merged"]
        else:
            cmd = [
                str(self.fastp_binary),
                "-i", str(self.read1),
                "-o", str(self.out1),
            ]

        # by default we trim the adapter lengths from the start of R1
        # and R2 reads based on length of `restriction_overhang` params.
        trim_front1 = len(self.data.params.restriction_overhang[0])
        trim_front2 = len(self.data.params.restriction_overhang[1])

        # HOwever, this can be overriden by setting `trim_reads` param.
        # to either a set of fixed positive integer lengths, or, if set
        # to a negative value, then NO TRIMMING will be done.
        if self.data.params.trim_reads[0]:
            if self.data.params.trim_reads[0] < 0:
                trim_front1 = 0
            else:
                trim_front1 = abs(self.data.params.trim_reads[0])
        if self.data.params.trim_reads[1]:
            if self.data.params.trim_reads[1] < 0:
                trim_front2 = 0
            else:
                trim_front2 = abs(self.data.params.trim_reads[1])

        cmd.extend([
            "-5",  # sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-3",  # sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
            # "-a", str(self.data.hackersonly.p5_adapter), # disable to allow auto-detection
            "-q", str(20 + self.data.hackers.phred_qscore_offset - 33),  # minqual
            "-l", str(self.data.params.filter_min_trim_len),  # minlen
            "-y", "-Y", "50",  # turns on and sets complexity filter to 50
            "-c",              # paired-end base correction
            "-w", "3",         # 3 proc. threads, 1 i/o thread.
            "--n_base_limit", str(self.data.params.max_low_qual_bases),
            "-j", str(self.json),
            "-h", str(self.html),
            "--trim_front1", str(trim_front1),
            "--trim_front2", str(trim_front2),
            # "--trim_tail1", str(abs(self.data.params.trim_reads[1])),
            # "--trim_tail2", str(abs(self.data.params.trim_reads[3])),
        ])

        # turn off filters if settings are lower than 2
        # hard coded fasta adapters file with Truseq and AAAAAAA
        extra_adapters = Path(__file__).parent / "adapters.fa"
        if self.data.params.filter_adapters == 2:
            cmd.extend(["--adapter_fasta", str(extra_adapters)])
        if self.data.params.filter_adapters == 1:
            cmd.extend("-A")
        if self.data.params.filter_adapters == 0:
            cmd.extend("-Q")
        self.command = cmd

    def run(self):
        """Run the fastp command on a subprocess"""
        print(" ".join(self.command))  # sends to logger.INFO on engines.
        with sps.Popen(self.command, stderr=sps.STDOUT, stdout=sps.PIPE) as proc:
            out = proc.communicate()
        if proc.returncode:
            logger.error(f"FASTP ERROR: {out}")
            raise IPyradError(out[0].decode())

    def parse_stats_from_json(self):
        """Get stats from the fastp JSON file."""
        with open(self.json, 'r', encoding="utf-8") as indata:
            jdata = json.loads(indata.read())
        return (jdata, [(str(self.out1), str(self.out2))])


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG", log_file="/tmp/test.log")

    # SE DATA
    # TEST = ip.load_json("/tmp/TEST.json")
    # TEST.run("2", force=False, quiet=True)

    # simulated PE DATA
    TEST = ip.load_json("/tmp/TEST5.json")
    TEST.run("2", force=True, quiet=False)
    print(TEST.stats)

    # # EMPIRICAL
    # TEST = ip.load_json("/tmp/PEDIC.json")
    # TEST.run("2", force=True, quiet=True)

    # TEST = ip.Assembly("PEDIC")
    # TEST.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    # TEST.params.project_dir = "/tmp"
    # TEST.run('12', force=True, quiet=True)
    # TEST = ip.load_json("/tmp/PEDIC.json")
    # TEST.run("2", force=True, quiet=True)

    # # Empirical pairddrad test
    # tdata = ip.load_json("/tmp/test-amaranth.json")
    # tdata.params.filter_adapters = 2
    # tdata.params.filter_min_trim_len = 50
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats)
    # print(tdata.stats_dfs.s2)

    # # Simulated SE RAD test
    # tdata = ip.load_json("/tmp/test-simrad.json")
    # tdata.params.filter_adapters = 2
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats.head())
    # print(tdata.stats_dfs.s2.head())

    # Simulated test that included merged assemblies
    # tdata = ip.load_json("/tmp/test-simrad.json")
    # tdata = ip.merge('merge-s2-test', [tdata, tdata])
    # tdata.params.filter_adapters = 2
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats.head())
    # print(tdata.stats_dfs.s2.head())


    # test trimming of PE-ddRAD fastq files
    # tdata = ip.load_json("/tmp/test-emppairgbs.json")
    # tdata.params.filter_adapters = 2
    # tdata.params.filter_min_trim_len = 50
    # tdata.params.trim_reads = (0, 0, 1, 0)
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats.T)
    # print(tdata.stats_dfs.s2.T)


    # # test trimming of PE-ddRAD fastq files
    # tdata = ip.load_json("/tmp/test-simpairddrad.json")
    # tdata.params.filter_adapters = 2
    # tdata.params.filter_min_trim_len = 50
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats.head())
    # print(tdata.stats_dfs.s2.head())


    # test by loading SE fastq files
    # tdata = ip.load_json("/tmp/test-simrad.json")
    # tdata.params.filter_adapters = 2
    # tdata.run("2", auto=True, force=True)
    # print(tdata.stats.head())


    # test on SE RAD data
    # tdata = ip.load_json("/home/deren/Documents/ipyrad/sandbox/oak-test.json")
    # tdata = idata.branch("test")
    # tdata.run("2", force=True, auto=True)

    # test on PE ddRAD data
    # idata = ip.load_json("/home/deren/Documents/ipyrad/sandbox/oak-test.json")
    # tdata = idata.branch("test")
    # tdata.run("2", force=True, auto=True)

