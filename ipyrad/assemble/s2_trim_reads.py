#!/usr/bin/env python

"""
Trim read (pairs) for adapters and low quality bases using 'fastp'.
This is processes the demultiplexed data files, and is intended that
the files have not been preprocessed by already. 

Further read processing to merge pairs or derep takes place in Step3.
"""

import os
import sys
import json
import subprocess as sps
import pandas as pd
from loguru import logger
from ipyrad.assemble.utils import IPyradError
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.base_step import BaseStep
from ipyrad.core.schema import Stats2


class Step2(BaseStep):
    """
    Run Step2 read trimming and filtering using fastp.
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, step=2, quiet=quiet, force=force)
        self.ctuple_dict = {}

        # threading: default 4-threads per job.
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view(
            self.ipyclient.ids[::4]
        )

    def run(self):
        """
        Run the component subfuctions to concatenate input files, 
        run the fastp binary, and cleanup tmp files.
        """
        self.concat_multiple_raws()
        self.distribute_fastp()

        # json must run first
        self.write_json_file()
        self.write_stats_file()


    def concat_multiple_raws(self):
        """
        Don't bother to parallelize this since its rarely run and when
        it is it's pretty I/O limited.

        Concatenate multiple raw files into a single 'concat' file. 
        This is needed when data from multiple step1 (e.g., demux) 
        runs are 'merged' into a single assembly object, or when a 
        demux run included merged technical replicates.
        """
        for sname in self.samples:
            sample = self.samples[sname]

            # check if user has moved any files since step 1
            for ftuple in sample.files.fastqs:
                for fastq in ftuple:
                    if fastq:
                        if not os.path.exists(fastq):
                            raise IOError(
                                f"fastq files from step 1 are missing: {fastq}")

            # if only one set of files then no concatenation is necessary
            if len(sample.files.fastqs) == 1:
                self.ctuple_dict[sname] = self.samples[sname].files.fastqs[0]

            # otherwise, we need to concatenate multiple fastq files.
            else:
                logger.debug(
                    "concatenating inputs: {}".format(sample.files.fastqs))

                # cat command works for gzip or not. Index 0 of tuples is R1s.
                cmd1 = ["cat"] + [i[0] for i in sample.files.fastqs]

                isgzip = ".gz"
                if not sample.files.fastqs[0][0].endswith(".gz"):
                    isgzip = ""

                # write to new concat handle
                conc1 = os.path.join(self.tmpdir, 
                    sample.name + "_R1_concat.fastq"
                )
                if sample.files.fastqs[0][0].endswith(".gz"):
                    conc1 += ".gz"

                # call and write to open outfile
                with open(conc1, 'w') as cout1:
                    proc1 = sps.Popen(
                        cmd1, stderr=sps.STDOUT, stdout=cout1, close_fds=True)
                    res1 = proc1.communicate()[0]

                # check for errors
                if proc1.returncode:
                    raise IPyradError("error in: {}, {}".format(cmd1, res1))

                # Only set conc2 if R2 actually exists
                conc2 = 0
                if "pair" in self.data.params.datatype:
            
                    # out _R2 filehandle
                    conc2 = os.path.join(
                        self.tmpdir,
                        sample.name + "_R2_concat.fq{}".format(isgzip)
                    )

                    # write concat results directly to _concat outfile.
                    cmd2 = ["cat"] + [i[1] for i in sample.files.fastqs]
                    with open(conc2, 'w') as cout2:
                        proc2 = sps.Popen(
                            cmd2, stderr=sps.STDOUT, stdout=cout2, close_fds=True)
                        res2 = proc2.communicate()[0]
            
                    # check for errors
                    if proc2.returncode:
                        raise IPyradError(
                            "Error concatenating fastq files. Make sure all "
                            "these files exist: {}\nError message: {}"
                            .format(cmd2, res2))

                # store new file handles
                self.ctuple_dict[sname] = (conc1, conc2)
                logger.debug('concat: {}'.format(self.ctuple_dict[sname]))


    def distribute_fastp(self):
        """
        Distribute N / 2 fastp jobs on N processors.
        """
        # send samples to cutadapt filtering
        logger.debug("submitting remote fastp jobs")

        def trim_reads(data, sname, read1, read2, workdir):
            """
            Remote function for applying fastp read trimming.
            """
            tool = ReadTrimming(
                data=data,
                sname=sname,
                read1=read1,
                read2=(read2 if read2 else None),
                workdir=workdir,
            )
            tool.run()
            return tool.parse_stats_from_json()

        # submit jobs to run
        jobs = {}
        for sname in self.samples:
            args = (
                self.data,
                sname, 
                self.ctuple_dict[sname][0], 
                self.ctuple_dict[sname][1],
                self.stepdir,
            )
            # logger.debug(" ".join(ReadTrimming(*args).command))
            rasync = self.lbview.apply(trim_reads, *args)
            jobs[sname] = rasync

        # wait for all to finish
        message = "processing reads"
        prog = AssemblyProgressBar(jobs, message, step=2, quiet=self.quiet)
        prog.block()
        prog.check()
        
        # check for errors
        for sname in prog.results:

            # enter stats on the sample
            jdata, filepaths = prog.results[sname]

            # store data
            self.samples[sname].files.edits = filepaths

            # store results
            stats = Stats2(
                reads_raw=jdata['summary']['before_filtering']['total_reads'],
                reads_filtered_by_Ns=jdata['filtering_result']['too_many_N_reads'],
                reads_filtered_by_low_quality=jdata['filtering_result']['low_quality_reads'],
                reads_filtered_by_low_complexity=jdata['filtering_result']['low_complexity_reads'],
                reads_filtered_by_minlen=jdata['filtering_result']['too_short_reads'],
                mean_len_R1_before_trimming=jdata['summary']['before_filtering']['read1_mean_length'],
                mean_len_R1_after_trimming=jdata['summary']['after_filtering']['read1_mean_length'],
                mean_len_R2_before_trimming=(
                    jdata['summary']['before_filtering'].get("read2_mean_length", None)),
                mean_len_R2_after_trimming=(
                    jdata['summary']['after_filtering'].get("read2_mean_length", None)),
                reads_passed_filter=jdata['summary']['after_filtering']['total_reads'],
            )
            self.samples[sname].stats_s2 = stats
            # for key in [
            #     "reads_raw", 
            #     "reads_filtered_by_Ns", 
            #     "reads_filtered_by_minlen", 
            #     "reads_filtered_by_low_complexity", 
            #     "reads_filtered_by_low_quality", 
            #     "reads_passed_filter"
            #     ]:
            #     sample.stats_dfs.s2[key] = sample.stats_dfs.s2[key] / 2


    def write_json_file(self):
        """
        Writes samples to the JSON file.        
        """
        # only advance the STATE for samples that were successful
        for sname in self.samples:
            if not self.samples[sname].stats_s2.reads_passed_filter:
                logger.warning(
                    f"sample {sname} has 0 reads after filtering "
                    "and will be excluded.")
            else:
                self.samples[sname].state = 2

        # put samples into Assembly and save updated JSON
        self.data.samples = self.samples
        self.data.save_json()


    def write_stats_file(self):
        """
        Write a easily readable tabular stats output file.
        """
        statsdf = pd.DataFrame(
            index=sorted(self.data.samples),
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
        for sname in self.data.samples:
            statsdict = self.data.samples[sname].stats_s2.dict()
            for i in statsdf.columns:
                if statsdict[i]:
                    statsdf.loc[sname, i] = statsdict[i]

        handle = os.path.join(self.stepdir, 's2_read_trim_stats.txt')
        with open(handle, 'w') as outfile:
            statsdf.fillna(value=0).to_string(outfile)
        logger.info("\n" + 
            statsdf.loc[:, ["reads_raw", "reads_passed_filter"]].to_string()
        )



class ReadTrimming:
    """
    Simple read trimming with fastp 
    https://github.com/OpenGene/fastp
    """
    def __init__(self, data, sname, read1, read2, workdir):

        # input args
        self.data = data
        self.read1 = os.path.realpath(os.path.expanduser(read1))
        self.read2 = (
            os.path.realpath(os.path.expanduser(read2)) if read2
            else None
        )
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.is_paired = self.read2 is not None
        self.fastp_binary = os.path.join(sys.prefix, "bin", "fastp")

        # output file paths (do not bother gzipping since these are tmp files)
        self.out1 = os.path.join(self.workdir, sname + ".trimmed_R1.fastq.gz")
        self.out2 = (
            os.path.join(self.workdir, sname + ".trimmed_R2.fastq.gz")
            if self.is_paired else ""
        )

        # paths to stats files
        self.json = os.path.join(self.workdir, f'{sname}.fastp.json')
        self.html = os.path.join(self.workdir, f'{sname}.fastp.html')

        # get the command
        self.command = []
        self.build_command()


    def build_command(self):
        """
        Calls fastp in subprocess and writes tmpfiles to workdir.
        0 = no trim or quality check.
        1 = only trimming
        2 = trimming and quality checks.

        (-A turns off adapter trimming)
        (-Q turns off quality filtering)
        """
        if self.is_paired:
            cmd = [
                self.fastp_binary, 
                "-i", self.read1,
                "-I", self.read2,
                "-o", self.out1,
                "-O", self.out2,
                "--detect_adapter_for_pe"
            ]
            # if we merged in step2 this would prevent branching at 
            # step 3 to try both denovo and reference, so not doing it.
            # if self.data.params.assembly_method:
                # cmd += ['-m', '--merged_out', self.out1 + ".merged"]

        else:
            cmd = [
                self.fastp_binary,
                "-i", self.read1,
                "-o", self.out1,
            ]

        cmd.extend([
            "-5",  # sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-3",  # sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
            # "-a", str(self.data.hackersonly.p5_adapter), # allow auto-detection
            "-q", str(20 + self.data.params.phred_qscore_offset - 33),  # minqual
            "-l", str(self.data.params.filter_min_trim_len),  # minlen
            "-y", "-Y", "50",  # complexity filter
            "-c",              # paired-end base correction
            "-w", "3",         # 3 proc. threads, 1 i/o thread.
            "--n_base_limit", str(self.data.params.max_low_qual_bases),
            "-j", self.json,
            "-h", self.html,
            "--trim_front1", str(self.data.params.trim_reads[0]),
            "--trim_tail1", str(self.data.params.trim_reads[1]),
            "--trim_front2", str(self.data.params.trim_reads[2]),
            "--trim_tail2", str(self.data.params.trim_reads[3]),
        ])

        # turn off filters if settings are lower than 2
        if self.data.params.filter_adapters == 1:
            cmd.extend("-A")
        if self.data.params.filter_adapters == 0:
            cmd.extend("-Q")
        self.command = cmd


    def run(self):
        """
        Run the command on a subprocess
        """
        print(" ".join(self.command))
        proc = sps.Popen(self.command, stderr=sps.STDOUT, stdout=sps.PIPE)
        out = proc.communicate()
        if proc.returncode:
            logger.error("FASTP ERROR: {}".format(out))
            raise IPyradError(out[0].decode())


    def parse_stats_from_json(self):
        """
        Get stats from JSON file.
        """
        with open(self.json, 'r') as indata:
            jdata = json.loads(indata.read())
        return (jdata, [(self.out1, self.out2)])



if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG", stderr=False, logfile="/tmp/test.log")

    # SE DATA
    TEST = ip.load_json("/tmp/TESTB.json")
    TEST.run("2", force=False, quiet=True)

    # PE DATA
    # TEST = ip.load_json("/tmp/TEST5.json")
    # TEST.run("2", force=True, quiet=True)

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

