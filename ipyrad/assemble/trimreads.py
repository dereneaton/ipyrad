#!/usr/bin/env python

"""
Trim read (pairs) for adapters and low quality bases using 'fastp'.
This is processes the demultiplexed data files, and is intended that
the files have not been preprocessed by already. 

Further read processing to merge pairs or derep takes place in Step3.
"""


from __future__ import print_function

import os
import sys
import json
import subprocess as sps
from loguru import logger
from ipyrad.assemble.utils import IPyradError, AssemblyProgressBar




class Step2:
    """
    Run Step2 read trimming and filtering using fastp.
    """
    def __init__(self, data, force, ipyclient):

        # store input params
        self.data = data
        self.force = force
        self.ipyclient = ipyclient

        # threading: default 4-threads per job.
        self.lbview = self.ipyclient.load_balanced_view(
            self.ipyclient.ids[::4]
        )

        # to be filled
        self.samples = []

        # run setup funcs
        self.print_headers()
        self.get_subsamples()
        self.setup_dirs()
        self.check_binaries()



    def print_headers(self):
        """
        Print the CLI header 
        """
        if self.data._cli:
            self.data._print(
                "\n{}Step 2: Filtering and trimming reads"
                .format(self.data._spacer)
            )


    def get_subsamples(self):
        """
        Select samples ready for this step 
        """
        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 2")

        # filter samples by state
        state1 = self.data.stats.index[self.data.stats.state == 1]
        statex = self.data.stats.index[self.data.stats.state > 1]

        # build list to run for samples being forced
        if self.force:
            subsamples = list(self.data.samples.values())
        else:
            # tell user which samples have already completed step 2
            if statex.any():
                print("skipping samples already finished step 2:\n{}"
                      .format(statex.tolist()))
            # run all samples in state 1
            subsamples = [self.data.samples[i] for i in state1]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_raw:
                checked_samples.append(sample)
            else:
                print("skipping {}; no reads found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 3")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_raw,
            reverse=True,
        )
        self.samples = checked_samples


    def setup_dirs(self):
        """
        create the output dir for this step
        """
        self.data.dirs.edits = os.path.join(
            self.data.params.project_dir,
            "{}_edits".format(self.data.name))
        if not os.path.exists(self.data.dirs.edits):
            os.makedirs(self.data.dirs.edits)        


    def check_binaries(self):
        """
        Check that fastp binary can be found
        """
        binary = os.path.join(sys.prefix, "bin", "fastp")
        if not os.path.exists(binary):
            raise IPyradError("'fastp' not found. See conda installation")


    def run(self):
        """
        Run the component subfuctions to concatenate input files, 
        run the fastp binary, and cleanup tmp files.
        """
        self.concat_multiple_raws()
        self.distribute_fastp()

        # save the assembly        
        self.data.save()


    def concat_multiple_raws(self):
        """
        concatenate multiple raw files into a single 'concat' file.
        """
        for sample in self.samples:
            if len(sample.files.fastqs) == 1:
                sample.files.concat = sample.files.fastqs

            else:
                logger.debug(
                    "concatenating inputs: {}".format(sample.files.fastqs))

                # cat command works for gzip or not. Index 0 of tuples is R1s.
                cmd1 = ["cat"] + [i[0] for i in sample.files.fastqs]

                isgzip = ".gz"
                if not sample.files.fastqs[0][0].endswith(".gz"):
                    isgzip = ""

                # write to new concat handle
                conc1 = os.path.join(
                    self.data.dirs.edits, 
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
                        self.data.dirs.edits, 
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
                sample.files.concat = [(conc1, conc2)]
                logger.debug('concat: {}'.format(sample.files.concat))


    def distribute_fastp(self):
        """
        Distribute N / 2 fastp jobs on N processors.
        """
        # send samples to cutadapt filtering
        logger.debug("submitting remote fastp jobs")
        jobs = {}
        for sample in self.samples:
            rasync = self.lbview.apply(trim_reads, *(self.data, sample))
            jobs[sample.name] = rasync

        # wait for all to finish
        printstr = ("processing reads    ", "s2")
        prog = AssemblyProgressBar(jobs, None, printstr, self.data)
        prog.block()
        prog.check()
        
        # check for errors
        for sample in self.samples:

            # enter stats on the sample
            sname = sample.name
            jdata, filepaths = jobs[sname].get()

            # parse from json to sample stats
            sample.stats_dfs.s2["reads_raw"] = int(
                jdata['summary']['before_filtering']['total_reads'])
            sample.stats_dfs.s2["reads_filtered_by_Ns"] = int(
                jdata['filtering_result']['too_many_N_reads'])
            sample.stats_dfs.s2['reads_filtered_by_low_quality'] = int(
                jdata['filtering_result']['low_quality_reads'])
            sample.stats_dfs.s2['reads_filtered_by_low_complexity'] = int(
                jdata['filtering_result']['low_complexity_reads'])
            sample.stats_dfs.s2['reads_filtered_by_minlen'] = int(
                jdata['filtering_result']['too_short_reads'])
            sample.stats_dfs.s2['avg_len_R1_before_trimming'] = int(
                jdata['summary']['before_filtering']['read1_mean_length'])
            sample.stats_dfs.s2['avg_len_R1_after_trimming'] = int(
                jdata['summary']['after_filtering']['read1_mean_length'])
            sample.stats_dfs.s2["reads_passed_filter"] = int(
                jdata['summary']['after_filtering']['total_reads'])

            # R2 optional
            sample.stats_dfs.s2['avg_len_R2_before_trimming'] = 0
            sample.stats_dfs.s2['avg_len_R2_after_trimming'] = 0
            if "read2_mean_length" in jdata['summary']['before_filtering']:
                sample.stats_dfs.s2['avg_len_R2_before_trimming'] = int(
                    jdata['summary']['before_filtering']['read2_mean_length'])
                sample.stats_dfs.s2['avg_len_R2_after_trimming'] = int(
                    jdata['summary']['after_filtering']['read2_mean_length'])

                for key in [
                    "reads_raw", 
                    "reads_filtered_by_Ns", 
                    "reads_filtered_by_minlen", 
                    "reads_filtered_by_low_complexity", 
                    "reads_filtered_by_low_quality", 
                    "reads_passed_filter"
                    ]:
                    sample.stats_dfs.s2[key] = sample.stats_dfs.s2[key] / 2

            sample.stats_dfs.s2 = sample.stats_dfs.s2.astype(int)
            logger.info(sample.stats_dfs.s2)

            # save stats to Assembly
            if sample.stats_dfs.s2.reads_passed_filter:
                sample.stats.state = 2
                sample.stats.reads_passed_filter = (
                    sample.stats_dfs.s2.reads_passed_filter)
                sample.files.edits = filepaths

            else:
                logger.warning("No reads passed filtering in {}".format(sname))

        # get step2 stats dataframe for Sample objects and stats filepath
        self.data.stats_dfs.s2 = self.data._build_stat("s2")
        self.data.stats_files.s2 = os.path.join(
            self.data.dirs.edits, 
            's2_rawedit_stats.txt')

        # write stats for all samples
        with open(self.data.stats_files.s2, 'w') as outfile:
            (
                self.data.stats_dfs.s2.fillna(value=0)
                .astype(int)
                .to_string(outfile)
            )


    def assembly_cleanup(self):
        """
        discard tmp concat files created prior to trimming. These are
        only created if you merged assemblies prior to step2 such that
        a sample had multiple input fastq files/pairs.
        """
        logger.warning("TODO: add concat cleanup")




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
        self.paired = read2 != ""
        self.fastp_binary = os.path.join(sys.prefix, "bin", "fastp")

        # output file paths (do not bother gzipping since these are tmp files)
        self.out1 = os.path.join(self.workdir, sname + ".trimmed_R1_.fastq.gz")
        self.out2 = os.path.join(self.workdir, sname + ".trimmed_R2_.fastq.gz")

        # paths to stats files
        self.json = os.path.join(self.workdir, '{}.json'.format(sname))
        self.html = os.path.join(self.workdir, '{}.html'.format(sname))


    def trim_reads(self):
        """
        Calls fastp in subprocess and writes tmpfiles to workdir.
        0 = no trim or quality check.
        1 = only trimming
        2 = trimming and quality checks.

        (-A turns off adapter trimming)
        (-Q turns off quality filtering)
        """
        if self.paired:
            cmd = [
                self.fastp_binary, 
                "-i", self.read1,
                "-I", self.read2,
                "-o", self.out1,
                "-O", self.out2,
            ]

        else:
            cmd = [
                self.fastp_binary,
                "-i", self.read1,
                "-o", self.out1,
            ]

        cmd.extend([
            "-5",  # sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-3",  # sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
            "-a", str(self.data.hackersonly.p5_adapter),
            "-q", str(20 + self.data.params.phred_Qscore_offset - 33),
            "-l", str(self.data.params.filter_min_trim_len),
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

        # turn off filters
        if self.data.params.filter_adapters == 1:
            cmd.extend("-A")
        if self.data.params.filter_adapters == 0:
            cmd.extend("-Q")

        # logger record
        logger.debug(" ".join(cmd))

        # run the command
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
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





def trim_reads(data, sample):
    """
    Remove function for applying fastp read trimming.
    """
    tool = ReadTrimming(
        data=data,
        sname=sample.name,
        read1=sample.files.concat[0][0],
        read2=sample.files.concat[0][1],
        workdir=data.dirs.edits,
    )
    tool.trim_reads()
    return tool.parse_stats_from_json()




if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    CURDIR = os.path.dirname(__file__)
    SIM_PREFIX = os.path.join(CURDIR, "../../tests/ipsimdata/")

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
    tdata = ip.load_json("/tmp/test-simrad.json")    
    tdata = ip.merge('merge-s2-test', [tdata, tdata])
    tdata.params.filter_adapters = 2
    tdata.run("2", auto=True, force=True)
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

