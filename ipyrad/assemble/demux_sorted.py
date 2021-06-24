#!/usr/bin/env python

"""
Loads files and counts nreads for sorted_fastq_path data input
"""

import os
import glob
import itertools
import subprocess
import pandas as pd
from loguru import logger
from ipyrad.assemble.utils import IPyradError, BADCHARS
from ipyrad.core.schema import SampleSchema, Stats1
from ipyrad.core.progress_bar import AssemblyProgressBar



class FileLinker:
    """
    Loads Samples from file names and check sample names for bad chars.
    """
    def __init__(self, step):
        # store input
        self.step = step
        self.quiet = step.quiet
        self.data = step.data
        self.fastqs = glob.glob(self.data.params.sorted_fastq_path)

        # to be filled
        self.ftuples = {}

        # parallel distributor
        self.lbview = step.ipyclient.load_balanced_view()


    def run(self):
        """
        checks files and then counts reads in each one.
        """
        self.check_file_suffixes()
        self.parse_sample_names_and_pairs()
        self.check_sample_names()

        # distributes jobs to parallel
        self.remote_run_counter()

        # save json file and csv stats file.
        self.write_json_and_stats()


    def check_file_suffixes(self):
        """
        Check file endings for unsupported types.
        """
        # Assert files are not .bz2 format
        if any([i.endswith(".bz2") for i in self.fastqs]):
            raise IPyradError((
                "Found bz2 formatted files in 'sorted_fastq_path': {} "
                "ipyrad does not support bz2 files. The only supported "
                "formats for samples are .gz, .fastq, and .fq. The easiest "
                "thing to do is probably go into your sorted_fastq_path "
                "directory and issue this command `bunzip2 *`. You "
                "will probably also need to update your params file to "
                "reflect the fact that sample raw files now probably end "
                "with .fq or .fastq.")
                .format(self.fastqs))

        # filter out any files without proper file endings. Raise if None left
        endings = ("gz", "fastq", "fq")
        fastqs = []
        for i in self.fastqs:
            if i.split(".")[-1] not in endings:
                logger.warning(f"skipping {i}, file suffix not recognized.")
            else:
                fastqs.append(i)
        self.fastqs = fastqs
        if not self.fastqs:
            raise IPyradError((
                "No fastq files found in 'sorted_fastq_path': {}\n"
                "Check that file names match the required convention for "
                "paired datatype, i.e., paired file names should be "
                "identical save for _R1_ and _R2_ (note the underscores "
                "before AND after R*).")
                .format(self.data.params.sorted_fastq_path))


    def parse_sample_names_and_pairs(self):
        """
        check for PE matching and extract sample names to ftuples
        """
        def drop_from_right(fname, idx):
            """
            used in pair name parsing to sub out _ delim chunks
            
            Example:
            ---------
            if idx = 1 then:
                name_prefix_001_R1_002.fastq.gz           # 1.
                ['name_prefix', '001', 'R1', '002']       # 2.
                ['name_prefix', '001', '002']             # 3.
            """
            fname = os.path.basename(fname).rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            chunks = fname.split("_")
            sublist = [j for i, j in enumerate(chunks[::-1]) if i != idx][::-1]
            return "_".join([i for i in sublist if i]).rstrip("_")

        # check if file names end with _ before the suffix and split
        # on two underscores, else split on last one.
        bases = sorted([
            os.path.basename(i).rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            for i in self.fastqs
        ])
        logger.info(bases)

        # link pairs into tuples
        if self.data.is_pair:

            # try these in order until ngroups == nfiles / 2
            idx = 0
            while 1:
                try:
                    # get groups up to an underscore delimiter
                    groups = itertools.groupby(
                        sorted(self.fastqs),
                        key=lambda x: drop_from_right(x, idx),
                    )
                    groups = {i: list(j) for i, j in groups}
                    logger.debug(f"{idx} {groups}")
                
                    assert len(groups) == len(self.fastqs) / 2
                    assert all(len(j) == 2 for i, j in groups.items())
                    logger.debug(f"using '_' rsplit = {idx}")
                    break
                except Exception:
                    pass
                # increase counter up to 5 _ back from end, then raise 
                idx += 1
                if idx > 5:
                    raise IPyradError(
                        "Cannot parse paired file names. File names must have "
                        "matching name prefix followed by _1 _2, _R1 _R2, "
                        "or _R1_ _R2_ followed by any subsequent suffix. "
                        f"Your filenames look like this: {self.fastqs}"
                        )

            # apply splitter to the full path names
            groups = itertools.groupby(
                sorted(self.fastqs),
                key=lambda x: drop_from_right(x, idx),
            )

            for fname, files in groups:
                fname = os.path.basename(fname)
                files = sorted(files)
                logger.debug(f"detected paired files: {files}")
                self.ftuples[fname] = (files[0], files[1])

            # file checks
            if not self.ftuples:
                raise IPyradError(
                    "No paired fastq files found. File names must have "
                    "matching name prefix followed by _1 _2, _R1 _R2, "
                    "or _R1_ _R2_.")

        # data are not paired, create empty tuple pair
        else:
            # print warning if _R2_ is in names when not paired
            endings = ("_R2", "_2", "_R2")
            warning = []
            for base in bases:
                if any([base.endswith(i) for i in endings]):
                    warning.append(base)
            if warning:
                message = (
                    "Input file names look suspiciously like paired-end"
                    "data but you selected single end. If so, you should "
                    "set the parameter 'datatype' to a paired option (e.g., "
                    "pairddrad or pairgbs) and re-run step 1, which will "
                    "require using the force flag (-f) to overwrite "
                    "existing data.\n{}".format(",".join(warning)))
                logger.warning(message)

            for i in self.fastqs:
                fname = os.path.basename(i.rstrip('.gz').rstrip('.fastq').rstrip('.fq'))
                self.ftuples[fname] = (i, "")

        logger.debug(f"ftuples: {self.ftuples}")


    def check_sample_names(self):
        """
        Do not allow bad characters in names.
        """
        snames = sorted(self.ftuples)
        for sname in snames:
            if any([i in sname for i in BADCHARS]):
                newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(
                    f"changing name {sname} to {newname} (hard characters).")
                self.ftuples[newname] = self.ftuples.pop(sname)


    def remote_run_counter(self):
        """
        Read in fastq files, count nreads for stats, and create Samples.
        """
        # submit jobs to run on client
        rasyncs = {}
        for sname in self.ftuples:
            rasyncs[sname] = self.lbview.apply(
                zbuf_count_lines, self.ftuples[sname][0]
            )

        # show progress bar until all jobs complete
        message = "loading reads"
        prog = AssemblyProgressBar(rasyncs, message, step=1, quiet=self.quiet)
        prog.block()
        prog.check()

        # collect results as they finish and show progress bar
        for sname in prog.results:
            # get nreads
            nreads = prog.results[sname]
            if not nreads:
                logger.warning(
                    f"sample {sname} has 0 reads and will be excluded.")
            else:
                # create a new Sample w/ stats and files
                stats_s1 = Stats1(reads_raw=nreads)
                sample = SampleSchema(name=sname, stats_s1=stats_s1)
                sample.files.fastqs = [self.ftuples[sname]]
                self.data.samples[sname] = sample
                logger.debug(f"new sample {sname}: {sample.stats_s1.reads_raw}")


    def write_json_and_stats(self):
        """
        Write updates to the project JSON files and copy stats to a 
        human-readable stats file.
        """
        logger.info(f"created {len(self.data.samples)} new samples")
        self.data.save_json()

        handle = os.path.join(self.step.stepdir, 's1_demultiplex_stats.txt')
        statsdf = pd.DataFrame(
            index=sorted(self.data.samples),
            columns=["reads_raw"],
            data=[
                self.data.samples[i].stats_s1.reads_raw 
                for i in sorted(self.data.samples)
            ])
        with open(handle, 'w') as outfile:
            statsdf.fillna(value=0).astype(int).to_string(outfile)
            # statsdf.to_csv(handle, sep="\t")


def zbuf_count_lines(filename):
    """
    Fast line counter using unix utils
    """
    # gunzip -c is the same as zcat but supported on more systems
    if filename.endswith(".gz"):
        cmd1 = ["gunzip", "-c", filename]
    else:
        cmd1 = ["cat", filename]

    # counts lines from gunzip stream
    cmd2 = ["wc", "-l"]

    # start process one and pipe to process 2
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE)
    res = proc2.communicate()[0]
    if proc2.returncode:
        raise IPyradError("error zbuf_count_lines {}:".format(res))
    nreads = int(res.split()[0]) / 4
    return nreads


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    DATA = ip.Assembly("TEST")
    DATA.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    DATA.params.project_dir = "/tmp"
    DATA.run('1', force=True, quiet=True)
    print(DATA.stats)

