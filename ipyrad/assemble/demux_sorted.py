#!/usr/bin/env python

"""Loads files and counts nreads for sorted_fastq_path data input
"""

from typing import List, Dict, Tuple, TypeVar
from pathlib import Path
import itertools
import subprocess
import pandas as pd
from loguru import logger
from ipyrad.assemble.utils import IPyradError, BADCHARS
from ipyrad.schema import Sample, Stats1
from ipyrad.core.progress import AssemblyProgressBar

logger = logger.bind(name="ipyrad")
Step1 = TypeVar("Step1")
Assembly = TypeVar("Assembly")


class FileLinker:
    """Loads Samples from file names and check sample names for bad chars.
    """
    def __init__(self, step: Step1):
        self.data: Assembly = step.data
        self.quiet: bool = step.quiet
        self.lbview = step.ipyclient.load_balanced_view()

        self.fastqs: List[Path] = []
        self.ftuples: Dict[str,Tuple[str,str]] = {}

    def run(self):
        """Checks files and then counts reads in each one."""
        self._get_fastq_files()
        self._check_file_suffixes()
        self._parse_sample_names_and_pairs()
        self._check_sample_names()
        self._remote_run_counter()
        self._write_stats()
        self.data.save_json()

    def _get_fastq_files(self):
        """Sets .fastqs with all Paths matching the regular expression path."""
        fastq_path = self.data.params.sorted_fastq_path
        self.fastqs = list(fastq_path.parent.glob(fastq_path.name))

    def _check_file_suffixes(self):
        """Check file endings for unsupported types."""
        # Assert files are not .bz2 format
        if any(i.suffix == ".bz2" for i in self.fastqs):
            raise IPyradError(BZ2_ERROR)

        # filter out any files without proper file endings. Raise if None left
        endings = (".gz", ".fastq", ".fq")
        keep = []
        for path in self.fastqs:
            if path.suffix not in endings:
                logger.warning(f"skipping {path}, file suffix not recognized.")
            else:
                keep.append(path)

        # update fastqs to be files with good suffix, raise if no files.
        self.fastqs = keep
        if not self.fastqs:
            raise IPyradError((
                "No fastq files found in 'sorted_fastq_path': {}\n"
                "Check that file names match the required convention for "
                "paired datatype, i.e., paired file names should be "
                "identical save for 'R1' and 'R2' in names."
                .format(self.data.params.sorted_fastq_path))
            )

    def _parse_sample_names_and_pairs(self):
        """Check for PE matching and extract sample names to ftuples.

        A fairly complicated process to flexibly find matching PE files.
        """
        def drop_from_right(fname: Path, idx: int) -> str:
            """Used in pair name parsing to sub out _ delim chunks
            
            Example
            -------
            if idx = 1 then:
                name_prefix_001_R1_002.fastq.gz           # 1.
                ['name_prefix', '001', 'R1', '002']       # 2.
                ['name_prefix', '001', '002']             # 3.
            """
            # fname = os.path.basename(fname).rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            fname = fname.name.rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            chunks = fname.split("_")
            sublist = [j for i, j in enumerate(chunks[::-1]) if i != idx][::-1]
            return "_".join([i for i in sublist if i]).rstrip("_")

        # check if file names end with _ before the suffix and split
        # on two underscores, else split on last one.
        bases = sorted([
            i.name.rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            for i in self.fastqs
        ])
        logger.debug(f"names: {bases}")

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

            for fname, paths in groups:
                paths = sorted(paths)
                logger.debug(f"detected paired files: {paths}")
                self.ftuples[fname] = (paths[0], paths[1])

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
                if any(base.endswith(i) for i in endings):
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
                # fname = os.path.basename(i.rstrip('.gz').rstrip('.fastq').rstrip('.fq'))
                fname = i.name.rstrip('.gz').rstrip('.fastq').rstrip('.fq')
                self.ftuples[fname] = (i, "")
        logger.debug(f"paired paths: {self.ftuples}")

    def _check_sample_names(self):
        """Do not allow bad characters in names."""
        snames = sorted(self.ftuples)
        for sname in snames:
            if any(i in sname for i in BADCHARS):
                newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(
                    f"changing name {sname} to {newname} (hard characters).")
                self.ftuples[newname] = self.ftuples.pop(sname)

    def _remote_run_counter(self):
        """Read in fastq files, count nreads for stats, and create Samples."""
        # submit jobs to run on client
        jobs = {}
        for sname, ftup in self.ftuples.items():
            jobs[sname] = self.lbview.apply(zbuf_count_lines, ftup[0])

        # show progress bar until all jobs complete
        msg = "loading reads"
        prog = AssemblyProgressBar(jobs, msg, step=1, quiet=self.quiet)
        prog.block()
        prog.check()

        # collect results as they finish and show progress bar
        for sname, nreads in prog.results.items():
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

    def _write_stats(self):
        """Write stats to human-readable stats file."""
        logger.info(f"created {len(self.data.samples)} new samples")
        statsdf = pd.DataFrame(
            index=sorted(self.data.samples),
            columns=["reads_raw"],
            data=[
                self.data.samples[i].stats_s1.reads_raw 
                for i in sorted(self.data.samples)
            ])
        handle = self.data.stepdir / 's1_demultiplex_stats.txt'
        with open(handle, 'w', encoding="utf-8") as out:
            statsdf.fillna(value=0).astype(int).to_string(out)
            # statsdf.to_csv(handle, sep="\t")


def zbuf_count_lines(filename):
    """Fast line counter using unix utils."""
    # gunzip -c is the same as zcat but supported on more systems
    if filename.suffix == ".gz":
        cmd1 = ["gunzip", "-c", filename]
    else:
        cmd1 = ["cat", filename]

    # counts lines from gunzip stream
    cmd2 = ["wc", "-l"]

    # start process one and pipe to process 2
    with subprocess.Popen(cmd1, stdout=subprocess.PIPE) as proc1:
        with subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE) as proc2:
            res = proc2.communicate()[0]
    if proc2.returncode:
        raise IPyradError("error zbuf_count_lines {}:".format(res))
    nreads = int(res.split()[0]) / 4
    return nreads


BZ2_ERROR = """
Found bz2 formatted files in 'sorted_fastq_path'.
ipyrad does not support bz2 files. The only supported
formats for samples are .gz, .fastq, and .fq. The easiest
thing to do is probably go into your sorted_fastq_path
directory and issue this command `bunzip2 *`. You
will probably also need to update your params file to
reflect the fact that sample raw files now probably end
with .fq or .fastq.
"""


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    DATA = ip.Assembly("TEST")
    DATA.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    DATA.params.project_dir = "/tmp"
    DATA.run('1', force=True, quiet=True)
    print(DATA.stats)
