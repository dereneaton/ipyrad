#!/usr/bin/env python

"""sra tools wrapper to download archived sequence data.

This requires that you have 'sra-tools' > v.2.11 installed.
Unfortunately, the conda installer for this is very unreliable,

Example
-------
>>> tool = ipa.sratools(accessions=["SRP065788"])
>>> metadata = tool.fetch_runinfo(fields=[1, 4, 6, 28, 29, 30])
>>> print(metadata)
>>> tool.run()
"""

from typing import Sequence, Optional
import sys
import time
import subprocess as sps
from io import StringIO
from pathlib import Path

# third party
import requests
import pandas as pd
from loguru import logger
from ipyrad.core.cluster import Cluster
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import IPyradError

logger = logger.bind(name="ipa")

# raise warning if missing imports
MISSING_IMPORTS = """
To use the ipa.sratools module you must install the sra toolkit to a
version >2.11 (2.8 will not work!). You can do this either by visting
the ncbi download page: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
or by using conda. Unfortunately the conda recipe for this is quite
unstable and often won't update to a recent enough version, but you
can try with this command:

  conda install 'sra-tools>=2.11' -c conda-forge -c bioconda

"""


ACCESSION_ID = """
Accession ID must be either a Run or Study accession, i.e.,
it must have one the following prefixes:
  Study: SRR, ERR, DRR
  Project: SRP, ERP, DRP
"""


class SRA:
    """Tool for downloading data from SRA, ERA, or DRA archives.

    Example
    -------
    >>> tool = ipa.sratools(accessions="SRP065788")
    >>> metadata = tool.fetch_runinfo(fields=[1, 4, 6, 28, 29, 30])
    >>> tool.run()

    Parameters
    ----------
    acccessions: Sequence[str]
        One or more accessions IDs from SRA, ERA, etc.
    workdir: str
        A path str to the directory where outputs should be written.
    name_fields: int or str
        Provide the index (1-indexed) of the name fields to be used
        as a prefix for fastq output files. The default is (1,30),
        which is the accession + SampleName fields. Use
        sra.fetch_fields to see all available fields and their
        indices. If multiple are listed then they will be joined
        by a "_" character. For example (29,30) would yield something
        like: latin-name_sample-name (e.g., mus_musculus-NR10123).
    dry_run: bool
        If True then a table of file names that _would_ be downloaded
        will be shown, but the actual files will note be downloaded.
    split_pairs: bool or None
        If True then pairs are split, if False they are not split, if
        None then we will auto-detect if paired or not and split pairs
        when detected. Forcing splitting can be helpful when the
        metadata was not set properly.
    binary_path: bool or None
        If None then binaries are looked for in the conda path. Else,
        user can override by entering a path to a `fastq-dump` binary.
    gzip: bool
        Gzip compress fastq files.
    """
    def __init__(
        self,
        accessions: Sequence[str],
        workdir: str = "sra-fastq-data",
        name_fields: Sequence[int] = (1, 30),
        name_separator: str = "_",
        dry_run: bool = False,
        split_pairs: Optional[bool] = None,
        binary_path: Optional[str] = None,
        ):

        # store attributes
        self.accessions = accessions
        self.workdir = Path(workdir).expanduser().absolute()
        self.name_fields = name_fields
        self.name_separator = name_separator
        self.dry_run = dry_run
        self.split_pairs = split_pairs
        self.binary_path = (
            Path(sys.prefix) / "bin" if binary_path is None else
            Path(binary_path).expanduser().absolute().parent)
        """: Prefix directory where fastq-dump and vdb-config binaries are"""
        self.sra_tmpdir = None
        """: Auto-set by check_binaries function"""
        self.is_sample = False
        self.is_project = False

        # if accession is a list then make it comma separated string
        if isinstance(self.accessions, (list, tuple)):
            pass
        if isinstance(self.accessions, (str)):
            self.accessions = [self.accessions]

        # get type
        if any(i in self.accessions[0] for i in ["SRR", "ERR", "DRR"]):
            self.is_sample = True
        elif any(i in self.accessions[0] for i in ["SRP", "ERP", "DRP"]):
            self.is_project = True
        else:
            raise IPyradError(ACCESSION_ID)

        # make sure required software if installed
        self.check_binaries()
        self.check_vdb()

    def check_binaries(self):
        """Find the fastq-dump binary in user $PATH."""
        binary = self.binary_path / "fastq-dump"
        cmd = ['which', str(binary)]
        with sps.Popen(cmd, stdout=sps.PIPE) as proc:
            comm = proc.communicate()[0].decode()
            if not comm:
                raise IPyradError(MISSING_IMPORTS)

    def check_vdb(self):
        """Check that user has run vdb-config, and report the tmpdir
        """
        # try calling fastq-dump and catch not-configured error
        binary = self.binary_path / "vdb-config"
        cmd = [str(binary), "-h"]
        with sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE) as proc:
            out, _ = proc.communicate()
            if "has not been configured" in out.decode():
                raise IPyradError(
                "To run sratools you must first configure the toolkit "
                "by running 'vdb-config -i' (a requirement of the sra-tools "
                "software).\n\nWhen doing so take note of the 'location of "
                "user-repository' setting which is the location where tmp "
                "files will be stored while downloading. Set this to a "
                "place with sufficient disk space. ipa.sratools will "
                "automatically remove any tmp files in this dir after it is "
                "finished using it."
            )
        try:
            vdb_path = Path("~/.ncbi/user-settings.mkfg").expanduser()
            with open(vdb_path, 'r', encoding="utf-8") as indata:
                lines = indata.readlines()
                for line in lines:
                    if line.startswith("/repository/user/default-path"):
                        self.sra_tmpdir = line.strip().split()[-1].strip('"')
        except Exception:
            pass
        if self.sra_tmpdir:
            self.sra_tmpdir = Path(self.sra_tmpdir).absolute()
            logger.info(f"vdb-config path (tmpdir) is {self.sra_tmpdir}")

    def run(self, cores: int=None):
        """Download accessions as fastq files into a designated workdir.

        Parameters
        ----------
        cores: int
            Number of cores for downloading files in parallel.
        """
        # ensure outdir exists
        self.workdir.mkdir(exist_ok=True)
        with Cluster(cores=cores) as ipyclient:

            # parallelize downloads
            lbview = ipyclient.load_balanced_view()

            # get run info and sort so largest samples are on top
            rundf = self.fetch_runinfo(list(range(31)))
            rundf = rundf.sort_values(
                by="spots",
                ascending=False,
            ).reset_index(drop=True)

            # make empty Accession field
            rundf["Accession"] = ""

            # choose spacer to replace spaces in names as different from name_sep
            otherspacer = "_" if self.name_separator != "_" else "-"

            # select names for downloaded .sra files
            if self.name_fields:

                # indices of runinfo fields for names
                fields = [i - 1 for i in fields_checker(self.name_fields)]

                # set new accession name
                for row in rundf.index:
                    rundf.loc[row, "Accession"] = (
                        self.name_separator.join(
                            [rundf.iloc[row, i] for i in fields]
                            )
                        ).replace(" ", otherspacer)

            # backup default naming scheme
            else:
                if rundf.SampleName.value_counts().max() > 1:
                    # set new accession name
                    for row in rundf.index:
                        rundf.loc[row, "Accession"] = (
                            self.name_separator.join(
                                [rundf.iloc[row, i] for i in [30, 1]]
                                )
                            )
                else:
                    rundf.Accession = rundf.SampleName

            # test run to see file names and location without download
            if self.dry_run:
                logger.info(
                    f"The following files will be written to: {self.workdir}\n"
                    f"{rundf.Accession}\n"
                )
                return

            # send download jobs
            msg = "downloading fastq data"
            jobs = {}
            for sidx in rundf.index:
                acc = rundf.Accession[sidx]
                srr = rundf.Run[sidx]
                jobs[sidx] = lbview.apply(
                    self._call_fastq_dump_on_srr, *(acc, srr, self.split_pairs)
                )
            prog = AssemblyProgressBar(jobs, msg, f"{len(jobs)} samples", False)
            prog.block()
            prog.check()

            # final report
            logger.info(f"\n{len(jobs)} fastq files downloaded to {self.workdir}")

    @property
    def fetch_fields(self):
        """The column names (fields) in an SRA Run Table."""
        fields = pd.DataFrame(COLNAMES, columns=["field"])
        fields.index += 1
        return fields

    @property
    def fields(self):
        """The column names (fields) in an SRA Run Table"""
        fields = pd.DataFrame(COLNAMES, columns=["field"])
        fields.index += 1
        return fields

    def fetch_runinfo(self, fields=None):
        """
        Query the RunInfo for a Sample or Run, returned as a DataFrame.
        The fields can be subselected. See <self>.fields for options.
        """
        logger.info("Fetching project data...")

        if fields is None:
            fields = list(range(31))
        fields = fields_checker(fields)

        sra_ids = []
        for accession in self.accessions:
            # SRA IDs from the SRP
            res = requests.get(
                url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={
                    "db": "sra",
                    "term": accession,
                    "tool": "ipyrad",
                    "email": "ipyrad@gmail.com",
                    "retmax": 1000,
                    },
                )
            sra_ids += [i[4:-5] for i in res.text.split() if "<Id>" in i]
            if not sra_ids:
                raise IPyradError(
                    "No SRA samples found in {}"
                    .format(self.accessions))
            time.sleep(3)

        # Get SRA Runinfo in batches of 20 at a time b/c who knows...
        blocks = []
        for block in range(0, len(sra_ids), 20):

            res = requests.get(
                url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={
                    "db": "sra",
                    "id": ",".join(sra_ids[block:block + 20]),
                    "tool": "ipyrad",
                    "email": "ipyrad@gmail.com",
                    "rettype": "runinfo",
                    "retmode": "text",
                    },
                )
            time.sleep(3)
            rundf = pd.read_csv(StringIO(res.text.strip()))
            blocks.append(rundf)

        rundf = pd.concat(blocks)
        rundf.reset_index(drop=True, inplace=True)
        return rundf.iloc[:, [i - 1 for i in fields]]


    def _call_fastq_dump_on_srr(self, acc, srr, paired):
        """
        calls fastq-dump on SRRs, relabels fastqs by their accession
        names, and writes them to the workdir. Saves temp sra files
        in the designated vdb config tmp folder and immediately
        removes them.
        """
        # build command for fastq-dumping
        fd_cmd = [
            str(self.binary_path / "fastq-dump"), srr,
            "--accession", acc,
            "--outdir", self.workdir,
            # "--disable-multithreading",
            ]
        if paired:
            fd_cmd += ["--split-spot"]

        # call fq dump command
        with sps.Popen(fd_cmd, stderr=sps.STDOUT, stdout=sps.PIPE) as proc:
            out, _ = proc.communicate()
            if proc.returncode:
                raise IPyradError(out.decode())

        # delete the temp sra file from the place
        sra_file = self.sra_tmpdir / srr + ".sra"
        sra_file.unlink(missing_ok=True)

def download_file(url, outname):
    """NOTE the stream=True parameter"""
    res = requests.get(url, stream=True)
    with open(outname, 'wb') as fout:
        for chunk in res.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                fout.write(chunk)
    return outname

def fields_checker(fields):
    """returns a fields argument formatted as a list of strings.
    and doesn't allow zero.
    """
    # make sure fields will work
    if isinstance(fields, int):
        fields = str(fields)

    if isinstance(fields, str):
        if "," in fields:
            fields = [str(i) for i in fields.split(",")]
        else:
            fields = [str(fields)]

    elif isinstance(fields, (tuple, list)):
        fields = [str(i) for i in fields]

    else:
        raise IPyradError("fields not properly formatted")

    ## do not allow zero in fields
    fields = [int(i) for i in fields if i != '0']
    return fields



## OMG FASTERQ-DUMP IS THE WORST DON"T TRY THIS!!!
# def fasterq_dump_file(path):
#     "Call fasterq-dump multi-threaded"

#     # get fastq conversion path
#     path = os.path.realpath(path)
#     fastqpath = path.rsplit(".sra", 1)[0] + ".fastq"

#     # call fasterq dump in a subprocess to write to .fastq
#     cmd = [
#         "fasterq-dump", path,
#         "-o", fastqpath,
#         "-t", os.path.join(tempfile.gettempdir(), "scratch"),
#         "--split-files",
#     ]
#     print("\n\n" + " ".join(cmd) + "\n\n")
#     proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
#     res = proc.communicate()

#     # check for errors
#     if proc.returncode:
#         raise IPyradError("error in fasterq-dump:\n{}\n{}"
#             .format(" ".join(cmd), res[0].decode())
#             )

#     # rename files in ipyrad format of "_R1_, _R2_"


#     # remove .sra file
#     print(path, fastqpath)
#     # os.remove(path)


FAILED_DOWNLOAD = """
Warning: One or more files failed to finish downloading or converting to fastq.
To avoid corruption the file was removed. Try downloading again to get
any missing files. The following samples were affected:
{}
"""


COLNAMES = [
    'Run',
    'ReleaseDate',
    'LoadDate',
    'spots',
    'bases',
    'spots_with_mates',
    'avgLength',
    'size_MB',
    'AssemblyName',
    'download_path',
    'Experiment',
    'LibraryName',
    'LibraryStrategy',
    'LibrarySelection',
    'LibrarySource',
    'LibraryLayout',
    'InsertSize',
    'InsertDev',
    'Platform',
    'Model',
    'SRAStudy',
    'BioProject',
    'Study_Pubmed_id',
    'ProjectID',
    'Sample',
    'BioSample',
    'SampleType',
    'TaxID',
    'ScientificName',
    'SampleName',
    'g1k_pop_code',
    'source',
    'g1k_analysis_group',
    'Subject_ID',
    'Sex',
    'Disease',
    'Tumor',
    'Affection_Status',
    'Analyte_Type',
    'Histological_Type',
    'Body_Site',
    'CenterName',
    'Submission',
    'dbgap_study_accession',
    'Consent',
    'RunHash',
    'ReadHash',
]


if __name__ == "__main__":

    # TEST
    # conda install 'sra-tools>=2.11' -c conda-forge -c bioconda
    import ipyrad.analysis as ipa
    ipa.set_log_level("DEBUG")

    # tool = ipa.sratools(accessions="SRP065788")
    tool = ipa.sratools(accessions=["SRR2895756", "SRR2895757"])
    metadata = tool.fetch_runinfo(fields=[1, 4, 6, 28, 29, 30])
    tool.run()
    # ipa.sratools(...)
