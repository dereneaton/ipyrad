#!/usr/bin/env python

"sra tools wrapper to download archived seq data"

# py2/3
from __future__ import print_function
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

# standard
import os
import time
import requests
import subprocess as sps

# third party
import pandas as pd
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import progressbar


## raise warning if missing imports
MISSING_IMPORTS = """
To use the ipa.sratools module you must install the sra-tools
software, which you can do with the following conda command. 

  conda install bioconda::sra-tools
"""

ACCESSION_ID = """
Accession ID must be either a Run or Study accession, i.e., 
it must have one the following prefixes:
  Study: SRR, ERR, DRR
  Project: SRP, ERP, DRP
"""


class SRA(object):
    """ ipyrad.analysis SRA download object"""
    def __init__(
        self, 
        accessions,
        workdir="sra-fastq-data",
        ):

        # store attributes
        self.accessions = accessions
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        self.is_sample = False
        self.is_project = False
        self._oldtmpdir = None

        ## cluster attributes
        self.ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }

        # if accession is a list then make it comma separated string
        if isinstance(self.accessions, (list, tuple)):
            pass
        if isinstance(self.accessions, (str)):
            self.accessions = [self.accessions]

        # get type
        if any([i in self.accessions[0] for i in ["SRR", "ERR", "DRR"]]):        
            self.is_sample = True
        elif any([i in self.accessions[0] for i in ["SRP", "ERP", "DRP"]]):
            self.is_project = True
        else:
            raise IPyradError(ACCESSION_ID)

        # make sure required software if installed
        self.check_binaries()


    def check_binaries(self):
        # check imports
        for binary in ['fastq-dump', 'fasterq-dump']:
            proc = sps.Popen(['which', binary], stdout=sps.PIPE)
            comm = proc.communicate()[0]
            if not comm:
                raise IPyradError(MISSING_IMPORTS)


    def run(
        self, 
        ipyclient=None,
        force=False, 
        name_fields=(1, 30), 
        name_separator="_", 
        dry_run=False, 
        split_pairs=None,
        gzip=False,
        show_cluster=False,
        auto=False,
        ):
        """
        Download the accessions into a the designated workdir. 

        Parameters
        ----------
        ipyclient: (ipyparallel.Client)
            If provided, work will be distributed across a parallel
            client, otherwise download will be run on a single core.

        force: (bool)
            If force=True then existing files with the same name
            will be overwritten. 

        name_fields: (int, str):
            Provide the index of the name fields to be used as a prefix
            for fastq output files. The default is 30, which is the 
            SampleName field. Use sra.fetch_fields to see all available
            fields and their indices. A likely alternative is 1 (Run). 
            If multiple are listed then they will be joined by a "_" 
            character. For example (29,30) would yield something like:
            latin-name_sample-name (e.g., mus_musculus-NR10123).

        dry_run: (bool)
            If True then a table of file names that _would_ be downloaded
            will be shown, but the actual files will note be downloaded.

        split_pairs: (bool or None)
            If True then pairs are split, if False they are not split, if None
            then we will auto-detect if paired or not.

        gzip: bool
            Gzip compress fastq files.

        auto: bool
            Automatically launch new ipcluster for parallelization and 
            shutdown when finished. See <object>.ipcluster for settings.
        """
        # ensure output directory, also used as tmpdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # distribute job wrapped in ipcluster cleanup
        pool = Parallel()
        # parallelize_run(
        #     self,
        #     run_kwargs={
        #         "ipyclient": ipyclient,
        #         "force": force,
        #         "name_fields": name_fields,
        #         "name_separator": name_separator,
        #         "dry_run": dry_run,
        #         "split_pairs": split_pairs,
        #         "gzip": gzip,
        #     },
        #     show_cluster=show_cluster, 
        #     auto=auto,
        # )


    def _run(
        self, 
        force, 
        ipyclient, 
        name_fields, 
        name_separator, 
        dry_run,
        split_pairs, 
        gzip):
        "Download files and fastq-dump them to workdir"

        # get run info and sort so largest samples are on top
        df = self.fetch_runinfo(list(range(31)), quiet=True)
        df = df.sort_values(by="spots", ascending=False)

        # parallelize downloads
        if ipyclient:
            lbview = ipyclient.load_balanced_view()

        # make empty Accession field
        df["Accession"] = ""

        # choose spacer to replace spaces in names as different from name_sep
        otherspacer = ("_" if name_separator != "_" else "-")

        # select names for downloaded .sra files
        if name_fields:

            # indices of runinfo fields for names
            fields = [i - 1 for i in fields_checker(name_fields)]

            # set new accession name
            for row in df.index:
                df.loc[row, "Accession"] = (
                    name_separator.join(
                        [df.iloc[row, i] for i in fields]
                        )
                    ).replace(" ", otherspacer)

        # backup default naming scheme
        else:
            if df.SampleName.value_counts().max() > 1:
                # set new accession name
                for row in df.index:
                    df.loc[row, "Accession"] = (
                        name_separator.join(
                            [df.iloc[row, i] for i in [30, 1]]
                            )
                        )
            else:
                df.Accession = df.SampleName       

        # test run to see file names and location without download
        if dry_run:
            print("\rThe following files will be written to: {}\n"
                .format(self.workdir))
            print("{}\n".format(df.Accession))
            return

        # send download jobs
        nfinished = 0
        ntotal = int(df.shape[0]) * 2
        start = time.time()
        message = "downloading/converting fastq data files"
        download_asyncs = {}
        for sidx in df.index:
            progressbar(nfinished, ntotal, start, message)
            acc = df.Accession[sidx]
            url = df.download_path[sidx]
            out = os.path.join(self.workdir, acc) + ".sra"
            out = os.path.realpath(os.path.expanduser(out))

            if ipyclient:
                download_asyncs[acc] = lbview.apply(download_file, *(url, out)) 
            else:
                download_asyncs[acc] = download_file(url, out)
                nfinished += 1
            time.sleep(1.1)

        # continue until all jobs finish
        while 1:
            # track progress and break
            progressbar(nfinished, ntotal, start, message)
            if nfinished == ntotal:
                print("")
                break

            # submit conversion job on finished downloads
            running = list(download_asyncs.keys())
            for key in running:
                if ipyclient:                
                    job = download_asyncs[key]
                    if job.ready():
                        if job.successful():
                            nfinished += 1  

                            # submit new job
                            srr = job.get()
                            paired = bool(split_pairs)
                            args = (srr, paired, gzip)
                            self._call_fastq_dump_on_SRRs(*args)
                            download_asyncs.pop(key)
                            nfinished += 1
                        else:
                            raise IPyradError(job.get())
                else:
                    srr = download_asyncs[key]
                    paired = bool(split_pairs)
                    args = (srr, paired, gzip)
                    self._call_fastq_dump_on_SRRs(*args)
                    download_asyncs.pop(key)
                    nfinished += 1
                
        # final report
        self._report(int(ntotal / 2))


    def _report(self, N):
        print("\n{} fastq files downloaded to {}".format(N, self.workdir))


    @property
    def fetch_fields(self):
        "The column names (fields) in an SRA Run Table."        
        fields = pd.DataFrame(COLNAMES, columns=["field"])
        fields.index += 1
        return fields

    @property
    def fields(self):
        "The column names (fields) in an SRA Run Table"
        fields = pd.DataFrame(COLNAMES, columns=["field"])
        fields.index += 1
        return fields


    def fetch_runinfo(self, fields=None, quiet=False):
        """
        Query the RunInfo for a Sample or Run, returned as a DataFrame. 
        The fields can be subselected. See <self>.fields for options.
        """
        if not quiet: 
            print("\rFetching project data...", end="")

        if fields is None:
            fields = list(range(30))
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
                    "email": "de2356@columbia.edu",
                    "retmax": 1000,
                    },
                )
            sra_ids += [i[4:-5] for i in res.text.split() if "<Id>" in i]
            if not sra_ids:
                raise IPyradError("No SRA samples found in {}"
                    .format(self.accessions))
            time.sleep(3)

        # SRA Runinfo for each ID
        res = requests.get(
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", 
            params={
                "db": "sra",
                "id": ",".join(sra_ids),
                "tool": "ipyrad", 
                "email": "de2356@columbia.edu",
                "rettype": "runinfo", 
                "retmode": "text",                
                },
            )
        df = pd.read_csv(StringIO(res.text.strip()))
        return df.iloc[:, [i - 1 for i in fields]]



    def _call_fastq_dump_on_SRRs(self, srr, paired, gzip):
        """
        calls fastq-dump on SRRs, relabels fastqs by their accession
        names, and writes them to the workdir. Saves temp sra files
        in the designated tmp folder and immediately removes them.
        """

        # build outname
        outname = os.path.split(srr)[-1]
        outname = outname.rsplit(".sra")[0]

        ## build command for fastq-dumping
        fd_cmd = [
            "fastq-dump", srr,
            "--accession", outname,
            "--outdir", self.workdir, 
            # "--disable-multithreading",
            ]
        if gzip:
            fd_cmd += ["--gzip"]
        if paired:
            fd_cmd += ["--split-files"]

        ## call fq dump command
        proc = sps.Popen(fd_cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc.communicate()

        if proc.returncode:
            raise IPyradError(o.decode())

        ## delete the temp sra file from the place 
        if os.path.exists(srr):
            os.remove(srr)



def download_file(url, outname):
    " NOTE the stream=True parameter"
    res = requests.get(url, stream=True)
    with open(outname, 'wb') as f:
        for chunk in res.iter_content(chunk_size=1024): 
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
    return outname



def fields_checker(fields):
    """
    returns a fields argument formatted as a list of strings.
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
To avoid corruption the file was file was removed. Try downloading again to get
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
