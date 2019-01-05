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
import tempfile
import requests
import subprocess as sps

# third party
import pandas as pd
import ipyparallel as ipp
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import progressbar


## raise warning if missing imports
MISSING_IMPORTS = """
To use the ipa.sratools module you must install two additional 
libraries which can be done with the following conda command. 

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
        accession,
        workdir="sra-fastq-data",
        ):

        # store attributes
        self.accession = accession
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        self.is_sample = False
        self.is_project = False
        self._oldtmpdir = None

        ## cluster attributes
        self._ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }

        ## 
        if any([i in self.accession for i in ["SRR", "ERR", "DRR"]]):        
            self.is_sample = True
        elif any([i in self.accession for i in ["SRP", "ERP", "DRP"]]):
            self.is_project = True
        else:
            raise IPyradError(ACCESSION_ID)

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
        ipyclient,
        force=False, 
        name_fields=(1, 30), 
        name_separator="_", 
        dry_run=False, 
        split_pairs=None,
        gzip=False,
        ):
        """
        Download the accessions into a the designated workdir. 

        Parameters
        ----------
        force: (bool)
            If force=True then existing files with the same name
            will be overwritten. 

        ipyclient: (ipyparallel.Client)
            If provided, work will be distributed across a parallel
            client, otherwise download will be run on a single core.

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
        """

        # ensure output directory, also used as tmpdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        try:
            ## register ipyclient for cleanup
            if ipyclient:
                self._ipcluster["pids"] = {}
                for eid in ipyclient.ids:
                    engine = ipyclient[eid]
                    if not engine.outstanding:
                        pid = engine.apply(os.getpid).get()
                        self._ipcluster["pids"][eid] = pid               

            ## submit jobs to engines or local 
            # self._submit_jobs(
            self._distribute_jobs(
                force=force, 
                ipyclient=ipyclient, 
                name_fields=name_fields, 
                name_separator=name_separator,
                dry_run=dry_run,
                split_pairs=split_pairs,
                gzip=gzip,
                )

        ## exceptions to catch, cleanup and handle ipyclient interrupts
        except KeyboardInterrupt:
            print("\nkeyboard interrupt...")


        finally:
            ## cleanup ipcluster shutdown
            if ipyclient:
                ## send SIGINT (2) to all engines still running tasks
                try:
                    ipyclient.abort()
                    time.sleep(0.5)
                    for engine_id, pid in self._ipcluster["pids"].items():
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                            print('killed', pid)
                        time.sleep(0.1)
                except ipp.NoEnginesRegistered:
                    pass

                # clean memory space
                if not ipyclient.outstanding:
                    ipyclient.purge_everything()

                # uh oh, kill everything, something bad happened
                else:
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                    print("\nwarning: ipcluster shutdown and must be restarted")


    def _distribute_jobs(
        self, 
        force, 
        ipyclient, 
        name_fields, 
        name_separator, 
        dry_run,
        split_pairs, 
        gzip):
        "Download files and fasterq-dump them to workdir"

        # get run info and sort so largest samples are on top
        df = self.fetch_runinfo(list(range(31)), quiet=True)
        df = df.sort_values(by="spots", ascending=False)

        # parallelize downloads
        if ipyclient:
            lbview = ipyclient.load_balanced_view()

        # make empty Accession field
        df["Accession"] = ""

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
                    )

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
        start = time.time()
        message = "downloading/converting fastq data files"
        progressbar(0, 1, start, message)
        download_asyncs = {}
        for sidx in df.index:
            acc = df.Accession[sidx]
            url = df.download_path[sidx]
            out = os.path.join(self.workdir, acc) + ".sra"
            out = os.path.realpath(os.path.expanduser(out))
            download_asyncs[acc] = lbview.apply(download_file, *(url, out))
            time.sleep(1.1)

        # collect results and send to fasterq-dump one at a time
        ntotal = len(download_asyncs) * 2
        nfinished = 0

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

        # final report
        self._report(int(ntotal / 2))



    # def _submit_jobs(self, 
    #     force, 
    #     ipyclient, 
    #     name_fields, 
    #     name_separator, 
    #     dry_run, 
    #     split_pairs,
    #     ):
    #     """
    #     Download the accessions into a the designated workdir. 
    #     If file already exists it will only be overwritten if 
    #     force=True. Temporary files are removed. 
    #     """

    #     ## get Run data with default fields (1,4,6,30)
    #     df = self.fetch_runinfo(list(range(31)), quiet=True)
    #     sys.stdout.flush()

    #     ## if not ipyclient then use multiprocessing?
    #     if ipyclient:
    #         lb = ipyclient.load_balanced_view()

    #     ## if Run has samples with same name (replicates) then 
    #     ## we need to include the accessions in the file names
    #     if name_fields:
    #         ## indexing requires -1 ints
    #         fields = [int(i) - 1 for i in fields_checker(name_fields)]
    #         ## make accession names, no spaces allowed
    #         df['Accession'] = pd.Series(df[df.columns[fields[0]]], index=df.index)
    #         for field in fields[1:]:
    #             df.Accession += name_separator + df[df.columns[field]]
    #         df.Accession = [i.replace(" ", "_") for i in df.Accession]    
    #         ## check that names are unique
    #         if not df.Accession.shape[0] == df.Accession.unique().shape[0]:
    #             raise IPyradError("names are not unique:\n{}"
    #                               .format(df.Accession))

    #     ## backup default naming scheme
    #     else:
    #         if len(set(df.SampleName)) != len(df.SampleName):
    #             accs = ("{}-{}".format(i, j) for (i, j) in zip(df.SampleName, df.Run))
    #             df.Accession = accs
    #         else:
    #             df.Accession = df.SampleName

    #     if dry_run:
    #         print("\rThe following files will be written to: {}".format(self.workdir))
    #         print("{}\n".format(df.Accession))
    #     else:
    #         ## iterate over and download
    #         asyncs = []
    #         for idx in df.index:

    #             ## get args for this run
    #             srr = df.Run[idx]
    #             outname = df.Accession[idx]
    #             fpath = os.path.join(self.workdir, outname + ".fastq.gz")

    #             ## get paired arg
    #             if not split_pairs:
    #                 paired = int(df.spots_with_mates[idx])
    #             else:
    #                 paired = bool(split_pairs)

    #             ## skip if exists and not force
    #             skip = False
    #             if force:
    #                 if os.path.exists(fpath):
    #                     os.remove(fpath)
    #             else:
    #                 if os.path.exists(fpath):                
    #                     skip = True
    #                     sys.stdout.flush()
    #                     print("[skip] file already exists: {}".format(fpath))

    #             ## single job progress bar
    #             tidx = df.Accession.shape[0]
    #             #if not ipyclient:
                    
    #             ## submit job to run
    #             if not skip:
    #                 args = (self, srr, outname, paired)
    #                 if ipyclient:
    #                     rasync = lb.apply_async(call_fastq_dump_on_SRRs, *args)
    #                     asyncs.append(rasync)
    #                 else:
    #                     print("Downloading file {}/{}: {}".format(idx + 1, tidx, fpath))
    #                     call_fastq_dump_on_SRRs(*args)
    #                     sys.stdout.flush()

    #         ## progress bar while blocking parallel
    #         if ipyclient:
    #             tots = df.Accession.shape[0]
    #             printstr = ("Downloading fastq files", "")
    #             start = time.time()
    #             while 1:
    #                 ready = sum([i.ready() for i in asyncs])
    #                 progressbar(tots, ready, start, printstr)
    #                 time.sleep(0.1)
    #                 if tots == ready:
    #                     print("")
    #                     break
    #             self._report(tots)

    #             ## check for fails
    #             for rasync in asyncs:
    #                 if not rasync.successful():
    #                     raise IPyradError(rasync.result())



    def _report(self, N):
        print("\n{} fastq files downloaded to {}".format(N, self.workdir))


    @property
    def fetch_fields(self):
        fields = pd.DataFrame(COLNAMES, columns=["field"])
        fields.index += 1
        return fields


    def fetch_runinfo(self, fields=None, quiet=False):
        "Requests based utils for fetching SRA IDs and URLs"
        # 
        if not quiet: 
            print("\rFetching project data...", end="")

        if fields is None:
            fields = list(range(30))
        fields = fields_checker(fields)

        # SRA IDs from the SRP 
        res = requests.get(
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", 
            params={
                "db": "sra",
                "term": self.accession,
                "tool": "ipyrad", 
                "email": "de2356@columbia.edu",
                },
            )
        sra_ids = [i[4:-5] for i in res.text.split() if "<Id>" in i]
        if not sra_ids:
            raise IPyradError("No SRA samples found in {}"
                .format(self.accession))

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
        time.sleep(3)
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



    # def _set_vdbconfig_path(self):

    #     ## get original path
    #     proc = sps.Popen(
    #         ['vdb-config', '-p'], 
    #         stderr=sps.STDOUT, stdout=sps.PIPE)
    #     o, e = proc.communicate()
    #     self._oldtmpdir = o.decode().split("root>")[1][:-2]

    #     ## set new temp dir 
    #     proc = sps.Popen(
    #         ['vdb-config', '-s', 
    #         'repository/user/main/public/root='+self.workdir], 
    #         stderr=sps.STDOUT, stdout=sps.PIPE)
    #     o, e = proc.communicate()
    #     #print('setting tmpdir to {}'.format(self.workdir))


    # def _restore_vdbconfig_path(self):
    #     ## set temp dir 
    #     if not self._oldtmpdir:
    #         self._oldtmpdir = os.path.join(os.path.expanduser("~"), "ncbi")
    #     proc = sps.Popen([
    #         'vdb-config', '-s', 
    #         'repository/user/main/public/root=' + self._oldtmpdir],
    #         stderr=sps.STDOUT, stdout=sps.PIPE)
    #     o, e = proc.communicate()
    #     #print('restoring tmpdir to {}'.format(self._oldtmpdir))



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



## OMG FASTERQ-DUMP IS THE WORST DON"T TRY THIS...
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
