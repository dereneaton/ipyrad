#!/usr/bin/env python

"sra tools wrapper to download archived seq data"

# py2/3
from __future__ import print_function

# standard
import os
import sys
import time
import glob
import shutil
import datetime
import multiprocessing
import subprocess as sps

# third party
import pandas as pd
import ipyparallel as ipp
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import Params, progressbar


# TODO: register this as an app with ncbi and avoid sratools altogether


## raise warning if missing imports
MISSING_IMPORTS = """
To use the ipa.sratools module you must install two additional 
libraries which can be done with the following conda command. 

  conda install -c bioconda sra-tools entrez-direct
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
        for binary in ['fastq-dump', 'esearch']:
            proc = sps.Popen(['which', binary], stdout=sps.PIPE)
            comm = proc.communicate()[0]
            if not comm:
                raise IPyradError(MISSING_IMPORTS)


    def run(self, 
        force=False, 
        ipyclient=None, 
        name_fields=30, 
        name_separator="_", 
        dry_run=False):
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
        """

        ## temporarily set directory for tmpfiles used by fastq-dump
        ## if this fails then just skip it.
        try:
            ## ensure output directory, also used as tmpdir
            if not os.path.exists(self.workdir):
                os.makedirs(self.workdir)

            ## get original directory for sra files 
            ## probably /home/ncbi/public/sra by default.
            self._set_vdbconfig_path()

            ## register ipyclient for cleanup
            if ipyclient:
                self._ipcluster["pids"] = {}
                for eid in ipyclient.ids:
                    engine = ipyclient[eid]
                    if not engine.outstanding:
                        pid = engine.apply(os.getpid).get()
                        self._ipcluster["pids"][eid] = pid               

            ## submit jobs to engines or local 
            self._submit_jobs(
                force=force, 
                ipyclient=ipyclient, 
                name_fields=name_fields, 
                name_separator=name_separator,
                dry_run=dry_run,
                )

        ## exceptions to catch, cleanup and handle ipyclient interrupts
        except KeyboardInterrupt:
            print("keyboard interrupt...")

        finally:
            ## reset working sra path
            self._restore_vdbconfig_path()

            ## if it made a new sra directory then it should be empty when 
            ## we are finished if all .sra files were removed. If so, then
            ## let's also remove the dir. if not empty, leave it.
            sradir = os.path.join(self.workdir, "sra")
            if os.path.exists(sradir) and (not os.listdir(sradir)):
                shutil.rmtree(sradir)
            else:
                ## print warning
                try:
                    if os.path.exists(sradir) and os.listdir(sradir):
                        print(FAILED_DOWNLOAD.format(os.listdir(sradir)))
                except OSError as inst:
                    ## If sra dir doesn't even exist something is broken.
                    raise IPyradError("Download failed. Exiting.")

                ## remove fastq file matching to cached sra file
                for srr in os.listdir(sradir):
                    isrr = srr.split(".")[0]
                    ipath = os.path.join(self.workdir, "*_{}*.gz".format(isrr))
                    ifile = glob.glob(ipath)[0]
                    if os.path.exists(ifile):
                        os.remove(ifile)
                ## remove cache of sra files
                if os.path.exists(sradir):
                    shutil.rmtree(sradir)

            ## cleanup ipcluster shutdown
            if ipyclient:
                ## send SIGINT (2) to all engines still running tasks
                try:
                    ipyclient.abort()
                    time.sleep(0.5)
                    for engine_id, pid in self._ipcluster["pids"].items():
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                        time.sleep(0.1)
                except ipp.NoEnginesRegistered:
                    pass
                ## clean memory space
                if not ipyclient.outstanding:
                    ipyclient.purge_everything()
                ## uh oh, kill everything, something bad happened
                else:
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                    print("\nwarning: ipcluster shutdown and must be restarted")
                    


    def _submit_jobs(self, 
        force, 
        ipyclient, 
        name_fields, 
        name_separator, 
        dry_run):
        """
        Download the accessions into a the designated workdir. 
        If file already exists it will only be overwritten if 
        force=True. Temporary files are removed. 
        """

        ## get Run data with default fields (1,4,6,30)
        df = self.fetch_runinfo(list(range(31)), quiet=True)
        sys.stdout.flush()

        ## if not ipyclient then use multiprocessing
        if ipyclient:
            lb = ipyclient.load_balanced_view()

        ## if Run has samples with same name (replicates) then 
        ## we need to include the accessions in the file names
        if name_fields:
            ## indexing requires -1 ints
            fields = [int(i) - 1 for i in fields_checker(name_fields)]
            ## make accession names, no spaces allowed
            df['Accession'] = pd.Series(df[df.columns[fields[0]]], index=df.index)
            for field in fields[1:]:
                df.Accession += name_separator + df[df.columns[field]]
            df.Accession = [i.replace(" ", "_") for i in df.Accession]    
            ## check that names are unique
            if not df.Accession.shape[0] == df.Accession.unique().shape[0]:
                raise IPyradError("names are not unique:\n{}"
                                  .format(df.Accession))

        ## backup default naming scheme
        else:
            if len(set(df.SampleName)) != len(df.SampleName):
                accs = (i+"-"+j for i, j in zip(df.SampleName, df.Run))
                df.Accession = accs
            else:
                df.Accession = df.SampleName

        if dry_run:
            print("\rThe following files will be written to: {}".format(self.workdir))
            print("{}\n".format(df.Accession))
        else:
            ## iterate over and download
            asyncs = []
            for idx in df.index:

                ## get args for this run
                srr = df.Run[idx]
                outname = df.Accession[idx]
                paired = df.spots_with_mates.values.astype(int).nonzero()[0].any()
                fpath = os.path.join(self.workdir, outname + ".fastq.gz")

                ## skip if exists and not force
                skip = False
                if force:
                    if os.path.exists(fpath):
                        os.remove(fpath)
                else:
                    if os.path.exists(fpath):                
                        skip = True
                        sys.stdout.flush()
                        print("[skip] file already exists: {}".format(fpath))

                ## single job progress bar
                tidx = df.Accession.shape[0]
                #if not ipyclient:
                    
                ## submit job to run
                if not skip:
                    args = (self, srr, outname, paired)
                    if ipyclient:
                        rasync = lb.apply_async(call_fastq_dump_on_SRRs, *args)
                        asyncs.append(rasync)
                    else:
                        print("Downloading file {}/{}: {}".format(idx + 1, tidx, fpath))
                        call_fastq_dump_on_SRRs(*args)
                        sys.stdout.flush()

            ## progress bar while blocking parallel
            if ipyclient:
                tots = df.Accession.shape[0]
                printstr = " Downloading fastq files | {} | "
                start = time.time()
                while 1:
                    elapsed = datetime.timedelta(seconds=int(time.time() - start))
                    ready = sum([i.ready() for i in asyncs])
                    progressbar(tots, ready, printstr.format(elapsed), spacer="")
                    time.sleep(0.1)
                    if tots == ready:
                        print("")
                        break
                self._report(tots)

                ## check for fails
                for rasync in asyncs:
                    if not rasync.successful():
                        raise IPyradError(rasync.result())



    def _report(self, N):
        print("{} fastq files downloaded to {}".format(N, self.workdir))


    @property
    def fetch_fields(self):
        fields = pd.DataFrame(
            data=[COLNAMES, list(range(1, len(COLNAMES) + 1))]
        ).T
        fields.columns = ['field', 'index']
        return fields


    def fetch_runinfo(self, fields=None, quiet=False):
        """
        Call esearch to grep SRR info for a project (SRP). Use the command
        sra.fetch_fields to see available fields to be fetched. This function
        returns a DataFrame with runinfo for the selected fields.

        Parameters:
        -----------
        Fields: (tuple or list)
            The default fields returned are 1-30. You can enter a list 
            or tuple of fewer numbers to select fewer fields. Example, 
            (1,4,6,29,30) returns a neat dataframe with Run IDs, 
            Number of reads (SE and PE), ScientificName, and SampleName. 
        """
        if not quiet:
            print("\rFetching project data...", end="")

        ## if no entry then fetch (nearly) all fields.
        if fields is None:  
            fields = range(30)
        fields = fields_checker(fields)

        ## command strings
        es_cmd = [
            "esearch", 
            "-db", "sra", 
            "-query", self.accession,
        ]

        ef_cmd = [
            "efetch", 
            "--format", "runinfo",
        ]

        cut_cmd = [
            "cut", 
            "-d", ",", 
            "-f", ",".join(fields),
        ]

        ## pipe commands together
        proc1 = sps.Popen(es_cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc2 = sps.Popen(ef_cmd, stdin=proc1.stdout, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc3 = sps.Popen(cut_cmd, stdin=proc2.stdout, stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc3.communicate()
        proc2.stdout.close()
        proc1.stdout.close()
        
        if o:
            vals = o.decode().strip().split("\n")
            names = vals[0].split(",")
            items = [i.split(",") for i in vals[1:] if i not in ["", vals[0]]]
            return pd.DataFrame(items, columns=names)
        else:
            raise IPyradError("no samples found in {}".format(self.accession))


    def _set_vdbconfig_path(self):

        ## get original path
        proc = sps.Popen(
            ['vdb-config', '-p'], 
            stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc.communicate()
        self._oldtmpdir = o.split("root>")[1][:-2]

        ## set new temp dir 
        proc = sps.Popen(
            ['vdb-config', '-s', 
            'repository/user/main/public/root='+self.workdir], 
            stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc.communicate()
        #print('setting tmpdir to {}'.format(self.workdir))


    def _restore_vdbconfig_path(self):
        ## set temp dir 
        if not self._oldtmpdir:
            self._oldtmpdir = os.path.join(os.path.expanduser("~"), "ncbi")
        proc = sps.Popen(
            ['vdb-config', '-s', 
            'repository/user/main/public/root='+self._oldtmpdir],
            stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc.communicate()
        #print('restoring tmpdir to {}'.format(self._oldtmpdir))



def call_fastq_dump_on_SRRs(self, srr, outname, paired):
    """
    calls fastq-dump on SRRs, relabels fastqs by their accession
    names, and writes them to the workdir. Saves temp sra files
    in the designated tmp folder and immediately removes them.
    """

    ## build command for fastq-dumping
    fd_cmd = [
        "fastq-dump", srr,
        "--accession", outname,
        "--outdir", self.workdir, 
        "--gzip",
        ]
    if paired:
        fd_cmd += ["--split-files"]

    ## call fq dump command
    proc = sps.Popen(fd_cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
    o, e = proc.communicate()

    ## delete the stupid temp sra file from the place 
    ## that it is very hard-coded to be written to, and 
    ## LEFT IN, for some crazy reason.
    srafile = os.path.join(self.workdir, "sra", srr+".sra")
    if os.path.exists(srafile):
        os.remove(srafile)



def fields_checker(fields):
    """
    returns a fields argument formatted as a list of strings.
    and doesn't allow zero.
    """
    ## make sure fields will work
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
    fields = [i for i in fields if i != '0']

    return fields


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
