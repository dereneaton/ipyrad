#!/usr/bin/env ipython2

""" D-statistic calculations """
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142
# pylint: disable=R0915
# pylint: disable=R0914
# pylint: disable=R0912


from __future__ import print_function
import os
import sys
import subprocess as sps
from ipyrad.assemble.util import IPyradWarningExit, progressbar


## raise warning if missing imports
MISSING_IMPORTS = """
You are missing required packages to use ipa.sratools. 
First run the following two conda install commands:

  conda install -c bioconda sra-tools
  conda install -c bioconda entrez-direct
"""

ACCESSION_ID = """
Accession ID must be either a Run or Study accession, i.e., 
it must have on the following prefixes:
  Study: SRR, ERR, DRR
  Project: SRP, ERP, DRP
"""


class SRA(object):
    """ ipyrad.analysis SRA download object"""
    def __init__(self, 
        accession,
        workdir,
        ):

        ## check imports
        for binary in ['fastq-dump', 'esearch']:
            if not sps.call("type " + binary, 
                        shell=True,
                        stdout=sps.PIPE,
                        stderr=sps.PIPE) == 0:
                raise IPyradWarningExit(MISSING_IMPORTS)

        ## store attributes
        self.accession = accession
        self.workdir = os.path.abspath(os.path.expanduser(workdir))
        self.is_sample = False
        self.is_project = False

        ## 
        if any([i in self.accession for i in ["SRR", "ERR", "DRR"]]):        
            self.is_sample = True
        elif any([i in self.accession for i in ["SRP", "ERP", "DRP"]]):
            self.is_project = True
        else:
            raise IPyradWarningExit(ACCESSION_ID)



    def run(self, force=False):
        """
        Download the accessions into a the designated workdir. 
        If file already exists it will only be overwritten if 
        force=True. Temporary files are removed. 
        """

        ## ensure output directory
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        ## download files
        if self.is_project:
            print("\rFetching project data...", end="")
            srrs, accs = self.fetch_runinfo()
            sys.stdout.flush()
            fn = 0
            for srr, acc in zip(srrs, accs):
                print("\rDownloading file {} of {}: {}.fastq.gz"\
                    .format(fn+1, len(srr), acc), end="")
                self._accession = srr

                skip = False
                fpath = os.path.join(self.workdir, acc+".fastq.gz")
                if force:
                    os.remove(fpath)
                else:
                    if os.path.exists(fpath):
                        skip = True
                        sys.stdout.flush()
                        print(" - skip - already exists in workdir")
                if not skip:
                    self._call_fastq_dump_on_SRRs(rename=acc)
                sys.stdout.flush()
                fn += 1
            self._report(fn)

        else:
            self._accession = self.accession
            self._call_fastq_dump_on_SRRs()
            self.report(1)


    def _report(self, N):
        print("{} fastq files downloaded to {}".format(N, self.workdir))


    def fetch_runinfo(self):
        """
        Call esearch to grep SRR info for a project (SRP). 
        Returns two lists: SRRs and ACCs. 
         """
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
            "-f", "1,30",
        ]

        grep_cmd = [
            "grep", "SRR"
        ]

        ## pipe commands together
        proc1 = sps.Popen(es_cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc2 = sps.Popen(ef_cmd, stdin=proc1.stdout, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc3 = sps.Popen(cut_cmd, stdin=proc2.stdout, stderr=sps.STDOUT, stdout=sps.PIPE)
        proc4 = sps.Popen(grep_cmd, stdin=proc3.stdout, stderr=sps.STDOUT, stdout=sps.PIPE)
        o, e = proc4.communicate()
        
        if o:
            srrlist = o.strip().split("\n")
            SRRs, ACCs = zip(*[i.split(",") for i in srrlist])
            return SRRs, ACCs 
        else:
            raise IPyradWarningExit("no samples found in {}".format(accession))



    def _call_fastq_dump_on_SRRs(self, rename=None):
        """
        calls fastq-dump on SRRs, relabels fastqs by their accession
        names, and writes them to the workdir.
        """

        ## build command
        fd_cmd = [
            "fastq-dump", self._accession, 
            "--accession", rename, 
            "--outdir", self.workdir, 
            "--gzip",
            ]

        ## call fq dump command
        proc = sps.Popen(fd_cmd)
        proc.communicate()

        ## delete the stupid temp sra file from the place 
        ## that it is very hard-coded to be written to.
        srafile = os.path.join(
            os.path.expanduser("~"), 
            "ncbi", 
            "public", 
            "sra", 
            self._accession+".sra")
        if os.path.exists(srafile):
            os.remove(srafile)


