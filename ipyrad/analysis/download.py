#!/usr/bin/env python

"""
A generic tool for streaming download of online files 
using the requests library, with option to gunzip.
"""

import os
import gzip
import subprocess
import requests
from loguru import logger


class Download:
    """
    Download a large file by streaming in chunks using the requests 
    library. The file from 'url' is downloaded to the location 'path'. 
    The file will not be downloaded if it already exists at this path, 
    unless you overwrite using the force option. You can optionally
    decompress the downloaded file with gunzip after downloading.

    Parameters
    ===========
    url (str):
        A valid URL destination of the file you wish to download.
    path (str):
        A valid path on your computer to use as the file name. If
        the directory does not exist we will try to create it. If
        the file is gzipped then the path should end with .gz.
    gunzip (bool):
        If the file ends with .gz and is gzipped then this will 
        create a copy that is decompressed without the .gz ending.
    force (bool):
        Overwrite existing file with the same name as 'path'.
    """
    def __init__(
        self, 
        url: str, 
        path: str, 
        gunzip: bool=False, 
        ):
        self.url = url
        self.path = path
        self.gunzip = gunzip
        self.gunzip_name = self.path.split(".gz")[0]


    def run(self, force):
        """
        Runs the tool to download the file as a stream to the 
        specified path and optionally decompress.
        """
        self.download(force)
        self.gunzip_file()
        # self.clean_data()


    def download(self, force):
        """
        Call the chunked download with requests
        """
        # attempt to create outdir if it doesn't exist
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

        # only run if the reference doesn't already exist
        if (not os.path.exists(self.path)) or (not force):
    
            # open a stream to url and write to file 1Mb at a time.
            res = requests.get(self.url, stream=True)
            with open(self.path, 'wb') as out:
                for chunk in res.iter_content(chunk_size=1024*1024):
                    if chunk:
                        out.write(chunk)
            logger.info(f"successful download: {self.path}")
        else:
            logger.warning(f"file already exists: {self.path}")


    def gunzip_file(self):
        """
        Make a decompressed copy of the file
        """
        if self.gunzip:
            try:
                if not os.path.exists(self.gunzip_name):
                    logger.info(
                        f'decompressing gzipped file to: {self.gunzip_name}')
                    with open(self.gunzip_name, 'w') as out:
                        with gzip.open(self.path, 'r') as stream:
                            while 1:
                                chunk = stream.read(1028 * 10)
                                if not chunk:
                                    break
                                enc = chunk.decode()
                                out.write(enc)
                else:
                    logger.warning(
                        f"decompressed file already exists: {self.gunzip_name}")
            except Exception:
                logger.error("error: couldn't gunzip file.")


    def clean_data(self):
        """
        NOT IMPLEMENTED

        Removes multi-allele and INDEL containing SNPs and all meta-data.
        The input must be bgzip, which requires the input to be non-zipped.

        If you get the error "[E::hts_idx_push] Chromosome blocks not continuous
        tbx_index_build failed: /tmp/wildlabdanio.vcf.gz" this means that 
        your VCF is malformed for use by tabix, and is not an ipa error.

        # this is the pipeline that is called:
          bgzip data.vcf
          tabix data.vcf.gz
          bcftools view -m2 -M2 -i'CIGAR="1X" & QUAL>30' data.vcf.gz -Ou | 
              bcftools annotate -x FORMAT,INFO > data.cleaned.vcf
          bgzip data.cleaned.vcf
        """

        # Check for required software...



        # compress the VCF as bgzip
        proc = subprocess.Popen(['bgzip', self.gunzip_name])
        proc.communicate()

        # tabix index the bgzip file
        proc = subprocess.Popen(["tabix", self.gunzip_name + ".gz"])
        proc.communicate()

        # 