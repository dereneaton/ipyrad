#!/usr/bin/env python

"""
A generic tool for streaming download of online files 
using the requests library, with option to gunzip.
"""

import os
import gzip
import subprocess
import requests



class Download:
    """
    Download a large file by streaming in chunks using the requests library.
    The file from 'url' is downloaded to the location 'path'. The file will
    not be downloaded if it already exists at this path, unless you overwrite
    using the force option.

    Parameters
    ===========
    url (str):
        A valid URL destination of the file you wish to download.
    path (str):
        A valid path on your computer to use as the file name.
    gunzip (bool):
        If the file ends with .gz and is gzipped then this will 
        write a copy that is decompressed without the .gz ending.
    force (bool):
        Overwrite existing file with the same name as 'path'.

    Returns
    ==========
    None
    """
    def __init__(self, url, path, gunzip=False, force=False):
        self.url = url
        self.path = path
        self.force = force
        self.gunzip = gunzip
        self.gunzip_name = self.path.split(".gz")[0]

        # download the data
        self.download()
        self.gunzip_file()
        # self.clean_data()


    def download(self):
        "call the chunked download with requests"
        # only run if the reference doesn't already exist
        if (not os.path.exists(self.path)) and (not self.force):
    
            # open a stream to url and write to file 1Mb at a time.
            res = requests.get(self.url, stream=True)
            with open(self.path, 'wb') as out:
                for chunk in res.iter_content(chunk_size=1024*1024):
                    if chunk:
                        out.write(chunk)
            print("successful download: {}".format(self.path))
        else:
            print("file already exists: {}".format(self.path))


    def gunzip_file(self):
        "make a decompressed copy of the file"
        if self.gunzip:
            try:
                if not os.path.exists(self.gunzip_name):
                    print('decompressing gzipped file')
                    with open(self.gunzip_name, 'w') as out:
                        with gzip.open(self.path, 'r') as stream:
                            while 1:
                                chunk = stream.read(1028 * 10)
                                if not chunk:
                                    break
                                enc = chunk.decode()
                                out.write(enc)
                else:
                    print("decompressed file already exists: {}"
                        .format(self.gunzip_name))
            except Exception:
                print("error: couldn't gunzip file.")


    def clean_data(self):
        """
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