#!/usr/bin/env python

"""Params and Hackers schemas for type checking and serialization.

Pydantic Models are similar to dataclases but they also include
*type validation*, meaning that if you try to set an attribute to
the wrong type it will raise an error.
"""

# pylint: disable=no-self-argument, no-name-in-module, no-self-use

import glob
import os
from enum import Enum
from typing import List, Tuple
from pydantic import BaseModel, Field, validator


class HackersSchema(BaseModel):
    """Hackers dict class for 'extra' options in ipyrad assembly."""
    random_seed: int = 42
    max_fragment_length: int = Field(50, help="limits maxlen of loci by tiling")
    max_inner_mate_distance: int = Field(None, help="used to delimit loci in ref assemblies")
    p5_adapter: str = "AGATCGGAAGAGC"
    p5_adapters_extra: List[str] = None
    p3_adapter: str = "AGATCGGAAGAGC"
    p3_adapters_extra: List[str] = None
    query_cov: float = None
    bwa_args: str = None
    demultiplex_on_i7_tags: bool = False
    declone_PCR_duplicates: bool = False
    merge_technical_replicates: bool = True
    exclude_reference: bool = True
    trim_loci_min_sites: int = 4

class AssemblyMethod(str, Enum):
    """supported assembly method categories"""
    DENOVO = "denovo"
    REFERENCE = "reference"

class DataType(str, Enum):
    """supported datatypes categories"""
    RAD = "rad"
    DDRAD = "ddrad"
    PAIRDDRAD = "pairddrad"
    PAIR3RAD = "pair3rad"
    GBS = "gbs"
    PAIRGBS = "pairgbs"
    _2BRAD = "2brad"

class ParamsSchema(BaseModel):
    """Assembly object parameters for ipyrad assembly."""
    assembly_name: str = Field(allow_mutation=False)
    project_dir: str = os.path.realpath("./")
    raw_fastq_path: str = None
    barcodes_path: str = None
    sorted_fastq_path: str = None
    assembly_method: AssemblyMethod = "denovo"
    reference_sequence: str = None
    datatype: DataType = "rad"
    restriction_overhang: Tuple[str, str] = ("TGCAG", "")
    max_low_qual_bases: int = 5
    phred_qscore_offset: int = 33
    min_depth_statistical: int = 6
    min_depth_majrule: int = 6
    max_depth: int = 10000
    clust_threshold: float = 0.85
    max_barcode_mismatch: int = 0
    filter_adapters: int = 2
    filter_min_trim_len: int = 35
    max_alleles_consens: int = 2
    max_n_consens: float = 0.05
    max_h_consens: float = 0.05
    min_samples_locus: int = 4
    max_snps_locus: float = 0.2
    max_indels_locus: int = 8
    max_shared_h_locus: float = 0.5
    trim_reads: Tuple[int, int, int, int] = (0, 0, 0, 0)
    trim_loci: Tuple[int, int, int, int] = (0, 0, 0, 0)
    output_formats: List[str] = ("p", "s", "l")
    pop_assign_file: str = None
    reference_as_filter: str = None

    class Config:
        """This is required to use allow_mutation=False on name.
        it enables type checking validation when using settattr in API.
        """
        validate_assignment = True

    def __str__(self):
        return self.json(indent=2)

    def __repr__(self):
        return self.json(indent=2)        

    ##################################################################
    # Below here, custom validator funcs in addition to type checking.
    ##################################################################

    @validator('assembly_name')
    def _name_validator(cls, value):
        """Names cannot have whitespace. Other strange characters are 
        simply replaced with a warning message printed.
        """
        if ' ' in value:
            raise ValueError('assembly_name cannot contain spaces')
        return value

    @validator('project_dir')
    def _proj_validator(cls, value):
        """Project_dir cannot have whitespace, is expanded, and created."""
        if ' ' in value:
            raise ValueError('project_dir cannot contain spaces')
        value = os.path.realpath(os.path.expanduser(value))
        if not os.path.exists(value):
            os.makedirs(value, exist_ok=True)
        return value

    @validator('raw_fastq_path')
    def _raw_validator(cls, value):
        """Expand path to check that some files match the regex string.
        If post-merge then no check.
        """
        if value and ("Merged:" not in value):
            value = os.path.realpath(os.path.expanduser(value))
            if os.path.isdir(value):
                raise ValueError(
                    "raw_fastq_path should not be a directory. To select "
                    "multiple files in a directory use a wildcard selector "
                    "e.g., './data/*.fastq.gz'.")
            if not glob.glob(value):
                raise ValueError(
                    "no files match the input string for raw_fastq_path. "
                    f"You entered: {value}")
        return value

    @validator('barcodes_path')
    def _barcode_validator(cls, value):
        """Check that barcodes file exists if a path was entered.
        If post-merge then no check.
        """
        if value and ("Merged:" not in value):
            # expand name
            fullpath = glob.glob(os.path.realpath(os.path.expanduser(value)))
            if not fullpath:
                raise ValueError(
                    f"barcodes_path file not found. You entered {value}")
            value = fullpath[0]
        return value

    @validator('sorted_fastq_path')
    def _sorted_validator(cls, value):
        """Expand path to check that some files match the regex string.
        If post-merge then no check.
        """
        if not value:
            return value
        if value and ("Merged:" not in value):
            value = os.path.realpath(os.path.expanduser(value))
            if os.path.isdir(value):
                raise ValueError(
                    "sorted_fastq_path should not be a directory. To select "
                    "multiple files in a directory use a wildcard selector "
                    "e.g., './data/*.fastq.gz'.")
            if not glob.glob(value):
                raise ValueError(
                    "no files match the input string for sorted_fastq_path. "
                    f"You entered: {value}")
        return value

    @validator("reference_sequence")
    def _reference_validator(cls, value):
        """Checks that reference file exists and expands path."""
        if value:
            fullpath = os.path.realpath(os.path.expanduser(value))
            if not glob.glob(fullpath):
                raise ValueError(
                    f"reference_sequence file not found. You entered {value}")
            value = glob.glob(fullpath)[0]
            if value.endswith(".gz"):
                raise ValueError("reference_sequence must be decompressed.")
        return value

    @validator("restriction_overhang")
    def _enzyme_validator(cls, value):
        """Check that the restriction enzymes do not contain bad chars?"""
        # TODO:
        return value

    @validator("datatype")
    def _datatype_validator(cls, value):
        """
        Return the Enum value (string)
        """
        return value.value

    @validator("assembly_method")
    def _method_validator(cls, value):
        """
        Return the Enum value (string)
        """
        return value.value        


    ### TODO...


    @validator("reference_as_filter")
    def _reference_filter_validator(cls, value):
        """Checks that reference file exists and expands path."""
        if value:
            fullpath = os.path.realpath(os.path.expanduser(value))
            if not glob.glob(fullpath):
                raise ValueError(
                    f"reference_as_filter file not found. You entered {value}")
            value = glob.glob(fullpath)[0]
            if value.endswith(".gz"):
                raise ValueError("reference_as_filter must be decompressed.")
        return value



if __name__ == "__main__":

    p = ParamsSchema(assembly_name="TEST", project_dir="./")
    print(p)
