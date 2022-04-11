#!/usr/bin/env python

"""Params and Hackers schemas for type checking and serialization.

Pydantic Models are similar to dataclases but they also include
*type validation*, meaning that if you try to set an attribute to
the wrong type it will raise an error.
"""

# pylint: disable=no-self-argument, no-name-in-module, no-self-use

import glob
from pathlib import Path
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
    phred_qscore_offset: int = 33    

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
    project_dir: Path = Path.cwd()
    raw_fastq_path: Path = None
    barcodes_path: Path = None
    sorted_fastq_path: Path = None
    assembly_method: AssemblyMethod = "denovo"
    reference_sequence: Path = None
    datatype: DataType = "rad"
    restriction_overhang: Tuple[str, str] = ("TGCAG", "")
    max_low_qual_bases: int = 5
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
    pop_assign_file: Path = None
    reference_as_filter: Path = None

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
    def _dir_validator(cls, value):
        """Project_dir cannot have whitespace, is expanded, and created."""
        if ' ' in value.name:
            raise ValueError('project_dir cannot contain spaces')
        value = value.expanduser().resolve()
        if value.exists():
            value.mkdir(exist_ok=True)
        return value

    @validator('raw_fastq_path', 'sorted_fastq_path', 'barcodes_path')
    def _path_validator(cls, value):
        """If a path was entered then it must match >=1 files."""
        if not value:
            return value
        if "Merged:" in value.name:
            return value
        value = value.expanduser().resolve()
        if value.is_dir():
            raise ValueError(
                "You entered a dir path where you must enter a file path.\n"
                "To select multiple files in a dir you can use regular \n"
                "expressions, e.g., './data/*.fastq.gz'.")
        if not glob.glob(str(value)):
            raise ValueError(f"no files match the input string: {value}")
        return value

    @validator("reference_sequence", "reference_as_filter")
    def _reference_validator(cls, value):
        """Checks that reference file exists and expands path."""
        if value:
            value = value.expanduser().resolve()
            if not value.exists():
                raise ValueError(f"no files match the input string: {value}")                
            if value.suffix == ".gz":
                raise ValueError(f"reference {value} must be decompressed.")
        return value

    @validator("restriction_overhang")
    def _enzyme_validator(cls, value):
        """Check that the restriction enzymes do not contain bad chars?"""
        # TODO:
        return value

    @validator("datatype")
    def _datatype_validator(cls, value):
        """Return the Enum value (string)"""
        return value.value

    @validator("assembly_method")
    def _method_validator(cls, value):
        """Return the Enum value (string)."""
        return value.value        


    ### TODO...


if __name__ == "__main__":

    p = ParamsSchema(assembly_name="TEST", project_dir="./")
    p.reference_as_filter = "../../tests/ipsimdata/gbs_example_genome.fa"
    print(p)
