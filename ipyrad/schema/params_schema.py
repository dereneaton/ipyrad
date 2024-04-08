#!/usr/bin/env python

"""Params and Hackers schemas for type checking and serialization.

Pydantic Models are similar to dataclases but they also include
*type validation*, meaning that if you try to set an attribute to
the wrong type it will raise an error. By using type validation the
pydantic models can be easily serialized to JSON and then reloaded
as the appropriate data types.

Some old parameters were deprecated because we now simply estimate
their values from the data. Examples:

hackers.px_adapters = Always use the Illumina adapters prefix
raw_fastq_path, barcodes_path = Options in demux
demultiplex_on_i7_tags = Option in demux
merge_technical_replicates = Option in demux
datatype = infer if data are paired, infer if revcomp during clustering.
assembly_method = infer based on presence of reference sequence.
restriction_overhang = infer from reads and use, or validate against user input
"""

# pylint: disable=no-self-argument, no-name-in-module, no-self-use

from pathlib import Path
from typing import List, Tuple, Dict
from pydantic import BaseModel, Field, field_validator
from loguru import logger

logger = logger.bind(name="ipyrad")


class Params(BaseModel):
    """Assembly object parameters for ipyrad assembly."""
    assembly_name: str = Field(frozen=True)  # allow_mutation=False)
    project_dir: Path = Field(default_factory=Path.cwd)
    fastq_paths: List[Path] | Path = Field(default_factory=list)
    reference_sequence: Path | None = None
    # filtering/trimming options
    filter_adapters: int = 2
    filter_min_trim_len: int = 35
    max_low_qual_bases: int = 5
    trim_reads: List[int] = [0, 0]
    phred_qscore_offset: int = 33
    # cluster building options
    reference_as_filter: Path | None = None
    clust_threshold: float = 0.85
    min_depth_statistical: int = 6
    min_depth_majrule: int = 6
    max_depth: int = 10_000
    # consens options
    max_alleles_consens: int = 2
    max_n_consens: float = 0.05
    max_h_consens: float = 0.05
    # locus building options
    min_samples_locus: int = 4
    max_snps_locus: float = 0.2
    max_indels_locus: int = 8
    max_shared_h_locus: float = 0.5
    trim_loci: List[int] = [0, 0, 0, 0]
    output_formats: List[str] = ["p", "s", "l"]
    pop_assign_file: Path | None = None
    # potpourri options
    exclude_reference: bool = True
    trim_locus_edges_min_sites: int = 4
    random_seed: int = 42
    max_fragment_length: int = Field(400, help="limits maxlen of loci by tiling. Updated automatically.")
    max_inner_mate_distance: int = Field(500, help="used to delimit loci in ref assemblies")
    declone_PCR_duplicates: bool = False
    technical_replicates: Dict[str, List[str]] = Field(default_factory=dict)# how should users enter this to params file?
    # NOT NECESSARY, WE CAN ESTIMATE FROM THE DATA
    restriction_overhang: List[str] = ["", ""]
    # datatype: DataType = "rad"
    # assembly_method: AssemblyMethod = "denovo"
    # cluster_query_cov: float

    class Config:
        """This is required to use allow_mutation=False on name.
        it enables type checking validation when using settattr in API.
        """
        validate_assignment = True

    def __str__(self):
        return self.model_dump_json(indent=2)

    def __repr__(self):
        return self.model_dump_json(indent=2)

    ##################################################################
    # Below here, custom validator funcs in addition to type checking.
    # These can be modified to support older params files when options
    # change, or to add extra validations.
    #
    # API: This is run during ip.load_json
    # CLI: This is run during ip.load_json w/ existing JSON file and
    # then run again to update the Assembly.ParamsSchema by setting
    # new values from the params file.
    ##################################################################

    @field_validator('assembly_name')
    @classmethod
    def _name_validator(cls, value: str) -> str:
        """Names cannot have whitespace. Other strange characters are
        replaced with a warning message printed. This becomes immutable.
        """
        if ' ' in value:
            raise ValueError('assembly_name cannot contain spaces')
        return value

    @field_validator('project_dir')
    @classmethod
    def _dir_validator(cls, value: Path) -> Path:
        """Project_dir cannot have whitespace, is expanded, and created."""
        if ' ' in value.name:
            raise ValueError('project_dir cannot contain spaces')
        value = value.expanduser().resolve()
        value.mkdir(exist_ok=True)
        return value

    @field_validator('fastq_paths')
    @classmethod
    def _path_validator(cls, value: str | Path | List[str | Path]) -> List[Path]:
        """If a path was entered then it must match >=1 files.

        Users are allowed to enter None to leave this field blank, in
        which case it only raises an exception if you run steps 1 or 2
        when these files are needed. However, if a user enters an
        invalid path where no files exist this will always raise an
        error during JSON or params file parsing.
        """
        # allow it to be empty so user can set it later in the API.
        if not value:
            return value
        # ensure input is a List[Path] or List[str]
        if not isinstance(value, list):
            value = [value]
        # iterate over list checking that each can be expanded to match
        # one more files. If so, save path to list.
        paths = []
        for path in value:
            # ensure it is Path object
            path = Path(path)

            # merged sample indicator
            # ... TODO

            # get full path
            path = path.expanduser().resolve()
            if path.is_dir():
                raise ValueError(
                    "You entered a dir path where you must enter a file path.\n"
                    "To select multiple files in a dir you can use regular \n"
                    "expressions, e.g., './data/*.fastq.gz'.")

            # raise warning if user entered a path but nothing matches it.
            if not list(path.parent.glob(path.name)):
                logger.warning(f"no files match the fastq_paths parameter: {path}")
            paths.append(path)
        return paths

    @field_validator("reference_sequence", "reference_as_filter")
    @classmethod
    def _reference_validator(cls, value: Path) -> Path:
        """Checks that reference file exists and expands path."""
        if value:
            value = value.expanduser().resolve()
            if not value.exists():
                raise ValueError(f"no files match the input string: {value}")
            if value.suffix == ".gz":
                raise ValueError(f"reference {value} must be decompressed.")
        return value

    @field_validator("restriction_overhang")
    @classmethod
    def _enzyme_validator(cls, value):
        """Check that the restriction enzymes do not contain bad chars?"""
        # parse paramsfile str
        if isinstance(value, str):
            if "," in value:
                value = value.strip().split(",")
                print(value)
        return value

    @field_validator("trim_reads")
    @classmethod
    def _trim_reads_validator(cls, value) -> List[int]:
        """Default is (0, 0). User can set to any integers. Negative
        values have special meaning of
        """
        # handle older format with 4 values by selecting the first,last
        if len(value) == 4:
            value = [value[0], value[-1]]
        if len(value) == 1:
            return [value, 0]
        assert len(value) == 2, "trim_reads must contain 2 values, e.g., (0, 0)"
        return list(value)

    @field_validator("trim_loci")
    @classmethod
    def _trim_loci_validator(cls, value) -> Tuple[int, int, int, int]:
        """Return a tuple[int,int]."""
        if value is None:
            return (0, 0, 0, 0)
        assert len(value) == 4, "trim_loci must contain 4 values, e.g., (0, 0, 0, 0)"
        return value

    # TODO...


if __name__ == "__main__":

    import ipyrad as ip

    p = Params(assembly_name="TEST", project_dir="./")
    p.min_depth_majrule = '2'
    p.restriction_overhang = "TGC,".split(',')
    p.project_dir = "/tmp"
    p.reference_as_filter = "../../tests/ipsimdata/gbs_example_genome.fa"
    p.fastq_paths = "../../pedtest/small_tmp_R1.*.gz"
    # p.fastq_paths = ["../../pedtest/small_tmp_R1.*.gz"]
    # p.fastq_paths = [Path('a'), Path('b')]  # "[PosixPath('/home/deren/Documents/ipyrad/pedtest/small_tmp_R1.*.gz')]"
    p.trim_reads = (0, 0)

    # p.max_indels_locus = "0 0".split(",")
    dump = p.model_dump_json()
    print(dump)
    # Params.model_validate_json(dump)

