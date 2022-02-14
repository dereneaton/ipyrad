#!/usr/bin/env python

"""JSON schema for storing project w/ Params and Sample info.

Pydantic Models are similar to dataclases but they also include
*type validation*, meaning that if you try to set an attribute to
the wrong type it will raise an error.
"""

# pylint: disable=no-self-argument, no-name-in-module

from typing import List, Tuple, Dict, Optional
from pathlib import Path
from pkg_resources import get_distribution
from pydantic import BaseModel, Field
from ipyrad.core.params_schema import ParamsSchema, HackersSchema


class Stats(BaseModel):
    state: int = None
    reads_raw: int = None
    reads_passed_filter: int = None
    reads_merged: int = None
    reads_mapped_to_ref_prop: int = None
    reads_mapped_to_ref_prop: float = None
    cluster_total: int = None
    clusters_hidepth: int = None
    hetero_est: float = None
    error_est: float = None
    reads_consens: int = None
    loci: int = None

class Stats1(BaseModel):
    reads_raw: int = 0

class Stats2(BaseModel):
    reads_raw: int = 0
    reads_filtered_by_Ns: int = 0
    reads_filtered_by_low_quality: int = 0
    reads_filtered_by_low_complexity: int = 0
    reads_filtered_by_minlen: int = 0
    mean_len_R1_before_trimming: float = 0.
    mean_len_R1_after_trimming: float = 0.
    mean_len_R2_before_trimming: Optional[float] = None
    mean_len_R2_after_trimming: Optional[float] = None
    reads_passed_filter: int = 0

class Stats3(BaseModel):
    merged_pairs: int = 0
    min_depth_maj_during_step3: int = 0
    min_depth_stat_during_step3: int = 0
    reads_mapped_to_ref: int = None
    reads_mapped_to_ref_prop: float = None
    clusters_total: int = None
    clusters_hidepth: int = None
    mean_depth_total: float = 0.
    mean_depth_mj: float = 0.
    mean_depth_stat: float = 0.
    std_depth_total: float = 0.
    std_depth_mj: float = 0.
    std_depth_stat: float = 0.
    filtered_bad_align: int = 0
    deduplicated_reads: int = 0
    deduplicated_reads_prop: int = 0
    max_hidepth_cluster_length: int = 0
    mean_hidepth_cluster_length: float = 0.
    std_hidepth_cluster_length: float = 0.
    depths_histogram: List[int] = None

class Stats4(BaseModel):
    hetero_est: float = 0.
    error_est: float = 0.
    min_depth_stat_during_step4: int = 0

class Stats5(BaseModel):
    cluster_total: int = None
    consensus_total: int = None
    filtered_by_depth: int = None
    filtered_by_max_h: int = None
    filtered_by_max_alleles: int = None
    filtered_by_max_n: int = None
    heterozygosity: float = None
    nsites: int = None
    nhetero: int = None
    min_depth_maj_during_step5: int = None
    min_depth_stat_during_step5: int = None

class Stats7(BaseModel):
    nloci: int = 0

class SampleFiles(BaseModel):
    """A dict-like class for storing file paths to samples."""
    fastqs: List[Tuple[str, str]] = None
    edits: List[Tuple[str, str]] = None
    mapped_reads: Tuple[str, str] = None
    unmapped_reads: Tuple[str, str] = None
    clusters: str = None
    consens: str = None
    depths: str = None
    database: str = None

class SampleSchema(BaseModel):
    """A dict-like class for representing individual Samples info."""
    name: str
    state: int = 1
    files: SampleFiles = SampleFiles()
    stats_s1: Stats1 = None
    stats_s2: Stats2 = None
    stats_s3: Stats3 = None
    stats_s4: Stats4 = None
    stats_s5: Stats5 = None
    stats_s7: Stats7 = None

class FilterStats(BaseModel):
    """A dict-like class with stats for filters applied during step 7."""
    nloci_before_filtering: int = 0
    filtered_by_rm_duplicates: int = 0
    filtered_by_min_sample_cov: int = 0
    filtered_by_max_indels: int = 0
    filtered_by_max_snps: int = 0
    filtered_by_max_shared_h: int = 0
    nloci_after_filtering: int = 0

class AssemblyStats(BaseModel):
    """A dict-like class with stats for completed assembly."""
    filters: FilterStats = FilterStats()
    nsnps: int = 0
    nloci: int = 0
    nbases: int = 0
    nsamples: int = 0
    sample_cov: Dict[str, int] = Field(default_factory=dict)
    locus_cov: Dict[int, int] = Field(default_factory=dict)
    var_sites: Dict[int, int] = Field(default_factory=dict)
    var_props: Dict[float, int] = Field(default_factory=dict)
    pis_sites: Dict[int, int] = Field(default_factory=dict)
    pis_props: Dict[float, int] = Field(default_factory=dict)

class Project(BaseModel):
    """Top-level schema for serializing Assembly objects to JSON and back.

    When Assembly objects calls `.save_json()` the object is serialized
    to JSON and written to disk. In each Step of the assembly in `.run`
    a Project is created from the Assembly and its Sample objects at
    that point to update it and re-serialize to JSON on disk.
    """
    version: str = str(get_distribution('ipyrad')).split()[1]
    params: ParamsSchema
    hackers: HackersSchema
    samples: Dict[str, SampleSchema]
    assembly_stats: AssemblyStats = None
    outfiles: Dict[str, Path] = Field(default_factory=dict)
    populations: Dict[str, Tuple[List[str], int]] = Field(default_factory=dict)

    def __str__(self):
        return self.json(indent=2)

    def __repr__(self):
        return self.json(indent=2)


if __name__ == "__main__":

    # user input params
    params = ParamsSchema(assembly_name="test").dict()
    hackers = HackersSchema()

    # step 1 creates a null sample
    samp1 = SampleSchema(name="A").dict()

    # collect stats for this step
    files = [('A_r1', 'A_r2')]
    stats = Stats1(reads_raw = 10000)

    # load back to schema
    samp1['files']['fastqs'] = files
    samp1['stats_s1'] = stats
    samp1 = SampleSchema(**samp1)
    print(samp1.json(indent=2, exclude_none=True))

    samples = [samp1, samp1]

    # step1 ends with Assembly creation.
    data = Project(
        params=ParamsSchema(**params),
        hackers=hackers,
        samples={sample.name: sample for sample in samples},
    )

    print(data.json(indent=2, exclude_none=True))
