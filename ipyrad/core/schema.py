#!/usr/bin/env python

"""
JSON schema for storing project w/ Params and Sample info.
"""

# pylint: disable=no-self-argument, no-name-in-module

from typing import List, Tuple, Dict, Optional
from pkg_resources import get_distribution
from pydantic import BaseModel, Field
from ipyrad.core.params_schema import ParamsSchema, HackersSchema


class Stats(BaseModel):
    state: int = Field(None)
    reads_raw: int = Field(None)
    reads_passed_filter: int = Field(None)
    reads_merged: int = Field(None)
    reads_mapped_to_ref_prop: int = Field(None)
    reads_mapped_to_ref_prop: float = Field(None)
    cluster_total: int = Field(None)
    clusters_hidepth: int = Field(None)
    hetero_est: float = Field(None)
    error_est: float = Field(None)
    reads_consens: int = Field(None)
    loci: int = Field(None)

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
    reads_mapped_to_ref: int = Field(None)
    reads_mapped_to_ref_prop: float = Field(None)
    clusters_total: int = Field(None)
    clusters_hidepth: int = Field(None)   
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
    depths_histogram: List[int] = Field(None)

class Stats4(BaseModel):
    hetero_est: float = 0.
    error_est: float = 0.
    min_depth_stat_during_step4: int = 0

class Stats5(BaseModel):
    cluster_total: int
    consensus_total: int
    filtered_by_depth: int
    filtered_by_max_h: int
    filtered_by_max_alleles: int
    filtered_by_max_n: int
    heterozygosity: float
    nsites: int
    nhetero: int
    min_depth_maj_during_step5: int = 0    
    min_depth_stat_during_step5: int = 0

class Stats7(BaseModel):
    nloci: int

class SampleFiles(BaseModel):
    fastqs: List[Tuple[str, str]] = Field(None)
    edits: List[Tuple[str, str]] = Field(None)
    mapped_reads: Tuple[str, str] = Field(None)
    unmapped_reads: Tuple[str, str] = Field(None)
    clusters: str = Field(None)
    consens: str = Field(None)
    depths: str = Field(None)
    database: str = Field(None)

class SampleSchema(BaseModel):
    name: str
    state: int = 1
    files: SampleFiles = SampleFiles()
    stats_s1: Stats1 = Field(None)
    stats_s2: Stats2 = Field(None)
    stats_s3: Stats3 = Field(None)
    stats_s4: Stats4 = Field(None)
    stats_s5: Stats5 = Field(None)
    stats_s7: Stats7 = Field(None)
    # stats: Stats = Field(None)

class FilterStats(BaseModel):
    nloci_before_filtering: int
    filtered_by_rm_duplicates: int
    filtered_by_min_sample_cov: int
    filtered_by_max_indels: int    
    filtered_by_max_snps: int
    filtered_by_max_shared_h: int
    nloci_after_filtering: int

class AssemblyStats(BaseModel):
    filters: FilterStats
    nsnps: int
    nloci: int
    nbases: int
    nsamples: int
    sample_cov: Dict[str, int]
    locus_cov: Dict[int, int]
    var_sites: Dict[int, int]
    var_props: Dict[float, int]
    pis_sites: Dict[int, int]
    pis_props: Dict[float, int]

class Project(BaseModel):
    version: str = str(get_distribution('ipyrad')).split()[1]
    params: ParamsSchema
    hackers: HackersSchema
    samples: Dict[str,SampleSchema]
    assembly_stats: AssemblyStats = Field(None)
    outfiles: Dict[str,str] = Field(None)

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
