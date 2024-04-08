#!/usr/bin/env python

"""...


"""

from typing import Dict, Tuple, List
from pathlib import Path
from importlib.metadata import distribution
from pydantic import BaseModel, Field
from ipyrad.schema.sample_schema import Sample
from ipyrad.schema.params_schema import Params

VERSION = distribution("ipyrad").version


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
    nsites: int = 0
    nsamples: int = 0
    sample_cov: Dict[str, int] = Field(default_factory=dict)
    locus_cov: Dict[int, int] = Field(default_factory=dict)
    var_sites: Dict[int, int] = Field(default_factory=dict)
    var_props: Dict[float, int] = Field(default_factory=dict)
    pis_sites: Dict[int, int] = Field(default_factory=dict)
    pis_props: Dict[float, int] = Field(default_factory=dict)


class StatsFiles(BaseModel):
    s1: Path = None
    s2: Path = None
    s3: Path = None
    s4: Path = None
    s5: Path = None
    s6: Path = None
    s7: Path = None

    def _clear_old_results(self, step: int) -> None:
        """Remove any existing results after the current state step."""
        for step in range(step + 1, 8):
            setattr(self, f"s{step}", None)


class Project(BaseModel):
    """Top-level schema for serializing Assembly objects to JSON and back.

    When Assembly objects calls `.save_json()` the object is serialized
    to JSON and written to disk. In each Step of the assembly in `.run`
    a Project is created from the Assembly and its Sample objects at
    that point to update it and re-serialize to JSON on disk.
    """
    version: str = VERSION
    params: Params
    samples: Dict[str, Sample]
    assembly_stats: AssemblyStats = None
    outfiles: Dict[str, Path] = Field(default_factory=dict)
    populations: Dict[str, Tuple[List[str], int]] = Field(default_factory=dict)
    stats_files: StatsFiles = StatsFiles()

    def __str__(self):
        return self.model_dump_json(indent=2)

    def __repr__(self):
        return self.model_dump_json(indent=2)


if __name__ == "__main__":

    sample1 = Sample(name="A")
    sample1.files.fastqs = [('a', 'b')]

    sample2 = Sample(name="B")
    sample2.files.fastqs = [('a', 'b')]

    samples = [sample1, sample2]
    params = Params(assembly_name="test").model_dump()

    # step1 ends with Assembly creation.
    data = Project(
        params=Params(**params),
        samples={sample.name: sample for sample in samples},
    )

    print(data)
    print(data.model_dump())
    data.save_json()
