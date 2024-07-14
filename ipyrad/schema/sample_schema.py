#!/usr/bin/env python

"""Serializable JSON schema for Sample object

PLANNED JSON SCHEMA:
---------------------
{
    name: 'a'',
    state: 2,
    files: {
        fastqs: [a_R1.fastq, a_R2.fastq],
        trims: [a_R1.trimmed.fastq, b_R2.trimmed.fastq],
        clusters: ...
        ...,
    },
    stats_s1: {
        ...,
    },
    stats_s2: {
        ...,
    },
    ...
}
"""

from typing import List, Tuple
from pathlib import Path
from pydantic import BaseModel

__all__ = [
    "Stats1",
    "Stats2",
    "Stats3",
    "Stats4",
    "Stats5",
    "Stats6",
    "Stats7",
    "Sample",
    "SampleFiles",
]


class Stats1(BaseModel):
    reads_raw: int = 0
    reads_filtered_by_Ns: int = 0
    reads_filtered_by_low_quality: int = 0
    reads_filtered_by_low_complexity: int = 0
    reads_filtered_by_minlen: int = 0
    adapter_trimmed_reads: int = 0
    adapter_trimmed_bases: int = 0
    mean_len_R1_before_trimming: float = 0.
    mean_len_R1_after_trimming: float = 0.
    mean_len_R2_before_trimming: float = 0.
    mean_len_R2_after_trimming: float = 0.
    reads_passed_filter: int = 0


class Stats2(BaseModel):
    """Results of clustering/mapping"""
    merged_pairs: int = 0
    merged_pairs_prop: float = 0
    reads_mapped_to_ref: int = None
    reads_mapped_to_ref_prop: float = None
    reads_mapped_to_ref_filter: int = None
    reads_mapped_to_ref_filter_prop: float = None
    dereplicated_reads: int = None
    cluster_seeds: int = None


class Stats3(BaseModel):
    """Results of cluster building."""
    clusters: int = 0
    cluster_depth_mean: float = 0.
    cluster_depth_median: float = 0.
    cluster_depth_std: float = 0.
    cluster_depth_histogram: List[Tuple[int, int]] = None
    filtered_bad_alignment: int = 0
    # pcr_duplicates: int = 0
    # pcr_duplicates_prop: float = 0

    # min_depth_maj_during_step3: int = 0
    # min_depth_stat_during_step3: int = 0
    # max_hidepth_cluster_length: int = 0
    # mean_hidepth_cluster_length: float = 0.
    # std_hidepth_cluster_length: float = 0.
    # depths_histogram: List[int] = None
    # clusters_total: int = None
    # clusters_hidepth: int = None
    # mean_depth_total: float = 0.
    # mean_depth_mj: float = 0.
    # mean_depth_stat: float = 0.
    # std_depth_total: float = 0.
    # std_depth_mj: float = 0.
    # std_depth_stat: float = 0.
    # deduplicated_reads: int = 0
    # deduplicated_reads_prop: int = 0


class Stats4(BaseModel):
    hetero_est: float = 0.
    error_est: float = 0.
    min_depth_stat_during_step4: int = 0


class Stats5(BaseModel):
    clusters_total: int = None
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


class Stats6(BaseModel):
    nloci_prefiltered: int = 0
    nsamples: int = 0


class Stats7(BaseModel):
    nloci: int = 0


class SampleFiles(BaseModel):
    fastqs: List[Tuple[Path, Path]] = None
    """: The fastq data files for this Sample."""
    trimmed: Tuple[Path, Path] = None
    """: pairs of fastq files that have been filtered."""
    clustmap: Tuple[Path, Path] = None
    """: denovo clustering: (derep, clustmap) tuple where dereps are R1nnnnR2."""
    bamfile: Path = None
    """: reference clustering: this is R1/R2 pairs mapped to ref."""
    clusters: Path = None
    """: fasta-like clusters."""
    consens: Path = None
    """: fasta consensus allele calls."""
    depths: Path = None
    """: hdf5 database storing site depths."""
    database: Path = None
    """: fasta alignments of loci across samples, not yet filtered."""


class Sample(BaseModel):
    """A dict-like class for representing individual Samples info."""
    name: str
    state: int = 1
    files: SampleFiles = SampleFiles()
    stats_s1: Stats1 = None
    stats_s2: Stats2 = None
    stats_s3: Stats3 = None
    stats_s4: Stats4 = None
    stats_s5: Stats5 = None
    stats_s6: Stats6 = None
    stats_s7: Stats7 = None

    def __str__(self):
        return self.model_dump_json(indent=2)


if __name__ == "__main__":

    import ipyrad as ip
    sample = Sample(name="A")
    sample.files.fastqs = [(Path('a1'), Path('a2'))]
    sample.files.fastqs = [(Path('/a1/x.fastq'), "null")]
    sample.stats_s1 = Stats1(reads_raw=1000)
    json = sample.model_dump_json(indent=2, exclude_none=True)
    print(json)
    with open("/tmp/_test.json", 'w', encoding="utf-8") as out:
        out.write(sample.model_dump_json(indent=2, exclude_none=True))
    s = Sample.model_validate_json(Path("/tmp/_test.json").read_text())
    print(s)
    # save assembly with this sample and re-load it.
    # params = ip.schema.Params(assembly_name="test").model_dump()
    # data = ip.schema.Project(
    #     params=ip.schema.Params(**params),
    #     samples={"A": sampleA.model_dump()},
    # )
    # data.save_json()
    # ip.load_json(...)

    # sample1 = Sample(name="A")
    # sample1.files.fastqs = ['a', 'b']

    # sample2 = Sample(name="B")
    # sample2.files.fastqs = ['a', 'b']

    # samples = [sample1, sample2]
    # params = Params(assembly_name="test").model_dump()

    # # step1 ends with Assembly creation.
    # data = Project(
    #     params=Params(**params),
    #     samples={sample.name: sample for sample in samples},
    # )

    # print(data)
