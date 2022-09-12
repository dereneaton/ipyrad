#!/usr/bin/env python

"""Utilities for step3 shared by both denovo and reference assemblies.

"""

from typing import Tuple, TypeVar, Iterator
import sys
import gzip
from pathlib import Path
from subprocess import Popen, PIPE, DEVNULL, STDOUT
from loguru import logger
import numpy as np
import pandas as pd

from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import IPyradError

logger = logger.bind(name="ipyrad")
Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")
BIN = Path(sys.prefix) / "bin"
BIN_BWA = str(BIN / "bwa")
BIN_SAMTOOLS = str(BIN / "samtools")  # indexing 
BIN_VSEARCH = str(BIN / "vsearch")    # dereplicating


class ClustMapBase:
    """Reference mapping assembly pipeline."""
    def __init__(self, step: "Step3"):
        # inherit attrs from step
        self.data = step.data
        self.samples = step.samples
        self.quiet = step.quiet
        self.lbview = step.lbview
        self.thview = step.thview

    def run(self):
        """Different methods apply to reference vs denovo."""
        # reference: 
        # [unique] .index_references
        # [shared] .concat_trimmed_files_from_assembly_merge
        # [shared] .mapping_to_reference_filter        
        # [unique] .join_pairs_for_derep
        # [shared] .tag_for_decloning
        # [shared] .dereplicate
        # [shared] .tag_back_to_header
        # [unique] .mapping

        # denovo: 
        # [unique] .index_reference_as_filter
        # [shared] .concat_trimmed_files_from_assembly_merge
        # [shared] .mapping_to_reference_filter            
        # [unique] .pair_merge_overlaps_with_vsearch
        # [unique] .pair_join_unmerged_end_to_end
        # [shared] .tag_for_decloning
        # [shared] .dereplicate
        # [shared] .tag_back_to_header
        # [unique] .cluster

    def concat_trimmed_files_from_assembly_merge(self):
        """Concat when assemblies were merged before step3.

        # i: edits/{sname}_edits_[R1,R2].fastq
        # o: tmpdir/{sname}_concat_edits_[R1,R2].fastq
        """
        # skip if no merge happened.
        if not any(len(self.samples[i].files.edits) > 1 for i in self.samples):
            return
        # else run the job
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            if len(self.samples[sname].files.edits) > 1:
                jobs[sname] = self.lbview.apply(concat_multiple_edits, *args)
        if jobs:
            msg = "concatenating merged assembly inputs"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def mapping_to_reference_filter(self):
        """Map reads to filter reference to get unmapped fastq.

        # i1: edits/{}_edits_R[1,2].fastq
        # i0: tmpdir/{}_concat_edits_R[1,2].fastq
        # o: tmpdir/{}_unmapped_R[1,2].fastq
        """
        if not self.data.params.reference_as_filter:
            return        
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.thview.apply(mapping_reads_minus, *args)
        if jobs:
            msg = "mapping to reference as filter"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

        # store stats
        for sname, val in prog.results.items():
            sample = self.data.samples[sname]
            sample.stats_s3.reads_mapped_to_ref_filter = val[0]
            sample.stats_s3.reads_mapped_to_ref_filter_prop = val[1]
            logger.debug(
                f"{sname} proportion refmap filtered reads: {val[1]:.4f}")            

    def decloning_transfer_tags_inline(self):
        """Moves the i5 tags from the index to start of reads for decloning.

        # i: tmpdir/{}_merged.fa
        # o: tmpdir/{}_decloned.fa
        """
        if (not self.data.is_pair) or (not self.data.hackers.declone_PCR_duplicates):
            return
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(tag_seq_for_decloning, *args)
        if jobs:
            msg = "tagging reads for decloning    "
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check()

    def dereplicate(self):
        """Dereplicate sequences (read pairs are merged).

        # i3: edits/{}_edits.fastq                 # se data
        # i2: tmpdir/{}_concat_edits.fastq         # se assembly merge
        # i1: tmpdir/{}_merged.fa                  # pe data
        # i0: tmpdir/{}_decloned.fa                # pe w/ declone
        # o: tmpdir/{}_derep.fa
        """
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.thview.apply(dereplicate, *args)
        msg = "dereplicating"
        prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
        prog.block()
        prog.check()

    def decloning_transfer_tags_to_header(self):
        """Moves the i5 tags from inline back to the header.

        # i: tmpdir/{}_derep.fa
        # o: tmpdir/{}_derep_tag.fa
        """
        if (not self.data.is_pair) or (not self.data.hackers.declone_PCR_duplicates):
            return
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(retag_header_after_derep_for_decloning, *args)
        if jobs:
            msg = "moving tags to header"
            prog = AssemblyProgressBar(jobs, msg, step=3, quiet=self.quiet)
            prog.block()
            prog.check() 

    def calculate_sample_stats(self):
        """Send samples to calc depths on remote, and then enter stats
        to sample objects non-parallel.
        """
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            jobs[sname] = self.lbview.apply(set_sample_stats, *args)

        # progress bar
        msg = "calculating stats"
        prog = AssemblyProgressBar(jobs, msg, 3, self.quiet)
        prog.block()
        prog.check()

        # store/save
        self.data.samples = prog.results
        self.data.save_json()

        # write stats text file
        handle = self.data.stepdir / "s3_cluster_stats.txt"
        with open(handle, 'w', encoding="utf-8") as out:
            table = pd.DataFrame({
                i: self.data.samples[i].stats_s3.dict() 
                for i in self.data.samples
            }).T
            table.sort_index(inplace=True)
            table.to_string(out)
            cols=["clusters_total", "clusters_hidepth", "mean_depth_stat"]
            logger.info("\n" + table.loc[:, cols].to_string())

    # def calculate_sample_stats(self):
    #     """Send samples to calc depths on remote, and then enter stats
    #     to sample objects non-parallel.
    #     """
    #     jobs = {}
    #     for sname in self.samples:
    #         jobs[sname] = self.lbview.apply(
    #             set_sample_stats,
    #             *(self.data, self.samples[sname])
    #         )

    #     # progress bar
    #     msg = "calculating stats"
    #     prog = AssemblyProgressBar(jobs, msg, 3, self.quiet)
    #     prog.block()
    #     prog.check()

    #     # store/save
    #     self.data.samples = prog.results  # {i: j for (i, j) in prog.results.items()}
    #     self.data.save_json()

    #     # write stats text file
    #     statsdf = pd.DataFrame(
    #         index=sorted(self.data.samples),
    #         columns=[
    #             "clusters_total",
    #             "clusters_hidepth",
    #             "mean_depth_mj",
    #             "mean_depth_stat",
    #             "reads_mapped_to_ref_prop",
    #         ],
    #     )
    #     for sname in self.data.samples:
    #         statsdict = self.data.samples[sname].stats_s3.dict()
    #         for i in statsdf.columns:
    #             statsdf.loc[sname, i] = statsdict[i]

    #     handle = os.path.join(self.data.stepdir, 's3_cluster_stats.txt')
    #     with open(handle, 'w', encoding="utf-8") as outfile:
    #         statsdf.to_string(outfile)
    #     logger.info("\n" + statsdf.to_string())


def concat_multiple_edits(data: Assembly, sample: Sample) -> None:
    """Concat files if multiple Assemblies were merged between steps 2-3.

    Create a temporary concatenated file for multiple edits input
    files, which arises when Assemblies were merged between steps
    2 and 3.
    """
    # define output files
    concat1 = data.tmpdir / f"{sample.name}_concat_R1.fastq.gz"
    concat2 = data.tmpdir / f"{sample.name}_concat_R2.fastq.gz"

    read1s = [i[0] for i in sample.files.edits]
    if len(read1s) > 1:
        cmd = ['cat'] + read1s
        with open(concat1, 'w', encoding="utf-8") as cout:
            with Popen(
                cmd, stderr=PIPE, stdout=cout, close_fds=True) as proc:
                res = proc.communicate()[1]
                if proc.returncode:
                    raise IPyradError(f"cmd: {' '.join(cmd)}\nerror: {res}")

    read2s = [i[1] for i in sample.files.edits if i[1]]
    if len(read2s) > 1:
        cmd = ['cat'] + read2s
        with open(concat2, 'w', encoding="utf-8") as cout:
            with Popen(cmd, stderr=PIPE, stdout=cout, close_fds=True) as proc:
                res = proc.communicate()[1]
                if proc.returncode:
                    raise IPyradError(f"cmd: {' '.join(cmd)}\nerror: {res}")


def mapping_reads_minus(data: Assembly, sample: Sample) -> Tuple[int, float]:
    """Map reads to the reference-filter fasta.

    Mapped reads are discarded and unmapped reads are kept as data.
    """
    reference = Path(data.params.reference_as_filter)
    if not reference.exists():
        raise IPyradError(f"reference_filter sequence not found: {reference}")

    # input reads are concat if present else trims
    read1 = data.tmpdir / f"{sample.name}_concat_R1.fastq.gz"
    read2 = data.tmpdir / f"{sample.name}_concat_R2.fastq.gz"
    read1 = read1 if read1.exists() else Path(sample.files.edits[0][0])
    read2 = read2 if read2.exists() else Path(sample.files.edits[0][1])

    # setup cmd1 (mapping w/ bwa)
    nthreads = max(1, data.ipcluster["threads"])
    cmd1 = [BIN_BWA, "mem", "-t", str(nthreads), "-M", str(reference)]
    cmd1 += [str(read1), str(read2) if read2 else ""]
    if data.hackers.bwa_args:
        for arg in data.hackers.bwa_args.split()[::-1]:
            cmd1.insert(2, arg)
    cmd1 += ['-o', str(data.tmpdir / f"{sample.name}_ref_filter.sam")]
    print(" ".join(cmd1))

    # run cmd1
    with Popen(cmd1, stderr=PIPE, stdout=DEVNULL) as proc1:
        error1 = proc1.communicate()[1].decode()
        if proc1.returncode:
            raise IPyradError(f"cmd: {' '.join(cmd1)}\nerror: {error1}")

    # setup cmd2 (sam to bam)
    cmd2 = [BIN_SAMTOOLS, 'view', '-b', '-F', '0x904']
    cmd2 += ['-U', str(data.tmpdir / f"{sample.name}_unmapped.bam")]
    cmd2 += [str(data.tmpdir / f"{sample.name}_ref_filter.sam")]
    print(' '.join(cmd2))

    # run cmd2
    with Popen(cmd2, stderr=PIPE, stdout=DEVNULL) as proc2:
        error2 = proc2.communicate()[1].decode()
        if proc2.returncode:
            raise IPyradError(f"cmd: {' '.join(cmd2)}\nerror: {error2}")

    # setup cmd3 (bam to fastq unmapped)
    cmd3 = [BIN_SAMTOOLS, 'fastq', '-v', '45']
    cmd3 += ['-1', str(data.tmpdir / f"{sample.name}_unmapped_R1.fastq")]
    cmd3 += ['-2', str(data.tmpdir / f"{sample.name}_unmapped_R2.fastq")]
    cmd3 += [str(data.tmpdir / f"{sample.name}_unmapped.bam")]
    print(' '.join(cmd3))

    # run cmd3
    with Popen(cmd3, stderr=PIPE, stdout=DEVNULL) as proc3:
        error3 = proc3.communicate()[1].decode()
        if proc3.returncode:
            raise IPyradError(f"cmd: {' '.join(cmd3)}\nerror: {error3}")

    # return the number of filtered reads
    unmapped = data.tmpdir / data.tmpdir / f"{sample.name}_unmapped_R1.fastq"
    with open(unmapped, 'r', encoding="utf-8") as ion:
        n_unmapped = int(sum(1 for i in ion) / 4)
    n_filtered = sample.stats_s2.reads_passed_filter - n_unmapped
    n_filtered_prop = n_filtered / sample.stats_s2.reads_passed_filter
    return n_filtered, n_filtered_prop


def tag_seq_for_decloning(data: Assembly, sample: Sample) -> None:
    """Tag reads w/ i5 inline prior to dereplicating.

    Pull i5 tag from Illumina index and insert into the fastq name
    header so we can use it later to declone even after the fastq
    indices are gone. THIS IS DONE BEFORE DEREPLICATION, so that
    identical reads are not collapsed if they have different i5s.

    3rad uses random adapters to identify pcr duplicates. For removing
    pcr dups later we need to insert the tag into the sequence here
    so they will not be dereplicated together.

    >>>                                                    i7       i5
    >>> >NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA
    >>> AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    the i5 is inserted to the read with "F" quality score

    >>> ********
    >>> >NB551405:60:H7T2GAFXY:1:11101:24455:4008 1:N:0:TATCGGTC+CAAGACAA
    >>> CAAGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    """
    # The input files are different depending on ref or not.
    if data.is_ref:
        tmpin = data.tmpdir / f"{sample.name}_joined.fastq"
    else:
        tmpin = data.tmpdir / f"{sample.name}_merged.fa"
    tmpout = data.tmpdir / f"{sample.name}_decloned.fa"

    # Remove adapters from head of sequence and write out
    # tmp_outfile is now the input file for the next step
    # first vsearch derep discards the qscore so we iterate pairs
    with open(tmpout, 'w', encoding="utf-8") as out:
        with open(tmpin, 'r', encoding="utf-8") as infile:

            # iterate over 2 lines a time
            duo = zip(*[infile] * 2)

            # a list to store until writing
            tmp = []
            for idx, read in enumerate(duo):

                # extract i5 if it exists else use empty string
                try:
                    indices = read[0].split(":")[-1]
                    index = indices.split("+")[1].strip()
                    # assert int(index), "failed to parse i5 index"
                except (IndexError, AssertionError):
                    index = ""

                # add i5 to the 5' end of the sequence and insert the
                # length of the i5 index (usually 8) into the header.
                newread = f"{read[0]}{index}{read[1]}"
                tmp.append(newread)

                # Write the data in chunks
                if not idx % 50_000:
                    out.write("".join(tmp))
                    tmp = []
            if tmp:
                out.write("".join(tmp))
    print(f"tagged with inline i5s for decloning: {sample.name}")
    tmpin.unlink()


def dereplicate(data: Assembly, sample: Sample) -> None:
    """Dereplicate reads and sort so that most replicated are on top.

    Paired data are dereplicated as joined reads.
    """
    # select input file with following precedence:
    # i3: edits/{}_edits.fastq                 # se data
    # i2: tmpdir/{}_concat_edits.fastq         # se assembly merge
    # i1: tmpdir/{}_merged.fa                  # pe data
    # i0: tmpdir/{}_decloned.fa                # pe w/ declone
    # o: tmpdir/{}_derep.fa
    infiles = [
        Path(sample.files.edits[0][0]),
        data.tmpdir / f"{sample.name}_concat_edits.fastq",
        data.tmpdir / f"{sample.name}_merged.fa",
        data.tmpdir / f"{sample.name}_decloned.fa",
    ]
    infiles = [i for i in infiles if i.exists()]
    infile = infiles[-1]

    # datatypes options
    strand = "plus"
    if data.params.datatype in ['gbs', 'pairgbs', '2brad']:
        strand = "both"

    # do dereplication with vsearch
    cmd = [
        BIN_VSEARCH,
        "--fastx_uniques", str(infile),
        "--strand", strand,
        "--fastaout", str(data.tmpdir / f"{sample.name}_derep.fa"),
        "--fasta_width", str(0),
        "--minseqlength", str(data.params.filter_min_trim_len),
        "--sizeout",
        "--relabel_md5",
        "--quiet",
        # "--threads", str(4),
        # "--fastq_qmax", "1000",
    ]

    # decompress argument (IF ZLIB is missing this will not work!!)
    # zlib is part of the conda installation.
    if infile.suffix == ".gz":
        cmd.append("--gzip_decompress")

    # build PIPEd job
    print(" ".join(cmd))  # engine sends to logger.info
    with Popen(cmd, stderr=STDOUT, stdout=PIPE, close_fds=True) as proc:
        errmsg = proc.communicate()[0]
        if proc.returncode:
            raise IPyradError(errmsg.decode())


def retag_header_after_derep_for_decloning(data: Assembly, sample: Sample) -> None:
    """Move inline i5 tag from start of seq back to header

    Example
    -------
    # input data derep file format
    >594732b799a25eb9b8ab4925f3a9a992;size=8
    GGGGGGGGATCGGAAGCACATACTATAATAAGGGGTAGGGTTTTATTGGCAGCAT
    >a9154e2d5348a59230c5ecd19e0afdf6;size=6
    GGGGGGGGATCGGTGCATTCCCCCAAGGGTGTCCTAAAGTTCCTCCACCAAACTG
    ******** <- i5 tag

    # output data derep_tag file
    >GGGGGGGG_594732b799a25eb9b8ab4925f3a9a992;size=8
    ATCGGAAGCACATACTATAATAAGGGGTAGGGTTTTATTGGCAGCATATTCAATC
    >GGGGGGGG_a9154e2d5348a59230c5ecd19e0afdf6;size=6
    ATCGGTGCATTCCCCCAAGGGTGTCCTAAAGTTCCTCCACCAAACTGTAGTACAG
    """
    # paired reads are merged or joined in the merged file
    # tmpin = data.tmpdir / f"{sample.name}_joined.fastq"
    tmpin = data.tmpdir / f"{sample.name}_derep.fa"
    tmpout = data.tmpdir / f"{sample.name}_derep_tag.fa"

    with open(tmpout, 'w', encoding="utf-8") as out:
        with open(tmpin, 'r', encoding="utf-8") as infile:

            # iterate over 2 lines a time
            duo = zip(*[infile] * 2)

            # a list to store until writing
            tmp = []
            for idx, read in enumerate(duo):

                # extract i5 from inline (assumes len=8 !!!)
                barcode, seq = read[1][:8], read[1][8:]

                # add i5 to the 5' end of the sequence and insert the
                # length of the i5 index (usually 8) into the header.
                newread = f">{barcode}_{read[0][1:]}{seq}"
                tmp.append(newread)

                # Write the data in chunks
                if not idx % 50_000:
                    out.write("".join(tmp))
                    tmp = []
            if tmp:
                out.write("".join(tmp))
    print(f"moved i5 tags to headers after derep: {sample.name}")


def iter_clusters(clustfile: Path, gzipped: bool=False) -> Iterator[str]:
    """Yields clusters between //\n// separators."""
    open_func = gzip.open if gzipped else open
    with open_func(clustfile, 'rt', encoding="utf-8") as clustio:
        data = []
        pairs = zip(*[iter(clustio)] * 2)
        for name, seq in pairs:
            if name[0] == ">":
                data.extend([name, seq])
            else:
                yield data
                data = []


def set_sample_stats(data: Assembly, sample: Sample) -> Sample:
    """Sets step3 stats to Samples from clusters files.

    This is used for both denovo and reference assemblies.
    Iterates over clustS files to count data, returns maxlen and
    depths arrays for each sample.
    """
    clustfile = data.stepdir / f"{sample.name}.clusters.gz"
    sample.state = 3
    sample.files.clusters = clustfile

    # get depths and lens distributions from clusters file.
    depths = []   # read depth: sum of 'sizes'
    clens = []    # lengths of clusters
    for clust in iter_clusters(sample.files.clusters, gzipped=True):
        names = clust[::2]
        sizes = [int(i.split(";")[-2][5:]) for i in names]
        depths.append(sum(sizes))
        clens.append(len(clust[1].strip()))
    clens, depths = np.array(clens), np.array(depths)

    # sample does not advance state to 3
    if not depths.size:
        sample.stats_s3.clusters_total = 0
        sample.stats_s3.clusters_hidepth = 0
        return sample

    # create mindepth masks
    maj_mask = depths >= data.params.min_depth_majrule
    hid_mask = depths >= data.params.min_depth_statistical

    # store length stats
    hilens = clens[hid_mask]
    sample.stats_s3.max_hidepth_cluster_length = int(hilens.max())
    sample.stats_s3.mean_hidepth_cluster_length = float(hilens.mean())
    sample.stats_s3.std_hidepth_cluster_length = float(hilens.std())

    # store n clusters stats
    sample.stats_s3.clusters_total = int(depths.size)
    sample.stats_s3.clusters_hidepth = int(depths[hid_mask].size)

    # store depths histogram as a dict. Limit to first 25 bins
    bars, _ = np.histogram(depths, bins=range(1, 26))
    sample.stats_s3.depths_histogram = [int(i) for i in bars]
    sample.stats_s3.mean_depth_total = float(depths.mean())
    sample.stats_s3.mean_depth_mj = float(depths[maj_mask].mean())
    sample.stats_s3.mean_depth_stat = float(depths[hid_mask].mean())
    sample.stats_s3.std_depth_total = float(depths.std())
    sample.stats_s3.std_depth_mj = float(depths[maj_mask].std())
    sample.stats_s3.std_depth_stat = float(depths[hid_mask].std())

    return sample
