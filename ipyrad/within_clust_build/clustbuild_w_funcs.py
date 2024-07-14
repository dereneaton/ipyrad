#!/usr/bin/env python

"""Functions associated with clustbuild_w.py

"""

from typing import Iterator, List, Tuple, Dict, TypeVar
import sys
import gzip
import itertools
import subprocess
from pathlib import Path
from collections import Counter
from loguru import logger
import numpy as np
from ipyrad.core.utils import comp
from ipyrad.core.exceptions import IPyradError
from ipyrad.schema import Sample

# testing shows spacers improve alignment edges for variable length
# input sequences and make it simple to handle empty halves of paired
# sequences. Padding is easily removed in later steps.
SPACER = "NNNNNNNNNN"
CHUNKSIZE = 2_000
BIN = Path(sys.prefix) / "bin"
BIN_MUSCLE = str(BIN / "muscle")
Assembly = TypeVar("Assembly")
logger = logger.bind(name="ipyrad")


def fill_match_dicts(matches: Path) -> Tuple[Dict, Dict]:
    """Return dictionaries for building clusters.

    hits_to_seed = {match_name: seed_name, ...}
    seeds_to_nhits = {seed_name: int size of seed's cluster, ...}
    """
    hits_to_seeds = {}
    seeds_to_nhits = Counter()
    # iterate over all [hit, seed, ...] rows.
    with open(matches, 'r', encoding="utf-8") as infile:
        for line in infile:
            hit, seed, _, _, strand, _ = line.strip().split()
            # store {hit: (seed, strand), ...}
            hits_to_seeds[hit] = (seed, strand)
            # store expected length of each cluster
            seeds_to_nhits[seed] += 1
    # add one more hit to each cluster for the seed itself
    for seed in seeds_to_nhits:
        seeds_to_nhits[seed] += 1
    return hits_to_seeds, seeds_to_nhits


def iter_derep(derep: Path) -> Iterator[Tuple[str, str]]:
    with open(derep, 'r', encoding="utf-8") as infile:
        for header, sequence in zip(infile, infile):
            yield header[1:].strip(), sequence.strip()


# def iter_align_muscle_chunks(sample: Sample, chunksize: int) -> Tuple[bool, str]:
#     """Generator of (bool, cluster) tuples for writing unaligned tmp files.

#     Separates singletons that don't need to be aligned from clusters
#     that do need muscle alignment. Yields clusters in chunks to be
#     written to tmpfiles.
#     """
#     # store singletons that do not need to be aligned
#     singletons = []
#     to_align = []

#     # iterate over unaligned clusters
#     for clust in iter_build_clusters(sample):

#         # if cluster is singleton then yield right away
#         if len(clust) == 1:
#             singletons.append(clust)

#         # else, send it to be aligned asynchronously
#         else:
#             to_align.append(clust)

#         # yield to_align chunk once it reaches CHUNKSIZE
#         if len(to_align) > CHUNKSIZE:
#             yield True, to_align
#             to_align = []

#     # finally, yield the unaligned chunk
#     yield False, singletons


def iter_build_clusters(sample: Sample, max_size: int = 50) -> Iterator[Tuple[bool, str]]:
    """Generator of unaligned clusters in fasta format.

    Parameters
    ----------
    max_size: int
        The max number of unique sequences to write to a cluster. This
        sets an upper limit to keep muscle alignment reasonably fast. If
        there are many singletons that represent errors to an already
        high depth cluster then it is preferred to just exclude them.
        Derep reads are already sorted by depth,...
    """
    clusters = {}
    hits_to_seeds, seeds_to_nhits = fill_match_dicts(sample.files.clustmap[1])

    # iterate over all sequences in the derep file
    for name, seq in iter_derep(sample.files.clustmap[0]):

        # it is a match to the centroid of a diff cluster seed
        if name in hits_to_seeds:
            seed, orient = hits_to_seeds[name]
            if orient == "-":
                seq = comp(seq)[::-1]
                newname = f">{name};-"
            else:
                newname = f">{name};+"

            # hit is same size as seed, but hit came up first.
            if seed not in clusters:
                clusters[seed] = [newname, seq]
            else:
                clusters[seed].extend([newname, seq])

        # it is a seed to the cluster centroid
        else:
            seed = name
            newname = f">{name};*"
            if seed not in clusters:
                clusters[seed] = [newname, seq]
            else:
                clusters[seed] = [newname, seq] + clusters[seed]

        # singleton (just seed, no hits)
        if seed not in seeds_to_nhits:
            yield False, "\n".join(clusters[seed][:max_size * 2])
            del clusters[seed]

        # yield cluster if it is complete
        else:
            if len(clusters[seed]) == seeds_to_nhits[seed] * 2:
                yield True, "\n".join(clusters[seed][:max_size * 2])
                del clusters[seed]

    # everything should hvae been yielded.
    assert not clusters


def write_tmp_clusters(tmpdir: Path, sample: Sample) -> None:
    r"""Write unaligned clusters to tmpdir as fasta with delim=//\n//.

    Gets clusters from iter_build_clusters generator func. Clusters
    are separated by //\\n//\\n because this acts as a unique delimiter
    by the muscle aligner function.
    """
    unaligned = []
    aligned = []
    chunk_idx = 1
    cidx = 0
    uidx = 0
    for to_align, clust in iter_build_clusters(sample):
        # print(clust)
        if to_align:
            unaligned.append(clust)
            cidx += 1
        else:
            aligned.append(clust)
            uidx += 1

        # occasionally write/dump stored clusters to file and clear mem
        if unaligned:
            if not cidx % CHUNKSIZE:
                unaligned_handle = tmpdir / f"{sample.name}_unaligned_{chunk_idx}.fa"
                with open(unaligned_handle, 'w', encoding="utf-8") as clustio:
                    clustio.write("\n//\n//\n".join(unaligned) + "\n//\n//\n")
                unaligned = []
                chunk_idx += 1

    # write any remaining unaligned
    if unaligned:
        unaligned_handle = tmpdir / f"{sample.name}_unaligned_{chunk_idx}.fa"
        with open(unaligned_handle, 'w', encoding="utf-8") as clustio:
            clustio.write("\n//\n//\n".join(unaligned) + "\n//\n//\n")

    # write singletons that don't need to be aligned
    if aligned:
        aligned_handle = tmpdir / f"{sample.name}_aligned_0.fa"
        with open(aligned_handle, 'w', encoding="utf-8") as clustio:
            clustio.write("\n//\n//\n".join(aligned) + "\n//\n//\n")
    logger.info(f"sample {sample.name} uidx={uidx} cidx={cidx}")


def iter_clusters(clusters: Path) -> Iterator[List[Tuple[str, str]]]:
    """Yields clusters between //\n// separators to feed to aligner."""
    xopen = gzip.open if clusters.suffix == ".gz" else open
    with xopen(clusters, 'rt') as clustio:
        data = []
        pairs = zip(clustio, clustio)
        for line1, line2 in pairs:
            # line1, line2 = line1.decode(), line2.decode()
            if line1[0] == ">":
                data.append((line1, line2))
            else:
                yield data
                data = []


def revise_mixed_cluster(clust: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Return cluster with only merged or non-merged readpairs, but not mixed.

    Which ever type is more common is returned, and the other type is
    discarded. This type of cluster can be common when the read pairs
    overlap only by a few base pairs at their ends, and read merging by
    vsearch is imperfect.
    """
    size_w_n = 0
    size_wo_n = 0
    for seq in clust:
        # get nreps if within-cluster, else 1 if across-cluster
        try:
            reps = int(seq[0].split(";")[1].split("=")[-1])
        except ValueError:
            reps = 1
        if "nnnn" in seq[1]:
            size_w_n += reps
        else:
            size_wo_n += reps

    reclust = []
    for seq in clust:
        header = seq[0]
        if size_wo_n >= size_w_n:
            if "nnnn" not in seq[1]:
                reclust.append((header, seq[1]))
        else:
            if "nnnn" in seq[1]:
                reclust.append((header, seq[1]))
    clust = reclust

    # in case seed was removed, relabel with *, +, -
    header = clust[0][0].strip()
    if header[-1] != "*":
        clust[0] = (header[:-1] + "*\n", clust[0][1])
        for lidx, seq in enumerate(clust[1:]):
            if ">" in seq[0]:
                if header[-1] == "-":
                    name = seq[0].strip()
                    if name[-1] == "-":
                        clust[lidx] = (name[:-1] + "+\n", seq[1])
                    else:
                        clust[lidx] = (name[:-1] + "-\n", seq[1])
    return clust


def iter_muscle_alignments(handle: Path) -> Iterator[Tuple[List[str], List[str]]]:
    """Generator of clusters aligned by muscle.

    Opens a subprocess running muscle and awaiting input on stdin which
    is then processed and read from stdout as a result. This cuts down
    on the overhead of opening and closing many subprocesses.
    """
    # muscle command to return alignment and then a spacer
    # muscle -quiet  -threads 1 -align /dev/fd/63 -output /dev/stdout
    cmd = [
        BIN_MUSCLE, "-quiet", "-threads", "1",
        "-align", "<(printf '{}')",
        "-output", "/dev/stdout",
        ";", "echo", "@@\n",
    ]
    cmd = " ".join(cmd)

    # Popen kwargs to pass a sequence alignment in
    kwargs = dict(
        stdout=subprocess.PIPE, stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT, bufsize=0,
    )

    # open a subprocess for each paired side
    with subprocess.Popen(
        'bash', **kwargs) as proc1, subprocess.Popen(
            'bash', **kwargs) as proc2:

        # iterate over clusters
        for clust in iter_clusters(handle):

            # choose block below based on frequency of 'nnnn' insert
            with_insert = sum('nnnn' in i[1] for i in clust)
            # logger.debug([with_insert, len(clust)])

            ###############################################################
            # pairs separator is variable in presence. If so, choose
            # the larger set: either the ones with 'nnnn' or without.
            if 0 < with_insert < len(clust):

                clust = revise_mixed_cluster(clust)
                # if reduced to a single uniq read then don't align it.
                # else, it will align in one of the two modes below.
                if len(clust) == 1:
                    ali = [clust[0][0], f"{SPACER}{clust[0][1].strip()}{SPACER}\n"]
                    yield ali, []
                    continue

            ###############################################################
            # NOT PAIRED (no pairs present, do single alignment)
            if not any("nnnn" in i[1] for i in clust):
                # insert a spacer around the sequence
                for idx, seq in enumerate(clust):
                    clust[idx] = (seq[0], f"{SPACER}{seq[1].strip()}{SPACER}\n")

                # build the cmd and send as bytes to running subprocess on stdin
                sequence = "".join("".join(i) for i in clust)
                # logger.warning(f"AAA\n{sequence}")
                cmd_with_seq = cmd.format(sequence)
                proc1.stdin.write(cmd_with_seq.encode())

                # stream bytes output from running subprocess stdout
                ali = [i.decode() for i in iter(proc1.stdout.readline, b'@@\n')]
                if ">" not in ali[0]:
                    raise IPyradError(
                        f"error in muscle alignment: {''.join(ali)}")

                # yield aligned cluster (R1 only)
                yield ali, []
                continue

            ##############################################################
            # PAIRED: pairs present. Split paired data on 'nnnn' separator.
            # Attaches 'edge blocks' "NNNNNN" to improve alignment at
            # edges of clusters. These don't need to be removed since
            # they will be trimmed in later step. Note: this handles the
            # case of an empty alignment half (e.g., nnnnATATAT...) by
            # providing the padding Ns still as input for the missing half.
            read1s = []
            read2s = []
            for seq in clust:
                pos = seq[1].find("nnnn")
                read1s.append((seq[0], f"{SPACER}{seq[1].strip()[:pos]}{SPACER}\n"))
                read2s.append((seq[0], f"{SPACER}{seq[1].strip()[pos + 4:]}{SPACER}\n"))

            # align reads simultaneously on separate subprocesses
            sequence = "".join("".join(i) for i in read1s)
            # logger.warning(f"BBB1\n{sequence}")
            cmd_with_seq = cmd.format(sequence)
            proc1.stdin.write(cmd_with_seq.encode())
            sequence = "".join("".join(i) for i in read2s)
            # logger.warning(f"BBB2\n{sequence}")
            cmd_with_seq = cmd.format(sequence)
            proc2.stdin.write(cmd_with_seq.encode())

            # collect alignments, they will be re-paired in next func
            ali1 = [i.decode() for i in iter(proc1.stdout.readline, b'@@\n')]
            ali2 = [i.decode() for i in iter(proc2.stdout.readline, b'@@\n')]
            for ali in (ali1, ali2):
                if ">" not in ali[0]:
                    raise IPyradError(
                        f"error in muscle alignment: {''.join(ali)}")
            yield ali1, ali2


def get_cluster_size(header: str) -> int:
    return int(header.split(";")[-2][5:])


def iter_muscle_alignments_formatted(handle: Path) -> Iterator[str]:
    """Generator of aligned, formatted, and sorted clusters.

    Iterates over alignments from iter_muscle_alignments and rejoins
    paired halves and formats as a string for writing.
    """
    # iterate over aligned chunks
    for ali1, ali2 in iter_muscle_alignments(handle):
        # get dict mapping {headers: sequences w/o newline breaks}
        head_to_seq1 = {}
        for line in ali1:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq1:
                    head_to_seq1[key] += line.strip()
                else:
                    head_to_seq1[key] = line.strip()

        # get dict mapping {headers: sequences w/o newline breaks}
        head_to_seq2 = {}
        for line in ali2:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq2:
                    head_to_seq2[key] += line.strip()
                else:
                    head_to_seq2[key] = line.strip()

        # sort the first reads by size and seed
        try:
            keys = sorted(head_to_seq1, key=get_cluster_size, reverse=True)
            seed = [i for i in keys if i[-1] == "*"][0]
            seed = keys.pop(keys.index(seed))
            order = [seed] + keys
        except IndexError:
            # report to logger
            print(f"@@DEBUG: ERROR IN CLUST SORT:\n\n{keys}\n\n{ali1}\n\n{ali2}")
            continue

        # sort read2s by their names in the same order.
        alignment = []
        if head_to_seq2:
            for key in order:
                alignment.append(
                    f"{key}\n{head_to_seq1[key]}nnnn{head_to_seq2[key]}")
        else:
            for key in order:
                alignment.append(f"{key}\n{head_to_seq1[key]}")
        yield "\n".join(alignment)


def write_alignments(handle: Path) -> None:
    """Writes alignments to tmpfiles.

    Iterates over chunks of aligned loci from generator and writes to a
    an output file CHUNKSIZE alignments at a time.
    """
    tmpout = handle.parent / handle.name.replace("_unaligned_", "_aligned_")
    chunk = []
    with open(tmpout, 'w', encoding="utf-8") as out:
        for idx, alignment in enumerate(iter_muscle_alignments_formatted(handle)):
            chunk.append(alignment)
            # if not idx % 1_000:
            #     out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")
            #     chunk = []
        if chunk:
            out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")
    print(f"@@INFO: finished alignment: {tmpout}", flush=True)


def reconcat(data: Assembly, sample: Sample) -> Tuple[int, float, float, float]:
    """Concatenate _aligned_ tmp chunks into a single .clusters.gz file.
    """
    chunks = list(data.tmpdir.glob(f"{sample.name}_aligned_*.fa"))

    # concatenate finished clusters
    clustfile = data.stepdir / f"{sample.name}.clusters.gz"
    with gzip.open(clustfile, 'wt') as out:
        depths = Counter()
        for fname in chunks:
            chunkdat = []
            for clust in iter_clusters(fname):
                # store depths info
                depth = sum(get_cluster_size(item[0]) for item in clust)
                depths[depth] += 1

                # store for writing
                chunkdat.append("".join("".join(i) for i in clust))
            # write to file
            out.write("\n".join(chunkdat))
            chunkdat = []

        # calculate stats
        nclusters = sum(depths.values())
        vals = list(itertools.chain(*[[key] * value for (key, value) in depths.items()]))
        mean = np.mean(vals)
        median = np.median(vals)
        std = np.std(vals)
        print(
            f"@@DEBUG: {sample.name} nclusters={nclusters}; "
            f"cluster_depths=(mean: {mean:.2f}, median: {median:.2f}, std: {std:.2f})"
        )
    return (nclusters, mean, median, std)


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")
    from ipyrad.core.load_json import load_json

    # data = load_json("/tmp/pairgbs_merge.json")
    # sample = data.samples["1A_0"]
    # print(sample.files.clustmap)

    data = load_json("../../pedtest/test-c85.json")
    sample = data.samples["kansuensis-DE366"]
    print(sample.files)
    # data.tmpdir = Path("/tmp")
    write_tmp_clusters(data.tmpdir, sample)
    # # for ufile in sorted(unaligned_files):
    # write_alignments(unaligned_files[2])
    # aligned_files = Path("/tmp/").glob("kansuensis-DE366_aligned_*.fa")
    # fnames = list(aligned_files)
    # reconcat(fnames)


    # for j, k in iter_muscle_alignments(ii):
    #     # pass
    #     try:
    #         print(f"{''.join(j)}")
    #     except TypeError:
    #         print("ERROR\n", j)
    # print(''.join(a))
    # print(''.join(b))
    # iclust = iter_clusters(
    # clust = next(iclust)


    # clust = "".join("".join(i) for i in clust)

    # cmd = [
    #     BIN_MUSCLE, "-quiet", "-threads", "1",
    #     "-align", "<(printf '{}')",
    #     "-output", "/dev/stdout",
    #     ";", "echo", "@@\n",
    # ]
    # cmd = " ".join(cmd)
    # byte_cmd = cmd.format(clust)
    # print(byte_cmd)


    # for cidx, clust in enumerate(iter_build_clusters(sample)):
    #     pass
    #     print(clust)
    #     print(f"CLUST={cidx}\n")
