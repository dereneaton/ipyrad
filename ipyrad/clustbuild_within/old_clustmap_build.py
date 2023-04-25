#!/usr/bin/env python

"""

1. If clusters can be built with derep, tsv then do it in next step.
2. If not, build in s2, and make s3 _alignments_within

"""

from typing import Iterator, List, Tuple, Dict
import subprocess
from pathlib import Path
from collections import Counter
from ipyrad.core import BaseStep, Assembly
from ipyrad.schema import Sample
from ipyrad.core.utils import comp


SPACER = "N" * 10


class Step3(BaseStep):
    def __init__(self, data, force, ipyclient):
        super().__init__(data, force)
        self.ipyclient = ipyclient

    def run(self):
        self._build()                             # -> tmp/.clusters.txt
        self.declone_clusters()                   # -> tmp/.clusters_decloned.txt
        self.muscle_chunk()                       # -> tmp/.ali
        self.muscle_align_chunks()                # -> tmp/.alignment
        self.calculate_sample_stats()


def iter_build_clusters(data: Assembly, sample: Sample) -> Iterator[str]:
    """Builds fasta like .clust.gz outputs from .utemp and .htemp.

    The .htemp files contain the seeds, and the .utemp contains info
    on the hits to the seeds.
    """
    infiles = [
        data.tmpdir / f"{sample.name}_derep.fa",
        data.tmpdir / f"{sample.name}_derep_tag.fa",
    ]
    derepfile = [i for i in infiles if i.exists()][-1]

    # i/o vsearch files
    uhandle = data.tmpdir / f"{sample.name}_utemp.tsv"
    usorthandle = data.tmpdir / f"{sample.name}_utemp_sort.tsv"
    hhandle = data.tmpdir / f"{sample.name}_htemp.fa"

    # Sort uhandle so we can read through matches efficiently
    cmd = ["sort", "-k", "2", str(uhandle), "-o", str(usorthandle)]
    subprocess.run(cmd, check=True)

    # load ALL derep reads into a dictionary (this can be a few
    # GB of RAM) and is larger if names are larger and if there are
    # more unique sequences. We are grabbing two lines at a time.
    alldereps = {}
    with open(derepfile, 'r', encoding="utf-8") as derepio:
        dereps = zip(*[derepio] * 2)
        for namestr, seq in dereps:
            nnn, sss = [i.strip() for i in (namestr, seq)]
            alldereps[nnn[1:]] = sss

    # func to sort clusters by size (used below)
    sorter = lambda x: int(x.split(";size=")[-1].split("\n")[0][:-2])

    # store observed seeds (this could count up to >million in bad data sets)
    seedsseen = set()

    # Iterate through the usort file grabbing matches to build clusters
    with open(usorthandle, 'r', encoding="utf-8") as usortio:
        lastseed = None    # keep track of last seed
        fseqs = []         # store a single cluster

        # iterate over all lines in the usort file.
        for line in usortio:
            hit, seed, _, _, ori, _ = line.strip().split()

            # same seed, append match
            if seed != lastseed:
                seedsseen.add(seed)

                # store the last cluster (fseq), count it, and clear fseq
                if fseqs:
                    # sort to seed then hits by size
                    fsort = sorted(fseqs[1:], key=sorter, reverse=True)
                    fseqs = [fseqs[0]] + fsort
                    yield "\n".join(fseqs)
                    fseqs = []

                # store the new seed on top of fseq list
                fseqs.append(f">{seed};*\n{alldereps[seed]}")
                lastseed = seed

            # add match to the seed
            # revcomp if orientation is reversed (comp preserves nnnn)
            seq = alldereps[hit] if ori == "+" else comp(alldereps[hit])[::-1]

            # only save if not too many indels. Here indels does not mean
            # the number of dashes, but rather, the number of gap-openings.
            # This is hard-coded here as a filter for bad alignments.
            nindels = len([i for i in seq.split("-") if i])
            if nindels <= data.max_indels:
                fseqs.append(f">{hit};{ori}\n{seq}")

    # write whatever is left over to the clusts file
    if fseqs:
        yield "\n".join(fseqs)

    # now write the seeds that had no hits from the data in htemp,
    # but do not write seeds if they appeared in usort above.
    with open(hhandle, 'r', encoding="utf-8") as iotemp:
        nohits = zip(*[iter(iotemp)] * 2)
        for line in nohits:
            line = line[0].strip()
            if line[1:] not in seedsseen:
                yield f"{line};*\n{alldereps[line[1:]]}"

    # delete large dict
    del alldereps


def write_clusters(data: Assembly, sample: Sample) -> None:
    """Write built clusters to .clust.txt file in chunks.

    Gets clusters from iter_build_clusters generator func.
    """
    clusthandle = data.tmpdir / f"{sample.name}_clust.txt"
    with open(clusthandle, 'w', encoding="utf-8") as clustio:
        clusters = []
        for idx, clust in enumerate(iter_build_clusters(data, sample)):
            clusters.append(clust)

            # occasionally write/dump stored clusters to file and clear mem
            if not idx % 10_000:
                clustio.write("\n//\n//\n".join(clusters) + "\n//\n//\n")
                clusters = []
        if clusters:
            clustio.write("\n//\n//\n".join(clusters) + "\n//\n//\n")


def muscle_chunker(data: Assembly, sample: Sample) -> None:
    """Splits .clust files into chunks for Muscle alignment.

    Each file is split into at most 10 chunks. Each chunk will be
    run on a separate core, with the largest cluster files starting
    first. Because the largest clusters are at the beginning of the
    clusters file, the chunks contain varying number of clusters to
    try to equalize their speed.
    """
    # only chunk up denovo data, refdata has its own chunking method which
    # makes equal size chunks, instead of uneven chunks like in denovo
    clustfiles = [
        data.tmpdir / (sample.name + "_clust.txt"),
        data.tmpdir / (sample.name + "_clust_decloned.txt"),
    ]
    clustfile = [i for i in clustfiles if i.exists()][-1]

    # get the number of clusters and chunk into 20 pieces
    nloci = sum(1 for i in iter_clusters(clustfile))
    optim = (nloci // 20) + (nloci % 20)
    inc = optim // 10

    # cluster generator
    clustgen = iter_clusters(clustfile)

    # splitting loci so first file is smaller and last file is bigger
    for idx in range(10):

        # how big is this chunk?
        this_chunk_size = optim + (idx * inc)
        left = nloci - this_chunk_size

        # grab next chunks-worth of data, or all that's left.
        if idx == 9:
            grabchunk = list(islice(clustgen, int(1e9)))
        else:
            grabchunk = list(islice(clustgen, this_chunk_size))
            nloci = left

        # write the chunk to file
        if grabchunk:
            tmpfile = data.tmpdir / f"{sample.name}_chunk_{idx}.ali"
            with open(tmpfile, 'w', encoding="utf-8") as out:
                grabclust = ["".join(i) for i in grabchunk]
                out.write("//\n//\n".join(grabclust).strip() + "\n//\n//\n")


def iter_muscle_alignments(handle: Path, maxdepth: int=100) -> Iterator[List[str]]:
    """Align .ali tmp cluster files using mucle.

    Note
    ----
    There is still some trouble with interrupting this when running
    on ipp engines; it can leave muscle jobs running.

    Parameters
    ----------
    handle: Path
    maxdepth: int
        This the max number of unique sequences that will be aligned.
        By using this we COULD greatly improve aligning speed for some
        datasets. Currently we set it to 100 so it usually has no
        effect.
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

    # create two persistent bash shells
    with subprocess.Popen(
        'bash', **kwargs) as proc1, subprocess.Popen(
            'bash', **kwargs) as proc2:

        # iterate over clusters
        for clust in iter_clusters(handle):

            # skip aligning if it is a singleton, but format like below
            if len(clust) == 2:
                yield clust, []
                continue

            # limit to first N unique sequences (pre-sorted by depth)
            # (this MAY be a bit limiting if decloning pcr duplicates.)
            # Here I considered that we could combine sequences with the
            # same sequence (already done except for when unique i5 umi
            # barcodes are present), but it would mess up the i5 tags
            # unless we combine multiple tags into the same sequence header
            # which could be done but hasn't been tried yet. Would require
            # updates on the parser end later too.
            clust = clust[:maxdepth]

            # ---------------------------------------------------------
            # pair separator is not consistently found. Single align it.
            # NOTE: this often leads to problems later, with poor alignment
            # in the middle and sometimes no data on one side of the nnnn.
            # Therefore this option was deprecated for the option below where
            # any merged sequence (w/o nnnn) are removed and only pairs kept.
            # ---------------------------------------------------------
            # if not all('nnnn' in i for i in clust[1::2]):
            #     proc1.stdin.write(cmd.format("\n".join(clust)).encode())
            #     ali = [i.decode() for i in iter(proc1.stdout.readline, b'@@\n')]
            #     if ">" not in ali[0]:
            #         raise IPyradError(
            #             f"error in muscle alignment: {''.join(ali)}")
            #     yield ali, []


            # pairs present SOME of the time. Keep only one set, either the ones
            # with 'nnnn' or without, depending on which set is larger.
            if not all('nnnn' in i for i in clust[1::2]):
                size_w_n = 0
                size_wo_n = 0
                for line in clust:
                    # get nreps if within-cluster, else 1 if across-cluster
                    if ">" in line:
                        try:
                            reps = int(line.split(";")[1].split("=")[-1])
                        except ValueError:
                            reps = 1
                    else:
                        if "nnnn" in line:
                            size_w_n += reps
                        else:
                            size_wo_n += reps

                reclust = []
                for line in clust:
                    if ">" in line:
                        header = line
                    else:
                        if size_wo_n >= size_w_n:
                            if "nnnn" not in line:
                                reclust.append(header)
                                reclust.append(line)
                        else:
                            if "nnnn" in line:
                                reclust.append(header)
                                reclust.append(line)
                clust = reclust

                # in case seed was removed, relabel with *, +, -
                header = clust[0].strip()
                if header[-1] != "*":
                    clust[0] = header[:-1] + "*\n"
                    for lidx, line in enumerate(clust):
                        if lidx >= 2:
                            if ">" in line:
                                if header[-1] == "-":
                                    if line[-1] == "-":
                                        clust[lidx] = line[:-1] + "+\n"
                                    else:
                                        clust[lidx] = line[:-1] + "-\n"

                # if reduced to a single uniq read then don't align it.
                if len(clust) == 2:
                    yield clust, []
                    continue

            # no pairs present, do single alignment
            if not any('nnnn' in i for i in clust[1::2]):
                for idx, line in enumerate(clust):
                    if line[0] != ">":
                        clust[idx] = SPACER + line.strip() + SPACER
                proc1.stdin.write(cmd.format("\n".join(clust)).encode())
                ali = [i.decode() for i in iter(proc1.stdout.readline, b'@@\n')]
                if ">" not in ali[0]:
                    raise IPyradError(
                        f"error in muscle alignment: {''.join(ali)}")
                yield ali, []
                continue

            # pairs present. Split paired data on 'nnnn' separator.
            # Attaches 'edge blocks' "NNNNNN" to improve alignment at
            # edges of clusters. These don't need to be removed since
            # they will be trimmed in step 5.
            headers = []
            read1s = []
            read2s = []
            for line in clust:
                if line[0] == ">":
                    headers.append(line)
                else:
                    pos = line.find("nnnn")
                    read1s.append(SPACER + line[:pos] + SPACER)
                    read2s.append(SPACER + line[pos + 4:].strip() + SPACER)
            clust1 = ["\n".join(i) for i in zip(headers, read1s)]
            clust2 = ["\n".join(i) for i in zip(headers, read2s)]

            # align reads simultaneously on separate subprocesses
            proc1.stdin.write(cmd.format("\n".join(clust1)).encode())
            proc2.stdin.write(cmd.format("\n".join(clust2)).encode())

            # collect alignments, they will be re-paired in next func
            ali1 = [i.decode() for i in iter(proc1.stdout.readline, b'@@\n')]
            ali2 = [i.decode() for i in iter(proc2.stdout.readline, b'@@\n')]
            for ali in (ali1, ali2):
                if ">" not in ali[0]:
                    raise IPyradError(
                        f"error in muscle alignment: {''.join(ali)}")
            yield ali1, ali2


def iter_alignment_format(handle: Path) -> Iterator[str]:
    """Generator to yield filtered, aligned, paired, and sorted clusters.

    Example for testing
    -------------------
    >>> tmpfile = "/tmp/test_tmp_clustmap/1A_0_chunk_0.ali"
    >>> align_gen = iter_alignment_format(tmpfile)
    >>> print(next(align_gen))
    """
    # lambda func to sort hits by size
    sorter = lambda x: int(x.split(";")[-2][5:])

    # iterate over aligned chunks
    for ali1, ali2 in iter_muscle_alignments(handle):
        # get dict mapping headers to sequences
        head_to_seq1 = {}
        for line in ali1:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq1:
                    head_to_seq1[key] += line.strip()
                else:
                    head_to_seq1[key] = line.strip()
        head_to_seq2 = {}
        for line in ali2:
            if line[0] == ">":
                key = line.strip()
            else:
                if key in head_to_seq2:
                    head_to_seq2[key] += line.strip()
                else:
                    head_to_seq2[key] = line.strip()

        # sort the first reads
        try:
            keys = sorted(head_to_seq1, key=sorter, reverse=True)
            seed = [i for i in keys if i[-1] == "*"][0]
            seed = keys.pop(keys.index(seed))
            order = [seed] + keys
        except IndexError:
            print(f"TRY:\n\n{keys}\n\n{ali1}\n\n{ali2}")
            continue

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

    Iterates over chunks of aligned loci from generator and
    writes to a concatenated output file.
    """
    print(f"aligning: {handle}") # engine sends to logger.info
    tmpout = handle.with_suffix(".alignment")
    chunk = []
    with open(tmpout, 'w', encoding="utf-8") as out:
        for idx, alignment in enumerate(iter_alignment_format(handle)):
            chunk.append(alignment)
            if not idx % 1_000:
                out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")
                chunk = []
        if chunk:
            out.write("\n//\n//\n".join(chunk) + "\n//\n//\n")


def reconcat(data: Assembly, sample: Sample) -> None:
    """Concatenate .aligned chunks into a single .clusters.gz file.
    """
    chunks = list(data.tmpdir.glob(f"{sample.name}_chunk_[0-9].alignment"))
    chunks.sort(key=lambda x: int(x.name.rsplit("_", 1)[-1][:-10]))

    # concatenate finished clusters
    clustfile = data.stepdir / f"{sample.name}.clusters.gz"
    with gzip.open(clustfile, 'w') as out:
        for fname in chunks:
            with open(fname, 'r', encoding="utf-8") as infile:
                dat = infile.read().strip()
                if dat:
                    out.write(dat.encode() + b"\n")


# def iter_seeded_clusters(derep) -> Iterator[Tuple[str, str]]:



def fill_match_dicts(matches) -> Tuple[Dict, Dict]:
    """Return dictionaries ...

    TODO: could filter on max_gaps here prior to alignment.
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


def iter_derep(derep) -> Iterator[Tuple[str, str]]:
    with open(derep, 'r', encoding="utf-8") as infile:
        for header, sequence in zip(infile, infile):
            yield header[1:].strip(), sequence.strip()


def iter_built_clusters(sample) -> Iterator[str]:
    """

    TODO: add '>' to it.
    """
    clusters = {}
    hits_to_seeds, seeds_to_nhits = fill_match_dicts(sample.files.clustmap[1])

    # iterate over all sequences in the derep file
    for name, seq in iter_derep(sample.files.clustmap[0]):

        # it is a match
        if name in hits_to_seeds:
            seed, orient = hits_to_seeds[name]
            if orient == "-":
                seq = comp(seq)[::-1]
                newname = name + ";-"
            else:
                newname = name + ";+"
            if seed not in clusters:
                clusters[seed] = [newname, seq]
            else:
                clusters[seed].extend([newname, seq])

        # it is a seed
        else:
            seed = name
            newname = name + ";*"
            if seed not in clusters:
                clusters[seed] = [newname, seq]
            else:
                clusters[seed] = [newname, seq] + clusters[seed]

        # singletone (just seed, no hits)
        if seed not in seeds_to_nhits:
            yield "\n".join(clusters[seed])
            del clusters[seed]

        # yield cluster if it is complete
        else:
            if len(clusters[seed]) == seeds_to_nhits[seed] * 2:
                yield "\n".join(clusters[seed])
                del clusters[seed]

    # everything should hvae been yielded.
    assert not clusters



if __name__ == "__main__":

    from ipyrad.core2.load_json import load_json
    data = load_json("/tmp/pairgbs_merge.json")
    sample = data.samples["1A_0"]
    print(sample.files.clustmap)

    for cidx, clust in enumerate(fill_2(sample)):
        pass
        print(clust)
        print(cidx, "\n")
