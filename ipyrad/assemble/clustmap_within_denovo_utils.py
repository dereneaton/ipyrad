#!/usr/bin/env python

"""Functions run on remote parallel cluster in clustmap.

"""

from typing import TypeVar, Iterator, List, Dict, Tuple
import sys
import gzip
from pathlib import Path
import subprocess
from subprocess import Popen, PIPE, STDOUT
from itertools import islice

from loguru import logger
from ipyrad.assemble.utils import IPyradError, comp
from ipyrad.assemble.clustmap_within_both import iter_clusters

# pylint: disable=too-many-branches, too-many-statements, consider-using-with, too-many-lines

logger = logger.bind(name="ipyrad")

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")

BIN = Path(sys.prefix) / "bin"
BIN_MUSCLE = str(BIN / "muscle")
BIN_VSEARCH = str(BIN / "vsearch")
BIN_BWA = str(BIN / "bwa")
BIN_SAMTOOLS = str(BIN / "samtools")
SPACER = "N" * 10

################################################
####
#### Prepare data for decloning and clustering.
####
################################################

def merge_pairs_with_vsearch(data: Assembly, sample: Sample) -> int:
    """Merge PE reads using vsearch to find overlap."""
    # input files (select only the top one)
    in1 = [
        data.tmpdir / f"{sample.name}_unmapped_R1.fastq",
        data.tmpdir / f"{sample.name}_concat_R1.fq.gz",
        sample.files.edits[0][0],
    ]
    in2 = [
        data.tmpdir / f"{sample.name}_unmapped_R2.fastq",
        data.tmpdir / f"{sample.name}_concat_R2.fq.gz",
        sample.files.edits[0][1],
    ]
    index = min([i for i, j in enumerate(in1) if Path(j).exists()])
    infile1 = in1[index]
    infile2 = in2[index]

    # define output files
    mergedfile = data.tmpdir / f"{sample.name}_merged.fa"
    nonmerged1 = data.tmpdir / f"{sample.name}_nonmerged_R1.fa"
    nonmerged2 = data.tmpdir / f"{sample.name}_nonmerged_R2.fa"

    # get the maxn and minlen values
    try:
        maxn = sum(data.params.max_low_qual_bases)
    except TypeError:
        maxn = data.params.max_low_qual_bases
    minlen = str(max(32, data.params.filter_min_trim_len))

    # vsearch merge can now take gzipped files (v.2.8)
    cmd = [
        BIN_VSEARCH,
        "--fastq_mergepairs", str(infile1),
        "--reverse", str(infile2),
        "--fastaout", str(mergedfile),
        "--fastaout_notmerged_fwd", str(nonmerged1),
        "--fastaout_notmerged_rev", str(nonmerged2),
        "--fasta_width", "0",
        "--fastq_minmergelen", str(minlen),
        "--fastq_maxns", str(maxn),
        "--fastq_minovlen", "10", # default value,
        "--fastq_maxdiffs", "4",  # <- fastp has done pe overlap correction
        # "--label_suffix", "_m1",
        "--fastq_qmax", "93",     # <- Set high to allow FASTQ+64
        "--threads", "2",
        "--fastq_allowmergestagger",
    ]
    print(" ".join(cmd)) # sends to logger.info
    with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
        res = proc.communicate()[0].decode()
        if proc.returncode:
            raise IPyradError(f"Error merge pairs:\n {cmd}\n{res}")

    # count number of lines / 2 from output
    with open(mergedfile, 'r', encoding="utf-8") as merged:
        nlines = sum(1 for line in merged)
    return int(nlines / 2)

def join_end_to_end(data: Assembly, sample: Sample) -> None:
    """Combine read pairs end to end for dereplication and/or clustering.
    """
    # output path (may already contain overlap merged pairs)
    joined = data.tmpdir / f"{sample.name}_merged.fa"

    # input paths (read pairs not overlap merged)
    nonmerged1 = data.tmpdir / f"{sample.name}_nonmerged_R1.fa"
    nonmerged2 = data.tmpdir / f"{sample.name}_nonmerged_R2.fa"

    # iterate over pairs and write append to out
    tmp = []
    with open(joined, 'a', encoding="utf-8") as out:
        # open both indata io
        with open(
            nonmerged1, 'r', encoding="utf-8") as in1, open(
            nonmerged2, 'r', encoding="utf-8") as in2:

            # generator to sample 4 lines
            quart1 = zip(*[in1] * 4)
            quart2 = zip(*[in2] * 4)

            # write out and fasta: >header\n{seq1}nnnn{seq2}
            for idx, (read1s, read2s) in enumerate(zip(quart1, quart2)):
                header = read1s[0][1:]
                read1 = read1s[1].strip()
                read2 = comp(read2s[1].strip())[::-1]
                newread = f">{header}{read1}nnnn{read2}"
                tmp.append(newread)

                # write on intervals
                if not idx % 50_000:
                    out.write("\n".join(tmp) + "\n")
                    tmp = []

        # write remaining to outfile
        if tmp:
            out.write("\n".join(tmp))

        # remove temporary files
        nonmerged1.unlink()
        nonmerged2.unlink()

################################################
####
#### Cluster reads and build clusters
####
################################################

def cluster(data: Assembly, sample: Sample) -> None:
    """Calls vsearch for clustering with cluster_smallmem.

    cov varies by data type, values were chosen
    based on experience, but could be edited by users.
    """
    # get dereplicated reads for denovo+reference or denovo-reference
    handles = [
        data.tmpdir / f"{sample.name}_derep.fa",
        data.tmpdir / f"{sample.name}_derep_tag.fa",
    ]
    derephandle = [i for i in handles if i.exists()][-1]
    assert derephandle, "bad derep handle"

    # create handles for the outfiles
    uhandle = data.tmpdir / f"{sample.name}_utemp.tsv"
    temphandle = data.tmpdir / f"{sample.name}_htemp.fa"

    # datatype specific optimization
    # minsl: the percentage of the seed that must be matched
    #    smaller values for RAD/ddRAD where we might want to combine, say 50bp
    #    reads and 100bp reads in the same analysis.
    # query_cov: the percentage of the query sequence that must match seed
    #    smaller values are needed for gbs where only the tips might overlap
    #    larger values for pairgbs where they should overlap near completely
    #    small minsl and high query cov allows trimmed reads to match to untrim
    #    seed for rad/ddrad/pairddrad.
    strand = "plus"
    cov = 0.5
    minsl = 0.5
    if data.params.datatype in ["gbs", "2brad"]:
        strand = "both"
        cov = 0.5
        minsl = 0.5
    elif data.params.datatype == 'pairgbs':
        strand = "both"
        cov = 0.75
        minsl = 0.75

    # If this value is not null (which is the default) then override query cov
    if data.hackers.query_cov:
        cov = str(data.hackers.query_cov)
        assert float(cov) <= 1, "query_cov must be <= 1.0"

    # get call string
    # we used to do --cluster_smallmem w/ --usersort to sort by depth
    # but now there is an option --cluster_size to accomplish this,
    # however, --cluster_fast sorts by length, and this MAY be better
    # when most sequences are similar depths (e.g., 2-3).
    cmd = [
        BIN_VSEARCH,
        "--cluster_fast", str(derephandle), #
        # "--cluster_smallmem", str(derephandle),
        # "--usersort",  # paired w/ cluster_smallmem  sorted by depth
        "--strand", strand,
        "--query_cov", str(cov),
        "--id", str(data.params.clust_threshold),
        "--minsl", str(minsl),
        "--userout", str(uhandle),
        "--userfields", "query+target+id+gaps+qstrand+qcov",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        "--threads", str(max(1, data.ipcluster['threads'])),
        "--notmatched", str(temphandle),
        "--fasta_width", "0",
        "--minseqlength", str(data.params.filter_min_trim_len),
        # "-fastq_qmax", "100",
        # "-fulldp",  # no longer relevant, it is now always on.
    ]

    # run vsearch
    print(" ".join(cmd))  # engine sends to logger.info
    with subprocess.Popen(cmd,
        stderr=subprocess.STDOUT, stdout=subprocess.PIPE, close_fds=True,
        ) as proc:
        res = proc.communicate()[0]

    # check for errors
    if proc.returncode:
        raise IPyradError(f"cmd {cmd}: {res}")

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
            hit, seed, _, indels, ori, _ = line.strip().split()

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

            # only save if not too many indels
            if int(indels) <= data.max_indels:
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

################################################
####
#### Align clusters
####
################################################


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

            # pairs present SOME of the time. Remove unpaired reads.
            if not all('nnnn' in i for i in clust[1::2]):
                reclust = []
                for line in clust:
                    if ">" in line:
                        header = line
                    else:
                        if "nnnn" in line:
                            reclust.append(header)
                            reclust.append(line)
                clust = reclust

                # if reduced to a single uniq read then don't align it.
                if len(clust) == 2:
                    header = clust[0].strip()
                    if header[-1] != "*":
                        clust[0] = header[:-1] + "*\n"
                    yield clust, []
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


###################################################



#############################
##
## DEPRECATED
##
#############################


def alignment_indel_filter(clust: Dict[str,str], max_internal_indels: int) -> bool:
    """Return True if internal indels exceeds the allowed maximum."""
    return any(
        i.strip("-").count("-") > max_internal_indels
        for i in clust.values()
    )


def gbs_trim(align1):
    """Trimming method applied to gbs and pairgbs denovo assemblies only.

    No reads can go past the left of the seed, or right of the least extended
    reverse complement match. Example below. m is a match. u is an area where
    lots of mismatches typically occur. The cut sites are shown.

    Original locus*
    Seed           TGCAG************************************-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm-----------------------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm------------------------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuuuuuuuu
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCAuuuuuuuu

    Trimmed locus*
    Seed           TGCAG************************************---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm---------------
    Forward-match  TGCAGmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm----------
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    Revcomp-match  ---------------mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmCTGCA
    Revcomp-match  --------------------------------mmmmmmmmmmmmmmmmmm
    Revcomp-match  ------------------------mmmmmmmmmmmmmmmmmmmmmmmmmm
    """
    # get the left and rigth most occurrences of the cut site
    leftmost = rightmost = None

    # create dict mapping {seqname: seq}
    ddd = {k: v for k, v in [j.rsplit("\n", 1) for j in align1]}

    # get the seed sequence which contains "*" in header
    seed = [i for i in ddd.keys() if i.rsplit(";")[-1][0] == "*"][0]

    # position of the leftmost sequence that is not a "-" char.
    leftmost = [i != "-" for i in ddd[seed]].index(True)

    # which sequences are revcomp matched to the seed
    revs = [i for i in ddd.keys() if i.rsplit(";")[-1][0] == "-"]

    # ...
    if revs:
        subright = max([
            [i != "-" for i in seq[::-1]].index(True)
            for seq in [ddd[i] for i in revs]
        ])
    else:
        subright = 0
    rightmost = len(ddd[seed]) - subright

    # if locus got clobbered then print place-holder NNN
    names, seqs = zip(*[i.rsplit("\n", 1) for i in align1])
    if rightmost > leftmost:
        newalign1 = [n + "\n" + i[leftmost:rightmost]
                     for i, n in zip(seqs, names)]
    else:
        newalign1 = [n + "\nNNN" for i, n in zip(seqs, names)]
    return newalign1


if __name__ == "__main__":
    pass
