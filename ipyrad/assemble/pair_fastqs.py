#!/usr/bin/env python

"""Group files into pairs based on file name matching patterns.

The `get_filenames_to_paired_fastqs` function is used in both demux
and step1 of assembly to return full Paths to fastq files as tuples
of either (R1, None) or (R1, R2).
"""

from typing import List, Tuple, Dict, Union
from pathlib import Path
import itertools


def get_paths_list_from_fastq_str(fastq_paths: Union[Path, List[Path]]) -> List[Path]:
    """Expand fastq_paths str argument into a list of Path objects.

    """
    expanded = []
    # ensure paths is a List[Path] but where the Path elements may be
    # regex path names that have not yet been expanded.

    # ensure it is a list
    if isinstance(fastq_paths, (str, Path)):
        fastq_paths = [fastq_paths]

    # ensure each is a Path object
    paths = []
    for path in fastq_paths:
        if isinstance(path, str):
            path = Path(path)
        paths.append(path)

    # for each Path in paths list expand into a list of Paths
    for path in paths:
        # expand a regex operator to possibly match multiple files
        # such as paired-end files.
        try:
            fastqs = list(path.parent.glob(path.name))
            assert fastqs
        except (ValueError, AssertionError):
            msg = f"No fastq data match input: {path}"
            raise ValueError(msg)
        expanded.extend([Path(i).expanduser().resolve() for i in fastqs])
    return expanded


def drop_from_right(path: Path, delim: str = "_", idx: int = 0) -> str:
    """Return a name with an underscore separated portion removed.

    This is used to find matching paired files by iteratively removing
    subsections of filenames to find pairs that match when the sections
    are removed (e.g., usually _R1 _R2 or _1 _2).

    Sample names are returned _without_ trailing suffixes, to prevent
    having to reparse for the sample name in calling functions.

    Example
    -------
    >>> path = Path("name_prefix_001_R1_002.fastq.gz")
    >>> drop_from_right(path, "_", 0)
    >>> # "name_prefix_001_R1
    >>> drop_from_right(path, "_", 1)
    >>> # "name_prefix_001_002
    """
    # save and remove suffixes (it seems this method is needed with .x.y
    # suffixes compared to using Path.stem or similar.)
    suffixes = path.suffixes

    # Allow periods ('.') in sample names. `path.suffixes` splits on '.',
    # so a sample name that includes a period will get mangled (#557).
    # We make an assumption that the 'true' suffixes (e.g. .fastq.gz, or .fasta)
    # will not contain underscores, and that the first false suffix that
    # precedes these will include an '_' in the content of the mate pair
    # indicator ('_R1_' or '_1_' or '_R2_' for example. Allows this to be
    # a legal sample name: 1A_0.1_R2_.fastq.gz
    while any(["_" in x for x in suffixes]):
        suffixes.pop(0)

    while path.suffix in suffixes:
        path = path.with_suffix('')

    # break file name on delimiter and get chunks in reverse order
    chunks = path.name.split(delim)[::-1]

    # get chunks minus the index from the right
    sublist = [j for i, j in enumerate(chunks) if i != idx][::-1]
    path = path.parent / "_".join([i for i in sublist if i]).rstrip(delim)
    # Removed 7/17/24 iao
    # Related to #557, path.with_suffix will _overwrite_ anything it considers
    # a current suffix (anything with a '.'), so we need to protect against it.
    # suffixes = path.suffixes + suffixes
    # path = path.with_suffix("".join(suffixes))
    return path


def get_fastq_tuples_dict_from_paths_list(fastqs: List[Path]) -> Dict[str, Tuple[Path, Path]]:
    """Return {sample_name_str: tuple_of_paired_fastqs} dict.

    If technical replicates are being merged then a sample can
    be assigned multiple pairs of paired fastqs files. Paired
    file names may differ by having the following diffs, which
    may occur anywhere in the file name.
    >>> '_1', '_2'
    >>> '_R1', '_R2'
    >>> '_R1_', '_R2_'

    This func works by splitting file names by "_" and "." examining
    each section in turn, starting from the end, to find one that
    appears to match the paired naming conventions above.
    """
    # dict to return
    snames_to_fastq_tuples = {}

    # look for PE matching file names.
    idx = 0
    paired = False
    while 1:

        # group names with matching names when _ section is removed.
        try:
            groups = itertools.groupby(
                sorted(fastqs),
                key=lambda x: drop_from_right(x, "_", idx)
            )
            assert groups
            # sort the tuples to (R1, R2) and remove suffixes from names
            sorted_tuple_groups = {}
            for (name, gtup) in groups:
                gtup = sorted(gtup)
                path = Path(name)
                # Removed 7/7/24 iao #557
                #suffixes = path.suffixes
                #while path.suffix in suffixes:
                #    path = path.with_suffix('')
                sorted_tuple_groups[path.name] = gtup
            groups = sorted_tuple_groups

            # groups = {i.with_suffix("").name: sorted(j) for (i, j) in groups}
            assert len(groups) == len(fastqs) / 2
            assert all(len(j) == 2 for i, j in groups.items())
            paired = True
            break

        # if >5 tries and pairs not found then either data are not
        # paired, or the file names format is strange. We will raise
        # a warning later if the data seem to be PE but it was not
        # detected.
        except (AssertionError, ValueError):
            idx += 1
            if idx > 4:
                print(
                    "  No PE fastq pairs detected based on filenames, "
                    "assuming SE data."
                )
                break

    # store paired tuples to each basename
    if paired:
        for name, paths in groups.items():
            snames_to_fastq_tuples[name] = tuple(
                [i.expanduser().resolve() for i in paths]
            )

    # store (R1, '') tuples for each basename, and report a warning
    # if names are suspiciously paired looking.
    else:
        for path in fastqs:
            subpath = path.with_suffix("")
            while subpath.suffix:
                subpath = subpath.with_suffix("")
            snames_to_fastq_tuples[subpath.name] = (path.expanduser().resolve(), "")

            # warning if the data appear to include R2s
            if any(i in str(path.name) for i in ("_R2_", "_2.", "_R2.")):
                print(
                    f"  fastq file name ({path.name}) has a filename "
                    "that suggests it may be an R2 read, but its paired "
                    "R1 file could not be found. "
                    "Paired files should have matching names except "
                    "for _1 _2, _R1 _R2, or any of these followed by a "
                    "'_' or '.'."
                )
    return snames_to_fastq_tuples


if __name__ == "__main__":

    import ipyrad as ip

    FASTQ_PATH = Path("../../pedtest/small_tmp_R*.gz")
    FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)
    print(FASTQS)
    print(pairs)

    FASTQ_PATH = Path("../../pedtest/NEW_fastqs/*fastq.gz")
    FASTQS = get_paths_list_from_fastq_str(FASTQ_PATH)
    pairs = get_fastq_tuples_dict_from_paths_list(FASTQS)
    # print(pairs)