#!/usr/bin/env python

"""Group files into pairs based on file name matching patterns.

The `get_filenames_to_paired_fastqs` function is used in both demux
and step1 of assembly to return full Paths to fastq files as tuples
of either (R1, None) or (R1, R2).
"""

from typing import List, Tuple, Dict, Union, Sequence
from pathlib import Path
from collections import Counter
import itertools

__all__ = [
    "get_paths_list_from_fastq_str",
    "get_single_clean_names_dict",
    "get_paired_clean_names_dict",
]


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
            raise ValueError(f"No fastq data match input: {path}")
        expanded.extend([Path(i).expanduser().resolve() for i in fastqs])
    return expanded


def get_shortest_shared_prefix(name1: str, name2: str) -> str:
    """Return the shortest name prefix shared by two strings.
    """
    # get first char where names do not match
    for num, (i, j) in enumerate(zip(name1, name2)):
        if i != j:
            break
    name = name1[:num]

    # strip ugly endings
    suffs = [".", "_", "_R"]
    while any(name.endswith(i) for i in suffs):
        for suff in suffs:
            name = name.rstrip(suff)
    return name


def get_stripped_length(path: Path) -> str:
    """Return the length of single suffix stripped name.
    """
    # strip ugly endings
    name = str(path.name)
    suffs = [".", "_", "_R1", ".gz", ".fq", ".fastq", "_1"]
    while any(name.endswith(i) for i in suffs):
        for suff in suffs:
            # rstrip(suff) treats suff as a list of characters to remove
            # and not as a suffix, so we need to use removesuffix(suff)
            name = name.removesuffix(suff)
    return len(path.name) - len(name)


def get_single_clean_names_dict(fastqs: List[Path]) -> Dict[str, Tuple[Path, Path]]:
    """Return ordered sample names from a list of path names.

    We strip off suffix endings from names. Some names may actually
    have these in their names, such as '_1', so we try to the most
    common suffix length to strip from all names that is composed of
    only suffix-like characters.
    """
    snames_to_fastq_tuples = {}
    stripped_lens = [get_stripped_length(i) for i in fastqs]
    # choose smaller trim len if a tie
    lens = Counter(stripped_lens)
    maxo = lens.most_common()[0][1]
    strip_len = min([i[0] for i in lens.most_common() if i[1] == maxo])
    # strip same length from all names
    names = [i.name[:-strip_len] for i in fastqs]
    return {i: (j, '') for i, j in zip(names, fastqs)}


def _sub_num_for_paired(path: Path, idx: int) -> str:
    """Substitute a 1 into a position in a path name.
    
    This is used to find unique pairs of path names differing by 1/2.
    """
    name = list(str(path))  # list(path.name)
    try:
        name[-idx] = "1"
    except IndexError:
        pass
    return "".join(name)


def get_paired_clean_names_dict(fastqs: List[Path]) -> Dict[str, Tuple[Path, Path]]:
    """Return dict with {name: (path1, path2)} based on file name matching.

    """
    snames_to_fastq_tuples = {}
    paths = sorted(fastqs)
    idx = 1

    # loop subbing '1' into name positions to find pairs differing by 1/2
    # all hits should happen about 5-15 positions back from end.
    while len(snames_to_fastq_tuples) < len(fastqs) / 2:
        groups = itertools.groupby(sorted(paths), lambda x: _sub_num_for_paired(x, idx))
        for key, group in groups:
            group = list(group)

            # a unique pair was found
            if key and len(group) == 2:
                r1, r2 = group
                name = get_shortest_shared_prefix(r1.name, r2.name)
                snames_to_fastq_tuples[name] = (r1, r2)
                paths.remove(r1)
                paths.remove(r2)
        idx += 1

        # give up finding pairs on large idx
        if idx > 100:
            return snames_to_fastq_tuples

    # check if pairs were detected and if so, break
    # print(len(snames_to_fastq_tuples), len(fastqs))
    return snames_to_fastq_tuples


def get_fastq_tuples_dict_from_paths_list(fastqs: List[Path], paired: bool) -> Dict[str, Tuple[Path, Path]]:
    """Return dict {clean_name: (Path, Path)} from input fastq paths.

    """

    if paired:
        # get as dict {str: (Path, Path)}
        # pe data
        ndict = get_paired_clean_names_dict(fastqs)
        if len(ndict) != len(fastqs) / 2:
            raise ValueError("Can't recognize paired fastq files. Check the format of fastq path names")
        return {i: ndict[i] for i in sorted(ndict)}
    else:
        # se data
        ndict = get_single_clean_names_dict(fastqs)
        return {i: ndict[i] for i in sorted(ndict)}


if __name__ == "__main__":

    # import ipyrad as ip

    # # se test
    # path = Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/*_R1*")
    # ndict = get_fastq_tuples_dict_from_paths_list(path)
    # for i in ndict:
    #     print(i, ndict[i])

    # # pe test
    # path = Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/*_R*")
    # ndict = get_fastq_tuples_dict_from_paths_list(path)
    # for i in ndict:
    #     print(i, ndict[i])

    # pretend files with '.' in names
    paths = [
        Path("/mnt/dir/genus.species_1.1.fq.gz"),
        Path("/mnt/dir/genus.species_1.2.fq.gz"),       
        Path("/mnt/dir/genus.species_2.1.fq.gz"),
        Path("/mnt/dir/genus.species_2.2.fq.gz"),
        Path("/mnt/dir/genus_species_sub2.1.fastq"),
        Path("/mnt/dir/genus_species_sub2.2.fastq"),
        Path("/mnt/dir/genus_species_x._R1_.fastq"),
        Path("/mnt/dir/genus_species_x._R2_.fastq"),        
    ]
    print(get_paired_clean_names_dict(paths))

    # expect: {genus.species_1: (Path, ''), genus.species_2: (Path, '')}
    paths = [
        Path("/mnt/dir/genus.species_1.1.fastq.gz"),
        Path("/mnt/dir/genus.species_2.1.fastq.gz"),
    ]
    print(get_single_clean_names_dict(paths))    
