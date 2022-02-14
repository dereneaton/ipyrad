#!/usr/bin/env python

"""Some utilities used in demux.py for demultiplexing.
"""

from typing import Dict, Tuple, List, TypeVar, Iterator
import io
import glob
import gzip
import itertools
from pathlib import Path
import subprocess as sps
from dataclasses import dataclass, field

from loguru import logger
import pandas as pd
from ipyrad.assemble.utils import IPyradError, AMBIGS, BADCHARS
from ipyrad.core.progress_bar import AssemblyProgressBar


Assembly = TypeVar("Assembly")
Client = TypeVar("Client")
CHUNKSIZE = 20_000

logger = logger.bind(name="ipyrad")

# pylint: disable=consider-using-with, too-many-nested-blocks


@dataclass
class SimpleDemux:
    data: Assembly
    quiet: bool
    ipyclient: Client

    # to be filled.
    fastq_paths: List[Path] = None
    """: List of Paths to fastq files, unpaired."""
    barcodes_path: Path = None
    """: Path to the barcodes file."""
    names_to_barcodes: Dict[str, str] = None
    """: A map of barcode strings to sample names, pre-expanded by off-by-N."""
    filenames_to_fastqs: Dict[str, List[Tuple[str,str]]] = None
    """: Dict mapping file short names to tuples of paired fastqs."""
    cutters: List[Tuple[str,str]] = None
    """: List of enzyme overhang sites as two tuples with two strings each."""
    barcodes_to_names: Dict[str, str] = None
    """: Dict of all acceptable barcodes (e.g., off-by-1) mapped to sample names."""

    # longbar: int = None
    file_stats: Dict[str, List] = None

    def run(self):
        """Run subfunctions of demultiplexing workflow."""
        self._get_file_paths()
        self._get_filenames_to_fastqs()
        self._get_names_to_barcodes()
        self._replace_bad_name_chars()        
        self._get_cutters_expanded()
        self._get_barcodes_to_names()
        # self._distribute_barmatching_jobs()

    def _get_file_paths(self) -> None:
        """Get fastq and barcodes file paths as Path objects."""
        self.fastq_paths = glob.glob(str(self.data.params.raw_fastq_path))
        if not self.fastq_paths:
            raise IPyradError(
                f"No fastq data found in {self.data.params.raw_fastq_path}")
        self.fastq_paths = [Path(i) for i in self.fastq_paths]

        # if regular expression pick first barcode path.
        self.barcodes_path = glob.glob(str(self.data.params.barcodes_path))
        if not self.barcodes_path:
            raise IPyradError(
                f"No barcodes file found at {self.data.params.barcodes_path}")
        self.barcodes_path = Path(self.barcodes_path[0])

    def _get_names_to_barcodes(self) -> None:
        """Fill .names_to_barcodes dict w/ info from barcodes file.

        This logs a WARNING if technical replicates are detected to
        make sure the user is aware of how they are being handled.
        """
        # parse the tabular barcodes file on whitespace. Expects 
        # there to be no header. There will be >=2 columns, >2 if 
        # combinatorial barcodes. 
        bardata = pd.read_csv(
            self.barcodes_path, header=None, delim_whitespace=True,
            ).dropna()
        if bardata.shape[1] == 2:
            bardata.columns = ["sample", "barcode1"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
        else:
            bardata = bardata.iloc[:, :2]
            bardata.columns = ["sample", "barcode1", "barcode2"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
            bardata["barcode2"] = bardata["barcode2"].str.upper()

        # check for replicate sample names in the barcodes file. These
        # are allowed, since a single sample can be sequenced multiple
        # times on the same plate with different barcodes attached, 
        # representing technical replicates. THere is a hackers option
        # for whether to combine tech reps, or keep as diff samples.
        if bardata.value_counts().max() > 1:
            # get duplicated names
            duplicated = (bardata.value_counts() > 1).index

            # warn that dups are present AND WILL BE merged.
            if self.data.hackers.merge_technical_replicates:
                logger.warning(
                    "Technical replicates are present (samples with same name "
                    "in the barcodes file) and will be merged into one sample.\n"
                    "To not merge replicate samples set "
                    "`hackers.merge_technical_replicates = False), which will "
                    "instead rename the samples with a '-technical-replicate-x' "
                    "suffix.\nYou can change this in the 'hacker' settings "
                    "in the project JSON file.")
            # warn that dups are present and WILL NOT be merged.
            else:
                logger.warning(
                    "Technical replicates are present (samples with same name "
                    "in the barcodes file) and WILL NOT be merged.\n"
                    "To instead merge replicate samples set "
                    "`hackers.merge_technical_replicates = True), which will "
                    "combine reads from the same sample name into a single "
                    "sample.\nYou can change this in the 'hacker' settings "
                    "in the project JSON file.")

            # either way, relabel the samples for now, and may or may not merge later.
            for dup in duplicated:
                ridxs = bardata[bardata.sample == dup]
                if ridxs.shape[0] > 1:
                    for idx, index in enumerate(ridxs.index):
                        newname = f"{dup}-technical-replicate-{idx}"
                        bardata.loc[index, 'sample'] = newname

        # make sure barcodes are valid characters.
        if not all(bardata["barcode1"].apply(set("RKSYWMCATG").issuperset)):
            raise IPyradError(
                "Barcodes file contains unexpected characters in the "
                "barcode sequences suggesting it is not correctly "
                "formatted. See documentation.")

        # convert bardata to a dictionary {sample: barcode}.
        # if combinatorial barcodes are present then combine them.
        if "barcode2" in bardata.columns:
            assert self.data.is_pair, (
                "only paired datatypes can make use of combinatorial barcodes.")
            self.names_to_barcodes = dict(
                zip(bardata["sample"], bardata["barcode1"] + bardata["barcode2"]))
        else:
            self.names_to_barcodes = dict(
                zip(bardata["sample"], bardata["barcode1"]))

        # report to logger
        logger.debug(f"barcodes map:\n{bardata}")

    def _get_filenames_to_fastqs(self) -> None:
        """Fill `names_to_fastqs` with paired fastq files.

        If technical replicates are being merged then a sample can 
        be assigned multiple pairs of paired fastqs files. Paired
        file names may differ by having the following diffs, which
        may occur anywhere in the file name.
        >>> '_1', '_2'
        >>> '_R1', '_R2'
        >>> '_R1_', '_R2_'

        This func works by splitting file names by "_" and examining
        each section in turn, starting from the end, to find one that
        appears to match the paired naming conventions above.
        """
        self.filenames_to_fastqs = {}
        # if data are not paired then there is nothing to look for.
        if not self.data.is_pair:
            bad_paths = []
            endings = ("_R2_", "_2.", "_R2.")
            for path in self.fastq_paths:
                if any(path.suffix == i for i in endings):
                    bad_paths.append(str(path))
            if bad_paths:
                logger.warning(
                    "Fastq file names looks suspicously like PE data, even "
                    "though you selected a SE data type. Consider changing "
                    "the datatype param to a paired method like 'pairddrad'.\n"
                    f"Example filenames: {bad_paths}"
                )
            for path in self.fastq_paths:
                name = path.with_suffix("").name
                self.filenames_to_fastqs[name] = (path, "")
            return None

        # data type is PE data.
        idx = 0
        while 1:
            try:
                # group names with matching names when _ section is removed.
                groups = itertools.groupby(
                    self.fastq_paths,
                    key=lambda x: drop_from_right(x, "_", idx))
                assert groups
                groups = {i.with_suffix("").name: list(j) for (i, j) in groups}                
                assert len(groups) == len(self.fastq_paths) / 2
                assert all(len(j) == 2 for i, j in groups.items())
                logger.debug(f"found PE matches: {groups}")
                break
            except AssertionError as inst:
                # if >5 tries and pairs not found then raise an error.
                idx += 1
                if idx > 4:
                    raise IPyradError(
                        "Cannot parse paired file names. File names must have "
                        "matching name prefix followed by _1 _2, _R1 _R2, "
                        "or _R1_ _R2_ followed by any subsequent suffix. "
                        f"Your data files look like this: {self.fastq_paths}"
                        ) from inst

        # store as outputs
        for name, paths in groups.items():
            self.filenames_to_fastqs[name] = tuple(paths)
        logger.debug(self.filenames_to_fastqs)
        return None

    def _replace_bad_name_chars(self) -> None:
        """Replaces bad characters in names in .names_to_barcodes."""
        names = list(self.names_to_barcodes)
        for name in names:
            if any(i in name for i in BADCHARS):
                newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(
                    f"changing name {name} to {newname} (bad characters).")
                self.names_to_barcodes[newname] = self.names_to_barcodes.pop(name)

    def _get_cutters_expanded(self) -> None:
        """Fills self.cutters with both resolutions if IUPAC ambig present.
        
        - ('TGCAG', '') -> [('TGCAG', ''), ('', '')]
        - ('TWGC', '') -> [('TAGC', 'TTGC'), ('', '')]
        - ('TWGC', 'AATT') -> [('TAGC', 'TTGC'), ('AATT', '')]
        """
        self.cutters = []
        for cutter in self.data.params.restriction_overhang:
            if not any(i in 'RKSYWM' for i in cutter):
                self.cutters.append((cutter, ""))
            else:
                cutter1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cutter]
                cutter2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cutter]
                self.cutters.append((cutter1, cutter2))

    def _get_barcodes_to_names(self) -> None:
        """Fills .barcodes_to_names with all acceptable barcodes -> name."""
        self.barcodes_to_names = {}
        bases = set("ACGTN")

        # store perfect match to barcode
        for name, barcode in self.names_to_barcodes.items():
            self.barcodes_to_names[barcode] = name

        # store off-by-n matches
        for _ in range(self.data.params.max_barcode_mismatch):
            # get current pairs so we don't change dict while iter'ing
            current_bars = list(self.barcodes_to_names.items())
            for barcode, name in current_bars:
                for idx, base in enumerate(barcode):
                    # get bases that are not this one including N
                    diffs = bases.difference(base)
                    # iter over barcodes that can be made by substituting
                    for diff in diffs:
                        lbarcode = list(barcode)
                        lbarcode[idx] = diff
                        sbarcode = "".join(lbarcode)

                        match = self.barcodes_to_names.get(sbarcode, None)
                        if match is None:
                            self.barcodes_to_names[sbarcode] = name
                        else:
                            if match != name:
                                logger.warning(
                                f"\nSample: {name} ({barcode}) is within "
                                f"{self.data.params.max_barcode_mismatch} "
                                f"base changes of sample ({self.barcodes_to_names[barcode]}). "
                                "Ambiguous barcodes that match to both samples "
                                "will arbitrarily be assigned to the first "
                                "sample. If you do not like this then lower "
                                "the value of max_barcode_mismatch and rerun "
                                "(recommended).")
        logger.debug(f"barcodes_to_names: {self.barcodes_to_names}")

    def _distribute_barmatching_jobs(self) -> None:
        """Send barcode matching jobs to remote engines."""
        lbview = self.ipyclient.load_balanced_view()
        jobs = {}
        for fidx, fname in enumerate(self.filenames_to_fastqs):
            fastqs = self.filenames_to_fastqs[fname],
            args = (self.data, fastqs, self.barcodes_to_names, self.cutters, fidx)
            # jobs[fname] = lbview.apply(barmatch, args)
        msg = "demultiplexing reads"
        prog = AssemblyProgressBar(jobs, msg, step=1, quiet=self.quiet)
        prog.update()
        prog.block()
        prog.check()
        # jobs write outputs to tmpfiles.

    def _collect_stats(self):
        pass


@dataclass
class BarMatching:
    """Base class for barcode matching. 

    See subclasses which have different versions of the function
    `_iter_matched_barcode` to find barcode matches based on i7, 
    combinatorial, or single inline barcodes. The subclasses all 
    share the functions of this class, which includes iterating
    over the fastq(s), storing stats, and writing to tmp files.
    """
    data: Assembly
    """: Assembly object with param settings."""
    fastqs: Tuple[str, str]
    """: A tuple with paired R1 and R2 fastq files."""
    barcodes_to_names: Dict[str, str]
    """: Dict matching barcodes to sample names."""
    cutters: Tuple[str, str]
    """: List of Tuples of RE overhangs."""
    fidx: int
    """: File index."""

    # stats counters
    barcode_misses: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that don't match."""
    barcode_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that match."""
    sample_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record number of hits per sample."""

    def _iter_fastq_reads(self):
        """Yields fastq quartets of lines from fastqs (gzip OK)."""
        # create first read iterator for paired data    
        opener = gzip.open if self.fastqs[0].suffix == ".gz" else io.open
        ofile1 = opener(self.fastqs[0], 'rt', encoding="utf-8")
        quart1 = zip(ofile1, ofile1, ofile1, ofile1)

        # create second read iterator for paired data
        if self.fastqs[1]:
            ofile2 = opener(self.fastqs[1], 'rt', encoding="utf-8")
            quart2 = zip(ofile2, ofile2, ofile2, ofile2)
        else:
            quart2 = iter(int, 1)

        # yield from iterators as 4 items as a time (fastq)
        for read1, read2 in zip(quart1, quart2):
            yield read1, read2

    def _iter_matched_barcode(self):
        """SUBCLASSES REPLACE THIS FUNCTION."""
        raise NotImplementedError("See subclasses.")

    def _iter_matched_chunks(self):
        """Stores matched reads until N then writes to file."""
        read1s = {}
        read2s = {}
        nstored = 0

        # iterate over matched reads
        for read1, read2, match in self._iter_matched_barcode():

            # store r1 as 4-line string
            fastq1 = "".join(read1)
            if match in read1s:
                read1s[match].append(fastq1)
            else:
                read1s[match] = [fastq1]

            # store r2 as 4-line string
            if read2:
                fastq2 = "".join(read2)
                if match in read2s:
                    read2s[match].append(fastq2)
                else:
                    read2s[match] = [fastq2]

            # write to file when size is big enough and reset.
            nstored += 1
            if nstored > CHUNKSIZE:
                yield read1s, read2s
                read1s = {}
                read2s = {}
                nstored = 0

        # write final chunk if data
        yield read1s, read2s

    def run(self):
        """Iterate over all lines matching barcodes and recording stats, 
        and write the matched reads to unique files in chunks.

        Write chunks to tmp files for each sample w/ data.
        Opens a file handle that is unique to this process/sample
        """
        for read1s, read2s in self._iter_matched_chunks():
            for name in read1s:

                # write to R1 chunk file.                
                path1 = self.data.tmpdir / f"{name}_R1.tmp{self.fidx}.fastq"
                data = read1s[name]
                with open(path1, 'a', encoding="utf-8") as out:
                    out.write("\n".join(read1s))
                    logger.debug(f"wrote demuliplex chunks to {path1}")

                # write to R2 chunk file.
                if read2s:
                    path2 = self.data.tmpdir / f"{name}_R2.tmp{self.fidx}.fastq"
                    data = read2s[name]
                    with open(path2, 'a', encoding="utf-8") as out:
                        out.write("\n".join(data))
                        logger.debug(f"wrote demuliplex chunks to {path2}")


@dataclass
class BarMatchingI7(BarMatching):
    """Subclass of Barmatching that matches barcode in i7 header.

    Example 3RAD R1 file with i7 tag in header
    ------------------------------------------
    >>> # asterisk part is the i7 --->                  ********
    >>> @NB551405:60:H7T2GAFXY:4:21612:8472:20380 1:N:0:TATCGGTC+ACCAGGGA
    >>> ATCGGTATGCTGGAGGTGGTGGTGGTGGAGGTGGACGTTACAAGGGTTCTGGTGGTAGCCGATCAG...
    >>> +
    >>> EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE...
    """
    def _iter_matched_barcode(self) -> Iterator[Tuple[str, str, str]]:
        """Find barcode in read and check for match.
        
        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():
            # pull barcode from header
            barcode = read1[0].strip().rsplit(":", 1)[-1].split("+")[0]
            # look for match
            match = self.barcodes_to_names.get(barcode)

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1
                yield read1, read2, match
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1                


@dataclass
class BarMatchingCombinatorial(BarMatching):
    """Subclass of Barmatching for combinatorial inline barcodes.

    Example R1 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= restriction overhang.
    >>> ********-----
    >>> @E00526:227:H53YNCCX2:8:1202:7710:23354 1:N:0:
    >>> CTGCAACTATCGGAGCGAATGAAAC........GACTCAACATAACGGGTCTGATCATTGAG
    >>> +
    >>> AA<FFJJJJJJJJJJJJJJJJJJJJ........JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

    Example R2 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= restriction overhang.
    >>> ********----
    >>> @E00526:227:H53YNCCX2:8:1202:7446:23354 2:N:0:CGAACTGT+ACAACAGT
    >>> ATGCTGTCGATCCCAACCACCACGC........TTTTTTTCTATCTCAACTATTTACAACAA
    >>> +
    >>> AAFFFJJJJJJJJJFJFJJJJJJ-F........AFJ<JFJJJJAJFFAA-F<A-AAF-AFFJ
    """
    def _iter_matched_barcode(self):
        """Find barcode in read and check for match.
        
        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():
            
            # pull barcode from start of R1 (barcode1 + RE1 overhang) 
            match_r1 = None
            
            # pull barcode from start of R2 (barcode2 + RE2 overhang) 
            match_r2 = None

            # look for matches
            match = self.barcodes_to_names.get(barcode)

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1
                yield read1, read2, match, barcode
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1                



    # # COMBINATORIAL BARCODES (BCODE1+BCODE2)
    # if '3rad' in self.data.params.datatype:
    #     barcode1 = find3radbcode(self.cutters, self.longbar[0], read1)
    #     barcode2 = find3radbcode(self.cutters, self.longbar[2], read2)
    #     barcode = barcode1 + "+" + barcode2

    # # USE BARCODE PARSER: length or splitting
    # else:
    #     # Parse barcode. Uses the parsing function selected above.
    #     barcode = self.demux(self.cutters, read1, self.longbar)
    # yield read1, read2, barcode





def drop_from_right(path: Path, delim: str = "_", idx: int = 0) -> str:
    """Return a name with an underscore separated portion removed.

    This is used within `get_paired_fastqs` to find matching pairs
    when R1 and R2 are removed from file names.

    Example
    -------
    >>> path = Path("name_prefix_001_R1_002.fastq.gz")
    >>> drop_from_right(path, "_", 1)
    >>> # "name_prefix_001_002.fastq.gz"
    """
    # save and remove suffixes
    suffixes = path.suffixes
    while path.suffix in suffixes:
        path = path.with_suffix('')

    # break file name on delimiter and get chunks in reverse order
    chunks = path.name.split(delim)[::-1]

    # get chunks minus the index from the right
    sublist = [j for i, j in enumerate(chunks) if i != idx][::-1]
    path = path.parent / "_".join([i for i in sublist if i]).rstrip(delim)
    path = path.with_suffix("".join(suffixes))
    return path


if __name__ == "__main__":

    import ipyrad as ip
    from ipyrad.assemble.s1_demux import Step1
    ip.set_log_level("DEBUG")

    DATA = ip.Assembly("TEST1")
    DATA.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_R1*.gz"    
    DATA.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes.txt"
    DATA.params.project_dir = "/tmp"
    DATA.params.max_barcode_mismatch = 0
    # DATA.run('1', force=True, quiet=True)

    DATA.params.raw_fastq_path = "../../tests/ipsimdata/pairddrad_example_*.gz"    
    DATA.params.barcodes_path = "../../tests/ipsimdata/pairddrad_example_barcodes.txt"
    DATA.params.datatype = "pairddrad"

    # import glob
    # fastqs = [Path(i) for i in glob.glob(str(DATA.params.raw_fastq_path))]
    # print(fastqs)
    # print(drop_from_right(fastqs[0], "_", 1))


    # TEST i7 demux.
    DATA = ip.Assembly("TEST_i7")
    DATA.params.raw_fastq_path = "../../sandbox/radcamp/SMALL_RAW_R*.fastq"
    DATA.params.barcodes_path = "../../sandbox/radcamp/SMALL_i7_barcodes.txt"
    DATA.params.project_dir = "/tmp"
    DATA.params.max_barcode_mismatch = 1
    DATA.hackers.demultiplex_on_i7_tags = True

    fastqs = [Path(i) for i in glob.glob(str(DATA.params.raw_fastq_path))]
    print(fastqs)
    print(drop_from_right(fastqs[0], "_", 0))

    # with ip.Cluster(4) as ipyclient:
    #     step = Step1(DATA, force=True, quiet=False, ipyclient=ipyclient)
    #     tool = SimpleDemux(step.data, quiet=False, ipyclient=step.ipyclient)
    #     tool.run()

    #     print(tool.filenames_to_fastqs)

        # barm = BarMatchingI7(
        #     tool.data, 
        #     list(tool.filenames_to_fastqs.values())[0],
        #     tool.barcodes_to_names,
        #     tool.cutters,
        #     0,
        # )
        # barm.run()
        # for i in barm._iter_fastq_reads():
            # print(i)
        # self.data, self.barcodes_to_names, 
        # self.filenames_to_fastqs[fname],
        # self.cutters, self.barcodes_to_names, fidx)
