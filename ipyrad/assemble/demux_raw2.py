#!/usr/bin/env python

"""Some utilities used in demux.py for demultiplexing.
"""

from typing import Dict, Tuple, List, TypeVar, Iterator
import gzip
import itertools
import subprocess
from pathlib import Path
from collections import Counter
from dataclasses import dataclass

from loguru import logger
import pandas as pd
from pandas.errors import ParserError

from ipyrad.core.schema import SampleSchema
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.demux_raw_barmatching import *
from ipyrad.assemble.utils import IPyradError, AMBIGS, BADCHARS


Assembly = TypeVar("Assembly")
Client = TypeVar("Client")
BASES = set("ACGTN")
logger = logger.bind(name="ipyrad")

# pylint: disable=consider-using-with, too-many-nested-blocks


@dataclass
class SimpleDemux:
    data: Assembly
    ipyclient: Client
    quiet: bool

    # attributes to be filled.
    fastq_paths: List[Path] = None
    """: List of Paths to fastq files, unpaired."""
    barcodes_path: Path = None
    """: Path to the barcodes file."""
    names_to_barcodes: Dict[str, Tuple[str, str]] = None
    """: A map of barcode strings to sample names, pre-expanded by off-by-N."""
    filenames_to_fastqs: Dict[str, List[Tuple[str, str]]] = None
    """: Dict mapping file short names to tuples of paired fastqs."""
    cuts1: List[str] = None
    """: List of enzyme overhang sites to match on read1s."""
    cuts2: List[str] = None
    """: List of enzyme overhang sites to match on read2s."""
    barcodes_to_names: Dict[str, str] = None
    """: Dict of all acceptable barcodes (e.g., off-by-1) mapped to sample names."""
    file_stats: Dict[str, List] = None
    """: Store stats per raw data file (pair)."""

    def __post_init__(self):
        """Run subfunctions to setup object."""
        self._get_file_paths()
        self._get_filenames_to_fastqs()
        self._get_names_to_barcodes()
        self._replace_bad_name_chars()
        self._get_cutters_expanded()
        self._get_barcodes_to_names_map()

    def run(self):
        """Run each file (pair) on separatre engine."""
        self._distribute_remote_jobs()
        self._write_stats()

    def _get_file_paths(self) -> None:
        """Get fastq and barcodes file paths as Path objects."""
        raws = self.data.params.raw_fastq_path
        self.fastq_paths = list(raws.parent.glob(raws.name))
        if not self.fastq_paths:
            raise IPyradError(
                f"No fastq data found in {self.data.params.raw_fastq_path}")
        self.fastq_paths = [Path(i) for i in self.fastq_paths]

        # if regular expression pick first barcode path.
        bars = self.data.params.barcodes_path
        self.barcodes_path = list(bars.parent.glob(bars.name))
        if not self.barcodes_path:
            raise IPyradError(
                f"No barcodes file found at {self.data.params.barcodes_path}")
        self.barcodes_path = Path(self.barcodes_path[0])

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
                groups = {i.with_suffix("").name: sorted(j) for (i, j) in groups}
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

    def _get_names_to_barcodes(self) -> None:
        """Fill .names_to_barcodes dict w/ info from barcodes file.

        This logs a WARNING if technical replicates are detected to
        make sure the user is aware of how they are being handled.
        """
        # parse the tabular barcodes file on whitespace. Expects
        # there to be no header. There will be >=2 columns, >2 if
        # combinatorial barcodes.
        try:
            bardata = pd.read_csv(
                self.barcodes_path, header=None, delim_whitespace=True,
            ).dropna()
        except ParserError as err:
            raise IPyradError(
                "Failed to parse barcodes file. Check that your sample\n"
                "names do not include spaces (invalid)."
            ) from err

        # the dataframe COULD have >3 columns, in which case we will
        # discard any extra columns to keep at most 3.
        bardata = bardata.iloc[:, :3]

        # set names on barcodes dataframe
        if bardata.shape[1] == 2:
            bardata.columns = ["sample", "barcode1"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
        else:
            bardata.columns = ["sample", "barcode1", "barcode2"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
            bardata["barcode2"] = bardata["barcode2"].str.upper()

        # check for replicate sample names in the barcodes file. These
        # are allowed, since a single sample can be sequenced multiple
        # times on the same plate with different barcodes attached,
        # representing technical replicates. THere is a hackers option
        # for whether to combine tech reps, or keep as diff samples.
        if bardata['sample'].value_counts().max() > 1:
            # get duplicated names
            duplicated = (bardata['sample'].value_counts() > 1).index

            # warn that dups are present AND WILL BE merged.
            if self.data.hackers.merge_technical_replicates:
                logger.warning(
                    "\nTechnical replicates are present (samples with same name "
                    "in the barcodes file)\nand will be merged into one sample. "
                    "To instead NOT merge replicate samples set\n"
                    "`hackers.merge_technical_replicates = False), which will "
                    "instead rename \nsamples with a '-technical-replicate-x' "
                    "suffix.\nYou can change this in the 'hacker' settings "
                    "in the project JSON file.\n")
            # warn that dups are present and WILL NOT be merged.
            else:
                logger.warning(
                    "\nTechnical replicates are present (samples with same name "
                    "in the barcodes file)\nand WILL NOT be merged. To instead"
                    "merge replicate samples set \n"
                    "`hackers.merge_technical_replicates = True "
                    "which will combine reads from the same sample name\n"
                    "into a single sample.\n"
                    "You can change this in 'hacker' settings of the project JSON file.")

            # either way, relabel the samples for now, and may or may not merge later.
            for dup in duplicated:
                ridxs = bardata[bardata['sample'] == dup]
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
            self.names_to_barcodes = dict(zip(
                bardata["sample"], zip(bardata["barcode1"], bardata["barcode2"])
            ))
        else:
            self.names_to_barcodes = dict(zip(
                bardata["sample"], ((i, "") for i in bardata["barcode1"])
            ))
        # report to logger
        logger.debug(f"barcodes map:\n{bardata}")

    def _replace_bad_name_chars(self) -> None:
        """Replaces bad characters in names in .names_to_barcodes."""
        names = list(self.names_to_barcodes)
        for name in names:
            if any(i in name for i in BADCHARS):
                newname = name
                for badchar in BADCHARS:
                    newname = newname.replace(badchar, "_")
                # newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(f"changing name {name} to {newname} (bad characters).")
                self.names_to_barcodes[newname] = self.names_to_barcodes.pop(name)

    def _get_cutters_expanded(self) -> None:
        """Fill `.cuts1` and `.cuts2` with ordered list of resolutions.

        Sequences will be searched for cut sites starting with the
        entered value and then proceeding to allow off-by-n matches.
        The first tested sequences will be the user entered value,
        with IUPAC resolved, followed by off-by-1 matches.
        """
        cuts1 = self.data.params.restriction_overhang[0]
        if any(i in 'RKSYWM' for i in cuts1):
            res1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cuts1]
            res2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cuts1]
            cuts1 = [res1, res2]
        else:
            cuts1 = [cuts1]
        self.cuts1 = cuts1 + list(set(itertools.chain(*[mutate(i) for i in cuts1])))

        cuts2 = self.data.params.restriction_overhang[1]
        if any(i in 'RKSYWM' for i in cuts2):
            res1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cuts2]
            res2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cuts2]
            cuts2 = [res1, res2]
        else:
            cuts2 = [cuts2]
        self.cuts2 = cuts2 + list(set(itertools.chain(*[mutate(i) for i in cuts2])))
        # logger.info(self.cuts1)
        # logger.info(self.cuts2)

        # convert all str to bytes
        self.cuts1 = [i.encode() for i in self.cuts1]
        self.cuts2 = [i.encode() for i in self.cuts2]

    def _get_barcodes_to_names_map(self) -> None:
        """Fills .barcodes_to_names with all acceptable barcodes: name.

        This updates the .barcodes_to_names from {str: Tuple[str,str]}
        to {str: str}.
        """
        # store perfect match to barcodes
        self.barcodes_to_names = {}

        # finished if no mismatch is allowed.
        if not self.data.params.max_barcode_mismatch:
            for name, barcode in self.names_to_barcodes.items():
                # convert tuple to string with _ separator
                barc = (
                    f"{barcode[0]}" if not barcode[1] else
                    f"{barcode[0]}_{barcode[1]}"
                )
                self.barcodes_to_names[barc.encode()] = name
            return

        # iterate over barcodes: names
        for name, barcode in self.names_to_barcodes.items():

            # get generators of off-by-n barcodes
            if self.data.params.max_barcode_mismatch == 1:
                gen1 = mutate(barcode[0])
                gen2 = mutate(barcode[1])
            else:
                gen1 = itertools.chain(*[(mutate(i)) for i in mutate(barcode[0])])
                gen2 = itertools.chain(*[(mutate(i)) for i in mutate(barcode[1])])
            bars1 = set(gen1)
            bars2 = set(gen2)

            # if only one barcode
            if not bars2:
                barcgen = iter(bars1)
            else:
                barcgen = (f"{i}_{j}" for (i, j) in itertools.product(bars1, bars2))

            warning = False
            for barc in barcgen:
                barc = barc.encode()
                if barc not in self.barcodes_to_names:
                    self.barcodes_to_names[barc] = name
                else:
                    logger.warning(
                        f"\nSample: {name} ({barc}) is within "
                        f"{self.data.params.max_barcode_mismatch} "
                        f"base changes of sample ({self.barcodes_to_names[barc]}).")
                    warning = True
            if warning:
                logger.warning(
                    "Ambiguous barcodes that match to multiple samples "
                    "will arbitrarily be assigned to the first sample.\n"
                    "If you do not like this then lower the value of "
                    "max_barcode_mismatch and rerun (recommended).")

    def _distribute_remote_jobs(self) -> None:
        """Send barcode matching jobs to remote engines."""

        # limit the max number of engines to ... 4?
        lbview = self.ipyclient.load_balanced_view(targets=self.ipyclient.ids[:4])

        # barmatching
        jobs = {}
        for fidx, fname in enumerate(self.filenames_to_fastqs):
            fastqs = self.filenames_to_fastqs[fname]
            args = (self.data, fastqs, self.barcodes_to_names, self.cuts1, self.cuts2, fidx)
            jobs[fname] = lbview.apply(barmatch, *args)
        msg = "demultiplexing reads"
        prog1 = AssemblyProgressBar(jobs, msg, step=1, quiet=self.quiet)
        prog1.update()
        prog1.block()
        prog1.check()

        # concatenating tmpfiles
        jobs = {}
        for name in self.names_to_barcodes:
            # if reps: submit only one job per set of technical replicates
            if self.data.hackers.merge_technical_replicates:
                name = name.split("-technical-replicate-")[0]
            if name not in jobs:
                jobs[name] = lbview.apply(concatenate_tmpfiles, *(self.data, name))
        msg = "concatenating chunked files"
        prog2 = AssemblyProgressBar(jobs, msg, step=1, quiet=self.quiet)
        prog2.update()
        prog2.block()
        prog2.check()

        # create samples and store fastq paths
        for name, fastqs in prog2.results.items():
            sample = SampleSchema(name=name)
            sample.files.fastqs.append(tuple(str(i) for i in fastqs))
            self.data.samples[name] = sample

        # store stats for writing the verbose output file
        self.file_stats = prog1.results

        # write stats to the sample
        for _, stats in self.file_stats.items():
            for sname, hits in stats[2].items():
                if self.data.hackers.merge_technical_replicates:
                    sname = sname.split("-technical-replicate-")[0]
                self.data.samples[sname].stats_s1.reads_raw += hits

        # set state to 1 on counts after counts concat above
        names = list(self.data.samples)
        for name in names:
            sample = self.data.samples[name]
            if sample.stats_s1.reads_raw:
                sample.state = 1
            else:
                self.data.samples.pop(name)
                logger.warning(f"sample {name} has 0 reads and will be excluded.")
        logger.info(f"created {len(self.data.samples)} new samples")
        self.data.save_json()

    def _write_stats(self):
        """Write to {project_dir}/`s1_demultiplex_stats.txt`.

        The stats file includes the number of reads per sample as well
        as information about demultiplexing in terms of nreads per file
        and the barcodes that were found.
        """
        # open the stats file for writing.
        stats_file = self.data.stepdir / "s1_demultiplex_stats.txt"
        outfile = open(stats_file, 'w', encoding="utf-8")

        # write the per-file stats
        outfile.write("# Raw file statistics\n######################\n")
        file_df = pd.DataFrame(
            index=sorted(self.file_stats),
            columns=["total_reads", "cut_found", "bar_matched"],
        )
        for key in sorted(self.file_stats):
            stats = self.file_stats[key]
            not_cut = sum(stats[0].values())
            matched = sum(stats[1].values())
            total = not_cut + matched
            file_df.loc[key, :] = total, total - not_cut, matched
        outfile.write(file_df.to_string() + "\n\n")

        # write sample nreads stats ----------------------------------
        outfile.write("# Sample demux statistics\n######################\n")
        sample_df = pd.DataFrame(
            index=sorted(self.data.samples),
            columns=["reads_raw"],
            data=[
                self.data.samples[i].stats_s1.reads_raw
                for i in sorted(self.data.samples)
            ],
        )
        outfile.write(sample_df.to_string() + "\n\n")
        logger.info("\n" + sample_df.to_string())

        # write verbose barcode information --------------------------
        outfile.write("# Barcode detection statistics\n######################\n")

        # record matches
        data = []
        bar_obs = Counter()
        for key in self.file_stats:
            bar_obs.update(self.file_stats[key][1])
        sorted_bar_obs = sorted(bar_obs, key=lambda x: bar_obs[x], reverse=True)
        for name, truebar in self.names_to_barcodes.items():
            for foundbar in sorted_bar_obs:
                if name == self.barcodes_to_names[foundbar]:
                    count = bar_obs[foundbar]
                    if count:
                        if "_" in foundbar:
                            foundbar = tuple(foundbar.split("_"))
                        else:
                            truebar = truebar[0]
                        data.append([name, truebar, foundbar, count])

        # record misses
        bad_bars = Counter()
        for key in sorted(self.file_stats):
            bad_bars.update(self.file_stats[key][0])
        bad_bar_obs = sorted(bad_bars, key=lambda x: bad_bars[x], reverse=True)
        for badbar in bad_bar_obs:
            count = bad_bars[badbar]
            if "_" in badbar:
                badbar = tuple(badbar.split("_"))
            data.append(["no_match", "", badbar, count])
        barcodes_df = pd.DataFrame(
            index=[i[0] for i in data],
            columns=["true_bar", "observed_bar", "N_records"],
            data=[i[1:] for i in data],
        )
        outfile.write(barcodes_df.to_string() + "\n")


######################################################################
######################################################################
##
##  Functions run on remote engines
##
######################################################################
######################################################################


def barmatch(data, fastqs, barcodes_to_names, cuts1, cuts2, fidx):
    """Starts barmatch process using the appropriate subclass."""
    if data.hackers.demultiplex_on_i7_tags:
        barmatcher = BarMatchingI7(
            data, fastqs, barcodes_to_names, cuts1, cuts2, fidx)
    elif b"_" in list(barcodes_to_names)[0]:
        barmatcher = BarMatchingCombinatorialInline(
            data, fastqs, barcodes_to_names, cuts1, cuts2, fidx)
    else:
        barmatcher = BarMatchingSingleInline(
            data, fastqs, barcodes_to_names, cuts1, cuts2, fidx)
    barmatcher.run()
    return barmatcher.barcode_misses, barmatcher.barcode_hits, barmatcher.sample_hits


def concatenate_tmpfiles(data: "Assembly", name: str) -> Tuple[str,str]:
    """write tmpfiles to stepdir."""
    r1s = list(data.tmpdir.glob(f"{name}_R1.tmp*.fastq.gz"))
    r2s = list(data.tmpdir.glob(f"{name}_R2.tmp*.fastq.gz"))

    # if not data was written for this sample return empties
    if not r1s:
        # printed msgs on engines are logged to INFO
        return ("", "")

    # get final output names
    bits = r1s[0].name.rsplit(".", 3)
    r1out = data.stepdir / ".".join(bits[:1] + bits[2:])
    r2out = ""
    if r2s:
        bits = r2s[0].name.rsplit(".", 3)
        r2out = data.stepdir / ".".join(bits[:1] + bits[2:])

    # if only one file then just rename and move to final dir
    if len(r1s) == 1:
        for tmpfile, outname in zip(r1s + r2s, [r1out, r2out]):
            tmpfile.rename(data.stepdir / outname)

    # multiple files: concatenate.
    else:
        with gzip.open(r1out, 'wb') as out:
            cmd = ["cat"] + r1s
            subprocess.run(cmd, check=True, stdout=out)
        if r2s:
            with gzip.open(r2out, 'wb') as out:
                cmd = ["cat"] + r2s
                subprocess.run(cmd, check=True, stdout=out)
    return (r1out, r2out)


def mutate(barcode: str) -> Iterator[str]:
    """Mutate a sequence by 1 base (ACGT). Used for barcode mismatch."""
    for pos, _ in enumerate(barcode):
        for sub in BASES:
            newbar = list(barcode)
            newbar[pos] = sub
            yield "".join(newbar)


def drop_from_right(path: Path, delim: str = "_", idx: int = 0) -> str:
    """Return a name with an underscore separated portion removed.

    This is used within `_get_filenames_to_fastqs` to find matching
    pairs when R1 and R2 are removed from file names.

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
    ip.set_log_level("DEBUG")


    DATA = ip.Assembly("TEST1")
    DATA.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_R1*.gz"
    DATA.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes.txt"
    DATA.params.project_dir = "/tmp"
    DATA.params.max_barcode_mismatch = 0
    DATA.run('1', force=True, quiet=True)
    print(DATA.stats)

    # DATA.params.raw_fastq_path = "../../tests/ipsimdata/pairddrad_example_*.gz"
    # DATA.params.barcodes_path = "../../tests/ipsimdata/pairddrad_example_barcodes.txt"
    # DATA.params.datatype = "pairddrad"

    # # import glob
    # # fastqs = [Path(i) for i in glob.glob(str(DATA.params.raw_fastq_path))]
    # # print(fastqs)
    # # print(drop_from_right(fastqs[0], "_", 1))


    # # TEST i7 demux.
    # DATA = ip.Assembly("TEST_i7")
    # DATA.params.raw_fastq_path = "../../sandbox/radcamp/SMALL_RAW_R*.fastq"
    # DATA.params.barcodes_path = "../../sandbox/radcamp/SMALL_i7_barcodes.txt"
    # DATA.params.project_dir = "/tmp"
    # DATA.params.max_barcode_mismatch = 1
    # DATA.hackers.demultiplex_on_i7_tags = True

    # FASTQS = [Path(i) for i in glob.glob(str(DATA.params.raw_fastq_path))]
    # print(FASTQS)
    # print(drop_from_right(FASTQS[0], "_", 0))

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
