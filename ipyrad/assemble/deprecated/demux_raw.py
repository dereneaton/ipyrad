#!/usr/bin/env python

"""Demultiplexing raw_fastq_path input data given inline or i7 barcodes.

Demultiplexing uses few cores since it is mostly i/o limited.
"""

import glob
from pathlib import Path
import itertools
from collections import Counter
from typing import Dict, List, Tuple
import pandas as pd
from loguru import logger
from ipyrad.core.schema import SampleSchema, Stats1
from ipyrad.assemble.utils import IPyradError, BADCHARS, AMBIGS
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.demux_raw_utils import BarMatch, collate_files

logger = logger.bind(name="ipyrad")

class SimpleDemux:
    def __init__(self, step):

        # store the step object inputs
        self.data = step.data
        self.quiet = step.quiet

        # get only every 4th available engine since i/o limits
        self.lbview = step.ipyclient.load_balanced_view(
            targets=step.ipyclient.ids[::4])

        # attrs to be filled (see type hints)
        self.fastqs: List[Path]
        self.ftuples: Dict[str,List[Tuple[str,str]]] = {}
        self.barcodes: Dict[str,str] = {}
        self.matchdict: Dict[str,str] = {}
        self.cutters: List[str] = []
        self.longbar: int = None
        self.file_stats: Dict[str,List] = {}

        # parse the barcodes file with inline and/or i7 tags
        # e.g., (i7 demuxing or merge technical replicates)
        self._get_file_paths()
        self._link_barcodes()
        self.update_barcodes_for_i7()
        self.check_sample_names()
        self.parse_file_names_and_pairs()
        self.expand_ambiguous_cutters()
        self.get_barcode_dict()


    def _get_file_paths(self):
        """Check .fastqs files are present and format as .ftuples."""
        fastqs = glob.glob(str(self.data.params.raw_fastq_path))
        if not fastqs:
            raise IPyradError(
                "No data found in {}. Fix path to data files."
                .format(self.data.params.raw_fastq_path))
        self.fastqs = [Path(i) for i in fastqs]

        # if regular expression pick first barcode path.
        barcode_path = glob.glob(str(self.data.params.barcodes_path))
        if not barcode_path:
            raise IPyradError("Barcodes file not found. You entered: "
                f"{self.data.params.barcodes_path}")
        self.barcode_path = Path(barcode_path[0])

    def _link_barcodes(self):
        """Parse Sample barcodes to dict using pandas.
        """
        # read in the file
        bdf = pd.read_csv(self.barcode_path, header=None, delim_whitespace=True)
        bdf = bdf.dropna()

        # make sure bars are upper case
        bdf[1] = bdf[1].str.upper()

        # if replicates are present
        if bdf[0].value_counts().max() > 1:

            # get duplicated names
            repeated = (bdf[0].value_counts() > 1).index

            # adds label to ONLY replicates which will stay permanently
            if not self.data.hackers.merge_technical_replicates:

                # print a warning about dups if the data are not demultiplexed
                logger.warning(
                    "Technical replicates (same name) are present and will not"
                    " be merged (hackers.merge_technical_replicates = False), "
                    "but instead renamed with a '-technical-replicate-x' "
                    "suffix. You can change this in the 'hacker' settings "
                    "in the project JSON file.")

            # adds tmp label to ALL sample which will be removed by end of s1
            else:
                # print a warning about dups if the data are not demultiplexed
                logger.warning(
                    "Technical replicates (same name) are present and will be "
                    "merged (hackers.merge_technical_replicates = True). "
                    "You can change this under the 'hacker' settings "
                    "in the project JSON file.")

            # either way we relabel them for now, but merge or not later
            for rep in repeated:
                farr = bdf[bdf[0] == rep]
                if farr.shape[0] > 1:
                    for idx, index in enumerate(farr.index):
                        bdf.loc[index, 0] = (
                            "{}-technical-replicate-{}".format(rep, idx))

        # make sure barcode chars are all proper
        if not all(bdf[1].apply(set("RKSYWMCATG").issuperset)):
            raise IPyradError(
                "Barcodes file contains unexpected characters in the "
                "barcode sequences suggesting it is not correctly "
                "formatted. See documentation.")

        # store barcodes as a dict
        self.barcodes = dict(zip(bdf[0], bdf[1]))

        # change dtype from 3RAD to pairddrad if i7 demuxing to avoid
        # looking for a second barcode.
        if self.data.hackers.demultiplex_on_i7_tags:
            self.data.params.datatype = "pairddrad"
            logger.info("setting datatype to pairrdrad for i7 demux.")

        # 3rad/seqcap use multiplexed barcodes
        if "3rad" in self.data.params.datatype:
            if not bdf.shape[1] == 3:
                raise IPyradError(
                    "pair3rad datatype should have two barcodes per sample.")

            # We'll concatenate them with a plus and split them later
            bdf[2] = bdf[2].str.upper()
            self.barcodes = dict(zip(bdf[0], bdf[1] + "+" + bdf[2]))
        else:
            self.barcodes = dict(zip(bdf[0], bdf[1]))
        logger.debug(f"barcodes: {self.barcodes}")

    def update_barcodes_for_i7(self):
        """updates longbar for the length of ..."""
        # Handle 3rad multi-barcodes. Gets len of the first one.
        blens = [len(i.split("+")[0]) for i in self.barcodes.values()]
        if len(set(blens)) == 1:
            self.longbar = (blens[0], 'same')
        else:
            self.longbar = (max(blens), 'diff')
        logger.debug("blens {}".format(blens))
        logger.debug("longbar {}".format(self.longbar))

        # i7 tags there will be only one barcode, so this overrides "datatype"
        # so that if you are using pair3rad if doesn't cause problems.
        # For pair3rad we need to add the length info for barcodes_R2
        if not self.data.hackers.demultiplex_on_i7_tags:
            if "3rad" in self.data.params.datatype:
                blens = [
                    len(i.split("+")[1]) for i in self.barcodes.values()
                ]
                self.longbar = (self.longbar[0], self.longbar[1], max(blens))

    def check_sample_names(self):
        """Replace sample names with bad characters."""
        snames = sorted(self.barcodes)
        for sname in snames:
            if any(i in sname for i in BADCHARS):
                newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(
                    f"changing name {sname} to {newname} (hard characters).")
                self.barcodes[newname] = self.barcodes.pop(sname)

    def parse_file_names_and_pairs(self):
        """Check for PE matching and extract sample names to ftuples."""
        def drop_from_right(chunks, idx):
            """
            used in pair name parsing to sub out _ delim chunks

            Example:
            ---------
            if idx = 1 then:
                name_prefix_001_R1_002.fastq.gz           # 1.
                ['name_prefix', '001', 'R1', '002']       # 2.
                ['name_prefix', '001', '002']             # 3.
            """
            sublist = [j for i, j in enumerate(chunks[::-1]) if i != idx][::-1]
            return "_".join([i for i in sublist if i]).rstrip("_")

        # check if file names end with _ before the suffix and split
        # on two underscores, else split on last one.
        bases = sorted([
            i.name.rstrip(".gz").rstrip(".fastq").rstrip(".fq")
            for i in self.fastqs
        ])

        # link pairs into tuples
        if self.data.is_pair:

            # try these in order until ngroups == nfiles / 2
            idx = 0
            while 1:
                try:
                    # get groups up to an underscore delimiter
                    groups = itertools.groupby(
                        bases,
                        key=lambda x: drop_from_right(x.split("_"), idx),
                    )
                    groups = {i: list(j) for i, j in groups}
                    assert len(groups) == len(self.fastqs) / 2
                    assert all(len(j) == 2 for i, j in groups.items())
                    logger.debug(f"using '_' rsplit = {idx}")
                    break
                except Exception:
                    pass
                # increase counter up to 5 _ back from end, then raise
                idx += 1
                if idx > 5:
                    raise IPyradError(
                        "Cannot parse paired file names. File names must have "
                        "matching name prefix followed by _1 _2, _R1 _R2, "
                        "or _R1_ _R2_ followed by any subsequent suffix. "
                        f"Your data files look like this: {self.fastqs}"
                        )

            # apply splitter to the full path names
            groups = itertools.groupby(
                sorted(self.fastqs),
                key=lambda x: drop_from_right(x.split("_"), idx),
            )

            for fname, files in groups:
                fname = fname.name
                files = sorted(files)
                logger.debug(f"detected paired files: {files}")
                self.ftuples[fname] = (files[0], files[1])

            # file checks
            if not self.ftuples:
                raise IPyradError(
                    "No paired fastq files found. File names must have "
                    "matching name prefix followed by _1 _2, _R1 _R2, "
                    "or _R1_ _R2_.")

        # data are not paired, create empty tuple pair
        else:
            # print warning if _R2_ is in names when not paired
            endings = ("_R2", "_2", "_R2")
            warning = []
            for base in bases:
                if any(base.endswith(i) for i in endings):
                    warning.append(base)
            if warning:
                message = (
                    "Input file names look suspiciously like paired-end "
                    "data but you selected single end. If so, you should "
                    "set the parameter 'datatype' to a paired option (e.g., "
                    "pairddrad or pairgbs) and re-run step 1, which will "
                    "require using the force flag (-f) to overwrite "
                    "existing data.\n{}".format(self.fastqs))
                logger.warning(message)

            for path in self.fastqs:
                fname = path.name.rstrip('.gz').rstrip('.fastq').rstrip('.fq')
                self.ftuples[fname] = (path, "")
        logger.debug(f"ftuples: {self.ftuples}")

    def expand_ambiguous_cutters(self):
        """Checks sample names and replaces bad chars in dict with _
        And returns a list of both resolutions of cut site 1 for ambigs.
        # (TGCAG, ) ==> [TGCAG, ]
        # (TWGC, ) ==> [TAGC, TTGC]
        # (TWGC, AATT) ==> [TAGC, TTGC]
        """
        # expand ambigs
        for cutter in self.data.params.restriction_overhang:
            if not any(i in "RKSYWM" for i in cutter):
                self.cutters.append([cutter, ""])
            else:
                cutter1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cutter]
                cutter2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cutter]
                self.cutters.append([cutter1, cutter2])

        assert self.cutters, (
            "Must enter a restriction_overhang for demultiplexing.")

    def get_barcode_dict(self):
        """Build full inverse barcodes dictionary."""
        bases = set("CATGN")
        poss = set()

        # do perfect matches
        for sname, barcode in self.barcodes.items():

            # store {barcode: name} mapping
            self.matchdict[barcode] = sname

            # record that this barcodes has been seen
            poss.add(barcode)

            # get n-off barcodes
            if self.data.params.max_barcode_mismatch:

                # iterate over bases in the barcode
                for idx1, base in enumerate(barcode):

                    # get bases that are not this one
                    diffs = bases.difference(base)

                    # iter over the ways this barcode has other bases as this pos.
                    for diff in diffs:
                        lbar = list(barcode)
                        lbar[idx1] = diff
                        tbar1 = "".join(lbar)

                        # if this new barcode has not been observed store it.
                        if tbar1 not in poss:
                            self.matchdict[tbar1] = sname
                            poss.add(tbar1)

                        # if it has been seen in another taxon, problem.
                        else:
                            logger.warning(
                                f"\nSample: {sname} ({barcode}) is within "
                                f"{self.data.params.max_barcode_mismatch} "
                                f"base changes of sample ({self.matchdict.get(barcode)}). "
                                "Ambiguous barcodes that match to both samples "
                                "will arbitrarily be assigned to the first "
                                "sample. If you do not like this then lower "
                                "the value of max_barcode_mismatch and rerun "
                                "(recommended).")

                    # if allowing two base difference things get big
                    # for each modified bar, allow one modification to other bases
                    if self.data.params.max_barcode_mismatch > 1:
                        for idx2, _ in enumerate(tbar1):
                            # skip the base that is already modified
                            if idx2 != idx1:
                                for diff in bases.difference(tbar1[idx2]):
                                    ltbar = list(tbar1)
                                    ltbar[idx2] = diff
                                    tbar2 = "".join(ltbar)
                                    if tbar2 not in poss:
                                        self.matchdict[tbar2] = sname
                                        poss.add(tbar2)
                                    else:
                                        if self.matchdict.get(tbar2) != sname:
                                            logger.warning(
                                f"\nSample: {sname} ({barcode}) is within "
                                f"{self.data.params.max_barcode_mismatch} "
                                f"base changes of sample ({self.matchdict.get(tbar2)}). "
                                "Ambiguous barcodes that match to both samples "
                                "will arbitrarily be assigned to the first "
                                "sample. If you do not like this then lower "
                                "the value of max_barcode_mismatch and rerun "
                                "(recommended).")
        logger.debug("cutters: {}".format(self.cutters))
        logger.debug("matchdict: {}...".format(str(self.matchdict)[:50]))

    def run(self):
        """Run complete steps"""
        self.remote_run_barmatch()
        self.concatenate_chunks()
        self.write_json()
        self.write_stats()

    def remote_run_barmatch(self):
        """Submit each file to be demux'd on a separate engine."""
        # remote function
        def barmatch(args):
            tool = BarMatch(*args)
            stats = tool.run()
            return stats

        # submit jobs
        rasyncs = {}
        for fidx, fname in enumerate(self.ftuples):
            args = (
                self.data,
                self.barcodes,
                self.ftuples[fname],
                self.longbar,
                self.cutters,
                self.matchdict,
                fidx
            )
            rasyncs[fname] = self.lbview.apply(barmatch, args)

        # progress bar info
        message = "sorting reads"
        prog = AssemblyProgressBar(rasyncs, message, step=1, quiet=self.quiet)
        prog.update()
        prog.block()
        prog.check()

        # collect and store results on per-file stats as jobs finish
        for fname in prog.results:
            self.file_stats[fname] = prog.results[fname]
        del rasyncs

    def concatenate_chunks(self):
        """
        If multiple chunk files match to the same sample name but with
        different barcodes (i.e., they are technical replicates) then this
        will assign all the files to the same sample name file and gzip it.
        """
        # get all the files
        ftmps = self.data.tmpdir.glob("tmp_*.fastq")

        # a dict to assign tmp files to names/reads
        r1dict = {}
        r2dict = {}

        # assign to name keys
        for ftmp in ftmps:
            ftmp_name = ftmp.name.split("_", 1)[1]
            tname, orient, _ = ftmp_name.rsplit("_", 2)

            # only merge technical reps if settings
            if self.data.hackers.merge_technical_replicates:
                if "-technical-replicate-" in tname:
                    sname = tname.rsplit("-technical-replicate", 1)[0]
                else:
                    sname = tname
            else:
                sname = tname

            # assign r1 and r2 to paired dicts
            if orient == "R1":
                if sname in r1dict:
                    r1dict[sname].append(ftmp)
                else:
                    r1dict[sname] = [ftmp]
            else:
                if sname in r2dict:
                    r2dict[sname].append(ftmp)
                else:
                    r2dict[sname] = [ftmp]

        # submit jobs to remote engines
        rasyncs = {}
        for sname in r1dict:
            tmp1s = sorted(r1dict.get(sname))
            tmp2s = sorted(r2dict.get(sname) if sname in r2dict else [])
            args = (self.data, sname, tmp1s, tmp2s)
            rasyncs[sname] = self.lbview.apply(collate_files, *args)

        # collate files progress bar
        message = "writing/compressing"
        prog = AssemblyProgressBar(rasyncs, message, 2, self.quiet)
        prog.block()
        prog.check()

    def write_json(self):
        """Create samples and save the project JSON file with new file
        paths and stats for samples.
        """
        # combine results across computed chunks
        sample_counter = Counter()
        for key in sorted(self.file_stats):
            sample_counter.update(self.file_stats[key][2])

        # combine results among technical replicates and create samples
        self.data.samples = {}
        for tname in sample_counter:

            # group replicates if in settings
            if self.data.hackers.merge_technical_replicates:
                sname = tname.split("-technical-replicate-")[0]
            else:
                sname = tname

            # create new sample
            if sname not in self.data.samples:
                stats_s1 = Stats1(reads_raw=sample_counter[tname])
                sample = SampleSchema(name=sname, stats_s1=stats_s1)
                sample.files.fastqs = [(
                    self.data.stepdir / f"{sname}_R1.fastq.gz",
                    self.data.stepdir / f"{sname}_R2.fastq.gz" if self.data.is_pair else "",
                )]
                self.data.samples[sname] = sample
            else:
                self.data.samples[sname].stats_s1.reads_raw += sample_counter[tname]

        # drop any samples with zero data and give warning message.
        for sname, sample in self.data.samples:
            if not sample.stats_s1.reads_raw:
                logger.warning(f"sample {sname} has 0 reads and will be excluded.")
        self.data.samples = {
            i: sample for i, sample in self.data.samples.items()
            if sample.stats_s1.reads_raw
        }

        # save project json
        logger.info(f"created {len(self.data.samples)} new samples")
        self.data.save_json()

    def write_stats(self):
        """
        The stats file includes the number of reads per sample as well
        as information about demultiplexing in terms of nreads per file
        and the barcodes that were found.
        """
        # open the stats file for writing.
        stats_file = self.data.stepdir / "s1_demultiplex_stats.txt"
        outfile = open(stats_file, 'w', encoding="utf-8")

        # write the per-file stats
        outfile.write("# Raw file statistics\n# -------------------\n")
        file_df = pd.DataFrame(
            index=[i.name for i in sorted(self.file_stats)],
            columns=["total_reads", "cut_found", "bar_matched"],
        )
        for key in sorted(self.file_stats):
            stats = self.file_stats[key]
            not_cut = stats[4].get("_") if "_" in stats[4] else 0
            file_df.loc[key, :] = stats[0], stats[0] - not_cut, stats[1]
        outfile.write(file_df.to_string() + "\n\n")

        # write sample nreads stats ----------------------------------
        outfile.write("# Sample demux statistics\n# -----------------------\n")
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
        outfile.write("# Barcode detection statistics\n# ----------------------------\n")
        bar_obs = Counter()
        for key in sorted(self.file_stats):
            bar_obs.update(self.file_stats[key][3])
        sorted_bar_obs = sorted(bar_obs, key=lambda x: bar_obs[x])
        data = []
        for tname in sorted(self.barcodes):
            true_barcode = self.barcodes[tname]
            for barcode in sorted_bar_obs:
                if self.matchdict[barcode] == tname:
                    count = bar_obs[barcode]
                    if count:
                        data.append([tname, true_barcode, barcode, count])
        bad_bars = Counter()
        for key in sorted(self.file_stats):
            bad_bars.update(self.file_stats[key][4])
        bad_bar_obs = sorted(bad_bars, key=lambda x: bad_bars[x], reverse=True)
        for barcode in bad_bar_obs:
            data.append(["no_match", "", barcode, bad_bars[barcode]])
        barcodes_df = pd.DataFrame(
            index=[i[0] for i in data],
            columns=["true_bar", "observed_bar", "N_records"],
            data=[i[1:] for i in data],
        )
        outfile.write(barcodes_df.to_string() + "\n")


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("INFO")

    # DEMUX ON I7 TAGS FOR PAIR3RAD TEST
