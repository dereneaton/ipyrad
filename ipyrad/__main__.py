#!/usr/bin/env python

"""Command line interface.

Examples
--------
>>> ipyrad demux --data RAW.fastq --barcodes BARS.csv --...
>>> ipyrad json -n NEW --project_dir /tmp --view
>>> ipyrad assemble -j NEW.json -s 123
>>> ipyrad branch -p params-NEW.txt -n NEW2 --clust-threshold 0.9 --...
>>> ipyrad stats -p params-NEW.txt

>>> ipyrad demux --data RAW.fastq --barcodes BARS.csv -m 1 -o /tmp/TEST
>>> ipyrad new -n NEW --project_dir /tmp --view
>>> ipyrad assemble -j NEW.json -s 123
>>> ipyrad branch -p params-NEW.txt -n NEW2 --clust-threshold 0.9 --...
>>> ipyrad stats -p params-NEW.txt
"""
import argparse
from pathlib import Path
import ipyrad as ip
from loguru import logger

logger = logger.bind(name="ipyrad")
PARAMS = ip.schema.Params.model_fields

# VERSION = str(get_distribution('ipyrad')).split()[1]
VERSION = str(ip.__version__)
HEADER = f"""
-------------------------------------------------------------
 ipyrad [v.{VERSION}]
 Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = " ipyrad command line tool. Select a positional subcommand:"
EPILOG = """\
Note
----
Each subcommand has its own additional help screen, e.g.,:
>>> ipyrad demux -h

Examples
--------
>>> # demux: demultiplexing data to samples by index or barcode
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv -c 10 -o ./demux
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --logger DEBUG
>>> ipyrad demux -d RAWS1/*.fastq.gz RAWS2/*.fastqs.gz -b BARCODES.tsv

>>> # new: create a new assembly file and set parameters
>>> ipyrad new -n NAME --fastq_paths ./demux/*.gz
>>> ipyrad new -n NAME --fastq_paths ./demux/*.gz --project_dir ./ ...

>>> # assemble: run assembly steps 1-7 (filtering, clustering, variants, formatting)
>>> ipyrad assemble -j NAME.json -s 1
>>> ipyrad assemble -j NAME.json -s 123 -c 40 -t 4
>>> ipyrad assemble -j NAME.json -s 456 -f --logger DEBUG

>>> # branch: create a copy/branch of an assembly to re-run w/ diff params
>>> ipyrad branch -j NAME.json -n NEW --cluster_threshold 0.9
>>> ipyrad assemble -j NEW.json -s 7 -f

>>> # merge: create new assembly containing union of samples from >=2 others.
>>> ipyrad merge -n NEW -j NAME1.json NAME2.json
>>> ipyrad assemble -j NEW.json -s 67 -f

>>> # print: show assembly parameters, stats, or file paths.
>>> ipyrad print -j NAME.json --stats
>>> ipyrad print -j NAME.json --params
"""

DEMUX_EPILOG = """\
Examples
--------
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv -c 10 -o ./demux
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --logger DEBUG
>>> ipyrad demux -d RAWS1/*.fastq.gz RAWS2/*.fastqs.gz -b BARCODES.tsv
"""

ASSEMBLE_EPILOG = """\
Examples
--------
>>> ipyrad assemble -j NAME.json -s 1
>>> ipyrad assemble -j NAME.json -s 123 -c 10 -t 2
>>> ipyrad assemble -j NAME.json -s 56 -c 40 -f --logger INFO
"""

BRANCH_EPILOG = """\
Examples
--------
>>> ipyrad branch -j NAME.json -n NEW
>>> ipyrad branch -j NAME.json -n NEW -s A B C D E
>>> ipyrad branch -j NAME.json -n NEW -e F G H
"""

MERGE_EPILOG = """\
Examples
--------
>>> ipyrad merge -n NEW -j NAME1.json NAME2.json
"""

STATS_EPILOG = """\
Examples
--------
>>> ipyrad print -j NAME.json --stats
"""

NEW_EPILOG = """\
Examples
--------
>>> ipyrad new -n name
>>> ipyrad new -n clust90 --clust_threshold 0.9
>>> ipyrad new -n clust90min4 --clust_threshold 0.9 --min_samples_locus 4

"""


def setup_demux_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad demux` subcommand parser.

    """
    demux = subparsers.add_parser(
        "demux",
        description=HEADER + "\n" + " ipyrad demux: demultiplex reads by index/barcode",
        help="Demultiplex unsorted data by index or barcode.",
        epilog=DEMUX_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    demux.add_argument(
        "-o", metavar="outpath", type=Path, default="./demux_fastqs",
        help="Path to use for output directory, e.g., './demux_fastqs'. Error if dir exists.",
    )
    demux.add_argument(
        "-d", metavar="data", type=Path, required=True, nargs="*",
        help="One or more paths to select fastq data file using regex (e.g., './data/*.fastq.gz')",
    )
    demux.add_argument(
        "-b", metavar="barcodes", type=Path, required=True,
        help=(
            "Path to a barcodes file (CSV, TSV, or whitespace delimited) "
            "where each line contains (name, barcode1, [optional barcode2]).")
    )
    demux.add_argument(
        "-re1", metavar="re1", type=str,
        help=(
            "The overhang on read1 from restriction digestion and ligation. "
            "If no value is entered then it will be estimated from the sequences.")
    )
    demux.add_argument(
        "-re2", metavar="re2", type=str,
        help=(
            "The overhang on read2 from restriction digestion and ligation. "
            "If no value is entered then it will be estimated from the sequences.")
    )
    demux.add_argument(
        "-m", "--max-mismatch", metavar="max_mismatch", type=int, default=0,
        help=(
            "The maximum number of allowed mismatches between true and "
            "oberved barcodes (Default=0).")
    )
    demux.add_argument(
        "-c", "--cores", metavar="cores", type=int, default=4,
        help="Number of cores for parallelization",
    )
    demux.add_argument(
        "--chunksize", type=float, default=int(1e7),
        help=(
            "Number of reads to process in between writing to disk. "
            "Larger values process faster, but use more RAM.")
    )
    demux.add_argument(
        "--i7", action="store_true",
        help="Demultiplex on i7 index instead of inline barcodes."
    )
    demux.add_argument(
        "--merge-technical-replicates", action="store_true",
        help=(
            "If a sample name appears >1 time in barcodes file the default "
            "behavior is to append '-technical-replicate-x' to each sample "
            "name, unless this option is turned on, in which case the data "
            "will be merged into a single sample.")
    )
    demux.add_argument(
        "--logger", type=str, nargs="*", default=("INFO", None),
        help=(
            "Logging info entered as one value for LOGLEVEL, or two values "
            "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    )
    # TOO RISKY perhaps, make the user remove existing dir themselves?
    # demux.add_argument(
    #     "--force", "-f", action="store_true",
    #     help=(
    #         "Force overwrite. Allows overwriting an existing directory of "
    #         "demultiplexed fastq files. (Be careful).")
    # )


def setup_assemble_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad assemble` subcommand parser.

    """
    assemble = subparsers.add_parser(
        "assemble",
        description=HEADER + "\n" + " ipyrad assemble: run assembly steps",
        help="Assemble data from fastqs to orthologs (steps 1-7).",
        epilog=ASSEMBLE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    assemble.add_argument(
        "-j", metavar="json", type=Path, required=True,
        help="JSON assembly file path ({assembly_name}.json)."
    )
    assemble.add_argument(
        "-s", metavar="steps", type=str, required=True,
        help="Steps of the assembly to run, e.g., 1 or 23 or 1234567."
    )
    assemble.add_argument(
        "-c", metavar="cores", type=int, default=4,
        help="Number of cores for parallelization.",
    )
    assemble.add_argument(
        "-t", metavar="threads", type=int, default=2,
        help="Number of threads per-core for multi-threaded steps.",
    )
    assemble.add_argument(
        "--force", "-f", action="store_true",
        help="Force overwrite existing results of an assembly step.",
    )
    assemble.add_argument(
        "--logger", type=str, nargs="*",
        help=(
            "Logging info entered as one value for LOGLEVEL, or two values "
            "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    )


def setup_new_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad new` subcommand parser.

    """
    new = subparsers.add_parser(
        "new",
        description=HEADER + "\n" + " ipyrad new: create a new named JSON Assembly file",
        help="Create a new named JSON Assembly file.",
        epilog=NEW_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # new.add_argument(
    #     "-n", metavar="name", type=str, required=True,
    #     help="Name of Assembly to use as prefix for JSON and output files."
    # )
    new.add_argument(
        "--force", "-f", action="store_true",
        help="Force overwrite of existing JSON file.",
    )
    new.add_argument(
        "--logger", type=str, nargs="*", default=("INFO", None),
        help=(
            "Set logging level as a single value to log to STDOUT (e.g., INFO) "
            "or as two values to log to a file (e.g., INFO log.txt).")
    )
    params = new.add_argument_group("params", "args to set parameters")
    params.add_argument(
        "--assembly_name", "-n", type=str, required=True,
        help="Name prefix used to for JSON file and outputs.")
    params.add_argument(
        "--project_dir", type=Path,
        help="The filepath at which to store the JSON file and output directories.")
    params.add_argument(
        "--fastq_paths", type=Path, nargs="*", required=True,
        help="One or more filepaths containing one or more fastq files (e.g., ./data/*.gz).")
    params.add_argument(
        "--reference_sequence", type=Path,
        help="Path to a reference genome fasta file. If None then denovo assembly is performed.")
    params.add_argument(
        "--filter_adapters", type=int,
        help="Set quality filtering level: 0=minimal, 1=quality-only, 2=quality+trim. Default=2")
    params.add_argument(
        "--filter_min_trim_len", type=int,
        help="Set minimum length of trimmed reads. Default=35.")
    params.add_argument(
        "--max_low_qual_bases", type=int,
        help="Maximum number of Ns allowed in reads. Default=5.")
    params.add_argument(
        "--trim_reads", type=int, nargs=2,
        help="Argument to perform global trimming at start and end of reads. Default=0 0.")
    params.add_argument(
        "--phred_qscore_offset", type=int,
        help="Phred quality score offset. Default=33.")
    params.add_argument(
        "--reference_as_filter", type=Path,
        help="Path to a reference genome fasta file used as a filter to remove mapped reads as contaminants.")
    params.add_argument(
        "--clust_threshold", type=float,
        help="Clustering threshold used to identify sequence homology in denovo assembly. Default=0.85")
    params.add_argument(
        "--min_depth_statistical", type=int,
        help="Minimum depth at which to make site variant calls. Default=6")
    params.add_argument(
        "--min_depth_majrule", type=int,
        help="Minimum depth at which to make majority-rule variant calls, if lower than min_depth_statistical. Default=6")
    params.add_argument(
        "--max_depth", type=int,
        help="Maximum depth at which to make site variant calls. Default=inf")
    params.add_argument(
        "--max_alleles_consens", type=int,
        help="Maximum number of alleles within a sample above it is excluded. Default=2")
    params.add_argument(
        "--max_n_consens", type=float,
        help="Maximum number or proportion of uncalled (N) bases in a consensus allele. Default=0.05")
    params.add_argument(
        "--max_h_consens", type=float,
        help="Maximum number or proportion of heterozygous bases in a consensus allele. Default=0.05")
    params.add_argument(
        "--min_samples_locus", type=int,
        help="Minimum number of samples that must have data at a locus to retain it. Default=4.")
    params.add_argument(
        "--max_snps_locus", type=float,
        help="Maximum proportion of sites that are variable in a locus, above which it is excluded. Default=0.2.")
    params.add_argument(
        "--max_indels_locus", type=float,
        help="Maximum number of sites that are indels in a locus, above which it is excluded. Default=8.")
    params.add_argument(
        "--max_shared_h_locus", type=float,
        help="Maximum number of sites that are indels in a locus, above which it is excluded. Default=8.")
    # ...


def setup_print_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad stats` subcommand parser.

    """
    pprint = subparsers.add_parser(
        "print",
        description=HEADER + "\n" + " ipyrad print: show summary of Assembly parameters, stats or filepaths.",
        help="Show summary of Assembly params, stats or filepaths.",
        epilog=STATS_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    pprint.add_argument(
        "--json", "-j", metavar="json", type=Path, required=True,
        help="path to {{name}}.json Assembly file."
    )
    pprint.add_argument(
        "--stats", action="store_true",
        help="Print Assembly stats."
    )
    pprint.add_argument(
        "--params", action="store_true",
        help="Print Assembly parameter settings."
    )


def setup_branch_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad branch` subcommand parser.

    """
    branch = subparsers.add_parser(
        "branch",
        help="Create a new branch (copy) of an assembly (JSON).",
        epilog=BRANCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    branch.add_argument(
        "-p", metavar="params", type=Path, required=True,
        help="path to params-{{name}}.txt used to select the assembly to branch.",
    )
    branch.add_argument(
        "-n", metavar="name", type=str, required=True,
        help="name for new assembly, creates ./{{name}}.json and params-{{name}}.txt."
    )
    branch_subsample = branch.add_mutually_exclusive_group()
    branch_subsample.add_argument(
        "-s", metavar="subsample", type=str, nargs="*",
        help="Name of one or more samples to keep (cannot use --exclude at same time as this option)."
    )
    branch_subsample.add_argument(
        "-e", metavar="exclude", type=str, nargs="*",
        help="Name of one or more samples to exclude (cannot use --subsample at same time as this option)."
    )
    branch.add_argument(
        "--force", "-f", action="store_true",
        help="force overwrite even if new named branch JSON already exists."
    )
    branch.add_argument(
        "--logger", type=str, nargs="*", default=("INFO", None),
        help=(
            "Logging info entered as one value for LOGLEVEL, or two values "
            "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    )


def setup_merge_subparser(subparsers: argparse._SubParsersAction) -> argparse._SubParsersAction:
    """Add `ipyrad merge` subcommand parser.

    """
    merge = subparsers.add_parser(
        "merge",
        help="Merge >=2 assemblies into one assembly containing their samples together.",
        epilog=MERGE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    merge.add_argument(
        "-p", metavar="params", type=Path, required=True, nargs="+",
        help="path to params-{{name}}.txt used to select assemblies to merge.",
    )
    merge.add_argument(
        "-n", metavar="name", type=str, required=True,
        help="name for new assembly, creates ./{{name}}.json and params-{{name}}.txt."
    )
    merge.add_argument(
        "--force", "-f", action="store_true",
        help="force overwrite even if new named merge JSON already exists."
    )


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        prog="ipyrad",
        description=HEADER + "\n" + DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action='version', version=f"ipyrad {ip.__version__}")
    subparsers = parser.add_subparsers(help="sub-commands", dest="subcommand")

    # add subcommands
    setup_demux_subparser(subparsers)
    setup_new_subparser(subparsers)
    setup_assemble_subparser(subparsers)
    setup_print_subparser(subparsers)
    setup_branch_subparser(subparsers)
    setup_merge_subparser(subparsers)
    return parser


def test():
    parser = setup_parsers()

    # cmd = "demux -h"
    cmd = "demux -d ./a.txt -b ./b.txt -o ./fastqs --i7 -m 1"
    cmd_list = cmd.split()
    args = parser.parse_args()
    print(args)
    print(parser.parse_args(cmd_list))

    cmd = "assemble -p params.txt -s 12 -c 4 -t 2"
    cmd_list = cmd.split()
    print(parser.parse_args(cmd_list))

    cmd = "branch -p params.txt -n NAME -s X Y Z"
    cmd_list = cmd.split()
    print(parser.parse_args(cmd_list))

    cmd = "merge -p params1.txt params2.txt -n NAME"
    cmd_list = cmd.split()
    print(parser.parse_args(cmd_list))

    cmd = "-h"
    cmd_list = cmd.split()
    print(parser.parse_args(cmd_list))

    # parser.parse_args(["root", "--help"])


def main():
    """Parse user CLI args and perform actions."""
    parser = setup_parsers()
    args = parser.parse_args()
    # args, unknown = parser.parse_known_args()

    # set logging ---------------------------------------------------
    if hasattr(args, "logger") and args.logger:
        if len(args.logger) > 1:
            ip.set_log_level(args.logger[0], args.logger[1])
        else:
            ip.set_log_level(args.logger[0])

    # parse unknown args --------------------------------------------
    # unkown is parsed from any additional "--key value" args
    # and will be converted to {key: value} and checked for format
    # and that the keys must be param keys.
    # kwargs = {}
    # idx = 0
    # if unknown:
    #     while idx < len(unknown):
    #         # get the next arg
    #         key = unknown[idx]
    #         idx += 1

    #         # if it doesn't start with '--' its an error
    #         if not key.startswith("--"):
    #             logger.warning(f"CLI arg '{key}' not recognized, skipped.")
    #             idx += 1

    #         # it does start with '--', now it must be followed by a value
    #         else:
    #             try:
    #                 val = unknown[idx]
    #                 key = key.lstrip("--")
    #                 if key in PARAMS:
    #                     kwargs[key] = val
    #                 else:
    #                     raise IndexError()
    #                 logger.debug(f"parsed CLI arg to param '{key}': {val}")
    #             except IndexError:
    #                 logger.warning(f"CLI arg '{key}' not recognized, skipped.")
    #             idx += 1

    # demultiplexing job --------------------------------------------
    if args.subcommand == "demux":
        from ipyrad.demux.demux import Demux
        Demux(
            fastq_paths=args.d,
            barcodes_path=args.b,
            outpath=args.o,
            re1=args.re1,
            re2=args.re2,
            cores=args.cores,
            max_barcode_mismatch=args.max_mismatch,
            merge_technical_replicates=args.merge_technical_replicates,
            chunksize=args.chunksize,
            i7=args.i7,
        ).run()

    # new params generate -------------------------------------------
    if args.subcommand == "new":
        tmp = ip.Assembly(args.assembly_name)
        set_params = {}
        for key, val in vars(args).items():
            if val is not None:
                if key not in ("subcommand", "force", "logger", "assembly_name"):
                    set_params[key] = val
        logger.info(f"setting {1 + len(set_params)} parameters")
        for key, val in set_params.items():
            setattr(tmp.params, key, val)
        # tmp.write_params(force=args.force)
        if tmp.json_file.exists() and not args.force:
            logger.error(f"JSON file ({tmp.json_file}) exists. Use --force to overwrite.")
            raise SystemExit(1)
        else:
            tmp.save_json()
        raise SystemExit(0)

    # assembly job --------------------------------------------------
    if args.subcommand == "assemble":
        # get Assembly from JSON
        data = ip.load_json(args.j)
        data.run(args.s, force=args.force, cores=args.c, threads=args.t)
        raise SystemExit(0)

    # branch job
    if args.subcommand == "branch":
        # load from json and create new copy
        data1 = ip.load_json(args.json)
        data2 = data1.branch(args.assembly_name)

        # set new parameters on the copy
        # ...

        # subsample/exclude samples on the copy
        # ...

        # save new JSON
        # data2.save_json()

    # merge job
    if args.subcommand == "merge":
        pass

    # stats job
    if args.subcommand == "print":
        tmp = ip.load_json(args.json)
        if args.params:
            logger.info(f"Printing parameter settings for Assembly ({tmp.name}):\n{tmp.params}")
        if args.stats:
            logger.info(f"Printing stats summary for Assembly ({tmp.name}):\n{tmp.stats}")


if __name__ == "__main__":

    # main()
    test()
    # cmd = "ipyrad demux --type inline"
    # parser = setup_parsers()
    # print(parser)