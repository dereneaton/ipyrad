#!/usr/bin/env python

"""
Examples
--------
>>> ipyrad demux --fastqs RAW.fastq --barcodes BARS.csv --...

>>> ipyrad new -n NEW
>>> ipyrad assemble -p params-NEW.txt -s 123
>>> ipyrad branch -p params-NEW.txt -n NEW2 --clust-threshold 0.9 --...
>>>
"""

import argparse
from pathlib import Path
from pkg_resources import get_distribution
import ipyrad as ip


VERSION = str(get_distribution('ipyrad')).split()[1]
HEADER = f"""
-------------------------------------------------------------
 ipyrad [v.{VERSION}]
 Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""


DESCRIPTION = "ipyrad command line tool. Select a subcommand."
EPILOG = """\
Examples
--------
# demux: demultiplexing data to samples by index or barcode
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv -c 10 -o ./demux
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
>>> ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --logger DEBUG
>>> ipyrad demux -d RAWS1/*.fastq.gz RAWS2/*.fastqs.gz -b BARCODES.tsv

# new: create a new assembly file as params-{name}.txt
>>> ipyrad new -n NAME

# assemble: filtering, clustering/mapping, variant calling, and formatting
>>> ipyrad assemble -p params -s 1
>>> ipyrad assemble -p params -s 123 -c 40 -t 4
>>> ipyrad assemble -p params -s 456 -f --logger DEBUG

# branch: create a copy/branch of an assembly to rerun w/ diff params
>>> ipyrad branch -p params-NEW1.txt -n NEW2
>>> ipyrad assemble -p params-NEW2.txt -s 7 -f

# merge: create new assembly containing union of samples from >=2 others.
>>> ipyrad merge -p params-NEW1.txt params-NEW2.txt -n NEW3
>>> ipyrad assemble -p params-NEW3.txt -s 67 -f

# stats: print summary of stats for an assembly
>>> ipyrad stats -p params-new1.txt
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
>>> ipyrad assemble -p params-new.txt -s 1
>>> ipyrad assemble -p params-new.txt -s 123 -c 10 -t 2
>>> ipyrad assemble -p params-new.txt -s 56 -c 40 -f --logger INFO
"""

BRANCH_EPILOG = """\
Examples
--------
>>> ipyrad branch -p params-new1.txt -n new2
>>> ipyrad branch -p params-new1.txt -n new2 -s A B C D E
>>> ipyrad branch -p params-new1.txt -n new2 -e F G H
"""

MERGE_EPILOG = """\
Examples
--------
>>> ipyrad merge -p params-new1.txt -p params-new2.txt -n new3
"""

STATS_EPILOG = """\
Examples
--------
>>> ipyrad stats -p params-new1.txt
"""


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        prog="ipyrad",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action='version', version=f"toytree {ip.__version__}")
    subparsers = parser.add_subparsers(help="sub-commands", dest="subcommand")

    ###################################################################
    # DEMUX
    demux = subparsers.add_parser(
        "demux",
        help="demultiplex unsorted data by index or barcode",
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
        help="One or more paths to select fastq data file using regex (e.g., './data/*.fastq.gz'",
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
        "-m", metavar="max_mismatch", type=int, default=0,
        help=(
            "The maximum number of allowed mismatches between true and "
            "oberved barcodes (Default=0).")
    )
    demux.add_argument(
        "-c", metavar="cores", type=int, default=4,
        help="Number of cores for parallelization",
    )
    demux.add_argument(
        "--chunksize", type=int, default=int(1e7),
        help=(
            "Number of reads to process in between writing to disk. "
            "Larger values process faster, but use more RAM.")
    )
    demux.add_argument(
        "--i7", action="store_true",
        help="demultiplex on i7 index instead of inline barcodes."
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
    demux.add_argument(
        "--force", "-f", action="store_true",
        help=(
            "Force overwrite. Allows overwriting an existing directory of "
            "demultiplexed fastq files. (Be careful).")
    )

    ###################################################################
    ###################################################################
    # ASSEMBLE
    assemble = subparsers.add_parser(
        "assemble",
        help="assemble data from fastqs to orthologs (steps 1-7).",
        epilog=ASSEMBLE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    assemble.add_argument(
        "-p", metavar="params", type=Path, required=True,
        help="Params file path (params-{name}.txt) for an assembly."
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

    ###################################################################
    ###################################################################
    # BRANCH
    branch = subparsers.add_parser(
        "branch",
        help="create a new branch (copy) of an assembly (JSON).",
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

    ###################################################################
    ###################################################################
    # MERGE
    merge = subparsers.add_parser(
        "merge",
        help="merge >=2 assemblies into one assembly containing their samples together.",
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

    ###################################################################
    # STATS
    stats = subparsers.add_parser(
        "stats",
        help="print summary of statistics for an assembly.",
        epilog=STATS_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    stats.add_argument(
        "-p", metavar="params", type=Path, required=True, nargs="+",
        help="path to params-{{name}}.txt to who stats for.",
    )
    return parser


def test():
    parser = setup_parsers()

    cmd = "demux -d ./a.txt -b ./b.txt -o ./fastqs"
    cmd_list = cmd.split()
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
    parser = setup_parsers()
    # example = ["demux", "-d", "../../pedtest/small_tmp_R*.gz", "-b", "../../pedtest/barcodes-true-plate1.csv", "--logger", "DEBUG", "-m", "1", "-o", "/tmp/TEST"])
    args = parser.parse_args()

    if args.logger:
        if len(args.logger) > 1:
            ip.set_log_level(args.logger[0], args.logger[1])
        else:
            ip.set_log_level(args.logger[0])

    # demultiplexing job.
    if args.demux:
        from ipyrad.demux.demux import Demux
        Demux(
            fastq_paths=args.d,
            barcodes_path=args.b,
            outpath=args.o,
            re1=args.re1,
            re2=args.re2,
            cores=args.c,
            merge_technical_replicates=args.m,
            chunksize=args.chunksize,
        ).run()


if __name__ == "__main__":

    main()
    # test()
    # cmd = "ipyrad demux --type inline"
    # parser = setup_parsers()
    # print(parser)