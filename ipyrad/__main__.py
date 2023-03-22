#!/usr/bin/env python

"""ipyrad command line interface (CLI).

Examples
--------
>>> ipyrad -n test
>>> ipyrad -p params-test.txt -s 1 -r
>>> ipyrad -p params-test.txt -s 2 -c 8 -t 2
>>> ipyrad -p params-test.txt -s 3 --logger INFO ipyrad_log.txt
>>> ipyrad -p params-test.txt -s 45 -c 100 --MPI
>>> ipyrad -p params-test.txt -b newtest
>>> ipyrad -p params-newtest.txt -s 6 --ipcluster default
>>> ipyrad -m combined params-test.txt params-newtest.txt
>>> ipyrad -p params-combined.txt -s 7 -f

Other convenience functions
---------------------------
>>> ipyrad --download SRP021469 sra-fastqs/
"""

from typing import Dict
from pathlib import Path
import argparse
import sys
import os

from pkg_resources import get_distribution
import ipyparallel as ipp
from loguru import logger

import ipyrad as ip
from ipyrad.core.params_schema import ParamsSchema
from ipyrad.assemble.utils import IPyradError, IPyradExit
from ipyrad.core.cluster import get_num_cpus

##############################################################
##############################################################
#    GLOBALS                                                 #
##############################################################
##############################################################

LIST_TYPE_PARAMS = [
    "restriction_overhang",
    "trim_reads",
    "trim_loci",
    "output_formats",
]
VERSION = str(get_distribution('ipyrad')).split()[1]
HEADER = f"""
 -------------------------------------------------------------
  ipyrad [v.{VERSION}]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------\
"""

EPILOG = """
  * Example command-line usage:
    ipyrad -n data                       # create new file called params-data.txt
    ipyrad -p params-data.txt -s 123     # run only steps 1-3 of assembly.
    ipyrad -p params-data.txt -s 3 -f    # run step 3, overwrite existing data.

  * HPC parallelization across 32 cores
    ipyrad -p params-data.txt -s 3 -c 32 --MPI

  * Print results summary
    ipyrad -p params-data.txt -r

  * Branch/Merging Assemblies
    ipyrad -p params-data.txt -b newdata
    ipyrad -m newdata params-1.txt params-2.txt [params-3.txt, ...]

  * Subsample taxa during branching
    ipyrad -p params-data.txt -b newdata taxaKeepList.txt

  * Download sequence data from SRA into directory 'sra-fastqs/'
    ipyrad --download SRP021469 sra-fastqs/

  * Documentation: http://ipyrad.readthedocs.io
"""

REQUIRED_ACTION = """
    Must provide action argument along with -p argument for params file.
    e.g., ipyrad -p params-test.txt -r              # shows results
    e.g., ipyrad -p params-test.txt -s 12           # runs steps 1 & 2
    e.g., ipyrad -p params-test.txt -b newbranch    # branch this assembly
"""

PARSER = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=EPILOG,
)

# add arguments
PARSER.add_argument(
    '-v', '--version', action='version',
    version=str(get_distribution('ipyrad')))
PARSER.add_argument(
    '-r', "--results", action='store_true',
    help="show results summary for Assembly in params.txt and exit")
PARSER.add_argument(
    '-f', "--force", action='store_true',
    help="force overwrite of existing data")
PARSER.add_argument(
    '-q', "--quiet", action='store_true',
    help="do not print to stderror or stdout.")
PARSER.add_argument(
    '-n', dest="new", type=str, default=None,
    help="create new file 'params-{new}.txt' in current directory")
PARSER.add_argument(
    '-p', dest="params", type=str, default=None,
    help="path to params file for Assembly: params-{assembly_name}.txt")
PARSER.add_argument(
    '-s', dest="steps", type=str, default=None,
    help="Set of assembly steps to run, e.g., -s 123")
PARSER.add_argument(
    '-b', dest="branch", type=str, default=None, nargs="*",
    help="create new branch of Assembly as params-{branch}.txt, and " + \
    "can be used to drop samples from Assembly.")
PARSER.add_argument(
    '-m', dest="merge", default=None, nargs="*",
    help="merge multiple Assemblies into one joint Assembly, and " + \
    "can be used to merge Samples into one Sample.")
PARSER.add_argument(
    "-c", metavar="cores", dest="cores",
    type=int, default=0,
    help="number of CPU cores to use (Default=0=All)")
PARSER.add_argument(
    "-t", metavar="threading", dest="threads",
    type=int, default=2,
    help="tune threading of multi-threaded binaries (Default=2)")
PARSER.add_argument(
    "--logger", action="store_true",
    help="print info to a logfile in ./ipyrad_log.txt.")
PARSER.add_argument(
    "--ipcluster", dest="ipcluster",
    type=str, nargs="?", const="default",
    help="connect to running ipcluster, enter profile name (default='default'")
PARSER.add_argument(
    "--download", dest="download", type=str,
    nargs="*", default=None,  # const="default",
    help="download fastq files by accession (e.g., SRP or SRR)")


class CLI:
    """Command line organization."""
    def __init__(self):
        self.args = None
        """: store argparse object."""
        self.parsedict: Dict[str, str] = None
        """: store loaded params file data as a str dict."""
        self.data: ip.Assembly = None
        """: store loaded Assembly object."""

    def _parse_argparse_to_args(self) -> None:
        """Fill CLI args. If no args then print help message.

        This is called by parse_params_file()
        """
        if len(sys.argv) == 1:
            PARSER.print_help()
        self.args = PARSER.parse_args()

    def _check_required_args(self) -> None:
        """User must enter -p or -n as an argument"""
        nargs = [self.args.params, self.args.new, self.args.download, self.args.merge]
        if not any(nargs):
            raise IPyradExit("\n" + "\n".join([
                "  ipyrad command must include either -p or -n ",
                "  run 'ipyrad -h' for further command line instructions\n",
            ]))

    def _parse_params(self) -> None:
        """Fill parsedict with params from paramsfile."""
        # check that params.txt file is correctly formatted.
        if not self.args.params:
            raise IPyradError("\n  No params file found\n")
        if not Path(self.args.params).exists():
            raise IPyradError("\n  No params file found\n")
        with open(self.args.params, 'r', encoding="utf-8") as paramsin:
            lines = paramsin.readlines()

        # get values from the params file lines
        vals = [i.split("##")[0].strip() for i in lines[1:] if i.strip()]

        # get keys in order from a tmp assembly
        keys = [i[0] for i in ip.Assembly('null').params]
        self.parsedict = {str(i): j for (i, j) in zip(keys, vals)}

    def _enable_logger(self):
        """ set logger to debugging """
        info = [
            f"ipyrad: {ip.__version__}",
            f"system: {os.uname()}",
            f"env: {sys.prefix}",
            f"python: {sys.executable}",
            f"params: {self.data.params}"
        ]
        logger.debug(info)

    def _flagnew(self):
        """Write a new default (but named) params file and exit 0."""
        tmp = ip.Assembly(self.args.new)

        # write the new params file
        tmp.write_params(force=self.args.force)

        # print log to screen
        curdir = os.path.realpath(os.path.curdir)
        print(f"New file 'params-{self.args.new}.txt' created in {curdir}\n")

    def _flagparams(self):
        """If params but no action then raise 1, else set steps to empty."""
        if self.args.params:
            if not any([self.args.branch, self.args.results, self.args.steps]):
                raise IPyradExit(REQUIRED_ACTION)
            if not self.args.steps:
                self.args.steps = ""

    def _flagdownload(self):
        """Download data and exit 0. If SRA install required exit 1."""
        if len(self.args.download) == 1:
            downloaddir = "sra-fastqs"
        else:
            downloaddir = self.args.download[1]

        # call the analysis tool
        import ipyrad.analysis as ipa
        sra = ipa.sratools(
            accessions=self.args.download[0],
            workdir=downloaddir,
        )

        # get run info and print spacer after
        df = sra.fetch_runinfo((1, 4, 6, 29, 30))
        print("")

        # rename spots for prettier printing and send to screen
        df.rename(columns={"spots_with_mates": "mates"}, inplace=True)
        print(df)

        # run download with default name_fields and separator
        sra.run(
            name_fields=(30, 1),
            name_separator="_",
            force=self.args.force,
            auto=True,
            show_cluster=True,
        )
        print("")

    def merge_assemblies(self):
        """Merge 2 or more assemblies in a new assembly.

        merge all given assemblies into a new assembly. Copies the params
        from the first passed in extant assembly. this function is called
        with the ipyrad -m flag. You must pass it at least 3 values, the first
        is a new assembly name (a new `param-newname.txt` will be created).
        The second and third args must be params files for currently existing
        assemblies. Any args beyond the third must also be params file for
        extant assemblies.
        """
        print(HEADER)
        print(f"\n  Merging assemblies: {self.args.merge[1:]}")

        # Make sure there are the right number of args
        if len(self.args.merge) < 3:
            sys.exit(_WRONG_NUM_CLI_MERGE)

        # Make sure the first arg isn't a params file, someone could do it...
        newname = self.args.merge[0]
        if os.path.exists(newname) and ("params-" in newname):
            sys.exit(_WRONG_ORDER_CLI_MERGE)

        # Make sure first arg will create a param file that doesn't exist
        if os.path.exists("params-" + newname + ".txt") and (not self.args.force):
            sys.exit(_NAME_EXISTS_MERGE.format("params-" + newname + ".txt"))

        # Make sure the rest of the args are params files that already exist
        assemblies_to_merge = self.args.merge[1:]
        for assembly in assemblies_to_merge:
            if not os.path.exists(assembly):
                sys.exit(_DOES_NOT_EXIST_MERGE.format(assembly))

        # Get assemblies for each of the passed in params files.
        # We're recycling some of the machinery for loading assemblies here
        assemblies = []
        for params_file in self.args.merge[1:]:
            self.args.params = params_file
            self.parse_params()
            self.get_assembly()
            assemblies.append(self.data)

        # Do the merge and try to save it.
        merged_assembly = ip.merge(newname, assemblies)
        merged_assembly.write_params(force=self.args.force)
        merged_assembly.save_json()
        print("\n  Merging succeeded. New params file for merged assembly:")
        print("\n    params-{}.txt\n".format(newname))
        IPyradExit(0)

    def get_assembly(self) -> None:
        """Load Assembly from the <project>/<name>.json file.

        Then updates json settings with params from the user-friendly
        params file, and return the updated Assembly object.
        """
        # is this an old params (pre v.1.0) formatted file?


        # load from json
        project_dir = Path(self.parsedict['project_dir'])
        assembly_name = self.parsedict['assembly_name']
        json_file = (Path(project_dir) / assembly_name).with_suffix('.json')

        # Create new Assembly instead of loading if NEW
        if self.args.steps:
            # starting a new assembly
            if '1' in self.args.steps:
                if self.args.force:
                    data = ip.Assembly(assembly_name)
                else:
                    if json_file.exists():
                        raise IPyradExit(
                            "Assembly already exists, use force to overwrite")
                    data = ip.Assembly(assembly_name)
            else:
                data = ip.load_json(json_file)
        else:
            data = ip.load_json(json_file)

        # Update json assembly with params in paramsfile in case they changed
        paramsdict = data.params.dict()
        for key, val in self.parsedict.items():
            if val:
                if key in LIST_TYPE_PARAMS: # if key "," in val:
                    val = [i.strip() for i in val.split(",")]
                paramsdict[key] = val

        # this can convert str inputs to the proper type, unless the
        # type is a container, which the above part converts to List
        data.params = ParamsSchema(**paramsdict)
        self.data = data

    def show_stats(self) -> None:
        """Load assembly or dies, and print stats to screen."""
        # Be nice if the user includes the extension.
        project_dir = Path(self.parsedict['project_dir'])
        assembly_name = str(self.parsedict['assembly_name'])
        json_file = project_dir / assembly_name
        if not json_file.suffix == ".json":
            json_file = json_file.with_suffix(".json")

        if not json_file.exists():
            raise IPyradExit(f"Cannot find assembly {json_file}")

        # load the assembly
        data = ip.load_json(json_file)
        print(
            "\nSummary stats of Assembly {}".format(data.name) +
            "\n------------------------------------------------")

        if not data.stats.empty:
            print(data.stats)
            print(
                "\n\nFull stats files (TODO... refactoring.)" +
                "\n------------------------------------------------")

            fullcurdir = os.path.realpath(os.path.curdir)
            for i in range(1, 8):
                #enumerate(sorted(data.stats_files)):
                key = "s" + str(i)
                try:
                    val = data.stats_files[key]
                    val = val.replace(fullcurdir, ".")
                    print("step {}: {}".format(i, val))
                except (KeyError, AttributeError):
                    print("step {}: None".format(i))
            print("\n")
        else:
            print("No stats to display")

    def branch_assembly(self) -> None:
        """Create branch/copy of current assembly and write params file."""
        # ensure workdir exists so there is somewhere to write to.
        self.data.params.project_dir.mkdir(exist_ok=True)

        # get arguments to branch command
        bargs = self.args.branch

        # get new name and trim off any suffices
        newname = Path(bargs[0]).stem

        # look for subsample arguments
        if len(bargs) > 1:

            # parse str or file of names to include/drop
            subargs = bargs[1:]

            # is there a '-' indicating to drop
            remove = 0
            if subargs[0] == "-":
                remove = 1
                subargs = subargs[1:]

            # is sample list a file?
            if Path(subargs[0]).exists():
                with open(subargs[0], 'r', encoding="utf-8") as infile:
                    subsamples = [
                        i.split()[0] for i in infile.readlines() if i.strip()
                    ]
            else:
                subsamples = subargs

            # check subsample drop names
            fails = [i for i in subsamples if i not in self.data.samples.keys()]
            if any(fails):
                raise IPyradError(
                    "\n  Failed: unrecognized names, check spelling:\n  {}"
                    .format("\n  ".join([str(i) for i in fails]))
                )

            # if drop then get subtract list
            if remove:
                print(f"  dropping {len(subsamples)} samples")
                subsamples = list(set(self.data.samples.keys()) - set(subsamples))

            # If the arg after the new param name is a file that exists
            new_data = self.data.branch(newname, subsamples)

        # keeping all samples
        else:
            new_data = self.data.branch(newname, None)

        print("  creating a new branch called '{}' with {} Samples"
              .format(new_data.name, len(new_data.samples)))

        print("  writing new params file to {}"
              .format("params-" + new_data.name + ".txt\n"))
        new_data.write_params(force=self.args.force)
        new_data.save_json()

    def parse_params_file(self) -> None:
        """Parse and store params."""
        self._parse_argparse_to_args()
        self._check_required_args()

        # if args.debug turn on the debugger
        if self.args.logger:
            ip.set_log_level("DEBUG", log_file="./ipyrad_log.txt")

        # run flags that are not step/run commands: -n, -m, --download
        # if run, these all end with a sys.exit
        if self.args.new:
            self._flagnew()
            raise IPyradExit(0)

        # check for download argument
        if self.args.download:
            self._flagdownload()
            raise IPyradExit(0)

        # check for merge of branches
        if self.args.merge:
            self.merge_assemblies()
            raise IPyradExit(0)

        # check that -p is accompanied by an action (-s, -r, -b)
        self._flagparams()

        # fill parsedict with params from paramsfile
        self._parse_params()

        # check for branching; needs to occur after Assembly has been loaded.
        if self.args.branch:
            self.get_assembly()
            self.branch_assembly()
            raise IPyradExit(0)

    def run(self):
        """Main function to connect a cluster and run assembly steps"""
        self.parse_params_file()
        self.get_assembly()

        # run assembly steps
        if self.args.steps:
            if not self.args.quiet:
                print(HEADER)

            # set CLI ipcluster terms
            self.data.ipcluster["cores"] = (
                self.args.cores if self.args.cores else get_num_cpus())
            self.data.ipcluster["threads"] = self.args.threads

            # if user entered information to connect to an already
            # running ipyclient then connect to it. TODO!
            # This usage is more likely to running MPI, probably.
            if self.args.ipcluster:
                ipyclient = ipp.Client(profile=self.args.ipcluster)
                self.data.ipcluster["cores"] = len(ipyclient)
                self.data.run(
                    self.args.steps,
                    force=self.args.force,
                    quiet=self.args.quiet,
                    ipyclient=ipyclient,
                )

            # else, start one here. This is the most common usage.
            else:
                # ncores = self.data.ipcluster.get("cores", get_num_cpus())
                self.data.run(
                    self.args.steps,
                    force=self.args.force,
                    quiet=self.args.quiet,
                    cores=self.data.ipcluster["cores"],
                )

        # show results summary
        if self.args.results:
            self.show_stats()


def cli():
    """Command line utility function."""
    CLI().run()


_WRONG_NUM_CLI_MERGE = """
  Error: Attempting to merge assemblies but wrong number of args.
  The format for the merging command is:

  ipyrad -m new_assembly_name params-1.txt params-2.txt
"""

_STEP_1_ASSEMBLY_EXISTS = """
    Looks like you're trying to re-run step 1 on an Assembly that already
    exists. If you created a branch after step 1 you can simply proceed to step
    2. If you wish to re-run step 1 you may use force to overwrite.
"""

_WRONG_ORDER_CLI_MERGE = """
  Looks like you're trying to pass in a params file that
  already exists. The first arg for -m should be the name of
  the new assembly (this will create a new params file for you)
"""

_NAME_EXISTS_MERGE = """
   New Assembly name already exists: {}
   The merge Assembly name must be unique
"""

_DOES_NOT_EXIST_MERGE = """
  All arguments after the first for -m need to be params files
  for assemblies with Samples that are all in state 1. This
  Assembly param file does not exist: {}
"""


if __name__ == "__main__":
    cli()
