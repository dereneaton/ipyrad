#!/usr/bin/env python

""" the main CLI for calling ipyrad """

from __future__ import print_function, division  

# from ipyrad.core.parallel import register_ipcluster
from ipyrad.assemble.utils import IPyradError, detect_cpus
from pkg_resources import get_distribution
from distutils.version import LooseVersion

import ipyparallel as ipp
import ipyrad as ip
import argparse
import requests
import sys
import os

from .core.Parallel import Parallel


class CLI:
    def __init__(self):

        # turn off interactive flag
        ip.__interactive__ = 0

        # the parser object will be used to fill args, parsedict, and load data
        self.args = None
        self.parsedict = None
        self.data = None
        self.parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=EPILOG)

        # check for a new version of ipyrad and warn user of upgrade
        #self.check_for_new_version()

        # parse the command line to get CLI args
        self._parse_command_line()
        self.args = self.parser.parse_args()

        # bail if no args for 'params' or 'new'
        self.check_args()

        # if args.debug turn on the debugger
        if self.args.debug:
            ip.__interactive__ = 1

        # run flags that are not step/run commands: -n, -m, --download
        # if run, these all end with a sys.exit
        if self.args.new:
            self._flagnew()
            sys.exit(0)

        # check for download argument
        if self.args.download:
            self._flagdownload()
            sys.exit(0)

        # check that -p is accompanied by an action (-s, -r, -m)
        self._flagparams()

        # check for merge of branches
        if self.args.merge:
            self.merge_assemblies()
            sys.exit(0)
        
        # fill parsedict with params from paramsfile
        self.parse_params()

        # functions below here involve an Assembly object.
        self.get_assembly()
        if self.args.branch:
            self.branch_assembly()
            sys.exit(0)

        self.run()



    def check_for_new_version(self):
        "Test if there's a newer version and nag the user to upgrade."

        try:
            # get response and parse it
            url = "https://anaconda.org/ipyrad/ipyrad"
            response = requests.get(url)
            htmldat = response.text.split("\n")

            # pull version from html
            curversion = (next(
                (x for x in htmldat if "subheader" in x), None)
                .split(">")[1]
                .split("<")[0])

            # check version against current
            if LooseVersion(ip.__version__) < LooseVersion(curversion):
                print(VERSION_UPDATE.format(curversion))
        
        # Let this fail silently
        except Exception:
            pass


    def check_args(self):
        if not any(i in ["params", "new"] for i in vars(self.args).keys()):
            print("""
        Bad arguments: ipyrad command must include at least one of the 
        following args: -p or -n. 
        """)
            self.parser.print_help()
            sys.exit(1)


    def parse_params(self):
        "Parse the params file args, create and return Assembly object."

        # check that params.txt file is correctly formatted.
        if not self.args.params:
            raise IPyradError("\n  No params file found\n")
        elif not os.path.exists(self.args.params):
            raise IPyradError("\n  No params file found\n")
        else:
            with open(self.args.params) as paramsin:
                lines = paramsin.readlines()

        # get values from the params file lines
        vals = [i.split("##")[0].strip() for i in lines[1:] if i.strip()]

        # get keys in order from a tmp assembly
        keys = [i[1:] for i in ip.Assembly('null', quiet=True).params._keys]
        
        # store as a dict
        self.parsedict = {str(i): j for (i, j) in zip(keys, vals)}


    def _parse_command_line(self):
        """ Parse CLI args."""

        # if no args then return help message
        if len(sys.argv) == 1:
            self.parser.print_help()

        ## add arguments 
        self.parser.add_argument(
            '-v', '--version', 
            action='version', 
            version=str(get_distribution('ipyrad'))
        )
        self.parser.add_argument('-r', "--results", action='store_true',
            help="show results summary for Assembly in params.txt and exit")

        self.parser.add_argument('-f', "--force", action='store_true',
            help="force overwrite of existing data")

        self.parser.add_argument('-q', "--quiet", action='store_true',
            help="do not print to stderror or stdout.")

        self.parser.add_argument('-d', "--debug", action='store_true',
            help="print lots more info to ipyrad_log.txt.")

        self.parser.add_argument('-n', dest="new", type=str, default=None, 
            help="create new file 'params-{new}.txt' in current directory")

        self.parser.add_argument('-p', dest="params", type=str, default=None,
            help="path to params file for Assembly: params-{assembly_name}.txt")

        self.parser.add_argument('-s', dest="steps", type=str, default=None,
            help="Set of assembly steps to run, e.g., -s 123")

        self.parser.add_argument('-b', dest="branch", type=str, default=None, 
            nargs="*",
            help="create new branch of Assembly as params-{branch}.txt, and " + \
            "can be used to drop samples from Assembly.")

        self.parser.add_argument('-m', dest="merge", default=None, nargs="*",
            help="merge multiple Assemblies into one joint Assembly, and " + \
            "can be used to merge Samples into one Sample.")

        self.parser.add_argument("-c", metavar="cores", dest="cores",
            type=int, default=0,
            help="number of CPU cores to use (Default=0=All)")

        self.parser.add_argument("-t", metavar="threading", dest="threads",
            type=int, default=2,
            help="tune threading of multi-threaded binaries (Default=2)")

        self.parser.add_argument("--MPI", action='store_true',
            help="connect to parallel CPUs across multiple nodes")

        self.parser.add_argument("--ipcluster", dest="ipcluster",
            type=str, nargs="?", const="default",
            help="connect to running ipcluster, enter profile name or profile='default'")

        self.parser.add_argument("--download", dest="download", type=str, 
        nargs="*", default=None,  # const="default",
            help="download fastq files by accession (e.g., SRP or SRR)")


    # def _hardlog_cli(self):

    #     # Log the current version. End run around the LOGGER
    #     # so it'll always print regardless of log level.
    #     with open(ip.__debugfile__, 'a') as logfile:
    #         logfile.write(HEADER)
    #         logfile.write("\n  Begin run: {}".format(
    #             time.strftime("%Y-%m-%d %H:%M")))
    #         logfile.write("\n  Using args {}".format(vars(self.args)))
    #         logfile.write("\n  Platform info: {}".format(os.uname()))


    # def _set_logger(self):

    #     ## Turn the debug output written to ipyrad_log.txt up to 11!
    #     ## Clean up the old one first, it's cleaner to do this here than
    #     ## at the end (exceptions, etc)
    #     #if os.path.exists(ip.__debugflag__):
    #     #    os.remove(ip.__debugflag__)

    #     if self.args.debug:
    #         print("  ** Enabling debug mode ** ")
    #         ip.set_logger_level("DEBUG")
    #         #ip._debug_on()
    #         #atexit.register(ip._debug_off)

    #     ## Only blank the log file if we're actually going to run a new
    #     ## assembly. This used to be in __init__, but had the side effect
    #     ## of occasionally blanking the log file in an undesirable fashion
    #     ## for instance if you run a long assembly and it crashes and
    #     ## then you run `-r` and it blanks the log, it's crazymaking.
    #     if os.path.exists(ip.__debugfile__):
    #         if os.path.getsize(ip.__debugfile__) > 50000000:
    #             with open(ip.__debugfile__, 'w') as clear:
    #                 clear.write("file reset")


    def _flagnew(self):
        # Create a tmp assembly, call write_params to make default params.txt
        tmpassembly = ip.Assembly(
            self.args.new, 
            quiet=True, 
            cli=True, 
            force=self.args.force,
        )

        # write the new params file
        tmpassembly.write_params(
            "params-{}.txt".format(self.args.new), 
            force=self.args.force,
        )

        # print log to screen
        print("\n  New file 'params-{}.txt' created in {}\n".format(
            self.args.new, 
            os.path.realpath(os.path.curdir))
        )


    def _flagparams(self):
        # if params then must provide action argument with it
        if self.args.params:
            if not any([self.args.branch, self.args.results, self.args.steps]):
                print("""
        Must provide action argument along with -p argument for params file. 
        e.g., ipyrad -p params-test.txt -r              # shows results
        e.g., ipyrad -p params-test.txt -s 12           # runs steps 1 & 2
        e.g., ipyrad -p params-test.txt -b newbranch    # branch this assembly
        """)
                sys.exit(1)

            if not self.args.steps:
                self.args.steps = ""


    def _flagdownload(self):
        ## if download data do it and then exit. Runs single core in CLI. 
        if len(self.args.download) == 1:
            downloaddir = "sra-fastqs"
        else:
            downloaddir = self.args.download[1]

        # call the analysis tool
        import ipyrad.analysis as ipa
        sra = ipa.sratools(
            accession=self.args.download[0],
            workdir=downloaddir,
            force=self.args.force,
        )

        ## get run info and print spacer after
        df = sra.fetch_runinfo((1, 4, 6, 29, 30))
        print("")

        ## rename spots for prettier printing and send to screen
        df.rename(columns={"spots_with_mates": "mates"}, inplace=True)
        print(df)

        ## run download with default name_fields and separator
        sra.run(
            name_fields=(30, 1), 
            name_separator="_", 
            force=self.args.force)
        print("")



    def merge_assemblies(self):
        """ 
        merge all given assemblies into a new assembly. Copies the params
        from the first passed in extant assembly. this function is called 
        with the ipyrad -m flag. You must pass it at least 3 values, the first
        is a new assembly name (a new `param-newname.txt` will be created).
        The second and third args must be params files for currently existing
        assemblies. Any args beyond the third must also be params file for
        extant assemblies.
        """
        print(HEADER)
        print("\n  Merging assemblies: {}".format(self.args.merge[1:]))

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

        # Do the merge
        merged_assembly = ip.merge(newname, assemblies)

        ## Write out the merged assembly params file and report success
        merged_assembly.write_params("params-{}.txt".format(newname), 
            force=self.args.force)

        print("\n  Merging succeeded. New params file for merged assembly:")
        print("\n    params-{}.txt\n".format(newname))
        sys.exit(1)


    def get_assembly(self):
        """ 
        loads assembly or creates a new one and set its params from 
        parsedict. Does not launch ipcluster. 
        """      
        # Be nice if the user includes the extension.
        project_dir = self.parsedict['project_dir']
        assembly_name = self.parsedict['assembly_name']
        json_file = os.path.join(project_dir, assembly_name)
        if not json_file.endswith(".json"):
            json_file += ".json"

        # make sure the working directory exists.
        if not os.path.exists(project_dir):
            os.mkdir(project_dir)

        # Create new Assembly instead of loading if NEW 
        if self.args.steps:
            # starting a new assembly
            if '1' in self.args.steps:
                if self.args.force:
                    data = ip.Assembly(assembly_name, cli=True)
                else:
                    if os.path.exists(json_file):
                        raise IPyradError(
                        "Assembly already exists, use force to overwrite")
                    else:
                        data = ip.Assembly(assembly_name, cli=True)
            else:
                data = ip.load_json(json_file, cli=True)
        else:
            data = ip.load_json(json_file, cli=True)

        # Update json assembly with params in paramsfile in case they changed
        for key, param in self.parsedict.items():
            if key not in ["assembly_name"]:
                data.set_params(key, param)

        # store it.
        self.data = data


    def show_stats(self):
        "loads assembly or dies, and print stats to screen"

        # Be nice if the user includes the extension.
        project_dir = self.parsedict['project_dir']
        assembly_name = self.parsedict['assembly_name']
        json_file = os.path.join(project_dir, assembly_name)
        if not json_file.endswith(".json"):
            json_file += ".json"

        if not os.path.exists(json_file):
            raise IPyradError(
                "Cannot find assembly {}".format(json_file))

        # load the assembly 
        data = ip.load_json(json_file, quiet=True, cli=True)
        print(
            "\nSummary stats of Assembly {}".format(data.name) +
            "\n------------------------------------------------")
        
        if not data.stats.empty:
            print(data.stats)
            print("\n\nFull stats files" + 
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


    def branch_assembly(self):
        """ 
        Load the passed in assembly and create a branch. Copy it
        to a new assembly, and also write out the appropriate params.txt
        """

        # get arguments to branch command
        bargs = self.args.branch

        # get new name
        newname = bargs[0]

        # trim .txt if it was accidentally added
        if newname.endswith(".txt"):
            newname = newname[:-4]

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
            if os.path.exists(subargs[0]):
                with open(subargs[0], 'r') as infile:
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
                        .format("\n  ".join([i for i in fails])))
            
            # if drop then get subtract list
            if remove:
                print("  dropping {} samples".format(len(subsamples)))
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

        new_data.write_params(
            "params-" + new_data.name + ".txt",
            force=self.args.force)


    def run(self):
        """ main function """

        if self.args.steps:
            # print header
            print(HEADER)

            # set to print headers
            self.data._headers = 1

            # set CLI ipcluster terms
            self.data.ipcluster["threads"] = self.args.threads

            # if ipyclient is running (and matched profile) then use that one
            if self.args.ipcluster:
                ipyclient = ipp.Client(profile=self.args.ipcluster)
                self.data.ipcluster["cores"] = len(ipyclient)

            # if not then we need to register and launch an ipcluster instance
            else:
                # set CLI ipcluster terms
                ipyclient = None
                self.data.ipcluster["cores"] = (
                    self.args.cores if self.args.cores else detect_cpus())
                self.data.ipcluster["engines"] = "Local"
                if self.args.MPI:
                    self.data.ipcluster["engines"] = "MPI"
                    if not self.args.cores:
                        raise IPyradError("must provide -c argument with --MPI")

            # get pool object
            pool = Parallel(
                tool=self.data, 
                rkwargs={"steps": self.args.steps, "force": self.args.force},
                ipyclient=ipyclient,
                show_cluster=True,
                auto=True,
                )
            pool.wrap_run()

            # # run assembly steps
            # steps = list(self.args.steps)
            # self.data.run(
            #     steps=steps, 
            #     force=self.args.force, 
            #     show_cluster=1, 
            #     ipyclient=ipyclient,
            #     )
        
        # show results summary                 
        if self.args.results:
            self.show_stats()
          


def sratools_download(SRP, workdir='SRA_fastqs', force=False):
    import ipyrad.analysis as ipa
    sra = ipa.sratools(accession=SRP, workdir=workdir)

    ## get run info and print spacer after
    df = sra.fetch_runinfo((1, 4, 6, 29, 30))
    print("")

    ## rename spots for prettier printing and send to screen
    df.rename(columns={"spots_with_mates": "mates"}, inplace=True)
    print(df)

    ## run download with default name_fields and separator
    sra.run(name_fields=(30, 1), name_separator="_", force=force)
    print("")




HEADER = """
 -------------------------------------------------------------
  ipyrad [v.{}]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------\
 """.format(ip.__version__)


EPILOG = """
  * Example command-line usage: 
    ipyrad -n data                       ## create new file called params-data.txt 
    ipyrad -p params-data.txt            ## run ipyrad with settings in params file
    ipyrad -p params-data.txt -s 123     ## run only steps 1-3 of assembly.
    ipyrad -p params-data.txt -s 3 -f    ## run step 3, overwrite existing data.

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

VERSION_UPDATE = """
  A new version of ipyrad is available (v.{}). To upgrade run:
  $ conda install -c ipyrad ipyrad
"""

_MPI_CORES_ERROR = """
  If using --MPI you must supply the number of cores with the -c argument.
"""

_WRONG_NUM_CLI_MERGE = """
  Error: Attempting to merge assemblies but wrong number of args.
  The format for the merging command is:
  
  ipyrad -m new_assembly_name params-1.txt params-2.txt
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



def main():
    CLI()


if __name__ == "__main__": 
    main()
