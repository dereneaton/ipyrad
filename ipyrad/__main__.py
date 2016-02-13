#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.parallel import ipcontroller_init
from ipyrad.assemble.util import IPyradError
import pkg_resources
import ipyrad as ip
import argparse
import logging as _logging
import sys
import os

# pylint: disable=W0212

LOGGER = _logging.getLogger(__name__)


def parse_params(args):
    """ Parse the params file args, create and return Assembly object."""
    ## check that params.txt file is correctly formatted.
    try:
        with open(args.params) as paramsin:
            plines = paramsin.readlines()
    except IOError as _:
        sys.exit("No params file found")

    ## check header: big version changes can be distinguished by the header
    assert len(plines[0].split()[0]) == 6, \
    "params file is not compatible with ipyrad v.{}.".format(ip.__version__) \
    + "Create a new params file with: ipyrad -n"

    ## check length
    ## assert len(plines) > 30, "params file error. Format not recognized."

    ## make into a dict
    items = [i.split("##")[0].strip() for i in plines[1:]]

    #keys = [i.split("]")[-2][-1] for i in plines[1:]]
    keys = range(len(plines)-1)
    parsedict = {str(i):j for i, j in zip(keys, items)}

    return parsedict



def showstats(parsedict):
    """ loads assembly or dies, and print stats to screen """

    project_dir = parsedict['1']
    ## Be nice if somebody also puts in the file extension
    assembly_name = parsedict['0']  ##.split(".assembly")[0]
    my_assembly = os.path.join(project_dir, assembly_name)

    ## If the project_dir doesn't exist don't even bother trying harder.
    if not os.path.isdir(project_dir):
        msg = """
    Trying to print stats for a project dir ({}) that doesn't exist. You must 
    run steps before you can show stats.
    """.format(project_dir)
        sys.exit(msg)

    if not assembly_name:
        msg = """
    Assembly name is not set in params.txt, meaning it was either changed or
    erased since the Assembly was started. Please restore the original name. 
    You can find the name of your Assembly in the "project dir": {}.
    """.format(project_dir)
        raise ip.assemble.util.IPyradParamsError(msg)

    try:
        print(my_assembly)
        if os.path.exists(my_assembly+".assembly"):
            data = ip.load.load_assembly(my_assembly, quiet=True)
        else:
            data = ip.load_json(my_assembly, quiet=True)

        print("Summary stats of Assembly {}".format(data.name) \
             +"\n------------------------------------------------")
        if not data.stats.empty:
            print(data.stats)
            print("\nFull stats files\n---------------------")

            fullcurdir = os.path.realpath(os.path.curdir)
            for key in sorted(data.statsfiles):
                val = data.statsfiles[key]
                val = val.replace(fullcurdir, ".")                
                print(key+":", val)
                print("\n----")
        else:
            print("No stats to display")

    except AssertionError as _:
        raise IPyradError("""
    Error: No Assembly file found at {}. 
    Check parameter settings for [project_dir]/[assembly_name]
    """).format(my_assembly)



def branch_assembly(args, parsedict):
    """ Load the passed in assembly and create a branch. Copy it
        to a new assembly, and also write out the appropriate params.txt
    """
    ## Get the current assembly
    data = getassembly(args, parsedict)
    new_data = data.branch(args.branch)

    print("  Creating a branch of assembly {} called {}".\
        format(data.name, new_data.name))

    print("  Writing new params file to {}"\
          .format("params-"+new_data.name+".txt"))
    new_data.write_params("params-"+new_data.name+".txt")



def getassembly(args, parsedict):
    """ loads assembly or creates a new one and set its params from 
    parsedict. Does not launch ipcluster. 
    """

    ## Creating an assembly with a full path in the name will "work"
    ## but it is potentially dangerous, so here we have assembly_name
    ## and assembly_file, name is used for creating new in cwd, file is
    ## used for loading existing.
    ##
    ## Be nice if the user includes the extension.
    project_dir = ip.core.assembly.expander(parsedict['1'])
    assembly_name = parsedict['0'] #.split(".assembly")[0]
    assembly_file = os.path.join(project_dir, assembly_name)

    ## Assembly creation will handle error checking  on
    ## the format of the assembly_name

    ## make sure the working directory exists.
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    ## Get cwd so we can pop back out after creating the new assembly
    ## We have to push and pop cwd so the assembly object gets created
    ## inside the working directory.
    cwd = os.path.realpath(os.path.curdir)
    os.chdir(project_dir)

    ## if forcing and doing step 1 then do not load existing Assembly
    if '1' in args.steps:
        ## create a new assembly object
        data = ip.Assembly(assembly_name)

    else:
        ## go back to cwd since existing will be loaded from its full path
        os.chdir(cwd)

        ## try loading an existing one
        try:
            if os.path.exists(assembly_file+".assembly"):
                ## json file takes precedence if both exist
                if os.path.exists(assembly_file+".json"):
                    data = ip.load_json(assembly_file)
                else:
                    data = ip.load.load_assembly(assembly_file)                    
            else:
                data = ip.load_json(assembly_file)

        ## if not found then create a new one
        except AssertionError:
            LOGGER.info("No current assembly found.")
            print("  No assembly found at: {}".format(assembly_file))

    ## ensure we are back where we belong in original cur dir
    os.chdir(cwd)

    ## for entering some params...
    for param in parsedict:
        ## trap assignment of assembly_name since it is immutable.
        if param == str(0):
            ## only pass to set_params if user tried to change assembly_name
            ## it will raise an Exit error
            if parsedict[param] != data.name:
                data.set_params(param, parsedict[param])
        else:
            ## all other params should be handled by set_params
            data.set_params(param, parsedict[param])

    return data



def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    ipyrad -n                            ## create new params.txt file.
    ipyrad -p params.txt                 ## run ipyrad with settings in params.txt.
    ipyrad -p params.txt -s 123          ## run only steps 1, 2 and 3 of assembly.
    ipyrad -p params.txt -s 4567         ## run steps 4, 5, 6 and 7 of assembly.
    ipyrad -p params.txt -s 3 -f         ## run step 3, overwrite existing data.

  * Preview mode:
    ipyrad -p params.txt -s --preview    ## Run quickly in preview mode to test

  * HPC parallelization
    ipyrad -p params.txt -s 3 --MPI      ## access cores across multiple nodes

  * Results summary quick view
    ipyrad -p params.txt -r              ## print summary stats to screen for params

  * Run quietly
    ipyrad -p params.txt -q              ## Run quietly, do not write output to std out

  * Branch assembly
    ipyrad -p params.txt -b new-params   ## Branch assembly for testing different
                                         ## parameter values w/o having to redo steps

  * Documentation: http://ipyrad.readthedocs.org/en/latest/
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('ipyrad')))

    parser.add_argument('-r', "--results", action='store_true',
        help="show summary of results for Assembly in params.txt and exit")

    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")

    #parser.add_argument('-q', "--quiet", action='store_true',
    #    help="do not print to stderror and stdout.")

    parser.add_argument('-n', metavar='new', dest="new", type=str, 
        default=None, 
        help="create params file as 'new-params.txt' in current directory")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str, default=None,
        help="path to params.txt file")

    parser.add_argument('-b', metavar='branch', dest="branch",
        type=str, default=None,
        help="create new branch of assembly designated by the -p flag")

    parser.add_argument('-s', metavar="steps", dest="steps",
        type=str, default="1234567",
        help="subset of assembly steps to perform (Default=1234567)")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=4,
        help="number of CPU cores to use (Default=4)")

    parser.add_argument("--MPI", action='store_true',
        help="connect to parallel CPUs across multiple nodes using MPI")

    parser.add_argument("--preview", action='store_true',
        help="run ipyrad in preview mode. Subset the input file so it'll run"\
            + "quickly so you can verify everything is working")


    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    if not any([args.params, args.results, args.new]):
        parser.print_help()
        sys.exit(1)

    return args


def main():
    """ main function """
    ## turn off traceback for the CLI
    ip.__interactive__ = 0

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## create new paramsfile if -n
    if args.new:

        ## Create a tmp assembly and call write_params to write out
        ## default params.txt file
        try:
            tmpassembly = ip.core.assembly.Assembly(args.new, quiet=True)
            tmpassembly.write_params("params-{}.txt".format(args.new),
                                     force=args.force)
        except Exception as inst:
            print(inst)
            sys.exit(2)

        print("\n  New file `params-{}.txt` created in: {}\n"\
              .format(args.new, os.path.realpath(os.path.curdir)))

        sys.exit(2)

    ## if showing results or branching, do not do steps and do not print header
    if args.results or args.branch:
        args.steps = ""
        print("")
    else:
        header = \
    "\n --------------------------------------------------"+\
    "\n  ipyrad [v.{}]".format(ip.__version__)+\
    "\n  Interactive assembly and analysis of RADseq data"+\
    "\n --------------------------------------------------"
        print(header)
    
    ## create new Assembly or load existing Assembly, quit if args.results
    if args.params:
        parsedict = parse_params(args)

        if args.results:
            showstats(parsedict)

        elif args.branch:
            branch_assembly(args, parsedict)

        else:
            ## run Assembly steps
            if args.steps:

                ## launch or load assembly with custom profile/pid
                data = getassembly(args, parsedict)

                ## store ipcluster info
                data._ipcluster["cores"] = args.cores

                if args.MPI:
                    data._ipcluster["engines"] = "MPI"
                else:
                    data._ipcluster["engines"] = "Local"

                ## launch ipcluster and register for later destruction
                data = ipcontroller_init(data)

                ## set to print headers
                data._headers = 1

                ## run assembly steps
                steps = list(args.steps)
                data.run(steps=steps, force=args.force, preview=args.preview)



if __name__ == "__main__": 
    main()


