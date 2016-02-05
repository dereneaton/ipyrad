#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.parallel import ipcontroller_init
from ipyrad.assemble.util import IPyradError
import pkg_resources
import ipyrad as ip
import argparse
import logging
import sys
import os
import socket
import ipyparallel as ipp
from collections import Counter

# pylint: disable=W0212

LOGGER = logging.getLogger(__name__)


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
    assembly_name = parsedict['0'].split(".assembly")[0]
    my_assembly = os.path.join(project_dir, assembly_name)

    ## If the project_dir doesn't exist don't even bother trying harder.
    if not os.path.isdir(project_dir):
        msg = "\n\nTrying to print stats for a project dir ({}) that doesn't exist.\n"\
            + "You must run steps before you can show stats.".format(project_dir)
        sys.exit(msg)

    if not assembly_name:
        msg = "\n\nAssembly name is not set in params.txt. This means somebody\n"\
            + "changed it or erased it, which is bad. Please restore the\n"\
            + "original name. You can find the name of your assembly in the\n"\
            + "project dir: {}.".format(project_dir)
        raise ip.assemble.util.IPyradParamsError(msg)

    try:
        data = ip.load.load_assembly(my_assembly,
                                    quiet=True, 
                                    launch=False)

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
        else:
            print("No stats to display")

    except AssertionError as _:
        raise IPyradError("""
    Error: No Assembly file found at {}. 
    Check parameter settings for [project_dir]/[assembly_name]
    """).format(my_assembly)



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
    assembly_name = parsedict['0'].split(".assembly")[0]
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

    ## if forcing or doing step 1 then do not load existing Assembly
    if args.force and '1' in args.steps:
        ## create a new assembly object
        data = ip.Assembly(assembly_name)
    else:
        ## try loading an existing one
        try:
            #print("Loading - {}".format(assembly_name))
            data = ip.load.load_assembly(assembly_name, launch=False)

        ## if not found then create a new one
        except AssertionError:
            LOGGER.info("No current assembly found, create new - {}".\
                        format(assembly_file))
            data = ip.Assembly(assembly_name)

    ## pop directory
    os.chdir(cwd)

    ## for entering some params...
    for param in parsedict:
        ## If the param isn't empty, and also if it isn't assembly_name
        ## trap assignment of assembly_name since it is immutable.
        if parsedict[param] and param != str(0):
            try:
                data.set_params(param, parsedict[param])
            except Exception as inst:
                print(inst)
                print("Bad parameter in the params file - param {} value {}".\
                      format(param, parsedict[param]))
                raise

    return data


def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    ipyrad -n                       ## create new params.txt file.
    ipyrad -p params.txt            ## run ipyrad with settings in params.txt.
    ipyrad -p params.txt -s 123     ## run only steps 1, 2 and 3 of assembly.
    ipyrad -p params.txt -s 4567    ## run steps 4, 5, 6 and 7 of assembly.
    ipyrad -p params.txt -s 3 -f    ## run step 3, overwrite existing data.

  * HPC parallelization options
    ipyrad -p params.txt -s 3 -c 64 --MPI   ## Access 64 cores using MPI.

  * Results summary quick view
    ipyrad -p params.txt -r         ## print summary stats to screen for params.

  * Documentation: http://ipyrad.readthedocs.org/en/latest/
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('ipyrad')))

    parser.add_argument('-r', "--results", action='store_true',
        help="show summary of results for Assembly in params.txt and exit")

    parser.add_argument('-n', "--new", action='store_true',
        help="create new default params.txt file in current directory")

    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")

    #parser.add_argument('-q', "--quiet", action='store_true',
    #    help="do not print to stderror and stdout.")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str, default=None,
        help="path to params.txt file")

    parser.add_argument('-s', metavar="steps", dest="steps",
        type=str, default="1234567",
        help="subset of assembly steps to perform. Default=1234567")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=4,
        help="number of CPU cores to use")

    parser.add_argument("--MPI", action='store_true',
        help="connect to parallel CPUs across multiple nodes using MPI")

    parser.add_argument("--preview", action='store_true',
        help="Run ipyrad in preview mode. Subset the input file so it'll run"\
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
            tmpassembly = ip.core.assembly.Assembly("my_new_assembly", quiet=True)
            tmpassembly.write_params("params.txt", force=args.force)
        except Exception as inst:
            print(inst)
            print("\nUse --force to overwrite\n")
            sys.exit(2)

        print("New file `params.txt` created in {}".\
               format(os.path.realpath(os.path.curdir)))

        sys.exit(2)

    ## if showing results, do not do any steps and do not print header
    if args.results:
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

        else:
            ## run Assembly steps
            if args.steps:
                ## launch ipcluster and register for later destruction
                if args.MPI:
                    controller = "MPI"
                else:
                    controller = "Local"

                ## launch or load assembly with custom profile/pid
                data = getassembly(args, parsedict)

                ## might want to use custom profiles instead of ...
                data._ipclusterid = ipcontroller_init(nproc=args.cores,
                                                      controller=controller)
                ## set to print headers
                data._headers = 1

                ## print the connection info
                # ipyclient = ipp.Client(cluster_id=data._ipclusterid)
                # dview = ipyclient[:]
                # ccx = Counter(dview.apply_sync(socket.gethostname))
                # for key in ccx:
                #     print("  Node: {}: {} cores".format(key, ccx[key]))
                #     print("")

                ## run assembly steps
                steps = list(args.steps)
                data.run(steps=steps, force=args.force, preview=args.preview)



if __name__ == "__main__": 
    main()


