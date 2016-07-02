#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.parallel import ipcontroller_init
from ipyrad.assemble.util import IPyradError, IPyradWarningExit, detect_cpus
import pkg_resources
import ipyrad as ip
import argparse
import logging
import atexit
import sys
import os

# pylint: disable=W0212

LOGGER = logging.getLogger(__name__)


def parse_params(args):
    """ Parse the params file args, create and return Assembly object."""
    ## check that params.txt file is correctly formatted.
    try:
        with open(args.params) as paramsin:
            plines = paramsin.readlines()
    except IOError as _:
        sys.exit("  No params file found")

    ## check header: big version changes can be distinguished by the header
    if not len(plines[0].split()[0]) == 7:
        raise IPyradWarningExit("""
    Error: file '{}' is not compatible with ipyrad v.{}.
    Please create an updated params file with the -n argument. 
    """.format(args.params, ip.__version__))

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
    assembly_name = parsedict['0']
    my_assembly = os.path.join(project_dir, assembly_name)

    ## If the project_dir doesn't exist don't even bother trying harder.
    if not os.path.isdir(project_dir):
        msg = """
    Trying to print stats for Assembly ({}) that doesn't exist. You must 
    first run steps before you can show results.
    """.format(project_dir)
        sys.exit(msg)

    if not assembly_name:
        msg = """
    Assembly name is not set in params.txt, meaning it was either changed or
    erased since the Assembly was started. Please restore the original name. 
    You can find the name of your Assembly in the "project dir": {}.
    """.format(project_dir)
        raise IPyradError(msg)

    data = ip.load_json(my_assembly, quiet=True)

    print("\nSummary stats of Assembly {}".format(data.name) \
         +"\n------------------------------------------------")
    
    if not data.stats.empty:
        print(data.stats)
        print("\n\nFull stats files"\
         +"\n------------------------------------------------")

        fullcurdir = os.path.realpath(os.path.curdir)
        for i in range(1, 8):
            #enumerate(sorted(data.stats_files)):
            key = "s"+str(i)
            try:
                val = data.stats_files[key]
                val = val.replace(fullcurdir, ".")                
                print("step {}: {}".format(i, val))
            except (KeyError, AttributeError):
                print("step {}: None".format(i))
        print("\n")
    else:
        print("No stats to display")



def branch_assembly(args, parsedict):
    """ 
    Load the passed in assembly and create a branch. Copy it
    to a new assembly, and also write out the appropriate params.txt
    """
    ## Get the current assembly
    bargs = args.branch

    if len(bargs) > 1:
        subsamples = bargs[1:]
    else:
        subsamples = []

    data = getassembly(args, parsedict)

    ## If the arg after the new param name is a file that exists
    ## then we'll just use that instead
    if os.path.exists(bargs[1]):
        new_data = data.branch(bargs[0], infile=bargs[1])
    else:
        new_data = data.branch(bargs[0], subsamples)

    print("  Creating a new branch called '{}' with {} Samples".\
        format(new_data.name, len(new_data.samples)))

    print("  Writing new params file to {}"\
          .format("params-"+new_data.name+".txt"))
    new_data.write_params("params-"+new_data.name+".txt", force=args.force)



def merge_assemblies(args):
    """ merge all given assemblies into a new assembly. Copies the params
    from the first passed in extant assembly. this function is called 
    with the ipyrad -m flag. You must pass it at least 3 values, the first
    is a new assembly name (a new `param-newname.txt` will be created).
    The second and third args must be params files for currently existing
    assemblies. Any args beyond the third must also be params file for
    extant assemblies.
    """
    print("\n  Merging assemblies: {}".format(args.merge[1:]))

    ## Make sure there are the right number of args
    if len(args.merge) < 3:
        msg = "\n  Attempting to merge assemblies but wrong number of args. "\
            + "The format for the merging command is:"\
            + "\n\n    ipyrad -m new_assembly_name params-1.txt params-2.txt"
        sys.exit(msg)
    ## Make sure the first arg isn't a params file, i could see someone
    ## someone trying this
    newname = args.merge[0]
    if os.path.exists(newname) and "params-" in newname:
        msg = "\n  Looks like you're trying to pass in a params file that "\
            + "already exists. The first arg\n  for -m should be the name of "\
            + "the new assembly (this will create a new params file for you)."
        sys.exit(msg) 
    ## Make sure first arg will create a param file that doesn't already exist
    if os.path.exists("params-" + newname + ".txt") and not args.force:
        msg = "\n  First argument for ipyrad -m should be the name of the new"\
            + " assembly. This will create\n  a new params file called params-"\
            + newname + ".txt, but this file already exists."
        sys.exit(msg)
    ## Make sure the rest of the args are params files that already exist
    assemblies_to_merge = args.merge[1:]
    for assembly in assemblies_to_merge:
        if not os.path.exists(assembly):
            msg = "\n  All arguments after the first for -m need to be params"\
                + " files for assemblies with\n  samples that are all in state 1."\
                + " This assembly param file does not exist: " + assembly
            sys.exit(msg)

    assemblies = []
    ## Get assemblies for each of the passed in params files.
    ## We're recycling some of the machinery for loading assemblies here
    for params_file in args.merge[1:]:
        args.params = params_file
        parsedict = parse_params(args)
        assemblies.append(getassembly(args, parsedict))

    ## Do the merge
    merged_assembly = ip.merge(newname, assemblies)

    ## Set the values for some params that don't make sense inside
    ## merged assemblies
    merged_names = ", ".join([x.name for x in assemblies])
    merged_assembly.paramsdict["raw_fastq_path"] = "Merged: " + merged_names
    merged_assembly.paramsdict["barcodes_path"] = "Merged: " + merged_names
    merged_assembly.paramsdict["sorted_fastq_path"] = "Merged: " + merged_names

    ## Write out the merged assembly params file and report success
    merged_assembly.write_params("params-{}.txt".format(newname), force=args.force)

    print("\n  Merging succeeded. New params file for merged assembly:")
    print("\n    params-{}.txt\n".format(newname))


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
    assembly_name = parsedict['0']
    assembly_file = os.path.join(project_dir, assembly_name)

    ## Assembly creation will handle error checking  on
    ## the format of the assembly_name

    ## make sure the working directory exists.
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    try:
        ## If 1 and force then go ahead and create a new assembly
        if ('1' in args.steps) and args.force:
            data = ip.Assembly(assembly_name)
        else:
            data = ip.load_json(assembly_file)

    except IPyradWarningExit as inst:
        ## if no assembly is found then go ahead and make one
        if '1' not in args.steps:
            raise IPyradWarningExit("""
    Error: Steps >1 ({}) requested but no current assembly found - {}
    """.format(args.steps, assembly_file))
        else:
            ## create a new assembly object
            data = ip.Assembly(assembly_name)

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
    """ Parse CLI args."""

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    ipyrad -n data                       ## create new file called params-data.txt 
    ipyrad -p params-data.txt            ## run ipyrad with settings in params file
    ipyrad -p params-data.txt -s 123     ## run only steps 1-3 of assembly.
    ipyrad -p params-data.txt -s 3 -f    ## run step 3, overwrite existing data.

  * HPC parallelization across multiple nodes
    ipyrad -p params-data.txt -s 3 --MPI    

  * Results summary quick view
    ipyrad -p params-data.txt -r 

  * Branch/Merging Assemblies
    ipyrad -p params-data.txt -b newdata  
    ipyrad -m newdata params-1.txt params2.txt [params-3.txt, ...]

  * Preview mode (subsamples data)
    ipyrad -p params-data.txt -s --preview  

  * Documentation: http://ipyrad.readthedocs.io
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('ipyrad')))

    parser.add_argument('-r', "--results", action='store_true',
        help="show results summary for Assembly in params.txt and exit")

    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")

    parser.add_argument('-q', "--quiet", action='store_true',
        help="do not print to stderror or stdout.")

    parser.add_argument('-d', "--debug", action='store_true',
        help="print lots more info to ipyrad_log.txt.")

    parser.add_argument('-n', metavar='new', dest="new", type=str, 
        default=None, 
        help="create new file 'params-{new}.txt' in current directory")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str, default=None,
        help="path to params file for Assembly: params-{assembly_name}.txt")

    parser.add_argument('-b', metavar='branch', dest="branch",
        type=str, default=None, nargs="*",
        help="create a new branch of the Assembly as params-{branch}.txt")

    parser.add_argument('-m', metavar='merge', dest="merge",
        default=None, nargs="*",
        help="merge all assemblies provided into a new assembly")

    parser.add_argument('-s', metavar="steps", dest="steps",
        type=str, default=None,
        help="Set of assembly steps to perform, e.g., -s 123 (Default=None)")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=0,
        help="number of CPU cores to use (Default=0=All)")

    parser.add_argument("--MPI", action='store_true',
        help="connect to parallel CPUs across multiple nodes")

    parser.add_argument("--preview", action='store_true',
        help="run ipyrad in preview mode. Subset the input file so it'll run"\
            + "quickly so you can verify everything is working")


    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    if not any(x in ["params", "new"] for x in vars(args).keys()):
        print("Bad arguments: ipyrad command must include at least one of"\
                +"`-p` or `-n`\n")
        parser.print_help()
        sys.exit(1)

    return args



def main():
    """ main function """
    ## turn off traceback for the CLI
    ip.__interactive__ = 0

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## Turn the debug output written to ipyrad_log.txt up to 11!
    ## Clean up the old one first, it's cleaner to do this here than
    ## at the end (exceptions, etc)
    if os.path.exists(ip.__debugflag__):
        os.remove(ip.__debugflag__)

    if args.debug:
        print("\n  ** Enabling debug mode ** ")
        ip.debug_on()
        atexit.register(ip.debug_off)        

    ## create new paramsfile if -n
    if args.new:
        ## Create a tmp assembly and call write_params to write out
        ## default params.txt file
        try:
            tmpassembly = ip.Assembly(args.new, quiet=True)
            tmpassembly.write_params("params-{}.txt".format(args.new), 
                                     force=args.force)
        except Exception as inst:
            print(inst)
            sys.exit(2)

        print("\n    New file 'params-{}.txt' created in {}\n".\
               format(args.new, os.path.realpath(os.path.curdir)))
        sys.exit(2)


    ## if params then must provide action argument with it
    if args.params:
        if not any([args.branch, args.results, args.steps]):
            print("""
    Must provide action argument along with -p argument for params file. 
    e.g., ipyrad -p params-test.txt -r              ## shows results
    e.g., ipyrad -p params-test.txt -s 12           ## runs steps 1 & 2
    e.g., ipyrad -p params-test.txt -b newbranch    ## branch this assembly
    """)
            sys.exit(2)

    if not args.params:
        if not any([args.branch, args.results, args.steps]):
            print("""
    Must provide params file for branching, doing steps or getting results.
    e.g., ipyrad -p params-test.txt -r              ## shows results
    e.g., ipyrad -p params-test.txt -s 12           ## runs steps 1 & 2
    e.g., ipyrad -p params-test.txt -b newbranch    ## branch this assembly
    """)

    ## if branching, or merging do not allow steps in same command
    ## print spacer
    if any([args.branch, args.merge]):        
        args.steps = ""    
        print("")    

    ## always print the header when doing steps
    header = \
    "\n --------------------------------------------------"+\
    "\n  ipyrad [v.{}]".format(ip.__version__)+\
    "\n  Interactive assembly and analysis of RADseq data"+\
    "\n --------------------------------------------------"

    ## if merging just do the merge and exit
    if args.merge:
        print(header)
        merge_assemblies(args)
        sys.exit(1)

    ## create new Assembly or load existing Assembly, quit if args.results
    elif args.params:
        parsedict = parse_params(args)

        if args.branch:
            branch_assembly(args, parsedict)

        elif args.steps:
            ## print header
            print(header)

            ## run Assembly steps
            ## launch or load assembly with custom profile/pid
            data = getassembly(args, parsedict)

            ## if cores was entered, limit cores to this number
            ## otherwise use all available cores. By default _ipcluster[cores] 
            ## is set to detect_cpus in Assembly.__init__)
            if args.cores:
                data.cpus = args.cores

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

        if args.results:
            showstats(parsedict)


if __name__ == "__main__": 
    main()
