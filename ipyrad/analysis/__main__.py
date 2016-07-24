#!/usr/bin/env python2
""" the main CLI for calling svd4tet """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.parallel import ipcontroller_init
from ipyrad.assemble.util import IPyradWarningExit
import pkg_resources
import ipyrad as ip
import numpy as np
import argparse
import logging
import sys
import os
import atexit

import ipyrad.analysis as ipa

# pylint: disable=W0212
# pylint: disable=C0301

LOGGER = logging.getLogger(__name__)


def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
  * Example command-line usage ---------------------------------------------- 

  * Read in sequence/SNP data file, provide linkage, and output name. 
     svd4tet -s data.phy                          ## use full sequence data 
     svd4tet -s data.snps.phy -n test2            ## use SNPs and name test2
     svd4tet -s data.snps.phy -l data.snps.map    ## use one SNP from each locus

  * Load saved/checkpointed analysis from json file, or force restart. 
     svd4tet -j test.json                    ## reads and writes with name 'test'
     svd4tet -j test.json -f                 ## reads in test, forces overwrite

  * Sampling modes: 'equal' can infer or accept a guide tree to sample more
    efficiently; quartets can be sampled randomly or default is to sample all. 
     svd4tet -s data.snps -m all                         ## sample all quartets
     svd4tet -s data.snps -m random -q 1e6 -x 123        ## random 1M randomly
     svd4tet -s data.snps -m equal -q 1e6 -t guide.tre   ## sample 1M across tree

  * HPC parallelization
      svd4tet -s data.phy              ## uses all cores on one machine
      svd4tet -s data.phy -c 20 --MPI  ## access 20 cores across multiple nodes

  * Documentation: http://ipyrad.readthedocs.org/en/latest/
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('ipyrad')))

    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")

    #parser.add_argument('-q', "--quiet", action='store_true',
    #    help="do not print to stderror or stdout.")

    parser.add_argument('-s', metavar="seq", dest="seq",
        type=str, default=None,
        help="path to input phylip file (SNPs of full sequence file)")

    parser.add_argument('-j', metavar='json', dest="json",
        type=str, default=None,
        help="load checkpointed/saved analysis from JSON file.")

    parser.add_argument('-m', metavar="method", dest="method",
        type=str, default="all",
        help="method for sampling quartets (all, random, or equal)")

    parser.add_argument('-q', metavar="nquartets", dest="nquartets",
        type=int, default=0,
        help="number of quartets to sample (if not -m all)")

    parser.add_argument('-b', metavar="boots", dest="boots",
        type=int, default=0,
        help="number of non-parametric bootstrap replicates")

    parser.add_argument('-l', metavar="map_file", dest="map",
        type=str, default=None,
        help="map file of snp linkages (e.g., ipyrad .snps.map)")

    parser.add_argument('-r', metavar="resolve", dest='resolve', 
        type=int, default=1, 
        help="randomly resolve heterozygous sites (default=1)")

    parser.add_argument('-n', metavar="name", dest="name",
        type=str, default="svd",
        help="output name prefix (default: 'svd')")

    parser.add_argument('-o', metavar="outdir", dest="outdir",
        type=str, default="./analysis_quartets",
        help="output directory (default: creates ./analysis_quartets)")

    parser.add_argument('-t', metavar="starting_tree", dest="tree",
        type=str, default=None,
        help="newick file starting tree for equal splits sampling")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=0,
        help="number of CPU cores to use (default = 0 = Use all)")

    parser.add_argument("-x", metavar="random_seed", dest="rseed",
        type=int, default=None,
        help="random seed for quartet sampling and/or bootstrapping")    

    parser.add_argument('-d', "--debug", action='store_true',
        help="print lots more info to ipyrad_log.txt.")

    parser.add_argument("--MPI", action='store_true',
        help="connect to parallel CPUs across multiple nodes")


    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    ## RAISE errors right away for some bad argument combinations:
    ## if 'random' require nquarts argument
    if args.method == 'random':
        if not args.nquartets:
            raise IPyradWarningExit(\
            "  Number of quartets (-n) is required with method = random")

    ## if 'equal' method require starting tree and nquarts
    if args.method == 'equal':
        if not args.nquartets:
            raise IPyradWarningExit(\
            "  Number of quartets (-n) is required with method = random")
        if not args.tree:
            raise IPyradWarningExit(\
        "  A starting tree (-t newick file) is required with method = equal")

    if not any(x in ["seq", "json"] for x in vars(args).keys()):
        print("""
    Bad arguments: svd4tet command must include at least one of (-s or -j) 
    """)
        parser.print_help()
        sys.exit(1)

    return args




def main():
    """ main function """
    ## not in ipython
    ip.__interactive__ = 0

    header = \
    "\n ---------------------------------------------------------------------"+\
    "\n  Analysis tools for ipyrad [v.{}]".format(ip.__version__)+\
    "\n  svd4tet -- fast quartet and tree inference "+\
    "\n ---------------------------------------------------------------------"
    print(header)

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## set random seed
    np.random.seed(args.rseed)

    ## debugger----------------------------------------
    if os.path.exists(ip.__debugflag__):
        os.remove(ip.__debugflag__)
    if args.debug:
        print("\n  ** Enabling debug mode ** ")
        ip.debug_on()
        atexit.register(ip.debug_off)      

    ## if JSON, load existing Quartet analysis -----------------------
    if args.json:
        data = ipa.svd4tet.load_json(args.json)

    ## else create a new tmp assembly for the seqarray-----------------
    else:
        ## create new Quartet class Object if it doesn't exist
        newjson = os.path.join(args.outdir, args.name+'.svd.json')
        if (not os.path.exists(newjson)) or args.force:
            data = ipa.svd4tet.Quartet(args.name, args.outdir, args.method)

        else:
            sys.exit("""
    Error: svd4tet analysis '{}' already exists in {} 
    Use the force argument (-f) to overwrite old analysis files, or,
    Use the JSON argument (-j {}/{}.svd.json) 
    to continue analysis of '{}' from last checkpoint.
    """.format(args.name, args.outdir, args.outdir, args.name, args.name))

        ## if more arguments for method 
        if args.nquartets:
            data.nquartets = int(args.nquartets)
        if args.boots:
            data.nboots = int(args.boots)

        ## store input files
        if args.map: 
            data.files.mapfile = args.map
        if args.tree:
            data.files.treefile = args.tree

        ## clear any existing hdf5 databases
        oldfiles = [data.h5in, data.h5out, data.files.qdump] + data.trees.values()
        for oldfile in oldfiles:
            if oldfile:
                if os.path.exists(oldfile):
                    os.remove(oldfile)

        ## get seqfile and names from seqfile
        data.files.seqfile = args.seq
        data.resolve_ambigs = args.resolve
        data.init_seqarray()
        data.parse_names()


    ## Use ipcluster info passed to command-line this time
    data._ipcluster["cores"] = args.cores
    if args.cores:
        data.cpus = int(args.cores)
    if args.MPI:
        data._ipcluster["engines"] = "MPI"
    else:
        data._ipcluster["engines"] = "Local"

    ## launch ipcluster and register for later destruction.
    data = ipcontroller_init(data)

    ## run svd4tet main function within a wrapper. The wrapper creates an 
    ## ipyclient view and appends to the list of arguments to run 'run'. 
    data.run(force=args.force)



if __name__ == "__main__": 
    main()


