#!/usr/bin/env python2
""" the main CLI for calling svd4tet """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.parallel import ipcontroller_init
from ipyrad.assemble.util import IPyradWarningExit
import pkg_resources
import ipyrad as ip
import argparse
import logging
import sys
import os
import atexit

import ipyrad.analysis as ipa

# pylint: disable=W0212


LOGGER = logging.getLogger(__name__)


def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 

  * Read in sequence/snp data file (-s) and provide output name (-o).
  * Files should be either phylip formated (whether they contain full
  * data or SNPs) or they should be in qusnps format (from ipyrad).
    svd4tet -s data.phy -o test1     ## use full sequence data 
    svd4tet -s data.snps -o test2    ## use all SNPs
    svd4tet -s data.usnps -o test2   ## use one SNP from each locus

  * Load existing checkpointed analysis from json file
    svd4tet -j test.json             ## reads and writes with name 'test'
    svd4tet -j test.json -f          ## reads in test, force (-f) overwrites

  * Sampling modes (-m option) all, random, or equal. Equal requires guide tree.
    svd4tet -s data.snps -m all -o test           ## sample all quartets
    svd4tet -s data.snps -m random -n 1e6 -o test ## random sample 1M quartets
    svd4tet -s data.snps -m equal -n 500000       ## sample 500K quartets evenly
            -t guide.tre -o test_equal            ## across splits of guide tree

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

    parser.add_argument('-q', "--quiet", action='store_true',
        help="do not print to stderror or stdout.")

    parser.add_argument('-s', metavar="seq", dest="seq",
        type=str, default=None,
        help="path to input phylip file (full seqs or SNPs)")

    parser.add_argument('-j', metavar='json', dest="json",
        type=str, default=None,
        help="load existing saved analysis from json file.")

    parser.add_argument('-m', metavar="method", dest="method",
        type=str, default="all",
        help="method for sampling quartets (all, random, or equal)")

    parser.add_argument('-n', metavar="nquartets", dest="nquartets",
        type=int, default=0,
        help="number of quartets to sample (if not -m all)")

    parser.add_argument('-b', metavar="boots", dest="boots",
        type=int, default=0,
        help="number of non-parametric bootstrap replicates")

    parser.add_argument('-t', metavar="starting_tree", dest="tree",
        type=str, default=None,
        help="newick file starting tree for equal splits sampling")

    parser.add_argument('-o', metavar="output_prefix", dest="output",
        type=str, default=None,
        help="output name prefix")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=0,
        help="number of CPU cores to use (default = 0 = Use all)")

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
        print("  Bad arguments: svd4tet command must include at least one of"\
             +"`-s` or `-j` to begin\n")
        parser.print_help()
        sys.exit(1)

    return args



def main():
    """ main function """
    ## not in ipython
    ip.__interactive__ = 0

    header = \
    "\n --------------------------------------------------"+\
    "\n  Analysis tools for ipyrad [v.{}]".format(ip.__version__)+\
    "\n  svd4tet -- fast quartet and tree inference "+\
    "\n --------------------------------------------------"
    print(header)

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## debugger
    if os.path.exists(ip.__debugflag__):
        os.remove(ip.__debugflag__)
    if args.debug:
        print("\n  ** Enabling debug mode ** ")
        ip.debug_on()
        atexit.register(ip.debug_off)      

    ## if JSON, load it
    if args.json:
        data = ip.load_json(args.json)
        data.outfiles.svdinput = data.outfiles.svdinput

    ## else create a tmp assembly for the seqarray
    else:
        if not args.output:
            raise IPyradWarningExit("  -o output_prefix required")
        if not args.seq:
            raise IPyradWarningExit("  -s sequence file required")
        ## create new JSON (Assembly) object
        data = ip.Assembly(args.output, quiet=True)
        data.outfiles.svdinput = args.seq
        data.set_params(1, "./")

        ## parse samples from the sequence file
        names = []
        with iter(open(args.seq, 'r')) as infile:
            infile.next().strip().split()
            while 1:
                try:
                    names.append(infile.next().split()[0])
                except StopIteration:
                    break
        ## store as Samples in Assembly
        data.samples = {name:ip.Sample(name) for name in names}

    ## store ipcluster info
    data._ipcluster["cores"] = args.cores
    if args.cores:
        data.cpus = args.cores

    if args.MPI:
        data._ipcluster["engines"] = "MPI"
    else:
        data._ipcluster["engines"] = "Local"

    ## launch ipcluster and register for later destruction
    data = ipcontroller_init(data)

    ## run svd4tet
    args = [data, args.boots, args.method, args.nquartets, args.force]
    data._clientwrapper(ipa.svd4tet.run, args, 45)


if __name__ == "__main__": 
    main()


