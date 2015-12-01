#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.newparamsfile import write_params
import pkg_resources
import ipyrad as ip
import argparse
import logging
import sys
import os

LOGGER = logging.getLogger(__name__)


def parse_params(args):
    """ Parse the params file args, create and return Assembly object."""
    ## check that params.txt file is correctly formatted.
    with open(args.params) as paramsin:
        plines = paramsin.readlines()

    ## check header: big version changes can be distinguished by the header
    assert len(plines[0].split()[0]) == 6, \
    "params file is not compatible with ipyrad v.{}.".format(ip.__version__) \
    + "Create a new params file with: ipyrad -n"

    ## check length
    assert len(plines) > 30, "params file error. Format not recognized."

    ## make into a dict
    items = [i.split("##")[0].strip() for i in plines[1:]]
    keys = range(1, 30)
    parsedict = {str(i):j for i, j in zip(keys, items)}

    return parsedict



def showstats(parsedict):
    """ loads assembly or dies, and print stats to screen """
    try:
        data = ip.load_assembly(
                    name=os.path.join(parsedict['1'], 
                                      parsedict['14']),
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

    except AssertionError as inst:
        sys.exit("Error: No Assembly file found at {}. ".format(inst)\
        +"\nCheck parameter settings for [working_dir]/[prefix_name]\n")



def getassembly(args, parsedict):
    """ loads assembly or creates a new one and set its params from dict"""
    ## if '1' then do not load existing Assembly
    if '1' in args.steps:
        if args.force:
            data = ip.Assembly(parsedict['14'])
        else:
            try:
                data = ip.load_assembly(os.path.join(parsedict['1'], 
                                                     parsedict['14']))
            except AssertionError:
                data = ip.Assembly(parsedict['14'])


    ## otherwise look for existing
    else:
        try:
            data = ip.load_assembly(os.path.join(parsedict['1'], 
                                                 parsedict['14']))
        except AssertionError:
            data = ip.Assembly(parsedict['14'])

    ## for entering some params...
    for param in parsedict:
        if parsedict[param]:
            data.set_params(param, parsedict[param])

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


    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## used to restrict some CLI workflow vs. interactive features
    ## does it need to be global?
    ip.__interactive__ = 1
    LOGGER.info("Running *CLI*")

    ## create new paramsfile if -n
    if args.new:
        write_params(ip.__version__)
        print("New file `params.txt` created in {}".\
               format(os.path.realpath(os.path.curdir)))
        sys.exit(2)

    ## if showing results, do not do step 1
    if args.results:
        args.steps = ""
        print("")
    else:
        header = \
    "\n --------------------------------------------------"+\
    "\n  ipyrad [v.{}]".format(ip.__version__)+\
    "\n  Interactive assembly and analysis of RADseq data"+\
    "\n --------------------------------------------------\n"
        print(header)


    ## create new Assembly or load existing Assembly, quit if args.results
    if args.params:
        parsedict = parse_params(args)

        if args.results:
            showstats(parsedict)

        else:
            data = getassembly(args, parsedict)
            ## For now print the params. 
            for key, item in data.paramsdict.items():
                LOGGER.debug("%s, %s", key, item)

            ## run Assembly steps
            steps = list(args.steps)
            data.run(steps=steps, force=args.force)

    print("")


if __name__ == "__main__": 
    main()


