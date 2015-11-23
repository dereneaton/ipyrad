#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

from ipyrad.core.newparamsfile import write_params
import pkg_resources
import ipyrad as ip
import argparse
import sys
import os



def parse_params(params):
    """ Parse the params file args, create and return Assembly object."""
    ## check that params.txt file is correctly formatted.
    with open(params) as paramsin:
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

    ## create a default Assembly object
    print('parsedict:\n', parsedict)
    data = ip.Assembly(parsedict['14'])

    ## set_params for all keys in parsedict 
    for param in parsedict:
        data.set_params(param, parsedict[param])

    return data


def parse_params():
    print("not yet")

def parse_command_line():
    """ Parse CLI. Only three options now. """
    ## -p  choose params file
    ## -s  subselect steps
    ## -n  create new params file
    print("not yet")

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    ipyrad -n                      ## create new params.txt file.
    ipyrad -p params.txt           ## run ipyrad with settings in params.txt.
    ipyrad -p params.txt -s 123    ## run only steps 1, 2 and 3 of assembly.
    ipyrad -p params.txt -s 4567   ## run steps 4, 5, 6 and 7 of assembly.

  * See documentation for writing advanced ipyrad scripts.
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('ipyrad')))

    parser.add_argument('-n', "--new", action='store_true',
        help="create new default params.txt file in current directory.")

    #parser.add_argument('-q', "--quiet", action='store_true',
    #    help="do not print to stderror and stdout.")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str, default="params.txt",
        help="path to params.txt file.")

    parser.add_argument('-s', metavar="steps", dest="steps",
        type=str, default="1234567",
        help="subset of assembly steps to perform. Default=1234567")

    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()
    return args



def main():
    """ main function """
    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## create new paramsfile if -n
    if args.new:
        write_params(ip.__version__)
        print("New file `params.txt` created in {}".\
               format(os.path.realpath(os.path.curdir)))
        sys.exit(2)

    ## parse params file
    data = parse_params(args.params)

    ## For now print the params. 
    print("")
    for key, item in data.paramsdict.items():
        print("{:<30} {:<20}".format(key, item))

    ## run assembly steps
    #data.run(args.steps)


if __name__ == "__main__": 
    main()
