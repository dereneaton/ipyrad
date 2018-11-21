#!/usr/bin/env python2

import os
import sys
import glob

def make(nexdir, outfile):

    nexs = glob.glob(nexdir + "/*.nexus")
    with open(outfile, 'w') as out:
        for i, nex in enumerate(nexs):
            with open(nex, 'r') as infile:
                lines = infile.readlines()[5:-2]
                lines = [line.replace("?", "N") for line in lines]
                out.write("".join(lines))
                out.write("//\t|{}|\n".format(i))
        

if __name__ == "__main__":
    print("  nexus2loci is experimental and very beta, so don't be surprised if it doesn't work.")
    if not len(sys.argv) == 3:
        msg = """
  nexus2loci requires 2 arguments: The directory of nexus files and the name of the output .loci
    usage: python nexus2loci.py <nexus_dir> <output.loci>
"""
        sys.exit(msg)
    nexdir = sys.argv[1]
    outfile = sys.argv[2]
    make(nexdir, outfile)
