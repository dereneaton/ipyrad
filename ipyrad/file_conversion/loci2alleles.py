#!/usr/bin/env python2

""" pyrad .loci to .alleles format conversion. This is a very simple 
conversion because the formats are very similar. More or less directly ripped 
from pyrad.alignable.makehaplos. """

## import libraries
import sys, os
from ipyrad.assemble.util import *


def make(data, samples):
    """ reads in .loci and builds alleles from case characters """
    
    #read in loci file
    outfile = open(os.path.join(data.dirs.outfiles, data.name+".alleles"), 'w')
    lines = open(os.path.join(data.dirs.outfiles, data.name+".loci"), 'r')

    ## Get the longest sample name for pretty printing
    longname = max(len(x) for x in data.samples.keys())

    ## Padding between name and sequence in output file. This should be the 
    ## same as write_outfiles.write_tmp_loci.name_padding
    name_padding = 5
    writing = []
    loc = 0
    for line in lines:
        if ">" in line:
            a, b = line.split(" ")[0],line.split(" ")[-1]
            a1, a2 = breakalleles(b.strip())

            ## Format the output string. the "-2" below accounts for the additional
            ## 2 characters added to the sample name that don't get added to the
            ## snpsites line, so you gotta bump this line back 2 to make it
            ## line up right.
            writing.append(a+"_0"+" "*(longname-len(a)-2+name_padding)+a1)
            writing.append(a+"_1"+" "*(longname-len(a)-2+name_padding)+a2)
        else:
            writing.append(line.strip())
        loc += 1

        ## print every 10K loci "
        if not loc % 10000:
            outfile.write("\n".join(writing)+"\n")
            writing = []

    outfile.write("\n".join(writing))
    outfile.close()


if __name__ == "__main__":
    pass
