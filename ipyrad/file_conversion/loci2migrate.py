#!/usr/bin/env python2

import numpy as np
import os
import sys
import gzip
from collections import OrderedDict, Counter

import logging
LOGGER = logging.getLogger(__name__)

def make( data, samples ):

    ## TODO: Fix migrate format
    print("Migrate format is still in development")
    return

    outfile  =  open(os.path.join(data.dirs.outfiles, data.name+".migrate"), 'w')
    infile =  open(os.path.join( data.dirs.outfiles, data.name+".loci" ), 'r' )

    ## TODO: Allow for subsampling the output by honoring the "samples" list passed in.

    ## Make sure we have population assignments for this format
    try:
        taxa = data.populations
    except AttributeError:
        LOGGER.error( "Migrate file output requires population assignments \
                        and this data.populations is empty. Make sure you have \
                        set the 'pop_assign_file' parameter, and make sure the \
                        path is correct and the format is right." )
        return

    ## TODO: Hax. pyRAD v3 used to allow specifying minimum coverage per group
    ## This is a hackish version where we just use mindepth_statistical. Maybe
    ## it's best not to filter
    ## minhits = [ data.paramsdict["mindepth_statistical"] ] * len(taxa)

    ## Hard coding minhits to 2 individuals per group.
    ## TODO: Fixme.
    minhits = [ 2 ] * len(taxa)

    print "\t    data set reduced for group coverage minimums"
    for i,j in zip(taxa,minhits):
        print "\t   ",i, taxa[i], "minimum=",j


    ## filter data to only the loci that have data
    ## for at least N individuals in each pop
    keep = []
    MINS = zip(taxa.keys(), minhits)

    ## read in data to sample names
    loci  = infile.read().strip().split("|")[:-1]
    for loc in loci:
        samps = [i.split()[0].replace(">","") for i in loc.split("\n") if ">" in i]
        ## filter for coverage
        GG = []
        for group,mins in MINS:
            GG.append( sum([i in samps for i in taxa[group]]) >= int(mins) )
        if all(GG):
            keep.append(loc)

    ## print data to file
    print >>outfile, len(taxa), len(keep), "( npops nloci for data set", data.name+".loci",")"
    
    ## print all data for each population at a time
    done = 0
    for group in taxa:
        ## print a list of lengths of each locus
        if not done:
            loclens = [len(loc.split("\n")[1].split()[-1].replace("x","n").replace("n","")) for loc in keep]
            print >>outfile, " ".join(map(str,loclens))
            done += 1

        ## print a list of number of individuals in each locus
        indslist = []
        for loc in keep:
            samps = [i.split()[0].replace(">","") for i in loc.split("\n") if ">" in i]
            inds = sum([i in samps for i in taxa[group]])
            indslist.append(inds)
        print >>outfile, " ".join(map(str,indslist)), group

        ## print sample id, spaces, and sequence data
        #for loc in range(len(keep)):
        for loc in range(len(keep)):
            seqs = [i.split()[-1] for i in keep[loc].split("\n") if \
                    i.split()[0].replace(">","") in taxa[group]]
            for i in range(len(seqs)):
                print >>outfile, group[0:8]+"_"+str(i)+\
                      (" "*(10-len(group[0:8]+"_"+str(i))))+seqs[i].replace("x","n").replace("n","")
            
    outfile.close()


if __name__ == "__main__":
    make( data, samples )
