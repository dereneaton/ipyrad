#!/usr/bin/env python2
""" a function to produce migrate-n input files from an ipyrad .loci file"""

import os
import numpy as np

## shorthand
OPJ = os.path.join



def loci2migrate(name, locifile, popdict, mindict=1):
    """  
    A function to build an input file for the program migrate from an ipyrad 
    .loci file, and a dictionary grouping Samples into populations. 

    Parameters:
    -----------
    name: (str)
       The name prefix for the migrate formatted output file.
    locifile: (str)
       The path to the .loci file produced by ipyrad. 
    popdict: (dict)
       A Python dictionary grouping Samples into Populations. 
    
    Examples:
    ---------
    You can create the population dictionary by hand, and pass in the path 
    to your .loci file as a string. 
       >> popdict = {'A': ['a', 'b', 'c'], 'B': ['d', 'e', 'f']}
       >> loci2migrate("outfile.migrate", "./mydata.loci", popdict)

    Or, if you load your ipyrad.Assembly object from it's JSON file, you can
    access the loci file path and population information from there directly. 
       >> data = ip.load_json("mydata.json")
       >> loci2migrate("outfile.migrate", data.outfiles.loci, data.populations)

    """

    ## I/O
    outfile = open(name+".migrate", 'w')
    infile = open(locifile, 'r')

    ## minhits dictionary can be an int (all same) or a dictionary (set each)
    if isinstance(mindict, int):
        mindict = {pop: mindict for pop in popdict}
    else:
        mindict = mindict

    ## filter data to only the loci that have data for mindict setting
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

    ## example run
    LOCIFILE = "./test.loci"
    POPDICT = {"A": ['a'], "B": ['b']}
    MINDICT = 1
    loci2migrate("test", LOCIFILE, POPDICT, MINDICT)



