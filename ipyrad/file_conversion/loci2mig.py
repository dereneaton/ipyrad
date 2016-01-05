#!/usr/bin/env python2

import numpy as np
import sys
import gzip
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
try:
    from collections import Counter
except ImportError:
    from counter import Counter
#import alignable

def make(WORK, outname, taxadict, minhits, seed):

    ## outfile
    outfile = open(WORK+"/outfiles/"+outname+".migrate", 'w')

    ## cleanup taxadict
    taxa = OrderedDict()
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            a = samp.split("/")[-1].replace(".consens.gz","")
            taxa[group].append(a)

    print "\t    data set reduced for group coverage minimums"
    for i,j in zip(taxa,minhits):
        print "\t   ",i, taxa[i], "minimum=",j

    #print taxadict.keys()

    ## filter data to only the loci that have data
    ## for at least N individuals in each pop
    keep = []
    MINS = zip(taxa.keys(), minhits)

    ## read in data to sample names
    loci  = open(WORK+"/outfiles/"+outname+".loci",'r').read().strip().split("|")[:-1]
    for loc in loci:
        samps = [i.split()[0].replace(">","") for i in loc.split("\n") if ">" in i]
        ## filter for coverage
        GG = []
        for group,mins in MINS:
            GG.append( sum([i in samps for i in taxa[group]]) >= int(mins) )
        if all(GG):
            keep.append(loc)

    ## print data to file
    print >>outfile, len(taxa), len(keep), "( npops nloci for data set", outname+".loci",")"
    
    ## print all data for each population at a time
    done = 0
    for group in taxadict:
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


# WORK = "/home/deren/Dropbox/Public/PyRAD_TUTORIALS/tutorial_RAD"
# outname = "c85m4p3"

# pops = ['pop1','pop2','pop3']
# samps = [ ["1A0","1B0","1C0","1D0"],
#           ["2E0","2F0","2G0","2H0"],
#           ["3I0","3J0","3K0","3L0"] ]

# taxadict = OrderedDict(zip(pops,samps))
# minhits = [4,4,4]
# seed = 112233

if __name__ == "__main__":
    make(WORK, outname, taxadict, minhits, seed)
