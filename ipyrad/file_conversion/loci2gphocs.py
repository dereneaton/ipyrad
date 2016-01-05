#!/usr/bin/env python2

# pyrad .loci to gphocs format conversion
#
# This is a very simple conversion because the formats are very similar.
#
# Isaac Overcast
# March 21, 2015

## import libraries
import sys, os

def make(WORK, outname):
    
    #read in loci file
    infile = open(WORK+"outfiles/"+outname+".loci")
    outfile = open(WORK+"outfiles/"+outname+".gphocs",'w')

    ## parse the loci
    ## Each set of reads at a locus is appended with a line
    ## beginning with // and ending with |x, where x in the locus id. 
    ## so after this call 'loci' will contain an array
    ## of sets of each read per locus.
    loci = infile.read().strip().split("//")[:-1]
    
    # Print the header, the number of loci in this file
    print >>outfile, len(loci)#-1
    
    # iterate through each locus, print out the header for each locus:
    # <locus_name> <n_samples> <locus_length>
    # Then print the data for each sample in this format:
    # <individual_name> <sequence>
    for i, loc in enumerate(loci):
        # Separate out each sequence within the loc block. 'sequences'
        # will now be a list strings containing name/sequence pairs.
        # We select each line in the locus string that starts with ">"
        names = [line.split()[0] for line in loc.strip().split("\n") if ">" in line]
        sequences = [line.split()[-1] for line in loc.strip().split("\n") if ">" in line] 

        # Strips off 'nnnn' separator for paired data
        # replaces '-' with 'N'
        editsequences = [seq.replace("n","").replace('-','N') for seq in sequences]
        sequence_length = len(editsequences[0])

        # get length of longest name and add 4 spaces
        longname = max(map(len,names))+4
        
        # Print out the header for this locus
        print >>outfile, 'locus'+ str(i), len(sequences), str( sequence_length )
    
        # Iterate through each sequence read at this locus and write it to the file.
        for name,sequence in zip(names,sequences):
            # Clean up the sequence data to make gphocs happy. Only accepts UPPER
            # case chars for bases, and only accepts 'N' for missing data. Also,
            # the .loci format prepends a '>' on to the individual names, so we have
            # to clean this up which is what the [1:] is doing.
            print >>outfile, name[1:]+" "*(longname-len(name))+sequence

if __name__ == "__main__":
    make(WORK, outfile)
