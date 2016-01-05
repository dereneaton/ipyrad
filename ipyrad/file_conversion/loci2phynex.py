#!/usr/bin/env python2

import numpy as np
import sys
import os
import glob


def update(idict, count, WORK, outname):
    """ updates dictionary with the next .5M reads 
    from the super long string phylip file. Makes
    for faster reading. """

    data = iter(open(WORK+"outfiles/"+outname+".phy"))
    ntax, nchar = data.next().strip().split()

    ## read in max N bp at a time                                                                            
    for line in data:
        tax, seq = line.strip().split()
        idict[tax] = idict[tax][100000:]
        idict[tax] += seq[count:count+100000]
    del line

    return idict
    


def makephy(WORK, outname, names, longname):
    """ builds phy output. If large files writes 50000 loci 
        at a time to tmp files and rebuilds at the end"""

    " order names "
    names = list(names)
    names.sort()
    
    " read in loci file "
    locus = iter(open(WORK+"outfiles/"+outname+".loci", 'rb'))

    " dict for saving the full matrix "
    fdict = {name:[] for name in names}

    " list for saving locus number and locus range for partitions "
    partitions = []
    loc_number = 1
    initial_pos = 1

    " remove empty column sites and append edited seqs to dict F "
    done = 0
    nloci = 0
    nbases = 0
    while nloci < 50000: #not done:
        seqs = []
        #arrayed = np.array([])
        anames = []
        while 1:
            ## get next locus
            try:
                samp = locus.next()
            except StopIteration:
                done = 1
                break
            if "//" in samp:
                nloci += 1
                break
            else:
                try:
                    name, seq = samp.split()
                except ValueError:
                    print samp
                anames.append(name[1:])
                seqs.append(seq.strip())
        ## reset
        arrayed = np.array([list(i) for i in seqs])
        if done:
            break
        ## create mask for columns that are empty or 
        ## that are paired-end separators (compatible w/ pyrad v2 and v3)
        #mask = [i for i in range(len(arrayed.T)) if np.any([
        ## still surely a better way to vectorize this...
        mask = [i for i in arrayed.T if any([j not in list("-Nn") for j in i])]
        masked = np.dstack(mask)[0]

        ## partition information
        loc_name = "p"+str(nloci)
        loc_range = str(initial_pos) + "-" +\
                    str(len(masked[0]) + initial_pos -1)
        initial_pos += len(masked[0])
        partitions.append(loc_name+"="+loc_range)

        ## uncomment to print block info (used to partition by locus)
        #blockend += minray
        #print blockend,
        #print loc
        #print arrayed

        ## append data to dict
        for name in names:
            if name in anames:
                #fdict[name].append(arrayed[anames.index(name), mask].tostring())
                fdict[name].append(masked[anames.index(name),:].tostring())
            else:
                fdict[name].append("N"*masked.shape[1])
                #fdict[name].append("N"*len(arrayed[0, mask]))
        ## add len to total length
        nbases += len(fdict[name][-1])

        ## after x iterations tmp pickle fdict?
        if not nloci % 1e4:
            ## concat strings
            for name in fdict:
                with open(os.path.join(WORK, "tmp", 
                    "{}_{}.tmp".format(name, nloci)), 'wb') as wout:
                    wout.write("".join(fdict[name]))
            del fdict
            fdict = {name:[] for name in names}

    ## print out .PHY file, if really big, pull form multiple tmp pickle
    superout = open(WORK+"outfiles/"+outname+".phy", 'wb')
    print >>superout, len(names), nbases
    if nloci < 1e4:
        for name in names:
            print >>superout, name+(" "*((longname+3)-\
                              len(name)))+"".join(fdict[name])
    else:
        for name in names:
            superout.write("{}{}{}".format(
                            name,
                            " "*((longname+3)-len(name)),
                            "".join(fdict[name])))
            tmpfiles = glob.glob(os.path.join(WORK, "tmp", name+"*.tmp"))
            tmpfiles.sort()
            for tmpf in tmpfiles:
                with open(tmpf, 'rb') as tmpin:
                    superout.write(tmpin.read())
            superout.write("\n")
    superout.close()
    raxml_part_out = open(WORK+"outfiles/"+outname+".phy.partitions", 'w')
    for partition in partitions:
        print >>raxml_part_out, "DNA, %s" % (partition)
    raxml_part_out.close()

    return partitions


def makenex(WORK, outname, names, longname, partitions):
    """ PRINT NEXUS """

    " make nexus output "
    data   = iter(open(WORK+"outfiles/"+outname+".phy"))
    nexout = open(WORK+"outfiles/"+outname+".nex", 'wb')

    ntax, nchar = data.next().strip().split(" ")

    print >>nexout, "#NEXUS"
    print >>nexout, "BEGIN DATA;"
    print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax,nchar)
    print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
    print >>nexout, "  MATRIX"

    idict = {}

    ## read in max 1M bp at a time
    for line in data:
        tax, seq = line.strip().split()
        idict[tax] = seq[0:100000]
    del line

    nameorder = idict.keys()
    nameorder.sort()

    n=0
    tempn=0
    sz = 100
    while n < len(seq):
        for tax in nameorder:
            print >>nexout, "  "+tax+" "*\
                             ((longname-len(tax))+3)+\
                             idict[tax][tempn:tempn+sz]
        n += sz
        tempn += sz
        print >>nexout, ""

        if not n % 100000:
            #print idict[tax][tempn:tempn+sz]
            idict = update(idict, n, WORK, outname)
            tempn -= 100000
            
    print >>nexout, ';'
    print >>nexout, 'END;'
    
    ### partitions info
    print >>nexout, "BEGIN SETS;"
    for partition in partitions:
        print >>nexout, "  CHARSET %s;" % (partition)
    print >>nexout, "END;"

    nexout.close()
        

def make(WORK, outfile, names, longname, formats):
    partitions = makephy(WORK, outfile, names, longname)
    makenex(WORK, outfile, names, longname, partitions)
    

if __name__ == "__main__":
    make()
