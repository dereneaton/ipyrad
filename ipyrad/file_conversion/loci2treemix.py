#!/usr/bin/env python2

import numpy as np
import sys
import gzip
from collections import OrderedDict, Counter
from ipyrad.assemble.util import *

def make(WORK, outname, taxadict, minhits):

    ## output files
    outfile = gzip.open(WORK+"/outfiles/"+outname+".treemix.gz",'w')

    ## cleanup taxadict to just sample names
    taxa = OrderedDict()
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            a = samp.split("/")[-1].replace(".consens.gz","")
            taxa[group].append(a)

    print "\t    data set reduced for group coverage minimums"        
    for i,j in zip(taxa,minhits):
        print "\t   ",i, taxa[i], 'minimum=',j
    
    ## read in data from unlinked_snps to sample names
    infile = open(WORK.rstrip("/")+"/outfiles/"+outname+".unlinked_snps",'r')
    dat = infile.readlines()
    nsamp,nsnps = dat[0].strip().split(" ")
    nsamp = int(nsamp)
    nsnps = int(nsnps)
    NDATA = np.empty([int(nsamp),int(nsnps)],dtype='object')
    excludes = 0

    ## read SNP matrix into a numpy.array
    for line in range(len(dat[1:])):
        a,b = dat[1:][line].split()
        NDATA[line] = list(b)
    sites = np.transpose(NDATA)

    ## unpack ambiguity bases and find two most common alleles
    ## at every SNP site, save to a list
    alleles = []
    for site in sites:
        ds = []
        for s in site:
            if s in list("RKSYWM"):
                ds.append(unstruct(s)[0])
                ds.append(unstruct(s)[1])
            else:
                ds.append(s)
                ds.append(s)
        snp = [s for s in ds if s not in ["N",'-']]
        a = Counter(snp).most_common(3)
        alleles.append([a[0][0],a[1][0]])

    ## create a dictionary mapping sample names to SNPs    
    SNPS = OrderedDict()
    for line in dat[1:]:
        a,b = line.split()
        SNPS[a] = b

    ## reduce Taxa dict to only samples that are in the unlinkedsnps alignment
    for key in taxa:
        replacement = []
        for val in taxa[key]:
            if val in SNPS.keys():
                replacement.append(val)
        taxa[key] = replacement

    ## create a dictionary with empty lists for each taxon 
    FREQ = OrderedDict()
    for tax in taxa:
        FREQ[tax] = []

    ## fill the FREQ dictionary with SNPs for all 
    ## samples in that taxon
    keeps = []
    for snp in range(int(nsnps)):
        GG = []
        ## if snp meets minhits requirement
        for tax,mins in zip(taxa,minhits):
            GG.append( sum([SNPS[i][snp] not in ["N","-"] for i in taxa[tax]]) >= int(mins))
        if all(GG):
            keeps.append(snp)


    for keep in keeps:
        for tax in FREQ:
            bunch = []
            for i in taxa[tax]:
                bunch.append(unstruct(SNPS[i][keep])[0])
                bunch.append(unstruct(SNPS[i][keep])[1])
                #print tax, i, SNPS[i][keep], bunch
            FREQ[tax].append("".join(bunch))

    ## header
    print >>outfile, " ".join(FREQ.keys())

    ## data to file
    for i,j in enumerate(keeps):
        a1 = alleles[j][0]
        a2 = alleles[j][1]
        H = [str(FREQ[tax][i].count(a1))+","+str(FREQ[tax][i].count(a2)) for tax in FREQ]
        HH = " ".join(H)

        ## exclude non-biallelic SNPs
        if " 0,0 " not in HH:
            ## exclude invariable sites given this sampling
            if not all([zz.split(",")[1] in '0' for zz in H]):
                print >>outfile, " ".join(H)
        else:
            excludes += 1

    outfile.close()
    

if __name__ == "__main__":
    make(WORK, outname, taxadict)
