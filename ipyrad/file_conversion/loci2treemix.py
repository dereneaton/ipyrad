#!/usr/bin/env python2

import numpy as np
import sys
import gzip
from collections import OrderedDict, Counter
from ipyrad.assemble.util import *
from ipyrad.file_conversion import loci2SNP

def make( data, samples ):

    ## output files
    outfile  =  gzip.open(os.path.join(data.dirs.outfiles, data.name+".treemix.gz"), 'w')

    try:
        infile =  open(os.path.join( data.dirs.outfiles, data.name+".usnps" ), 'r' )
    except FileNotFoundError:
        LOGGER.info( "unlinked_snps file doesn't exist. Try creating it." )
        try:
            ## Try generating the .unlinked_snps file
            data.paramsdict["output_formats"] = "usnps,"+data.paramsdict["output_formats"]
            loci2SNP( data, samples )
        except Exception:
            LOGGER.error("Treemix file conversion requires .unlinked_snps file, which does not exist. \
                        Make sure the param `output_formats` includes at least `usnps,treemix` \
                        and rerun step7()")
            return

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

    ## Hard coding minhits to 2 individuals per group. Fixme.
    minhits = [ 2 ] * len(taxa)

    print "\t    data set reduced for group coverage minimums"        
    for i,j in zip(taxa,minhits):
        print "\t   ",i, taxa[i], 'minimum=',j
    
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
