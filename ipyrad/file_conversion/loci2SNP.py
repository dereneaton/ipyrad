#!/usr/bin/env python2

import numpy as np
import sys
import gzip
from collections import Counter
from itertools import chain
from ipyrad.assemble.util import *

def make( data, samples ):

    ploidy = data.paramsdict["ploidy"]
    names = map( lambda x: x.name, samples )
    longname = max(map(len,names))
    inloci   =  open(os.path.join( data.dirs.outfiles, data.name+".loci" ), 'r').read()

    ## Potential outfiles
    snpsout = os.path.join( data.dirs.outfiles, data.name+".snps" )
    usnpsout = os.path.join( data.dirs.outfiles, data.name+".usnps" )
    structout = os.path.join( data.dirs.outfiles, data.name+".str" )
    genoout = os.path.join( data.dirs.outfiles, data.name+".snps.geno" )
    ugenoout = os.path.join( data.dirs.outfiles, data.name+".usnps.geno" )

    ## Output file for writing some stats
    statsfile= os.path.join( data.dirs.outfiles, data.name+".snps.stats" )

    ## The output formats to write
    formats = data.paramsdict["output_formats"]

    seed = data._hackersonly["random_seed"]
    np.random.seed(int(seed))

    " output .snps and .unlinked_snps"
    S = {}      ## snp dict
    Si = {}     ## unlinked snp dict
    for name in list(names):
        S[name] = []
        Si[name] = []

    " record bi-allelic snps"
    bis = 0

    " for each locus select out the SNPs"
    for loc in inloci.strip().split("|")[:-1]:
        pis = ""
        ns = []
        ss = []
        cov = {}  ## record coverage for each SNP
        for line in loc.split("\n"):
            if ">" in line:
                ns.append(line.split()[0].replace(">",""))
                ss.append(line.split()[-1])
            else:
                pis = [i[0] for i in enumerate(line) if i[1] in list('*-')]
                LOGGER.debug("pis - {}".format(pis))                

        " assign snps to S, and record coverage for usnps"
        for tax in S:
            if tax in ns:
                if pis:
                    for snpsite in pis:
                        snpsite -= (longname+5)
                        S[tax].append(ss[ns.index(tax)][snpsite])
                        LOGGER.debug("pis2 - snpsite {} tax {} ss {}".format(snpsite, tax, ss[ns.index(tax)]))
                        if snpsite not in cov:
                            cov[snpsite] = 1
                        else:
                            cov[snpsite] += 1
                        "downweight selection of gap sites "
                        if ss[ns.index(tax)][snpsite] != '-':
                           cov[snpsite] += 1
            else:
                if pis:
                    for snpsite in pis:
                        S[tax].append("N")
                    Si[tax].append("N")

        " randomly select among snps w/ greatest coverage for unlinked snp "
        maxlist = []
        for j,k in cov.items():
            if k == max(cov.values()):
                maxlist.append(j)

        " Is bi-allelic ? "
        bisnps = []
        for i in maxlist:
            if len(set([ss[ns.index(tax)][i] for tax in S if tax in ns])) < 3:
                bisnps.append(i)

        #rando = pis[np.random.randint(len(pis))]
        #rando -= (longname+5)
        if bisnps:
            rando = bisnps[np.random.randint(len(bisnps))]
        elif maxlist:
            rando = maxlist[np.random.randint(len(maxlist))]
        tbi = 0
        for tax in S:
            if tax in ns:
                if pis:
                    " if none are bi-allelic "
                    if not bisnps:
                        tbi += 1
                    Si[tax].append(ss[ns.index(tax)][rando])
            if pis:
                " add spacer between loci "                
                S[tax].append(" ")
            else:
                " invariable locus "
                S[tax].append("_ ")
        bis += tbi
    " names "
    SF = list(S.keys())
    SF.sort()

    ## Write linked snps format
    if "snps" in formats:
        with open(snpsout, 'w') as outfile:
            print >>outfile, "## %s taxa, %s loci, %s snps" % (len(S),
                                                               len("".join(S.values()[0]).split(" "))-1,
                                                               len("".join(S[SF[0]]).replace(" ","")))
            for i in SF:
                print >>outfile, i+(" "*(longname-len(i)+3))+"".join(S[i])

    ## Write unlinked snps format
    if "usnps" in formats:
        with open(usnpsout, 'w') as outfile:
            print >>outfile, len(Si),len("".join(Si.values()[0]))
            for i in SF:
                print >>outfile, i+(" "*(longname-len(i)+3))+"".join(Si[i])

        with open(statsfile, 'a') as statsout:
            print >>statsout, "sampled unlinked SNPs=",len(Si.values()[0])
            print >>statsout, "sampled unlinked bi-allelic SNPs=",len(Si.values()[0])-bis

    ## Write STRUCTURE format
    if "str" in formats:
        with open(structout, 'w') as outfile:
        
            B = {'A': '0',
                 'T': '1',
                 'G': '2',
                 'C': '3',
                 'N': '-9',
                 '-': '-9'}
            if ploidy > 1:
                for line in SF:
                    print >>outfile, line+(" "*(longname-len(line)+3))+\
                          "\t"*6+"\t".join([B[unstruct(j)[0]] for j in Si[line]])
                    print >>outfile, line+(" "*(longname-len(line)+3))+\
                          "\t"*6+"\t".join([B[unstruct(j)[1]] for j in Si[line]])
            else:
                for line in SF:
                    print >>outfile, line+(" "*(longname-len(line)+3))+\
                          "\t"*6+"\t".join([B[unstruct(j)[1]] for j in Si[line]])

    ## Do linked and unlinked snps in .geno format
    if "geno" in formats:
        with open(ugenoout, 'w') as outfile:

            for i in range(len(Si.values()[0])):
                getref = 0
                ref = "N"
                while ref == "N":
                    ref = unstruct(Si[SF[getref]][i])[0]
                    getref += 1
                SNProw = "".join(map(str,[unstruct(Si[j][i]).count(ref) if Si[j][i] != "N" \
                                          else "9" for j in SF]))
                ## print ref,SNProw
                if len(set(SNProw)) > 1:
                    print >>outfile, SNProw 

        with open(genoout, 'w') as outfile:
            for i in range(len(S.values()[0])):
                if S[SF[0]][i].strip("_").strip():
                    getref = 0
                    ref = "N"
                    while ref == "N":
                        #print i, S[SF[0]][i]
                        ref = unstruct(S[SF[getref]][i])[0]
                        getref += 1
                        SNProw = "".join(map(str,[unstruct(S[j][i]).count(ref) if \
                                                  S[j][i] != "N" else "9" for j in SF]))
                    ## print ref,SNProw
                    if len(set(SNProw)) > 1:
                        print >>outfile, SNProw 


if __name__ == "__main__":
    make( data, samples )
