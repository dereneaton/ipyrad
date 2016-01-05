#!/usr/bin/env python2

import time
import numpy as np
#import alignable


def make(WORK, version, outname, mindepth, names):
    outfile  =  open(WORK+"/outfiles/"+outname+".vcf", 'w')
    inloci   =  WORK+"/outfiles/"+outname+".loci"
    names = list(names)
    names.sort()

    print >>outfile, "##fileformat=VCFv4.1"
    print >>outfile, "##fileDate="+time.strftime("%Y%m%d")
    print >>outfile, "##source=pyRAD.v."+str(version)
    print >>outfile, "##reference=common_allele_at_each_locus"
    print >>outfile, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    print >>outfile, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    print >>outfile, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    print >>outfile, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"
    print >>outfile, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print >>outfile, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    print >>outfile, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    print >>outfile, "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO    ","FORMAT"]+list(names))

    loci = open(inloci).read().split("|")[:-1]
    snps = 0
    vcflist = []
    for locusnumber in range(len(loci)):
        samps = [i.split()[0][1:] for i in loci[locusnumber].strip().split("\n") if ">" in i]
        loc = np.array([tuple(i.split()[-1]) for i in loci[locusnumber].strip().split("\n") if ">" in i])
        NS = str(len(loc))
        DP = str(mindepth)
        for base in range(len(loc.T)):
            col = []
            site = list(loc.T[base])
            site = list("".join(site).replace("-","").replace("N",""))
            if site:
                for bb in site:
                    if bb in list("RKYSWM"):
                        col += alignable.unstruct(bb)[0]
                        col += alignable.unstruct(bb)[1]
                    else:
                        col += bb
                REF = alignable.most_common([i for i in col if i not in list("-RKYSWMN")])
                ALT = set([i for i in col if (i in list("ATGC-N")) and (i!=REF)])
                if ALT:
                    snps += 1
                    GENO = [REF]+list(ALT)
                    GENOS = []
                    for samp in names:
                        if samp in samps:
                            idx = samps.index(samp)
                            f = alignable.unstruct(loc.T[base][idx])
                            if ('-' in f) or ('N' in f):
                                GENOS.append("./.")
                            else:
                                GENOS.append(str(GENO.index(f[0]))+"|"+str(GENO.index(f[1])))
                        else:
                            GENOS.append("./.")
                    vcflist.append("\t".join([`locusnumber+1`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
                                              ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS))
        if not locusnumber % 1000:
            outfile.write( "\n".join(vcflist)+"\n" )
            vcflist = []
                                              
                    #print >>outfile, "\t".join([`locusnumber+1`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
                    #                            ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS)
    

    outfile.write( "\n".join(vcflist) )
    outfile.close()

if __name__ == "__main__":
    make(WORK, version, outname, mindepth, names)
