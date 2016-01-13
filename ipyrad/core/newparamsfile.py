#!/usr/bin/env python2.7

""" creates a new params.txt file with default entries """

from __future__ import print_function


def write_params(version):
    """ write default parameter settings to file params.txt. """
    paramstr = \
"""------ ipyrad params file (v.{}) ---------------------------------------
test                       ## [{}] Assembly name (unique prefix for file names)
./                         ## [{}] Working dir (made in curdir if not present)
./null_*.fastq.gz          ## [{}] Location of raw non-demultiplexed fastq files
./null_barcodes.txt        ## [{}] Location of barcodes file
                           ## [{}] Location of demultiplexed/sorted fastq files
                           ## [{}] Location of reference genome file
denovo                     ## [{}] Assembly method (denovo, ref., hyb., subref.)
rad                        ## [{}] Datatype (see docs): rad, gbs, ddrad, etc.
(TGCAG,)                   ## [{}] Restriction overhang (cut1,) or (cut1, cut2)
5                          ## [{}] Max low quality base calls (Q<20) in a read
6                          ## [{}] Mindepth for statistical base calls
6                          ## [{}] Mindepth for majority rule base calls
0.85                       ## [{}] clustering threshold for de novo assembly
4                          ## [{}] min number of Samples with data at a locus
.25                        ## [{}] max number/prop Samples heterozyg. at a site
33                         ## [{}] phred Q score offset (only alternative=64)
1                          ## [{}] max number mismatches in barcodes
0                          ## [{}] filter for adapters/primers (1 or 2=stricter)
32                         ## [{}] minimum length of reads after adapter trim
2                          ## [{}] max_alleles_consens: 1=haploid, 2=diploid, >2=(see docs)
1000                       ## [{}] max cluster depth within samples
5                          ## [{}] maxNs (uncalled bases) in consensus (R1, R2)
5                          ## [{}] maxHs (heterozygotes) in consensus (R1, R2)
(100,100)                  ## [{}] maxSNPs in a locus (R1, R2)
(5,100)                    ## [{}] maxIndels in a locus (R1, R2)
(1,2,2,1)                  ## [{}] trim overhang (see docs) (R1>, <R1, R2>, <R2)
(0,0)                      ## [{}] edit cut-sites (R1, R2) (see docs)
lpn                        ## [{}] output formats (see docs)
------ optional: list group/clade assignments below this line (see docs) ---
""".format(version, *range(28))

    with open("params.txt", 'w') as paramsfile:
        print(paramstr, file=paramsfile)

if __name__ == "__main__":
    write_params("test")

