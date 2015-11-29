#!/usr/bin/env python2.7

""" creates a new params.txt file with default entries """

from __future__ import print_function


def write_params(version):
    """ write default parameter settings to file params.txt. """
    paramstr = \
"""------ ipyrad params file (v.%s) ---------------------------------------
./                         ## [1] Working directory
./*.fastq.gz               ## [2] Location of raw non-demultiplexed fastq files
./*_barcodes.txt           ## [3] Location of barcodes file.
./fastqs/*.fastq.gz        ## [4] Location of demultiplexed/sorted fastq files
TGCAG                      ## [5] Restriction overhang (cut1) or (cut1, cut2)
5                          ## [6] max low quality base calls (Q<20) in a read
4                          ## [7] N engines (threads) per job
6                          ## [8] Mindepth for statistical base calls
6                          ## [9] Mindepth for majority rule base calls
rad                        ## [10] datatype (see docs): rad, gbs, ddrad, etc.
0.85                       ## [11] clustering threshold for de novo assembly
4                          ## [12] min number of Samples with data at a locus
.25                        ## [13] max number/prop Samples heterozyg. at a site
ipyrad_test                ## [14] prefix name for saved output files
33                         ## [15] phred Q score offset (only alternative=64)
1                          ## [16] max number mismatches in barcodes
0                          ## [17] filter for adapters/primers (1 or 2=stricter)
32                         ## [18] minimum length of reads after adapter trim
2                          ## [19] ploidy: 1=haploid, 2=diploid, >2=(see docs)
1000                       ## [20] max cluster depth within samples
5                          ## [21] max Ns (uncalled bases) in consensus reads
5                          ## [22] max Hs (heterozygotes) in consensus reads
100,100                    ## [23] max SNPs in a locus (first,second for pairs)
5,100                      ## [24] max indels in a locus ("")
1,2,2,1                    ## [25] trim overhang (see docs)
0                          ## [26] hierarchical clustering (prob deprecated...)
denovo                     ## [27] clustering method (denovo, reference, hybrid)
./*.fa                     ## [28] reference genome file
lpn                        ## [29] output formats (see docs)
------ optional: list group/clade assignments below this line (see docs) ---
""" % version

    with open("params.txt", 'w') as paramsfile:
        print(paramstr, file=paramsfile)

if __name__ == "__main__":
    write_params("test")

