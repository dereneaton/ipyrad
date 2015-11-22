#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+ 
import optparse


def main():
    import ipyrad as ip
    from ipyrad.core.assembly import Assembly
    """ a toy func """
    ## create a default Assembly object to get the cluster running
    TEST = Assembly("test-refseq")

    ##
    #TEST = ip.load_assembly("/tmp/ipyrad-test/test-refseq.assembly")
    TEST.set_params(27, "hybrid")
    TEST.set_params(1, "/tmp/ipyrad-test")
    TEST.set_params(2, "./tests/data/sim_rad_test_R1_.fastq.gz")
    TEST.set_params(3, "./tests/data/sim_rad_test_barcodes.txt")
    TEST.set_params(28, "/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa" )
    TEST.step1()
    TEST.step2()
    TEST.step3(["1A_0"], preview=True, force=True)
    print(TEST.stats)

def parse_params():
    print("not yet")

def parse_command_line():
    """ Parse CLI. Only three options now. """
    ## -p  choose params file
    ## -s  subselect steps
    ## -n  create new params file
    print("not yet")


if __name__ == "__main__": 
    
    ## parse params file input
    parse_params()
    main()
    ## 
