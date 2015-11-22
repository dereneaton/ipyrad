#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+ 


def main():
    import ipyrad as ip
    from ipyrad.core.assembly import Assembly
    """ a toy func """
    xpoint = 3
    print(xpoint, "CLI not ready yet")
    #TEST = Assembly("test-refseq")
    TEST = ip.load_assembly("/tmp/ipyrad-test/test-refseq.assembly")
    TEST.set_params(27, "hybrid")
    #TEST.set_params(1, "/tmp/ipyrad-test")
    #TEST.set_params(2, "./tests/data/sim_rad_test_R1_.fastq.gz")
    #TEST.set_params(3, "./tests/data/sim_rad_test_barcodes.txt")
    #TEST.set_params(28, "/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa" )
    #TEST.step1()
    #TEST.step2()
    TEST.step3(["1A_0"], preview=True, force=True)
    print(TEST.stats)

if __name__ == "__main__": 
    main()
