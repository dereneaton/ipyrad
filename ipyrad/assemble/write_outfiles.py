#!/usr/bin/env python2.7

""" Apply filters and write output files """

from __future__ import print_function

from ipyrad.file_conversion import *

import logging
LOGGER = logging.getLogger(__name__)

## List of all possible output formats.
output_formats = ['alleles', 'phy', 'nex', 'snps', 'vcf', 'usnps',
                     'str', 'geno', 'treemix', 'migrate', 'gphocs']

def run( data, samples, force, ipyclient ):
    """ Check all samples requested have been clustered (state=6), 
    make output directory, then create the requested outfiles.
    """
    if any([i.stats.state <=5 for i in samples]):
        print("  Step 7: Not all samples are aligned.")
        print("  Here are states for all the samples requested:")
        for i in samples:
            print("\t{}\t=\t{}".format( i.name, str(i.stats.state ) ) )
        print("  All samples should be in state 6 for writing outfiles. Try rerunning step6()")

    ## prepare dirs
    data.dirs.outfiles = os.path.join(data.dirs.working, data.name+"_outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    ## Make the .loci file from the vcf generated in step6()
#    loci_from_vcf( data, samples, force )

    ## Make all requested outfiles
    make_outfiles( data, samples, force )

def make_outfiles( data, samples, force ):
    """ Get desired formats from paramsdict and write files to outfiles directory """

    for filetype in data.paramsdict["output_formats"]:
        print( "Doing - ", filetype )


if __name__ == "__main__":
    import ipyrad as ip
#    TESTFILE = "/tmp/ipyrad-test/test-refseq.assembly"
#    TESTER = ip.load.load_assembly(TESTFILE)
    TESTER = ip.core.assembly.Assembly( "test" )
    TESTER.set_params( "output_formats", "vcf,snps" )
    TESTER.get_params()
    TESTER.set_params( "output_formats", "*" )
    TESTER.get_params()
    TESTER.step7()
