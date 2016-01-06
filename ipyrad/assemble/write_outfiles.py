#!/usr/bin/env python2.7

""" Apply filters and write output files """

from __future__ import print_function

import h5py
from ipyrad.file_conversion import *

import logging
LOGGER = logging.getLogger(__name__)

## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized.
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
    data.dirs.outfiles = os.path.join(data.dirs.working, "outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    ## Make the .loci file from the vcf generated in step6()
    loci_from_unfilteredvcf( data, samples, force )

    ## Make all requested outfiles
    make_outfiles( data, samples, force )



def loci_from_unfilteredvcf( data, samples, force ):
    """ Read in the unfiltered vcf and supercatg from step6. Apply filters for coverage,
    heterozygosity, number of snps, etc. Write out .loci to output directory """

    unfiltered_vcf = os.path.join(data.dirs.consens, data.name+".vcf")

    supercatg = h5py.File(data.database, 'r')

    # Do filtering: max_shared_heterozygosity, minsamp, maxSNP, etc.

    ## Write out .loci
    locifile = os.path.join( data.dirs.outfiles, data.name+".loci" )



def make_outfiles( data, samples, force ):
    """ Get desired formats from paramsdict and write files to outfiles directory """

    ## Read in the input .loci file that gets transformed into all output formats
    locifile = os.path.join( data.dirs.outfiles, data.name+".loci" )


    for filetype in data.paramsdict["output_formats"]:
        LOGGER.info( "Doing - ", filetype )

        # phy & nex come from loci2phynex
        if filetype in ["phy", "nex"]:
            filetype = "phynex"
        # All these file types come from loci2SNP
        elif filetype in ["snps", "usnps", "str", "geno"]:
            filetype = "SNP"

        ## Everything else has its own loci2*.py conversion file.
        ## Get the correct module for this filetype
        ## globals() here gets the module name, and then we can call
        ## the .make() function. This is a little tricky.
        format_module = globals()[ "loci2"+filetype ]

        ## do the call to make the new file format
        format_module.make( data, samples )

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
