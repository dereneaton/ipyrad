#!/usr/bin/env python2

import logging
LOGGER = logging.getLogger(__name__)

def make( data, samples ):
    """ Convert vcf from step6 to .loci format to facilitate downstream format conversion """
    invcffile   =  os.path.join( data.dirs.consens, data.name+".vcf" )
    outlocifile  =  os.path.join( data.dirs.outfiles, data.name+".loci" )

    import( invcffile, outlocifile )


def import( vcffile, locifile ):
    """ Function for importing a vcf file into loci format. Arguments
        are the input vcffile and the loci file to write out. """

    try:
        ## Get names of all individuals in the vcf
        with open( invcffile, 'r' ) as invcf:
            for line in invcf:
                if line.split()[0] == "#CHROM":
                    ## This is maybe a little clever. The names in the vcf are everything after
                    ## the "FORMAT" column, so find that index, then slice everything after it.
                    names_col = line.split().index( "FORMAT" ) + 1
                    names = line.split()[ names_col:]
                    LOGGER.debug( "Got names - %s", names )
                    break

    ## Get the column to start reading at
