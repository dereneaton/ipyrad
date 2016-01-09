#!/usr/bin/env python2.7

""" Apply filters and write output files """

from __future__ import print_function

import h5py
import os
from ipyrad.file_conversion import *

import logging
LOGGER = logging.getLogger(__name__)

## List of all possible output formats. This is global because it's
## referenced by assembly.py and also paramsinfo. Easier to have it
## centralized.
OUTPUT_FORMATS = ['alleles', 'phy', 'nex', 'snps', 'vcf', 'usnps',
                  'str', 'geno', 'treemix', 'migrate', 'gphocs']


def run(data, samples, force, ipyclient):
    """ Check all samples requested have been clustered (state=6), 
    make output directory, then create the requested outfiles.
    """
    if any([i.stats.state <= 5 for i in samples]):
        print("  Step 7: Not all samples are aligned.")
        print("  Here are states for all the samples requested:")
        for i in samples:
            print("\t{}\t=\t{}".format(i.name, str(i.stats.state)))
        print("  All samples should be in state 6 for writing outfiles. "\
               +"Try rerunning step6()")
        ## TODO: Bail out here? Probably not good to write out data
        ## if all the requested samples don't exist.

    ## prepare dirs
    data.dirs.outfiles = os.path.join(data.dirs.working, "outfiles")
    if not os.path.exists(data.dirs.outfiles):
        os.mkdir(data.dirs.outfiles)

    LOGGER.info("Applying filters")
    ## Apply filters to supercatg and superhdf5 and write vcf
    filter_all_clusters( data, samples, ipyclient )

    LOGGER.info("Make .loci from filtered .vcf")
    ## Make .loci from the filtered vcf
    vcf2loci.make(data, samples, force)

    LOGGER.info("Convert .loci to all requested output file formats")
    ## Make all requested outfiles from the filtered .loci file
    make_outfiles(data, samples, force)


def filter_all_clusters( data, samples, ipyclient ):
    """ Read in the catclust.gz aligned clusters and the HDF5 supercatg
    database. Run through and filter each cluster

    Modelled on cluster_within.multi_muscle_align
    """
    ## create loadbalanced ipyclient
    lbview = ipyclient.load_balanced_view()

    ## Find optim cluster size
    ## TODO: Improve this, it's naive.
    optim = np.mean(data.stats["reads_consens"])/10

    ## Split dataset into chunks
    ## tmpnames holds the list of tmpfile names so we can rejoin them at the end
    tmpnames = []
    try:
        clustfile = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
        with gzip.open(clustfile, 'rb') as clustio:

            ## write optim clusters to each tmp file
            inclusts = iter(clustio.read().strip().split("//\n//\n"))
            grabchunk = list(itertools.islice(inclusts, optim))
            while grabchunk:
                with tempfile.NamedTemporaryFile('w+b',
                                                 delete=False,
                                                 dir=data.dirs.outfiles,
                                                 prefix=data.name+"_",
                                                 suffix='.chunk') as out:
                    out.write("//\n//\n".join(grabchunk))
                tmpnames.append(out.name)
                grabchunk = list(itertools.islice(inclusts, optim))

        ## Distribute chunks to parallel filter processes

        ## create job queue
        submitted_args = []
        for fname in tmpnames:
            submitted_args.append([data, samples, fname])

        ## run filter_stacks on all tmp files            
        results = lbview.map_async(filter_stacks, submitted_args)
        results.get()

        ## Concatenate filtered tmp files into final vcf/loci
        outvcf = os.path.join(data.dirs.outfiles, data.name+".full.vcf")
        make_vcfheader(data, samples, outvcf)

        with open(outvcf, 'wb') as out:
            for fname in tmpnames:
                with open(fname) as infile:
                    out.write(infile.read())

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        ## still delete tmpfiles if job was interrupted
        for fname in tmpnames:
            if os.path.exists(fname):
                os.remove(fname)
        del lbview

def filter_stacks(data, samples, fname):
    """ Filter one chunk of stacks and write out .tmp vcf/loci files
    This function runs in parallel, reads in a chunk of the stacks of reads,
    applies user specified filters, and writes out a tmp vcf style file of the results.
    """

    ## Read in chunk file
    try:
        with open(fname) as infile:
            loci = infile.read().split("//\n//\n")

    ## Apply various filters
    loci = filter_excludes(data, loci)
    loci = filter_minsamp(data, loci)
    loci = filter_maxSNP(data, loci)
    loci = filter_maxhet(data, loci)
    loci = filter_maxindels(data, loci)

    ## Write out .tmp vcf

def filter_excludes(data, loci):
    """ Remove excludes and outgroups
    """

def filter_minsamp(data, loci):
    """ Filter minimum # of samples per locus
    """
    # data.paramsdict["minsamp"]

    return loci

def filter_maxSNP(data, loci):
    """ Filter max # of SNPs per locus
    """
    # data.paramsdict["max_SNPs_locus"]

    return loci

def filter_maxhet(data, loci):
    """ Filter max shared heterozygosity per locus
    """

    # data.paramsdict["max_shared_heterozygosity"]

    return loci

def filter_maxindels(data, loci):
    """ Filter max # of indels per locus
    """
    # data.paramsdict["max_Indels_locus"]

    return loci


def loci_from_unfilteredvcf(data, samples, force):
    """ Read in the unfiltered vcf and supercatg from step6. Apply filters for
    coverage, heterozygosity, number of snps, etc. Write out .loci to output 
    directory """
    ## get unfiltered vcf handle
    unfiltered_vcf = os.path.join(data.dirs.consens, data.name+".vcf")

    supercatg = h5py.File(data.database, 'r')

    # Do filtering: max_shared_heterozygosity, minsamp, maxSNP, etc.
    finalfilter(data, samples, supercatg, unfiltered_vcf)

    ## Write out .loci
    locifile = os.path.join(data.dirs.outfiles, data.name+".loci")



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
