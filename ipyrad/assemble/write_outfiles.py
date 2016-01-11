#!/usr/bin/env python2.7

""" Apply filters and write output files. The basic body plan of this
code is as follows:
  * Read in the final aligned clusters file
  * Make sizeable chunks of loci
  * Distribute these chunks to parallel filters
  * Combine output of parallel filters into final loci file
  * Write out the output in full vcf format
"""

from __future__ import print_function

import h5py
import os
import gzip
import tempfile
import itertools
import numpy as np
from collections import Counter
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
    #vcf2loci.make(data, samples)

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

        ## Gather tmp loci files
        outloci = os.path.join(data.dirs.outfiles, data.name+".loci")
        with open(outloci, 'wb') as out:
            for fname in tmpnames:
                with open(fname.replace("chunk", "loci")) as infile:
                    out.write(infile.read())

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        ## still delete tmpfiles if job was interrupted
        ## fname - catclust.gz tmp files
        ## loc_name - tmp loci chunks
        for fname in tmpnames:
            if os.path.exists(fname):
                os.remove(fname)
            loc_name = fname.replace("chunk", "loci")
            if os.path.exists(loc_name):
                os.remove(loc_name)
        del lbview

def filter_stacks(args):
    """ Filter one chunk of stacks and write out .tmp vcf/loci files
    This function runs in parallel, reads in a chunk of the stacks of reads,
    applies user specified filters, and writes out a tmp vcf style file of the results.

    The design of the filtering steps intentionally sacrifices some performance
    for an increase in readability, and extensibility. Calling multiple filter
    functions ends up running through the sequences per stack several times, 
    but I felt this design made more sense, and also will easily allow us to
    add more filters in the future.
    """
    data, samples, fname = args

    ##mbda x: "_".join(x[0]. Read in chunk file
    ## info about filtered loci do we want?
    try:
        with open(fname) as infile:
            loci = infile.read().split("//\n//\n")

        ## Filter out the all the extraneous stuff after the sample name
        ## Locus names are of the form "samplename_4532", so we split on the
        ## underscore and return everything but the last element
        ## There's a better way to handle this by testing substring of data.sample names
        ##
        ## Converts each locus into a list of piars for name and sequence
        for i, loc in enumerate(loci):
            seqs = itertools.izip(*[iter(loc.split())]*2)
            loci[i] = map(lambda x: ("_".join(x[0].split("_")[:-1]), x[1]), seqs)
            ## Old way
            #loci[i] = "".join(map(lambda x: "_".join(x[0].split("_")[:-1]) + "\n" + x[1] + "\n", seqs))

        ## Apply various filters
        loci = filter_excludes(data, loci)
        #LOGGER.debug("loci out after excludes - {}".format(loci))
        loci = filter_duplicates(data, loci)
        #LOGGER.debug("loci out after dupes - {}".format(loci))
        loci = filter_minsamp(data, loci)
        #LOGGER.debug("loci out after minsamp - {}".format(loci[0]))
        loci = filter_maxhet(data, loci)
        #LOGGER.debug("loci out after maxhet - {}".format(loci))
        loci = filter_maxindels(data, loci)
        #LOGGER.debug("loci out after indels - {}".format(loci))
        loci = filter_maxSNP(data, loci)
        #LOGGER.debug("loci out after maxSNP - {}".format(loci[0]))

    except Exception as e:
        ## Do something real here
        LOGGER.warn(e)
        raise

    ## Write out .tmp loci
    write_tmp_loci(data, loci, fname)

    ## Write out .tmp vcf
    #write_tmp_vcf(data, loci, fname)

def filter_excludes(data, loci):
    """ Remove excludes and outgroups
    """
    ## Get all the samples to exclude. The 'or' conditional is a hack to account
    ## for the fact that paramsdict may either be an empty string or a list of
    ## names to exclude. List and "" don't append right, so if you do this 'or'
    ## and either of the lists is empty you'll get ["", "1B_0"]
    excludes = (data.paramsdict["excludes"] or [""]) \
                + (data.paramsdict["outgroups"] or [""])
    LOGGER.info("Excluding these individuals - {}".format(excludes))

    count = 0
    for i, loc in enumerate(loci):

        ## Get the count of excludes for logging/stats
        count += len(filter(lambda x: x[0] in excludes, loc))
 
        ## Actually filter the little buggers
        loci[i] = filter(lambda x: x[0] not in excludes, loc)

    LOGGER.info("Filterered exclude/outgroup sequences - {}".format(count))
    return loci

def filter_duplicates(data, loci):
    """ Don't allow multiple sequences from the same individual in a locus.
    """
    count = 0
    for i, loc in enumerate(loci):
        names = [x[0] for x in loc]

        ## set() dereplicates duplicate elements, so if len() is longer
        ## than set(len()) it means there was a dupe
        if len(names) > len(set(names)):
            loc = []
            count +=1 

    loci = filter(lambda x: x != [], loci)

    LOGGER.info("Filterered duplicates - {}".format(count))
    return loci


def filter_minsamp(data, loci):
    """ Filter minimum # of samples per locus
    """
    minsamp = data.paramsdict["minsamp"]
    
    count = 0
    for i, loc in enumerate(loci):
        nseqs = len([x[0] for x in loc])
        if nseqs < minsamp:
            loci[i] = []
            count += 1

    ## Filter out the empty lists left in the wake
    loci = filter(lambda x: x != [], loci)

    LOGGER.info("Filterered minsamp - {}".format(count))

    return loci


def filter_maxSNP(data, loci):
    """ Filter max # of SNPs per locus. Do R1 and R2 separately if PE. Also generate
    and append the snpsite line for the .loci format.
    """
    maxSNPS1 = data.paramsdict["max_SNPs_locus"]

    count = 0
    for i, loc in enumerate(loci):

## PE is broken for now, or it just does it as one big locus
##        if "pair" in data.paramsdict["datatype"]:
##            ## Set the values of max SNPs for R1 and R2 potentially separately. If only 
##            ## one value of max_SNPs_locus is set in paramsdict then use this value for R1 and R2
##            try:
##                maxSNPS2 = data.paramsdict["max_SNPS_locus"][1]
##            except AttributeError as e:
##                maxSNPS2 = maxSNPS1
##
##          loc2 = loc.split("SSSS")
##            nsnps, snpsites = count_snps(loc)

        ## Just pass in the list of sequences
        ## return is the # of snps and a .loci formatted string of snp positions
        seqs = [x[1] for x in loc]
        nsnps, snpsites = count_snps(seqs)

        LOGGER.debug("snpsites - {}".format(snpsites))
        if nsnps > maxSNPS1:
            loci[i] = []
            count += 1

        else:
            loci[i].append(snpsites)

    ## Filter out the empty lists left in the wake
    loci = filter(lambda x: x != [], loci)

    LOGGER.info("Filterered maxSNPs - {}".format(count))
    return loci


def filter_maxhet(data, loci):
    """ Filter max shared heterozygosity per locus
    """
    ## TODO: Figure out how to handle the 'SSSS' PE separator,
    ## cuz right now this is going to filter out all non merged
    ## PE reads....

    ## maxHet can be either an absolute value or a proportion.
    ## Since loci can be of different lengths we have to set
    ## it differently for each locus
    maxHet = data.paramsdict["max_shared_heterozygosity"]

    count = 0
    for i, loc in enumerate(loci):

        ## Pass the list of sequences to count_polymorphic_sites,
        ## which returns a list of counts per sequence.
        seqs = [x[1] for x in loc]
        nhets = count_shared_polymorphisms(seqs)

        ## Set the right value for maxHet
        if isinstance(maxHet, float):
            maxHet = int(maxHet * len(loc))

        ## Filter if any sequence has too many heterozygous sites
        if max(nhets) > maxHet:
            count += 1
            loci[i] = []

    ## Filter out the empty lists left in the wake
    loci = filter(lambda x: x != [], loci)

    LOGGER.info("Filterered max shared heterozygosity- {}".format(count))
    return loci

def filter_maxindels(data, loci):
    """ Filter max # of indels per locus
    """

    ## TODO: Fix this to evaluate PE reads individually
    maxIndel = data.paramsdict["max_Indels_locus"]

    count = 0
    for i, loc in enumerate(loci):
        ## Get just the sequences
        seqs = [x[1] for x in loci]
        ## Make counters for each sequence
        counts = [Counter(i) for i in seqs]
        ## If any sequence has more than the allowable # of indels then toss it
        if any([x["-"] > maxIndel for x in counts]):
            loci[i] = []
            count += 1

    ## Filter out the empty lists left in the wake
    loci = filter(lambda x: x != [], loci)

    LOGGER.info("Filterered max indels- {}".format(count))
    return loci

## This isn't being used
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

## Utility subfunctions
def count_shared_polymorphisms(seqs):
    """ Count the number of shared polymorphic sites at every base in a locus. If any 
    base contains too may shared polymorphisms (paramsdict[max_shared_heterozygosity])
    then we'll throw out the whole locus.
    Input is a list of sequences, output is a list of ints representing counts
    of shared polyorphisms per base
    """
    ## Transpose the lost of seqs, so we now have a list of bases at each locus.
    ## Stacks of sequences should always be the same length, but if they aren't
    ## we'll protect by filling with N's, rather than truncating, which is what zip()
    ## would do.
    ## Makes a list of counters for each stack of bases of this locus
    stacks = [Counter(i) for i in itertools.izip_longest(*seqs, fillvalue="N")]

    ## This is maybe 20% too clever, but i wanted to do it with list comprehension,
    ## so here we are. Read this as "At each base get the max number of counts of any
    ## ambiguity code". Below y.items is a counter for each position, and x[0] is the counter
    ## key, x[1] is the count for that key. Another stupid thing is that max() in python 2.7
    ## doesn't handle empty lists, so there's this trick max(mylist or [0]), empty list
    ## evaluates as false, so it effectively defaults max({}) = 0.
    max_counts = [max([x[1] for x in y.items() if x[0] in "RYSWKM"] or [0]) for y in stacks]

    return max_counts

def count_snps(seqs):
    count = 0

    ## As long as we're ripping through looking for snps, keep track of the snp array
    ## to output to the .loci file
    snpsite = [" "]*len(seqs[0])

    ## Transpose to make stacks per base, same as count_shared_polymorphisms
    ## so see that function for more explicit documentation.
    stacks = [Counter(i) for i in itertools.izip_longest(*seqs, fillvalue="N")]

    ## At each stack, if the length of the counter is > 1 then its a snp
    for i, stack in enumerate(stacks):
        if len(stack) > 1:
            LOGGER.debug("Got a snp {}".format(stack))
            ## Count the number of keys that occur once, these are autapomorphies
            autapomorphies = stack.values().count(1)

            LOGGER.debug("Autapomorphies - {}".format(autapomorphies))
            ## If the number of autapomporphies is one less than the total length
            ## of the counter then this is still a parsimony uninformative site.                
            if autapomorphies == (len(stack) - 1):
                snpsite[i] = "-"
            ## If not then it is a parsimony informative site
            else:
                snpsite[i] = "*"

    ## When you're done, rip through the snpsite string and count snp markers "*" & "-"
    nsnps = sum(1 for x in snpsite if x in "*-")

    ## Just return the snpsite line
    snps = ("//", "".join(snpsite)+"|")
    #seqs.append(("//", "".join(snpsite)+"|"))

    return nsnps, snps

## File output subfunctions
def write_tmp_loci(data, loci, fname):
    """ Write out the filtered chunk to a tmp file which will be collated by the
    top level run() thread"""

    ## Get longest sample name for pretty printing
    longname_len = max(len(x) for x in data.samples.keys())
    ## Padding distance between name and seq this could be a hackers only param
    name_padding = 5

    with open(fname.replace("chunk","loci"), 'w') as outfile:
        for loc in loci:
            for seq in loc:

                ## Backwards compatibility with .loci format which doesn't have 
                ## leading '>' on the snpsites line
                if seq[0] == "//":
                    name = seq[0] 
                else:
                    name = ">" + seq[0]

                name +=  " " * (longname_len - len(name)+ name_padding)
                outfile.write(name + seq[1] +"\n")

def make_vcfheader(data, samples, outvcf):
    LOGGER.debug("Entering make_vcfheader()")

if __name__ == "__main__":
    import ipyrad as ip
    TESTFILE = "/tmp/ipyrad-test/test-refseq.assembly"
    TESTER = ip.load.load_assembly(TESTFILE)
#    TESTER = ip.core.assembly.Assembly( "/tmp/ipyrad-test-pair" )
    TESTER.set_params( "output_formats", "vcf,snps" )
    TESTER.get_params()
    TESTER.set_params( "output_formats", "*" )
    TESTER.get_params()
    TESTER.step7()
