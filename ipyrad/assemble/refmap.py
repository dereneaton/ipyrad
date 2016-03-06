#!/usr/bin/env python2.7

""" 
Aligns reads to a reference genome. Returns unaligned reads to the de novo
pipeline. Aligned reads get concatenated to the *clustS.gz files at the 
end of step3
"""

from __future__ import print_function

import os
import re
import sys
import gzip
import copy
import shutil
import tempfile
import itertools
import subprocess
import numpy as np
import ipyrad
from .util import *
from ipyrad.assemble.rawedit import comp

import logging
LOGGER = logging.getLogger(__name__)

# pylint: disable=W0142


def index_reference_sequence(data, force=False):
    """ Attempt to index the reference sequence. This is a little naive
    in that it'll actually _try_ do to the reference every time, but it's
    quick about giving up if it detects the indices already exist. You could
    also test for existence of both index files, but i'm choosing to just let
    smalt do that for us ;) """

    refseq_file = data.paramsdict['reference_sequence']

    #TODO: Here test if the indices exist already
    # These are smalt specific index files. We don't ever reference
    # them directly except here to make sure they exist, so we don't need
    # to keep them around.
    index_sma = refseq_file+".sma"
    index_smi = refseq_file+".smi"

    if os.path.isfile(index_sma) and os.path.isfile(index_smi):
        if force:
            print("    Force reindexing of reference sequence")
        else:
            print("    Reference sequence index exists")
            return

    msg = "    **************************************************************\n"\
        + "    Indexing reference sequence. This only needs to be done once.\n"\
        + "    This should take less than 10 minutes.\n"\
        + "    **************************************************************"

    LOGGER.info(msg)
    if data._headers:
        print(msg)

    cmd = ipyrad.bins.smalt\
        + " index "\
        + " -k "+ str(data._hackersonly["smalt_index_wordlen"])\
        + " " + refseq_file + " " + refseq_file

    LOGGER.debug(cmd)
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

    if data._headers:
        print("    Done indexing reference sequence")



def mapreads(args):
    """ Attempt to map reads to reference sequence. This reads in the 
    samples.files.edits .fasta files, and attempts to map each read to the 
    reference sequences. Unmapped reads are dropped right back in the de novo 
    pipeline. Reads that map successfully are processed and pushed downstream 
    and joined with the rest of the data post musle_align. The read mapping 
    produces a sam file, unmapped reads are pulled out of this and dropped back 
    in place of the edits (.fasta) file. The raw edits .fasta file is moved to 
    .<sample>.fasta to hide it in case the mapping screws up and we need to 
    roll-back. Mapped reads stay in the sam file and are pulled out of the 
    pileup later."""

    ## get args
    data, sample, noreverse, nthreads = args
    LOGGER.debug("Entering mapreads(): %s %s %s", \
                                    sample.name, noreverse, nthreads)

    ## Test edits actually exist
    if sample.files.edits == []:
        LOGGER.debug("Sample edits empty. Rerun step2(force=True)")
        LOGGER.debug("Sample files - %s", sample.files)
        sys.exit("Sample edits empty. Rerun step2(force=True)")

    ## Set the edited fastq to align for this individual, we set this here
    ## so we can overwrite it with a truncated version if we're in preview mode.
    ## Set the default sample_fastq file to align as the entire original
    ##
    ## refmap_init copies the file paths for the full fastq files to
    ## files.edits_refmap_clean. It makes new temp files at sample.files.edits
    ## because this is where we'll dump the unmapped reads for downstream denovo.
    ## files.edits gets put back to point to the original fastq files at the end
    ## of split_among_processors in cluster_within.
    sample_fastq = [sample.files.edits_refmap_clean[0][0]]

    ## If pair append R2 to sample_fastq
    if 'pair' in data.paramsdict["datatype"]:
        sample_fastq.append(sample.files.edits_refmap_clean[0][1])

    ## Files we'll use during reference sequence mapping
    ##
    ## samhandle - Raw output from smalt reference mapping
    ## unmapped_bamhandle - bam file containing only unmapped reads
    ## mapped_bamhandle - bam file with only successfully mapped reads
    ##    In the case of paired end, we only consider pairs where both
    ##    reads map successfully
    ## sorted_bamhandle - Sorted bam file
    ##    This is a precursor for mpileup (bam must be sorted for piluep to work)
    ##    I don't think this is strictly necessary any more....
    ##    TODO: Delete sorting?
    ## unmapped_fastq_handle - The final output file of this function
    ##    writes out unmapped reads to the 'edits' directory as .fq
    ##    which is what downstream analysis expects
    samhandle = os.path.join(data.dirs.refmapping, sample.name+".sam")
    sample.files.unmapped_reads = unmapped_bamhandle = \
                os.path.join(data.dirs.refmapping, sample.name+"-unmapped.bam")
    sample.files.mapped_reads = mapped_bamhandle = \
                os.path.join(data.dirs.refmapping, sample.name+"-mapped.bam")
    
    sorted_mapped_bamhandle = sample.files["mapped_reads"]

    unmapped_fastq_handle = sample.files.edits[0][0]

    ## If PE set the output file path for R2
    if 'pair' in data.paramsdict["datatype"]:
        unmapped_fastq_handle_R2 = sample.files.edits[0][1]


    ## get smalt call string
    ##  -f sam       : Output as sam format, tried :clip: to hard mask output but it fsck
    ##                 and shreds the unmapped reads (outputs empty fq)
    ##  -l [pe,mp,pp]: If paired end select the orientation of each read
    ##  -n #         : Number of threads to use
    ##  -x           : Perform a more exhaustive search
    ##  -c #         : fraction of the query read length that must be covered
    ##  -o           : output file
    ##               : Reference sequence
    ##               : Input file(s), in a list. One for forward and one for reverse reads.

    ## Wrap the whole thing in a try so we can clean up if something fails or user cancelse
    ## This is pretty monolithic, it could get broken up into subroutines
    try:

        ## We are assuming Illumina paired end reads for the orientation of mate pairs.
        ## (orientation: ---> <----).
        if 'pair' in data.paramsdict["datatype"]:
            pairtype = " -l pe "
        else:
            pairtype = " "
    
        cmd = ipyrad.bins.smalt+\
            " map -f sam -n " + str(nthreads) +\
            pairtype+\
            " -x -c " + str(data.paramsdict['clust_threshold'])+\
            " -o " + samhandle +\
            " " + data.paramsdict['reference_sequence'] +\
            " " + " ".join( sample_fastq )
    
        LOGGER.debug( "%s", cmd )
        try:
            subprocess.check_call(cmd, shell=True,
                                 stderr=subprocess.STDOUT,
                                 stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as inst :
            ## Handle error in outside try statement
            raise
    
        ## Get the mapped and unmapped reads from the sam. For PE both reads must map
        ## successfully in order to qualify.
        ##   1) Get mapped reads and convert to bam
        ##   2) Sort them and save the path to the bam to sample.files.mapped_reads
        ## The mapped reads are synced back into the pipeline downstream, after
        ## muscle aligns the umapped reads.
        ##
        ## samtools view arguments
        ##   -b = write to .bam
        ##   -F = Select all reads that DON'T have this flag. 
        ##         0x4 (segment unmapped)
        ##   -U = Write out all reads that don't pass the -F filter (all unmapped reads
        ##        go to this file.
        ##        <TODO>: This is deeply hackish right now, it will need some
        ##                serious thinking to make this work for PE, etc.
        sam_filter_flag = " -F 0x4 "
    
        if 'pair' in data.paramsdict["datatype"]:
            ## Additionally for PE only output read pairs that both align
            sam_filter_flag += " -f 0x2 "
    
        cmd = ipyrad.bins.samtools+\
            " view -b"+\
                sam_filter_flag+\
                " -U " + unmapped_bamhandle+\
                " " + samhandle+\
                " > " + mapped_bamhandle
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ## Step 2 sort mapped bam
        ##   -T = Temporary file name, this is required by samtools, ignore it
        ##        Here we hack it to be samhandle.tmp cuz samtools cleans it up
        ##   -O = Output file format, in this case bam
        ##   -o = Output file name
        cmd = ipyrad.bins.samtools+\
            " sort -T "+samhandle+".tmp" +\
            " -O bam "+mapped_bamhandle+\
            " -o "+sorted_mapped_bamhandle
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ## Step 3 index mapped reads
        ## Samtools pileup needs the bam to be indexed
        ## No arguments, a very simple function. It writes the index to 
        ## a default location
        cmd = ipyrad.bins.samtools+\
            " index " + mapped_bamhandle
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    

        ##############################################
        ## Do unmapped
        ## Output the unmapped reads to the original
        ## sample.edits fq path.
        ##############################################
    
        outfiles = [unmapped_fastq_handle]
        if 'pair' in data.paramsdict["datatype"]:
            outfiles.append(unmapped_fastq_handle_R2)
            outflags = " -1 " + outfiles[0]+\
                       " -2 " + outfiles[1]
        else:
            outflags = " -0 " + outfiles[0]
    
        cmd = ipyrad.bins.samtools+\
            " bam2fq " + outflags+\
                " " + unmapped_bamhandle
    
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ## This is the end of processing for each sample. Stats
        ## are appended to the sample for mapped and unmapped reads 
        ## during cluster_within.cleanup()

    except subprocess.CalledProcessError:
            LOGGER.error("Error in reference mapping. Try copy/pasting and running this "+\
                            "command by hand:\n\t%s", cmd)
            ## TODO: Maybe make this not bail out hard?
            sys.exit("Error in reference mapping \n{}\n{}\n{}."\
                     .format(inst, subprocess.STDOUT, cmd))
            raise
    except KeyboardInterrupt:
            LOGGER.error(  "Handling user initiated shutdown.")
            sys.exit("Handling user initiated shutdown.")
    except Exception as inst:
            LOGGER.error(  "Error in reference mapping - {}".format(inst))
            raise
    finally:
        pass




def finalize_aligned_reads( data, sample, ipyclient ):
    """ Run bedtools to get all overlapping regions. Then split this
        list of regions and distribute it across a bunch of parallel threads.

        1) Run bedtools merge to get a list of all contiguous blocks of bases
        in the reference seqeunce where one or more of our reads overlap.
        The output will look like this:
            1       45230754        45230783
            1       74956568        74956596
            ...
            1       116202035       116202060

        The structure here is copied from multi_muscle_align. Each thread will
        run on 100 regions.
    """

    lbview = ipyclient.load_balanced_view()

    try:

        ## Regions is a giant list of 5-tuples, of which we're only really 
        ## interested in the first three: chrom, start, and end position.
        regions = bedtools_merge(data, sample)
        
        if len(regions) > 0:
            ## Empty array to hold our chunks
            tmp_chunks = []

            all_regions = iter(regions.strip().split("\n"))
            grabchunk = list(itertools.islice(all_regions, 100))
            while grabchunk:
                tmp_chunks.append(grabchunk)
                grabchunk = list(itertools.islice(all_regions, 100))

            submitted_args = []
            for chunk in tmp_chunks:
                # LOGGER.info("Whats this: %s", [data, sample, chunk])
                submitted_args.append([data, sample, chunk])

            ## run get_aligned_reads on all region chunks            
            results = lbview.map_async(get_aligned_reads, submitted_args)
            results.get()
        else:
            LOGGER.info("No reads mapped to reference sequence.")

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        pass



def get_aligned_reads( args ):
    """Pull aligned reads out of sorted mapped bam files and
    append them to the clustS.gz file so the fall into downstream analysis
    Here's the exact order of ops:
    
    1) Coming into this function we have sample.files.mapped_reads 
        as a sorted bam file, and a passed in list of regions to evaluate.
    2) Get all reads overlapping with each individual region.
    3) Write the aligned sequences out to a fasta file (vsearch doesn't 
       handle piping data) (IT WILL VERY SOON, actually).
    4) Call vsearch to derep reads
    5) Append to the clustS.gz file.

    ## The old way was a bit more laborious, and also didn't work great. Plus
    ## it had the disadvantage of not including bases that don't align to 
    ## the reference (they would get hard masked by the pileup command). 
    ## Get rid of all this when you get sick of looking at it...
    2) Samtools pileup in each region individually. This will build stacks
       of reads within each contiguous region.
    3) Decompile the mpileup format into fasta to get the raw sequences within
       each stack.
    """

    #reload(ipyrad.assemble.refmap)
    #import ipyrad.assemble.refmap
    data, sample, regions = args
    ## Keep track of all the derep'd fasta files per stack, we concatenate them
    ## all to the end of the clustS.gz file at the very end of the process
    derep_fasta_files = []
    aligned_seq_files = []
    reads_merged = 0

    # Wrap this in a try so we can clean up if it fails.

    try:
        ## For each identified region, build the pileup and write out the fasta
        for line in regions:
            LOGGER.info("line %s in region %s", line, regions)

            # Blank lines returned from bedtools make the rest of the 
            # pipeline angry. Filter them.
            if line == "":
                continue

            chrom, region_start, region_end = line.strip().split()[0:3]

            # Here aligned seqs is a list of files 1 for SE or 2 for PE
            args = [data, sample, chrom, region_start, region_end]
            aligned_seqs = bam_region_to_fasta(*args)

            ## If bam_region_to_fasta fails for some reason it'll return [], 
            ## in which case skip the rest of this. Normally happens if reads
            ## map successfully, but too far apart.
            if not aligned_seqs:
                continue

            # This whole block is getting routed around at this point. I'm not 
            # deleting it cuz i'm precious about it, and cuz if we ever decide 
            # to use pileup to call snps it could be useful. The next two 
            # functions generate pileups per region and then backtransform
            # pileup to aligned fasta.
            # pileup_file, read_labels = bam_to_pileup( data, sample, chrom, 
            # region_start, region_end )
            ## Test the return from bam_to_pileup. If mindepth isn't satisfied 
            # this function will return
            ## an empty string and empty list. In this case just bail out on 
            # this pileup, bad data.
            #if pileup_file == "":
            #    LOGGER.debug( "not enough depth at - %s, %s, %s", chrom, 
                #region_start, region_end )
            #    continue
            #aligned_seqs = mpileup_to_fasta( data, sample, pileup_file )

            #aligned_fasta_file = write_aligned_seqs_to_file( data, 
                #sample, aligned_seqs, read_labels )

            ## merge fastq pairs
            if 'pair' in data.paramsdict['datatype']:
                ## merge pairs that overlap and combine non-overlapping
                ## pairs into one merged file. merge_pairs takes the unmerged
                ## files list as an argument because we're reusing this code 
                ## in the refmap pipeline, trying to generalize.
                LOGGER.debug("Merging pairs - %s", sample.files)
                mergefile, nmerged = refmap_merge_pair(data, sample, aligned_seqs)
                
                ## Update the total number of merged pairs
                reads_merged += nmerged
                #sample.stats.reads_merged += nmerged
                sample.merged = 1
                aligned_fasta = mergefile
            else:
                ## If SE we don't need to merge, and the aligned fasta are 
                ## just the first element of the list returned above
                aligned_fasta = aligned_seqs[0]

            ## Post refmap derep brings mapped reads into alignment with the
            ## format expected downstream (clust.gz)
            derep_fasta = derep_and_sort(data, sample, aligned_fasta)

            ## Derep_fasta_files are merged for PE
            derep_fasta_files.append(derep_fasta)
            aligned_seq_files.append(aligned_seqs)
        append_clusters(data, sample, derep_fasta_files)

    except Exception as inst:
        LOGGER.warn("Caught exception inside get_aligned_reads - {}".format(inst))
        LOGGER.debug("Current contents of:\nderep_fasta_files {}\n"\
           "aligned_seq_files {}".format(derep_fasta_files, aligned_seq_files))
        raise
    finally:
        LOGGER.debug("Total merged reads for {} - {}".format(sample.name, reads_merged))
        sample.stats.reads_merged = reads_merged
        # Clean up all the tmp files
        # Be a little careful. Don't remove files if they're already gone :-/
        for i in derep_fasta_files:
            if os.path.isfile(i):
                os.remove(i)
        for j in aligned_seq_files:
            if os.path.isfile(j[0]):
                os.remove(j[0])
            if 'pair' in data.paramsdict['datatype']:
                if os.path.isfile(j[1]):
                    os.remove(j[1])


def refmap_merge_pair(data, sample, aligned_seqs):
    """This formats the data for the call to util/merge_pairs().
    merge_pairs() arguments and expectations changed enough to where
    it would be painful to perform these gymanstics in get_aligned_reads.
    This util function is a workaround."""
    sample.files.edits = [(aligned_seqs[0], aligned_seqs[1])]

    sample = merge_pairs(data, sample)

    return sample.files.merged, sample.stats.reads_merged
    

def bedtools_merge(data, sample):
    """ Get all contiguous genomic regions with one or more overlapping
    reads. This is the shell command we'll eventually run

        bedtools bamtobed -i 1A_0.sorted.bam | bedtools merge [-d 100]
            -i <input_bam>  :   specifies the input file to bed'ize
            -d <int>        :   For PE set max distance between reads
    """
    LOGGER.debug("Entering bedtools_merge: %s", sample.name)

    if 'pair' in data.paramsdict["datatype"]:
        ## check mean insert size for this sample
        ## and update hackersonly.max_inner_mate_distance
        ## if need be. This value controls how far apart 
        ## mate pairs can be to still be considered for 
        ## bedtools merging downstream
        cmd = ipyrad.bins.samtools+\
            " stats " + sample.files.mapped_reads + " | grep SN"
        LOGGER.debug("%s", cmd)
        ret = subprocess.check_output(cmd, shell=True,
                             stderr=subprocess.STDOUT)

        avg_insert = 0
        stdv_insert = 0
        avg_len = 0
        for line in ret.split("\n"):
            if "insert size average" in line:
                avg_insert = float(line.split(":")[-1].strip())
            elif "insert size standard deviation" in line:
                stdv_insert = float(line.split(":")[-1].strip())
            elif "average length" in line:
                avg_len = float(line.split(":")[-1].strip())

        LOGGER.debug("avg {} stdv {} avg_len {}".format(avg_insert, stdv_insert, avg_len))

        ## If all values return successfully set the max inner mate distance to 
        if all([avg_insert, stdv_insert, avg_len]):
            data._hackersonly["max_inner_mate_distance"] = avg_insert + (3 * stdv_insert) - (2 * avg_len)
        else:
            ## If something fsck then set a relatively conservative distance
            data._hackersonly["max_inner_mate_distance"] = 200
        LOGGER.debug("inner mate distance for {} - {}".format(sample.name,\
                    data._hackersonly["max_inner_mate_distance"]))

    ## Set the -d flag to tell bedtools how far apart to allow mate pairs.
    if 'pair' in data.paramsdict["datatype"]:
        bedtools_dflag = " -d " + str(data._hackersonly["max_inner_mate_distance"])
    else:
        bedtools_dflag = " "

    cmd = ipyrad.bins.bedtools+\
        " bamtobed "+\
        " -i " + sample.files.mapped_reads+\
        " | " +\
        ipyrad.bins.bedtools +\
        " merge "+\
        bedtools_dflag
    LOGGER.debug("%s", cmd)
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    ## Report the number of regions we're returning
    LOGGER.debug("bedtools_merge: Got # regions: %s", str(len(result.strip().split("\n"))))
    return result



def bam_region_to_fasta( data, sample, chrom, region_start, region_end):
    """ Take the chromosome position, and start and end bases and output sequences
    of all reads that overlap these sites. This is the command we're building:

        samtools view -b 1A_sorted.bam 1:116202035-116202060 | samtools bam2fq <options> -
                -b      : output bam format
                -0      : For SE, output all reads to this file
                -1/-2   : For PE, output first and second reads to different files
                -       : Tell samtools to read in from the pipe

    Write out the bam file, then use samtools bam2fq to write out the reads
    to individual files. Return the file name(s) for each read.
    """
    LOGGER.debug( "Entering bam_region_to_fasta: %s %s %s %s", sample.name, chrom, region_start, region_end )

    ## Return it as a list so we can handle SE and PE without a bunch of monkey business.
    outfiles = []

    tmp_outfile = sample.files.mapped_reads+"-"+str(chrom)+"-"+str(region_start)

    try:
        ## Make the samtools view command to output bam in the region of interest
        view_cmd = ipyrad.bins.samtools+\
            " view "+\
            " -b "+\
            sample.files.mapped_reads+\
            " " + chrom + ":" + region_start + "-" + region_end

        ## Set output files and flags for PE/SE
        outfiles = [ tmp_outfile+"R1.fq" ]
        if 'pair' in data.paramsdict["datatype"]:
            outfiles.append( tmp_outfile+"R2.fq" )
            outflags =  " -1 " + outfiles[0]+\
                        " -2 " + outfiles[1]
        else:
            outflags = " -0 " + outfiles[0]

        bam2fq_cmd = ipyrad.bins.samtools+\
            " bam2fq " + outflags + " - "

        cmd = " | ".join( (view_cmd, bam2fq_cmd) )
        LOGGER.debug( "%s", cmd )
        subprocess.call( cmd , shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)

    except Exception as inst:
        LOGGER.warn( inst )

    ## post-flight checks for consistency
    if 'pair' in data.paramsdict["datatype"]:
        len1 = sum(1 for line in open(outfiles[0]))
        len2 = sum(1 for line in open(outfiles[1]))
        if not len1 == len2:
            LOGGER.warn("Read files in region {}:{}-{} different lengths, probably"+\
                        "distance between reads exceeds data._hackersonly[`max_inner_mate_distance`]={}".format(\
                        chrom, region_start, region_end, data._hackersonly["max_inner_mate_distance"]))
            outfiles = []

    return outfiles



## This is the "old" way, where i was writing fasta files by hand rather than using
## bam2fq like makes more sense...
def bam_region_to_fasta2(data, sample, chrom, region_start, region_end):
    """ Take the chromosome position, and start and end bases and output 
    sequences of all reads that overlap these sites. This is the command we're 
    building:

        samtools view 1A_sorted.bam 1:116202035-116202060

    We also have to track the names of each read name (QNAME) in this region, 
    so we can reconstruct the fasta downstream. This command will output all
    overlapping reads in headerless sam format. QNAME is the first field and
    sequence data is the 9th.
    """
    LOGGER.debug("Entering bam_region_to_fasta: %s %s %s %s", \
                 sample.name, chrom, region_start, region_end)

    cmd = ipyrad.bins.samtools+\
        " view "+\
        sample.files.mapped_reads+\
        " " + chrom + ":" + region_start + "-" + region_end
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)

    # Use list comprehension to pull out the zeroth/ninth element of each row
    # TODO: Easy optimization here, here's probably a "smart" way to do this in one line.
    read_labels = [x.split("\t")[0] for x in result.strip("\n").split("\n")]
    sequence_data = [x.split("\t")[9] for x in result.strip("\n").split("\n")]

    return sequence_data, read_labels



def bam_to_pileup(data, sample, chrom, region_start, region_end):
    """ Take the chromosome position, and start and end bases and output a 
    pileup of all reads that overlap these sites. This is the command we're 
    building:

        samtools mpileup -f MusChr1.fa \
                         -r 1:116202035-116202060 \
                         -o out.pileup 1A_0.sorted.bam

    We also have to track the names of each read name (QNAME) in this region, 
    so we can reconstruct the fasta downstream. This command will output all
    overlapping reads in headerless sam format. QNAME is the first field.

        samtools view 1A_sorted.bam 1:116202035-116202060

    NB: This function is not currently used in the pipeline. It works good tho,
        so I'm keeping it around in case we need it in the future.
    """
    LOGGER.debug("Entering bam_to_pileup: %s %s %s %s", \
                 sample.name, chrom, region_start, region_end)

    ## make output directory for pileups
    ## These aren't really strictly necessary to keep around, but for
    ## dev and debugging it's good to be able to see what's going on.
    ## This dir should probably be cleaned up after execution.  
    ## This is actually stupid here, but it's just for testing
    ## TODO: Remove the function to keep pileups before shipping.
    ##
    ## Just replace all this shit with this:
    ##     pileup_file = data.dirs.refmapping+"/"+sample.name+"-"+\
    ##     region_start+".pileup"
    pileup_dir = os.path.join(data.dirs.refmapping, "pileups")
    if not os.path.exists(pileup_dir):
        os.makedirs(pileup_dir)

    pileup_file = pileup_dir+"/"+sample.name+"-"+region_start+".pileup"

    cmd = ipyrad.bins.samtools+\
        " view "+\
        sample.files.mapped_reads+\
        " " + chrom + ":" + region_start + "-" + region_end
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)

    # Use list comprehension to pull out the first element of each row,
    read_labels = [x.split("\t")[0] for x in result.strip("\n").split("\n")]
    
    # Test the number of reads overlapping in this region. If it's strictly
    # less than mindepth parameter, then just toss it out, save us some time.
    # I fuckin hate multiple returns, but in this case i'll make an exception.
    # I still might move this to the end, just on principle.
    #
    # We'll leave low depth clusters in for now. If you decide to get
    # rid of them just uncomment these two lines.
    # 
    #if len(read_labels) < data.paramsdict["mindepth_majrule"]:
    #    return( "", [] )

    cmd = ipyrad.bins.samtools+\
        " mpileup "+\
        " -f " + data.paramsdict['reference_sequence']+\
        " -r " + chrom+":"+region_start+"-"+region_end+\
        " -o " + pileup_file+\
        " " + sample.files.mapped_reads
    #    LOGGER.debug( "%s", cmd )
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)

    return pileup_file, read_labels




def mpileup_to_fasta(data, sample, pileup_file):
    """ Takes a pileup file and decompiles it to fasta. It is currently 
    "working" but there are some bugs. If you want to actually use this it'll 
    need some tlc.
    
    NB: This function is not currently used in the pipeline. It's a good idea 
    tho, so I'm keeping it around in case we need it in the future.
    """

    LOGGER.debug("Entering mpileup_to_fasta: %s %s", sample.name, pileup_file)

    with open(pileup_file, 'r') as pfile:
        pileup = []
        for line in pfile:
            dat = np.array(line.strip().split())

            # Handle the case where the number of aligned reads is 0, in which case there
            # is no 4th field for the sequence info. If this line has 0 reads, just skip it.
            if len(dat) == 4:
                continue
            refbase, count, seqs = dat[[2,3,4]]
            pileup.append([refbase, count, seqs])
    
        # TODO: Fix this, max len may not be the len of the max depth in the pileup
        nreads = max([int(x[1]) for x in pileup])
        
        # Remove the ^ and the mapping qual score indicating the start of a new read segment
        # because it's not useful for us. We will track the construction of the alignments
        # in a structured np array.
        # Fill in all missing values with N, it'll make things less messy down the road
        for line in pileup:
            line[2] = re.sub(r"\^.", "", line[2])
            line[2] = line[2] + ((nreads - len(line[2])) * "N")
    
        # Make an np output array of strings of length
        # TODO: Fix this, the max len of any sequence may NOT be the length of 
        # the pileup (because insertions lengthen the sequence).
        maxseqlen = "S"+str(len(pileup))
        #seqs = np.zeros(len(pileup), dtype=[("seqs", np.str_, nreads)])
        seqs = np.zeros(nreads, dtype=maxseqlen)
    
        # This is a hack to keep track of the reads that haven't been
        # terminated by a $. This is an artifact of the wacky way pileup
        # files represent sequences. Basically, as we go through we'll assign
        # bases to each read, when a read terminates we remove that read
        # from the positional array.
        incomplete_reads = np.array(range(0,nreads))
        
        # Now run through the list of mapped reads at each base and do the dirty work.
        # Iterate base by base converting one base for each seqeunce. Replace missing
        # values with N
        for ind, line in enumerate(pileup):
    
            # Fetch the reference base
            refbase = line[0]
            
            # Fetch the line of bases per sequence for this base
            bases = line[2]
    
            # While we're going through each line we have to keep track of all
            # the reads that terminate during that line (as indicated by the $)
            # After we finish reading each line, then we'll delete all the
            # completed_reads from the incomplete_reads index array.
            completed_reads = np.array(-1)
    
            # Track the number of $ so far in this line, so we can accurately update the
            # completed position index. This is brutal, unforgivable hax.
            n_dollarsigns = 0
    
            # We need to use an iterator here because I need to next()
            # the iterator a couple times for example with "+" insertions
            # to get the bases to insert, or for "-" to delete the bases
            # to remove.
            b = iter(bases)
            
            # Now do the REAL dirty work
            # There's _gotta_ be a better way to do this...
            #
            # NB: This code is making a _strong_ assumption that there will
            # never been double-digit insertions/deletions. If there is a
            # double digit indel it'll probably behave in mysterious ways.
            for pos, base in enumerate(b):
                ret = ""
                # Now do it real
                if base in "Nn":
                    ret = base
                #if base == "n":
                #    ret = "n"
                elif base == ".":
                    # Dot matches reference sequence on the forward strand
                    ret = refbase
                elif base == ",":
                    # Comma matches reference sequence on the reverse strand
                    ret = refbase
                elif base in "ACGTacgt":
                    # Capital letter indicates difference at base on the 
                    # forward strand.
                    # Lower case indicates difference on reverse strand
                    ret = base
                elif base == "-":
                    # Minus indicates a deletion of one or more bases
                    # In pileup format deletions look like this "-1a".
                    # If you see a minus you need to consume the next n+1
                    # elements of the iterator, where n = the # after the -.
                    # Get the number of deletions here
                    ndeletions = int(next(b))
                    # consume the deleted bases
                    for _ in range(ndeletions):
                        _ = next(b)
                    ret = ""
                elif base == "*":
                    # * indicates a deleted base so return nothing
                    ret = ""
                elif base == "+":
                    # Plus indicates an insertion wrt the reference
                    # Insertions look like this "+2ga", the number of insertions
                    # followed by the bases.
                    ninsertions = int(next(b))
                    # gather the inserted bases
                    ret = ""
                    for _ in range(ninsertions):
                        ret += next(b).upper()
                elif base == "$":
                    # Return the N's to pad the remaining sequence
                    # TODO: Handle padding of sequences in a better way.
                    # It could be done here, but it seems messy.
                    #ret = "N" * (len(seqs[pos]))
                    #print(len(seqs[pos]))
                    #print(ret)
                    # Update the count of dollar signs per line
                    n_dollarsigns += 1
                    completed_reads = np.append(completed_reads, pos-n_dollarsigns)
                else:
                    # TODO Handle '<' & '>'
                    ret = "wat"
                if pos < len(incomplete_reads):
                    seqs[incomplete_reads[pos]] = seqs[incomplete_reads[pos]]+ret
            incomplete_reads = np.delete(incomplete_reads, completed_reads)
        ## Done processing one line of the pileup
    return seqs



def write_aligned_seqs_to_file(data, sample, aligned_seqs, read_labels):
    """ Because vsearch doesn't handle named pipes, or piping at all
    we need to write the aligned sequences per read out to a fasta
    formatted file. 
    """
    with tempfile.NamedTemporaryFile('w', delete=False,
                                          dir=data.dirs.refmapping,
                                          prefix=sample.name+"_",
                                          suffix='.fa') as out:
        for i, line in enumerate(aligned_seqs):
            out.write(">"+read_labels[i]+"\n")
            out.write(line+"\n")
    return out.name



def derep_and_sort(data, sample, aligned_fasta_file):

    #TODO: Delete this. I'm setting filter_min_trim_len
    # to 10 for testing, should delete this prior to shipping
    #data.set_params(18, 10)

    # This is how it's done in cluster_within
    if "gbs" in data.paramsdict["datatype"]:                                                                                                         
        reverse = " -strand both "
    else:
        reverse = " "

    #outfile = os.path.join(data.dirs.refmapping, 
    #                       aligned_fasta_file.split("/")[-1]+".map_derep.fa")
    outfile = tempfile.NamedTemporaryFile(dir=data.dirs.refmapping, 
                                          prefix="derep_and_sort", 
                                          suffix=".fa", delete=False)

    # Stacks are never going to be too big, so faking this to set
    # nthreads to 1, maybe fix this to use the real value if it'll 
    # help performance
    nthreads = 1
    ## do dereplication with vsearch
    cmd = ipyrad.bins.vsearch+\
          " -derep_fulllength "+ aligned_fasta_file+\
          reverse+\
          " -output "+outfile.name+\
          " -sizeout "+\
          " -threads "+str(nthreads)+\
          " -fasta_width 0"+\
          " --minseqlength "+ str(data.paramsdict["filter_min_trim_len"])

    ## run vsearch
    LOGGER.debug("%s", cmd)
    try:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(cmd)
        LOGGER.error(inst)
        sys.exit("Error in vsearch: \n{}\n{}\n{}."\
                 .format(inst, subprocess.STDOUT, cmd))

    return outfile.name



def append_clusters(data, sample, derep_fasta_files):
    """ Append derep'd mapped fasta stacks to the clust.gz file.
    This goes back into the pipeline _before_ the call to muscle
    for alignment.
    """
    LOGGER.debug("Entering append_clusters.")

    ## get clustfile
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clust.gz")

    ## Set the write mode for opening clusters file. If we are
    ## doing "reference" assembly this means we're throwing
    ## out the denovo reads and only keeping the reference
    ## mapped sequences, so we set the flag to 'wb' and overwrite
    ## otherwise use 'ab' to append.
    if data.paramsdict["assembly_method"] == "reference":
        LOGGER.debug("Doing \"reference\" assembly")
        write_flag = 'wb'
    else:
        write_flag = 'ab'


    ## A little bit of monkey business here to get the expected
    ## format right. Downstream expects name lines to end with
    ## * if it's the most abundant read, and + if it's anything else
    ##
    ## TODO: refmapping is not currently checking for a max *
    ## of snps per locus. Probably should fix that.
    with gzip.open(sample.files.clusters, write_flag) as out:
        for fname in derep_fasta_files:
            LOGGER.info("derep_fasta_files %s", fname)
            # We need to update and accumulate all the seqs before
            # we write out to the file or the ipp threads will step
            # all over each other
            seqs = []
            with open(fname) as infile:
                for i, duo in enumerate(itertools.izip(*[iter(infile)]*2)):
                    if i == 0:
                        name = duo[0].strip()+"*"
                    else:
                        name = duo[0].strip()+"+"
                    seqs.append(name+"\n"+duo[1])
#                    LOGGER.info("".join(seqs))
            out.write("".join(seqs))
            out.write("//\n//\n")



def refmap_init(data, sample):
    """Set the mapped and unmapped reads files for this sample
    """
    LOGGER.debug("Entering refmap_init - {}".format(sample.name))
    sample.files.unmapped_reads = os.path.join(\
                             data.dirs.refmapping, sample.name+"-unmapped.bam")
    sample.files.mapped_reads = os.path.join(\
                             data.dirs.refmapping, sample.name+"-mapped.bam")

    ## Save original edits paths to a new files reference in sample
    sample.files["edits_refmap_clean"] = [copy.deepcopy(sample.files.edits[0])]

    ## Build a new tuple for refmap unmapped reads .edits files
    refmap_unmapped_edits = [sample.files.edits[0][0]+".refmap"]
    if "pair" in data.paramsdict["datatype"]:
        refmap_unmapped_edits.append(sample.files.edits[0][1]+".refmap")

    ## Set the edits to point to the unmapped reads
    ## This gets cleaned up at the end of cluster_within
    sample.files.edits = [tuple(refmap_unmapped_edits)]
    return sample



def refmap_stats(data, sample):
    """ Get the number of mapped and unmapped reads for a sample
    and update sample.stats """
    cmd = ipyrad.bins.samtools+\
    " flagstat "+sample.files.unmapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_unmapped_reads"] = int(result.split()[0])

    cmd = ipyrad.bins.samtools+\
    " flagstat "+sample.files.mapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_mapped_reads"] = int(result.split()[0])


def refmap_cleanup(data, sample):
    """ Clean up any loose ends here. Nasty files laying around, etc.
    Also, importantly, recover the files.edits files we stepped on earlier
    when we dropped the unmapped reads back on top of edits and hid
    the originals.
    """
    LOGGER.debug("Entering refmap_cleanup - {}".format(sample.name))
    ## If edits and edits_preview_bak are the same then something borked
    ## so don't delete any files
    if sample.files.edits == sample.files.edits_refmap_clean:
        sample.files.pop("edits_refmap_clean", None)
        return

    ## Remove the unmapped fastq files
    for f in sample.files.edits[0]:
        if(os.path.exists(f)):
            os.remove( f )

    ## Restore original paths to full fastq files
    sample.files.edits = sample.files.edits_refmap_clean
    ## Remove the tmp file reference. The second arg defines what to return
    ## if the key doesn't exist.
    sample.files.pop("edits_refmap_clean", None)


if __name__ == "__main__":
    from ipyrad.core.assembly import Assembly

    shutil.copy("/tmp/wat", "/tmp/watt")
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.set_params(28, '/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa')
    DATA.get_params()
    print(DATA.log)
    #DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
