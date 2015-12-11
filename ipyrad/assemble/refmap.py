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
from shutil import copy2
import tempfile
import itertools
import subprocess
import numpy as np
from ipyrad.assemble.rawedit import comp

import logging
LOGGER = logging.getLogger(__name__)

def mapreads(args):
    """ Attempt to map reads to reference sequence. This reads in the 
    samples.files.edits .fasta files, and attempts to map each read to the 
    reference sequences. Unmapped reads are dropped right back in the de novo 
    pipeline. Reads that map successfully are processed and pushed downstream 
    and joined with the rest of the data post musle_align. The read mapping 
    produces a sam file, unmapped reads are pulled out of this and dropped back 
    in place of the edits (.fasta) file. The raw edits .fasta file is moved to 
    .<sample>.fasta to hide it in case the mapping screws up and we need to roll-back.
    Mapped reads stay in the sam file and are pulled out of the pileup later."""

    ## get args
    data, sample, preview, noreverse, nthreads = args
    LOGGER.debug("Entering mapreads(): %s %s %s %s", sample.name, preview, noreverse, nthreads )

    ## preview
    if preview:
        print("preview: in run_full, using", nthreads)

    ## Files we'll use during reference sequence mapping
    ##
    ## samhandle - Raw output from smalt reference mapping
    ## unmapped_bamhandle - bam file containing only unmapped reads
    ## mapped_bamhandle - bam file with only successfully mapped reads
    ##    In the case of paired end, we only consider pairs where both
    ##    reads map successfully
    ## sorted_*_bamhandle - Sorted bam files for mapped and unmapped
    ##    This is a precursor to bam2fq which requires sorted bams
    ## unmapped_fastq_handle - The final output file of this function
    ##    writes out unmapped reads to the 'edits' directory as .fq
    ##    which is what downstream analysis expects
    samhandle = os.path.join(data.dirs.refmapping, sample.name+".sam")
    unmapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-unmapped.bam")
    mapped_bamhandle = os.path.join(data.dirs.refmapping, sample.name+"-mapped.bam")
    sorted_unmapped_bamhandle = sample.files["unmapped_reads"]
    sorted_mapped_bamhandle = sample.files["mapped_reads"]

    ## TODO, figure out why edits is a tuple? There could be multiple edits files, yes?
    ## but by the time we get to step three they are all collapsed in to one big file
    ## which is the first element of the tuple, yes?

    ##
    ## TODO: This is hackish, we are overwriting the fastq that contains all reads
    ## might want to preserve the full .fq file in case of a later 'force' command
    unmapped_fastq_handle = fastq_file = sample.files.edits[0][0]

    ## Check for a hidden .fastq file, this is the original fq, pre-ref mapping.
    ## if it exists, recover it. 
    fastq_dotfile = data.dirs.edits + "/." + fastq_file.split("/")[-1]

    # TODO: Why is this fucking broken? If I import copy2 at the top of the file
    # it tells me "global name shutil not found". UGNN!
    from shutil import copy2
    if os.path.isfile( fastq_dotfile ): 
        copy2( fastq_dotfile, fastq_file )

    # Preserve the original full .fq file as a hidden "dot" file.
    copy2( unmapped_fastq_handle, fastq_dotfile )

################
# TODO: Paired end isn't handled yet, but it needs to be.
    ## datatype variables
    if data.paramsdict["datatype"] in ['gbs']:
        reverse = " -strand both "
        cov = " -query_cov .35 "
    elif data.paramsdict["datatype"] in ['pairgbs', 'merged']:
        reverse = "  -strand both "
        cov = " -query_cov .60 "
    else:  ## rad, ddrad, ddradmerge
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override reverse clustering option
    if noreverse:
            reverse = " -leftjust "
            print(noreverse, "not performing reverse complement clustering")
################


    ## get call string
    cmd = data.smalt+\
        " map -f sam -n " + str(nthreads) +\
        " -o " + samhandle +\
        " " + data.paramsdict['reference_sequence'] +\
        " " + sample.files.edits[0][0]

    LOGGER.debug( "%s", cmd )
    ## run smalt
    if preview:
        ## make this some kind of wait command that kills after a few mins
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    else:
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)

    ## Get the reads that map successfully. For PE both reads must map
    ## successfully in order to qualify.
    ##   1) Get mapped reads and convert to bam
    ##   2) Sort them and save the path to the bam to sample.files.mapped_reads
    ## The mapped reads are synced back into the pipeline downstream, after
    ## muscle aligns the umapped reads.
    ##
    ## samtools view arguments
    ##   -b = write to .bam
    ##   -F = Select all reads that DON'T have this flag (0x4 means unmapped)
    ##   -U = Write out all reads that don't pass the -F filter (all unmapped reads
    ##        go to this file.
    ##        <TODO>: This is deeply hackish right now, it will need some
    ##                serious thinking to make this work for PE, etc.
    cmd = data.samtools+\
        " view -b -F 0x4 "+\
            " -U " + unmapped_bamhandle + \
            " " + samhandle+\
            " > " + mapped_bamhandle
    LOGGER.debug( "%s", cmd )
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## Step 2 sort mapped bam
    ##   -T = Temporary file name, this is required by samtools, you can ignore it
    ##        Here we just hack it to be samhandle.tmp cuz samtools will clean it up
    ##   -O = Output file format, in this case bam
    ##   -o = Output file name
    cmd = data.samtools+\
        " sort -T "+samhandle+".tmp" +\
        " -O bam "+mapped_bamhandle+\
        " -o "+sorted_mapped_bamhandle
    LOGGER.debug( "%s", cmd )
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ##############################################
    ## Do unmapped
    ##############################################

    ## Fetch up the unmapped reads and write them to the edits .fasta file
    ## In order to do this we have to go through a little process:
    ##   1) Get the unmapped reads and convert to bam
    ##   2) Sort them 
    ##   3) Dump bam to fastq

    ## Step 1 get unmapped reads with samtools view
    ##   -b = write out to .bam
    ##   -f = Select only reads with associated flags set (0x4 means unmapped)
#    cmd = data.samtools+\
#        " view -b -f 0x4 "+samhandle+\
#            " > " + unmapped_bamhandle
#    subprocess.call(cmd, shell=True,
#                         stderr=subprocess.STDOUT,
#                         stdout=subprocess.PIPE)

    ## Step 2 sort unmapped bam
    ##   -T = Temporary file name, this is required by samtools, you can ignore it
    ##        Here we just hack it to be samhandle.tmp cuz samtools will clean it up
    ##   -O = Output file format, in this case bam
    ##   -o = Output file name
    cmd = data.samtools+\
        " sort -T "+samhandle+".tmp" +\
        " -O bam "+unmapped_bamhandle+\
        " -o "+sorted_unmapped_bamhandle
    LOGGER.debug( "%s", cmd )
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## Step 3 Dump sorted bam to fastq
    ## No args for this one. Pipe it through gzip, because the downstream
    ## analysis expects gzipped fq
    cmd = data.samtools+\
        " bam2fq "+sorted_unmapped_bamhandle+\
        " | gzip > "+unmapped_fastq_handle
    LOGGER.debug( "%s", cmd )
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

    ## bs for experimentation. is samtools/smalt fscking the output?
    ## TODO: get rid of this when mapped reads get returned to the 
    ## pipeline right. or whenever, it's not useful.
#    cmd = data.samtools+\
#        " bam2fq "+sorted_mapped_bamhandle+\
#        " | gzip >> "+unmapped_fastq_handle
#    subprocess.call(cmd, shell=True,
#                         stderr=subprocess.STDOUT,
#                         stdout=subprocess.PIPE)

    ## This is the end of processing for each sample. Stats
    ## are appended to the sample for mapped and unmapped reads 
    ## during cluster_within.cleanup()

def getalignedreads( data, sample ):
    """Pull aligned reads out of sorted mapped bam files and
    append them to the clustS.gz file so the fall into downstream analysis """

    mapped_fastq_handle = "/tmp/wat.fq"
    ## Build the samtools bam2fq command to push bam out
    cmd = data.samtools+\
        " bam2fq "+sample.files["mapped_reads"]+\
        " > "+mapped_fastq_handle
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)

def bamtopileup( data, sample):
    LOGGER.debug( "Entering bamtopileup: %s", sample.name )

def pileuptofasta( data, sample):
        # Read in the mpileup file for one locus
        #TESTFILE = "/Volumes/WorkDrive/ipyrad/refhacking/pero.1locus.mpileup"
    #    TESTFILE = "/Volumes/WorkDrive/ipyrad/refhacking/perom/pero.1locus.pileup"
    #TESTFILE = "/Volumes/WorkDrive/ipyrad/refhacking/1A.pileup"
    with open( TESTFILE, 'r' ) as pfile:
        pileup = []
        for line in pfile:
            dat = np.array(line.strip().split())
            # Handle the case where the number of aligned reads is 0, in which case there
            # is no 4th field for the sequence info. Maybe just pass overit.
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
                    # Capital letter indicates difference at base on the forward strand
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
                    for i in range(ndeletions):
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
                    ret=""
                    for i in range(ninsertions):
                        ret+=next(b).upper()
                elif base == "$":
                    # Return the N's to pad the remaining sequence
                    ret = "N" * (len(seqs[pos]))
                    print(len(seqs[pos]))
                    print(ret)
                    # Update the count of dollar signs per line
                    n_dollarsigns +=1
                    completed_reads = np.append( completed_reads, pos-n_dollarsigns)
                else:
                    # TODO Handle '<' & '>'
                    ret = "wat"
                if pos < len(incomplete_reads):
                    seqs[incomplete_reads[pos]] = seqs[incomplete_reads[pos]]+ret
            if( completed_reads.size > 1 ):
                print("completed-", completed_reads)
                print("incomplete-", incomplete_reads)
                print("\n")
            incomplete_reads = np.delete( incomplete_reads, completed_reads )
        print(incomplete_reads)
    print( seqs )  

def refmap_stats( data, sample ):
    """ Get the number of mapped and unmapped reads for a sample
    and update sample.stats """
    cmd = data.samtools+\
    " flagstat "+sample.files.unmapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_unmapped_reads"] = int(result.split()[0])

    cmd = data.samtools+\
    " flagstat "+sample.files.mapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_mapped_reads"] = int(result.split()[0])


if __name__ == "__main__":
    from ipyrad.core.assembly import Assembly
    import shutil
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
    pass
