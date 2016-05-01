#!/usr/bin/env python2.7

""" 
Aligns reads to a reference genome. Returns unaligned reads to the de novo
pipeline. Aligned reads get concatenated to the *clustS.gz files at the 
end of step3
"""

from __future__ import print_function

import os
import gzip
import shutil
import subprocess
import ipyrad
from util import *
from ipyrad.assemble.rawedit import comp

import logging
LOGGER = logging.getLogger(__name__)

# pylint: disable=W0142
# pylint: disable=W0212



def index_reference_sequence(data, force=False):
    """ 
    Attempt to index the reference sequence. This is a little naive
    in that it'll actually _try_ do to the reference every time, but it's
    quick about giving up if it detects the indices already exist. You could
    also test for existence of both index files, but i'm choosing to just let
    smalt do that for us ;) 
    """

    refseq_file = data.paramsdict['reference_sequence']

    # These are smalt specific index files. We don't ever reference
    # them directly except here to make sure they exist, so we don't need
    # to keep them around.
    index_sma = refseq_file+".sma"
    index_smi = refseq_file+".smi"
    # samtools specific index
    index_fai = refseq_file+".fai"

    if all([os.path.isfile(i) for i in [index_sma, index_smi, index_fai]]):
        if force:
            print("    Force reindexing of reference sequence")
        else:
            print("    Reference sequence index exists")
            return

    msg = """\
    *************************************************************
    Indexing reference sequence. 
    This only needs to be done once, and takes just a few minutes
    ************************************************************* """
    if data._headers:
        print(msg)

    ## smalt index for mapping
    cmd = ipyrad.bins.smalt\
        + " index "\
        + " -k "+ str(data._hackersonly["smalt_index_wordlen"])\
        + " " + refseq_file + " " + refseq_file
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

    ## simple samtools index for grabbing ref seqs
    cmd = ipyrad.bins.samtools\
        + " faidx "\
        + " " + refseq_file
    subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)        

    if data._headers:
        print("    Done indexing reference sequence")



def mapreads(args):
    """ 
    Attempt to map reads to reference sequence. This reads in the 
    samples.files.edits .fasta files, and attempts to map each read to the 
    reference sequences. Unmapped reads are dropped right back in the de novo 
    pipeline. Reads that map successfully are processed and pushed downstream 
    and joined with the rest of the data post musle_align. 

    Mapped reads end up in a sam file.
    Unmapped reads are pulled out and put in place of the edits (.fasta) file. 

    The raw edits .fasta file is moved to .<sample>.fasta to hide it in case 
    the mapping screws up and we need to roll-back. Mapped reads stay in the 
    sam file and are pulled out of the pileup later.
    """


    ## get args
    data, sample, noreverse, nthreads = args
    LOGGER.info("Entering mapreads(): %s %s %s", \
                                    sample.name, noreverse, nthreads)

    ## Test that edits filenames actually exist, meaning the edits files exist
    if sample.files.refmap_edits == []:
        raise IPyradWarningExit("""\
    Sample edits files empty. Rerun step2 with force argument.
    Sample files - {}""".format(sample.files))

    ## Set the edited fastq to align for this individual, we set this here
    ## so we can overwrite it with a truncated version if we're in preview mode.
    ## Set the default sample_fastq file to align as the entire original

    ## refmap_init copies the file paths for the full fastq files to
    ## files.edits_refmap_clean. It makes new temp files at sample.files.edits
    ## because this is where we'll dump the unmapped reads for downstream 
    ## denovo files.edits gets put back to point to the original fastq files 
    ## at the end of apply_jobs in cluster_within.
    #sample_fastq = [sample.files.edits[0][0]]

    ## This is the (perhaps merged) input reads file from upstream with the
    ## derep'd reads ready for mapping
    input_reads = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")

    ## Files we'll use during reference sequence mapping
    ## * samhandle - Raw output from smalt reference mapping
    ## * unmapped_bamhandle - bam file containing only unmapped reads
    ## * mapped_bamhandle - bam file with only successfully mapped reads
    ##    In the case of paired end, we only consider pairs where both
    ##    reads map successfully
    ## * unmapped_fastq_handle - Samtools writes out unmapped R1/R2
    ##    reads to the 'edits' directory as .fq. Downstream expects
    ##    these reads to be merged.
    ## * unmapped_fastq_merged - the final merged unmapped file that lives here:
    ##    os.path.join(data.dirs.edits, sample.name+"_refmap_derep.fastq")
    sam = os.path.join(data.dirs.refmapping, sample.name+".sam")
    unmapped_bam = sample.files.unmapped_reads = \
                os.path.join(data.dirs.refmapping, sample.name+"-unmapped.bam")
    mapped_bam = sample.files.mapped_reads = \
                os.path.join(data.dirs.refmapping, sample.name+"-mapped.bam")
    sorted_mapped_bam = sample.files["mapped_reads"]

    ## These are the file(s) that samtools bam2fq will write out
    ## which are then subsequently merged. These files won't be used
    ## after the merging
    unmapped_fastq_r1 = sample.files.refmap_edits[0][0]
    ## If PE set the output file path for R2
    if 'pair' in data.paramsdict["datatype"]:
        unmapped_fastq_r2 = sample.files.refmap_edits[0][1]

    ## Input files initially passed to smalt for read mapping. If doing
    ## SE then you're good to go, just use the derep fastq from upstream.
    ## If doing PE you have to read in the derep file and split it back
    ## into two file R1, and R1. Since we need a file to actually
    ## write to for these, and since we don't need to keep them around
    ## we'll use the refmap_edits files to tem hold the split reads
    ## pre-mapping. This will obviously get written over when samtools
    ## dumps the unmapped reads, but it doesn't matter.
    if 'pair' in data.paramsdict["datatype"]:
        ## TODO: func to split?
        sample_fastq = [sample.files.refmap_edits[0][0],\
                        sample.files.refmap_edits[0][1]]
        split_merged_reads(sample_fastq, input_reads)
    else:
        sample_fastq = [input_reads]

    ## Final output files is merged derep'd refmap'd reads
    ## This will be the only file for unmapped reads used downstream
    unmapped_merged_handle = os.path.join(data.dirs.edits, sample.name+"-refmap_derep.fastq")

    ## get smalt call string
    ##  -f sam       : Output as sam format, tried :clip: to hard mask output 
    ##                 but it shreds the unmapped reads (outputs empty fq)
    ##  -l [pe,mp,pp]: If paired end select the orientation of each read
    ##  -n #         : Number of threads to use
    ##  -x           : Perform a more exhaustive search
    ##  -y #         : proportion matched to reference (sequence similarity)
    ##  -o           : output file
    ##               : Reference sequence
    ##               : Input file(s), in a list. One for R1 and one for R2

    ##  -c #         : proportion of the query read length that must be covered


    ## We are assuming Illumina paired end reads for the orientation 
    ## of mate pairs (orientation: ---> <----).
    if 'pair' in data.paramsdict["datatype"]:
        pairtype = " -l pe "
    else:
        pairtype = " "

    ## wrap in a try statement
    try:
        cmd = ipyrad.bins.smalt+\
            " map -f sam -n " + str(nthreads) +\
            pairtype+\
            " -x -y " + str(data.paramsdict['clust_threshold'])+\
            " -o " + sam +\
            " " + data.paramsdict['reference_sequence'] +\
            " " + " ".join(sample_fastq)
        LOGGER.debug(cmd)
        subprocess.check_call(cmd, shell=True,
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    
        ## Get the mapped and unmapped reads from the sam. For PE both reads 
        ## must map successfully in order to qualify.
        ##   1) Get mapped reads and convert to bam
        ##   2) Sort them and save the path to the bam to sample.files.mapped_
        ## reads. The mapped reads are synced back into the pipeline 
        ## downstream, after muscle aligns the umapped reads.
        ##
        ## samtools view arguments
        ##   -b = write to .bam
        ##   -F = Select all reads that DON'T have this flag. 
        ##         0x4 (segment unmapped)
        ##   -U = Write out all reads that don't pass the -F filter 
        ##        (all unmapped reads go to this file).
        sam_filter_flag = " -F 0x4 "
    
        if 'pair' in data.paramsdict["datatype"]:
            ## Additionally for PE only output read pairs that both align
            sam_filter_flag += " -f 0x2 "
    
        cmd = ipyrad.bins.samtools+\
            " view -b"+\
                sam_filter_flag+\
                " -U " + unmapped_bam+\
                " " + sam+\
                " > " + mapped_bam
        LOGGER.info("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ## Step 2 sort mapped bam
        ##   -T = Temporary file name, this is required by samtools, ignore it
        ##        Here we hack it to be samhandle.tmp cuz samtools cleans it up
        ##   -O = Output file format, in this case bam
        ##   -o = Output file name
        cmd = ipyrad.bins.samtools+\
            " sort -T "+sam+".tmp" +\
            " -O bam "+mapped_bam+\
            " -o "+sorted_mapped_bam
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ## Step 3 index mapped reads
        ## Samtools pileup needs the bam to be indexed
        ## No arguments, a very simple function. It writes the index to 
        ## a default location
        cmd = ipyrad.bins.samtools+\
            " index " + mapped_bam
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    
        ##############################################
        ## Do unmapped
        ## Output the unmapped reads to the original
        ## sample.edits fq path.
        ##############################################
        outfiles = [unmapped_fastq_r1]
        if 'pair' in data.paramsdict["datatype"]:
            outfiles.append(unmapped_fastq_r2)
            outflags = " -1 " + outfiles[0]+\
                       " -2 " + outfiles[1]
        else:
            outflags = " -0 " + outfiles[0]

        cmd = ipyrad.bins.samtools+\
            " bam2fq " + outflags+\
                " " + unmapped_bam

        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)

        ## Finally, merge the unmapped reads, which is what cluster()
        ## expects. If SE, just rename the outfile. In the end
        ## <sample>-refmap_derep.fq will be the final output
        if 'pair' in data.paramsdict["datatype"]:
            LOGGER.info("Merging unmapped reads {} {}".format(outfiles[0],
                                                              outfiles[1]))
            ## merge_pairs wants the files to merge in this stupid format,
            ## also the '1' at the end means "really merge" don't just join w/ nnnn
            merge_pairs(data, [(outfiles[0], outfiles[1])], unmapped_merged_handle, 1)
        else:
            LOGGER.debug("Renaming unmapped reads file from {} to {}"\
                        .format(outfiles[0], unmapped_merged_handle))
            os.rename(outfiles[0], unmapped_merged_handle)


    except subprocess.CalledProcessError as inst:
        ## Handle error outside try statement
        LOGGER.error(inst)
        raise 



def ref_muscle_chunker(args):
    """ 
    Run bedtools to get all overlapping regions. Then derep the list of reads
    split this list of 
    regions and distribute it across a bunch of parallel threads.

    1) Run bedtools merge to get a list of all contiguous blocks of bases
    in the reference seqeunce where one or more of our reads overlap.
    The output will look like this:
            1       45230754        45230783
            1       74956568        74956596
            ...
            1       116202035       116202060

    """
    ## parse args
    data, sample = args
    LOGGER.info('entering ref_muscle_chunker')

    ## Regions is a giant list of 5-tuples, of which we're only really 
    ## interested in the first three: chrom, start, and end position.
    regions = bedtools_merge(data, sample)

    if len(regions) > 0:
        ## Empty array to hold our chunks
        get_overlapping_reads([data, sample, regions])
        #for region in regions.strip().split("\n"):
        #    LOGGER.info("region %s", region)

    else:
        msg = "No reads mapped to reference sequence - {}".format(sample.name)
        LOGGER.warn(msg)


def get_overlapping_reads(args):
    """
    Pull aligned reads out of sorted mapped bam files and append them to the 
    clust.gz file so the fall into downstream (muscle alignment) analysis. 
    Here's the exact order of ops:
    
    1) Coming into this function we have sample.files.mapped_reads 
        as a sorted bam file, and a passed in list of regions to evaluate.
    2) Get all reads overlapping with each individual region.
    3) Write the aligned sequences out to a fasta file (vsearch doesn't 
        handle piping data) (IT WILL VERY SOON, actually).
    4) Call vsearch to derep reads
    5) Append to the clust.gz file.
    """

    data, sample, regions = args
    locus_list = []
    reads_merged = 0

    ## Set the write mode for opening clusters file.
    ## 1) if "reference" then only keep refmapped, so use 'wb' to overwrite 
    ## 2) if "denovo+reference" then 'ab' adds to end of denovo clust file
    if data.paramsdict["assembly_method"] == "denovo+reference":
        write_flag = 'ab'
    else:
        write_flag = 'wb'
    sample.files.clusters = os.path.join(
                                data.dirs.clusts, sample.name+".clust.gz")
    outfile = gzip.open(sample.files.clusters, write_flag)

    ## If we are appending reference reads to the denovo then you need to 
    ## prepend the first stack with a //\n// separator
    if data.paramsdict["assembly_method"] == "denovo+reference":
        outfile.write("\n//\n//\n")

    # Wrap this in a try so we can clean up if it fails.
    try:
        ## For each identified region, build the pileup and write out the fasta
        for line in regions.strip().split("\n"):
            #LOGGER.info("line %s ", line)

            # Blank lines returned from bedtools make the rest of the 
            # pipeline angry. Filter them.
            if line == "":
                continue

            chrom, region_start, region_end = line.strip().split()[0:3]


            ## bam_region_to_fasta returns a chunk of fasta sequence
            args = [data, sample, chrom, region_start, region_end]
            clust = bam_region_to_fasta(*args)

            ## If bam_region_to_fasta fails for some reason it'll return [], 
            ## in which case skip the rest of this. Normally happens if reads
            ## map successfully, but too far apart.
            if not clust:
                continue

            ## Derep_fasta_files are merged for PE
            #derep_fasta_files.append(derep_fasta)
            locus_list.append(clust)

            ## write chunk of 1000 loci and clear list to minimize memory
            if not len(locus_list) % 1000:
                outfile.write("\n//\n//\n".join(locus_list))
                outfile.write("\n//\n//\n")
                locus_list = []
        ## write remaining
        outfile.write("\n//\n//\n".join(locus_list))
        outfile.write("\n//\n//\n")

        ## close handle
        outfile.close()


    except Exception as inst:
        LOGGER.error("Exception inside get_aligned_reads - {}".format(inst))
        raise

    finally:
        LOGGER.info("Total merged reads for {} - {}"\
                     .format(sample.name, reads_merged))
        sample.stats.reads_merged = reads_merged

        # Clean up all the tmp files
        # Be a little careful. Don't remove files if they're already gone :-/
        # for i in derep_fasta_files:
        #     if os.path.isfile(i):
        #         os.remove(i)
        # for j in aligned_seq_files:
        #     if os.path.isfile(j[0]):
        #         os.remove(j[0])
        #     if 'pair' in data.paramsdict['datatype']:
        #         if os.path.isfile(j[1]):
        #             os.remove(j[1])


def split_merged_reads(sample_fastq, input_reads):
    """
    Reads in the merged derep file from upstream and splits
    back into two files for R1 and R2. Arguments are:
        - sample_fastq: a list of the two file paths to write out
        - input_reads: the path to the input merged reads
    """
    R1_file = sample_fastq[0]
    R2_file = sample_fastq[1]
    R1 = open(R1_file, 'w')
    R2 = open(R2_file, 'w')

    with open(input_reads, 'r') as infile: 
        ## Read in the infile two lines at a time, seq name and sequence
        duo = itertools.izip(*[iter(infile)]*2)
        while 1:
            try:
                itera = duo.next()
            except StopIteration:
                break
            ## R1 needs a newline, but R2 inherits it from the original file
            R1.write(itera[0]+itera[1].split("nnnn")[0]+"\n")
            R2.write(itera[0]+itera[1].split("nnnn")[1])
    R1.close()
    R2.close()


def bedtools_merge(data, sample):
    """ 
    Get all contiguous genomic regions with one or more overlapping
    reads. This is the shell command we'll eventually run

        bedtools bamtobed -i 1A_0.sorted.bam | bedtools merge [-d 100]
            -i <input_bam>  :   specifies the input file to bed'ize
            -d <int>        :   For PE set max distance between reads
    """
    LOGGER.info("Entering bedtools_merge: %s", sample.name)

    if 'pair' in data.paramsdict["datatype"]:
        ## check mean insert size for this sample
        ## and update hackersonly.max_inner_mate_distance
        ## if need be. This value controls how far apart 
        ## mate pairs can be to still be considered for 
        ## bedtools merging downstream
        cmd = ipyrad.bins.samtools+\
            " stats " + sample.files.mapped_reads + " | grep SN"
        LOGGER.debug("%s", cmd)
        ret = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

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

        LOGGER.debug("avg {} stdv {} avg_len {}"\
                     .format(avg_insert, stdv_insert, avg_len))

        ## If all values return successfully set the max inner mate distance to 
        if all([avg_insert, stdv_insert, avg_len]):
            hack = avg_insert + (3 * stdv_insert) - (2 * avg_len)
            data._hackersonly["max_inner_mate_distance"] = hack
        else:
            ## If something fsck then set a relatively conservative distance
            data._hackersonly["max_inner_mate_distance"] = 200
        LOGGER.debug("inner mate distance for {} - {}".format(sample.name,\
                    data._hackersonly["max_inner_mate_distance"]))

    ## Set the -d flag to tell bedtools how far apart to allow mate pairs.
    if 'pair' in data.paramsdict["datatype"]:
        bed_dflag = " -d " + str(data._hackersonly["max_inner_mate_distance"])
    else:
        bed_dflag = " "

    cmd = ipyrad.bins.bedtools+\
        " bamtobed "+\
        " -i " + sample.files.mapped_reads+\
        " | " +\
        ipyrad.bins.bedtools +\
        " merge "+\
        bed_dflag
    LOGGER.debug("%s", cmd)
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    ## Report the number of regions we're returning
    LOGGER.info("bedtools_merge: Got # regions: %s", 
                 str(len(result.strip().split("\n"))))
    return result



def bam_region_to_fasta(data, sample, chrom, region_start, region_end):
    """ 
    Take the chromosome position, and start and end bases and 
    return sequences of all reads that overlap these sites. 
    This is the command we're building:

    samtools view -b 1A_sorted.bam 1:116202035-116202060 | \
             samtools bam2fq <options> -

            -b      : output bam format
            -0      : For SE, output all reads to this file
            -1/-2   : For PE, output first and second reads to different files
            -       : Tell samtools to read in from the pipe

    Write out the sam output and parse it to return as fasta for clust file. 
    We also grab the reference sequence with a @REF header to aid in alignment.
    This is removed post-alignment. 
    """

    LOGGER.info("Grabbing bam_region_to_fasta: %s %s %s %s", 
                sample.name, chrom, region_start, region_end)

    try:
        ## Make the samtools view command output bam in the region of interest
        bamf = os.path.join(data.dirs.refmapping, sample.name+"-mapped.bam")

        view_ref = ipyrad.bins.samtools+\
                " faidx "+\
                data.paramsdict["reference_sequence"]+\
                " {}:{}-{}".format(chrom, str(int(region_start)+1), region_end)

        ## Reach out and grap the reference sequence that overlaps this region
        ## we'll paste this in at the top of each stack to aid alignment
        ref = subprocess.check_output(view_ref, shell=True)
        ## parse sam to fasta. Save ref location to name
        name, seq = ref.strip().split("\n", 1)
        seq = "".join(seq.split("\n"))
        fasta = [name[1:]+"_REF;+"+"\n"+seq]

        ## if PE then you have to merge the reads here
        if "pair" in data.paramsdict["datatype"]:
            ## TODO: Can we include the reference sequence in the PE clust.gz?
            ## if it's longer than the merged pairs it makes hella indels
            ## Drop the reference sequence for now...
            fasta = []

            ## Set output files and flags for PE/SE
            ## Create temporary files for R1, R2 and merged
            prefix = os.path.join(data.dirs.refmapping, sample.name \
                            + chrom + region_start + region_end)
            R1 = prefix+"-R1"
            R2 = prefix+"-R2"
            merged = prefix+"-merged"
            outflags = " -1 " + R1 +\
                        " -2 " + R2
            pe_view_cmd = ipyrad.bins.samtools + \
                    " view -b "+\
                    bamf + \
                    " {}:{}-{}".format(chrom, region_start, region_end)
            bam2fq_cmd = ipyrad.bins.samtools+" bam2fq " + outflags + " - "
            cmd = " | ".join((pe_view_cmd, bam2fq_cmd))

            LOGGER.info(cmd)
            try:
                subprocess.check_output(cmd, shell=True)
                nmerged = merge_pairs(data, [(R1, R2)], merged, 1)

                infile = open(merged)
                quatro = itertools.izip(*[iter(infile)]*4)
                while 1:
                    try:
                        itera = quatro.next()
                    except StopIteration:
                        break
                    ## TODO: figure out a real way to get orientation for PE
                    orient = "+"
                    label = itera[0].split(";")[0] + ";" + "{}:{}-{}".format(\
                            chrom, str(int(region_start)+1), region_end) \
                            + ";" + itera[0].split(";")[1] + ";" + orient
                    fasta.append(label+"\n"+itera[1].strip())
            except Exception as inst:
                LOGGER.error("Exception in bam_region_to_fasta doing PE - {}".format(inst))
                raise
            ## always clean up the temp files
            finally:
                files = []
                #files = [R1, R2, merged]
                for f in files:
                    if os.path.exists(f):
                        os.remove(f)
        else:
            ## Kind of dumb to do PE and SE so different, but this
            ## is probably faster for SE since you don't have to
            ## write to a bunch of files
            se_view_cmd = ipyrad.bins.samtools+\
                        " view "+\
                        bamf+\
                        " {}:{}-{}".format(chrom, region_start, region_end)
            sam = subprocess.check_output(se_view_cmd, shell=True)

            ## do not join seqs that 
            for idx, line in enumerate(sam.strip().split("\n")):
                bits = line.split("\t")

                ## Read in the 2nd field (FLAGS), convert to binary
                ## and test if the 7th bit is set which indicates revcomp
                orient = "+"
                if int('{0:012b}'.format(int(bits[1]))[7]):
                    orient = "-"
                ## Rip insert the mapping position between the seq label and
                ## the vsearch derep size
                label = bits[0].split(";")[0] + ";" + "{}:{}-{}".format(chrom, \
                        str(int(region_start)+1), region_end) + ";" + \
                        bits[0].split(";")[1] + ";"

                fasta.append(label+orient+"\n"+bits[9])

        return "\n".join(fasta)

    except Exception as inst:
        LOGGER.error("Caught exception in bam_region_to_fasta - {}".format(inst))
        raise


def refmap_init(args):
    """
    Set the mapped and unmapped reads files for this sample
    """
    #LOGGER.info("""\
    #    refmap_init | sample %s | file %s \
    #    """, sample.name, sample.files.edits)
    data, sample = args
    sample.files.unmapped_reads = os.path.join(\
                             data.dirs.refmapping, sample.name+"-unmapped.bam")
    sample.files.mapped_reads = os.path.join(\
                             data.dirs.refmapping, sample.name+"-mapped.bam")

    ## Build a new tuple for refmap unmapped reads .edits files
    refmap_unmapped_edits = [sample.files.edits[0][0]+".refmap"]
    if "pair" in data.paramsdict["datatype"]:
        refmap_unmapped_edits.append(sample.files.edits[0][1]+".refmap")
    sample.files.refmap_edits = [tuple(refmap_unmapped_edits)]
    #LOGGER.info("after %s", sample.files.refmap_edits)

    return sample


def refmap_stats(data, sample):
    """ 
    Get the number of mapped and unmapped reads for a sample
    and update sample.stats 
    """
    cmd = ipyrad.bins.samtools+" flagstat "+sample.files.unmapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_unmapped_reads"] = int(result.split()[0])

    cmd = ipyrad.bins.samtools+" flagstat "+sample.files.mapped_reads
    result = subprocess.check_output(cmd, shell=True,
                                          stderr=subprocess.STDOUT)
    sample.stats["refseq_mapped_reads"] = int(result.split()[0])



def refmap_cleanup(data, sample):
    """ 
    Clean up any loose ends here. Nasty files laying around, etc.
    Also, importantly, recover the files.edits files we stepped on earlier
    when we dropped the unmapped reads back on top of edits and hid
    the originals.
    """
    LOGGER.info("Entering refmap_cleanup - {}".format(sample.name))
    ## If edits and edits_preview_bak are the same then something borked
    ## so don't delete any files
    if sample.files.edits == sample.files.edits_refmap_clean:
        sample.files.pop("edits_refmap_clean", None)
        return

    ## Remove the unmapped fastq files
    for sfile in sample.files.edits[0]:
        if os.path.exists(sfile):
            LOGGER.info("removing %s", sample.files.edits[0])
            os.remove(sfile)

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
    #DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)
