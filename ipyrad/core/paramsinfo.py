#!/usr/bin/env python2.7

""" Return explanation and options for each parameter. 
    ip.get_params_info(1) or ip.get_params_info("working_directory") 
    return the same result. If not argument, a summary of the available
    parameters and their numbered references is returned. 
"""

from __future__ import print_function
from collections import OrderedDict


def get_params_info(param=""):
    """ Returns detailed information for the numbered parameter. 
        Further information is available in the tutorial."""

    paramsinfo = OrderedDict([
    ("1", """
        (1) working_directory ------------------------------------------------
        Path to the working directory where all data files will be saved. 
        This parameter affects all steps of assembly (1-7). 
        Examples: 
        ----------------------------------------------------------------------
        data.setparams(1) = "/home/user/rad_analysis/"   ## full path example
        data.setparams(1) = "./"                         ## current work. dir.
        data.setparams('working_directory') = "./"       ## verbose
        ----------------------------------------------------------------------
        """),

    ("2", """
        (2) raw_fastq_path ---------------------------------------------------
        The directory or files (selected with * wildcard selector) in which 
        FASTQ data files reside. Files can be gzipped. This parameter affects 
        only step 1 of assembly. Examples:
        ----------------------------------------------------------------------
        data.setparams(2) = "/home/user/rad_analysis/*.fastq  ## full path
        data.setparams(2) = "raw/*.fq"                        ## relative path
        data.setparams(2) = "raw/*.fastq.gz"                  ## gzipped data 
        data.setparams("raw_fastq_path") = "raw/*.fastq.gz"   ## verbose 
        ----------------------------------------------------------------------
        """),

    ("3", """
        (3) barcodes_path ----------------------------------------------------
        Path to the barcodes file used in step 1 of assembly for 
        demultiplexing. If data are already demultiplexed this can be left 
        blank. This parameter affects only step 1 of assembly. Examples:
        ----------------------------------------------------------------------
        data.setparams(3) = "/home/user/rad_analysis/barcodes.txt ## full path
        data.setparams(3) = "./barcodes.txt                   ## relative path
        data.setparams("barcodes_path") = "./barcodes.txt"    ## verbose
        ----------------------------------------------------------------------
        """),

    ("4", """
        (4) sorted_fastq_path ------------------------------------------------
        Path to demultiplexed fastq data. If left blank, this is assigned
        automatically to <data.name>_fastq/ within the working directory. If your
        data are already demultiplexed then you must enter the location of your  
        data here. Wildcard selectors can be used to select a subsample of files 
        within a directory, else all files are selected in the directory.
        This parameter affects only step 2 of assembly. 
        Examples:
        ----------------------------------------------------------------------
        data.setparams(4) = "/home/user/data/*.fastq    ## set data location
        data.setparams(4) = ""                ## defaults to working directory
        data.setparams(4) = "./"                  ## uses current directory
        data.setparams("sorted_fastq_path") = ""  ## Use default
        ----------------------------------------------------------------------
        """),

    ("5", """
        (5) datatype ---------------------------------------------------------
        Options: rad, gbs, ddrad, pairddrad, pairgbs, merged.
        This parameter affects all steps of assembly (1-7).         
        Examples:
        ----------------------------------------------------------------------
        data.setparams(7) = 'rad'                     ## rad data type
        data.setparams(7) = 'gbs'                     ## gbs data type
        data.setparams(7) = 'pairddrad'               ## gbs data type        
        data.setparams(7) = 'merged'                  ## merged data type
        data.setparams("datatype") = 'ddrad'          ## verbose
        ----------------------------------------------------------------------
        """),

    ("6", """
        (6) restriction_overhang ---------------------------------------------
        A tuple containing one or two restriction overhangs. Single digest 
        RADseq with sonication requires only one overhange, all other data 
        types should have two. The first is used for detecting barcodes, the 
        second is not required, but is used in filtering, and is needed for 
        removal from short DNA fragments. Use .preview() methods (see 
        documentation) to ensure that restriction overhangs are entered 
        correctly. This parameter affects steps 1,2,4,5, and 7 of assembly. 
        Examples:
        ----------------------------------------------------------------------
        data.setparams(8) = ("TGCAG", "")           ## default rad (PstI)
        data.setparams(8) = ("CWGC", "CWGC")        ## gbs or pairgbs (ApeKI)
        data.setparams(8) = ("CAGT", "AATT")        ## ddrad (ApeKI, MSI)
        data.setparams(8) = ("CAGT", "AATT")        ## pairddrad (ApeKI, MSI)        
        data.setparams("restriction_overhang") = ("CAGT", "AATT")   ## verbose
        ----------------------------------------------------------------------
        """),

    ("7", """
        (7) mindepths --------------------------------------------------------
        A tuple containing two values, the mindepth for statistical base calls
        based a binomial probability with H and E estimated from the data, and
        the mindepth for majority-rule base calls. Base calls are made at >= 
        the value entered. For most reasonable estimates of E and H, 
        statistical base calls cannot be made below 5 or 6, and will instead 
        be called N. It may often be advantageous to use a low value for 
        majrule calls to preserve most data during assembly within-samples, 
        so that more data is clustered between samples. Low depth data can be 
        filtered out later from the final data set if needed. 
        The parameter affects steps 5 and 7 of assembly. 
        Examples:
        ----------------------------------------------------------------------
        data.setparams(9) = (6, 6)    ## only stat base calls down to depth=6
        data.setparams(9) = (10, 5)   ## stat calls above 9, majrule from 9-5.
        data.setparams(9) = (10, 1)   ## stat calls above 9, majrule from 9-1.
        data.setparams(mindepths) = (6, 1)    ## verbose
        ----------------------------------------------------------------------
        """),

    ("10", """
        (10) clust_threshold -------------------------------------------------
        Clustering threshold. 
        Examples:
        ----------------------------------------------------------------------
        data.setparams(10) = .85          ## clustering similarity threshold
        data.setparams(10) = .90          ## clustering similarity threshold
        data.setparams(10) = .95          ## very high values not recommended 
        data.setparams("clust_threshold") = .83  ## verbose
        ----------------------------------------------------------------------
        """),

    ("11", """
        (11) minsamp ---------------------------------------------------------
        Minimum number of samples a locus must be shared across to be included
        in the exported data set following filtering for sequencing depth, 
        paralogs, ...
        Examples
        ----------------------------------------------------------------------
        data.setparams(11) = 4            ## min 4; most inclusive phylo data 
        data.setparams(11) = 20           ## min 20; less data, less missing
        data.setparams(11) = 1            ## min 1; most data, most missing
        data.setparams("minsamp") = 4     ## verbose
        ----------------------------------------------------------------------
        """),
    ("12", """
        (12) max_shared_heterozygosity ---------------------------------------
        ...
        ----------------------------------------------------------------------
        data.setparams(12) = .25          ## set as proportion of samples
        data.setparams(12) = 4            ## set as number of samples
        data.setparams(12) = 9999         ## set arbitrarily high
        data.setparams("max_shared_heterozygosity) = 4      ## verbose
        ----------------------------------------------------------------------
        """),

    ("13", """
        (13) prefix_outname --------------------------------------------------

        ----------------------------------------------------------------------
        data.setparams(13) = test          ## set a name
        data.setparams(13) = c85d4m8p4     ## set a name of parameters values
        data.setparams("prefix_outname") = c85d4m8p4   ## verbose
        ---------------------------------------------------------------------- 
        """),
    ("27", """
        (27) assembly_method -------------------------------------------------
        A string specifying the desired assembly method. There are three 
        available options for assembly method:
            denovo    -   Denovo assembly is the classic pyrad method, and
                          it is the <default> unless otherwise specified.
                          Denovo will cluster and align all reads from scratch
            reference -   Reference assembly will map and align reads to the
                          provided reference sequence, which must be specified
                          in parameter 28 (reference_sequence). Strict refer-
                          ence assembly will throw out all unmapped reads, 
                          which could be a significant proportion detpending
                          on the distance between your reference and study
                          species. Most times you'll want to use 'hybrid'.
            hybrid    -   Hybrid assembly will attempt to map and align reads
                          to the provided reference sequence which must be set
                          in parameter 28. It will also denovo assemble all
                          unmapped reads, and then merge assembled mapped and 
                          unmapped for downstream analysis. This is what you'll
                          want most of the time if you're passing in a refer-
                          ence sequence.
        ----------------------------------------------------------------------
        data.setparams(27) = denovo        ## set a name
        data.setparams(27) = hybrid        ## set a name of parameters values
        data.setparams("assembly_method") = reference   ## verbose
        ---------------------------------------------------------------------- 
        """),
    ("28", """
        (28) reference_sequence ----------------------------------------------
        The path to the reference sequence you desire to map your reads to.
        The reference may be either fasta or gzipped fasta. It should be a 
        complete reference sequence, including all chromosomes, scaffolds, and
        contigs in one huge file (most reference sequences available will be
        in this format, especially non-model references). The first time you 
        attempt to use this sequence it will be indexed (we are using smalt 
        for reference mapping). This is a time intensive process so expect the 
        first run to take some time, certainly more than ten minutes, but less 
        than an hour. If you desire to index the reference yourself you can do 
        this, but best not to unless you really care about smalt indexing 
        settings. We chose conservative defaults that have worked well for us 
        on other projects. 

        A word on the format of the path (this is important). The path may
        either be a full path (desirable) or a path relative to the directory
        you are running ipyrad from (supported but be careful of the path).
        ----------------------------------------------------------------------
        data.setparams(28) = /home/wat/data/reference.fa  ## set a full path
        data.setparams(28) = ./data/reference.fa.gz       ## set a relative path
        data.setparams("reference_sequence") = ./data/reference.fa   ## verbose
        ---------------------------------------------------------------------- 
        """),

    ])

    if param == "*":
        for key in paramsinfo:
            print(paramsinfo[str(key)])
    elif param:
        try: 
            print(paramsinfo[str(param)])
        except (KeyError, ValueError) as err:
            print("\tKey name/number not recognized", err)
    else:
        print("Enter a name or number for explanation of the parameter\n")
        for key in paramsinfo:
            print(paramsinfo[str(key)].split("\n")[1][2:-10])


if __name__ == "__main__":
    pass
