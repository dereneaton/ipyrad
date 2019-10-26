.. include:: global.rst

.. _parameters:

Assembly: Parameters
====================
The parameters contained in a params file affect the actions that are performed
during each step of an ipyrad assembly. The defaults that we chose are fairly
reasonable values for most assemblies, however, you will always need to modify
at least a few of them (for example, to indicate the location of your data),
and often times you will want to modify many of the parameters.
The ability to easily assemble your data set under a range of parameter
settings is one of the main features of *ipyrad*.

Below is an explanation of each parameter setting, the steps of the assembly
that it effects, and example entries for the parameter into a params.txt file.


.. _params_file:
Parameters (Params File)
------------------------
The parameter input file, which typically includes ``params.txt`` in its name, 
can be created with the ``-n`` option from the ipyrad command line. This file 
lists all of the :ref:`parameter settings<paramater settings>` 
necessary to complete an assembly. A description of how to create and use a 
params file can be found in the :ref:`introductory tutorial<tutorial_intro_cli>`. 

.. code-block:: bash
    
    # create a ipyrad params file from the command line
    >>> ipyrad -n test 


.. code-block:: parsed-literal

    ------- ipyrad params file (v.0.9.1-dev)----------------------------------------
    test                           ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
    /home/deren/ipyrad/tests       ## [1] [project_dir]: Project dir (made in curdir if not present)
                                   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                                   ## [3] [barcodes_path]: Location of barcodes file
                                   ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
    denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                                   ## [6] [reference_sequence]: Location of reference sequence file
    rad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
    TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
    5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
    33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
    6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
    6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
    10000                          ## [13] [maxdepth]: Max cluster depth within samples
    0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
    0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
    0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
    35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
    2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
    0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
    0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
    4                              ## [21] [min_samples_locus]: Min # samples per locus for output
    0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
    8                              ## [23] [max_Indels_locus]: Max # of indels per locus
    0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
    0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
    0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
    p, s, l                        ## [27] [output_formats]: Output formats (see docs)
                                   ## [28] [pop_assign_file]: Path to population assignment file    



.. _assembly_name:

0. Assembly name
-----------------
The Assembly name is used as the prefix for all output files. It should be a
unique identifier for the assembly, meaning the set of parameters you are using
for the current data set. When I assemble multiple data with different parameter
combinations I usually either name them consecutively (e.g., data1, data2), or
with names indicating their parameter combinations (e.g., data_clust90,
data_clust85). The Assembly name cannot be changed after an Assembly is created
with the ``-n`` flag, but a new Assembly with a different name can be created
by branching the Assembly (see :ref:`branching workflow<branching_workflow>`).

Affected steps: 1-7
Example: new Assemblies are created with the -n or -b options to ipyrad:

.. code-block:: bash

    >>> ipyrad -n data1                       ## create a new assembly named data1
    >>> ipyrad -p params-data1.txt -b data2   ## create a branch assembly named data2


.. _project_dir:

1. Project dir
--------------
A project directory can be used to group together multiple related Assemblies.
A new directory will be created at the given path if it does not already exist.
A good name for Project_dir will generally be the name of the organism being studied.
The project dir path should generally not be changed after an analysis is initiated,
unless the entire directory is moved to a different location/machine.

Affected steps: 1-7
Example entries into params.txt:

.. code-block:: bash

    /home/deren/ipyrad/tests/finches   ## [1] create/use project dir called finches
    finches                            ## [1] create/use project dir called finches


.. _raw_fastq_path:

2. Raw fastq path
-----------------
This is a path to the location of raw (non-demultiplexed) fastq data files. If
your data are already demultiplexed then this should be left blank. The input
files can be gzip compressed (i.e., have name-endings with .gz). If you enter
a path for raw data files then you should also enter a path to a
:ref:`barcodes file<barcodes file>`.
To select multiple files, or all files in a directory, use a wildcard character (*).

Affected steps = 1
Example entries into params.txt:

.. code-block:: bash

    /home/deren/ipyrad/tests/data/*.fastq.gz     ## [2] select all gzip data files
    ~/ipyrad/tests/data/*.fastq.gz               ## [2] select all gzip data files
    ./ipsimdata/rad_example*.fastq.gz            ## [2] select files w/ `rad_example` in name


.. _barcodes_path:

3. Barcodes path
----------------
This is a path to the location of a barcodes_file_. This is used in step1
for demuliplexing, and can also be used in step2 to improve the detection of
adapter/primer sequences that should be filtered out. If your data are already
demultiplexed the barcodes path can be left blank.

Affected steps = 1-2.
Example entries into params.txt:

.. code-block:: python

    /home/deren/ipsimdata/rad_example_barcodes.txt    ## [3] select barcode file
    ./ipsimdata/rad_example_barcodes.txt              ## [3] select barcode file


.. _sorted_fastq_path:

4. Sorted fastq path
--------------------
This is a path to the location of sorted fastq data. If your data are already
demultiplexed then this is the location from which data will be loaded when
you run step 1. A wildcard character can be used to select multiple
files in directory.

Affected steps = 1
Example entries into params.txt:

.. code-block:: python

    /home/deren/ipyrad/tests/ipsimdata/*.fastq.gz    ## [4] select all gzip data files
    ~/ipyrad/tests/ipsimdata/*.fastq                 ## [4] select all fastq data files
    ./ipsimdata/rad_example*.fastq.gz                ## [4] select files w/ `rad_example` in name


.. _assembly_method:

5. Assembly method
--------------------
There are four :ref:`Assembly_methods<Assembly_methods>` options in ipyrad:
denovo, reference, denovo+reference, and denovo-reference.
The latter three all require a reference sequence file (param #6) in fasta
format. See the :ref:`tutorials` for an example.

Affected steps = 3, 6
Example entries into params.txt:

.. code-block:: python

    denovo                            ## [5] denovo assembly
    reference                         ## [5] reference assembly
    denovo+reference                  ## [5] reference addition assembly
    denovo-reference                  ## [5] reference subtraction assembly


.. _reference_sequence:

6. Reference sequence
---------------------
The reference sequence file should be in fasta format. It does
not need to be a complete nuclear genome, but could also be any
other type of data that you wish to map RAD data to; for example
plastome or transcriptome data.

.. code-block:: python

    ~/ipyrad/tests/ipsimdata/rad_example_genome.fa   ## [6] select fasta file
    ./data/finch_full_genome.fasta                   ## [6] select fasta file


.. _datatype:

7. Datatype
------------
There are now many forms of restriction-site associated DNA library preparation
methods and thus many differently named data types. Currently, we categorize
these into :ref:`six data types <data_types>`. Follow the link
to deteremine the appropriate category for your data type.


.. code-block:: python

    rad                       ## [7] rad data type (1 cutter, sonication)
    pairddrad                 ## [7] paired ddrad type (2 different cutters)
    pairgbs                   ## [7] paired gbs type (1 cutter cuts both ends)



.. _restriction_overhang:

8. Restriction_overhang
-----------------------
The restriction overhang is used during demultiplexing (step1) and also to detect
and filter out adapters/primers (in step2), if the `filter_adapters` parameter
is turned on. Identifying the correct sequence to enter for the restriction_overhang
can be tricky. You do not enter the restriction recognition sequence, but rather
the portion of this sequence that is left attached to the sequenced read after
digestion. For example, the enzyme PstI has the following palindromic sequence,
with `^` indicating the cut position.

.. code-block:: python

    5'...C TGCA^G...'3
    3'...G^ACGT C...'5

Digestion with this enzyme results in DNA fragments with the sequence ``CTGCA``
adjacent to the cut site, which when sequenced results in the reverse complement
``TGCAG`` as the restriction overhang at the beginning of each read.
The easiest way to identify the restriction overhang is simply to look at
the raw (or demultiplexed) data files yourself. The restriction overhang
will be the (mostly) invariant sequence that occurs at
the very beginning of each read if the data are already demultiplexed, or right
after the barcode in the case of non-demultplexed data. Use the command
below to peek at the first few lines of your fastQ files to find the invariant
sequence.

.. code-block:: bash

    ## gunzip decompresses the file,
    ## the flag -c means print to screen,
    ## and `less` tells it to only print the first 100 lines
    zless 100 my_R1_input_file.fastq.gz

This will print something like the following. You can see that each of the
lines of sequence data begins with `TGCAG` followed by variable sequence data.
For data that used two-cutters (ddrad), you will likely not be able to see the second
cutter overhang for single-end reads, but if your data are paired-end, then
the `_R2_` files will begin with the second `restriction_overhang`. The second
restriction_overhang is only used to detect adapters/primers if the
`filter_adapters`_ parameter is set > 1. The second `restriction_overhang`
can optionally be left blank.

.. parsed-literal::

    @HWI-ST609:152:D0RDLACXX:2:2202:18249:93964 1:N:0:
    TGCAGCAGCAAGTGCTATTCGTACAGTCATCGATCAGGGTATGCAACGAGCAGAAGTCATGATAAAGGGTCCCGGTCTAGGAAGAGACGCAGCATTA
    +
    BDFFHHHHHJJJHIJJJJJIJIGJJJJJJJJJJJJJJJJDGGIJJJJJHIIIJJJJHIIHIGIGHHHHFFFFEDDDDDACCDCDDDDDDDDDBBDC:
    @HWI-ST609:152:D0RDLACXX:2:2202:18428:93944 1:N:0:
    TGCAGGATATATAAAGAATATACCAATCCTAAGGATCCATAGATTTAATTGTGGATCCAACAATAGAAACATCGGCTCAACCCTTTTAGTAAAAGAT
    +
    ADEFGHFHGIJJJJJJIJJIIJJJIJJIJGIJJJJJJJJIJJIJJJJIIIGGHIEGHJJJJJJG@@CG@@DDHHEFF>?A@;>CCCDC?CDDCCDDC
    @HWI-ST609:152:D0RDLACXX:2:2202:18489:93970 1:N:0:
    TGCAGCCATTATGTGGCATAGGGGTTACATCCCGTACAAAAGTTAATAGTATACCACTCCTACGAATAGCTCGTAATGCTGCGTCTCTTCCTAGACC
    +
    BDFFHHHHHJJJJIJJIJJJJJJJHIJJJJJJJJHIIJJJJIFHIJJJJFGIIJFHIJJIJJIFBHHFFDFEBACEDCDDDDBBBDCCCDDDCDDC:
    @HWI-ST609:152:D0RDLACXX:2:2202:18714:93960 1:N:0:
    TGCAGCATCTGGAAATTATGGGGTTATTTCACAGAAGCTGGAATCTCTTGGGCAATTTCACAGAATCTGGGAATATCTGGGGTAAATCTGCAAGATC
    +
    BDEFHHHHHJJJIJJJJJJJJJJCGIIJJGHJJJJJJJJJJIJJIJJJIHIJJJJJJJHHIIJJJJJJJHGHHHFEFFFDEDABDDFDDEDDDDDDA
    @HWI-ST609:152:D0RDLACXX:2:2202:20484:93905 1:N:0:
    TGCAGAGGGGGATTTTCTGGAGTTCTGAGCATGGACTCGTCCCGTTTGTGCTGCTCGAACACTGACGTTACTTCGTTGATCCCTATGGACTTGGTCA
    +
    ADEEHHHHGJJHIIJJJJJJIJFHJIIIJJJJIIIJJJJEHHIIGHIGGHGHHHHFFEEEDDDD;CBBDDDACC@DC?<CDDCCCCCCA@CCD>>@:
    @HWI-ST609:152:D0RDLACXX:2:2202:20938:93852 1:N:0:
    TGCAGAAGCTGGAGATTCTGGGGCAGCTTTGCAGCAAGCTGAAAATTCTGGGGGTCGATCTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

In some cases restriction enzymes can bind to more than one specific sequence,
for example ApoI will bind to AATTY (i.e. AATTC and AATTT). If you used an 
enzyme with reduced specificity you can include ambiguity codes in the 
restriction overhang sequence.

Affected steps = 1,2. Example entries to params.txt file:

.. code-block:: bash

    TGCAG                     ## [8] single cutter (e.g., rad, gbs)
    TGCAG, AATT               ## [8] double digest (e.g., ddrad, pairddrad)
    CWGC                      ## [8] single cutter w/ degenerate base


NB: 3RAD and SeqCap data can use up to 4 restriction enzymes. If you have this
kind of data, simply list all the restriction overhangs for all your cutters.

.. code-block:: bash

    CTAGA, CTAGC, AATTC               ## [8] 3rad data (multiple cutters)


.. _max_low_qual_bases:

9. max_low_qual_bases
---------------------
During step 2 bases are trimmed from the 3' end of reads when the quality score
is consistently below 20 (which can be modified by modifying phred_Qscore_offset_). 
However, your reads may still contain some number of ambiguous (N) sites that 
were not trimmed based on quality scores, and these will affect the efficiency
and accuracy of clustering downstream. This parameter sets the upper limit on the 
number of Ns allowed in reads. The default value for
`max_low_qual_bases` is 5. I would generally recommend against increasing 
this value greatly. 

Affected steps = 2. Example entries to params.txt:

.. code-block:: bash

    0                      ## [9] allow zero low quality bases in a read
    5                      ## [9] allow up to five low quality bases in a read


.. _phred_Qscore_offset:

10. Phred_Qscore_offset
------------------------
Bases are trimmed from the 3' end of reads if their quality scores is below
this 20. The default offset for quality scores is 33. Some 
older data use a qscore offset of 64, but this is increasingly rare. You
can toggle the offset number to change the threshold for trimming. 
For example, reducing the offset from 33 to 23 is equivalent to changing the 
minimum quality score from 20 to 10, which is approximately 
95% probability of a correct base call.

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    33                 ## [10] default offset of 33, converts to min score=20
    43                 ## [10] offset increased by 10, converts to min score=30
    64                 ## [10] offset used by older data, converts to min score=20.


.. _mindepth_statistical:

11. mindepth_statistical
-------------------------
This is the minimum depth at which statistical base calls will be made during
step 5 consensus base calling. By default this is set to 6, which for most
reasonable error rates estimates is approximately the minimum depth at which a
heterozygous base call can be distinguished from a sequencing error.

Affected steps = 4, 5. Example entries to params.txt

.. parsed-literal::

    6                 ## [11] set mindepth statistical to 6
    10                ## [11] set to 10


.. _mindepth_majrule:

12. mindepth_majrule
---------------------
This is the minimum depth at which majority rule base calls are made during
step 5 consensus base calling. By default this is set to the same value as
mindepth_statistical, such that only statistical base calls are made. This
value must be <= mindepth_statistical. If lower, then sites with coverage
>= mindepth_majrule and < mindepth_statistical will make majority rule calls.
If your data set is very low coverage such that many clusters are excluded due
to low sequencing depth then lowering mindepth_majrule can be an effective way
to increase the amount of usable information in your data set. However, you
should be aware the majority rule consensus base calls will underestimate
heterozygosity.

Affected steps = 4, 5. Example entries to params.txt:

.. parsed-literal::

    6                 ## [12] set to relatively high value similar to mindepth_stat
    2                 ## [12] set below the statistical limit for base calls.


.. _maxdepth:
13. maxdepth
-------------
Sequencing coverage is often highly uneven among due to differences in the
rate at which fragments are amplified during library preparation, the extent
to which varies across different library prep methods. Moreover, repetitive
regions of the genome may appear highly similar and thus cluster as high depth
clusters. Setting a maxdepth helps to remove the latter problem, but at the
expense of potentially removing good clusters that simply were sequenced
to high depth. The default maxdepth is set quite high (10,000), but you may
change it as you see fit.

Affected steps = 4, 5. Example entries to params.txt:

.. parsed-literal::

    10000             ## [13] maxdepth above which clusters are excluded.



.. _clust_threshold:

14. clust_threshold
--------------------
This the level of sequence similarity at which two sequences are identified
as being homologous, and thus cluster together. The value should be entered
as a decimal (e.g., 0.90). We do not recommend using values higher than 0.95,
as homologous sequences may not cluster together at such high threshold due
to the presence of Ns, indels, sequencing errors, or polymorphisms.

Affected steps = 3, 6. Example entries to params.txt:

.. parsed-literal::

    0.90              ## [14] clust threshold set to 90%
    0.85              ## [14] clust threshold set to 85%


.. _max_barcodes_mismatch:

15. max_barcodes_mismatch
---------------
The maximum number of allowed mismatches between the barcodes in the barcodes
file and those found in the sequenced reads. Default is 0. Barcodes usually differ
by a minimum of 2 bases, so I would not generally recommend using a value >2.

Affected steps = 1. Example entries to params.txt:

.. parsed-literal::

    0              ## [15] allow no mismatches
    1              ## [15] allow 1 mismatched base


.. _filter_adapters:

16. filter_adapters
--------------------
It is important to remove Illumina adapters from your data if present. 
Depending on the fidelity of the size selection procedure implemented during
library preparation there is often at least some small proportion of sequences
in which the read length is longer than the actual DNA fragment, such that the
primer/adapter sequence ends up in the read. This occurs more commonly in 
double-digest (GBS, ddRAD) data sets that use a common cutter, and can be 
especially problematic for GBS data sets, in which short fragments are sequenced 
from either end. The `filter_adapters` parameter has three settings (0, 1, or 2). 
**If 0**, then reads are only removed if they contain more Ns than allowed by the
`max_low_qual_bases` parameter. **If 1**, then reads are trimmed to the first base
which has a Qscore < 20 (on either read for paired data), and also removed if there
are too many Ns. **If 2**, then reads are searched for the common Illumina adapter, 
plus the reverse complement of the second cut site (if present), plus the barcode
(if present), and this part of the read is trimmed. This filter is applied using 
code from the software `cutadapt`, which allows for errors within the adapter 
sequence. 

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    0                ## [16] No adapter filtering
    1                ## [16] filter based on quality scores
    2                ## [16] strict filter for adapters


.. _filter_min_trim_len:

17. filter_min_trim_len
------------------------
During step 2 if `filter_adapters` is > 0 reads may be trimmed to a shorter length
if they are either low quality or contain Illumina adapter sequences. 
By default ipyrad will keep trimmed reads down to a minimum length of 35bp. 
If you want to set a higher limit you can do so here.

Affected steps = 2. Example entries to params.txt

.. parsed-literal::

    50                ## [17] minimum trimmed seqlen of 50
    75                ## [17] minimum trimmed seqlen of 75


.. _max_alleles_consens:

18. max_alleles_consens:
-------------------------
This is the maximum number of unique alleles allowed in (**individual**) consens
reads after accounting for sequencing errors. Default=2, which is fitting for
diploids.
At this setting any locus which has a sample with more than 2 alleles detected
will be  excluded/filtered out. If set to max_alleles_consens = 1 (haploid)
then error-rate and heterozygosity are estimated with H fixed to 0.0 in step 4,
and base calls are made with the estimated error rate, and any consensus reads
with more than 1 allele present are excluded.
If max_alleles_consens is set > 2 then more alleles are allowed, however,
heterozygous base calls are still made under the assumption of diploidy
i.e., hetero allele frequency=50%.

Affected steps = 4, 7. Example entries to params.txt

.. parsed-literal::

    2                ## [18] diploid base calls, exclude if >2 alleles
    1                ## [18] haploid base calls, exclude if >1 allele
    4                ## [18] diploid-base calls, exclude if >4 alleles

.. _max_Ns_consens:


19. max_Ns_consens:
--------------------
The maximum fraction of uncalled bases allowed in consens seqs. If a
base call cannot be made confidently (statistically) then it is called
as ambiguous (N). You do not want to allow too many Ns in consensus reads
or it will affect their ability to cluster with consensus reads from other
Samples, and it may represent a poor alignment. Default is 0.05.

Affected steps = 5. Example entries to params.txt

.. parsed-literal::

    0.1                ## [19] allow max of 10% Ns in a consensus seq


.. _max_Hs_consens:

20. max_Hs_consens:
--------------------
The maximum fraction of heterozygous bases allowed in consens seqs.
This filter helps to remove poor alignments which will tend to have an
excess of Hs. Default is 0.05.

Affected steps = 5. Example entries to params.txt

.. parsed-literal::

    0.05                ## [20] allow max of 5% Hs in a consensus seq


.. _min_samples_locus:

21. min_samples_locus
----------------------
The minimum number of Samples that must have data at a given locus for it to
be retained in the final data set. If you enter a number equal to the full
number of samples in your data set then it will return only loci that have
data shared across all samples. Whereas if you enter a lower value, like 4, 
it will return a more sparse matrix, including any loci for which at least 
four samples contain data. This parameter is overridden if a min_samples values 
are entered in the :ref:`popfile<pop_assign_file>`. Default value is 4.

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    4                ## [21] create a min4 assembly
    12               ## [21] create a min12 assembly


.. _max_SNPs_locus:

22. max_SNPs_locus
-------------------
Maximum number of SNPs allowed in a final locus.
This can remove potential effects of poor alignments in repetitive regions
in a final data set by excluding loci with more than N snps. Setting lower
values is likely only helpful for extra filtering of very messy data sets.
The default is 0.2

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    0.2               ## [22] allow max of 20% SNPs per locus.
    0.05              ## [22] allow max of 5% SNPs locus.


.. _max_Indels_locus:

23. max_Indels_locus
---------------------
The maximum number of Indels allowed in a final locus. This helps to filter
out poor final alignments, particularly for paired-end data.
The default is 8.

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    5                ## [23] allow max of 5 indels per locus.


.. _max_shared_Hs_locus:

24. max_shared_Hs_locus
------------------------
Maximum number (or proportion) of shared polymorphic sites in a locus.
This option is used to detect potential paralogs, as a shared heterozygous
site across many samples likely represents clustering of paralogs with a
fixed difference rather than a true heterozygous site. Default is 0.5.

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    0.25             ## [24] allow hetero site to occur across max of 25% of Samples


.. _trim_reads:

25. trim_reads
--------------
Sometimes you can look at your fastq data files and see that there was a problem
with the sequencing such that the cut site which should occur at the beginning of
your reads is either offset by one or more bases, or contains many errors. You 
can trim off N bases from the beginning or end of R1 and R2 reads during step 2 
by setting the number of bases here. This could similarly be used to trim all 
reads to a uniform length (though uniform read lengths are not required in ipyrad).

Affected steps = 2. Example entries to params.txt

.. parsed-literal::

    0, 0, 0, 0       ## [25] does nothing
    5, 0, 0, 0       ## [25] trims first 5 bases from R1s
    5, -5, 0, 0      ## [25] trims first 5 bases and last five based from R1s
    5, 80, 0, 0      ## [25] trims first 5 bases from R1s and trims maxlen to 80
    5, 75, 5, 75     ## [25] trims first 5 from R1 and R1, and maxlen to 75.


.. _trim_loci:

26. trim_loci
------------------
Trim N bases from the edges of final aligned loci. This can be useful in denovo
data sets in particular, where the 3' edge of reads is less well aligned than the
5' edge, and thus error rates are sometimes higher at the ends of reads. 

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    0, 0, 0, 0     ## [26] no locus edge trimming
    0, 5, 0, 0     ## [26] trims first 5 bases from R1s in aligned locus
    0, 5, 5, 0     ## [26] trims last 5 bases from R1s and first 5 from R2s


.. _output_formats:

27. output_formats
-------------------
Disk space is cheap, and these are quick to make, so by default we make all
formats. More are coming (alleles, treemix, migrate-n, finestructure). 
The short list of available options is below but see 
:ref:`output formats<full_output_formats>` section for full descriptions of 
the available formats. 

.. code-block:: python

    p: PHYLIP (Full dataset)
    s: PHYLIP (SNPs only)
    u: PHYLIP (One SNP per locus)
    n: NEXUS
    k: STRUCTURE 
    g: EIGENSTRAT .geno
    G: G-PhoCS
    v: VCF (SNPs only)

Affected steps = 7. Example entries to params.txt

.. parsed-literal::

    *                     ## [27] Make all output datatypes
    n, v, g               ## [27] Only write out nexus, vcf and geno formats
    u,k                   ## [27] Only write out unlinked snps in phylip, and structure



.. _pop_assign_file:

28. pop_assign_file
--------------------
Population assignment file for creating population output files, or 
assigning min_samples_locus value to each population. Enter a path to 
the file. (see below for details of the file).

Affected step: 7. Example entries to params.txt

.. parsed-literal::

    /home/user/ipyrad/popfile.txt        ## [28] example...

The population assignment file should be formatted as a plain-txt, whitespace
delimited list of individuals and population assignments. Care should be taken
with spelling and capitalization. Each line should contain a sample name followed
by a population name to which that sample is assigned. One or more additional 
lines should be included that start with one or more "#" characters. These 
special lines tell ipyrad how many samples must have data within each population
for the locus to be retained in the final assembly, and thus assign different 
min_samples_locus values to each population. This will override the global
min_samples_locus value. 

See the example below.

.. parsed-literal::

    Sample1 pop1
    Sample2 pop1
    Sample3 pop1
    Sample4 pop2
    Sample5 pop2

    # pop1:2 pop2:2
