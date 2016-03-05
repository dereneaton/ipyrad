.. include:: global.rst

.. _parameters:

Assembly Parameters
====================

ipyrad_ performs a series of :ref:`7 steps <7_steps>` 
to assemble RAD-seq data sets by sorting
filtering, clustering, and formatting the data into output files. During each 
step a number of parameters are used which affect how that step is performed. 
The defaults that we chose are fairly reasonable values for most assemblies, 
however, you will always need to modify at least a few of them (for example, to 
indicate the location of your data), and often times you will want to modify 
many of the parameters. The ability to easily assemble your data set under a range
of parameter settings is one of the main features of ipyrad_. 

Below is an explanation of each parameter setting, the steps of the assembly 
that it effects, and example entries for the parameter into a params.txt file.


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
by copying (branching) the Assembly 
(see :ref:`branching workflow<branching_workflow>`).

Affected steps: 1-7  

Example entries into params.txt:  

.. code-block:: python

    data1                      ## [0] name the Assembly data1
    clust90_minsamp4           ## [0] name based on some param settings


.. _project_dir:

1. Project dir
--------------
The Project directory is the location where a group of Assemblies which share
data files will be saved. This can be either a name or a path. If it is a path
then the a new directory will be created at the given path if it does not already
exist. If a name is entered then a new directory with that name will be created in the
current directory if it does not already exist. A good name for Project_dir will
generally be the name of the organism being studied.

Affected steps: 1-7  

Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/finches   ## [1] create/use project dir called finches
    finches                            ## [1] create/use project dir called finches


.. _raw_fastq_path:

2. Raw fastq path
-----------------
This is a path to the location of raw (non-demultiplexed) fastq data files. If
your data are already demultiplexed then this should be left blank. The input
files can be gzip compressed (i.e., have name-endings with .gz). If you enter
a path for raw data files then you should also have a path to a 
:ref:`barcodes file<barcodes file>` for parameter #3 (`barcodes path`_). 
To select multiple files, or all files in a directory, use a wildcard character (*).

Affected steps = ``1``. Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/*.fastq.gz  ## [2] select all gzip data files
    ~/ipyrad/tests/data/*.fastq.gz            ## [2] select all gzip data files
    ./data/sim_rad*.fastq.gz                  ## [2] select files w/ `sim_rad` in name


.. _barcodes_path:

3. Barcodes path
----------------
This is a path to the location of the barcodes_file_. This is used in step1
for demuliplexing, and can also be used in step2 to improve the detection of
adapter/primer sequences that should be filtered out. If your data are already
demultiplexed the barcodes path can be left blank. 

Affected steps = ``1, 2``. Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/sim_barcodes.txt  ## [3] select barcode file
    ~/tests/data/sim_barcodes.txt                   ## [3] select barcode file


.. _sorted_fastq_path:

4. Sorted fastq path
--------------------
This is a path to the location of sorted fastq data. If your data are already
demultiplexed then this is the location that will be used in step 1 to load 
the data into ipyrad. A wildcard character can be used to select multiple 
files in directory. 

Affected steps = ``1``. Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/*.fastq.gz  ## [4] select all gzip data files
    ~/ipyrad/tests/data/*.fastq               ## [4] select all fastq data files
    ./data/sim_rad*.fastq.gz                  ## [4] select files w/ `sim_rad` in name


.. _assembly_method:

5. Assembly method
--------------------
There are four :ref:`Assembly_methods<Assembly_methods>` options in ipyrad_: 
denovo, reference, denovo+reference, and denovo-reference. 
The latter three all require a reference sequence file (param #6) in fasta 
format. See the :ref:`tutorials` for an example. 

Affected steps: ``3, 6``. Example entries into params.txt:  

.. code-block:: python

    denovo                            ## [5] denovo assembly
    reference                         ## [5] reference assembly
    denovo+reference                  ## [5] reference addition assembly
    denovo-reference                  ## [5] reference subtraction assembly


.. _reference_sequence:

6. Reference sequence
---------------------
The reference sequence file should be in fasta format. It does 
not need to be a complete nuclear genome, but can also be useful for 
filtering/selecting plastome or transcriptome data, 
or for other uses as well. 

.. code-block:: python

    ~/ipyrad/tests/data/sim_mt_genome.fasta   ## [6] select fasta file
    ./data/finch_full_genome.fasta            ## [6] select fasta file


.. _datatype:

7. Datatype
------------
There are now many forms of restriction-site associated DNA library preparation
methods and thus many differently named data types. Currently, we categorize 
these into :ref:`six data types <Supported data types>`. Follow the link
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
    gunzip -c my_R1_input_file.fastq.gz | less 100

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


Affected steps = 1,2. Example entries to params.txt file:

.. code-block:: bash

    TGCAG                     ## [8] single cutter (e.g., rad, gbs)
    TGCAG, AATT               ## [8] double digest (e.g., ddrad, pairddrad) 
    CWGC                      ## [8] single cutter w/ degenerate base




.. _max_low_qual_bases:

9. max_low_qual_bases
---------------------
During step 2 low quality base calls are converted to Ns, and reads with more
than some number of low quality base calls are excluded. The default value for 
`max_low_qual_bases` is 4. Allowing too many Ns in sequences will reduce the 
likelihood that homologous reads cluster together, so I do not recommend 
increasing this value above the number of allowed differences between reads
based on the `clust_threshold` and sequence length. 

Affected steps = 2. Example entries to params.txt:

.. code-block:: bash
    
    0                      ## [9] allow zero low quality bases in a read
    4                      ## [9] allow up to four low quality bases in a read


.. _phred_Qscore_offset:

10. Phred_Qscore_offset
------------------------
The threshold at which a base call is considered a low quality base call 
during step 2 filtering is determined by the `phred_Qscore_offset`. The 
default offset is 33, which is equivalent to a minimum qscore of 20. Some 
older data use a qscore offset of 64. You can toggle the offset number 
to change the threshold for low qual bases.

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    33                 ## [10] default offset of 33, converts to min score=20
    23                 ## [10] offset reduced by 10, converts to min score=30
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
This the level of sequene similarity at which two sequences are identified
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
file and those found in the sequenced reads. Default is 1. Barcodes usually differ
by a minimum of 2 bases, so I would not generally recommend using a value >2. 

Affected steps = 1. Example entries to params.txt:

.. parsed-literal::

    0              ## [15] allow no mismatches
    1              ## [15] allow 1 mismatched base


.. _filter_adapters:

16. filter_adapters
--------------------
Depending on the fidelity of the size selection procedure implemented during 
library preparation there are usually at least some proportion of sequences
in which the read length is longer than the actual DNA fragment, such that the 
primer/adapter sequence ends up in the read. This can be a problem and should
be filtered out. It occurs more commonly in double-digest data sets that use 
a common cutter, and can be espeically problematic for gbs data sets, in 
which short fragments are sequenced from either end.
If filter_adapters is set to 0 then no check is performed for sequence adapters.
If it is set to 1 then step 2 performs a fuzzy check for the reverse 
complement match of the cut site followed by the beginning of the adapter 
sequence. If it is set to 2 is performs an even fuzzier match to catch 
adapters even if there are sequencing errors or Ns in them, but has a higher
rate of false positive matches. 

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    0                ## [16] No adapter filtering
    1                ## [16] filter for adapters
    2                ## [16] strict filter for adapters


.. _filter_min_trim_len:

17. filter_min_trim_len
------------------------




