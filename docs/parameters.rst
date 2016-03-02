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
data_clust85). The Assembly name cannot be changed after an Assembly is created, 
but a new Assembly with a different name can be created by copying (branching)
the Assembly (see :ref:`branching workflow<branching_workflow>`).

Affected steps: 1-7  

Example entries into params.txt:  

.. code-block:: python

    data1                      ## [0] name the Assembly data1
    clust90_minsamp4           ## [0] name the Assembly based on some param settings


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

    TGCAG                     ## single cutter (e.g., rad, gbs)
    TGCAG, AATT               ## double digest (e.g., ddrad, pairddrad) 
    CWGC                      ## single cutter w/ degenerate base




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
    
    0                      ## allow zero low quality bases in a read
    4                      ## allow up to four low quality bases in a read


.. _phred_Qscore_offset:

10. Phred_Qscore_offset
------------------------
The threshold at which a base call is considered a low quality base call 
during step 2 filtering is determined by the `phred_Qscore_offset`. The 
default offset is 33, which is equivalent to a minimum qscore of 20. Some 
older data use a qscore offset of 64. You can toggle the offset number 
to change the threshold for low qual bases.

Affected steps = 2. Example entries to param.txt:

.. code-block:: bash

    33                 ## default offset of 33, converts to min score=20
    23                 ## offset reduced by 10, converts to min score=30
    64                 ## offset used by older data, converts to min score=20.


.. _mindepth_statistical:

11. Mindepth_statistical
-------------------------



.. _mindepth_majrule:

12. Mindepth_majrule
---------------------

.. _maxdepth:
13. Maxdepth
-------------
...

.. _clust_threshold:
14. Clust_threshold
--------------------
...




