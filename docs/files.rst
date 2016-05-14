
.. include:: global.rst  

.. _files:

Input data/files
=================
ipyrad_ can be used to assemble any type of data that is generated using a 
restriction digest method (RAD, ddRAD, GBS) or related amplification-based 
process (e.g., NextRAD, RApture), both of which yield data that is anchored
on at least one side so that reads are expected to align fairly closely. 
ipyrad is not intended for constructing long contigs from many partially 
overlapping (i.e., shotgun) sequences, however, ipyrad can accomodate paired-end
data with or without merged reads, and works with reads of various lengths, so 
that older data can be easily combined with new data of different lengths. 

Depending on how and where your sequence data were generated you may receive the
data in one giant file, or in many smaller files. The files may contain data
from all of your individuals mixed up together, in which case they need
to be demultiplexed based on barcodes or indices; or the data may 
already be demultiplexed, in which case each of your data files corresponds to 
a different Sample. 


__`multiplexed (raw) sequence files`__
If your data are not yet sorted among individuals/samples then you will need 
to have their barcode information organized into a 
:ref:`barcodes file<barcodes_file>`. Sample names are taken from the barcodes 
file. The raw data file(s) should be entered in the ``raw_fastq_path`` parameter
of the params file. 

demultiplexed (sorted) sequence files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If your data are already sorted then you simply have to enter the path to the 
data files in the ``sorted_fastq_path`` parameter.
The :ref:`cookbook recipes <cookbook_recipes>` section provides more complex
methods for combining data from multiple sequencing runs into the same 
individual, or for using multiple barcodes file.


.. note:: 

    It's worth paying careful attention to file names before starting
    an analysis since these names, and any included typos, will be perpetuated 
    through all the resulting data files. Do not include spaces in file names.


.. _file_names:
Input file names
-----------------
If your data are not yet demultiplexed then Sample names will come from the 
:ref:`barcodes files<barcodes_file>`, as shown below. 
Otherwise, if data files are already sorted among Samples (demultiplexed) 
then Sample names will be extracted from the file names. 
The file names should not have any spaces in them. 
If you are using a paired-end data type then the rules for file names are a bit 
more strict than for single-end data. Every read1 file must contain the string 
``_R1_`` in it, and every R2 file must match exactly to the name of the R1 file
except that it has ``_R2_``. See the tutorials for an example. 



.. _barcodes_file:

Barcodes file
--------------
The barcodes file is a simple table linking barcodes to samples. 
Barcodes can be of varying lengths. 
Each line should have one name and then one barcode, separated by a tab or 
space. The names that you enter in the barcodes file are the names 
that will end up in your output files, so it is useful to check for 
typos or other errors, or to shorten the names as you see fit before 
running step1. Do not include any spaces in Sample names. 

.. parsed-literal:: 
    sample1     ACAGG
    sample2     ATTCA  
    sample3     CGGCATA  
    sample4     AAGAACA  


.. _params_file:
Params file
------------
The parameter input file, which typically includes ``params.txt`` in its name, 
can be created with the ``-n`` option from the ipyrad command line. This file 
lists all of the :ref:`parameter settings<paramater settings>` 
necessary to complete an assembly. 
A description of how to create and use a parmas file can be found in the 
:ref:`introductory tutorial<tutorial_intro_cli>`. 



.. _data_types:
Supported data types
--------------------

There is increasingly a large variety of ways to generate reduced representation 
genomic data sets using either restriction digestion or primer sets, 
many of which can be assembled in ipyrad_. Because it is difficult to keep up with 
all of the names, we use our own terminology, described below, to group together
data types that can be analyzed using the same bioinformatic methods. 
If you have a data type that is not described below and you're not sure if it 
can be analyzed in ipyrad_ :ref:`let us know here<gitter>`. 


rad 
^^^^
This category includes data types which use a single cutter to generate 
DNA fragments for sequencing based on a single cut site. 

e.g., RAD-seq, NextRAD


ddrad 
^^^^^
This category includes data types which select fragments that were digested
by two different restriction enzymes which cut the fragment on either end. 
During assembly this type of data is analyzed differently from the **rad** data
type by more stringent filtering that looks for occurrences of the second 
(usually more common) cutter. 

e.g., double-digest RAD-seq


gbs
^^^
This category includes any data type which selects fragments that were digested
by a **single enzyme** that cuts both ends of DNA fragments. This data type requires 
reverse-complement clustering because the forward vs reverse adapters can attach
to either end of each fragment, and thus when shorter fragments are sequenced 
from either end the resulting reads often overlap partially or completely. 
When analyzing GBS data we strongly recommend using a stringent setting for 
the `filters_adapters` parameter.

e.g., genotyping-by-sequencing (Elshire et al.), EZ-RAD (Toonin et al.)


pairddrad
^^^^^^^^^
This category is for paired-end data from fragments that were generated 
through restriction digestion using two different enzymes. During step 3 the 
paired-reads will be tested for :ref:`paired read merging<paired_read_merging>`
if they overlap partially. 

e.g., double-digest RAD-seq (w/ paired-end sequencing)


pairgbs
^^^^^^^
This category is for paired-end data from fragments that were generated by 
digestion with a **single enzyme** that cuts both ends of the fragment. 
Because the forward adapter might bind to either end of these fragments,
approximately half of the matches are expected to be reverse-complemented 
with perfect overlap. Paired reads are checked for merging before clustering/mapping.


e.g., genotyping-by-sequencing, EZ-RAD, (w/ paired-end sequencing)


2brad
^^^^^^
This category is for a special class of reads sequenced fragments generated using
a type IIb restriction enzyme. The reads are usually very short in length, and 
are treated slightly differently in steps 2 and 7. 
