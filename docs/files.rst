
.. include:: global.rst  

.. _files:

Input data files
================
ipyrad_ can be used to assemble any type of data that is generated using a 
restriction digest method (RAD, ddRAD, GBS) or related amplification-based 
process (e.g., NextRAD, RApture), both of which yield data that is anchored
on at least one side so that reads are expected to align fairly closely. 
ipyrad is not intended for constructing long contigs from many partially 
overlapping (i.e., shotgun) sequences. It can, however, accomodate paired-end
reads, and has methods for detecting and merging overlaps. ipyrad can also
combine reads of various lengths, so that older data is easily combined 
with newer data of different lengths. 

Depending how and where your sequence data were generated you may receive 
data as one giant file, or in many smaller files. The files may contain data
from all of your individuals mixed up together, or as separate files for each 
Sample. If they are mixed up then the data need to be demultiplexed based on 
barcodes or indices. Input files to ipyrad can be either demultiplexed or not:

**multiplexed (raw) sequence files** -- If your data are not yet sorted 
among individuals/samples then you will need to have barcode 
information organized into a :ref:`barcodes file<barcodes_file>`. 
Sample names are taken from the barcodes file. 
The raw data file path(s) will need to be entered in the 
``raw_fastq_path`` parameter of the params file. 

**demultiplexed (sorted) sequence files** -- If your data are already sorted 
then you simply have to enter the path to the data files in the 
``sorted_fastq_path`` parameter. The :ref:`cookbook recipes <cookbook_recipes>`
section provides more complex methods for combining data from multiple 
sequencing runs into the same individual, or for using multiple barcodes file.


Should you pre-filter data?
---------------------------
ipyrad has methods for filtering based on Q-scores, and for detecting and 
removing barcodes+adapters. For paired-end data ipyrad uses vsearch to test
for merged reads following the similar protocol as the program PEAR. The 
merged and non-merged reads are combined into a single downstream analysis. 


.. _file_names:
File names
-----------
If demultiplexing, then Sample names will be extracted from
the :ref:`barcodes files<barcodes_file>`. Whereas if your data are already 
sorted demultiplexed then Sample names are extracted from the file names 
directly. Do not include spaces in file names. For paired-end data we need
to be able to identify which R1 and R2 files go together, and so we require that
every read1 file name contains the string ``_R1_`` (*with underscores before
and after*), and every R2 file name must match exactly the R1 file
except that it has ``_R2_`` in place of ``_R1_``. 
See the example data for an example. 

.. note:: 

    Pay careful attention to file names at the very beginning of an analysis since 
    these names, and any included typos, will be perpetuated through all the 
    resulting data files. Do not include spaces in file names.


.. _barcodes_file:
Barcodes file
--------------
The barcodes file is a simple table linking barcodes to samples. 
Barcodes can be of varying lengths. 
Each line should have one name and then one barcode, separated by whitespace
(a tab or spaces). 

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
necessary to complete an assembly. A description of how to create and use a 
params file can be found in the :ref:`introductory tutorial<tutorial_intro_cli>`. 


.. _data_types:
Supported data types
--------------------
There is increasingly a large variety of ways to generate reduced representation 
genomic data sets using either restriction digestion or primer sets, and ipyrad 
aims to be flexible enough to handle all of these types. Because it is difficult
to keep up with all of the names, we use our own terminology, described below, 
to group together data types that can be analyzed using the same bioinformatic 
methods. If you have a data type that is not described below and you're not sure
if it can be analyzed in ipyrad_ :ref:`let us know here<gitter>`. 


**rad** -- This category includes data types which use a single cutter to generate 
DNA fragments for sequencing based on a single cut site. *e.g., RAD-seq, NextRAD*. 


**ddrad** -- This category is very similar data types which select fragments that were digested
by two different restriction enzymes which cut the fragment on either end. 
During assembly this type of data is analyzed differently from the **rad** data
type by more stringent filtering that looks for occurrences of the second 
(usually more common) cutter. *e.g., double-digest RAD-seq*. 


**gbs** -- This category includes any data type which selects fragments that were digested
by a **single enzyme** that cuts both ends of DNA fragments. This data type requires 
reverse-complement clustering because the forward vs reverse adapters can attach
to either end of each fragment, and thus when shorter fragments are sequenced 
from either end the resulting reads often overlap partially or completely. 
When analyzing GBS data we strongly recommend using a stringent setting for 
the `filters_adapters` parameter. *e.g., genotyping-by-sequencing 
(Elshire et al.), EZ-RAD (Toonin et al.)*. 


**pairddrad** -- This category is for paired-end data from fragments that were 
generated through restriction digestion using two different enzymes. 
During step 3 the paired-reads will be tested for :ref:`paired 
read merging<paired_read_merging>` if they overlap partially. Because two 
different cutters are used reverse-complement clustering is not necessary. 
*e.g., double-digest RAD-seq (w/ paired-end sequencing)*. 


**pairgbs** -- This category is for paired-end data from fragments that were
generated by digestion with a **single enzyme** that cuts both ends of the
fragment. Because the forward adapter might bind to either end of these
fragments,approximately half of the matches are expected to be reverse-
complemented with perfect overlap. Paired reads are checked for merging before 
clustering/mapping. *e.g., genotyping-by-sequencing, EZ-RAD, (w/ paired-end sequencing)*. 


**2brad** -- This category is for a special class of sequenced fragments 
generated using a type IIb restriction enzyme. The reads are usually very short
in length, and are treated slightly differently in steps 2 and 7. 
(We are looking for people to do more testing of this method on empirical data).

**pair3rad** -- This category is for 3Rad/RadCap data that uses multiplexed
barcodes. 3Rad/RadCap can use up to **four restriction enzymes**, and also
uses a suite of custom adapters to control for PCR duplicates. This data is
always paired end, since one barcode is ligated to each read. PCR clones
are removed in step 3, after merging but before dereplication. The **pair3rad**
`datatype` is used for both 3Rad and RadCap types because these datatypes
only differ in how they are generated, not how they are demultiplexed and
filtered. *See Glenn et al 2016, and Hoffberg et al 2016*
