
.. include:: global.rst  


Getting Started: Demultiplexing
================================
Demultiplexing is the process of sorting sequenced reads into separate files for each sample in a sequenced run. You may have received your data already demultiplexed, with a separate file for each sample. If so, then you can proceed to the next section. If your data are not yet sorted into separate files then you will need to perform *demultiplexing* during step 1 of the ipyrad assembly.


Multiplexing and Multiple Libraries
-----------------------------------
If your data are not yet sorted among individuals/samples then you will need to have barcode/index information organized into a 
:ref:`barcodes file<barcodes_file>` to sort data to separate files for each sample. *ipyrad* has several options for demultiplexing by internal barcodes or external i7 indices, and for combining samples from many different sequencing runs together into a single analysis, or splitting them into separate analyses, as well as for merging data from multiple sequenced lanes into the same sample names (e.g., technical replicates). See the Demultiplexing section for simple examples, and the Cookbook section for further detailed examples.


If demultiplexing then Sample names will be extracted from
the :ref:`barcodes files<barcodes_file>`, whereas if your data are already demultiplexed then Sample names are extracted from the file names directly. Do not include spaces in file names. For paired-end data we need to be able to identify which R1 and R2 files go together, and so we require that every read1 file name contains the string ``_R1_`` (*with underscores before and after*), and every R2 file name must match exactly the R1 file except that it has ``_R2_`` in place of ``_R1_``. See the tutorial data files for an example. 

.. note:: 

    Pay careful attention to file names at the very beginning of an analysis since these names, and any included typos, will be perpetuated through all the resulting data files. Do not include spaces in file names.


.. _sample_names:

Sample Names
-------------
When demultiplexing Sample names will be extracted from
the :ref:`barcodes files<barcodes_file>` whereas if your data are already demultiplexed then Sample names are extracted from file names 
directly. Do not include spaces in file names. For paired-end data we need to be able to identify which R1 and R2 files go together, and so we require that every read1 file name contains the string ``_R1_`` (*with underscores before and after*), and every R2 file name must match exactly the R1 file except that it has ``_R2_`` in place of ``_R1_``. 
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



Combinatorial indexing
-----------------------
To perform combinatorial indexing you will need to enter two barcodes for 
each sample name. These should be ordered so that the barcode on read1 is 
first and the barcode on read2 second. A simple way to ensure that barcodes
are attached to your reads in the way that you expect is to look at the raw 
data files (e.g., use the command line tool `less`) and check for the 
barcode sequences. 


.. parsed-literal:: 
    sample1     ACAGG	TTCCG
    sample2     ATTCA	CCGGAA
    sample3     CGGCAT	GAGTCC
    sample4     AAGAAC	CACCG


i7 indexing
-------------
External barcodes/indexes can also be attached external to the sequenced read
on the Illumina adapters. This is often used to combine multiple plates together
onto a single sequencing run. You can find the i7 index in the header line of 
each read in a fastq file. *ipyrad* can demultiplex using i7 indices if you 
turn on a special flag. An example of how to do this using the *ipyrad* API 
is available in the cookbook section.


.. parsed-literal:: 
    lib1     CCGGAA
    lib2     AATTCC


Combining multiple libraries
----------------------------
With *ipyrad* it is very easy to combine multiple sequenced libraries into a
single assembly. This is accomplished by demultiplexing each lane of data 
separately and then combining the libraries using *merging*. See the merging
section for details and examples in the cookbook section.
