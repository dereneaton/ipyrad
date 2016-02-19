
.. _files:


Input data/files
=================
ipyrad_ can be used to assemble any kind of data that is generated using a 
restriction digest method (RAD-seq) or 


.. _data_types:
Supported data types
^^^^^^^^^^^^^^^^^^^^^
There are increasingly a huge range of ways to generate reduced representation 
genomic data sets using either restriction digestion or a set of primers or baits, 
many of which can be assembled in ipyrad. Because it is difficult to keep up with 
all of the names, we use our own terminology, described below, to group together
data types that can be analyzed using the same methods. If you have a data type
that isn't described below and your're not sure if it can be analyzed in ipyrad
:ref:`let us know here<gitter>`. 

is some confusion in the literature about the names of 


**RAD-seq**
RAD

**dd-RAD**
double-digest RAD-seq describes . This includes what double-digest GBS. 

**GBS**
This includes EZ-RAD, ... and others.

**paired-GBS**
This includes paired-end EZ-RAD...

**paired-GBS**
...


.. _input_files:
FASTQ input files
^^^^^^^^^^^^^^^^^^^^
Depending on how and where your sequence data are generated you may receive the
data in a single giant file, or in many smaller files. The files may contain data
from all of your individuals still mixed up together, in which case the data need
to be demultiplexed based on their attached barcodes or index; or your data may 
already be demultiplexed, in which case each of your data files corresponds to 
a different sample. 

** multiplexed (raw) sequence files **  
If your data are not yet sorted among individuals/samples then you will need 
to have their barcodes information organized into a barcodes_file_. Sample names 
will be taken from the barcodes file. 

** demultiplexed (sorted) sequence files **  
If your data are already sorted then you simply have to enter the path to the 
data files in the ``sorted_fastq_path`` parameter of ipyrad. 
Sample names come from file names. 

.. note:: The file names matter.


.. _file_names:
Input file names
^^^^^^^^^^^^^^^^^
If you are using a paired-end data type then the rules for filenames are much 
more strict than for single-end data. Every read1 file must contain the string 
``_R1_`` in it, and every R2 file must match exactly to the name of the R1 file
except that it has ``_R2_``. 

* Link to recipe: renaming Samples different from file names.
* Link to recipe: combining multiple input files for individual Samples. 


.. _barcodes_file:

Barcodes file
^^^^^^^^^^^^^^
The barcodes file is a simple table linking barcodes to samples. 
Barcodes can be of varying lengths. 
Each line should have one name and then one barcode, separated by a tab or 
space. The names that you enter in the barcodes file are the names 
that will end up in your output files, so it is useful to check for 
typos or other errors, or to shorten the names as you see fit before 
running step1. 

.. parsed-literal:: 
    sample1     ACAGG
    sample2     ATTCA  
    sample3     CGGCATA  
    sample4     AAGAACA  


.. _params_file:
Params file
^^^^^^^^^^^^
The parameter input file, which will typically include ``params.txt`` in its name, 
can be created with the ``-n`` option from the ipyrad command line. This file 
lists all of the :ref:`paramater settings` necessary to run an assembly.
A description of how to create and use a parmas file can be found in the 
:ref:`introductory tutorial<tutorial_intro_cli>`. 