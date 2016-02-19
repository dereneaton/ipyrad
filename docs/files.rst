
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
FASTQ files -- raw or sorted data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Data can be already demultiplexed (sorted among individuals/samples) or it might
be in the form of 'multiplexed' data, in which case you will need the barcode
information to sort Samples out. ...


.. _file_names:
Input file names
^^^^^^^^^^^^^^^^^
blah blah ``_R1_`` and ``_R2_``. Sample names come from file names. 

* Link to recipe: renaming Samples different from file names.
* Link to recipe: combining multiple input files for individual Samples. 


.. _barcodes_file:
Barcodes file
^^^^^^^^^^^^^^
The barcodes file is a simple table linking barcodes to samples. 
Barcodes can be of varying lengths. 
Each line should have one name and then one barcode, separated by a tab or 
space.


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