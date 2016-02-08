.. include:: global.rst

.. _parameters:

Assembly Parameters
====================

To assemble data in ipyrad_ a series of 7 functions (steps_) are performed to
sort, filter, cluster, and format the data into output files. Just how these 
different steps function is controlled by a set of parameters that can be 
modified by users. 

The default settings for each parameter are fairly reasonable values for most
assemblies, but you will always need to modify some of them (for example, to 
indicate the location of your data), and often times you will want to modify 
many others as well. The ability to easily assemble your data set under a range
of parameter settings is one of the strengths of ipyrad. 

Below is an explanation of each parameter setting and which steps of the 
assembly it effects. 


0. Assembly name
-----------------
The Assembly name is used as the prefix for all output files. It should be a 
unique identifier for the Assembly, meaning the set of parameters you are using 
for the current data set. When I assemble multiple data with different parameter
combinations I usually either name them consecutively (e.g., data1, data2), or 
with names indicating their parameters combinations (e.g., data_clust90, 
data_clust85). The Assembly name cannot be changed after an Assembly is created, 
but a new Assembly with a different name can be created by copying (branching)
the Assembly (see Branching_). 

.. code-block:: python

data1                      ## [0] name the Assembly data1
clust90_minsamp4           ## [0] name the Assembly based on parameter settings


1. Project dir
----------------
The Project directory is the location where a group of Assemblies which share
data files will be saved. This can be either a name or a path. If it is a path
then the a new directory will be created at the given path if it does not already
exist. If it is a name then a new directory with that name will be created in the
current directory if it does not already exist. A good name for Project_dir will
generally be the name of the organism being studied (e.g., white_crowned_sparrow). 

.. code-block:: python

/home/deren/ipyrad/tests/finches   ## [1] create/use project dir called finches
finches                            ## [1] create/use project dir called finches


2. Location of non-demultiplexed data
--------------------------------------
This is a path to the location of raw (non-demultiplexed) fastq data files. The 
files can be gzip compressed (i.e., have name-endings .fastq.gz). If you enter
a path for raw data files then you should also have a path to a barcodes file
for parameter 3. To select multiple files, or all files in a directory, use a
wildcard character (*), as in the examples:

.. code-block:: python

/home/deren/ipyrad/tests/data/*.fastq.gz  ## [2] select all gzip data files
~/ipyrad/tests/data/*.fastq.gz            ## [2] select all gzip data files
./data/sim_rad*.fastq.gz                  ## [2] select `sim_rad` data files


3. Location of barcodes file
-----------------------------
This is a path to the location of the barcodes_file_. 



