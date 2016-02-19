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
of parameter settings is one of the strengths of ipyrad_. 

Below is an explanation of each parameter setting, the steps of the assembly 
that it effects, and example entries for the parameter into a params.txt file.

.. _assembly_name:

0. Assembly name
-----------------
The Assembly name is used as the prefix for all output files. It should be a 
unique identifier for the Assembly, meaning the set of parameters you are using 
for the current data set. When I assemble multiple data with different parameter
combinations I usually either name them consecutively (e.g., data1, data2), or 
with names indicating their parameters combinations (e.g., data_clust90, 
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
exist. If it is a name then a new directory with that name will be created in the
current directory if it does not already exist. A good name for Project_dir will
generally be the name of the organism being studied (e.g., white_crowned_sparrow). 

Affected steps: 1-7  

Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/finches   ## [1] create/use project dir called finches
    finches                            ## [1] create/use project dir called finches


.. _raw_fastq_path:

2. Raw fastq path
-----------------
This is a path to the location of raw (non-demultiplexed) fastq data files. The 
files can be gzip compressed (i.e., have name-endings .fastq.gz). If you enter
a path for raw data files then you should also have a path to a barcodes file
for parameter 3. To select multiple files, or all files in a directory, use a
wildcard character (*).

Affected steps: 1  

Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/*.fastq.gz  ## [2] select all gzip data files
    ~/ipyrad/tests/data/*.fastq.gz            ## [2] select all gzip data files
    ./data/sim_rad*.fastq.gz                  ## [2] select `sim_rad` data files

.. _barcodes_path:

3. Barcodes path
----------------
This is a path to the location of the barcodes_file_. 

Affected steps: 1   

Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/sim_barcodes.txt  ## [3] select barcode file
    ~/tests/data/sim_barcodes.txt                   ## [3] select barcode file


.. _sorted_fastq_path:

4. Sorted fastq path
--------------------
This is a path to the location of sorted fastq data. If your data are already
demultiplexed then this is the location that will be used in step 1 to load 
the data into ipyrad. 

Affected steps: 1  

Example entries into params.txt:  

.. code-block:: python

    /home/deren/ipyrad/tests/data/*.fastq.gz  ## [4] select all gzip data files
    ~/ipyrad/tests/data/*.fastq               ## [4] select all fastq data files
    ./data/sim_rad*.fastq.gz                  ## [4] select `sim_rad` data files

.. _assembly_method:

5. Assembly method
--------------------
There are four :ref:`Assembly_methods<Assembly_methods>` options in ipyrad_: 
denovo, reference, reference_add, and
reference_sub. The latter three all require a reference sequence file (param 6).

Affected steps: 3,6
Example entries into params.txt:  

.. code-block:: python

    denovo                            ## [5] denovo assembly
    reference                         ## [5] reference assembly
    reference_add                     ## [5] reference addition assembly
    reference_sub                     ## [5] reference subtraction assembly

.. _reference_sequence:

6. Reference sequence
---------------------
...

.. _datatype:

7. Datatype
------------
...


.. _restriction_overhang:

8. Restriction_overhang
-----------------------
...


.. _max_low_qual_bases:

9. max_low_qual_bases
---------------------
...

.. _phred_Qscore_offset:

10. Phred_Qscore_offset
------------------------
...

.. _mindepth_statistical:

11. Mindepth_statistical
-------------------------
...

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




