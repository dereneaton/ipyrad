.. include:: global.rst

.. _outline:

Assembly Outline
================

The typical workflow to move from fastq formatted sequence data to assembled 
comparative data sets in ipyrad_ involves 7 sequential steps. Separating the
assembly into steps allows for a modular workflow where jobs can easily be 
restarted from any step, and in which multiple final assembled data sets can 
be created by applying different parameter settings to it at different stages
of the assembly. 


Seven Steps
------------

1. Demultiplexing / Loading fastq files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 1 involves loading sequence data into a named "Assembly" and sorting the 
sequences among a number of "Samples" (individuals). If the data are not yet 
demultiplexed then step 1 uses information from a barcodes_file_ to assign 
sequences to Samples. If the data are already demultiplexed then step 1 simply
reads the data in to count how many reads are assigned to each Sample. 

2. Filtering / Editing reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...

3. Clustering / Mapping reads within Samples and alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...

4. Joint estimation of heterozygosity and error rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...

5. Consensus base calling and filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...

6. Clustering / Mapping reads among Samples and alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...

7. Filtering and formatting output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
...


Schematic Example
------------------

Basic workflow
^^^^^^^^^^^^^^^

**Example CLI basic workflow**

.. code-block:: bash

    ## create a params.txt file and rename it data1, and then use a text editor
    ## to edit the parameter settings in data1-params.txt
    ipyrad -n data1

    ## run steps 1-7 using the settings in data1-params.txt
    ipyrad -p data1-params.txt


**Example API basic workflow**

.. code-block:: python

    ## import ipyrad 
    import ipyrad as ip

    ## create an Assembly and modify some parameter settings
    data1 = ip.Assembly("data1")
    data1.set_params("project_dir", "example")
    data1.set_params("raw_fastq_path", "data/*.fastq")
    data1.set_params("barcodes_path", "barcodes.txt")   

    ## run steps 1-7
    data1.run()


### INSERT IMAGE HERE




Branching workflow
-------------------


**Example CLI workflow**


.. code-block:: bash

    ## create a params.txt file and rename it data1, and then use a text editor
    ## to edit the parameter settings in data1-params.txt
    ipyrad -n data1

    ## run steps 1-2 using the default settings
    ipyrad -p data1-params.txt -s 12

    ## branch to create a 'copy' of this assembly named data2
    ipyrad -p data1-params.txt -b data2

    ## edit data2-params.txt to a different parameter settings in a text editor,
    ## for example, change the clustering threshold from 0.85 to 0.90

    ## now run the remaining steps (3-7) on each data set
    ipyrad -p data1-params.txt -s 34567
    ipyrad -p data2-params.txt -s 34567


**Example API workflow**


.. code-block:: python

    ## import ipyrad 
    import ipyrad as ip

    ## create an Assembly and modify some parameter settings
    data1 = ip.Assembly("data1")
    data1.set_params("project_dir", "example")
    data1.set_params("raw_fastq_path", "data/*.fastq")
    data1.set_params("barcodes_path", "barcodes.txt")   

    ## run steps 1-2
    data1.run("12")

    ## create a new branch of this Assembly named data2
    ## and change some parameter settings 
    data2 = data1.copy("data2")
    data2.set_params("clust_threshold", 0.90)

    ## run steps 3-7 for the two Assemblies
    data1.run("34567")
    data2.run("34567")



### INSERT IMAGE HERE


