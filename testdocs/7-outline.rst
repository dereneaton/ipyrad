
.. include:: global.rst

.. _outline:


Assembly: Seven Steps
=====================
The goal of the *assembly* process is to convet raw or sorted fastq data into assembled loci that can be formatted for downstream analyses in phylogenetic or population genetic inference software. In *ipyrad* we have purposefully atomized this process into :ref:`seven sequential steps <seven_steps>` to create a modular workflow that can be easily restarted if interrupted, and can be :ref:`branched <branching_workflow>` at 
different points to create assemblies under different combinations of parameter settings. 


Basic Assembly Workflow
------------------------
The simplest use of ipyrad is to assemble a data set under a single set of parameters defined in a params file. Step 1 loads/assigns data to each sample; steps 2-5 process data for each sample; step 6 identifies orthologs across samples; and step 7 filters the orthologs and writes formatted files for downstream analyses. 

.. image:: ./images/steps.png

The code to run a basic workflow is quite simple:

.. code-block:: bash
    
    # create an initial Assembly params file
    >>> ipyrad -n data1 

    # enter values into the params file using a text editor
    ## ... editing params-data1.txt

    # select a params file (-p) and steps to run (-s) for this assembly
    >>> ipyrad -p params-data1.txt -s 1234567



.. include:: global.rst




Advanced Branching workflow
----------------------------
A more effective way to use *ipyrad* can be to create branching
assemblies in which multiple data sets are assembled under different parameter 
settings. The schematic below shows an example where an assembly 
is branched at step3. The new branch will inherit file paths and statistics 
from the first Assembly, but can then apply different parameters going forward.
Branching does not create hard copies of existing data files, and so is not 
an "expensive" action in terms of disk space or time. We suggest it be used 
quite liberally whenever applying a new set of parameters. 

.. image:: images/steps_branching.png


The code to run a branching workflow is only a bit more complex than the basic
workflow. You can find more branching examples in the 
:ref:`advanced tutorial<tutorial_advanced_cli>` and 
:ref:`cookbook<cookbook>` sections. 

.. code-block:: bash
    
    ## create an initial Assembly and params file, here called 'data1'
    >>> ipyrad -n data1 

    ## edit the params file for data1 with your text editor
    ## ... editing params-data1.txt

    ## run steps 1-2 with the params file
    >>> ipyrad -p params-data1.txt -s 12

    ## create a new branch of 'data1' before step3, here called 'data2'.
    >>> ipyrad -p params-data1.txt -b data2

    ## edit the params file for data2 using a text editor
    ## ... editing params-data2.txt

    ## run steps 3-7 for both assemblies
    >>> ipyrad -p params-data1.txt -s 34567
    >>> ipyrad -p params-data2.txt -s 34567










.. _seven_steps:


Seven Steps
-------------


1. Demultiplexing / Loading fastq files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 1 loads sequence data into a named :ref:`Assembly<Assembly>` and assigns
reads to :ref:`Samples<Samples>` (individuals). If the data are not yet 
demultiplexed then step 1 uses information from a :ref:`barcodes file<barcodes_file>` 
to demultiplex the data, otherwise, it simply reads the data for each Sample. 

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step1:  

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`raw_fastq_path<raw_fastq_path>`  
* :ref:`barcodes_path<barcodes_path>`  
* :ref:`sorted_fastq_path<sorted_fastq_path>`  
* :ref:`*datatype<datatype>`  
* :ref:`restriction_overhang<restriction_overhang>`  
* :ref:`max_barcode_mismatch<max_barcode_mismatch>`  



2. Filtering / Editing reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 2 uses the quality score recorded in the fastQ data files to filter low 
quality base calls. Sites with a score below a set value are changed into “N”s, 
and reads with more than the number of allowed “N”s are discarded. The threshold
for inclusion is set with the :ref:`phred_Qscore_offset<phred_Qscore_offset>` 
parameter. An optional filter can be applied to remove adapters/primers
(see :ref:`filter_adapters<filter_adapters>`), and there is an 
optional filter to clean up the edges of poor quality reads
(see :ref:`edit_cutsites<edit_cutsites>`).

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step2: 

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`barcodes_path<barcodes_path>`  
* :ref:`*datatype<datatype>`  
* :ref:`restriction_overhang<restriction_overhang>`  
* :ref:`max_low_qual_bases<max_low_qual_bases>`  
* :ref:`filter_adapters<filter_adapters>`  
* :ref:`filter_min_trim_len<filter_min_trim_len>`  
* :ref:`edit_cut_sites<edit_cut_sites>`  


3. Clustering / Mapping reads within Samples and alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 3 first dereplicates the sequences from step 2, recording the number of 
times each unique read is observed. If the data are paired-end, it then uses
vsearch_ to merge paired reads which overlap. The resulting data are 
then either de novo clustered (using vsearch_) or mapped to a reference 
genome (using bwa_ and bedtools_), depending on the selected assembly method.
In either case, reads are matched together on the basis of sequence similarity
and the resulting clusters are aligned using muscle_. 

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step3: 

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`*assembly_method<assembly_method>`  
* :ref:`*datatype<datatype>`  
* :ref:`*clust_threshold<clust_threshold>`  


4. Joint estimation of heterozygosity and error rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step4 jointly estimates sequencing error rate and heterozygosity based on counts
of site patterns across clustered reads. 
These estimates are used in step5 for consensus base calling. If the 
max_alleles_consens is set to 1 (haploid) then heterozygosity is fixed to 0 and 
only error rate is estimated. For all other settings of max_alleles_consens 
a diploid model is used (i.e., two alleles are expected to occur equally). 

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step4: 

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`*datatype<datatype>`  
* :ref:`*restriction_overhang<restriction_overhang>`  
* :ref:`*max_alleles_consens<max_alleles_consens>`  


5. Consensus base calling and filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step5 estimates consensus allele sequences from clustered reads given the estimated
parameters from step 4 and a binomial model. During this step we filter for maximum
number of undetermined sites (Ns) per locus (max_Ns_consens). The number of alleles
at each locus is recorded, but a filter for max_alleles is not applied until step7. 
Read depth information is also stored at this step for the VCF output in step7. 

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step5:  

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`*datatype<datatype>`  
* :ref:`*max_Ns_consens<max_Ns_consens>`  


6. Clustering / Mapping reads among Samples and alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step6 clusters consensus sequences across Samples using the same assembly method 
as in step 3. One allele is randomly sampled before clustering so that ambiguous
characters have a lesser effect on clustering, but the resulting data retain
information for heterozygotes. The clustered sequences are then aligned using 
muscle_.

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step6: 

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`*datatype<datatype>`  


7. Filtering and formatting output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step7 applies filters to the final alignments and saves the final data in a 
number of possible :ref:`output formats<full_output_formats>`. This step is most 
often repeated at several different settings for the parameter 
:ref:`min_samples_locus` to create different assemblies with different 
proportions of missing data (see :ref:`branching_workflow`). 

The following :ref:`parameters<parameters>` are *potentially*
used or required (\*) for step7:   

* :ref:`*assembly_name<assembly_name>`  
* :ref:`*project_dir<project_dir>`  
* :ref:`*datatype<datatype>`  
* :ref:`min_samples_locus<min_samples_locus>`  
* :ref:`max_Indels_locus<max_Indels_locus>`  
* :ref:`max_shared_Hs_locus<max_shared_Hs_locus>`  
* :ref:`max_alleles_consens<max_alleles_consens>`  
* :ref:`trim_overhang<trim_overhang>`  
* :ref:`output_formats<output_formats>`  
* :ref:`pop_assign_file<pop_assign_file>`  





**Example CLI branching workflow**


.. code-block:: bash

    ## create a params.txt file and rename it data1, and then use a text editor
    ## to edit the parameter settings in data1-params.txt
    ipyrad -n data1

    ## run steps 1-2 using the default settings
    ipyrad -p params-data1.txt -s 12

    ## branch to create a 'copy' of this assembly named data2
    ipyrad -p params-data1.txt -b data2

    ## edit data2-params.txt to a different parameter settings in a text editor,
    ## for example, change the clustering threshold from 0.85 to 0.90

    ## now run the remaining steps (3-7) on each data set
    ipyrad -p params-data1.txt -s 34567
    ipyrad -p params-data2.txt -s 34567


**Example Python API branching workflow**


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
    data2 = data1.branch("data2")
    data2.set_params("clust_threshold", 0.90)

    ## run steps 3-7 for the two Assemblies
    data1.run("34567")
    data2.run("34567")
