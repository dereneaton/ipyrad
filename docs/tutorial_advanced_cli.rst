

.. _tutorial_advanced_cli:


Advanced tutorial -- CLI
========================
This is the advanced tutorial for the command line interface to ipyrad. 
In this tutorial we will introduce two new methods that were not
used in the introductory tutorial, but which provide a lot of exciting new 
functionality to ipyrad. The first is ``branching``, which is used to 
efficiently assemble multiple data sets under a range of parameter settings, 
and the second is ``reference mapping``, which is a way to leverage information
from any kind of reference genomic data (e.g., full genome, transcriptome, 
plastome, etc). 



Branching Assemblies
~~~~~~~~~~~~~~~~~~~~
If you've already been through the introductory tutorial you'll remember that 
a typical ipyrad analysis runs through seven sequential steps to take data 
from its raw state into finished output files with aligned data. 
After finishing one assembly, it is common that we might want to create a 
second assembly of our data under a different set of parameters; 
say by changing the ``clust_threshold`` from 0.85 to 0.90, or changing 
``min_samples_locus`` from 4 to 20. 

If we were to restart our analysis from the very beginning that would be really 
inefficient. Thus, we could just change a few parameters in the params file 
and re-run the existing assembly to achieve new results, but that's not quite 
ideal either as it would require the user to rename/move the results files, 
to avoid overwriting existing files, and they would lose a record of the parameter
setting that the first assembly was created under. This was the motivation for the 
branching assembly method which is aimed at allowing efficient re-use of existing
data files, while also keeping good records of which parameters were used for 
different resulting assemblies. 

At its core, branching creates a copy of an Assembly object (the object that is
saved as a ``.json`` file by ipyrad) such that a new Assembly inherits all of 
the information from it's parent Assembly, including filenames, samplenames, 
and assembly statistics. The branching process always requires a 
new :ref:`assembly_name<assembly_name>`, which is important in allowing the 
new branch to re-use old files and statistics, but create new files separate from
its parent branch going forward which are saved with a the new assembly_name prefix.


Reference Sequence Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~
ipyrad_ offers three :ref:`assembly_methods<assembly_methods>` which can utilize
a reference sequence file. The first method, ``reference``, maps all RAD sequences
to the reference file in order to determine homology among sequences. The second
method, ``denovo+reference``, uses the reference first to identify homology 
based on anything that matches to the reference, and then the remaining reads are
all dumped into the standard ``denovo`` ipyrad pipeline and clustered together. 
The final method, ``denovo-reference``, removes any reads which match to the 
reference and only retains non-matching sequences to be used in a denovo analysis. 




Running ipyrad CLI on a cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation<installation>` section, ipyrad_ is very 
easy to use on a cluster because as long as it is installed using conda_ it does
not require the user to load any external modules or software to run on an HPC
cluster. Really, there is only **one** extra argument that you need to remember
to use which is the ``--MPI`` argument. This ensures that processing cores which
are split across different nodes of a cluster can all see the same data. Using 
ipyrad_ with the --MPI flag on an HPC machine should allow users 



Getting started
~~~~~~~~~~~~~~~
First let's grab the :ref:`simulated data sets from here<simulated_data>` to 
use for this analysis. We'll use the data set called ``sim_rad_test``. 

.. code:: bash

    ## create a new Assembly under the name data1. This will 
    ## create a new params file called params-data1.txt
    ipyrad -n data1

Then use a text editor to make the following changes to ``params-data1.txt``. 
All other parameters can be left at their default values.

.. parsed-literal::

    ## ./advanced_tutorial                  ## [1] [project_dir]
    ## ./simdata/sim_rad_test.fastq.gz      ## [2] [raw_fastq_path]
    ## ./simdata/sim_rad_test_barcodes.txt  ## [3] [barcodes_path]

Let's then run steps 1-2 to demultiplex the data and run it through the filtering
step. 

.. code:: bash
    ipyrad -p params-data1.txt -s 12

Now look in your working directory. You'll notice that ipyrad_ created a new 
directory here called ``advanced_tutorial``. This is our project directory. If 
we had set a full path for ``project_dir`` in the params file (i.e., one like 
``/home/user/myfolder/`` instead of ``./myfolder``) the new directory would be 
created in the full path location. If the project_dir already exists that is fine.
You'll see that in the project_dir ipyrad_ has created two new directories which 
have names prefixed by the name of our assembly, ``data1``. 

.. code:: bash
    ls ./advanced_tutorial

.. parsed-literal::
    data1_edits   data1_fastqs   data1.json

The other saved file is a json file which you can look at with a text editor if 
you wish. It's used by ipyrad_ to store information about your Assembly. You 
should generally not edit this file by hand. 


Branching example
~~~~~~~~~~~~~~~~~~
...





Assembly and Sample objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assembly and Sample objects are used by *ipyrad* to access data stored
on disk and to manipulate it. Each biological sample in a data set is
represented in a Sample object, and a set of Samples is stored inside an
Assembly object. The Assembly object has functions to assemble the data,
and stores a log of all steps performed and the resulting statistics of
those steps. Assembly objects can be copied or merged to allow branching
events where different parameters can subsequently be applied to
different Assemblies going forward. Examples of this are shown below.

To create an Assembly object call ``ip.Assembly()`` and pass a name for
the data set. An Assembly object does not initially contain Samples,
they will be created either by linking fastq files to the Assembly
object if data are already demultiplexed, or by running ``step1()`` to
demultiplex raw data files, as shown below.

.. code:: python

    ## create an Assembly object called data1. 
    data1 = ip.Assembly("data1")
    
    ## The object will be saved to disk using its assigned name
    print "Assembly object named", data1.name


.. parsed-literal::

    Assembly object named data1


Modifying assembly parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the parameter settings are linked to an Assembly object, which
has a set of default parameters when it is created. These can be viewed
using the ``get_params()`` function. To get more detailed information
about all parameters use ``ip.get_params_info()`` or to select a single
parameter use ``ip.get_params_info(3)``. Assembly objects have a
function ``set_params()`` that can be used to modify parameters.

.. code:: python

    ## modify parameters for this Assembly object
    data1.set_params(1, "./test_rad")
    data1.set_params(2, "./data/sim_rad_test_R1_.fastq.gz")
    data1.set_params(3, "./data/sim_rad_test_barcodes.txt")
    #data1.set_params(2, "~/Dropbox/UO_C353_1.fastq.part-aa.gz")
    #data1.set_params(3, "/home/deren/Dropbox/Viburnum_revised.barcodes")
    data1.set_params(7, 3)
    data1.set_params(10, 'rad')
    
    ## print the new parameters to screen
    data1.get_params()


.. parsed-literal::

      1   project_dir                   ./test_rad                                   
      2   raw_fastq_path                ./data/sim_rad_test_R1_.fastq.gz             
      3   barcodes_path                 ./data/sim_rad_test_barcodes.txt             
      4   sorted_fastq_path                                                          
      5   restriction_overhang          ('TGCAG', '')                                
      6   max_low_qual_bases            5                                            
      7   N_processors                  3                                            
      8   mindepth_statistical          6                                            
      9   mindepth_majrule              6                                            
      10  datatype                      rad                                          
      11  clust_threshold               0.85                                         
      12  minsamp                       4                                            
      13  max_shared_heterozygosity     0.25                                         
      14  prefix_outname                data1                                        
      15  phred_Qscore_offset           33                                           
      16  max_barcode_mismatch          1                                            
      17  filter_adapters               0                                            
      18  filter_min_trim_len           35                                           
      19  ploidy                        2                                            
      20  max_stack_size                1000                                         
      21  max_Ns_consens                5                                            
      22  max_Hs_consens                8                                            
      23  max_SNPs_locus                (100, 100)                                   
      24  max_Indels_locus              (5, 99)                                      
      25  trim_overhang                 (1, 2, 2, 1)                                 
      26  hierarchical_clustering       0                                            


Starting data
~~~~~~~~~~~~~

If the data are already demultiplexed then fastq files can be linked
directly to the Data object, which in turn will create Sample objects
for each fastq file (or pair of fastq files for paired data). The files
may be gzip compressed. If the data are not demultiplexed then you will
have to run the step1 function below to demultiplex the raw data.

.. code:: python

    ## This would link fastq files from the 'sorted_fastq_path' if present
    ## Here it does nothing b/c there are no files in the sorted_fastq_path
    data1.link_fastqs()


.. parsed-literal::

    0 new Samples created in data1.
    0 fastq files linked to Samples.


Step 1: Demultiplex the raw data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This uses the barcodes information to demultiplex reads in data files
found in the 'raw\_fastq\_path'. It will create a Sample object for each
sample that will be stored in the Assembly object.

.. code:: python

    ## run step 1 to demultiplex the data
    data1.step1()
    
    ## print the results for each Sample in data1
    print data1.stats.head()


.. parsed-literal::

          state  reads_raw  reads_filtered  clusters_total  clusters_kept  
    1A_0      1      20099             NaN             NaN            NaN   
    1B_0      1      19977             NaN             NaN            NaN   
    1C_0      1      20114             NaN             NaN            NaN   
    1D_0      1      19895             NaN             NaN            NaN   
    2E_0      1      19928             NaN             NaN            NaN   
    
          hetero_est  error_est  reads_consens  
    1A_0         NaN        NaN            NaN  
    1B_0         NaN        NaN            NaN  
    1C_0         NaN        NaN            NaN  
    1D_0         NaN        NaN            NaN  
    2E_0         NaN        NaN            NaN  


Step 2: Filter reads
~~~~~~~~~~~~~~~~~~~~

If for some reason we wanted to execute on just a subsample of our data,
we could do this by selecting only certain samples to call the ``step2``
function on. Because ``step2`` is a function of ``data``, it will always
execute with the parameters that are linked to ``data``.

.. code:: python

    %%time
    ## example of ways to run step 2 to filter and trim reads
    #data1.step2("1B_0")                 ## run on a single sample
    #data1.step2(["1B_0", "1C_0"])       ## run on one or more samples
    data1.step2(force=True)              ## run on all samples, skipping finished ones
    
    ## print the results
    print data1.stats.head()

Step 3: clustering within-samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's imagine at this point that we are interested in clustering our
data at two different clustering thresholds. We will try 0.90 and 0.85.
First we need to make a copy the Assembly object. This will inherit the
locations of the data linked in the first object, but diverge in any
future applications to the object. Thus, they can share the same working
directory, and will inherit shared files, but create divergently linked
files within this directory. You can view the directories linked to an
Assembly object with the ``.dirs`` argument, shown below. The
prefix\_outname (param 14) of the new object is automatically set to the
Assembly object name.

.. code:: python

    ## run step 3 to cluster reads within samples using vsearch
    #data1.step3(['2E_0'], force=True, preview=True)  # ["2H_0", "2G_0"])
    data1.step3(force=True)
    ## print the results
    print data1.stats.head()

Branching Assembly objects
~~~~~~~~~~~~~~~~~~~~~~~~~~

And you can see below that the two Assembly objects are now working with
several shared directories (working, fastq, edits) but with different
clust directories (clust\_0.85 and clust\_0.9).

.. code:: python

    ## create a branch of our Assembly object
    data2 = data1.branch(newname="data2")
    
    ## set clustering threshold to 0.90
    data2.set_params(11, 0.90)
    
    ## look at inherited parameters
    data2.get_params()

.. code:: python

    ## run step 3 to cluster reads within samples using vsearch
    data2.step3(force=True)  # ["2H_0", "2G_0"])
    
    ## print the results
    print data2.stats

.. code:: python

    print "data1 directories:"
    for (i,j) in data1.dirs.items():
        print "{}\t{}".format(i, j)
        
    print "\ndata2 directories:"
    for (i,j) in data2.dirs.items():
        print "{}\t{}".format(i, j)

.. code:: python

    ## TODO, just make a [name]_stats directory in [work] for each data obj
    data1.statsfiles


Saving stats outputs
~~~~~~~~~~~~~~~~~~~~

.. code:: python

    data1.stats.to_csv("data1_results.csv", sep="\t")
    data1.stats.to_latex("data1_results.tex")

Example of plotting with *ipyrad*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a a few simple plotting functions in *ipyrad* useful for
visualizing results. These are in the module ``ipyrad.plotting``. Below
is an interactive plot for visualizing the distributions of coverages
across the 12 samples in the test data set.

.. code:: python

    import ipyrad.plotting as iplot
    
    ## plot for one or more selected samples
    iplot.depthplot(data1, ["1A_0", "1B_0"])
    
    ## plot for all samples in data1
    #iplot.depthplot(data1)
    
    ## save plot as pdf and html
    iplot.depthplot(data1, outprefix="testfig")

Step 4: Joint estimation of heterozygosity and error rate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import ipyrad as ip
    data1 = ip.load_assembly("test_rad/data1")

.. code:: python

    ## run step 4
    data1.step4("1A_0", force=True)
    
    ## print the results
    print data1.stats

Step 5: Consensus base calls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #import ipyrad as ip
    
    ## reload autosaved data. In case you quit and came back 
    #data1 = ip.load_dataobj("test_rad/data1.assembly")

.. code:: python

    ## run step 5
    data1.step5()
    
    ## print the results
    print data1.stats

.. code:: python

    data1.samples["1A_0"].stats

Quick parameter explanations are always on-hand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ip.get_params_info(10)

Log history
~~~~~~~~~~~

A common problem after struggling through an analysis is that you find
you've completely forgotten what parameters you used at what point, and
when you changed them. The log history time stamps all calls to
``set_params()``, as well as calls to ``step`` methods. It also records
copies/branching of data objects.

.. code:: python

    for i in data1.log:
        print i

Saving Assembly objects
~~~~~~~~~~~~~~~~~~~~~~~

Assembly objects can be saved and loaded so that interactive analyses
can be started, stopped, and returned to quite easily. The format of
these saved files is a serialized 'dill' object used by Python.
Individual Sample objects are saved within Assembly objects. These
objects to not contain the actual sequence data, but only link to it,
and so are not very large. The information contained includes parameters
and the log of Assembly objects, and the statistics and state of Sample
objects. Assembly objects are autosaved each time an assembly ``step``
function is called, but you can also create your own checkpoints with
the ``save`` command.

.. code:: python

    ## save assembly object
    #ip.save_assembly("data1.p")
    
    ## load assembly object
    #data = ip.load_assembly("data1.p")
    #print data.name
