


.. include:: global.rst


.. _API:


The ipyrad API
===============
The API_ (application program interface) for ipyrad_ is a way of directly 
accessing the nuts and bolts of ipyrad_ using Python_, and in particular
we have written it to work interactively in IPython_. 
This has a number of advantages over the CLI_ (command line interface) 
in that there is greater flexibility for creating highly 
complex branching assemblies, or for applying ipyrad_ in a 
non-standard way. Perhaps the best feature, though, is that you can 
perform entire analyses within Jupyter :ref:`notebooks<notebooks>` to create 
documented reproducible code for your analyses. 

Running ipyrad API in parallel
==============================
If you are using the API then you must have an ipcluster instance
started in order to parallelize your code. This can be started locally
by opening a separate terminal and running (``ipcluster start -n=10``)
to start 10 engines. Or, to run your code on a remote cluster set up
your ipcluster instance following `this tutorial <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__.

Cookbooks -- coming soon
----------------
We plan to add example cookbook recipes of IPython code that 
can used to perform non-standard procedures in ipyrad that
are not easily implemented in the CLI. 


Two main functions of the API
------------------------------

* Assembly -- perform all assembly steps available in the CLI but with several additional options.
* Analysis -- analyze and compare the size and distribution of data sets, create plots, calculate population genetic statistics, and perform phylogenetic
analyses. 


Getting started with IPython/Jupyter notebooks
-----------------------------------------------
Our goal in developing the ipyrad API is not only to get people writing Python 
scripts, but also to encourage the use of an exciting new tool called
Jupyter notebooks, an excellent tool for reproducible science. 

The envisioned usage of the ipyrad Python API is to run test assemblies within
a Jupyter notebook on a local computer.  
Once you've tested that your assembly looks good, and that 
your selected parameters seems appropriate you can then take your script and 
submit a long running job to a larger computing cluster. 

Beyond the very easy way that the CLI provides for connecting to HPC 
clusters, the API provides even more flexibility for doing very 
advanced computation for users that are fluent with the 
ipyparallel package. 


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


.. _cookbook_recipes_API:

Cookbook recipes - API 
-----------------------

.. toctree::
   :maxdepth: 1

   ipyrad_scripts.rst
   pedicularis.rst
   HPC_script.rst



old api stuff
-----------------
Begin by creating an Assembly class object, this is equivalent 
to running the (-n) argument from the command line. It initializes
an the Assembly object that will be used to store the assembly 
parameters that will be used. 

.. code:: python

   ## initalize an Assembly class object
   data = ip.Assemble("test")
   
The assembly object has a number of attributes and functions
that you can access interactively in IPython by using tab-completion. 
Simply type the name of your object (here named 'data'), followed
by tab to see all of the options. 

.. code:: python

    data.<tab>     
              data.barcodes         data.dirs             data.name              
              data.branch           data.files            data.outfiles          
              data.build_stat       data.get_params       data.paramsdict
              data.clust_database   data.link_fastqs      data.populations       
              data.database         data.link_populations data.run               


Two important functions of the Assembly include ``get_params()`` and
``set_params()``, which are used to view and modify the parameter 
settings, respectively. 

.. code:: python

    data.get_params()


.. parsed-literal::

  0   assembly_name               test                                         
  1   project_dir                 ./                                           
  2   raw_fastq_path                                                           
  3   barcodes_path                                                            
  4   sorted_fastq_path                                                        
  5   assembly_method             denovo                                       
  6   reference_sequence                                                       
  7   datatype                    rad                                          
  8   restriction_overhang        ('TGCAG', '')                                
  9   max_low_qual_bases          5                                            
  10  phred_Qscore_offset         33                                           
  11  mindepth_statistical        6                                            
  12  mindepth_majrule            6                                            
  13  maxdepth                    10000                                        
  14  clust_threshold             0.85                                         
  15  max_barcode_mismatch        0                                            
  16  filter_adapters             0                                            
  17  filter_min_trim_len         35                                           
  18  max_alleles_consens         2                                            
  19  max_Ns_consens              (5, 5)                                       
  20  max_Hs_consens              (8, 8)                                       
  21  min_samples_locus           4                                            
  22  max_SNPs_locus              (20, 20)                                     
  23  max_Indels_locus            (8, 8)                                       
  24  max_shared_Hs_locus         0.5                                          
  25  edit_cutsites               (0, 0)                                       
  26  trim_overhang               (4, 4, 4, 4)                                 
  27  output_formats              *                                            
  28  pop_assign_file                       


You can change parameters by providing either the name or index number 
of the parameter you want to change and the new parameter setting
to ``set_params``. Below we also show what happens if you enter an invalid
parameter.


.. code:: python

    ## change a few parameters
    data.set_params("project_dir", "iptest")
    data.set_params("raw_fastq_path", "ipsimdata/rad_example_R1_.fastq.gz")
    data.set_params("barcodes_path", "ipsimdata/rad_example_barcodes.txt")

    ## if you enter a parameter setting that is invalid an error will raise
    data.set_params("max_alleles_consens", "four")


.. parsed-literal::

    IPyradError:     
        Error setting parameter 'max_alleles_consens'
        invalid literal for int() with base 10: 'four'
        You entered: four


The next important function is the ``run()`` command, which is used
to run steps of the assembly. 

.. code:: python 

    ## run a few steps of the assembly
    data.run("12")


Let's imagine at this point that we are interested in clustering our
data at two different clustering thresholds. We will try 0.90 and 0.85.
First we need to make a copy the Assembly object. This will inherit the
locations of the data linked in the first object, but diverge in any
future applications to the object. Thus, they can share the same working
directory, and will inherit shared files, but create divergently linked
files within this directory. You can view the directories linked to an
Assembly object with the ``.dirs`` argument, shown below. 


.. code:: python

    ## run step 3 to cluster reads within samples
    data.run(3) 

    ## print the results
    print data1.stats

    ## create a new branch
    data2 = data.branch("c85")
    data2.set_params("clust_threshold", 0.85)

    ## run steps on new assembly
    data2.run(3, force=True)



And you can see below that the two Assembly objects are now working with
several shared directories (working, fastq, edits) but with different
clust directories (clust\_0.85 and clust\_0.9).

.. code:: python

    print "data1 directories:"
    for (i,j) in data.dirs.items():
        print "{}\t{}".format(i, j)
        
    print "\ndata2 directories:"
    for (i,j) in data2.dirs.items():
        print "{}\t{}".format(i, j)

