
User guide to the *ipyrad* API
==============================

Welcome! This tutorial will introduce you to the basic and advanced
features of working with the *ipyrad* API to assemble RADseq data in
Python. The API offers many advantages over the command-line interface,
but requires a little more work up front to learn the necessary tools
for using it. This includes knowing some very rudimentary Python, and
setting up a Jupyter notebook.

.. code:: ipython2

    import ipyrad as ip

Getting started with Jupyter notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial is an example of a `Jupyter
Notebook <http://jupyter.org/Jupyter>`__. If you've installed *ipyrad*
then you already have jupyter installed as well, which you can start
from the command-line (type ``jupyter-notebook``) to launch an
interactive notebook like this one. For some background on how jupyter
notebooks work I would recommend searching on google, or watching this
`YouTube video <https://www.youtube.com/watch?v=HW29067qVWk&t=47s>`__.
Once you have the hang of it, follow along with this code in your own
notebook.

Connecting your notebook to a cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have a separate tutorial about using Jupyter notebooks and connecting
Jupyter notebooks to a computing cluster (see
`here <http://ipyrad.readthedocs.io/HPC_Tunnel.html#>`__). For this
notebook I will assume that you are running this code in a Jupyter
notebook, and that you have an *ipcluster* instance running either
locally or remotely on a cluster. If an *ipcluster* instance is running
on its default settings (default profile) then *ipyrad* will
automatically use all available cores started on that cluster instance.

Import Python libraries
~~~~~~~~~~~~~~~~~~~~~~~

The only library we need to import is *ipyrad*. The *import* command is
usually the first code called in a Python document to load any necessary
packages. In the code below, we use a convenient trick in Python to tell
it that we want to refer to *ipyrad* simply as *ip*. This saves us a
little space since we might type the name many times. Below that, we use
the print statement to print the version number of *ipyrad*. This is
good practice to keep a record of which software version we are using.

.. code:: ipython2

    ## this is a comment, it is not executed, but the code below it is.
    import ipyrad as ip
    
    ## here we print the version
    print ip.__version__


.. parsed-literal::

    0.7.8


The *ipyrad* API data structures
--------------------------------

There are two main objects in *ipyrad*: Assembly class objects and
Sample class objects. And in fact, most users will only ever interact
with the Assembly class objects, since Sample objects are stored inside
of the Assembly objects, and the Assembly objects have functions, such
as merge, and branch, that are designed for manipulating and exchanging
Samples between different Assemblies.

Assembly Class objects
~~~~~~~~~~~~~~~~~~~~~~

Assembly objects are a unique data structure that ipyrad uses to store
and organize information about how to Assemble RAD-seq data. It contains
functions that can be applied to data, such as clustering, and aligning
sequences. And it stores information about which settings (prarmeters)
to use for assembly functions, and which Samples the functions should be
applied to. You can think of it mostly as a container that has a set of
rules associated with it.

To create a new Assembly object use the ``ip.Assembly()`` function and
pass it the name of your new Assembly. Creating an object in this way
has exactly the same effect as using the **-n {name}** argument in the
*ipyrad* command line tool, except in the API instead of creating a
params.txt file, we store the new Assembly information in a Python
variable. This can be named anything you want. Below I name the variable
*data1* so it is easy to remember that the Assembly name is also data1.

.. code:: ipython2

    ## create an Assembly object named data1. 
    data1 = ip.Assembly("data1")



.. parsed-literal::

    New Assembly: data1


Setting parameters
~~~~~~~~~~~~~~~~~~

You now have a Assembly object with a default set of parameters
associated with it, analogous to the params file in the command line
tool. You can view and modify these parameters using two arguments to
the Assembly object, ``set_params()`` and ``get_params()``.

.. code:: ipython2

    ## setting/modifying parameters for this Assembly object
    data1.set_params('project_dir', "pedicularis")
    data1.set_params('sorted_fastq_path', "./example_empirical_rad/*.gz")
    data1.set_params('filter_adapters', 2)
    data1.set_params('datatype', 'rad')
    
    ## prints the parameters to the screen
    data1.get_params()


.. parsed-literal::

    0   assembly_name               data1                                        
    1   project_dir                 ./pedicularis                                
    2   raw_fastq_path                                                           
    3   barcodes_path                                                            
    4   sorted_fastq_path           ./example_empirical_rad/*.gz                 
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
    16  filter_adapters             2                                            
    17  filter_min_trim_len         35                                           
    18  max_alleles_consens         2                                            
    19  max_Ns_consens              (5, 5)                                       
    20  max_Hs_consens              (8, 8)                                       
    21  min_samples_locus           4                                            
    22  max_SNPs_locus              (20, 20)                                     
    23  max_Indels_locus            (8, 8)                                       
    24  max_shared_Hs_locus         0.5                                          
    25  trim_reads                  (0, 0, 0, 0)                                 
    26  trim_loci                   (0, 0, 0, 0)                                 
    27  output_formats              ['p', 's', 'v']                              
    28  pop_assign_file                                                          


Instantaneous parameter (and error) checking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A nice feature of the ``set_params()`` function in the *ipyrad* API is
that it checks your parameter settings at the time that you change them
to make sure that they are compatible. By contrast, the *ipyrad* CLI
does not check params until you try to run a step function. Below you
can see that an error is raised when we try to set the
"clust\_threshold" parameters with an integer, since it requires the
value to be a float (decimal). It's hard to catch every possible error,
but we've tried to catch many of the most common errors in parameter
settings.

.. code:: ipython2

    ## this is expected to raise an error, since the clust_threshold 
    ## parameter cannot be 2.0
    data1.set_params("clust_threshold", 2.0)


::


    ---------------------------------------------------------------------------

    IPyradError                               Traceback (most recent call last)

    <ipython-input-5-feb2fc2a340d> in <module>()
          1 ## this is expected to raise an error, since the clust_threshold
          2 ## parameter cannot be 2.0
    ----> 3 data1.set_params("clust_threshold", 2.0)
    

    /home/deren/Documents/ipyrad/ipyrad/core/assembly.pyc in set_params(self, param, newvalue)
        797         except Exception as inst:
        798             raise IPyradWarningExit(BAD_PARAMETER\
    --> 799                                     .format(param, inst, newvalue))
        800 
        801 


    /home/deren/Documents/ipyrad/ipyrad/assemble/util.pyc in __init__(self, *args, **kwargs)
         50     def __init__(self, *args, **kwargs):
         51         if ipyrad.__interactive__:
    ---> 52             raise IPyradError(*args, **kwargs)
         53         else:
         54             SystemExit.__init__(self, *args, **kwargs)


    IPyradError:     Error setting parameter 'clust_threshold'
        clust_threshold must be a decimal value between 0 and 1.
        You entered: 2.0
        


Attributes of Assembly objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assembly objects have many attributes which you can access to learn more
about your Assembly. To see the full list of options you can type the
name of your Assembly variable, followed by a '.', and then press . This
will use tab-completion to list all of the available options. Below I
print a few examples.

.. code:: ipython2

    print data1.name


.. parsed-literal::

    data1


.. code:: ipython2

    ## another example attribute listing directories
    ## associated with this object. Most are empty b/c
    ## we haven't started creating files yet. But you 
    ## can see that it shows the fastq directory. 
    print data1.dirs


.. parsed-literal::

    fastqs : 
    edits : 
    clusts : 
    consens : 
    outfiles : 
    


Sample Class objects
~~~~~~~~~~~~~~~~~~~~

Sample Class objects correspond to individual samples in your study.
They store the file paths pointing to the data that is saved on disk,
and they store statistics about the results of each step of the
Assembly. Sample class objects are stored inside Assembly class objects,
and can be added, removed, or merged with other Sample class objects
between differnt Assemblies.

Creating Samples
~~~~~~~~~~~~~~~~

Samples are created during step 1 of the ipyrad Assembly. This involves
either demultiplexing raw data files or loading data files that are
already demultiplexed. For this example we are loading demultiplexed
data files. Because we've already entered the path to our data files in
``sorted_fastq_path`` of our Asssembly object, we can go ahead and run
step 1 to create Sample objects that are linked to the data files.

.. code:: ipython2

    ## run step 1 to create Samples objects
    data1.run("1")



.. parsed-literal::

    Assembly: data1
    [####################] 100%  loading reads         | 0:00:11 | s1 | 


The ``.run()`` command
~~~~~~~~~~~~~~~~~~~~~~

The run function is equivalent to the ``-s`` argument in the ipyrad
command line tool, and tell ipyrad which steps (1-7) of the assembly to
run. If a step has already been run on your samples they will be skipped
and it will print a warning. You can enforce overwriting the existing
data using the force flag.

.. code:: ipython2

    ## The force flag allows you to re-run a step that is already finished
    data1.run("1", force=True)


.. parsed-literal::

    Assembly: data1
    [####################] 100%  loading reads         | 0:00:11 | s1 | 


The run command will automatically parallelize work across all cores of
a running ipcluster instance (remember, you should have started this
outside of notebook. Or you can start it now.) If ipcluster is running
on the default profile then ipyrad will detect and use it when the run
command is called. However, if you start an ipcluster instance with a
specific profile name then you will need to connect to it using the
ipyparallel library and then pass the connection client object to
ipyrad. I'll show an example of that here.

.. code:: ipython2

    ## this is the explicit way to connect to ipcluster
    import ipyparallel
    
    ## connect to a running ipcluster instance
    ipyclient = ipyparallel.Client()
    
    ## if you used a named profile then enter that
    ipyclient = ipyparallel.Client(profile="default")

.. code:: ipython2

    ## call the run function of ipyrad and pass it the ipyclient
    ## process that you want the work distributed on.
    data1.run("1", ipyclient=ipyclient, force=True)


.. parsed-literal::

    Assembly: data1
    [####################] 100%  loading reads         | 0:00:10 | s1 | 


Samples stored in an Assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can see below that after step 1 has been run there will be a
collection of Sample objects stored in an Assembly that can be accessed
from the attribute ``.Samples``. They are stored as a dictionary in
which the keys are Sample names and the values of the dictionary are the
Sample objects.

.. code:: ipython2

    ## Sample objects stored as a dictionary
    data1.samples




.. parsed-literal::

    {'29154_superba': <ipyrad.core.sample.Sample at 0x7f98b1c5fc50>,
     '30556_thamno': <ipyrad.core.sample.Sample at 0x7f98b1c5fd10>,
     '30686_cyathophylla': <ipyrad.core.sample.Sample at 0x7f98b1c48650>,
     '32082_przewalskii': <ipyrad.core.sample.Sample at 0x7f98b1c39fd0>,
     '33413_thamno': <ipyrad.core.sample.Sample at 0x7f98b1c39dd0>,
     '33588_przewalskii': <ipyrad.core.sample.Sample at 0x7f98b1c22090>,
     '35236_rex': <ipyrad.core.sample.Sample at 0x7f98b1c22290>,
     '35855_rex': <ipyrad.core.sample.Sample at 0x7f98b1c48b50>,
     '38362_rex': <ipyrad.core.sample.Sample at 0x7f98b1c5f850>,
     '39618_rex': <ipyrad.core.sample.Sample at 0x7f98b1c93090>,
     '40578_rex': <ipyrad.core.sample.Sample at 0x7f98b1c37750>,
     '41478_cyathophylloides': <ipyrad.core.sample.Sample at 0x7f98b1bfff10>,
     '41954_cyathophylloides': <ipyrad.core.sample.Sample at 0x7f98b1c6abd0>}



The progress bar
~~~~~~~~~~~~~~~~

As you can see running a step of the analysis prints a progress bar
similar to what you would see in the *ipyrad* command line tool. There
are some differences, however. It shows on the far right "s1" to
indicate that this was step 1 of the assembly, and it does not print
information about our cluster setup (e.g., number of nodes and cores).
This was a stylistic choice to provide a cleaner output for analyses
inside Jupyter notebooks. You can view the cluster information when
running the step functions by adding the argument ``show_cluster=True``.

.. code:: ipython2

    ## run step 1 to create Samples objects
    data1.run("1", show_cluster=True, force=True)



.. parsed-literal::

    host compute node: [4 cores] on oud
    Assembly: data1
    [####################] 100%  loading reads         | 0:00:09 | s1 | 


Viewing results of Assembly steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Results for each step are stored in Sample class objects, however,
Assembly class objects have functions available for summarizing the
stats of all Sample class objects that they contain, which provides an
easy way to view results. This includes ``.stats`` attribute, and the
``.stats_dfs`` attributes for each step.

.. code:: ipython2

    ## print full stats summary
    print data1.stats


.. parsed-literal::

                            state  reads_raw
    29154_superba               1     696994
    30556_thamno                1    1452316
    30686_cyathophylla          1    1253109
    32082_przewalskii           1     964244
    33413_thamno                1     636625
    33588_przewalskii           1    1002923
    35236_rex                   1    1803858
    35855_rex                   1    1409843
    38362_rex                   1    1391175
    39618_rex                   1     822263
    40578_rex                   1    1707942
    41478_cyathophylloides      1    2199740
    41954_cyathophylloides      1    2199613


.. code:: ipython2

    ## print full stats for step 1 (in this case it's the same but for other
    ## steps the stats_dfs often contains more information.)
    print data1.stats_dfs.s1


.. parsed-literal::

                            reads_raw
    29154_superba              696994
    30556_thamno              1452316
    30686_cyathophylla        1253109
    32082_przewalskii          964244
    33413_thamno               636625
    33588_przewalskii         1002923
    35236_rex                 1803858
    35855_rex                 1409843
    38362_rex                 1391175
    39618_rex                  822263
    40578_rex                 1707942
    41478_cyathophylloides    2199740
    41954_cyathophylloides    2199613


Branching to subsample taxa
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Branching in the *ipyrad* API works the same as in the CLI, but in many
ways is easier to use because you can access attributes of the Assembly
objects much more easily, such as when you want to provide a list of
Sample names in order to subsample (exclude samples) during the
branching process. Below is an example.

.. code:: ipython2

    ## access all Sample names in data1
    allsamples = data1.samples.keys()
    print "Samples in data1:\n", "\n".join(subsamples)


.. parsed-literal::

    Samples in data1:
    30686_cyathophylla
    33413_thamno
    30556_thamno
    32082_przewalskii
    29154_superba
    41478_cyathophylloides
    40578_rex
    35855_rex
    33588_przewalskii
    39618_rex
    38362_rex
    35236_rex
    41954_cyathophylloides


.. code:: ipython2

    ## Drop the two samples from this list that have "prz" in their names.
    ## This is an easy programmatic way to remove the outgroup samples
    ## in this data set.
    subs = [i for i in allsamples if "prz" not in i]
    
    ## use branching to create new Assembly named 'data2'
    ## with only Samples whose name is in the subs list
    data2 = data1.branch("data2", subsamples=subs)
    print "Samples in data2:\n", "\n".join(data2.samples)


.. parsed-literal::

    Samples in data2:
    30686_cyathophylla
    33413_thamno
    41478_cyathophylloides
    29154_superba
    40578_rex
    35855_rex
    30556_thamno
    39618_rex
    38362_rex
    35236_rex
    41954_cyathophylloides


Branching to iterate over parameter settings
--------------------------------------------

This is the real bread and butter of the *ipyrad* API.

You can write simple for-loops using Python code to apply a range of
parameter settings to different branched assemblies. Furthermore, using
branching this can be done in a way that greatly reduces the amount of
computation needed to produce multiple data sets. Essentially, branching
allows you to recycle intermediate states that are shared between
branched Assemblies. This is particularly useful when assemblies differ
by only one or few parameters that are applied late in the assembly
process. To set up efficient branching code in this way requires some
prior knowledge about when (which step) each parameter is applied in
ipyrad. That information is available in the documentation
(http://ipyrad.readthedocs.io/parameters.html).

When setting up for-loop routines like the one below it may be helpful
to break the script up among multiple cells of a Jupyter notebook so
that you can easily restart from one step or another. It may also be
useful to subsample your data set to a small number of samples to test
the code first, and if all goes well, then proceed with your full data
set.

An example to create many assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example below we will create 8 complete Assemblies which vary in
three different parameter combinations (filter\_setting,
clust\_threshold, and min\_sample).

.. code:: ipython2

    ## Start by creating an initial assembly, setting the path to your data, 
    ## and running step1. I set a project-dir so that all of our data sets
    ## will be grouped into a single directory called 'branch-test'.
    data = ip.Assembly("base")
    data.set_params("project_dir", "branch-test")
    data.set_params("raw_fastq_path", "./ipsimdata/rad_example_R1_.fastq.gz")
    data.set_params("barcodes_path", "./ipsimdata/rad_example_barcodes.txt")
    
    ## step 1: load in the data
    data.run('1')


.. parsed-literal::

    New Assembly: base
    Assembly: base
    [####################] 100%  sorting reads         | 0:00:02 | s1 | 
    [####################] 100%  writing/compressing   | 0:00:00 | s1 | 


.. code:: ipython2

    ## let's create a dictionary to hold the finished assemblies
    adict = {}
    
    ## iterate over parameters settings creating a new named assembly
    for filter_setting in [1, 2]:
        ## create a new name for the assembly and branch
        newname = data.name + "_f{}".format(filter_setting)
        child1 = data.branch(newname)
        child1.set_params("filter_adapters", filter_setting)
        child1.run("2")
        
        ## iterate over clust thresholds
        for clust_threshold in ['0.85', '0.90']:
            newname = child1.name + "_c{}".format(clust_threshold[2:])
            child2 = child1.branch(newname)
            child2.set_params("clust_threshold", clust_threshold)
            child2.run("3456")
            
            ## iterate over min_sample coverage
            for min_samples_locus in [4, 12]:
                newname = child2.name + "_m{}".format(min_samples_locus)
                child3 = child2.branch(newname)
                child3.set_params("min_samples_locus", min_samples_locus)
                child3.run("7")
                
                ## store the complete assembly in the dictionary by its name
                ## so it is easy for us to access and retrieve, since we wrote
                ## over the variable name 'child' during the loop. You can do
                ## this using dictionaries, lists, etc., or, as you'll see below,
                ## we can use the 'load_json()' command to load a finished assembly
                ## from its saved file object.
                adict[newname] = child3


.. parsed-literal::

    Assembly: base_f1
    [####################] 100%  processing reads      | 0:00:03 | s2 | 
    Assembly: base_f1_c85
    [####################] 100%  dereplicating         | 0:00:00 | s3 | 
    [####################] 100%  clustering            | 0:00:00 | s3 | 
    [####################] 100%  building clusters     | 0:00:00 | s3 | 
    [####################] 100%  chunking              | 0:00:00 | s3 | 
    [####################] 100%  aligning              | 0:00:11 | s3 | 
    [####################] 100%  concatenating         | 0:00:00 | s3 | 
    [####################] 100%  inferring [H, E]      | 0:00:03 | s4 | 
    [####################] 100%  calculating depths    | 0:00:00 | s5 | 
    [####################] 100%  chunking clusters     | 0:00:00 | s5 | 
    [####################] 100%  consens calling       | 0:00:14 | s5 | 
    [####################] 100%  concat/shuffle input  | 0:00:00 | s6 | 
    [####################] 100%  clustering across     | 0:00:01 | s6 | 
    [####################] 100%  building clusters     | 0:00:00 | s6 | 
    [####################] 100%  aligning clusters     | 0:00:04 | s6 | 
    [####################] 100%  database indels       | 0:00:00 | s6 | 
    [####################] 100%  indexing clusters     | 0:00:01 | s6 | 
    [####################] 100%  building database     | 0:00:00 | s6 | 
    Assembly: base_f1_c85_m4
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:00 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f1_c85_m4_outfiles
    
    Assembly: base_f1_c85_m12
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:01 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f1_c85_m12_outfiles
    
    Assembly: base_f1_c90
    [####################] 100%  dereplicating         | 0:00:00 | s3 | 
    [####################] 100%  clustering            | 0:00:00 | s3 | 
    [####################] 100%  building clusters     | 0:00:00 | s3 | 
    [####################] 100%  chunking              | 0:00:00 | s3 | 
    [####################] 100%  aligning              | 0:00:11 | s3 | 
    [####################] 100%  concatenating         | 0:00:00 | s3 | 
    [####################] 100%  inferring [H, E]      | 0:00:03 | s4 | 
    [####################] 100%  calculating depths    | 0:00:00 | s5 | 
    [####################] 100%  chunking clusters     | 0:00:00 | s5 | 
    [####################] 100%  consens calling       | 0:00:14 | s5 | 
    [####################] 100%  concat/shuffle input  | 0:00:00 | s6 | 
    [####################] 100%  clustering across     | 0:00:01 | s6 | 
    [####################] 100%  building clusters     | 0:00:00 | s6 | 
    [####################] 100%  aligning clusters     | 0:00:04 | s6 | 
    [####################] 100%  database indels       | 0:00:00 | s6 | 
    [####################] 100%  indexing clusters     | 0:00:01 | s6 | 
    [####################] 100%  building database     | 0:00:00 | s6 | 
    Assembly: base_f1_c90_m4
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:00 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f1_c90_m4_outfiles
    
    Assembly: base_f1_c90_m12
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:00 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f1_c90_m12_outfiles
    
    Assembly: base_f2
    [####################] 100%  processing reads      | 0:00:04 | s2 | 
    Assembly: base_f2_c85
    [####################] 100%  dereplicating         | 0:00:00 | s3 | 
    [####################] 100%  clustering            | 0:00:00 | s3 | 
    [####################] 100%  building clusters     | 0:00:00 | s3 | 
    [####################] 100%  chunking              | 0:00:00 | s3 | 
    [####################] 100%  aligning              | 0:00:11 | s3 | 
    [####################] 100%  concatenating         | 0:00:00 | s3 | 
    [####################] 100%  inferring [H, E]      | 0:00:03 | s4 | 
    [####################] 100%  calculating depths    | 0:00:00 | s5 | 
    [####################] 100%  chunking clusters     | 0:00:00 | s5 | 
    [####################] 100%  consens calling       | 0:00:15 | s5 | 
    [####################] 100%  concat/shuffle input  | 0:00:00 | s6 | 
    [####################] 100%  clustering across     | 0:00:01 | s6 | 
    [####################] 100%  building clusters     | 0:00:00 | s6 | 
    [####################] 100%  aligning clusters     | 0:00:04 | s6 | 
    [####################] 100%  database indels       | 0:00:00 | s6 | 
    [####################] 100%  indexing clusters     | 0:00:01 | s6 | 
    [####################] 100%  building database     | 0:00:00 | s6 | 
    Assembly: base_f2_c85_m4
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:01 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f2_c85_m4_outfiles
    
    Assembly: base_f2_c85_m12
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:01 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f2_c85_m12_outfiles
    
    Assembly: base_f2_c90
    [####################] 100%  dereplicating         | 0:00:00 | s3 | 
    [####################] 100%  clustering            | 0:00:00 | s3 | 
    [####################] 100%  building clusters     | 0:00:00 | s3 | 
    [####################] 100%  chunking              | 0:00:00 | s3 | 
    [####################] 100%  aligning              | 0:00:11 | s3 | 
    [####################] 100%  concatenating         | 0:00:00 | s3 | 
    [####################] 100%  inferring [H, E]      | 0:00:03 | s4 | 
    [####################] 100%  calculating depths    | 0:00:00 | s5 | 
    [####################] 100%  chunking clusters     | 0:00:00 | s5 | 
    [####################] 100%  consens calling       | 0:00:16 | s5 | 
    [####################] 100%  concat/shuffle input  | 0:00:00 | s6 | 
    [####################] 100%  clustering across     | 0:00:01 | s6 | 
    [####################] 100%  building clusters     | 0:00:00 | s6 | 
    [####################] 100%  aligning clusters     | 0:00:05 | s6 | 
    [####################] 100%  database indels       | 0:00:00 | s6 | 
    [####################] 100%  indexing clusters     | 0:00:01 | s6 | 
    [####################] 100%  building database     | 0:00:00 | s6 | 
    Assembly: base_f2_c90_m4
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:01 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f2_c90_m4_outfiles
    
    Assembly: base_f2_c90_m12
    [####################] 100%  filtering loci        | 0:00:00 | s7 | 
    [####################] 100%  building loci/stats   | 0:00:00 | s7 | 
    [####################] 100%  building vcf file     | 0:00:01 | s7 | 
    [####################] 100%  writing vcf file      | 0:00:00 | s7 | 
    [####################] 100%  building arrays       | 0:00:01 | s7 | 
    [####################] 100%  writing outfiles      | 0:00:00 | s7 | 
    Outfiles written to: ~/Documents/ipyrad/tests/branch-test/base_f2_c90_m12_outfiles
    


Saving Assembly objects
~~~~~~~~~~~~~~~~~~~~~~~

Assembly objects (and the Sample objects they contain) are automatically
saved each time that you use the ``.run()`` function. However, you can
also save by calling the ``.save()`` function of an Assembly object.
This updates the JSON file. Additionally, Assembly objects have a
function called ``.write_params()`` which can be invoked to create a
params file for use by the *ipyrad* command line tool.

.. code:: ipython2

    ## save assembly object (also auto-saves after every run() command)
    data1.save()
    
    ## load assembly object
    data1 = ip.load_json("pedicularis/data1.json")
    
    ## write params file for use by the CLI
    data1.write_params()
