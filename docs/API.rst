


.. include:: global.rst


.. _API:


The ipyrad API
===============
The API_ (application program interface) for ipyrad_ is a way of directly 
accessing the nuts and bolts of ipyrad_ using Python_. 
This has a number of advantages over the CLI_ in that there is a lot more 
flexibility for creating highly complex branching assemblies, or for applying
ipyrad_ in a non-standard way. It's best feature, though, is that you can 
perform entire analyses within Jupyter :ref:`notebooks<notebooks>` to create 
documented reproducible code for your analyses. 


Why use the API? 
----------------
The API provides a more flexible framework for writing code to perform
complex branching assemblies than the CLI, and it also provides a more 
method of running remote code over very large computing clusters. 
Because it is interactive you can easily access the results and statistics 
from each step and use these in downstream analyses to visualize, analyze, 
and compare assemblies. 


Two main functions of the API
------------------------------

* Assembly -- perform all the assembly steps available in the CLI.
* Analysis -- analyze and compare the size and distribution of data sets, 
create plots, calculate population genetic statistics, and perform phylogenetic
analyses. 

Most users will start by using the CLI, but may eventually find that the
API provides functionality that is necessary for advanced usage. For example, 
the API can be used to merge two Assemblies so that Samples from different 
Assemblies can be clustered together. This can be useful when combining data 
from a previous study with newly collected data. 


Getting started with IPython/Jupyter notebooks
-----------------------------------------------
Our goal with using the ipyrad API is not only to get people writing Python 
scripts, but also to encourage the use of the a really exciting new tool called
Jupyter notebooks...
There may be a slight learning curve, however, for users who have no prior
experience with Python scripting. 

The envisioned usage of the ipyrad Python API is to run test assemblies within
a Jupyter notebook on a local computer using the **preview mode** method to 
execute quickly. Once the 



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



.. _cookbook_recipes:

Cookbook recipes - API 
-----------------------


.. toctree::
   :maxdepth: 1

   ipyrad_scripts.rst
   pedicularis.rst
   HPC_script.rst






old api stuff
-----------------

All of the parameter settings are linked to an Assembly object, which
has a set of default parameters when it is created. These can be viewed
using the ``get_params()`` function. To get more detailed information
about all parameters use ``ip.get_params_info()`` or to select a single
parameter use ``ip.get_params_info(3)``. Assembly objects have a
function ``set_params()`` that can be used to modify parameters.


If the data are already demultiplexed then fastq files can be linked
directly to the Data object, which in turn will create Sample objects
for each fastq file (or pair of fastq files for paired data). The files
may be gzip compressed. If the data are not demultiplexed then you will
have to run the step1 function below to demultiplex the raw data.

If for some reason we wanted to execute on just a subsample of our data,
we could do this by selecting only certain samples to call the ``step2``
function on. Because ``step2`` is a function of ``data``, it will always
execute with the parameters that are linked to ``data``.

.. code:: python

    ## example of ways to run step 2 to filter and trim reads
    #data1.step2("1B_0")                 ## run on a single sample
    #data1.step2(["1B_0", "1C_0"])       ## run on one or more samples
    data1.step2(force=True)              ## run on all samples, skipping finished ones
    
    ## print the results
    print data1.stats.head()


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


.. code:: python

    data1.stats.to_csv("data1_results.csv", sep="\t")
    data1.stats.to_latex("data1_results.tex")


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


.. code:: python

    import ipyrad as ip
    data1 = ip.load_assembly("test_rad/data1")

.. code:: python

    ## run step 4
    data1.step4("1A_0", force=True)
    
    ## print the results
    print data1.stats


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

.. code:: python

    ip.get_params_info(10)


