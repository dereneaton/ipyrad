
.. include:: global.rst 

.. _tutorial_advanced_cli:


Advanced tutorial -- CLI
========================
This is the advanced tutorial for the command line interface to ipyrad. 
In this tutorial we will introduce two new methods that were not
used in the introductory tutorial, but which provide some exciting new 
functionality. The first is ``branching``, which is used to 
efficiently assemble multiple data sets under a range of parameter settings, 
and the second is ``reference mapping``, which is a way to leverage information
from reference genomic data (e.g., full genome, transcriptome, 
plastome, etc) during assembly. 


Branching Assemblies
~~~~~~~~~~~~~~~~~~~~
If you've already been through the introductory tutorial you'll remember that 
a typical ipyrad analysis runs through seven sequential steps to take data 
from its raw state to finished output files of aligned data. 
After finishing one assembly, it is common that we might want to create a 
second assembly of our data under a different set of parameters; 
say by changing the ``clust_threshold`` from 0.85 to 0.90, or changing 
``min_samples_locus`` from 4 to 20. 

If we were to restart our analysis from the very beginning that would be really 
inefficient. So one way to go about this would be to change a few parameters in 
the params file to try to re-run an existing assembly by re-using some of the
existing data files. This approach is a little tricky, since the user would need
to know which files to rename/move, and it has the problem that previous results
files and parameters could be overwritten so that you lose information about how
the first data set was assembled. Simplifying this process is the motivation 
behind the branching assembly process in ipyrad, which does all of this renaming 
business for you, allowing efficient re-use of existing data files, while also 
keeping separate records (params files) of which parameters were used in each assembly. 

At its core, branching creates a copy of an Assembly object (the object that is
saved as a ``.json`` file by ipyrad) such that the new Assembly inherits all of 
the information from it's parent Assembly, including filenames, samplenames, 
and assembly statistics. The branching process requires a 
new :ref:`assembly_name<assembly_name>`, which is important so that all new files
created along this branch will be saved with a unique filename prefix. 
We'll show an example of a branching process below, but first we need to 
describe reference mapping, since for our example we will be creating two 
branches which are assembled using different ``assembly_methods``. 


Reference Sequence Mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~
ipyrad_ offers four :ref:`assembly methods<assembly_methods>`, three of which 
can utilize a reference sequence file. The first method, called ``reference``, 
maps RAD sequences to a reference file to determine homology and discards all
sequences which do not match to it. The second method, ``denovo+reference``, 
uses the reference first to identify homology, but then the remaining unmatched
sequences are all dumped into the standard ``denovo`` ipyrad pipeline to be clustered.
In essence, the reference file is simply used to assist the denovo assembly, and 
to add additional information. The final method, ``denovo-reference``, 
removes any reads which match to the reference and retains only non-matching 
sequences to be used in a denovo analysis. In other words, it allows the use 
of a reference sequence file as a filter to remove reads which match to it. You
can imagine how this would be useful for removing contaminants, plastome data, 
symbiont-host data, or coding/non-coding regions.  


Preview Mode
~~~~~~~~~~~~~
...


Running ipyrad CLI on a cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As explained in the :ref:`installation<installation>` section, ipyrad_ is very 
easy to use on an HPC cluster because as long as it is installed using conda_ it
does not require the user to load any external modules or software. 
Really, there is only **one** extra argument that you need to remember
to use which is the ``--MPI`` argument. This ensures that processing cores which
are split across different nodes of a cluster can all see the same data. Using 
ipyrad_ with the --MPI flag on an HPC machine should allow users to split jobs
across dozens or hundreds of cores to assemble data sets very rapidly. As an 
example, a large phylogenetic-scale RAD-seq data set analyzed on my desktop 
computer with 10 cores took ~2 days, while on an HPC system with 64 cores it 
took only ~16 hours. More detailed speed comparisons are in the works. 



Getting started
~~~~~~~~~~~~~~~
Let's first download the example simulated data sets for ipyrad_. Copy and paste
the code below into a terminal. This will create a new directory called 
``ipsimdata/`` in your current directory containing all of the necessary files.

.. code:: bash

    ## The curl command needs a capital O, not a zero.
    curl -LkO https://github.com/dereneaton/ipyrad/blob/master/tests/ipsimdata.tar.gz
    tar -xvzf ipsimdata.tar.gz


If you look in the ``ipsimdata/`` directory you'll see there are a number of example
data sets. For this tutorial we'll be using one called ``sim_rad_test``. Let's 
start by creating a new Assembly, and then we'll edit the params file to 
tell it how to find the input data files for this data set.

.. code:: bash

    ## creates a new Assembly named data1
    ipyrad -n data1


.. parsed-literal::

    New file params-data1.txt created in /home/deren/Documents/ipyrad


As you can see, this created a new  params file for our Assembly. We need to 
edit this file since it contains only default values. Use any text editor to 
open the params file ``params-data1.txt`` and enter the values 
below for parameters 1, 2, and 3. All other parameters can be left at their 
default values for now. This tells ipyrad that we are going to use the name
``iptutorial`` as our project_dir (where output files will be created), and 
that the input data and barcodes file are located in ``ipsimdata/``.

.. parsed-literal::

    ## ./iptutorial                              ## [1] [project_dir] ...
    ## ./ipsimdata/sim_rad_test_R1_.fastq.gz     ## [2] [raw_fastq_path] ...
    ## ./ipsimdata/sim_rad_test_barcodes.txt     ## [3] [barcodes_path] ...


Now we're ready to start the assembly. Let's begin by running just steps 1 and 2
to demultiplex and filter the sequence data. This will create a bunch of new 
files in the ``iptutorial/`` directory. 

.. code:: bash 

    ipyrad -p params-data1.txt -s 12


.. parsed-literal::

    Step1: Demultiplexing fastq data to Samples.
      Saving Assembly.
    Step2: Filtering reads 
      Saving Assembly.


Inside ``iptutorial`` you'll see that ipyrad_ has created two subdirectories 
with names prefixed by the assembly_name ``data1``. The other saved file is a 
``.json`` file,  which you can look at with a text editor if you wish. 
It's used by ipyrad_ to store information about your Assembly. 
You'll notice that ipyrad_ prints "Saving Assembly" quite often. 
This allows the assembly to be restarted easily from any point 
if it ever interrupted. In general, you should not mess with the .json file, 
since editing it by hand could cause errors in your assembly. 

.. code:: bash
    ls ./iptutorial

.. parsed-literal::
    data1_edits/   data1_fastqs/   data1.json


Branching example
~~~~~~~~~~~~~~~~~
For this example we will branch our Assembly before running step3 so that we can
see the results when the data are asembled with different assembly_methods. Our
existing assembly ``iptest1`` is using the denovo method. Let's create a branch
called ``iptest2`` which will use reference assembly. First we need to run the 
branch command, then we'll edit the new params file to change the assembly_method
and add the reference sequence file. 


.. code:: bash

    ## create a new branch of the Assembly iptest1
    ipyrad -p params-iptest1.txt -b iptest2
    
.. parsed-literal::

    New file params-iptest2.txt created in /home/deren/Documents/ipyrad


And make the following edits to ``params-iptest2.txt``:

.. parsed-literal::

    ## reference                               ## [5] [assembly_method] ...
    ## ./ipsimdata/sim_mt_genome.fa            ## [6] [reference_sequence] ...


Now we can run steps 3-7 on these two assemblies each using their own params 
file and each will create its own output files and saved results. It's worth 
noting that we could also have finished the first assembly all the way through 
and then gone back and made a new branch of it later. The only difference in that
case is that the second assembly would have saved file locations for the outputs
from steps 3-7, and so when you try to run it it will give a warning that the 
assembly has already finished steps 3-7. You can override this warning by passing
the ``-f`` flag, or ``--force``, which tells ipyrad that you think you know what
you're doing and want it to go ahead. 

.. code:: bash
   
   ## assemble the first data set denovo
   ipyrad -p params-iptest1.txt -s 34567

   ## assemble the second data set using reference mapping
   ipyrad -p params-iptest2.txt -s 34567


.. parsed-literal::

    --------------------------------------------------
     ipyrad [v.0.1.47]
     Interactive assembly and analysis of RADseq data
    --------------------------------------------------
     loading Assembly: iptest1 [/home/deren/Documents/ipyrad/iptest1.json]
     ipyparallel setup: Local connection to 4 Engines

     Step3: Clustering/Mapping reads
       Saving Assembly.
     Step4: Joint estimation of error rate and heterozygosity
       Saving Assembly.   
     Step5: Consensus base calling
       Diploid base calls and paralog filter (max haplos = 2)
       error rate (mean, std):  0.00075, 0.00002
       heterozyg. (mean, std):  0.00196, 0.00018
       Saving Assembly.
     Step6: Clustering across 12 samples at 0.85 similarity
       Saving Assembly.
     Step7: Filtering and creating output files 
       Saving Assembly.

     loading Assembly: iptest2 [/home/deren/Documents/ipyrad/iptest2.json]

     Step3: Clustering/Mapping reads
       Saving Assembly.
     Step4: Joint estimation of error rate and heterozygosity
       Saving Assembly.   
     Step5: Consensus base calling
       Diploid base calls and paralog filter (max haplos = 2)
       error rate (mean, std):  0.00075, 0.00002
       heterozyg. (mean, std):  0.00196, 0.00018
       Saving Assembly.
     Step6: Clustering across 12 samples at 0.85 similarity
       Saving Assembly.
     Step7: Filtering and creating output files 
       Saving Assembly.






old api tutorial
~~~~~~~~~~~~~~

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


