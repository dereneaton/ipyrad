
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

It would be wholy inefficient to restart from the beginning for each assembly 
that uses different parameter settings. A better way would be to re-use existing
data files and only rerun steps downstream from where parameter changes have
an effect. This approach is a little tricky, since the user would need
to know which files to rename/move to avoid existing results files and parameter
information from being overwritten and lost. 

The motivation behind the branching assembly process in ipyrad is to simplify 
this process. ipyrad does all of this renaming business for you, and creates
new named files in a way the retains records of the existing assemblies and 
effectively re-uses existing data files. 

At its core, branching creates a copy of an Assembly object (the object that is
saved as a ``.json`` file by ipyrad) such that the new Assembly inherits all of 
the information from it's parent Assembly, including filenames, samplenames, 
and assembly statistics. The branching process requires a 
new :ref:`assembly_name<assembly_name>`, which is important so that all new files
created along this branch will be saved with a unique filename prefix. 
We'll show an example of a branching process below, but first we need to 
describe reference mapping, since for our example we will be creating two 
branches which are assembled using different 
:ref:`assembly methods<assembly_methods>`. 


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
computer with 12 cores took ~2 days, while on an HPC system with access to
48 cores it took only ~12 hours. More detailed speed comparisons are in the 
works. For most steps of ipyrad the speed improvement is linear with 
the number of cores. 



Getting started
~~~~~~~~~~~~~~~
Let's first download the example simulated data sets for ipyrad_. Copy and paste
the code below into a terminal. This will create a new directory called 
``ipsimdata/`` in your current directory containing all of the necessary files.

.. code:: bash

    ## The curl command needs a capital O, not a zero.
    curl -LkO https://github.com/dereneaton/ipyrad/raw/master/tests/ipsimdata.tar.gz
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

    ## enter these lines into the params-data1.txt file
    ./iptutorial                              ## [1] [project_dir] ...
    ./ipsimdata/sim_rad_test_R1_.fastq.gz     ## [2] [raw_fastq_path] ...
    ./ipsimdata/sim_rad_test_barcodes.txt     ## [3] [barcodes_path] ...


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
existing assembly ``data1`` is using the denovo method. Let's create a branch
called ``data2`` which will use reference assembly. First we need to run the 
branch command, then we'll edit the new params file to change the assembly_method
and add the reference sequence file. 


.. code:: bash

    ## create a new branch of the Assembly iptest1
    ipyrad -p params-data1.txt -b data2
    
.. parsed-literal::

    New file params-data2.txt created in /home/deren/Documents/ipyrad


And make the following edits to ``params-data2.txt``:

.. parsed-literal::

    ## reference                               ## [5] [assembly_method] ...
    ## ./ipsimdata/sim_mt_genome.fa            ## [6] [reference_sequence] ...


Now we can run steps 3-7 on these two assemblies each using their own params 
file and each will create its own output files and saved results. 

.. code:: bash
   
   ## assemble the first data set denovo
   ipyrad -p params-data1.txt -s 34567

   ## assemble the second data set using reference mapping
   ipyrad -p params-data2.txt -s 34567


.. parsed-literal::

    --------------------------------------------------
     ipyrad [v.0.1.73]
     Interactive assembly and analysis of RADseq data
    --------------------------------------------------
     loading Assembly: data1 [/home/deren/Documents/ipyrad/data1.json]
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

     loading Assembly: data2 [/home/deren/Documents/ipyrad/data2.json]

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


Now let's suppose we're interested in the effect of missing data on our assemblies
and we want to assemble each data set with a different ``min_samples_locus`` 
setting. Maybe at 4, 8, and 12 (ignore the fact that the example data set 
has no missing data, and so this has no practical effect). It's worth 
noting that we can branch assemblies after an analysis has finished as well. 
The only difference is that the new assembly will think that it has already 
finished all of the steps, and so if we ask it to run them again it will instead
want to skip over them. You can override this behavior by passing the ``-f`` flag, 
or ``--force``, which tells ipyrad that you want it to run the step even though
it's already finished it. The two assemblies we finished were both assembled at
the default value of 4 for ``min_samples_locus``, so below I set up code to 
branch and then run step7 on each of these assemblies with a new setting of 8 or 12. 

.. code:: bash
   
   ## branch data1 to make min8 and min12 data sets
   ipyrad -p params-data1.txt -b data1-min8
   ipyrad -p params-data1.txt -b data1-min12

   ## use a text editor to set min_samples_locus to the new value (4 or 8) in each

   ## branch data2 to make min8 and min12 data sets
   ipyrad -p params-data2.txt -b data2-min8
   ipyrad -p params-data2.txt -b data2-min12

   ## use a text editor to set min_samples_locus to the new value in each

   ## run step7 on using the new min_samples_locus settings
   ipyrad -p params-data1-min8.txt -s 7
   ipyrad -p params-data1-min12.txt -s 7
   ipyrad -p params-data2-min8.txt -s 7
   ipyrad -p params-data2-min12.txt -s 7


Now if we look in our project_dir ``iptutorial/`` we see that the fastq/ 
and edits/ directories were created using just the first assembly ``data1``, 
while the clust/ and consens/ directories were created for both ``data1`` and
``data2``, since both completed steps 3-6. Finally, you can see that each 
assembly has its own ``outfiles/`` directory with the results of step7. 

.. code:: bash

   ## use ls -l to view inside the project directory as a list
   ls -l iptutorial/

I show the file tree structure a bit more clearly below:

.. parsed-literal::  

   iptutorial/
       data1.json
       data1_fastqs/
       data1_edits/
       data1_clust_0.85/
       data1_consens/
       data1_outfiles/
       data1_min8.json
       data1_min8_outfiles/
       data1_min12.json 
       data1_min12_outfiles/
       data2.json
       data2_clust_0.85/
       data2_consens/
       data2_outfiles/
       data2_min8.json
       data2_min8_outfiles/
       data2.json
       data2_min12_outfiles/


In your working directory you will have the four params files which 
have the full set of parameters used in each of your assemblies. 
This makes for a good reproducible workflow, and can be referenced later
as a reminder of the parameters used for each data set. 


What's next
~~~~~~~~~~~
Check out the :ref:`example empirical data sets<pedicularis_cli>`
to see how this process looks
when you run it on a relatively quick set of real data. 


Writing ipyrad scripts
~~~~~~~~~~~~~~~~~~~~~~
From the code above you may have noticed that the only thing stopping you from
being able to write one long script that creates a whole range of assemblies is 
when you have to edit the new params files by hand. We've purposefully avoided 
creating an ipyrad command to change parameters on the fly, since this would 
make it so that the params file are not a good record of the parameter set used
throughout an entire assembly. 

However, if you're a very programmatic type of person who would prefer that 
all of your branching and parameter changing could take place within a single
script you'll want to check out the :ref:`ipyrad API<API>`, which provides a 
more elegant pure Python way to edit parameters in your code while 
assembling data. 


