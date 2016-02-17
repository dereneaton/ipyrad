
ipyrad command line tutorial
=========================

This is the full tutorial for the command line interface for ipyrad. In this
tutorial we'll walk through the entire assembly and analysis process. This is 
meant as a broad introduction to familiarize users with the general workflow,
and some of the parameters and terminology. For simplicity we'll use 
single-end RAD-Seq as the example data, but the core concepts will apply
to assembly of other data types (GBS and paired-end). 

If you are new to RADseq analyses, this tutorial will provide a simple overview 
of how to execute ipyrad, what the data files look like, and how to check that 
your analysis is working, and the expected output formats.

Each cell in this tutorial beginning with the header (%%bash) indicates that the 
code should be executed in a command line shell, for example by copying and 
pasting the text into your terminal (but excluding the %%bash header). All 
lines in code cells beginning with ## are comments and should not be copied
and executed.

Getting Started
~~~~~~~~~~~~~~~

If you haven't already installed ipyrad go here first: :ref:`Installation <installation>`

We provide a very small sample data set that we recommend using for this tutorial.
Full datasets can take days and days to run, whereas with the simulated data
you could complete the whole tutorial in an afternoon. 

First make a new directory and fetch & extract the test data.
.. code-block:: bash
    mkdir ipyrad-test
    cd ipyrad-test
    curl -O https://github.com/dereneaton/ipyrad/blob/master/tests/ipyrad_tutorial_data.tgz
    tar -xvzf ipyrad_tutorial_data.tgz

You should now see a folder in your current directory called `data`. This 
directory contains two files we'll be using:
    - sim_rad_test_R1_.fastq.gz - Illumina fastQ formatted reads (gzip compressed)
    - sim_rad_test_barcodes.txt - Mapping of barcodes to sample IDs

It also contains many other simulated datasets, as well as a simulated 
reference genome, so you can experiment with other datatypes after you get
comfortable with RAD.

Create a new parameters file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ipyrad uses a text file to hold all the parameters for a given assembly. 
Start by creating a new parameters file with the `-n` flag. This flag
requires you to pass in a name for your assembly. In the example we use 
`ipyrad-test` but the name can be anything at all. Once you start 
analysing your own data you might call your parameters file something 
more informative, like the name of your organism.

.. code:: bash
    ipyrad -n ipyrad-test

This will create a file in the current directory called `params-ipyrad-test.txt`.
The params file lists on each line one parameter followed by a ## mark, 
and then a short description of the parameter. Lets take a look at it.

.. code:: bash
    cat params-ipyrad-test.txt

.. parsed-literal::
    ------ ipyrad params file (v.0.1.47)--------------------------------------------
    ipyrad-test                    ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
    ./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
                                   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                                   ## [3] [barcodes_path]: Location of barcodes file
                                   ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
    denovo                         ## [5] [assembly_method]: Assembly method (denovo, hybrid, reference_only, denovo_only)
                                   ## [6] [reference_sequence]: Location of reference sequence file
    rad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
    TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
    5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
    33                             ## [10] [phred_Qscore_offset]: phred Q score offset (only alternative=64)
    6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
    6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
    1000                           ## [13] [maxdepth]: Max cluster depth within samples
    0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
    1                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
    0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
    35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
    2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
    5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
    8, 8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
    4                              ## [21] [min_samples_locus]: Min # samples per locus for output
    100, 100                       ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
    5, 99                          ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
    0.25                           ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
    0, 0                           ## [25] [edit_cutsites]: Edit cut-sites (R1, R2) (see docs)
    1, 2, 2, 1                     ## [26] [trim_overhang]: Trim overhang (see docs) (R1>, <R1, R2>, <R2)
    *                              ## [27] [output_formats]: Output formats (see docs)
                                   ## [28] [pop_assign_file]: Path to population assignment file
                                   ## [29] [excludes]: Samples to be excluded from final output files
                                   ## [30] [outgroups]: Outgroup individuals. Excluded from final output files

In general the defaults are sensible, and we won't mess with them for now, but there
are a few parameters we *must* change. We need to set the path to the raw data we 
want to analyse, and we need to set the path to the barcodes file.

In your favorite text editor open params-ipyrad-test.txt and change these two lines
to look like this, and then save it:

.. parsed-literal::
    ./data/sim_rad_test_R1_.fastq.gz         ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
    ./data/sim_rad_test_barcodes.txt         ## [3] [barcodes_path]: Location of barcodes file

Input data format
~~~~~~~~~~~~~~~~~
Before we get started let's take a look at what the raw data look like.

Your input data will be in fastQ format, usually ending in .fq, .fastq, or
.fq.gz, .fastq.gz. Your data could be split among multiple files, or all 
within a single file (de-multiplexing goes much faster if they happen to 
be split into multiple files). The file/s may be compressed with gzip so 
that they have a .gz ending, but they do not need to be. The location of 
these files should be entered on line 2 of the params file. Below are 
the first three reads in the example file.

.. code:: bash
    ## For your personal edification here is what this is doing:
    ##  gzip -c: Tells gzip to unzip the file and write the contents to the screen
    ##  head -n 12: Grabs the first 12 lines of the fastq file. Fastq files
    ##      have 4 lines per read, so the value of `-n` should be a multiple of 4
    ##  cut -c 1-90: Trim the length of each line to 90 characters
    ##      we don't really need to see the whole sequence we're just trying
    ##      to get an idea.

    gzip -c ./data/sim_rad_test_R1_.fastq.gz | head -n 12 | cut -c 1-90

And here's the output:

.. parsed-literal::
    @lane1_fakedata0_R1_0 1:N:0:
    TTTTAATGCAGTGAGTGGCCATGCAATATATATTTACGGGCGCATAGAGACCCTCAAGACTGCCAACCGGGTGAATCACTATTTGCTTAG
    +
    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    @lane1_fakedata0_R1_1 1:N:0:
    TTTTAATGCAGTGAGTGGCCATGCAATATATATTTACGGGCGCATAGAGACCCTCAAGACTGCCAACCGGGTGAATCACTATTTGCTTAG
    +
    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    @lane1_fakedata0_R1_2 1:N:0:
    TTTTAATGCAGTGAGTGGCCATGCAATATATATTTACGGGCGCATAGAGACCCTCAAGACTGCCAACCGGGTGAATCACTATTTGCTTAG
    +
    BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB


Each read takes four lines. The first is the name of the read (its 
location on the plate). The second line contains the sequence data. 
The third line is a spacer. And the fourth line the quality scores 
for the base calls. In this case arbitrarily high since the data 
were simulated.

These are 100 bp single-end reads prepared as RADseq. The first 
six bases form the barcode and the next five bases (TGCAG) the 
restriction site overhang. All following bases make up the sequence 
data.

Step 1: Demultiplex the raw data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 1 reads in the barcodes file and the raw data. It scans through
the raw data and sorts each read based on the mapping of samples to 
barcodes. At the end of this step we'll have a new directory in our project_dir
called `ipyrad-test_fastqs`. Inside this directory will be individual
fastq.gz files for each sample.

**NB:** You'll notice the name of this output directory bears a strong
resemblence to the name of the assembly we chose at the time
of the params file creation. Assembling rad-seq type sequence
data requires a lot of different steps, and these steps generate a 
_LOT_ of intermediary files. ipyrad organizes these files into 
directories, and it prepends the name of your assembly to each
directory with data that belongs to it. One result of this is that
you can have multiple assemblies of the same raw data with different
parameter settings and you don't have to manage all the files
yourself! (See :ref:`Branching assemblies <advanced_CLI>` for more
info). Another result is that *you should not rename or move any
of the directories inside your project directory*, unless you know
what you're doing or you don't mind if your assembly breaks. 

Lets take a look at the barcodes file for the simulated data. You'll 
see sample names (left) and their barcodes (right) each on a 
separate line with a tab between them.
.. code:: bash
    cat ./data/sim_rad_test_barcodes.txt

.. parsed-literal::
    1A_0    CATCAT
    1B_0    AGTGAT
    1C_0    ATGGTA
    1D_0    GTGGGA
    2E_0    AGGGAA
    2F_0    AAAGTG
    2G_0    GATATA
    2H_0    GAGGAG
    3I_0    GGGATT
    3J_0    TAATTA
    3K_0    TGAGGG
    3L_0    ATATTA

Now lets run step 1!

.. code:: bash
    ## -p indicates the params file we wish to use
    ## -s indicates the step to run
    ipyrad -p params-ipyrad-test.txt -s 1

.. parsed-literal::
 --------------------------------------------------
  ipyrad [v.0.1.47]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  New Assembly: ipyrad-test
  ipyparallel setup: Local connection to 4 Engines

  Step1: Demultiplexing fastq data to Samples.
    Saving Assembly.

There are 4 main parts to this step:
- Create a new assembly. Since this is our first time running
any steps we need to initialize our assembly.
- Start the parallel cluster. ipyrad uses a parallelization 
library called ipyparallel. Every time we start a step we 
fire up the parallel clients. This makes your assemblies go
**smokin'** fast.
- Actually do the demuliplexing.
- Save the state of the assembly.

Have a look at the results of this step in the `ipyrad-test_fastqs`
output directory:
.. code:: bash
   ls ipyrad-test_fastqs 

.. parsed-literal::
    1A_0_R1_.fastq.gz        1D_0_R1_.fastq.gz        2G_0_R1_.fastq.gz        3J_0_R1_.fastq.gz        s1_demultiplex_stats.txt
    1B_0_R1_.fastq.gz        2E_0_R1_.fastq.gz        2H_0_R1_.fastq.gz        3K_0_R1_.fastq.gz
    1C_0_R1_.fastq.gz        2F_0_R1_.fastq.gz        3I_0_R1_.fastq.gz        3L_0_R1_.fastq.gz

A more informative metric of success might be the number
of raw reads demultiplexed for each sample. Fortunately 
ipyrad tracks the state of all your steps in your current 
assembly, so at any time you can ask for results by 
invoking the `-r` flag.

.. code:: bash
    ## -r fetches informative results from currently 
    ##      executed steps
    ipyrad -p params-ipyrad-test.txt -r

.. parsed-literal::
    Summary stats of Assembly ipyrad-test
    ------------------------------------------------
          reads_raw  state
    1A_0      20099      1
    1B_0      19977      1
    1C_0      20114      1
    1D_0      19895      1
    2E_0      19928      1
    2F_0      19934      1
    2G_0      20026      1
    2H_0      19936      1
    3I_0      20084      1
    3J_0      20011      1
    3K_0      20117      1
    3L_0      19901      1

If you want to get even **more** info ipyrad tracks all kinds of
wacky stats and saves them to a file inside the directories it
creates for each step. For instance to see full stats for step 1:

.. code:: bash                                                                                                                                       
    cat ./ipyrad-test_fastqs/s1_demultiplex_stats.txt

And you'll see a ton of fun stuff I won't copy here in the interest
of conserving space. Please go look for yourself if you're interested.

Step 2: Filter reads
~~~~~~~~~~~~~~~~~~~~
This step filters reads based on quality scores, and can be used to 
detect Illumina adapters in your reads, which is sometimes a problem 
with homebrew type library preparations. Here the filter is set to the 
default value of 0, meaning it filters only based on quality scores of 
base calls. The filtered files are written to a new directory called 
`ipyrad-test_edits`.

.. code:: bash
    ipyrad -p params-ipyrad-test.txt -s 2

.. parsed-literal::
 --------------------------------------------------
  ipyrad [v.0.1.47]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  loading Assembly: ipyrad-test [/private/tmp/ipyrad-test/ipyrad-test.json]
  ipyparallel setup: Local connection to 4 Engines

  Step2: Filtering reads 
    Saving Assembly.

Again, you can look at the results output by this step and also some 
handy stats tracked for this assembly.
.. code:: bash
    ## View the output of step 2
    ls ipyrad-test_edits

.. parsed-literal::                                                                                                                                  
    1A_0_R1_.fastq       1C_0_R1_.fastq       2E_0_R1_.fastq       2G_0_R1_.fastq       3I_0_R1_.fastq       3K_0_R1_.fastq       s2_rawedit_stats.txt
    1B_0_R1_.fastq       1D_0_R1_.fastq       2F_0_R1_.fastq       2H_0_R1_.fastq       3J_0_R1_.fastq       3L_0_R1_.fastq

.. code:: bash
    ## Get current stats including # raw reads and # reads
    ## after filtering.
    ipyrad -p params-ipyrad-test.txt -r

.. parsed-literal::
    Summary stats of Assembly ipyrad-test
    ------------------------------------------------
          reads_filtered  reads_raw  state
    1A_0           20099      20099      2
    1B_0           19977      19977      2
    1C_0           20114      20114      2
    1D_0           19895      19895      2
    2E_0           19928      19928      2
    2F_0           19934      19934      2
    2G_0           20026      20026      2
    2H_0           19936      19936      2
    3I_0           20084      20084      2
    3J_0           20011      20011      2
    3K_0           20117      20117      2
    3L_0           19901      19901      2

You might also take a gander at the filtered reads:
.. code:: bash
    head -n 12 ./ipyrad-test_fastqs/1A_0_R1_.fastq


Step 3: clustering within-samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 3 de-replicates and then clusters reads within each sample 
by the set clustering threshold and writes the clusters to new 
files in a directory called `ipyrad-test_clust_0.85`.

The true name of this output directory will be dictated by the value
you set for the `clust_threshold` parameter in the params file. 
.. parsed-literal::
    0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
You can see the default value is 0.85, so our default directory is 
named accordingly. This value dictates the percentage of sequence
similarity that reads must have in order to be considered reads
at the same locus. You'll more than likely want to experiment
with this value, but 0.85 is a reliable default, balancing
over-splitting of loci vs over-lumping. Don't mess with this
until you feel comfortable with the overall workflow, and also
until you've learned about :ref:`Branching assemblies <advanced_CLI>`.

Now lets run step 3:
.. code:: bash
    ipyrad -p params-ipyrad-test.txt -s 3

.. parsed-literal::
 --------------------------------------------------
  ipyrad [v.0.1.47]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  loading Assembly: ipyrad-test [/private/tmp/ipyrad-test/ipyrad-test.json]
  ipyparallel setup: Local connection to 4 Engines

  Step3: Clustering/Mapping reads
    Saving Assembly.

And we can examine the output:

.. code:: bash
    ipyrad -p params-ipyrad-test.txt -r

.. parsed-literal::                                                                                                                                  
    Summary stats of Assembly ipyrad-test
    ------------------------------------------------
          clusters_hidepth  clusters_total  reads_filtered  reads_raw  state
    1A_0              1000            1000           20099      20099      3
    1B_0              1000            1000           19977      19977      3
    1C_0              1000            1000           20114      20114      3
    1D_0              1000            1000           19895      19895      3
    2E_0              1000            1000           19928      19928      3
    2F_0              1000            1000           19934      19934      3
    2G_0              1000            1000           20026      20026      3
    2H_0              1000            1000           19936      19936      3
    3I_0              1000            1000           20084      20084      3
    3J_0              1000            1000           20011      20011      3
    3K_0              1000            1000           20117      20117      3
    3L_0              1000            1000           19901      19901      3


Step 4: Joint estimation of heterozygosity and error rate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 5: Consensus base calls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
