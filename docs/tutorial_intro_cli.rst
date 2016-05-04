

.. include:: global.rst 

.. _tutorial_intro_cli:


Introductory tutorial - CLI
============================

This is the full introductory tutorial for the command line interface to ipyrad. 
Here we will walk through an entire assembly process. The goal is
to become familiarized with the general workflow, terminology, data files, and 
parameter settings in ipyrad. We will use a single-end RAD-seq data set
as an example, but the core concepts apply to other data types as well 
(e.g., GBS and paired-end). Follow along by copy/pasting the 
code-blocks into a command line terminal. 

.. note:: 

    If you haven't yet installed ipyrad go here first: 
    :ref:`Installation <installation>`


Getting the data
~~~~~~~~~~~~~~~~~
The example data set for this tutorial can be assembled in just a few minutes 
on a typical laptop computer. Use the commands below to download and extract
the data. This will create a new directory called ``ipsimdata/`` located 
in your current directory. 

.. code-block:: bash

    ## The curl command needs a capital O, not a zero
    >>> curl -LkO https://github.com/dereneaton/ipyrad/raw/master/tests/ipsimdata.tar.gz
    >>> tar -xvzf ipsimdata.tar.gz


Use the command `ls` to look inside this directory. You'll see that
it contains many different files representing different test data sets. 

.. code-block:: bash  

    ## the command ls shows you the files inside a directory 
    ## the -l flag shows it as a list.
    >>> ls -l ipsimdata/

.. parsed-literal::
    -rw-rw-r-- 1 deren deren    144 May  3 15:56 gbs_example_barcodes.txt
    -rw-rw-r-- 1 deren deren  26792 May  3 15:58 gbs_example_genome.fa
    -rw-rw-r-- 1 deren deren 840593 May  3 15:57 gbs_example_R1_.fastq.gz
    -rw-rw-r-- 1 deren deren    144 May  3 15:57 pairddrad_example_barcodes.txt
    -rw-rw-r-- 1 deren deren 840599 May  3 15:57 pairddrad_example_R1_.fastq.gz
    -rw-rw-r-- 1 deren deren 849834 May  3 15:57 pairddrad_example_R2_.fastq.gz
    -rw-rw-r-- 1 deren deren    144 May  3 15:57 pairddrad_wmerge_example_barcodes.txt
    -rw-rw-r-- 1 deren deren  36958 May  3 15:58 pairddrad_wmerge_example_genome.fa
    -rw-rw-r-- 1 deren deren 838658 May  3 15:57 pairddrad_wmerge_example_R1_.fastq.gz
    -rw-rw-r-- 1 deren deren 847002 May  3 15:57 pairddrad_wmerge_example_R2_.fastq.gz
    -rw-rw-r-- 1 deren deren    144 May  3 15:57 pairgbs_example_barcodes.txt
    -rw-rw-r-- 1 deren deren 840597 May  3 15:57 pairgbs_example_R1_.fastq.gz
    -rw-rw-r-- 1 deren deren 847587 May  3 15:57 pairgbs_example_R2_.fastq.gz
    -rw-rw-r-- 1 deren deren    144 May  3 15:57 pairgbs_wmerge_example_barcodes.txt
    -rw-rw-r-- 1 deren deren  36958 May  3 15:58 pairgbs_wmerge_example_genome.fa
    -rw-rw-r-- 1 deren deren 838656 May  3 15:57 pairgbs_wmerge_example_R1_.fastq.gz
    -rw-rw-r-- 1 deren deren 844596 May  3 15:57 pairgbs_wmerge_example_R2_.fastq.gz
    -rw-rw-r-- 1 deren deren    144 May  3 15:56 rad_example_barcodes.txt
    -rw-rw-r-- 1 deren deren  26792 May  3 15:58 rad_example_genome.fa
    -rw-rw-r-- 1 deren deren 840593 May  3 15:56 rad_example_R1_.fastq.gz


For this introductory tutorial we will use just two files from 
this directory. The file ``rad_example_R1_.fastq.gz`` contains Illumina 
fastQ formatted reads and is gzip compressed. This is a typical format for raw 
data. The other file, ``rad_example_barcodes.txt``, is a tab-separated table 
matching barcodes to sample IDs. 


Input files
~~~~~~~~~~~
Before we get started let's take a look at what the raw data looks like. Your 
input data will be in fastQ format, usually ending in ``.fq``, ``.fastq``,
``.fq.gz``, or ``.fastq.gz``. It can be split among multiple files, or all 
within a single file. In this tutorial the data are not yet demultiplexed 
(sorted into separate files for each sample), and so we will start by demultiplexing
the data files. For this, we enter the location of our fastQ data files 
on line 2 of the params file ('raw_fastq_path'). If the data were already 
demultiplexed we would instead enter the location of the data files on line 
4 of the params file ("sorted_fastq_path"). 
Below are the first three reads in the example file.

.. code-block:: bash

    ## For your personal edification here is what this is doing:
    ##  gzip -c: Tells gzip to unzip the file and write the contents to the screen
    ##  head -n 12: Grabs the first 12 lines of the fastq file. 

    >>> gunzip -c ./ipsimdata/sim_rad_test_R1_.fastq.gz | head -n 12 


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

These are 100 bp single-end reads prepared as RAD-seq. The first 
six bases form the barcode (TTTTAA) and the next five bases (TGCAG) the 
restriction site overhang. All following bases make up the sequence 
data.

Lets also take a look at the barcodes file for the simulated data. You'll 
see sample names (left) and their barcodes (right) each on a 
separate line with a tab between them.

.. code-block:: bash

    >>> cat ./ipsimdata/sim_rad_test_barcodes.txt

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



Create an ipyrad params file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ipyrad uses a simple text file to hold all the parameters for a given assembly. 
Start by creating a new params file using the ``-n`` flag, followed
by a name for your assembly. In the example we use the name ``iptest``, 
but the name can be anything at all. Once you start analysing your own data 
you might call your params file something more informative, like the name 
of your organism. We will refer to this as the "assembly_name". 

.. code-block:: bash  

    >>> ipyrad -n iptest


.. parsed-literal::
    New file params-iptest.txt created in /home/deren/Documents/ipyrad/tests


This will create a file in the current directory called ``params-iptest.txt``.
The params file lists on each line one parameter followed by a ## mark, 
then the name of the parameter, and then a short description of its 
purpose. Take a look at it by using the unix command 'cat' (or you can
use any text editor you like).

.. code-block:: bash  

    >>> cat params-iptest.txt


.. parsed-literal::
    ------ ipyrad params file (v.0.2.0)---------------------------------------------
    iptest                         ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
    ./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
                                   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                                   ## [3] [barcodes_path]: Location of barcodes file
                                   ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
    denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                                   ## [6] [reference_sequence]: Location of reference sequence file
    rad                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
    TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
    5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
    33                             ## [10] [phred_Qscore_offset]: phred Q score offset (only alternative=64)
    6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
    6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
    10000                          ## [13] [maxdepth]: Max cluster depth within samples
    0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
    1                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
    0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
    35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
    2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
    5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
    8, 8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
    4                              ## [21] [min_samples_locus]: Min # samples per locus for output
    50, 50                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
    8, 8                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
    0.25                           ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
    0, 0                           ## [25] [edit_cutsites]: Edit cut-sites (R1, R2) (see docs)
    1, 2, 2, 1                     ## [26] [trim_overhang]: Trim overhang (see docs) (R1>, <R1, R2>, <R2)
    *                              ## [27] [output_formats]: Output formats (see docs)
                                   ## [28] [pop_assign_file]: Path to population assignment file
                                   ## [29] [excludes]: Samples to be excluded from final output files
                                   ## [30] [outgroups]: Outgroup individuals. Excluded from final output


In general the default parameter values are sensible, and we won't 
mess with them for now, but there are a few parameters we *must* change. 
We need to set the path to the raw data we want to analyse, and we need 
to set the path to the barcodes file.

In your favorite text editor open ``params-iptest.txt`` and change these two lines
to look like this, and then save it. Be careful of typos, if you enter the path
incorrectly ipyrad will raise an error and tell you that it can't find your 
data files:

.. parsed-literal::
    ./ipsimdata/rad_example_R1_.fastq.gz       ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
    ./ipsimdata/rad_example_barcodes.txt       ## [3] [barcodes_path]: Location of barcodes file


Step 1: Demultiplex the raw data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now we will start assembling the data with ipyrad. 
Step 1 reads in the barcodes file and the raw data. It scans through
the raw data and sorts each read based on the mapping of samples to 
barcodes. At the end of this step we'll have a new directory in our project_dir
called ``iptest_fastqs/``. Inside this directory will be individual
fastq.gz files for each sample.

.. note:: 

    You'll notice the name of this output directory bears a strong
    resemblence to the name of the assembly we chose at the time
    of the params file creation. Assembling rad-seq type sequence
    data requires a lot of different steps, and these steps generate a 
    **lot** of intermediary files. ipyrad organizes these files into 
    directories, and it prepends the name of your assembly to each
    directory with data that belongs to it. One result of this is that
    you can have multiple assemblies of the same raw data with different
    parameter settings and you don't have to manage all the files
    yourself! (See :ref:`Branching assemblies <tutorial_advanced_CLI>` for more
    info). Another result is that **you should not rename or move any
    of the directories inside your project directory**, unless you know
    what you're doing or you don't mind if your assembly breaks. 

Now lets run step 1! For the simulated data this will take just a few seconds.

.. code-block:: bash

    ## -p indicates the params file we wish to use
    ## -s indicates the step to run (in this case 1)
    >>> ipyrad -p params-iptest.txt -s 1

.. parsed-literal::
  --------------------------------------------------
   ipyrad [v.0.2.0]
   Interactive assembly and analysis of RADseq data
  --------------------------------------------------
   New Assembly: iptest
   ipyparallel setup: Local connection to 4 Engines
 
   Step1: Demultiplexing fastq data to Samples
   [####################] 100%  sorting reads         | 0:00:04 
   [####################] 100%  writing files         | 0:00:00 
   Saving Assembly.


There are 4 main parts to this step: (1) It creates a new Assembly called iptest, 
since this is our first time running any steps for the named assembly; (2) It 
launches a number of parallel Engines, by default this is the number of available
CPUs on your machine; (3) It performs the step functions, in this case it sorts
the data and writes the outputs; and (4) It saves the Assembly. 

Have a look at the results of this step in the ``iptest_fastqs/``
output directory:

.. code-block:: bash

   >>> ls iptest_fastqs 

.. parsed-literal::
    1A_0_R1_.fastq.gz        1D_0_R1_.fastq.gz        2G_0_R1_.fastq.gz        3J_0_R1_.fastq.gz        s1_demultiplex_stats.txt
    1B_0_R1_.fastq.gz        2E_0_R1_.fastq.gz        2H_0_R1_.fastq.gz        3K_0_R1_.fastq.gz
    1C_0_R1_.fastq.gz        2F_0_R1_.fastq.gz        3I_0_R1_.fastq.gz        3L_0_R1_.fastq.gz

A more informative metric of success might be the number of raw reads 
demultiplexed for each sample. Fortunately ipyrad tracks the state of all your 
steps in your current assembly, so at any time you can ask for results by 
invoking the ``-r`` flag.

.. code-block:: bash

    ## -r fetches informative results from currently executed steps
    >>> ipyrad -p params-iptest.txt -r


.. parsed-literal::
    Summary stats of Assembly iptest
    ------------------------------------------------
          state  reads_raw
    1A_0      1      20144
    1B_0      1      20024
    1C_0      1      20055
    1D_0      1      19927
    2E_0      1      19936
    2F_0      1      20094
    2G_0      1      19930
    2H_0      1      20048
    3I_0      1      19952
    3J_0      1      20164
    3K_0      1      19993
    3L_0      1      20035
    
    
    Full stats files
    ------------------------------------------------
    step 1: ./iptest_fastqs/s1_demultiplex_stats.txt
    step 2: None
    step 3: None
    step 4: None
    step 5: None
    step 6: None
    step 7: None


If you want to get even **more** info ipyrad tracks all kinds of
wacky stats and saves them to a file inside the directories it
creates for each step. For instance to see full stats for step 1:

.. code-block:: bash

    >>> cat ./iptest_fastqs/s1_demultiplex_stats.txt

And you'll see a ton of fun stuff I won't copy here in the interest
of conserving space (you'll see more still on real empirical data versus 
the simulated data here). Please go look for yourself if you're interested.


Step 2: Filter reads
~~~~~~~~~~~~~~~~~~~~
This step filters reads based on quality scores, and can be used to 
detect Illumina adapters in your reads, which is sometimes a problem 
with homebrew type library preparations. Here the filter is set to the 
default value of 0 (zero), meaning it filters only based on quality scores of 
base calls. The filtered files are written to a new directory called 
``iptest_edits/``.

.. code-block:: bash

    >>> ipyrad -p params-iptest.txt -s 2

.. parsed-literal::
  --------------------------------------------------
   ipyrad [v.0.2.0]
   Interactive assembly and analysis of RADseq data
  --------------------------------------------------
   loading Assembly: iptest
   from saved path: ~/Documents/ipyrad/tests/iptest.json
   ipyparallel setup: Local connection to 4 Engines
 
   Step2: Filtering reads 
   [####################] 100%  processing reads      | 0:00:33 
   Saving Assembly.


Again, you can look at the results output by this step and also some 
handy stats tracked for this assembly.

.. code-block:: bash

    ## View the output of step 2
    >>> ls iptest_edits/

.. parsed-literal::                                                                                                                                  
    1A_0_R1_.fastq       1C_0_R1_.fastq       2E_0_R1_.fastq       2G_0_R1_.fastq       3I_0_R1_.fastq       3K_0_R1_.fastq       s2_rawedit_stats.txt
    1B_0_R1_.fastq       1D_0_R1_.fastq       2F_0_R1_.fastq       2H_0_R1_.fastq       3J_0_R1_.fastq       3L_0_R1_.fastq

.. code-block:: bash

    ## Get current stats including # raw reads and # reads
    ## after filtering.
    >>> ipyrad -p params-iptest.txt -r

.. parsed-literal::
    Summary stats of Assembly iptest
    ------------------------------------------------
          state  reads_raw  reads_filtered
    1A_0      2      20144           20144
    1B_0      2      20024           20024
    1C_0      2      20055           20055
    1D_0      2      19927           19927
    2E_0      2      19936           19936
    2F_0      2      20094           20094
    2G_0      2      19930           19930
    2H_0      2      20048           20048
    3I_0      2      19952           19952
    3J_0      2      20164           20164
    3K_0      2      19993           19993
    3L_0      2      20035           20035
    
    
    Full stats files
    ------------------------------------------------
    step 1: ./iptest_fastqs/s1_demultiplex_stats.txt
    step 2: ./iptest_edits/s2_rawedit_stats.txt
    step 3: None
    step 4: None
    step 5: None    
    step 6: None
    step 7: None


You might also take a gander at the filtered reads:

.. parsed-literal::
    >>> head -n 12 ./iptest_edits/1A_0_R1_.fastq


Step 3: clustering within-samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 3 de-replicates and then clusters reads within each sample 
by the set clustering threshold and then writes the clusters to new 
files in a directory called ``iptest_clust_0.85/``. Intuitively
we are trying to identify all the reads that map to the same locus
within each sample. The clustering threshold specifies the minimum 
percentage of sequence similarity below which we will consider two 
reads to have come from different loci.

The true name of this output directory will be dictated by the value
you set for the ``clust_threshold`` parameter in the params file. 

.. parsed-literal::
    0.85             ## [14] [clust_threshold]: proportion identical for clustering

You can see the default value is 0.85, so our default directory is 
named accordingly. This value dictates the percentage of sequence
similarity that reads must have in order to be considered reads
at the same locus. You may want to experiment
with this value, but 0.85-0.90 is a fairly reliable range, balancing
over-splitting of loci vs over-lumping. Don't mess with this
until you feel comfortable with the overall workflow, and also
until you've learned about :ref:`Branching assemblies <branching_workflow>`.

Later you will learn how to incorporate information from a reference 
genome (if you have one) to improve clustering at this step. For now, bide your
time (but see :ref:`Reference sequence mapping <tutorial_advanced_cli>` if 
you're impatient).

Now lets run step 3:

.. code-block:: bash

    >>> ipyrad -p params-iptest.txt -s 3

.. parsed-literal::
  --------------------------------------------------
   ipyrad [v.0.2.0]
   Interactive assembly and analysis of RADseq data
  --------------------------------------------------
   loading Assembly: iptest
   from saved path: ~/Documents/ipyrad/tests/iptest.json
   ipyparallel setup: Local connection to 4 Engines
 
   Step3: Clustering/Mapping reads
   [####################] 100%  dereplicating         | 0:00:01 
   [####################] 100%  clustering            | 0:00:03 
   [####################] 100%  chunking              | 0:00:00 
   [####################] 100%  aligning              | 0:00:44 
   [####################] 100%  concatenating         | 0:00:00 
   Saving Assembly.


Again we can examine the results. The stats output tells you how many clusters 
were found, and the number of clusters that pass the mindepth thresholds. 
We'll go into more detail about mindepth settings in some of the advanced tutorials
but the important thing to know is that by default step 3 will filter out clusters
that only have a handful of reads since we have little power to make confident
base calls in low depth clusters. 

.. code-block:: bash

    >>> ipyrad -p params-iptest.txt -r

.. parsed-literal:: 
    Summary stats of Assembly iptest
    ------------------------------------------------
          state  reads_raw  reads_filtered  clusters_total  clusters_hidepth
    1A_0      3      20144           20144            1000              1000
    1B_0      3      20024           20024            1000              1000
    1C_0      3      20055           20055            1000              1000
    1D_0      3      19927           19927            1000              1000
    2E_0      3      19936           19936            1000              1000
    2F_0      3      20094           20094            1000              1000
    2G_0      3      19930           19930            1000              1000
    2H_0      3      20048           20048            1000              1000
    3I_0      3      19952           19952            1000              1000
    3J_0      3      20164           20164            1000              1000
    3K_0      3      19993           19993            1000              1000
    3L_0      3      20035           20035            1000              1000
    

    Full stats files
    ------------------------------------------------
    step 1: ./iptest_fastqs/s1_demultiplex_stats.txt
    step 2: ./iptest_edits/s2_rawedit_stats.txt
    step 3: ./iptest_clust_0.85/s3_cluster_stats.txt
    step 4: None
    step 5: None
    step 6: None
    step 7: None


The aligned clusters found during this step are now located in 
``./iptest_clust_0.85/``. You can get a feel for what this looks like
by examining a portion of one of the files using the command below.

.. code-block:: bash                                                                                    

    ## Same as above, gunzip -c means print to the screen and 
    ## `head -n 28` means just show me the first 28 lines. If 
    ## you're interested in what more of the loci look like
    ## you can increase the number of lines you ask head for,
    ## e.g. ... | head -n 100
    >>> gunzip -c iptest_clust_0.85/1A_0.clustS.gz | head -n 28

Reads that are sufficiently similar (based on the above sequence similarity 
threshold) are grouped together in clusters separated by "//". For the first
cluster below there is clearly one allele (homozygote) and one read with a 
(simulated) sequencing error. This is apparent in the 'size=' field of the two
reads for this cluster. For the second cluster it seems there are two alleles 
(heterozygote), and a read with a sequencing error. For the third 
cluster it's a bit harder to say. Is this a homozygote with lots of sequencing
errors, or a heterozygote with few reads for one of the alleles? Thankfully,
untangling this mess is what steps 4 and 5 are all about.

.. parsed-literal::
    >1A_0_1164_r1;size=16;*0
    TGCAGCTATTGCGACAAAAACACGACGGCTTCCGTGGGCACTAGCGTAATTCGCTGAGCCGGCGTAACAGAAGGAGTGCACTGCCACGTGCCCG
    >1A_0_1174_r1;size=1;+1
    TGCAGCTATTGCGACAAAAACACGACGGCTTCCGTGGGCACTAGCGTAATTCGCTGAGCCGGCGTAACAGAAGGAGTGCACTGCCACATGCCCG
    //
    //
    >1A_0_4137_r1;size=10;*0
    TGCAGGGTCGCCGGCAACTCAGCATTTTAACTCCGCGGGTTACACGTGCGGAGGCCTACTGGCTATCATTTTTAGGGTGCATTTGGTCGGCTGG
    >1A_0_4130_r1;size=6;+1
    TGCAGGGTCGCCGGCAACTCAGCATTTTAACTCCGCGGGTTACACGTGTGGAGGCCTACTGGCTATCATTTTTAGGGTGCATTTGGTCGGCTGG
    >1A_0_4131_r1;size=1;+2
    TGCAGGGTCGCCGGCAACTCAGCATTTTAACTCCGCGGGTTACACGTGTCGAGGCCTACTGGCTATCATTTTTAGGGTGCATTTGGTCGGCTGG
    //
    //
    >1A_0_6246_r1;size=15;*0
    TGCAGATACAAAAGCTTGCCCACTAAGTTGTGTGATCACTGTCTTATTACGGTGGCCTCCTTCAAGCTTCGAACGAGTTGTGGATCGGTAGGCT
    >1A_0_6259_r1;size=1;+1
    TGCAGATACAAAAGCTTGCCCACTAAGTTGTGTGATCACTGTCTTATTACGGTGGCCTCCTTCAAGCTTCGAACGAGTTGTGGATCGGGAGGCT
    >1A_0_6264_r1;size=1;+2
    TGCAGATACAAAAGCTTGCCCACTAAGTTGTGTGATCACTGTCTTATTACGGTGGCCTCCTACAAGCTTCGAACGAGTTGTGGATCGGTAGGCT
    >1A_0_6268_r1;size=1;+3
    TGCAGATTCAAAAGCTTGCCCACTAAGTTGTGTGATCACTGTCTTATTACGGTGGCCTCCTTCAAGCTTCGAACGAGTTGTGGATCGGTAGGCT
    //
    //

Step 4: Joint estimation of heterozygosity and error rate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Jointly estimate sequencing error rate and heterozygosity to disentangle
which reads are "real" and which are sequencing error. We need to know
which reads are "real" because in diploid organisms there are a maximum of 2
alleles at any given locus. If we look at the raw data and there are 5 or 
ten different "alleles", and 2 of them are very high frequency, and the rest 
are singletons then this gives us evidence that the 2 high frequency alleles 
are good reads and the rest are probably not. This step is pretty 
straightforward, and pretty fast. Run it thusly:

.. code-block:: bash

    >>> ipyrad -p params-iptest.txt -s 4

.. parsed-literal::
  --------------------------------------------------
   ipyrad [v.0.2.0]
   Interactive assembly and analysis of RADseq data
  --------------------------------------------------
   loading Assembly: iptest
   from saved path: ~/Documents/ipyrad/tests/iptest.json
   ipyparallel setup: Local connection to 4 Engines
 
   Step4: Joint estimation of error rate and heterozygosity
   [####################] 100%  inferring [H, E]      | 0:01:09 
   Saving Assembly.


This step does not produce new output files, only a stats file with the 
estimated heterozygosity and error rate parameters. 
You can also invoke the ``-r`` flag to see the estimated values.

.. code-block:: bash

    >>> ipyrad -p params-iptest.txt -r

.. parsed-literal::
    Summary stats of Assembly iptest
    ------------------------------------------------
          state  reads_raw  reads_filtered  clusters_total  clusters_hidepth
    1A_0      4      20144           20144            1000              1000
    1B_0      4      20024           20024            1000              1000
    1C_0      4      20055           20055            1000              1000
    1D_0      4      19927           19927            1000              1000
    2E_0      4      19936           19936            1000              1000
    2F_0      4      20094           20094            1000              1000
    2G_0      4      19930           19930            1000              1000
    2H_0      4      20048           20048            1000              1000
    3I_0      4      19952           19952            1000              1000
    3J_0      4      20164           20164            1000              1000
    3K_0      4      19993           19993            1000              1000
    3L_0      4      20035           20035            1000              1000
    
          hetero_est  error_est  
    1A_0    0.001866   0.000739  
    1B_0    0.001988   0.000728  
    1C_0    0.001833   0.000733  
    1D_0    0.001719   0.000752  
    2E_0    0.002034   0.000726  
    2F_0    0.001866   0.000748  
    2G_0    0.001899   0.000733  
    2H_0    0.002261   0.000767  
    3I_0    0.001886   0.000776  
    3J_0    0.001718   0.000764  
    3K_0    0.002046   0.000744  
    3L_0    0.001796   0.000719  
    
    
    Full stats files
    ------------------------------------------------
    step 1: ./iptest_fastqs/s1_demultiplex_stats.txt
    step 2: ./iptest_edits/s2_rawedit_stats.txt
    step 3: ./iptest_clust_0.85/s3_cluster_stats.txt
    step 4: ./iptest_clust_0.85/s4_joint_estimate.txt
    step 5: None
    step 6: None
    step 7: None



Step 5: Consensus base calls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 5 uses the inferred error rate and heterozygosity to call the consensus
of sequences within each cluster. Here we are identifying what we believe
to be the real haplotypes at each locus within each sample.

.. code-block:: bash

    ipyrad -p params-iptest.txt -s 5

.. parsed-literal::

    --------------------------------------------------
     ipyrad [v.0.1.73]
     Interactive assembly and analysis of RADseq data
    --------------------------------------------------
     loading Assembly: iptest [~/Documents/ipyrad/tests/iptest.json]
     ipyparallel setup: Local connection to 4 Engines

     Step5: Consensus base calling 
       Diploid base calls and paralog filter (max haplos = 2)
       error rate (mean, std):  0.00075, 0.00002
       heterozyg. (mean, std):  0.00198, 0.00018
       Saving Assembly.

Again we can ask for the results:

.. code-block:: bash

    ipyrad -p params-iptest.txt -r

And here the important information is the number of ``reads_consens``. This is 
the number of "good" reads within each sample that we'll send on to the next step.
As you'll see in examples with empirical data, this is often a step where many
reads are filtered out of the data set. If not data were filtered, then the 
number of reads_consens should be equal to the number of clusters_hidepth.

.. parsed-literal::

    Summary stats of Assembly iptest
        ------------------------------------------------
              state  reads_raw  reads_filtered  clusters_total  clusters_hidepth
    1A_0      5      20099           20099            1000              1000
    1B_0      5      19977           19977            1000              1000
    1C_0      5      20114           20114            1000              1000
    1D_0      5      19895           19895            1000              1000
    2E_0      5      19928           19928            1000              1000
    2F_0      5      19934           19934            1000              1000
    2G_0      5      20026           20026            1000              1000
    2H_0      5      19936           19936            1000              1000
    3I_0      5      20084           20084            1000              1000
    3J_0      5      20011           20011            1000              1000
    3K_0      5      20117           20117            1000              1000
    3L_0      5      19901           19901            1000              1000

          hetero_est  error_est  reads_consens
    1A_0    0.002223   0.000756           1000
    1B_0    0.001910   0.000775           1000
    1C_0    0.002260   0.000751           1000
    1D_0    0.001876   0.000731           1000
    2E_0    0.001809   0.000770           1000
    2F_0    0.002103   0.000725           1000
    2G_0    0.001910   0.000707           1000
    2H_0    0.002215   0.000755           1000
    3I_0    0.001877   0.000784           1000
    3J_0    0.001698   0.000741           1000
    3K_0    0.001821   0.000767           1000
    3L_0    0.002003   0.000755           1000


    Full stats files
    ------------------------------------------------
    step 1: ./iptest_fastqs/s1_demultiplex_stats.txt
    step 2: ./iptest_edits/s2_rawedit_stats.txt
    step 3: ./iptest_clust_0.85/s3_cluster_stats.txt
    step 4: ./iptest_clust_0.85/s4_joint_estimate.txt
    step 5: ./iptest_consens/s5_consens_stats.txt
    step 6: None
    step 7: None


Step 6: Cluster across samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 6 clusters consensus sequences across samples. Now that we have good 
estimates for haplotypes within samples we can try to identify similar sequences
at each locus between samples. We use the same clustering threshold as step 3
to identify sequences between samples that are probably sampled from the 
same locus, based on sequence similarity.

.. code-block:: bash

    ipyrad -p params-iptest.txt -s 6

.. parsed-literal::

    --------------------------------------------------
     ipyrad [v.0.1.73]
     Interactive assembly and analysis of RADseq data
    --------------------------------------------------
     loading Assembly: iptest [~/Documents/ipyrad/tests/iptest.json]
     ipyparallel setup: Local connection to 4 Engines

     Step6: Clustering across 12 samples at 0.85 similarity

This step differs from previous steps in that we are no longer applying the
steps to each Sample individually, but instead we are applying a function
to all of the Samples collectively. Our end result is that we now have a 
map telling us which loci cluster together from which Samples. This step 
does not produce a stats output. It might be more enlightening to consider 
the output of step 6 by examining the file that contains the aligned 
reads clustered across samples:

.. code-block:: bash

    gunzip -c iptest_consens/iptest_catclust.gz | head -n 30 | less

The final output of step 6 is a file in ``iptest_consens/`` called 
``iptest_catclust.gz``. There is actually also a database file called
``iptest_test.hdf5`` which records not only the clustered sequences but also
information about read depths at each site, but that file is not 
human-readable. The file ``iptest_catclust.gz`` on the other hand provides 
a good sanity check for examining the aligned reads across all samples. 
Executing the above command you'll see the output below which 
shows all the reads that align at one particular locus. You'll see the 
sample name of each read followed by the sequence of the read at that locus
for that sample. If you wish to examine more loci you can increase the number
of lines you want to view by increasing the value you pass to ``head`` in
the above command (e.g. ``... | head -n 300 | less``

.. parsed-literal::

    1C_0_691
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGAGTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    3L_0_597
    TGCAGGGTGGGTKGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTGTAATCGAGTATTAGCGCGGAAGC
    2E_0_339
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    2F_0_994
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    3K_0_941
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    1B_0_543
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    3J_0_357
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    2H_0_106
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    3I_0_202
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTACTAGCGCGGAAGC
    2G_0_575
    TGCAGSGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGKTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    1D_0_744
    TGCAGGGTGGGTGGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    1A_0_502
    TGCAGGGTGGGTTGTGTTATTTAACATCCAATGCTTAAAGTTTCGATTAGGGGCCTGTTACCGTAGAGTTTTAATCGAGTATTAGCGCGGAAGC
    //
    //

Step 7: Filter and write output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The final step is to filter the data and write output files in many 
convenient file formats. First we apply filters for maximum number of 
indels per locus, max heterozygosity per locus, max number of snps 
per locus, and minimum number of samples per locus. All these filters
are configurable in the params file and you are encouraged to explore 
different settings, but the defaults are quite good and quite conservative.

After running step 7 like so:

.. code-block:: bash

    ipyrad -p params-iptest.txt -s 7

A new directory is created called ``iptest_outfiles``. This directory contains
all the output files specified in the params file. The default is to 
create all supported output files which include 
.phy, .nex, .geno, .str, as well as many others (forthcoming). Explore each of
these files as you wish.

Congratulations! You've completed your first toy assembly. Now you can try applying
what you've learned to assemble your own real data. Please consult the docs for many
of the more powerful features of ipyrad including reference sequence mapping, 
assembly branching, and post-processing analysis including svdquartets and 
many population genetic summary statistics.
