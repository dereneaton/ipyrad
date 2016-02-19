

.. _tutorial_intro_cli:


Introductory tutorial - CLI
============================

This is the full tutorial for the command line interface for ipyrad. In this
tutorial we'll walk through the entire assembly and analysis process. This is 
meant as a broad introduction to familiarize users with the general workflow,
and some of the parameters and terminology. For simplicity we'll use 
single-end RAD-Seq as the example data, but the core concepts will apply
to assembly of other data types (GBS and paired-end). 

If you are new to RADseq analyses, this tutorial will provide a simple overview 
of how to execute ipyrad, what the data files look like, how to check that 
your analysis is working, and what the final output formats will be. You can 
follow along by copy/pasting the code-blocks into a command line terminal. 


Getting Started
~~~~~~~~~~~~~~~

.. _note::

    If you haven't already installed ipyrad go here first: 
    :ref:`Installation <installation>`

.. _note::

    If you haven't already installed ipyrad go here first


.. _warning::
    How do I make warnings?

.. _DANGER::
    Is it all about making it in CAPS?

.. _DANGER::

    or Is it about making it in CAPS and spaces?


We provide a very small sample data set that we recommend using for this tutorial.
Full datasets can take several hours to several days to complete, 
whereas with the simulated data you can complete an assembly in just 
a few minutes. 

First download and extract a set of example data from the web using the command 
below. This will create a directory called ``ipsimdata/`` in your current directory
containing a number of test data sets. 

.. code-block:: bash

    ## The curl command needs a capital O, not a zero
    curl -O https://github.com/dereneaton/ipyrad/blob/master/tests/ipsimdata.tar.gz
    tar -xvzf ipsimdata.tar.gz


This directory contains many simulated datasets, as well as a simulated 
reference genome that we will use in other tutorials. For this introductory
tutorial we will use just the following two files from this directory:

    - ``sim_rad_test_R1_.fastq.gz`` - Illumina fastQ formatted reads (gzip compressed)
    - ``sim_rad_test_barcodes.txt`` - Mapping of barcodes to sample IDs


Create a new parameters file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ipyrad uses a text file to hold all the parameters for a given assembly. 
Start by creating a new parameters file with the ``-n`` flag. This flag
requires you to pass in a name for your assembly. In the example we use 
``ipyrad-test`` but the name can be anything at all. Once you start 
analysing your own data you might call your parameters file something 
more informative, like the name of your organism.

.. code-block:: bash

    ipyrad -n ipyrad-test

This will create a file in the current directory called ``params-ipyrad-test.txt``.
The params file lists on each line one parameter followed by a ## mark, 
then the name of the parameter, and  then a short description of its 
purpose. Lets take a look at it.

.. code-block:: bash

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

In your favorite text editor open ``params-ipyrad-test.txt`` and change these two lines
to look like this, and then save it:

.. parsed-literal::
    ./data/sim_rad_test_R1_.fastq.gz         ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
    ./data/sim_rad_test_barcodes.txt         ## [3] [barcodes_path]: Location of barcodes file

Input data format
~~~~~~~~~~~~~~~~~
Before we get started let's take a look at what the raw data looks like.

Your input data will be in fastQ format, usually ending in ``.fq``, ``.fastq``,
``.fq.gz``, or ``.fastq.gz``. Your data could be split among multiple files, or all 
within a single file (de-multiplexing goes much faster if they happen to 
be split into multiple files). The file/s may be compressed with gzip so 
that they have a .gz ending, but they do not need to be. The location of 
these files should be entered on line 2 of the params file. Below are 
the first three reads in the example file.

.. code-block:: bash

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
called ``ipyrad-test_fastqs``. Inside this directory will be individual
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
info). Another result is that **you should not rename or move any
of the directories inside your project directory**, unless you know
what you're doing or you don't mind if your assembly breaks. 

Lets take a look at the barcodes file for the simulated data. You'll 
see sample names (left) and their barcodes (right) each on a 
separate line with a tab between them.

.. code-block:: bash

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

Now lets run step 1! For the simulated data this will take < 1 minute.

.. code-block:: bash

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
    - Create a new assembly. Since this is our first time running any steps we need to initialize our assembly.
    - Start the parallel cluster. ipyrad uses a parallelization library called ipyparallel. Every time we start a step we fire up the parallel clients. This makes your assemblies go **smokin'** fast.
    - Actually do the demuliplexing.
    - Save the state of the assembly.

Have a look at the results of this step in the ``ipyrad-test_fastqs``
output directory:

.. code-block:: bash

   ls ipyrad-test_fastqs 

.. parsed-literal::
    1A_0_R1_.fastq.gz        1D_0_R1_.fastq.gz        2G_0_R1_.fastq.gz        3J_0_R1_.fastq.gz        s1_demultiplex_stats.txt
    1B_0_R1_.fastq.gz        2E_0_R1_.fastq.gz        2H_0_R1_.fastq.gz        3K_0_R1_.fastq.gz
    1C_0_R1_.fastq.gz        2F_0_R1_.fastq.gz        3I_0_R1_.fastq.gz        3L_0_R1_.fastq.gz

A more informative metric of success might be the number
of raw reads demultiplexed for each sample. Fortunately 
ipyrad tracks the state of all your steps in your current 
assembly, so at any time you can ask for results by 
invoking the ``-r`` flag.

.. code-block:: bash

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

.. code-block:: bash

    cat ./ipyrad-test_fastqs/s1_demultiplex_stats.txt

And you'll see a ton of fun stuff I won't copy here in the interest
of conserving space. Please go look for yourself if you're interested.

Step 2: Filter reads
~~~~~~~~~~~~~~~~~~~~
This step filters reads based on quality scores, and can be used to 
detect Illumina adapters in your reads, which is sometimes a problem 
with homebrew type library preparations. Here the filter is set to the 
default value of 0 (zero), meaning it filters only based on quality scores of 
base calls. The filtered files are written to a new directory called 
``ipyrad-test_edits``.

.. code-block:: bash

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

.. code-block:: bash

    ## View the output of step 2
    ls ipyrad-test_edits

.. parsed-literal::                                                                                                                                  
    1A_0_R1_.fastq       1C_0_R1_.fastq       2E_0_R1_.fastq       2G_0_R1_.fastq       3I_0_R1_.fastq       3K_0_R1_.fastq       s2_rawedit_stats.txt
    1B_0_R1_.fastq       1D_0_R1_.fastq       2F_0_R1_.fastq       2H_0_R1_.fastq       3J_0_R1_.fastq       3L_0_R1_.fastq

.. code-block:: bash

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
.. code-block:: bash

    head -n 12 ./ipyrad-test_fastqs/1A_0_R1_.fastq


Step 3: clustering within-samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 3 de-replicates and then clusters reads within each sample 
by the set clustering threshold and then writes the clusters to new 
files in a directory called ``ipyrad-test_clust_0.85``. Intuitively
we are trying to identify all the reads that map to the same locus
within each sample. The clustering threshold specifies the minimum 
percentage of sequence similarity below which we will consider two 
reads to have come from different loci.

The true name of this output directory will be dictated by the value
you set for the ``clust_threshold`` parameter in the params file. 

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

Later you will learn how to incorporate information from a reference 
genome to improve clustering at this this step. For now, bide your
time (but see :ref:`Reference sequence mapping <advanced_CLI>` if 
you're impatient).

Now lets run step 3:

.. code-block:: bash

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

Again we can examine the results. The stats output tells you how many clusters 
were found, and the number of clusters that pass the mindepth thresholds. 
We'll go into more detail about mindepth settings in some of the advanced tutorials
but for now all you need to know is that by default step 3 will filter out clusters
that only have a handful of reads on the assumption that these are probably
all mostly due to sequencing error.

.. code-block:: bash

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

Again, the final output of step 3 is dereplicated, clustered files for each sample 
in ``./ipryad-test_clust_0.85/``. You can get a feel for what this looks like
by examining a portion of one of the files.

.. code-block:: bash                                                                                                                                 

    ## Same as above, gunzip -c means print to the screen and 
    ## `head -n 28` means just show me the first 28 lines. If 
    ## you're interested in what more of the loci look like
    ## you can increase the number of lines you ask head for,
    ## e.g. ... | head -n 100
    gunzip -c ipyrad-test_clust_0.85/1A_0.clustS.gz | head -n 28

Reads that are sufficiently similar (based on the above sequence similarity 
threshold) are grouped together in clusters separated by "//". For the first
cluster below there is clearly one allele (homozygote) and one read with a 
(simulated) sequencing error. For the second cluster it seems there are two alleles 
(heterozygote), and a couple reads with sequencing errors. For the third 
cluster it's a bit harder to say. Is this a homozygote with lots of sequencing
errors, or a heterozygote with few reads for one of the alleles?

Thankfully, untangling this mess is what step 4 is all about.

.. parsed-literal::
    >1A_0_1164_r1;size=16;*0
    TGCAGCTATTGCGACAAAAACACGACGGCTTCCGTGGGCACTAGCGTAATTCGCTGAGCCGGCGTAACAGAAGGAGTGCACTGCCACGTGCCCG
    >1A_0_1174_r1;size=1;+1
    TGCAGCTATTGCGACAAAAACACGACGGCTTCCGTGGGCACTAGCGTAATTCGCTGAGCCGGCGTAACAGAAGGAGTGCACTGCCACATGCCCG
    //
    //
    >1A0_8280_r1;size=10;
    TGCAGCGTATATGATCAGAACCGGGTGAGTGGGTACCGCGAACCGAAAGGCATCGAAAGTTTAGCGCAGCACTAATCTCA
    >1A0_8290_r1;size=8;+
    TGCAGCGTATATGATCAGAACCGGGTGAGTGGGTACCGCGAACCGAAAGGCACCGAAAGTTTAGCGCAGCACTAATCTCA
    >1A0_8297_r1;size=1;+
    TGCAGCGTATATGATCAGAACCGGGTGAGTGGGAACCGCGAACCGAAAGGCACCGAAAGTTTAGCGCAGCACTAATCTCA
    >1A0_8292_r1;size=1;+
    TGCAGCCTATATGATCAGAACCGGGTGAGTGGGTACCGCGAACCGAAAGGCACCGAAAGTTTAGCGCAGCACTAATCTCA
    //
    //
    >1A_0_2982_r1;size=17;*0
    TGCAGACGTGGAGTAACCGGCGGCCTTTAGTCTTAGTAGTGTCCGGGGTACCCGTTGGTTTGTCGTAGTGAGTTCGGTAGGCAAACTTCTGGCC
    >1A_0_2983_r1;size=1;+1
    TGCAGACGTGGAGTATCCGGCGGCCTTTAGTCTTAGTAGTGTCCGGGGTACCCGTTGGTTTGTCGTAGTGAGTTCGGTAGGCAAACTTCTGGCC
    >1A_0_2985_r1;size=1;+2
    TGCAGACGTGGAGTAACCGGCGGCCTTTAGTCTAAGTAGTGTCCGGGGTACCCGTTGGTTTGTCGTAGTGAGTTCGGTAGGCAAACTTCTGGCC
    >1A_0_2988_r1;size=1;+3
    TGCAGACGAGGAGTAACCGGCGGCCTTTAGTCTTAGTAGTGTCCGGGGTACCCGTTGGTTTGTCGTAGTGAGTTCGGTAGGCAAACTTCTGGCC
    >1A_0_3002_r1;size=1;+4
    TGCAGACGTGGAGCAACCGGCGGCCTTTAGTCTTAGTAGTGTCCGGGGTACCCGTTGGTTTGTCGTAGTGAGTTCGGTAGGCAAACTTCTGGCC
    //
    //


Step 4: Joint estimation of heterozygosity and error rate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Jointly estimate sequencing error rate and heterozygosity to help us figure
out which reads are "real" and which are sequencing error. We need to know
which reads are "real" because in diploid organisms there are a maximum of 2
alleles at any given locus. If we look at the raw data and there are 5 or 
ten different "alleles", and 2 of them are very high frequency, and the rest 
are singletons then this gives us evidence that the 2 high frequency alleles 
are good reads and the rest are probably junk. This step is pretty straightforward, 
and pretty fast. Run it thusly:

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -s 4

.. parsed-literal::
 --------------------------------------------------                                                                                                  
  ipyrad [v.0.1.47]                                                                                                                                  
  Interactive assembly and analysis of RADseq data                                                                                                   
 --------------------------------------------------                                                                                                  
  loading Assembly: ipyrad-test [/private/tmp/ipyrad-test/ipyrad-test.json]                                                                          
  ipyparallel setup: Local connection to 4 Engines                                                                                                   
                                                                                                                                                     
  Step4: Joint estimation of error rate and heterozygosity                                                                                           
    Saving Assembly.

In terms of results, there isn't as much to look at as in previous steps, though
you can invoke the ``-r`` flag to see the estimated heterozygosity and error
rate per sample.

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -r

.. parsed-literal::
    Summary stats of Assembly ipyrad-test
    ------------------------------------------------
          clusters_hidepth  clusters_total  error_est  hetero_est  reads_filtered
    1A_0              1000            1000   0.000757    0.002212           20099
    1B_0              1000            1000   0.000774    0.001883           19977
    1C_0              1000            1000   0.000745    0.002223           20114
    1D_0              1000            1000   0.000734    0.001894           19895
    2E_0              1000            1000   0.000778    0.001800           19928
    2F_0              1000            1000   0.000728    0.002082           19934
    2G_0              1000            1000   0.000707    0.001825           20026
    2H_0              1000            1000   0.000756    0.002190           19936
    3I_0              1000            1000   0.000778    0.001848           20084
    3J_0              1000            1000   0.000739    0.001705           20011 
    3K_0              1000            1000   0.000768    0.001857           20117
    3L_0              1000            1000   0.000756    0.001979           19901 


Step 5: Consensus base calls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 5 uses the inferred error rate and heterozygosity to call the consensus
of sequences within each cluster. Here we are identifying what we believe
to be the real haplotypes at each locus within each sample.

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -s 5

.. parsed-literal::                                                                                                                                  
 --------------------------------------------------                                                                                                  
  ipyrad [v.0.1.47]                                                                                                                                  
  Interactive assembly and analysis of RADseq data                                                                                                   
 --------------------------------------------------                                                                                                  
  loading Assembly: ipyrad-test [/private/tmp/ipyrad-test/ipyrad-test.json]                                                                          
  ipyparallel setup: Local connection to 4 Engines                                                                                                   
                                                                                                                                                     
  Step5: Consensus base calling                                                                                                                      
    Diploid base calls and paralog filter (max haplos = 2)                                                                                           
    error rate (mean, std):  0.00075, 0.00002                                                                                                        
    heterozyg. (mean, std):  0.00196, 0.00018                                                                                                        
    Saving Assembly. 

Again we can ask for the results:

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -r

And here the important information is the number of ``reads_consens``. This is 
the number of "good" reads within each sample that we'll send on to the next step.

.. parsed-literal::
          clusters_hidepth  clusters_total  error_est  hetero_est  reads_consens
    1A_0              1000            1000   0.000757    0.002212           1000
    1B_0              1000            1000   0.000774    0.001883           1000
    1C_0              1000            1000   0.000745    0.002223           1000
    1D_0              1000            1000   0.000734    0.001894           1000
    2E_0              1000            1000   0.000778    0.001800           1000
    2F_0              1000            1000   0.000728    0.002082           1000
    2G_0              1000            1000   0.000707    0.001825           1000
    2H_0              1000            1000   0.000756    0.002190           1000
    3I_0              1000            1000   0.000778    0.001848           1000
    3J_0              1000            1000   0.000739    0.001705           1000
    3K_0              1000            1000   0.000768    0.001857           1000
    3L_0              1000            1000   0.000756    0.001979           1000

Step 6: Cluster across samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 6 clusters consensus sequences across samples. Now that we have good 
estimates for haplotypes within samples we can try to identify similar sequences
at each locus between samples. We use the same clustering threshold as step 3
to identify sequences between samples that are probably sampled from the same locus,
based on sequence similarity.

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -s 6

.. parsed-literal::
 --------------------------------------------------
  ipyrad [v.0.1.47]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  loading Assembly: ipyrad-test [/private/tmp/ipyrad-test/ipyrad-test.json]
  ipyparallel setup: Local connection to 4 Engines

  Step6: Clustering across 12 samples at 0.85 similarity
    Saving Assembly.

Since in general the stats for results of each step are sample based, the 
output of  ``-r`` at this point is less useful. You can still try it though.

.. code-block:: bash

    ipyrad -p params-ipyrad-test.txt -r

It might be more enlightening to consider the output of step 6 by examining
the file that contains the reads clustered across samples:

.. code-block:: bash

    gunzip -c ipyrad-test_consens/ipyrad-test_catclust.gz | head -n 30 | less

The final output of step 6 is a file in ``ipyrad-test_consens`` called 
``ipyrad-test_catclust.gz``. This file contains all aligned reads across
all samples. Executing the above command you'll see the output below which 
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

    ipyrad -p params-ipyrad-test.txt -s 7

A new directory is created called ``ipyrad-test_outfiles``. This directory contains
all the output files specified in the params file. The default is to 
create all supported output files which include .phy, .nex, .geno, .treemix, .str, as
well as many others.

Congratulations! You've completed your first toy assembly. Now you can try applying
what you've learned to assemble your own real data. Please consult the docs for many
of the more powerful features of ipyrad including reference sequence mapping, 
assembly branching, and post-processing analysis including svdquartets and 
many population genetic summary statistics.
