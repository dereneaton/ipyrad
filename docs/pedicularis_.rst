
.. include:: global.rst

.. _pedicularis_cli:


Eaton & Ree (2013) data set
===========================
In this tutorial we demonstrate a denovo assembly for an empirical data set 
to give a general idea of the typical results you might expect to recover 
and typical run times. This example is run on a 4-core laptop with 8GB RAM, 
to show that you do not need a super computer to assemble your data. However, 
using more cores will improve the speed of ipyrad approximately 
linearly, so if you have access to a large cluster I recommend using it. 

Here we use the 13 taxa *Pedicularis* data set from **Eaton and Ree (2013)**. 
This data set is composed of single-end 75bp reads for a 
RAD-seq library prepared with the PstI enzyme for 13 individuals.  
The original paper is available open-access: 
(:ref:`eaton_and_ree<eaton_and_ree>`). 


Download the empirical example data set (*Pedicularis*)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These data are archived on the NCBI sequence read archive (SRA) under 
accession id SRP021469. For convenience, the data are also hosted at a 
public Dropbox link which is a bit easier to access. Run the code below to 
download and decompress the fastq data files, which will save them into a 
directory called ``example_empirical_data/``. The compressed file size is 
approximately 1.1GB.

.. code:: bash

    ## curl grabs the data from a public dropbox url
    ## the curl command uses an upper-case o argument, not a zero.
    >>> curl -LkO https://dl.dropboxusercontent.com/u/2538935/example_empirical_rad.tar.gz
    
    ## the tar command decompresses the data directory
    >>> tar -xvzf example_empirical_rad.tar.gz


Setup a base params file
~~~~~~~~~~~~~~~~~~~~~~~~
We start by using the ``-n`` argument to create a new named Assembly. 
I use the name ``base`` to indicate that this is the base assembly from 
which we will later create several branches.

.. code:: bash

    >>> ipyrad -n "base"

.. parsed-literal::
    New file 'params-base.txt' created in /home/deren/Downloads


The data come to us already demultiplexed so we are going to simply set the 
**sorted\_fastq\_path** to tell ipyrad the location of the data files, 
and also set a **project\_dir**, which will group all of our analyses into 
a single directory. For the latter I use the name of our study organism, 
"pedicularis". 

.. parsed-literal::
    ## Use your text editor to enter the following values:
    ## The wildcard (*) tells ipyrad to select all files ending in .gz
    pedicularis                       ## [1] [project_dir] ...
    example_empirical_rad/*.gz        ## [4] [sorted_fastq_path] ...

For now we'll leave the remaining parameters at their default values.


Step 1: Load the fastq data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You always start an ipyrad assembly by running step 1. 
When the data location is entered as a **sorted_fastq_path** (param 4), as 
opposed to the **raw_fastq_path** (param 2), step 1 simply counts the 
number of reads for each Sample and parses the file names to 
extract names for each Sample. For example, the file ``29154_superba.fastq.gz`` 
will be assigned to Sample ``29154_superba``. Now, run step 1 (-s 1) and 
tell ipyrad to print the results when it is finished (-r). 

.. code:: bash

    >>> ipyrad -p params-base.txt -s 1 -r


.. parsed-literal:: 
 --------------------------------------------------
  ipyrad [v.0.3.15]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  loading Assembly: base
  from saved path: ~/Downloads/pedicularis/base.json
  New Assembly: base
  ipyparallel setup: Local connection to 4 Engines

  Step1: Linking sorted fastq data to Samples
    Linking to demultiplexed fastq files in:
      /home/deren/Downloads/example_empirical_rad/*.gz
    13 new Samples created in base.
    13 fastq files linked to 13 new Samples.
  Saving Assembly.

  Summary stats of Assembly base
  ------------------------------------------------
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


Finish the assembly
~~~~~~~~~~~~~~~~~~~
Because the point of this tutorial is to demonstrate run times and 
statistics, I will leave the rest of the parameters at their
defaults and simply run all remaining steps. Further below I will 
explain in more detail the stats files for each step and what the values mean. 
To fully assemble this data set on a 4-core laptop takes about 5 hours. With 
access to 24 cores it would take less than 1 hour. 


.. code:: bash

    ## run steps 2-7
    >>> ipyrad -p params-base.txt -s 234567


.. parsed-literal::
 --------------------------------------------------
  ipyrad [v.0.3.15]
  Interactive assembly and analysis of RADseq data
 --------------------------------------------------
  loading Assembly: base
  from saved path: ~/Downloads/pedicularis/base.json
  ipyparallel setup: Local connection to 4 Engines

  Step2: Filtering reads 
  [####################] 100%  processing reads      | 0:37:50 

  Step3: Clustering/Mapping reads
  [####################] 100%  dereplicating         | 0:00:31 
  [####################] 100%  clustering            | 0:19:07
  [####################] 100%  chunking              | 0:00:01 
  [####################] 100%  aligning              | 2:53:23 
  [####################] 100%  concatenating         | 0:00:14 

  Step4: Joint estimation of error rate and heterozygosity
  [####################] 100%  inferring [H, E]      | 0:08:17 

  Step5: Consensus base calling 
  Mean error  [0.00339 sd=0.00168]
  Mean hetero [0.01908 sd=0.00297]
  [####################] 100%  consensus calling     | 0:18:17 

  Step6: Clustering across 13 samples at 0.85 similarity
  [####################] 100%  concat/shuffle input  | 0:00:05 
  [####################] 100%  clustering across     | 0:03:17 
  [####################] 100%  building clusters     | 0:00:43 
  [####################] 100%  aligning clusters     | 0:07:33 
  [####################] 100%  indexing clusters     | 0:24:42 
  [####################] 100%  building database     | 0:06:25 

  Step7: Filter and write output files for 13 Samples
  [####################] 100%  filtering loci        | 0:03:43 
  [####################] 100%  building loci/stats   | 0:00:05 
  [####################] 100%  building vcf file     | 0:04:41 
  [####################] 100%  writing outfiles      | 0:00:12 
  Outfiles written to: ~/Downloads/pedicularis/base_outfiles


  Summary stats of Assembly base
  ------------------------------------------------
                          state  reads_raw  reads_filtered  clusters_total
  29154_superba               6     696994          644514          121823   
  30556_thamno                6    1452316         1358296          191762   
  30686_cyathophylla          6    1253109         1109175          212955   
  32082_przewalskii           6     964244          887289          136218   
  33413_thamno                6     636625          570041          153729   
  33588_przewalskii           6    1002923          929348          143611   
  35236_rex                   6    1803858         1668629          382823   
  35855_rex                   6    1409843         1307294          159765   
  38362_rex                   6    1391175         1294460          121620   
  39618_rex                   6     822263          756623          134339   
  40578_rex                   6    1707942         1591828          202486   
  41478_cyathophylloides      6    2199740         2056389          156560   
  41954_cyathophylloides      6    2199613         2008385          270071   
  
                          clusters_hidepth  hetero_est  error_est  reads_consens  
  29154_superba                      32503    0.021164   0.003767          30934  
  30556_thamno                       49918    0.020225   0.003297          47289  
  30686_cyathophylla                 49563    0.017679   0.002749          46957  
  32082_przewalskii                  38720    0.022536   0.004018          36590  
  33413_thamno                       27742    0.021708   0.002996          26259  
  33588_przewalskii                  42843    0.022167   0.003103          40701  
  35236_rex                          51861    0.020122   0.002004          49786  
  35855_rex                          53107    0.020848   0.006000          49974  
  38362_rex                          50231    0.014790   0.002639          48381  
  39618_rex                          40723    0.019533   0.003938          38877  
  40578_rex                          53346    0.018743   0.001154          51263  
  41478_cyathophylloides             52494    0.013695   0.001281          50426  
  41954_cyathophylloides             70699    0.014845   0.007090          67582    



Take a look at the stats summary 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each assembly that finishes step 7 will create a stats.txt output summary
in the 'assembly_name'_outfiles/ directory. This includes information about 
which filters removed data from the assembly, how many loci were recovered
per sample, how many samples had data for each locus, and how many variable
sites are in the assembled data. 


.. code:: python

    cat ./pedicularis/base_outfiles/base_stats.txt


.. parsed-literal::

  ## The number of loci caught by each filter.
  ## ipyrad API location: [assembly].statsfiles.s7_filters
  
                              locus_filtering
  total_prefiltered_loci                84221
  filtered_by_rm_duplicates              1304
  filtered_by_max_indels                 1750
  filtered_by_max_snps                      0
  filtered_by_max_shared_het             2865
  filtered_by_min_sample                37299
  filtered_by_max_alleles                9112
  total_filtered_loci                   38656


  ## The number of loci recovered for each Sample.
  ## ipyrad API location: [assembly].stats_dfs.s7_samples

                          sample_coverage
  29154_superba                     19269
  30556_thamno                      29914
  30686_cyathophylla                24779
  32082_przewalskii                 13600
  33413_thamno                      16453
  33588_przewalskii                 15988
  35236_rex                         31414
  35855_rex                         31159
  38362_rex                         31640
  39618_rex                         25565
  40578_rex                         32010
  41478_cyathophylloides            29699
  41954_cyathophylloides            26647


  ## The number of loci for which N taxa have data.
  ## ipyrad API location: [assembly].stats_dfs.s7_loci
  
      locus_coverage  sum_coverage
  1                0             0
  2                0             0
  3                0             0
  4             4892          4892
  5             3611          8503
  6             3198         11701
  7             2929         14630
  8             3117         17747
  9             4130         21877
  10            4840         26717
  11            5300         32017
  12            4290         36307
  13            2349         38656


  ## The distribution of SNPs (var and pis) across loci.
  ## var = all variable sites (pis + autapomorphies)
  ## pis = parsimony informative site (minor allele in >1 sample)
  ## ipyrad API location: [assembly].stats_dfs.s7_snps

       var  sum_var   pis  sum_pis
  0   1707        0  9461        0
  1   3344     3344  8974     8974
  2   4400    12144  6801    22576
  3   4681    26187  4522    36142
  4   4527    44295  3071    48426
  5   3984    64215  2010    58476
  6   3396    84591  1339    66510
  7   2857   104590   906    72852
  8   2328   123214   625    77852
  9   1863   139981   427    81695
  10  1537   155351   244    84135
  11  1188   168419   135    85620
  12   910   179339    60    86340
  13   606   187217    34    86782
  14   417   193055    22    87090
  15   318   197825    12    87270
  16   208   201153     8    87398
  17   135   203448     4    87466
  18   100   205248     1    87484
  19    59   206369     0    87484
  20    38   207129     0    87484
  21    31   207780     0    87484
  22    12   208044     0    87484
  23     5   208159     0    87484
  24     1   208183     0    87484
  25     4   208283     0    87484




Take a peek at the .loci output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the first places I look when an assembly finishes. It provides a 
clean view of the data with variable sites (-) and parsimony informative
SNPs (*) highlighted. Use the unix commands **less** or **head** to look at this
file briefly. Each locus is labelled with a number corresponding to the locus 
order before filters are applied in step 7. If you branch this assembly and 
run step 7 again with a different set of parameters you may recover fewer 
or more total loci. 


.. code:: bash

  ## head -n 50 prints just the first 50 lines of the file to stdout
  head -n 50 pedicularis/base_outfiles/base.loci

.. parsed-literal::
  
    30686_cyathophylla         ATGCAGGAG-ATCAAACATA-CGACAAAGAAATAAATTTATTGGATTCGTAGTGGATTTATTCTTTTTT
    32082_przewalskii          ATGCAGGTATWTCMWWCATATCGACAAAGAAATAAAKTTTTTGGATTCGTAGTGGATTTATTCTTYTTT
    33413_thamno               ATGCAGGNGTATCAAACATATTGACAAAGAAATAAATTTACTGGATTCGTAGTGGNTTTATTCTTTTTT
    33588_przewalskii          ATGCAGGTATATCAAACATATCGACAAAGAAATAAATTTTTTGGATTCGTAGTTTATTCTTTCTTTTTT
    35236_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGTGGATTTATTCTTTTTT
    35855_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTAYTGGATTCGTAGTGGATTTATTCKTTTTT
    38362_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGTGGGTTTATTCTTTTTT
    40578_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGTGGATTTATTCGTTTTT
    41478_cyathophylloides     ATGCAGGAG-AAGAAACATA-CGACAAAGAAATAAATTTATTGGATTCGTAGTGGATTTATTCTTTTTT
    41954_cyathophylloides     ATGCAGGAG-AAGAAACATA-CGACANNGAAATAAATTTATTGGATTCGTAGTGGATTTATTCTTTTTT
    //                                *\ * -* *---     -              -  *\ *            ---  --   * -   |1|
    29154_superba              GTTCTGGAGTTGTTCAGGGTACTGTTGTGCAGCCATTGCAGCAAATTGAGCCGCAACATCAGCATCAGC
    30686_cyathophylla         GTTCTGGAGTTGTTCAGGGTAGTGTTGTGCAGCCATTGCAGCAAATTGACCCGCAACATCAGCATCAGC
    32082_przewalskii          GTTCTGGAGTTGTTCAGGGTAGTGTTGTGCAGCCATTGCAGCAGATTGAGCCG---CATCAGGATGGTC
    33588_przewalskii          GTTCTGGAGTTGTTCAGGGTAGTGTTGTGCAGCCATTGCAGCAGATTGAGCCG---CATCAGGATGGTC
    35855_rex                  GTTCTGGGGTTGTTCAAGGTAGTGTTGTGCAACCATTGCAGCAGATTGAGCCGCAACATCAGCATCAGT
    41478_cyathophylloides     GTTCTGGAGTTGTTCAGGGCAGTGTTGTGCAGCCATTGCAGCAAATTGAGCCGCAACATCAGCATCAGC
    41954_cyathophylloides     GTTCTGGAGTTGTTCAGGGCAGTGTTGTGCAGCCATTGCAGCAAATTGAGCCGCAACATCAGCATCAGC
    //                                -        -  * -         -           *     -            *  *\ *\ *-|4|
    30556_thamno               GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCATGTT
    30686_cyathophylla         GAGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTGTTTATACACTTTCCCTTCCTGTT
    32082_przewalskii          GTGATCCATTATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    33413_thamno               GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    33588_przewalskii          GTGATCCATTATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    35236_rex                  GTGATCAATAATTCACTTTCCCGCTTCACRCKAACAACACASACTATTTATACACTTTCCCTTCCTGTT
    35855_rex                  GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    38362_rex                  GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    39618_rex                  GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    40578_rex                  GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    41478_cyathophylloides     GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    41954_cyathophylloides     GTGATCAATAATTCACTTTCCCGCTTCACGCTAACAACACAGACTATTTATACACTTTCCCTTCCTGTT
    //                          -    *  *                   - -         -   -                  -    |6|
    33413_thamno               TCAGATACGAGAGAAA-AAAATACACATAATTGAAGCTGTATAGAAAAAGAATAATATGTGAGAAAGGA
    33588_przewalskii          ACAGATACGAGAGAAATACATTACACATAGTTGAAGCTGTATAGAAAAGGAATGATATGTGACAAAGGA
    35855_rex                  TCAGATACGAGAGAAA-AAAATACACATAATTGAAGCTGTATAGAAAAAGAATAATATGTGAGAAGGGA
    39618_rex                  TCAGATACGAGAGAAA-AAAATACACATAATTGAAGCTGTATAGAAAAAGAATAATATGTGAGAAAGGA
    //                         -                 - -        -                  -    -        -  -   |15|
    35236_rex                  ATGGAGAGTTCATAGCAATTCCTGATGAGTTTTTAATACAAAAAACAGGAGATCCTTTCAAACTTATCT
    35855_rex                  ATGGAGAGTTCATAGCAATTCCTGATGAGTTTTTAATACAAAAAACAGGAGATCCTTTCAAACTTATCT
    40578_rex                  ATGGAGAGTTCATAGCAATTCCTGATGAGTTTTTAATACAAAAAACAGGAGATCCTTTCAAACTTATCT
    41478_cyathophylloides     ATGGAGAATTCATAGCAGTTCCTTATGAGTTTTTAATACAAAAAACAGGCGATCCTTTCAAACTTATCT
    41954_cyathophylloides     ATGGAGAATTCATAGCAGTTCCTTATGAGTTTTTAATACAAAAAACAGGCGATCCTTTCAAACTTATCT
    //                                *         *     *                         *                   |17|
    30556_thamno               AAAACCCAAAGGAAATAGAGGCATTGAATCATGTATCAAGGGAGAAATATTTGTTTGATCACTAAAATA
    35236_rex                  AAAACCCAAAGGAAATAGAGGCATTGAATCATGTATCAAGGGAGAAATATTTGTTTGATCACCAAAATA
    35855_rex                  AAAACCCAAAGGAAATTGAGGCATTGANTCATGTNTCAAGGGAGAAATATTTGTGTGATCACCAAAATA
    39618_rex                  AAAACCCAAAGGAAATTGAGGCATTGAATCATGTATCAAGGGAGAAATATTTGTGTGATCACCAAAATA
    41478_cyathophylloides     AAAACCCAAAGGAAATAGAGGCATTGAATCATGTATCAAGGGAGAAATATTTGTTTGATCACCTTGTTT
    41954_cyathophylloides     AAAACCCAAAGGAAATAGAGGCATTGAATCATGTATCAAGGGAGAAATATTTGTTTGATCACCTTGTTT
    //                                         *                                     *       -*\ *\ *\ * *|18|


peek at the .phy files
~~~~~~~~~~~~~~~~~~~~~~
This is the concatenated sequence file of all loci in the data set. It is typically
used in phylogenetic analyses, like in the program *raxml*. This super matrix is
13 taxon deep by 2.68 Mb long. 


.. code:: bash

    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/base_outfiles/base.phy


.. parsed-literal::

  13 2684753
  29154_superba              NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  30556_thamno               NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  30686_cyathophylla         ATGCAGGAGNATCAAACATANCGACAAAGAAATAAATTTATTGGATTCGTAGT
  32082_przewalskii          ATGCAGGTATWTCMWWCATATCGACAAAGAAATAAAKTTTTTGGATTCGTAGT
  33413_thamno               ATGCAGGNGTATCAAACATATTGACAAAGAAATAAATTTACTGGATTCGTAGT
  33588_przewalskii          ATGCAGGTATATCAAACATATCGACAAAGAAATAAATTTTTTGGATTCGTAGT
  35236_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGT
  35855_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTAYTGGATTCGTAGT
  38362_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGT
  39618_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  40578_rex                  ATGCAGGTGTATCAAACATATCGACAAAGAAATAAATTTATTGGATTCGTAGT
  41478_cyathophylloides     ATGCAGGAGNAAGAAACATANCGACAAAGAAATAAATTTATTGGATTCGTAGT
  41954_cyathophylloides     ATGCAGGAGNAAGAAACATANCGACANNGAAATAAATTTATTGGATTCGTAGT



peek at the .snps.phy file
~~~~~~~~~~~~~~~~~~~~~
This is similar to the phylip file format, but only variable site columns are 
included. All SNPs are in the file, in contrast to the .u.snps.phy file, which 
randomly selects only a single SNP per locus. 


.. code:: bash

    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/base_outfiles/base.snps.phy


.. parsed-literal::

  13 208283
  29154_superba              NNNNNNNNNNNNNNNNNNNAGTCGAGCCAGCNNNNNNNNNNNNNNNNNNNNNN
  30556_thamno               NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAAGTGAANNNNNNNNNNNNAT
  30686_cyathophylla         AGATCAAACTATGGATATTAGTGGACCCAGCAAAGTGGCNNNNNNNNNNNNNN
  32082_przewalskii          TAWTCMWWCKTTGGATATYAGTGGGGGGGTCTCTGTGACNNNNNNNNNNNNNN
  33413_thamno               NGATCAAATTACGGNTATTNNNNNNNNNNNNTAAGTGACTAAAAAGANNNNNN
  33588_przewalskii          TAATCAAACTTTTTACTTTAGTGGGGGGGTCTCTGTGACACTGGGCANNNNNN
  35236_rex                  TGATCAAACTATGGATATTNNNNNNNNNNNNTAARKSACNNNNNNNNGAGAAT
  35855_rex                  TGATCAAACTAYGGATAKTGATGAGGCCAGTTAAGTGACTAAAAAGGGAGATG
  38362_rex                  TGATCAAACTATGGGTATTNNNNNNNNNNNNTAAGTGACNNNNNNNNNNNNNN
  39618_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTAAGTGACTAAAAAGANNNNTG
  40578_rex                  TGATCAAACTATGGATAGTNNNNNNNNNNNNTAAGTGACNNNNNNNNGAGANN
  41478_cyathophylloides     AGAAGAAACTATGGATATTAGCGGAGCCAGCTAAGTGACNNNNNNNNAGTCAT
  41954_cyathophylloides     AGAAGAAACTATGGATATTAGCGGAGCCAGCTAAGTGACNNNNNNNNAGTCAT


peek at the .vcf.gz file
~~~~~~~~~~~~~~~~~~~~~
The VCF output for ipyrad contains the full sequence information for all samples
as well as the sequencing depth information for all base calls that were made. 
This file should be easily parsable if users want to extract information or 
modify it so that this file can be used in other software such as GATK. We 
are working on developing our own population-aware genotype caller that will
correct low-depth base calls at this stage. Stay tuned. 


.. code:: bash

    ## gunzip -c decompresses the file and passes it to the pipe (|)
    ## head -n 50 reads data from the pipe and show the first 50 lines. 
    ## and we pipe this to 'cut', which shows only the first 80 rows of data
    ## for easier viewing. 

    gunzip -c pedicularis/base_outfiles/base.vcf.gz | head -n 50 | cut -c 1-80


.. parsed-literal::

  ##fileformat=VCFv4.0
  ##fileDate=2016/06/29
  ##source=ipyrad_v.0.3.15
  ##reference=pseudo-reference (most common base at site)
  ##phasing=unphased
  ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
  RAD_1_  0 . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  1 . T . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  2 . G . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 0/0:0,0,0,11  0/0:0,0,0,6 0/0:0,0,0,7 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  3 . C . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:6,0,0,0 0/0:11,0,0,0  0/0:6,0,0,0 0/0:7,0,0,0 0/0:20,0,0,0  0/0:15,0,0,0  0/0:24,0,0,0  ./.:0,0,0,0 0/0:11,0,0,0  0/0:13,0,0,0  0/0:6,0,0,0
  RAD_1_  4 . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  5 . G . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 0/0:0,0,0,11  0/0:0,0,0,6 0/0:0,0,0,7 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  6 . G . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 0/0:0,0,0,11  0/0:0,0,0,6 0/0:0,0,0,7 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  7 . T A 13  PASS  NS=10;DP=112  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 1/1:0,6,0,0 0/0:0,0,6,0 ./.:0,0,5,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,14,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  1/1:0,13,0,0  1/1:0,6,0,0
  RAD_1_  8 . G A 13  PASS  NS=10;DP=114  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 1/1:0,6,0,0 0/0:0,0,0,6 1/1:0,7,0,0 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  9 . T . 13  PASS  NS=7;DP=94  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  ./.:0,0,0,0 ./.:0,0,0,0
  RAD_1_  10  . A T 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 1/0:0,6,5,0 0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  11  . T A 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  1/1:0,13,0,0  1/1:0,6,0,0
  RAD_1_  12  . C G 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:6,0,0,0 0/0:11,0,0,0  0/0:6,0,0,0 0/0:7,0,0,0 0/0:20,0,0,0  0/0:15,0,0,0  0/0:24,0,0,0  ./.:0,0,0,0 0/0:11,0,0,0  1/1:0,0,0,13  1/1:0,0,0,6
  RAD_1_  13  . A C 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 1/0:5,6,0,0 0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  14  . A T 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 1/0:0,6,5,0 0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  15  . A T 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 1/0:0,6,5,0 0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  16  . C . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:6,0,0,0 0/0:11,0,0,0  0/0:6,0,0,0 0/0:7,0,0,0 0/0:20,0,0,0  0/0:15,0,0,0  0/0:24,0,0,0  ./.:0,0,0,0 0/0:11,0,0,0  0/0:13,0,0,0  0/0:6,0,0,0
  RAD_1_  17  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  18  . T . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  19  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  20  . T . 13  PASS  NS=7;DP=94  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  ./.:0,0,0,0 ./.:0,0,0,0
  RAD_1_  21  . C T 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:6,0,0,0 0/0:11,0,0,0  1/1:0,0,6,0 0/0:7,0,0,0 0/0:20,0,0,0  0/0:15,0,0,0  0/0:24,0,0,0  ./.:0,0,0,0 0/0:11,0,0,0  0/0:13,0,0,0  0/0:6,0,0,0
  RAD_1_  22  . G . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 0/0:0,0,0,11  0/0:0,0,0,6 0/0:0,0,0,7 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  23  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  24  . C . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:6,0,0,0 0/0:11,0,0,0  0/0:6,0,0,0 0/0:7,0,0,0 0/0:20,0,0,0  0/0:15,0,0,0  0/0:24,0,0,0  ./.:0,0,0,0 0/0:11,0,0,0  0/0:13,0,0,0  0/0:6,0,0,0
  RAD_1_  25  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  26  . A . 13  PASS  NS=10;DP=118  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  ./.:0,5,0,0
  RAD_1_  27  . A . 13  PASS  NS=10;DP=118  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  ./.:0,5,0,0
  RAD_1_  28  . G . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,0,6 0/0:0,0,0,11  0/0:0,0,0,6 0/0:0,0,0,7 0/0:0,0,0,20  0/0:0,0,0,15  0/0:0,0,0,24  ./.:0,0,0,0 0/0:0,0,0,11  0/0:0,0,0,13  0/0:0,0,0,6
  RAD_1_  29  . A . 13  PASS  NS=10;DP=118  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,12,0,0  0/0:0,6,0,0
  RAD_1_  30  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  31  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  32  . T . 13  PASS  NS=10;DP=117  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,14,0  0/0:0,0,23,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  33  . A . 13  PASS  NS=10;DP=118  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,14,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  34  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  35  . A . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 0/0:0,11,0,0  0/0:0,6,0,0 0/0:0,7,0,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0
  RAD_1_  36  . T G 13  PASS  NS=10;DP=118  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 1/0:0,0,6,5 0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,14,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  37  . T . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  38  . T . 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,0,6,0 0/0:0,0,11,0  0/0:0,0,6,0 0/0:0,0,7,0 0/0:0,0,20,0  0/0:0,0,15,0  0/0:0,0,24,0  ./.:0,0,0,0 0/0:0,0,11,0  0/0:0,0,13,0  0/0:0,0,6,0
  RAD_1_  39  . A T 13  PASS  NS=10;DP=119  GT:CATG ./.:0,0,0,0 ./.:0,0,0,0 0/0:0,6,0,0 1/1:0,0,11,0  0/0:0,6,0,0 1/1:0,0,7,0 0/0:0,20,0,0  0/0:0,15,0,0  0/0:0,24,0,0  ./.:0,0,0,0 0/0:0,11,0,0  0/0:0,13,0,0  0/0:0,6,0,0



.. downstream analyses
.. ~~~~~~~~~~~~~~~~~~~
.. We are in the process of developing many downstream tools for analyzing 
.. RAD-seq data. The first of which is the program svd4tet, which is installed 
.. alongside ipyrad during the conda installation. 
.. This program implements the svdquartets algorithm of Chifmann & Kubatko (2014). 
.. It includes alternative methods for sampling quartets over
.. very large trees to heuristically reduce the total number of quartets needed 
.. in order to resolve a large tree. For this demonstration we'll simply run 
.. the default option which is to infer all quartets. 
.. The summary tree that it spits out is the same tree inferred in Eaton & Ree 
.. (2013), and if we plot the tree with support values you will see that we find 
.. lower support across the edges of the tree that are known to be involved in 
.. introgression. This part of the tutorial is under construction -- more to come. 

.. .. code:: bash

..     ## run svd4tet on the unlinked snps file (-s) (.u.snps.phy suffix)
..     ## and give it the output prefix (-o) 'pedictree'
..     svd4tet -s pedicularis/base_outfiles/base.u.snps.phy -o pedictree

.. .. parsed-literal::
    
