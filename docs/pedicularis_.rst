
.. include:: global.rst

.. _pedicularis_cli:


Eaton & Ree (2013) single-end RAD data set
==========================================

**This tutorial is in the process of being updated.**  

Here we demonstrate a denovo assembly for an empirical RAD data set to 
give a general idea of the typical results you might expect to recover 
and typical run times. This example was run on a 4-core laptop with 8GB RAM, 
and takes about 2.25 to run completely, showing that you do not need
a super computer to assemble many data sets. However, using more cores 
will improve the speed of ipyrad approximately linearly, so if you have 
access to a large cluster go ahead and use it.

We will use the 13 taxa *Pedicularis* data set from **Eaton and Ree (2013)**.
This data set is composed of single-end 75bp reads from a RAD-seq library 
prepared with the PstI enzyme for 13 individuals. The original paper is 
available open-access: (:ref:`Eaton & Ree<eatonAndRee>`). 

This data set also serves as an example for several of our **cookbook recipes** that 
demonstrate downstream analyses methods. So after you finish this assembly head over
there to check out fun ways to analyze the resulting data files.


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


Setup a params file
~~~~~~~~~~~~~~~~~~~~~~~~
Always start by using the ``-n {name}`` argument to create a new named Assembly. 
I'll use the name ``base`` to indicate that this is the base assembly from which
we will later create several branches. 

.. code:: bash

    >>> ipyrad -n "base"

.. parsed-literal::
    New file 'params-base.txt' created in /home/deren/tests


The data come to us already demultiplexed so we are going to simply set the 
**sorted\_fastq\_path** to tell ipyrad the location of the data files. You can 
select multiple files at once using regular expressions, in this example we
use an asterisk (`*.gz`) to select all files in the directory ending in .gz. We also
set a **project\_dir**, which is useful for grouping all our results into a single 
directory. For this we'll use the name of our study organism, "pedicularis". 
Take note when entering the values below into your params file that they 
correspond to parameters 1 and 4, respectively.

.. parsed-literal::
    ## Use your text editor to enter the following values:
    ## The wildcard (*) tells ipyrad to select all files ending in .gz
    pedicularis                       ## [1] [project_dir] ...
    example_empirical_rad/*.gz        ## [4] [sorted_fastq_path] ...  


For now we'll leave the remaining parameters at their default values. 
.. See some of the cookbook recipes to see this data set assembled with other params. 


Step 1: Load the fastq data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Start an ipyrad assembly by running step 1. When the data location is entered 
as a **sorted_fastq_path** (param 4), as opposed to the **raw_fastq_path** 
(param 2), step 1 simply counts the number of reads for each Sample and 
parses the file names to extract names for each Sample. For example, the 
file ``29154_superba.fastq.gz`` will be assigned to Sample ``29154_superba``.
We use the `-s` argument followed by `1` to tell ipyrad to run step 1. We also
pass it the `-r` argument so that it will print a results summary when finished. 

.. code:: bash

    >>> ipyrad -p params-base.txt -s 1 -r


.. parsed-literal:: 
 -------------------------------------------------------------
  ipyrad [v.0.4.1]
  Interactive assembly and analysis of RAD-seq data
 -------------------------------------------------------------
  New Assembly: base
  local compute node: [4 cores] on oud

  Step 1: Loading sorted fastq data to Samples
  Loading data from: /home/deren/Downloads/example_empirical_rad/*.gz
  [####################] 100%  loading reads         | 0:00:13
  13 fastq files loaded to 13 Samples.

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


Run the remaining assembly steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Because the point of this tutorial is to demonstrate run times and 
statistics, I will leave the rest of the parameters at their
defaults and simply run all remaining steps. Further below I will 
explain in more detail the stats files for each step and what the values mean. 
To fully assemble this data set on a 4-core laptop takes about 2.25 hours. 
With access to 24 cores it would take only about 20 minutes.


.. code:: bash

    ## run steps 2-7
    >>> ipyrad -p params-base.txt -s 234567


.. parsed-literal::
  --------------------------------------------------------------------
  ipyrad [v.0.4.1]
  Interactive assembly and analysis of RAD-seq data
 --------------------------------------------------------------------
  loading Assembly: pedic
  from saved path: ~/Downloads/pedicularis/pedic.json
  local compute node: [4 cores] on oud

  Step 2: Filtering reads 
  [####################] 100%  processing reads      | 0:03:50 

  Step 3: Clustering/Mapping reads
  [####################] 100%  dereplicating         | 0:00:31 
  [####################] 100%  clustering            | 0:26:46 
  [####################] 100%  building clusters     | 0:00:43 
  [####################] 100%  chunking              | 0:00:01 
  [####################] 100%  aligning              | 1:00:40 
  [####################] 100%  concatenating         | 0:00:12 

  Step4: Joint estimation of error rate and heterozygosity
  [####################] 100%  inferring [H, E]      | 0:09:17 

  Step 5: Consensus base calling 
  Mean error  [0.00309 sd=0.00112]
  Mean hetero [0.01651 sd=0.00282]
  [####################] 100%  consensus calling     | 0:18:17 

  Step 6: Clustering across 13 samples at 0.9 similarity
  [####################] 100%  concat/shuffle input  | 0:00:06  
  [####################] 100%  clustering across     | 0:05:15  
  [####################] 100%  building clusters     | 0:00:06  
  [####################] 100%  aligning clusters     | 0:08:25  
  [####################] 100%  database indels       | 0:00:10  
  [####################] 100%  indexing clusters     | 0:00:43  
  [####################] 100%  building database     | 0:00:01 

  Step 7: Filter and write output files for 13 Samples
  [####################] 100%  filtering loci        | 0:00:13  
  [####################] 100%  building loci/stats   | 0:00:01  
  [####################] 100%  building vcf file     | 0:01:37  
  [####################] 100%  writing vcf file      | 0:00:01  
  [####################] 100%  writing outfiles      | 0:00:03  
  Outfiles written to: ~/Downloads/pedicularis/base_outfiles


  Summary stats of Assembly base
  ------------------------------------------------
                          state  reads_raw  reads_passed_filter  clusters_total  
  29154_superba               6     696994               689987          136851 
  30556_thamno                6    1452316              1440303          215671   
  30686_cyathophylla          6    1253109              1206943          245142
  32082_przewalskii           6     964244               955477          154639
  33413_thamno                6     636625               626082          174937
  33588_przewalskii           6    1002923               993868          161186
  35236_rex                   6    1803858              1787337          422427
  35855_rex                   6    1409843              1397062          187048
  38362_rex                   6    1391175              1379600          136588
  39618_rex                   6     822263               813983          149946
  40578_rex                   6    1707942              1695515          225191 
  41478_cyathophylloides      6    2199740              2185350          176495
  41954_cyathophylloides      6    2199613              2176179          306151

                          clusters_hidepth  hetero_est  error_est  reads_consens  
  29154_superba                      35021    0.017512   0.002765          34598  
  30556_thamno                       52659    0.019254   0.003597          51726  
  30686_cyathophylla                 54136    0.015769   0.003694          53301  
  32082_przewalskii                  42375    0.017746   0.004038          41839  
  33413_thamno                       30889    0.018300   0.002378          30511  
  33588_przewalskii                  46501    0.017635   0.002902          45901  
  35236_rex                          54736    0.015610   0.001742          54197  
  35855_rex                          56583    0.022157   0.005108          55519  
  38362_rex                          53160    0.013540   0.001944          52623  
  39618_rex                          43924    0.016187   0.003047          43381  
  40578_rex                          56662    0.016197   0.001921          56050  
  41478_cyathophylloides             55294    0.010843   0.002080          54746  
  41954_cyathophylloides             76585    0.013932   0.004930          75645  


  Full stats files
  ------------------------------------------------
  step 1: None
  step 2: ./pedicularis/linktest_edits/s2_rawedit_stats.txt
  step 3: ./pedicularis/linktest_clust_0.9/s3_cluster_stats.txt
  step 4: ./pedicularis/linktest_clust_0.9/s4_joint_estimate.txt
  step 5: ./pedicularis/linktest_consens/s5_consens_stats.txt
  step 6: ./pedicularis/linktest_consens/s6_cluster_stats.txt
  step 7: ./pedicularis/linktest_outfiles/linktest_stats.txt


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

                               total_filters  applied_order  retained_loci
   total_prefiltered_loci              96461              0          96461
   filtered_by_rm_duplicates            3794           3794          92667
   filtered_by_max_indels                132            116          92551
   filtered_by_max_snps                   38              9          92542
   filtered_by_max_shared_het            864            658          91884
   filtered_by_min_sample              47725          45906          45978
   filtered_by_max_alleles             10441           4546          41432
   total_filtered_loci                 41432              0          41432


   ## The number of loci recovered for each Sample.
   ## ipyrad API location: [assembly].stats_dfs.s7_samples

                                     sample_coverage
   29154_superba                     21015
   30556_thamno                      31804
   30686_cyathophylla                26616
   32082_przewalskii                 13314
   33413_thamno                      18633
   33588_przewalskii                 15523
   35236_rex                         33236
   35855_rex                         33372
   38362_rex                         33763
   39618_rex                         28073
   40578_rex                         34129
   41478_cyathophylloides            30896
   41954_cyathophylloides            28266


   ## The distribution of SNPs (var and pis) across loci.
   ## var = all variable sites (pis + autapomorphies)
   ## pis = parsimony informative site (minor allele in >1 sample)
   ## ipyrad API location: [assembly].stats_dfs.s7_snps

        var  sum_var    pis  sum_pis
   0   2899        0  12974        0
   1   5312     5312  11192    11192
   2   6670    18652   7533    26258
   3   6570    38362   4498    39752
   4   5700    61162   2639    50308
   5   4567    83997   1303    56823
   6   3519   105111    730    61203
   7   2440   122191    352    63667
   8   1571   134759    132    64723
   9    979   143570     47    65146
   10   544   149010     22    65366
   11   322   152552      6    65432
   12   160   154472      2    65456
   13    97   155733      2    65482
   14    42   156321      0    65482
   15    18   156591      0    65482
   16    13   156799      0    65482
   17     6   156901      0    65482
   18     1   156919      0    65482
   19     2   156957      0    65482



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

.. code:: bash 

  29154_superba              TCTGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTTTCGATCTCAGGCG
  30556_thamno               TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCTCAGGCG
  30686_cyathophylla         TCCAGTCCCGCGGGTGATCAAGGCCCCACCACCGCATCTCACATTCTCGATCTCAGGCG
  33413_thamno               TCCGGTCCTTCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCTCAGGCG
  35236_rex                  TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTMGATCTCAGGCG
  35855_rex                  TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCTCAGGCG
  38362_rex                  TCCGGTCCTTCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCTCAGGCG
  40578_rex                  TCCGGTCCYKCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTCGATCTCAGGCG
  41478_cyathophylloides     TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTATCGATCTCAGGCG
  41954_cyathophylloides     TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTATCGATCTCAGGCG
  //                           --    **                         -         * *           |1|
  29154_superba              TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTAGTATTGTGAAATATATGCTTAAA
  30556_thamno               TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  30686_cyathophylla         TAAAAGCGAGTCACATCTAATGATCTANAATCTGTGGTATTGTGAAATATATGCTTAAA
  33413_thamno               TAAAAGCAAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  35236_rex                  TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  35855_rex                  TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  38362_rex                  TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTCAAA
  39618_rex                  TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTCAAA
  40578_rex                  TAAAAGCGAGTCACATCTAATGATCTAAAMTCTGTGGTATTGTGAAATATATGCTTAAA
  41478_cyathophylloides     TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  41954_cyathophylloides     TAAAAGCGAGTCACATCTAATGATCTAAAATCTGTGGTATTGTGAAATATATGCTTAAA
  //                                -                     -     -                   *   |3|
  29154_superba              AATGGGTTGTTCCATGGATAACAACTCCGTTTTATRCCAAATACTGTGACACGCACRCA
  32082_przewalskii          AATGGGTTGTTCCATGGTTAACAACTCCGTTTTATGCCAACTACTGCGACACACACGCA
  33588_przewalskii          AATGGGTTGTTCCATGGTTAACAACTCCGTTTTATGCCAACTACTGCGACACGCACGCA
  41478_cyathophylloides     AATGGGTTGTTCCATGGATAACAACTCCGTTTTATGCCAAATACTGTGACACGCACGCA
  41954_cyathophylloides     AATGGGTTGTTCCATGGATAACAACTCCGTTTTATGCCAAATACTGTGACACGCACGCA
  //                                          *                 -    *     *     -   -  |5|
  29154_superba              AGCCGATTCGGTCGCGAGCAGCGATATTTTGTTTCCCCTCAAAATCTTCACAATCTCTA
  30686_cyathophylla         AGCCGATTTGGTTGCGAGCAGCGATATTTTGTTTCCCCTCAAAATCTTCACAATCTCCG
  35236_rex                  AGCYGATTTGGTCGCGAGCAGCGATGTTTTGTTTCCCCTCAAAATCTTCATAATCTCTA
  38362_rex                  AGCYGATTTGGTYGCGAGCAGCGATRTTTTGYTTCCCCTCAAAATCTTCAYAATCTCYR
  41478_cyathophylloides     AGCCGATTTGGTTGCGAGCAGCGATATTTTGTTTCCCCTCAAAATCTTCACAATCTCCA
  41954_cyathophylloides     AGCCGATTTGGTTGCGAGCAGCGATATTTTGTTTCCCCTCAAAATCTTCACAATCTCCA
  //                            *    -   *            *     -                  *      **|7|
  29154_superba              TCGACGCCATGTATGACTGTTCAAAATATCAAATGTACT-ATTACNACCACCCTTTTTT
  30686_cyathophylla         TCGACGCCATGTATGACTGTTCAAAATATCAAATGTACTAATTACCACCACCCTTTTTT
  38362_rex                  TCGACGCCATGTATGACTGTTCAAAATATCAAATGTACT-ATTACCACCACCCTTTTTT
  40578_rex                  TCGACGCCATNTATGACTGTTCAAAATATCAAACGTACT-ATTACCACCACCCTTTTTT
  //                                                          -                         |14|
  29154_superba              ATCGATCATTTCGCCTCACAGTTGCTGGGTGCAGAAAAANNTCTTCATCTGATTCAGGT
  30556_thamno               ATCGATCATTTCTTCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  30686_cyathophylla         ATCGATCATTTCGCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  32082_przewalskii          ATCGATCATTTCGCCTCACAGTTGCTGGATGCAGAAAAAATTCTTCATCTGATTCAGGT
  33413_thamno               ATCGATCATTTCTNCTCACAGTTGCTGGGTNCAGAAAA---------------------
  33588_przewalskii          ATCGATCATTTCGCCTCACAGTTGCTGGATGCAGAAAAAATTCTTCATCTGATTCAGGT
  35236_rex                  ATCGATCATTTCTCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  35855_rex                  ATCGATCATTTCTCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  38362_rex                  ATCGATCATTTCTCCTCACAGTTGCTGGGTGCAAAAAAAATTCTTCATCTGATTCAGGT
  39618_rex                  ATCGATCATTTCTCCTCACAGTTGCTGGGTGCAAAAAAAATTCTTCATCTGATTCAGGT
  40578_rex                  ATCGATCATTTCTCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  41478_cyathophylloides     ATCGATCATTTCGCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  41954_cyathophylloides     ATCGATCATTTCGCCTCACAGTTGCTGGGTGCAGAAAAAATTCTTCATCTGATTCAGGT
  //                                     *-              *    *                         |16|


peek at the .phy files
~~~~~~~~~~~~~~~~~~~~~~
This is the concatenated sequence file of all loci in the data set. It is typically
used in phylogenetic analyses, like in the program *raxml*. This super matrix is
13 taxon deep by 2.44 Mbp long. 


.. code:: bash

  ## cut -c 1-80 prints only the first 80 characters of the file
  cut -c 1-80 pedicularis/base_outfiles/base.phy


.. parsed-literal::

   13 2444589
   29154_superba              TCTGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTTTCGATCT
   30556_thamno               TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCT
   30686_cyathophylla         TCCAGTCCCGCGGGTGATCAAGGCCCCACCACCGCATCTCACATTCTCGATCT
   32082_przewalskii          NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
   33413_thamno               TCCGGTCCTTCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCT
   33588_przewalskii          NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
   35236_rex                  TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTMGATCT
   35855_rex                  TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCT
   38362_rex                  TCCGGTCCTTCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTAGATCT
   39618_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
   40578_rex                  TCCGGTCCYKCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTCTCGATCT
   41478_cyathophylloides     TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTATCGATCT
   41954_cyathophylloides     TCCGGTCCCGCGGGTGATCAAGGCCCCACCACCGCGTCTCACATTATCGATCT


peek at the .snps.phy file
~~~~~~~~~~~~~~~~~~~~~
This is similar to the phylip file format, but only variable site columns are 
included. All SNPs are in the file, in contrast to the .u.snps.phy file, which 
randomly selects only a single SNP per locus. 


.. code:: bash

    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/base_outfiles/base.snps.phy


.. parsed-literal::

   13 156957
   29154_superba              TGCGGTCGAATARATGRCCCATCTATGCGGCATGCACACTMGCTCGGTACAAT
   30556_thamno               CGCGGCAGAGTNNNNNNNNNNNNNNNTTGGNNNNYGYGCTAGTTTGGTACAAT
   30686_cyathophylla         CACGACCGAGTNNNNNNCTTATCCGTGCGGTATGCGCANNNNNNNNNTANNNN
   32082_przewalskii          NNNNNNNNNNNTGCCAGNNNNNNNNNGCAGNNNNNNNNNNNNNNNNNCTNNNN
   33413_thamno               CGTTGCAAAGTNNNNNNNNNNNNNNNTNGGNNNNYGCGNNNNNNNNNNNNNNN
   33588_przewalskii          NNNNNNNNNNNTGCCGGNNNNNNNNNGCAGNNNNNNNNNNNNNNNNNCTNNNN
   35236_rex                  CGCGGCMGAGTNNNNNNYTCGTTTANTCGGNNNNYGCGCTAACATGGTACAAG
   35855_rex                  CGCGGCAGAGTNNNNNNNNNNNNNNNTCGGNNNNYGCGTCAGCTTAGTACAAT
   38362_rex                  CGTTGCAGAGCNNNNNNYTYRYYYRTTCGANNNNNNNNTTAGCTTGATACGAT
   39618_rex                  NNNNNNNGAGCNNNNNNNNNNNNNNNTCGANNNNNNNNTTAGCTTGATACGTT
   40578_rex                  CGYKGCCGMGTNNNNNNNNNNNNNNCTCGGNNNNYGYGNNNNNNNNNTAMAAT
   41478_cyathophylloides     CGCGGACGAGTAGATGGCTTATCCANGCGGTGATNNNNNNNNNNNNNTACAAT
   41954_cyathophylloides     CGCGGACGAGTAGATGGCTTATCCANGCGGTGATNNNNNNNNNNNNNTACAAT


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
   ##fileDate=2016/10/04
   ##source=ipyrad_v.0.4.1
   ##reference=pseudo-reference (most common base at site)
   ##phasing=unphased
   ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
   ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   ##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	29154_superba	30556_thamno	30686_cyathophylla
   2	1	.	T	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,9,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	2	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	3	.	C	T	13	PASS	NS=10;DP=158	GT:CATG	1/1:0,0,9,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	4	.	G	A	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	1/1:0,9,0,0	
   2	5	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	6	.	T	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,9,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	7	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	8	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	9	.	C	T	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	10	.	G	T	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	11	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	12	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	13	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	14	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	15	.	T	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,8,1	0/0:0,0,15,0	0/0:0,0,9,0	
   2	16	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	17	.	A	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,9,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	18	.	T	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,9,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	19	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	20	.	A	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,9,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	21	.	A	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,9,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	22	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	23	.	G	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,0,0,9	0/0:0,0,0,15	0/0:0,0,0,9	
   2	24	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	25	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	26	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	27	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	28	.	A	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,9,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	29	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	30	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	31	.	A	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:0,9,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	32	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	33	.	C	.	13	PASS	NS=10;DP=158	GT:CATG	0/0:9,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	34	.	G	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,0,0,8	0/0:0,0,0,15	0/0:0,0,0,9	
   2	35	.	C	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:8,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	36	.	G	A	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,0,0,8	0/0:0,0,0,15	1/1:0,9,0,0	
   2	37	.	T	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,0,8,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	38	.	C	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:8,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	39	.	T	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,0,8,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	40	.	C	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:8,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	41	.	A	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,8,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	42	.	C	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:8,0,0,0	0/0:15,0,0,0	0/0:9,0,0,0	
   2	43	.	A	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,8,0,0	0/0:0,15,0,0	0/0:0,9,0,0	
   2	44	.	T	.	13	PASS	NS=10;DP=157	GT:CATG	0/0:0,0,8,0	0/0:0,0,15,0	0/0:0,0,9,0	
   2	45	.	T	.	13	PASS	NS=10;DP=156	GT:CATG	0/0:0,0,8,0	0/0:0,0,15,0	0/0:0,0,8,0	
   2	46	.	C	A,T	13	PASS	NS=10;DP=156	GT:CATG	2/2:0,0,8,0	0/0:15,0,0,0	0/0:8,0,0,0
   2	47	.	T	.	13	PASS	NS=10;DP=156	GT:CATG	0/0:0,0,8,0	0/0:0,0,15,0	0/0:0,0,8,0	
   2	48	.	C	A	13	PASS	NS=10;DP=156	GT:CATG	0/0:8,0,0,0	1/1:0,15,0,0	0/0:8,0,0,0	
   2	49	.	G	.	13	PASS	NS=10;DP=156	GT:CATG	0/0:0,0,0,8	0/0:0,0,0,15	0/0:0,0,0,8	
   2	50	.	A	.	13	PASS	NS=10;DP=156	GT:CATG	0/0:0,8,0,0	0/0:0,15,0,0	0/0:0,8,0,0	



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
    
