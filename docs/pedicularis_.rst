
.. include:: global.rst

.. _pedicularis_cli:


Eaton & Ree (2013) single-end RAD data set
==========================================

Here we demonstrate a *denovo* assembly for an empirical RAD data set to 
give a general idea of the results you might expect to recover. 
This example was run on a 4-core laptop with 8GB RAM, and takes about 2.25 hours
to run completely, showing that you do not need a super computer to assemble 
many data sets. However, using more cores will improve the speed of ipyrad 
approximately linearly, so if you have access to a large cluster go ahead and use it.

We will use the 13 taxa *Pedicularis* data set from Eaton and Ree (2013) 
(`open access link <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739883/pdf/syt032.pdf>`__).
This data set is composed of single-end 75bp reads from a RAD-seq library 
prepared with the PstI enzyme. This data set also serves as an example for several 
of our `analysis cookbooks <http://ipyrad.readthedocs.io/analysis.html>`__ 
that demonstrate methods for analyzing RAD-seq results. So after you finish this
assembly head over there to check out fun ways to analyze the data. 


Download the data set (*Pedicularis*)
---------------------------------------
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
-------------------
Always start by using the ``-n {name}`` argument to create a new named Assembly. 
I'll use the name ``base`` to indicate this is the base assembly from which
we will later create several branches. 

.. code:: bash

    >>> ipyrad -n pedicularis

This will print the message:

.. parsed-literal::
    New file 'params-pedicularis.txt' created in /home/deren/Documents/ipyrad/tests


In this case, the data come to us already demultiplexed so we are going to simply set the 
**sorted\_fastq\_path** to tell ipyrad the location of our data files. You can 
select multiple files at once using regular expressions, in this example we
use an asterisk (`*.gz`) to select all files in the directory ending in *.gz*. We also
set a **project\_dir**, which is useful for grouping all our results into a single 
directory. For this we'll use name the project directory "analysis-ipyrad". 
If this folder doesn't exist then ipyrad will create it. Take note when entering 
the values below into your params file that they properly correspond to 
parameters 1 and 4, respectively.

.. parsed-literal::
    ## Use your text editor to enter the following values:
    ## The wildcard (*) tells ipyrad to select all files ending in .gz
    analysis-ipyrad                   ## [1] [project_dir] ...
    example_empirical_rad/*.gz        ## [4] [sorted_fastq_path] ...  

We'll add a few additional options as well to: filter for adapters (param 16); 
trim the 3' edge of R1 aligned loci by 5bp (param 26; this is optional, but helps 
to remove poorly aligned 3' edges); and produce all output formats (param 27).

.. parsed-literal::
    ## enter the following params as well
    2                                ## [16] [filter_adapters] ...
    0, 5, 0, 0                       ## [26] [trim_loci] ...
    *                                ## [27] [output_formats] ...

We'll leave the remaining parameters at their default values. 


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
    ipyrad [v.0.5.15]
    Interactive assembly and analysis of RAD-seq data
   -------------------------------------------------------------
    New Assembly: pedicularis
    host compute node: [20 cores] on tinus

    Step 1: Loading sorted fastq data to Samples
    [####################] 100%  loading reads         | 0:00:11  
    13 fastq files loaded to 13 Samples.


  Summary stats of Assembly pedicularis
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
The example here was run on a 20-core workstation and can finish in ~20 minutes.


.. code:: bash

    ## run steps 2-7
    >>> ipyrad -p params-base.txt -s 234567


.. parsed-literal::
   -------------------------------------------------------------
    ipyrad [v.0.5.15]
    Interactive assembly and analysis of RAD-seq data
   -------------------------------------------------------------
    loading Assembly: pedicularis
    from saved path: ~/Documents/ipyrad/tests/analysis-ipyrad/pedicularis.json
    host compute node: [20 cores] on tinus

    Step 2: Filtering reads 
    [####################] 100%  processing reads      | 0:01:21  

    Step 3: Clustering/Mapping reads
    [####################] 100%  dereplicating         | 0:00:09  
    [####################] 100%  clustering            | 0:05:02  
    [####################] 100%  building clusters     | 0:00:30  
    [####################] 100%  chunking              | 0:00:05  
    [####################] 100%  aligning              | 0:03:27  
    [####################] 100%  concatenating         | 0:00:17  

    Step 4: Joint estimation of error rate and heterozygosity
    [####################] 100%  inferring [H, E]      | 0:01:17  

    Step 5: Consensus base calling 
    Mean error  [0.00283 sd=0.00081]
    Mean hetero [0.01563 sd=0.00238]
    [####################] 100%  calculating depths    | 0:00:05  
    [####################] 100%  chunking clusters     | 0:00:07  
    [####################] 100%  consens calling       | 0:03:12  

    Step 6: Clustering at 0.85 similarity across 13 samples
    [####################] 100%  concat/shuffle input  | 0:00:06  
    [####################] 100%  clustering across     | 0:03:16  
    [####################] 100%  building clusters     | 0:00:06  
    [####################] 100%  aligning clusters     | 0:01:14  
    [####################] 100%  database indels       | 0:00:15  
    [####################] 100%  indexing clusters     | 0:00:09  
    [####################] 100%  building database     | 0:00:30  

    Step 7: Filter and write output files for 13 Samples
    [####################] 100%  filtering loci        | 0:00:06  
    [####################] 100%  building loci/stats   | 0:00:01  
    [####################] 100%  building vcf file     | 0:00:08  
    [####################] 100%  writing vcf file      | 0:00:00  
    [####################] 100%  building arrays       | 0:00:04  
    [####################] 100%  writing outfiles      | 0:01:48  
    Outfiles written to: ~/Documents/ipyrad/tests/analysis-ipyrad/pedicularis_outfiles


  Summary stats of Assembly pedicularis
  ------------------------------------------------
                          state  reads_raw  reads_passed_filter  clusters_total  \
  29154_superba               6     696994               689996          130735   
  30556_thamno                6    1452316              1440314          199587   
  30686_cyathophylla          6    1253109              1206947          233183   
  32082_przewalskii           6     964244               955480          146566   
  33413_thamno                6     636625               626084          169514   
  33588_przewalskii           6    1002923               993873          153089   
  35236_rex                   6    1803858              1787366          410136   
  35855_rex                   6    1409843              1397068          169357   
  38362_rex                   6    1391175              1379626          128389   
  39618_rex                   6     822263               813990          142844   
  40578_rex                   6    1707942              1695523          215721   
  41478_cyathophylloides      6    2199740              2185364          166229   
  41954_cyathophylloides      6    2199613              2176210          293120   

                          clusters_hidepth  hetero_est  error_est  reads_consens  
  29154_superba                      34539    0.015084   0.002612          32913  
  30556_thamno                       51736    0.016421   0.003716          48957  
  30686_cyathophylla                 53357    0.014842   0.003001          50649  
  32082_przewalskii                  41518    0.018446   0.002874          39315  
  33413_thamno                       30913    0.017537   0.002662          29417  
  33588_przewalskii                  45282    0.018394   0.002772          42987  
  35236_rex                          53678    0.015655   0.001939          51485  
  35855_rex                          55421    0.019357   0.003986          52107  
  38362_rex                          51863    0.012369   0.002065          49989  
  39618_rex                          43044    0.014691   0.002916          41122  
  40578_rex                          55350    0.015747   0.002098          53177  
  41478_cyathophylloides             53965    0.012430   0.001714          51816  
  41954_cyathophylloides             73857    0.012264   0.004415          70662  


  Full stats files
  ------------------------------------------------
  step 1: ./analysis-ipyrad/pedicularis_s1_demultiplex_stats.txt
  step 2: ./analysis-ipyrad/pedicularis_edits/s2_rawedit_stats.txt
  step 3: ./analysis-ipyrad/pedicularis_clust_0.85/s3_cluster_stats.txt
  step 4: ./analysis-ipyrad/pedicularis_clust_0.85/s4_joint_estimate.txt
  step 5: ./analysis-ipyrad/pedicularis_consens/s5_consens_stats.txt
  step 6: ./analysis-ipyrad/pedicularis_consens/s6_cluster_stats.txt
  step 7: ./analysis-ipyrad/pedicularis_outfiles/pedicularis_stats.txt



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
  total_prefiltered_loci              88341              0          88341
  filtered_by_rm_duplicates            2566           2566          85775
  filtered_by_max_indels                518            518          85257
  filtered_by_max_snps                  212            121          85136
  filtered_by_max_shared_het            946            908          84228
  filtered_by_min_sample              39170          38942          45286
  filtered_by_max_alleles             10196           5101          40185
  total_filtered_loci                 40185              0          40185


  ## The number of loci recovered for each Sample.
  ## ipyrad API location: [assembly].stats_dfs.s7_samples

                          sample_coverage
  29154_superba                     20755
  30556_thamno                      30996
  30686_cyathophylla                26288
  32082_przewalskii                 14496
  33413_thamno                      18214
  33588_przewalskii                 16846
  35236_rex                         32353
  35855_rex                         32397
  38362_rex                         32795
  39618_rex                         27194
  40578_rex                         33154
  41478_cyathophylloides            30667
  41954_cyathophylloides            27961


  ## The number of loci for which N taxa have data.
  ## ipyrad API location: [assembly].stats_dfs.s7_loci

      locus_coverage  sum_coverage
  1                0             0
  2                0             0
  3                0             0
  4             5136          5136
  5             3702          8838
  6             3311         12149
  7             2942         15091
  8             3028         18119
  9             4014         22133
  10            4904         27037
  11            5486         32523
  12            4740         37263
  13            2922         40185


  ## The distribution of SNPs (var and pis) across loci.
  ## var = all variable sites (pis + autapomorphies)
  ## pis = parsimony informative site (minor allele in >1 sample)
  ## ipyrad API location: [assembly].stats_dfs.s7_snps

       var  sum_var    pis  sum_pis
  0   2107        0  10483        0
  1   3878     3878   9695     9695
  2   5048    13974   7088    23871
  3   5365    30069   4765    38166
  4   4921    49753   3084    50502
  5   4330    71403   1960    60302
  6   3532    92595   1260    67862
  7   2975   113420    819    73595
  8   2253   131444    489    77507
  9   1743   147131    270    79937
  10  1331   160441    144    81377
  11   948   170869     73    82180
  12   665   178849     26    82492
  13   388   183893     19    82739
  14   271   187687      9    82865
  15   178   190357      0    82865
  16   111   192133      1    82881
  17    65   193238      0    82881
  18    39   193940      0    82881
  19    27   194453      0    82881
  20    10   194653      0    82881



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
  13 2577585
  29154_superba              AATGATGGTGGTACACATATTAATTACAATTTGGACAACGGCGGCTTTGTTCA
  30556_thamno               ACAGATGGTGGTACACATGTCAATTACAATTTGGATAACGGCGGNNNNNNNNN
  30686_cyathophylla         AATGATGGTGGTACACATATTAATTACAATTTGGACAACGGCGGCTTTGTTCA
  32082_przewalskii          NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  33413_thamno               AGTGATGGTGGTACACATGTCNANTACAATTTGGACAACGGCGGCTTTGTTCN
  33588_przewalskii          NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  35236_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  35855_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  38362_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  39618_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  40578_rex                  AATGATGGTGGTACACATATYAATTACAAYTTGGAYAACGGCGGCTTTGTTCA
  41478_cyathophylloides     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  41954_cyathophylloides     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN


peek at the .snps.phy file
~~~~~~~~~~~~~~~~~~~~~
This is similar to the phylip file format, but only variable site columns are 
included. All SNPs are in the file, in contrast to the .u.snps.phy file, which 
randomly selects only a single SNP per locus. 


.. code:: bash

    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/base_outfiles/base.snps.phy


.. parsed-literal::
  13 194653
  29154_superba              ATATTCAAACTATTCAAAGTAACTGATGAAAYCTAGGGGAKCAGTTCGCGTGC
  30556_thamno               CAGCTTAAATTATNNGGCGCAACCGGAGAAANNNNNNNNNGAAGGTTACATNN
  30686_cyathophylla         ATATTCAAACTATAANNNNNAACTGATGAAACTTGTCGGNGCAGGTTACATGC
  32082_przewalskii          NNNNNNNNNNNNNACNNNNNCANNNNNNNNNNNNNNNNNCNNNNGTCGCGTNN
  33413_thamno               GTGCTCAAATTAANNNNNNNAGCCGAAGAAACCCGGCATNGCAKGTTANANNN
  33588_przewalskii          NNNNNNNNNNNNNACNNNNNAANNNNNNNNNNNNNNNNNCNNNNGTCGCGYNN
  35236_rex                  NNNNNNAAATTGTNNGGCGTAACCAAAGAAANNNNNNNNNGCAGGTTAAATNN
  35855_rex                  NNNNNNACATTATNNNNNNNAATCGAAGAAANNNNNNNNNNNNNGTTACATGC
  38362_rex                  NNNNNNGAATTATNNNNNNNAACCAAAGAAACCCG-CCGNNNNNGTTACATNN
  39618_rex                  NNNNNNGAATTATNNNNNNNAACCAAAGAAACCCG-CCGNNNNNGKTACATNN
  40578_rex                  ATAYYYRAATYATNNGGCKTAACCGAAGAGGNNNNNNNNNNNNNGTTACATRC
  41478_cyathophylloides     NNNNNNATTTTATACNNNNNAACTGATTGAACCTAGGGGAGCGGGTTACATGT
  41954_cyathophylloides     NNNNNNATTTTATACNNNNNAANNNNNNNNNCCTAGGGGAGCGGGTTACATGT



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

    head -n 50 pedicularis/base_outfiles/base.vcf | cut -c 1-80


.. parsed-literal::
  ##fileformat=VCFv4.0
  ##fileDate=2017/02/14
  ##source=ipyrad_v.0.5.15
  ##reference=pseudo-reference (most common base at site)
  ##phasing=unphased
  ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Dat
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  ##FORMAT=<ID=CATG,Number=1,Type=String,Description="Base Counts (CATG)">
  #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  29154_superba 30556_thamno  30
  locus_1 2 . A C,G 13  PASS  NS=5;DP=49  GT:DP:CATG  0/0:9:0,9,0,0 1/1:7:7,0,0,0
  locus_1 3 . T A 13  PASS  NS=5;DP=49  GT:DP:CATG  0/0:9:0,0,9,0 1/1:7:0,7,0,0 0
  locus_1 19  . A G 13  PASS  NS=5;DP=49  GT:DP:CATG  0/0:9:0,9,0,0 1/1:7:0,0,0,7 
  locus_1 21  . C T 13  PASS  NS=5;DP=49  GT:DP:CATG  1/1:9:0,0,9,0 0/0:7:7,0,0,0 
  locus_1 30  . T C 13  PASS  NS=5;DP=49  GT:DP:CATG  0/0:9:0,0,9,0 0/0:7:0,0,7,0 
  locus_1 36  . C T 13  PASS  NS=5;DP=49  GT:DP:CATG  0/0:9:9,0,0,0 1/1:7:0,0,7,0 
  locus_2 15  . A G 13  PASS  NS=11;DP=210  GT:DP:CATG  0/0:12:0,12,0,0 0/0:24:0,2
  locus_2 16  . A T,C 13  PASS  NS=11;DP=210  GT:DP:CATG  0/0:12:0,12,0,0 0/0:24:0
  locus_2 18  . A T 13  PASS  NS=11;DP=210  GT:DP:CATG  0/0:12:0,12,0,0 0/0:24:0,2
  locus_2 20  . T C 13  PASS  NS=11;DP=210  GT:DP:CATG  1/1:12:12,0,0,0 0/0:24:0,0
  locus_2 29  . T C 13  PASS  NS=11;DP=209  GT:DP:CATG  0/0:12:0,0,12,0 0/0:23:0,0
  locus_2 30  . A G 13  PASS  NS=11;DP=209  GT:DP:CATG  0/0:12:0,12,0,0 0/0:23:0,2
  locus_2 47  . T A 13  PASS  NS=11;DP=210  GT:DP:CATG  0/0:12:0,0,12,0 0/0:24:0,0
  locus_3 46  . A T 13  PASS  NS=6;DP=69  GT:DP:CATG  1/1:10:0,0,10,0 ./.:0:0,0,0,
  locus_3 62  . C A 13  PASS  NS=6;DP=68  GT:DP:CATG  0/0:10:10,0,0,0 ./.:0:0,0,0,
  locus_6 11  . G A 13  PASS  NS=4;DP=67  GT:DP:CATG  1/1:11:0,11,0,0 0/0:7:0,0,0,
  locus_6 29  . G A 13  PASS  NS=4;DP=67  GT:DP:CATG  1/1:11:0,11,0,0 0/0:7:0,0,0,
  locus_6 34  . C A 13  PASS  NS=4;DP=67  GT:DP:CATG  1/1:11:0,11,0,0 0/0:7:7,0,0,
  locus_6 35  . G T 13  PASS  NS=4;DP=67  GT:DP:CATG  0/0:11:0,0,0,11 0/0:7:0,0,0,
  locus_6 40  . T C 13  PASS  NS=4;DP=67  GT:DP:CATG  0/0:11:0,0,11,0 1/1:7:7,0,0,
  locus_9 19  . A C 13  PASS  NS=13;DP=224  GT:DP:CATG  0/0:12:0,12,0,0 0/0:22:0,2
  locus_9 25  . A G 13  PASS  NS=13;DP=224  GT:DP:CATG  0/0:12:0,12,0,0 0/0:22:0,2
  locus_11  4 . C T 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:11,0,0,0 0/0:17:17,
  locus_11  13  . C T 13  PASS  NS=10;DP=137  GT:DP:CATG  1/1:11:0,1,10,0 0/0:17:17
  locus_11  21  . G A 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,0,0,11 0/0:17:0,
  locus_11  23  . A G 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,11,0,0 1/1:17:0,
  locus_11  24  . A T 13  PASS  NS=10;DP=137  GT:DP:CATG  1/1:11:0,0,11,0 0/0:17:0,
  locus_11  38  . G T 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,0,0,11 0/0:17:0,
  locus_11  42  . A G 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,11,0,0 0/0:17:0,
  locus_11  54  . A G 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,11,0,0 0/0:17:0,
  locus_11  55  . A G 13  PASS  NS=10;DP=137  GT:DP:CATG  0/0:11:0,11,0,0 0/0:17:0,
  locus_12  6 . C T 13  PASS  NS=7;DP=94  GT:DP:CATG  1/0:7:4,0,3,0 ./.:0:0,0,0,0 
  locus_12  20  . C T 13  PASS  NS=7;DP=94  GT:DP:CATG  0/0:7:7,0,0,0 ./.:0:0,0,0,0
  locus_12  31  . T C 13  PASS  NS=7;DP=94  GT:DP:CATG  0/0:7:0,0,7,0 ./.:0:0,0,0,0
  locus_12  33  . G A 13  PASS  NS=7;DP=94  GT:DP:CATG  1/1:7:0,7,0,0 ./.:0:0,0,0,0
  locus_12  37  . G T 13  PASS  NS=5;DP=52  GT:DP:CATG  0/0:7:0,0,0,7 ./.:0:0,0,0,0
  locus_12  43  . C G 13  PASS  NS=7;DP=94  GT:DP:CATG  1/1:7:0,0,0,7 ./.:0:0,0,0,0
  locus_12  45  . G C,A 13  PASS  NS=7;DP=94  GT:DP:CATG  0/0:7:0,0,0,7 ./.:0:0,0,0
  locus_12  47  . G T 13  PASS  NS=7;DP=94  GT:DP:CATG  0/0:7:0,0,0,7 ./.:0:0,0,0,0




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
    
