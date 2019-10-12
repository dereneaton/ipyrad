
.. include:: global.rst

.. _pedicularis_cli:


Eaton & Ree (2013) single-end RAD data set
==========================================

Here we demonstrate a *denovo* assembly for an empirical RAD data set to 
give a general idea of the results you might expect to recover. 
This example was run on a 8-core laptop with 16GB RAM, and takes 
about 1 hour to run completely. 

We will use the 13 taxa *Pedicularis* data set from Eaton and Ree (2013) 
(`open access link <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739883/pdf/syt032.pdf>`__).
This data set is composed of single-end 75bp reads from a RAD-seq library 
prepared with the PstI enzyme. This data set also serves as an example for several 
of our `analysis tools <https://ipyrad.readthedocs.io/en/latest/API-analysis/>`__ 
to demonstrate methods for analyzing RAD-seq results. So after you finish this
assembly head over there to check out fun ways to analyze the data. 


Download the data set (*Pedicularis*)
---------------------------------------
These data are archived on the NCBI sequence read archive (SRA) under 
accession id SRP021469. We've written a convenient wrapper for sra-tools 
that allows ipyrad to download data from SRA, SRP, ERA, etc., 
IDs easily (`more info here <https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-sratools.html>`__). Run the code below to download and decompress the fastq data files, 
which will save them into a directory called ``example_empirical_data/``, 
or whatever you wish to name it. The directory will be created if it doesn't
already exist. The compressed file size is approximately 1.1GB.

.. code:: bash

  # first we need to download an additional tool
  >>> conda install sra-tools -c bioconda

  # then, use ipyrad to download the fastq data from the SRA database
  >>> ipyrad --download SRP021469 example_empirical_data/


The `--download` function will print the following output:

.. code:: parsed-literal

  Fetching project data...
             Run    spots  mates                ScientificName              SampleName
  0   SRR1754715   696994      0           Pedicularis superba           29154_superba
  1   SRR1754720  1452316      0       Pedicularis thamnophila            30556_thamno
  2   SRR1754730  1253109      0      Pedicularis cyathophylla      30686_cyathophylla
  3   SRR1754729   964244      0       Pedicularis przewalskii       32082_przewalskii
  4   SRR1754728   636625      0       Pedicularis thamnophila            33413_thamno
  5   SRR1754727  1002923      0       Pedicularis przewalskii       33588_przewalskii
  6   SRR1754731  1803858      0               Pedicularis rex               35236_rex
  7   SRR1754726  1409843      0               Pedicularis rex               35855_rex
  8   SRR1754725  1391175      0               Pedicularis rex               38362_rex
  9   SRR1754723   822263      0               Pedicularis rex               39618_rex
  10  SRR1754724  1707942      0               Pedicularis rex               40578_rex
  11  SRR1754722  2199740      0  Pedicularis cyathophylloides  41478_cyathophylloides
  12  SRR1754721  2199613      0  Pedicularis cyathophylloides  41954_cyathophylloides
  Parallel connection | latituba: 8 cores
  [####################] 100% 0:01:43 | downloading/extracting fastq data 

  13 fastq files downloaded to /home/deren/Documents/ipyrad/sandbox/pedicularis/example_empirical_data



Setup a params file
-------------------
Always start an ipyrad assembly by using the ``-n {name}`` argument to 
create a new named Assembly. I'll use the name ``pedicularis`` to 
indicate taxa being assembled.

.. code:: bash

  >>> ipyrad -n pedicularis

This will print the message:

.. code:: parsed-literal

  New file 'params-pedicularis.txt' created in /home/deren/Documents/ipyrad/sandbox


In this example, the data come to us already demultiplexed so we are going to 
simply set the **sorted\_fastq\_path** to tell ipyrad the location of our data 
files. Published data sets will typically be available in this way, already 
demultiplexed. You can select multiple files at once using regular expressions, 
in this example we use an asterisk (`*.gz`) to select all files in the 
directory ending in *.gz*. We also set a **project\_dir**, which is useful for
grouping all our results into a single directory. For this we'll use name 
the project directory "analysis-ipyrad". If this folder doesn't exist then 
ipyrad will create it. Take note when entering the values below into your 
params file that they properly correspond to parameters 1 and 4, respectively.
Use any text editor to edit the params file.

.. code:: parsed-literal
    # Use your text editor to enter the following values:
    # The wildcard (*) tells ipyrad to select all files ending in .gz
    analysis-ipyrad                    ## [1] [project_dir] ...
    example_empirical_data/*.fastq     ## [4] [sorted_fastq_path] ...  

We'll add a few additional options as well to: filter for adapters (param 16); 
trim the 3' edge of R1 aligned loci by 5bp (param 26; this is optional, but helps 
to remove poorly aligned 3' edges); and produce all output formats (param 27).

.. code:: parsed-literal
    # enter the following params as well
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

    >>> ipyrad -p params-pedicularis.txt -s 1 -r


.. code:: parsed-literal

 -------------------------------------------------------------
  ipyrad [v.0.9.14]
  Interactive assembly and analysis of RAD-seq data
 ------------------------------------------------------------- 
  Parallel connection | latituba: 8 cores
  
  Step 1: Loading sorted fastq data to Samples
  [####################] 100% 0:00:06 | loading reads          
  13 fastq files loaded to 13 Samples.

  Parallel connection closed.

  Summary stats of Assembly pedicularis
  ------------------------------------------------
                                     state  reads_raw
  29154_superba_SRR1754715               1     696994
  30556_thamno_SRR1754720                1    1452316
  30686_cyathophylla_SRR1754730          1    1253109
  32082_przewalskii_SRR1754729           1     964244
  33413_thamno_SRR1754728                1     636625
  33588_przewalskii_SRR1754727           1    1002923
  35236_rex_SRR1754731                   1    1803858
  35855_rex_SRR1754726                   1    1409843
  38362_rex_SRR1754725                   1    1391175
  39618_rex_SRR1754723                   1     822263
  40578_rex_SRR1754724                   1    1707942
  41478_cyathophylloides_SRR1754722      1    2199740
  41954_cyathophylloides_SRR1754721      1    2199613


  Full stats files
  ------------------------------------------------
  step 1: ./analysis-ipyrad/pedicularis_s1_demultiplex_stats.txt
  step 2: None
  step 3: None
  step 4: None
  step 5: None
  step 6: None
  step 7: None


Run the remaining assembly steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Because the point of this tutorial is to demonstrate run times and 
statistics, I will leave the rest of the parameters at their
defaults and simply run all remaining steps. Further below I will 
explain in more detail the stats files for each step and the meaning
of the stats values. 

.. code:: bash

    ## run steps 2-7
    >>> ipyrad -p params-pedicularis.txt -s 234567 -r


.. code:: parsed-literal

 -------------------------------------------------------------
  ipyrad [v.0.9.14]
  Interactive assembly and analysis of RAD-seq data
 ------------------------------------------------------------- 
  Parallel connection | latituba: 8 cores
  
  Step 2: Filtering and trimming reads
  [####################] 100% 0:01:59 | processing reads     
  
  Step 3: Clustering/Mapping reads within samples
  [####################] 100% 0:00:12 | dereplicating          
  [####################] 100% 0:07:34 | clustering/mapping     
  [####################] 100% 0:00:01 | building clusters      
  [####################] 100% 0:00:00 | chunking clusters      
  [####################] 100% 0:12:11 | aligning clusters      
  [####################] 100% 0:00:04 | concat clusters        
  [####################] 100% 0:00:03 | calc cluster stats     
  
  Step 4: Joint estimation of error rate and heterozygosity
  [####################] 100% 0:00:31 | inferring [H, E]       
  
  Step 5: Consensus base/allele calling 
  Mean error  [0.00314 sd=0.00090]
  Mean hetero [0.02171 sd=0.00385]
  [####################] 100% 0:00:03 | calculating depths     
  [####################] 100% 0:00:05 | chunking clusters      
  [####################] 100% 0:09:12 | consens calling        
  [####################] 100% 0:00:06 | indexing alleles       
  
  Step 6: Clustering/Mapping across samples 
  [####################] 100% 0:00:03 | concatenating inputs   
  [####################] 100% 0:02:06 | clustering across    
  [####################] 100% 0:00:02 | building clusters      
  [####################] 100% 0:01:31 | aligning clusters      
  
  Step 7: Filtering and formatting output files 
  [####################] 100% 0:00:14 | applying filters       
  [####################] 100% 0:00:09 | building arrays        
  [####################] 100% 0:00:10 | writing conversions    
  [####################] 100% 0:01:27 | indexing vcf depths    
  [####################] 100% 0:00:24 | writing vcf output     

  Parallel connection closed.

  Summary stats of Assembly pedicularis
  ------------------------------------------------
                                     state  reads_raw  ...  error_est  reads_consens
  29154_superba_SRR1754715               6     696994  ...   0.003211          29903
  30556_thamno_SRR1754720                6    1452316  ...   0.003184          43870
  30686_cyathophylla_SRR1754730          6    1253109  ...   0.003297          45856
  32082_przewalskii_SRR1754729           6     964244  ...   0.003079          34733
  33413_thamno_SRR1754728                6     636625  ...   0.003317          26228
  33588_przewalskii_SRR1754727           6    1002923  ...   0.003267          38137
  35236_rex_SRR1754731                   6    1803858  ...   0.002206          46683
  35855_rex_SRR1754726                   6    1409843  ...   0.004316          46234
  38362_rex_SRR1754725                   6    1391175  ...   0.002350          46081
  39618_rex_SRR1754723                   6     822263  ...   0.003636          37259
  40578_rex_SRR1754724                   6    1707942  ...   0.002229          48255
  41478_cyathophylloides_SRR1754722      6    2199740  ...   0.001721          47976
  41954_cyathophylloides_SRR1754721      6    2199613  ...   0.005028          64654

  [13 rows x 8 columns]


  Full stats files
  ------------------------------------------------
  step 1: ./analysis-ipyrad/pedicularis_s1_demultiplex_stats.txt
  step 2: ./analysis-ipyrad/pedicularis_edits/s2_rawedit_stats.txt
  step 3: ./analysis-ipyrad/pedicularis_clust_0.85/s3_cluster_stats.txt
  step 4: ./analysis-ipyrad/pedicularis_clust_0.85/s4_joint_estimate.txt
  step 5: ./analysis-ipyrad/pedicularis_consens/s5_consens_stats.txt
  step 6: ./analysis-ipyrad/pedicularis_across/pedicularis_clust_database.fa
  step 7: ./analysis-ipyrad/pedicularis_outfiles/pedicularis_stats.txt



Take a look at the stats summary 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each assembly that finishes step 7 will create a stats.txt output summary
in the 'assembly_name'_outfiles/ directory. This includes information about 
which filters removed data from the assembly, how many loci were recovered
per sample, how many samples had data for each locus, and how many variable
sites are in the assembled data. 

.. code:: bash

    >>> cat ./analysis-ipyrad/pedicularis_outfiles/pedicularis_stats.txt

.. code:: parsed-literal

  ## The number of loci caught by each filter.
  ## ipyrad API location: [assembly].stats_dfs.s7_filters

                              total_filters  applied_order  retained_loci
  total_prefiltered_loci                  0              0          80481
  filtered_by_rm_duplicates             828            828          79653
  filtered_by_max_indels               1290           1290          78363
  filtered_by_max_SNPs                  946            914          77449
  filtered_by_max_shared_het            718            699          76750
  filtered_by_min_sample              35889          35672          41078
  total_filtered_loci                 39671          39403          41078


  ## The number of loci recovered for each Sample.
  ## ipyrad API location: [assembly].stats_dfs.s7_samples

                                     sample_coverage
  29154_superba_SRR1754715                     21095
  30556_thamno_SRR1754720                      31418
  30686_cyathophylla_SRR1754730                26754
  32082_przewalskii_SRR1754729                 14507
  33413_thamno_SRR1754728                      18504
  33588_przewalskii_SRR1754727                 16928
  35236_rex_SRR1754731                         32588
  35855_rex_SRR1754726                         32462
  38362_rex_SRR1754725                         33372
  39618_rex_SRR1754723                         27673
  40578_rex_SRR1754724                         33260
  41478_cyathophylloides_SRR1754722            31255
  41954_cyathophylloides_SRR1754721            28381


  ## The number of loci for which N taxa have data.
  ## ipyrad API location: [assembly].stats_dfs.s7_loci

      locus_coverage  sum_coverage
  1                0             0
  2                0             0
  3                0             0
  4             5347          5347
  5             3809          9156
  6             3488         12644
  7             3096         15740
  8             3330         19070
  9             4217         23287
  10            5102         28389
  11            5422         33811
  12            4562         38373
  13            2705         41078


  The distribution of SNPs (var and pis) per locus.
  ## var = Number of loci with n variable sites (pis + autapomorphies)
  ## pis = Number of loci with n parsimony informative site (minor allele in >1 sample)
  ## ipyrad API location: [assembly].stats_dfs.s7_snps
  ## The "reference" sample is included if present unless 'exclude_reference=True'

       var  sum_var    pis  sum_pis
  0   1806        0  10232        0
  1   3528     3528   9935     9935
  2   4758    13044   7601    25137
  3   5297    28935   5109    40464
  4   5183    49667   3212    53312
  5   4650    72917   2029    63457
  6   3893    96275   1225    70807
  7   3340   119655    757    76106
  8   2529   139887    456    79754
  9   1967   157590    302    82472
  10  1480   172390    132    83792
  11  1145   184985     68    84540
  12   796   194537     16    84732
  13   541   201570      4    84784
  14   151   203684      0    84784
  15    13   203879      0    84784
  16     0   203879      0    84784
  17     0   203879      0    84784
  18     0   203879      0    84784
  19     1   203898      0    84784


  ## Final Sample stats summary
                                     state  reads_raw  reads_passed_filter  clusters_total  clusters_hidepth  hetero_est  error_est  reads_consens  loci_in_assembly
  29154_superba_SRR1754715               7     696994               689996          126896             34145    0.024641   0.003211          29903             21095
  30556_thamno_SRR1754720                7    1452316              1440314          192920             50491    0.022635   0.003184          43870             31418
  30686_cyathophylla_SRR1754730          7    1253109              1206947          225144             52464    0.020622   0.003297          45856             26754
  32082_przewalskii_SRR1754729           7     964244               955480          142366             41046    0.027211   0.003079          34733             14507
  33413_thamno_SRR1754728                7     636625               626084          165338             30754    0.024820   0.003317          26228             18504
  33588_przewalskii_SRR1754727           7    1002923               993873          148920             44642    0.025917   0.003267          38137             16928
  35236_rex_SRR1754731                   7    1803858              1787366          401906             52694    0.019709   0.002206          46683             32588
  35855_rex_SRR1754726                   7    1409843              1397068          164312             54484    0.025071   0.004316          46234             32462
  38362_rex_SRR1754725                   7    1391175              1379626          124417             51061    0.016379   0.002350          46081             33372
  39618_rex_SRR1754723                   7     822263               813990          138973             42451    0.022817   0.003636          37259             27673
  40578_rex_SRR1754724                   7    1707942              1695523          210842             54539    0.019760   0.002229          48255             33260
  41478_cyathophylloides_SRR1754722      7    2199740              2185364          162093             53191    0.015180   0.001721          47976             31255
  41954_cyathophylloides_SRR1754721      7    2199613              2176210          286667             72791    0.017415   0.005028          64654             28381


  ## Alignment matrix statistics:
  snps matrix size: (13, 203898), 35.26% missing sites.
  sequence matrix size: (13, 2840602), 35.60% missing sites.


Take a peek at the .loci output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the first place I look when an assembly finishes. It provides a 
clean view of the data with variable sites (-) and parsimony informative
SNPs (*) highlighted. Use the unix commands **less** or **head** to look at this
file briefly. Each locus is labelled with a number corresponding to the locus 
order before filters are applied in step 7. If you branch this assembly and 
run step 7 again with a different set of parameters you may recover fewer 
or more total loci. 


.. code:: bash

  ## head -n 50 prints just the first 50 lines of the file to stdout
  >>> head -n 50 analysis-ipyrad/pedicularis_outfiles/pedicularis.loci

.. code:: bash 

  29154_superba_SRR1754715              TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGGGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  30556_thamno_SRR1754720               TAGGGTGGGTCKCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  30686_cyathophylla_SRR1754730         TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGGGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  32082_przewalskii_SRR1754729          TAGGGTGGGTCTCGTTCAAGGTATTCGAACAAGAGGGTACCCTGCGAACTTCCAAATTCACCCTCNTCG
  33413_thamno_SRR1754728               TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  33588_przewalskii_SRR1754727          TAGGGTGGGTCTCNTTCAAGGTATTCGAACAASAGGGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  35236_rex_SRR1754731                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  35855_rex_SRR1754726                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  38362_rex_SRR1754725                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  39618_rex_SRR1754723                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCG
  40578_rex_SRR1754724                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCCTGCGAACTTCCAAATTCACCCTCATCA
  41478_cyathophylloides_SRR1754722     TAGGGTGGGTCTCGTTCAAGGTATTCGAACAATAGGGTACCCAGCGAACTTCCAAATTCACCCTCATCG
  41954_cyathophylloides_SRR1754721     TAGGGTGGGTCTCGTTCAAGGTATTCGAACAATAGGGTACCCAGCGAACTTCCAAATTCACCCTCATCG
  //                                               -                    *  *      *                         -|0|
  29154_superba_SRR1754715              ACGACGTCTCTCCCCGAGCCGGCTATCAGGAGACGGATTTTCGAGATGGGGGGTCGTTTTGCTGTTTGT
  30556_thamno_SRR1754720               ACGACGTCTCTCCCCGAGCCGGCTATCAGGAGACGGATTTTCGAGATGGGGGGTCGTTTTGCTGTTTGT
  33413_thamno_SRR1754728               ACGACGTCTCTCCCCGAGMCGGCTATCAGGAGACGGATTTTCGAGATGGGGGGTCGTTTTGCTGTTTGT
  40578_rex_SRR1754724                  ACGACGTCTCTCCCCGAGCCGGCTATCAGGAGACGGATTTTCGAGATGGGGGGTCGTTTTGCTGTTTGT
  //                                                      -                                                  |1|
  29154_superba_SRR1754715              CTTGGCACTGAATTAGCAGAACTTCAACAATTAAGTCTCCAGTATAATTGAATTYGATTTAATTTAATT
  30686_cyathophylla_SRR1754730         CTTGGCACTGAATTAGCAGAACTTCAACAATTAAGTCTCCAGTATAACTGAATTTGATTTAATTTAATT
  35236_rex_SRR1754731                  CTTGGCACTGAAGTAGCAGAACTTCAACAATTAAGTCTCCAGTATAATTGAATTTGATTTAATTTAATG
  35855_rex_SRR1754726                  CTTGGCACTGAAGTAGCAGAACTTCAACAATTAAGTCTCCGGTATAATTGAATTTGATTTAATTTAATG
  38362_rex_SRR1754725                  CTTGGCACTGAAGTAGCAGAACTTCAACAATTAAGTCTCAGGTATAATTGAATTTGATTTAATTTAATT
  40578_rex_SRR1754724                  CTTGGCACTCAAGTAGCAGAACTTCAACAATTAAGTCTCCGGTATAATTGAATTTGATTTAATTTAATG
  41478_cyathophylloides_SRR1754722     CTTGGCACTGAAGTAGCAGAACTTCAACAATTAAGTCTCCAGTATAATTGAATTTGATTTAATTTAATT
  41954_cyathophylloides_SRR1754721     CTTGGCACTGAAGTAGCAGAACTTCAACAATTAAGTCTCCAGTATAATTGAATTTGATTTAATTTAATT
  //                                             -  *                          -*      -      -             *|2|
  29154_superba_SRR1754715              GATCCTGAAATGACAASAAACATAACANGGGGGTAATTTTTTGTAATTAT---CCCTTAGA-TAAACTATACA
  33413_thamno_SRR1754728               GATCCTGAAACGACAACAAACATAACACGGGGGTAATYTTTTGTAATTAT---CCCTTMGA-TAAACTATACA
  35236_rex_SRR1754731                  GATCCTGAAACGACAACAAACATAACACGGGGGTAATTTTTTGTAATTAT---CCCTTAGA-TAAACTATACA
  38362_rex_SRR1754725                  GATCCTGAAACGACAACAAACATAACACGGGGGTAATTTTTTGTAATTAT---CCCTTAGA-TAAACTATACA
  39618_rex_SRR1754723                  GATCCTGAAACGACAACAAACATAACACGGGGGTAATTTTTTGTAATTAT---CCCTTAGA-TAAACTATACA
  40578_rex_SRR1754724                  GATCCTGAAACGACAACAAACATAACAYGGGGGTAATTTTTTGTAATTAY---CCCTTAGA-TAAACTATACA
  41478_cyathophylloides_SRR1754722     GATCCTGAAATGACAACAAACATAACAGGGGGGTAATTTTTTGTAATTATCCCCCCTTAGATTAAACTA----
  41954_cyathophylloides_SRR1754721     GATCCTGAAATGACAACAAACATAACAGGGGGGTAATTTTTTGTAATTATCCCCCCTTAGATTAAACT-----
  //                                              *     -          *         -           -        -              |3|
  29154_superba_SRR1754715              AAAACAGGATGAGTGCATATCTCTCGTTCTAACTACTGCAATGCTAGGNAAATAAAATACAGACTAAAA
  30686_cyathophylla_SRR1754730         AAAACAGGATGAGTGCATATCTCTCGTTCTAACTACTGCAATGCTAGGTAAATAAAATACAGACTAAAA
  32082_przewalskii_SRR1754729          AAAACAGGATGAGTGCATATCTCTCGTACTAACTACTGCAATGCTAGGTAAATAAAATACAGACTAAAA
  33588_przewalskii_SRR1754727          AAAACAGGATGAGTGCATATCTCTCGTACTAACTACTGCAATGCTAGGTAAATAAAATACAGACTAAAA
  41478_cyathophylloides_SRR1754722     AAAACAGGATGAGTGCATATCTCTCGTTTTAACTACTGCAATGCTAGGTAAATAAAATAGAGACTAAAA
  41954_cyathophylloides_SRR1754721     AAAACAGGATGAGTGCATATCTCTCGTTTTAACTACTGCAATGCTAGGTAAATAAAATAGAGACTAAAA
  //                                                               **                              *         |4|




peek at the .phy files
~~~~~~~~~~~~~~~~~~~~~~
This is the concatenated sequence file of all loci in the data set. It is typically
used in phylogenetic analyses, like in the program *raxml*. This super matrix is
13 taxon deep by 2.44 Mbp long. 


.. code:: bash

  ## cut -c 1-80 prints only the first 80 characters of the file
  >>> cut -c 1-80 analysis-ipyrad/pedicularis_outfiles/pedicularis.phy


.. code:: parsed-literal

  13 2840602
  29154_superba_SRR1754715              TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGGGTACCC
  30556_thamno_SRR1754720               TAGGGTGGGTCKCGTTCAAGGTATTCGAACAACAGAGTACCC
  30686_cyathophylla_SRR1754730         TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGGGTACCC
  32082_przewalskii_SRR1754729          TAGGGTGGGTCTCGTTCAAGGTATTCGAACAAGAGGGTACCC
  33413_thamno_SRR1754728               TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  33588_przewalskii_SRR1754727          TAGGGTGGGTCTCNTTCAAGGTATTCGAACAASAGGGTACCC
  35236_rex_SRR1754731                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  35855_rex_SRR1754726                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  38362_rex_SRR1754725                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  39618_rex_SRR1754723                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  40578_rex_SRR1754724                  TAGGGTGGGTCTCGTTCAAGGTATTCGAACAACAGAGTACCC
  41478_cyathophylloides_SRR1754722     TAGGGTGGGTCTCGTTCAAGGTATTCGAACAATAGGGTACCC
  41954_cyathophylloides_SRR1754721     TAGGGTGGGTCTCGTTCAAGGTATTCGAACAATAGGGTACCC


peek at the .snps.phy file
~~~~~~~~~~~~~~~~~~~~~
This is similar to the phylip file format, but only variable site columns are 
included. All SNPs are in the file, in contrast to the .u.snps.phy file, which 
randomly selects only a single SNP per locus. 


.. code:: bash

    ## cut -c 1-80 prints only the first 80 characters of the file
    >>> cut -c 1-80 analysis-ipyrad/pedicularis_outfiles/pedicularis.snps.phy


.. code:: parsed-literal

  13 203898
  29154_superba_SRR1754715              TCGTGCGTCATYTTSNTTATCCGAYYACTGTGTAAGCCGGGG
  30556_thamno_SRR1754720               KCATGCNNNNNNNNNNNNNNNNGACCGNNGTGTAAGCCGGGG
  30686_cyathophylla_SRR1754730         TCGTGNGTCACTTNNNNNNTCCGACCACTGTGTAAGCCAGGG
  32082_przewalskii_SRR1754729          TGGTGNNNNNNNNNNNNNNACCNNNNNNNNNNNNNNNTGCAA
  33413_thamno_SRR1754728               TCATGMNNNNNNNCCCYTMNNNGACCGNNGTGTAAGCNNNNN
  33588_przewalskii_SRR1754727          TSGTGNNNNNNNNNNNNNNACCATCCGNNACATAAGCTGCAA
  35236_rex_SRR1754731                  TCATGNGGCATTGCCCTTANNNGACCGYTGTGTRRGCNNNNN
  35855_rex_SRR1754726                  TCATGNGGCGTTGNNNNNNNNNGACCGCTGTGTAAGCNNNNN
  38362_rex_SRR1754725                  TCATGNGGAGTTTCCCTTANNNGACCGNNTTGTAAGANNNNN
  39618_rex_SRR1754723                  TCATGNNNNNNNNCCCTTANNNGACCGNNNNNNNNNNNNNNN
  40578_rex_SRR1754724                  TCATACCGCGTTGCCYTYANNNGACCGCWGTGAAARCNNNNN
  41478_cyathophylloides_SRR1754722     TTGAGNGGCATTTTCGTTATTGGACCGCTGTGTAAGCCGGGG
  41954_cyathophylloides_SRR1754721     TTGAGNGGCATTTTCGTTATTGGACCGCTGTGTAAGCCGGGG


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

    >>> head -n 50 analysis-ipyrad/pedicularis_outfiles/pedicularis.vcf | cut -c 1-80


.. parsed-literal::
  ##fileformat=VCFv4.0
  ##fileDate=2017/02/14
  ##source=ipyrad_v.0.7.28
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

