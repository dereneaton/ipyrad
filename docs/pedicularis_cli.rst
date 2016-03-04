
.. include:: global.rst

.. _pedicularis_cli:


Empirical example (*Pedicularis*) - CLI
========================================
For this tutorial we will assemble a single-end RAD-seq data set of
13 individuals from the *Cyathophora* clade of the angiosperm genus 
*Pedicularis*, originally published by **Eaton and Ree (2013)** 
(:ref:`link to open access article 
<http://sysbio.oxfordjournals.org/content/62/5/689.full>`). All of the code 
on this page uses the CLI, and thus should be executed in a terminal. 


Download the fastq files
~~~~~~~~~~~~~~~~~~~~~~~~
The data are hosted online at the NCBI sequence read archive (SRA) under 
accession id SRP021469. For convenience, I've also hosted the data at a 
publicly available dropbox link which we will use to download the data here, 
since it's a bit easier. Run the code below to download and decompress 
the fastq files. They will be saved in a directory called 
``example_empirical_data/`` in your current directory. 
The total size is approximately 1.1GB.

.. code:: bash

    ## curl grabs the data from a public dropbox url
    ## the curl command uses an upper-case o argument, not a zero.
    curl -LskO https://dl.dropboxusercontent.com/u/2538935/example_empirical_rad.tar.gz
    
    ## the tar command decompresses the data directory
    tar -xvf example_empirical_rad.tar.gz


Starting an ipyrad analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Start by using the ``-n`` argument to ipyrad followed by a name
for your assembly. This creates a parameter input file (params-name.txt) 
which includes the Assembly name. I'll use the name ``base`` to start, 
to indicate that this is the base assembly from which we will later 
create new branches.

.. code:: bash

    ipyrad -n "base"


.. parsed-literal::

    New file `params-base.txt` created in /home/deren/Downloads



Edit the params file
~~~~~~~~~~~~~~~~~~~~

For this data set the data are already demultiplexed so we are going to
set the **sorted\_fastq\_path** to tell it the location of the fastQ
data files. I also change the **project\_dir** to "pedicularis". All
other parameters are left at their default values for now.

.. code:: python

    ## use your text editor to change the following params 
    ## use a wildcard character to select all 13 gzipped files.
    
    pedicularis                    ## [1] [project_dir] ...
    example_empirical_rad/*.gz     ## [4] [sorted_fastq_path] ...


Run step1 to load in the fastq data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ## Now run step 1 of the assembly 
    ## the -p flag tells ipyrad which assembly to use (params-base.txt)
    ## the -s flag tells ipyrad which step to run (1)
    
    ipyrad -p params-base.txt -s 1 


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      New Assembly: base
      ipyparallel setup: Local connection to 4 Engines
    
      Step1: Linking sorted fastq data to Samples
    
        Linking to demultiplexed fastq files in:
          /home/deren/Downloads/example_empirical_rad/*.gz
        13 new Samples created in `base`.
        13 fastq files linked to 13 new Samples.
        Saving Assembly.


We can use the -r flag to see the results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ipyrad -p params-base.txt -r



.. parsed-literal::

    
    
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
    
    
    Full stats files
    ------------------------------------------------
    step 1: None
    step 2: None
    step 3: None
    step 4: None
    step 5: None
    step 6: None
    step 7: None
    
    


Next we run step2 to filter the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assembling this complete data set takes several hours depending on how
many processors are available. Using four cores it can finish about 1.5
hours if we subsample the data set using the --preview method in ipyrad.
If run during step2 this function subsamples 100K reads from each
sample. As you can see above, this is only about 5-10% of the total
reads. This time I add the results call into the same cell so that it
runs step 2 and then prints the results.

.. code:: python

    %%bash
    ipyrad -p params-base.txt -s 2 --preview
    ipyrad -p params-base.txt -r


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step2: Filtering reads 
        Running preview mode: subselecting maximum of 100000 reads per sample    
        Saving Assembly.


.. code:: python

    %%bash
    ipyrad -p params-base.txt -r


.. parsed-literal::

    
    
    Summary stats of Assembly base
    ------------------------------------------------
                            state  reads_raw  reads_filtered
    29154_superba               2     696994           92448
    30556_thamno                2    1452316           93666
    30686_cyathophylla          2    1253109           89122
    32082_przewalskii           2     964244           92016
    33413_thamno                2     636625           89428
    33588_przewalskii           2    1002923           92418
    35236_rex                   2    1803858           92807
    35855_rex                   2    1409843           92883
    38362_rex                   2    1391175           93363
    39618_rex                   2     822263           92096
    40578_rex                   2    1707942           93386
    41478_cyathophylloides      2    2199740           93846
    41954_cyathophylloides      2    2199613           91756
    
    
    Full stats files
    ------------------------------------------------
    step 1: None
    step 2: ./pedicularis/base_edits/s2_rawedit_stats.txt
    step 3: None
    step 4: None
    step 5: None
    step 6: None
    step 7: None
    
    


Run step 3 (clustering and aligning)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is generally one of the longest running steps, depending on how
many unique clusters (loci) there are in each sample. Using more
processors will allow it run much faster. On my laptop with 4 cores this
step finishes in approximately 30 minutes. From the results you can see
that there are many clusters found in each sample (clusters\_total), but
many fewer (~2%) that were recovered at high depth (clusters\_hidepth).
The coverage would of course be much better if we did not subsample the
data set in step2 using --preview mode.

.. code:: python

    %%bash
    ipyrad -p params-base.txt -s 3
    ipyrad -p params-base.txt -r


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step3: Clustering/Mapping reads
        Saving Assembly.
    
    
    Summary stats of Assembly base
    ------------------------------------------------
                            state  reads_raw  reads_filtered  clusters_total  \
    29154_superba               3     696994           92448           45531   
    30556_thamno                3    1452316           93666           45745   
    30686_cyathophylla          3    1253109           89122           50306   
    32082_przewalskii           3     964244           92016           44242   
    33413_thamno                3     636625           89428           52053   
    33588_przewalskii           3    1002923           92418           46674   
    35236_rex                   3    1803858           92807           57801   
    35855_rex                   3    1409843           92883           45139   
    38362_rex                   3    1391175           93363           41580   
    39618_rex                   3     822263           92096           47295   
    40578_rex                   3    1707942           93386           45295   
    41478_cyathophylloides      3    2199740           93846           41965   
    41954_cyathophylloides      3    2199613           91756           47735   
    
                            clusters_hidepth  
    29154_superba                        978  
    30556_thamno                         987  
    30686_cyathophylla                   757  
    32082_przewalskii                    686  
    33413_thamno                         728  
    33588_przewalskii                    904  
    35236_rex                            767  
    35855_rex                           1106  
    38362_rex                           1140  
    39618_rex                           1258  
    40578_rex                            832  
    41478_cyathophylloides               992  
    41954_cyathophylloides              1307  
    
    
    Full stats files
    ------------------------------------------------
    step 1: None
    step 2: ./pedicularis/base_edits/s2_rawedit_stats.txt
    step 3: ./pedicularis/base_clust_0.85/s3_cluster_stats.txt
    step 4: None
    step 5: None
    step 6: None
    step 7: None
    
    


Run Step 4 (joint estimation of error rate & heterozygosity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step runs pretty fast. It should finish in about 10 minutes. As you
can see in the results the error rate is about 10X the heterozygosity
estimate. The latter does not vary significantly across samples. With
data of greater depth the estimates will be more accurate.

.. code:: python

    %%bash
    ipyrad -p params-base.txt -s 4 
    ipyrad -p params-base.txt -r


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step4: Joint estimation of error rate and heterozygosity
        Saving Assembly.
    
    
    Summary stats of Assembly base
    ------------------------------------------------
                            state  reads_raw  reads_filtered  clusters_total  \
    29154_superba               4     696994           92448           45531   
    30556_thamno                4    1452316           93666           45745   
    30686_cyathophylla          4    1253109           89122           50306   
    32082_przewalskii           4     964244           92016           44242   
    33413_thamno                4     636625           89428           52053   
    33588_przewalskii           4    1002923           92418           46674   
    35236_rex                   4    1803858           92807           57801   
    35855_rex                   4    1409843           92883           45139   
    38362_rex                   4    1391175           93363           41580   
    39618_rex                   4     822263           92096           47295   
    40578_rex                   4    1707942           93386           45295   
    41478_cyathophylloides      4    2199740           93846           41965   
    41954_cyathophylloides      4    2199613           91756           47735   
    
                            clusters_hidepth  hetero_est  error_est  
    29154_superba                        978    0.038530   0.006630  
    30556_thamno                         987    0.038266   0.006009  
    30686_cyathophylla                   757    0.044680   0.004627  
    32082_przewalskii                    686    0.046796   0.007077  
    33413_thamno                         728    0.041466   0.004528  
    33588_przewalskii                    904    0.041445   0.011253  
    35236_rex                            767    0.042423   0.005119  
    35855_rex                           1106    0.035123   0.012086  
    38362_rex                           1140    0.041206   0.004702  
    39618_rex                           1258    0.040696   0.009077  
    40578_rex                            832    0.045177   0.002789  
    41478_cyathophylloides               992    0.041085   0.004468  
    41954_cyathophylloides              1307    0.032387   0.013090  
    
    
    Full stats files
    ------------------------------------------------
    step 1: None
    step 2: ./pedicularis/base_edits/s2_rawedit_stats.txt
    step 3: ./pedicularis/base_clust_0.85/s3_cluster_stats.txt
    step 4: ./pedicularis/base_clust_0.85/s4_joint_estimate.txt
    step 5: None
    step 6: None
    step 7: None
    
    


Run step 5 (consensus base calls)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is another step that can be computationally intensive. Here it
takes about 15 minutes on 4 cores. Although many clusters are filtered
out at this step (especially due to low depth) their information is
retained for the VCF output later so that the coverage/depth of excluded
reads can be examined.

.. code:: python

    %%bash
    ipyrad -p params-base.txt -s 5
    ipyrad -p params-base.txt -r



.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step5: Consensus base calling 
        Diploid base calls and paralog filter (max haplos = 2)
        error rate (mean, std):  0.00703, 0.00331
        heterozyg. (mean, std):  0.04071, 0.00396
        Saving Assembly.
    
    
    Summary stats of Assembly base
    ------------------------------------------------
                            state  reads_raw  reads_filtered  clusters_total  \
    29154_superba               5     696994           92448           45531   
    30556_thamno                5    1452316           93666           45745   
    30686_cyathophylla          5    1253109           89122           50306   
    32082_przewalskii           5     964244           92016           44242   
    33413_thamno                5     636625           89428           52053   
    33588_przewalskii           5    1002923           92418           46674   
    35236_rex                   5    1803858           92807           57801   
    35855_rex                   5    1409843           92883           45139   
    38362_rex                   5    1391175           93363           41580   
    39618_rex                   5     822263           92096           47295   
    40578_rex                   5    1707942           93386           45295   
    41478_cyathophylloides      5    2199740           93846           41965   
    41954_cyathophylloides      5    2199613           91756           47735   
    
                            clusters_hidepth  hetero_est  error_est  reads_consens  
    29154_superba                        978    0.038530   0.006630            821  
    30556_thamno                         987    0.038266   0.006009            810  
    30686_cyathophylla                   757    0.044680   0.004627            606  
    32082_przewalskii                    686    0.046796   0.007077            523  
    33413_thamno                         728    0.041466   0.004528            597  
    33588_przewalskii                    904    0.041445   0.011253            709  
    35236_rex                            767    0.042423   0.005119            629  
    35855_rex                           1106    0.035123   0.012086            844  
    38362_rex                           1140    0.041206   0.004702            943  
    39618_rex                           1258    0.040696   0.009077           1011  
    40578_rex                            832    0.045177   0.002789            689  
    41478_cyathophylloides               992    0.041085   0.004468            872  
    41954_cyathophylloides              1307    0.032387   0.013090            983  
    
    
    Full stats files
    ------------------------------------------------
    step 1: None
    step 2: ./pedicularis/base_edits/s2_rawedit_stats.txt
    step 3: ./pedicularis/base_clust_0.85/s3_cluster_stats.txt
    step 4: ./pedicularis/base_clust_0.85/s4_joint_estimate.txt
    step 5: ./pedicularis/base_consens/s5_consens_stats.txt
    step 6: None
    step 7: None
    
    


Step 6 (clustering and aligning across samples)
-----------------------------------------------

This step clusters consensus loci across Samples using the same
threshold for sequence similarity as used in step3.

.. code:: python

    %%bash
    ipyrad -p params-base.txt -s 6


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step6: Clustering across 13 samples at 0.85 similarity
        Saving Assembly.


Branch the assembly
-------------------

Here we will branch the assembly to create different assemblies that we
will use as our final outputs. The main parameter we will focus on is
the ``min_samples_locus``, which is the minimum number of samples that
must have data at a locus for the locus to be retained in the data set.
We create a ``min4``, ``min8``, and ``min12`` data sets.

.. code:: python

    %%bash
    
    ipyrad -p params-base.txt -b min4
    ipyrad -p params-base.txt -b min8
    ipyrad -p params-base.txt -b min12



.. parsed-literal::

    
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      Creating a branch of assembly base called min4
      Writing new params file to params-min4.txt
    
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      Creating a branch of assembly base called min8
      Writing new params file to params-min8.txt
    
      loading Assembly: base [~/Downloads/pedicularis/base.json]
      Creating a branch of assembly base called min12
      Writing new params file to params-min12.txt


Change the parameter settings in params.txt for each assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ## Enter the changes below into the params files using a text editor
    
    ## params-min4.txt
    4       ## [21] [min_samples_locus] ...
    
    ## params-min8.txt
    8       ## [21] [min_samples_locus] ...
    
    ## params-min12.txt
    12      ## [21] [min_samples_locus] ...


Step 7 (final filtering and create output files)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filter and create output files for the three assemblies with different
values for the parameter ``min_samples_locus``.

.. code:: python

    %%bash
    ## assemble the three data sets
    ## take note: The force argument is required
    ## if the base assembly already created output files 
    ## for step 7. 
    ipyrad -p params-min4.txt -s 7 -f
    ipyrad -p params-min8.txt -s 7 -f
    ipyrad -p params-min12.txt -s 7 -f


.. parsed-literal::

    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: min4 [~/Downloads/pedicularis/min4.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step7: Filter and write output files for 13 Samples.
        Outfiles written to: ~/Downloads/pedicularis/min4_outfiles
        Saving Assembly.
    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: min8 [~/Downloads/pedicularis/min8.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step7: Filter and write output files for 13 Samples.
        Outfiles written to: ~/Downloads/pedicularis/min8_outfiles
        Saving Assembly.
    
     --------------------------------------------------
      ipyrad [v.0.1.70]
      Interactive assembly and analysis of RADseq data
     --------------------------------------------------
      loading Assembly: min12 [~/Downloads/pedicularis/min12.json]
      ipyparallel setup: Local connection to 4 Engines
    
      Step7: Filter and write output files for 13 Samples.
        Outfiles written to: ~/Downloads/pedicularis/min12_outfiles
        Saving Assembly.


.. code:: python

    cat ./pedicularis/min4_outfiles/min4_stats.txt


.. parsed-literal::

    
    
    ## The number of loci caught by each filter.
    ## ipyrad API location: [assembly].statsfiles.s7_filters
    
                               locus_filtering
    total_prefiltered_loci                1206
    filtered_by_rm_duplicates                2
    filtered_by_max_indels                 159
    filtered_by_max_snps                     0
    filtered_by_max_hetero                  15
    filtered_by_min_sample                 921
    filtered_by_edge_trim                    0
    total_filtered_loci                    221
    
    
    ## The number of loci recovered for each Sample.
    ## ipyrad API location: [assembly].stats_dfs.s7_samples
    
                            sample_coverage
    29154_superba                       151
    30556_thamno                        120
    30686_cyathophylla                  118
    32082_przewalskii                   132
    33413_thamno                        155
    33588_przewalskii                    89
    35236_rex                           154
    35855_rex                           145
    38362_rex                           154
    39618_rex                           156
    40578_rex                            90
    41478_cyathophylloides              156
    41954_cyathophylloides               97
    
    
    ## The number of loci for which N taxa have data.
    ## ipyrad API location: [assembly].stats_dfs.s7_loci
    
        locus_coverage  sum_coverage
    1              NaN             0
    2              NaN             0
    3              NaN             0
    4               65            65
    5               36           101
    6               16           117
    7               10           127
    8               10           137
    9                6           143
    10               2           145
    11              12           157
    12               7           164
    13              57           221
    
    
    ## The distribution of SNPs (var and pis) across loci.
    ## pis = parsimony informative site (minor allele in >1 sample)
    ## var = all variable sites (pis + autapomorphies)
    ## ipyrad API location: [assembly].stats_dfs.s7_snps
    
        var  sum_var   pis  sum_pis
    0   260      260  1140     1140
    1   130      390    45     1185
    2   130      520    12     1197
    3   133      653     3     1200
    4    90      743     4     1204
    5   110      853     0     1204
    6    77      930     0     1204
    7    70     1000     2     1206
    8    64     1064     0     1206
    9    32     1096     0     1206
    10   31     1127     0     1206
    11   30     1157     0     1206
    12   16     1173     0     1206
    13   14     1187     0     1206
    14    5     1192     0     1206
    15    6     1198     0     1206
    16    2     1200     0     1206
    17    4     1204     0     1206
    18    2     1206     0     1206

Take a peek at the .loci (easily human-readable) output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ## head -n 50 prints just the first 50 lines of the file to stdout
    head -n 50 pedicularis/min4_outfiles/min4.loci


.. parsed-literal::

    29154_superba              CCTTGGTSACCTTMGCWCCWGAYGGRTCCTTCTTCTCCACACTCTTKATRACACCAACAGCAACAGTC
    32082_przewalskii          CCTTGGTSACCTTRGCWCCWGAYGG-TCCTTCTTCTCCACACTCTTGATRACACCAACAGCAACAGTC
    35236_rex                  CCTTGGTCACCTTAGCACCTGATGG-TCCTTCTTCTCCACACTCTTGATGACACCAACAGCAACAGTC
    38362_rex                  CCTTGGTCACCTTAGCACCTGATGGRTCCTTCTTCTCCACACTCTTGATGACACCAACAGCAAC-GTC
    //                                -     -  -  -  -  -                    -  -                  |
    33413_thamno               TAGACAACCAGTGCCTTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    33588_przewalskii          GAGACAACCAGTGCCTTCTTGTCTATCAGCCTCACACCTGTCTTCGGTACTTTCGGTACTTAGAAGCA
    35855_rex                  TAGACAACCAGTGCCGTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    38362_rex                  TAGACAACCAGTGCCTTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    39618_rex                  TAGACAACCAGTGCCTTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    40578_rex                  TAGACAACCAGTGCCTTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    41478_cyathophylloides     TAGACAACCAGTGCCGTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    41954_cyathophylloides     TAGACAACCAGTGCCGTCTTGTCTATCAGTCTCACACCTGTCTTCGGTACTTGCGGTACTTAGAAGCA
    //                         -              *             -                      -               |
    29154_superba              AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    30556_thamno               AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    30686_cyathophylla         AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    32082_przewalskii          AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    33413_thamno               AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    33588_przewalskii          AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    35236_rex                  AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    35855_rex                  AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    38362_rex                  AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    39618_rex                  AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    40578_rex                  AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    41478_cyathophylloides     AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGC-CTAGCTAACGCGCCC
    41954_cyathophylloides     AGCAAGCGAAGAAAACGTAAGGGCGCGCGTTAGCACTCCTGCAAGAAAACGGCNCTAGCTAACGCGCCC
    //                                                                                              |
    30686_cyathophylla         TAGCAATAAATGCAAGAATATTTACTTCCATAATTTCGTCGGTTTTTTAATTCGCAATAACTCGGGAT
    32082_przewalskii          TAGCAATAAATGCAAGAATATTGACTTCCATAATTTCGTCGGTTTTTTAATTCGCAATAACTCGGGAT
    33588_przewalskii          TAGCAATAAATGCAAGAATATTGACTTCCATAATTTCGTCGGTTTTTTAATTCGCAATAACTCGGGAT
    35236_rex                  TAGCAATAAATGCAAGAATATTKACTTCCATAATTTCGTCKGTTTTTTAATTCGCAATAACTCGGGAT
    38362_rex                  TAGCAATAAATGCAAGAATATTTACTTCCATAATTTCGTCTGTTTTTTAATTCGCAATAACTCGGGAT
    39618_rex                  TAGCAATAAATGCAAGAATATTTACTTCCATAATTTCGTCTGTTTTTTAATTCGCAATAACTCGGGAT
    40578_rex                  TAGCAATAAATGCAAGAATATTTACTTCCATAATTTCGTCTGTTTTTTAATTCGCAATAACTCGGGAT
    41478_cyathophylloides     TAGCAATAAATGCAAGAATATTT-CTTCCATAATTTCGTCGGTTTTTTAATTCGCAATAACTCGGGAT
    //                                               *                 *                           |
    35236_rex                  CTCTAGGTGGAGCTCCAGCTGGGTCTGAACCAGATCCTCCGTAAKCGGATCATCATGTGCGAGTTGAC
    35855_rex                  CTCTAGGTGGAGCTCCAGCTGGGTCTGAACCAGATCCTCCGTAAGCGGATCATCATGTGCGAGTGGAC
    38362_rex                  CTCTAGGTGGAGCTCCAGCTGGGTCTGAACCAGATCCTCCGTAAGCGGATCATCATGTGCGAGTTGAC
    39618_rex                  CTCTAGGTGGAGCTCCAGCTGGGTCTGAACCAGATCCTCCGTAAGCGGATCATCATGTGCGAGTTGAC
    //                                                                     -                   -   |
    30556_thamno               C-TTCTGATTAATCTG-AAATTGTAATCAAATGAAATYAAACAGCAAAAACAATGACTSGATAAACTA
    33413_thamno               CTTTCTGWTTAATCTGMAAATTGTAATCAAATGAAATCAAACARCAAAAACAATGACTYGAYAAWCYR
    35236_rex                  C-TTCTGATTAATCTG-AAATTGTAATCAAATGAAATCAAACA-CAAAAACAATGACT-GATAAACTA
    41478_cyathophylloides     CTTTCTGATTAATCTGCAAATTGTAATCAAATGAAATCAAACAGCAAAAACAATAACTTGATAAAATA
    //                                -        -                    -     -          -   -  -  ----|
    30556_thamno               GAAAGATWT-AYTGTAGACGTAWTKGATCRSAGGWKGAGGTGATGWATCATAWTCAT-ATCAGAGGAG
    38362_rex                  GAAAGATTTCACTGTAGACGTAATGGATCAGAGGTTGAGGTGATGRATCATAATCATGATCAGAGGWG
    39618_rex                  GAAAGATTTCACTGTAGACGTAWWGGATCMSAGGWKGAGGTGATGRATCATAATCATKATCAGAGGAG


peek at the .phy files
~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/min4_outfiles/min4.phy


.. parsed-literal::

    13 15034
    29154_superba              CCTTGGTSACCTTMGCWCCWGAYGGRTCCTTCTTCTCCACACTCTTKATRACA
    30556_thamno               NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    30686_cyathophylla         NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    32082_przewalskii          CCTTGGTSACCTTRGCWCCWGAYGGNTCCTTCTTCTCCACACTCTTGATRACA
    33413_thamno               NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    33588_przewalskii          NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    35236_rex                  CCTTGGTCACCTTAGCACCTGATGGNTCCTTCTTCTCCACACTCTTGATGACA
    35855_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    38362_rex                  CCTTGGTCACCTTAGCACCTGATGGRTCCTTCTTCTCCACACTCTTGATGACA
    39618_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    40578_rex                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    41478_cyathophylloides     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    41954_cyathophylloides     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN


peek at the .snp file
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/min4_outfiles/min4.snp


.. parsed-literal::

    13 711
    29154_superba              SMWWYRKRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGYATNYGAMGT
    30556_thamno               NNNNNNNNNNNNNNNNA-YGGSTACTACWYWTKRSWKWW-TAGTAT-NNNNGT
    30686_cyathophylla         NNNNNNNNNNNNTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCRACTT
    32082_przewalskii          SRWWY-GRNNNNGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATNNNNNGG
    33413_thamno               NNNNNNNNTTTGNNNNWMCRGYYWCYRYNNNNNNNNNNNNNNNNNNNNNNNGT
    33588_przewalskii          NNNNNNNNGTCTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNRTRKNNNNNGG
    35236_rex                  CAATT-GGNNNNKKKTA-C-G-TACTACNNNNNNNNNNNNNNG-ATMNNNNGT
    35855_rex                  NNNNNNNNTGTGNNGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGRCGT
    38362_rex                  CAATTRGGTTTGTTGTNNNNNNNNNNNNTCATGAGTTRAGTWNNNNMCGACGT
    39618_rex                  NNNNNNNNTTTGTTGTNNNNNNNNNNNNTCWWGMSWKRAKTAGTATMCGRCGT
    40578_rex                  NNNNNNNNTTTGTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATNNNNNGT
    41478_cyathophylloides     NNNNNNNNTGTGTGNNACCGATTAATACWKATKRSWKRAGYAGTATNCRACGT
    41954_cyathophylloides     NNNNNNNNTGTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGT


peek at the .snp file for the min12 assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    %%bash
    ## cut -c 1-80 prints only the first 80 characters of the file
    cut -c 1-80 pedicularis/min12_outfiles/min12.snp


.. parsed-literal::

    13 44
    29154_superba              GTTGACGCGGTTCTCCCCCGCGCAGTGCAGTCCGCTAANNAGTT
    30556_thamno               GTTGACGCGGTTCTCACCCGCGCAGCGCGGTACACTAAGAAATT
    30686_cyathophylla         TTTGACGCGGTTCTCACCCGCGCAGCGCGGTCCACTAAGAAATT
    32082_przewalskii          GGTGATGCAACCCCTACCCGCGTGGCGTARYCAGTTARGAGACC
    33413_thamno               GTTGACGCGGTTTTCACCCGCGCAGCGCGGTACACTAAGAAATT
    33588_przewalskii          GGTGATGCAACCCCTACCCGCGTG-TC-ARYCAGTTARGAGACC
    35236_rex                  GTTGACGCGGTTCTCACCCGCGCAGCGCGRYACACTWAKMAATT
    35855_rex                  GTKRMYKMGGTTCTCACCCGCGCAGCGCGGTACACTAAGANNNT
    38362_rex                  GTTGACGCGGTTCTCACCCGCGCAGCGCGGTCCACTAAGAAATT
    39618_rex                  GTTGACGCGGTTCTCAYYYRYRCAGCGCGGTCCACTAAGAAATT
    40578_rex                  GTTGACGCGGTTCTCACCCGCGCAGCGCGGTACACTAAGAAATT
    41478_cyathophylloides     GTTGACGCGGTTCTCACCCGCGCAATGCAGTCCGCCAAGAAATT
    41954_cyathophylloides     GTTGACGCGGTTCTCACCCGCGCAATGCAGTCCGCCAAGAAATT



downstream...
~~~~~~~~~~~~~

