
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
    13 new Samples created in `base`.
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
To fully assemble this data set on a 4-core laptop takes about 7 hours. With 
access to 24 cores it would take closer to 1 hour. 


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




Take a look at the stats summary 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each assembly that finishes step 7 will create a stats.txt output summary
in the 'assembly_name'_outfiles/ directory. This includes information about 
which filters removed data from the assembly, how many loci were recovered
per sample, how many samples had data for each locus, and how many variable
sites are in the assembled data. 


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


Take a peek at the .loci output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the first places I look when an assembly finishes. It provides a 
clean view of the data with variable sites (-) and parsimony informative
SNPs (*) highlighted. Use the unix commands **less** or **head** to look at this
file briefly.

.. code:: bash

    ## head -n 50 prints just the first 50 lines of the file to stdout
    head -n 50 pedicularis/min4_outfiles/min4.loci



peek at the .phy files
~~~~~~~~~~~~~~~~~~~~~~
This is the concatenated sequence file of all loci in the data set. It is typically
used in phylogenetic analyses, like in the program *raxml*. 


.. code:: bash

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
This is similar to the phylip file format, but only variable site columns are 
included. All SNPs are the file, in contrast to the .usnps file, which selects
only a single SNP per locus. 


.. code:: bash

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
Similar to above but you can see that there is much less missing data (Ns). 

.. code:: bash

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


downstream analyses
~~~~~~~~~~~~~~~~~~~
We are in the process of developing many downstream tools for analyzing 
RAD-seq data. The first of which is the program svd4tet, which is installed 
alongside ipyrad during the conda installation. 
This program implements the svdquartets algorithm of Chifmann & Kubatko (2014). 
It includes alternative methods for heuristically inferring quartets over
very large trees, but for this demonstration we'll simply run the default 
option to infer all quartets. As you can see, the summary tree that it spits
out is the same tree inferred in Eaton & Ree (2013), and if we plot the 
tree with support values you will see that we find lower support across the 
edges of the tree that are known to be involved in introgression. This 
part of the tutorial is under construction -- more to come. 

.. code:: bash

    ## run svd4tet on the unlinked snps file (-s) (.u.snps.phy suffix)
    ## and give it the output prefix (-o) 'pedictree'
    svd4tet -s pedicularis/base_outfiles/base.u.snps.phy -o pedictree

.. parsed-literal::
    
