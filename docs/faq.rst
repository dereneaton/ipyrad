
.. _faq:  
Frequenty Asked Questions
=========================

Tuning Assembly Parameters
==========================
This is a collection of papers related to tuning different assembly
parameters and how these may or may not impact downstream analysis. 
Just because we include a paper here doesn't mean we endorse the results, 
only that we think it's worth considering. Where necessary, we 
highlight specific papers and interesting/important results below.

* Huang & Knowles 2016 - `Unforeseen Consequences of Excluding Missing Data from Next-Generation Sequences: Simulation Study of RAD Sequences <https://academic.oup.com/sysbio/article/65/3/357/2468879>`__
    Especially useful. TL;DR Don't over-filter missing data. In practice this means allowing for more permissive `min_samples_locus` parameter settings. Sometimes people have a tendency to treat RAD data as if it were a giant multi-locus Sanger dataset. They want to see a "complete data matrix", so use very high `min_samples_locus` settings in order to reduce the amount of missing data. This is bad, and this paper shows us why.
* Shafer et al. 2016 - `Bioinformatic processing of RAD-seq data dramatically impacts downstream population genetic inference <https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12700>`__
    
* Linck & Battey 2018 - `Minor allele frequency thresholds dramatically affect population structure inference with genomic datasets <https://www.biorxiv.org/content/biorxiv/early/2018/10/21/188623.full.pdf>`__
    Don't over-filter your data. Filtering too strictly on MAF reduces power to detect population structure. Given that low-frequency variants retain much of the signal of recent history this result seems sensible.

Troubleshooting Procedures
==========================

Troubleshooting ipyparallel issues
----------------------------------
Sometimes ipyrad can have trouble talking to the ipyparallel
cluster on HPC systems. First we'll get an interactive shell
on an HPC compute node (YMMV with the `qsub -I` here, you might
need to specify the queue and allocate specific resource).

.. code-block:: bash

    qsub -I

.. code-block:: bash

    ipcluster start --n 4 --daemonize

Then type `ipython` to open an ipython session.

.. code-block:: python

    import ipyparallel as ipp

    rc = ipp.Client()
    rc[:]

The result should look something like this:
.. parsed-literal::

    Out[1]: <DirectView [0, 1, 2, 3]>

.. code-block:: python

    import ipyparallel as ipp

    rc = ipp.Client(profile="default")
    rc[:]

.. code-block:: python

    import ipyrad as ip

    ## location of your json file here
    data = ip.load_json("dir/path.json")

    print data._ipcluster

.. code-block:: python

    data = ip.Assembly('test')

    data.set_params("raw_fastq_path", "path_to_data/\*.gz")
    data.set_params("barcodes_path", "path_to_barcode.txt")

    data.run('1')

    print data.stats
    print data._ipcluster

.. parsed-literal::

    {'profile': 'default', 'engines': 'Local', 'quiet': 0, 'cluster_id': '', 'timeout': 120, 'cores': 48}

.. code-block:: python

    data.write_params('params-test.txt')

Don't forget to stop the ipcluster when you are done.

.. code-block:: bash

    ipcluster stop

Running ipyrad on HPC that restricts write-access to /home on compute nodes
---------------------------------------------------------------------------

Some clusters forbid writing to `/home` on the compute nodes. It guarantees that users 
only write to scratch drives or high performance high volume disk, and not the user 
home directory (which is probably high latency/low volume). They have write access on 
login, just not inside batch jobs. This manifests in weird ways, it's hard to debug,
but you can fix it by adding an `export` inside your batch script.

.. code-block:: bash

    export HOME=/<path>/<to>/<some>/<writable>/<dir>

In this way, `ipcluster` and `ipyrad` will both look in `$HOME` for the `.ipython` directory.

ipyrad crashes during dereplication in step 3
---------------------------------------------

.. parsed-literal::

    ERROR sample [XYZ] failed in step [derep_concat_split]; error: EngineError(Engine '68e79bbc-0aae-4c91-83ec-97530e257387' died while running task u'fdef6e55-dcb9-47cb-b4e6-f0d2b591b4af')

If step 3 crashes during dereplication you may see an error like above. Step 3
can take quite a lot of memory if your data do not de-replicate very efficiently.
Meaning that the sample which failed may contain a lot of singleton reads. 

You can take advantage of the following steps during step 2 to better filter your 
data so that it will be cleaner, and thus dereplicate more efficiently. This will
in turn greatly speed up the step3 clustering and aligning steps. 

* Use the "filter_adapters" = 2 argument in ipyrad which will search for and remove Illumina adapters. 
* Consider trimming edges of the reads with the "trim_reads" option. An argument like (5, 75, 5, 75) would trim the first five bases of R1 and R2 reads, and trim all reads to a max length of 75bp. Trimming to a fixed length helps if your read qualities are variable, because the reads may be trimmed to variable lengths. 
* Try running on a computer with more memory, or requesting more memory if on a cluster.

Collisions with other local python/conda installs
-------------------------------------------------

.. parsed-literal::

    Failed at nopython (nopython frontend)
    UntypedAttributeError: Unknown attribute "any" of type Module(<module 'numpy' from...

In some instances if you already have conda/python installed the local environment
variable PYTHONPATH will be set, causing python to use versions of modules 
outside the miniconda path set during ipyrad installation. This error can be fixed by 
blanking the PYTHONPATH variable during execution (as below), or by adding the export
to your ~/.bashrc file.

.. code-block:: bash

    export PYTHONPATH=""; ipyrad -p params.txt -s 1

Why doesn't ipyrad handle PE original RAD?
------------------------------------------
Paired-End RAD protocol is tricky to denovo assemble. Because of the sonication step R2 
doesn't line up nicely. ipyrad makes strong assumptions about how r1 and r2 align, 
assumptions which are met by PE gbs and ddrad, but which are not met by original RAD. 
This doesn't matter (as much) if you have a reference genome, but if you don't have a 
reference it's a nightmare... dDocent has a PE-RAD mode, but I haven't evaluated it. 
I know that people have also used stacks (because stacks treats r1 andr2 as independent 
loci). If people ask me how to denovo assemble with PE-RAD in ipyrad I tell them to 
just assemble it as SE and ignore R2.

Why doesn't ipyrad write out the .alleles format with phased alleles like pyrad used to?
----------------------------------------------------------------------------------------
We're hoping to provide something similar eventually, the problem with the pyrad alleles 
file is that the alleles are only phased correctly when we enforce that reads must align 
almost completely, i.e., they are not staggered in their overlap. So the alleles are 
correct for RAD data, because the reads match up perfectly on their left side, however, 
staggered overlaps are common in other data sets that use very common cutters, like 
ezRAD and some GBS, and especially so when R1 and R2 reads merge. So we needed to change 
to an alternative way of coding the alleles so that we can store both phased and unphased 
alleles, and its just taking a while to do. So for now we are only providing unphased 
alleles, although we do save the estimated number of alleles for each locus. This 
information is kind of hidden under the hood at the moment though.

Why is my assembly taking FOREVER to run?
-----------------------------------------
There have been a few questions recently about long running jobs (e.g., >150 hours), which 
in my experience should be quite rare when many processors are being used. In general, 
I would guess that libraries which take this long to run are probably overloaded with 
singleton reads, meaning reads are not clustering well within or across samples. This 
can happen for two main reasons: (1) Your data set actually consists of a ton of 
singleton reads, which is often the case in libraries that use very common cutters like 
ezRAD; or (2) Your data needs to be filtered better, because low quality ends and 
adapter contamination are causing the reads to not cluster.

If you have a lot of quality issues or if your assemby is taking a long time to cluster 
here are some ways to filter more aggressively, which should improve runtime and the
quality of the assembly:

* Set filter_adapters to 2 (stringent=trims Illumina adapters)
* Set phred_Qscore_offset to 43 (more aggressive trimming of low quality bases from 3' end of reads
* Hard trim the first or last N bases from raw reads by setting e.g., trim_reads to (5, 5, 0, 0)
* Add additional 'adapter sequences' to be filtered (any contaminant can be searched for, I have added long A-repeats in one library where this appeared common). This can be done easily in the API, but requires editing the JSON file for the CLI.

I still don't understand the `max_alleles_consens` parameter
------------------------------------------------------------
In step 5 base calls are made with a diploid model using the parameters estimated in
step 4. The only special case in when `max_alleles_consens` = 1, in which case the step 4
heterozygosity estimate will be fixed to zero and the error rate will suck up all of the 
variation within sites, and then the step 5 base calls will be haploid calls. For all 
other values of `max_alleles_consens`, base calls are made using the diploid model using 
the H and E values estimated in step 4. **After site base calls are made** ipyrad then counts 
the number of alleles in each cluster. This value is now simply stored in step 5 for use 
later in step 7 to filter loci, under the assumption that if a locus has paralogs in one 
sample then it probably has them in other samples but there just wasn't enough variation to 
detect them.

Why does it look like ipyrad is only using 1/2 the cores I assign, and what does the `-t` flag do?
--------------------------------------------------------------------------------------------------
Most steps of ipyrad perform parallelization by multiprocessing, meaning that jobs are 
split into smaller bits and distributed among all of the available cores. However, some 
parts of the analysis also use multithreading, where a single function is performed over 
multiple cores. More complicated, parts like step3 perform several multithreaded jobs in 
parallel using multiprocessing... you still with me? The -c argument is the total number 
of cores that are available, while the -t argument allows more fine-tuned control of how 
the multithreaded functions will be distributed among those cores. For example, the 
default with 40 cores and -t=2 would be to start 20 2-threaded vsearch jobs. There are 
some parts of the code that cannot proceed until other parts finish, so at some points 
the code may run while using fewer than the total number of cores available, which is 
likely what you are seeing in step 3. Basically, it will not start the aligning step 
until all of the samples have finished clustering. It's all fairly complicated, but we 
generally try to keep everything working as efficiently as possible. If you have just 
one or two samples that are much bigger (have more data) than the rest, and they are 
taking much longer to cluster, then you may see a speed improvement by increasing the 
threading argument (e.g., -t 4).

How to fix the GLIBC error
--------------------------
If you ever see something that looks like this `/lib64/libc.so.6: version `GLIBC_2.14' not found`
it's probably because you are on a cluster and it's using an old version of GLIBC. To
fix this you need to recompile whatever binary isn't working on your crappy old machine.
Easiest way to do this is a conda local build and install. Using `bpp` as the example:

.. parsed-literal:: 

    git clone https://github.com/dereneaton/ipyrad.git
    conda build ipyrad/conda.recipe/bpp/
    conda install --use-local bpp

How do I interpret the `distribution of SNPs (var and pis) per locus` in the *_stats.txt output file
----------------------------------------------------------------------------------------------------
Here is an example of the first few lines of this block in the stats file:

.. parsed-literal:: 

    bash    var  sum_var    pis  sum_pis
    0    661        0  10090        0
    1   1660     1660   5070     5070
    2   2493     6646   1732     8534
    3   2801    15049    483     9983
    4   2683    25781    147    10571
    5   2347    37516     59    10866
    6   1740    47956     17    10968
    7   1245    56671      7    11017

**pis** is exactly what you think, it's the count of loci with *n* parsimony informative sites. So row 0 is loci with no pis, row 1 is loci with 1 pis, and so on.

**sum_pis** keeps a running total of the counts for all pis across all loci up to that point, which is why the sum looks weird, but i assure you its fine. For the row that records 3 pis per site, you see the # pis = 483 and 483 * 3 + 8534 = 9983.

**var** is a little trickier and here's where the docs are a little goofy. This keeps track of the number of loci with n variable sites including autapomorphies and pis within each locus. So row 0 is all totally monomorphic loci. row 1 is all loci with *either* one pis or one autapomorphy. Row 2 is all loci with *either* two pis, or two autapomorphies, *OR* one of each, and so on.

**sum_var** is calculated identical to **sum_pis**, so it does look weird but it's right.

The reason the counts in, for example, row 1 do not appear to agree for var and pis is because the value of row 1 for pis *includes all* loci with only one pis irrespective of the number of autapomorphies, whereas the value for var records all loci with *only one* of either of these. 

How to fix the `IOError(Unable to create file IOError(Unable to create file...` error
-------------------------------------------------------------------------------------
The HDF5_USE_FILE_LOCKING error is caused by the fact that your cluster filesystem is NFS (or some other network based filesystem). You can disable hdf5 file locking by setting an environment variable `export  HDF5_USE_FILE_LOCKING=FALSE`. See here for more info:

http://hdf-forum.184993.n3.nabble.com/HDF5-files-on-NFS-td4029577.html

Why am I getting the 'empty varcounts' error during step 7?
-----------------------------------------------------------
Occasionally during step 7 you will see this error:

.. code-block::
    Exception: empty varcounts array. This could be because no samples                                                                                                    
    passed filtering, or it could be because you have overzealous filtering.                                                                                              
    Check the values for `trim_loci` and make sure you are not trimming the                                                                                               
    edge too far.

This can actually be caused by a couple of different problems that all result in the same behavior, namely that you are filtering out *all* loci.

**trim_loci** It's true that if you set this parameter too aggressively all loci will be trimmed completely and thus there will be no data to output.

**min_samples_locs** Another way you can eliminate all data is by setting this parameter too high. Try dropping it way down, to like 3, then rerunning to get a better idea of what an appropriate value would be based on sample depths.

**pop_assign_file** A third way you can get this error is related to the previous one. The last line of the pop_assign_file is used for specifying min_sample per population for writing a locus. If you mis-specify the values for the pops in this line then it's possible to filter out all your data and thus obtain the above error.

How do I fix this error: "OSError: /lib64/libpthread.so.0: version `GLIBC_2.12' not found"?
-------------------------------------------------------------------------------------------
This error crops up if you are running ipyrad on a cluster that has an older version of GLIBC. The way to work around this is to install specific versions of some of the requirements that are compiled for the older version. Thanks to Edgardo M. Ortiz for this solution.


First clean up your current environment:

.. code-block:: bash

    module unload python2
    rm -rf miniconda2 .conda

    bash Miniconda2-latest-Linux-x86_64.sh
    source ~/.bashrc

then install the old version of llvmlite (and optionally the old versions of pyzmq and ipyparallel if necessary):

.. code-block:: bash

    conda install llvmlite=0.22

    conda install pyzmq=16
    conda install ipyparallel=5.2

and finally reinstall ipyrad:

.. code-block:: bash
    conda install -c ipyrad ipyrad
    conda install toytree -c eaton-lab

optional:

.. code-block:: bash
    conda clean --all

Why am i getting the 'ERROR   R1 and R2 files are not the same length.' during step 1?
--------------------------------------------------------------------------------------
This is almost certainly a disk space issue. Please be sure you have _plenty_ of disk space on whatever drive you're doing your assembly on. Running out of disk can cause weird problems that seem to defy logic, and that are a headache to debug (like this one). Check your disk space: `df -h`

Why does the number of pis recovered in the output stats change when I change the value of `max_snp_locus`?
-----------------------------------------------------------------------------------------------------------
While it does seem that the # of pis shouldn't change under varying `max_snp_locus` thresholds, it is in fact not true. This is because the setting is for max __SNP__ per locus, not max __PIS__. So for example if you have `max_snp_locus` set to 5, and you have a locus with 5 singleton snps and one doubleton snp (which is parsimony informative), then this locus would be filtered out. However if you set `max_snp_locus` to 10, then this locus would be included and the 'pis' counter would be incremented by 1. In this way you can see that the number of PIS recovered will change because of variation in this parameter setting.

Can ipyrad assemble MIG-seq data?
---------------------------------
MIG-seq (multiplexed ISSR genotyping by sequencing) is a method proposed by Suyama and Matsuki (2015), which involves targeting variable regions between simple sequence repeats (SSR). The method produces data that is somewhat analogous to ddRAD, in that you have the variable region which is flanked on either side by sequences that are known to be repeated randomly and at some appreciable frequency throughout the genome. Check out the `figure from the manuscript <https://www.nature.com/articles/srep16963/figures/1>`__. Anyway.... yes, ipyrad can assemble this kind of data, though there are some tricks. Primarily we recommend higher values of `filter_min_trim_len` and `clust_threshold`. If sequenced on a desktop NGS platform (Ion Torrent PGM, MiSeq) it also helps to reduce both `mindepth` params to recover more clusters.

Why are my ipcluster engines dying silently on cluster compute notes?
---------------------------------------------------------------------
This is a nasty bug that's bitten me more than once. If you are having trouble with cluster engines running jobs and then dying silently it may be because the cluster is headless and the engines are trying to interact with a GUI backend. This causes nasty things to happen. Here are a couple links that provide workable solutions:

https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/o0pnE9PEqA0

https://github.com/ipython/ipyparallel/issues/213

Why are my ipyrad.analysis.structure runs taking so long/not doing anything?
----------------------------------------------------------------------------
See the previous FAQ answer. It's typical for HPC cluster systems to be configured without a GUI backend. Unfortunately ipyparallel and this particular GUI-less environment have a hard time interacting (for complicated reasons). We have derived a workaround that allows the parallelization to function. You should execute the following commands in a terminal on your cluster head node.

VERY IMPORTANT: This environment variable needs to be set in both .bashrc and .profile so that it is picked up when you run ipyparallel in either the head node of the cluster or on compute nodes.

.. code-block:: bash
    $ echo "# Prevent ipyparallel engines from dying in a headless environment" >> ~/.bashrc
    $ echo "export QT_QPA_PLATFORM=offscreen" >> ~/.bashrc
    $ echo "export QT_QPA_PLATFORM=offscreen" >> ~/.profile
    $ source ~/.bashrc
    $ source ~/.profile
    $ env | grep QT

Why is my structure analysis crashing when it looks like it should be working?
------------------------------------------------------------------------------
When running structure, specifically in the `get_clumpp_table` call, you might be told that "No files ready for XXX-K-2 in </your/structure/folder>", when in fact there are files ready. Well it turns out that CLUMPP has a 100 character file name limit, and it'll crash with names longer than this. The ipyrad.analysis.structure functions use absolute paths to specify file names, so it's not hard to see how this 100 character limit could be violated. Try moving your structure analysis to a place higher in the file system hierarchy. Baffling!

Problems with SRATools analysis package
---------------------------------------
Occasionlly with the sratools package you might have some trouble with downloading. It could look something like this `Exception in run() - index 29 is out of bounds for axis 0 with size 1`. This is a problem with your `esearch` install, which does not have https support built in. You can verify with with the command `esearch -db sra -query SRP021469`, which should give you an https protocol support error. You can easily fix this by installing `ssleay`: `conda install -c bioconda perl-net-ssleay`. Thanks to @ksil91.

ValueError in step 7
--------------------
During step 7 if you see something like this `error in filter_stacks on chunk 0: ValueError(zero-size array to reduction operation minimum which has no identity)` it means that one of your filtering parameters is filtering out all the loci. This is bad, obviously, and it's probably because one of your filtering parameters is too strict. Take a look at a couple of the samples in the *_consens directory to make sure they are reasonable, then try adjusting your filtering parameters based on how the consensus reads look.
