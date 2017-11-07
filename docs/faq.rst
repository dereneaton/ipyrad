
.. _faq:  


FAQ
===

* wat

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

```
git clone https://github.com/dereneaton/ipyrad.git
conda build ipyrad/conda.recipe/bpp/
conda install --use-local bpp
```
