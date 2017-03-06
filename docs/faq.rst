
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

Why doesn't ipyrad handle PE original RAD
-----------------------------------------
Paired-End RAD protocol is tricky to denovo assemble. Because of the sonication step R2 
doesn't line up nicely. ipyrad makes strong assumptions about how r1 and r2 align, 
assumptions which are met by PE gbs and ddrad, but which are not met by original RAD. 
This doesn't matter (as much) if you have a reference genome, but if you don't have a 
reference it's a nightmare... dDocent has a PE-RAD mode, but I haven't evaluated it. 
I know that people have also used stacks (because stacks treats r1 andr2 as independent 
loci). If people ask me how to denovo assemble with PE-RAD in ipyrad I tell them to 
just assemble it as SE and ignore R2.
