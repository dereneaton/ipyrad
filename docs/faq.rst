
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
on an HPC compute node.

.. code-block:: bash
    qsub -I

.. code-block:: bash
    ipcluster start --n 4 --daemonize

Then type `ipython` top open an ipython session.

.. code-block:: python

    import ipyparallel as ipp

    rc = ipp.Client()
    rc[:]

The result should look something like this:
.. parsed_literal::

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

.. parsed_literal::

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

Other random problems with solutions
------------------------------------
.. parsed_literal::

    Failed at nopython (nopython frontend)
    UntypedAttributeError: Unknown attribute "any" of type Module(<module 'numpy' from...

In some instances if you already have conda/python installed the local environment
variable PYTHONPATH will be set, causing python to use versions of modules 
outside the miniconda path set during ipyrad installation. This error can be fixed by 
blanking the PYTHONPATH variable during execution (as below), or by adding the export
to your ~/.bashrc file.

.. code-block:: bash

    export PYTHONPATH=""; ipyrad -p params.txt -s 1

