
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
