.. include:: global.rst  

.. _analysis:


Analysis tools
==============

ipyrad_ includes a suite of analysis tools for comparing and visualizing
the size and completeness of assembled RAD-seq data sets, and also for 
calculating population-genetic statistics or performing several comparative
genomic analyses. Many of these tools are still in development and will 
be expanded with time.

Some of the tools are available as separate command line programs (e.g., svd4tet)
while others are accessed through the ipyrad_ Python API (e.g., plotting functions).
For the latter, a bit of familiarity with Python programming will generally help. 
However, we are working to develop a suite of examples in the form of 
:ref:`cookbook recipes<cookbook_recipes>` which should make it easy to use.


Loading an Assembly into the API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If your data set was assembled using the ipyrad_ CLI then all the information 
for the entire assembly can be loaded from the `.json` file which will be
located in the `project dir`. In the example below we load an assembly named
`my_assembly` from a project dir called `testdir`. 


.. code-block:: python

    ## import ipyrad 
    import ipyrad as ip

    ## load in the complete assembly
    data1 = ip.load_json("testdir/my_assembly.json")

    ## downstream analyses and plotting examples
    ## coming soon...


Loading other data into the API (coming soon)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Alternatively, if you want to analyze a data set that was not created in 
ipyrad_, but using some other software, it can be loaded in from either a 
.loci or .vcf formatted file. The analysis options are much more limited in 
this case (e.g., a lot of information needed for plotting is missing).


.. .. code-block:: python

..     ## If your data set was assembled in an old version of pyrad
..     ## the data can be loaded from a .loci file
..     data1 = ip.load.load_loci("testdir/old_assembly.loci")

..     ## or, if your data was assembled elsewhere it can be loaded from a vcf
..     data1 = ip.load.load_vcf("testdir/somedata.vcf")



Plotting (in development)
^^^^^^^^
.. toctree::
   :maxdepth: 2

   plotting.rst


Introgression analyses (in development)
^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   dstats.rst


Population genetic analyses (in development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   popgen.rst


SVD4tet -- species tree inference (in development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Infer species trees from quartets inferred from SNP data. We have developed
an implementation of the SVDquartets algorithm of 
:ref:`Chiffman and Kubatko (2014)<svdquartets>` similar to the one 
currently implemented in Paup*. However, our implementation differs 
in several respects, with a focus on making use of RAD-seq data, and of 
efficiently sampling quartets over both small and very large trees (>200 tips). 
The analysis can be massively parallelized using MPI, and allows checkpointing
so that analyses can be stopped and restarted at a later time. 
We include three sampling schemes: 'full sampling', 'random sampling', and 
a novel approach for sampling from very large trees called 'even sampling'. 

Quartet joining is performed in the qQMC software both with and without
quartet weights (see :ref:`Avni et al. 2014<weighted_quartets>`)

The command line program svd4tet is installed alongside ipyrad and can be 
called from the command line. 

Link to a notebook tutorial for svd4tet and also to further cookbook recipes. 