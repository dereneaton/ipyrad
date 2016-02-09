.. include:: global.rst  

.. _analysis:


Analysis tools
==============

ipyrad_ includes a suite of analysis tools for comparing and visualizing
the size and completeness of assembled RAD-seq data sets, and also for 
calculating population-genetic statistics or performing several comparative
genomic analyses. These tools are in development and will be expanded with time.

The analysis tools are accessed through the ipyrad_ Python API. Thus, a bit of 
familiarity with Python programming will generally help to make these analyses
run more easily. However, we provide a suite of examples in the form of 
:ref:`cookcook recipes<cookbook_recipes>` which should make it easy to use.



Loading an Assembly into the API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If your data set was assembled using the ipyrad_ CLI then all the information 
for the entire assembly can be loaded from the `.assembly` file which will be
located in the `project dir`. In the example below we load an assembly named
`my_assembly` from a project dir called `testdir`. 


.. code-block:: python

    ## import ipyrad 
    import ipyrad as ip

    ## load in the complete assembly
    data1 = ip.load.load_assembly("testdir/my_assembly.json")


Loading other data into the API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Alternatively, if you want to analyze a data set that was not created in 
ipyrad_, but using some other software, it can be loaded in from either a 
.loci or .vcf formatted file. The analysis options are much more limited in 
this case (e.g., a lot of information needed for plotting is missing).


.. code-block:: python

    ## If your data set was assembled in an old version of pyrad
    ## the data can be loaded from a .loci file
    data1 = ip.load.load_loci("testdir/old_assembly.loci")

    ## or, if your data was assembled elsewhere it can be loaded from a vcf
    data1 = ip.load.load_vcf("testdir/somedata.vcf")



Plotting
^^^^^^^^
.. toctree::
   :maxdepth: 2

   plotting.rst


Introgression analyses
^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   dstats.rst


Population genetic analyses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   popgen.rst




