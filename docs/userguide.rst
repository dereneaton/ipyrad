
.. include:: global.rst  

.. _userguide:


.. _tutorials:

Tutorials
=========

The ipyrad_ command line interface (CLI) is accessed through a terminal. Use 
the -help (-h) flag to print a help screen with a description of the main 
arguments to the CLI. Detailed instructions are available through the tutorials
below. 

.. code-block:: bash

    >>> ipyrad -h


.. _introductory_tutorial:
Introductory tutorial
---------------------
To learn the basics this is the best place to start. We run through an example
simulated single-end RAD-seq data set to give detailed descriptions of files 
and statistics produced by each step of an assembly. 

.. toctree::
   :maxdepth: 1

   tutorial_intro_cli.rst


.. _advanced_tutorial:
Advanced tutorials
----------------------
An introduction to advanced methods in ipyrad, including branching methods 
for assembling data sets under a range of parameter settings, and an example of 
both denovo and reference assembly methods. 

.. toctree::
   :maxdepth: 1

   tutorial_advanced_cli.rst


.. _empirical_example:

Empirical example
-----------------
The following tutorials show example assemblies with publicly available 
empirical data sets. The first analyzes a small RAD-seq data set from 
Eaton and Ree (2013). 


.. toctree::
    :maxdepth: 1

    pedicularis_.rst
    xxxpedicularis_cli.rst



Cookbook recipes
-----------------
Many more tutorials coming soon. 


.. toctree::
    :maxdepth: 1

    HPC_script.rst
    pedicularis_api.rst
    viburnum.rst
