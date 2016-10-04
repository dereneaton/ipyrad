
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
Introductory tutorials
---------------------
Start here to learn the basics. We run through an example simulated single-end 
RAD-seq data set and give detailed descriptions of files and statistics 
produced by each step of an assembly. 

.. toctree::
   :maxdepth: 1

   tutorial_intro_cli.rst


.. .. _advanced_tutorial:
.. Advanced tutorials
.. ----------------------
Next, try some advanced methods, like using branching to
assemble data sets under a range of parameter settings. Here we assemble a
data set using both denovo and reference assembly methods.

.. toctree::
   :maxdepth: 1

   tutorial_advanced_cli.rst


.. _empirical_example:

Empirical examples
-----------------
The following tutorials demonstrate assemblies of publicly available 
empirical data sets representing different data types. 
The first analysis (Eaton and Ree, 2013) can be assembled very quickly, 
and is re-used in our analysis cookbook recipes, below. The others include
tips for optimizing ipyrad for use with that data type. 

.. toctree::
    :maxdepth: 1

    pedicularis_.rst
    pairddrad_.rst
    pairgbs_.rst
    xxxpedicularis_cli.rst


Cookbook recipes
-----------------
Many more tutorials coming soon. 


.. toctree::
    :maxdepth: 1

    HPC_script.rst
    pedicularis_api.rst
    viburnum.rst
    cleaning_up_pairs.rst

