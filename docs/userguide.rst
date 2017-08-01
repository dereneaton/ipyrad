
.. include:: global.rst  

.. _userguide:


.. _tutorials:

Tutorials (running ipyrad)
==========================

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
produced by each step of an assembly. Next, try some advanced methods, 
like using branching to assemble data sets under a range of parameter settings, 
and assemble data with respect to a reference genome.

   + `Introductory tutorial to the ipyrad CLI <tutorial_intro_cli.html>`__ (command line)  
   + `Advanced tutorial to the ipyrad CLI <tutorial_advanced_cli.html>`__ (command line)  
   + `Running on ipyrad on HPC <HPC_script.html>`__  
   + `Introductory tutorial to the ipyrad API <API_user-guide.html>`__  


.. _empirical_example:

Empirical examples
-----------------
The following tutorials demonstrate assemblies of publicly available 
empirical data sets representing different data types. The first analysis 
(Eaton and Ree, 2013) can be assembled very quickly, and is re-used in our 
analysis cookbook recipes, below. The others include tips for optimizing 
ipyrad for use with that data type. 

   + `Pedicularis CLI <pedicularis_.html>`__ (command line)  
   + `Pedicularis API <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-1-pedicularis.ipynb>`__ (run in jupyter-notebook)  
   + `Finch API <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-2-Finches.ipynb>`__ (run in jupyter-notebook)  
   + `see more ipyrad API cookbooks here <http://ipyrad.readthedocs.io/analysis.html>`__  

    .. pedicularis_api.rst
    .. viburnum.rst
    .. cleaning_up_pairs.rst
    .. comment out

