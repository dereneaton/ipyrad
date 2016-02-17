
.. include:: global.rst  

.. _userguide:

User Guide
============
The most quick and dirty way to get started with ipyrad is to use the CLI
as described in the :ref:`CLI Quick Guide <quickguide_CLI>`. To get started
with writing scripts for assembly and/or analyses using the ipyrad API see the 


.. _CLI:

Tutorial - CLI
---------------------------------------
The ipyrad_ command line interface (CLI_) is designed to be easy to use and will
be familiar to users of its predecessor program pyrad_. First things first, once
ipyrad_ is installed open a terminal and type:

.. code-block:: bash

    ipyrad -h

and the following help screen will appear:

.. code-block:: none

    usage: ipyrad [-h] [-v] [-r] [-n] [-f] [-p params] [-s steps] [-c cores]
                  [--MPI] [--preview]

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit
      -r, --results  show summary of results for Assembly in params.txt and exit
      -f, --force    force overwrite of existing data
      -n new         create new params file with name 'new' in current directory
      -p params      path to params file (e.g., new-params.txt)
      -s steps       subset of assembly steps to perform (Default=1234567)
      -c cores       number of CPU cores to use (Default=all)
      --MPI          use MPI to connect to parallel CPUs across multiple nodes
      --preview      run ipyrad in preview mode. Subset the input file so it'll
                     run quickly so you can verify everything is working

      * Example command-line usage: 
        ipyrad -n test                     ## create new file test-params.txt.
        ipyrad -p test-params.txt          ## run ipyrad with settings in params.
        ipyrad -p test-params.txt -s 123   ## run only steps 1, 2 and 3 of assembly.
        ipyrad -p test-params.txt -s 4567  ## run steps 4, 5, 6 and 7 of assembly.
        ipyrad -p test-params.txt -s 3 -f  ## run step 3, overwrite existing data.

      * HPC parallelization
        ipyrad -p test-params.txt -s 3 --MPI   ## access cores across multiple nodes

      * Results summary quick view
        ipyrad -p test-params.txt -r         ## print summary stats to screen for params.

      * Documentation: http://ipyrad.readthedocs.org/en/latest/




Example Tutorials
------------------
The following tutorials show an example run through for an entire data set of 
each common data type, and explains some vagaries unique to each data type. 

* :ref:`Full Tutorial (RAD-Seq) <full_tutorial_CLI.rst>`

* :ref:`Basic RAD <ipyrad_scripts.rst>`
* :ref:`Basic ddRAD <quickguide_CLI.rst>`
* Basic GBS  
* Basic paired ddRAD  
* Basic paired GBS  

.. toctree::                                                                                                                                         
   :maxdepth: 2                                                                                                                                      
                                                                                                                                                     
   advanced_CLI.rst  




.. _API:

Tutorial - API (Python)
-------------------------------------------
The API_ (application program interface) for ipyrad_ is a way of directly 
accessing the nuts and bolts of ipyrad_ using Python_. 
This has a number of advantages over the CLI_ in that there is a lot more 
flexibility for creating highly complex branching assemblies, or for applying
ipyrad_ in a non-standard way. It's best feature, though, is that you can 
perform entire analyses within Jupyter :ref:`notebooks<notebooks>` to create 
documented reproducible code for your analyses. 

There may be a slight learning curve, however, for users who have no prior
experience with Python scripting. 

Python it provides a much more flexible framework for writing code to perform
complex branching assemblies than the CLI can provide. Because it is interactive
you can more easily access the results and statistics from each step. There are
two main functions of the API: 

* Assembly -- write scripts to perform the assembly 
* Analysis -- analyze and compare the size and distribution of data sets

The envisioned usage of the ipyrad Python API is to run test assemblies within
a Jupyter notebook on a local computer using the **preview mode** method to 
execute quickly. Once the 


Cookbook recipes
----------------

:ref:`API Quick Guide <quickguide_API>`. 
:ref:`ipyrad_scripts.rst <ipyrad_scripts.rst>`
:ref:`test_rad.rst <test_rad.rst>`


.. toctree::
   :maxdepth: 2

   quickguide_CLI.rst
   quickguide_API.rst
   ipyrad_scripts.rst
   test_rad.rst

