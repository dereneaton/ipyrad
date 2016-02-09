

.. _userguide:

User Guide
============
The most quick and dirty way to get started with ipyrad is to use the CLI
as described in the :ref:`CLI Quick Guide <quickguide_CLI>`. To get started
with writing scripts for assembly and/or analyses using the ipyrad API see the 


.. _CLI:

Tutorial - CLI
---------------------------------------
The ipyrad_ CLI_ (command line interface) is the easiest way to use ipyrad_, and 
will be familiar to users of its predecessor program pyrad_. 

.. code-block:: python

    ipyrad -h



Example Tutorials
------------------
The following tutorials show an example run through for an entire data set of 
each common data type, and explains some vagaries unique to each data type. 

* :ref:`Basic RAD <ipyrad_scripts.rst>`
* :ref:`Basic ddRAD <quickguide_CLI.rst>`
* Basic GBS  
* Basic paired ddRAD  
* Basic paired GBS  





.. _API:

Tutorial - API (Python)
-------------------------------------------
The ipyrad_ API_ (application program interface) is a way of directly accessing 
the Python_ functions that power ipyrad_. For users who are comfortable using 
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

