.. ipyrad documentation master file, created by
   sphinx-quickstart on Sun Sep  6 00:44:22 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ipyrad: assembly and analysis of RADseq data sets
==================================================

Welcome!
========

ipyrad_ is an interactive toolkit for assembly and analysis of genomic RADseq data sets.
Our goal is to support all restriction-site associated data types (e.g., RAD, ddRAD, GBS;
see Data_types_). 


How is it different from pyrad?
-------------------------------

ipyrad_ is a complete re-write of pyrad_ built with a very different philosophy in mind. 
While it retains the easy-to-use command line interface (CLI_) that will be familiar to pyrad_ users,
the real power of ipyrad_ comes from its implementation through a Python API, which allows users to 
write scripts that detail complex assemblies able to construct multiple data sets under multiple 
parameter settings. Other improvements include: 

    - 3 modes of assembly: *de novo*, reference alignment, or hybrid (*de novo* & reference).  
    - Parallel implementation using ipyparallel_ which utilizes MPI allowing greater use of HPC clusters.   
    - Better checkpointing. If your job is ever interrupted you should be able to simply restart the
      script and continue from where it left off.  
    - Faster code (speed comparisons forthcoming with publication).  
    - Write highly reproducible documented code with Jupyter Notebooks (see Notebook_workflow_).    
    - No external installations: vsearch, muscle and all other dependencies are installed with ipyrad_ 
      using conda (see Installation_).   


Documentation
=============

.. toctree::
   :maxdepth: 2

   Ethos.rst
   Installation.rst
   Command_line_interface.rst
   ipyrad_scripts.rst
   Notebook_workflow.rst
   Tutorials.rst
   Data_types.rst
   Citing.rst
   License.rst
   Contributions.rst
   Dependencies.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

