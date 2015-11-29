
.. _features:


Features
========


What does ipyrad do?
----------------------
:ref:`ipyrad <ipyrad>` can be used to assemble RADseq data sets using 
:ref:`*de novo* assembly <denovo>`, :ref:`reference mapping assembly <reference>`, 
or `:ref:`hybrid assembly <hybrid>` -- a combination of the two approaches. 
Assembled data sets can be output in a large variety of `formats`_, facilitating downstream genomic analyses for both population genetic and phylogenetic 
studies. It also includes methods for visualizing data and results, inferring population genetic statistics, and inferring genomic introgression.


How is it different from pyrad?
-------------------------------
:ref:`ipyrad <ipyrad>` is a complete re-write of :ref:`pyrad <pyrad>` with 
an expanded focus on speed and flexibility. While we continue in the minimalist 
:ref:`ethos <ethos>` of :ref:`pyrad <pyrad>` which emphasizes a simple installation procedure and ease-of-use, :ref:`ipyrad <ipyrad>` offers a Python API which allows users to interactively access data and results 

through simple Python scripts. We continue to support a command line interface (CLI_) that will be familiar to legacy pyrad_ users, but the real power of ipyrad_ comes from its implementation as a Python module which allows users to design complex
assemblies that construct multiple data sets under multiple sets
of parameter settings; to directly access assembly statistics; to plot assembly results;
and to perform interactive downstream analyses.


New Features
------------
Major new features and improvements include:

    - New assembly methods: *de novo*, reference alignment, or hybrid (*de novo* & reference).
    - Parallel implementation using ipyparallel_ which utilizes MPI allowing greater use of HPC clusters.
    - Better checkpointing. If your job is ever interrupted you should be able to simply restart the
      script and continue from where it left off.
    - Faster code (speed comparisons forthcoming with publication).
    - Write highly reproducible documented code with Jupyter Notebooks (see Notebook_workflow_).
    - No external installations: vsearch, muscle and all other dependencies are installed with ipyrad_
      using conda (see Installation_).

