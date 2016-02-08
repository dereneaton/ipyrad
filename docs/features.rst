
.. include:: global.rst  

.. _features:


Features
========

What does ipyrad do?
--------------------
ipyrad_ is a toolbox for assembly and analysis of RAD-seq type genomic data sets. 
Notably, it has four assembly_methods_ by which to assemble data: denovo, reference, 
reference addition, and reference subtraction. Assembled data sets are created
in variety of output_formats_, facilitating downstream genomic analyses for both 
population genetic and phylogenetic studies. ipyrad_ also includes methods for 
visualizing_ and analyzing_ data and results. 

How is it different from pyrad?
-------------------------------
ipyrad_ is a complete re-write of pyrad_ with an expanded focus on speed and 
flexibility. While we continue in the minimalist ethos_ of pyrad_, which emphasized
a simple installation procedure and ease-of-use, ipyrad_ offers many new 
features, particularly through its Python API_. However, we also continue to 
support the command line interface CLI_ that will be familiar pyrad_ users.


Major New Features in ipyrad
----------------------
* New assembly methods: :ref:`de novo <denovo>`, :ref:`reference mapping <refalign>`, :ref:`reference_sub <reference_sub>`, :ref:`reference_add <reference_add>`.
* New parallel implementation that can utilize MPI, PBS, and LSF to full utilize HPC clusters.
* Better checkpointing. Interrupted jobs can be easily restarted to continue from where they left off.
* Faster code (speed comparisons forthcoming with publication).
* Write highly reproducible documented code with Jupyter Notebooks (see :ref:`Notebooks <notebooks>`).
* No external installations: all other dependencies are installed with conda (see :ref:`installation <installation>`).


Coming Soon
-----------
* Introgression analyses (ABBA-BABA tests) 
* Population genetic statistics 

