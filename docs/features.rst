
.. include:: global.rst  

.. _features:


Features
========

What does ipyrad do?
--------------------
ipyrad_ is a toolbox for assembly and analysis of RAD-seq type genomic data sets. 
Notably, it has four :ref:`assembly methods<assembly_methods>` 
by which to assemble data: denovo, reference, 
reference addition, and reference subtraction. Assembled data sets are created
in a variety of :ref:`output formats<full_output_formats>`, facilitating downstream 
genomic analyses for both population genetic and phylogenetic studies. ipyrad_ 
also includes methods for visualizing and analyzing data and results. 

How is it different from pyrad?
-------------------------------
ipyrad_ is a complete re-write of pyrad_ with an expanded focus on speed and 
flexibility. While we continue in the minimalist :ref:`ethos<ethos>` of pyrad_, 
which emphasized a simple installation procedure and ease-of-use, ipyrad_ offers
many new features, and is now easily extensible to new data types and models 
through its Python :ref:`API<API>`. 


Major New Features in ipyrad
----------------------
* New :ref:`assembly methods<assembly_methods>`: *de novo* and reference-based methods.
* Improved :ref:`checkpointing<checkpointing>`. Interrupted jobs are easily restarted.
* Much faster code (speed comparisons forthcoming with publication).
* State of the art parallel implementation (ipyparallel) for running on computing clusters. 
* Easy installation: all dependencies are included during :ref:`installation<installation>`.


Coming Soon
-----------
* Downstream analysis tools
* Quartet-based species tree inference (_tetrad_ program)
* Introgression analyses (ABBA-BABA tests) 
* Population genetic statistics 

