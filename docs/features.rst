
.. _features:


Features
========


What does ipyrad do?
----------------------
:ref:`ipyrad <ipyrad>` can be used to assemble RADseq data sets using 
:ref:`*de novo* assembly <denovo>`, :ref:`reference mapping assembly <reference>`, 
or `:ref:`hybrid assembly <hybrid>` -- a combination of the two approaches. 
Assembled data sets can be output in a large variety of :ref:`formats <formats>`, facilitating downstream genomic analyses for both population genetic and phylogenetic 
studies. It also includes methods for visualizing data and results, inferring population genetic statistics, and inferring genomic introgression.


How is it different from pyrad?
-------------------------------
:ref:`ipyrad <ipyrad>` is a complete re-write of :ref:`pyrad <pyrad>` with 
an expanded focus on speed and flexibility. While we continue in the minimalist 
:ref:`ethos <ethos>` of :ref:`pyrad <pyrad>` which emphasizes a simple installation procedure and ease-of-use, :ref:`ipyrad <ipyrad>` offers many new features through 
its Python API. However, we also continue to support the command line interface (:ref:`CLI <CLI>`) that will be familiar to :ref:`pyrad <pyrad>` users.


Major New Features in ipyrad
----------------------
* New assembly methods: :ref:`de novo <denovo>`, :ref:`reference mapping <refalign>`, or :ref:`hybrid <hybrid>`.
* New parallel implementation that can utilize MPI, PBS, and LSF to full utilize HPC clusters.
* Better checkpointing. Interrupted jobs can be easily restarted to continue from where they left off.
* Faster code (speed comparisons forthcoming with publication).
* Write highly reproducible documented code with Jupyter Notebooks (see :ref:`Notebooks <notebooks>`).
* No external installations: all other dependencies are installed with conda (see :ref:`installation <installation>`).


Coming Soon
-----------
* Introgression analyses (ABBA-BABA tests) 
* Population genetic statistics 