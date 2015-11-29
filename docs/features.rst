
.. _features:


Features
========


What does ipyrad do?
----------------------
:ref:`ipyrad <ipyrad>` can be used to assemble RADseq data sets using 
:ref:`*de novo* assembly <denovo>`, :ref:`reference mapping assembly <reference>`, 
or `:ref:`hybrid assembly <hybrid>` -- a combination of the two approaches. 
Assembled data sets can be output in a large variety of `formats`, facilitating downstream genomic analyses for both population genetic and phylogenetic 
studies. It also includes methods for visualizing data and results, inferring population genetic statistics, and inferring genomic introgression.


How is it different from pyrad?
-------------------------------
:ref:`ipyrad <ipyrad>` is a complete re-write of :ref:`pyrad <pyrad>` with 
an expanded focus on speed and flexibility. While we continue in the minimalist 
:ref:`ethos <ethos>` of :ref:`pyrad <pyrad>` which emphasizes a simple installation procedure and ease-of-use, :ref:`ipyrad <ipyrad>` offers an additional Python API
with which to interactively access data and results, as well as new methods for 
assembly and visualization of results, and more rigorous code testing. 

We continue to support the command line interface (:ref:`CLI <CLI>`) that 
is familiar to :ref:`pyrad <pyrad>` users, but the real strength 
:ref:`ipyrad <ipyrad>` comes from writing Python scripts. We provide many 
example scripts with the hope of promoting the use of complex but easily 
understandable code to construct multiple data sets under multiple sets
of parameter settings. 


Major New Features in ipyrad
----------------------

* New assembly methods: *de novo*, reference alignment, or hybrid (*de novo* & reference).
* Parallel implementation using :ref:`ipyparallel <ipyparallel>` which can utilize MPI, PSB, and LSF to make greater use of HPC clusters.
* Better checkpointing. If your job is ever interrupted you should be able to simply restart the script and continue from where it left off.
* Faster code (speed comparisons forthcoming with publication).
* Write highly reproducible documented code with Jupyter Notebooks (see :ref:`Notebook_workflow <notebooks>`).
* No external installations: vsearch, muscle and all other dependencies are installed  using conda (see :ref:`installation <installation>`).


Coming Soon
-----------
* Introgression analyses (ABBA-BABA tests) 
* Population genetic statistics 