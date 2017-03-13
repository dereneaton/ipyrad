.. include:: global.rst

.. _analysis:


Analysis tools
==============

ipyrad_ includes a suite of analysis tools for comparing and visualizing
the size and completeness of assembled RAD-seq data sets, for
calculating population-genetic statistics, and for performing comparative
genomic analyses. In addition, there of course countless other tools available
for the downstream analysis of RAD-seq data. In this section of the documentation
we hope to provide many examples to guide users through such analyses. Many
examples are shared in the form of Jupyter Notebooks, which are a useful tool
for doing reproducible science, and our first tutorial provides a crash course
if using Jupyter with ipyrad. You can follow along and use these tutorials 
without Jupyter as well. 


Using Jupyter notebooks 
^^^^^^^^^^^^^^^^^^^^^^^^^
This is an optional tool to use with ipyrad.
+ `Intro to Jupyter Notebooks Video <https://www.youtube.com/watch?v=HW29067qVWk&t=47s>`__  
+ `Jupyter Git and ipyrad <...>`__  
+ `setup SSH Tunneling to a HPC Cluster <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__


*ipyrad* API Cookbooks 
^^^^^^^^^^^^^^^^^^^^^^^^^^
These notebooks show example usage of the ipyrad API.

+ Intro to the ipyrad API <tutorial_API.html>`__
+ `Quantify and plot shared RAD data with ipyrad <visualize.html>`__
+ `Example RAD assembly with ipyrad API <http://nbviewer.jupyter.org/github/dereneaton/ficus-rad/blob/master/Ficus_Jander_assembly.ipynb>`__
+ `Example PE GBS assembly with ipyrad API <http://nbviewer.jupyter.org/github/dereneaton/pedicularis-WB-GBS/blob/master/nb1-WB-assembly.ipynb>`__

Downstream analysis tools
^^^^^^^^^^^^^^^^^^^^^^^^^^
These notebook show how to do parallelized downstream analyses in Jupyter-notebooks.

+ `phylogenetic inference with RAxML <http://ipyrad.readthedocs.io/raxml.html>`__
+ `concordance tree inference with BUCKy (parallel) <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-bucky.ipynb>`__
+ `STRUCTURE analysis (parallel) <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-structure-pedicularis.ipynb>`__
+ `quartet species tree inference with tetrad <http://ipyrad.readthedocs.io/tetrad.html>`__
+ ABBA BABA test for introgression <cookbook-dstats>`__
+ PCA analysis of genetic structure <cookbook-PCA>`__
