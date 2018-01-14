.. include:: global.rst

.. _analysis:


The ipyrad analysis toolkit
===========================

ipyrad_ includes a suite of analysis tools that are designed to make it 
easy to run inference programs (e.g., STRUCTURE, BUCKy, BPP) on the data
in an efficient way by sampling distributions of loci or SNPs from your RAD
data, grouping individuals into populations, filtering for missing data, 
and parallelizing computation. 

In this section of the documentation we have a number of example analyses
in the form of Jupyter notebooks, which is a useful tool for doing 
reproducible science. In fact, ipyrad has been designed since its inception 
for the goal of working seamlessly within jupyter. Check out the tutorials 
below on using Jupyter notebooks, and using ipyrad in notebooks. 
Then check out the analysis tools notebooks. 


Jupyter notebooks 
^^^^^^^^^^^^^^^^^^
This is an optional tool to use with ipyrad, but one that we strongly 
recommend learning. See the video and link below to learn about notebooks and
how to run them locally or on an HPC cluster. 

+ `Intro to Jupyter Notebooks <https://www.youtube.com/watch?v=HW29067qVWk&t=47s>`__ (Video)
+ `Running jupyter on a HPC cluster <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__ 
.. + `Intro to ipyparallel (the parallel client used by ipyrad) <http://ipyrad.readthedocs.io/ipyparallel.html>`__
.. + `Using ipyrad with jupyter and git <http://ipyrad.readthedocs.io/Git_Sync.html>`__>`__


*ipyrad* API full example notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These notebooks show example usage of the ipyrad API.

+ `Pedicularis API <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-1-pedicularis.ipynb>`__ (run in jupyter-notebook)
+ `Finch API <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-2-Finches.ipynb>`__ (run in jupyter-notebook)  


*ipyrad* API analysis tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^
These notebook show how to do parallelized downstream analyses in Jupyter-notebooks, and to generate advanced input files for many programs using the ipyrad analysis tools. 

+ `TETRAD quartet species tree inference <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-quartet-species-tree.ipynb>`__
+ `RAxML concatenation tree inference <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-raxml-pedicularis.ipynb>`__
+ `BPP species tree and delimitation <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-bpp-species-delimitation.ipynb>`__
+ `TREEMIX admixture graph inference <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-treemix-pedicularis.ipynb>`__
+ `BUCKY concordance tree inference <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-bucky.ipynb>`__
+ `STRUCTURE population structure <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-structure-pedicularis.ipynb>`__
+ `STRUCTURE with pop assignments <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/example-structure-popdata.ipynb>`__
+ `ABBA BABA admixture <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-abba-baba.ipynb>`__


*command line programs*
^^^^^^^^^^^^^^^^^^^^^^^
These pages discuss further information about some command-line analysis tools 
that are frequently used with RAD-seq data.  

+ `RAxML phylogenetic inference CLI <http://ipyrad.readthedocs.io/raxml.html>`__
+ `TETRAD command line <http://ipyrad.readthedocs.io/tetrad.html>`__   


Other jupyter notebooks (ipyrad in the wild)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
+ You can contribute here. Let us know.
