


.. include:: global.rst


.. _API:


The basics
----------
The API_ (application program interface) for ipyrad_ is a way of directly 
accessing the nuts and bolts of ipyrad_ using Python_. 
This has a number of advantages over the CLI_ in that there is a lot more 
flexibility for creating highly complex branching assemblies, or for applying
ipyrad_ in a non-standard way. It's best feature, though, is that you can 
perform entire analyses within Jupyter :ref:`notebooks<notebooks>` to create 
documented reproducible code for your analyses. 

Why use the API? 
^^^^^^^^^^^^^^^^^^^
The API provides a much more flexible framework for writing code to perform
complex branching assemblies than the CLI can provide. Because it is interactive
you can more easily access the results and statistics from each step. There are
two main functions of the API: 

* Assembly -- write scripts to perform the assembly 
* Analysis -- analyze and compare the size and distribution of data sets, 
create plots, calculate population genetic statistics, and perform phylogenetic
analyses. 

Before we jump into describing the API usage 

If you are already a pro with using Python then skip the next few sections
to get to the meat of using the :ref:`ipyrad API<ipyrad_API>`. 


Getting started with IPython/Jupyter notebooks
-----------------------------------------------
Our goal with using the ipyrad API is not only to get people writing Python 
scripts, but also to encourage the use of the a really exciting new tool called
Jupyter notebooks...
There may be a slight learning curve, however, for users who have no prior
experience with Python scripting. 

The envisioned usage of the ipyrad Python API is to run test assemblies within
a Jupyter notebook on a local computer using the **preview mode** method to 
execute quickly. Once the 


SVD4tet -- species tree inference
------------------------------------------------


Cookbook recipes - API 
-----------------------




.. toctree::
   :maxdepth: 1

   ipyrad_scripts.rst
   test_rad.rst

