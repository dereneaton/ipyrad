
.. include:: global.rst  

.. _userguide:


Outline
=======
The most quick and dirty way to get started with ipyrad is to use the CLI
as described in the :ref:`Quick Guide<quickguide_CLI>`. A much more detailed
walk-through of a RAD-seq assembly with explanations for each step in available
in the :ref:`introductory tutorial<tutorial_intro_cli>`. And once you're 
comfortable with that you should explore the :ref:`advanced tutorial<tutorial_advanced_cli>`
to learn how to really take advantage of some tricks that ipyrad offers for 
really efficiently assembling data. We also provide introductory guides for 
several common :ref:`data types<data types>` which have a few differences each
that users should be aware of. 


.. _quickguide_CLI:

Quickguide -- command line interface
-----------------------------------
The ipyrad_ command line interface (CLI_) is designed to be easy to use and will
be familiar to users of its predecessor program pyrad_. First things first, once
ipyrad_ is installed open a terminal and type:

.. code-block:: bash

    ipyrad -h

and a help screen will appear with a short description of the arguments to the 
command line. 



.. _CLI:

Tutorials -- command line interace
-----------------------------------
The following tutorials show an example run through for an entire data set of 
each common data type, and explains some vagaries unique to each data type. 

.. toctree::
   :maxdepth: 2

   tutorial_intro_cli.rst
   tutorial_advanced_cli.rst




.. _API:

Tutorials - API 
-------------------------------------------
The API_ (application program interface) for ipyrad_ is a way of directly 
accessing the nuts and bolts of ipyrad_ using Python_. 
This has a number of advantages over the CLI_ in that there is a lot more 
flexibility for creating highly complex branching assemblies, or for applying
ipyrad_ in a non-standard way. It's best feature, though, is that you can 
perform entire analyses within Jupyter :ref:`notebooks<notebooks>` to create 
documented reproducible code for your analyses. 

There may be a slight learning curve, however, for users who have no prior
experience with Python scripting. 

Python it provides a much more flexible framework for writing code to perform
complex branching assemblies than the CLI can provide. Because it is interactive
you can more easily access the results and statistics from each step. There are
two main functions of the API: 

* Assembly -- write scripts to perform the assembly 
* Analysis -- analyze and compare the size and distribution of data sets

The envisioned usage of the ipyrad Python API is to run test assemblies within
a Jupyter notebook on a local computer using the **preview mode** method to 
execute quickly. Once the 


Cookbook recipes - API 
-----------------------


.. toctree::
   :maxdepth: 1

   ipyrad_scripts.rst
   test_rad.rst

