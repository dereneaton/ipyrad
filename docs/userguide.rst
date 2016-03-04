
.. include:: global.rst  

.. _userguide:


.. _tutorials:

Tutorials
=========
The most quick and dirty way to get started with ipyrad is to use the 
command line interface (CLI). A detailed walk-through of a RAD-seq assembly 
with explanations for each step, and the resulting data files that it creates, 
is available in the :ref:`introductory tutorial<tutorial_intro_cli>`. Once you're 
comfortable with that you should explore the 
:ref:`advanced tutorial<tutorial_advanced_cli>`
to learn how to take advantage of some tricks that ipyrad offers for 
more efficiently assembling data. In the coming weeks we plan to add more
tutorials and cookbook recipes for advanced usage. 


The command line interface (CLI)
---------------------------------
The ipyrad_ command line interface (CLI_) is designed to be easy to use and will
be familiar to users of its predecessor program pyrad_. Once ipyrad_ is installed
you can start using the CLI_ by simply opening a terminal and typing the 
following:

.. code-block:: bash

    ipyrad -h

This prints a help screen with a short description of the main arguments to the
CLI. There is actually a second way to use ipyrad_ separate from the CLI_ which
is to write scripts using the Python :ref:`API<API>`. Because the CLI_ in 
generally easier to use we focus on this for the main tutorials. 


.. _CLI:
Introductory tutorials 
----------------------
The following tutorials show an example run through of an entire data set 
to demonstrate either basic or advanced principles of using ipyrad. 
I recommend that anyone start by reading through the Introductory tutorial. 

.. toctree::
   :maxdepth: 1

   tutorial_intro_cli.rst
   tutorial_advanced_cli.rst
   tutorial_preview_mode.rst
   tutorial_intro_gbs.rst
   tutorial_paired_gbs.rst


.. _empirical_examples:

Examples with empirical data
-----------------------------
The following tutorials show example assemblies with publicly available 
empirical data sets. The first analyzes a small RAD-seq data set from 
Eaton and Ree (2013) using preview-mode in the CLI. 
The second tutorial analyzes the same data set using the API, 
and includes downstream visualization of data sharing among 
samples, and analysis of the results using the svd4tet species tree approach, 
and ABBA-BABA tests for introgression.


.. toctree::
    :maxdepth: 1
    pedicularis_cli.rst
    pedicularis_api.rst
    viburnum.rst


