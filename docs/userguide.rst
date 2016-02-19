
.. include:: global.rst  

.. _userguide:


Tutorials
=========
The most quick and dirty way to get started with ipyrad is to use the CLI
as described in the :ref:`Quick Guide<quickguide_CLI>`. A much more detailed
walk-through of a RAD-seq assembly with explanations for each step in available
in the :ref:`introductory tutorial<tutorial_intro_cli>`. And once you're 
comfortable with that you should explore the :ref:`advanced tutorial<tutorial_advanced_cli>`
to learn how to really take advantage of some tricks that ipyrad offers for 
really efficiently assembling data. We also provide introductory guides for 
several common :ref:`data types<data types>` which have a few differences each
that users should be aware of. 


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
is to write scripts using the Python API_, which is described further below. 
Because the CLI_ in generally easier the primary tutorials focus on this 
interface. 


.. _CLI:

Introductory tutorials 
----------------------
The following tutorials show an example run through for an entire data set of 
each common data type, and explains some vagaries unique to each data type. I 
recommend that everyone starts by reading through the Introductory tutorial. 

.. toctree::
   :maxdepth: 1

   tutorial_intro_cli.rst ``(start here)``
   tutorial_advanced_cli.rst
   tutorial_intro_gbs.rst
   tutorial_paired_gbs.rst


