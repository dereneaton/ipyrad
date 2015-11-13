.. ipyrad documentation master file, created by
   sphinx-quickstart on Sun Sep  6 00:44:22 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: global.rst  


*ipyrad*: interactive assembly and analysis of RADseq data sets
-------------------------------------------------
Welcome to `*ipyrad*`_, an interactive toolkit for assembly and analysis of 
restriction-site associated genomic data sets (e.g., RAD, ddRAD, GBS) with 
the following goals:

- Provide an easy-to-use and intuitive workflow to convert raw data to formatted output data files.
- Provide a framework for producing complex assemblies within a simple and reproducible framework.
- Provide visualization and checks on the quality of data assemblies. 
- Provide interactive access to assembled data and statistics for downstream analyses and plotting.

Read more about our broader goals behind *ipyrad* here_. 

.. here_ Ethos.rst


What does _ipyrad_ do?
----------------------
`*ipyrad*`_ can be used to assemble RADseq data sets using either `*de novo* assembly`_, 
`reference mapping assembly`_, or `*hybrid assembly*`_ -- a combination
of the two approaches. Assembled data sets are output in a huge variety of 
`output formats`_ facilitating downstream genomic analyses for both population
genetic analyses as well as for phylogenetic studies of highly divergent species. 
*ipyrad* also includes methods for visualizing data and results, inferring 
population genetic statistics, and inferring genomic introgression. 


How is it different from *pyrad*?
-------------------------------

ipyrad_ is a complete re-write of pyrad_ with an expanded focus on speed and flexibility.
While we continue in the minimalist ethos of pyrad_ which emphasized a simple 
installation procedure and ease-of-use, ipyrad_ offers an additional interactive 
interface with which to access data and results through simple Python scripts. 
We continue to support a command line interface (CLI_) that will be familiar
to legacy pyrad_ users, but the real power of ipyrad_ comes from its 
implementation as a Python module which allows users to design complex 
assemblies that construct multiple data sets under multiple sets
of parameter settings; to directly access assembly statistics; to plot assembly results;
and to perform interactive downstream analyses. 


Major new features and improvements include: 

    - New assembly methods: *de novo*, reference alignment, or hybrid (*de novo* & reference).  
    - Parallel implementation using ipyparallel_ which utilizes MPI allowing greater use of HPC clusters.   
    - Better checkpointing. If your job is ever interrupted you should be able to simply restart the
      script and continue from where it left off.  
    - Faster code (speed comparisons forthcoming with publication).  
    - Write highly reproducible documented code with Jupyter Notebooks (see Notebook_workflow_).    
    - No external installations: vsearch, muscle and all other dependencies are installed with ipyrad_ 
      using conda (see Installation_).   


Documentation
-------------

.. toctree::
   :maxdepth: 2

   Ethos.rst
   Features.rst
   Installation.rst
   Quick-guide.rst
   Assembly.rst
   Tutorials.rst
   Citing.rst
   License.rst
   Contributions.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

