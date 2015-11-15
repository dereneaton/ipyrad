
.. _ Ethos:


How is it different from *pyrad*?
-------------------------------
ipyrad_ is a complete re-write of pyrad_ with an expanded focus on speed and flexibility.
While we continue in the minimalist ethos of pyrad_ which emphasized a simple
installation procedure and ease-of-use, ipyrad_ differs in offering an additional 
interactive interface through which to access data and results with simple Python scripts.
We continue to support a command line interface (CLI_) that will be familiar
to legacy pyrad_ users, but the real power of ipyrad_ comes from its
implementation as a Python module which allows users to design complex
assemblies that construct multiple data sets under multiple sets
of parameter settings; to directly access assembly statistics; to plot assembly results;
and to perform interactive downstream analyses.


Features
--------
Major new features and improvements include:

    - New assembly methods: *de novo*, reference alignment, 
         or hybrid (*de novo* & reference).
    - Parallel implementation using ipyparallel_ which utilizes MPI 
         allowing greater use of HPC clusters.
    - Better checkpointing. If your job is ever interrupted you should 
         be able to simply restart the
      script and continue from where it left off.
    - Faster code (speed comparisons forthcoming with publication).
    - Write highly reproducible documented code with Jupyter Notebooks (see Notebook_workflow_).
    - No external installations: vsearch, muscle and all other 
         dependencies are installed with ipyrad_ using conda (see Installation_).


.. include:: global.rst
