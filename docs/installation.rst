
.. _installation:

Installation
============

.. toctree::
   :maxdepth: 2

   macports-installation.rst

Using Pip / Easy Install
------------------------

The easiest way to install *ipyrad* and all of its dependencies is 
to use *conda*, which is the current standard for installing Python
packages. Follow the very simple instructions to install *Anaconda*
or *Miniconda* for Python2.7 here_. Once installed you can use *conda* 
to install additional Python packages including *ipyrad* with commands 
like below. 

.. _here: http://conda.pydata.org/docs/install/quick.html

    $ conda update 
    $ conda install ipyrad  


In contrast to its predecessor (*pyrad*), *ipryad* makes use of many more
Python resources which are easily bundled with the installation when *conda*
is used. These include the following: 

- Numpy -- Scientific processing  
- Scipy -- Scientific processing  
- Pandas -- Used for manipulating data frames  
- Sphinx -- Used for building documentation  
- IPython -- Interactivity  
- ipyparallel -- Parallel, threading, MPI support
- jupyter -- Creating reproducible notebooks  
- H5py -- Database structure for large data sets  
- Dill -- Store pickle objects of complex Classes  
- Toyplot -- [optional]...   

Installing on HPC
-----------------
One of the greatest strengths of using *conda* for installation is that it 
creates a Python package directory in your home directory called ``~/anaconda/``
where new packages are installed, and because they are not stored in a 
system-wide directory you do not need administrator privileges to install 
new packages.  








