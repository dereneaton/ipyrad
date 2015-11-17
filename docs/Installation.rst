
.. include:: global.rst

.. _installation:


Installation with conda
-----------------------

The easiest way to install *ipyrad* and all of its dependencies is 
to use conda_, which is a command line program for installing Python
packages. If you do not have *conda* installed, follow these 
instructions_ to install either *Anaconda* or *Miniconda* 
for Python2.7. If you're working on an HPC system you can install
*conda* in your home directory even without administrative privileges. 

Either installation wil create a command line program called *conda* 
which can be used to install Python packages. The main difference between the
two is that *Anaconda* will also install a large suite of commonly used Python 
packages along with it, whereas *Miniconda* is a bare bones version that
includes only the framework to install new packages. Unless you're really hard
up for disk space I recommend installing *Anaconda*. 

To install *ipyrad* using *conda* simply type the following into a terminal ::

    $ conda upgrade                 ## updates conda 
    $ conda install ipyrad          ## installs the latest release

If you wish to install a specific version of ipyrad, or to upgrade to the 
latest release from an older version, you could use one of the following commands::

    $ conda install ipyrad=0.7.0    ## install ipyrad v.0.7.0
    $ conda update ipyrad           ## update to the latest


Dependencies
------------
All required dependencies for *ipyrad* should be installed along with 
it when using *conda*. This will include the following Python packages, 
some of which install additional dependencies of their own. 

- Numpy -- Scientific processing  
- Scipy -- Scientific processing  
- Pandas -- Used for manipulating data frames  
- Sphinx -- Used for building documentation  
- IPython2 -- Interactive version of Python 2.7
- ipyparallel -- Parallel, threading, MPI support
- jupyter -- Creating reproducible notebooks  
- Cython -- C bindings for Python 
- H5py -- Database and HDF5 headers
- Dill -- Store pickle objects of complex Classes  
- Toyplot -- [optional].


HPC installation
----------------
One of the benefits of using *conda* for installation is that it 
creates a Python package directory in your home directory called 
``~/anaconda/`` (or ``~/miniconda/``) where new packages are installed.
Make sure you follow the installation instructions_ so that Python
scripts will look in this directory by default. 
Because these Python packages are not stored in a system-wide 
directory you will not need administrator privileges to install 
new packages, nor will you have to load these modules from the 
system before using them. 

TODO: is MPI an exception to this?








