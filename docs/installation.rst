
.. include:: global.rst

.. _installation:


Installation
============

The easiest way to install ipyrad_ and all of its dependencies is with conda_,
a command line program for installing Python packages. If you do not have 
conda_ installed, follow these instructions_ to install either *Anaconda* or 
*Miniconda* for Python2.7. 

If you're working on an :ref:`HPC <HPC_installation>` system you can install 
*conda* into your local directory (e.g., /home/user), for which you will not 
need administrative privileges. This way you can install and access ipyrad_ and 
all its dependencies locally. Specific HPC directions are here 
(HPC_installation_).

The only difference between *Anaconda* and *Miniconda* is that *Anaconda* 
installs a large suite of commonly used Python packages along with the base
installer, whereas *Miniconda* is a bare bones version that includes only 
the framework for installing new packages. Unless you're really hard
up for disk space I recommend installing *Anaconda*. 

To install ipyrad_ using *conda* simply type the following into a terminal

.. code-block:: bash
::
    $ conda update conda                 ## updates conda 
    $ conda install -c ipyrad ipyrad     ## installs the latest release

If you wish to install a specific version of ipyrad, or to upgrade to the 
latest release from an older version, you could use one of the following commands::

.. code-block:: bash
::
    $ conda install -c ipyrad ipyrad=0.1.40    ## install ipyrad v.0.1.40
    $ conda update -c ipyrad ipyrad            ## update to the latest


Dependencies
------------
conda will install all of the following required dependencies during the
installation of ipyrad. 

Python Packages:  

* Numpy -- Scientific processing  
* Scipy -- Scientific processing  
* Pandas -- Used for manipulating data frames  
* Sphinx -- Used for building documentation  
* IPython2 -- Interactive version of Python 2.7  
* ipyparallel -- Parallel, threading, MPI support  
* jupyter -- Creating reproducible notebooks  
* Cython -- C bindings for Python  
* H5py -- Database and HDF5 headers  
* Dill -- Store pickle objects of complex Classes   
* Toyplot -- [optional]  


Executables:  

* vsearch  
* muscle  
* smalt  
* samtools  
* hdf5  
* mpirun, mpiexec  


.. _HPC_installation:

HPC installation
----------------
One of the benefits of using *conda* for installation is that it 
creates a Python package directory where you install it. So if you install into
your home directory it creates ``~/anaconda/`` (or ``~/miniconda/``). It will

.. code-block:: bash

    ## download miniconda  
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

    ## Install miniconda. Follow the directions, by default it will propose installing
    ## to your home directory, which should be fine, e.g., `/home/user/miniconda2`
    ## When asked yes/no whether to append the miniconda directory to $PATH, say yes.  
    bash Miniconda-latest-Linux-x86_64.sh

    ## You could now quit and reconnect, or just run the following command 
    ## which reloads .bashrc so that miniconda will now be in your path. This
    ## way the conda program can be found and run by calling conda.  
    source ~/.bashrc

    ## upgrade conda  
    conda upgrade conda

    ## install ipyrad  
    conda install -c ipyrad ipyrad


Because these Python packages and executables are not stored in a system-wide 
directory you will not need administrator privileges to install them, nor will
you have to load these modules from the system before using them. It is 
important to note that if you let conda append the anaconda directory to your 
$PATH, which you should let it do, that this will become the default location
for these executables. If you want to hide anaconda so that your system defaults
to the system software then you will have to remove the anaconda directory from
your PATH variable in the file `~/.bashrc`, which is in your home directory.

