
.. include:: global.rst

.. _installation:

Installation
============

Conda install
-------------

The easiest way to install ipyrad_ and all of its dependencies is with conda_,
a command line program for installing Python packages. If you do not have 
conda_ installed, follow these instructions_ to install either *Anaconda* or 
*Miniconda* for Python2.7. 

The only difference between *Anaconda* and *Miniconda* is that *Anaconda* 
installs a large suite of commonly used Python packages along with the base
installer, whereas *Miniconda* is a bare bones version that includes only 
the framework for installing new packages. Unless you're hard up for disk space
I recommend installing *Anaconda*. 

If you're working on an :ref:`HPC <HPC_installation>` system you can install 
conda_ into a locally owned directory (e.g., /home/user) without need for 
administrative privileges. This is useful because it then allows you to install 
and access ipyrad_ and all its dependencies (other Python modules and 
executables) locally, without needing to load them from the system-wide 
software. More detailed HPC directions are here: HPC_installation_.

Once conda_ is installed, ipyrad_ can be installed easily by entering the 
following into a terminal:

.. code-block:: bash  

    conda update conda                 ## updates conda 
    conda install -c ipyrad ipyrad     ## installs the latest release

If you wish to install a specific version of ipyrad_, or to upgrade to the 
latest release from an older version, you could use one of the following 
commands::

.. code-block:: bash  

    conda install -c ipyrad ipyrad=0.1.40    ## install specific version
    conda update -c ipyrad ipyrad            ## update to the latest


Dependencies
------------
The conda_ installation will install the following required dependencies:

**Python Packages**:  

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


**Executables**:  

* vsearch -- used for de novo clustering
* muscle -- used for sequence alignment
* smalt -- used for reference mapping
* samtools -- used for reference mapping
* bedtools -- used for reference mapping
* hdf5 -- used for large array storage/access
* mpirun, mpiexec -- used for parallelization


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

