
.. include:: global.rst

.. _installation:

Installation
============

Conda install
-------------

The easiest way to install ipyrad_ and all of its dependencies is with conda_,
a command line program for installing Python packages. If you do not have 
conda_ installed, you can find detailed installation 
:ref:`here <...>`, or simply
follow the directions outlined below. You will need to install either anaconda_ 
or miniconda_ for Python2.7. The only difference between the two is that 
anaconda_ installs a large suite of commonly used Python packages along with the 
base installer, whereas miniconda_ installs only a bare bones version that 
includes just the framework for installing new packages. 

First, we need to download the conda_ installer, I'll use miniconda_ for my example. 
There are separate installers for Linux and Mac, so choose only the one that
is correct for your system. If you are working on an :ref:`HPC<HPC_installation>`
cluster it is most likely a Linux machine. Skip this section if you already 
have conda installed.

.. code-block:: bash

    ## EITHER download Miniconda for Linux 
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh    

    ## OR download Miniconda for Mac
    ## wget is not available by default on Mac so you have to use curl. The
    ## -O flag is a capital o not a zero
    curl -O https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh

    ## Install miniconda. Follow the directions, by default it will propose installing
    ## to your home directory, which should be fine, e.g., `/home/user/miniconda2`
    ## When asked yes/no whether to append the miniconda directory to $PATH, say yes.  
    ## The example here is for the Linux installer.
    bash Miniconda-latest-Linux-x86_64.sh

    ## You could now quit and reopen the terminal, or just run the following command 
    ## which reloads ~/.bashrc so that miniconda will now be in your path. This
    ## way the conda program can be found from the terminal by simply typing conda.  
    source ~/.bashrc

    ## test that conda is installed by printing info about conda
    conda info


During installation conda_ will ask if it can append the newly created 
miniconda/ (or anaconda/) directory to your $PATH, say yes. What this does
is add a line to your **~/.bashrc** file so that the anaconda directory becomes 
the default location to search for Python and Python modules, and also so that 
it can find executables in this directory. If you find that setting this path 
interferes with any of your other software you can always comment out the 
line from **~/.bashrc**. However, the whole point of conda_ is to create 
unique environments in which software packages are protected from conflicting
with each other, so if you run into problems it can likely be fixed using conda_. 

Once conda_ is installed, ipyrad_ can be installed by typing the following 
command into a terminal:

.. code-block:: bash  

    conda update conda                 ## updates conda 
    conda install -c ipyrad ipyrad     ## installs the latest release

If you wish to install a specific version of ipyrad_, or to upgrade from an 
older version to the most recent release, you could use one of the following 
commands:

.. code-block:: bash  

    conda install -c ipyrad ipyrad=0.1.40    ## install specific version
    conda update -c ipyrad ipyrad            ## update to the latest


.. _HPC_installation: 

HPC installation
-----------------

If you're working on an :ref:`HPC <HPC_installation>` system you can install 
conda_ into a locally owned directory (e.g., /home/user) without need for 
administrative privileges. This is useful because it then allows you to install 
and access ipyrad_ and all its dependencies (other Python modules and 
executables) locally, without needing to load them from the system-wide 
software. 


Included Dependencies
------------
The conda_ installation will install the following required dependencies:

**Python Packages**:  

* Numpy -- Scientific processing  
* Scipy -- Scientific processing  
* Pandas -- Used for manipulating data frames  
* Sphinx -- Used for building documentation  
* IPython2 -- Interactive version of Python 2.7  
* ipyparallel_ -- Parallel, threading, MPI support  
* jupyter -- Creating reproducible notebooks  
* Cython -- C bindings for Python  
* H5py -- Database and HDF5 headers  
* Dill -- Store pickle objects of complex Classes   
* Toyplot -- Plotting 


**Executables**:  

* vsearch -- used for de novo clustering
* muscle_ -- used for sequence alignment
* smalt -- used for reference mapping
* samtools -- used for reference mapping
* bedtools -- used for reference mapping
* hdf5 -- used for large array storage/access
* mpirun, mpiexec -- used for parallelization

