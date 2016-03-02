
.. include:: global.rst

.. _installation:

Installation
============
We put significant effort into making the installation process for ipyrad as easy 
as possible. Please let us know if any 


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


Mac install
^^^^^^^^^^^^

.. code-block:: bash

    ## The curl command is used to download the installer from the web. 
    ## Take note that the -O flag is a capital o not a zero.
    curl -O https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh

    ## Install miniconda. By default it will propose installing to your 
    ## home directory, which should be fine, e.g., `/home/user/miniconda2`
    ## When asked yes/no to append the miniconda directory to $PATH, say yes.
    bash Miniconda-latest-MacOSX-x86_64.sh

    ## Now either quit and reopen the terminal, or run the following command 
    ## to reload your ~/.bash_profile so that miniconda will be in your path.
    ## This is necessary so that the conda program can be found from the 
    ## terminal by simply typing conda. If a ~/.bash_profile does not exist 
    ## it might alternatively be named ~/.bashrc.
    source ~/.bash_profile

    ## test that conda is installed. This will print info about your conda install.
    conda info


Linux install
^^^^^^^^^^^^^^

    ## The curl command is used to download the installer from the web. Take note
    ## that the -O flag is a capital o not a zero.
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

    ## Install miniconda. Follow the directions, by default it will propose installing
    ## to your home directory, which should be fine, e.g., `/home/user/miniconda2`
    ## When asked yes/no whether to append the miniconda directory to $PATH, say yes.  
    bash Miniconda-latest-Linux-x86_64.sh

    ## You could now quit and reopen the terminal, or just run the following command 
    ## which reloads your ~/.bash_profile so that miniconda will now be in your path.
    ## This is necessary so that the conda program can be found from the terminal by 
    ## simply typing conda. 
    source ~/.bashrc

    ## test that conda is installed. This will print info about your conda install.
    conda info


During installation conda_ will ask if it can append the newly created 
miniconda/ (or anaconda/) directory to your $PATH, say yes. What this does
is add a line to your **~/.bashrc** file so that the anaconda directory becomes 
the default location to search for Python and Python modules, and also so that 
it can find executables in this directory. If you find that setting this path 
interferes with any of your other software you can always comment out the appended
line from **~/.bashrc**. However, the whole point of conda_ is to create 
unique environments in which software packages are protected from conflicting
with each other, so if you run into problems it can likely be fixed using conda_. 


Install ipyrad
^^^^^^^^^^^^^^^
Once conda_ is installed, ipyrad_ can be installed by typing the following 
command into a terminal. This sometimes takes a few minutes to check all of the
dependencies before the installation finishes, so be patient. Make sure you 
do not forget the -c flag. This tells conda that the ipyrad package is located
in a channel called ipyrad.

.. code-block:: bash  

    conda update conda                 ## updates conda 
    conda install -c ipyrad ipyrad     ## installs the latest release

If you wish to install a specific version of ipyrad_, or to upgrade from an 
older version to the most recent release, you can use one of the following 
commands:

.. code-block:: bash  

    conda install -c ipyrad ipyrad=0.1.60    ## install specific version
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

