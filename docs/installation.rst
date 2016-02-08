
.. include:: global.rst

.. _installation:

Installation
============

Conda install
-------------

The easiest way to install ipyrad_ and all of its dependencies is with conda_,
a command line program for installing Python packages. If you do not have 
conda_ installed, you can find detailed installation instructions_here_, or simply
follow the directions outlined below. You will need to install either `anaconda` or 
`miniconda` for Python2.7. The only difference between the two is that 
`anaconda` installs a large suite of commonly used Python packages along with the 
base installer, whereas `miniconda` is a bare bones version that includes only 
the framework for installing new packages. 


First, we need to download the conda installer, I'll use `miniconda` for my example. 
There are separate installers for Linux and Mac, so choose only the one that
is correct for your system. If you are working on an HPC cluster it is most \
likely a Linux machine. 

.. code-block:: bash

    ## EITHER download Miniconda for Linux 
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh    

    ## OR download Miniconda for Mac
    wget https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh

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
`miniconda/` (or `anaconda/`) directory to your `$PATH`, say yes. What this does
is add a line to your `~/.bashrc` file so that the anaconda directory becomes 
the default location to search for Python and Python modules, and also so that 
it can find executables in this directory. If you find that setting this path 
interferes with any of your other software you can always comment out the 
line from `~/.bashrc`. But, if conflicts arise my advice would be rather than 
working around conda_, to instead learn more about it (instructions_), 
since it is a powerful tool for creating multiple environments with different 
software packages that do not conflict. 

If you're working on an :ref:`HPC <HPC_installation>` system you can install 
conda_ into a locally owned directory (e.g., /home/user) without need for 
administrative privileges. This is useful because it then allows you to install 
and access ipyrad_ and all its dependencies (other Python modules and 
executables) locally, without needing to load them from the system-wide 
software. 

Once conda_ is installed, ipyrad_ can be installed easily by entering the 
following into a terminal:

.. code-block:: bash  

    conda update conda                 ## updates conda 
    conda install -c ipyrad ipyrad     ## installs the latest release

If you wish to install a specific version of ipyrad_, or to upgrade to the 
latest release from an older version, you could use one of the following 
commands:

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
* Toyplot -- Plotting 


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
One of the benefits of using conda_ for installation is that it 
creates a Python package directory where you install it. So if you install into
your home directory it creates `~/anaconda/` (or `~/miniconda/`). 


Because these Python packages and executables are not stored in a system-wide 
directory you will not need administrator privileges to install them, nor will
you have to load these modules from the system before using them. You should find
that you can simply type the name of the software and load it, even on jobs 
submitted using `qsub` or `sbatch`. 

As an example, you can play around with ipyrad within an IPython terminal:  

.. code-block:: bash  

    ## open an ipython shell in the terminal
    ipython

In the IPython session load ipyrad. See API_ usage for more details:  

.. code-block:: python  

    ## import ipyrad under its common shortname
    import ipyrad as ip

    ## create a test Assembly object
    data = ip.Assembly("test")

    ## print the default parameters
    data.get_params()


It is important to note that if you let conda append the conda directory to your 
$PATH, which you should let it do, that this will become the default location
for these executables. If you want to hide anaconda so that your system defaults
to the system software then you will have to remove the anaconda directory from
your PATH variable in the file `~/.bashrc`, which is in your home directory.


