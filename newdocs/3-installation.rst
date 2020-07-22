
.. _installation: 

Installation
============

ipyrad can be installed using pip or conda. We strongly recommend the conda version. If you are not familiar with conda then please check out our long-form installation instructions `below <longform_>`__ to start by installing the conda package manager.


Conda install
-------------
ipyrad is available for Python >=2.7 and >=3.5.

.. code:: bash

    conda install ipyrad -c conda-forge -c bioconda


Recommended additional packages
-------------------------------
The ipyrad API provides a powerful interface to using ipyrad for assembling and analyzing data inside of jupyter notebooks, a tool for reproducible science. Many of our downstream analysis tutorials are focused on using the API in jupyter. You can install jupyter easily with the conda command below. In addition, if you wish to distribute jobs over multiple nodes of a HPC cluster then you must install the additional tool mpi4py (more details in Parallelization section).

.. code:: bash

    conda install notebook -c conda-forge
    conda install mpi4py -c conda-forge


We generally recommend setting conda-forge as your default conda channel 
as this reduces the likelihood that you will run into incompatibilities 
later if you install some software dependencies with 
or without it. 

.. code:: bash
    conda config --add channels conda-forge
    conda config --set channel_priority strict


Alternative: install from GitHub
--------------------------------
You can alternatively install ipyrad from its source code on GitHub. This is not recommended unless you're involved in development. 

.. code::bash
    
    # install external requirements first (e.g., using conda)
    conda install ipyrad -c conda-forge -c ipyrad 
    conda install mpi4py notebook -c conda-forge

    # clone the master branch from repo
    git clone -b master https://github.com/dereneaton/ipyrad

    # cd into source and install w/ pip (notice final . in command)
    # this local ipyrad copy will take precedence over the conda copy.
    cd ./ipyrad
    pip install -e .

or alternatively (for version 0.9.56, for example):

.. code::bash

    pip install git+https://github.com/dereneaton/ipyrad.git@0.9.56

Details: dependencies:
----------------------
The following Python packages are installed as dependencies of ipyrad:

    - numpy
    - scipy
    - pandas
    - h5py
    - mpi4py
    - numba
    - ipyparallel
    - pysam
    - cutadapt
    - requests


.. _longform:


Details: Long-form instructions
-------------------------------
We put significant effort into making the installation process for ipyrad as easy as possible, whether you are working on your own desktop computer, or remotely on a large computing cluster. Simply copy and paste a few lines of code below and you will be ready to go.

The easiest way to install ipyrad and all of its dependencies is with conda, a command line program for installing Python packages. Follow
these instructions to first install conda for Python 2 or 3 on your system (the code below is for Python3 since this is now recommended).

Conda comes in two flavors, anaconda and miniconda. The only difference between the two is that anaconda installs a large suite of commonly used Python packages along with the base installer, whereas miniconda installs only a bare bones version that includes just the framework for installing new packages. I recommend miniconda, and that's what we'll use here. 

The code below includes a line that will download the conda installer. **Make sure you follow either the Linux or Mac instructions**, whichever is appropriate for your system. If you are working on an HPC cluster it is almost certainly Linux.

While conda is installing it will ask you to answer **yes** to a few questions. This includes whether it can append the newly created miniconda/ (or anaconda/) directory to your $PATH, say **yes**. What this does is add a line to your **~/.bashrc** (or **~/.bash_profile on Mac**) file so that the software in your conda directory can be automatically found by the systems whenever you login. 


Mac install instructions for *conda*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # The curl command is used to download the installer from the web.
    # Take note that the -O flag is a capital o not a zero.
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

    # Install miniconda into $HOME/miniconda3
    #  * Type 'yes' to agree to the license
    #  * Press Enter to use the default install directory
    #  * Type 'yes' to initialize the conda install
    bash Miniconda3-latest-Linux-x86_64.sh

    # Refresh your terminal session to see conda
    bash

    # test that conda is installed. Will print info about your conda install.
    conda info

Linux install instructions for conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Fetch the miniconda installer with wget
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # Install miniconda into $HOME/miniconda3
    #  * Type 'yes' to agree to the license
    #  * Press Enter to use the default install directory
    #  * Type 'yes' to initialize the conda install
    bash Miniconda3-latest-Linux-x86_64.sh

    # Refresh your terminal session to see conda
    bash

    # test that conda is installed. Will print info about your conda install.
    conda info


.. _HPC_installation:

Details: ipyrad on HPC
^^^^^^^^^^^^^^^^^^^^^^
If you're working on an HPC cluster we still recommend that you follow the 
instructions above to install your own local miniconda directory that you can
use to install local software into. However, you can alternatively ask your 
administrator to install ipyrad into a system-wide conda distribution (and
a specific conda environment) which you and many other users can then use. The 
drawback of this approach is that if you want to upgrade or install additional
software tools you need to ask your administrator and this will likely cause delays.
