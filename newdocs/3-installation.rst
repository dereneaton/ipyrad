
.. _installation: 

Installation
============

ipyrad can be installed using pip or conda. We strongly recommend the conda version. If you are not familiar with conda then please check out our long-form installation instructions `below <longform_>`__ to start by installing the conda package manager.


Conda install
-------------
ipyrad is available for Python >=3.5 (3.10 is recommended).

.. code:: bash

    conda install ipyrad -c conda-forge -c bioconda


Recommended additional conda settings 
-------------------------------------

We generally recommend installing ipyrad into a clean conda environment,
which (after conda is installed) can be achieved like this:

.. code:: bash

    conda create -n ipyrad
    conda activate ipyrad
    conda install ipyrad -c conda-forge -c bioconda

We also recommend setting conda-forge as your default conda channel 
as this reduces the likelihood that you will run into incompatibilities 
later if you install some software dependencies with 
or without it. 

.. code:: bash

    conda config --add channels conda-forge
    conda config --set channel_priority strict

Fixing conda install stuck on 'Solving Environment'
---------------------------------------------------
Recently (2023-ish) conda has had some problems resolving dependencies and can get stuck solving
them 'forever'. A workaround, for the moment, is to install and use the `libmamba-solver`.

.. code:: bash

    conda update -n base conda
    conda install -n base conda-libmamba-solver
    conda config --set solver libmamba

Installing ipyrad on Mac M1/M2
------------------------------
As of mid-2023 the bioconda channel does not have support for ARM architecture (Mac M1/M2).
A workaround is to install the Mac Intel version of conda (x86_64) using the Terminal App
running with Rosetta.

Check if you have Rosetta installed, and if not then install it. If you have Rosetta
installed, or if the `pgrep oahd` command below returns an integer (whole number value)
then you can skip this step:

.. code:: bash

    # Check if Rosetta is installed
    pgrep oahd

    # If this does not return any information then run the Rosetta installer
    softwareupdate --install-rosetta --agree-to-license

With Rosetta installed you can now run the Terminal App in x86_64 compatibility mode.
Open Applications->Utilities->Terminal (Or spotlight search Terminal):

.. code:: bash

    # Activate the terminal in x86_64 mode    
    env /usr/bin/arch -x86_64 /bin/zsh ---login

    # Fetch the miniconda Mac x86 installer and run it
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

    # Accept the license, the default install location, and say 'yes' to initialization
    sh Miniconda3-latest-MacOSX-x86-64.sh

    # Activate the new conda install by running `zsh`. The prompt will change to show
    # the (base) conda environment is now activated
    zsh

    # Create and activate a new conda env for ipyrad
    conda create -n ipyrad
    conda activate ipyrad

    # Set the default solver to libmamba
    conda config --set solver libmamba

    # Install ipyrad
    conda install -c conda-forge -c bioconda ipyrad

Add the following two lines to the end of .zshrc to make it easier to switch the Terminal
App to using x86 vs arm:

.. code:: bash

    alias x86="$env /usr/bin/arch -x86_64 /bin/zsh ---login"
    alias arm="$env /usr/bin/arch -arm64 /bin/zsh ---login"

After editing the .zshrc run `source ~/.zshrc`. Now you can type `x86` to switch to Rosetta
x86 mode in the Terminal. Check the architecture by typing `arch`, it will say `i386` if it
is in x86_64 compatibility mode.

One gotcha can be if you have previously installed conda on your system with osx-arm64 your
config might want to install arm64 packages. You can disable this by removing this option from
the conda config

.. code:: bash

    # Use conda config --show to inspect 'subdirs', if osx-arm64 shows up then remove it
    conda config --remove subdirs osx-arm64

Alternative: install from GitHub
--------------------------------
You can alternatively install ipyrad from its source code on GitHub. This is not recommended unless you're involved in development. 

.. code:: bash
    
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

.. code:: bash

    pip install git+https://github.com/dereneaton/ipyrad.git@0.9.56

Details: dependencies:
----------------------
The following Python packages are installed as dependencies of ipyrad:

    - numpy
    - scipy
    - pandas
    - h5py
    - mpi4py
    - notebook
    - numba
    - ipyparallel
    - pysam
    - cutadapt
    - requests
    - muscle
    - samtools
    - bedtools
    - bwa
    - vsearch

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
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

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
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

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
