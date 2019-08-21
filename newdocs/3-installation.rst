
.. _installation: 

Installation
============

ipyrad can be installed using pip or conda. We strongly recommend the conda version. If you are not familiar with conda then please check out our long-form installation instructions `below <longform_>`__ to start by installing the conda package manager.


Conda install
-------------
ipyrad is available for Python >=2.7 and >=3.5.

.. code:: bash

	conda install ipyrad -c bioconda


Recommended additional packages
-------------------------------
The ipyrad API provides a powerful interface to using ipyrad for assembling and analyzing data inside of jupyter notebooks, a tool for reproducible science. Many of our downstream analysis tutorials are focused on using the API in jupyter. You can install jupyter easily with the conda command below. In addition, if you wish to distribute jobs over multiple nodes of a HPC cluster then you must install the additional tool mpi4py (more details in Parallelization section).

.. code:: bash

	conda install notebook -c conda-forge
    conda install mpi4py -c conda-froge


Alternative: install from GitHub
--------------------------------
You can alternatively install ipyrad from its source code on GitHub. This is not recommended unless you're involved in development. 

.. code::bash
	
	# install external requirements first (e.g., using conda)
	conda install vsearch muscle bedtools bwa samtools pysam -c bioconda
    conda install mpi4py notebook -c conda-forge

	# clone the master branch from repo
	git clone -b master https://github.com/dereneaton/ipyrad

	# cd into source and install w/ pip (notice final . in command)
	cd ./ipyrad
	pip install -e .


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

Conda comes in two flavors, anaconda_ and miniconda_. The only difference between the two is that anaconda_ installs a large suite of commonly used Python packages along with the base installer, whereas miniconda_ installs only a bare bones version that includes just the framework for installing new packages. I recommend miniconda, and that's what we'll use here. 

The code below includes a line that will download the conda_ installer. **Make sure you follow either the Linux or Mac instructions**, whichever is appropriate for your system. If you are working on an HPC cluster it is almost certainly Linux.

While conda is installing it will ask you to answer **yes** to a few questions. This includes whether it can append the newly created miniconda/ (or anaconda/) directory to your $PATH, say **yes**. What this does is add a line to your **~/.bashrc** (or **~/.bash_profile on Mac**) file so that the software in your conda directory can be automatically found by the systems whenever you login. 


Mac install instructions for *conda*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # The curl command is used to download the installer from the web.
    # Take note that the -O flag is a capital o not a zero.
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

    # Install miniconda into $HOME/miniconda3
    bash Miniconda3-latest-MacOSX-x86_64.sh -b 

    # Make it so miniconda is always in your PATH when you open a terminal.
    echo 'PATH=$HOME/miniconda3/bin/:$PATH' >> ~/.bash_profile
    source ~/.bash_profile

    # test that conda is installed. Will print info about your conda install.
    conda info

    # Now run the following command to reload your ~/.bash_profile so that 
    # miniconda will be in your path. This is necessary so that the conda 
    # program can be found from the terminal by simply typing conda. If a 
    # ~/.bash_profile does not exist it might alternatively be named ~/.bashrc.
    source ~/.bash_profile

    # test that conda is installed. This will print info about your conda install.
    conda info


Linux install instructions for conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # The curl command is used to download the installer from the web. Take note
    # that the -O flag is a capital o not a zero.
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # Install miniconda into $HOME/miniconda3
    bash Miniconda3-latest-Linux-x86_64.sh -b 

    # Make it so miniconda is always in your PATH when you open a terminal.
    echo 'PATH=$HOME/miniconda3/bin/:$PATH' >> ~/.bashrc
    source ~/.bashrc

    # test that conda is installed. Will print info about your conda install.
    conda info


.. _HPC_installation:

Details: ipyrad on HPC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you're working on an HPC cluster you should still follow the exact same instructions above to install conda_ into your local home directory (e.g., /home/user). This does not require administrative privileges. In fact, the whole point is to create a local repository for software that you control yourself, separate from the system-wide software. This is 
what we recommend, however, if there is a system-wide version of ipyrad 
installed then you can use that on HPC as well. 
