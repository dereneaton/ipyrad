

.. _installation: 


Installation
============

ipyrad can be installed using pip or conda. We strongly recommend the conda version. If you are not familiar with conda then please check out our long-form installation instructions `below <longform_>`__.


Conda install
-------------
ipyrad is available for Python >=2.7 and >=3.5.

.. code:: bash

	conda install ipyrad -c bioconda


Recommended additional packages
-------------------------------
If you plan to use the ipyrad API in IPython or Jupyter (which we do recommend, especially for downstream analyses) then you should also install IPython and Jupyter, which you can do with the following command.

.. code:: bash

	conda install notebook


Optional: install from GitHub
-----------------------------
You can alternatively install ipyrad from its source code on GitHub. This is not recommended unless you're involved in development. 

... code::bash
	
	# install external requirements first (e.g., using conda)
	conda install vsearch muscle bedtools bwa samtools mpi4py -c bioconda -c conda-forge

	# clone the master branch from repo
	git clone -b master https://github.com/dereneaton/ipyrad

	# cd into source and install w/ pip (notice final . in command)
	cd ./ipyrad
	pip install .


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


Details: install Conda
----------------------
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

    # Install miniconda. By default it will install into your home directory.
    # e.g., `/home/user/miniconda3`
    bash Miniconda3-latest-MacOSX-x86_64.sh -b 

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

    # Install miniconda. Follow the directions, by default it will propose installing
    # to your home directory, which should be fine, e.g., `/home/user/miniconda3`
    # When asked yes/no whether to append the miniconda directory to $PATH, say yes.
    bash Miniconda3-latest-Linux-x86_64.sh -b 

    # You could now quit and reopen the terminal, or just run the following command
    # which reloads your ~/.bashrc so that miniconda will now be in your path.
    # This is necessary so that the conda program can be found from the terminal by
    # simply typing conda.
    source ~/.bashrc

    # test that conda is installed. This will print info about your conda install.
    conda info


.. _HPC_installation:

Details: using ipyrad on a HPC cluster
--------------------------------------
If you're working on an HPC cluster you should still follow the exact same instructions above to install conda_ into your local home directory (e.g., /home/user). This does not require administrative privileges. In fact, the whole point is to create a local repository for software that you control yourself, separate from the system-wide software. 

This is useful because it then allows you to install and access ipyrad_ and all its dependencies (other Python modules and executables), and to update them yourself. Lot's of useful software is available on conda, which you can find and install by googling conda and the software name. 
