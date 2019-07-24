

.. _installation: 


Installation
============

ipyrad can be installed using pip or conda. We strongly recommend the conda 
version. If you are not familiar with conda then please check out our long-form
installation instructions `below <longform_>`__.


Conda install
-------------
ipyrad is available for Python >=2.7 and >=3.5.

.. code:: bash

	conda install ipyrad -c bioconda -c conda-forge


Recommended additional packages
-------------------------------
If you plan to use the ipyrad API in IPython or Jupyter (which we do recommend,
especially for downstream analyses) then you should also install IPython and 
Jupyter, which you can do with the following command.

.. code:: bash

	conda install notebook


Alternative: install from GitHub
--------------------------------
You can alternatively install ipyrad from its source code on GitHub. This 
is not recommended unless you're involved in development. 

... code::bash
	
	# install external requirements first (e.g., using conda)
	conda install vsearch muscle bedtools bwa samtools mpi4py -c bioconda -c conda-forge

	# clone the master branch from repo
	git clone -b master https://github.com/dereneaton/ipyrad

	# cd into source and install w/ pip (notice final . in command)
	cd ./ipyrad
	pip install .


Dependencies:
-------------
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


Long-form instructions
----------------------
We put significant effort into making the installation process for ipyrad
as easy as possible, whether you are working on your own desktop computer, or
remotely on a large computing cluster. Simply copy and paste a few lines of
code below and you will be ready to go.


Conda install
-------------
The easiest way to install ipyrad and all of its dependencies is with conda,
a command line program for installing Python packages. Follow
these instructions to first install conda for Python 2 or 3 on 
your system (the code below is for Python3 since this is now recommended).

Conda comes in two flavors, anaconda_ and miniconda_. The only difference
between the two is that anaconda_ installs a large suite of commonly used
Python packages along with the base installer, whereas miniconda_ installs
only a bare bones version that includes just the framework for installing
new packages. I recommend miniconda, and that's what we'll use here. 

The code below includes a line that will download the conda_ installer. 
**Make sure you follow either the Linux or Mac instructions**, whichever is 
appropriate for your system. If you are working on an HPC cluster it is 
almost certainly Linux.

While conda is installing it will ask you to answer **yes** to a few questions. 
This includes whether it can append the newly created miniconda/ (or anaconda/) 
directory to your $PATH, say **yes**. What this does is add a line to your 
**~/.bashrc** (or **~/.bash_profile on Mac**) file so that the software in your
conda directory can be automatically found by the systems whenever you login. 


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



ipyrad install
--------------
Once conda_ is installed, ipyrad can be installed by typing the following
command into a terminal. This sometimes takes a few minutes to check all of the
dependencies before the installation finishes, so be patient. Make sure you
do not forget the ``-c ipyrad`` flag. This tells conda that the ipyrad package
is located in a channel called ipyrad.

.. code-block:: bash

    conda update conda                                # updates conda
    conda install ipyrad -b bioconda -c conda-forge   # installs the latest release

If you wish to install a specific version of ipyrad_, or to upgrade from an
older version to the most recent release, you can use one of the following
commands:

.. code-block:: bash

    conda install -c bioconda ipyrad=0.9.10      ## install specific version
    conda update -c bioconda ipyrad              ## update to the latest


.. _HPC_installation:

How does this work on a HPC cluster?
------------------------------------
If you're working on an HPC cluster you should still follow the exact same 
instructions above to install conda_ into your local home directory 
(e.g., /home/user). This does not require administrative privileges. In fact, 
the whole point is to create a local repository for software that you control
yourself, separate from the system-wide software. 

This is useful because it then allows you to install and access ipyrad_ and all 
its dependencies (other Python modules and executables), and to update them 
yourself. Lot's of useful software is available on conda, which you can find 
and install by googling conda and the software name. 


How do I ignore or remove conda?
---------------------------------
Conda is super easy to ignore or remove if you ever find that it is not working
for you. Conda itself, as well as all of the software that it installs is 
located in the miniconda/ directory, and so you *can* remove all of it by 
removing that directory. I would advise, however, that a much simpler way to 
switch on/off conda software would be to simply comment out the line in your 
``~/.bashrc`` file that appends miniconda/ to your PATH. Then run 
``source ~/.bashrc`` and your system will completely ignore the conda software. 
Likewise, you can uncomment the line, re-source the file, and your conda 
software will be back. Conda is hugely popular, but it is also quite new, 
and actively under development, which has caused some issues with compatibility
when major updates have arisen over the last 1-2 years. If you have a quite old
conda distribution (pre v.4) that is giving you troubles when you try to update
software I would recommend removing it and reinstalling. You can then reinstall
all of your conda software quite easily.
