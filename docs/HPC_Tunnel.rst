
.. include:: global.rst

.. _HPCscript:


Running jupyter-notebooks locally and remotely (on HPC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *ipyrad* API was specifically designed for use inside 
`jupyter-notebooks <http://jupyter.org>`_, a tool for reproducible science.
This section of the documentation is about how to start and run jupyter
notebooks, which you can then use to run your ipyrad analyses using
the ipyrad API. For instructions on how to use the ipyrad API 
(after you have a notebook started) go here: (:ref:`ipyrad API <API>`__). 
An example of a complete notebook showing assembly and analysis of 
a RAD data set with the ipyrad API can be found here:
(`Pedicularis API <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-1-pedicularis.ipynb>`__).


(:ref:`full API example <http://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/tests/cookbook-empirical-API-1-pedicularis.ipynb>`).  

Jupyter notebooks allow you to run interactive code that can be 
documented with embedded Markdown (words and fancy text) 
to create a shareable and executable document. 
Running *ipyrad* interactively in a notebook 
is easy to do on a laptop or workstation, and slightly more difficult
to run an HPC cluster, but after reading this tutorial you will 
hopefully find it easy to do. If this is your
first time using jupyter it will be easiest to start by trying 
it on your laptop first before trying to use jupyter on a cluster. 
In the case of running on a cluster our example below include an 
example job submission script for SLURM, but other job 
submission systems should be similar. 


The following tools are used in this section:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
+ ipyrad (used for RAD-seq assembly)  
+ jupyter-notebook (an environment in which we run Python code)  
+ ipcluster (used to parallelize code within a notebook)  
+ ssh (used to connect to a notebook running on HPC)  


Starting a jupyter-notebook **locally**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To start a jupyter-notebook on your local computer (e.g., laptop)
execute the command below in a terminal. This will start a local 
notebook server and open a window in your default web-browser. 
Leave the server running the terminal. You will not need to 
touch that again until you want to stop the notebook server.
You can now interact with the notebook server through your 
web-browser. You should see a page showing the files and folders
in your directory where you started the notebook. In the upper
right you will see a tab where you can select <new> and then 
<Python 2> to start a new Python notebook.

.. code-block:: bash

    jupyter-notebook


Starting a jupyter-notebook **remotely**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Because jupyter works by sending and receiving information 
(i.e., it's a server) it is easy to run a jupyter notebook through 
your browser even if the notebook server is running on a remote computer 
that is far away, for example on a computing cluster. Start by 
assigning a password to your notebook server which will give it 
added security. 

.. code-block:: bash

    ## Run this on the remote mahcine (i.e., the cluster)
    ## It will ask you to enter a password which will be 
    ## encrypted and stored for use when connecting.
    jupyter-notebook password


.. code-block:: bash

    ## Run this on the remote machine (i.e., the cluster)
    jupyter-notebook --no-browser --ip=$(hostname -i) --port=9999  


Once the notebook starts it will print some information including 
the IP address of the machine your are connected to (this will something
like 10.115.0.25), and the port number that it is using (this will
probably be 9999 if that is what you entered above, however, if 9999
is already in use then it will select a different port number, so check
the output). You will need these to pieces of information, the IP-address
and the port number, for the next command. Replace the values that are 
between brackets with appropriate values. 


.. code-block:: bash

    ## Run this on your local machine (i.e., your laptop)
    ssh -N -L <port>:<ip-address>:<port> <user>@<login>


.. code-block:: bash

    ## This would be an example with real values entered:
    ssh -N -L 9999:10.115.0.25:9999 deren@hpc.columbia.edu  



Starting jupyter through a batch script:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tldr; short video tutorial.

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/hjBJw1fY5Uo" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
    <br>


Step 1: Submit a batch script to start jupyter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Copy and paste the code block below into a text editor and save the script as 
``slurm_jupyter.sbatch``. The #SBATCH section of the script may need to be edited 
slightly to conform to your cluster. The stdout (output) of the job will be 
printed to a log file named ``jupyter-log-%J.txt``, where %J will be replaced 
by the job ID number. We'll need to look at the log file once the job starts
to get information about how to connect to the jupyter server that we've started.


Single Node setup:
This example would connect to one node with 20 cores available. 

.. code-block:: bash

    #!/bin/bash
    #SBATCH --partition general
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 4:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name tunnel
    #SBATCH --output jupyter-log-%J.txt

    ## get tunneling info
    XDG_RUNTIME_DIR=""
    ipnport=$(shuf -i8000-9999 -n1)
    ipnip=$(hostname -i)

    ## print tunneling instructions to jupyter-log-{jobid}.txt 
    echo -e "
        Copy/Paste this in your local terminal to ssh tunnel with remote 
        -----------------------------------------------------------------
        ssh -N -L $ipnport:$ipnip:$ipnport user@host                     
        -----------------------------------------------------------------

        Then open a browser on your local machine to the following address
        ------------------------------------------------------------------
        localhost:$ipnport  (prefix w/ https:// if using password)
        ------------------------------------------------------------------
        "

    ## start an ipcluster instance and launch jupyter server
    jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip


Now submit the sbatch script to the cluster to reserve the node and 
start the jupyter notebook server running on it. The notebook server 
will continue running until it hits the walltime limit, or you stop it.

.. code-block:: bash

    ## submit the job script to your cluster job scheduler
    sbatch slurm_jupyter.sbatch



Step 2: Connecting to the jupyter server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After submitting your sbatch script to the queue you can check to see if
it has started with the ``squeue -u {username}`` command. 
Once it starts information will be printed to the log file which 
we named ``jupyter-log-{jobid}.txt``. Use the command ``less`` 
to look at this file and you should see something like below. 

.. code-block:: yaml

     Copy/paste this in your local terminal to ssh tunnel with remote
     ---------------------------------------------------------------- 
     ssh -N -L 8193:xx.yyy.zzz:8193 user@host
     ---------------------------------------------------------------
 
     Then open a browser on your local machine to the following address
     ------------------------------------------------------------------
     localhost:8193  (prefix w/ https:// if using password)
     ------------------------------------------------------------------

Follow the instructions and paste the `ssh` code block into a terminal on your 
local machine (e.g., laptop). This creates the SSH tunnel from your local 
machine to the port on the cluster where the jupyter server is running. 
As long as the SSH tunnel is open you will be able to interact with the 
jupyter-notebook through your browser. You can close the SSH tunnel at any time 
and the notebook will continue running on the cluster. You can 
re-connect to it later by re-opening the tunnel with the same SSH command.


.. code-block:: bash

    ## This would be an example with real values entered:
    ssh -N -L 8193:10.115.0.25:8193 deren@hpc.columbia.edu  


Security/tokens
~~~~~~~~~~~~~~~~
If you did not create a password earlier, then when you connect to 
the jupyter-notebook server it will ask you for a password/token. 
You can find an automatically generated token in your jupyter-log 
file near the bottom. It is the long string printed after the word 
`token`. Copy just that portion and paste it in the token cell. I 
find it easier to use password. See the jupyter documentation for how
to setup further security. 


Using jupyter
~~~~~~~~~~~~~~
Once connected you can open an existing notebook or create a new one. 
The notebooks are physically located on your cluster, meaning all of your data 
and results will be saved there. I usually keep notebooks associated with 
different projects in different directories, where each directory is also a 
github repo, which makes them easy to share. When running ipyrad I usually set 
the "project_dir" be a location in the scratch directory of the cluster, since
it is faster for reading/writing large files. 


Using ipcluster on a multi-node MPI setup:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the example above we started a notebook on a node with 20 cores available. 
Once connected, the first I would do is typically to start an ipcluster instance
running in a terminal so that I can use it to parallelize computations
(see our `ipyparallel tutorial`__). If you want to connect to multiple nodes, 
however, then it is better to start the ipcluster instance separately in its 
own separate job submission script. Here is an example. Importantly, we will 
tell ipcluster to use a specific `--profile` name, in this case named `MPI60`, 
to indicate that we're connecting to 60 cores using MPI. When we connect
to the client later we will need to provide the profile name. I name this file
``slurm_ipcluster_MPI.sbatch``. 

For this setup we also add a command to load the MPI module. You will probably
need to modify ``module load OpenMPI`` to whatever the appropriate module 
name is for MPI on your system. If you do not know what this is then look
it up or ask the system administrator. 


.. code-block:: bash

    #!/bin/bash
    #SBATCH --partition general
    #SBATCH --nodes 3
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 30-00:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name MPI60
    #SBATCH --output ipcluster-log-%J.txt

    ## set the profile name here
    profile="MPI60"

    ## Start an ipcluster instance. This server will run until killed.
    module load OpenMPI
    sleep 10
    ipcluster start --n=60 --engines=MPI --ip='*' --profile=$profile



Connecting to the ipcluster instance in Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now when you are in the jupyter notebook you can connect to this ipcluster
instance -- which is running as a completely separate job on your cluster -- 
with the following simple Python code. The object ``ipyclient`` can then
be used to distribute your computation on the remote cluster. When you
run ipyrad pass the ipyclient object to tell it this is the cluster you want
computation to occur on. The results of your computation will still be 
printed in your jupyter notebook.


.. code-block:: python
    
    import ipyrad as ip
    import ipyparallel as ipp

    ## connect to the client
    ipyclient = ipp.Client(profile="MPI60")

    ## print how many engines are connected
    print(len(ipyclient), 'cores')

    ## or, use ipyrad to print cluster info
    ip.cluster_info(ipyclient)


.. code-block:: yaml

    60 cores
    host compute node: [20 cores] on c14n02.farnam.hpc.yale.edu
    host compute node: [20 cores] on c14n03.farnam.hpc.yale.edu
    host compute node: [20 cores] on c14n04.farnam.hpc.yale.edu


When running the ipyrad API you would distribute work by passing the
ipyclient object in the ipyclient argument. See the ipyrad API for more
information. 

.. code-block:: python

    ## run step 3 of ipyrad assembly across 60 cores of the cluster    
    data.run(steps='3', ipyclient=ipyclient)



The slurm_jupyter.sbatch script explained
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
So what is the sbatch script above doing? 
The ``XDG_RUNTIME_DIR`` command is a little obscure, it simply fixes a 
bug where SLURM otherwise sets this variable to something that
is incompatible with jupyter. The ``ipnport`` is a random number between 8000-9999
that selects which port we will use to send data on. The ``ipnip`` is the ip 
address of the login node that we are tunneling through. The ``echo`` commands 
simply print the tunneling information to the log file. 

In the multi-node ipcluster script we use a the ``module load``
command to load the system-wide MPI software. Then we call ipcluster
with arguments to find cores across all available nodes using MPI, and
we provide a name (profile) for this cluster so it will be easy 
to connect to.


Restarting ipcluster
~~~~~~~~~~~~~~~~~~~~~
Once the connection is established you can later stop and restart ``ipcluster`` 
if you run into a problem with the parallel engines, for example, you might 
have a stalled job on one of the engines. The easiest way to do this is to stop 
the ``ipcluster`` instance by starting a new terminal from the jupyter dashboard, 
by selecting [new]/[terminal] on the right side, and then following
the commands below to restart ``ipcluster``. If you are using a multi-node
setup then you will need to resubmit the ipcluster job through a script in 
order to connect to multiple computers again. 

.. code-block:: bash

    ## stop the running ipcluster instance
    ipcluster stop

    ## start a new ipcluster instance viewing all nodes
    ipcluster start


Terminating the connection
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To stop a running jupyter notebook just cancel the job on your cluster's queue, 
or if working locally, just press control-c in the terminal window. If you 
disconnect from a remote notebook and later reconnect you can continue 
using the notebook without needed to restart it by going to the menu and 
select kernel reconnect. If progress bars were printing output while you 
were disconnected it may not show up, but the job will have kept running.
The loss of progress bars is a shortcoming that will likely be 
fixed in the near future. 

