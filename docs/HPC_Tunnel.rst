
.. _HPCscript:

SSH Tunnel Jupyter notebook to an HPC cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *ipyrad* API was designed with the intention of being used inside Jupyter 
notebooks, which provide a convenient way of executing interactive code, and
documenting it with embedded Markdown, to create highly reproducible workflows.
Running ipyrad interactively in this way is simple to do on a normal 
laptop or workstation. Simply call `jupyter-notebook` from a terminal
to open a notebook in your browser, and setup a local cluster using
`ipcluster`. This is explained in the standard *ipyrad* API tutorial. 

Running Jupyter notebooks on a remote cluster is only slightly more difficult, 
but hugely advantageous, because of course you have access to massively more 
computing power. This tutorial explains how to setup a SSH Tunnel so that you 
can execute code in a Jupyter notebook on your local computer (i.e., your laptop) 
but have the actual heavy computation be executed remotely on compute nodes of 
your cluster. We have instructions below for both TORQUE (i.e., qsub based systems)
as well as SLURM (i.e., sbatch based systems). 


Step 1: Launch jupyter-notebook and ipcluster on a compute node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For SLURM I save the following script as `slurm_launch_jupyter_cluster.sh`, and 
submit it with the sbatch command. This will need to be edited slightly to 
conform to your cluster, by changing the name of the partition (queue), 
changing the number of nodes to whatever works best for you, and setting 
the walltime limit. Fortunately, for this cluster the walltime
limit is 30 days, so I can submit the script and keep a connection to 32 cores
for that amount of time. 

.. code-block:: bash

    #!/bin/bash
    # set the number of nodes and processes per node
    #SBATCH --partition general
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 8
    #SBATCH --exclusive
    #SBATCH --time 30:00:00
    #SBATCH --job-name jupyter-ipcluster
    #SBATCH --output jupyterlog.txt

    ## launch ipcluster engines across available cpus
    ipcluster start --n=32 --engines=MPI --ip=* --daemonize

    ## necessary fix for a slurm bug that otherwise crashes jupyter
    XDG_RUNTIME_DIR=""

    ## open a tunnel between compute and login nodes on port 8181
    NOTEBOOKPORT=8181
    ssh -N -f -R $NOTEBOOKPORT:localhost:$NOTEBOOKPORT $SLURM_SUBMIT_HOST

    ## launch notebook
    jupyter-notebook --no-browser --port=$NOTEBOOKPORT


Alternatively, if you are using TORQUE, then submit the following script using 
the qsub command instead, which I save under the name 
`torque_launch_jupyter_cluster.sh`. 

.. code-block:: bash

    #!/bin/bash
    # set the number of nodes and processes per node    
    #PBS -l nodes=4:ppn=8
    #PBS -walltime=30:00:00
    #PBS -j oe
    #PBS -N jupyter-cluster
    #PBS -q general

    ## launch ipcluster engines across all available cpus
    ipcluster start --n=32 --engines=MPI --ip=* --daemonize

    ## open a tunnel between compute and login nodes on port 8181
    NOTEBOOKPORT=8181
    ssh -N -f -R $NOTEBOOKPORT:localhost:$NOTEBOOKPORT $PBS_O_HOST

    ## launch notebook
    jupyter-notebook --no-browser --port=$NOTEBOOKPORT


Step 2: Open ssh connection to your cluster from local
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is similar to the normal way you would login to your HPC cluster, except that
you tell it to forward all information it receives on port 8181 to your 
local port 8181. Also change the login credentials to your name and host. 

.. code-block:: bash
    
    ssh -N -L localhost:8181:localhost:8181 user@hpc_login_node.edu


Step 3: Tunnel from local computer to notebook 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now we simply open a browser to http://localhost:8181

You should see the Jupyter notebook view of your filesystem on the HPC cluster. 
You can open an existing notebook, or start a new one. The notebooks are located
on your cluster, meaning all of your data and results will be saved there. I 
like to store my notebooks inside directories that are each separate git repos
in my home directory, and to store all of my big data in a scratch directory. 
You can see an example like that :ref:`here<http://nbviewer.jupyter.org/github/dereneaton/pedicularis-WB-GBS/blob/master/nb-WB-Pedicularis.ipynb>`. This way, the notebook records all of 
the code you execute in your notebook which can be saved to your git repo, 
but all of the giant data is still saved in scratch. 


Connecting multiple notebook at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to run jobs in multiple notebooks simultaneously then you should open
a second port, rather than run both notebooks through the same port, otherwise
they will be sharing the same ipcluster instance, and thus fight over the 
available engines. Instead start a second ipcluster by submitting a second 
submission script to your cluster to launch a different ipcluster instance and 
jupyter-notebook. Make sure you designate a *different* port number. 
You can use any port number between 8000-9000. 


Terminating the connection
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To disconnect the jupyter notebook and ipcluster running remotely simply kill/cancel
the running job on your cluster. To terminate the SSH connection from your local 
machine that is viewing an open port, you can simply close/cancel the ssh connection
running in a terminal. If you have it running in the background and can't find the
running ssh job, you can run the following to find whatever is looking into your
open port (e.g., 8181). Then simply call 'kill' to terminate that process id. 

.. code-block:: bash

    ## which PID is using port 8181?
    user@login$ lsof -ti:8181

    ## let's say it returned pid=31189. To kill it do the following:
    user@login$ kill 31189

