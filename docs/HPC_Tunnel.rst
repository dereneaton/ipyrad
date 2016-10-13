
.. _HPCscript:

Tunnel Jupyter notebook to an HPC cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *ipyrad* API was designed with the intention of being used inside Jupyter 
notebooks, which provide a convenient way of executing interactive code, and
documenting it with embedded with Markdown, to create highly reproducible workflows.
Running ipyrad interactively in this way is simple to do on a normal 
laptop or workstation. Simply call `jupyter-notebook` from a terminal
to open a notebook in your browser. You can setup your local cluster setup using
`ipcluster`. 

Running Jupyter notebooks is more difficult, however, when you want to execute
code remotely on a HPC cluster. This tutorial explains how to setup an SSH Tunnel 
so that you can execute code through a Jupyter notebook on your local computer 
(i.e., your laptop) while having the actual heavy computation be executed remotely
on compute nodes of your cluster. 


#### Step 1: Set up SSH tunneling  
The instructions below are for connecting to multiple compute nodes interactively, 
however, you could similarly execute the code on compute nodes through a 
submission script. The main thing we are trying to do is to tunnel information 
using SSH from the compute node to our local machine, using the login node as 
an intermediate. 

.. .. parsed-literal::

.. code-block:: bash

    ## use SSH to connect to the login node with your credentials
    user@local$ ssh user@hpc_login_node.edu  

    ## use the 'screen' command to allow connecting/disconnecting from login
    user@login$ screen

    ## from the (screened) login node connect to a compute node interactively
    user@login$ qsub -I -l nodes=4:ppn=8 

    ## from the compute node launch ipcluster to connect all accessable nodes
    user@compute$ ipcluster start --n=32 --engines=MPI --ip=* --daemonize --profile=ipyrad

    ## from the compute node start a Jupyter-notebook sending data to port 8888
    user@compute$ jupyter-notebook --no-browser --port=8888  

    ## take note of your compute node hostname
    user@compute$ <- it will be in place of 'compute' here.


#### Step 2: Go back to login node (use ctrl-d to disconnect from 'screen')  

.. code-block:: bash

    ## tell login node to tunnel data over port 8888 from compute node
    user@login$ ssh -N -f -L localhost:8888:localhost:8888 user@compute

    ## tell your local machine where to look for the Jupyter notebook
    user@local$ ssh -N -f -L localhost:8887:localhost:8888 user@hpc_login_node.edu



#### Step 3: Open a browser to http://localhost:8887  

You should see a Jupyter notebook view of the files in your home directory 
on the HPC cluster. You can open an existing notebook, or start a new one. The notebooks
are located on your cluster, meaning all of your data and results will be saved there. 
I like to store my notebook in a git repo in my home directory, and to store all of 
my data that I am working on at the time in a scratch directory. This way, all of the
code you execute in your notebook can be saved to your git repo, and you basically have
a perfect supplementary materials document for your study. When finished with my assembly, 
I might upload the final outfiles to somewhere more permanent than the scratch dir, until they can be archived for publication. They may be too big for github, in which case Zenodo is another good choice. 


#### Other SSH notes:  

If you try to form the SSH tunnel but the port is already being used, first check whether it is one of your previous jobs using the port, in case you may want to close it. You can see which PID is associated with the port by using the command:  

.. code-block:: bash

    ## which PID is using port 8888?
    user@login$ lsof -ti:8888

    ## let's say it returned pid=311. To kill it do the following:
    user@login$ kill 311

