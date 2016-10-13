
.. _HPCscript:

Tunnel Jupyter notebook to an HPC cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ipyrad API was designed with the intention of being used within a Jupyter 
notebook, which provides a convenient way of executing code embedded between
Markdown cells to create highly reproducible workflows. Running ipyrad 
interactively in this way is simple to do on a normal laptop or workstation, 
but more difficult to run remotely on a HPC cluster. This tutorial explains 
how to setup an SSH Tunnel so that you can execute code through a Jupyter notebook
on your local computer (i.e., your laptop) while having the actual heavy 
computation executing remotely on compute nodes of a HPC cluster. 


#### Step 1: Set up SSH tunneling
The instructions below are for connecting to multiple compute nodes interactively, 
however, you could similarly execute the interactive part as a submission script.
The main thing we are trying to do is to Tunnel information using SSH from the 
compute node to our local machine, using the login node as an intermediate. 

.. parsed-literal::

    ## use SSH to connect to the login node with your credentials
    user@local$ ssh user@hpc_login_node.edu  

    ## use the 'screen' command to allow connecting/disconnecting from login
    user@login$ screen

    ## from the (screened) login node connect to a compute node interactively
    user@login$ qsub -I -l nodes=4:ppn=8 

    ## from the compute node start a Jupyter-notebook sending data to port 8888
    user@compute$ jupyter-notebook --no-browser --port=8888 &

    ## from the compute node launch ipcluster engines on all connected compute nodes
    ## is not running interactively, do not use the --daemonize argument
    user@compute$ ipcluster start --n=32 --engines=MPI --ip=* --daemonize


#### Step 2: Go back to login node (use ctrl-d to disconnect from 'screen')

    ## tell login node to tunnel data over port 8888 from compute node
    user@login$ ssh -N -f -L localhost:8888:localhost:8888 user@compute

	## tell your local machine where to look for the Jupyter notebook
	user@local$ ssh -N -f -L localhost:8887:localhost:8888 user@hpc_login_node.edu



#### Step 3: Open a browser to http://localhost:8887
You should see a Jupyter notebook view of the files in your home directory 
on the HPC cluster. 
