
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
submit it with the `sbatch` command. The code may need to be edited slightly to 
conform to your cluster, which may require setting the name of the partition (queue), 
or changing the number of nodes and walltime to within the limits that are
allowed. 

.. code-block:: bash

    #!/bin/bash
    # set the number of nodes and processes per node
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 30-00:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name jptr
    #SBATCH --output jupyterlog.txt

    ## pick a port number between 8000-9999
    NOTEBOOKPORT=8888
    
    ## load MPI [the exact command may vary on your cluster!]
    module load MPI/OpenMPI

    ## launch ipcluster engines across available cpus
    ipcluster start --n=80 --engines=MPI --ip=* --daemonize

    ## necessary fix for a slurm bug that otherwise crashes jupyter
    XDG_RUNTIME_DIR=""

    ## open a tunnel between compute and login nodes on port 8181
    ssh -N -X -f -R $NOTEBOOKPORT:localhost:$NOTEBOOKPORT $SLURM_SUBMIT_HOST

    ## launch notebook
    jupyter-notebook --no-browser --port=$NOTEBOOKPORT


On the login node of your cluster submit the submission script. 

.. code-block:: bash

    user@login-node$ sbatch slurm_launch_jupyter_cluster.sh


Alternatively, if you are using TORQUE, then submit the following script using 
the qsub command instead, which I save under the name `torque_launch_jupyter_cluster.sh`. 

.. code-block:: bash

    #!/bin/bash
    # set the number of nodes and processes per node    
    #PBS -l nodes=4:ppn=8
    #PBS -walltime=30:00:00
    #PBS -j oe
    #PBS -o jupyterlog.txt
    #PBS -N jptr
    #PBS -q general

    ## choose a port number between 8000-9999
    NOTEBOOKPORT=8888

    ## launch ipcluster engines across all available cpus
    ipcluster start --n=32 --engines=MPI --ip=* --daemonize

    ## open a tunnel between compute and login nodes on port 8181
    ssh -N -X -f -R $NOTEBOOKPORT:localhost:$NOTEBOOKPORT $PBS_O_HOST

    ## launch notebook
    jupyter-notebook --no-browser --port=$NOTEBOOKPORT


.. code-block:: bash

    ## submit the submission script
    user@login-node$ qsub torque_launch_jupyter_cluster.sh


Step 2: Open ssh connection to your cluster from local
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is similar to the normal way you would login to your HPC cluster, except that
you tell it to forward all information it receives on port xxxx to your 
local port xxxx. Also change the login credentials to your name and host. 
If you forget which port you entered in your submission script you can 
check the log file on your cluster, which we named jupyterlog.txt.

.. code-block:: bash
    
    user@local$ ssh -N -L 8888:localhost:8888 user@hpc_login_node.edu


Step 3: Tunnel from local computer to notebook 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now simply open a browser to http://localhost:8181

You should see the Jupyter notebook view of your filesystem on the HPC cluster. 
You can open an existing notebook, or start a new one. The notebooks are 
physically located on your cluster, meaning all of your data and results will be 
saved there. I usually sync my working directories in which notebooks reside 
using github, which makes them easy to share. I usually set the "project_dir"
parameter in ipyrad to be in a scratch directory. 
You can see an example of this type of setup here:
:ref:`here<http://nbviewer.jupyter.org/github/dereneaton/pedicularis-WB-GBS/blob/master/nb-WB-Pedicularis.ipynb>`. 
This way, the notebook records all of the code you execute in your notebook 
which can be saved to your git repo, but all of the giant data is 
saved in scratch. 


Connecting multiple notebook at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to run multiple notebooks simultaneously you can do so from 
a single port, by simply opening new notebooks from the Homepage. 
If you started an ipcluster instance in your submission script, then 
all notebooks can access this instance. If you would rather divide up between
multiple notebooks you can do so by opening a terminal from the Jupyter
homepage, and running the command `ipcluster stop` to stop the instance 
that is running. Then you can start two separate ipcluster instances in 
the terminal by assigning each a different number of clusters (-n=X) and 
assigning them different IDs (cluster-id=X). In your notebooks you then 
have to tell your Assemblies which ipcluster instance to connect to by 
assigning a 'cluster_id' dictionary variable. For example, 
`Assembly._ipcluster["cluster_id"] = "ip-8888". 


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
    user@local$ lsof -ti:8181

    ## let's say it returned pid=31189. To kill it do the following:
    user@local$ kill 31189

    ## when I SSH locally I see the error `channel 2: open failed: connect failed: Connection refused`:
    Check to make sure you are entering the correct port number. If you did and you still see this message,
    try running the SSH script from a new terminal, and try connecting to localhost:8181 in a new browser window.
    If that still doesn't work, try a different port number, the one you chose may already be in use. 



