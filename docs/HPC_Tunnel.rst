
.. _HPCscript:

SSH Tunnel Jupyter notebook to an HPC cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *ipyrad* API was specifically designed for use inside jupyter-notebooks, 
a tool for reproducible science. Notebooks allow you to run interactive code 
that can documented with embedded Markdown to create shareable and executable
documentation. Running *ipyrad* interactively in a notebook is easy to do on 
a laptop or workstation. Just type `jupyter-notebook` into a terminal
and a notebook dashboard will automatically launch in your default browser.
See our [introductory tutorial on the ipyrad API]. 

Running jupyter-notebooks on a remote HPC cluster is only slightly more 
difficult, but hugely advantageous, because you have access to massively 
more computing power. This tutorial explains how to setup a SSH Tunnel so that
you can execute code in a Notebook on your local computer (i.e., your laptop)
but have the heavy computation occurring remotely on compute nodes of your 
cluster. Instructions below are for the SLURM (sbatch) job submission 
system, we also have an [example using TORQUE (qsub) here]. 


Step 1: Submit a batch script to launch a jupyter-notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Copy and paste the codeblock below into a text editor and save the script as 
`slurm_jupyter.sh`. The code may need to be edited slightly to conform to your 
cluster, which may require setting the name of the partition (queue), or 
changing the number of nodes and walltime to within the limits that are 
allowed. 

.. code-block:: bash

    #!/bin/bash
    #SBATCH --partition general
    #SBATCH --nodes 3
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 30-00:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name jptr60
    #SBATCH --output jupyter-log-%J.txt

    ## a required bugfix for slurm/jupyter
    XDG_RUNTIME_DIR=""

    ## selects a random port number 
    ipnport=$(shuf -i8000-9999 -n1)

    ## gets ip of compute node host
    ipnip=$(hostname -i)

    ## prints tunneling instructions to ipyrad-log file
    echo -e "\n\n   Copy/Paste this in your local terminal to ssh tunnel with remote "
    echo        "   ------------------------------------------------------------------"
    echo        "   ssh -N -L $ipnport:$ipnip:$ipnport $USER@$SLURM_SUBMIT_HOST "
    echo        "   ------------------------------------------------------------------"
    echo -e "\n\n   Then open a browser on your local machine to the following address"
    echo        "   ------------------------------------------------------------------"
    echo        "   localhost:$ipnport"
    echo -e     "   ------------------------------------------------------------------\n\n"
    sleep 1

    ## start an ipcluster instance here to init MPI
    ipcluster start --n=60 --engines=MPI --ip=* --daemonize

    ## start notebook on remote host 
    jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip


When executed, this script does the following: 

(i) parses the #SBATCH info to request resources, here asking for 3 nodes 
with 20 cores each, to not allow others to connect to our nodes (exclusive),
to let us stay connected for 30 days, to have at least 4GB RAM per CPU, to
name the job 'jptr60' (corresponding to the fact that we requested 60 CPUs)
and to create an output file called 'jupyter-log-%J.txt', where %J will be 
replaced by the job ID number. 

(ii) we get some info for tunneling jupyter. The XDG_RUNTIME_DIR argument
fixes a bug where SLURM otherwise sets this variable to something that is 
incompatible with jupyter. The `ipnport` argument selects a random port number 
between 8000-9999. The `ipnip` is the ip address of the login node we are 
connected to. The `echo` command prints information to the log file about how 
to connect to our notebook once it has started. 

(iii) We start an ipcluster instance. This is the parallel client that ipyrad 
will use to send and receive data across multiple nodes or machines. There are 
many ways to set this up, see the ipyparallel docs for tips. The arguments I 
pass to it here should generally work for any system, though. 

(iiii) Launch a jupyter-notebook server from this ip address exporting data to 
the port number that we specified. 


Now that we have our sbatch script prepared, simply submit the job to the queue
using the `sbatch` command. 

.. code-block:: bash

    user@login-node$ sbatch slurm_jupyter.sh

You can check the queue to see if the job has started using the `squeue` command. 
Once it has started information will be printed to the log file, which will be 
named `jupyter-log-{jobid}.txt`. Use the command `less` to look at this file and
you should see something like below. 


.. code-block:: yaml

     Copy/Paste this in your local terminal to ssh tunnel with remote 
     ---------------------------------------------------------------- 
     ssh -N -L 8193:xx.yyy.zzz:8193 user@remote.hpc.edu
     ---------------------------------------------------------------
 
 
     Then open a browser on your local machine to the following address
     ------------------------------------------------------------------
     localhost:8193
     ------------------------------------------------------------------

Follow the instructions from the logfile and paste the `ssh` code block into 
a terminal on your local machine (e.g., laptop). This creates the SSH tunnel
from your local machine to the remote compute node on your cluster. As long
as the SSH tunnel is open you should be able to view the Jupyter-notebook in 
your browser by going to the localhost address listed. You can close the SSH
tunnel at any time and your code will continue to run on the Jupyter-notebook, 
and you can re-connect later by re-opening the tunnel with the same SSH command.


Security
~~~~~~~~
When you connect to your jupyter-notebook server in your browser you will likely
be asked for a password/token. You can find the token in your jupyter-log file
near the bottom. It is the long string printed after the word `token`. 


Using jupyter
~~~~~~~~~~~~~~
Once connected, you can open an existing notebook, or start a new one. The notebooks are 
physically located on your cluster, meaning all of your data and results will be 
saved there. I usually sync my working directories in which notebooks reside 
using github, which makes them easy to share. I usually set the "project_dir"
parameter in ipyrad to be in a scratch directory. 
You can see an example of this type of setup here:
:ref:`here<http://nbviewer.jupyter.org/github/dereneaton/pedicularis-WB-GBS/blob/master/nb-WB-Pedicularis.ipynb>`. 
This way, the notebook records all of the code you execute in your notebook 
which can be saved to your git repo, while all of the giant data is 
saved in scratch. 


Restarting ipcluster
~~~~~~~~~~~~~~~~~~~~~



Connecting multiple notebook at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to run multiple notebooks simultaneously you can do so from 
a single port, by simply opening new notebooks from the Dashboard. 
If you started an ipcluster instance in your submission script, then 
all notebooks can access this instance. If you would rather divide the cores 
so only some of them are available to each notebook the easiest way to do this
is to start a new separate ipcluster instance for each. To do this, connect to 
a terminal from your Jupyter dashboard by clicking [New] and then [Terminal]. 
Then stop your existing ipcluster instance by running `ipcluster stop`. 
Now you can start a new distince `ipcluster` instances in 
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


Troubleshooting
~~~~~~~~~~~~~~~
+ I see the error `channel X: open failed: connect failed: Connection refused`  

Check to make sure you are entering the correct port number. If you did and you still see this message,
I find the problem is most easily fixed by closing the terminal on your local
machine and opening a new one. For some reason this seems to reset something 
that allows the connection to work again. 


