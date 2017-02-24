
.. include:: global.rst

.. _HPCscript:

Run jupyter-notebook on an HPC cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *ipyrad* API was specifically designed for use inside 
`jupyter-notebooks <http://jupyter.org>`_,  
a tool for reproducible science. 
Notebooks allow you to run interactive code that can be documented with 
embedded Markdown to create a shareable and executable document.
Running *ipyrad* interactively in a notebook is easy to do on 
a laptop or workstation. Simply type ``jupyter-notebook`` into a terminal
and a notebook dashboard will open in your default browser.
For more information see our [introductory tutorial on the ipyrad API] (coming soon). 

Running jupyter-notebooks on a remote HPC cluster is only slightly more 
difficult, but hugely advantageous, because you have access to massively 
more computing power. This tutorial explains how to start a notebook server
on your HPC cluster, and connect to it from your local computer (i.e., your laptop), 
so that you can interact with the notebook in your browser but still have 
the heavy computation occurring remotely on a cluster. 
Instructions below are for the SLURM (sbatch) job submission 
system, we have additional examples available for [TORQUE (qsub) submission 
scripts] (soon) and [others] (soon). 


tldr; Video tutorial
~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="//www.youtube.com/embed/dQw4w9WgXcQ" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
    <br>

Step 1: Submit a batch script to launch a notebook server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    echo -e "\n"
    echo    "  Paste ssh command in a terminal on local host (i.e., laptop)"
    echo    "  ------------------------------------------------------------"
    echo -e "  ssh -N -L $ipnport:$ipnip:$ipnport $USER@$SLURM_SUBMIT_HOST\n"
    echo    "  Open this address in a browser on local host; see token below"
    echo    "  ------------------------------------------------------------"
    echo -e "  localhost:$ipnport                                      \n\n"

    ## start an ipcluster instance and launch jupyter server
    ipcluster start --daemonize
    jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip


Multi-node MPI setup:
For this setup you will have to replace ``module load OpenMPI`` with the 
appropriate module command to load MPI on your system. If you do not know what
this is then look it up for your cluster or ask the system administrator. 

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

    ## get tunneling info
    XDG_RUNTIME_DIR=""
    ipnport=$(shuf -i8000-9999 -n1)
    ipnip=$(hostname -i)

    ## print tunneling instructions to jupyter-log-{jobid}.txt 
    echo -e "\n"
    echo    "  Paste ssh command in a terminal on local host (i.e., laptop)"
    echo    "  -------------------------------------------------------------"
    echo -e "  ssh -N -L $ipnport:$ipnip:$ipnport $USER@$SLURM_SUBMIT_HOST\n"
    echo    "  Open this address in a browser on local host; see token below"
    echo    "  -------------------------------------------------------------"
    echo -e "  localhost:$ipnport                                      \n\n"

    ## initiate MPI & start ipcluster engines using MPI
    module load OpenMPI
    ipcluster start --n=60 --engines=MPI --ip=* --daemonize

    ## start notebook on remote host 
    jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip


If you want to know the details of what this script is doing jump down to 
the section titled 
:ref:`The slurm_jupyter.sbatch script explained<The slurm_jupyter.sbatch script explained>`. 
For now, you can simply submit it to the cluster queue using the sbatch command.

.. code-block:: bash

    user@login-node$ sbatch slurm_jupyter.sbatch



Connecting to the jupyter server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
After submitting your sbatch script to the queue you can check to see if
it has started with the ``squeue`` command. Once it starts information will 
be printed to the log file, ``jupyter-log-{jobid}.txt``. Use the command
``less`` to look at this file and you should see something like below. 

.. code-block:: yaml

     Paste ssh command in a terminal on local host (i.e., laptop)
     ---------------------------------------------------------------- 
     ssh -N -L 8193:xx.yyy.zzz:8193 user@remote.hpc.edu
     ---------------------------------------------------------------
 
     Open this address in a browser on local host; see token below
     ------------------------------------------------------------------
     http://localhost:8193
     ------------------------------------------------------------------

Follow the instructions and paste the `ssh` code block into a terminal on your 
local machine (e.g., laptop). This creates the SSH tunnel from your local 
machine to the port on the cluster where the jupyter server is running. 
As long as the SSH tunnel is open you will be able to interact with the 
jupyter-notebook through your browser. You can close the SSH tunnel at any time 
the notebook will continue running on the cluster. You can also re-connect to it 
later by re-opening the tunnel with the same SSH command.

Security/tokens
~~~~~~~~~~~~~~~~
When you connect to the jupyter-notebook server it will likely ask for a 
password/token. You can find an automatically generated token in your 
jupyter-log file near the bottom. It is the long string printed after the word 
`token`. Copy just that portion and paste it in the token cell. 

Using jupyter
~~~~~~~~~~~~~~
Once connected, you can open any existing notebook, or create a new one. 
The notebooks are physically located on your cluster, meaning all of your data 
and results will be saved there. I usually keep notebooks associated with 
different projects in different directories, where each directory is also a 
github repo, which makes them easy to share. When running ipyrad I usually set 
the "project_dir" be a location in the scratch directory of the cluster, since
it is faster for reading/writing large files. 
You can see an example of this type of setup using the ipyrad API here
(`API empirical notebook <http://nbviewer.jupyter.org/github/dereneaton/pedicularis-WB-GBS/blob/master/nb-WB-Pedicularis.ipynb>`_).


The slurm_jupyter.sbatch script explained
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
What is the sbatch script doing? The ``XDG_RUNTIME_DIR`` command is a little obscure 
and simply fixes a bug where SLURM otherwise sets this variable to something that
is incompatible with jupyter. The ``ipnport`` is a random number between 8000-9999
that selects which port we will use to send data on. The ``ipnip`` is the ip 
address of the login node that we are tunneling through. The ``echo`` commands 
simply print the tunneling information to the log file. In the multi-node 
script there is an additional argument for loading the MPI module, which isn't
necessary on all clusters, some initiate MPI automatically when we run ipcluster, 
but the safest bet is to always load your system MPI module. The final two 
commands are the most important. The first starts an ``ipcluster`` 
instance which will ensure that we can connect to all of the requested CPUs. 
There are many ways to start the parallel client (see the ipyparallel docs), 
but the arguments we used should generally work for most systems.
The final command starts the ``jupyter-notebook`` server, telling it
to forward data to the port that we specified, from the IP address we specified. 


Restarting ipcluster
~~~~~~~~~~~~~~~~~~~~~
Once the connection is established you can later stop and restart ``ipcluster`` 
if you run into a problem with the parallel engines, for example, you might 
have a stalled job on one of the engines. The easiest way to do this is to stop 
the ``ipcluster`` instance by starting a new terminal from the jupyter dashboard, 
by selecting [new]/[terminal] on the right side, and then following
the commands below to restart ``ipcluster``. This does not *always* work for 
reconnecting to multiple nodes over MPI, that will depend on your system, but 
it should always work for restarting a local connection, i.e., no extra 
arguments are passed to ipcluster. 

.. code-block:: bash

    ## stop the running ipcluster instance
    ipcluster stop

    ## start a new ipcluster instance viewing all nodes
    ipcluster start


Connecting multiple notebook at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to run multiple notebooks simultaneously in different tabs and 
have each of them access a different subset of your engines that are available
you can do so using the ``cluster-id`` argument to ipcluster. If you do this 
you will need to tell ipyrad that you are using a non-default ``cluster-id`` 
by setting it in the ipcluster info for your Assembly object (in the JSON 
file for CLI, or in the attribute for the API). 


Terminating the connection
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Just cancel the job on your clusters job queue. You can close the local connections
at any time and reconnect to them later. Remember, the serving is running on the 
cluster. 