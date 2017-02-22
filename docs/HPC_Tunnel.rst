
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
the heavy computation occurring remotely on the cluster. 
Instructions below are for the SLURM (sbatch) job submission 
system, we have [examples using TORQUE (qsub) submission scripts available as well]. 


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
``slurm_jupyter.sbatch``. The #SBATCH section of the script will need to be edited 
slightly to conform to your cluster, which may require setting the name of the 
partition (queue), or changing the number of nodes and walltime limits. In this 
script we are requesting 60 cores (3 nodes, 20 cores per node). The stdout (output)
of the job will be printed to the log file named ``jupyter-log-%J.txt``, where 
%J will be replaced by the job ID number. 

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


What is this script doing? The ``XDG_RUNTIME_DIR`` command is a little obscure 
and simply fixes a bug where SLURM otherwise sets this variable to something that
is incompatible with jupyter. The ``ipnport`` is a random number between 8000-9999
that will be used to send data. The ``ipnip`` is the ip address of the login 
node that we are connected to. The ``echo`` commands simply print this 
information to the log file so we will know how to connect to our notebook 
server once it has started. 

The final two commands are the most important. The first starts an ``ipcluster`` 
instance which will ensure that we can connect to all of the requested CPUs. 
There are many ways to start this parallel client (see the ipyparallel docs), 
but the arguments we used above should generally work for most systems.
The final command starts the ``jupyter-notebook`` server, telling it
to forward data to the port that we specified. 

Now you can submit the script to the queue using the ``sbatch`` command:

.. code-block:: bash

    user@login-node$ sbatch slurm_jupyter.sbatch

You can check the queue to see if the job has started using the ``squeue`` command. 
Once it has started information will be printed to the log file, which is
named ``jupyter-log-{jobid}.txt``. Use the command ``less`` to look at this file and
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
tunnel at any time and your code will continue to run on the Jupyter-notebook. 
You could re-connect later to the same notebook by re-opening the tunnel with 
the same SSH command.

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


Restarting ipcluster
~~~~~~~~~~~~~~~~~~~~~
It is necessary to start the ``ipcluster`` instance in our sbatch script in order
to initialize a connection to all of the avialable CPUs. However, once the connection
has been established we can later stop and restart ``ipcluster`` however we wish
and it will continue to find the same engines. Sometimes if an error arises and 
you want to kill the ipcluster engines the easiest way is to stop the ``ipcluster``
instance. You can do this by starting a new terminal from the jupyter dashboard, 
and selecting [new]/[terminal] on the right side. In the terminal run the following
commands to restart ``ipcluster``. You can close the tab if you wish but the 
terminal will remain running on the remote system. You can use ``ctrl-c`` to
stop the ipcluster instance after you restart it once in this way. 

.. code-block:: bash

    ## stop the running ipcluster instance
    ipcluster stop

    ## start a new ipcluster instance viewing all nodes
    ipcluster start --n=60 --engines=MPI --ip=*


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
To close/disconnect the jupyter-notebook and ipcluster instance simply kill/cancel
the job running on your cluster. To terminate the SSH connection from your local 
machine that is viewing an open port, you can simply close the ssh connection
running in a terminal. 