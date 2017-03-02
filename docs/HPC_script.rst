



.. _HPCscript:

Running ipyrad on a cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^

High performance computing (HPC) clusters are accessible to most 
users with a University affiliation. The HPC software set up often varies 
among different institutions, so we provide here a general workflow 
for running ipyrad on TORQUE or SLURM submission software.

We designed ipyrad from the beginning to work seamlessly on HPC setups. This 
begins with the fact that the software is installed locally (in your miniconda directory)
and thus does not require you to locate or load any other software on the HPC system. 
Secondly, ipyrad can distribute work across as many cores and nodes as you have 
access to with just a few simple commands. 


Basics of HPC
----------------
First you should connect to the cluster from your local machine using SSH with 
your login credentials.

.. parsed-literal::

    ## connect to the login node from your computer
    user@local$ ssh user@address_of_HPC_cluster

You will then be connected to the ``login``, which is *not* where you should run 
ipyrad or any other software. It is just a landing point from which you can access
your files on the system, and from where you should submit jobs to processing nodes.
Jobs can be submitted using either the *sbatch* or *qsub* commands, depending on your
system. Here we show two ways to run ipyrad, either interactively, or by submitting job
scripts. If possible, I prefer interactive jobs because you can watch the 
progress more easily, but either way is fine. 

.. note::

    After you master using HPC with simple submission scripts, check out our
    `tutorials on using ipyrad in Jupyter Notebooks with SSH tunneling <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__, it's kind of advanced stuff, but a really effective way to run interactive 
    code once you have it working. 


General notes on running ipyrad on HPC
---------------------------------------
When running ipyrad on a cluster you should make sure to tell it explicitly how many 
cores you plan to connect to with the ``-c`` flag, and also use the ``--MPI`` flag
to allow it to efficiently access cores across multiple nodes. 


Submitting job scripts
----------------------
Below is an example *qsub* script, which is used to submit jobs on a TORQUE system, 
as well as an example *sbatch* script that is used to submit jobs on a SLURM system. 
As you can see they're pretty similar. 


.. parsed-literal::

    #!/bin/sh
    # set the number of nodes and processes per node
    #PBS -N ipyrad
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=1:ppn=20

    ## change into your home dir, or a specific place from there
    cd $PBS_O_WORKDIR/myanalysis/

    ## call some ipyrad commands
    ipyrad -p params-demux.txt -s 1 -c 20  
    ipyrad -p params-test.txt -b newbranch  
    ipyrad -p params-newbranch.txt -s 2345 -c 20  
    ipyrad -p params-newbranch.txt -r   


.. parsed-literal::
    ## submit the qsub script
    user@login$ qsub qsub_script.sh


And here is an example *sbatch* script:

.. parsed-literal::

    #!/bin/bash
    # set the number of nodes and processes per node
    #SBATCH --partition general
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 16
    #SBATCH --exclusive
    #SBATCH --time 30-00:00:00
    #SBATCH --mem-per-cpu 2000
    #SBATCH --job-name ipyrad
    #SBATCH --output ipyrad_output.txt

    ## change into whichever directory you which to run code from
    cd $HOME/myanalysis

    ## call ipyrad on your params file
    module load OpenMPI
    ipyrad -p params-test.txt -s 1234567 -c 64 --MPI

.. parsed-literal::
    ## submit the qsub script
    user@login$ sbatch slurm_script.sh


Running interactive jobs
------------------------
For testing purposes it is best to login interactively to a compute node. 
This can be done on TORQUE with the *-I* argument. Sometimes you have to 
provide additional arguments such as the name of the queue you are connecting to.
This information should be available from your institution. Gaining access to 
the node may be instant, or it may take hours depending on the size of your 
cluster and how many users are active.

.. parsed-literal::
    ## connect to a compute node interactively
    user@login$ qsub -I 

You could similarly provide all of the typical qsub arguments with this command: 

.. parsed-literal::
    ## ask for 64 cores across 8 nodes from queue 'fas_general' 
    ## and request 24 hours of wall time.
    user@login$ qsub -I -l nodes=8:ppn=8 -l walltime=24:00:00 -q "fas_general"
    
.. parsed-literal::

    ## On SLURM systems the command is somewhat different.
    user@login$ srun -A lfg_lab -p serial -t 120:00:00 -N 1 -n 5 --pty --mem-per-cpu=6000 /bin/bash



Optional: Controlling ipcluster by hand
------------------------------------
ipyrad uses a program called *ipcluster* to control parallelization, most of which 
occurs behind the scenes for the user. However, it is possible to gain more 
fine-tuned control of the connection to parallel CPUs by starting the ipcluster
instance yourself, and using the `--ipcluster` argument to ipyrad to tell it to 
find your running ipcluster instance. 

This has proved useful on a few HPC clusters where compute nodes spin up 
very slowly, and ipyrad would quit after a few minutes if it didn't find the 
connected CPUs it was looking for. To work around this the user can spin up
ipcluster with the arguments listed below, then add in a sleep command to tell
the system to wait a minute, and then the ipyrad command. 

.. parsed-literal::

    ## Login in to an interactive node
    user@login$ qsub -I -l nodes=1:ppn=20 -l walltime=24:00:00
    
    ## Now that you are on the compute node, start an ipcluster instance 
    user@compute$ ipcluster start --n 20 --daemonize

    ## Wait for ipcluster. Sleeping for 60 seconds should be sufficient.
    sleep 60

    ## Then run ipyrad like normal but with --ipcluster so it knows to look for 
    ## your specific ipcluster instance.
    user@compute$ ipyrad -p params-test.txt -s 2 --ipcluster



