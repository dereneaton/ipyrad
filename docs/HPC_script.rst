



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
The first step to doing work on a HPC cluster is to connect from your local 
terminal to the HPC cluster using SSH.

.. parsed-literal::

    ## connect to the login node from your computer
    >>> ssh user@address_of_HPC_cluster

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
cores you plan to connect to with the `-c` flag, and also use the `--MPI` flag
to allow it to efficiently access cores across multiple nodes. 


Submitting job scripts
----------------------
Below is an example *qsub* script, which is used to submit jobs on a TORQUE system, 
as well as an example *sbatch* script that is used to submit jobs on a SLURM system. 
As you can see they're pretty similar. 


.. parsed-literal::

    #!/bin/sh

    #PBS -N ipyrad
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=2:ppn=8

    ## change into your home dir, or a specific place from there
    cd $PBS_O_WORKDIR/myanalysis/

    ## call some ipyrad commands
    ipyrad -p params-demux.txt -s 1 -c 16 --MPI
    ipyrad -p params-test.txt -b newbranch
    ipyrad -p params-newbranch.txt -s 2345 -c 16 --MPI
    ipyrad -p params-newbranch.txt -r 


.. parsed-literal::
    ## submit the qsub script
    >>> qsub qsub_script.sh


And here is an example *sbatch* script:

.. parsed-literal::

    #!/bin/sh

    #PBS -N ipyrad-test
    #PBS -j oe
    #PBS -m ae
    #http://ipyrad.readthedocs.io/HPC_Tunnel.htmlPBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=8:ppn=8

    ## change into whichever directory you which to run code from
    cd $PBS_O_WORKDIR

    ## call ipyrad on your params file
    ipyrad -p params-test.txt -s 1 -c 64 --MPI


.. parsed-literal::

    ## submit the qsub script
    >>> qsub qsub_script.sh



Running interactive jobs
------------------------
From there you can then try logging in interactively to a 
computing node by using the command ``qsub`` with the ``-I``
argument. Sometimes you have to provide additional
arguments such as the name of the queue you are connecting to.
This information should be available from your institution.
Gaining access to the node may be instant, or it may take hours
depending on the size of your cluster and how many users are 
active.

.. parsed-literal::

    ## connect to an node interactively
    >>> qsub -I 


You can use all of the typical qsub arguments to connect
to multiple CPUs, or specific queues. You can think of a 
computing cluster as many smaller computers linked together. 
To access CPUs that are not only from a single computer, 
but from many of those computers you have to provide 
a few extra arguments. A single computer, or ``node``, 
will typically have 4, 8, 16, or 32 cores. The command below
will access 64 cores spread across 8 8-core nodes. 


.. parsed-literal::

    ## ask for 64 cores across 8 nodes from queue 'fas_general' 
    ## and request 24 hours of wall time.
    >>> qsub -I -l nodes=8:ppn=8 -l walltime=24:00:00 -q "fas_general"
    
.. parsed-literal::

    ## On SLURM systems the command is somewhat different.
    >>> srun -A lfg_lab -p serial -t 120:00:00 -N 1 -n 5 --pty --mem-per-cpu=6000 /bin/bash



Optional: Running ipcluster by hand
------------------------------------

On some hpc compute nodes ipcluster does not spin up fast enough
and ipyrad times out. To work around this it is possible to start
ipcluster by hand, wait for it to fully fire up, then connect to
it with the CLI. ipyrad has an argument `--ipcluster`, which when 
enabled will tell it to skip trying to create an ipcluster 
instance and to instead connect to the existing ipcluster 
instance with `profile=ipryad`. It's up to you to remember to run 
`ipcluster stop --profile=ipyrad` when you're done. The 
`--daemonize` flag tells ipcluster to run in the background.

.. parsed-literal::
    ## Get an interactive shell on a compute node
    $ qsub -I -l nodes=8:ppn=8 -l walltime=24:00:00
    
    ## start an ipcluster instance with --profile=ipyrad
    $ ipcluster start --n 48 --profile=ipyrad --daemonize

    ## run ipyrad with --ipcluster flag so it knows to look for 
    ## that specific ipcluster instance
    $ ipyrad -p params-test.txt -s 2 --ipcluster

In the event that you want to run ipcluster by hand _and_ take
advantage of MPI you need to add a couple more arguments.

.. parsed-literal::
    ## start ipcluster with MPI enabled
    $ ipcluster start --n 48 --profile=ipyrad --daemonize --ip=* --engines=MPI

    ## run ipyrad and talk to the MPI enabled ipcluster you just started
    $ ipyrad -p params-test.txt -s 2 --ipcluster
