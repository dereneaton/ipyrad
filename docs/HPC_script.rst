

.. _HPCscript:

Running ipyrad on a cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^

High performance computing (HPC) clusters are accessible to most 
users with a University affiliation and can massively increase the speed
at which your analyses are run. HPC systems often vary from among different 
institutions and we have aimed to make ipyrad's parallel workflow flexible 
enough to work across most (hopefully any) system. Our examples below 
are for TORQUE or SLURM submission software, but again, the general idea
should hold for other systems as well.

We designed ipyrad from the beginning to work seamlessly on HPC systems. 
This begins with the fact that the software is installed locally 
(in your miniconda directory) and thus does not require you to locate or 
load many other types of software on the HPC system. 
Secondly, ipyrad can distribute work across as many cores and nodes as you 
have access to with just a single extra command. 


Basics of HPC
----------------
To begin connect to your cluster's login node using SSH with your login credentials.

.. parsed-literal::

    user@local$ ssh user@address_of_HPC_cluster

Once you are connected to the login node you can begin to submit jobs using 
submission scripts. You should *not* run your jobs directly on the login node.
It is simply a landing point from which to to access files on the system, and 
to submit jobs to processing nodes. You can now submit your submission script
with the appropriate command for you system. For a SLURM system the command 
is `sbatch`, like below. 


.. parsed-literal::
    user@login$ sbatch slurm_script.sh


.. note::

    After you master using HPC with simple submission scripts, check out our
    `tutorials on using ipyrad in Jupyter Notebooks with SSH tunneling <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__, it's kind of advanced stuff, but a really effective way to run interactive code once you have it working. 


Part I: Using multiple cores on a single node
---------------------------------------------
Running ipyrad on a single node is very simple since it is essentially the 
same process as running it on your own computer. As a best practice 
for organizing ipyrad files on your system, I typically run ipyrad from the 
location where my params file is located, and I use the params file itself 
to set the explicit path to my data files. This involves setting the path to
your raw data (usually somewhere in your scratch directory) and I also set 
the `project_dir` to be in my scratch directory. The project directory will 
be created once the analysis starts if it does not already exist.
This kind of setup allows you to keep your params files in one convenient 
place (I use a local directory called ipyrad-analyses/) and to have each 
params file write to its own separate output location. 

Below we show two example submission scripts to run a single-node job 
on either a TORQUE (qsub) or SLURM (sbatch) system:


Example *qsub* script:

.. parsed-literal::

    #!/bin/sh
    #PBS -N ipyrad
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=1:ppn=20

    ## change into the directory with your params file
    cd $PBS_O_WORKDIR/ipyrad-analyses/

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
    #SBATCH --partition general
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 7-00:00:00
    #SBATCH --mem-per-cpu 2000
    #SBATCH --job-name ipyrad
    #SBATCH --output ipyrad_output.txt

    ## change into the directory where your params file resides
    cd $HOME/ipyrad-analyses/

    ## call ipyrad on your params file
    ipyrad -p params-test.txt -s 1234567 -c 20

.. parsed-literal::
    ## submit the qsub script
    user@login$ sbatch slurm_script.sh



Part II: Running multi-node jobs
--------------------------------
Accessing cores from multiple nodes (essentially multiple computers) 
requires that you use the `--MPI` flag to turn on the *message passing interface*
and that you also tell ipyrad explicitly how many cores you are planning to 
connect to with the `-c` flag. For MPI, this is the one case where you do 
need to load software that is external to ipyrad using a `module load` command. 
Use the appropriate command on your system (usually something like 
`module load MPI`, but check to be sure what it is exactly on your system). 
Examples below:

Example *qsub* script:

.. parsed-literal::

    #!/bin/sh
    #PBS -N ipyrad
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=4:ppn=20

    ## load MPI
    module load MPI

    ## change into your home dir, or a specific place from there
    cd $PBS_O_WORKDIR/ipyrad-analyses/

    ## call some ipyrad commands 
    ipyrad -p params-demux.txt -s 1 -c 80 --MPI
    ipyrad -p params-test.txt -b newbranch  
    ipyrad -p params-newbranch.txt -s 2345 -c 80 --MPI
    ipyrad -p params-newbranch.txt -r 


.. parsed-literal::
    ## submit the qsub script
    user@login$ qsub qsub_script.sh


And here is an example *sbatch* script:

.. parsed-literal::

    #!/bin/bash
    #SBATCH --partition general
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 7-00:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name ipyrad
    #SBATCH --output ipyrad_output.txt

    ## change into the directory where your params file resides
    cd $HOME/ipyrad-analyses/

    ## call ipyrad on your params file
    ipyrad -p params-test.txt -s 1234567 -c 80 --MPI

.. parsed-literal::
    ## submit the qsub script
    user@login$ sbatch slurm_script.sh



Running interactive jobs
------------------------
For testing purposes it is best to login interactively to a compute node. 
This can be done on TORQUE with the -I argument. Sometimes you have to 
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
    user@login$ qsub -I -l nodes=8:ppn=8 -l walltime=24:00:00 -q "queue_name"
    
.. parsed-literal::
    ## On SLURM systems the command is somewhat ugly.
    user@login$ srun -p general -t 120:00:00 -N 1 -n 5 --pty --mem-per-cpu=4000 /bin/bash



Optional: Controlling ipcluster by hand
------------------------------------
ipyrad uses a program called *ipcluster* (from the ipyparallel Python module)
to control parallelization, most of which occurs behind the scenes for the user.
However, it is possible to gain more fine-tuned control of the connection to 
parallel CPUs by starting the ipcluster instance yourself, and using the 
`--ipcluster` argument to ipyrad to tell it to find your running ipcluster 
instance. 

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
    user@compute$ sleep 60

    ## Then run ipyrad like normal but with --ipcluster so it knows to look for 
    ## your specific ipcluster instance.
    user@compute$ ipyrad -p params-test.txt -s 2 --ipcluster



