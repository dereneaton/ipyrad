



.. _HPCscript:

Running ipyrad on a cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^

High performance computing (HPC) clusters are accessible to most 
users with University affiliations, providing cheap or free access
to dozens or hundreds of computing cores, which is key for 
assembling genomic data set quickly. The HPC software set up often varies 
between different institutions but we provide here a general workflow 
for running ipyrad on TORQUE systems, which are those that 
use the ``qsub`` submission format. 

When you login to an HPC cluster you will be connected to a
``login`` node, and from there you can submit ``jobs`` which 
are sent to processing nodes. That's where the heavy computing happens.
Here we show two ways to run ipyrad, interactively, and by submitting job
scripts. If it's available, interactive jobs are easier to run, 
but job scripts are usually better suited for very long running jobs. 
I suggest that you start by trying interactive mode, since it
is much better for trouble shooting.

We've designed ipyrad to be very easy to run on HPC setups. 
Because the software is installed locally (in your miniconda directory)
you can simply call ipyrad the same as you would on your local
machine. No need to ask you administrator to install the software
globally, and no need to load software modules before running 
ipyrad. 

The first step to running ipyrad on an HPC is to connect 
from your local terminal to 
the HPC cluster, usually involving an `ssh` script. 

.. parsed-literal::

    ## connect to the login node from your computer
    >>> ssh user@address_of_HPC_cluster


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
    >>> qsub -I -l nodes=8:ppn=8:walltime=24:00:00 -q "fas_general"
    

Submitting job scripts
----------------------
Because it sometimes takes a while to connect to an interactive
node, it is common practice to instead submit scripts that will
be run whenever the node becomes available. Here is an example
qsub script that could make with a text editor and then 
save with a name such as ``qsub_script.sh``:

.. parsed-literal::

    #!/bin/sh

    #PBS -N ipyrad-test
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=1:ppn=8

    ## change into whichever directory you will run code from
    cd $PBS_O_WORKDIR

    ## call ipyrad on your params file
    ipyrad -p params-test.txt -s 1 


Then you submit this script to a queue with a qsub command such as:

.. parsed-literal::

    ## submit the qsub script
    >>> qsub qsub_script.sh


Here is another example script that would connect to more processors
spread across multiple nodes:

.. parsed-literal::

    #!/bin/sh

    #PBS -N ipyrad-test
    #PBS -j oe
    #PBS -m ae
    #PBS -M youremail@institution.edu
    #PBS -q queue_name
    #PBS -l nodes=8:ppn=8

    ## change into whichever directory you which to run code from
    cd $PBS_O_WORKDIR

    ## call ipyrad on your params file
    ipyrad -p params-test.txt -s 1 -c 64 --MPI


.. parsed-literal::

    ## submit the qsub script
    >>> qsub qsub_script.sh




