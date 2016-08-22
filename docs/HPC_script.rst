



.. _HPCscript:

Command-line interface
^^^^^^^^^^^^^^^^^^^^^^

High performance computing (HPC) clusters are accessible to most 
users with University affiliations, providing cheap or free access
to dozens or hundreds of computing cores allowing for very fast
assembly of genomic data. HPC set ups often vary between 
different institutions but we provide a general recommended 
workflow here that works for most common setups.  

When you login to an HPC cluster you are usually connected to 
``login`` node, and from there you can submit ``jobs`` which 
are sent to processing nodes. That's where the heavy computing happens.
Here we show two ways to run ipyrad, interactively, and through job
submission. If it's available, interactive jobs are easier to run, 
but very long running jobs are better suited for submitting as 
scripts. 

The first step is to use connect from your local terminal to 
the HPC cluster, usually involving an `ssh` script. 

.. parsed-literal::

    ## connect to the login node from your computer
    ssh user@address_of_HPC_cluster

From there you can then try logging in interactively to a 
computing node by using the command `qsub` with the ``-I`` 
argument. Sometimes you have to provide other default 
arguments such as the name of the queue you are connecting to. 
This information should be available from your institution. 

.. parsed-literal::

    ## try logging into an interactive node, depending
    ## on the queue this may be instant, or may take hours
    qsub -I 


You can use all of the typical qsub arguments here to connect
to more CPUs, or specific queues. 

    ## ask for 64 cores across 8 nodes from queue 'fas_general' 
    ## and request 24 hours of wall time.
    qsub -I -l nodes=8:ppn=8:walltime=24:00:00 -q "fas_devel"
    

Because it sometimes takes a while to connect to an interactive
node, it is common practice to instead submit scripts that will
be run whenever the node becomes available. Here is an example
qsub script which is saved as ``qsub_script.sh``:

.. parsed-literal::

#!/bin/sh

#PBS -N ipyrad-test
#PBS -j oe
#PBS -m ae
#PBS -M deren.eaton@yale.edu
#PBS -q fas_normal
#PBS -l nodes=1:ppn=8

## change into whichever directory you which to run code from
cd $PBS_O_WORKDIR

## call ipyrad on your params file
ipyrad -p params-test.txt -s 1 


Then you submit this script to a queue with a qsub command such as:

.. parsed-literal::

    ## submit the qsub script
    qsub qsub_script.sh


Here is another example script that would connect to more processors
spread across multiple nodes:

.. parsed-literal::

#!/bin/sh

#PBS -N ipyrad-test
#PBS -j oe
#PBS -m ae
#PBS -M deren.eaton@yale.edu
#PBS -q fas_normal
#PBS -l nodes=8:ppn=8

## change into whichever directory you which to run code from
cd $PBS_O_WORKDIR

## call ipyrad on your params file
ipyrad -p params-test.txt -s 1 -c 64 --MPI





