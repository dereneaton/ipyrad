.. include:: global.rst  

.. _analysis:


Maximum-likelihood phylogenetic inference
=========================================
The program RAxML is a standard tool for phylogenetic inference
popular for its speed and ease of use. It offers a huge
variety of analysis methods, and you could spend hours exploring its 
massive documentation. Here I list a few tips for working with 
large concatenated RAD-seq alignments (.phy output file) and analyzing
them on a HPC cluster. 


Installing raxml on a cluster
-----------------------------
Many old versions of raxml are floating around the internet and your cluster
almost certainly does not have a recent version installed. So the best thing
to do is install it yourself (you do not need administrative privileges for this.)
We will install the MPI-enabled version of raxml in a local directory, which
you can remove at any time if you like. 

.. code:: bash  

    ## (optional) create directories to store the software. 
    ## I use ~/local/src to store source code and 
    ## I use ~/local/bin to store binaries.
    mkdir -p ~/local/src ~/local/bin
    
    ## cd to where you want to store the raxml source code. 
    cd ~/local/src

    ## use git to clone the raxml github repo into your src dir
    git clone https://github.com/stamatak/standard-RAxML.git

    ## now cd into the raxml directory
    cd standard-RAxML.git

    ## compile the AVX2.MPI version of raxml
    make -f Makefile.AVX2.MPI.gcc

    ## (optional) copy the binary to your binaries dir
    cp raxml-MPI-AVX2 ~/local/bin

    


Running raxml 
----------------------------------------
Below I show how to submit a job script to run raxml, but first let's 
go though some of the options.

.. code:: bash

    raxmlHPC-MPI-AVX2 -f a \                   ## do rapid-bootstrapping & full search
                      -m GTRGAMMA \            ## use GTRGAMMA model
                      -N 100 \                 ## 100 searches from parsimony start trees
                      -x 12345 \               ## bootstrap random seed 
                      -p 54321 \               ## parsimony random seed
                      -n outname \             ## a name for your output files
                      -w outdir \              ## a directory for your output files
                      -s inputfile.phy \       ## your sequence alignment
                      -o outgroup1,outgroup2   ## set your outgroups!


Should I use the GTRCAT?
------------------------
GTRCAT is a speed improvement for modeling rate variation under the GTRGAMMA model. 
It is particularly designed for modeling rate heterogeneity across very large trees
(e.g., hundreds of taxa), and is not recommended for smaller trees. In fact the raxml
docs state in **bold font** that using it for less than 50 taxa is a bad idea. If your 
tree has >100 taxa then I would say go for it.


Setting an outgroup
-------------------
If you have prior information about which clade is the outgroup I recommend 
setting this value in the command string. List each sampled taxon
in the outgroup comma-separated. Setting the outgroups now makes your life easier.
If you do not set the outgroup but try to re-root your tree later the node labels 
indicating bootstrap support values can easily become misplaced. 


Running parallel raxml on a cluster
-----------------------------------
The method we are using will distribute 100 replicate analyses across all of the
cores you are connected to (including across multiple nodes) using MPI (a 
way of sharing information between computers). But we need to make sure we tell
the program explicitly how we many cores will be available. If you have ipyrad 
installed then you will already have MPI installed (just type mpiexec), but your
system probably has a version installed as well.

Below is an example SLURM (sbatch) submission script, you can make something similar
but slightly different for other systems such as TORQUE (qsub). Save the file 
with a name like *raxml-script.sh*. 

SLURM (sbatch) example.      
.. code:: python
    #!/bin/bash
    # set the number of nodes and processes per node
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 8
    #SBATCH --exclusive
    #SBATCH --time 10-00:00:00
    #SBATCH --mem-per-cpu 2000
    #SBATCH --job-name raxml-0
    #SBATCH --output raxml-0

    ## make sure you're in your home directory
    cd $HOME

    ## call mpiexec and raxml, use -np for number of cores.
    ~/miniconda/bin/mpiexec -np 32 ~/local/bin/raxml-MPI-AVX2 \
                      -f a \                   
                      -m GTRGAMMA \            
                      -N 100 \                 
                      -x 12345 \               
                      -p 54321 \              
                      -n raxml-0 \      
                      -w raxml_runs/ \         
                      -s test/test_outfiles/test.phy \    
                      -o outg1,outg2


Then submit the job to run on your cluster. Sometimes you will have to 
add additional arguments to the submission script, such as the name of the
queue that you are submitting to. 

.. code:: python 
    sbatch raxml_script.sh


Plotting the trees
------------------
There are many ways to do this. I prefer using the *ape* package in R. 
You can make a simple tree plot with the code below. 


.. code:: python

    ## load the ape library
    library(ape)

    ## read in the tree file
    tre <- read.tree("raxout/RAxML_bipartitions.name.tre")

    ## ladderize the tree (makes it prettier)
    ltre <- ladderize(tre)

    ## plot the tree with bootstrap support on node labels 
    plot(ltre, cex=0.7)
    nodelabels(ltre$node.label, cex=0.7)

    

    
