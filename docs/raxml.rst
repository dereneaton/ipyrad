.. include:: global.rst  

.. _analysis:


Maximum-likelihood phylogenetic inference
=========================================
The program RAxML is a standard tool for phylogenetic inference
popular for its speed and ease of use. It offers a huge
variety of analysis methods, and you could spend hours exploring its 
massive documentation. Optimizing it for use on a cluster can be a 
bit a tricky. Here I list a few tips for working with 
large concatenated RAD-seq alignments (.phy output file) and analyzing
them on a HPC cluster based on information from the 
`v.8.0 raxml documentation <https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwj4opaa1IjRAhVB1oMKHajPAMAQFggcMAA&url=https%3A%2F%2Fbioinformatics.oxfordjournals.org%2Fcontent%2Fsuppl%2F2014%2F01%2F18%2Fbtu033.DC1%2FNewManual.pdf&usg=AFQjCNH_8fbJI7fBU6yVL74UFKRzZhftFg&sig2=3GfktJYcAdFcSxRWs0TgFw>`__


Installing raxml on a cluster
-----------------------------
There are many versions of raxml, and it is updated frequently, so the 
version that is lying around on your cluster may very well be outdated, 
or not the version that is best for you. You can ask your administrator
to install the latest version, or install it yourself *locally*
(you do not need administrative privileges for this.)
The code below installs three versions, the PTHREADS (threaded version), 
MPI (can use processors from different nodes), and Hybrid
(a mix of the first two). This installation will put the executables
in a local directory called `~/local/bin/`.

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

    ## compile the AVX2.PTHREADS version of raxml
    rm *.o
    make -f Makefile.AVX2.PTHREADS.gcc

    ## compile the hybrid version
    rm *.o
    make -f Makefile.AVX2.HYBRID.gcc

    ## (optional) copy the binary to your binaries dir
    cp raxml-MPI-AVX2 ~/local/bin
    cp raxml-PTHREADS-AVX2 ~/local/bin
    cp raxml-HYBRID-AVX2 ~/local/bin


Why multiple versions?
---------------------------
If you only plan to use a single compute node on your cluster then you should 
just use the PTHREADS (threaded) version, as this will most efficiently make use
of the cores on that node. The MPI version is needed to make use of cores spread
across multiple nodes, however, it only offers a subset of the functions that
are available in the threaded version. Mostly it is used for distributing many
independent bootstrap analyses. The HYBRID approach makes use of MPI to distribute
threaded jobs across different compute nodes. This can be the most efficient method
but can also be a bit tricky to get working.


Running raxml (threaded) on a single node
------------------------------------------
This code run the (-f a) method, which performs the *standard hill-climbing
algorithm* to find the best scoring ML tree *and* it performs a rapid bootstrap
analysis. We tell it how many bootstraps with the -N option.

.. code:: bash

    ## this is an example call to run raxml tree inference w/ bootstrapping
    raxmlHPC-PTHREADS-AVX2 -f a \              ## do rapid-bootstrapping & full search
                      -T 20 \                  ## number of threads available
                      -m GTRGAMMA \            ## use GTRGAMMA model
                      -N 100 \                 ## 100 searches from parsimony start trees
                      -x 12345 \               ## bootstrap random seed 
                      -p 54321 \               ## parsimony random seed
                      -n outname \             ## a name for your output files
                      -w outdir \              ## a directory for your output files
                      -s inputfile.phy \       ## your sequence alignment
                      -o outgroup1,outgroup2   ## set your outgroups!



Running raxml (HYBRID) across multiple nodes
--------------------------------------------
The HYBRID version of raxml is best used for large-scale bootstrapping when you 
have access to many cores spread across multiple compute nodes. 
Because this version uses MPI you must call an MPI executable 
(e.g., mpiexec or mpirun) before the command to specify the number of nodes 
and then -T to specify the number of threads per node. It is best that you are
connected to many cores with the same number of cores. 


.. code:: bash

    ## this is an example call to run raxml tree inference w/ bootstrapping
    mpiexec -np 4 raxmlHPC-HYBRID-AVX2 -f a \    ## do rapid-bootstrapping & full search
                      -T 20 \                    ## number of threads available
                      -m GTRGAMMA \              ## use GTRGAMMA model
                      -N 100 \                   ## 100 searches from parsimony start trees
                      -x 12345 \                 ## bootstrap random seed 
                      -p 54321 \                 ## parsimony random seed
                      -n outname \               ## a name for your output files
                      -w outdir \                ## a directory for your output files
                      -s inputfile.phy \         ## your sequence alignment
                      -o outgroup1,outgroup2     ## set your outgroups!


Should I use the GTRCAT model?
------------------------------
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


Submitting jobs to run on a cluster
-----------------------------------
The method we are using will distribute 100 replicate analyses across all of the
cores you are connected to (including across multiple nodes) using MPI. 
But we need to make sure we tell the program explicitly how we many cores 
will be available. If you have ipyrad installed then you will already have 
MPI installed (just type mpiexec), but your system probably has a version 
installed as well.

Below is an example SLURM (sbatch) submission script, you can make something similar
but slightly different for other systems such as TORQUE (qsub). Save the file 
with a name like *raxml-script.sh*. 

.. code:: python

    #!/bin/bash
    # set the number of nodes and processes per node
    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 20
    #SBATCH --exclusive
    #SBATCH --time 10-00:00:00
    #SBATCH --mem-per-cpu 4000
    #SBATCH --job-name raxml-0
    #SBATCH --output raxml-0

    ## make sure you're in your home directory
    cd $HOME

    ## you can load a system-wide MPI module if available
    #module load MPI/OpenMPI

    ## call mpiexec and raxml, use -np for number of cores.
    ~/miniconda/bin/mpiexec -np 4 ~/local/bin/raxml-HYBRID-AVX2 \
                      -T 20 \
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

    

    
