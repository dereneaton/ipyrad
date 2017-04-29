.. include:: global.rst  

.. _analysis:


Maximum-likelihood phylogenetic inference
=========================================
RAxML is a standard tool for phylogenetic inference popular for its speed and 
ease of use, and is among the most commonly used software for analyzing RAD-seq 
alignments. While it is easy to run a multi-threaded version of RAxML, which 
can take advantage of many threads on a single machine, it is a bit more difficult
to optimize a run that is parallized over many connected machines on 
a HPC cluster, which requires the MPI version. Below we list a few common 
commands for analyzing large RAD-seq alignments in RAxML. This is not an 
exhaustive tutorial, just a quick quide. 

More information about RAxML can be found in the  
`v.8.0 raxml documentation <https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwj4opaa1IjRAhVB1oMKHajPAMAQFggcMAA&url=https%3A%2F%2Fbioinformatics.oxfordjournals.org%2Fcontent%2Fsuppl%2F2014%2F01%2F18%2Fbtu033.DC1%2FNewManual.pdf&usg=AFQjCNH_8fbJI7fBU6yVL74UFKRzZhftFg&sig2=3GfktJYcAdFcSxRWs0TgFw>`__ 
and on the google group `raxml forum <https://groups.google.com/forum/#!topic/raxml/>`__. 


Installing raxml on a cluster
-----------------------------
There are many versions of raxml available and the one on your system may not 
be up to date. You can ask your administrator to install the latest version, or 
install it yourself *locally* (you do not need administrative privileges to do 
so.) I usually recommend using `conda`, which makes it quite easy to install: 

.. code:: bash

    ## one way of installing raxml is with conda
    conda install raxml -c bioconda


However, you will probably be able to get a bit faster performance if you 
build raxml from source on your machine, since conda does not yet handle 
well checking for various threading/compiling options. 

As an alternative to using conda, or your HPC system's version of raxml
(if there is one), the code below can be used to install raxml locally. 
I will assume that your machine has AVX2 mode available. 
This installation will put the executables in a local directory called 
`~/local/bin/` which you will want to add to your $PATH 
(add it to your .bashrc file).

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
independent bootstrap analyses, which later need to be combined using the 
options in raxml to do this. The HYBRID approach makes use of MPI to distribute
threaded jobs across different compute nodes, but I think it's pretty 
tricky to get working right.


Running raxml (threaded) on a single node
------------------------------------------
This is the most common code I use to analyses RAD-seq alignments. It runs 
the (-f a) method, which performs the *standard hill-climbing algorithm* to 
find the best scoring ML tree *and* it performs a rapid bootstrap
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
in the outgroup comma-separated. Setting the outgroups makes your life easier.
If you do not set the outgroup but try to re-root your tree later the node labels 
indicating bootstrap support values can easily become misplaced in many 
tree plotting programs. 


Submitting an MPI HYBRID job to run on a cluster (experimental)
---------------------------------------------------------------
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
    module load OpenMPI

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
There are many ways to do this. I wrote a Python program called toytree which
I prefer, but another popular alternative is the *ape* package in R. 

To make a simple tree plot in python use the code below. 

.. code:: python

    ## load modules
    import toyplot
    import toytree

    ## draw a tree
    tre = toytree.tree("raxout/RAxML_bipartitions.name.tre")
    canvas, axes = tre.draw(
        width=400,
        node_labels=tre.get_node_values('support'),
        )

    ## save the tree
    toyplot.html.render(canvas, "mytree.html")
    import toyplot.pdf
    toyplot.pdf.render(canvas, "mytree.pdf")



To make a simple tree plot in R use the code below. 

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

    

    
