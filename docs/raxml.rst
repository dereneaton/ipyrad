.. include:: global.rst  

.. _analysis:


Maximum-likelihood phylogenetic inference
=========================================

To program RAxML has become a standard tool for phylogenetic inference
due to its speed and easy to use command-line interface. There is a huge
variety of types of analyses that can be done in RAxML, many of which you
can see by using the help command ``-h``. A good place to start for any
phylogenetic question is by performing a ML inference on a concatenated
supermatrix of the data. The code below can be used to perform rapid 
bootstrap resampling in addition to a full tree inference. 

.. code:: python

    raxmlHPC-PTHREADS-AVX -f a \ 
                          -m GTRGAMMA \
                          -N 100 \
                          -x 12345 \
                          -p 54321 \
                          -T 20 \
                          -n outname \
                          -w outdir \
                          -s inputfile.phy \
                          -o outgroup1,outgroup2,outgroup3  


Which model do I use?
---------------------
It is fairly standard to use the GTR+GAMMA model in RAxML due to its fast speed. 
I tend to use the GTR-GAMMA model over GTRCAT even though it is slower. I believe
GTRCAT performs poorly if you do not have a very large number of taxa in the data set. 

  
Setting an outgroup
-------------------
If you have prior information about which clade is the outgroup I recommend 
setting this value in the command string. It should include each sampled taxon
in the outgroup comma-separated. If you do not set the outgroup and then 
try to re-root the tree later the node labels indicating bootstrap support
values will be misplaced. 


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

    

    
