
Cookbook for running BUCKy in parallel in a Jupyter notebook
------------------------------------------------------------

This notebook uses the *Pedicularis* example data set from the first
empirical ipyrad tutorial. Here I show how to run BUCKy on a large set
of loci parsed from the output file with the ``.loci`` ending. All code
in this notebook is Python. You can simply follow along and execute this
same code in a Jupyter notebook of your own.

Software requirements for this notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All required software can be installed through conda by running the
commented out code below in a terminal.

.. code:: ipython2

    ## conda install -c BioBuilds mrbayes
    ## conda install -c ipyrad ipyrad
    ## conda install -c ipyrad bucky

.. code:: ipython2

    ## import Python libraries
    import ipyparallel as ipp
    import subprocess as sps
    import ipyrad as ip
    import glob
    import os
    import ipyrad.file_conversion as ifc

Cluster setup
~~~~~~~~~~~~~

To execute code in parallel we will use the ``ipyparallel`` Python
library. A quick guide to starting a parallel cluster locally can be
found `here <link>`__, and instructions for setting up a remote cluster
on a HPC cluster is available
`here <http://ipyrad.readthedocs.io/HPC_Tunnel.html>`__. In either case,
this notebook assumes you have started an ``ipcluster`` instance that
this notebook can find, which the cell below will test.

.. code:: ipython2

    ## look for running ipcluster instance, and create load-balancer
    ipyclient = ipp.Client()
    lbview = ipyclient.load_balanced_view()
    print "{} engines found".format(len(ipyclient))


.. parsed-literal::

    4 engines found


.. code:: ipython2

    %%px
    ## push imports to parallel engines
    import subprocess as sps
    import glob
    import os

Input/output organization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    ## enter the data file for your analysis here
    LOCIFILE = "/home/deren/Documents/ipyrad/tests/branch-test/base_outfiles/base.loci"
    WORKDIR = "analysis-bucky"

.. code:: ipython2

    ## This ensures the file paths are Full Paths (not relative) 
    LOCIFILE = os.path.realpath(LOCIFILE)
    WORKDIR = os.path.realpath(WORKDIR)
    print "infile is:", LOCIFILE
    print "outdir is:", WORKDIR


.. parsed-literal::

    infile is: /home/deren/Documents/ipyrad/tests/branch-test/base_outfiles/base.loci
    outdir is: /home/deren/Documents/ipyrad/tests/analysis-bucky


Set up some tests
~~~~~~~~~~~~~~~~~

List the names of the samples you wish to include in your analysis.
BUCKy begins to perform less well when the number of tips is >10 or so,
so you might want to try focus your analysis on subsampled sets of taxa.
Here we select just 9 of the 13 samples in the data set, with just one
representative of each species or subspecies.

.. code:: ipython2

    ## make a list of sample names you wish to include in your BUCKy analysis 
    SUBSAMPLES = [
        "29154_superba", 
        "30686_cyathophylla", 
        "41478_cyathophylloides", 
        "33413_thamno", 
        "30556_thamno",
        "35236_rex",
        "40578_rex", 
        "38362_rex", 
        "33588_przewalskii",
    ]

Sample loci and write NEXUS files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    ## create a name for this particular data set
    NAME = "example"
    
    ## create nexus files for this data set
    ifc.loci2multinex(name=NAME, 
                      locifile=LOCIFILE, 
                      subsamples=SUBSAMPLES, 
                      minSNPs=2, 
                      outdir=WORKDIR)


.. parsed-literal::

    wrote 709 nexus files to /home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-example


An example nexus file
~~~~~~~~~~~~~~~~~~~~~

Nexus files are written to a new directory called ``bucky-{name}``,
where name is the name entered into the ``loci2multinex()`` function. If
you entered a ``outdir`` argument as well then this new directory will
be made as a subdirectory inside that outdir. Above we used
name="example" and outdir=WORKDIR, which created files in the directory
shown above.

.. code:: ipython2

    ## get RUNDIR relative to WORKDIR to ensure it is a Full Path
    RUNDIR = os.path.join(WORKDIR, "bucky-{}".format(NAME))
    
    ## print an example nexus file
    with open(os.path.join(RUNDIR, "1.nex")) as nex:
        print nex.read()


.. parsed-literal::

    #NEXUS
    begin data;
    dimensions ntax=9 nchar=66;
    format datatype=dna interleave=yes gap=- missing=N;
    matrix
    30686_cyathophylla      CTTGGCAGGTGGCAGTTCGTTGCTGTTATATGCTGTAAGAAAAT-AAAAAAAAATCACCTGTTTAG
    33413_thamno            CTTGGCAGGTGGCAGTTTGTTGCTGTTTTATGCTGTAAGAAAAT--AAAAAAAACCACCTGTTTAG
    30556_thamno            CTTNGCAGGTGGCAGTTTGTTGCTGTTTTATGCTGTAAGAAAAT-NAAAAAAAATCACCTGTTTAG
    33588_przewalskii       CTTGGCAGGTGGCAGTTCGTTGCTGAAATATGCTGTAAGAAAAT-AAAGAAAAATCATTT-TTTGG
    29154_superba           CTTGGCAGTTGGCATTTCGTTGCTGTTATATGCTGTAAGAAAAT-AAAAAAAAATCACCTGTTTAA
    40578_rex               CTTGGCAGGTGGCAGTTTGTTGCTGTTTTATGCTGTAAGAAAAT--AAAAAAAATCACCTGTTTAG
    41478_cyathophylloides  CTTGGCAGGTGGCAGTTCGTTGCTGTTATATGCTGTAAGAAAATAAAAAAAAAATCACCTGTTTAG
    38362_rex               CTTGGCAGGTGGCAGTTTGTTGCTGTTTTATGCTGTAAGAAAATAAAAAAAAAATCACCTGTTTAG
    35236_rex               CTTGGCAGGTGGCAGTTTGTTGCTGTTTTATGCTGTAAGAAAAT--AAAAAAAATCACCTGTTTAG
    
        ;
    
    begin mrbayes;
    set autoclose=yes nowarn=yes;
    lset nst=6 rates=gamma;
    mcmc ngen=2000000 samplefreq=1000 printfreq=2000000;
    sump burnin=1000000;
    sumt burnin=1000000;
    end;
    


.. code:: ipython2

    ## get all nexus files for this data set
    nexfiles = glob.glob(os.path.join(RUNDIR, "*.nex"))

A Python function to call ``mrbayes``, ``mbsum`` and ``bucky``.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    def mrbayes(infile):
        ## double check file path
        infile = os.path.realpath(infile)
        if not os.path.exists(infile):
            raise Exception("infile not found; try using a fullpath")
            
        ## call mrbayes
        cmd = ['mb', infile]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        stdout = proc.communicate()
        
        ## check for errors
        if proc.returncode:
            return stdout

.. code:: ipython2

    def mbsum(dirs):
        trees1 = glob.glob(os.path.join(dirs, "*.run1.t"))
        trees2 = glob.glob(os.path.join(dirs, "*.run2.t"))
        tidx = 0
        for tidx in xrange(len(trees1)):
            cmd = ["mbsum", 
                   "-n", "0", 
                   "-o", os.path.join(dirs, str(tidx))+".in", 
                   trees1[tidx], 
                   trees2[tidx]]
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            proc.communicate()
        print "summed {} trees in: {}".format(tidx, dirs)

.. code:: ipython2

    def bucky(outname, indir, alpha, nchains, nreps, niter):
        ## check paths
        if not os.path.exists(indir):
            raise Exception("infiles not found; try using a fullpath")
        
        ## call bucky 
        infiles = os.path.join(indir, "*.in")
        cmd = ["bucky", 
               "-a", str(alpha),
               "-c", str(nchains),
               "-k", str(nreps),
               "-n", str(int(niter)), 
               "-o", outname, 
               infiles]
        
        cmd = " ".join(cmd)
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE, shell=True)
        stdout = proc.communicate()
        if proc.returncode:
            return " ".join(cmd), stdout

Run mrbayes on all nexus files in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is important that the lists contain the full paths to the files.

.. code:: ipython2

    ## send jobs to the parallel engines
    asyncs = []
    for nexfile in nexfiles:
        async = lbview.apply(mrbayes, nexfile)
        asyncs.append(async)

Track progress of the mrbayes runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to check the progress interactively then execute the cell
below, which will tell you how many jobs have finished. The cell below
that uses a wait() statement to block progress until all of the mrbayes
jobs are finished.

.. code:: ipython2

    ready =  [i for i in asyncs if i.ready()]
    failed = [i for i in ready if not i.successful()]
    
    ## print progress
    print "mrbayes batch runs:"
    print "{} jobs submitted".format(len(asyncs))
    print "{} jobs finished".format(len(ready))
    
    ## print errors, if any.
    if any(failed):
        print failed[0].exception()
        print failes[0].result()


.. parsed-literal::

    mrbayes batch runs:
    722 jobs submitted
    35 jobs finished


.. code:: ipython2

    ## waits until all mrbayes runs are finished
    ipyclient.wait()




.. parsed-literal::

    True



Summarize the mrbayes posteriors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    ## run mbsum on each directory of tree files
    mbsum(RUNDIR1)
    mbsum(RUNDIR2)


.. parsed-literal::

    summed 9 trees in: /home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13
    summed 0 trees in: /home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp9


Run BUCKy to infer concordance factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    nchains = 4
    nreps = 4
    niter = 1e6
    alphas = [0.1, 1, 10]
    
    ## submit jobs to run at several values of alpha
    bsyncs = []
    for alpha in alphas:
        outname = os.path.join(RUNDIR, "bucky-{}".format(alpha))
        args = (outname, RUNDIR, alpha, nchains, nreps, niter)
        async = lbview.apply(bucky, *args)
        bsyncs.append(async)

Track progress of Bucky runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    ready =  [i for i in bsyncs if i.ready()]
    failed = [i for i in ready if not i.successful()]
    print "bucky batch runs:"
    print "{} jobs submitted".format(len(bsyncs))
    print "{} jobs finished".format(len(ready))
    if len(ready) == len(bsyncs):
        ## print errors, if any.
        if any(failed):
            print failed[0].exception()



.. parsed-literal::

    bucky batch runs:
    3 jobs submitted
    0 jobs finished


.. code:: ipython2

    ipyclient.wait()




.. parsed-literal::

    True



Results
~~~~~~~

Look at individual results files for final stats.

.. code:: ipython2

    results = glob.glob(os.path.join(RUNDIR, "bucky-*.concordance"))


.. code:: ipython2

    results




.. parsed-literal::

    ['/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-1.txt.concordance',
     '/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-0.1.concordance',
     '/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-0.1.txt.concordance',
     '/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-1.concordance',
     '/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-10.txt.concordance',
     '/home/deren/Documents/ipyrad/tests/analysis-bucky/bucky-samp13/bucky-10.concordance']


