

==============================
Analysis: Overview
==============================

The ipyrad-analysis toolkit is a Python interface for running a suite of 
evolutionary analysis tools with convenient features for filtering for
missing data, grouping individuals into populations, dropping samples, 
and more. This makes it easy to run many different analyses from a single
assembled dataset where filtering is applied at the analysis stage. 
These tools are designed to be implemented in jupyter notebooks, 
creating a reproducible document for all of your analyses. 


.. image:: ./images/ipa-flow.svg


Common Syntax
-------------
All of these tools share a common syntax making them easy to use without 
having to worry about creating different input files, or learning new 
file formats. 

.. code-block:: python

   # the analysis tools are a subpackage of ipyrad
   import ipyrad.analysis as ipa

   # a large suite of tools are available 
   tool = ipa.structure(data="./outfiles/data.snps.hdf5")

   # all tools share a common syntax for distributing work in parallel.
   tool.run()



imap sample filtering
---------------------
All tools allow you to enter a dictionary object to select which samples 
will be included in the analyses. Any sample not present as a value in 
the dictionary will be excluded. For any tools that use populations assignments
samples can be grouped into populations by using key:value pairs like below.

Common usage:

.. code-block:: python

	# create an imap dictionary with only the samples you want to include
	IMAP = {'pop1': ['name1', 'name2', 'name3'], 'pop2':  ['name4', 'name5', 'name6']}

	# init tool with data, selected samples in imap, and filtering option
	tool = ipa.window_extracter(data="test.seqs.hdf5", imap=IMAP, mincov=0.9)

	# this will generate results for only the 6 selected samples with filters applied.
	tool.run()


Another commonly useful trick:

.. code-block:: python

	# create a list of sample names you wish to exclude
	EXCLUDE = ['A', 'B', 'C']

	# init a simple tool that will be used only to access the sample names
	null = ipa.snps_extracter(data='test.snps.hdf5')

	# create an imap with names from the null tool, but skipping those in EXCLUDE
	IMAP = {'keep': [i for i in null.names if i not in EXCLUDE]}

	# now init the tool you wish to run using this imap
	tool = ipa.pca(data='test.snps.hdf5', imap=IMAP)



Site/SNP filtering 
------------------
The ipyrad analysis tools are designed specifically to help you deal with missing
data in your analyses. This is done first by allowing you to apply **filtering** 
when you start your analyses. 


**mincov (minimum total sample coverage)**: The `mincov` option is used to filter 
sites for missing data. It sets the minimum threshold for the number of samples that 
must have data at a site for it to be retained in the analysis. This can be set
either as an integer (e.g., `2` = at least 2 samples must have data), or as a proportion
(e.g., `0.2` = at least 20% of samples must have data).


**minmap (minimum imap sample coverage)**: The `minmap` option acts the same as 
the `mincov` parameter, however, it applies to every population in the imap
dictionary separately. If any population in the imap dictionary has a sample 
coverage below the minimum for that population in the `minmap` dictionary then
the site will be removed. This allows you to ensure there is equal coverage
across different populations or species. For example, requiring that every
site in your analyses has data from at least 5 individuals in each of 5 different
populations. This can be set either as an integer (e.g., `2` = at least 2 samples
must have data in every population in `imap`), or it can be set using a proportion
(e.g., `0.2` = at least 20% of samples must have data in each population); 
or, different values can be set for every population in imap using a dictionary
(e.g., `minmap={'a': 0.2, 'b': 0.5, 'c': 0.1}`).

**minmaf (minimum minor allele frequency)**: The `minmaf` option applies a 
minor allele frequency filter. This can remove singleton variants from your
dataset which may be more likely to represent errors, or to be of little 
power to certain analyses. *There are many programs for which this type of
filter is inappropriate*, and we usually do not provide it as an option
if we think it should not be used. It can be useful for analyses like
PCA and STRUCTURE to reduce noise. 



Simple and massive parallelization
----------------------------------
The ipyrad-analysis tools use the same parallelization framework as
ipyrad, based on the `ipyparallel <https://ipyparallel.rtfd.io>`_
package. This allows for flexible parallelization, whether working
on a laptop or a massive computing cluster. Advanced users can 
start their own ipcluster instances and connect to them in each tool.
But generally users can acheive any desired outcome using our simple
automated workflow within the `.run()` function of each tool.


In this example the treeslider tool will distribute hundreds of 
tree inference jobs in parallel, by creating a queue and starting
each job as the previous one finishes so that it does not use
more resources that are available. Let's imagine you are on a workstation or 
cluster with 80 cores available. If you simply called `.run(auto=True)`
the default behavior will be to use all 80 cores. If instead you want
it to use only 40 cores then you can set this will 'cores' in 
the `.ipcluster` dict of the tool object. If the tool that is being called 
uses multi-threading then you can also specify the number of threads
per job in the `ipcluster` dict. In this case it will run 4 simultaneous 
raxml jobs each using 10 threads each. If you do not set the `threads` 
option then it will use a default option set on each tool 
(usually nthreads=ncores).

.. code-block:: python

	# init a tool
	tool = ipa.treeslider(data="data.seqs.hdf5")

	# set the parallelization strategy
	tool.ipcluster['cores'] = 40
	tool.ipcluster['threads'] = 10

	# run on a computer with at least 40 cores
	tool.run(auto=True)


Simulation Tutorials
--------------------
In addition to empirical examples of each tutorial we also provide 
notebooks that test each tool on simulated data. These notebooks begin
with a few blocks of code to setup the simulation scenario using the 
`ipcoal <https://ipcoal.rtfd.io>`_ Python package. When you click on any 
of the links in the tutorials section the notebook will open in a page
called nbviewer, and in the upper right of the page there is a circular
icon that says "open in binder" if you hover over it. This will open
an interactive cloud-based version of the notebook. This is a great 
way to explore the tools yourself without having to install anything.
You can modify the simulation scenario to match your empirical data
and use this to setup expectations for your results, or explore biases
that can be introduced by ILS, or missing data, or sampling methods.
