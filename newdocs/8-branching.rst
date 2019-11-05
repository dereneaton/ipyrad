
.. _branching_workflow:

Assembly: Branching and Merging
===============================

Branching can be used to create a new assembly with a different name and 
to which you can apply new parameters and create new downstream files ...

Merging can be used to combine samples from multiple libraries.


.. _dropping_samples:

Drop samples by branching
--------------------------
Another point where branching is useful is for adding or dropping
samples from an assembly, either to analyze a subset of samples 
separately from others, or to exclude samples with low coverage. 
The `branching` and `merging` fuctions in ipyrad make this easy. 
By requiring a branching process in order to drop samples from an assembly ipyrad inherently forces you to retain the parent assembly as a copy. This provides a nice fail safe so that you can mess around with your new branched assembly without affecting it's pre-branched parent assembly. 

Examples using the ipyrad CLI

.. code-block:: bash

    ## branch and only keep 3 samples from assembly data1
    >>> ipyrad -n data1 -b data2 1A0 1B0 1C0

    ## and/or, branch and only exclude 3 samples from assembly data1
    >>> ipyrad -n data1 -b data3 - 1A0 1B0 1C0


Examples using the ipyrad Python API 

.. code-block:: bash

    ## branch and only keep 3 samples from assembly data1
    >>> data1.branch("data2", subsamples=["1A0", "1B0", "1C0"])

    ## and/or, branch and only exclude 3 samples from assembly data1
    >>> keep_list = [i for i in data1.samples.keys() if i not in ["1A0", "1B0", "1C0"]]
    >>> data1.branch("data3", subsamples=keep_list)


Merge Samples or Libraries
---------------------------

There are a number of ways to enter your data into *ipyrad* and we've
tried to make it as easy as possible to combine data from multiple
libraries and multiple plates in a simple and straightforward way. Here
we demonstrate a number of ways to demultiplex and load data under
different scenarios:

1. `One Library One Lane of
   sequencing <#one-library-one-lane-of-sequencing>`__
2. `One Library Multiple lanes of
   sequencing <#one-library-multiple-lanes-of-sequencing>`__
3. `Multiple libraries Multiple lanes of
   sequencing <#multiple-libraries-multiple-lanes-of-sequencing>`__
4. `Separate multiple libraries from one lane of
   sequencing <#separate-multiple-libraries-from-one-lane-of-sequencing>`__
5. `Alternative: Doing all of this with the API instead of the
   CLI <#alternative:-using-the-ipyrad-api-to-do-these-things>`__

1. One library One Lane of sequencing
-------------------------------------

First create a new Assembly, here we'll call it ``demux1``. Then use a
text-editor to edit the params file to enter the ``raw_fastq_path`` and
the ``barcodes_path`` information that is needed to demultiplex the
data. To automate the process of editing the params file I use the
command-line program ``sed`` here to substitute in the new values.

.. code:: python

    ## create a new assembly
    ipyrad -n demux1


.. parsed-literal::

    
      New file 'params-demux1.txt' created in /home/deren/Documents/ipyrad/tests
    


.. code:: python

    ## edit the params file to enter your raw_fastq_path and barcodes path
    sed -i '/\[2] /c\ipsimdata/rad_example_R1_.fastq.gz  ## [2] ' params-demux1.txt
    sed -i '/\[3] /c\ipsimdata/rad_example_barcodes.txt  ## [3] ' params-demux1.txt


.. parsed-literal::

    

.. code:: python

    ## run step 1 to demultiplex the data
    ipyrad -p params-demux1.txt -s 1 


.. parsed-literal::

    
     -------------------------------------------------------------
      ipyrad [v.0.5.15]
      Interactive assembly and analysis of RAD-seq data
     -------------------------------------------------------------
      loading Assembly: demux1
      from saved path: ~/Documents/ipyrad/tests/demux1.json
      New Assembly: demux1
      host compute node: [40 cores] on tinus
    
      Step 1: Demultiplexing fastq data to Samples
    
      [####################] 100%  sorting reads         | 0:00:06  
      [####################] 100%  writing/compressing   | 0:00:00  
    


The demultiplexed data is now located in the directory
``<project_dir>/<assembly_name>/``, which in this case is in
``./demux1_fastqs/``. The Assembly ``demux1`` knows the location of the
data, and so from here you can proceed in either of two ways. (1) You
simply continue on to step 2 using this Assembly object (demux1), or (2)
You create a new 'branch' of this Assembly, which will start by reading
in the ``sorted_fastq_data``. The latter is sometimes more clear in that
you keep separate the demultiplexing steps from the assembly steps. It
does not make a difference in this example, where we have only one
library and one lane of data, but as you will see in the examples below,
that it is sometimes easier to create multiple separate demux libraries
that are then merged into a single Object for assembling.

.. code:: python

    ## option 1: continue to assemble this data set
    ipyrad -p params-demux1 -s 234567

.. code:: python

    ## OR, option 2: create a new Assembly and enter path to the demux data
    ipyrad -n New
    
    ## enter path to the 'sorted_fastq_data' in params
    sed -i '/\[4] /c\./demux1_fastq/*.gz  ## [2] ' params-New.txt
    
    ## assemble this data set 
    ipyrad -p params-New.txt -s 1234567

2. One Library Multiple Lanes of Sequencing
-------------------------------------------

There are two options for how to join multiple lanes of sequence data
that are from the same library (i.e., there is only one barcodes file).
(1) The simplest way is to simply put the multiple raw fastq data files
into the same directory and select them all when entering the
``raw_fastq_path`` using a wildcard selector (e.g., "\*.fastq.gz"). (2)
The second way is to create two separate demux Assemblies and the merge
them, which I demonstrate below. Because the two demultiplexed lanes
each use the same barcodes file the Samples will have identical names.
*ipyrad* will recognize this during merging and read both input files
for each Sample in step 2.

.. code:: python

    ## create demux Assembly object for lane 1 
    ipyrad -n lane1raws 


.. parsed-literal::

    
      New file 'params-lane1raws.txt' created in /home/deren/Documents/ipyrad/tests
    


.. code:: python

    ## create demux Assembly object for lane 2 
    ipyrad -n lane2raws


.. parsed-literal::

    
      New file 'params-lane2raws.txt' created in /home/deren/Documents/ipyrad/tests
    


.. code:: python

    ## edit the params file for lane1 to enter its raw_fastq_path and barcodes file
    sed -i '/\[2] /c\ipsimdata/rad_example_R1_.fastq.gz  ## [2] ' params-lane1raws.txt
    sed -i '/\[3] /c\ipsimdata/rad_example_barcodes.txt  ## [3] ' params-lane1raws.txt
    
    ## edit the params file for lane2 to enter its raw_fastq_path and barcodes file
    sed -i '/\[2] /c\ipsimdata/rad_example_R1_.fastq.gz  ## [2] ' params-lane2raws.txt
    sed -i '/\[3] /c\ipsimdata/rad_example_barcodes.txt  ## [3] ' params-lane2raws.txt


.. parsed-literal::

    

.. code:: python

    ## demultiplex lane1
    ipyrad -p params-lane1raws.txt -s 1 


.. parsed-literal::

    
     -------------------------------------------------------------
      ipyrad [v.0.5.15]
      Interactive assembly and analysis of RAD-seq data
     -------------------------------------------------------------
      New Assembly: lane1raws
      host compute node: [40 cores] on tinus
    
      Step 1: Demultiplexing fastq data to Samples
    
      [####################] 100%  sorting reads         | 0:00:06  
      [####################] 100%  writing/compressing   | 0:00:01  
    


.. code:: python

    ## demultiplex lane2
    ipyrad -p params-lane2raws.txt -s 1 


.. parsed-literal::

    
     -------------------------------------------------------------
      ipyrad [v.0.5.15]
      Interactive assembly and analysis of RAD-seq data
     -------------------------------------------------------------
      New Assembly: lane2raws
      host compute node: [40 cores] on tinus
    
      Step 1: Demultiplexing fastq data to Samples
    
      [####################] 100%  sorting reads         | 0:00:06  
      [####################] 100%  writing/compressing   | 0:00:00  
    


.. code:: python

    ## merge the two lanes into one Assembly named both
    ipyrad -m both params-lane1raws.txt params-lane2raws.txt


.. parsed-literal::

    
    
     -------------------------------------------------------------
      ipyrad [v.0.5.15]
      Interactive assembly and analysis of RAD-seq data
     -------------------------------------------------------------
    
      Merging assemblies: ['params-lane1raws.txt', 'params-lane2raws.txt']
      loading Assembly: lane1raws
      from saved path: ~/Documents/ipyrad/tests/lane1raws.json
      loading Assembly: lane2raws
      from saved path: ~/Documents/ipyrad/tests/lane2raws.json
    
      Merging succeeded. New params file for merged assembly:
    
        params-both.txt
    


.. code:: python

    ## print merged stats of new Assembly
    ipyrad -p params-both.txt -r 


.. parsed-literal::

    
    Summary stats of Assembly both
    ------------------------------------------------
          state  reads_raw
    1A_0      1      39724
    1B_0      1      40086
    1C_0      1      40272
    1D_0      1      39932
    2E_0      1      40034
    2F_0      1      39866
    2G_0      1      40060
    2H_0      1      40398
    3I_0      1      39770
    3J_0      1      39644
    3K_0      1      39930
    3L_0      1      40016
    
    
    Full stats files
    ------------------------------------------------
    step 1: ./lane1raws_fastqs/s1_demultiplex_stats.txt
    step 2: None
    step 3: None
    step 4: None
    step 5: None
    step 6: None
    step 7: None
    
    


.. code:: python

    ## run remaining steps on the merged assembly
    ipyrad -p params-both.txt -s 234567

3. Multiple Libraries Multiple Lanes of Sequencing
--------------------------------------------------

The recommended way to combine multiple lanes of data is the same as we
just demonstrated above, however, in this case because the Samples in
each Object come from a different library, they will have different
names. Imagine that each lane of sequencing contains a library with 48
Samples in it. In the example above (One library multiple lanes) the
Samples would be combined so that you have 48 Samples, and each Sample
has data from two fastq files. Alternatively, the merging in this
example would combine the two libraries that contain different Samples
into a single data set with 96 Samples, where each Sample has one lane
of data.

4. Separate Multiple Libraries from One Lane of Sequencing
----------------------------------------------------------

.. code:: python

    ## create new Assembly named lib1
    ipyrad -n lib1 
    
    ## enter raw_fastq_path and barcodes_path into params
    sed -i '/\[2] /c\ipsimdata/rad_example_R1_.fastq.gz  ## [2] ' params-lib1.txt
    sed -i '/\[3] /c\ipsimdata/rad_example_barcodes.txt  ## [3] ' params-lib1.txt
    
    ## demultiplex the lane of data
    ipyrad -p params-lib1.txt -s 1 
    
    ## create a new branch with only the Samples for project 1
    ipyrad -p params-lib1.txt -b project1 1A_0 1B_0 1C_0 1D_0 
    
    ## create a another branch with only the Samples for project 2
    ipyrad -p params-lib1.txt -b project2 2E_0 2F_0 2G_0 2H_0 

.. code:: python

    ## assemble project 1 
    ipyrad -p params-project1 -s 234567

.. code:: python

    ## assemble project 2
    ipyrad -p params-project2 -s 234567

5. Alternative: Using the *ipyrad* API to do these things
---------------------------------------------------------

Using the *ipyrad* API is an alternative to using the
command-line-interface (CLI) above. As you can see below, writing code
with the Python API can be much simpler and more elegant. We recommend
using the API inside a Jupyter-notebook.

.. code:: python

    ## import ipyrad
    import ipyrad as ip

.. code:: python

    ## one lane one library
    data1 = ip.Assembly("data1")
    data1.set_params("raw_fastq_path", "ipsimdata/rad_example_R1_.fastq.gz")
    data1.set_params("barcodes_path", "ipsimdata/rad_example_barcodes.txt")
    data.run("123467")

.. code:: python

    ## one library multiple lanes
    lib1lane1 = ip.Assembly("lib1lane1")
    lib1lane1.set_params("raw_fastq_path", "ipsimdata/rad_example_R1_.fastq.gz")
    lib1lane1.set_params("barcodes_path", "ipsimdata/rad_example_barcodes.txt")
    lib1lane1.run("1")
    
    lib1lane2 = ip.Assembly("lib1lane2")
    lib1lane2.set_params("raw_fastq_path", "ipsimdata/rad_example_R1_.fastq.gz")
    lib1lane2.set_params("barcodes_path", "ipsimdata/rad_example_barcodes.txt")
    lib1lane2.run("1")
    
    merged = ip.merge("lib1-2lanes", [lib1lane1, lib1lane2])
    merged.run("234567")

.. code:: python

    ## multiple libraries multiple lanes
    lib1lane1 = ip.Assembly("lib1lane1")
    lib1lane1.set_params("raw_fastq_path", "ipsimdata/lib1_lane1_R1_.fastq.gz")
    lib1lane1.set_params("barcodes_path", "ipsimdata/lib1_barcodes.txt")
    lib1lane1.run("1")
    
    lib1lane2 = ip.Assembly("lib1lane2")
    lib1lane2.set_params("raw_fastq_path", "ipsimdata/lib1_lane2.fastq.gz")
    lib1lane2.set_params("barcodes_path", "ipsimdata/lib1_barcodes.txt")
    lib1lane2.run("1")
    
    lib2lane1 = ip.Assembly("lib1lane1")
    lib2lane1.set_params("raw_fastq_path", "ipsimdata/lib2_lane1.fastq.gz")
    lib2lane1.set_params("barcodes_path", "ipsimdata/lib2_barcodes.txt")
    lib2lane1.run("1")
    
    lib2lane2 = ip.Assembly("lib1lane2")
    lib2lane2.set_params("raw_fastq_path", "ipsimdata/lib2_lane2_.fastq.gz")
    lib2lane2.set_params("barcodes_path", "ipsimdata/lib2_barcodes.txt")
    lib2lane2.run("1")
    
    fulldata = ip.merge("fulldata", [lib1lane1, lib1lane2, lib2lane1, lib2lane2])
    fulldata.run("234567")

.. code:: python

    ## splitting a library into different project
    project1 = ["sample1", "sample2", "sample3"]
    project2 = ["sample4", "sample5", "sample6"]
    
    proj1 = fulldata.branch("proj1", subsamples=project1)
    proj2 = fulldata.branch("proj2", subsamples=project2)
    
    proj1.run("234567", force=True)
    proj2.run("234567", force=True)

.. code:: python

    ## print stats of project 1
    print proj1.stats


For advanced examples see the CookBook section. 
