
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

TODO...

For advanced examples see the CookBook section. 