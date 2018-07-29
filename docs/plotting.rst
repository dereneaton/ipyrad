.. include:: global.rst  

.. _plotting:


Plotting
========
The ipyrad_ plotting library has a number of basic function written with the 
Python module Toyplot_. For each of these I will demonstrate a simple usage 
of the available function, as well as how to access the raw data in case you
want to create your own plotting functions using Python or R. An important
note about using the ipyrad_ plotting functions is that they have to be loaded
from a separate module, like below. 

.. code-block:: python  

    ## import ipyrad and the ipyrad plotting library
    import ipyrad as ip
    import ipyrad.plotting as ipp


Depth plots
^^^^^^^^^^^
Sequencing depth varies greatly among RAD-seq data sets depending on the library
construction method and the sequencing effort, and it can vary across samples
as well. Any Samples that have completed step 3 of an ipyrad assembly will have
depth information available which we can access through the API to analyze and 
visualize. 

.. meth:: depthplot

    describe the positional arguments to depthplot here...
    

.. code-block:: python

    ## load the assembly that is passed step 3
    data1 = ip.load_json("tests/data1.json")

    ## this example has 12 Samples
    data1.samples

    ## let's look at the depth data for Sample '1A_0'
    data1.samples["1A_0"].depths

    ## the data are stored as a numpy array, so there are many
    ## different operations you can perform to analyze it. 
    data1.samples["1A_0"].depths.mean()

    ## using these data you can a create plots using any plotting library.
    ## We provide a simple create a plot, see cookbook for further options to depthplot
    canvas = ipp.depthplot(data1, ...)

    ## save the plot as a pdf
    canvas.render_pdf("depthfig.pdf")



