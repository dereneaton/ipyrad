.. include:: global.rst  

.. _plotting:


Plotting
========
...


Depth plots
^^^^^^^^^^^
Sequencing depth varies greatly among RAD-seq data sets depending on the library
construction method and the sequencing effort, and it can vary across samples
as well. Any Samples that have completed step 3 of an ipyrad assembly will have
depth information available which we can access through the API to analyze and 
visualize. 

.. code-block:: python

    ## import ipyrad
    import ipyrad as ip

    ## load the assembly that is passed step 3
    data1 = ip.load.load_assembly("tests/data1")

    ## this example has 12 Samples
    data1.samples

    ## let's look at the depth data for Sample '1A_0'
    data1.samples["1A_0"].depths

    ## the data are stored as a numpy array, so there are many
    ## different operations you can perform to analyze it. 
    data1.samples["1A_0"].depths.mean()

    ## using these data you can create plot on your own using 
    ## various plotting libraries. We provide some basic plotting
    ## functions in the ipyrad.plotting library
    import ipyrad.plotting as ipp

    ## create a plot, see cookbook for further options to depthplot
    canvas = ipp.depthplot(data1, ...)

    ## save the plot as a pdf
    canvas.render_pdf("depthfig.pdf")



