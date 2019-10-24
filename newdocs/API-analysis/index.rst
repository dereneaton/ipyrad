.. toytree documentation master file, created by
   sphinx-quickstart on Tue Oct 16 20:42:08 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==========================
API: ipyrad analysis tools
==========================

The ipyrad-analysis toolkit is a Python interface for taking the output
files produced in a ipyrad assembly and running a suite of evolutionary
analysis tools with convenient features for filtering for missing data, 
grouping individuals into populations, dropping samples, and more. 

All of these tools share a common syntax making them easy to use without 
having to worry about creating different input files, or learn new 
file formats. They are designed for use within Jupyter notebooks, a tool 
for reproducible science. See the examples below. 


.. code-block:: python

   # the analysis tools are a subpackage of ipyrad
   import ipyrad as ipa

   # a large suite of tools are available 
   tool = ipa.structure(data="./outfiles/data.snps.hdf5")

   # all tools share a common syntax for setting params
   # and distributing work in parallel.
   tool.run()


.. toctree::
   :titlesonly:

   ipa.vcf_to_hdf5 <cookbook-vcf2hdf5.ipynb>
   ipa.treemix <cookbook-treemix.ipynb>
   ipa.pca <cookbook-pca.ipynb>
   ipa.raxml <cookbook-raxml.ipynb>
   ipa.tetrad <cookbook-tetrad.ipynb>
   ipa.structure <cookbook-structure.ipynb>
   ipa.sratools <cookbook-sratools.ipynb>
   ipa.baba <cookbook-abba-baba.ipynb>
   ipa.bucky <cookbook-bucky.ipynb>
   ipa.bpp <cookbook-bpp.ipynb>
   ipa.window_extacter <cookbook-window_extracter.ipynb>
   ipa.treeslider 1 <cookbook-treeslider.ipynb>
   ipa.treeslider 2 <cookbook-treeslider-downstream.ipynb>
   ipa.distance <cookbook-genetic-distance.ipynb>
   ipa.mrbayes <cookbook-mrbayes.ipynb>
