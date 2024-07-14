

====================================
Analysis: Tutorials
====================================

This page contains links to all ipyrad-analysis tutorials demonstrated
with both empirical and simulated examples. Please take note of the 
appropriate citations for any tools used.



ipa.vcf_to_hdf5 (format/filter)
-------------------------------
**vcf_to_hdf5 is used to convert VCF formatted SNP data to HDF5 format**.
This allows you to convert data produced by other bioinformatic tools (e.g., 
stacks, GATK) for use in ipyrad-analysis. The HDF5 file format 
is useful because it is super fast and efficient for loading data, 
which allows us to load it, filter it, and convert it to the appropriate 
format for each analysis tool, so you can easily apply different filters 
for different tools and never worry about formatting data files, 
or overloading your RAM.

:Empirical tutorial: `vcf_to_hdf5 conversion on RAD and WGS datasets <cookbook-vcf_to_hdf5-empirical.html>`_
:Simulated tutorial: `vcf_to_hdf5 conversion on simulated VCF with missing data <cookbook-vcf_to_hdf5-empirical.html>`_
:Citation: ipyrad-analysis



.. Convert sequences to HDF5
.. -------------------------
.. If you assembled your data using ipyrad then you will not need this tool 
.. since HDF5 files are automatically produced as outputs. This tool allows you 
.. to convert a sequence alignment from NEXUS format (e.g., produced by some 
.. external tool) to the HDF5 database format (seqs.hdf5) used by ipyrad-analysis. 
.. The seqs.hdf5 file contains the sequence data and information about which 
.. sequences are from which RAD loci, and in the case of reference-mapped data, 
.. the genomic positions of the loci. Several ipyrad-analysis tools can use this 
.. information to combine information (concatenate) nearby loci, such 
.. as window-extracter and treeslider. The NEXUS file must contain 
.. a partition block to assign sequences to loci.
.. + `ipa.nex2hdf5 tutorial on empirical xyxy data <...>`_



ipa.window_extracter (format/filter)
------------------------------------
**window_extracter** is used to filter and subsample sequence alignments**.
You can apply filters for missing data. 
Write filtered alignment to phylip or nexus format for downstream analyses.

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-structure-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-structure-ipcoal.ipynb>`_. 
:Citation: `ipyrad-analysis <...>`_





ipa.pca (pop structure)
------------------------------------
**This tool is used to implement PCA and other dimensionality reduction methods
to infer and visualize population structure**. You can apply filters or imputation
to reduce the impact of missing data and we provide built-in plotting tools.
In addiction to principal components analysis (PCA) you can implement 
t-SNE or UMAP.

:Empirical tutorial: `PCA/tSNE/UMAP in 7 species (35 samples) oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-pca-empirical.ipynb>`_
:Simulated tutorial: `Setup and run PCA on your own simulated scenario) <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-pca-ipcoal.ipynb>`_
:Citation to external software: `ipyrad-analysis`, `scikit-learn`




ipa.structure (pop structure)
------------------------------
**STRUCTURE is used to infer population structure.** The ipyrad-analysis tool 
acts as a simple wrapper for calling raxml from a jupyter notebook
to perform maximum likelihood tree inference.

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-structure-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-structure-ipcoal.ipynb>`_. 
:Citation of external software: `Pritchard <...>`_



ipa.baba (admixture inference)
------------------------------
This tutorial shows how to use the ipa.baba tool to analyze SNP data to test
hypotheses of admixture on a species tree hypothesis. You can implement 
4-taxon (abba-baba or D-statistic) tests, or 5-taxon 
(partitioned D-statistics) tests. Results are returned in a pandas DataFrame, 
which can be save to CSV, and plotting tool for summarizing results.

:Empirical tutorial: `Example analysis on 7 species of oaks <...>`_
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/newdocs/API-analysis/cookbook-baba-ipcoal.ipynb>`_.  
:Citations of methods: `Durand <...>`_, `Eaton <...>`_


ipa.treemix (admixture inference)
---------------------------------
Treemix infers a distance-based admixture graph from F-statistics inferred
from SNP frequency data. This method assumes that all variants were present
in a common ancestor and sorted among descendant populations by drift or gene
flow, and so it is not suitable for deeper-scale phylogenetic relationships.
However, it can be a useful exploratory analysis tool. The ipyrad-analysis tool
filters data, groups individuals into populations to format input files, and 
has plotting tools for visualizing results.

:Empirical tutorial: `Example analysis on 7 species of oaks <...>`_
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/newdocs/API-analysis/cookbook-treemix.ipynb>`_. 
:Citation of external software: `Pickerell <...>`_


ipa.raxml (gene trees)
----------------------
**RAxML is a tool for phylogenetic inference.** The ipyrad-analysis tool 
acts as a simple wrapper for calling raxml from a jupyter notebook
to perform maximum likelihood tree inference. This reduces clutter when 
calling RAxML many hundreds of times using the same options. We
recommend to use ipa.window_extracter to optimize alignments (.phy) 
for missing data prior to raxml inference.

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-raxml-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-raxml-ipcoal.ipynb>`_. 
:Citation of external software: `Stamatakis <...>`_



ipa.treeslider (gene trees)
---------------------------
**Treeslider is used to perform phylogenetic inference across whole genomes**.
It can apply to genomic windows containing RAD loci mapped to a reference genome, 
or to all individual anonymous loci in a denovo assembly. It allows you to 
automate the application of `window_extracter` to apply common filtering options
to all windows, followed by an inference tool with its options (e.g., raxml, 
mrbayes). 

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-treeslider-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-treeslider-ipcoal.ipynb>`_. 
:Citation of external software: `ipyrad-analysis <...>`_



ipa.clade_weights (gene trees)
------------------------------
**clade_weights** is an ipa tool used to summarize and visualize spatial 
patterns of gene trees along chromosomes. To use this method your data must
be aligned to a relatively well-assembled reference genome. The tool is
modeled on the TWISST method, which should be cited. 

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-treeslider-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-treeslider-ipcoal.ipynb>`_. 
:Citation of external software: `ipyrad-analysis <...>`_, TWISST



ipa.mrbayes (gene trees)
---------------------------
**Mrbayes is a Bayesian phylogenetic inference tool**.
most commonly used to inferring dated phylogenies**.
...

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-mb-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-mb-ipcoal.ipynb>`_. 
:Citation of external software: `ipyrad-analysis <...>`_




ipa.bpp (species trees)
---------------------------
**BP&P is a Bayesian phylogenetic inference tool for species trees**.
...

:Empirical tutorial: `Example 00 analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-bpp-00-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-bpp-00-ipcoal.ipynb>`_. 

:Empirical tutorial: `Example 01 analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-bpp-00-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-bpp-00-ipcoal.ipynb>`_. 

:Citation of external software: `ipyrad-analysis <...>`_





ipa.astral (species tree)
-------------------------
**Astral is used to infer a species tree from a collection of gene trees**. 
The recommended workflow in ipyrad-analysis is to first run ipa.treeslider to 
get a collection of hundreds or thousands of trees and then to input them to 
ipa.astral. This is demonstrated in the tutorials.
...

:old tutorial: `ASTRAL species tree inference <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/newdocs/API-analysis/cookbook-astral.ipynb>`_
:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-astral-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-astral-ipcoal.ipynb>`_. 
:Citation of external software: `ipyrad-analysis <...>`_



ipa.tetrad (species tree)
--------------------------
**Tetrad is a tool for inferring a species tree from SNP data**. It follows the
algorithm implemented in SVDQuartets to infer quartet trees from phylogenetic
invariant patterns in SNP data. You can also run this tool externally as a command
line tool (see `https://github.com/eaton-lab/tetrad <https://github.com/eaton-lab/tetrad>`_), 
the ipyrad-analysis wrapper is simply a convenience for running this tool in a 
jupyter notebook.

:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-tetrad-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-tetrad-ipcoal.ipynb>`_. 
:Citation of external software: `Eaton 2015`, `Kubatko`.





ipa.snaq (network inference)
----------------------------
...
:old tutorial: `SNAQ phylogenetic network inference <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/newdocs/API-analysis/cookbook-snaq.ipynb>`_
:Empirical tutorial: `Example analysis on 7 species of oaks <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-snaq-empirical.ipynb>`_. 
:Simulated tutorial: `Setup your own tests on a simulated scenario <https://nbviewer.jupyter.org/github/dereneaton/ipyrad/blob/master/testdocs/analysis/cookbook-snaq-ipcoal.ipynb>`_. 
:Citation of external software: `ipyrad-analysis <...>`_







ipa.digest_genomes (assembly)
------------------------------
**This tool is used to perform *in silico* digestion of a fasta genome file
to produce RADseq-like fastq Illumina reads**. You can use this to extract
RAD data from published genomes that you can then include in your assembly.
This is particularly useful for including additional outgroup samples.









..  toctree
..    :maxdepth:1

..    ipa.vcf_to_hdf5 <cookbook-vcf2hdf5.ipynb>
..    ipa.treemix <cookbook-treemix.ipynb>
..    ipa.pca <cookbook-pca.ipynb>
..    ipa.raxml <cookbook-raxml.ipynb>
..    ipa.mrbayes <cookbook-mrbayes.ipynb>
..    ipa.tetrad <cookbook-tetrad.ipynb>
..    ipa.structure <cookbook-structure.ipynb>
..    ipa.sratools <cookbook-sratools.ipynb>
..    ipa.baba <cookbook-abba-baba.ipynb>
..    ipa.bucky <cookbook-bucky.ipynb>
..    ipa.window_extracter <cookbook-window_extracter.ipynb>
..    ipa.digest_genomes <cookbook-digest_genomes.ipynb>
..    ipa.treeslider <cookbook-treeslider.ipynb>

.. 
   ipa.bpp <cookbook-bpp.ipynb>
   ipa.treeslider 1 <cookbook-treeslider.ipynb>
   ipa.treeslider 2 <cookbook-treeslider-downstream.ipynb>
   ipa.distance <cookbook-genetic-distance.ipynb>
