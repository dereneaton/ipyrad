
.. _assembly_guidelines:  

Guidelines for RADSeq Assemblies
================================

What parameters have the biggest effect on accuracy of genotypes? What parameter testing should be done with a new taxon or dataset?
------------------------------------------------------------------------------------------------------------------------------------
Step 1: `trim_reads` - This isn't a parameter *per se* but it is a parameter
that can be used in conjunction with `fastqc` in order to remove regions with
excessive error from the raw data. Distal ends of reads (especially paired-end
data) can start to accummulate lots of low quality bases. This will increase
the error rate of basecalling, decrease the performance, and increase the
false positive rate. It's **always** good to look at your data with `fastqc`
before you start an assembly.

Step 3: `clust_threshold` - This is probably the single most important parameter
in the whole process, controlling the sequence similarity threshold above which
two reads/loci are considered to be from the same genomic position. Values that
are too low will cause *over-lumping* (non-homologous regions clustering) which
will inflate heterozygosity and error rate. Values that are too high will cause
*over-splitting* (homologous regions failing to cluster) which will depress
heterozygosity and error rate. The optimization procedure of Mastretta-Yanes et al
2015 is straightforward to adapt from stacks to ipyrad and I've seen people use
this to good effect. One can also evaluate the effect of this parameter on the
`hetero_est` and `error_est` from the step 4 stats file. Error rate should be on
the order of 0.001 and heterozygosity should be on the order of 0.01. Excessively
high or low values indicate `clust_threshold` misspecification.

What are the relevant differences between the RADseq pipelines for someone choosing to do their first RADseq analysis?
----------------------------------------------------------------------------------------------------------------------
Key strengths of `ipyrad`:

* Simple, easy to use CLI interface
* Numerous statistics and reporting on assembly quality at each step of the process
* Powerful API mode for jupyter notebook assemblies and reproducibility
* Massively Parallel (w/ MPI across HPC compute nodes)
* Built-in analysis tools for running and plotting downstream popgen and phylogenetic analyses (PCA, STRUCTURE, RAxML, etc)
* Documentation: https://ipyrad.readthedocs.io and extensive tutorials: https://radcamp.github.io/
* Support: https://gitter.im/dereneaton/ipyrad

How can choices made during genotyping or filtering affect downstream analyses and summary statistics (e.g., Fst, demographic analyses, phylogenetics)?
-------------------------------------------------------------------------------------------------------------------------------------------------------
Misspecified `clust_threshold`: This will have obvious consequences for
over/under-estimating popgen sumstats (Fst, pi, etc). Similarly with demographic
analysis, if heterozygosity is inflated/deflated because of misspecified clust_threshold
then estimates of Ne, divergence time, and migration rates will all be skewed.

Missing data and `min_samples_locus`: This gets tricky because you can wind up
in hot water one way or the other depending on how permissive/conservative the
setting of this parameter is. If `min_samples_locus` is very high, then it will
remove many/most loci with missing data. This will do two things: 1) it will bias
against low-frequency variants (see Huang & Knowles 2016) and 2) it will also
bias toward more conserved regions which will have consequences downstream. If
`min_samples_locus` is low (*recommended*) then the data matrix will contain
some amount of missing data which is okay, but which will need to be properly
handled downstream (see below).

Some guidelines to follow when filtering your data.
---------------------------------------------------
The greatest crime to commit against RADSeq data is to over-filter on
missingness. Do not do this! The `min_samples_locus` parameter is the most
often abused parameter of all. Do not treat RADSeq data as if it's a massive
multi-locus dataset. Do not expect a complete data matrix. Embrace the
uncertainty and propagate it downstream.

When filtering data think in terms of **orders of magnitude**. Look at this
example from real data, where is the vast majority of filtering happening:

.. code-block:: bash

                                total_filters  applied_order  retained_loci
    total_prefiltered_loci              87851              0          87851
    filtered_by_rm_duplicates            2045           2045          85806
    filtered_by_max_indels               1170           1170          84636
    filtered_by_max_snps                  888            135          84501
    filtered_by_max_shared_het           2775           2427          82074
    filtered_by_min_sample              39613          38851          43223
    filtered_by_max_alleles             10210           3793          39430
    total_filtered_loci                 39430              0          39430

The **right** way to handle RADSeq data is to set a very permissive `min_samples_locus`
and then **deal** with the missingness downstream. A good example of how to deal with
missing data in a principled way is illustrated in the ipyrad analysis PCA tutorial:
https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-pca.html#Subsampling-with-replication.
Similar procedures can be applied to other analysis methods.

* Other useful guidelines in the ipyrad faq: https://ipyrad.readthedocs.io/en/latest/faq.html

What are some trade-offs between reference-based and de novo assembly, for example when the reference is of poor quality?
-------------------------------------------------------------------------------------------------------------------------
Reference assembly:
* Faster
* Provides genomic coordinates which can be useful (annotations)
* Perhaps more effective paralog filtering (?)
* Advanced end-users only (requires special care and attention to details)

de novo assembly:
* Much more reliably produces "reasonable" looking assemblies

Other useful ideas
------------------
* As long as `clust_threshold` is in the ballpark it is okay.
* Allele dropout and MAF filtering don't matter that much.

Step 1: `max_barcode_mismatch` - In practice you shouldn't need to allow any barcode
mismatches during step 1 demultiplexing, the amount of reads recovered by allowing
1 or more mismatches tends to be a fraction of a percent, so not worth the effort. This
parameter is only really useful for triaging very low quality data where each read
recovered really counts.

Step 2: `filter_adapters` - This should **always** be enabled. There's very little
performance penalty for filtering adapters, and the downstream consequences of not
filtering can be severe.

Step 3: `reference_as_filter` - ipyrad allows to provide a reference sequence to
**remove** mapped sequences prior to step 3 clustering. This can be **extremely**
useful for removing contaminants (e.g. microbial contaminants or mtDNA/cpDNA bycatch).

Step 7: When combining multiple plates or data from multiple lanes/runs/experiments
it's very important to evaluate the results for **batch effects.** For example, if
one combines samples from libraries that use two different restriction enzymes or
two different size selection windows then the samples will share much more data *within*
libraries than *among* libraries, regardless of the true relatedness of the samples.
ipyrad implements a `sharing` analysis tool to evaluate the extent of locus sharing
across samples. An example notebook is here:
https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-sharing.html

Step 7: `pop_assign_file` - This parameter can be used to specify a file which maps
samples to *a priori* designated populations which can then be used for setting
`min_samples` values per population. This can be useful if, for example, you have
suspected populations or sampling sites and you want to retain a minimum number of
samples per pop/site for downstream analysis.
