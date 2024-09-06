
.. _release_notes:

Release Notes
=============

0.9.102
-------
**September 06, 2024**

- fix SE loading to allow mixed .fq and .fastq file names

0.9.101
-------
**September 03, 2024**

- Fix bug in SE sorted_fastqs

0.9.100
-------
**August 28, 2024**

- Use the PE/SE information in the `datatype` param to determine whether to search for paired files.
- hackers mask cut sites in s3
- hackers mask cut sites in s3
- new hackers option 'mask_restriction_sites' for s7

0.9.99
------
**August 16, 2024**

- Step 7 uses uint32 for edges for ref assemblies #571

0.9.98
------
**August 15, 2024**

0.9.97
------
**August 12, 2024**

- fix for pd whitespace warning
- new name parser
- Documenting a bit better the tweaks in pair_fastqs to allow . in sample name
- pair_fastqs.py - remove some of the shuffling with suffixes
- allow . in sample names for step 1 loading from fastqs #557
- verify even number of files for PE data when linking fastqs
- fix for pd warnings
- update for toy3 compat and pd warning
- fix pandas warnings
- Update faq.rst - Add link to Hemstrom et al 2024

0.9.96
------
**June 18, 2024**

- Fix consens_se filter_mindepth to better handle haploid data (removes spurious het calls in haploid)
- ipa.treemix test for toytree < v3
- Removed very old comment about alleles format in the faq

0.9.94
------
**February 02, 2024**

- Fix pandas delim_whitespace deprecation
- fix popgen merge
- ipa.popgen minor edits
- update assembly method docs to remove deprecated methods

0.9.93
------
**August 16, 2023**

- add mac arm install instructions
- Allow popgen.py to handle missing data a bit better"
- Update faq.rst
- test muscle binary version for #504
- analysis.sratools allow sample names to be all integers

0.9.92
------
**May 31, 2023**

- analysis.popgen: Allow missing data without crashing
- Fixed the popgen sumstats analysis tool and the cookbook

0.9.91
------
**May 27, 2023**

- Fix numpy deprecation of np.int and np.bool
- Update environment.yml
- Update README.rst
- Update README.rst
- README - Add binder badge
- Allow mrbayes.py to pass in the binary path

0.9.90
------
**April 03, 2023**

- Fixing a bunch of AttributeErrors on h5py attribute decode attempts
- Pin muscle < v5.0 in meta.yaml

0.9.89
------
**March 30, 2023**

- Fix for #502
- Update install instructions  to point to libmamba-solver

0.9.88
------
**March 12, 2023**

- Update faq.rst
- Update faq.rst
- Solving retro-compatibility in windows (#496)

0.9.87
------
**November 28, 2022**

- Fix oops in clustmap.py #494

0.9.86
------
**November 18, 2022**

- step3: Filter by min MAPQ in samtools view #494
- Allow project_dir directories to create recursively
- mod rtfd reqs.txt to unpin jinja2 version
- Add rtfd yml file

0.9.85
------
**September 15, 2022**

- Fix #493 - ValueError in step 6 if all R2s empty
- Update 3-installation.rst
- Update README.rst
- Update faq.rst
- Handle call to decode() in case older version of h5py file.
- Clean up tempdir after step 5
- Update faq.rst

0.9.84
------
**February 05, 2022**

- fix pd 1.4.0 issue with setting options #475

0.9.83
------
**February 01, 2022**

- ipyrad.core.sample Set dtype of stats pd.Series as object to silence a pd warning.
- Test vsearch version to handle derep_fulllength deprecation #469
- structure error reporting (#471)

0.9.82
------
**January 21, 2022**

- pin vsearch version <=2.19 #469
- analysis.snps_extracter: allow passing in vcf file (assuming RADSeq data).
- Add check for executable of structure binary
- kmeans imputation niters fix (#459)
- cookbook-bpp: Update tauprior to bpp 4.x inverse gamma standard
- cookbook-bpp remove randomize_order flag from run() as this doesn't exist anymore

0.9.81
------
**June 16, 2021**

- Step5 add filter_max_alleles() function to honor the max_alleles_consens parameter.

0.9.80
------
**June 09, 2021**

- samtools defaults to using the  flag to generate .csi index files, to allow for chrom size > 500Mb

0.9.79
------
**June 09, 2021**

- Step 2 - Do NOT take the complement for the p5 (R2) adapter sequence, otherwise the adapter trimming doesn't work.

0.9.78
------
**May 28, 2021**

- #444: fix bytes to string comp oops

0.9.77
------
**May 11, 2021**

0.9.76
------
**May 11, 2021**

0.9.75
------
**May 11, 2021**

0.9.74
------
**May 08, 2021**

0.9.73
------
**May 07, 2021**

0.9.72
------
**May 07, 2021**

- Docs: add some assembly guidelines
- Docs: update sharing and popgen sumstats
- Step1.FileLinker fix oops with trapping bz2 formatted input files.

0.9.71
------
**April 15, 2021**

- demultiplex: Properly handle 3rad barcodes of different length

0.9.70
------
**April 01, 2021**

- Allow snps_extractor to handle snps.hdf5 files with names not encoded as bytes
- Fixing mismatch of #SBATCH and command parameters (#440)

0.9.69
------
**March 19, 2021**

- ipa.structure: handle edge case with running more jobs with different K values

0.9.68
------
**February 23, 2021**

- Handle i7 demux to strip trailing newline from barcode
- demultiplex.py: Allow very short input fq files.
- Fix project_dir inconsistency in merged assemblies.
- Raise an error if setting a bad parameter in API mode #354

0.9.67
------
**February 21, 2021**

- Don't blank extra p3/p5 adapters in step2.check_adapters if filter_adapters != 3.

0.9.66
------
**February 18, 2021**

- analysis.popgen: _collect_results() function to compile and organize results from the engines. It's a little ugly.
- analysis.popgen Change default minmap value to 4
- analysis.popgen Chunking of loci and parallelization of analysis working.
- analysis.popgen._fst_full() implemented. Had to slightly tweak the code that was there, and it works a bit different than the other sumstats, but whatever for now it works.
- analysis.popgen implemented Dxy
- analysis.popgen._process_locus_pops() to prepare data for pairwise summary statistics
- analysis.popgen: Stats returned as dict instead of tuple, and remote processes dump full stats to pickle files.
- Add a _process_locus function to remove a bunch of redundancy, and split consens bases to account for diploidy in all the calcs
- analysis.popgen: If no imap then all samples are considered from the same population.
- analysis.popgen: refactor _Watterson to accept sequence length and return raw and per base theta values.
- analysis.popgen: Refactor into 'public' sumstat methods which can be called on array-like locus data, and 'semi-private' methods (e.g. _pi, or _Watterson) that are much more efficient but accept very specific information. The Processor.run() function will do some housekeeping per locus and call the 'semi-private' methods, for efficiency. Public methods are for testing.
- analysis.popgen Tajima's D function
- analysis.popgen Watterson's theta function
- clustmap: Handle bam files with very large chromosomes (#435)
- clustmap_across: Catch samtools indexing failure for large bam files and try again with .csi indexing with the  flag.
- Adopting the Processor() structure to enable parallelization. Also implemented nucleotide diversity.
- analysis.popgen add parallelization support for calculating all sumstats per locus.
- analysis.popgen _fst() is working for snp data
- Add an option to LocusExtracter.get_locus to return a pd.DataFrame indexed by sample names
- analysis.popgen getting it into current analysis tools format and starting to flesh it out.
- CLI honors the -q flag now
- analysis.sharing allow subsampling inside draw() to prevent recalculating data
- analysis.shared allow sorting by mean pairwise sharing or missing
- Allow analysis.sharing to reorder samples, and add a progress bar
- Add new pairwise snp/missingness sharing analysis tool.
- add locus sharing/missingness sharing analysis tool skel
- pin numpy >=1.15, a known working version to address #429
- Fix weird race condition with branching and pop_assign_files. See #430.

0.9.65
------
**January 22, 2021**

- Fix #429 weird masking bug in older versions of numpy.
- Add refs to the analysis.popgen tool.

0.9.64
------
**January 16, 2021**

- replaced core.Assembly.database which actually wasn't doing anything with snps_database and seqs_database to retain info about the hdf5 files in the assembly object
- fix empirical api structure params format
- Allow structure to accept vcf files and auto-convert to hdf5
- fix oops in i7 demux cookbook

0.9.63
------
**December 17, 2020**

- Fix off-by-one error in nexus output
- update struct testdocs
- Add Tajima's D denominator equation to the popgen analysis tool, because I coded it before and it's a nightmare.
- use quiet in lex
- plot posteriors range limits removed
- Actually fix pca cookbook
- fix malformatted pca cookbook
- raxml w/ gamma rates

0.9.62
------
**November 03, 2020**

- updating mb cookbook
- bpp transform parser bugfix
- add bpp dev nb
- allow non-blocking calls again in bpp
- bpp result split on delimiter
- allow empty node dist in bpp
- add option to extract as string/bytes
- snps_extracter: report maf filter as unique filtered
- snps extracter reports n unlinked snps
- tmp dir for new cookbook docs
- allow maf filter in structure
- node_dists check for tips in bpp plot
- new plotting tools for bpp results
- added download util
- add minmaf filter to snps_extracter, pca
- additionall bpp summary tools

0.9.61
------
**October 15, 2020**

- Fix nasty nasty bug in step 7 for issue #422

0.9.60
------
**October 13, 2020**

- Fix mapping minus reads bug step 3

0.9.59
------
**September 20, 2020**

- In structure._call_structure() self is not in scope, so now we pass in the full path to the structure binary.

0.9.58
------
**September 10, 2020**

- fix oops in window_extracter
- Allow scaffold idxs to be int
- bugfix to skip scaffs with no hits even when end=None
- end arg offset bugfix, may have affected size of windows in treeslider
- Add handler for malformed loci in baba.py
- fix to allow digest of single end RAD
- Changed astral annotation to default (#417)
- Update dependency documentation
- baba2 draw fix
- return canvas in draw
- Fix an oops iff hackersonly.declone_PCR_duplicates && reference assembly
- Merge pull request #409 from camayal/master
- path fix for conda envs
- path fix for conda envs
- update for py38
- Update faq.rst
- pca docs
- Change exit code for successful merging from 1 to 0.
- Remove the png outfile from pca because it wants ghostscript package installed and apparently toytree doesn't have that as a dependency. Annoying.
- Update README.rst
- Allow pca to write png as well.
- Added fade option for blocks and tooltips
- Merge pull request #1 from camayal/camayal-baba2-work
- Changes in Drawing class

0.9.57
------
**August 01, 2020**

- Add functionality to pca to allow adding easy text titles to the plots with the  param
- Document writing pca figure to a file
- Allow pca.draw() to write out pdf and svg
- force option to remove cftable in snaq
- fix ntaxa in phy header for lex
- mb ipcoal doc up
- get locus phy uppers
- changed default nboots in baba
- Set pandas max_colwidth=250 to allow for very long sample names.
- Fix oops in step 7 where trim_loci 3' R1 wasn't being used.
- Allow imap pops to be specified as np.ndarray in snps_extracter
- fix path snaq log
- Fix branching oops
- snaq working
- snaq
- network updated
- testing network analysis
- ast conda bin path fix
- ast conda bin path fix
- ast conda bin

0.9.56
------
**June 29, 2020**

- revise installation docs to add conda-forge recommendation to reduce conflict errors
- update README install instructions

0.9.55
------
**June 26, 2020**

- bpp with ipcoal demoup
- update tmx tool for toytree v2
- conversion funcs in snps_extracter
- binder add ipcoal and CF as first env
- ipcoal terremix docs
- tool and docs for genotype freq conversion output
- tool and docs for genotype freq conversion output
- baba ipcoal notebook
- update 2brad docs slightly
- Add force flag to force overwrite vcf_to_hdf5 (prevent redundant conversion)
- better prior plots and transformers
- mb more prior options and fixed tree support
- pca allow setting colors
- pca cna set opacity
- pca cna set opacity
- pca cna set opacity
- pca cna set opacity
- fix to allow custom colors in pca
- docs
- fix py2 compat
- wex ints
- fix api oops in window extracter docs for scaffold_idxs
- wex start end as ints
- re-supporting single or multiple locs/chrom in wex
- option to not subsample SNPs
- add mrbayes keepdir arg to organize files
- simplified wex
- extracter filters invariant sites when subsampled by IMAP
- added pca panel plot
- subsample loci jit for snps_extractor of linked SNPS
- major baba2 update, to replace baba eventually
- axes label fix and figt cleanup
- handle int chrom names

0.9.54
------
**May 31, 2020**

- off-by-one to ref pos in s3 applied again here

0.9.53
------
**May 19, 2020**

- Fix off by 1 error in step 3 for PE data.
- Fix toytree documentation in baba cookbook
- Fix py2 compat by removing trailing commas from function argument lists in a couple of anaysis tools.
- Fix oops in handling errors during convert_outputs

0.9.52
------
**May 09, 2020**

- Fix nasty off-by-one error in reference positions
- multiple default clock models
- multiple default clock models
- multiple default clock models
- multiple default clock models
- multiple default clock models
- multiple default clock models
- multiple default clock models
- wex concat name drop fix
- ts tmpdir renamed bootsdir
- umap learn conda instructions
- ts: added dryrun method
- wex: remove print debug statement
- Fix baba cookbook docs
- add umap option
- added pseudocode for a further imputer in prog
- prettier bpp plot
- pca analysis tool passes through quiet flag to subfunctions
- warning about missing ref only up with no ref1 or ref2
- merge fix
- improving ipabpp summary funcs
- ensure conda ipcluster bin on stop
- bpp prior checks, new ctl build for 4.0, parsing results funcs
- Add a helpful message if merging assemblies with technical replicates beyond step 3.
- missing import
- Handle empty imap population in snps_extractor
- binary fix
- syntaxerr on quiet
- hide toytree dep in ast
- ast better error message
- ip assemble shows cluster on run by default
- show_cluster func now listens to param arg
- big update for bpp 4.0, uses lex
- wex and ts both use idxs in param name now
- simple astral run tool
- lex: imap/minmap filtering fix
- wex: imap/minmap filtering fix
- fixed warning message
- default minmap to 0 if imap and minmap empty
- hide toyplot dependency
- simple option to keep treefiles in treeslider
- under the hood mods to pca draw func to make it more atomic
- Update faq.rst
- Set default filter_adapters parameter to 2
- raise warning if ref+ or ref- and method not ref
- notes on window extracter

0.9.51
------
**April 17, 2020**

- 1 index POS in vcf output
- minmap default is 0
- bugfix: apply imapdrop only when imap
- faster extraction and mincov after minmap in lex
- mincov applies after minmap in wex
- scaff arg entered later in cov tool
- rmincov added to ts
- option to keep all tmp files in treeslider
- major fix to names sorting in wex
- names offset by scaff length in cov plot
- set default inner mate to 500 and use it unless user changes to None, in which case we estimate from reads
- tmp working baba update
- added locus extracter
- option to keep all files in treeslider
- added cov plot tool

0.9.50
------
**April 05, 2020**

- Actually fix FASTQ+64 problem. Max fastq_qmax is 126, so this is set to 93 now (93+33=126)

0.9.49
------
**April 02, 2020**

- Allow high fastq_qmax in pair merging to allow FASTQ+64 data

0.9.48
------
**April 01, 2020**

- Record refseq mapped/unmapped for both SE & PE
- wextract minmap+consred minmap default added
- treeslider default args typed
- tested working wextracter
- baba merge
- new dict for translation
- updating bpp for 4.0

0.9.47
------
**March 24, 2020**

- Fix snpstring length oops in .alleles outputs so they line up right.

0.9.46
------
**March 24, 2020**

- Fix pd.as_matrix() call which is deprecated.
- Force pca.draw() to honor the length of the color list if it is sufficiently long to color all samples in the imap, or at least use the length of the color list to set the value of the  variable.
- Fix oops in baba.py for importing msprime. Pushed it to Sim.__init__, since if you want to do baba, and don't care about sims, then you shouldn't have to install msprime.
- h5py warning fix
- use _pnames to use filtered names in run()

0.9.45
------
**March 08, 2020**

- Allow more flexibility in sorted fastqs directory (DO NOT DELETE if it points to projectdir + _fastqs)
- window extracter fix for multiple loci w/ reduce

0.9.44
------
**March 04, 2020**

- Fix the treemix output so it actually generates. Took WAYYYY longer than I thought it would.
- Update faq.rst
- Update faq.rst

0.9.43
------
**February 26, 2020**

- Fix off by 2 error in minsamp when using reference sequence
- window extacter working for denovo loci
- Cleaning up a TON of sphinx warnings from the docs and fixing a bunch of docs issues.
- fix oops in baba.py (import sys)

0.9.42
------
**February 19, 2020**

- Fix oops in step 6 which was leaving bad sample names hanging after alignment.

0.9.41
------
**February 18, 2020**

- Set s6.data.ncpus value when routing around hierarchical clustering for ref based assemblies.
- disable hierarchical clustering until further testing
- split samples evenly among cgroups for hierarch clust
- digest genomes uses qual score B instead of b

0.9.40
------
**February 16, 2020**

- subsample loci func added
- counts rm duplicates in denovo and works with step6 skipping alignment of loci w/ dups
- denovo paired aligned separately again
- fastq qmax error in merge denovo fixed

0.9.39
------
**February 15, 2020**

- Why can't i figure out how to comment out this plotting code right? wtf!

0.9.38
------
**February 15, 2020**

- commented out the import of the baba_plot plotting function and the baba.plot() method as these are broken rn, and also the plotting/baba_plotting routine tries to access toyplot in a way that breaks the conda build since toyplot isn't a strict requirement. We could fix this in the future, but i'm tring to get the bioconda package to build successfully rn.

0.9.37
------
**February 15, 2020**

- fix import checking for baba_panel_plot.py

0.9.36
------
**February 15, 2020**

- Handle external imports in the baba module in the same way as the other analysis tools to fix the broken bioconda build.
- Add a pops file to the ipsimdata.tar.gz because it's always useful.
-  "Updating ipyrad/__init__.py to version - 0.9.35

0.9.35
------
**February 12, 2020**

0.9.35
------
**February 12, 2020**

- Fix a bug in step 5 handling of RemoteError during indexing alleles.
- Report debug traceback for all crashes, not just API. This is essentially making the debug flag useless in v.0.9

0.9.34
------
**February 09, 2020**

- Roll back baba code to 0.7 version which doesn't use the current analysis format, but which still works. Saved ongoing baba code as baba.v0.9.py

0.9.33
------
**February 06, 2020**

- Fix major oops in consens_se which failed step 5 every time. Bad!
- In step 6 use the sample.files.consens info, rather than data.dirs to allow for merging assemblies after step 5 where data.dirs is invalid/empty.

0.9.32
------
**February 04, 2020**

- #392 allow scaffold names to be int
- Add sensible error handling if only a few samples fail step 5.
- add docs to clustmap_across
- fix for name re-ordering in window-extracter with multiple regions selected
- added comments
- added comments
- added sys
- Actually handle failed samples in step 2.
- fix for new h5py warning
- fix for new sklearn warning

0.9.31
------
**January 19, 2020**

- Fix error in bucky (progressbar hell).
- Add error handling in a couple cases if run() hasn't been called, e.g. before draw, and also add the pcs() function as a convenience.
- Removed support for legacy argument format from bpp.py and updated the docs.
- Allow PCA() to import data as vcf.
- Add support for importing VCF into PCA

0.9.30
------
**January 16, 2020**

- Fix whoops with bucky progressbars

0.9.29
------
**January 15, 2020**

- Fix bucky progressbar calls.
- Fixed progressbar calls in bucky.py
- Conda install instructions were wrong.
- add future print import to fasttree.py
- Add releasenotes to the TOC

0.9.28
------
**January 12, 2020**

- Fix versioner.py to actually record releasenotes
- Fix releasenotes in versioner script
- Fix releasenotes

0.9.27
------
**January 12, 2020**

- Fix releasenotes

0.9.26
------
**January 01, 2020**

- During steps 5 & 6 honor the filter_min_trim_len parameter, which is useful in some cases (2brad).
- In step 3, force vsearch to honor the filter_min_trim_len param, otherwise it defaults to --minseqlength 32, which can be undesirable in some cases.

0.9.25
------
**December 31, 2019**

- concatedits files now write to the tmpdir, rather than edits (#378), also handle refmap samples with no reads that map to the reference, also change where edits files are pulling from during PE merging to allow for assembly merging after step 2. phew.
- digested genome bugfix - check each fbit 0 and -1
- digest genomes nscaffolds arg support
- docs cookbook updates
- new consens sampling function and support for window extracter to concatenate
- comment about zlib

0.9.24
------
**December 24, 2019**

- Fix IPyradError import

0.9.23
------
**December 24, 2019**

- Regress baba.py

0.9.22
------
**December 23, 2019**

- Add support for .ugeno file
- Add support for .ustr format
- Remove duplication of code in write_str()
- Fix docs for output formats
- Add back output formats documentation

0.9.21
------
**December 23, 2019**

- Fix stupid bug introduced by fe8c2dfc282e177a7c18f6e2e23ef84d284a9e3f

0.9.20
------
**December 18, 2019**

- Expose analysis.baba for testing
- fasterq-dump seems to be only avail on linux
- Fix bug in handling sample names in the pops file. re: #375.
- Allow faidict scaffold names to be int (cast to dtype=object)

0.9.19
------
**December 03, 2019**

- Fix step 6 with pop_assign_file
- fix for empty samples after align
- list missing as ./. in VCF (like we used to)

0.9.18
------
**November 23, 2019**

- Fix oops handling missing data in vcf to hdf5
- mb binary path bugfix
- treeslider mb bugfix
- treeslider mb working
- treemix support for conda env installations
- additional drawing options for pca
- raxml cookbook update
- tetrad notebook updated
- Fix oops in params.py checking for lowercase overhangs seqs
- Fix a nasty stupid bug setting the overhang sequence
- Add back the docs about merging
- Error checking in step 5.
- Forbid lowercase in overhang sequence

0.9.17
------
**November 04, 2019**

- Ooops. Allow popsfile w/o crashing, and allow populations to be integer values
- cookbooks added link to nb
- pca stores results as attr instead or returning

0.9.16
------
**October 31, 2019**

- commented fix of optim chunksize calc
- treeslider now working with mb
- toggle to write in nexus
- mb saves convergence stats as df
- single-end mapping infiles bugfix
- pca cookbook update
- added fasttree tool
- update cookbooks index
- cookbooks updated headers
- warning about denovo-ref to use new param
- clustmap keeps i5s and can do ref minus
- window extracter updated
- mb load existing results and bugfix result paths
- update treemix and mb docs
- Fix calculation of optim during step 6 aligning. 3-4x speedup on this medium sized simulated data i'm working on.
- Fix oops in how optim was being counted. Was counting using _unsorted_ seed handle, I switched it to use sorted and now it works more like expected
- Clean up clust.txt files after step 3 finishes
- Update docs and parameter descriptions to reflect new reality of several params
- Add handlers for denovo +/- reference.
- vcf tool docs
- tools docs update
- enable vcf_to_hdf5
- pca reps legend looks nice
- added replicate clouds to pca
- vcf to hdf5 converter tested empirically
- Add hils.py from the hotfix branch
- Pull from the correct repo inside meta.yaml
- VCF 9's fixed to be .
- Add back tetrad docs
- default hackers set to  yes merge tech reps, and cleanup
- behavior for duplicates in barcodes file
- bugfix: error reporting for barcodes within n
- find binary from env or user entered
- find ipcluster from conda env bin
- bugfix: allow demux i7s even if datatype=pair3rad
- add the notebook tunnel docs back
- pedicularis cli tutorial updated
- df index needed sorting
- Allow sample names to be integers. wtf, how did this never come up before?
- Add the McCartney-Melstad ref to the faq
- fix typo in bpp cookbook
- Fix bpp Params import
- docs nav bar cleanup
- testing binder w/o treemix
- Add advanced tutorial back (?), maybe as a placeholder.
- Remove references to smalt and replace with bwa. That's some old-ass junk!
- Fixed the versioner.py script and added the faq.rst to the newdocs

0.9.14
------
**October 05, 2019**

- binder update
- docs update
- add bpp docs
- indentation in docs
- add i7 demux cookbook
- mroe analysis cookbooks
- analysis cookbooks
- merged clustmap
- Fix a nasty bug with stats for assemblies where chunks end up empty after filtering
- Fix step 3 to allow some tmpchunks to be empty without raising an error during chunk aligning
- Fix a bug in bpp.py
- Fix a nasty error in jointestimate.stackarray() where some long reads were slipping in over the maxlen length and causing a broadcast error

0.9.13
------

- py2 bug: print missing as float
- py2 bug fix: database ordering
- allow iterable params object for py2 and 3
- Fix an edge case to protect against empty chunks post-filtering during step 7
- install docs update
- Fix CLI so merge works
- Fix the max_shared_Hs param description to agree with only having one value, rather than 1 value perper R1/R2
- Ooops. checked in a pdb.set_trace in write outfiles. sorry\!
- add deps to newdocs
- bug fix for a rare trim that leaves >=1 all-N rows. Filter it.
- documenting a hard coded backward compatibility in write_output.Processor()
- hdf5 formatting for window slider in both denovo and ref
- sratools up to date with CLI working too
- Don't pester about mpi4py if you're not actually using MPI (CLI mode)
- Allow for user to not input overhang sequences and jointestimate will just proceed with the edges included.
- chunked downloads bug fix

**<sunspots cause discontinuity in version history>**

0.7.30
------
**March 09, 2019**

- Fix pca for scikit 1.2.0 API and a few minor fixes.
- Update faq.rst
- Update faq.rst

0.7.29
------
**January 21, 2019**

- Fix nasty ValueError bug in step 7 (re: merged PE loci)
- Update faq.rst
- Update faq.rst
- Update faq.rst
- Adding more docs
- Starting list of papers related to assembly parameters
- Remove 'skip' flags from meta.yaml, because False is default now
- Add funcsigs dependency
- Fix baba.py so max locus length is autodetected from the data, instead of being fixed at 300
- Adding a nexus2loci.py conversion script which takes in a directory of nexus alignments and writes out a .loci file. This is as stupid as possible and it makes a lot of assumptions about the data, so don't be surprised if it doesn't work right.
- added missing dependency on cutadapt (#314)
- Add support for finding bins in a virtualenv environment installed with pip
- add missing requirement: dask[array] (#313)
- Update faq.rst
- Update faq.rst
- fix branching docs
- Fix a nasty bug in sra tools if you try to dl more than 50 or 60 samples.
- fix dox
- Fix references to load_assembly to point to load_json
- Removing docs of preview mode
- Purge references to preview mode. Clean up some deprecated code blocks in demux.
- Remove import of util.* from load, and include only the few things it needs, remove circular dependency.
- Add docs about structure parallel runs failing silently
- Removing the restriction on ipyparallel version to obtain the 'IPython cluster' tab in notebooks.
- Adding docs about engines that die silently on headless nodes
- Add title and save ability to pca.plot()
- Make pca.plot() less chatty
- Forbid nPCs < n samples
- Update ipyrad meta.yaml to specify ipyparallel, and scikit-allel version.
- Fix pis docs in faq
- Update full_tutorial_CLI.rst
- Update full_tutorial_CLI.rst
- Update full_tutorial_CLI.rst
- Update full_tutorial_CLI.rst
- Adding scikit-allel dependency for pca analysis tool
- Update cookbook-PCA-pedicularis.ipynb
- Fix a bug that was causing _link_fastqs to fail silently.
- fixing inconsistencies in the pedicularis CLI tutorial
- Big update to the PCA cookbook.

0.7.28
------
**June 18, 2018**

- Add functions for missingness, trim missing, and fill missing.
- Adding PCA cookbook
- pcs are now stored as pandas, also, you can specify ncomps

0.7.27
------
**June 15, 2018**

- Add distance plot, and pca.pcs to hold coordinates per sample
- remove some crust from pca.pywq

0.7.26
------
**June 14, 2018**

- Adding analysis.pca
- Allow passing in just a dict for specifying populations to _link_populations(), and assume all minsamps = 0
- Some of step 2 docs were outdated
- Fix stupid link
- Adding some docs about MIG-seq.
- Damn this cluster config mayhem is a mess.
- Fix faq re pyzmq
- adding docs about max_snp_locus settings
- Fix merge conflict
- Add docs to fix the GLIBC error
- Docs for r1/r2 not the same length

0.7.25
------
**May 17, 2018**

- nb showing fix for 6-7 branching
- nb showing fix for 6-7 branching
- fixed branching between 6-7 when using populations information
- suppress h5py warning
- Allow sample names to be numbers as well.

0.7.24
------
**May 03, 2018**

- Better handling of utf-8 in sample names by default.
- Add docs in the faq about the empty varcounts array
- Catch an exception in sratools raised by non-existant sra directory.
- Add HDF5 file locking fix to the faq.
- Add docs to peddrad notebook.
- Adding PE-ddRAD analysis notebook.
- Add the right imports error message to the structure analysis tool.

0.7.23
------
**February 21, 2018**

- some releasenotes fixes
- Fix filter_min_trim_len not honoring the setting in the params file.


0.7.22
------
**February 13, 2018**

- bug fix to bpp.py
- updated tetrad cookbook
- ipa: structure has max_var_multiple option, and documentation now includes it.
- update baba cookbook
- API user guide update
- bug fix: allow for 'n' character in reftrick
- ipa: can reload structure results, better API design for summarizing results, better documentation
- allow subsetting in baba plot, and bug fix for generate_tests dynamic func
- undo dumb commit
- added --download to the docs example

0.7.21
------
**January 23, 2018**

- Fix step 2 with imported fastq ungzipped.
- docs update
- update ipa structure notebook
- update ipyparallel tutorial
- update ipa structure notebook
- docs updates
- improved cleanup on sra tools
- updated bucky cookbook
- updated --help for sra download
- updated docs for sra download

0.7.20
------
**January 09, 2018**

- fixed gphocs output format
- A note to add a feature for the future.
- abba baba cookbook updated for new code
- updated baba plot to work better with updated toytree
- baba: added functions for parsing results of 5-taxon tests and improved plotting func.
- added notes
- added CLI command to do quick downloads from SRA. Useful for tutorials especially
- update bpp cookbook
- added functions to calculate Evanno K and to exlude reps based on convergence stats
- added funcs to bpp tool to load existing results and to parse results across replicates
- ipp jobs are submitted as other jobs finish so that RAM doesn't fill up with queued arrays

0.7.19
------
**November 16, 2017**

- bugfix; error was raised in no barcodes during step2 filtering for gbs data. Now just a warning is printed
- Fixed structure conda meta.yaml
- Fix ipcluster warning message.
- Adding to the faq explaining stats better
- new working meta.yaml
- trying alternatives with setup files for jupyter conda bug fix
- updating setup.py stuff to try to fix jupyter missing in conda install

0.7.18
------
**November 13, 2017**

- allow user to set bpp binary path if different from default 'bpp'
- skip concat edits of merged reads if merge file exists unless force flag is set
- added a progress bar tracker for reference indexing
- speed improvement to refmapping, only tests merge of read pairs if their mapped positions overlap
- update to docs
- update API userguide
- added twiist tool
- update bpp notebook
- tetrad bug fix for OSX users for setting thread limit
- added check for structure path in structure.py
- allow setting binary path and check for binary added to bpp.py
- Update requirements.txt
- Added to the faq how to fix the GLIBC error.
- Fix logging of superints shape.
- Test for samples in the populations file not in the assembly.

0.7.17
------
**October 28, 2017**

- Properly handle empty chunks during alignment. Very annoying.

0.7.16
------
**October 28, 2017**

- Fix SE reference bug causing lots of rm_duplicates.
- Lowered min_se_refmap_overlap and removed useless code to recalibrate it based on filter_min_trim_len.
- Actually fix conda package.
- aslkfljsdjsdffd i don't know how this shit works.
- Fixing build still.
- Fix typo in meta.yaml.

0.7.15
------
**October 01, 2017**

- Fix conda build issue.

0.7.14
------
**September 28, 2017**

- Fix orientation of R2 for pe refmap reads.
- better error reporting, and ensure * at top of stacks
- quickfix from last commit, keep first st seq after pop to seed in align
- edge trim in s7 cuts at 4 or minsamp
- added adapter-barcode order checking for cases where merged samples, and pegbs data is analyzed either as pe or forced into se.
- update to gbs edge trimming, stricter filtering on partial overlapping seqs
- Add a comment line to the pysam conda build to make it easier to build on systems with older glibc.
-  "Updating ipyrad/__init__.py to version - 0.7.13
- API style modifications to tetrad

0.7.13
------
**September 05, 2017**

- API style modifications to tetrad

0.7.13
------
**September 04, 2017**

- Add support for optional bwa flags in hackersonly.
- Force resetting the step 6 checkpointing if step 5 is re-run.
- fix for max_shared_Hs when a proportion instead of a fixed number. Now the proportion is applied to every locus based on teh number of samples in that locus, not the total N samples
- access barcode from assembly not sample unless multiple barcodes per sample. Simpler.
- added back in core throttling in demux step b/c it is IO limited
- fix to progress bar fsck, and fix to cluster location used in step4 that was breaking if assemblies were merged between 3 and 4
- step 6 clustering uses threading options from users for really large systems to avoid RAM limits
- fix for progress bar printing in tetrad, and to args entry when no tree or map file
- fix to default ncbi sratools path

0.7.12
------
**August 28, 2017**

- update ezrad notebook
- ezrad-test notebook up
- Update cookbook-empirical-API-1-pedicularis.ipynb
- big improvements to sratools ipa, now better fetch function, easier renaming, and wraps utility to reassign ncbi dump locations
- fix for bucky bug in error reporting
- wrote tetrad CLI to work with new tetrad object
- rewrite of tetrad, cleaner code big speed improvements
- allow more flexible name entry for paired data, i.e., allow _R1.fastq, or _1.fastq instead of only _R1_.fastq, etc.
- Fixed denovo+reference assembly method.
- update bpp cookbook
- update bpp cookbook
-  "Updating ipyrad/__init__.py to version - 0.7.11
- removed repeat printing of error statements
- added more warning and reports to bpp analysis tool

0.7.11
------
**August 14, 2017**

- removed repeat printing of error statements
- added more warning and reports to bpp analysis tool

0.7.11
------
**August 14, 2017**

- better error checking in bucky run commandipa tools
- added workdir default name to sra tools ipa tool
- improved error checking in step 6
- bugfix for VCF output where max of 2 alternative alleles were written although there could sometimes be 3

0.7.10
------
**August 08, 2017**

- fix misspelled force option in ipa bucky tool
- bpp ipa tool changed 'locifile' arg to 'data' but still support old arg, and removed 'seed' arg from run so that the only 'seed' arg is in paramsdict
- bugfix to not remove nex files in bucky ipa tool

0.7.9
-----
**August 07, 2017**

- cleaner shutdown of tetrad on interrupt. Bugfix to stats counter for quartets sampled value. Cleaner API access by grouoping attributes into params attr.
- cleanup rawedit to shutdown cleaner when interrupted
- modified run wrapper in assembly object to allow for cleaner shutdown of ipyclient engines
- bug fix so that randomize_order writes separate seqfiles for each rep in bpp analysis tool
- Adding error handling, prevent tmp files being cleaned up during DEBUG, and fix tmp-align files for PE refmap.
- Derep and cluster 2brad on both strands.
- Actually fix refmap PE merging.
- Fix merging for PE refmap.
- Add a switch to _not_ delete temp files if DEBUG is on. Helpful.
- 2 new merge functions for PE refmap. One is slowwwww, the other uses pipes, but doesn't work 100% yet.
- New hackersonly parameter to switch merging PE after refmap.
- bugfix to ipa baba plotting function for updated toyplot
- Reduce minovlen length for merging reference mapped PE reads.
- docs update
- docs update
- improved design of --ipcluster flag in tetrad CLI
- improved design of --ipcluster flag in tetrad CLI
- improved design of --ipcluster flag in ipyrad CLI

0.7.8
-----
**July 28, 2017**

- bpp randomize-order argument bugfix
- added .draw to treemix object
- update tuts

0.7.7
-----
**July 27, 2017**

- Proper support for demux 2brad.

0.7.6
-----
**July 27, 2017**

- Fix very nasty refmap SE bug.
- update tutorials -- added APIs
- update tutorials -- added APIs
- testing MBL slideshow
- API cookbooks updated
- cleanup of badnames in sratools

0.7.5
-----
**July 26, 2017**

- Added error handling in persistent_popen_align3
- Catch bad seeds in step 6 sub_build_clustbits().

0.7.4
-----
**July 26, 2017**

- Actually fix the step 6 boolean mask error.
- Fix for boolean mask array length bug in step 6.
- add -noss option for treemix ipa
- mods to tetrad and sratools ipa
- ensure ints not floats for high depth base counts
- sratools updates
- improvements to sratools
- added extra line ending to step7 final print statement
- add dask to environment.yaml
- added sratools to ipyrad.analysis

0.7.3
-----
**July 23, 2017**

- Better handling for restarting jobs in substeps of step 6.
- Fixed the fscking pysam conda-build scripts for osx.
- Add patch for pysam build on osx
- Fix for conda-build v3 breaking meta.yaml
- Using htslib internal to pysam and removing bcftools/htslib/samtools direct dependencies.
- Add force flag to force building clusters if utemp exists.
- conda recipe updates
- updateing conda recipe
- ensure stats are saved as floats
- fix to bug introduced just now to track progress during s6 clustering
- Fix an issue with merged assemblies and 3rad.
- fix for step 6 checkpoints for reference-based analyses
- conda recipe tweaking
- conda recipe updates
- fix to conda recipes
- update bucky cookbook
- added shareplot code
- bucky ipa update remove old files
- conda recipe updated
-  "Updating ipyrad/__init__.py to version - 0.7.2
- update conda recipe
- update pysam to correct version

0.7.2
-----
**July 10, 2017**

- update conda recipe
- update pysam to correct version (0.11.2.2)
- added bucky ipa code
- bucky cookbook up
- automatically merges technical replicates in demux
- check multiple barcodes in samples that were merge of technical replicates
- fix for alleles output error
- Added checkpointing/restarting from interrupt to step 6. 
- Added cli detection for better spacer printing. 
- bpp bug fixe to ensure full path names
- API user guide docs update.
- cookcook updates tetrad and treemix
- new _cli, _checkpoint, and _spacer attributes, and new 'across' dir for step 6
- load sets cli=False by default, and it saves checkpoint info
- allow profile with ipcluster
- treemix report if no data is written (i.e., all filtered)
- fix to allow setting nquartets again. 
- Better integration of API/CLI. 
- Bug fix to Tree drawing when no boots in tetrad. 
- tetrad fix for compatibility with new toytree rooting bug fix for saving features.
- cli is now an attribute of the Assembly object that is set to True by __main__ at runtime, otherwise 0.
- cluster_info() now prints instead of return
- rehaul of bucky ipa tools
- print cluster_info now skips busy engines
- unroot tetrad tree on complete

0.7.1
-----
**June 16, 2017**

- Actually handle SE reference sequence clustering.
- Prevent empty clust files from raising an error. Probably only impacts sim data.
- If debug the retain the bed regions per sample as a file in the refmap directory.
- updated tunnel docs
- HPC tunnel update
- support for parsing supervised structure analyses in ipa
- HPC tunnel docs update
- update analysis docs
- ipa.treemix params
- more params added to ipa.treemix
- cookbook update treemix
- fix to conda rec
- treemix ipa updates

0.7.0
-----
**June 15, 2017**

- put a temporary block on denovo+ref
- added treemix ipa funcs
- update conda recipe
- added notebook for structure with popdata
- updated tetrad notebook
- update bpp notebook
- fix missing newline in alleles
- ipa structure file clobber fix
- cleaner and more consistent API attr on ipa objects
- Added docs for the -t flag.
- fix in ipa.structure so replicate jobs to do not overwrite
- Fix bad link in docs.
- better method to find raxml binary in analysis tools
- consens bugfix for new ipmlementation
- ensure h5 files are closed after dask func
- fix to parse chrom pos info from new consens name format
- removed deprecated align funcs
- removed hardcoded path used in testing
- removed deprecated align funcs. Made it so build_clusters() does nothing for 'reference' method since there is a separate method in ref for chunking clusters
- some new simpler merge funcs
- make new ref funcs work with dag map
- new build funcs usign pysam

0.6.27
------
**June 03, 2017**

- Step 6 import fullcomp from util.

0.6.26
------
**June 01, 2017**

- Step 4 - Handle the case where no clusters have sufficient depth for statistical basecalling.

0.6.25
------
**May 30, 2017**

- Fix a bug in refmap that was retaining the reference sequence in the final clust file on rare occasions.

0.6.24
------
**May 25, 2017**

- Bug fix for "numpq" nameerror

0.6.23
------
**May 24, 2017**

- bug fix for numq error in s5

0.6.22
------
**May 22, 2017**

- Fixed bug in vcf output for reference mapped.

0.6.21
------
**May 19, 2017**

- Fix new chrom/pos mechanism to work for all assembly methods.
- Change chroms dtype to int64. Reference sequence CHROM is now 1-indexed. Anonymous loci are -1 indexed.
- Switch chroms dataset dtype to int64.
- Fix for alleles output.
- Fix nasty PE refmap merging issue.
- Fix massive bug in how unmapped reads are handled in refmap.
- added md5 names to derep and simplified code readability within pairmerging
- fix for binary finder
- added dask to conda recipe
- added dask dependency

0.6.20
------
**May 10, 2017**

- added dask dependency
- vcf building with full ref info
- bug fix to alleles output and support vcf chrompos storage in uint64
- simpler and slightly faster consens calls and lower memory and stores chrompos as uint64s
- chrompos now stored as uint64
- reducing memory load in race conditions for parallel cutadapt jobs
- Squash Cosmetic commit logs in releasenotes. Add more informative header in step 7 stats file.
- Trying to catch bad alignment for PE in step 6.

0.6.19
------
**May 04, 2017**

- Handle empty locus when building alleles file. Solves the ValueError "substring not found" during step 7.
- workshop notebook uploaded

0.6.18
------
**May 03, 2017**

- update to analysis tools
- accepted the local bpp notebook
- complete bpp notebook up
- notebook updates
- raxml docs
- raxml cookbook up
- docs update
- raxml docs updated
- links to miniconda updated
- fix for tetrad restarting bootstraps
- removed bitarray dependency
- adding restart checkpoints in step6

0.6.17
------
**April 26, 2017**

- support for alleles file in bpp tools
- align names in alleles output
- bugfix to name padding in .alleles output
- slight delay between jobs
- bpp store asyncs
- bpp store asyncs
- update bpp cookbook
- testing html
- testing html
- new filter_adapters=3 option adds filtering of poly-repeats
- conda recipe update for cutadapt w/o need of add-channel

0.6.16
------
**April 25, 2017**

- alleles output now supported
- Additional documentation for max_alleles_consens parameter.
- support alleles output, minor bug fixes for step6, much faster alignment step6
- lower default 'cov' value for vsearch within clustering in RAD/ddrad/pairddrad
- tetrad bug, use same ipyclient for consensus tree building
- store asyncs in the structure object
- allow passing in ipyclient explicitly in .run() in tetrad
- fix for time stamp issue in tetrad
- Better testing for existence of all R2 files for merged assemblies.
- notebook updates
- tunnel docs update
- updated HPC docs
- tetrad cookbook updated
- HPC docs update
- bpp cookbook good to go
- update tetrad notebook
- missing import

0.6.15
------
**April 18, 2017**

- Actually fix gphocs output.
- allow passing in ipyclient in API
- baba notebook update
- cleaner api for bpp object
- new analysis setup
- updated analysis tools without ete
- adding doc string

0.6.14
------
**April 13, 2017**

- Fixed CHROM/POS output for reference mapped loci.

0.6.13
------
**April 13, 2017**

- Fix gphocs output format.
- If the user removes the population assignment file blank out the data.populations dictionary.

0.6.12
------
**April 10, 2017**

- Prevent versioner from including merge commits in the release notes cuz they are annoying.
- Add the date of each version to the releasenotes docs, for convenience.
- Experimenting with adding date to releasenotes.rst
- added more attributres to tree
- change alpha to >=
- tip label and node label attributes added to tree
- tetrad ensure minrank is int
- fix structure obj removing old files
- lots of cleanup to baba code
- edit to analysis docs
- Handle pop assignment file w/o the min sample per pop line.
- merge conflict resolved
- bug fix for tuples in output formats json
- sim notebook started
- cookbook abba-baba updated
- tetrad cookbook api added
- added option to change line spacing on progress bar
- major overhaul to ipyrad.analysis and plotting
- option to buffer line spacing on cluster report
- Removed confusing punctuation in warning message
- Make vcf and loci output files agree about CHROM number per locus.
- Cosmetic change to debug output.
- Make the new debug info append instead of overwrite.
- Fix annoying bug with output_format param I introduced recently.
- Add platform info to default log output on startup.
- Actually write the error to the log file on cutadapt failure.
- Write the version and the args used to the log file for each run. This might be annoying, but it could be useful.
- bpp randomize option added to write
- adding bpp cookbook update
- updating analysis tools for new bpp baba and tree
- merge resolved
- analysis init update for new funcs
- apitest update
- abba cookbook update
- update bpp cookbook
- small edit to HPC docs
- tetrad formatting changing
- updated analysis tools cookbooks
- docs analysis page fix
- added header to bpp convert script

0.6.11
------
**March 27, 2017**

- Fix a bug in PE refmapping.
- Fix error reporting if when testing for existence of the clust_database file at beginning of step 7.
- Fix bug reading output formats from params file.
- Add docs for dealing with long running jobs due to quality issues.
- bug fix for output format empty
- structure cookbook update
- pushing analysis tools
- svg struct plot added
- structure cookbook updates
- struct image added for docs
- update structure cookbook for new code
- Actually fix the output_format default if blank.
- Set blank output formats in params to default to all formats.
- Add a filter flag for samtools to push secondary alignments to the unmapped file.
- rm old files
- shareplot code in progress
- work in progress baba code notebook
- a decent api intro but bland
- beginnings of a migrate script
- raxml docs updated, needs work still
- analysis docs page update
- structure parallel wrapper scripts up in analysis
- simplifying analysis imports
- cleanup top imports
- Adding support for G-PhoCS output format.
- Fix wacky reporting of mapped/unmapped reads for PE.
- Document why we don't write out the alleles format currently.
- module init headers
- added loci2cf script
- update structure notebook with conda recipes
- fileconversions updated
- loci2cf func added
- cookbook bucky docs up
- loci2multinex and bucky notebook updated
- BUCKy cookbook updated
- bucky conda recipe up
- fix to API access hint
- cleaner code by moving msgs to the end
- slight modification to paired adapter trimming code
- cleaner Class Object in baba
- minor change to cluster_info printing in API

0.6.10
------
- Filter reference mapped reads my mapq < 30, and handle the occasional malformed region string in bam_region_to_fasta.
- Handle PE muscle failing alignment.
- Cosmetic faq.rst
- Cosmetic faq.rst
- Cosmetic
- Cosmetic docs changes.
- Add docs for step 3 crashing bcz of lack of memory.
- Catch a bug in alignment that would crop up intermittently.
- removed the --profile={} tip from the docs
- Fix notebook requirement at runtime error.
- Fix formatting of output nexus file.

0.6.9
-----
- Changed the sign on the new hackersonly parameter min_SE_refmap_overlap.
- added a persistent_popen function for aligning, needs testing before implementing
- debugger in demux was printing way too much
- bugfix for empty lines in branching subsample file
- Add a janky version checker to nag the user.

0.6.8
-----
- Actually remove the reference sequence post alignment in step 3. This was BREAKING STUFF.
- updated notebook requirement in conda recipe
- Handle conda building pomo on different platforms.
- Ooops we broke the versioner.py script. Now it's fixed.
- conda recipe updates
- conda recipe updates
- conda recipe updates
- testing git lfs for storing example data

0.6.7
-----
- Fixed stats reported for filtered_by_depth during step 5.
- Add new hackersonly parameter min_SE_refmap_overlap and code to refmap.py to forbid merging SE reads that don't significantly overlap.
- Use preprocessing selectors for linux/osx for clumpp.
- Add url/md5 for mac binary to clumpp meta.yaml
- conda recipes update
- getting ipyrad to conda install on other envs
- updating versions for conda, rtd, setup.py
- moving conda recipes
- conda recipe dir structure
- bpp install bug fix
- bpp recipe fix
- conda recipes added
- Roll back change to revcomp reverse strand SE hits. Oops.
- fix merge conflect with debug messages.
- Fix a bug in refmap, and handle bad clusters in cluster_within.
- Actually revcomp SE - strand reads.
- updated HPC docs
- updated HPC docs
- updated HPC docs

0.6.6
-----
- bug fix in building_arrays where completely filtered array bits would raise index error -1
- tunnel docs updates
- method docs updated to say bwa
- some conda tips added
- fix for name parsing of non gzip files that was leaving an underscore
- Allow get_params using the param string as well as param index
- Update hpc docs to add the sleep command when firing up ipcluster manually.
- Fixed some formatting issues in the FAQ.rst.

0.6.5
-----
- Fixed 2 errors in steps 3 and 4.
-  "Updating ipyrad/__init__.py to version - 0.6.4
- left a debugging print statement in the code
- removed old bin
-  "Updating ipyrad/__init__.py to version - 0.6.4

0.6.4
-----
- left a debugging print statement in the code
- removed old bin
-  "Updating ipyrad/__init__.py to version - 0.6.4

0.6.4
-----

0.6.4
-----
- update to docs parameters
- bug fix for merging assemblies with a mix of same named and diff named samples

0.6.3
-----
- Fixed a bug i introduced to assembly. Autotroll.

0.6.2
-----
- Fix subtle bug with migration to trim_reads parameter.

0.6.1
-----
- Fixed malformed nexus output file.
- cookbook updates to docs
- updated cookbook structure pedicularis

0.6.0
-----
- trim reads default 0,0,0,0. Similar action to trim loci, but applied in step 2 to raws
- trim_reads default is 0,0
- raise default cov/minsl for gbs data to 0.5 from 0.33
- prettifying docs
- pedicularis docs update v6 way way faster
- updated tutorial
- fixing links in combining data docs
- updating tutorial for latest version/speed
- added docs for combining multiple plates
- added docs for combining multiple plates
- added docs for combining multiple plates
- Removed  from output formats defaults (it doesn't do anything)
- baba cookbooks [unfinished] up
- finally added osx QMC and fixed bug for same name and force flag rerun
- put back in a remove tmpdirs call
- removed a superfluous print statement
- bug fix to mapfile, now compatible with tetrad
- paramsinfo for new trimreads param
- branching fix for handling new param names and upgrading to them
- better handling of pairgbs no bcode trimming. Now handles --length arg
- better handling of KBD in demux. Faster compression.
- forgot sname var in cutadaptit_single
- Fix step 2 for PE reads crashing during cutatapt.
- Test for bz2 files in sorted_fastq_path and nag the user bcz we don't support this format.
- Step 1 create tmp file for estimating optim chunk size in project_dir not ./
- Add force flag to mapreads(), mostly to save time on rerunning if it crashes during finalize_mapping. Also fixed a nasty bug in refmapping.
- Added text to faq about why PE original RAD is hard to assemble, cuz people always ask.
- Better handling of loci w/ duplicate seqs per sample.
- Fix a bug that munged some names in branching.
- merge conflict
- modified for new trim param names
- support for new trim_loci param
- support for updated cutadapt
- bugfix for hackerdict modify of cov
- chrom only for paired data
- changed two parameter names (trims)
- tested out MPI checks
- cutadapt upgrade allow for --length option
- Moved log file reset from init to main to prevent -r from blanking the log >:{
- Moved log file reset from __init__ to __main__
- Don't bother aligning clusters with duplicates in step 6.
- baba update
- remove print statement left in code
- same fix to names parser, better.
- added comment ideas for chrompos in refmap
- bug fix, Sample names were being oversplit if they had '.' in them
- test labels, improved spacing, collapse_outgroups options added to baba plots
- Fix debug message in refmap and don't raise on failure to parse reference sequence.
- attempts to make better cleanup for interrupt in API
- some cleanup to calling steps 1,2 funcs
- speed testing demux code with single vs multicore
- moved setting of ['merged'] to replace filepath names to Assembly instead of main so that it also works for the API
- added a np dict-like arr to be used in baba, maybe in ref.
- baba plotting functions added
- Better handling of tmpdir in step 6.
- added baba cookbook
- only map chrom pos if in reference mode
- new batch and plotting functions
- trim .txt from new branch name if accidentally added to avoid Assembly name error
- added a name-checker to the branch-drop CLI command
- Fixed legend on Pedicularis manuscript analysis trees.
- Cosmetic change
- Adding manuscript analysis tree plotting for empirical PE ddRAD refmap assemblies.
- More or less complete manuscript analysis results.
- Actually fix vcf writing CHROM/POS information from refseq mapped reads.
- Handle monomorphic loci during vcf construction.
- removed deprecated subsample option from jointestimate
- --ipcluster method looks for default profile and cluster-id instance
- clode cleanup and faster haploid E inference
- simplified cluster info printing
- enforce ipyclient.shutdown at end of API run() if engine jobs are not stopped
- code cleanup. Trying to allow better KBD in step2
- lots of cleanup to DAG code. Now ok for individual samples to fail in step3, others will continue. Sorts clusters by derep before align chunking
- Allow assemblies w/o chrom/pos data in the hdf5 to continue using the old style vcf position numbering scheme.
- Don't print the error message about samples failing step 4 if no samples actually fail.
- Set a size= for reference sequence to sort it to the top of the chunk prior to muscle aligning.
- Allow samples with very few reads to gracefully fail step 4.
- Better error handling during reference mapping for PE.
- Fix error reporting in merge_pairs().
- Add CHROM/POS info to the output vcf file. The sorting order is a little wonky.
- Handle empty project_dir when running -r.
- a clean bighorse notebook run on 100 cores
- Fix minor merge conflict in ref_muscle_chunker.
- Use one persistant subprocess for finalizing mapped reads. Big speed-up. Also fix a stupid bug in estimating insert size.
- Better handling of errors in merge_pairs, and more careful cleanup on error.
- If /dev/shm exists, use it for finalizing mapped reads.
- Handle a case where one or the other of the PE reads is empty.
- cleaner print cpus func
- Adding a new dataset to the catg and clust hdf5 files to store CHROM and POS info for reference mapped reads.
- added cleanhorse notebook
- working on notebook
- cleanup up redundancy
- MUCH FASTER STEP 4 using numba array building and vectorized scipy
- MUCH FASTER MUSCLE ALIGNING. And a bug fix to a log reporter
- bug fix to error/log handler
- Finish manuscript refmap results analysis. Added a notebook for plotting trees from manuscript Pedicularis assembly.
- Better checking for special characters in assembly names, and more informative error message.
- added a test on big data
- broken notebook
- development notebook for baba
- working on shareplots
- testing caching numba funcs for faster run starts
- added optional import of subprocess32
- docs update
- progress on baba
- added option to add additional adapters to be filtered from paired data
- Adding pairwise fst to manuscript analysis results. Begin work on raxml for manuscript analysis results.
- Change a log message from info to warn that handles exceptions in rawedit.
- abba baba updated
- Fixed link in tetrad doc and cosmetic change to API docs.
- Add comments to results notebooks.
- Adding manuscript reference mapping results.
- Manuscript analysis reference sequence mapping horserace updates. Stacks mostly done. dDocent started.
- Adding ddRAD horserace nb.
- Better cleanup during refmap merge_pairs (#211).
- update for raxml-HYBRID
- update raxml docs
- cleanup old code
- update raxml docs
- updating raxml docs
- update to bucky cookbook

0.5.15
------
- bug fix to ensure chunk size of the tmparray in make-arrays is not greater than the total array size
- fix for vcf build chunk error 'all input arrays must have the same number of dimensions'. This was raised if no loci within a chunk passed filtering
- allow vcf build to die gracefully
- api cleanup

0.5.14
------
- updated docs for popfile
- fix for long endings on new outfile writing method
- Made max size of the log file bigger by a zero.
- Be nice and clean up a bunch of temporary files we'd been leaving around.
- Better handling for malformed R1/R2 filenames.
- api notebook update
- more verbose warning on ipcluster error
- allow setting ipcluster during Assembly instantiation
- improved populations parser, and cosmetic
- greatly reduced memory load with new func boss_make_arrays that builds the arrays into a h5 object on disk, and uses this to build the various output files. Also reduced disk load significantly by fixing the maxsnp variable bug which was making an empty array that was waay to big. Also added support for nexus file format. Still needs partition info to be added.
- CLI ipcluster cluster-id='ipyrad-cli-xxx' to more easily differentiate from API
- added note on threading
- API cleanup func names
- write outfiles h5 mem limit work around for build-arrays
- step 1 with sorted-fastq-path no longer creates empty fastq dirs

0.5.13
------
- API user guide updated
- Added ipyclient.close() to API run() to prevent 'too many files open' error.
- Bug fix for concatenation error in vcf chunk writer
- added smarter chunking of clusters to make for faster muscle alignments
- closed many subprocess handles with close_fds=True
- added closure for open file handle
- cleanup of API attributes and hidden funcs with underscores

0.5.12
------
- Refmap: actually fix clustering when there are no unmapped reads.
- Updated docs for parameters.

0.5.11
------
- Refmap: Handle case where all reads map to reference sequence (skip unmapped clustering).
- More refined handling of reference sequences with wacky characters in the chrom name like | and (. Who would do that?
- Raxml analysis code added to Analysis Tools: http://ipyrad.readthedocs.io/analysis.html
- HPC tunneling documentation updated with more troubleshooting
- Better handling of final alignments when they contain merged and unmerged sequences (#207)
- added finetune option to loci2bpp Analysis tools notebook.
- More improvements to manuscript analysis.
- Finished simulated analysis results and plotting.
- Improve communication if full raw path is wonky.
- Horserace is complete for simulated and empirical. Continued improvement to gathering results and plotting.

0.5.10
------
- Fix for 3Rad w/ only 2 cutters during filtering.
- Better handling for malformed 3rad barcodes file.

0.5.9
-----

0.5.8
-----
- improved progress bar
- merge fix
- notebook testing geno build
- Fix to memory handling on vcf build, can now handle thousands of taxa. Also, now saves filepaths to json and API object.
- progres on dstats package
- More progress on manuscript horserace. Analysis is done, now mostly working on gathering results.

0.5.7
-----
- Fix error handing during writing of vcf file.

0.5.6
-----
- notebook testing
- purge after each step to avoid memory spillover/buildup
- better handling of memory limits in vcf build. Now producing geno output files. Better error reporting when building output files
- added a global dict to util
- new smaller limit of chunk sizes in h5 to avoid memory limits
- analysis docs update
- Document weird non-writable home directory on cluster issues.
- docs update for filtering differences
- merge fix
- tetrad notebook edits
- dstat calc script editing
- Added code to copy barcodes during assembly merge. Barcodes are needed for all PE samples in step 2.

0.5.5
-----
- Better handling for PE with loci that have some merged and some unmerged reads.
- Allow other output formats to try to build if vcf fails.
- Fixed bug that was forcing creation of the vcf even if it wasn't requested.

0.5.4
-----
- More improved handling for low/no depth samples.
- Better handling for cleanup of samples with very few reads.

0.5.3
-----
- Catch sample names that don't match barcode names when importing demux'd pair data.
- Serious errors now print to ipyrad_log.txt by default.

0.5.2
-----
- Handle sample cleanup if the sample has no hidepth clusters.
- Fix for declone_3rad on merged reads.
- Better support for 3rad lining presorted fastqs.
- bucky cookbook updated
- dstat code updates
- bucky cookbook uploaded

0.5.1
-----
- added tetrad docs
- make tetrad work through API
- added tetrad notebook

0.5.0
-----
- Swap out smalt for bwa inside refmapping. Also removes reindexing of reference sequence on -f in step 3.
- fix for array error that was hitting in Ed's data, related to 2X count for merged reads. This is now removed.
- bug fix for 4/4 entries in vcf when -N at variable site.
- prettier printing of stats file

0.4.9
-----
- fix for array error that was hitting in Ed's data, related to 2X count for merged reads. This is now removed.
- bug fix for 4/4 entries in vcf when -N at variable site.
- prettier printing in s5 stats file
- hotfix for large array size bug introduced in 0.4.8


0.4.8
-----
- bug fix to measure array dims from mindepth settings, uses statistical for s4, and majrule for s5
- adding bwa binary for mac and linux
- improved N removal from edges of paired reads with variable lengths
- new parsing of output formats, and fewer defaults
- only snps in the vcf is new default. Added pair support but still need to decide on spacer default. New cleaner output-formats stored as a tuple
- small fix for better error catching
- new hidepth_min attr to save the mindepth setting at the time when it is used
- mindepth settings are now checked separately from other parameters before 'run' to see if they are incompatible. Avoids race between the two being compared individually in set-params.
- new functions in steps 3-5 to accomodate changes to mindepth settings so that clusters-hidepth can be dynamically recalculated
- fix to SSH tunnel docs
- hotfix for step5 sample save bug. pushed to soon

0.4.7
-----
- make compatible with changes to s6
- allow sample to fail s2 without crashing
- cleaner progress bar and enforced maxlen trimming of longer reads
- lowered maxlen addon, enforced maxlen trimming in singlecat
- updates to docs
- testing new maxlen calculation to better acommodate messy variable len paired data sets.
- update to docs about pre-filtering
- temporary fix for mem limit in step 6 until maxlen is more refined
- Fix bug in refmap.

0.4.6
-----
- Nicely clean up temp files if refmap merge fails.

0.4.5
-----
- Add docs for running ipcluster by hand w/ MPI enabled.
- Fix PE refmap bug #148
- Documenting PYTHONPATH bug that crops up occasionally.
- Adjusted fix to bgzip test.
- Fixed a bug w/ testing for bgzip reference sequence. Also add code to fix how PE ref is handled to address #148.
- fix for last fix
- fix for last push gzip
- collate with io.bufferedwriter is faster
- faster collating of files
- Continuing work on sim and empirical analysis.
- rev on barcode in step2 filter pairgbs
- faster readcounter for step1 and fullcomp on gbs filter=2 barcode in step2
- tunnel docs update
- working on a SSH tunnel doc page
- Handle OSError in the case that openpty() fails.

0.4.4
-----
- Handle blank lines at the top of the params file.

0.4.3
-----
- making smoother progress bar in write vcfs
- bugfix for jointestimate
- testing bugfixes to jointestimate
- default to no subsampling in jointestimate call
- testing bugfixes to jointestimate
- added hackersonly option for additional adapters to be filtered
- bug fix to joint H,E estimate for large data sets introduced in v.0.3.14 that was yielding inflated rates.
- fix for core count when using API
- Added plots of snp depth across loci, as well as loci counts per sample to results notebook.
- phylogenetic_invariants notebook up
- some notes on output formats plans
- removed leftjust arg b/c unnecessary and doesn't work well with left trimmed data

0.4.2
-----
- Merging for Samples at any state, with warning for higher level states. Prettier printing for API. Fix to default cores setting on API.
- fix for merged Assemblies/Samples for s2
- fix for merged Assemblies&Samples in s3
- removed limit on number of engines used during indexing
- Added ddocent to manuscript analysis.
- tutorial update
- in progress doc notebook
- parallel waits for all engines when engines are designated, up until timeout limit
- parallelized loading demux files, added threads to _ipcluster dict, removed print statement from save
- vcf header was missing
- added step number to progress bar when in interactive mode
- added warning message when filter=2 and no barcodes are present
- improved kill switch in step 1
- use select to improve cluster progress bar
- added a CLI option to fine-tune threading
- added dstat storage by default
- new default trim_overhang setting and function (0,0,0,0)
- fix for overzealous warning message on demultiplexing when allowing differences

0.4.1
-----
- Fixed reference before assignment error in step 2.

0.4.0
-----
- Cosmetic change
- new sim data and notebook up
- Added aftrRAD to the manuscript analysis horserace
- made merging reads compatible with gzipped files from step2
- modify help message
- made TESTS global var, made maparr bug fix to work with no map info
- More carefully save state after completion of each step.
- limit vsearch merging to 2 threads to improve parallel, but should eventually make match to cluster threading. Added removal of temp ungzipped files.
- more detailed Sample stats_df.s2 categories for paired data
- made merge command compatible with gzip outputs from step2
- simplified cutadapt code calls
- updates to simdata notebook
- merge conflict fix
- new stats categories for step2 results
- added adapter seqs to hackersdict
- much faster vcf building
- new step2 quality checks using cutadapt
- small changes to use stats from new s2 rewrite. Breaks backwards compatibility with older assemblies at step3
- massive rewrite of cluster across, faster indexing, way less memory overhead
- just added a pylint comment
- Adding cutadapt requirement for conda build
- Suppress numpy mean of empty slice warnings.
- Merged PR from StuntsPT. Fix to allow param restriction_overhang with only one enzyme to drop the trailing comma (,).
- Merge branch 'StuntsPT-master'
- Adding a FAQ to the docs, including some basic ipyparallel connection debugging steps.
- Adding documentation for the  CLI flag for attaching to already running cluster.
- Update docs to include more specifics about ambiguous bases in restriction overhang seqs.
- Get max of max_fragment_length for all assemblies during merge()
- Make gbs a special case for handling the restriction overhang.
- Changed the way single value tuples are handled.
- cleaning up releasenotes
- added networkx to meta.yaml build requirements
-  "Updating ipyrad/__init__.py to version - 0.3.42

0.3.42
------
- always prints cluster information when not using ipcluster[profile] = default
- broke and then fixed samtools sorting on mac (BAM->bam)
- better error message at command line
- cleaned code base, deleting deprecated funcs.
- revcomp function bug fix to preserve lower case pair splitter nnnn for pairgbs data
- Adding requirement for numba >= 0.28 to support
- Updating mac and linux vsearch to 2.0
- docs updates (pull request #186) from StuntsPT/master
- Added a troubleshooting note.
- wrapped long running proc jobs so they can be killed easily when engines are interrupted
- fix for API closing ipyclient view
- fix for piping in subprocess
- bug fix for missing subprocess module for zcat, and new simplified sps calls.
- merge fix
- allow for fuzzy match characters in barcode path
- new simulated data set
- uploaded cookbook for simulating data
- no longer register ipcluster to die at exit, but rather call shutdown explicitly for CLI in the finally call of run()
- massive code cleanup in refmapping, though mostly cosmetic. Simplified file paths and calls to subprocess.
- massive restructuring to organize engine jobs in a directed acyclic graph to designate dependencies to ipyparallel. Lot's of code cleanup for subprocess calls.
- fix for progress bar cutting short in step 6. And simplified some code calling tmpdir.
- Adding notebooks for ipyrad/pyrad/stacks simulated/emprical horserace.
- Better handling for mindepth_statistical/majrule. Enforce statistical >= majrule.
- Allow users with SE data to only enter a single value for edit_cutsites.
- Properly finalize building database progress bar during step 6, even if some samples fail.
- allow max_indels option for step 3 in API. Experimental.
- bug fix to indel filter counter. Now applies in step7 after ignoring terminal indels, only applies to internal indels
- much faster indexing using sorted arrays of matches from usort. Faster and more efficient build clusters func.
- rewrote build_clusters func to be much faster and avoid memory limits. Other code cleanup. Allow max_indel_within option, though only in API currently.
- numba update requirement to v.0.28

0.3.41
------
- Reverting a change that broke cluster_within

0.3.40
------
- Set vsearch to ignore max phred q score on merging pairs
- Added bitarray dependency to conda build

0.3.39
------
- Fix vsearch fastq max threshold arbitrarily high. Also remove debug crust.

0.3.38
------
- Handle samples with few reads, esp the case where there are no matches during clustering.
- Handle samples with few or no high depth reads. Just ignore them and inform the user.

0.3.37
------
- Fix to allow pipe character  in chrom names of reference sequences
- Tweak to calculation of inner mate distance (round up and cast to int)
- Refmap: fix calc inner mate distance PE, handle samples w/ inner mate distance > max, and handle special characters in ref seq chromosome names
- Add a test to forbid spaces in project directory paths
- Cosmetic docs fix
- Cosmetic fix to advanced CLI docs
- Added more explicit documentation about using the file to select samples during branching
- Clarifying docs for qscore offset in the default params file
- Cosmetic change to docs
- Rolling back changes to build_clusters
-  "Updating ipyrad/__init__.py to version - 0.3.36
- hotfix for edgar fix break
-  "Updating ipyrad/__init__.py to version - 0.3.36
- hotfix for edgar fix break

0.3.36
------
- hotfix for edgar fix break
-  "Updating ipyrad/__init__.py to version - 0.3.36
- hotfix for edgar fix break

0.3.36
------
- hotfix for edgar fix break

0.3.36
------
- hotfix for memory error in build_clusters, need to improve efficiency for super large numbers of hits
- more speed testing on tetrad
- merge conflict
- cleaner print stats for tetrad
- finer tuning of parallelization tetrad

0.3.35
------
- Handled bug with samtools and gzip formatted reference sequence
- Fixed a bug where CLI was not honoring -c flag
- debugging and speed tests
- added manuscript dir
- Update on Overleaf.
- Manuscript project created
- speed improvements to tetrad
- smarter/faster indexing in tetrad matrix filling and speed up from skipping over invariant sites
- finer tuning of bootstrap restart from checkpoint tetrad
- print bigger trees for tetrad
- fix to printing checkpoint info for tetrad
- bug fix for limiting n cores in tetrad
- made an extended majority rule consensus method for tetrad to avoid big import packages just for this.
- testing timeout parallel
- test notebook update
- adding consensus mj50 function

0.3.34
------
- new --ipcluster arg allows using a running ipcluster instance that has profile=ipyrad
- temporary explicit printing during ipcluster launch for debugging
- also make longer timeout in _ipcluster dict of Assembly object

0.3.33
------
- temporary explicit printing during ipcluster launch for debugging
-  "Updating ipyrad/__init__.py to version - 0.3.33
- also make longer timeout in _ipcluster dict of Assembly object

0.3.33
------
- also make longer timeout in _ipcluster dict of Assembly object

0.3.33
------
- increased timeout for ipcluster instance from 30 seconds to 90 seconds
- Added sample populations file format example
- quick api example up
- merge conflict
- removed chunksize=5000 option
- Update README.rst

0.3.32
------
- Fix optim chunk size bug in step 6 (very large datasets overflow hdf5 max chunksize 4GB limit)
- Doc update: Cleaned up the lists of parameters used during each step to reflect current reality.
- Fixed merge conflict in assembly.py
- Fix behavior in step 7 if requested samples and samples actually ready differ
- Removing references to deprecated params (excludes/outgroups)
- Simple error handling in the event no loci pass filtering
- changed tetrad default mode to MPI
- release notes update

0.3.31
------
- changed name of svd4tet to tetrad
- improved message gives info on node connections for MPI
- added a test script for continuous integration
- big cleanup to ipcluster (parallel) setup, better for API/CLI both
- modified tetrad ipcluster init to work the same as ipyrad's
- generalized ipcluster setup

0.3.30
------
- Changed behavior of step 7 to allow writing output for all samples that are ready. Allows the user to choose whether to continue or quit.
- Fixed very stupid error that was not accurately tracking max_fragment_length.
- Better error handling on malformed params file. Allows blank lines in params (prevents that gotcha).
- Cosmetic changes to step 7 interaction if samples are missing from db
- prettier splash
- edited splash length, added newclient arg to run
- testing MPI on HPC multiple nodes
- updating docs parameters

0.3.29
------
- Temp debug code in jointestimate for tracking a bug
- Step 5 - Fixed info message for printing sample names not in proper state. Cosmetic but confusing.

0.3.28
------
- Added statically linked binaries for all linux progs. Updated version for bedtools and samtools. Updated vsearch but did not change symlink (ipyrad will still use 1.10)
- Bugfix that threw a divide by zero error if no samples were actually ready for step 5

0.3.27
------
- Fixed a race condition where sometimes last_sample gets cleaned up before the current sample finishes, caused a KeyError. Very intermittent and annoying, but should work now

0.3.26
------
- fix merge conflict
- removed future changes to demultiplex, fixed 1M array size error
- added notes todo
- removed unnecessary imports
- removed backticks from printouts
- removed backticks from printouts
- removed unnecessary '\' from list of args
- code cleanup for svd4tet
- update to some error messages in svd4tet
- slight modification to -n printout
- updated analysis docs
- minor docs edits
- updated releasenotes

0.3.25
------
- better error message if sample names in barcodes file have spaces in them
- VCF now writes chr ('chromosomes' or 'RAD loci') as ints, since vcftools and other software hate strings apparently
- fix for concatenating multiple fastq files in step2
- fix for cluster stats output bug

0.3.24
------
- added nbconvert as a run dependency for the conda build

0.3.23
------
- svd4tet load func improved
- fixed bug with floating point numbers on weights. More speed improvements with fancy matrix tricks.
- added force support to svd4tet
- update releasenotes
- added stats storage to svd4tet
- loci bootstrap sampling implemented in svd4tet
- init_seqarray rearrangement for speed improvement to svd4tet
- removed svd and dstat storage attributes from Assembly Class
- added a plink map output file format for snps locations
- further minimized depth storage in JSON file. Only saved here for a quick summary plot. Full info is in the catg file if needed. Reduces bloat of JSON.
- huge rewrite of svd4tet with Quartet Class Object. Much more concise code
- big rearrangement to svd4tet CLI commands
- code cleanup


0.3.22
------
- only store cluster depth histogram info for bins with data. Removes hugely unnecessary bloat to the JSON file.
- fixed open closure
- massive speed improvement to svd4tet funcs with numba jit compiled C code
- added cores arg to svd4tet

0.3.21
------
- new defaults - lower maxSNPs and higher max_shared_Hs
- massive reworking with numba code for filtering. About 100X speed up.
- reworking numba code in svd4tet for speed
- added debugger to svd4tet
- numba compiling some funcs, and view superseqs as ints instead of strings gives big speedups
- fix to statcounter in demultiplex stats
- improvement to demultiplexing speed
- releasenotes update
- minor fix to advanced tutorial
- updated advanced tutorial
- forgot to rm tpdir when done
- testing s6

0.3.20
------
- bug fix for max_fragment_len errors for paired data and gbs
- fix for gbs data variable cluster sizes.
- prettier printing, does not explicitly say 'saving', but it's still doing it.
- numba update added to conda requirements
- Wrote some numba compiled funcs for speed in step6
- New numba compiled svd func can speed up svd4tet
- update to analysis tools docs

0.3.19
------
- fix for bug in edge trimming when assembly is branched after s6 clustering, but before s7 filtering

0.3.18
------
- Better error handling for alignment step, and now use only the consensus files for the samples being processed (instead of glob'ing every consens.gz in the working directory
- Fix a bug that catches when you don't pass in the -p flag for branching
- cleaning up the releasenotes

0.3.17
------
- removed the -i flag from the command line.
- fix for branching when no filename is provided.
- Fix so that step 6 cleans up as jobs finish. This fixes an error raised if a dummy job finishes too quick. 
- removed a redundant call to open the allhaps file
- Added a check to ensure R2 files _actually exist. Error out if not. Updated internal doc for link_fastq().
- tmp fix for svd4tet test function so we can put up this hotfix

0.3.16
------
- working on speed improvements for svd4tet. Assembly using purging cleanup when running API.
- fix for KeyError caused by cleanup finishing before singlecats in step6
- update to empirical tutorial

0.3.15
------
- write nexus format compatible with ape in svd4tet outputs.
- closing pipe was causing a stall in step6.

0.3.14
------
- merge conflict fix
- set subsample to 2000 high depth clusters. Much faster, minimal decrease in accuracy. Slightly faster code in s4.
- better memory handling. Parallelized better. Starts non-parallel cleanups while singlecats are running = things go faster.
- cluster was commented out in s6 for speed testing

0.3.13
------
- Replaced direct call to  with ipyrad.bins.vsearch
- Fixed reference to old style assembly method reference_sub
- Added ability to optionally pass in a flat file listing subsample names in a column.
- Set a conditional to make sure params file is passed in if doing -b, -r, or -s
- Softened the warning about overlapping barcodes, and added a bit more explanation
- Set default max barcode mismatch to 0

0.3.12
------
- Fixed infinite while loop inside __name_from_file

0.3.11
------
- Fixed commented call to cluster(), step 6 is working again
- Added a check to ensure barcodes contain only IUPAC characters
- Fixed demultiplex sorting progress bar
- append data.name to the tmp-chunks directory to prevent users from running multiple step1 and stepping on themselves
- Update README.rst
- Added force flag for merging CLI
- Bug in rawedit for merged assemblies
- much faster indel entry in step6
- chunks size optimization
- optimizing chunk size step6
- merge for lowmem fixes to step6
- decided against right anchoring method from rad muscle alignments. Improved step6 muscle align progress bar
- reducing memory load in step6
- debug merge fix
- improvement to debug flag. Much improved memory handling for demultiplexing

0.3.10
------
- versioner now actually commits the releasenotes.rst

0.3.9
-----
- Versioner now updates the docs/releasenotes.rst
- Eased back on the language in the performance expectations note
- fixed all links to output formats file
- blank page for recording different performance expectations

0.3.5
-----
- Added `-m` flag to allow merging assemblies in the CLI

0.2.6
-----
- Fix to SNP masking in the h5 data base so that stats counts match the number of snps in the output files. 


0.1.39
------
- Still in development


0.1.38
------
- Still in development. 
- Step7 stats are now building. Extra output files are not. 
- New better launcher for Clients in ipyparallel 5


0.1.37
------
- conda installation mostly working from ipyrad channel

