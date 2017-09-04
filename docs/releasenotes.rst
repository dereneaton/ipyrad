

.. _release_notes:

Release Notes
=============

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


