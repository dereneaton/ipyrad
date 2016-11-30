

.. _release_notes:

Release Notes
=============

0.5.9
-----

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
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6

0.5.6
-----
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6

0.5.6
-----
-  "Updating ipyrad/__init__.py to version - 0.5.6
-  "Updating ipyrad/__init__.py to version - 0.5.6

0.5.6
-----
-  "Updating ipyrad/__init__.py to version - 0.5.6

0.5.6
-----

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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- Better support for 3rad lining presorted fastqs.
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- bucky cookbook updated
- docs update
- docs update
- dstat updates
- docs update
- docs update
- docs update
- docs update
- docs update
- docs updates
- docs update
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- cosmetic
- default to no subsampling in jointestimate call
- testing bugfixes to jointestimate
- added hackersonly option for additional adapters to be filtered
- bug fix to joint H,E estimate for large data sets introduced in v.0.3.14 that was yielding inflated rates.
- fix for core count when using API
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- Added plots of snp depth across loci, as well as loci counts per sample to results notebook.
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- phylogenetic_invariants notebook up
- some notes on output formats plans
- removed leftjust arg b/c unnecessary and doesn't work well with left trimmed data

0.4.2
-----
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- Merging for Samples at any state, with warning for higher level states. Prettier printing for API. Fix to default cores setting on API.
- fix for merged Assemblies/Samples for s2
- fix for merged Assemblies&Samples in s3
- removed limit on number of engines used during indexing
- Added ddocent to manuscript analysis.
- tutorial update
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- in progress doc notebook
- parallel waits for all engines when engines are designated, up until timeout limit
- parallelized loading demux files, added threads to _ipcluster dict, removed print statement from save
- vcf header was missing
- added step number to progress bar when in interactive mode
- added warning message when filter=2 and no barcodes are present
- improved kill switch in step 1
- cosmetic changes
- cosmetic changes
- use select to improve cluster progress bar
- added a CLI option to fine-tune threading
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- added dstat storage by default
- new default trim_overhang setting and function (0,0,0,0)
- fix for overzealous warning message on demultiplexing when allowing differences

0.4.1
-----
- Fixed reference before assignment error in step 2.

0.4.0
-----
- Cosmetic change
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- new sim data and notebook up
- Added aftrRAD to the manuscript analysis horserace
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- made merging reads compatible with gzipped files from step2
- modify help message
- made TESTS global var, made maparr bug fix to work with no map info
- More carefully save state after completion of each step.
- limit vsearch merging to 2 threads to improve parallel, but should eventually make match to cluster threading. Added removal of temp ungzipped files.
- more detailed Sample stats_df.s2 categories for paired data
- made merge command compatible with gzip outputs from step2
- simplified cutadapt code calls
- cosmetic
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- Merge branch 'master' of https://github.com/StuntsPT/ipyrad into StuntsPT-master
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
- cosmetic changes to merge function in util
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- hotfix for memory error in build_clusters, need to improve efficiency for super large numbers of hits
- more speed testing on tetrad
- merge conflict
- cosmetic
- cosmetic
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- cosmetic changes
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
- cosmetic code cleanup
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
- cosmetic
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
- cosmetic code fix

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
- just cosmetic code cleanup.
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
- cosmetic changes
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
- cosmetic changes
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


