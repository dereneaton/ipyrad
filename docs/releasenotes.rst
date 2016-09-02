

.. _release_notes:

Release Notes
=============

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


