

.. _release_notes:

Release Notes
=============

0.3.25
------
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- better error message if sample names in barcodes file have spaces in them
- VCF now writes chr ('chromosomes' or 'RAD loci') as ints, since vcftools and other software hate strings apparently
- fix for concatenating multiple fastq files in step2
- fix for cluster stats output bug

0.3.24
------
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- added nbconvert as a run dependency for the conda build

0.3.23
------
- cosmetic code cleanup
- svd4tet load func improved
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
-  "Updating ipyrad/__init__.py to version - 0.3.22


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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
- cosmetic code fix

0.3.18
------
- Better error handling for alignment step, and now use only the consensus files for the samples being processed (instead of glob'ing every consens.gz in the working directory
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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
- Merge branch 'master' of https://github.com/dereneaton/ipyrad
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


