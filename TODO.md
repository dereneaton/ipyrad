
# Main CLI workflow
```bash
$ ipyrad -h
$ ipyrad demux -h
$ ipyrad new -h
$ ipyrad assemble -h
$ ipyrad download -h
$ ipyrad merge -h
$ ipyrad branch -h
$ ipyrad stats -h
```

# Checklist
- [x] create samples w/ .files attr and infer restriction overhang and barcodes.
- [x] perform fastp filtering on samples.
- [x] save stats_s1 results to new JSON.
- [x] param option to keep restriction overhang on sequences? Decided No, b/c it can be variable...
- [x] input: [data1/*.gz data2/*gz] or data1/.gz
- [x] out: trimmed/trimmed

### demux
- [x] .demux update for pandas deprecation of delim_whitespace arg.
- [ ] catch memory errors and log report before dying... Seems difficult.
- [ ] do not report zero reads found when merging technical replicates.

### merge
- [ ] --merge-samples option which drops those samples back to state=1
- [ ] only allow merge samples when they are in state = 1.
- [ ] merge at any step, but walk back sample states and log report it.
- [ ] only allow merge assemblies when...
	- 1 (before clustering)
	- 3 (after cluster building, before consens calling)
	- 5 (before clustering across)
	- The only time to merge samples with the same name from diff assemblies
	  into a single sample with data from both is between steps 1-2. Any other
	  time and an exception will be raised saying you cannot combine the samples.

### step 1 (trim/)
- [x] store .files.fastqs as List[Tuple[Path, Path | str]]
- [x] store .files.trimmed as Tuple[Path, Path]
- [x] allow merging samples at this stage?

### step 2 (within_clust/, within_clust_build/, within_refmap/)
- [x] what is diff between derep-fulllength and fastqx-uniques?: depr to support fasta/q
- [ ] test reference-minus filtering.
- [x] check threading of paired merging: it is 2 hard-coded.
- [x] perform dereplicating
- [ ] sort samples by raw reads
- [x] log report on each finished.
- [x] try using centroids
- [x] perform clustering.
- [ ] perform mapping.
- [ ] update base step message
- [x] add more logging info on who is being clustered/finished.
- [x] .derep should be \_derep...
- [ ] vsearch --fastq_join option w/ --join_padgap nnnn
- [ ] vsearch merge test w/ option --fastq_nostagger instead of allowmergestagger. Stagger are then discarded. Maybe conservative is better, since stagger part should already be trimmed. Examine diff in real data.
- [ ] increase default clust threshold? Run on ipsimdata.
- [x] in: trimmed list
- [x] out: (derep,matches)
- [ ] out: (bam, None)
- [ ] out: store number of unaligned clusters?

### step 3 
- [x] build clusters from clust files
- [x] perform alignments of clusters
- [x] max_size of alignments
- [ ] deprecate max_depth? is too low! use a multiplier? 
- [ ] consider not allowing left of seed start.
- [ ] in sublime KBD kill is still not catching... try in CLI. Why so hard to stop muscle?
- [ ] write stats.

### step 4
- [ ] ref: use pileup to call variants?
- [ ] denovo: infer H,E and call consens alleles
- [ ] keep sumdepth info in the header.

### step 5
- [ ] cluster/map across samples
- [ ] rm duplicates? options: drop lowest size, or, use NJ to split into separate.

### step 6
- [ ] build-align unfiltered clusters. Very fast for ref.

### step 7
- [ ] filtering and formatting.

### Analysis
- [ ] plink imputation for REF data.
- [ ] allow ipa.bpp to accept no tree arg.


## Multilane example --------------------------------------
ipyrad demux -o lib1 -...
ipyrad demux -o lib2 -...

ipyrad new -n lib1 --fastq_path lib1/
ipyrad new -n lib2 --fastq_path lib2/ --reference_sequence ?

<!-- step1 of assembly EXPECTS 2 files per sample, reports stats on each. -->
ipyrad assemble -p params-lib1 -s 1
ipyrad assemble -p params-lib2 -s 1 

<!-- before step2 is the only chance to merge SAMPLES -->
ipyrad merge -p params-lib1 params-lib2 -n libs12

<!-- step2-3 will cluster ... -->
ipyrad assemble -p params-libs12 -s 23

<!-- later steps can branch and/or merge samples into subsets for speed or param testing-->
ipyrad assemble -p params-libs12 -s 456

<!-- branch to try different parameters -->
ipyrad assemble -p params-libs12 -s 7 



### step 3
- step3 reference examine clusters with different trim_edges settings...
- refminus mapping strips i5 and so can't currently be combined with decloning.
- collapse same info derep'd reads (differ only by N or len) to speed align while retaining i5 tags?

### step 56
- debugging still in `consens_utils._iter_build_consens` for denovo pair to trim N around insert.
- persist alleles?
- store nalleles?
- do we still need 'concat_denovo_consen()' from consens_utils?

### step 6
- merging duplicates for denovo? or drop one?

### step 7
[x] rename reference to `!reference*` and do not allow sample names to start with ! or ' ' ... so we can simply sort names. Can't do this, bad name char. Easier to just special sort.
[ ] to_locus in s7.write_outputs_processor for reference is fully borken
[ ] snpsmap.tsv for reference data distinguish loc from scaff
[ ] use INFO to logger when using population information.
- Continue on denovo depths script for combining vcf depth info.
	- test with more trim and indels
	- test on empirical rad denovo
	- test on sim pairddrad denovo
	- test on sim pairgbs denovo
	- test on empirical pair3rad
[ ] develop alternative depths grepper for reference data.
[ ] save version to VCF output and stats w/o import error. maybe use package distr...

### IPA
- PCA: no imap crashes drawing
- PCA: imap but no minmap complains of missing minmap
- Coverage plot (look in isolation dir on oud)

### SRA
- update to sra-tools 3.0 with user download/install and optional provide binary instructions.

### CLI
- test merge
- logger: default WARNING to FILE. Options for logger to stdout w/ print suppressed.
- JSON & HDF5 and VCF2HDF5 CONVERTERS.
- Instead of error on dir exists and no-force, should we just skip to next step?

### DENOVO
- test with reference as filter to remove reads.


# Empirical test datasets

- denovo 
	- single RAD
		- ...
	- single GBS
		- ...
	- paired GBS
		- ...
	- paired 3RAD

- reference
	- single RAD
		- ...
	- single GBS
		- ...
	- paired GBS
		- ...
	- paired 3RAD
		- ...
- test VCF correctness.
