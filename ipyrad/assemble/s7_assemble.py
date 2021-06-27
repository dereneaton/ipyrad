#!/usr/bin/env python

"""
Filter loci and generate output files.
"""

# standard lib imports
import os
import glob
from collections import Counter

# third party imports
from loguru import logger
import numpy as np
import pandas as pd
import h5py
from numba import njit

# internal imports
from ipyrad.core.schema import Project, AssemblyStats, FilterStats, Stats7
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.utils import IPyradError, clustdealer, splitalleles
from ipyrad.assemble.utils import GETCONS, DCONS, chroms2ints

# helper classes ported to separate files.
from ipyrad.assemble.write_outputs_helpers import ChunkProcess
from ipyrad.assemble.write_outputs_converter import Converter
from ipyrad.assemble.write_outputs_vcf import FillVCF, build_vcf

# pylint: disable=too-many-branches, too-many-statements


class Step7(BaseStep):
    """
    Organization for step7 funcs.
    Users can branch assemblies between 6-7, but cannot merge.
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 7, quiet, force)        
        logger.info(f"samples: {sorted(self.samples)}")

        self.data.tmpdir = self.tmpdir
        self.data.stepdir = self.stepdir
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()

        # store step samples to data, and drop reference if not in hackers.
        self.data.samples = self.samples
        if self.data.hackers.exclude_reference and self.data.is_ref:
            self.data.samples.pop("reference")

        # info gathered from step6 results and ncpus
        self.chunksize = 0
        self.clust_database = ""
        self.check_database_files()
        self.get_chunksize()

        # store stats to JSON
        self.proj = Project.parse_file(self.data.json_file)

        # dict mapping of samples to padded names for loci file aligning.
        self.data.pnames, self.data.snppad = self.get_padded_names()

        # output file formats to produce ('l' is required).
        self.outformats = set(['l']).union(
            set(self.data.params.output_formats))

        # make new database files
        self.proj.outfiles = {
            'loci': os.path.join(self.stepdir, f"{self.data.name}.loci")}
        self.proj.outfiles['seqs_database'] = os.path.join(
            self.stepdir, f"{self.data.name}.seqs.hdf5")
        self.proj.outfiles['snps_database'] = os.path.join(
            self.stepdir, f"{self.data.name}.snps.hdf5")
        for base in self.outformats:
            self.proj.outfiles[OUT_SUFFIX[base]] = ""


    def check_database_files(self):
        """
        All samples must have the same clust database file meaning
        that they were run together in step6. No merging of assemblies
        can occur between steps 6-7. This checks that the names in the 
        database file are a superset of the names in self.samples
        """
        msg = (
            "It appears you likely merged assemblies between steps 6 "
            "and 7, which is not allowed. Samples must be merged into "
            "the same assembly before orthologs are identified in step 6. "
            "The simplest solution is to now re-run steps 6 on the merged "
            "assembly with the force=True flag before running step 7."
        )
        dbfiles = [
            self.samples[sname].files.database for sname 
            in self.samples if sname != "reference"
        ]
        if len(set(dbfiles)) > 1:
            raise IPyradError(
                f"Samples have different database files.\n{msg}")
        self.clust_database = dbfiles[0]


    def get_chunksize(self):
        """
        get nloci and ncpus to chunk and distribute work across processors
        """
        # this file is inherited from step 6 to allow step7 branching.
        with open(self.clust_database, 'r') as inloci:
            # skip header
            inloci.readline()
            # get nraw loci
            nraws = sum(1 for i in inloci if i == "//\n") // 2

        # chunk to approximately 4 chunks per core
        ncpus = len(self.ipyclient.ids)
        self.chunksize = sum([
            (nraws // (ncpus * 4)),
            (nraws % (ncpus * 4)),
        ])
        logger.debug(f"using chunksize {self.chunksize} for {nraws} loci")


    def get_padded_names(self):
        """
        Get padded names sorted alphanumeric with 'reference' top top.
        """
        # get longest name
        longname = max(len(i) for i in self.data.samples)
        pnames = {
            i: i.ljust(longname + 5) for i in self.data.samples
        }
        snppad = "//".ljust(longname + 5)
        return pnames, snppad


    def run(self):
        """
        All steps to complete step7 assembly
        """
        # split clusters into bits given nengines.
        self.split_clusters()

        # applies filtering and trimming to the aligned clusters
        # and writes processed chunks to the tmpdir, and returns stats.
        results = self.apply_filters_and_trimming()
        self.collect_stats(results)
        self.write_stats_files()
        self.write_data()
        self.remote_write_outfiles()

        # select additional formats to produce
        # self.store_file_handles()

        # write loci and alleles outputs (parallelized on 3 engines)
        # self.remote_build_arrays_and_write_loci()

        # send conversion jobs from array files to engines

        # send jobs to build vcf
        # throttle job to avoid memory errors based on catg size
        # if 'v' in self.formats:
            # self.remote_fill_depths()
            # self.remote_build_vcf()


    def split_clusters(self):
        """
        Splits the step6 clust_database into chunks to be processed
        in parallel by ChunkProcessor to apply filters.
        """
        with open(self.clust_database, 'rt') as clusters:
            # skip header
            clusters.readline()

            # build iterator
            pairdealer = zip(*[clusters] * 2)

            # grab a chunk of clusters
            idx = 0
            while 1:

                # if an engine is available pull off a chunk
                try:
                    done, chunk = clustdealer(pairdealer, self.chunksize)
                except IndexError as err:
                    msg = "clust_database formatting error in {}".format(chunk)
                    logger.error(msg)
                    raise IPyradError from err

                # write to tmpdir and increment counter
                if chunk:
                    chunkpath = os.path.join(self.tmpdir, f"chunk-{idx}")
                    with open(chunkpath, 'wt') as outfile:
                        outfile.write("//\n//\n".join(chunk))
                    idx += 1

                # break on final chunk
                if done:
                    break


    def apply_filters_and_trimming(self):
        """
        Calls process_chunk() function in parallel.
        """
        def process_chunk(data, chunksize, chunkfile):
            """
            init a ChunkProcessor, run it and collect results.
            """
            # process chunk writes to files and returns proc with features.
            proc = ChunkProcess(data, chunksize, chunkfile)
            proc.run()
            return proc.stats

        # organize chunks and submit to job engine
        chunks = glob.glob(os.path.join(self.tmpdir, "chunk-*"))
        chunks = sorted(chunks, key=lambda x: int(x.rsplit("-")[-1]))        
        jobs = {}
        for chunk in chunks:
            args = (self.data, self.chunksize, chunk)
            jobs[chunk] = self.lbview.apply(process_chunk, *args)
        msg = "applying filters"
        prog = AssemblyProgressBar(jobs, msg, 7, self.quiet)
        prog.block()
        prog.check()
        return prog.results


    def collect_stats(self, results):
        """
        Collect results from ChunkProcess and write stats file.
        """
        # join dictionaries into global stats
        nbases = 0
        locus_covs = Counter({})
        sample_covs = Counter({})
        var_sites = Counter({})
        pis_sites = Counter({})
        var_props = Counter({})
        pis_props = Counter({})
        for chunkfile in results:
            data = results[chunkfile]
            sample_covs.update(data['sample_cov'])
            locus_covs.update(data['locus_cov'])
            var_sites.update(data['var_sites'])
            pis_sites.update(data['pis_sites'])
            var_props.update(data['var_props'])
            pis_props.update(data['pis_props'])
            nbases += data['nbases']

        # reorder site dicts
        var_sites = {i: var_sites[i] for i in sorted(var_sites)}
        pis_sites = {i: pis_sites[i] for i in sorted(pis_sites)}

        # load all of the filters and concatenate
        filters = glob.glob(os.path.join(self.data.tmpdir, "chunk-*.csv"))
        filters = pd.concat([
            pd.read_csv(i, index_col=0) for i in filters],
        ).sum(axis=0)

        # store stats to Project
        self.proj.assembly_stats = AssemblyStats(
            sample_cov=sample_covs,
            locus_cov=locus_covs,
            var_sites=var_sites,
            var_props=var_props,
            pis_sites=pis_sites,
            pis_props=pis_props,
            nloci=sum(locus_covs.values()),
            nsnps=sum([i * var_sites[i] for i in var_sites]),
            nbases=nbases,
            nsamples=len(sample_covs),
            filters=FilterStats(
                nloci_before_filtering=sum(locus_covs.values()) + filters.sum(),
                nloci_after_filtering=sum(locus_covs.values()),
                filtered_by_rm_duplicates=filters.dups,
                filtered_by_min_sample_cov=filters.minsamp,
                filtered_by_max_indels=filters.maxind,
                filtered_by_max_snps=filters.maxvar,
                filtered_by_max_shared_h=filters.maxshared,
            ),
        )

        # store locus stats to Sample objects
        for sname in self.samples:
            self.data.samples[sname].stats_s7 = Stats7(nloci=sample_covs[sname])
            self.data.samples[sname].state = 7

        # write step7 json
        with open(self.data.json_file, 'w') as out:
            out.write(self.proj.json(indent=4, exclude_none=True))


    def write_stats_files(self):
        """
        Write the s7_stats file using results stored in JSON file.
        """
        outstats = open(os.path.join(self.data.stepdir, "s7_assembly_stats.txt"), 'wt')

        if not self.proj.assembly_stats.filters.nloci_after_filtering:
            raise IPyradError("no loci passed filtering")

        outstats.write(
            "## The number of loci before and after filtering, showing the "
            "in which filters were applied.\n\n")

        fdata = pd.DataFrame(
            index=list(self.proj.assembly_stats.filters.dict()),
            columns=["filtered", "retained_loci"],
        )
        fdata.loc["nloci_before_filtering", "filtered"] = 0
        fdata.loc["nloci_before_filtering", "retained_loci"] = (
            self.proj.assembly_stats.filters.nloci_before_filtering)
        fdata.loc["filtered_by_rm_duplicates", "filtered"] = (
            self.proj.assembly_stats.filters.filtered_by_rm_duplicates)
        fdata.loc["filtered_by_rm_duplicates", "retained_loci"] = (
            fdata.iloc[0, 1] - sum(fdata.iloc[:2, 0]))
        fdata.loc["filtered_by_min_sample_cov", "filtered"] = (
            self.proj.assembly_stats.filters.filtered_by_min_sample_cov)
        fdata.loc["filtered_by_min_sample_cov", "retained_loci"] = (
            fdata.iloc[0, 1] - sum(fdata.iloc[:3, 0]))
        fdata.loc["filtered_by_max_indels", "filtered"] = (
            self.proj.assembly_stats.filters.filtered_by_max_indels)
        fdata.loc["filtered_by_max_indels", "retained_loci"] = (
            fdata.iloc[0, 1] - sum(fdata.iloc[:4, 0]))
        fdata.loc["filtered_by_max_snps", "filtered"] = (
            self.proj.assembly_stats.filters.filtered_by_max_snps)
        fdata.loc["filtered_by_max_snps", "retained_loci"] = (
            fdata.iloc[0, 1] - sum(fdata.iloc[:5, 0]))
        fdata.loc["filtered_by_max_shared_h", "filtered"] = (
            self.proj.assembly_stats.filters.filtered_by_max_shared_h)
        fdata.loc["filtered_by_max_shared_h", "retained_loci"] = (
            fdata.iloc[0, 1] - sum(fdata.iloc[:6, 0]))
        fdata.loc["nloci_after_filtering"] = (
            0, self.proj.assembly_stats.filters.nloci_after_filtering)
        fdata.to_string(buf=outstats)
        

        outstats.write(
            "\n\n## The number of loci recovered for each sample.\n\n")
        pd.DataFrame(
            index=list(self.proj.assembly_stats.sample_cov),
            columns=["sample_coverage"],
            data=list(self.proj.assembly_stats.sample_cov.values()),
        ).to_string(buf=outstats)


        outstats.write(
            "\n\n## The number of loci for which N taxa have data.\n"
            "## The 'reference' sample is included if present unless 'exclude_reference=True'\n\n")
        ldata = pd.DataFrame(
            index=list(self.proj.assembly_stats.locus_cov),
            columns=["locus_coverage", "summed_locus_coverage"],
        )
        ldata.locus_coverage = list(self.proj.assembly_stats.locus_cov.values())
        for i in ldata.index:
            ldata.loc[i, 'summed_locus_coverage'] = sum(ldata.locus_coverage[:i])
        ldata.to_string(buf=outstats)

        
        outstats.write(
            "\n\n## The distribution of % variable sites and % pis sites per locus.\n"
            "## loci are binned into 0.1% intervals except the smallest bin for invariant loci.\n"
            "## var = any variable site (pis + autapomorphies).\n"
            "## pis = only parsimony informative variable sites (minor allele in >1 sample).\n\n")
        sdata = pd.DataFrame(
            index=list(self.proj.assembly_stats.var_props),
            columns=["var_sites", "pis_sites"],
        )
        sdata.var_sites = list(self.proj.assembly_stats.var_props.values())
        sdata.pis_sites = list(self.proj.assembly_stats.pis_props.values())
        sdata.to_string(buf=outstats)


        outstats.write(
            "\n\n## The distribution of N variable sites and N pis sites per locus.\n"
            "## var = any variable site (pis + autapomorphies).\n"
            "## pis = only parsimony informative variable sites (minor allele in >1 sample).\n\n")
        maxval = max([
            max(self.proj.assembly_stats.var_sites.keys()),
            max(self.proj.assembly_stats.pis_sites.keys()),
        ])            
        ssdata = pd.DataFrame(
            index=range(maxval + 1),
            columns=["var_sites", "pis_sites"],
        )
        ssdata.var_sites = [
            self.proj.assembly_stats.var_sites[i] 
            if i in self.proj.assembly_stats.var_sites
            else 0 
            for i in range(maxval + 1)
        ]
        ssdata.pis_sites = [
            self.proj.assembly_stats.pis_sites[i] 
            if i in self.proj.assembly_stats.pis_sites
            else 0
            for i in range(maxval + 1)
        ]        
        ssdata.to_string(buf=outstats)

        outstats.write(
            "\n\n## Final sample stats summary\n"
            "## See JSON file or assembly subfolders for detailed stats on each assembly step.\n\n"
        )
        self.data.stats.to_string(buf=outstats)

        outstats.write(
            "\n\n## Alignment matrix statistics:\n"
            "## See ipyrad-analysis toolkit for tools to subsample loci, \n"
            "## SNPs, or matrices with options/filters for missing data.\n\n"
            )
        outstats.close()


    def write_data(self):
        """
        Write the .loci file....
        """
        msg = "writing loci and database files"
        jobs = {0: self.lbview.apply(write_loci, self.data)}
        jobs = {1: self.lbview.apply(fill_seq_array, *(self.data, self.proj))}
        jobs = {2: self.lbview.apply(fill_snp_array, *(self.data, self.proj))}        
        prog = AssemblyProgressBar(jobs, msg, 7, self.quiet)
        prog.block()
        prog.check()


    def store_file_handles(self):
        """
        Fills self.data.outfiles with names of the files to be produced.
        """
        # always produce a .loci file + whatever they ask for.
        testformats = list(self.formats)
        for outf in testformats:

            # if it requires a pop file and they don't have one then skip
            # and write the warning to the expected file, to prevent an
            # annoying message every time if you don't have a pops file, but
            # still to be transparent about skipping some files. This caused
            # me some real amount of pain, like "why isnt' the treemix file
            # being created, fudkckkk!!!1" And then like 10 minutes later, oh
            # yeah, no pops file, fml. 3/2020 iao.
            if (outf in ("t", "m")) and (not self.data.populations):
                outfile = os.path.join(
                    self.data.dirs.outfiles,
                    self.data.name + OUT_SUFFIX[outf][0],
                )
                with open(outfile, 'w') as out:
                    out.write(POPULATION_REQUIRED.format(outf))

                # remove format from the set
                self.formats.discard(outf)
                continue

            # store handle to data object
            for ending in OUT_SUFFIX[outf]:

                # store 
                self.data.outfiles[ending[1:]] = os.path.join(
                    self.data.dirs.outfiles,
                    self.data.name + ending)           


    def remote_write_outfiles(self):
        """
        Calls convert_outputs() in parallel.
        """
        def convert_outputs(data, outf):
            """
            Call the Converter class functions to write formatted output files
            from the HDF5 database inputs.
            """
            try:
                Converter(data).run(outf)
            except Exception as inst:
                # Allow one file to fail without breaking all step 7
                msg = ("Error creating outfile: {}\n{}\t{}"
                    .format(OUT_SUFFIX[outf], type(inst).__name__, inst))
                logger.exception(msg)
                raise IPyradError(msg) from inst

        msg = "writing conversions"
        jobs = {}
        for outf in self.outformats:
            jobs[outf] = self.lbview.apply(convert_outputs, *(self.data, outf))

        # iterate until all chunks are processed
        prog = AssemblyProgressBar(jobs, msg, 7, self.quiet)
        prog.block()
        prog.check()

        # store results to project
        for key in prog.results:
            outfile = prog.results[key]
            outname = OUT_SUFFIX[key]
            self.proj.outfiles[outname] = outfile
        print(self.proj.outfiles)


    def remote_fill_depths(self):
        """
        Call fill_vcf_depths() in parallel.
        """
        printstr = ("indexing vcf depths ", "s7")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasyncs = {}
        for sample in self.data.samples.values():
            if not sample.name == "reference":
                rasyncs[sample.name] = self.lbview.apply(
                    fill_vcf_depths, *(self.data, self.nsnps, sample))
        # iterate until all chunks are processed
        prog.jobs = rasyncs
        prog.block()
        prog.check()


    def remote_build_vcf(self):
        """
        Calls build_vcf() in parallel. 
        """
        printstr = ("writing vcf output  ", "s7")        
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        rasync = self.lbview.apply(build_vcf, self.data)          
        prog.jobs = {0: rasync}
        prog.block()
        prog.check()



###############################################################


def fill_vcf_depths(data, nsnps, sample):
    """
    Get catg depths for this sample.
    """
    filler = FillVCF(data, nsnps, sample)
    filler.run()

    # write vcfd to file and cleanup
    vcfout = os.path.join(data.tmpdir, sample.name + ".depths.hdf5")
    with h5py.File(vcfout, 'w') as io5:
        io5.create_dataset(
            name="depths",
            data=filler.vcfd,
            )
    del filler



# ------------------------------------------------------------
# funcs parallelized on remote engines 
# -------------------------------------------------------------
def write_loci(data):
    """
    Write the .loci file from processed loci chunks. Tries to write
    allele files with phased diploid calls if present.
    """
    # get faidict mapping {int: scaffname}
    if data.is_ref:
        faidict = chroms2ints(data, True)

    # gather all loci bits
    locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
    sortbits = sorted(locibits, 
        key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

    # write to file while adding counters to the ordered loci
    outfile = os.path.join(data.stepdir, f"{data.name}.loci")
    outloci = open(outfile, 'w')

    idx = 0
    for bit in sortbits:
        # store until writing
        lchunk = []

        # LOCI ONLY: iterate through chunk files
        indata = open(bit, 'rt')
        for line in iter(indata):

            # write name, seq pairs
            if "|\n" not in line:
                lchunk.append(line)  # [:5] + line[5:].upper())
                
            # write snpstring and info
            else:
                snpstring, nidxs = line.rsplit("|", 2)[:2]
                if data.params.assembly_method == 'reference':
                    refpos = nidxs.split(",")[0]

                    # translate refpos chrom idx (1-indexed) to chrom name
                    cid, rid = refpos.split(":")
                    cid = faidict[int(cid) - 1]
                    lchunk.append(
                        "{}|{}:{}:{}|\n".format(snpstring, idx, cid, rid))
                else:
                    lchunk.append(
                        "{}|{}|\n".format(snpstring, idx))
                idx += 1
        # close bit handle            
        indata.close()
        outloci.write("".join(lchunk))
    outloci.close()



def pseudoref2ref(pseudoref, ref):
    """
    Reorder psuedoref (observed bases at snps sites) to have the ref allele
    listed first. On rare occasions when ref is 'N' then 
    Called in fill_snps_array.
    """
    # create new empty array
    npseudo = np.zeros(pseudoref.shape, dtype=np.uint8)

    # at all sites where pseudo 0 matches reference, leave it
    matched = np.where(pseudoref[:, 0] == ref)[0]
    npseudo[matched] = pseudoref[matched, :]

    # at other sites, shift order so ref is first
    notmatched = np.where(pseudoref[:, 0] != ref)[0]
    for row in notmatched:
        dat = list(pseudoref[row])

        # skips if ref allele is missing (N)
        try:
            # pop ref and insert it
            new = dat.pop(dat.index(ref[row]))    
            dat.insert(0, new)
            npseudo[row] = dat
        except ValueError:
            npseudo[row] = pseudoref[row]

    return npseudo



def fill_seq_array(data, proj):
    """
    Fills the HDF5 seqs array from loci chunks and stores phymap.
    This contains the full sequence data for all sites >mincov. 
    """
    # gather all loci bits
    locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
    sortbits = sorted(locibits, 
        key=lambda x: int(x.rsplit("-", 1)[-1][:-5]))

    # init/reset hdf5 database
    outfile = os.path.join(data.stepdir, f"{data.name}.seqs.hdf5")
    with h5py.File(outfile, 'w') as io5:

        # temporary array data sets 
        phy = io5.create_dataset(
            name="phy",
            shape=(proj.assembly_stats.nsamples, proj.assembly_stats.nbases), 
            dtype=np.uint8,
        )
        # temporary array data sets 
        phymap = io5.create_dataset(
            name="phymap",
            shape=(proj.assembly_stats.nloci, 5),
            dtype=np.uint64,
        )

        # alphanumeric name order except 'reference' on top.
        snames = sorted(data.samples)
        if "reference" in snames:
            snames.remove("reference")
            snames = ["reference"] + snames
        sidxs = {sample: i for (i, sample) in enumerate(snames)}

        # store attrs of the reference genome to the phymap
        if data.params.assembly_method == 'reference':
            io5["scaffold_lengths"] = get_fai_values(data, "length")
            io5["scaffold_names"] = get_fai_values(data, "scaffold").astype("S")
            phymap.attrs["reference"] = data.params.reference_sequence
        else:
            phymap.attrs["reference"] = "pseudoref"

        # store names and 
        phymap.attrs["phynames"] = snames
        phymap.attrs["columns"] = [
            b"chroms", b"phy0", b"phy1", b"pos0", b"pos1", 
        ]

        # skip the first row when finding invariant sites because reference
        # is present and always has data.
        skip_row = False
        if data.is_ref and data.is_pair and (not data.hackers.exclude_reference):
            skip_row = True

        # iterate through file
        gstart = 0
        start = end = 0
        maxsize = 100000
        tmploc = {}
        mapends = []
        mapchroms = []
        mappos0 = []
        mappos1 = []
        mapstart = 0
        locidx = 0

        # array to store until writing; TODO: Accomodate large files...
        tmparr = np.zeros((proj.assembly_stats.nsamples, maxsize + 50000), dtype=np.uint8)
        
        # iterate over chunkfiles
        for bit in sortbits:
            # iterate lines of file until locus endings
            indata = open(bit, 'rt')
            for line in iter(indata):
                
                # still filling locus until |\n
                if "|\n" not in line:

                    # if empty skip
                    try:
                        name, seq = line.split()
                        tmploc[name] = seq
                    except ValueError:
                        continue

                # locus is full, dump it
                else:
                    # convert seqs to an array
                    locidx += 1

                    # parse chrom:pos-pos                   
                    if data.is_ref:
                        lineend = line.split("|")[1]
                        chrom = int(lineend.split(":")[0])
                        pos0, pos1 = 0, 0                    
                        pos0, pos1 = (
                            int(i) for i in lineend
                            .split(":")[1]
                            .split(",")[0]
                            .split("-")
                        )

                    # seq ordered into array by snames as int8 (py2/3 checked)
                    loc = np.array([
                        list(tmploc[i].encode()) for i in snames if i in tmploc
                    ]).astype(np.int8)
                    
                    # drop the site that are all N or - (e.g., pair inserts)
                    if skip_row:
                        mask = np.all(loc[1:, :] == 78, axis=0)
                    else:
                        mask = np.all((loc == 45) | (loc == 78), axis=0)
                    loc = loc[:, np.invert(mask)]
                    
                    # store end position of locus for map
                    end = start + loc.shape[1]

                    # checked for py2/3 (keeping name order straight important)
                    lidx = 0
                    for name in snames:
                        if name in tmploc:
                            sidx = sidxs[name]
                            tmparr[sidx, start:end] = loc[lidx]
                            lidx += 1
                    mapends.append(gstart + end)

                    if data.is_ref:
                        mapchroms.append(chrom)
                        mappos0.append(pos0)
                        mappos1.append(pos1)
                    else:
                        mapchroms.append(locidx - 1)
                    
                    # reset locus
                    start = end
                    tmploc = {}
                    
                # dump tmparr when it gets large
                if end > maxsize:
                    
                    # trim right overflow from tmparr (end filled as 0s)
                    trim = np.where(tmparr != 0)[1]
                    if trim.size:
                        trim = trim.max() + 1
                    else:
                        trim = tmparr.shape[1]

                    # fill missing with 78 (N)
                    tmparr[tmparr == 0] = 78
                    
                    # dump tmparr to hdf5
                    phy[:, gstart:gstart + trim] = tmparr[:, :trim]                       
                    phymap[mapstart:locidx, 0] = mapchroms
                    phymap[mapstart:locidx, 2] = mapends
                    if data.is_ref:
                        phymap[mapstart:locidx, 3] = mappos0
                        phymap[mapstart:locidx, 4] = mappos1                       
                    mapstart = locidx
                    mapends = []
                    mapchroms = []
                    mappos0 = []
                    mappos1 = []
                    
                    # reset
                    tmparr = np.zeros(
                        (proj.assembly_stats.nsamples, maxsize + 50000),
                        dtype=np.uint8,
                    )
                    gstart += trim
                    start = end = 0

            # close bit handle
            indata.close()

        if start == 0 and end == 0:
            # The last chunk fell exactly on the maxsize boundary so it has
            # already been trimmed and dumped to the phy. In this case the for
            # loop on the line iterator has fallen through (no more data) so
            # there is no final chunk to trim.
            pass
        else:
            # trim final chunk tmparr to size
            trim = np.where(tmparr != 0)[1]
            if trim.size:
                trim = trim.max() + 1
            else:
                trim = tmparr.shape[1]

            # fill missing with 78 (N)
            tmparr[tmparr == 0] = 78

            # dump tmparr and maplist to hdf5. Because we dropped sites that are
            # all N or - the length of phy data can be less than nbases and so
            # there can be 000 at the end of phy. This is ok, we trim it when
            # writing phylip file, but good to be aware it's there for other things
            phy[:, gstart:gstart + trim] = tmparr[:, :trim]
            phymap[mapstart:locidx, 0] = mapchroms
            phymap[mapstart:locidx, 2] = mapends
            if data.is_ref:
                phymap[mapstart:locidx, 3] = mappos0
                phymap[mapstart:locidx, 4] = mappos1
        phymap[1:, 1] = phymap[:-1, 2]

        # fill 'scaffold' information for denovo data sets from the data
        if data.params.assembly_method != "reference": 
            # 1-index the "chromosomes" and store lengths as loclens
            phymap[:, 0] += 1  # <- does not cause problems for outfiles...
            io5["scaffold_names"] = (io5["phymap"][:, 0]).astype("S")
            io5["scaffold_lengths"] = io5["phymap"][:, 2] - io5["phymap"][:, 1]

        # write stats to the output file
        outstats = os.path.join(data.stepdir, "s7_assembly_stats.txt")
        with open(outstats, 'a') as outstats:
            trim = phymap[-1, 2]
            missmask = phy[:trim] == 78
            missmask += phy[:trim] == 45
            missing = 100 * (missmask.sum() / float(phy[:trim].size))
            outstr = (
                "sequence matrix size: ({}, {}), {:.2f}% missing sites.\n"
                .format(
                    len(snames),
                    trim,
                    max(0, missing),
                )
            )
            outstats.write(outstr)



def fill_snp_array(data, proj):
    """
    Fills the SNPS HDF5 database with variables sites.
    """
    # gather all loci bits
    locibits = glob.glob(os.path.join(data.tmpdir, "*.loci"))
    sortbits = sorted(
        locibits, 
        key=lambda x: int(x.rsplit("-", 1)[-1][:-5])
    )

    # open new database file handle
    outfile = os.path.join(data.stepdir, f"{data.name}.snps.hdf5")
    with h5py.File(outfile, 'w') as io5:

        # Database files for storing arrays on disk. 
        # Should optimize for slicing by rows if we run into slow writing, or 
        # it uses too much mem. For now letting h5py to auto-chunking.
        io5.create_dataset(
            name="snps",
            shape=(proj.assembly_stats.nsamples, proj.assembly_stats.nsnps),
            dtype=np.uint8,
        )
        # store snp locations:
        # (loc-counter, loc-snp-counter, loc-snp-pos, chrom, chrom-snp-pos)
        io5.create_dataset(
            name="snpsmap",
            shape=(proj.assembly_stats.nsnps, 5),
            dtype=np.uint32,
        )
        # store snp locations
        io5.create_dataset(
            name="pseudoref",
            shape=(proj.assembly_stats.nsnps, 4),
            dtype=np.uint8,
        )
        # store genotype calls (0/0, 0/1, 0/2, etc.)
        io5.create_dataset(
            name="genos",
            shape=(proj.assembly_stats.nsnps, proj.assembly_stats.nsamples, 2),
            dtype=np.uint8,
        )

        # alphanumeric name order except 'reference' on top.
        snames = sorted(data.samples)
        if "reference" in snames:
            snames.remove("reference")
            snames = ["reference"] + snames
        sidxs = {sample: i for (i, sample) in enumerate(snames)}

        # store sample names and snpmap columns names as attributes
        io5["snps"].attrs["names"] = [i.encode() for i in snames]
        io5["snpsmap"].attrs["columns"] = [
            b"locus", b"locidx", b"locpos", b"scaf", b"scafpos",
        ]

        # iterate through file
        start = end = 0
        tmploc = {}
        locidx = 1
        snpidx = 1            

        # array to store until writing
        tmparr = np.zeros(
            shape=(proj.assembly_stats.nsamples, proj.assembly_stats.nsnps), 
            dtype=np.uint8,
        )
        tmpmap = np.zeros((proj.assembly_stats.nsnps, 5), dtype=np.uint32)

        # iterate over chunkfiles
        for bit in sortbits:
            # iterate lines of file until locus endings
            indata = open(bit, 'rt')
            for line in iter(indata):

                # while still filling locus until |\n store name,seq in dict
                if "|\n" not in line:
                    try:
                        name, seq = line.split()
                        tmploc[name] = seq
                    except ValueError:
                        continue

                # locus is full, dump it
                else:
                    # convert seqs to an np.int8 array, checked py2/3
                    loc = np.array(
                        [list(tmploc[i].encode()) for i in snames if i in tmploc]
                    ).astype(np.int8)
                    snps, idxs, _ = line[len(data.snppad):].rsplit("|", 2)
                    snpsmask = np.array(list(snps)) != " "
                    snpsidx = np.where(snpsmask)[0]

                    # select only the SNP sites
                    snpsites = loc[:, snpsmask]

                    # store end position of locus for map
                    end = start + snpsites.shape[1]

                    # checked for py2/3 (keeping name order straight important)
                    lidx = 0
                    for name in snames:
                        if name in tmploc:
                            sidx = sidxs[name]
                            tmparr[sidx, start:end] = snpsites[lidx, :]
                            lidx += 1

                    # store snpsmap data 1-indexed with chroms info
                    if data.is_ref:
                        chrom, pos = idxs.split(",")[0].split(":")
                        start = int(pos.split("-")[0])
                        chromidx = int(chrom)
                        for isnp in range(snpsites.shape[1]):
                            isnpx = snpsidx[isnp]
                            tmpmap[snpidx - 1] = (
                                locidx, isnp, isnpx, chromidx, isnpx + start,
                            )
                            snpidx += 1

                    # store snpsmap data (snpidx is 1-indexed)
                    else:
                        for isnp in range(snpsites.shape[1]):
                            tmpmap[snpidx - 1] = (
                                locidx, isnp, snpsidx[isnp], 0, snpidx,
                            )
                            snpidx += 1
                    locidx += 1

                    # reset locus
                    start = end
                    tmploc = {}

            # close file handle
            indata.close()

        # fill missing with 78 (N)
        tmparr[tmparr == 0] = 78

        # dump tmparr and maplist to hdf5       
        io5['snps'][:] = tmparr[:]
        io5['snpsmap'][:] = tmpmap
        del tmparr

        # write stats output
        nsnps = proj.assembly_stats.nsnps
        outstats = os.path.join(data.stepdir, "s7_assembly_stats.txt")
        with open(outstats, 'a') as outstats:
            missmask = io5["snps"][:] == 78
            missmask += io5["snps"][:] == 45
            missing = 100 * (missmask.sum() / float(io5["snps"][:nsnps].size))
            outstr = (
                "snps matrix size: ({}, {}), {:.2f}% missing sites.\n"
                .format(
                    len(snames),
                    nsnps,
                    missing,
                )
            )
            outstats.write(outstr)

        # fill in the reference and geno arrays
        # convert snps to characters uppered to get most common as pseudoref
        snparr = io5["snps"][:].view("S1")
        snparr = np.char.upper(snparr).view(np.uint8)

        # store pseudo-ref (most common base)
        # with ambiguous bases resolved: (87, 78, 0, 0).
        if data.params.assembly_method != 'reference':
            io5['pseudoref'][:] = reftrick(snparr, GETCONS)

        else:
            ref = snparr[snames.index('reference')]   
            pseudoref = reftrick(snparr, GETCONS)
            io5['pseudoref'][:] = pseudoref2ref(pseudoref, ref)

        # fill for each taxon
        for sidx, _ in enumerate(snames):
            resos = [DCONS[i] for i in snparr[sidx, :]]

            # pseudoref version
            io5['genos'][:, sidx, :] = get_genos(
                np.array([i[0] for i in resos]), 
                np.array([i[1] for i in resos]),
                io5['pseudoref'][:]
            )


###############################################################


@njit
def reftrick(iseq, consdict):
    """
    Returns the most common base at each site in order.
    Called from fill_snp_array().
    """
    altrefs = np.zeros((iseq.shape[1], 4), dtype=np.uint8)
    altrefs[:, 1] = 46

    for col in range(iseq.shape[1]):
        # expand colums with ambigs and remove N-
        fcounts = np.zeros(111, dtype=np.int64)
        counts = np.bincount(iseq[:, col])  #, minlength=90)
        fcounts[:counts.shape[0]] = counts
        # set N and - to zero, wish numba supported minlen arg
        fcounts[78] = 0
        fcounts[45] = 0
        # add ambig counts to true bases
        for aidx in range(consdict.shape[0]):
            nbases = fcounts[consdict[aidx, 0]]
            for _ in range(nbases):
                fcounts[consdict[aidx, 1]] += 1
                fcounts[consdict[aidx, 2]] += 1
            fcounts[consdict[aidx, 0]] = 0

        # now get counts from the modified counts arr
        who = np.argmax(fcounts)
        altrefs[col, 0] = who
        fcounts[who] = 0

        # if an alt allele fill over the "." placeholder
        who = np.argmax(fcounts)
        if who:
            altrefs[col, 1] = who
            fcounts[who] = 0

            # if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 2] = who
            fcounts[who] = 0

            # if 3rd or 4th alleles observed then add to arr
            who = np.argmax(fcounts)
            altrefs[col, 3] = who

    return altrefs


@njit
def get_genos(f10, f01, pseudoref):
    """
    Returns genotype as 0/1/2/3 or 9 for missing.
    Called from fill_snp_array
    """
    res = np.zeros((f10.size, 2), dtype=np.uint8)
    for i in range(f10.size):
        match = np.where(f10[i] == pseudoref[i])[0]
        if match.size:
            res[i, 0] = match[0]
        else:
            res[i, 0] = 9
        match = np.where(f01[i] == pseudoref[i])[0]
        if match.size:
            res[i, 1] = match[0]
        else:
            res[i, 1] = 9
    return res



def get_fai_values(data, value):
    """
    Returns ... from the indexed reference genome as an array
    """
    reference_file = data.params.reference_sequence
    fai = pd.read_csv(
        reference_file + ".fai",   
        names=['scaffold', 'length', 'sumsize', 'a', 'b'],
        sep="\t",
    )
    return fai[value].values  


BADPOP_SAMPLES = """
There are sample names in the populations assignments that are not present in 
this assembly. This is likely due to a typo and should be corrected. The 
following sample names are in the pop assignments but not in this Assembly:
{}
"""
POPULATION_REQUIRED = """\
Warning: Skipping output format '{}'. Requires population assignments.
You can alternatively create this type of file using ipyrad-analysis 
after assembling your data.
"""

OUT_SUFFIX = {
    'l': ('loci',),
    'p': ('phy',),
    's': ('snps', 'snpsmap',),
    'n': ('nex',),
    'k': ('str',),
    # 'a': ('alleles',),
    'g': ('geno',),
    'G': ('gphocs',),
    'u': ('usnps', 'ustr', 'ugeno'),
    'v': ('vcf',),
    't': ('treemix',),
    'm': ('migrate',),
}


if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG", logfile="/tmp/test.log")

    TEST = ip.load_json("/tmp/TEST1.json")
    TEST.run("7", force=True, quiet=False)    

    # tdata = ip.load_json("/tmp/test-simpairddrad.json")
    # tdata.params.output_formats = "lpsnkaguvtm"
    # tdata.run("7", auto=True, force=True)
    # logger.info(tdata.stats.T)

    # tdata = ip.load_json("/tmp/test-amaranth.json")
    # tdata.run("7", auto=True, force=True)
    # print(tdata.stats)
    # print(tdata.stats_dfs.s5)

    # self.data.hackersonly.declone_PCR_duplicates:
    # tdata = ip.load_json("/tmp/test-amaranth-denovo.json")
    # tdata.ipcluster['cores'] = 4
    # tdata.run("7", auto=True, force=True)
    # logger.info(tdata.stats.T)


    ## compare .loci file to extracted loci from HDF5.
    # with h5py.File("./refdata_outfiles/refdata.seqs.hdf5", "r") as io5:
    #     names = (io5["phymap"].attrs['phynames'])
    #     scaffidx, start, end, pos0, pos1 = io5["phymap"][2]
    #     print(scaffidx, start, end, pos0, pos1)
    #     arr = io5["phy"][:, start:end]
    #     for idx in range(arr.shape[0]):
    #         print(
    #           f"{names[idx]}\t", 
    #           b"".join(arr[idx, :10].view("S1")).decode(), 
    #           b"".join(arr[idx, -10:].view("S1")).decode())