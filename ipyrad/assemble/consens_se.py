#!/usr/bin/env python

"""
Call consensus base calls on paired or single-end stacks/contigs
"""

# py2/3 compatible
from __future__ import print_function
from builtins import range
try:
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip

import os
import time
import gzip
import glob
import shutil
import warnings
import subprocess as sps
from collections import Counter

from loguru import logger
import numpy as np
import pandas as pd
import scipy.stats

import ipyrad as ip
from ipyrad.assemble.jointestimate import recal_hidepth
from ipyrad.assemble.utils import IPyradError, clustdealer, PRIORITY
from ipyrad.assemble.utils import AssemblyProgressBar

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py



TRANS = {
    (71, 65): 82,
    (71, 84): 75,
    (71, 67): 83,
    (84, 67): 89,
    (84, 65): 87,
    (67, 65): 77,
    (65, 67): 77,
    (65, 84): 87,
    (67, 84): 89,
    (67, 71): 83,
    (84, 71): 75,
    (65, 71): 82,
}


DCONS = {
    82: (71, 65),
    75: (71, 84),
    83: (71, 67),
    89: (84, 67),
    87: (84, 65),
    77: (67, 65),
    78: (9, 9),
    45: (9, 9),
    67: (67, 67),
    65: (65, 65),
    84: (84, 84),
    71: (71, 71),
}


class Step5:
    """
    Organized Step 5 functions for all datatype and methods
    """
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.print_headers()
        self.samples = self.get_subsamples()
        self.isref = bool("reference" in data.params.assembly_method)
        self.ipyclient = ipyclient
        self.lbview = ipyclient.load_balanced_view()
        self.setup_dirs()


    def print_headers(self):
        """
        print headers for the CLI
        """
        if self.data._cli:
            self.data._print(
                "\n{}Step 5: Consensus base/allele calling "
                .format(self.data._spacer)
            )



    def get_subsamples(self):
        """
        Apply state, ncluster, and force filters to select samples
        """
        # bail out if no samples ready
        if not hasattr(self.data.stats, "state"):
            raise IPyradError("No samples ready for step 5")

        # filter samples by state
        state3 = self.data.stats.index[self.data.stats.state < 4]
        state4 = self.data.stats.index[self.data.stats.state == 4]
        state5 = self.data.stats.index[self.data.stats.state > 4]

        # tell user which samples are not ready for step5
        if state3.any():
            print("skipping samples not in state==4:\n{}"
                  .format(state3.tolist()))

        if self.force:
            # run all samples above state 3
            subs = self.data.stats.index[self.data.stats.state > 3]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 5
            if state5.any():
                print("skipping samples already finished step 5:\n{}"
                      .format(state5.tolist()))

            # run all samples in state 4
            subsamples = [self.data.samples[i] for i in state4]

        # check that kept samples have clusters
        checked_samples = []
        for sample in subsamples:
            if sample.stats.clusters_hidepth:
                checked_samples.append(sample)
            else:
                print("skipping {}; no clusters found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 5")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.clusters_hidepth,
            reverse=True,
        )

        # if sample is already done skip
        if "hetero_est" not in self.data.stats:
            for sample in checked_samples:
                sample.stats.hetero_est = 0.001
                sample.stats.error_est = 0.0001

        if self.data._cli:
            print(u"  Mean error  [{:.5f} sd={:.5f}]".format(
                self.data.stats.error_est.mean(),
                self.data.stats.error_est.std(),
            ))
            print(u"  Mean hetero [{:.5f} sd={:.5f}]".format(
                self.data.stats.hetero_est.mean(),
                self.data.stats.hetero_est.std(),
            ))
        return checked_samples



    def setup_dirs(self):
        """
        setup directories, remove old tmp files
        """
        # final results dir
        self.data.dirs.consens = os.path.join(
            self.data.dirs.project, 
            "{}_consens".format(self.data.name))
        if not os.path.exists(self.data.dirs.consens):
            os.mkdir(self.data.dirs.consens)

        # tmpfile dir (zap it if it exists)
        self.data.tmpdir = os.path.join(
            self.data.dirs.project, 
            "{}-tmpdir".format(self.data.name))
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.mkdir(self.data.tmpdir)

        # assign output file handles for s6
        for sample in self.samples:
            if not self.isref:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.gz".format(sample.name))
            else:
                sample.files.consens = os.path.join(
                    self.data.dirs.consens, 
                    "{}.consens.bam".format(sample.name))

        # zap existing consens files for selected samples
        for sample in self.samples:
            if os.path.exists(sample.files.consens):
                os.remove(sample.files.consens)
            dfile = sample.files.consens.rsplit(".", 2)[0]
            sample.files.database = dfile + ".catg.hdf5"
            if os.path.exists(dfile):
                os.remove(dfile)

        # set up parallel client: allow user to throttle cpus
        self.lbview = self.ipyclient.load_balanced_view()
        if self.data.ipcluster["cores"]:
            self.ncpus = self.data.ipcluster["cores"]
        else:
            self.ncpus = len(self.ipyclient.ids)



    def run(self):
        """
        Run the main functions on the parallel client
        """
        self.remote_calculate_depths()
        self.remote_make_chunks()
        statsdicts = self.remote_process_chunks()
        self.remote_concatenate_chunks()
        self.data_store(statsdicts)
        self.data.save()



    def remote_calculate_depths(self):
        """
        Checks whether mindepth has changed and calc nclusters and maxlen
        """
        # send jobs to be processed on engines
        start = time.time()
        printstr = ("calculating depths  ", "s5")
        prog = AssemblyProgressBar({}, start, printstr, self.data)
        prog.update()

        jobs = {}
        maxlens = []
        for sample in self.samples:
            rasync = self.lbview.apply(recal_hidepth, *(self.data, sample))
            jobs[sample.name] = rasync

        # block until finished
        prog.jobs = jobs
        prog.block()
        prog.check()

        # check for failures and collect results
        for sample in self.samples:
            hidepth, maxlen, _, _ = jobs[sample.name].get()
            sample.stats["clusters_hidepth"] = hidepth
            sample.stats_dfs.s3["clusters_hidepth"] = hidepth
            maxlens.append(maxlen)
        
        # update hackersdict with max fragement length
        self.data.hackersonly.max_fragment_length = max(maxlens)
        logger.debug("max_fragment_length set to {}".format(max(maxlens)))



    def remote_make_chunks(self):
        """
        split clusters into chunks for parallel processing
        """
        printstr = ("chunking clusters   ", "s5")
        prog = AssemblyProgressBar({}, None, printstr, self.data)
        prog.update()

        # send off samples to be chunked
        jobs = {}
        for sample in self.samples:
            args = (self.data, sample, len(self.ipyclient))
            rasync = self.lbview.apply(make_chunks, *args)
            jobs[sample.name] = rasync

        # block until finished
        prog.jobs = jobs
        prog.block()
        prog.check()



    def remote_process_chunks(self):
        """
        Process the cluster chunks into arrays and consens or bam files
        """
        # send chunks to be processed
        start = time.time()
        printstr = ("consens calling     ", "s5")
        prog = AssemblyProgressBar({}, start, printstr, self.data)
        prog.update()

        # submit jobs (10 per sample === can be hundreds of jobs...)
        jobs = {sample.name: [] for sample in self.samples}
        for sample in self.samples:
            chunks = glob.glob(os.path.join(
                self.data.tmpdir,
                "{}.chunk.*".format(sample.name)))
            chunks.sort(key=lambda x: int(x.split('.')[-1]))

            # submit jobs
            for chunk in chunks:
                args = (self.data, sample, chunk, self.isref)
                rasync = self.lbview.apply(process_chunks, *args)
                jobs[sample.name].append(rasync)
               
        # track progress - just wait for all to finish before concat'ing
        allsyncs = chain(*[jobs[i] for i in jobs])
        prog.jobs = dict(enumerate(allsyncs))
        prog.update()
        prog.block()
        prog.check()

        # collect all results for a sample and store stats 
        statsdicts = {}
        for sample in self.samples:
            statsdicts[sample.name] = [i.get() for i in jobs[sample.name]]
        return statsdicts



    def remote_concatenate_chunks(self):
        """
        Concatenate chunks and relabel for joined chunks. This spends
        most of its time storing CATG data that will probably not be used,
        but is important for saving SNP depths info.
        """
        # concatenate and store catgs
        start = time.time()
        printstr = ("indexing alleles    ", "s5")
        prog = AssemblyProgressBar({1:1}, start, printstr, self.data)
        prog.update()

        # concat catgs for each sample
        asyncs1 = {}
        for sample in self.samples:
            args = (self.data, sample, self.isref)
            asyncs1[sample.name] = self.lbview.apply(concat_catgs, *args)

        # select the next job
        if self.isref:
            concat_job = concat_reference_consens
        else:
            concat_job = concat_denovo_consens

        # collect all results for a sample and store stats 
        asyncs2 = {}
        for sample in self.samples:
            args = (self.data, sample)
            rasync = self.lbview.apply(concat_job, *args)
            asyncs2[sample.name] = rasync
            
        # track progress of stats storage
        alljobs = list(asyncs1.values()) + list(asyncs2.values())
        prog.jobs = dict(enumerate(alljobs))
        prog.update()
        prog.block()
        prog.check()



    def data_store(self, statsdicts):
        """
        Store assembly object stats
        """       
        # store sample stats
        for sample in self.samples:
            store_sample_stats(self.data, sample, statsdicts[sample.name])

        # build Assembly stats
        self.data.stats_dfs.s5 = self.data._build_stat("s5")

        # write stats file
        self.data.stats_files.s5 = os.path.join(
            self.data.dirs.consens, 
            's5_consens_stats.txt')

        with open(self.data.stats_files.s5, 'w') as out:
            self.data.stats_dfs.s5.to_string(
                buf=out,
                formatters={
                    'clusters_total': '{:.0f}'.format,
                    'filtered_by_depth': '{:.0f}'.format,
                    'filtered_by_maxH': '{:.0f}'.format,
                    'filtered_by_maxN': '{:.0f}'.format,
                    'reads_consens': '{:.0f}'.format,
                    'nsites': '{:.0f}'.format,
                    'nhetero': '{:.0f}'.format,
                    'heterozygosity': '{:.5f}'.format
                })




def make_chunks(data, sample, ncpus):
    """
    Split job into bits and pass to the client
    """
    # counter for split job submission
    num = 0

    # set optim size for chunks in N clusters. The first few chunks take longer
    # because they contain larger clusters, so we create 4X as many chunks as
    # processors so that they are split more evenly.
    optim = int(
        (sample.stats.clusters_total // ncpus) + \
        (sample.stats.clusters_total % ncpus))

    # open to clusters
    with gzip.open(sample.files.clusters, 'rb') as clusters:
        # create iterator to sample 2 lines at a time
        pairdealer = izip(*[iter(clusters)] * 2)

        # Use iterator to sample til end of cluster
        done = 0
        while not done:
            # grab optim clusters and write to file.
            done, chunk = clustdealer(pairdealer, optim)
            chunk = [i.decode() for i in chunk]

            # make file handle
            chunkhandle = os.path.join(
                data.tmpdir,
                "{}.chunk.{}.{}".format(sample.name, optim, num * optim))

            # write to file
            if chunk:
                with open(chunkhandle, 'wt') as outchunk:
                    outchunk.write("//\n//\n".join(chunk) + "//\n//\n")
                num += 1
    logger.debug("chunked {} to {} files".format(sample.files.clusters, num))



def process_chunks(data, sample, chunkfile, isref):
    """
    Remote callable function to use Processor class instance
    """
    proc = Processor(data, sample, chunkfile, isref)
    proc.run()
    return proc.counters, proc.filters     



class Processor:
    """
    The consensus calling process that calls alleles and genotypes and
    applies filters to consens seqs.
    """
    def __init__(self, data, sample, chunkfile, isref):

        # input params
        self.data = data
        self.sample = sample
        self.chunkfile = chunkfile
        self.isref = isref

        # set max limits and params from data
        self.nalleles = 1
        self.tmpnum = int(self.chunkfile.split(".")[-1])
        self.optim = int(self.chunkfile.split(".")[-2])
        self.est_err = self.data.stats.error_est.mean()
        self.est_het = self.data.stats.hetero_est.mean()
        self.maxlen = self.data.hackersonly.max_fragment_length
        self.maxhet = self.data.params.max_Hs_consens
<<<<<<< HEAD
        self.maxn = self.data.params.max_Ns_consens
        self.maxa = self.data.params.max_alleles_consens
=======
        self.maxalleles = self.data.params.max_alleles_consens

>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7
        # not enforced for ref
        self.maxn = self.data.params.max_Ns_consens
        if self.isref:
            self.maxn = int(1e6)

        # target attrs to be filled
        self.counters = {}
        self.filters = {}
        self.storeseqs = {}

        # tmp attrs to be filled
        self.names = None
        self.seqs = None
        self.ref_position = None
        self.consens = None
        self.hidx = None
        self.nheteros = None

        # local copies to use to fill the arrays
        self.catarr = np.zeros((self.optim, self.maxlen, 4), dtype=np.uint32)
        self.nallel = np.zeros((self.optim, ), dtype=np.uint8)
        self.refarr = np.zeros((self.optim, 3), dtype=np.int64)

        # store data for stats counters.
        self.counters = {
            "name": self.tmpnum,
            "heteros": 0,
            "nsites": 0,
            "nconsens": 0,
        }

        # store data for what got filtered
        self.filters = {
            "depth": 0,
            "maxh": 0,
            "maxn": 0,
<<<<<<< HEAD
            "maxa": 0,
=======
            "maxalleles": 0,
>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7
        }

        # store data for writing
        self.storeseq = {}

        # if reference-mapped then parse the fai to get index number of chroms
        if self.isref:
            fai = pd.read_csv(
                self.data.params.reference_sequence + ".fai",
                names=['scaffold', 'size', 'sumsize', 'a', 'b'],
                sep="\t",
                dtype=object)
            self.faidict = {j: i for i, j in enumerate(fai.scaffold)}
            self.revdict = {j: i for i, j in self.faidict.items()}


    def run(self):
        """
        The main function to run the processor
        """
        self.process_chunk()
        self.write_chunk()


    def process_chunk(self):
        """
        iterate over clusters processing each sequentially.
        """
        logger.debug("consens calling on {}".format(self.chunkfile))

        # load 2 lines at a time to get <name> <sequence>
        inclust = open(self.chunkfile, 'rb')
        pairdealer = izip(*[iter(inclust)] * 2)
        done = 0

        # stream through the clusters
        while not done:

            # repeat pairdealer until an entire cluster is pulled in.
            # returns 1 if it reaches the end of the file.
            done, chunk = clustdealer(pairdealer, 1)

            # if chunk was returned then process it.
            if not chunk:
                continue

            # fills .name and .seqs attributes
            # seqs is uint8 with missing as 78 or 45                
            self.parse_cluster(chunk)

            # denovo only: mask repeats (drops lowcov sites w/ dashes) and
            # converts remaining dashes to Ns
            if not self.isref:
                # return 1 if entire seq was not masked
                if self.mask_repeats():
                    continue

            # return 1 if cov > mindepth_mj (or >max), filter[depth] + 1
            if self.filter_mindepth():
                continue

            # fills .consens with genotype calls and trims .seqs
            # consens is uint8 with all missing as Ns (78)
            triallele = self.new_build_consens()

            # simple 3rd-allele filter, filter[maxallele] + 1
            # Returns 1 if allele not allowed. See also allele filter.
            if self.filter_triallele(triallele):
                continue

            # fills .hidx and .nheteros
            self.get_heteros()

            # return 1 if too many heterozygote sites
            if self.filter_maxhetero():
                continue

            # return 1 if too many N or too short
            # maxN set arbitrarily high for REF to allow contig-building.
            if self.filter_max_n_min_len():
                continue

<<<<<<< HEAD
                            # return 1 if not too many N or too short 
                            if self.filter_maxN_minLen():
                                
                                # filter for max haplotypes...
                                # self.get_alleles()
                                # ...
                                if self.filter_max_alleles():
                                    # store result
                                    self.store_data()
=======
            # Infer haplotypes from reads at variable geno sites
            if self.nheteros > 1:
                if self.filter_alleles():
                    continue

            # store result
            self.store_data()
>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7

        # cleanup close handle
        inclust.close()



    def parse_cluster(self, chunk):
        """
        Read in cluster chunk to get .names & .seqs and ref position
        """
        # get names and seqs
        piece = chunk[0].decode().strip().split("\n")
        self.names = piece[0::2]
        seqs = piece[1::2]
        reps = [int(n.split(";")[-2][5:]) for n in self.names]

        # fill array with copies of each seq
        sarr = np.zeros((sum(reps), len(seqs[0])), dtype=np.uint8)
        idx = 0
        for ireps, iseq in zip(reps, seqs):
            for _ in range(ireps):
                sarr[idx] = np.array(list(iseq)).astype(bytes).view(np.uint8)
                idx += 1

        # trim array to allow maximum length
        self.seqs = sarr[:, :self.maxlen]

        # ref positions (chromint, pos_start, pos_end)
        self.ref_position = (-1, 0, 0)
        if self.isref:
            # parse position from name string
            rname = self.names[0].rsplit(";", 2)[0]
            chrom, posish = rname.rsplit(":")
            pos0, pos1 = posish.split("-")
            # drop `tag=ACGTACGT` if 3rad and declone_PCR_duplicates is set
            pos1 = pos1.split(";")[0]
            # pull idx from .fai reference dict
            chromint = self.faidict[chrom] + 1
            self.ref_position = (int(chromint), int(pos0), int(pos1))



    def mask_repeats(self):
        """
        Removes mask columns with low depth repeats from denovo clusters.
        """
        # get column counts of dashes
        idepths = np.sum(self.seqs == 45, axis=0).astype(float)

        # get proportion of bases that are dashes at each site
        props = idepths / self.seqs.shape[0] 

        # if proportion of dashes sites is more than 0.8?
        keep = props < 0.9

        # report to logger
        # if np.any(props >= 0.9):
        #     logger.debug(
        #         "repeat masked {}:\n{}"
        #         .format(np.where(~keep)[0], self.seqs[:, ~keep])
        #     )

        # drop sites from seqs
        self.seqs = self.seqs[:, keep]
        self.seqs[self.seqs == 45] = 78

        # apply filter in case ALL sites were dropped
        if not self.seqs.size:
            return 1
        return 0



    def filter_mindepth(self):
        """
        return 1 if READ depth < minimum
        """
        sumreps = self.seqs.shape[0]
        bool1 = sumreps >= self.data.params.mindepth_majrule
        bool2 = sumreps <= self.data.params.maxdepth
        if bool1 & bool2:       
            return 0

        # return that this cluster was filtered out
        self.filters['depth'] += 1
        return 1



    def new_build_consens(self):
        """
        Replace function below using int arrays. 
        Returns 1 if evidence of a triallele.
        """
        # get unphased genotype call of this sequence
        consens, triallele = new_base_caller(
            self.seqs, 
            self.data.params.mindepth_majrule, 
            self.data.params.mindepth_statistical, 
            self.est_het, 
            self.est_err,
        )

        # TODO: for denovo we could consider trimming Ns from internal
        # split near pair separator, applied here.
        # ...

        # trim Ns (uncalled genos) from the left and right
        trim = np.where(consens != 78)[0]
        ltrim = trim.min()
        rtrim = trim.max() + 1
        self.consens = consens[ltrim: rtrim]
        self.seqs = self.seqs[:, ltrim: rtrim]

        # update position for trimming
        self.ref_position = (
            self.ref_position[0], 
            self.ref_position[1] + ltrim,
            self.ref_position[1] + ltrim + rtrim,
        )
        return triallele



    def filter_triallele(self, triallele):
        """
        A simple filter on 3rd alleles if only allowing max 2
        """
        if triallele:
            if self.data.params.max_alleles_consens < 3:
                self.filters['maxalleles'] += 1
                return 1
        return 0



    def get_heteros(self):
        """
        Record the indices of heterozygous sites which will be used
        to identify the number of alleles at this locus. Easier to 
        redo this now after the seqs have been trimmed.

        # ambiguous sites
               R   K   S   Y   W   M
        array([82, 75, 83, 89, 87, 77], dtype=uint8)
        """
        hsites = np.any([
            self.consens == 82,
            self.consens == 75,
            self.consens == 83,
            self.consens == 89,
            self.consens == 87,
            self.consens == 77,                                                            
        ], axis=0)

        self.hidx = np.where(hsites)[0]
        self.nheteros = self.hidx.size



    def filter_maxhetero(self):
        """
        Return 1 if it PASSED the filter, else 0
        """
        if self.nheteros > (self.consens.size * self.maxhet):
            self.filters['maxh'] += 1
            return 1
        return 0



    def filter_max_n_min_len(self):
        """
        Return 1 if it PASSED the filter, else 0. The minLen filter
        here is treated the same as maxN since Ns compose the adjacent
        edges of the locus that were uncalled.
        """
        # trimmed too short in size then count as maxN
        if self.consens.size < self.data.params.filter_min_trim_len:
            self.filters['maxn'] += 1
            return 1
        if (self.consens == 78).sum() > (self.consens.size * self.maxn):
            self.filters['maxn'] += 1
            return 1
        return 0


<<<<<<< HEAD
    #def get_alleles(self):
    def filter_max_alleles(self):
=======

    def filter_alleles(self):
>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7
        """
        Infer the number of alleles from haplotypes.
        
        Easy case:
        AAAAAAAAAATAAAAAAAAAACAAAAAAAA
        AAAAAAAAAACAAAAAAAAAATAAAAAAAA

        T,C
        C,T

        Challenging case:
        AAAAAAAAAATAAAAAAAA-----------
        AAAAAAAAAACAAAAAAAA-----------
        -------------AAAAAAAACAAAAAAAA
        -------------AAAAAAAATAAAAAAAA

        T -
        C -
        - C
        - T
        ---------
        """
        # array of hetero sites (denovo: this can include Ns)
        harray = self.seqs[:, self.hidx]

        # remove varsites containing N
        harray = harray[~np.any(harray == 78, axis=1)]

        # remove alleles containing a site not called in bi-allele genos
        calls0 = np.array([DCONS[i][0] for i in self.consens[self.hidx]])
        calls1 = np.array([DCONS[i][1] for i in self.consens[self.hidx]])
        bicalls = (harray == calls0) | (harray == calls1)
        harray = harray[np.all(bicalls, axis=1), :]

        # get counts of each allele (e.g., AT:2, CG:2)
        ccx = Counter([tuple(i) for i in harray])

        # if max-alleles or fewer, then great, move on.
        if len(ccx) <= self.maxalleles:
            self.nalleles = len(ccx)
            return 0

        # too many alleles, BUT, try dropping v. low cov alleles
        alleles = []
        for allele in ccx:
            if (ccx[allele] / self.seqs.shape[0]) >= 0.1:
                alleles.append(allele)  # ccx[allele]
        
        # then filter if all were not dropped as lowcov
        arr = np.array(alleles)
        if not arr.size:
            self.filters['maxalleles'] += 1            
            self.nalleles = len(ccx)
            return 1

        # or if all bi-alleles are not still present
        try:
            if np.any(np.all(arr == arr[0], axis=1)):
                # logger.debug(
                    # "dropped lowcov alleles but remaining do not cover "
                    # "genotype calls, so filter is applied."
                    # )
                self.filters['maxalleles'] += 1
                self.nalleles = len(ccx)
                return 1
        except Exception as inst:
            logger.warning(ccx)
            logger.warning(arr)
            logger.exception(inst)
            raise IPyradError from inst

        # if passed filter-filter, how many alleles now?
        if len(alleles) <= self.maxalleles:
            self.nalleles = len(alleles)
            return 0

        # too many alleles
        self.filters['maxalleles'] += 1        
        self.nalleles = len(ccx)
        return 1

        # # require that remaining alleles contain hetero calls.

        # # if max-alleles or fewer, then great, move on.
        # if len(ccx) <= self.maxalleles:
        #     self.nalleles = len(ccx)
        #     return 0

        # # if max-alleles or fewer, then great, move on.
        # if len(alleles) <= self.maxalleles:
        #     self.nalleles = len(alleles)
        #     return 0

        #if self.nalleles == 2:
        #    self.consens = encode_alleles(self.consens, self.hidx, alleles)

<<<<<<< HEAD
        if self.nalleles > self.maxa:
            self.filters['maxa'] += 1
            return 0
        return 1
            # Do we need this still? iao 6/2021
            #if self.nalleles == 2:
            #    try:
            #        self.consens = storealleles(self.consens, self.hidx, alleles)
            #    except (IndexError, KeyError):
            #        pass
=======
>>>>>>> ff8f2462f57837696fa3c37046cbc7368d88b0d7


    def store_data(self):
        """
        Store the site count data for writing to HDF5, and stats.
        """
        # current counter
        cidx = self.counters["nconsens"]
        self.nallel[cidx] = self.nalleles
        self.refarr[cidx] = self.ref_position

        # store a reduced array with only CATG
        catg = np.array(
            [np.sum(self.seqs == i, axis=0) for i in (67, 65, 84, 71)],
            dtype=np.uint16,
        ).T

        # do not allow ints larger than 65535 (uint16)
        self.catarr[cidx, :catg.shape[0], :] = catg

        # store the seqdata and advance counters
        self.storeseq[cidx] = bytes(self.consens)
        self.counters["name"] += 1
        self.counters["nconsens"] += 1
        self.counters["heteros"] += self.nheteros



    def write_chunk(self):
        """
        Writes chunk of consens reads to disk, stores depths, alleles, and 
        chroms, and stores stats. For denovo data it writes consens chunk
        as a fasta file. For reference data it writes as a SAM file that is
        compatible to be converted to BAM (very stringent about cigars).
        """
        # the reference arr is bigger than the actual seqarr, so find the
        # size of the actually filled data array (end)
        end = np.where(np.all(self.refarr == 0, axis=1))[0]
        if np.any(end):
            end = end.min()
        else:
            end = self.refarr.shape[0]

        # write final consens string chunk
        consenshandle = os.path.join(
            self.data.tmpdir,
            "{}_tmpcons.{}.{}".format(self.sample.name, end, self.tmpnum)
        )

        # write chunk. 
        if self.storeseq:
            outfile = open(consenshandle, 'wt')

            # denovo just write the consens simple
            if not self.isref:
                seqlist = []
                for key in self.storeseq:
                    seq = self.storeseq[key].decode()
                    cstring = ">{}_{}\n{}".format(self.sample.name, key, seq)
                    seqlist.append(cstring)
                outfile.write("\n".join(seqlist))
                # constring = "\n".join(
                    # [">" + self.sample.name + "_" + str(key) + \
                    # "\n" + self.storeseq[key].decode() 
                    # for key in self.storeseq])
                # outfile.write(constring)

            # reference needs to store if the read is revcomp to reference
            else:
                # seqlist = []
                # for key in self.storeseq:
                    # pass
                constring = "\n".join(
                    ["{}:{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                    .format(
                        self.sample.name,
                        self.refarr[i][0],
                        self.refarr[i][1],
                        self.refarr[i][1] + len(self.storeseq[i]),
                        #self.refarr[i][2],
                        0,
                        self.revdict[self.refarr[i][0] - 1],
                        self.refarr[i][1],
                        0,
                        make_cigar(
                            np.array(list(self.storeseq[i].decode()))
                            ),
                        "*",
                        0,
                        #self.refarr[i][2] - self.refarr[i][1],
                        len(self.storeseq[i]),
                        self.storeseq[i].decode(),
                        "*",
                    ) for i in self.storeseq]
                )
                outfile.write(constring)

            outfile.close()

        # store reduced size arrays with indexes matching to keep indexes
        tmp5 = consenshandle.replace("_tmpcons.", "_tmpcats.")
        with h5py.File(tmp5, 'w') as io5:
            # reduce catarr to uint16 size for easier storage
            self.catarr[self.catarr > 65535] = 65535
            self.catarr = self.catarr.astype(np.uint16)
            io5.create_dataset(name="cats", data=self.catarr[:end])
            io5.create_dataset(name="alls", data=self.nallel[:end])
            io5.create_dataset(name="chroms", data=self.refarr[:end])
        del self.catarr
        del self.nallel
        del self.refarr

        # return stats and skip sites that are Ns (78)
        #self.counters['nsites'] = sum([len(i) for i in self.storeseq.values()])
        self.counters['nsites'] = sum(
            sum(1 if j != 78 else 0 for j in i) 
            for i in self.storeseq.values()
        )
        del self.storeseq



def new_base_caller(sarr, mindepth_mj, mindepth_st, est_het, est_err):
    """
    Call site genotypes from site counts. Uses scipy binomial.
    """
    # start with an array of Ns
    cons = np.zeros(sarr.shape[1], dtype=np.uint8)
    cons.fill(78)

    # record if evidence of a tri-allele
    triallele = 0

    # iterate over columns
    for cidx in range(sarr.shape[1]):

        # get col, nmask, dashmask, and bothmask
        col = sarr[:, cidx]
        nmask = col == 78
        dmask = col == 45
        bmask = nmask | dmask

        # if non-masked site depth is below mindepth majrule fill with N (78)
        if np.sum(~bmask) < mindepth_mj:
            cons[cidx] = 78
            continue
                
        # if not variable (this will include the 'nnnn' pair separator (110))
        if np.all(col == col[0]):
            cons[cidx] = col[0]
            continue

        # make statistical base calls on allele frequencies
        counts = np.bincount(col, minlength=79)
        counts[78] = 0

        # get allele freqs (first-most, second, third = p, q, r)
        pbase = np.argmax(counts)
        nump = counts[pbase]
        counts[pbase] = 0

        qbase = np.argmax(counts)
        numq = counts[qbase]
        counts[qbase] = 0

        rbase = np.argmax(counts)
        numr = counts[rbase]
        counts[rbase] = 0

        # if third allele occurs >X% then fill N and mark as paralog
        if (numr / (nump + numq + numr)) > 0.2:
            triallele = 1

        # based on biallelic depth
        bidepth = nump + numq
        if bidepth < mindepth_mj:
            cons[cidx] = 78
            continue

        # if depth is too high, reduce to sampled int
        if bidepth > 500:
            nump = int(500 * (nump / float(bidepth)))
            numq = int(500 * (numq / float(bidepth)))

        # make majority-rule call
        if bidepth < mindepth_st:
            if nump == numq:
                cons[cidx] = TRANS[(pbase, qbase)]
            else:
                cons[cidx] = pbase        

        # make statistical base call
        ishet, prob = get_binom(nump, numq, est_err, est_het)
        if prob < 0.95:
            cons[cidx] = 78
        else:
            if ishet:
                cons[cidx] = TRANS[(pbase, qbase)]
            else:
                cons[cidx] = pbase

    return cons, triallele



def concat_catgs(data, sample, isref):
    """
    Concat catgs into a single sample catg and remove tmp files
    """
    # collect tmpcat files written by write_chunks()
    tmpcats = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcats.*".format(sample.name))
    )
    tmpcats.sort(key=lambda x: int(x.split(".")[-1]))

    # get full nrows of the new h5 from the tmpcat filenames
    nrows = sum([int(i.rsplit(".", 2)[-2]) for i in tmpcats])

    # get shape info from the first cat, (optim, maxlen, 4)
    optim = min(nrows, 5000)
    with h5py.File(tmpcats[0], 'r') as io5:
        _, maxlen, _ = io5['cats'].shape

    # Check values of nrows and maxlen are > 0
    # This literally shouldn't happen, but it does, or has at least twice.
    # Related to issue #369
    if not all([nrows, maxlen]):
        raise IPyradError(
            "Error in concat_catgs both nrows and maxlen must be positive."
            "\nsample: {}\tnrows: {}\tmaxlen: {}"
            .format(sample.name, nrows, maxlen))

    # fill in the chunk array
    with h5py.File(sample.files.database, 'w') as ioh5:
        dcat = ioh5.create_dataset(
            name="catg",
            shape=(nrows, maxlen, 4),
            dtype=np.uint32,
            chunks=(optim, maxlen, 4),
            compression="gzip")
        dall = ioh5.create_dataset(
            name="nalleles", 
            shape=(nrows, ),
            dtype=np.uint8,
            chunks=(optim, ),
            compression="gzip")
        
        # only create chrom for reference-aligned data
        if isref:
            dchrom = ioh5.create_dataset(
                name="chroms",
                shape=(nrows, 3),
                dtype=np.int64,
                chunks=(optim, 3),
                compression="gzip")

        # Combine all those tmp cats into the big cat
        start = 0
        for icat in tmpcats:
            addon = int(icat.rsplit(".", 2)[-2])
            end = start + addon
            io5 = h5py.File(icat, 'r')
            dcat[start:end] = io5['cats'][:]
            dall[start:end] = io5['alls'][:]
            if isref:
                dchrom[start:end] = io5['chroms'][:]
            start = end
            io5.close()
            os.remove(icat)


def concat_denovo_consens(data, sample):
    """
    Concatenate consens bits into fasta file for denovo assemblies
    """
    # collect consens chunk files
    combs1 = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcons.*".format(sample.name)))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    # stream through file renaming consens reads and write out
    idx = 0
    with gzip.open(sample.files.consens, 'wt') as out:
        for fname in combs1:
            cstack = []
            #starter = int(fname.rsplit(".")[-1])
            with open(fname) as infile:
                for line in infile:
                    if line.startswith(">"):
                        name, num = line.rsplit("_", 1)
                        line = name + "_{}".format(idx)
                        idx += 1
                        #int(num) + starter)
                    cstack.append(line.strip())
                out.write("\n".join(cstack) + "\n")
            os.remove(fname)


def concat_reference_consens(data, sample):
    """
    Concatenates consens bits into SAM for reference assemblies
    """
    # write sequences to SAM file
    combs1 = glob.glob(os.path.join(
        data.tmpdir,
        "{}_tmpcons.*".format(sample.name)))
    combs1.sort(key=lambda x: int(x.split(".")[-1]))

    # open sam handle for writing to bam
    samfile = sample.files.consens.replace(".bam", ".sam")    
    with open(samfile, 'w') as outf:

        # parse fai file for writing headers
        fai = "{}.fai".format(data.params.reference_sequence)
        fad = pd.read_csv(fai, sep="\t", names=["SN", "LN", "POS", "N1", "N2"])
        headers = ["@HD\tVN:1.0\tSO:coordinate"]
        headers += [
            "@SQ\tSN:{}\tLN:{}".format(i, j)
            for (i, j) in zip(fad["SN"], fad["LN"])
        ]
        outf.write("\n".join(headers) + "\n")

        # write to file with sample names imputed to line up with catg array
        counter = 0
        for fname in combs1:
            with open(fname) as infile:
                # impute catg ordered seqnames 
                data = infile.readlines()
                fdata = []
                for line in data:
                    name, chrom, rest = line.rsplit(":", 2)
                    fdata.append(
                        "{}_{}:{}:{}".format(name, counter, chrom, rest)
                        )
                    counter += 1
                outf.write("".join(fdata) + "\n")

    # convert to bam
    cmd = [ip.bins.samtools, "view", "-Sb", samfile]
    with open(sample.files.consens, 'w') as outbam:
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=outbam)
        comm = proc.communicate()[1]
    if proc.returncode:
        raise IPyradError("error in samtools: {}".format(comm))
    else:
        pass
        #os.remove(samfile)
        #for fname in combs1:
        #    os.remove(fname)


def store_sample_stats(data, sample, statsdicts):
    """
    Not parallel, store the sample objects stats
    """
    # record results
    xcounters = {
        "nconsens": 0,
        "heteros": 0,
        "nsites": 0,
    }
    xfilters = {
        "depth": 0,
        "maxh": 0,
        "maxn": 0,
        "maxa": 0,
    }

    # merge finished consens stats
    for counters, filters in statsdicts:
        # sum individual counters
        for key in xcounters:
            xcounters[key] += counters[key]
        for key in xfilters:
            xfilters[key] += filters[key]

    # set Sample stats_dfs values
    if int(xcounters['nsites']):
        prop = int(xcounters["heteros"]) / float(xcounters['nsites'])
    else:
        prop = 0

    # store stats attributes to the sample
    sample.stats_dfs.s5.nsites = int(xcounters["nsites"])
    sample.stats_dfs.s5.nhetero = int(xcounters["heteros"])
    sample.stats_dfs.s5.filtered_by_depth = xfilters['depth']
    sample.stats_dfs.s5.filtered_by_maxH = xfilters['maxh']
    sample.stats_dfs.s5.filtered_by_maxN = xfilters['maxn']
    sample.stats_dfs.s5.filtered_by_maxAlleles = xfilters['maxa']
    sample.stats_dfs.s5.reads_consens = int(xcounters["nconsens"])
    sample.stats_dfs.s5.clusters_total = sample.stats_dfs.s3.clusters_total
    sample.stats_dfs.s5.heterozygosity = float(prop)

    # set the Sample stats summary value
    sample.stats.reads_consens = int(xcounters["nconsens"])

    # save state to Sample if successful
    if sample.stats.reads_consens:
        sample.stats.state = 5
    else:
        print("No clusters passed filtering in Sample: {}".format(sample.name))


def make_cigar(arr):
    """
    Writes a cigar string with locations of indels and lower case ambigs
    """
    # simplify data
    arr[np.char.islower(arr)] = '.'
    indel = np.bool_(arr == "-")
    ambig = np.bool_(arr == ".")
    arr[~(indel + ambig)] = "A"

    # counters
    cigar = ""
    mcount = 0
    tcount = 0
    lastbit = arr[0]
    for i, j in enumerate(arr):

        # write to cigarstring when state change
        if j != lastbit:
            if mcount:
                cigar += "{}{}".format(mcount, "M")
                mcount = 0
            else:
                cigar += "{}{}".format(tcount, CIGARDICT.get(lastbit))
                tcount = 0
            mcount = 0
            
        # increase counters
        if (j == '.' or j == '-'):
            tcount += 1
        else:
            mcount += 1
        lastbit = j
        
    # write final block
    if mcount:
        cigar += "{}{}".format(mcount, "M")
    if tcount:
        cigar += '{}{}'.format(tcount, CIGARDICT.get(lastbit))
    return cigar


CIGARDICT = {
    '-': "I",
    '.': "S",
}


# this is used in write_chunk for reference mapped data.
def make_allele_cigar(seq, on='.', letter='S'):
    iii = seq.split(on)
    cig = ""
    lastletter = ""
    isi = 0
    for i, j in enumerate(iii):    
        if j == "":
            isi += 1
               
        # store matches
        else:
            if isi:
                # skip if first indel base
                if cig:
                    isi += 1
                cig += "{}{}".format(isi, letter)
                isi = 0
                lastletter = letter
            
            if lastletter == "M":
                cig += "{}{}".format(1, letter)
            cig += "{}M".format(len(j))
            lastletter = "M"
    if isi:
        cig += "{}{}".format(isi, letter)
    return cig



# not currently used
def make_indel_cigar(seq, on='-', letter='I'):
    iii = seq.split(on)
    cig = ""
    lastletter = ""
    isi = 0
    for i, j in enumerate(iii):    
        if j == "":
            isi += 1
               
        # store matches
        else:
            if isi:
                # skip if first indel base
                if cig:
                    isi += 1
                cig += "{}{}".format(isi, letter)
                isi = 0
                lastletter = letter
            
            if lastletter == "M":
                cig += "{}{}".format(1, letter)
                
            acig = ''.join('.' if i.islower() else i for i in j)
            mcig = make_allele_cigar(acig)
            cig += mcig
            lastletter = "M"
    if isi:
        cig += "{}{}".format(isi, letter)
    return cig







def get_binom(base1, base2, est_err, est_het):
    """
    return probability of base call
    """
    prior_homo = (1. - est_het) / 2.
    prior_hete = est_het

    ## calculate probs
    bsum = base1 + base2
    hetprob = scipy.special.comb(bsum, base1) / (2. ** (bsum))
    homoa = scipy.stats.binom.pmf(base2, bsum, est_err)
    homob = scipy.stats.binom.pmf(base1, bsum, est_err)

    ## calculate probs
    hetprob *= prior_hete
    homoa *= prior_homo
    homob *= prior_homo

    ## final
    probabilities = [homoa, homob, hetprob]
    bestprob = max(probabilities) / float(sum(probabilities))

    ## return
    if hetprob > homoa:
        return True, bestprob
    return False, bestprob



# not currently used in reference assemblies
def mask_repeats(consens, arrayed):
    """
    Checks for interior Ns in consensus seqs and removes those that are at
    low depth, here defined as less than 1/3 of the average depth. The prop 1/3
    is chosen so that mindepth=6 requires 2 base calls that are not in [N,-].

    Python3 notes:
    consens and arrayed are both in bytes in entry. Consens is converted to
    unicode for operations, and both are returned as bytes.
    """
    ## default trim no edges
    consens[consens == b"-"] = b"N"
    consens = b"".join(consens)

    ## split for pairs
    try:
        cons1, cons2 = consens.split(b"nnnn")
        split = consens.index(b"nnnn")
        arr1 = arrayed[:, :split]
        arr2 = arrayed[:, split + 4:]
    except ValueError:
        cons1 = consens
        cons2 = ""
        arr1 = arrayed

    ## trim from left and right of cons1
    edges = [None, None]
    lcons = len(cons1)
    cons1 = cons1.lstrip(b"N")
    edges[0] = lcons - len(cons1)

    ## trim from right if nonzero
    lcons = len(cons1)
    cons1 = cons1.rstrip(b"N")
    if lcons - len(cons1):
        edges[1] = -1 * (lcons - len(cons1))

    ## trim same from arrayed
    arr1 = arr1[:, edges[0]:edges[1]]

    ## trim from left and right of cons2 if present
    if cons2:
        ## trim from left and right of cons1
        edges = [None, None]
        lcons = len(cons2)
        cons2 = cons2.lstrip(b"N")
        edges[0] = lcons - len(cons2)

        ## trim from right if nonzero
        lcons = len(cons2)
        cons2 = cons2.rstrip(b"N")
        if lcons - len(cons2):
            edges[1] = -1 * (lcons - len(cons2))

        ## trim same from arrayed
        arr2 = arr2[:, edges[0]:edges[1]]

        ## reconstitute pairs
        consens = cons1 + b"nnnn" + cons2
        consens = np.array(list(consens), dtype=bytes)
        sep = np.array(arr1.shape[0] * [list(b"nnnn")])
        arrayed = np.hstack([arr1, sep, arr2])

    ## if single-end...
    else:
        consens = np.array(list(cons1), dtype=bytes)
        arrayed = arr1

    ## get column counts of Ns and -s
    ndepths = np.sum(arrayed == b'N', axis=0)
    idepths = np.sum(arrayed == b'-', axis=0)

    ## get proportion of bases that are N- at each site
    nons = ((ndepths + idepths) / float(arrayed.shape[0])) >= 0.75
    ## boolean of whether base was called N
    isn = consens == b"N"
    ## make ridx
    ridx = nons * isn

    ## apply filter
    consens = consens[~ridx]
    arrayed = arrayed[:, ~ridx]

    return consens, arrayed



def nfilter4(consens, hidx, arrayed):
    """
    Applies max haplotypes filter returns pass and consens
    """
    # if less than two Hs then there is only one allele
    if len(hidx) < 2:
        return consens, 1

    # store base calls for hetero sites
    harray = arrayed[:, hidx]

    # remove any reads that have N or - base calls at hetero sites
    # these cannot be used when calling alleles currently.
    harray = harray[~np.any(harray == b"-", axis=1)]
    harray = harray[~np.any(harray == b"N", axis=1)]

    # get counts of each allele (e.g., AT:2, CG:2)
    ccx = Counter([tuple(i) for i in harray])

    ## Two possibilities we would like to distinguish, but we can't. Therefore,
    ## we just throw away low depth third alleles that are within seq. error.
    ## 1) a third base came up as a sequencing error but is not a unique allele
    ## 2) a third or more unique allele is there but at low frequency

    ## remove low freq alleles if more than 2, since they may reflect
    ## sequencing errors at hetero sites, making a third allele, or a new
    ## allelic combination that is not real.
    if len(ccx) > 2:
        totdepth = harray.shape[0]
        cutoff = max(1, totdepth // 10)
        alleles = [i for i in ccx if ccx[i] > cutoff]
    else:
        alleles = ccx.keys()

    ## how many high depth alleles?
    nalleles = len(alleles)

    ## if 2 alleles then save the phase using lowercase coding. 
    # todo: store for each allele whether it is phased or not.
    if nalleles == 2:
        try:
            consens = storealleles(consens, hidx, alleles)
        except (IndexError, KeyError):
            pass

        return consens, nalleles

    ## just return the info for later filtering
    else:
        return consens, nalleles



def encode_alleles(consens, hidx, alleles):
    """ 
    Store phased allele data for diploids 
    """
    ## find the first hetero site and choose the priority base
    ## example, if W: then priority base in A and not T. PRIORITY=(order: CATG)
    bigbase = PRIORITY[consens[hidx[0]]]

    ## find which allele has priority based on bigbase
    bigallele = [i for i in alleles if i[0] == bigbase][0]

    ## uplow other bases relative to this one and the priority list
    ## e.g., if there are two hetero sites (WY) and the two alleles are
    ## AT and TC, then since bigbase of (W) is A second hetero site should
    ## be stored as y, since the ordering is swapped in this case; the priority
    ## base (C versus T) is C, but C goes with the minor base at h site 1.
    #consens = list(consens)
    for hsite, pbase in zip(hidx[1:], bigallele[1:]):
        if PRIORITY[consens[hsite]] != pbase:
            consens[hsite] = consens[hsite].lower()

    ## return consens
    return consens






# DEPRECATED

# def build_consens_and_array(self):
#     """
#     Makes base calls and converts - to N and trims terminal Ns. Setting 
#     internal - to Ns makes handling the insert of paired reference mapped
#     data much better... but puts N's into denovo data where we might other
#     wise choose to drop those columns... Think more about this...
#     """
#     # get stacks of base counts
#     sseqs = [list(seq) for seq in self.seqs]
#     arrayed = np.concatenate(
#         [[seq] * rep for (seq, rep) in zip(sseqs, self.reps)]
#     ).astype(bytes)

#     # ! enforce maxlen limit !
#     self.arrayed = arrayed[:, :self.maxlen]
                
#     # get unphased consens sequence from arrayed
#     self.consens = base_caller(
#         self.arrayed, 
#         self.data.params.mindepth_majrule, 
#         self.data.params.mindepth_statistical,
#         self.esth, 
#         self.este,
#     )

#     # trim Ns from the left and right ends
#     mask = self.consens.copy()
#     mask[mask == b"-"] = b"N"
#     trim = np.where(mask != b"N")[0]

#     # bail out b/c no bases were called 
#     if not trim.size:
#         self.filters['depth'] += 1
#         return 0

#     ltrim, rtrim = trim.min(), trim.max()
#     self.consens = self.consens[ltrim:rtrim + 1]
#     self.arrayed = self.arrayed[:, ltrim:rtrim + 1]

#     # update position for trimming
#     self.ref_position = (
#         self.ref_position[0], 
#         self.ref_position[1] + ltrim,
#         self.ref_position[1] + ltrim + rtrim + 1,
#     )
#     return 1



# def base_caller(arrayed, mindepth_majrule, mindepth_statistical, est_het, est_err):
#     """
#     Call all sites in a locus array. Can't be jit'd yet b/c scipy
#     """
#     # an array to fill with consensus site calls
#     cons = np.zeros(arrayed.shape[1], dtype=np.uint8)
#     cons.fill(78)
#     arr = arrayed.view(np.uint8)

#     # iterate over columns
#     for col in range(arr.shape[1]):
#         # the site of focus
#         carr = arr[:, col]

#         # if site is all dash then fill it dash (45)
#         if np.all(carr == 45):
#             cons[col] = 45
            
#         # else mask all N and - sites for base call
#         else:
#             mask = carr == 45
#             mask += carr == 78
#             marr = carr[~mask]
            
#             # call N if no real bases, or below majrule.
#             if marr.shape[0] < mindepth_majrule:
#                 cons[col] = 78
                
#             # if not variable
#             elif np.all(marr == marr[0]):
#                 cons[col] = marr[0]

#             # estimate variable site call
#             else:
#                 # get allele freqs (first-most, second, third = p, q, r)
#                 counts = np.bincount(marr)

#                 pbase = np.argmax(counts)
#                 nump = counts[pbase]
#                 counts[pbase] = 0

#                 qbase = np.argmax(counts)
#                 numq = counts[qbase]
#                 counts[qbase] = 0

#                 ## based on biallelic depth
#                 bidepth = nump + numq
#                 if bidepth < mindepth_majrule:
#                     cons[col] = 78

#                 else:
#                     # if depth is too high, reduce to sampled int
#                     if bidepth > 500:
#                         base1 = int(500 * (nump / float(bidepth)))
#                         base2 = int(500 * (numq / float(bidepth)))
#                     else:
#                         base1 = nump
#                         base2 = numq

#                     # make statistical base call
#                     if bidepth >= mindepth_statistical:
#                         ishet, prob = get_binom(base1, base2, est_err, est_het)
#                         if prob < 0.95:
#                             cons[col] = 78
#                         else:
#                             if ishet:
#                                 cons[col] = TRANS[(pbase, qbase)]
#                             else:
#                                 cons[col] = pbase

#                     # make majrule base call
#                     else:
#                         if nump == numq:
#                             cons[col] = TRANS[(pbase, qbase)]
#                         else:
#                             cons[col] = pbase
#     return cons.view("S1")




if __name__ == "__main__":


    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    # self.data.hackersonly.declone_PCR_duplicates:
    tdata = ip.load_json("/tmp/test-amaranth-denovo.json")
    # tdata.ipcluster['cores'] = 4
    tdata.run("5", auto=True, force=True)
    logger.info(tdata.stats.T)
