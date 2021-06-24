#!/usr/bin/env python

"""
Call consensus alleles from clustered stacks within samples.
"""

from typing import Dict
import os
import gzip
import glob
import itertools
from collections import Counter
import numpy as np
import pandas as pd
from loguru import logger

from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.schema import Stats5
from ipyrad.assemble.base_step import BaseStep
from ipyrad.assemble.utils import IPyradError, clustdealer
from ipyrad.assemble.consens_utils import (
    Processor, 
    concat_catgs,
    concat_denovo_consens,
    concat_reference_consens,
)


class Step5(BaseStep):
    """
    Run Step3 clustering/mapping using vsearch or bwa
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data, 5, quiet, force)
        self.ipyclient = ipyclient
        self.lbview = self.ipyclient.load_balanced_view()
        self.max_frag = self.data.hackers.max_fragment_length
        self.nclusters_dict: Dict[str,int] = {}
        self.data.tmpdir = self.tmpdir
        self.data.stepdir = self.stepdir

        # sample paths to be used
        for sname in self.samples:
            self.samples[sname].files.depths = (
                os.path.join(self.stepdir, f"{sname}.catg.hdf5"))
            if self.data.is_ref:
                self.samples[sname].files.consens = (
                    os.path.join(self.stepdir, f"{sname}.bam"))
            else:
                self.samples[sname].files.consens = (
                    os.path.join(self.stepdir, f"{sname}.consens.gz"))

    def run(self):
        """
        Submit jobs to run either denovo, reference, or complex.
        """
        self.calculate_depths_and_max_frag()
        self.make_chunks()
        self.process_chunks()
        self.concatenate_chunks()
        self.data.save_json()


    def calculate_depths_and_max_frag(self):
        """
        Checks whether mindepth has changed and calc nclusters and maxlen
        """
        # send jobs to be processed on engines
        jobs = {}
        for sname in self.samples:
            args = (self.data, self.samples[sname])
            rasync = self.lbview.apply(recal_hidepth_majrule, *args)
            jobs[sname] = rasync
        msg = "calculating depths"
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()

        # check for failures and collect results
        for sname in prog.results:
            n_hidepth_clusts, max_frag = prog.results[sname]
            self.data.max_frag = int(max(self.max_frag, max_frag))
            self.nclusters_dict[sname] = n_hidepth_clusts
        logger.debug(
            f"max_fragment_length will be constrained to {self.data.max_frag}")
        logger.debug(
            f"high depth clusters for consensus calling: {self.nclusters_dict}")


    def make_chunks(self):
        """
        split clusters into chunks for parallel processing
        """
        msg = "chunking clusters"
        jobs = {}
        for sname in self.samples:
            args = (
                self.data, 
                self.samples[sname], 
                self.nclusters_dict[sname],
                len(self.ipyclient),
            )
            rasync = self.lbview.apply(make_chunks, *args)
            jobs[sname] = rasync
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()


    def process_chunks(self):
        """
        Process the cluster chunks into arrays and consens or bam files
        """
        def processor_wrap(data, sample, chunkfile):
            "wrapper for class function"
            proc = Processor(data, sample, chunkfile)
            proc.run()
            return proc.counters, proc.filters

        # get mean values across all samples
        self.data.error_est = np.mean([
            self.data.samples[i].stats_s4.error_est 
            for i in self.data.samples
        ])
        self.data.hetero_est = np.mean([
            self.data.samples[i].stats_s4.error_est 
            for i in self.data.samples
        ])

        # submit jobs (10 per sample === can be hundreds of jobs...)
        jobs = {sname: [] for sname in self.samples}
        for sname in self.samples:
            chunks = glob.glob(os.path.join(self.tmpdir, f"{sname}.chunk.*"))
            chunks.sort(key=lambda x: int(x.split('.')[-1]))
            for chunk in chunks:
                args = (self.data, self.samples[sname], chunk)
                jobs[sname].append(self.lbview.apply(processor_wrap, *args))

        # send chunks to be processed
        allsyncs = itertools.chain(*jobs.values())
        allsyncs = dict(enumerate(allsyncs))
        msg = "consens calling"
        prog = AssemblyProgressBar(allsyncs, msg, 5, self.quiet)
        prog.block()
        prog.check()

        # collect results from jobs since where they're grouped by sname
        for sname in jobs:
            rdicts = [i.get() for i in jobs[sname]]
            stat_dict = Counter()
            filter_dict = Counter()
            for tup in rdicts:
                stat_dict.update(tup[0])
                filter_dict.update(tup[1])

            # ONLY ADVANCE STATE IF CONSENS READS EXIST
            if stat_dict['nconsens']:
                self.samples[sname].state = 5

            # STORE STATS 
            self.samples[sname].stats_s5 = Stats5(
                cluster_total=self.nclusters_dict[sname],
                consensus_total=stat_dict['nconsens'],
                filtered_by_depth=filter_dict['depth'],
                filtered_by_max_h=filter_dict['maxh'],
                filtered_by_max_alleles=filter_dict['maxalleles'],
                filtered_by_max_n=filter_dict['maxn'],
                nsites=stat_dict['nsites'],
                nhetero=stat_dict['heteros'],
                heterozygosity=stat_dict['heteros'] / stat_dict['nsites'],
                min_depth_maj_during_step5=self.data.params.min_depth_majrule,
                min_depth_stat_during_step5=self.data.params.min_depth_statistical,
            )

        # WRITE TO STATS FILE and LOGGER
        stats_file = os.path.join(self.stepdir, "s5_stats_consensus.txt")
        statsdf = pd.DataFrame(
            index=sorted(self.samples),
            columns=list(self.samples[sname].stats_s5.dict()),
        )
        for sname in self.samples:
            istats = self.samples[sname].stats_s5.dict()
            for column in statsdf.columns:
                statsdf.loc[sname, column] = istats[column]
        with open(stats_file, 'w') as out:
            out.write(statsdf.to_string())
            logger.info("\n" + statsdf.iloc[:, :7].to_string())


    def concatenate_chunks(self):
        """
        Concatenate chunks and relabel for joined chunks. This spends
        most of its time storing CATG data that will probably not be used,
        but is important for saving SNP depths info.
        """
        # concat catgs for each sample
        jobs1 = {}
        jobs2 = {}
        for sname in self.samples:

            # submit h5 counts concat
            jobs1[sname] = self.lbview.apply(
                concat_catgs, *(self.data, self.samples[sname]))

            # submit seq file concat
            if not self.data.is_ref:
                jobs2[sname] = self.lbview.apply(
                    concat_denovo_consens, *(self.data, self.samples[sname]))
            else:
                jobs2[sname] = self.lbview.apply(
                    concat_reference_consens, *(self.data, self.samples[sname]))               

        jobs = list(jobs1.values()) + list(jobs2.values())
        jobs = dict(enumerate(jobs))
        msg = "indexing alleles"
        prog = AssemblyProgressBar(jobs, msg, 5, self.quiet)
        prog.block()
        prog.check()


def recal_hidepth_majrule(data, sample):
    """
    Returns the number of loci above the majrule threshold and the 
    max fragment length that will be allowed.
    """
    # if nothing changed then return max fragment length
    check1 = data.params.min_depth_majrule == sample.stats_s3.min_depth_maj_during_step3
    check2 = data.params.min_depth_statistical == sample.stats_s3.min_depth_stat_during_step3
    if check1 and check2:
        maxfrag = (
            4 + sample.stats_s3.mean_hidepth_cluster_length + 
            (2 * sample.stats_s3.std_hidepth_cluster_length)
        )
        return sample.stats_s3.clusters_hidepth, maxfrag

    # otherwise calculate depth again given the new mindepths settings.
    with gzip.open(sample.files.clusters, 'rt') as infile:
        pairdealer = zip(*[infile] * 2)

        # storage
        counts = []
        depths = []
        maxlen = []

        # start with cluster 0
        tdepth = 0
        tlen = 0
        tcount = 0

        # iterate until empty
        while 1:
            try:
                name, seq = next(pairdealer)
            except StopIteration:
                break

            # if at the end of a cluster
            if name == "//\n":
                depths.append(tdepth)
                maxlen.append(tlen)
                counts.append(tcount)
                tlen = 0
                tdepth = 0
                tcount = 0
            else:
                tdepth += int(name.strip().split("=")[-1][:-2])
                tlen = len(seq)
                tcount += 1

    # convert to arrays
    maxlens, depths = np.array(maxlen), np.array(depths)

    # get mask of clusters that are hidepth
    stat_mask = depths >= data.params.min_depth_statistical

    # get frag lenths of clusters that are hidepth
    lens_above_st = maxlens[stat_mask]

    # calculate frag length from hidepth lens
    try:       
        maxfrag = 4 + int(lens_above_st.mean() + (2. * lens_above_st.std()))
    except Exception as inst:
        raise IPyradError(
            "No clusts with depth sufficient for statistical basecalling. "
            f"I recommend you branch to drop this sample: {sample.name}"
            ) from inst
    return stat_mask.sum(), maxfrag


def make_chunks(data, sample, nclusters, ncpus):
    """
    Split job into bits and pass to the client
    """
    # counter for split job submission
    num = 0

    # set optim size for chunks in N clusters. The first few chunks take longer
    # because they contain larger clusters, so we create 4X as many chunks as
    # processors so that they are split more evenly.
    optim = int((nclusters // ncpus) + (nclusters % ncpus))

    # open to clusters
    with gzip.open(sample.files.clusters, 'rt') as clusters:
        pairdealer = zip(*[clusters] * 2)

        # Use iterator to sample til end of cluster
        done = 0
        while not done:
            # grab optim clusters and write to file.
            done, chunk = clustdealer(pairdealer, optim)
            chunkhandle = os.path.join(
                data.tmpdir, f"{sample.name}.chunk.{optim}.{num * optim}")

            # write to file
            if chunk:
                with open(chunkhandle, 'wt') as outchunk:
                    outchunk.write("//\n//\n".join(chunk) + "//\n//\n")
                num += 1


def store_sample_stats(sample, statsdicts):
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
        "maxalleles": 0,
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
    sample.stats_dfs.s5.filtered_by_max_alleles = int(xfilters['maxalleles'])
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




if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG", stderr=False, logfile="/tmp/test.log")
   
    TEST = ip.load_json("/tmp/TEST1.json")
    TEST.run("5", force=True, quiet=False)
