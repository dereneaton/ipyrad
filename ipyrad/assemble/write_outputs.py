#!/usr/bin/env python

from __future__ import print_function
from builtins import range

import os
import shutil
import numpy as np
from .util import IPyradError  # , splitalleles, BTS, AMBIGS, TRANSFULL

# suppress this terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class Step7:
    def __init__(self, data, ipyclient):
        self.data = data
        self.samples = self.get_subsamples()
        self.setup_dirs()

    def run(self):
        self.init_database()
        self.remote_fill_filters_and_edges()
        self.remote_build_loci_and_stats()
        self.remote_fill_depths()
        self.remote_fill_arrs()
        self.remote_build_vcf()
        self.remote_build_conversions()


    def get_subsamples(self):
        "get subsamples for this assembly. All must have been in step6"
        # filter samples by state
        state5 = self.data.stats.index[self.data.stats.state < 6]
        state6 = self.data.stats.index[self.data.stats.state == 6]
        state7 = self.data.stats.index[self.data.stats.state > 6]

        # tell user which samples are not ready for step5
        if state5.any():
            print("skipping samples not in state==6:\n{}"
                  .format(state5.tolist()))

        if self.force:
            # run all samples above state 5
            subs = self.data.stats.index[self.data.stats.state > 5]
            subsamples = [self.data.samples[i] for i in subs]

        else:
            # tell user which samples have already completed step 6
            if state7.any():
                raise IPyradError(
                    "some samples are already in state==7. If you wish to \n" \
                  + "create new outfiles for this assembly use the force \n" \
                  + "argument.")
            # run all samples in state 6
            subsamples = [self.data.samples[i] for i in state6]

        # check that kept samples were in the step6 database
        checked_samples = []
        for sample in subsamples:
            if sample.stats.reads_consens:
                checked_samples.append(sample)
            else:
                print("skipping {}; no consensus reads found.")
        if not any(checked_samples):
            raise IPyradError("no samples ready for step 6")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_consens,
            reverse=True,
        )
        return checked_samples


    def setup_dirs(self):
        "Create temp h5 db for storing filters and depth variants"
        # get stats from step6 h5 and create new h5
        outdir = os.path.join(
            self.data.paramsdict["project_dir"],
            "{}_outfiles".format(self.data.name))
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
        if not os.path.exists(outdir):
            os.mkdirs(outdir)
        #self.data.database = h5py.File(self.data.database, 'w')


    def remote_fill_filters_and_edges(self):
        # slice out loci 1000 at a time
        for hslice in range(0, 5000, 1000):
            fill_filters_and_edges(self.data, hslice)

    def remote_build_loci_and_stats(self):
        build_loci_and_stats(self.data)

    def build_vcf(self):
        pass

    def build_conversions(self):
        pass


def fill_filters_and_edges(data, hslice):
    """
    filter loci in this chunk and trim edges.
    """
    # open database for writing
    io5 = h5py.File(data.database, 'w')
    filters = io5["filters"]
    edges = io5["edges"]

    # load loci into an array (1000, nsamples, maxlen)
    arr = np.array()

    # 



    # close database file
    io5.close()


def apply_filters(data, hslice):
    """
    Grab a chunk of loci from the HDF5 database. Apply filters and fill the
    the filters boolean array.

    The design of the filtering steps intentionally sacrifices some performance
    for an increase in readability, and extensibility. Calling multiple filter
    functions ends up running through the sequences per stack several times,
    but I felt this design made more sense, and also will easily allow us to
    add more filters in the future.
    """
    LOGGER.info("Entering filter_stacks")

    ## open h5 handles
    io5 = h5py.File(data.clust_database, 'r')
    co5 = h5py.File(data.database, 'r')
    ## get a chunk (hslice) of loci for the selected samples (sidx)
    #superseqs = io5["seqs"][hslice[0]:hslice[1], sidx,]
    ## get an int view of the seq array
    #superints = io5["seqs"][hslice[0]:hslice[1], sidx, :].view(np.int8)

    ## we need to use upper to skip lowercase allele storage
    ## this slows down the rate of loading in data by a ton.
    superints = np.char.upper(
        io5["seqs"][hslice[0]:hslice[1], sidx, ]
        ).view(np.int8)
    LOGGER.info("superints shape {}".format(superints.shape))

    ## fill edge filter
    ## get edges of superseqs and supercats, since edges need to be trimmed
    ## before counting hets, snps, inds. Technically, this could edge trim
    ## clusters to the point that they are below the minlen, and so this
    ## also constitutes a filter, though one that is uncommon. For this
    ## reason we have another filter called edgfilter.
    splits = co5["edges"][hslice[0]:hslice[1], 4]
    edgfilter, edgearr = get_edges(data, superints, splits)
    del splits
    LOGGER.info('passed edges %s', hslice[0])

    ## minsamp coverages filtered from superseqs
    minfilter = filter_minsamp(data, superints)
    LOGGER.info('passed minfilt %s', hslice[0])

    ## maxhets per site column from superseqs after trimming edges
    hetfilter = filter_maxhet(data, superints, edgearr)
    LOGGER.info('passed minhet %s', hslice[0])

    ## ploidy filter
    pldfilter = io5["nalleles"][hslice[0]:hslice[1]].max(axis=1) > (
        data.paramsdict["max_alleles_consens"])

    ## indel filter, needs a fresh superints b/c get_edges does (-)->(N)
    indfilter = filter_indels(data, superints, edgearr)
    LOGGER.info('passed minind %s', hslice[0])

    ## Build the .loci snpstring as an array (snps)
    ## shape = (chunk, 1) dtype=S1, or should it be (chunk, 2) for [-,*] ?
    snpfilter, snpsarr = filter_maxsnp(data, superints, edgearr)

    LOGGER.info("edg %s", edgfilter.sum())
    LOGGER.info("min %s", minfilter.sum())
    LOGGER.info("het %s", hetfilter.sum())
    LOGGER.info("pld %s", pldfilter.sum())
    LOGGER.info("snp %s", snpfilter.sum())
    LOGGER.info("ind %s", indfilter.sum())

    ## SAVE FILTERS AND INFO TO DISK BY SLICE NUMBER (.0.tmp.h5)
    chunkdir = os.path.join(data.dirs.outfiles, data.name + "_tmpchunks")

    handle = os.path.join(chunkdir, "edgf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, edgfilter)

    handle = os.path.join(chunkdir, "minf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, minfilter)

    handle = os.path.join(chunkdir, "hetf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, hetfilter)

    handle = os.path.join(chunkdir, "snpf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, snpfilter)

    handle = os.path.join(chunkdir, "pldf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, pldfilter)

    handle = os.path.join(chunkdir, "indf.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, indfilter)

    handle = os.path.join(chunkdir, "snpsarr.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, snpsarr)

    handle = os.path.join(chunkdir, "edgearr.{}.npy".format(hslice[0]))
    with open(handle, 'wb') as out:
        np.save(out, edgearr)

    io5.close()
    co5.close()



def get_ref_variant_depths(data, samples, locids, site):
    "get variant site depths from individual catg files for each sample"

    # get catg arrays
    catgs = [
        os.path.join(
            data.dirs.consens,
            "{}.catg".format(sample.name)
        ) for sample in samples]

    # open all arrays for getting... would this work in there are say 5K samps?

