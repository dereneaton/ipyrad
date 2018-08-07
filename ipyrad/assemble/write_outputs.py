#!/usr/bin/env python

from __future__ import print_function
try:
    from builtins import range
    from itertools import izip, chain
except ImportError:
    from itertools import chain
    izip = zip


import os
import shutil
import numpy as np
from .utils import IPyradError, clustdealer  # , splitalleles, BTS, AMBIGS, TRANSFULL

# suppress this terrible h5 warning
import warnings
with warnings.catch_warnings(): 
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class Step7:
    def __init__(self, data, force, ipyclient):
        self.data = data
        self.force = force
        self.ipyclient = ipyclient
        self.samples = self.get_subsamples()
        self.setup_dirs()
        self.init_vars()
        self.init_database()

    def run(self):
        self.remote_fill_vars_filters_and_edges()
        self.remote_build_loci_and_stats()
        self.remote_fill_depths()
        self.remote_fill_arrs()
        self.remote_build_vcf()
        self.remote_build_conversions()

    ## init functions -----------------------------------
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
            raise IPyradError("no samples ready for step 7")

        # sort samples so the largest is first
        checked_samples.sort(
            key=lambda x: x.stats.reads_consens,
            reverse=True,
        )
        return checked_samples


    def setup_dirs(self):
        "Create temp h5 db for storing filters and depth variants"

        # make new output directory
        self.data.dirs.outfiles = os.path.join(
            self.data.paramsdict["project_dir"],
            "{}_outfiles".format(self.data.name),
            )
        if os.path.exists(self.data.dirs.outfiles):
            shutil.rmtree(self.data.dirs.outfiles)
        if not os.path.exists(self.data.dirs.outfiles):
            os.makedirs(self.data.dirs.outfiles)

        # make tmpdir directory
        self.data.tmpdir = os.path.join(
            self.data.paramsdict["project_dir"],
            "{}_outfiles".format(self.data.name),
            "tmpdir",
            )
        if os.path.exists(self.data.tmpdir):
            shutil.rmtree(self.data.tmpdir)
        if not os.path.exists(self.data.tmpdir):
            os.makedirs(self.data.tmpdir)

        # make new database file
        self.data.database = os.path.join(
            self.data.dirs.outfiles,
            self.data.name + ".hdf5",
            )


    def init_vars(self):
        "generate empty h5 database to be filled and setup chunk sizes"   

        # count number of loci
        self.rawloci = os.path.join(
            self.data.dirs.across, 
            self.data.name + "_raw_loci.fa")
        with open(self.rawloci, 'r') as inloci:
            self.nraws = sum(1 for i in inloci if i == "//\n") // 2

        # chunk to approximately 2 chunks per core
        self.ncpus = len(self.ipyclient.ids)
        self.chunks = ((self.nraws // (self.ncpus * 2)) + \
                       (self.nraws % (self.ncpus * 2)))


    def init_database(self):

        # open new database file handle
        with h5py.File(self.data.database, 'w') as io5:

            # attributes
            io5.attrs["samples"] = [i.name.encode() for i in self.samples]
            io5.attrs["filters"] = [b"duplicates", b"indels", b"alleles"]

            # arrays
            io5.create_dataset(
                name="edges", 
                shape=(self.nraws, 5),
                dtype=np.uint16,
                chunks=(self.chunks, 5),
                compression='gzip')
            io5.create_dataset(
                name="chroms", 
                shape=(self.nraws, 3), 
                dtype=np.int64, 
                chunks=(self.chunks, 3),
                compression="gzip")
            io5.create_dataset(
                name="filters", 
                shape=(self.nraws, 3), 
                dtype=np.int64, 
                chunks=(self.chunks, 3),
                compression="gzip")

            # superseqs array
            #io5.create_dataset()


    def iterator(self):
        "..."
        
        # data has to be entered in blocks
        clusters = open(self.rawloci, 'r')
        pairdealer = izip(*[iter(clusters)] * 2)

        # iterate over clusters
        done = 0
        iloc = 0
        cloc = 0
        chunkseqs = np.zeros(self.chunks, dtype="S1")
        chunkedge = np.zeros(self.chunks, dtype=np.uint16)

        while 1:
            try:
                done, chunk = clustdealer(pairdealer, 1)
            except IndexError:
                raise IPyradError("rawloci formatting error in %s", chunk)
            print(chunk)




            # if chunk is full put into superseqs and reset counter
            if cloc == self.chunks:
                superseqs[iloc - cloc:iloc] = chunkseqs
                splits[iloc - cloc:iloc] = chunkedge
                ## reset chunkseqs, chunkedge, cloc
                cloc = 0
                chunkseqs = np.zeros((chunksize, len(samples), maxlen), dtype="S1")
                chunkedge = np.zeros((chunksize), dtype=np.uint16)

            ## get seq and split it
            if chunk:
                try:
                    fill = np.zeros((len(samples), maxlen), dtype="S1")
                    fill.fill("N")
                    piece = chunk[0].decode().strip().split("\n")
                    names = piece[0::2]
                    seqs = np.array([list(i) for i in piece[1::2]])
                    
                    ## fill in the separator if it exists
                    separator = np.where(np.all(seqs == 'n', axis=0))[0]
                    if np.any(separator):
                        chunkedge[cloc] = separator.min()

                    # fill in the hits
                    # seqs will be (5,) IF the seqs are variable lengths, which 
                    # can happen if it had duplicaes AND there were indels, and 
                    # so the indels did not get aligned
                    try:
                        shlen = seqs.shape[1]
                    except IndexError as inst:
                        shlen = min([len(x) for x in seqs])

                    for name, seq in zip(names, seqs):
                        sidx = snames.index(name.rsplit("_", 1)[0])
                        #fill[sidx, :shlen] = seq[:maxlen]
                        fill[sidx, :shlen] = seq[:shlen]

                    ## PUT seqs INTO local ARRAY
                    chunkseqs[cloc] = fill

                except Exception as inst:
                    LOGGER.info(inst)
                    LOGGER.info("\nfill: %s\nshlen %s\nmaxlen %s", fill.shape, shlen, maxlen)
                    LOGGER.info("dupe chunk \n{}".format("\n".join(chunk)))

                ## increase counters if there was a chunk
                cloc += 1
                iloc += 1
            if done:
                break

        ## write final leftover chunk
        superseqs[iloc - cloc:, ] = chunkseqs[:cloc]
        splits[iloc - cloc:] = chunkedge[:cloc]

        ## close super
        io5.close()
        clusters.close()

        ## edges is filled with splits for paired data.
        LOGGER.info("done filling superseqs")

        ## close handle
        #os.remove(infile)







    ## run functions -------------------------------------
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

    # trim edges based on 



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


