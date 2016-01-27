#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0212

import os
import sys
import gzip
import h5py
import random
import shutil
import itertools
import subprocess
import numpy as np
import pandas as pd
from util import *
from ipyrad.assemble.cluster_within import muscle_call, parsemuscle


import logging
LOGGER = logging.getLogger(__name__)



def muscle_align_across(args):
    """ aligns clusters and fills a tmparray with indels"""
    ## parse args, names are used to order arrays by taxon names
    data, samples, chunk = args
    snames = [sample.name for sample in samples]

    ## data are already chunked, read in the whole thing
    infile = open(chunk, 'rb')
    clusts = infile.read().split("//\n//\n")[:-1]
    out = []

    ## array to store indel information 
    maxlen = data._hackersonly["max_fragment_length"]
    if 'pair' in data.paramsdict["datatype"]:
        maxlen *= 2
    indels = np.zeros((len(clusts), len(samples), maxlen), dtype=np.bool)

    ## iterate over clusters and align
    loc = 0
    for loc, clust in enumerate(clusts):
        stack = []
        lines = clust.strip().split("\n")
        names = [i.split()[0][1:] for i in lines]
        seqs = [i.split()[1] for i in lines]

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+";*"+str(i) for i, j in enumerate(names)]

        ## don't bother aligning singletons
        if len(names) <= 1:
            if names:
                stack = [names[0]+"\n"+seqs[0]]
        else:
            ## split seqs before align if PE. If 'nnnn' not found (single end 
            ## or merged reads) then `except` will pass it to SE alignment. 
            try:
                seqs1 = [i.split("nnnn")[0] for i in seqs] 
                seqs2 = [i.split("nnnn")[1] for i in seqs]

                string1 = muscle_call(data, names, seqs1)
                string2 = muscle_call(data, names, seqs2)
                anames, aseqs1 = parsemuscle(string1)
                anames, aseqs2 = parsemuscle(string2)

                ## resort so they're in same order
                aseqs = [i+"nnnn"+j for i, j in zip(aseqs1, aseqs2)]
                for i in xrange(len(anames)):
                    stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])
                    ## store the indels and separator regions as indels
                    locinds = np.array(list(aseqs[i])) == "-"
                    locinds += np.array(list(aseqs[i])) == "n"
                    sidx = [snames.index(anames[i].rsplit("_", 1)[0])]
                    indels[loc, sidx, :locinds.shape[0]] = locinds

            except IndexError:
                string1 = muscle_call(data, names, seqs)
                anames, aseqs = parsemuscle(string1)
                for i in xrange(len(anames)):
                    stack.append(anames[i].rsplit(';', 1)[0]+"\n"+aseqs[i])                    
                    ## store the indels
                    locinds = np.array(list(aseqs[i])) == "-"
                    sidx = snames.index(anames[i].rsplit("_", 1)[0])
                    indels[loc, sidx, :locinds.shape[0]] = locinds

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    infile.close()
    with open(chunk, 'wb') as outfile:
        outfile.write("\n//\n//\n".join(out)+"\n")#//\n//\n")

    return chunk, indels, loc+1



def multi_muscle_align(data, samples, nloci, clustbits, ipyclient):
    """ Splits the muscle alignment across nthreads processors, each runs on 
    clustbit clusters at a time. This is a kludge until I find how to write a 
    better wrapper for muscle. 
    """
    ## create loaded 
    lbview = ipyclient.load_balanced_view()

    try: 
        ## create job queue for clustbits
        submitted_args = []
        for fname in clustbits:
            submitted_args.append([data, samples, fname])

        ## run muscle on all tmp files            
        results = lbview.map_async(muscle_align_across, submitted_args)
        ## return tuple of (chunkname, indelarray, nloci)
        indeltups = results.get()
        ## sort into input order by chunk names
        indeltups.sort(key=lambda x: int(x[0].rsplit("_", 1)[1]))

        ## get dims for full indel array
        maxlen = data._hackersonly["max_fragment_length"]
        if 'pair' in data.paramsdict["datatype"]:
            maxlen *= 2

        ## INIT INDEL ARRAY
        ## build an indel array for ALL loci in cat.clust.gz
        ioh5 = h5py.File(os.path.join(
                            data.dirs.consens, data.name+".indels"), 'w')
        iset = ioh5.create_dataset("indels", (nloci, len(samples), maxlen),
                                   dtype=np.bool)
                                   #chunks=(nloci/10, len(samples), maxlen),
                                   #compression="gzip")

        ## enter all tmpindel arrays into full indel array
        for tup in indeltups:
            start = int(tup[0].rsplit("_", 1)[1])
            iset[start:start+tup[2]] = tup[1]
        ioh5.close()

        ## concatenate finished reads into a tmp file 
        ## TODO: Remove this output... no longer necessary. Tho it provides
        ## a decent sanity check.
        outhandle = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
        with gzip.open(outhandle, 'wb') as out:
            for fname in clustbits:
                with open(fname) as infile:
                    out.write(infile.read()+"//\n//\n")

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        ## still delete tmpfiles if job was interrupted. 
        for path in ["_cathaps.tmp", "_catcons.tmp", ".utemp", ".htemp"]:
            fname = os.path.join(data.dirs.consens, data.name+path)
            if os.path.exists(path):
                pass
                #os.remove(path)
        tmpdir = os.path.join(data.dirs.consens, "tmpaligns")
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)



def cluster(data, noreverse):
    """ calls vsearch for clustering across samples. """
    ## input and output file handles
    cathaplos = os.path.join(data.dirs.consens, data.name+"_cathaps.tmp")
    uhaplos = os.path.join(data.dirs.consens, data.name+".utemp")
    hhaplos = os.path.join(data.dirs.consens, data.name+".htemp")

    ## parameters that vary by datatype 
    ## (too low of cov values yield too many poor alignments)
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .60 "
    elif data.paramsdict["datatype"] == "pairgbs":
        reverse = " -strand both "
        cov = " -query_cov .90 "
    else:
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    ## get call string. Thread=0 means all. 
    ## old userfield: -userfields query+target+id+gaps+qstrand+qcov" \
    cmd = data.bins.vsearch+\
        " -cluster_smallmem "+cathaplos \
       +reverse \
       +cov \
       +" -id "+str(data.paramsdict["clust_threshold"]) \
       +" -userout "+uhaplos \
       +" -notmatched "+hhaplos \
       +" -userfields query+target+qstrand" \
       +" -maxaccepts 1" \
       +" -maxrejects 0" \
       +" -minsl 0.5" \
       +" -fasta_width 0" \
       +" -threads 0" \
       +" -fulldp " \
       +" -usersort " \
       +" -quiet"

    try:
        LOGGER.debug(cmd)
        subprocess.check_call(cmd, shell=True, 
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



def build_h5_array(data, samples, nloci):
    """ build full catgs file """
    ## catg array of prefiltered loci (4-dimensional) aye-aye!
    ## this can be multiprocessed!! just sum arrays at the end
    ## but one big array will overload memory, so need to be done in slices
    ## or better yet, using dask to automatically chunk it up...

    ## from dask.distributed import dask_client_from_ipclient
    ## dclient = dask_client_from_ipclient(ipclient)

    ## sort to ensure samples will be in alphabetical order, tho they should be.
    samples.sort(key=lambda x: x.name)

    ## initialize an hdf5 array of the super catg matrix with dims
    maxlen = data._hackersonly["max_fragment_length"]
    if 'pair' in data.paramsdict["datatype"]:
        maxlen *= 2

    ## open h5py handle
    ioh5 = h5py.File(data.database, 'w')

    ## INIT FULL CATG ARRAY
    ## store catgs with a .10 loci chunk size
    supercatg = ioh5.create_dataset("catgs", (nloci, len(samples), maxlen, 4),
                                    dtype=np.uint32)
                                    #chunks=(nloci/10, len(samples), maxlen, 4),
                                    #compression="gzip")
    supercatg.attrs["samples"] = [i.name for i in samples]
    supercatg.attrs["chunksize"] = 1000

    ## INIT FULL SEQS ARRAY
    ## array for clusters of consens seqs
    superseqs = ioh5.create_dataset("seqs", (nloci, len(samples), maxlen),
                                     dtype="|S1")
    superseqs.attrs["samples"] = [i.name for i in samples]
    superseqs.attrs["chunksize"] = 1000    

    ## INIT FULL SNPS ARRAY
    ## array for snp string, 2 cols, - and *
    snps = ioh5.create_dataset("snps", (nloci, maxlen, 2), dtype=np.bool)
    snps.attrs["names"] = ["-", "*"]

    ## INIT FULL FILTERS ARRAY
    ## array for filters that will be applied in step7
    filters = ioh5.create_dataset("filters", (nloci, 6), dtype=np.bool)
    filters.attrs["filters"] = ["duplicates", "max_indels", "max_snps", 
                                "max_hets", "min_samps", "bad_edges"]
    filters.attrs["chunks"] = 1000        

    ## INIT FULL EDGE ARRAY
    ## array for edgetrimming 
    edges = ioh5.create_dataset("edges", (nloci, 5), dtype=np.uint16)
    edges.attrs["names"] = ["R1_L", "R1_R", "R2_L", "R2_R", "sep"]

    ## RUN SINGLECAT, FILL FILTERS
    ## for each sample fill its own hdf5 array with catg data & indels. 
    ## maybe this can be parallelized. Can't right now since we pass it 
    ## an open file object (indels). Room for speed improvements, tho.
    ipath = os.path.join(data.dirs.consens, data.name+".indels")
    with h5py.File(ipath, 'r') as indh5:
        for sidx, sample in enumerate(samples):
            ## loads the array for one sample into memory. 
            ## TODO: Chunk to allow efficient reading along y axis. 
            indels = indh5["indels"][:, sidx, :]
            ## return which loci were filtered b/c of duplicates.
            dupfilter, indfilter = singlecat(data, sample, nloci, indels)
            filters[:, 0] += dupfilter
            filters[:, 1] += indfilter

    ## FILL SUPERCATG -- TODO: can still parallelize singlecat.
    ## combine indvidual hdf5 arrays into supercatg
    h5handles = []
    for sidx, sample in enumerate(samples):
        h5handle = os.path.join(data.dirs.consens, sample.name+".hdf5")
        h5handles.append(h5handle)
        ## open, read, and close individual data base
        with h5py.File(h5handle, 'r') as singleh5:
            icatg = singleh5[sample.name]
            supercatg[:, sidx, :, :] = icatg

    ## FILL SUPERSEQS
    fill_superseqs(data, samples, superseqs)

    ## close the big boy
    ioh5.close()

    ## clean up / remove individual catg files
    for handle in h5handles:
        os.remove(handle)
    ## remove indels array
    os.remove(ipath)

    ## set sample states
    for sample in samples:
        sample.stats.state = 6



def singlecat(data, sample, nloci, indels):
    """ 
    Orders catg data for each sample into the final locus order. This allows
    all of the individual catgs to simply be combined later. They are also in 
    the same order as the indels array, so indels are inserted from the indel
    array that in passed in. 
    """
    LOGGER.debug("singlecat: %s", sample.name)

    ## if an hdf5 file already exists delete it.
    h5handle = os.path.join(data.dirs.consens, sample.name+".hdf5")
    if os.path.exists(h5handle):
        os.remove(h5handle)

    ## INIT SINGLE CATG ARRAY
    ## create an h5 array for storing catg infor for this sample
    new_h5 = h5py.File(h5handle, 'w')
    maxlen = data._hackersonly["max_fragment_length"]
    if 'pair' in data.paramsdict["datatype"]:
        maxlen *= 2

    #nloci = sample.stats.clusters_hidepth
    icatg = new_h5.create_dataset(sample.name, (nloci, maxlen, 4), 
                                  dtype=np.uint32)#,
                                  #chunks=(nloci/10, maxlen, 4))#,
                                  #compression="gzip")

    ## LOAD IN STEP5 CATG ARRAY
    ## get catg from step5 for this sample, the shape is (nloci, maxlen)
    old_h5 = h5py.File(sample.files.database, 'r')
    catarr = old_h5["catg"][:]

    ## get utemp cluster hits as pandas data frame
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)
    ## add name and index columns to dataframe
    updf.loc[:, 3] = [i.rsplit("_", 1)[0] for i in updf[0]]
    ## create a column with only consens index
    updf.loc[:, 4] = [i.rsplit("_", 1)[1] for i in updf[0]]

    ## for each locus in the udic (groups) ask if sample is in it
    dupfilter = np.zeros(nloci, dtype=np.bool)
    indfilter = np.zeros(nloci, dtype=np.bool)
    udic = updf.groupby(by=1, sort=False)

    for iloc, seed in enumerate(udic.groups.iterkeys()):
        ipdf = udic.get_group(seed)
        ask = ipdf.where(ipdf[3] == sample.name).dropna()
        if ask.shape[0] == 1: 
            ## if multiple hits of a sample to a locus then it is not added
            ## to the locus, and instead the locus is masked for exclusion
            ## using the filters array.
            icatg[iloc] = catarr[int(ask[4]), :icatg.shape[1], :]
        elif ask.shape[0] > 1:
            ## store that this cluster failed b/c it had duplicate samples. 
            dupfilter[iloc] = True

    ## for each locus in which Sample was the seed
    seedmatch1 = (sample.name in i for i in udic.groups.keys())
    seedmatch2 = (i for i, j in enumerate(seedmatch1) if j)
    for iloc in seedmatch2:
        icatg[iloc] = catarr[iloc, :icatg.shape[1], :]

    ## close the old hdf5 connections
    old_h5.close()

    ## insert indels into new_h5 (icatg array) which now has loci in the same
    ## order as the final clusters from the utemp file
    for iloc in xrange(icatg.shape[0]):
        ## indels locations
        indidx = np.where(indels[iloc, :])[0]
        if np.any(indidx):
            ### apply indel filter 
            #ind1 = len(indidx1) <= data.paramsdict["max_Indels_locus"][0]
            #ind2 = len(indidx2) <= data.paramsdict["max_Indels_locus"][1]
            #if ind1 and ind2:
            #    indfilter[iloc] = True
            #    print(iloc, ind1, ind2, icatg.shape[0])
            if len(indidx) < sum(data.paramsdict["max_Indels_locus"]):
                indfilter[iloc] = True

            ## insert indels into catg array
            newrows = icatg[iloc].shape[0] + len(indidx)
            not_idx = np.array([k for k in range(newrows) if k not in indidx])
            ## create an empty new array with the right dims
            newfill = np.zeros((newrows, 4), dtype=icatg.dtype)
            ## fill with the old vals leaving the indels rows blank
            newfill[not_idx, :] = icatg[iloc]
            ## store new data into icatg
            icatg[iloc] = newfill[:icatg[iloc].shape[0]]

    ## close the new h5 that was written to
    new_h5.close()

    ## returns tmpfilter 
    return dupfilter, indfilter



def fill_superseqs(data, samples, superseqs):
    """ fills the superseqs array with seq data from cat.clust """
    ## samples are already sorted
    snames = [i.name for i in samples]
    ## get maxlen again
    maxlen = data._hackersonly["max_fragment_length"]
    if 'pair' in data.paramsdict["datatype"]:
        maxlen *= 2

    ## data has to be entered in blocks
    infile = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
    clusters = gzip.open(infile, 'r')
    pairdealer = itertools.izip(*[iter(clusters)]*2)

    ## iterate over clusters
    done = 0
    iloc = 0
    while not done:
        try:
            done, chunk = clustdealer(pairdealer, 1)
        except IndexError:
            raise IPyradError("clustfile formatting error in %s", chunk)    
        ## input array must be this shape
        if chunk:
            fill = np.zeros((len(samples), maxlen), dtype="|S1")
            piece = chunk[0].strip().split("\n")
            names = piece[0::2]
            seqs = np.array([list(i) for i in piece[1::2]])

            ## fill in the hits
            indices = range(len(snames))
            shlen = seqs.shape[1]
            for name, seq in zip(names, seqs):
                sidx = snames.index(name.rsplit("_", 1)[0])
                fill[sidx, :shlen] = seq
                indices.remove(sidx)

            ## fill in the misses
            for idx in indices:
                fill[idx] = np.array(["N"]*maxlen)

            ## PUT INTO ARRAY
            superseqs[iloc] = fill

        ## increase counter
        iloc += 1
    ## close handle
    clusters.close()



def build_reads_file(data):
    """ reconstitutes clusters from .utemp and htemp files and writes them 
    to chunked files for aligning in muscle. Return a dictionary with 
    seed:hits info from utemp file.
    """
    ## read in cluster hits as pandas data frame
    uhandle = os.path.join(data.dirs.consens, data.name+".utemp")
    updf = pd.read_table(uhandle, header=None)

    ## load full fasta file into a Dic
    conshandle = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    consdf = pd.read_table(conshandle, delim_whitespace=1, header=None)
    printstring = "{:<%s}    {}" % max([len(i) for i in set(consdf[0])])
    consdic = consdf.set_index(0)[1].to_dict()

    ## make an tmpout directory and a printstring for writing to file
    tmpdir = os.path.join(data.dirs.consens, "tmpaligns")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    ## groupby index 1 (seeds) 
    groups = updf.groupby(by=1, sort=False)

    ## a chunker for writing every N
    optim = 100
    #if len(groups) > 2000:
    #    optim = len(groups) // 10

    ## get seqs back from consdic
    clustbits = []
    locilist = []
    loci = 0

    ## iterate over seeds and add in hits seqs
    for seed in set(updf[1].values):
        ## get dataframe for this locus/group (gdf) 
        gdf = groups.get_group(seed)
        ## set seed name and sequence
        seedseq = consdic.get(">"+seed)
        names = [">"+seed]
        seqs = [seedseq]
        ## iterate over group dataframe (gdf) and add hits names and seqs
        ## revcomp the hit if not '+' in the df.
        for i in gdf.index:
            hit = gdf[0][i]
            names.append(">"+hit)
            if gdf[2][i] == "+":
                seqs.append(consdic[">"+hit])
            else:
                seqs.append(fullcomp(consdic[">"+hit][::-1]))

        ## append the newly created locus to the locus list
        locilist.append("\n".join([printstring.format(i, j) \
                        for i, j in zip(names, seqs)]))
        loci += 1

        ## if enough loci have finished write to file to clear the mem
        if not loci % optim:
            ## a file to write results to
            handle = os.path.join(tmpdir, "tmp_"+str(loci - optim))
            with open(handle, 'w') as tmpout:
                tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
                locilist = []
                clustbits.append(handle)
    ## write the final remaining to file
    if locilist:
        handle = os.path.join(tmpdir, "tmp_"+str(loci - len(locilist)))
        with open(handle, 'w') as tmpout:
            tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
            clustbits.append(handle)

    ## return stuff
    return clustbits, loci



def build_input_file(data, samples, outgroups, randomseed):
    """ Make the concatenated consens files to input to vsearch. 
    Orders reads by length and shuffles randomly within length classes"""

    ##  make a copy list that will not have outgroups excluded
    conshandles = list(itertools.chain(*[samp.files.consens \
                                         for samp in samples]))
    ## remove outgroup sequences, add back in later to bottom after shuffling
    ## outgroups could be put to end of sorted list
    conshandles.sort()
    assert conshandles, "no consensus files found"
                
    ## output file for consens seqs from all taxa in consfiles list
    allcons = os.path.join(data.dirs.consens, data.name+"_catcons.tmp")
    allhaps = open(allcons.replace("_catcons.tmp", "_cathaps.tmp"), 'w')

    ## combine cons files and write as pandas readable format to all file
    ## this is the file that will be read in later to build clusters
    printstr = "{:<%s}    {}" % 100   ## max len allowable name
    with open(allcons, 'wb') as consout:
        for qhandle in conshandles:
            with gzip.open(qhandle, 'r') as tmpin:
                data = tmpin.readlines()
                names = iter(data[::2])
                seqs = iter(data[1::2])
                consout.write("".join(printstr.format(i.strip(), j) \
                              for i, j in zip(names, seqs)))

    ## created version with haplos that is also shuffled within seqlen classes
    random.seed(randomseed)
    ## read back in cons as a data frame
    consdat = pd.read_table(allcons, delim_whitespace=1, header=None)
    ## make a new column with seqlens and then groupby seqlens
    consdat[2] = pd.Series([len(i) for i in consdat[1]], index=consdat.index)
    lengroups = consdat.groupby(by=2)
    lenclasses = sorted(set(consdat[2]), reverse=True)
    ## write all cons in pandas readable format
    for lenc in lenclasses:
        group = lengroups.get_group(lenc)
        ## shuffle the subgroup
        shuf = group.reindex(np.random.permutation(group.index))
        ## convert ambiguity codes into a sampled haplotype for any sample
        ## to use for clustering, but ambiguities are still saved in allcons
        writinghaplos = []
        for ind in shuf.index:
            ## TODO: Is this too much in memory for super huge data sets?
            ## may need to be chunked.
            writinghaplos.append("\n".join([shuf[0][ind], 
                                            splitalleles(shuf[1][ind])[0]]))
        allhaps.write("\n".join(writinghaplos)+"\n")
    allhaps.close()



def run(data, samples, noreverse, force, randomseed, ipyclient):
    """ subselect and pass args for across-sample clustering """

    ## clean the slate
    if os.path.exists(data.database):
        os.remove(data.database)

    ## make file with all samples reads, allow priority to shunt outgroups
    ## to the end of the file
    outgroups = []
    LOGGER.info("creating input files")
    build_input_file(data, samples, outgroups, randomseed)

    ## call vsearch
    LOGGER.info("clustering")    
    cluster(data, noreverse)

    ## build consens clusters and returns chunk handles to be aligned
    LOGGER.info("building consens clusters")        
    clustbits, nloci = build_reads_file(data)

    ## muscle align the consens reads and creates hdf5 indel array
    LOGGER.info("muscle alignment & building indel database")
    multi_muscle_align(data, samples, nloci, clustbits, ipyclient)

    ## builds the final HDF5 array which includes three main keys
    ## /catg -- contains all indiv catgs and has indels inserted
    ##   .attr['samples'] = [samples]
    ## /filters -- empty for now, will be filled in step 7
    ##   .attr['filters'] = [f1, f2, f3, f4, f5]
    ## /seqs -- contains the clustered sequence data as string arrays
    ##   .attr['samples'] = [samples]
    LOGGER.info("building full database")    
    ## calls singlecat func inside
    build_h5_array(data, samples, nloci)

    ## invarcats()
    ## invarcats()


if __name__ == "__main__":
    ## test...
    import ipyrad as ip

    ## get path to test dir/ 
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__)
               )
           )
       )

    ## run test on pairgbs data1
    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "Ron", "Ron"))
    TEST.step6(force=True)
    print(TEST.stats)

    ## run test on pairgbs data1
    # TEST = ip.load.load_assembly(os.path.join(\
    #                      ROOT, "tests", "test_pairgbs", "test_pairgbs"))
    # TEST.step6(force=True)
    # print(TEST.stats)

    ## run test on rad data1
    TEST = ip.load.load_assembly(os.path.join(\
                         ROOT, "tests", "test_rad", "data1"))
    TEST.step6(force=True)
    print(TEST.stats)

    ## load test data (pairgbs)
    # DATA = ip.load.load_assembly(\
    #     "/home/deren/Documents/ipyrad/tests/test_pairgbs/test_pairgbs.assembly")
    # ## run step 6
    # DATA.step6(force=True)


