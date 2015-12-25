#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101

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
from collections import OrderedDict
from .util import *
from ipyrad.assemble.cluster_within import muscle_call, parsemuscle


import logging
LOGGER = logging.getLogger(__name__)


def breakalleles(consensus):
    """ break ambiguity code consensus seqs into two alleles """
    allele1 = ""
    allele2 = ""
    bigbase = ""
    for base in consensus:
        if base in tuple("RKSYWM"):
            unhet1, unhet2 = unhetero(base)
            hetset = set([unhet1, unhet2])
            allele1 += uplow((unhet1, unhet2))
            allele2 += hetset.difference(uplow((unhet1, unhet2))).pop()
            if not bigbase:
                bigbase = uplow((allele1, allele2))
        elif base in tuple("rksywm"):
            unhet1, unhet2 = unhetero(base)
            hetset = set([unhet1, unhet2])
            allele2 += uplow((unhet1, unhet2))
            allele1 += hetset.difference(uplow((unhet1, unhet2))).pop()
        else:
            allele1 += base
            allele2 += base
    return allele1, allele2



def fullcomp(seq):
    """ returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences"""
    ## this is surely not the most efficient...
    seq = seq.replace("A", 'u')\
             .replace('T', 'v')\
             .replace('C', 'p')\
             .replace('G', 'z')\
             .replace('u', 'T')\
             .replace('v', 'A')\
             .replace('p', 'G')\
             .replace('z', 'C')

    seq = seq.replace('R', 'u')\
             .replace('Y', 'v')\
             .replace('K', 'p')\
             .replace('M', 'z')\
             .replace('u', 'Y')\
             .replace('v', 'R')\
             .replace('p', 'M')\
             .replace('z', 'K')

    seq = seq.replace('r', 'u')\
             .replace('y', 'v')\
             .replace('k', 'p')\
             .replace('m', 'z')\
             .replace('u', 'y')\
             .replace('v', 'r')\
             .replace('p', 'm')\
             .replace('z', 'k')
    return seq



def muscle_align_across(args):
    """ aligns reads, does split then aligning for paired reads """
    ## parse args
    data, samples, chunk = args
    snames = [sample.name for sample in samples]

    ## data are already chunked, read in the whole thing
    infile = open(chunk, 'rb')
    clusts = infile.read().split("//\n//\n")[:-1]
    out = []
    ## array to store indel information
    maxlen = [225 if 'pair' in data.paramsdict["datatype"] else 125][0]
    indels = np.zeros((len(clusts), len(samples), maxlen), dtype=np.int8)

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
            ## split seqs if paired end seqs
            if 'pair' in data.paramsdict["datatype"]:
                seqs1 = [i.split("ssss")[0] for i in seqs] 
                seqs2 = [i.split("ssss")[1] for i in seqs]

                string1 = muscle_call(data, names[:200], seqs1[:200])
                string2 = muscle_call(data, names[:200], seqs2[:200])
                anames, aseqs1 = parsemuscle(string1)
                anames, aseqs2 = parsemuscle(string2)

                ## resort so they're in same order
                aseqs = [i+"ssss"+j for i, j in zip(aseqs1, aseqs2)]
                somedic = OrderedDict()
                for i in range(len(anames)):
                    somedic[anames[i].rsplit(';', 1)[0]] = aseqs[i]
                    locinds = np.array(aseqs[i] == "-", dtype=np.int32)
                    indels[loc][snames.index(anames[i]\
                                .rsplit("_", 1)[0])][:maxlen] = locinds

            else:
                string1 = muscle_call(data, names[:200], seqs[:200])
                anames, aseqs = parsemuscle(string1)
                somedic = OrderedDict()
                for i in range(len(anames)):
                    somedic[anames[i].rsplit(";", 1)[0]] = aseqs[i]
                    locinds = np.array(aseqs[i] == "-", dtype=np.int32)
                    indels[loc][snames.index(anames[i]\
                                .rsplit("_", 1)[0])][:maxlen] = locinds

            for key in somedic.iterkeys():
                stack.append(key+"\n"+somedic[key])
            #for key in keys:
            #    if key[-1] == '-':
            #        ## reverse matches should have --- overhang
            #        if set(somedic[key][:4]) == {"-"}:
            #            stack.append(key+"\n"+somedic[key])
            #        else:
            #            pass ## excluded
            #    else:
            #        stack.append(key+"\n"+somedic[key])

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    infile.close()
    outfile = open(chunk, 'wb')#+"_tmpout_"+str(num))
    outfile.write("\n//\n//\n".join(out)+"\n")#//\n//\n")
    outfile.close()

    return chunk, indels, loc+1



def multi_muscle_align(data, samples, clustbits, ipyclient):
    """ Splits the muscle alignment across nthreads processors, each runs on 
    1000 clusters at a time. This is a kludge until I find how to write a 
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
        indeltups = results.get()
        ## concatenate indel arrays in correct order
        nloci = 1000
        maxlen = [225 if 'pair' in data.paramsdict["datatype"] else 125][0]
        ioh5 = h5py.File(os.path.join(
                            data.dirs.consens, data.name+"indels"), 'w')
        dset = ioh5.create_dataset("indels", (nloci, len(samples), maxlen),
                                   dtype='i4', 
                                   chunks=(nloci/10, len(samples), maxlen),
                                   compression="gzip")
        ## sort into input order
        indeltups.sort(key=lambda x: int(x[0].rsplit("_", 1)[1]))
        for tup in indeltups:
            #print(tup[0][-10:], tup[1].shape, tup[2])
            start = int(tup[0].rsplit("_", 1)[1])
            dset[start:start+tup[2]] = tup[1]
        ioh5.close()

        ## concatenate finished reads
        outhandle = os.path.join(data.dirs.consens, data.name+"_catclust.gz")
        with gzip.open(outhandle, 'wb') as out:
            for fname in clustbits:
                with open(fname) as infile:
                    out.write(infile.read()+"//\n//\n")
        #return nloci

    except Exception as inst:
        LOGGER.warn(inst)
        raise

    finally:
        ## still delete tmpfiles if job was interrupted
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
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .35 "
    elif data.paramsdict["datatype"] == "pairgbs":
        reverse = " -strand both "
        cov = " -query_cov .60 "
    else:
        reverse = " -leftjust "
        cov = " -query_cov .90 "

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        print(noreverse, "not performing reverse complement clustering")

    ## get call string
    cmd = data.bins.vsearch+\
        " -cluster_smallmem "+cathaplos \
       +reverse \
       +cov \
       +" -id "+str(data.paramsdict["clust_threshold"]) \
       +" -userout "+uhaplos \
       +" -notmatched "+hhaplos \
       +" -userfields query+target+id+gaps+qstrand+qcov" \
       +" -maxaccepts 1" \
       +" -maxrejects 0" \
       +" -minsl 0.5" \
       +" -fasta_width 0" \
       +" -threads 0" \
       +" -fulldp " \
       +" -usersort " \
       +" -quiet"

    try:
        subprocess.check_call(cmd, shell=True, 
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



def build_catg_file(data, samples, udic):
    """ build full catgs file """
    ## catg array of prefiltered loci (4-dimensional) aye-aye!
    ## this can be multiprocessed!! just sum arrays at the end
    ## but one big array will overload memory, so need to be done in slices
    ## or better yet, using dask to automatically chunk it up...

    ## from dask.distributed import dask_client_from_ipclient
    ## dclient = dask_client_from_ipclient(ipclient)

    ## for each catg sample fill its own hdf5
    for sample in samples:
        singlecat(data, sample, udic)

    ## initialize an hdf5 array of the super catg matrix
    maxlen = 125
    nloci = 1000
    ioh5 = h5py.File(data.database, 'w')
    ## probably have to do something better than .10 loci chunk size
    supercatg = ioh5.create_dataset("catgs", (nloci, len(samples), maxlen, 4),
                                    dtype='i4', 
                                    chunks=(nloci/10, len(samples), maxlen, 4),
                                    compression="gzip")

    ## sum individual hdf5 arrays into supercatg



def singlecat(data, sample, udic):
    """ run one samples """
    ## create 
    ioh5 = h5py.File(os.path.join(data.dirs.consens, sample.name+".hdf5"), 'w')
    ## probably have to do something better than .10 loci chunk size
    maxlen = 125
    nloci = 1000
    icatg = ioh5.create_dataset(sample.name, (nloci, maxlen, 4),
                                dtype='i4',) 
                                #chunks=(nloci/10, len(samples), maxlen, 4),
                                #compression="gzip")

    ## get catg for one sample
    catarr = np.load(sample.files.database)
    #catarr = np.load(os.path.join(data.dirs.consens, sample.name+".catg"))
    ## for each read in catarr, if it g ask which locus it is in

    ## gonna want to 'filter' the pdf, e.g., :
    ## dff.groupby('B').filter(lambda x: len(x) > 2)
    ## groups.get_group('3J_0_594').apply(lambda x: x[0].rsplit("_", 1)[0], axis=1)
    for iloc, seed in enumerate(udic.groups.iterkeys()):
        ipdf = udic.get_group(seed)

    ## for each locus in full data set
    itax = 10 ## get this
    nloci = 1000
    for iloc in xrange(nloci):
        getloc = udic.get_group('...')
        dset[iloc][itax] = catarr[getloc]

    print(dset.name, dset.shape)
    print(catarr.shape)




def build_reads_file(data):
    """ reconstitutes clusters from .utemp and htemp files and writes them 
    to chunked files for aligning in muscle. Return a dictionary with 
    seed:hits info from utemp file
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
    optim = 100

    ## groupby index 1 (seeds) 
    groups = updf.groupby(by=1, sort=False)

    ## get seqs back from consdic
    clustbits = []
    locilist = []
    loci = 0
    for seed in set(updf[1].values):
        ## get sub-DF
        gdf = groups.get_group(seed)
        ## set seed name and sequence
        seedseq = consdic.get(">"+seed)
        names = [">"+seed]
        seqs = [seedseq]
        ## iterate over gdf
        for i in gdf.index:
            hit = gdf[0][i]
            names.append(">"+hit)
            if gdf[4][i] == "+":
                seqs.append(consdic[">"+hit])
            else:
                seqs.append(fullcomp(hit[::-1]))
        locilist.append("\n".join([printstring.format(i, j) \
                        for i, j in zip(names, seqs)]))

        loci += 1
        if not loci % optim:
            ## a file to write results to
            handle = os.path.join(tmpdir, "tmp_"+str(loci))
            with open(handle, 'w') as tmpout:
                tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
                locilist = []
                clustbits.append(handle)
    if locilist:
        with open(os.path.join(tmpdir, "tmp_"+str(loci)), 'w') as tmpout:
            tmpout.write("\n//\n//\n".join(locilist)+"\n//\n//\n")
            clustbits.append(handle)

    return groups, clustbits



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
            writinghaplos.append(shuf[0][ind]+'\n'+\
                                 breakalleles(shuf[1][ind])[0])
        allhaps.write("\n".join(writinghaplos))
    allhaps.close()



def run(data, samples, noreverse, force, randomseed, ipyclient):
    """ subselect and pass args for across-sample clustering """

    ## make file with all samples reads, allow priority to shunt outgroups
    ## to the end of the file
    outgroups = []
    LOGGER.info("creating input files")
    build_input_file(data, samples, outgroups, randomseed)

    ## call vsearch
    LOGGER.info("clustering")    
    cluster(data, noreverse)

    ## build consens clusters and return hit dict
    LOGGER.info("building consens clusters")        
    ugroups, clustbits = build_reads_file(data)

    ## muscle align the consens reads and creates hdf5 indel array
    LOGGER.info("muscle alignment & building indel database")
    multi_muscle_align(data, samples, clustbits, ipyclient)

    ## build supercatg file and insert indels
    ## this can't do all loci at once, needs chunking @ < (1e6, 100) 
    LOGGER.info("... not yet building full database")    
    #build_catg_file(data, samples, ugroups)

    ## convert full catg into a vcf file
    LOGGER.info("...not yet converting to VCF")        

    ## invarcats()
    ## invarcats()


if __name__ == "__main__":
    ## test...
    import ipyrad as ip

    ## reload autosaved data. In case you quit and came back 
    DATA = ip.load.load_assembly(\
        "/home/deren/Documents/ipyrad/tests/test_rad/data1.assembly")

    ## run step 6
    DATA.step6(force=True)


