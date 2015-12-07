#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

from __future__ import print_function
# pylint: disable=E1101

import os
import sys
import gzip
import copy
import random
import itertools
import subprocess
import numpy as np
import cPickle as pickle
from ipyrad.assemble.consens_se import unhetero, uplow


def breakalleles(consensus):
    """ break ambiguity code consensus seqs
    into two alleles """
    a1 = ""
    a2 = ""
    bigbase = ""
    for base in consensus:
        if base in tuple("RKSYWM"):
            a,b = unhetero(base)
            d = set([a,b])
            a1 += uplow((a,b))
            a2 += d.difference(uplow((a,b))).pop()
            if not bigbase:
                bigbase = uplow((a,b))
        elif base in tuple("rksywm"):
            a,b = unhetero(base)
            d = set([a,b])
            a2 += uplow((a,b))
            a1 += d.difference(uplow((a,b))).pop()
        else:
            a1 += base
            a2 += base
    return a1,a2



def fullcomp(seq):
    """ returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences"""
    ## this is probably not the most efficient...
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



def cluster(data, noreverse, nthreads):
    """ calls vsearch for clustering across samples. """
    ## input and output file handles
    cathaplos = os.path.join(data.dirs.consens, "cat.haplos")
    uhaplos = os.path.join(data.dirs.consens, "cat.utemp")
    hhaplos = os.path.join(data.dirs.consens, "cat.htemp")    

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
    cmd = data.vsearch+\
        " -cluster_smallmem "+cathaplos \
       +reverse \
       +cov \
       +" -id "+str(data.paramsdict["clust_threshold"]) \
       +" -userout "+uhaplos \
       +" -userfields query+target+id+gaps+qstrand+qcov" \
       +" -maxaccepts 1" \
       +" -maxrejects 0" \
       +" -minsl 0.5" \
       +" -fulldp " \
       +" -threads "+str(nthreads) \
       +" -usersort "+\
       +" -notmatched "+hhaplos \
       +" -fasta_width 0" \
       +" -quiet"

    try:
        subprocess.check_call(cmd, shell=True, 
                                   stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



def oldcluster(data, handle, gid, quiet):
    """ clusters with vsearch across samples """

    if params["datatype"] == 'pairddrad':
        ## use first files for split clustering "
        if gid:
            ## hierarchical clustering save temps
            notmatched = " -notmatched "+params["work"]+\
                           "prefix/"+handle.split("/")[-1].\
                           replace(".firsts_", "._temp_")
            userout = " -userout "+params["work"]+"prefix/"+\
                      handle.split("/")[-1].replace(".firsts_", ".u_")
        else:
            notmatched = ""
            userout = " -userout "+handle.replace(".firsts_", ".u")
    else:
        ## use haplos files "
        if gid:
            ## hierarchical clustering save temps
            notmatched = " -notmatched "+params["work"]+"prefix/"+\
                         handle.split("/")[-1].replace(".haplos_", "._temp_")
            userout = " -userout "+params["work"]+"prefix/"+\
                        handle.split("/")[-1].replace(".haplos_", ".u_")
        else:
            notmatched = ""
            userout = " -userout "+handle.replace(".haplos_", ".u")

    if params["datatype"] in ['gbs', 'pairgbs', 'merged']:
        parg = " -strand both"
        cov = " -query_cov .60 "
    else:
        parg = " -leftjust "
        cov = " -query_cov .90 "
    if 'vsearch' not in params["vsearch"]:
        mask = ""
        threads = " -threads 1"
    else:
        mask = " -qmask "+params["mask"]
        threads = " -threads 8"
    if quiet:
        quietarg = " -quiet "
    else:
        quietarg = ""

    cmd = params["vsearch"]+\
        " -cluster_smallmem "+handle+\
        parg+\
        mask+\
        threads+\
        " -id "+params["wclust"]+\
        userout+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -fulldp"+\
        " -usersort"+\
        cov+\
        notmatched+\
        quietarg
    subprocess.call(cmd, shell=True)


def makeclust(params, handle, gid, minmatch):
    """ reconstitutes clusters from .u and temp files and
    stores consens file names to cluster numbers and dumps to a dict """

    ## locus counter
    olddict = {}
    locus = 1

    ## read in cluster hits and seeds files
    if not gid:
        userout = open(handle.replace(".haplos_", ".u"), 'r')
        outfile = gzip.open(handle.replace(".haplos_"+gid, \
                                           ".clust_"+gid+".gz"), 'w')
    else:
        userout = open(params["work"]+'prefix/'+handle.split("/")[-1].\
                                       replace(".haplos_", ".u_"), 'r')
        nomatch = open(params["work"]+'prefix/'+handle.split("/")[-1].\
                                       replace(".haplos_", "._temp_"), 'r')
        outfile = open(params["work"]+'prefix/'+handle.split("/")[-1].\
                                       replace(".haplos_", ".seed_"), 'w')

    ## load full fasta file into a Dic
    consdic = {}   ## D
    if params["datatype"] == 'pairddrad':
        if gid:
            infile = open(handle.replace(".haplos_"+gid, ".firsts_"+gid))
        else:
            infile = gzip.open(handle.replace(".haplos_"+gid,
                                              ".consens_"+gid+".gz"))
    else:
        infile = gzip.open(handle.replace(".haplos_"+gid, 
                                          ".consens_"+gid+".gz"))

    duo = itertools.izip(*[iter(infile)]*2)  ## L
    while 1:
        try: 
            name, seq = duo.next()
        except StopIteration:
            break
        consdic[name.strip()] = seq.strip()
    infile.close()

    ## load .u match info into a Dic
    udic = {}
    for line in [line.split("\t") for line in userout.readlines()]:
        if ">"+line[1] in udic:
            udic[">"+line[1]].append([">"+line[0], line[4]])
        else:
            udic[">"+line[1]] = [[">"+line[0], line[4]]]
    
    ## if tier 1 of hierarchical clustering "
    if gid:
        if int(minmatch) == 1:
            ## no reduction, write seeds only "
            singles = nomatch.read().split(">")[1:]
            for i in singles:
                i = i.split("\n")[0]+"\n"+"".join(i.split("\n")[1:]).upper()
                print >>outfile, ">"+i+"\n//"
                #print >>outfile, str(locus)
            del singles
        else:       
            for key, values in udic.items():
                ## reduction, only write seed if minimum hits reached
                if (len(values)+1) >= int(minmatch):
                    ## fix for if short seqs are excluded during clustering
                    if consdic.get(key):
                        seq = key+"\n"+consdic[key]+"\n"
                        seq += "//\n"
                        #seq += str(locus)+"\n"
                        outfile.write(seq)
    else:
        ## map sequences to clust file in order
        seq = ""
        for key, values in udic.items():
            if consdic.get(key):   
                ## fix for if short seqs are excluded during clustering
                seq = key+"\n"+consdic[key]+'\n'
                names = [i[0] for i in values]
                orient = [i[1] for i in values]
                for i in range(len(names)):
                    if consdic.get(names[i]):  
                        if orient[i] == "+":
                            seq += names[i]+'\n'+consdic[names[i]]+"\n"
                        else:
                            seq += names[i]+'\n'+\
                                   fullcomp(consdic[names[i]][::-1])+"\n"
                seq += str(locus)+"\n//\n"
                outfile.write(seq)

                ## store names in locus number
                ## with ">" char trimmed off
                #print 'names,'
                #print names, key
                olddict[locus] = [i[1:] for i in names+[key]]
                #print locus
                #print len(names+[key]), names+[key]
                locus += 1

    outfile.close()
    userout.close()
    if gid:
        nomatch.close()

    ## dump the locus2name map
    pickleout = gzip.open(handle.replace("haplos_", "loc2name.map"), 'wb')
    pickle.dump(olddict, pickleout)
    pickleout.close()


# def splitter(handle):
#     """ splits first reads from second reads for pairddrad data
#         for split clustering method """
#     infile = open(handle)
#     if os.path.exists(handle.replace(".haplos", ".firsts")):
#         os.remove(handle.replace(".haplos", ".firsts"))
        
#     orderfirsts = open(handle.replace(".haplos", ".firsts"), 'w')
#     dp = itertools.izip(*[iter(infile)]*2)
#     ff = []
#     cnts = 0
#     for d in dp:
#         n,s = d
#         ## checking fix to pairddrad splitting problem...
#         ## backwards compatible with pyrad v2
#         s1 = s.replace("X", "x").replace("x", "n").split("nn")[0]
#         ff.append(n+s1+"\n")
#         cnts += 1
#     orderfirsts.write("".join(ff))
#     orderfirsts.close()
#     return handle.replace(".haplos", ".firsts")



def makecons(data, samples, outgroups, randomseed):
    """ Make the concatenated consens files to input to vsearch. 
    Orders reads by length and shuffles randomly within length classes"""

    ##  make a copy list that will not have outgroups excluded
    conshandles = [data.samples[samp].files.consens for samp in samples]

    ## remove outgroup sequences, add back in later to bottom after shuffling "
    ##    outgroup = ""
                
    ## output file for consens seqs from all taxa in consfiles list
    allcons = gzip.open(data.dirs.consens, "cat_consens.tmp")

    ## combine cons files into haplos
    with open(allcons, 'wb') as consout:
        for qhandle in conshandles:
            with gzip.open(qhandle, 'r') as tmpin:
                consout.write(tmpin.read())

    ## shuffle sequences in vsearch
    cmd = data.vsearch \
        +" --shuffle "+allcons \
        +" --output "+allcons.replace("_consens.gz", "_shuf.tmp")
    subprocess.check_call(cmd, shell=True)

    ## order by length in vsearch
    cmd = data.vsearch \
        +" --sortbylength "+allcons.replace("_consens.gz", "_shuf.tmp") \
        +" --output "+allcons.replace("_consens.gz", "_sort.tmp")
    subprocess.check_call(cmd, shell=True)

    sys.exit()
    ## message to shell
    print("Step 6: clustering across {} samples at {} similarity").\
         format(len(samples), data.paramsdict["clust_threshold"]) 

    ## write output to .consens_.gz file

    ## convert ambiguity codes into a sampled haplotype for any sample
    ## to use for clustering, but save ambiguities for later

    ## output file
    outhaplos = open(outhandle, 'w')

    ## input file
    infile = gzip.open(params["work"]+"clust"+params["wclust"]+\
                       "/cat.consens_"+gid+".gz")
    lines = iter(infile.readlines())
    infile.close()
    
    ## write to haplo files in fasta format
    writinghaplos = []

    for line in lines:
        if ">" in line:
            writinghaplos.append(line.strip())
        else:
            allele = breakalleles(line)[0]
            writinghaplos.append(allele.strip())
    outhaplos.write("\n".join(writinghaplos))
    outhaplos.close()



def main(params, inlist, gid, group, minhit, quiet):
    """ setup to call functions """

    ## make outhandle name, ... why do I need gid?
    outhandle = params["work"]+"clust"+params["wclust"]+\
                "/cat.haplos_"+gid

    ## make file with all reads 
    makecons(params, inlist, outhandle, gid, minhit, quiet) 

    if params["datatype"] == 'pairddrad':
        ## split first from second reads
        splithandle = splitter(outhandle)
        cluster(params, splithandle, gid, quiet)
    else:
        cluster(params, outhandle, gid, quiet)

    ## make clusters with .haplos, .u, and .temp files"
    makeclust(params, outhandle, gid, minhit)

    ## align clusters
    ##

    ## load into hdf database with site counts..?
    ##




def run(data, samples, ipyclient, noreverse, force):
    """ subselect and pass args for across-sample clustering """

    ## make file with all samples reads, allow priority to shunt outgroups
    ## to the end of the file
    outgroups = []
    randomseed = np.random.randint(1, int(1e9))
    makecons(data, samples, outgroups, randomseed)

    ## 



def oldpre():
    if not params["hierarch"]:
        ## Regular clustering
        if "," in params["subset"]:
            inlist = [params["work"]+"clust"+params["wclust"]+\
                      "/"+i+".consens*" for i in params["subset"]\
                      .strip().split(",")]
        else:
            inlist = glob.glob(params["work"]+"clust"+params["wclust"]+\
                               "/"+params["subset"]+"*.consens*")

        excludes = params["exclude"].strip().split(",")
        names = [i.split("/")[-1].replace(".consens.gz", "") for i in inlist]
        fulllist = [i for i, j in zip(inlist, names) if j not in excludes]

        ## cluster consens files in inlist
        cluster_cons7_shuf.main(params, fulllist, "", "", "", quiet)
        if not quiet:
            sys.stderr.write("\n\tfinished clustering\n")


    else:
        ## Hierarchical clustering
        print gids
        print groups
        print minhits

        ## make a new dir/ for hierarchs
        if not os.path.exists(params["work"]+"prefix/"):
            os.makedirs(params["work"]+"prefix/")

        ## queue up jobs
        work_queue = multiprocessing.Queue()
        result_queue = multiprocessing.Queue()

        ## submit jobs
        for (gid, minhit) in zip(gids, minhits):
            inlist = groups[gid]
            work_queue.put([params, inlist, gid, minhit, "", quiet])
    #     inlist = Hgroups[Hgid]
    #     work_queue.put([vsearch, wclust, datatype, 
    #         outgroup, seed,
    #         Hgid, Hminhit, inlist,
    #         WORK, MASK, 1 ])
                        
        ## execute first tier jobs "    
        jobs = []
        for i in range(params["parallel"]):
            worker = Worker(work_queue,
                            result_queue,
                            cluster_cons7_shuf.main)
            jobs.append(worker)
            worker.start()
            for j in jobs:
                j.join()

        ## cluster second tier
        tier2clust.main(params)
        #  Hgids, seed, WORK, MASK)
        print "\n\tfinished clustering\n"


if __name__ == "__main__":
    main()

