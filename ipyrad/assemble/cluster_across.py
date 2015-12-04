#!/usr/bin/env python2

""" cluster across samples using vsearch with options for 
    hierarchical clustering """

import os
import sys
import itertools
import random
import subprocess
import gzip
import copy
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



def makecons(params, inlist, outhandle, gid, minhit, quiet):
    """ Make the concatenated consens files to imput to vsearch. 
    Orders reads by length and shuffles randomly within length classes"""

    ## make list of consens files but cats already made
    consfiles = [i for i in inlist if "/cat.cons" not in i]
    consfiles = [i for i in consfiles if "/cat.group" not in i]

    if not consfiles:
        sys.exit("no consens files found")

    ##  make a copy list that will not have outgroups excluded
    conswithouts = copy.copy(inlist)
    
    ##are files gzipped ?
    gzp = ""
    if any([i.endswith(".gz") for i in consfiles]):
        gzp = ".gz"

    ## remove previous files if present "
    check = params["work"]+'clust'+params["wclust"]+'/cat.consens_'+gid+gzp
    if os.path.exists(check):
        os.remove(check)
    check = params["work"]+'clust'+params["wclust"]+'/cat.group_'+gid+gzp
    if os.path.exists(check):
        os.remove(check)

    ## remove outgroup sequences, add back in later to bottom after shuffling "
    outgroup = ""
    if params["outgroup"]:
        outgroup = params["outgroup"].strip().split(",")
        if len(outgroup) > 1:
            for samp in outgroup:
                check = params["work"]+"clust"+params["wclust"]+\
                        "/"+samp+".consens"+gzp
                if check in consfiles:
                    consfiles.remove(check)
        else:
            outgroup = params["work"]+"clust"+params["wclust"]+"/"+\
                       params["outgroup"]+".consens"+gzp
            if outgroup in consfiles:
                consfiles.remove(outgroup)
                
    ## output file for consens seqs from all taxa in consfiles list
    out = gzip.open(params["work"]+'clust'+params["wclust"]+\
                    '/cat.group_'+gid+gzp, 'w')

    ## iterate over files and...
    for qhandle in consfiles:
        if gzp:
            consfile = gzip.open(qhandle)
        else:
            consfile = open(qhandle)
        duo = itertools.izip(*[iter(consfile)]*2)
        while 1:
            try: 
                itera = duo.next()
            except StopIteration: 
                break
            print >>out, itera[0].strip()+"    "+itera[1].strip()
        consfile.close()
    out.close()

    ## message to shell
    if not quiet:
        if gid:
            sys.stderr.write('\n\tstep 6: clustering across '+\
                             str(len(consfiles))+\
                             " samples at "+params["wclust"]+" similarity "+
                             "\n\tfor group ("+str(gid)+") retaining seeds "+\
                             "w/ minimum of "+str(minhit)+\
                             " hits\n\n")
        else:
            sys.stderr.write('\n\tstep 6: clustering across '+\
                             str(len(consfiles))+\
                             ' samples at '+params["wclust"]+\
                             ' similarity \n\n')

    ## make list of random number and data
    if params["seed"]:
        random.seed(params["seed"])

    ## open file for reading consensus reads grouped together in one file
    source = gzip.open(params["work"]+'clust'+params["wclust"]+\
                       '/cat.group_'+gid+".gz", 'r')

    ## generator to add a random number next to every sequence
    data = ((random.random(), line) for line in source)

    ## sort by the random number into a list (now stored in memory)
    randomized_data = sorted(data)
    source.close()
    
    ## order by size while retaining randomization within size classes
    splitlines = (line.split('    ') for _, line in randomized_data)
    equalspacers = iter("".join([i[0]+" "*(100-len(i[0])), i[1]]) \
                        for i in splitlines)
    orderedseqs = sorted(equalspacers, key=len, reverse=True)
    nitera = iter(["**".join([i.split(" ")[0], i.split(" ")[-1]]) \
                         for i in orderedseqs])

    ## write output to .consens_.gz file
    ## NB: could probably speed this up
    out = gzip.open(params["work"]+'clust'+params["wclust"]+\
                    '/cat.consens_'+gid+".gz", 'wb')
    while 1:
        try:
            name, seq = nitera.next().split("**")
        except StopIteration: 
            break
        print >>out, name+'\n'+seq.strip()
    
    ##  add outgroup taxa back onto end of file.
    ## append to existing consens_file
    if outgroup:
        if len(outgroup) > 1:
            for samp in outgroup:
                xoutg = params["work"]+"clust"+\
                        params["wclust"]+"/"+samp+".consens.gz"
                if xoutg in conswithouts:
                    oreads = gzip.open(xoutg)
                    duo = itertools.izip(*[iter(oreads)]*2)
                    while 1:
                        try: 
                            oitera = duo.next()
                        except StopIteration:
                            break
                        print >>out, oitera[0].strip()+"\n"+\
                                     oitera[1].strip()
                    oreads.close()
        elif len(outgroup) == 1:
            xoutg = params["work"]+"clust"+params["wclust"]+\
                    "/"+outgroup[0]+".consens.gz"
            if xoutg in conswithouts:
                oreads = gzip.open(xoutg)
                duo = itertools.izip(*[iter(oreads)]*2)
                while 1:
                    try:
                        oitera = duo.next()
                    except StopIteration:
                        break
                    print >>out, oitera[0].strip()+\
                                 "\n"+oitera[1].strip()
                oreads.close()
        else:
            pass
    out.close()        

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

    ## make 
    makecons(params, inlist, outhandle, gid, minhit, quiet) 

    if params["datatype"] == 'pairddrad':
        ## split first from second reads
        splithandle = splitter(outhandle)
        cluster(params, splithandle, gid, quiet)
    else:
        cluster(params, outhandle, gid, quiet)

    ## make clusters with .haplos, .u, and .temp files"
    makeclust(params, outhandle, gid, minhit)



def run(data, samples, ipyclient, preview, noreverse, force):
    """ subselect and pass args for across-sample clustering """

    ## list of samples to submit to queue

    if all([isinstance(i, str) for i in samples]):
        subsamples = []

    ## 
    outhandle = os.path.join(data.dirs.consens, "cat.haplos")
    ##
    makecons(data, subsamples, outhandle)



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

