#!/usr/bin/env python2.7

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

from __future__ import print_function
# pylint: disable=E1101
# pylint: disable=F0401
# pylint: disable=W0142
# pylint: disable=W0212
# pylint: disable=C0301

import os
import sys
import gzip
import itertools
import subprocess
import numpy as np
import ipyrad
import time
import datetime
import shlex

from collections import Counter
from refmap import *
from util import * 

import logging
LOGGER = logging.getLogger(__name__)



def sample_cleanup(data, sample):
    """ stats, cleanup, and link to samples """
    
    ## set cluster file handles
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")

    if 'pair' in data.paramsdict["datatype"]:
        sample.files.merged = os.path.join(data.dirs.edits,
                                           sample.name+"_merged_.fastq")
        ## record how many read pairs were merged
        with open(sample.files.merged, 'r') as tmpf:
            sample.stats.reads_merged = len(tmpf.readlines()) // 4

    ## get depth stats
    infile = gzip.open(sample.files.clusters)
    duo = itertools.izip(*[iter(infile)]*2)
    depth = []
    thisdepth = 0
    maxlen = 0
    while 1:
        try:
            itera = duo.next()
            seqlen = len(itera[1])
            itera = itera[0]
            if seqlen > maxlen:
                maxlen = seqlen
        except StopIteration:
            break
        if itera != "//\n":
            try:
                tdepth = int(itera.split(";")[-2][5:])
                ## double depth for merged reads
                if "_m1;s" in itera:
                    tdepth *= 2 
                thisdepth += tdepth

            except IndexError:
                ## TODO: if no clusts pass align filter this will raise
                LOGGER.debug("Here %s, %s", sample.name, itera)
                raise IPyradError("ERROR 63: bad cluster file: %s", sample.name)
        else:
            ## append and reset
            depth.append(thisdepth)
            thisdepth = 0
    infile.close()

    ## If our longest sequence is longer than the current max_fragment_length
    ## then update max_fragment_length
    ## Padding with 4 extra bases to account for pair + sequence separator.
    if maxlen > data._hackersonly["max_fragment_length"]:
        data._hackersonly["max_fragment_length"] = maxlen + 4

    if depth:
        ## make sense of stats
        depth = np.array(depth)
        keepmj = depth[depth >= data.paramsdict["mindepth_majrule"]]    
        keepstat = depth[depth >= data.paramsdict["mindepth_statistical"]]

        ## sample summary stat assignments
        sample.stats["state"] = 3
        sample.stats["clusters_total"] = len(depth)
        sample.stats["clusters_hidepth"] = \
                                    max([len(i) for i in (keepmj, keepstat)])
        ## store large list of depths as a counter dict
        sample.depths = dict(Counter(depth))

        ## get stats for depths
        mmj = depth[depth >= data.paramsdict["mindepth_majrule"]]
        mms = depth[depth >= data.paramsdict["mindepth_statistical"]]
       
        ## sample stat assignments
        sample.stats_dfs.s3["merged_pairs"] = sample.stats.reads_merged
        sample.stats_dfs.s3["clusters_total"] = depth.shape[0]
        sample.stats_dfs.s3["clusters_hidepth"] = int(sample.stats["clusters_hidepth"])
        sample.stats_dfs.s3["avg_depth_total"] = depth.mean()
        sample.stats_dfs.s3["avg_depth_mj"] = mmj.mean()
        sample.stats_dfs.s3["avg_depth_stat"] = mms.mean()

        sample.stats_dfs.s3["sd_depth_total"] = depth.std()
        sample.stats_dfs.s3["sd_depth_mj"] = mmj.std()
        sample.stats_dfs.s3["sd_depth_stat"] = mms.std()

    else:
        print("no clusters found for {}".format(sample.name))

    ## Get some stats from the bam files
    ## This is moderately hackish. samtools flagstat returns
    ## the number of reads in the bam file as the first element
    ## of the first line, this call makes this assumption.
    if not data.paramsdict["assembly_method"] == "denovo":
        refmap_stats(data, sample)



def muscle_align(args):
    """ aligns reads, does split then aligning for paired reads """
    ## parse args
    data, chunk = args

    LOGGER.debug("aligning chunk %s", chunk)

    ## data are already chunked, read in the whole thing. 
    ## bail out if there is not a chunk file, means there were few clusts
    try:
        infile = open(chunk, 'rb')
    except IOError:
        return 0
    clusts = infile.read().split("//\n//\n")
    out = []

    ## iterate over clusters and align
    for clust in clusts:
        stack = []
        lines = clust.split("\n")
        names = lines[::2]
        seqs = lines[1::2]

        ## append counter to end of names b/c muscle doesn't retain order
        names = [j+str(i) for i, j in enumerate(names)]        

        ## don't bother aligning singletons
        if len(names) <= 1:
            if names:
                stack = [names[0]+"\n"+seqs[0]]
        else:
            ## split seqs if paired end seqs
            try:
                seqs1 = [i.split("nnnn")[0] for i in seqs] 
                seqs2 = [i.split("nnnn")[1] for i in seqs]
                ## muscle align
                string1 = muscle_call(data, names[:200], seqs1[:200])
                string2 = muscle_call(data, names[:200], seqs2[:200])
                ## resort so they're in same order as original
                anames, aseqs1 = parsemuscle(string1)
                anames, aseqs2 = parsemuscle(string2)
                ## get leftlimit of seed, no hits can go left of this 
                ## this can save pairgbs from garbage
                idxs = [i for i, j in enumerate(aseqs1[0]) if j != "-"]
                leftlimit = min(0, idxs)
                aseqs1 = [i[leftlimit:] for i in aseqs1]
                ## preserve order in ordereddict
                aseqs = zip(aseqs1, aseqs2) 

                ## append to somedic if not bad align
                for i in range(len(anames)):
                    ## filter for max internal indels
                    intindels1 = aseqs[i][0].rstrip('-').lstrip('-').count('-')
                    intindels2 = aseqs[i][1].rstrip('-').lstrip('-').count('-')
                    #if intindels1 <= data.paramsdict["max_Indels_locus"][0] & \
                    #   intindels2 <= data.paramsdict["max_Indels_locus"][1]:
                    if (intindels1 <= 5) and (intindels2 <= 5):
                        stack.append("\n".join([anames[i], 
                                        aseqs[i][0]+"nnnn"+aseqs[i][1]]))
                        #somedic[anames[i]] = aseqs[i]]))
                    else:
                        LOGGER.info("""
                high indels: %s
                1, 2: (%s, %s)
                """, aseqs[i], intindels1, intindels2)

            except IndexError:
                string1 = muscle_call(data, names[:200], seqs[:200])
                anames, aseqs = parsemuscle(string1)
                ## Get left and right limits, no hits can go outside of this. 
                ## This can save gbs overlap data significantly. 
                if 'gbs' in data.paramsdict['datatype']:
                    ## left side filter is the left limit of the seed, unless 
                    ## there is a mixture of merged and non-merged reads. If 
                    ## any reads are merged then the left limit is the left
                    ## limit of a merged read.

                    # ## if merged
                    # if any(["_m1;s" in nam for nam in anames]):
                    #     LOGGER.info("""
                    #                 THIS IS A MERGE/NON_MERGE MATCH: 
                    #                 %s""", aseqs)
                    #     idxs = []
                    #     for nam, seq in zip(anames, aseqs):
                    #         if "_m1;s" in nam:
                    #             idxs.append(\
                    #           min([i for i, j in enumerate(seq) if j != "-"]))

                    ## much simpler if no merged reads
                    ## else:
                    idxs = [i for i, j in enumerate(aseqs[0]) if j != "-"]

                    ## apply left side filter
                    if idxs:
                        leftlimit = max(0, min(idxs))
                        aseqs = [i[leftlimit:] for i in aseqs]
                        LOGGER.info('leftlimit %s', leftlimit)

                    ## right side filter is the reverse seq that goes the least
                    ## far to the right. 
                    ## Get reverse seqs and their index (index, rseq)
                    revs = [(i, aseqs[i]) for i in range(len(aseqs)) if \
                            anames[i].split(";")[-1][0] == '-']
                    ## get right side filter and remove if there's a bad match.
                    idxs = []
                    for ridx, rseq in revs:
                        try:
                            idxs.append(max(\
                                [i for i, j in enumerate(rseq) if j != "-"]))
                        except ValueError as _:
                            LOGGER.debug("\
                            Found chunk that contains a locus that's all "\
                          +"indels. Throw it out and count it as filtered.")
                            ## Remove the seq name from the names list, and 
                            ## continue with the next iteration of the for loop,
                            ## effectively drops the rseq. Use list comprehension
                            ## to drop the idx'th element and then convert back to tuple
                            anames = tuple([x for i, x in enumerate(anames) if i != ridx])
                            aseqs = tuple([x for i, x in enumerate(aseqs) if i != ridx])
                            continue
                    if idxs:
                        rightlimit = min(idxs)
                        aseqs = [i[:rightlimit] for i in aseqs]
                        #LOGGER.info('rightlimit %s', rightlimit)                    

                for i in range(len(anames)):
                    ## filter for max internal indels 
                    intind = aseqs[i].rstrip('-').lstrip('-')
                    ind1 = intind.count('-') <= \
                                data.paramsdict["max_Indels_locus"][0]
                    ## could also filter by max indel inserts                                
                    #ind2 = len([i.split("-") for i in intind if i]) < 3
                    if ind1:
                        stack.append("\n".join([anames[i], aseqs[i]]))
                        #somedic[anames[i]] = aseqs[i]
                    else:
                        LOGGER.info("high indels: %s", aseqs[i])

            ## save dict into a list ready for printing
            #for key in somedic.keys():
            #    stack.append(key+"\n"+somedic[key])

        if stack:
            out.append("\n".join(stack))

    ## write to file after
    infile.close()
    outfile = open(chunk, 'wb')#+"_tmpout_"+str(num))
    outfile.write("\n//\n//\n".join(out)+"\n")#//\n//\n")
    outfile.close()



def parsemuscle(out):
    """ parse muscle string output into two sorted lists. Sorts them first."""
    lines = out[1:].split("\n>")
    names = [line.split("\n", 1)[0] for line in lines]
    seqs = [line.split("\n", 1)[1].replace("\n", "") for line in lines]
    tups = zip(names, seqs)
    ## who knew, zip(*) is the inverse of zip
    anames, aseqs = zip(*sorted(tups, 
                        key=lambda x: int(x[0].split(";")[-1][1:])))
    ## remove reference if present
    if "_REF;+0" in anames[0]:
        return anames[1:], aseqs[1:]
    else:
        return anames, aseqs



# def test_gaps(inputstr, gapopen):
#     """ 
#     Try 3 different values for gapopen to minimize misaligns
#     """
#     ##
#     args = shlex.split("{} -quiet -clw -gapopen {} -in -"\
#                        .format(ipyrad.bins.muscle, gapopen))
#     proc1 = subprocess.Popen(args,
#                           stdin=subprocess.PIPE,
#                           stdout=subprocess.PIPE,
#                           stderr=subprocess.STDOUT)
#     proc2 = subprocess.Popen(["tr", "-d", "-C", "\*"],
#                           stdin=proc1.stdout,
#                           stdout=subprocess.PIPE)
#     proc3 = subprocess.Popen(["wc", "-c"], 
#                           stdin=proc2.stdout,
#                           stdout=subprocess.PIPE)
#     ## pass sequences to p1 and get clw alignment returned
#     clwa = proc1.communicate(inputstr)[0]
#     star = proc2.communicate()
#     nstar = int(proc3.communicate()[0].strip())

#     LOGGER.info("clwa %s", clwa)
#     LOGGER.info("nstars %s", nstar)
#     ## 
#     return clwa, nstar




def muscle_call(data, names, seqs):
    """ Makes subprocess call to muscle. A little faster than before. 
    TODO: Need to make sure this works on super large strings and does not 
    overload the PIPE buffer.
    """

    ## get input string
    inputstr = "\n".join([">{}\n{}".format(i, j) for i, j in zip(names, seqs)])
    proc1 = subprocess.Popen([ipyrad.bins.muscle, 
        "-quiet", "-gapopen", "-1200"], 
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE)

    return proc1.communicate(inputstr)[0]

    # ## default test
    # clwa, nstar = test_gaps(inputstr, -600)

    # ## try two more stringent values
    # for gaps in [-1200, -2400]:
    #     tclwa, tnstar = test_gaps(inputstr, gaps)
    #     if tnstar < nstar:
    #         clwa = tclwa
    #     LOGGER.info("CLWA %s", clwa)

    # return clwa

    #inputstring = "\n".join(">"+i+"\n"+j for i, j in zip(names, seqs))
    #return subprocess.Popen(ipyrad.bins.muscle, 
    #                        stdin=subprocess.PIPE, 
    #                        stdout=subprocess.PIPE)\
    #                        .communicate(inputstring)[0]



def build_clusters(data, sample):
    """ combines information from .utemp and .htemp files 
    to create .clust files, which contain un-aligned clusters """

    ## derepfile 
    derepfile = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")

    ## vsearch results files
    ufile = os.path.join(data.dirs.clusts, sample.name+".utemp")
    htempfile = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## create an output file to write clusters to        
    clustfile = gzip.open(os.path.join(data.dirs.clusts,
                          sample.name+".clust.gz"), 'wb')
    sample.files["clusts"] = clustfile

    ## if .u files are present read them in as userout
    try:
        userout = open(ufile, 'rb').readlines()
    except IOError:
        inst = """
    No clusters (.utemp hits) found for {}. If you are running preview mode and
    the size of the truncated input file isn't big enough try increasing the 
    size of <your_assembly>._hackersonly[\"preview_truncate_length\"
    """.format(sample.name)
        LOGGER.warn(inst)
        raise IPyradError(inst)

    ## load derep reads into a dictionary
    hits = {}  
    ioderep = open(derepfile, 'rb')
    dereps = itertools.izip(*[iter(ioderep)]*2)
    for namestr, seq in dereps:
        _, count = namestr.strip()[:-1].split(";size=")
        hits[namestr] = [int(count), seq.strip()]
    ioderep.close()

    ## create dictionary of .utemp cluster hits
    udic = {} 
    for uline in userout:
        uitems = uline.strip().split("\t")
        ## if the seed is in udic
        if ">"+uitems[1]+"\n" in udic:
            ## append hit to udict
            udic[">"+uitems[1]+'\n'].append([">"+uitems[0]+"\n", 
                                            uitems[4],
                                            uitems[5].strip(),
                                            uitems[3]])
        else:
            ## write as seed
            udic[">"+uitems[1]+'\n'] = [[">"+uitems[0]+"\n", 
                                         uitems[4],
                                         uitems[5].strip(),
                                         uitems[3]]]

    ## map sequences to clust file in order
    seqslist = [] 
    for key, values in udic.iteritems():
        ## this is the seed. Store the left most non indel base for the seed. 
        ## Do not allow any other hits to go left of this (for pairddrad)
        seedhit = hits[key][1]
        seq = [key.strip()+"*\n"+seedhit]

        ## allow only N internal indels in hits to seed for within-sample clust
        ## prior to alignment. This improves alignments. 
        ## could be written a little more cleanly but this allows better debug.
        for i in xrange(len(values)):
            inserts = int(values[i][3])
            if values[i][1] == "+":
                fwdseq = hits[values[i][0]][1]
                if inserts < 6:
                    seq.append(values[i][0].strip()+"+\n"+fwdseq)
                else:
                    fwdseq = hits[values[i][0]][1]                    
                    LOGGER.debug("exc indbld: %s %s", inserts, fwdseq)
            ## flip to the right orientation 
            else:
                revseq = comp(hits[values[i][0]][1][::-1])
                if inserts < 6:
                    seq.append(values[i][0].strip()+"-\n"+revseq)
                else:
                    LOGGER.debug("exc indbld: %s %s", inserts, revseq)

        seqslist.append("\n".join(seq))
    clustfile.write("\n//\n//\n".join(seqslist)+"\n")

    ## make Dict. from seeds (_temp files) 
    iotemp = open(htempfile, 'rb')
    invars = itertools.izip(*[iter(iotemp)]*2)
    seedsdic = {k:v for (k, v) in invars}  
    iotemp.close()

    ## create a set for keys in I not in seedsdic
    set1 = set(seedsdic.keys())   ## htemp file (no hits) seeds
    set2 = set(udic.keys())       ## utemp file (with hits) seeds
    diff = set1.difference(set2)  ## seeds in 'temp not matched to in 'u
    if diff:
        for i in list(diff):
            clustfile.write("//\n//\n"+i.strip()+"*\n"+hits[i][1]+'\n')
    #clustfile.write("//\n//\n\n")
    clustfile.close()
    del dereps
    del userout
    del udic
    del seedsdic



def setup_dirs(data):
    """ sets up directories for step3 data """
    ## make output folder for clusters  
    data.dirs.clusts = os.path.join(
                os.path.realpath(data.paramsdict["project_dir"]),
                data.name+"_"+"clust_"+str(data.paramsdict["clust_threshold"]))
    if not os.path.exists(data.dirs.clusts):
        os.makedirs(data.dirs.clusts)

    ## If ref mapping, init samples and make the refmapping output directory.
    if not data.paramsdict["assembly_method"] == "denovo":
        ## make output directory for read mapping process
        data.dirs.refmapping = os.path.join(
                        os.path.realpath(data.paramsdict["project_dir"]),
                        data.name+"_refmapping")
        if not os.path.exists(data.dirs.refmapping):
            os.makedirs(data.dirs.refmapping)

    ## clust directory
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)



def null_func(args):
    """ 
    Takes a list of args and prints them to the logger. Used as a null func
    to skip steps that are only used in some assembly steps. 
    """
    LOGGER.info("null func skipping")#  %s", str(arglist))



def apply_jobs(data, samples, ipyclient, noreverse, force, preview):
    """ pass the samples to N engines to execute run_full on each.

    :param data: An Assembly object
    :param samples: one or more samples selected from data
    :param ipyclient: ipyparallel load_balanced_view client
    :param noreverse: toggle revcomp clustering despite datatype default
    :param force: force
    :param preview: run preview
    :param align_only: skips clustering/mapping and aligns existing files

    :returns: None
    """
    ## make directories
    setup_dirs(data)

    ## Create threaded_view of engines by grouping only ids that are threaded
    hostdict = get_threaded_view(ipyclient)
    threaded_views = {}
    for key, val in hostdict.items():
        ## e.g., client.load_balanced_view([1,3])
        threaded_views[key] = ipyclient.load_balanced_view(val)

    ## A single load-balanced view for FUNCs 3-4
    lbview = ipyclient.load_balanced_view()

    ## FUNC 1: derep ---------------------------------------------------
    res_derep = {}
    for sample in samples:
        args = [data, sample]
        res_derep[sample] = lbview.apply(derep_concat_split, args)

    ## FUNC 2: init ref -----------------------------------------------------
    ## pull out the seqs that match to the reference
    done_ref = 1
    res_refinit = {}
    mcfunc = null_func
    if "reference" in data.paramsdict["assembly_method"]:
        done_ref = 0
        ## if reference set func to mapping func
        mcfunc = refmap_init
    ## Initialize the mapped and unmapped file paths per sample
    for sample in samples:
        args = [data, sample]
        with lbview.temp_flags(after=res_derep[sample]):            
            res_refinit[sample] = lbview.apply(mcfunc, args)            
    #[res_refinit[i].get() for i in res_refinit]
    #return 1


    ## FUNC 3: mapreads ----------------------------------------------------
    ## Submit samples to reference map else null
    res_ref = {}
    mcfunc = null_func
    if "reference" in data.paramsdict["assembly_method"]:
        mcfunc = mapreads
    for sample in samples:
        args = [data, sample, noreverse, 1]
        with lbview.temp_flags(after=res_refinit[sample]):
            res_ref[sample] = lbview.apply(mcfunc, args)
    ## test just up to this point [comment this out when done]
    #[res_ref[i].get() for i in res_ref]
    #return 1


    ## FUNC 4: clustering ---------------------------------------------------
    mcfunc = null_func
    done_clust = 1
    if data.paramsdict["assembly_method"] in ["denovo", "denovo+reference"]: 
        done_clust = 0
        mcfunc = clust_and_build
    ## require that the sample successfully finished previous step
    res_clust = {}
    for sample in samples:
        args = [data, sample, noreverse, 1] 
        with lbview.temp_flags(after=res_ref[sample]):
            res_clust[sample] = lbview.apply(mcfunc, args)
    ## test just up to this point [comment this out when done]
    #[res_clust[i].get() for i in res_clust]
    #return 1


    ## FUNC 5: reference cleanup -------------------------------------------
    ## Pull in alignments from mapped bam files and write them to the clust.gz 
    ## to fold them back into the pipeline. If we are doing "denovo" then 
    ## don't call this, but less obvious, "denovo-reference" intentionally 
    ## doesn't call this to effectively discard reference mapped reads.
    mcfunc = null_func
    if data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]: 
        mcfunc = ref_muscle_chunker
    ## requires sample to have finished the previous step before running
    res_clean = {}
    for sample in samples:
        with lbview.temp_flags(after=res_clust[sample]):
            res_clean[sample] = lbview.apply(mcfunc, [data, sample])
    ## test just up to this point [comment this out when done]
    #[res_clean[i].get() for i in res_clean]
    #return 1


    ## FUNC 6: split up clusters into chunks -------------------------------
    #mcfunc = null_func
    #if data.paramsdict["assembly_method"] in ["denovo", "denovo+reference"]: 
    mcfunc = muscle_chunker
    res_chunk = {}
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')    
    for sample in samples:
        with lbview.temp_flags(after=res_clean[sample]):
            args = [data, sample, tmpdir]
            res_chunk[sample] = lbview.apply(mcfunc, args)
    ## test just up to this point [comment this out when done]
    #[res_chunk[i].get() for i in res_chunk]
    #return 1


    ## FUNC 7: align chunks -------------------------------------------------
    res_align = {sample:[] for sample in samples}
    for sample in samples:
        ## get all chunks for this sample
        with lbview.temp_flags(after=res_chunk[sample]):
            for i in range(10):
                chunk = os.path.join(tmpdir, sample.name+"_chunk_{}.ali".format(i))
                res_align[sample].append(lbview.apply(muscle_align, [data, chunk]))
    ## get all async results for aligners
    all_aligns = list(itertools.chain(*res_align.values()))    
    ## test just up to this point [comment this out when done]
    #[i.get() for i in all_aligns]
    #return 1


    ## FUNC 6: concat chunks -------------------------------------------------
    res_concat = {}
    for sample in samples:
        #LOGGER.info('resalign[sample] %s', res_align[sample])
        tmpids = list(itertools.chain(*[i.msg_ids for i in res_align[sample]]))
        #LOGGER.info('tmpids %s', tmpids)
        with lbview.temp_flags(after=tmpids):
            res_concat[sample] = lbview.apply(reconcat, [data, sample])
    ## test just up to this point [comment this out when done]
    #[res_concat[i].get() for i in res_concat]

    ## wait func
    tmpids = list(itertools.chain(*[i.msg_ids for i in res_concat.values()]))
    with lbview.temp_flags(after=tmpids):
        res = lbview.apply(time.sleep, 0.1)    

    ## print progress bars
    clusttotal = len(res_clust)
    all_refmap = len(res_ref)
    all_aligns = list(itertools.chain(*res_align.values()))
    aligntotal = len(all_aligns)

    while 1:
        if not res.ready():
            if not done_ref:
                ## prints a progress bar
                fref = sum([res_ref[i].ready() for i in res_ref])
                elapsed = datetime.timedelta(seconds=int(res.elapsed))                
                progressbar(all_refmap, fref, 
                            " mapping reads     | {}".format(elapsed))
                ## go to next print row when done
                if fref == all_refmap:
                    done_ref = 1
                    print("")

            elif not done_clust:
                ## prints a progress bar
                fclust = sum([res_clust[i].ready() for i in res_clust])
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(clusttotal, fclust, 
                            " clustering reads  | {}".format(elapsed))
                ## go to next print row when done
                if fclust == clusttotal:
                    done_clust = 1
                    print("")

            else:
                falign = sum([i.ready() for i in all_aligns])
                elapsed = datetime.timedelta(seconds=int(res.elapsed))
                progressbar(aligntotal, falign, 
                    " aligning clusters | {}".format(elapsed))
            sys.stdout.flush()
            time.sleep(1)

        else:
            elapsed = datetime.timedelta(seconds=int(res.elapsed))                            
            progressbar(20, 20,
                " aligning clusters | {}".format(elapsed))
            print("")
            break

    ## Cleanup -------------------------------------------------------
    for sample in samples:
        sample_cleanup(data, sample)

    data_cleanup(data)



def jobs_cleanup(data, samples, preview):
    """
    For preview/refmap restore original sample.files.edits paths and 
    clean up the tmp files.

    If we did refmapping return the samples.files.edits to their original
    condition. Have to do this before restoring preview files because
    the refmap edits backup will be a backup of the preview truncated 
    files. The order of these two functions matters.
    """

    if not "denovo" == data.paramsdict["assembly_method"]: 
        for sample in samples:
            refmap_cleanup(data, sample)

    if preview:
        for sample in samples:
            try:
                ## If edits and edits_preview_bak are the same then 
                ## something borked so don't delete any files
                if sample.files.edits == sample.files.edits_preview_bak:
                    sample.files.pop("edits_preview_bak", None)
                    continue

                for tmpf in sample.files.edits[0]:
                    if os.path.exists(tmpf):
                        os.remove(tmpf)

                ## Restore original paths to full fastq files
                sample.files.edits = sample.files.edits_preview_bak
                ## Remove the tmp file reference. The second arg defines 
                ## what to return if the key doesn't exist.
                sample.files.pop("edits_preview_bak", None)
            except Exception:
                pass



def derep_and_sort(data, sample):
    """ 
    Dereplicates reads and sorts so reads that were highly replicated are at
    the top, and singletons at bottom, writes output to derep file. Paired
    reads are dereplicated as one concatenated read and later split again. 
    """

    LOGGER.info("in the real derep; %s", sample.name)
    ## reverse complement clustering for some types    
    if "gbs" in data.paramsdict["datatype"]:
        reverse = " -strand both "
    else:
        reverse = " "

    ## do dereplication with vsearch
    cmd = ipyrad.bins.vsearch\
         +" -derep_fulllength "+sample.files.edits[0][0]\
         +reverse \
         +" -output "+os.path.join(data.dirs.edits, sample.name+"_derep.fastq")\
         +" -sizeout " \
         +" -threads 1 "\
         +" -fasta_width 0"
    LOGGER.info(cmd)

    ## run vsearch
    try:
        LOGGER.info(cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        LOGGER.error(inst)
        raise IPyradError(inst)



def data_cleanup(data):
    """ cleanup / statswriting function for Assembly obj """
    data.stats_dfs.s3 = data.build_stat("s3")
    data.stats_files.s3 = os.path.join(data.dirs.clusts, "s3_cluster_stats.txt")
    with open(data.stats_files.s3, 'w') as outfile:
        data.stats_dfs.s3.to_string(outfile)



def concat_edits(data, sample):
    """ 
    Concatenate if multiple edits files for a sample
    """
    LOGGER.debug("Entering concat_edits: %s", sample.name)
    ## if more than one tuple in the list
    if len(sample.files.edits) > 1:
        ## create temporary concat file
        tmphandle1 = os.path.join(data.dirs.edits,
                              "tmp1_"+sample.name+".concat")       
        with open(tmphandle1, 'wb') as tmp:
            for editstuple in sample.files.edits:
                with open(editstuple[0]) as inedit:
                    tmp.write(inedit)

        if 'pair' in data.paramsdict['datatype']:
            tmphandle2 = os.path.join(data.dirs.edits,
                                "tmp2_"+sample.name+".concat")       
            with open(tmphandle2, 'wb') as tmp:
                for editstuple in sample.files.edits:
                    with open(editstuple[1]) as inedit:
                        tmp.write(inedit)
            sample.files.edits = [(tmphandle1, tmphandle2)]
        else:
            sample.files.edits = [(tmphandle1, )]
    return sample



def cluster(data, sample, noreverse, nthreads):
    """ calls vsearch for clustering. cov varies by data type, 
    values were chosen based on experience, but could be edited by users """
    ## get files
    derephandle = os.path.join(data.dirs.edits, sample.name+"_derep.fastq")
    uhandle = os.path.join(data.dirs.clusts, sample.name+".utemp")
    temphandle = os.path.join(data.dirs.clusts, sample.name+".htemp")

    ## datatype specific optimization
    ## minsl: the percentage of the seed that must be matched
    ##    smaller values for RAD/ddRAD where we might want to combine, say 50bp 
    ##    reads and 100bp reads in the same analysis. 
    ## query_cov: the percentage of the query sequence that must match seed
    ##    smaller values are needed for gbs where only the tips might overlap
    ##    larger values for pairgbs where they should overlap near completely
    ##    small minsl and high query cov allows trimmed reads to match to untrim
    ##    seed for rad/ddrad/pairddrad. 
    if data.paramsdict["datatype"] == "gbs":
        reverse = " -strand both "
        cov = " -query_cov .33 " 
        minsl = " 0.33"
    elif data.paramsdict["datatype"] == 'pairgbs':
        reverse = "  -strand both "
        cov = " -query_cov .75 " 
        minsl = " 0.75"  
    else:  ## rad, ddrad
        reverse = " -leftjust "
        cov = " -query_cov .90 "
        minsl = " 0.5"

    ## override query cov
    if data._hackersonly["query_cov"]:
        cov = " -query_cov "+str(data._hackersonly["query_cov"])
        assert data._hackersonly["query_cov"] <= 1, "query_cov must be <= 1.0"

    ## override reverse clustering option
    if noreverse:
        reverse = " -leftjust "
        LOGGER.warn(noreverse, "not performing reverse complement clustering")

    ## get call string
    cmd = ipyrad.bins.vsearch+\
        " -cluster_smallmem "+derephandle+\
        reverse+\
        cov+\
        " -id "+str(data.paramsdict["clust_threshold"])+\
        " -userout "+uhandle+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -minsl "+str(minsl)+\
        " -fulldp"+\
        " -threads "+str(nthreads)+\
        " -usersort "+\
        " -notmatched "+temphandle+\
        " -fasta_width 0"

    ## run vsearch
    try:
        LOGGER.debug("%s", cmd)
        subprocess.call(cmd, shell=True,
                             stderr=subprocess.STDOUT,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as inst:
        sys.exit("Error in vsearch: \n{}\n{}".format(inst, subprocess.STDOUT))



def muscle_chunker(args):
    """ 
    Splits the muscle alignment into chunks. Each runs on N clusters at a time.
    """
    ## parse args
    LOGGER.info("inside chunker")
    data, sample, tmpdir = args

    ## get the number of clusters
    clustfile = os.path.join(data.dirs.clusts, sample.name+".clust.gz")
    with gzip.open(clustfile, 'rb') as clustio:
        tclust = clustio.read().count("//")//2
        optim = (tclust//10) + (tclust%10)
        LOGGER.info("optim for align chunks: %s", optim)

    ## write optim clusters to each tmp file
    clustio = gzip.open(clustfile, 'rb')
    inclusts = iter(clustio.read().strip().split("//\n//\n"))
    grabchunk = list(itertools.islice(inclusts, optim))

    idx = 0
    while grabchunk:
        tmpfile = os.path.join(tmpdir, sample.name+"_chunk_{}.ali".format(idx)) 
        with open(tmpfile, 'wb') as out:
            out.write("//\n//\n".join(grabchunk))
        idx += 1
        grabchunk = list(itertools.islice(inclusts, optim))
    clustio.close()



def reconcat(args):
    """ aligns chunked seqs using muscle """

    ## parse args
    data, sample = args

    ## get chunks
    tmpdir = os.path.join(data.dirs.project, data.name+'-tmpalign')
    chunks = glob.glob(os.path.join(tmpdir, sample.name+"_chunk_*"))
   
    ## sort by chunk number
    chunks.sort(key=lambda x: int(x.rsplit("_", 1)[-1][:-4]))

    LOGGER.info("chunk %s", chunks)
    ## concatenate finished reads
    sample.files.clusters = os.path.join(data.dirs.clusts,
                                         sample.name+".clustS.gz")
    ## reconcats aligned clusters
    with gzip.open(sample.files.clusters, 'wb') as out:
        for fname in chunks:
            with open(fname) as infile:
                out.write(infile.read()+"//\n//\n")
            os.remove(fname)



def alignment_cleanup(tmpnames):
    ## still delete tmpfiles if job was interrupted
    for fname in tmpnames:
        if os.path.exists(fname):
            os.remove(fname)

    if os.path.exists(tmpdir):
        try:
            shutil.rmtree(tmpdir)
        except OSError as inst:
            ## In some instances nfs creates hidden dot files in directories
            ## that claim to be "busy" when you try to remove them. Don't
            ## kill the run if you can't remove this directory.
            LOGGER.warn("Failed to remove tmpdir {}".format(tmpdir))



def derep_concat_split(args):
    """ 
    Running on remote Engine. Refmaps, then merges, then dereplicates, 
    then denovo clusters reads. 
    """

    ## get args
    data, sample = args
    LOGGER.info("INSIDE derep %s", sample.name)

    ## concatenate edits files in case a Sample has multiple, and 
    ## returns a new Sample.files.edits with the concat file
    sample = concat_edits(data, sample)

    ## Denovo: merge or concat fastq pairs [sample.files.pairs]
    ## Reference: only concat fastq pairs  []
    ## Denovo + Reference: 
    if 'pair' in data.paramsdict['datatype']:
        ## merge pairs that overlap and combine non-overlapping
        ## pairs into one merged file. merge_pairs takes the unmerged
        ## files list as an argument because we're reusing this code 
        ## in the refmap pipeline, trying to generalize.
        LOGGER.debug("Merging pairs - %s", sample.files)
        merge = 1
        if data.paramsdict["assembly_method"] in ["reference", "denovo+reference"]:
            merge = 0
        sample = merge_pairs(data, sample, merge)
        sample.files.edits = [(sample.files.merged, )]
        LOGGER.debug("Merged file - {}".format(sample.files.merged))

    ## convert fastq to fasta, then derep and sort reads by their size
    derep_and_sort(data, sample)
    


def clust_and_build(args):
    """ 
    cluster and build clusters
    """
    ## parse args
    data, sample, noreverse, nthreads = args

    ## cluster derep fasta files in vsearch 
    cluster(data, sample, noreverse, nthreads)

    ## cluster_rebuild. Stop and print warning if no .utemp hits
    try:
        build_clusters(data, sample)
        ## record that it passed the clustfile build
        return 1
    except IPyradError as inst:
        print(inst)
        return 0



def run(data, samples, noreverse, force, preview, ipyclient):
    """ run the major functions for clustering within samples """

    ## list of samples to submit to queue
    subsamples = []

    ## if sample is already done skip
    for sample in samples:
        ## If sample not in state 2 don't try to cluster it.
        if sample.stats.state < 2:
            print("""\
    Sample not ready for clustering. First run step2 on sample: {}""".\
    format(sample.name))
            continue

        if not force:
            if sample.stats.state >= 3:
                print("""\
    Skipping {}; aleady clustered. Use force to re-cluster""".\
    format(sample.name))
            else:
                if sample.stats.reads_filtered:
                    subsamples.append(sample)
        else:
            ## force to overwrite
            if sample.stats.reads_filtered:            
                subsamples.append(sample)

    ## run subsamples 
    if not subsamples:
        print("  No Samples ready to be clustered. First run step2().")

    else:
        ## arguments to apply_jobs, inst catches exceptions
        args = [data, subsamples, ipyclient, noreverse, force, preview]
        inst = ""

        ## wrap job in try/finally to ensure cleanup
        try:
            apply_jobs(*args)
        except Exception as inst:
            LOGGER.warn(inst)
            raise
        finally:
            pass
            #jobs_cleanup(data, samples, preview)
            #if inst:
            #    raise inst


if __name__ == "__main__":
    ## test...

    ## reload autosaved data. In case you quit and came back 
    DATA = ipyrad.load_json("cli/cli.json")

    ## run step 6
    DATA.step3(force=True)

    # DATA = Assembly("test")
    # DATA.get_params()
    # DATA.set_params(1, "./")
    # DATA.set_params(28, '/Volumes/WorkDrive/ipyrad/refhacking/MusChr1.fa')
    # DATA.get_params()
    # print(DATA.log)
    # DATA.step3()
    #PARAMS = {}
    #FASTQS = []
    #QUIET = 0
    #run(PARAMS, FASTQS, QUIET)

    
