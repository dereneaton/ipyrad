#!/usr/bin/env python2.7

""" call consensus base calls on single-end data """

from __future__ import print_function
# pylint: disable=E1101

import ipyparallel as ipp
import scipy.stats
import scipy.misc
import itertools
import numpy
import time
import gzip
import os

from .demultiplex import zcat_make_temps
from ipyrad.assemble.worker import ObjDict
from collections import Counter

import logging
LOGGER = logging.getLogger(__name__)


def binomprobr(base1, base2, error, het):
    """
    given two bases are observed at a site
    n1 and n2, and the error rate e, the
    probability the site is truly aa,bb,ab
    is calculated using binomial distribution
    as in Li_et al 2009, 2011, and if
    coverage > 500, 500 dereplicated reads were
    randomly sampled.
    """
    ## major allele freq
    mjaf = base1/float(base1+base2)
    prior_homo = ((1.-het)/2.)
    prior_het = het

    ## get probabilities
    hetro = scipy.misc.comb(base1+base2, base1)/(2.**(base1+base2))
    homoa = scipy.stats.binom.pmf(base1, base1+base2, error)
    homob = scipy.stats.binom.pmf(base2, base1+base2, error)

    ## calculate probs
    homoa = homoa*prior_homo
    homob = homob*prior_homo
    hetro = hetero*prior_het

    ## return 
    probabilities = [homoa, homob, hetro]
    genotypes = ['aa', 'bb', 'ab']
    bestprob = max(probabilities)/float(sum(probabilities))

    return [bestprob, mjaf, genotypes[probabilities.index(max(probabilities))]]



def simpleconsensus(base1, base2):
    """
    majority consensus calling for sites
    with too low of coverage for
    statistical calling. Only used
    with 'lowcounts' option. Returns 
    the most common base. Returns consistent
    alphabetical order for ties.
    """
    #qQn = ['aa','bb','ab']
    maf = base1/(base1+base2)
    return [1.0, maf, 'aa']



def hetero(base1, base2):
    """
    returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites.
    """
    iupac = "N"
    trans = {('G', 'A'):"R",
             ('G', 'T'):"K",
             ('G', 'C'):"S",
             ('T', 'C'):"Y",
             ('T', 'A'):"W",
             ('C', 'A'):"M"}
    order1 = trans.get((base1, base2))
    order2 = trans.get((base2, base1))
    if order1:
        iupac = order1
    elif order2:
        iupac = order2
    else:
        ## one or both are Ns
        if [base1, base2].count("N") == 2:
            pass
        elif base1 == "N":
            iupac = base2
        elif base2 == "N":
            iupac = base1
        ## one or both are (-)s            
        else:
            if base1 == "-":
                iupac = base2
            elif base2 == "-":
                iupac = base1
            else:
                LOGGER.error("unrecognized sequence: %s %s", base1, base2)
    return iupac



def unhetero(amb):
    " returns bases from ambiguity code"
    amb = amb.upper()
    trans = {"R":("G", "A"),
             "K":("G", "T"),
             "S":("G", "C"),
             "Y":("T", "C"),
             "W":("T", "A"),
             "M":("C", "A")}
    return trans.get(amb)



def uplow(hsite):
    """ allele precedence used in assigning upper and lower
    case letters ot a consensus sequence to store the the phased
    allele pattern for diploids. 
    G > T > C > A """
    prec = {('G', 'A'):"G",
            ('A', 'G'):"G",
            ('G', 'T'):"G",
            ('T', 'G'):"G",
            ('G', 'C'):"G",
            ('C', 'G'):"G",
            ('T', 'C'):"T",
            ('C', 'T'):"T",
            ('T', 'A'):"T",
            ('A', 'T'):"T",
            ('C', 'A'):"C",
            ('A', 'C'):"C"}
    bigbase = prec.get(hsite)
    if not bigbase:
        bigbase = hsite[0]
    return bigbase
    


def findalleles(consens, heteros, counted_h, sloc):
    """ store phased allele data for diploids """
    ## find the first hetero site from the left
    alleles = [i[0] for i in counted_h.most_common(2)]
    bigbase = uplow(unhetero(consens[heteros[0]]))
    try: 
        query = [i for i in alleles if i[0] == bigbase][0]
    except IndexError:
        pass
        # #print "THIS ONE HERE"
        # print alleles, 'alleles'
        # print consens, 'consens'
        # print heteros, 'heteros'
        # print counted_h, 'counted_h'
        # print bigbase, 'bigbase'
        # for i in sloc:
        #     print "".join(i)

    for ind, hsite in enumerate(heteros[1:]):
        bigbase = uplow(unhetero(consens[hsite]))
        if query[ind+1] != bigbase:
            consens[hsite] = consens[hsite].lower()
    return consens



def removerepeat_Ns(shortcon, stacked):
    """ checks for interior Ns in consensus seqs
        remove those that arise next to *single repeats*
        of at least 3 bases on either side, which may be
        sequencing errors on deep coverage repeats """
    ## TODO: continue trying to improve removal of columns with 
    ## indel singletons, likely errors.

    ### what is the most common site depth
    #comdepth = Counter([sum(i.values()) for i in stacked]).most_common(1)[0][1]
    #print comdepth, 'comdepth'

    ### below this depth it is likely a sequencing repeat
    #mindepth = comdepth/10.
    #print mindepth, 'mindepth'

    nnsites = [i for i, j in enumerate(shortcon) if j == "N"]
    repeats = set()
    for nsite in nnsites:
        ## grab the three bases to the left of this site
        isvar = len(set(list(shortcon)[nsite-3:nsite]))
        if isvar < 2:
            ## could use mindepth here to increase stringency
            repeats.add(nsite)
        ## grab the three bases to the rigth of this n site
        isvar2 = len(set(list(shortcon)[nsite+1:nsite+4]))
        if isvar2 < 2:
            ## could use mindepth here to increase stringency            
            repeats.add(nsite)

    ## remove repeat sites from shortcon and stacked
    shortcon = "".join([j for (i, j) in enumerate(shortcon) if \
                                      i not in repeats])
    stacked = [j for (i, j) in enumerate(stacked) if \
                                      i not in repeats]
    return shortcon, stacked



def gbs_edgefilter(name, seq, leftjust, rights):
    """ filter edges in case of partial overlap in revcomp seqs. 
    Does not allow reverse reads to go left of the seed to prevent 
    resequencing of the barcode."""

    ## leftjust is seed's left most non missing base"
    if name.split(";")[-1][0] == "*":
        leftjust = seq.index([i for i in seq if i not in list("-N")][0])

    ## rightjust is the shortest reverse hit "
    if name.split(";")[-1][0] == "-":
        rights.append(max(-1, [seq.rindex(i) for i in seq \
                               if i in list("ATGC")]))
    return leftjust, rights



def consensus(args):
    """
    from a clust file handle, reads in all copies at a locus and sorts
    bases at each site, tests for errors at the site according to error 
    rate, calls consensus.
    """

    ## unpack args
    data, sample, tmpchunk, point = args

    ## read in cluster file 2 lines at a time
    infile = gzip.open(tmpchunk) #sample.files["clusters"])
    duo = itertools.izip(*[iter(infile)]*2)

    ## store read depth info for later output files
    datadict = {}

    ## counters 
    locus = 0
    minsamp_filtered = 0
    nheteros = 0

    ## iterate over clusters
    while 1:
        try: 
            first = duo.next()
        except StopIteration:
            break
        itera = [first[0], first[1]]
        fname = itera[0].split(";")[0]

        ## local containers and counters for this locus"
        locus += 1         ## recording n loci
        sloc = []          ## list for sequence data 
        nloc = []          ## list for names used for gbs filters

        ## grab seqs until end of cluster
        while itera[0] != "//\n":
            ## append sequence * number of dereps "
            nreps = int(itera[0].split(";")[-2].split("=")[1])
            for _ in xrange(nreps):
                sloc.append(tuple(itera[1].strip())) 
                nloc.append(itera[0])
            ## move on to the next sequence
            itera = duo.next()

        ## now that all seqs in this loc are read in 
        ## check that none overlap leftjust overhang if gbs
        if data.paramsdict["datatype"] in ['gbs', 'merged']:
            ## TODO: test these new changes to gbs filter
            ## edge filters
            leftjust = rightjust = None
            rights = []

            ## get leftjust and rights
            for i, j in zip(nloc, sloc):
                leftjust, rights = gbs_edgefilter(i, j, leftjust, rights)
            if rights:
                ## record in name that there was a reverse hit"
                fname = fname[:-2]+"c1"
                try: 
                    rightjust = min([min(i) for i in rights])
                except ValueError:
                    sloc = ""
            
            for seq in xrange(len(sloc)):
                sloc[seq] = sloc[seq][leftjust:]
                if rightjust:
                    sloc[seq] = sloc[seq][:rightjust+1]

        ## Apply depth filter
        if (len(sloc) >= data.paramsdict["mindepth_majrule"]) and \
           (len(sloc) <= data.paramsdict["max_stack_size"]):

            ## this loc passed the minsamp filter
            minsamp_filtered += 1

            ## get stacks of bases at each site
            arrayed = numpy.array(sloc)
            stacked = [Counter(seq) for seq in arrayed.T] 

            ## apply functions to list of sites in stacked
            ## filter by site for paralogs and make consens calls
            consens = [filter2(data, site) for site in stacked]

            ## filtered by locus for paralog
            if "@" not in consens:
                ## get hetero sites
                heteros = [i[0] for i in enumerate(consens) \
                                  if i[1] in list("RKSYWM")]

                ## filter for max number of hetero sites
                exceedmaxploid = 0
                if len(heteros) <= data.paramsdict["max_Hs_consens"]:
                    ## filter for more than x alleles given ploidy. Only 
                    ## relevant if locus is polymorphic at more than one site
                    if len(heteros) > 1:
                        consens, exceedmaxploid = filter3(data, consens,
                                                          heteros, sloc)

                    ## if the locus passed paralog filtering
                    if not exceedmaxploid:

                        consens = "".join(consens).replace("-", "N")
                        ## if a site is stripped then I need to remove the site
                        ## from the site counter (stacked)
                        shortconl = consens.lstrip("N")
                        if len(shortconl) < len(consens):
                            stacked = stacked[-len(shortconl):]
                        shortcon = consens.rstrip("N")
                        if len(shortcon) < len(shortconl):
                            stacked = stacked[:len(shortcon)]                            

                        ## this function which removes low coverage sites next 
                        ## to poly repeats that are likely sequencing errors 
                        ## also edits 'stacked'
                        shortcon, stacked = removerepeat_Ns(shortcon, stacked)

                        ## only allow maxN internal "N"s in a locus
                        if shortcon.count("N") <= int(
                                           data.paramsdict["max_Ns_consens"]):
                            ## minimum length for clustering in vsearch
                            if len(shortcon) >= 32:
                                ## keep for counter
                                nheteros += len(heteros)

                                ## store the consens seq
                                #consdic[fname] = shortcon

                                ## create dataobj w/ name fname
                                dataobj = ObjDict()
                                ## store qual and depth data
                                dataobj.seq = shortcon #[len(cut1):]
                                dataobj.Cs = [i["C"] for i in stacked] 
                                dataobj.As = [i["A"] for i in stacked] 
                                dataobj.Ts = [i["T"] for i in stacked] 
                                dataobj.Gs = [i["G"] for i in stacked]
                                #Cs = [i["C"] for i in stacked] 
                                #As = [i["A"] for i in stacked] 
                                #Ts = [i["T"] for i in stacked] 
                                #Gs = [i["G"] for i in stacked]
                                #dfconsens = pd.DataFrame([list(shortcon), 
                                #                          Cs, As, Ts, Gs])
                                tag = "_".join(fname.split("_")[-2:])
                                datadict[tag] = dataobj
                                #datadict[tag] = dfconsens
                        else:
                            pass #print "maxN filtered loc", locus
                    else:
                        pass #print "ploid filtered loc", locus
                else:
                    pass #print "maxH filtered loc", locus
            else:
                pass #print "third base filtered loc", locus
        else:
            pass #print "mindepth filtered loc", locus

    data.dirs.consens = os.path.join(data.dirs.clusts,
                                     "consens")

    #if not os.path.exists(os.path.join(...)):
    #    os.mkdir(consensdir)
    ## get filename
    consenshandle = "" #os.path.join([consensdir, sample.name+"consens.gz"])

    ## write to file
    with gzip.open(consenshandle, 'wb') as outfile:
        outfile.write("\n".join([">"+sample.name+"_"+obj+"\n"+\
                                 datadict[obj].seq for obj in datadict]))
        #for obj in datadict:
        #    outfile.write(">"+sample.name+"_"+obj+"\n"+datadict[obj].seq+"\n")

    ## count the number of polymorphic sites
    if 'ddrad' in data.paramsdict["datatype"]:
        if 'pair' in data.paramsdict["datatype"]:
            sub = 4 + len(data.paramsdict["restriction_overhang"][0])
        else:
            sub = len(data.paramsdict["restriction_overhang"][0])
    elif 'gbs' in data.paramsdict["datatype"]:
        if 'pair' in data.paramsdict["datatype"]:
            sub = 4 + len(data.paramsdict["restriction_overhang"][0])*2 
            #  (len(params["cut"])*2)
        else:
            sub = len(data.paramsdict["restriction_overhang"][0])
    elif data.paramsdict["datatype"] == "merged":
        sub = len(data.paramsdict["restriction_overhang"][0])*2 
    else:
        sub = len(data.paramsdict["restriction_overhang"][0])
    nsites = sum([len(datadict[i].seq)-sub for i in datadict])
    ldic = len(datadict)
    try: 
        poly = nheteros/float(nsites)
    except ZeroDivisionError:
        poly = 0.

    ## dump the quality score and depth info into a pickle
    #pickleout = gzip.open(handle.replace("clustS.gz", "bindata"), 'wb')
    #pickle.dump(datadict, pickleout)
    #pickleout.close()

    return [sample.name, locus, minsamp_filtered, ldic, nsites, nheteros, poly]



def basecaller(data, site, base1, base2):
    """ inputs data to binomprobr and gets alleles
        correctly oriented """

    ## make statistical base call
    if base1+base2 >= data.paramsdict["mindepth_statistical"]:
        prob, _, who = binomprobr(base1, base2, 
                                     data.stats.error_est.mean(),
                                     data.stats.hetero_est.mean())
    elif base1+base2 >= data.paramsdict["mindepth_majrule"]:
        prob, _, who = simpleconsensus(base1, base2)

    ## if the base could be called with 95% probability
    if float(prob) >= 0.95:
        if who not in "ab":
            ## site is homozygous
            cons = site.most_common(1)[0][0]
        else:
            ## site is heterozygous
            cons = hetero(*[i[0] for i in site.most_common(2)])
    else:
        cons = "N"
    return cons



def filter2(data, site):
    """ site filter, calls basecaller function after filtering to remove
    third base in diploids, and Ns. """

    ## remove Ns counts from site dict
    if "N" in site:
        site.pop("N")

    ## default string to return
    cons = ""

    ## depth of coverage
    depthofcoverage = sum(site.values())

    ## get the most common alleles
    if len(site.values()) > 0:
        base1 = site.most_common(1)[0][1]
    else:
        base1 = base2 = 0
    if len(site.values()) > 1:
        base2 = site.most_common(2)[1][1]
    else:
        base2 = 0
    if len(site.values()) > 2:
        base3 = site.most_common(3)[2][1]
    else:
        base3 = 0

    ## speed hack based on ploidy
    ## if diploid don't allow a third base >20%
    if base3:
        if data.paramsdict["ploidy"] == 2:
            if base3/float(depthofcoverage) > 0.20:
                #print base3
                #print depthofcoverage
                #print base3/float(depthofcoverage)
                #mark as a paralog
                cons = "@"

    if not cons:
        ## if site depth at this point after removing Ns and third bases
        ## is below the limit then no call
        if base1+base2 < data.paramsdict["mindepth_majrule"]:
            cons = "N"
        else:
            ## if depth > 500 reduce to randomly sampled 500 
            if base1+base2 >= 500: 
                firstfivehundred = numpy.array(tuple("A"*base1+"B"*base2))
                numpy.random.shuffle(firstfivehundred)
                base1 = list(firstfivehundred[:500]).count("A")
                base2 = list(firstfivehundred[:500]).count("B")

            ## make the base call using a method depending on depth
            cons = basecaller(data, site, base1, base2)
    return cons



def filter3(data, consens, heteros, sloc):
    """ filters the consensus locus for paralogs based 
    on the user supplied ploidy level"""
    ## TODO: add back in the findalleles function

    ## store matchings of hetero sites
    h_ordered = []

    ## get matchings from each read in sloc
    for read in sloc:
        read_ordered = ""
        for hsite in heteros:
            ssite = read[hsite]
            #if ssite in ["-", "N"]:
            #    read_ordered += "N"
            #else:
            read_ordered += ssite
        h_ordered.append(read_ordered)

    ## count how many times each allele arose
    ## filter out the alleles that have missing or Ns
    counted_h = Counter([i for i in h_ordered if "N" not in i])

    ## exclude if allele occurrence is within the sequencing error rate
    total = len(counted_h.values())
    checkkeys = counted_h.keys() 
    mmm = total*len(sloc[0])*data.stats.error_est.mean()

    ## TODO; does this work if majrule=1?
    for key in checkkeys:
        if counted_h[key] < max(2, round(mmm)):
            counted_h.pop(key)

    ## only allow as many alleles as expected given ploidy
    exceedmaxploid = 0
    if counted_h:
        if len(counted_h) > data.paramsdict["ploidy"]:
            exceedmaxploid = 1
            #print max(2, round(total*len(sloc[0])*params["E"]))
            #print counted_h, 'counted_h after exceedmaxploid'

        ## set the upper vs lower case for bases to save 
        ## the phased order of diploid alleles
        #if params["haplos"] == 2:
        #elif len(counted_h) == 2:
            ## can only phase if two alleles are detected 
            #print counted_h.keys()
            #if all([i[0] != i[1] for i in range(len(counted_h.iterkeys()))])
            #if all([i[0] != i[1] for i in counted_h.keys()]):
            #    print 'yes'
            ## TODO: correct alleles counting for .alleles output
            #if all([counted_h.keys()[0][i] != counted_h.keys()[1][i] \
            #        for i in range(len(counted_h.keys()[0]))]):
            #    consens = findalleles(consens, heteros, counted_h, sloc)
            ## else... could filter if alleles are not detected
            ## after removing low copy possible error alleles
    return consens, exceedmaxploid



# def up_std(handle, mindepth):
#     " function to calculate mean and SD of clustersize"
#     infile = gzip.open(handle)
#     duo = itertools.izip(*[iter(infile)]*2)
#     itera = duo.next()[0].strip()
#     depth = []
#     thisdepth = int(itera.split(";")[1].replace("size=", ""))
#     while 1:
#         try: 
#             itera = duo.next()[0].strip()
#         except StopIteration: 
#             break
#         if itera != "//":
#             thisdepth += int(itera.split(";")[1].replace("size=", ""))
#         else:
#             depth.append(thisdepth)
#             thisdepth = 0
#     infile.close()
#     keep = [i for i in depth if i >= mindepth]
#     if keep:
#         me = numpy.mean(keep)
#         std = numpy.std(keep)
#     else:
#         me = 0.0
#         std = 0.0
#     return me, std



def run_full(data, sample, ipyclient):
    """ split job into bits and passses to the client """

    ## counter for split job submission
    submitted = 0
    num = 0

    ## set optim size for chunks in N clusters
    optim = 1000
    if sample.stats.clusters_kept > 50000:
        optim = 2000
    if sample.stats.clusters_kept > 100000:
        optim = 3000

    ## break up the file into smaller tmp files for each engine
    chunkslist = []
    for clustfile in sample.files.clusters:
        args = [data, clustfile, num, optim]
        _, achunk = zcat_make_temps(args)
        chunkslist += achunk

    import sys
    sys.exit()
    ## send chunks across engines, will delete tmps if failed
    try:
        submitted_args = []
        for tmpchunk in chunkslist:
            ## used to increment names across processors
            point = num*optim  ## 10000
            args = [data, sample, tmpchunk, point]
            submitted_args.append(args)
            submitted += 1
            num += 1

        ## first launch of ipyclient
        def launch():
            """ launch ipyclient """
            lbview = ipyclient.load_balanced_view()
            return lbview

        ## launch within try statement in case engines aren't ready yet
        ## and try 30 one second sleep/wait cycles before giving up on engines
        tries = 30
        while tries:
            try:
                lbview = launch()
                tries = 0
            except ipp.NoEnginesRegistered:
                time.sleep(1)
                tries -= 1

        results = lbview.map_async(consens, submitted_args)

        try:
            results.get()
        except (TypeError, AttributeError, NameError):
            for key in ipyclient.history:
                if ipyclient.metadata[key].error:
                    LOGGER.error("step2 error: %s", 
                        ipyclient.metadata[key].error)
                    raise SystemExit
                if ipyclient.metadata[key].stdout:
                    LOGGER.error("step2 stdout:%s", 
                        ipyclient.metadata[key].stdout)
                    raise SystemExit            
        del lbview

    finally:
        ## if process failed at any point delete tmp files
        for tmpchunk in chunkslist:
            os.remove(tmpchunk)

    return submitted, results


        




def run(data, sample, ipyclient, preview=0, force=False):
    """ checks if the sample should be run and passes the args """

    ## if very few reads, raise warning
    if sample.stats.clusters_kept < 50:
        print("warning: Sample {} has very few reads ({}).".\
            format(sample.name, sample.stats.reads_raw))

    ## if sample is already done skip
    if not force:
        if sample.stats.state >= 5:
            print("Skipping {}. Consensus sequences already called. "\
                 +"Use force=True to overwrite".format(sample.name))
        else:
            submitted, results = run_full(data, sample, ipyclient)
            cleanup(data, sample, submitted, results)
    else:
        submitted, results = run_full(data, sample, ipyclient)
        cleanup(data, sample, submitted, results)



if __name__ == "__main__":
    import ipyrad as ip
    TESTER = ip.load_dataobj("test")
    TESTER.step5.run()
