#!/usr/bin/env python2.7

""" call consensus base calls on single-end data """

from __future__ import print_function
import multiprocessing
import glob
import itertools
import sys
import scipy.stats
import scipy.misc
import numpy
import os
import operator
import gzip
#import cPickle as pickle
from collections import Counter


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
    ab = scipy.misc.comb(base1+base2, base1)/(2.**(base1+base2))
    aa = scipy.stats.binom.pmf(base1, base1+base2, error)
    bb = scipy.stats.binom.pmf(base2, base1+base2, error)

    ## calculate probs
    aa = aa*prior_homo
    bb = bb*prior_homo
    ab = ab*prior_het

    ## return 
    probabilities = [aa, bb, ab]
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
    trans = {('G', 'A'):"R",
             ('G', 'T'):"K",
             ('G', 'C'):"S",
             ('T', 'C'):"Y",
             ('T', 'A'):"W",
             ('C', 'A'):"M"}
    order1 = trans.get((base1, base2))
    order2 = trans.get((base2, base1))
    if order1:
        return order1
    elif order2:
        return order2
    else:
        if [base1, base2].count("N") == 2:
            return "N"
        elif base1 == "N":
            return base2
        elif base2 == "N":
            return base1
        else:
            if base1 == "-":
                return base2
            elif base2 == "-":
                return base1
            #sys.exit("report this error: bad site ("+\
            #         base1+" "+base2+")")


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
        print "THIS ONE HERE"
        print alleles, 'alleles'
        print consens, 'consens'
        print heteros, 'heteros'
        print counted_h, 'counted_h'
        print bigbase, 'bigbase'
        for i in sloc:
            print "".join(i)

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



def gbs_edgefilter(itera, leftjust, rights):
    """ filter edges in case of overlap """

    ## leftjust is seed's left "
    if itera[0].strip()[-1] == ";":
        leftjust = itera[1].index([i for i in itera[1] \
                                   if i not in list("-N")][0])

    ## rightjust is the shortest reverse hit "
    if itera[0].strip()[-1] == "-":
        rights.append(max(-1, [itera[1].rindex(i) for i \
                          in itera[1] if i in list("ATGC")]))
    return leftjust, rights



def consensus(data, sample):
    """
    from a clust file handle,
    reads in all copies at a locus and sorts
    bases at each site, tests for errors at the
    site according to error rate, calls consensus
    """

    ## read in cluster file 2 lines at a time
    infile = gzip.open(sample.files["clusters"])
    duo = itertools.izip(*[iter(infile)]*2)

    ## dataobj for storing base call probabilities
    ## and read depth info for later output files
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
        fname = itera[0].strip().split(";")[0]

        ## edge filters
        leftjust = rightjust = None
        rights = []

        ## local containers and counters for this locus"
        locus += 1         ## recording n loci
        sloc = []          ## list for sequence data 

        ## grab seqs until end of cluster
        while itera[0] != "//\n":
            ## append sequence * number of dereps "
            nreps = int(itera[0].strip().split(";")[1].replace("size=", ""))
            for _ in xrange(nreps):
                sloc.append(tuple(itera[1].strip())) ##.replace("-", "N"))) 
                
            ## record left and right most index of seed and hits (for GBS) "
            if params["datatype"] in ['gbs', 'merged']:
                leftjust, rights = gbs_edgefilter(itera, leftjust, rights)

            ## get the next sequence
            itera = duo.next()

        ## now that all seqs in this loc are read in 
        ## check that none overlap leftjust overhang if gbs
        if params["datatype"] in ['gbs', 'merged']:
            if rights:
                ## record in name that there was a reverse hit"
                fname = "_".join(fname.split("_")[0:-1])+"_c1"
                try: 
                    rightjust = min([min(i) for i in rights])
                except ValueError:
                    sloc = ""
            
            for seq in xrange(len(sloc)):
                sloc[seq] = sloc[seq][leftjust:]
                if rightjust:
                    sloc[seq] = sloc[seq][:rightjust+1]

        ## Apply depth and paralog filters "
        if (len(sloc) >= min(params["lowcounts"], params["mindepth"])) and\
           (len(sloc) <= upper_sd):

            ## this loc passed the minsamp filter
            minsamp_filtered += 1

            ## get stacks of bases at each site
            arrayed = numpy.array(sloc)
            stacked = [Counter(seq) for seq in arrayed.T] 

            ## iterate over sites in stacked
            ## apply functions to list of sites in stacked
            consens = [filter2(params, site) for site in stacked]

            ## filtered by site for paralog
            if "@" not in consens:
                ## filtered by locus for paralog
                heteros = [i[0] for i in enumerate(consens) \
                                  if i[1] in list("RKSYWM")]

                ## filter for max number of hetero sites
                exceedmaxploid = 0
                if len(heteros) < params["maxH"]:
                    ## filter for more than x alleles given ploidy
                    ## only relevant if locus is polymorphic at 
                    ## more than one site
                    if len(heteros) > 1:
                        consens, exceedmaxploid = filter3(params,
                                                          consens,
                                                          heteros,
                                                          sloc)

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
                        ## to poly repeates that are likely sequencing errors 
                        ## also edits 'stacked'
                        shortcon, stacked = removerepeat_Ns(shortcon, stacked)

                        ## only allow maxN internal "N"s in a locus
                        if shortcon.count("N") <= int(params["maxN"]):
                            ## minimum length for clustering in vsearch
                            if len(shortcon) >= 32:
                                ## keep for counter
                                nheteros += len(heteros)

                                ## store the consens seq
                                #consdic[fname] = shortcon

                                ## create dataobj w/ name fname
                                dataobj = Consobj()
                                ## store qual and depth data
                                dataobj.seq = shortcon #[len(cut1):]
                                dataobj.Cs = [i["C"] for i in stacked] 
                                dataobj.As = [i["A"] for i in stacked] 
                                dataobj.Ts = [i["T"] for i in stacked] 
                                dataobj.Gs = [i["G"] for i in stacked] 
                                tag = "_".join(fname.split("_")[-2:])
                                datadict[tag] = dataobj
                        else:
                            pass
                            #print "maxN filtered loc", locus
                    else:
                        pass
                        #print "ploid filtered loc", locus
                else:
                    pass
                    #print "maxH filtered loc", locus
            else:
                pass
                #print "third base filtered loc", locus
        else:
            pass
            #print "mindepth filtered loc", locus

    outfile = gzip.open(handle.replace(".clustS", ".consens"), 'w+')
    for obj in datadict:
        outfile.write(">"+sampname+"_"+obj+"\n"+datadict[obj].seq+"\n")
    #for i in consdic.items():
    #    outfile.write(str(i[0])+'\n'+str(i[1])+"\n")
    outfile.close()
    if not quiet:
        sys.stderr.write(".")

    ## count the number of polymorphic sites
    if 'ddrad' in params["datatype"]:
        if 'pair' in params["datatype"]:
            sub = 4 + len(params["cut"].split(",")[0])
        else:
            sub = len(params["cut"].split(",")[0])
    elif 'gbs' in params["datatype"]:
        if 'pair' in params["datatype"]:
            sub = 4 + (len(params["cut"])*2)
        else:
            sub = len(params["cut"])
    elif params["datatype"] == "merged":
        sub = len(params["cut"])*2
    else:
        sub = len(params["cut"])
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

    return [handle.split("/")[-1], locus, minsamp_filtered, 
                                   ldic, nsites, nheteros, 
                                   round(poly, 7)]



def basecaller(params, site, base1, base2):
    """ inputs data to binomprobr and gets alleles
        correctly oriented """

    ## make statistical base call
    if base1+base2 >= params["mindepth"]:
        prob, mjaf, who = binomprobr(base1, base2, 
                                     float(params["E"]),
                                     float(params["H"]))
    elif base1+base2 >= params["lowcounts"]:
        prob, mjaf, who = simpleconsensus(base1, base2)

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


def filter2(params, site):
    """ site filter """

    ## remove Ns from site
    if "N" in site:
        site.pop("N")

    ## string to return
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
        if params["haplos"] == 2:
            if base3/float(depthofcoverage) > 0.20:
                #print base3
                #print depthofcoverage
                #print base3/float(depthofcoverage)
                #print site
                cons = "@"

    if not cons:
        ## if depth is below the limit then no call
        if base1+base2 < min(params["mindepth"], params["lowcounts"]):
            cons = "N"
        else:
            ## if depth > 500 reduce to 500 for base calling
            if base1+base2 >= 500: 
                firstfivehundred = numpy.array(tuple("A"*base1+"B"*base2))
                numpy.random.shuffle(firstfivehundred)
                base1 = list(firstfivehundred[:500]).count("A")
                base2 = list(firstfivehundred[:500]).count("B")

            ## make the base call using a method depending on depth
            cons = basecaller(params, site, base1, base2)
            ## get the correct base
    return cons



def filter3(params, consens, heteros, sloc):
    """ filters the consensus locus for paralogs based 
    on the user supplied ploidy level"""

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
    for key in checkkeys:
        if counted_h[key] < max(2, round(total*len(sloc[0])*params["E"])):
            counted_h.pop(key)

    ## only allow as many alleles as expected given ploidy
    exceedmaxploid = 0
    if counted_h:
        if len(counted_h) > params["haplos"]:
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



def up_std(handle, mindepth):
    " function to calculate mean and SD of clustersize"
    infile = gzip.open(handle)
    duo = itertools.izip(*[iter(infile)]*2)
    itera = duo.next()[0].strip()
    depth = []
    thisdepth = int(itera.split(";")[1].replace("size=", ""))
    while 1:
        try: 
            itera = duo.next()[0].strip()
        except StopIteration: 
            break
        if itera != "//":
            thisdepth += int(itera.split(";")[1].replace("size=", ""))
        else:
            depth.append(thisdepth)
            thisdepth = 0
    infile.close()
    keep = [i for i in depth if i >= mindepth]
    if keep:
        me = numpy.mean(keep)
        std = numpy.std(keep)
    else:
        me = 0.0
        std = 0.0
    return me, std



def run(data, samples):
    """ calls the main function """
    # load up work queue
    work_queue = multiprocessing.Queue()

    # iterate over files
    submitted = 0

    ## sort samples by clustsize
    samples.sort(key=lambda x: x.stats["clusters_kept"], reverse=True)

    ## put on work queue
    for sample in samples:
        ## only require state 3, param estimates are Assembly averages
        if sample.stats["state"] in [3, 4]:
            if sample.stats["clusters_kept"]:
                work_queue.put([data, sample])
                submitted += 1
        elif sample.stats["state"] >= 5:
            print("sample {} already passed state 5".format(sample.name))

    ## create a queue to pass to workers to store the results
    fake_queue = multiprocessing.Queue()

    try:
        ##  spawn workers
        jobs = []
        for _ in range(min(data.paramsdict["N_processors"], submitted)):
            work = worker.Worker(work_queue, fake_queue, consensus)
            work.start()
            jobs.append(work)
        for job in jobs:
            job.join()

    finally:
        ## get results and cleanup
        cleanup(data, samples)



if __name__ == "__main__":
    run()
