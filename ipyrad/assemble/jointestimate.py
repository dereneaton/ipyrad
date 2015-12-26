#!/usr/bin/env python2

""" jointly infers heterozygosity and error rate from stacked sequences """

from __future__ import print_function
import scipy.stats
import scipy.optimize
import numpy as np
import itertools
import os
import gzip
from collections import Counter

# pylint: disable=E1101


def get_freqs(stack):
    """ returns a list as frequencies for ATGC"""
    sump = sum([sum(cc.values()) for cc in stack])
    #sump = sum([sum(i) for i in site])
    totalcount = Counter()
    for stackcount in stack:
        totalcount += stackcount
    basefreqs = np.array([totalcount["A"],
                          totalcount["T"],
                          totalcount["G"],
                          totalcount["C"]])/float(sump)
    return basefreqs


def likelihood1(errors, base_frequencies, stacks):
    """probability homozygous"""
    ## make sure base_frequencies are in the right order
    #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    totals = np.array([stacks.sum(axis=1)]*4).T
    prob = scipy.stats.binom.pmf(totals-stacks, totals, errors)
    return np.sum(base_frequencies*prob, axis=1)



def likelihood2(errors, base_frequencies, stacks):
    """probability of heterozygous"""
    returns = np.zeros([len(stacks)])
    for idx, stackl in enumerate(stacks):
        spair = list(itertools.combinations(stackl, 2))
        bpair = list(itertools.combinations(base_frequencies, 2))
        one = 2.*np.product(bpair, axis=1)
        tot = stackl.sum() #np.sum(spair, axis=1)
        atwo = tot - np.array([i[0] for i in spair]) -\
                     np.array([i[1] for i in spair])
        two = scipy.stats.binom.pmf(atwo, tot, (2.*errors)/3.)
        three = scipy.stats.binom.pmf(np.array([i[0] for i in spair]),
                                      np.array([i[0]+i[1] for i in spair]),
                                      0.5)
        four = 1.-np.sum(base_frequencies**2)
        returns[idx] = np.sum(one*two*(three/four))
    return np.array(returns)



def get_diploid_lik(starting_params, base_frequencies, stacks, stackcounts):
    """ Log likelihood score given values [H,E] """
    hetero = starting_params[0]
    errors = starting_params[1]
    if (hetero <= 0.) or (errors <= 0.):
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = (1.-hetero)*likelihood1(errors, base_frequencies, stacks)
        lik2 = (hetero)*likelihood2(errors, base_frequencies, stacks)
        liks = lik1+lik2
        logliks = np.log(liks[liks > 0])*stackcounts[liks > 0]
        #logliks = np.log(liks)*stackcounts
        score = -logliks.sum()
    return score


def get_haploid_lik(errors, base_frequencies, tabled_stacks):
    """ Log likelihood score given values [E] """
    hetero = 0.
    listofliks = np.zeros([len(tabled_stacks)])
    if errors <= 0.:
        score = np.exp(100)
    else:
        for idx, uniqstack in enumerate(tabled_stacks):
            loglik = ((1.-hetero)*\
                     likelihood1(errors, base_frequencies, uniqstack))+\
                     (hetero*likelihood2(errors, base_frequencies, uniqstack))
            if loglik > 0:
                #listofliks.append(tabled_stacks[uniqstack]*np.log(loglik))
                listofliks[idx] = tabled_stacks[uniqstack]*np.log(loglik)
        score = -sum(listofliks)
    return score



def tabledstack(stack):
    """ makes a dictionary with counts of base counts [x,x,x,x],
        greatly speeds up Likelihood calculation"""
    countedstacks = []
    for stackcount in stack:
        ## convert to string for use with counter
        countedstacks.append(str([stackcount["A"],
                                  stackcount["T"],
                                  stackcount["G"],
                                  stackcount["C"]]))
    return Counter(countedstacks)



def countlist(data, sample, subsample=None):
    """ makes a list of lists of reads at each site """
    infile = gzip.open(sample.files.clusters)
    duo = itertools.izip(*[iter(infile)]*2)
    stacked = []

    ## speed hack for subsampling
    if subsample:
        until = int(subsample)
    else:
        until = int(1e6)

    while len(stacked) < until:
        try:
            itera = duo.next()
        except StopIteration:
            break
        #itera = [first[0], first[1]]
        thisstack = []
        while itera[0] != "//\n":
            nreps = int(itera[0].split(";")[1][5:])

            ## record left and right most for cutting if gbs merge data
            ## ....
            ## append sequence * number of dereps
            for _ in range(nreps):
                thisstack.append(tuple(itera[1].strip()))
            itera = duo.next()

        ## don't use sequence edges / restriction overhangs
        cutl1, cutl2 = [len(i) for i in data.paramsdict["restriction_overhang"]]
        if cutl2:
            cutl2 *= -1
        else:
            cutl2 = None
        sarray = np.array(thisstack)[:, cutl1:cutl2]

        ## enforce minimum depth for estimates
        if sarray.shape[0] >= data.paramsdict["mindepth_statistical"]:
            ## make list for each site in sequences
            res = [Counter(seq) for seq in sarray.T]
            ## exclude sites with indels
            stacked += [i for i in res if "-" not in i]
    return stacked



def toarray(uniqstack):
    """ converts string lists to arrays"""
    return np.array([int(i) for i in uniqstack[1:-1].split(',')])



def optim(args):
    """ func scipy optimize to find best parameters"""

    ## split args
    data, sample, subsample = args

    ## make a list of Counter objects for each site in each stack
    stacked = countlist(data, sample, subsample)

    ## get base frequencies
    base_frequencies = get_freqs(stacked)

    ## get tabled counts of base patterns
    tabled_stacks = tabledstack(stacked)

    ## put into array
    stacks = np.array([toarray(i) for i in tabled_stacks.keys()])
    stackcounts = np.array(tabled_stacks.values())
    del stacked

    ## if data are haploid fix H to 0
    if data.paramsdict["ploidy"] == 1:
        starting_params = [0.001]
        hetero = 0.
        errors = scipy.optimize.fmin(
                                get_haploid_lik,
                                starting_params,
                                (base_frequencies, 
                                    stacks,
                                    stackcounts),
                                disp=False,
                                full_output=False)
    else:
        starting_params = [0.01, 0.001]
        hetero, errors = scipy.optimize.fmin(get_diploid_lik,
                                             starting_params,
                                             (base_frequencies, 
                                                stacks,
                                                stackcounts),
                                             disp=False,
                                             full_output=False)
    ## store result in a tempfile
    handle = os.path.join(os.path.dirname(sample.files["clusters"]),
                          "tmp_"+sample.name+".het")
    with open(handle, 'wb') as out:
        out.write("{}\t{}\t{}".format(sample.name, hetero, errors))

    return [stacks, stackcounts]



def run(data, samples, subsample, force, ipyclient):
    """ calls the main functions """

    # if haploid data
    if data.paramsdict["ploidy"] == 1:
        print("Applying haploid-based test (infer E with H fixed to 0).")

    submitted_args = []
    ## if sample is already done skip
    for sample in samples:
        if not force:
            if sample.stats.state >= 4:
                print(sample.name+"already estimated. Use force=True "+\
                      "to overwrite.")
            elif sample.stats.state < 3:
                print(sample.name+"not clustered yet. Run step3() first.")
            elif sample.stats.clusters_hidepth < 100:
                print("skipping {}. Too few reads ({}). Use force=True "+\
                     "to override".format(sample.name, sample.stats.reads_raw))
            else:
                submitted_args.append([data, sample, subsample])
        else:
            if sample.stats.state < 3:
                print(sample.name+" not clustered yet. Run step3() first.")
            else:
                submitted_args.append([data, sample, subsample])

    ## if jobs then run
    if submitted_args:
        ## sort by cluster size
        submitted_args.sort(key=lambda x: x[1].stats.clusters_hidepth, 
                                                      reverse=True)
        lbview = ipyclient.load_balanced_view()
        results = lbview.map_async(optim, submitted_args)
        results.get()

        ## get results and remove temp files
        cleanup(data, samples)



def cleanup(data, samples):
    """ stores results to the Assembly object, writes to stats file, 
    and cleans up temp files """

    ## collect all samples that have finished
    for sampobj in samples:
        tempres = os.path.join(
                      os.path.dirname(
                          sampobj.files["clusters"]),
                          "tmp_"+sampobj.name+".het")

        ## if sample finished
        if os.path.exists(tempres):
            ## get stats
            _, est_h, est_e = open(tempres).read().split("\t")

            ## sample assignments
            sampobj.stats["state"] = 4
            sampobj.stats["hetero_est"] = est_h
            sampobj.stats["error_est"] = est_e
            data._stamp("s4 params estimated on "+sampobj.name)
            os.remove(tempres)
        else:
            pass
            #print("{} did not finish".format(sampobj.name))


    ## create a final stats file
    #statsfile = os.path.join(os.path.dir(""))
    #data.statsfiles["s4"] = statsfile


    #if not os.path.exists(data.paramsdict["working_directory"]+\
    #    "stats/Pi_E_estimate.txt"):
    #     outstats = open(params["work"]+"stats/Pi_E_estimate.txt", 'w')
    #     outstats.write("taxa\tH\tE\n")
    # else:
    #     outstats = open(params["work"]+"stats/Pi_E_estimate.txt", 'a')

    # ## remove stats temp files
    # for handle in funcfiles:
    #     end = handle.split("/")[-1].replace(".clustS.gz", "")
    #     tempfile = params["work"]+"stats/."+end+".temp"
    #     line = open(tempfile).readlines()
    #     outstats.write(line[0])
    #     os.remove(tempfile)
    # outstats.close()


if __name__ == "__main__":

    import ipyrad as ip

    ## get path to test dir/ 
    ROOT = os.path.realpath(
       os.path.dirname(
           os.path.dirname(
               os.path.dirname(__file__)
               )
           )
       )

    ## run test on RAD data1
    TEST = ip.load_assembly(os.path.join(ROOT, "tests", "test_rad", "data1"))
    TEST.step4(force=True)
    print(TEST.stats)

    ## run test on messy data set
    #TEST = ip.load_assembly(os.path.join(ROOT, "tests", "radmess", "data1"))

    ## check if results are correct

    ## cleanup

