#!/usr/bin/env python

"popgen tools"

from __future__ import print_function, division
from itertools import chain

import glob
import h5py
import itertools
import math
import numpy as np
import pandas as pd
import os
import pickle
import time
from collections import Counter
from functools import reduce
from ipyrad import Assembly
from ipyrad.assemble.utils import IPyradError
from itertools import combinations
from scipy.stats import entropy, hmean
from ..core.Parallel import Parallel
from .locus_extracter import LocusExtracter
from .utils import Params, ProgressBar
from ..assemble.utils import DCONS

_BAD_IMAP_ERROR = """
Samples in imap not in the hdf5 file: {}"
"""

_SKIP_SAMPLES_WARN = """
Skipping samples in hdf5 not present in imap: {}"
"""

class Popgen(object):
    """
    Analysis functions for calculating theta, Fst, Fis, thetaW, etc.

    Some functions follow Ferretti et al 2012 for calculating stats while
    accounting for missing data:

    Ferretti, L., Raineri, E., & Ramos-Onsins, S. (2012). Neutrality tests for
    sequences with missing data. Genetics, 191(4), 1397-1401.

    Another useful resource for calculating sumstats with missing data:

    Korunes, K., & Samuk, K. (2021). pixy: Unbiased estimation of nucleotide
    diversity and divergence in the presence of missing data. Molecular Ecology
    Resources.
    """

    def __init__(
        self,
        data,
        imap=None,
        minmap=None,
        mincov=4,
        workdir="analysis-popgen",
        quiet=False,
        ):
        
        # set attributes
        self.imap = (imap if imap else {})
        self.minmap = minmap #(minmap if minmap else {i: 4 for i in self.imap})
        self.mincov = mincov
        self.npops = (len(self.imap) if imap else 1)
        self.quiet = quiet
        self.nboots = 100

        # Data from the LocusExtracter
        self.seqs = ''

        # SNP data from snps.hdf5["snps"]
        # This is the way I started with, but maybe not the best approach
        # get rid of this when you're sick of looking at it once the seqs
        # approach is working.
        self.snps = pd.DataFrame()

        # i/o paths
        self.workdir=os.path.abspath(workdir)
        self.mapfile = ""
        self._check_files(data)
        self._check_samples()
        #self.maparr = np.zeros()

        # parallelization
        self.ipcluster = {
            "cluster_id": "",
            "profile": "default",
            "engines": "Local",
            "quiet": 0,
            "timeout": 60,
            "cores": 0,
            "threads": 1,
            "pids": {},
            }

        # params to pass to the processor
        self.params = Params()
        self.params["imap"] = self.imap
        self.params["workdir"] = self.workdir

        # results dataframes
        self.results = Params()

        # maybe these aren't necessary?
        # pairwise Fst between all populations
#        arrfst = np.zeros((self.npops, self.npops), dtype=np.uint64)
#        self.results.fst = pd.DataFrame(
#            arrfst
#            )

        # individual pi 
#        arrpi = np.zeros(self.npops, dtype=np.uint64)
#        self.results.pi = pd.DataFrame(
#            arrpi
#            )

        # population thetas
#        arrtheta = np.zeros(self.npops, dtype=np.uint64)
#        self.results.theta = pd.DataFrame(
#            arrtheta
#            )


    def _check_files(self, data):
        "check input files and file paths"

        if isinstance(data, Assembly):
            try:
                # Since v0.9.63-ish the snps_database is stored as an
                # assembly parameter. If it's not there try to construct
                # it from data.outfiles
                self.datafile = data.seqs_database
            except AttributeError:
                self.datafile = os.path.join(data.dirs.outfiles + 
                                            f"{data.name}.seqs.hdf5")
        else:
            self.datafile = data

        # check data file
        if os.path.exists(self.datafile):
            self.datafile = os.path.realpath(self.datafile)
        else:
            raise IPyradError(
                "data file does not exist. Check path: {}"
                .format(self.datafile))

        # check map file
        #if self.mapfile:
        #    self.mapfile = os.path.realpath(self.mapfile)

        # check workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # load the snp data
        # Not using the snp data internally at all. iao 5/2023
        #with h5py.File(self.snpfile, 'r') as io5:
        #    for idx, name in enumerate(io5["snps"].attrs["names"]):
        #        self.snps[name.decode("utf-8")] = io5["snps"][idx]
        #    # TODO This is temporary to keep _fst running with the snps data
        #    self.data = self.snps


    def _check_samples(self):
        "Read in list of sample names from the datafile"

        with h5py.File(self.datafile, 'r') as io5:
            try:
                self.samples = [i.decode() for i in io5["phymap"].attrs["phynames"]]
            except AttributeError:
                self.samples = [i for i in io5["phymap"].attrs["phynames"]]
        if self.imap:
            # Check agreement between samples in imap and hdf5 file
            imap_samps = list(chain(*self.imap.values()))
            in_imap_not_hdf5 = set(imap_samps).difference(self.samples)
            in_hdf5_not_imap = set(self.samples).difference(imap_samps)
            if in_imap_not_hdf5:
                # Error if you pass in a sample in imap not in hdf5
                raise IPyradError(_BAD_IMAP_ERROR.format(in_imap_not_hdf5))
            if in_hdf5_not_imap:
                if not self.quiet:
                    # Warn here because this is a valid way to remove samples
                    print(_SKIP_SAMPLES_WARN.format(in_hdf5_not_imap))
        else:
            # If no imap then all samples are from one population
            self.imap["pop1"] = self.samples


    def run(self, ipyclient=None, force=False, show_cluster=True, auto=False):
        """
        Submits popgen jobs to run on a cluster (ipyparallel Client). An
        ipyclient connection is optional. If no ipyclient then it runs
        serially on one core.

        Parameters:
        -----------
        ipyclient: (type=ipyparallel.Client); Default=None.
            If you started an ipyclient manually then you can
            connect to it and use it to distribute jobs here.

        force: (type=bool); Default=False.
            Force overwrite of existing output with the same name.

        show_cluster: (type=bool); Default=False.
            Print information about parallel connection.

        auto: (type=bool); Default=False.
            Let ipyrad automatically manage ipcluster start and shutdown.
            This will connect to all avaiable cores by default, but can
            be modified by changing the parameters of the .ipcluster dict
            associated with this tool.
        """
        pool = Parallel(
            tool=self,
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force},
            )
        pool.wrap_run()


    def _run(self, force=False, ipyclient=None):

        self.asyncs = []

        lbview = ipyclient.load_balanced_view()

        # apply locus extracter filtering
        self.lex = LocusExtracter(
            data=self.datafile,
            imap=self.imap,
            minmap=self.minmap,
            mincov=self.mincov,
#            minsnps=self.minsnps,
#            maxmissing=self.maxmissing,
#            minlen=self.minlen,
        )

        # Extract loci from seqs data
        self.lex.run(ipyclient=ipyclient, force=True, show_cluster=False)

        # For each locus, calculate all the sumstats of interest
        nloci = len(self.lex.loci)
        chunksize = int(np.round(nloci/len(ipyclient)))
        chunkstart = 0
        while True:
            # Test for the last chunk and pare it back to the correct size
            if nloci < (chunkstart + chunksize):
                chunksize = nloci - chunkstart
            # Carve off big chunks to pass out to the cluster engines
            chunk = [self.lex.get_locus(x, as_df=True) for x in range(chunkstart,
                                                                      chunkstart+chunksize)]
            args = [self.params, chunkstart, chunk]
            rasync = lbview.apply(_calc_sumstats, *args)
            self.asyncs.append(rasync)

            chunkstart += chunksize
            if chunkstart >= nloci: break

        # setup progress bar
        prog = ProgressBar(nloci, None, "Calculating sumstats for nloci {}".format(nloci))
        prog.finished = 0
        prog.update()

        # block until jobs are done with a progress bar.
        while 1:
            # break between checking progress
            prog.update()
            time.sleep(5)
            # calc finished jobs
            finished = [i.ready() for i in self.asyncs]
            if not all(finished):
                prog.finished = len(finished)
            else:
                # all jobs finished
                prog.finished = nloci
                prog.update()
                print("")
                break

        # Check failed jobs
        for job in self.asyncs:
            if not job.successful():
                job.get()

        self._collect_results()


    def _collect_results(self):
        "Collect results from Processor and write results files."

        # organize stats averaged across all loci into dataframes
        single_popstats = pd.DataFrame(
            columns=["mean_pi", "mean_Watterson", "mean_TajimasD"],
            index=self.imap.keys(),
        )

        avg_dxy = pd.DataFrame(
            columns=self.imap.keys(),
            index=self.imap.keys(),
        )
        avg_fst = pd.DataFrame(
            columns=self.imap.keys(),
            index=self.imap.keys(),
        )

        # load pickled dictionaries into a dict
        pickles = glob.glob(os.path.join(self.workdir, "*.p"))

        pdicts = {}
        for pkl in pickles:
            with open(pkl, 'rb') as inp:
                pdicts[int(pkl.rsplit("/", 1)[-1][:-2])] = pickle.load(inp)

        # Make one giant dictionary with all the results
        full_res = {}
        for d in [pdicts]: full_res.update(d)

        # Sort locus ids so output files will be in locus order
        pidx = sorted(full_res.keys())

        pi_dict = {}
        w_theta_dict = {}
        tajd_dict = {}

        dxy_dict = {}

        # Unpack the results from each chunk per summary statistic
        for idx in pidx:
            pi_dict.update(full_res[idx]["pi"])
            w_theta_dict.update(full_res[idx]["Watterson"])
            tajd_dict.update(full_res[idx]["TajimasD"])
            # Collect pi stats

        prog = ProgressBar(len(self.imap), None, "Collating sumstats for npops {}".format(len(self.imap)))
        prog.finished = 0
        prog.update()

        # Fill one dataframe per population for within population stats per locus
        popstats = {}
        site_pis = {}
        for pidx, pop in enumerate(self.imap):
            site_pis[pop] = {}
            popstats[pop] = pd.DataFrame([], columns=["pi",
                                                        "raw_pi",
                                                        "Watterson",
                                                        "raw_Watterson",
                                                        "TajimasD"], index=range(len(self.lex.loci)))
            for lidx in range(len(self.lex.loci)):
                popstats[pop]["pi"].loc[lidx] = pi_dict[lidx][pop]["pi_per_base"]
                popstats[pop]["raw_pi"].loc[lidx] = pi_dict[lidx][pop]["pi"]
                popstats[pop]["Watterson"].loc[lidx] = w_theta_dict[lidx][pop]["w_theta_per_base"]
                popstats[pop]["raw_Watterson"].loc[lidx] = w_theta_dict[lidx][pop]["w_theta"]
                popstats[pop]["TajimasD"].loc[lidx] = tajd_dict[lidx][pop]

                # Accumulate site pis per locus
                site_pis[pop][lidx] = pi_dict[lidx][pop]["site_pi"]

            with open(os.path.join(self.workdir, "{}-sites.pi".format(pop)), 'w') as sites_file:
                sites_file.write("locus\tsite\tpi\n")
                for lidx in site_pis[pop]:
                    for k, v in site_pis[pop][lidx].items():
                        sites_file.write("\t".join([str(x) for x in [lidx, k, v]]) + "\n")

            # Simple means for quick reporting
            single_popstats["mean_pi"][pop] = popstats[pop]["pi"].mean()
            single_popstats["mean_Watterson"][pop] = popstats[pop]["Watterson"].mean()
            single_popstats["mean_TajimasD"][pop] = popstats[pop]["TajimasD"].mean()

            prog.finished = pidx
            prog.update()
        self.results["within"] = single_popstats

        # Pull between population sumstats
        # Each full_res key contains a dictionary of dataframes showing pairwise
        # dxy among all populations. Extract all the dataframes and flatten them to a
        # list of dataframes
        dxys_per_locus = [y for x in full_res for y in full_res[x]["Dxy"].values()]
        # Add all dataframes together and divide by the number of dataframes to
        # get average dxy among pops
        dxys = reduce(lambda x, y: x.add(y, fill_value=0), dxys_per_locus)/len(dxys_per_locus)
        self.results["between"] = {}
        self.results["between"]["Dxy"] = dxys
 
        # Repeat the same procedure for Fst
        # There are actually 3 versions of Fst: Fst, Fst_adj, and Fst_Nm
        fsts_per_locus = [y["Fst"] for x in full_res for y in full_res[x]["Fst"].values()]
        fsts = reduce(lambda x, y: x.add(y, fill_value=0), fsts_per_locus)/len(fsts_per_locus)
        self.results["between"]["Fst"] = fsts
        # Repeat for Fst_adj
        fsts_per_locus = [y["Fst_adj"] for x in full_res for y in full_res[x]["Fst"].values()]
        fsts = reduce(lambda x, y: x.add(y, fill_value=0), fsts_per_locus)/len(fsts_per_locus)
        self.results["between"]["Fst_adj"] = fsts

        prog.finished = len(self.imap)
        prog.update()
        print("")


# ------------------------------------------------------------
# Classes initialized and run on remote engines.
# ------------------------------------------------------------
def _calc_sumstats(params, start_lidx, loci):
    # process chunk writes to files and returns proc with features.
    proc = Processor(params, start_lidx, loci)
    proc.run()

    with open(proc.outfile, 'wb') as outpickle:
        pickle.dump(proc.results, outpickle)


##############################################################

class Processor(object):
    def __init__(self, params, start_lidx, loci):
        """
        """

        # init data
        self.imap = params.imap
        self.outfile = os.path.join(params.workdir, "{}.p".format(start_lidx))
        self.start_lidx = start_lidx
        self.nloci = len(loci)
        self.loci = iter(loci)

        self.results = Params()
        self.results.pi = {}
        self.results.Watterson = {}
        self.results.TajimasD = {}
        self.results.Fst = {}
        self.results.Dxy = {}


    def run(self):
        lidx = self.start_lidx
        # iterate through loci in the chunk
        while 1:
            try:
                locus = next(self.loci)
                self.results.pi[lidx] = {}
                self.results.Watterson[lidx] = {}
                self.results.TajimasD[lidx] = {}

                # within pops stats
                for pop in self.imap:
                    # Carve off just the samples for this population
                    try:
                        # The locus may not have data for all samples in the population
                        # so `intersection` retains the sample names common to the locus
                        # index and the samples in the imap pop
                        cts, sidxs, length = self._process_locus(
                                                    locus.loc[locus.index.intersection(self.imap[pop])])
                    except KeyError:
                        raise Exception("Error in Processor.run() lidx: {}".format(lidx))

                    # Number of segregating sites
                    S = len(sidxs)
                    # Number of samples
                    n = len(locus)
                    pi_res = self._pi(cts, sidxs, length)
                    w_theta_res = self._Watterson(S=S, n=n, length=length)
                    tajD_res = self._TajimasD(S=S, n=n,
                                                pi=pi_res["pi"],
                                                w_theta=w_theta_res["w_theta"])
                    # store results
                    self.results.pi[lidx][pop] = pi_res
                    self.results.Watterson[lidx][pop] = w_theta_res
                    self.results.TajimasD[lidx][pop] = tajD_res

                # between pops stats
                Dxy_arr = pd.DataFrame(
                    data=np.zeros((len(self.imap), len(self.imap))),
                    index=self.imap.keys(),
                    columns=self.imap.keys(),
                )
                for pops in combinations(self.imap, 2):
                    try:
                        pop_cts, sidxs = self._process_locus_pops(locus, pops)
                        Dxy_res = self._dxy(*pop_cts.values(), len(locus))
                        Dxy_arr[pops[1]][pops[0]] = Dxy_res
                    except KeyError:
                        pass
                self.results.Dxy[lidx] = Dxy_arr

                try:
                    Fst_res = self._fst_full(locus)
                    self.results.Fst[lidx] = {}
                    self.results.Fst[lidx]["Fst"] = Fst_res[0]
                    self.results.Fst[lidx]["Fst_adj"] = Fst_res[1]
                    self.results.Fst[lidx]["Fst_Nm"] = Fst_res[2]
                except Exception:
                    # Does not work if any pops have no samples
                    # Doesn't work well with missing data
                    pass

                lidx += 1
            except StopIteration:
                break


    def _process_locus(self, locus):
        """
        Helper function to convert a locus into a list of counters per
        variable site, the indices of each variable site, and the length
        of the locus.
        """
        # Count numbers of unique bases per site
        # Don't consider indels (45) and N (78). Similar to how pixy
        # does it and should result in less biased pi estimates.
        cts = np.array(locus.apply(lambda bases:\
                        Counter(x for x in bases if x not in [45, 78])))
        # Only consider variable sites
        snps = np.array([len(x) for x in cts]) > 1
        # Indexes of variable sites
        sidxs = np.where(snps)[0]

        # Split consensus bases and disambiguate (convert diploid to haploid)
        # Only retain this info for variable sites
        var_cts = cts[snps]

        # An inner function to split alleles for base counts within
        # a given site. Will disambiguate ambiguous bases and sum all CATGs:
        # Counter({82: 5, 71: 4, 84: 3})
        # The Rs (82) gets split up into As (65) and Gs (71)
        # Counter({71: 13, 65: 5, 84: 6})
        def dcons(counter):
            return Counter(list(itertools.chain(*[DCONS[x]*ct for x, ct in counter.items()])))
        # Map this function across all counters per variable site.
        # There's probably a better way to do this, but this works.
        var_cts = list(map(dcons, var_cts))
        return var_cts, sidxs, len(cts)


    def _process_locus_pops(self, locus, pops):
        """
        A helper function to return counts at variable sites between two
        (or more) populations. There's probably a way to hack the
        _process_locus() function to do this, but it would get ugly fast.

        :param array-like locus: The locus to process.
        :param list pops: The two or more populations to consider. These
            should be valid keys in the self.imap dictionary.

        :returns tuple: A dict with one record per population containing
            a counter at all shared variable sites. The second element of the
            tuple is a sorted list of all shared variable sites.
        """
        pop_cts = {}
        pop_sidxs = {}

        for pop in pops:
            # Get counts for this pop
            cts = np.array(locus.loc[self.imap[pop]].apply(
                            lambda bases: Counter(x for x in bases if x not in [45, 78])))
            # Only consider variable sites
            snps = np.array([len(x) for x in cts]) > 1
            # store indexes of variable sites and counts for all bases
            pop_cts[pop] = cts
            pop_sidxs[pop] = np.where(snps)[0]

        # variable sites across all pops. Compile all the sidxs per population
        # make a set and sort it. A little ugly.
        sidxs = sorted(list(set((itertools.chain(*[pop_sidxs[x] for x in pop_sidxs])))))

        # Now go back through and subset the counts to retain all shared
        # variable sites. You have to do it this way because some sites may be
        # invariable within a population, but you still want counts while
        # accounting for missing data and indels.
        for pop in pops:
            pop_cts[pop] = pop_cts[pop][sidxs]

        # pop_cts here will be a dictionary with list of equal length for all
        # populations including counts per base at each sidxs site.
        return pop_cts, sidxs


    # Within population summary statistics
    # The semi-private (single "_") methods are driven by the run() method
    # and require super specific, individualized data which is parsed and
    # passed to them by run(). They will be difficult to call by hand, but
    # will be much more efficient in aggregate than the public methods below.
    def _pi(self, cts, sidxs, length):
        """
        :param array-like cts: A list of counters for each variable site.
        :param array-like sidxs: List of base positions of variable sites
        """
        pi = 0
        site_pi = {}
        for idx, sidx in enumerate(sidxs):
            site_pi[sidx] = 0
            # Enumerate the possible comparisons and for each
            # comparison calculate the number of pairwise differences,
            # summing over all sites in the sequence.
            for c in combinations(cts[idx].values(), 2):
                n = c[0] + c[1]
                n_comparisons = float(n) * (n - 1) / 2
                site_pi[sidx] += float(c[0]) * (n-c[0]) / n_comparisons
        if site_pi:
            # Average over the length of the whole sequence.
            pi = sum(site_pi.values())

        return {"pi":pi, "pi_per_base":pi/length,  "site_pi":site_pi}


    def _Watterson(self, S, n, length):
        w_theta = S/(n*(1/hmean(list(range(1, n+1)))))
        return {"w_theta":w_theta, "w_theta_per_base":w_theta/length}


    def _TajimasD(self, S, n, pi, w_theta):
        D = 0
        ddenom = 0
        if S > 0:
            d_num = pi - w_theta
            ddenom = self._TajimasD_denom(S, n)
        if ddenom != 0:
            D = d_num/ddenom
        return D


    def _TajimasD_denom(self, S, n):
        """
        Tajima's D denominator. I toiled over this to get it right and it is
        known to be working.

        This page has a nice worked example with values for each
        subfunction so you can check your equations:
        https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-
        quantitative-genomics-fall-2005/study-materials/tajimad1.pdf

        :param int N: The number of samples
        :param int S: The number of segregating sites.
        """
        b1 = (n+1)/float(3*(n-1))
        a1 = sum([1./x for x in range(1, n)])
        c1 = b1 - (1./a1)
        e1 = c1/a1
        a2 = sum([1./(x**2) for x in range(1, n)])
        b2 = (2.*(n**2 + n + 3))/(9*n*(n-1))
        c2 = b2 - (n+2)/(a1*n) + (a2/(a1**2))
        e2 = c2/(a1**2+a2)
        ddenom = math.sqrt(e1*S + e2*S*(S-1))
        return ddenom


    # The "public" methods are wholly independent and operate at the level
    # of a given locus. They are also substantially redundant, so they may be
    # called by hand, but the bulk of the work is done by the semi-private
    # (single "_") methods above.
    def pi(self, locus):
        """
        Calculate nucleotide diversity per site and also average per base.

        :param array-like locus: An np.array or pd.DataFrame of aligned loci
            as might be returned by calling LocusExtracter.get_locus(as_df=True).

        :return tuple: Returns a tuple with raw pi, pi_per_base (averaged across
            the length of the whole locus), and a dictionary containing values of
            pi per site, with keys as the base positions.
        """
        cts, sidxs, length = self._process_locus(locus)
        return self._pi(cts, sidxs, length)


    def Watterson(self, locus):
        """
        Calculate Watterson's theta and optionally average over sequence length.
    
        :param array-like locus: The DNA sequence(s) over which to
            calculate the statistic. This should be formatted in the same way
            as the result from a call to LocusExtracter.get_locus(), i.e. as
            an array or DataFrame with bases coded as int8 ascii values.
    
        :return tuple: The value of Watterson's estimator of theta, both the
            raw value and scaled to per base.
        """
        n = len(locus)
        cts, sidxs, length = self._process_locus(locus)
        return self._Watterson(S=len(sidxs), n=n, length=length)


    def TajimasD(self, locus):
        """
        Calculate Tajima's D for a given locus.

        :param array-like locus: Locus data in the same format as the other
            functions.

        :return float: Tajima's D calculated on the data for this locus.
        """
        n = len(locus)
        cts, sidxs, length = self._process_locus(locus)
        S = len(sidxs)
        pi = self._pi(cts, sidxs, length)["pi"]
        w_theta = self._Watterson(S, n, length)["w_theta"]
        return self._TajimasD(S, n, pi, w_theta)


    # Between population summary statistics
    def _dxy(self, cts_a, cts_b, length):
        """
        Calculate Dxy, the absolute sequence divergence, between two
        populations. The input are counts of each base within each population
        at sites that vary in either or both populations.

        :param list cts_a: A list of counters for each base that varies in
            either population for population 1.
        :param list cts_b: Same thing for population 2.
        :param int length: The length of the sequence to calculate Dxy over

        :return float: The raw Dxy value (unscaled by sequence length).
        """
        Dxy = 0
        ncomps = 0
        for cta, ctb in zip(cts_a, cts_b):
            ncomps += sum(list(cta.values())) *\
                        sum(list(ctb.values()))

            for ka, va in cta.items():
                for kb, vb in ctb.items():
                    if ka == kb: continue
                    Dxy += va*vb
        # Monomorphic locus
        if ncomps == 0: return 0
        return Dxy/ncomps/length


    def _fst_full(self, locus):
        """
        This function operates a bit differently than the others because it is
        largely code that Deren implemented. It works great, but it operates on
        the full locus sequence, rather than only the variable sites, and it
        does all pairwise comparisons at once rather than one at a time.

        The other thing to remember about this version is that it isn't
        currently  accounting for '-' and 'N' in the calculations, which we
        are explicitly doing with the other sumstats (e.g. pi/Dxy).

        Calculate population fixation index Fst using Hudson's estimator.
        Hudson et al. (1992) "Estimation of Levels of Gene Flow From DNA 
        Sequence Data", also returns Fstd with correction for the number of 
        subpopulations (using the number sampled, since the true number is 
        unknown) and number of migrants (Nm) derived from Li (1976b) as 
        described in the Hudson et al. (1992) paper. 

        See also Bhatia et al 2013 Estimating and interpreting FST: The impact
        of rare variants, though this formulation isn't implemented here.
        """       
        # init fst matrix df with named rows and cols
        npops = len(self.imap)
        farr = pd.DataFrame(
            data=np.zeros((npops, npops)),
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        darr = pd.DataFrame(
            data=np.zeros((npops, npops)),
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        narr = pd.DataFrame(
            data=np.zeros((npops, npops)),
            index=self.imap.keys(),
            columns=self.imap.keys(),
        ) 
        d = npops

        # iterate over pairs of pops and fill Fst values
        pairs = itertools.combinations(self.imap.keys(), 2)
        for (pop1, pop2) in pairs:
            pop1idx = self.imap[pop1]
            pop2idx = self.imap[pop2]
            popaidx = pop1idx + pop2idx

            within1 = list(itertools.combinations(pop1idx, 2))
            within2 = list(itertools.combinations(pop2idx, 2))
            withins = within1 + within2
            allpairs = itertools.combinations(popaidx, 2)
            betweens = itertools.filterfalse(
                lambda x: bool(x in withins),
                allpairs
            )
            # A list of T/F values for differences between each pair per site
            diff = [locus.loc[i] != locus.loc[j] for (i, j) in withins]
            # sum of all pairwise differences per site
            sums = np.sum(diff, axis=0)
            # average number of differences per site
            a = np.sum(sums) / sums.shape[0]

            # Same as above, but now pairwise for samples between sites
            diff = [locus.loc[i] != locus.loc[j] for (i, j) in betweens]
            sums = np.sum(diff, axis=0)
            b = np.sum(sums) / sums.shape[0]

            # Fst - Hudson 1992 Eq 3
            farr.loc[pop1, pop2] = 1 - (a / b)
            # Fst adjusted for known # of subpops - Hudson 1992 Eq 6
            darr.loc[pop1, pop2] = 1 - (a / (((1 / d) * a) \
                                        + (((d - 1) / d) * b)))
            # Nm adjusted for known # of subpops - Hudson 1992 Eq 7
            narr.loc[pop1, pop2] = (((d - 1) / d) * (1 / 2) * (a / (b - a)))

        return farr, darr, narr


    def _fis(self):
        "calculate population inbreeding Fis after filtering"
        pass


    # "Public" methods for the between populations stats
    def Dxy(self, locus, pops):
        """
        Calculate Dxy for a given locus.
        """
        pop_cts, sidxs = self._process_locus_pops(locus, pops)
        return self._dxy(*pop_cts.values(), len(locus))


    # utils
    
    def hill_number(abunds, order=0):
        """
        Get one hill humber from an array-like of values.
    
        :param array-like abunds: An `array-like` of counts of individuals per
            species.
        :param float order: The order of the Hill number to calculate.
    
        :return: The Hill number of order `order`.
        """
        abunds = np.array(abunds)
        ## Degenerate edge cases can cause all zero values, particulary for pi
        if not np.any(abunds):
            return 0
        if order == 0:
            return len(np.nonzero(abunds)[0])
        if order == 1:
            h1 = np.exp(entropy(abunds))
            return h1
        tot = float(np.sum(abunds))
        proportions = np.array(abunds[abunds > 0])/tot
        prop_order = proportions**order
        h2 = np.sum(prop_order)**(1./(1-order))
        return h2
    
    def hill_numbers(abunds, orders, granularity=None, do_negative=False):
        """
        Get all hill numbers from 0 to 'orders' from an array-like of abundances.
        If `granularity` is specified then calculate fractional Hill numbers
        such that the returned `numpy.array` contains `granularity` count of Hill
        numbers calculated and equally spaced between 0 and `orders`. The
        `do_negative` parameter will include both positive and negative values
        of Hill numbers up to order equal to `orders`.
    
        :param array-like abunds: An `array-like` of counts of individuals per
            species.
        :param int orders: The max order of Hill numbers to calculate.
        :param int granularity: The number of equally spaced fractional Hill
            numbers to calculate between 0 and `orders`. If not specified then
            returns only integer values of Hill numbers between 0 and `orders`.
        :param bool do_negative: Whether to calculate negative as well as positive
            Hill numbers between 0 and `orders`.
    
        :return numpy.array: An array of Hill numbers from 0 to (+/-) `orders`
            in intervals of defined by the `granularity` parameter.
        """
        ret = []
        min_order = 0
        if not granularity: granularity = orders + 1
        if do_negative:
            min_order = -orders
            granularity *= 2
        for order in np.linspace(min_order, orders, granularity):
            ret.append(hill_number(abunds, order))
        return np.array(ret)



