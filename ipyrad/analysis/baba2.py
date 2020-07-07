#!/usr/bin/env python

"D-statistic calculations"

# py2/3 compat
from __future__ import print_function, division
from builtins import range

import time
from collections import OrderedDict

# import scipy.stats as st  ## used for dfoil
from numba import njit
import pandas as pd
import numpy as np

# ipyrad tools
from .utils import ProgressBar
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError
from .snps_extracter import SNPsExtracter

# check for toytree and toyplot
try:
    import toytree
    import toyplot
    import scipy.stats as sc
except ImportError:
    pass


"""
TODO: 
- sliding window analysis (wrap around SNPsExtracter...)
"""


class Baba:
    """
    ipyrad.analysis Baba Class object.

    Parameters
    ----------
    data : str
        Path to the .snps.hdf5 input file.

    Functions
    ---------
    run_test: 

    run:

    generate_tests_from_tree:

    plot:

    """
    def __init__(
        self, 
        data=None,
        ):

        # store tests
        self.data = data

        # results storage
        self.results_table = None
        self.taxon_table = None

        # cluster attributes
        self.ipcluster = {
            "cluster_id": "", 
            "profile": "default",
            "engines": "Local", 
            "quiet": 0, 
            "timeout": 60, 
            "cores": 0, 
            "threads": 2,
            "pids": {},
            }



    def _get_pop_freqs(self, arr, imap):
        """
        Calculates frequencies of 'derived' alleles at each site in 
        each population group after choosing 'ancestral' allele 
        based on what is present in the outgroup.

        # nb: ma arrays do not support numba jit
        """
        barr = np.zeros((4, arr.shape[1]), dtype=np.float)

        # iterate over populations
        for pidx, pop in enumerate(("p1", "p2", "p3", "p4")):

            # get sample idxs
            sidx = [self.snex.names.index(i) for i in imap[pop]]

            # mask missing data
            marr = np.ma.array(data=arr[sidx, :], mask=arr[sidx, :] == 9)

            # proportion derived
            freq = (marr / 2).mean(axis=0).data
            barr[pidx] = freq

        # if p4 deriv freq is >50% then invert to make it 50% ancestor.
        flip = barr[-1] >= 0.5
        barr[:, flip] = 1 - barr[:, flip]
        return barr


    @staticmethod
    @njit
    def _get_dstat(barr):
        """
        Returns D-stat and return abba baba freq arrays
        """
        abba = (1 - barr[0]) * (barr[1]) * (barr[2]) * (1 - barr[3])
        baba = (barr[0]) * (1 - barr[1]) * (barr[2]) * (1 - barr[3])
        sabba = abba.sum()
        sbaba = baba.sum()
        dstat = (sabba - sbaba) / (sabba + sbaba)
        return dstat, sabba, sbaba



    def _get_boots(self, imap, nboots):
        """
        Returns array of bootstrap replicate D-stats
        """
        boots = np.zeros(nboots)
        for idx in range(nboots):
            arr = self.snex.subsample_loci(quiet=True)
            barr = self._get_pop_freqs(arr, imap)
            dhat, _, _ = self._get_dstat(barr)
            boots[idx] = dhat
        return boots



    def _get_test_results(self, imap, nboots):
        """
        Returns a row of results for a single test.
        """
        barr = self._get_pop_freqs(self.snex.snps, imap)
        dhat, abba, baba = self._get_dstat(barr)
        boots = self._get_boots(imap, nboots)
        zstat = abs(dhat) / boots.std()
        return dhat, boots.std(), zstat, abba, baba  # .sum(), baba.sum()



    def _get_pop_freqs_5(self, arr, imap):
        """
        5-taxon allels frequencies
        """
        barr = np.zeros((5, arr.shape[1]), dtype=float)
        for pidx, pop in enumerate(["p1", "p2", "p3_1", "p3_2", "p4"]):
            sidx = [self.snex.names.index(i) for i in imap[pop]]
            marr = np.ma.array(data=arr[sidx, :], mask=arr[sidx, :] == 9)
            freq = (marr / 2).mean(axis=0).data
            barr[pidx] = freq
        flip = barr[-1] >= 0.5
        barr[:, flip] = 1 - barr[:, flip]
        return barr



    @staticmethod
    @njit
    def _get_partitioned_dstats(barr):
        """
        Calculate partitioned D-statistics from 5 taxon tree.
        """     
        abbba = (1 - barr[0]) * (barr[1]) * (barr[2] * barr[3]) * (1 - barr[4])
        ababa = (1 - barr[0]) * (barr[1]) * ((1 - barr[2]) * barr[3]) * (1 - barr[4])
        abbaa = (1 - barr[0]) * (barr[1]) * (barr[2] * (1 - barr[3])) * (1 - barr[4])

        babba = (barr[0]) * (1 - barr[1]) * (barr[2] * barr[3]) * (1 - barr[4])
        baaba = (barr[0]) * (1 - barr[1]) * ((1 - barr[2]) * barr[3]) * (1 - barr[4])
        babaa = (barr[0]) * (1 - barr[1]) * (barr[2] * (1 - barr[3])) * (1 - barr[4])

        sabbba = abbba.sum()
        sbabba = babba.sum()
        sababa = ababa.sum()
        sbaaba = baaba.sum()
        sabbaa = abbaa.sum()
        sbabaa = babaa.sum()

        dstat12 = (sabbba - sbabba) / (sabbba + sbabba)
        dstat1 = (sabbaa - sbabaa) / (sabbaa + sbabaa)
        dstat2 = (sababa - sbaaba) / (sababa + sbaaba)
        return (
            dstat12, dstat1, dstat2, 
            sabbba, sbabba, 
            sabbaa, sbabaa,
            sababa, sbaaba, 
            )



    def _get_boots_5(self, imap, nboots):
        boots12 = np.zeros(nboots)
        boots1 = np.zeros(nboots)
        boots2 = np.zeros(nboots)
        for idx in range(nboots):
            arr = self.snex.subsample_loci(quiet=True)
            barr = self._get_pop_freqs_5(arr, imap)
            res = self._get_partitioned_dstats(barr)
            boots12[idx] = res[0]
            boots1[idx] = res[1]
            boots2[idx] = res[2]
        return boots12, boots1, boots2



    def _get_test_results_5(self, imap, nboots):
        """
        Returns a row of results for a single test
        """
        barr = self._get_pop_freqs_5(self.snex.snps, imap)
        res = self._get_partitioned_dstats(barr)
        boots = self._get_boots_5(imap, nboots)
        zstat12 = abs(res[0]) / boots[0].std()
        zstat1 = abs(res[1]) / boots[1].std()
        zstat2 = abs(res[2]) / boots[2].std()
        return res, boots, (zstat12, zstat1, zstat2)



    def run_partitioned_test(self, imap, minmap=None, nboots=100, quiet=False):
        """
        Returns Partitioned D-statistic results for a single test.
        """
        # parse genotypes for this subsapmle 
        self.snex = SNPsExtracter(
            self.data, imap=imap, minmap=minmap, mincov=4, quiet=quiet,
        )
        self.snex.parse_genos_from_hdf5()

        # get test results and arrange in dataframe
        res = self._get_test_results_5(imap, nboots)
        resdf = pd.DataFrame({
            "D12": res[0][0],
            "D1": res[0][1],
            "D2": res[0][2],
            "boot12std": res[1][0].std(),
            "boot1std": res[1][1].std(),
            "boot2std": res[1][2].std(),
            "Z12": res[2][0],
            "Z1": res[2][1],
            "Z2": res[2][2],                        
            "ABBBA": res[0][3],
            "BABBA": res[0][4],
            "ABBAA": res[0][5],
            "BABAA": res[0][6],
            "ABABA": res[0][7],
            "BAABA": res[0][8],
            "nSNPs": self.snex.snpsmap.shape[0],
            "nloci": len(set(self.snex.snpsmap[:, 0])),
        }, index=[0])
        return resdf



    def run_test(self, imap, minmap=None, nboots=100, quiet=False):
        """
        Return D-statistic results for a single test. 

        Parameters:
        -----------
        imap (dict):
            A test is defined using a dict mapping the keys p1, p2, p3, and p4
            to lists of sample names. When multiple samples represent a tip 
            the allele frequency is used. 

        minmap (dict):
            A dictionary containing the same keys as imap and with a float 
            or int as the value representing the minimum number (or proportion)
            of samples in the list that must have data at a SNP for it to be
            included in the analysis. This only applies if you have more than
            one sample per tip.

        nboots (int):
            The number of bootstrap samples used to calculate std dev. and 
            measure Z-score for significance testing.

        quiet (bool):
            Verbosity of SNP filtering.       

        Returns
        -------
        result (pandas.DataFrame):
            A DataFrame 
        """

        # check on imaps


        # check on minmaps


        # parse genotypes for this subsapmle 
        self.snex = SNPsExtracter(
            self.data, imap=imap, minmap=minmap, mincov=4, quiet=quiet,
        )
        self.snex.parse_genos_from_hdf5()

        # get test results and arrange in dataframe
        res = self._get_test_results(imap, nboots)
        res = pd.DataFrame({
            "D": res[0],
            "bootstd": res[1],
            "Z": res[2],
            "ABBA": res[3],
            "BABA": res[4],
            "nSNPs": self.snex.snpsmap.shape[0],
            "nloci": len(set(self.snex.snpsmap[:, 0])),
        }, index=[0])
        return res



    def _run(self, imaps, minmaps, nboots, ipyclient):

        # load-balancer
        lbview = ipyclient.load_balanced_view()

        # store the set of tests used here
        self.tests = imaps

        # expand minmaps
        if minmaps is None:
            minmaps = [
                {i: 1 for i in ["p1", "p2", "p3", "p4"]} 
                for j in range(len(imaps))
            ]

        # distribute job
        dfs = {}

        # distribute jobs
        rasyncs = {}
        idx = 0
        for imap, minmap in zip(imaps, minmaps):
            args = (self.data, imap, minmap, nboots, True)
            rasync = lbview.apply(remote_run, *args)
            rasyncs[idx] = rasync
            idx += 1

        # setup progress bar
        prog = ProgressBar(len(imaps), None, "abba-baba tests")
        prog.finished = 0
        prog.update()

        while 1:
            # check for completed
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                dfs[idx] = rasyncs[idx].get()
                prog.finished += 1
                del rasyncs[idx]

            # show progress
            prog.update()
            time.sleep(0.9)
            if not rasyncs:
                print("")
                break

        # concat results to df
        df = pd.concat([dfs[i] for i in range(len(imaps))], ignore_index=True)
        self.results_table = df

        # concat sample names to df strings
        self.taxon_table = pd.DataFrame(imaps).applymap(lambda x: ",".join(x))



    def run(self, imaps, minmaps=None, nboots=1000, auto=True, ipyclient=None, show_cluster=False):
        """
        Run a batch of dstat tests in parallel on a list of test dictionaries.
        The list of tests can either be set on the .tests attribute of the
        baba object, or auto-generated by calling .generate_tests_from_tree().

        Parameters:
        -----------
        auto (bool):
            Automatically start and stop parallel cluster using all cores
            available or with fine tuning by .ipcluster attribute params.

        force (bool):
            Overwrite existing results CSV file with this [workdir]/[name].csv

        ipyclient (ipyparallel.Client object):
            An ipyparallel client object to distribute jobs to a cluster. 
            This is an optional alternative to using 'auto=True'.

        show_cluster (bool):
            Verbose option to print information about n cores in cluster.
        """
        # distribute jobs in a wrapped cleaner function
        pool = Parallel(
            tool=self,
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={'imaps': imaps, 'minmaps': minmaps, 'nboots': nboots},
            )
        pool.wrap_run()

        # batch(self, ipyclient)
        # ## skip this for 5-part test results
        # if not isinstance(self.results_table, list):
        #     self.results_table.nloci = (
        #         np.nan_to_num(self.results_table.nloci).astype(int))



    def generate_tests_from_tree(
        self,
        tree,
        constraint_dict=None, 
        constraint_exact=False,
        return_idxs=False,
        quiet=False):
        """ 
        Returns a list of all possible 4-taxon tests on a tree (newick file). 
        The number of possible tests can be greatly reduced by setting 
        constraints on the taxon sampling using the constraint_dict arg. 

        Parameters:
        -----------
        constraint_dict: dict
            The constraint dict will limit the tests generated to only include
            the taxa listed in the dict. 

        constraint_exact: bool or list
            If constraint_exact is True then only samples meeting the exact 
            entries in the constraint_dict will be returned, as opposed to all
            subsets of those entries. If a list then different values can be 
            applied to [p1, p2, p3, p4]. For example, if the constraint_dict is
            {"p1": sample1, "p2": sample2, "p3": sample3, "p4": [sample4, sample5]},
            then with constraint_exact==False you get:

            sample1, sample2, sample3, sample4
            sample1, sample2, sample3, sample5
            sample1, sample2, sample3, [sample4, sample5]

            and with constraint_exact==True you get only:

            sample1, sample2, sample3, [sample4, sample5]
        """
        # init traversal extraction
        tests = TreeParser(tree, constraint_dict, constraint_exact).tests

        # print message success
        if not quiet:
            print("{} tests generated from tree".format(len(tests)))

        # convert tests to lists of names
        if not return_idxs:
            ntests = []
            for test in tests:
                tdict = {
                    "p1": tree.get_tip_labels(test[0]),
                    "p2": tree.get_tip_labels(test[1]),
                    "p3": tree.get_tip_labels(test[2]),
                    "p4": tree.get_tip_labels(test[3]),
                }
                ntests.append(tdict)
            tests = ntests
        else:
            tests = list(tests)

        # return the set of tests
        return tests



    def draw(self, tree, width=500, height=500, sort=False, prune=False, *args, **kwargs):
        """ 
        Draw a multi-panel figure with tree, tests, and results 

        Parameters:
        -----------
        width: int
            Width in pixels

        height: int
            Height in pixels

        prune: bool
            Prune the tree to only draw tips that are involved in tests.\
            
        sort: bool
            

        """

        # make the plot
        drawing = Drawing(self.results_table, self.taxon_table, tree, width, height, sort=sort, prune=prune)
        return (drawing.canvas, )




class Drawing:
    def __init__(self, res, tax, tree, width=500, height=500, sort=False, prune=False):
        
        self.tests = tax
        self.res = res
        self.ntests = res.shape[0]
        
      
        # if prune tree
        if prune:
            intree = set([])
            for cell in self.tests.values.flatten():
                for tax_ in cell.split(","):
                    intree.add(tax_)
            tree = tree.drop_tips(
                [i for i in tree.get_tip_labels() if i not in intree]
            )
        
        
        # define tree, original tree or prunned tree
        self.tree = tree
        
        
        if sort:
            # split to make cell into a list
            sindex = (
                self.tests
                .applymap(lambda x: x.split(","))
                .applymap(self.tree.get_mrca_idx_from_tip_labels)
                .sort_values(by=["p4", "p3", "p2", "p1"])
            ).index

            # rearrange tables by sindex
            self.tests = self.tests.loc[sindex]
            self.res = self.res.loc[sindex]
            self.tests.reset_index(drop=True, inplace=True)
            self.res.reset_index(drop=True, inplace=True)

        
        




        # canvas and axes components
        self.canvas = toyplot.Canvas(width, height)
        self.add_tree_to_canvas()
        self.add_zscores_to_canvas()
        self.add_histos_to_canvas()
        self.add_test_idxs_to_canvas()
        self.add_tip_names_to_canvas()
        self.add_tests_to_canvas()


    def add_tree_to_canvas(self):
        ax0 = self.canvas.cartesian(bounds=("50%", "90%", "5%", "19%"), show=False)
        self.tree.draw(
            axes=ax0, 
            ts='n', 
            layout='d', 
            tip_labels=False, 
            tip_labels_align=True, 
            xbaseline=0.5,
        )
        ax0.rectangle(
            0, self.tree.ntips, 
            0, self.tree.treenode.height, 
            style={"fill": "none"},
        )


    def add_test_idxs_to_canvas(self):
        # test names
        ax4 = self.canvas.cartesian(bounds=("91%", "95%", "21%", "80%"), show=False)
        ax4.rectangle(
            0, 1, 
            0, self.ntests + 1, 
            style={"fill": "none"})
        ax4.text(
            np.repeat(0, self.ntests),
            np.arange(self.ntests) + 1, 
            [str(i) for i in range(self.ntests)][::-1],
            style={"fill": "black", "text-anchor": "start"}
        )


    def add_tip_names_to_canvas(self):
        # tip names
        ax5 = self.canvas.cartesian(bounds=("50%", "90%", "80%", "97%"), show=False)
        ax5.rectangle(0, self.tree.ntips, 0, 1, style={"fill": "none"})
        ax5.text(
            np.arange(self.tree.ntips) + 0.5,
            np.repeat(0.9, self.tree.ntips),
            self.tree.get_tip_labels(),
            angle=-90,
            style={"fill": "black", "text-anchor": "start"},
            annotation=True,
        )


    def add_tests_to_canvas(self):
        # add tests bars to axes
        ax1 = self.canvas.cartesian(
            bounds=("50%", "90%", "21%", "80%"), 
            show=False,
            padding=0,
        )

        # spacer rect
        ax1.rectangle(
            0, self.tree.ntips, 
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1, 
            },
        )

        # coloring
        COLORS = toyplot.color.Palette()
        colors = [COLORS[0], COLORS[1], toyplot.color.black, COLORS[7]]
        TIPS = self.tree.get_tip_labels()

        # draw blocks
        for idx in range(self.ntests):

            # line tracing
            hidx = self.ntests - idx
            ax1.hlines(hidx, color=toyplot.color.black, style={"stroke-dasharray": "2,4"})

            # get test [name1, name2, name3]
            for cidx, pop in enumerate(["p1", "p2", "p3", "p4"]):
                test = self.tests.iloc[idx][pop]

                # get name indices [0, 2, 3]
                tidxs = sorted([TIPS.index(i) for i in test.split(",")])

                # draw blocks connecting index to next until no more.
                blocks = []
                block = [tidxs[0], tidxs[0]]
                for i in range(1, len(tidxs)):
                    if tidxs[i] - tidxs[i - 1] == 1:
                        block[-1] = tidxs[i]
                    else:
                        blocks.append(block)
                        block = [tidxs[i], tidxs[i]]
                blocks.append(block)
                blocks[-1][-1] = tidxs[-1]

                # draw them (left, right, top, bottom)
                for block in blocks:
                    ax1.rectangle(
                        a=block[0] + 0.25,
                        b=block[1] + 0.75,
                        c=hidx + 0.25, 
                        d=hidx - 0.25,
                        style={
                            "fill": colors[cidx],
                            "stroke": toyplot.color.black,
                            "stroke-width": 0.5,
                        },
                    )
        ax1.hlines(
            [0, self.ntests + 1], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5}
        )
        ax1.vlines(
            [0, self.tree.ntips], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},
        )        


    def add_zscores_to_canvas(self):
        # add zscores bars to axes
        ax2 = self.canvas.cartesian(
            bounds=("25%", "47%", "21%", "80%"), 
            yshow=False,
            padding=0,
        )

        # the longest bar space
        maxz = max(self.res.Z) + (max(self.res.Z) * .10)

        # spacer rect
        ax2.rectangle(
            -maxz, 0,
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1,
            },
        )

        # add data bars
        for idx in range(self.ntests):
            hidx = self.ntests - idx
            ax2.hlines(hidx, color='black', style={"stroke-dasharray": "2,4"})
            ax2.rectangle(
                0, -self.res.Z[idx],
                hidx - 0.25, hidx + 0.25, 
                color=toyplot.color.black,
            )


        # stylring 
        ax2.x.spine.show = False
        ax2.x.label.text = "Z-score"
        ax2.x.ticks.locator = toyplot.locator.Extended(5, only_inside=True)
        ax2.vlines(
            [ax2.x.domain.min, ax2.x.domain.max, 0, -maxz], 
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},
        )
        ax2.hlines(
            [0, self.ntests + 1],
            style={"stroke": toyplot.color.black, "stroke-width": 1.5},            
        )


    def add_histos_to_canvas(self):
        # add histograms to axes
        ax3 = self.canvas.cartesian(
            bounds=("5%", "22%", "21%", "80%"), 
            yshow=False, 
            padding=0,
        )

        zmin = min(self.res.D - 3.25 * self.res.bootstd[0])
        zmax = max(self.res.D + 3.25 * self.res.bootstd[0])

        # draw outline and fill
        ax3.rectangle(
            zmin, zmax,
            0, self.ntests + 1, 
            style={
                "fill": "grey", 
                "fill-opacity": 0.1, 
            },
        )

        # iterate over tests to add histos
        for idx in range(self.ntests):
            hidx = self.ntests - idx

            # get fill color
            if self.res.Z[idx] < 2.5:
                fill = toyplot.color.Palette()[7]
            else:
                if self.res.D[idx] > 0:
                    fill = toyplot.color.Palette()[1]
                else:
                    fill = toyplot.color.Palette()[0]

            # histogram fill
            points = np.linspace(zmin, zmax, 30)
            density = sc.norm.pdf(
                points, loc=self.res.D[idx], scale=self.res.bootstd[idx],
            )
            ax3.fill(
                points, density / density.max() * 0.7,
                baseline=np.repeat(hidx - 0.25, len(points)),
                style={
                    "stroke": 'black', 
                    "stroke-width": 0.5, 
                    "fill": fill},
            )

        # Z=0 indicator    
        ax3.vlines(
            0, 
            style={
                "stroke": "grey", 
                "stroke-dasharray": "2,4", 
                "stroke-width": 1,
            })

        ax3.vlines(
            [zmin, zmax],
            style={"stroke": "black", "stroke-width": 1.5},
        )
        ax3.hlines(
            [0, self.ntests + 1],
            style={"stroke": "black", "stroke-width": 1.5},
        )        

        # style axes
        ax3.x.label.text = "D-statistic"
        ax3.x.spine.show = False
        ax3.x.ticks.locator = toyplot.locator.Explicit(
            [zmin, 0.0, zmax],
            ["{:.1f}".format(i) for i in [zmin, 0.0, zmax]],
        )



class TreeParser:
    def __init__(self, tree, constraint_dict, constraint_exact):
        "Traverses tree to build test sets given constraint options."

        # store sets of four-taxon splits
        self.testset = set()
        self.hold = [0, 0, 0, 0]

        # tree to traverse
        self.tree = toytree.tree(tree)
        if not self.tree.is_rooted(): 
            raise IPyradError(
                "generate_tests_from_tree(): tree must be rooted and resolved")

        # store contraints
        self.cdict = OrderedDict((i, []) for i in ["p1", "p2", "p3", "p4"])
        # self.cdict = [(0, 0, 0, 0) for i in ]

        # constraints entered as a dict or tuple: (0, 1, 10, 13)
        if isinstance(constraint_dict, dict):
            for key, val in constraint_dict.items():
                if isinstance(val, int):
                    val = tree.get_tip_labels(val)
                self.cdict[key] = val

        elif isinstance(constraint_dict, (list, tuple, np.ndarray)):
            for cidx, pop in enumerate(["p1", "p2", "p3", "p4"]):
                const = constraint_dict[cidx]
                if isinstance(const, int):
                    self.cdict[pop] = (
                        tree.get_tip_labels(const)
                    )

        # constraint setting [True, True, False, False]
        self.xdict = constraint_exact
        if isinstance(self.xdict, bool):
            self.xdict = [self.xdict] * 4
        if isinstance(self.xdict, (tuple, list, np.ndarray)):
            if len(self.xdict) != len(self.cdict):
                print(self.xdict, self.cdict)
                raise Exception(
                    "constraint_exact must be bool or list of bools length N")
        self.xdict = np.array(self.xdict).astype(bool)

        # get tests
        self.loop(self.tree.treenode)

        # order and check redundancy
        tests = []
        coords = tree.get_node_coordinates(layout='d')
        for test in self.testset:
            stest = sorted(test[:2], key=lambda x: coords[x, 0])
            ntest = stest[0], stest[1], test[2], test[3]
            if ntest not in tests:
                tests.append(ntest)
        self.tests = tests


    def loop(self, node):  # , idx):
        "getting closer...."
        for topnode in node.traverse():
            for oparent in topnode.children:
                for onode in oparent.traverse():
                    if self.test_constraint(onode, 3):                       
                        self.hold[3] = onode.idx

                        node2 = oparent.get_sisters()[0]
                        for topnode2 in node2.traverse():
                            for oparent2 in topnode2.children:
                                for onode2 in oparent2.traverse():
                                    if self.test_constraint(onode2, 2):                       
                                        self.hold[2] = onode2.idx

                                        node3 = oparent2.get_sisters()[0]
                                        for topnode3 in node3.traverse():
                                            for oparent3 in topnode3.children:
                                                for onode3 in oparent3.traverse():
                                                    if self.test_constraint(onode3, 1):                       
                                                        self.hold[1] = onode3.idx

                                                        node4 = oparent3.get_sisters()[0]
                                                        for topnode4 in node4.traverse():
                                                            for onode4 in topnode4.traverse():
                                                                if self.test_constraint(onode4, 0):
                                                                    self.hold[0] = onode4.idx
                                                                    self.testset.add(tuple(self.hold))
                                                            # for oparent4 in topnode4.children:
                                                            #     for onode4 in oparent4.traverse():
                                                            #         if self.test_constraint(onode4, 0):
                                                            #             self.hold[0] = onode4.idx
                                                            #             self.testset.add(tuple(self.hold))


    def test_constraint(self, node, idx):
        names = set(node.get_leaf_names())
        const = set(list(self.cdict.values())[idx])
        if const:
            if self.xdict[idx]:
                if names == const:
                    return 1
                else:
                    return 0
            else:
                if len(names.intersection(const)) == len(names):
                    return 1
                else:
                    return 0        
        return 1        



def remote_run(data, imap, minmap, nboots, quiet):
    "to be called on ipengine"
    self = Baba(data)
    res = self.run_test(imap, minmap, nboots, quiet)
    return res


#######################################################################

if __name__ == "__main__":

    ## test input files
    LOCIFILE = "/home/deren/Dropbox/RADexplore/EmpVib/"\
              + "vib_half_64tip_c85d6m4p99.loci"

    # ## taxon list to parse from LOCIFILE
    TAXONLIST = ['acutifolium_DRY3_MEX_006',
                 'sulcatum_D9_MEX_003',
                 'jamesonii_D12_PWS_1636',
                 'triphyllum_D13_PWS_1783',
                 'dentatum_ELS4']

    ## calculate dstats