#!/usr/bin/env python

"""Calculate D-statistics from RAD loci using bootstrap resampling.


"""

from typing import Dict, List, Union, Optional

from loguru import logger
import pandas as pd
import numpy as np
import toytree
# from numba import njit
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.core.parallel import Cluster
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.snps_extracter import SNPsExtracter
from ipyrad.analysis.baba_drawing import Drawing

logger = logger.bind(name="ipa")


class Baba:
    """...

    Parameters
    ----------
    data: str
        A file path to an input .snps.hdf5 data file.
    """
    def __init__(self, data: str):
        self.data = data

        # attrs to be filled.
        self.results_table: pd.DataFrame=None
        """A pandas DataFrame with results from the last test ran."""
        self.taxon_table: pd.DataFrame=None
        """A pandas DataFrame with sample names from the last test run."""
        self._ext: 'ipa.snps_extracter'=None
        """An ipa.snps_extracter tool used to filter snps."""

    def run_test(
        self,
        imap: Dict[str,List[str]],
        minmap: Union[Dict[str, Union[int, float]], int, float],
        nboots: int=100,
        ):
        """Return a DataFrame results for a single D-statistic test.

        Parameters
        ----------
        imap: dict
            A test is defined using a dict mapping the keys p1, p2, p3,
            and p4 to lists of sample names. When multiple samples
            represent a tip the allele frequency is used.
        minmap: dict, int, or float.
            A dict mapping the same keys as in imap to int or float
            values representing the minimum number (or proportion)
            of samples in the list that must have data at a SNP for
            it to be included in the analysis. This only applies if
            you have more than one sample per tip, since the minmap
            will default to a minimum of 1 per group.
        nboots: int
            The number of bootstrap samples used to calculate std dev
            and measure Z-score for significance testing.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the number of observed site patterns and
            the test significance based on bootstrap resampling loci
            with replacement.
        """

        # check on imaps
        assert sorted(imap.keys()) == ["p1", "p2", "p3", "p4"], (
            "Malformed imap: must contain keys 'p1', 'p2', 'p3', and 'p4'.\n"
            f"You entered: {imap}"
        )

        # check on minmaps
        if isinstance(minmap, (int, float)):
            minmap = {i: minmap for i in imap}

        # parse genotypes for this subsapmle
        self._ext = SNPsExtracter(
            self.data,
            imap=imap,
            minmap=minmap,
            mincov=4,
        )
        self._ext.run(cores=1)

        # get test results and arrange in dataframe
        return self._get_test_results(imap, nboots)

    def run_partitioned_test(
        self,
        imap: Dict[str,List[str]],
        minmap: Dict[str,Union[int,float]]=None,
        nboots: int=100,
        ):
        """Return partitioned D-statistic results for a single test.

        Parameters
        ----------
        imap: dict
            ...
        minmap: dict
            ...
        nboots: int
            ...
        """
        # parse genotypes for this subsapmle
        self._ext = SNPsExtracter(
            self.data,
            imap=imap,
            minmap=minmap,
            mincov=4,
        )
        self._ext.run()
        return self._get_test_results_5(imap, nboots)

    def _get_test_results(self, imap, nboots):
        """Return a pd.Series of results for a single test."""
        barr = self._get_pop_freqs(self._ext.genos, imap)
        dhat, abba, baba = self._get_dstat(barr)
        boots = self._get_boots(imap, nboots)
        boots_std = boots.std()
        zstat = abs(dhat) / boots_std
        return pd.Series(
            name="D-statistic",
            dtype=float,
            data={
                "D": dhat,
                "D_bootstrap_std": boots_std,
                "Z": zstat,
                "ABBA": abba,
                "BABA": baba,
                "nSNPs": self._ext.snpsmap.shape[0],
                "nloci": len(set(self._ext.snpsmap[:, 0])),
            },
        )

    def _get_test_results_5(self, imap, nboots):
        """Returns a pd.Series of results for a single test."""
        barr = self._get_pop_freqs_5(self._ext.genos, imap)
        res = self._get_partitioned_dstats(barr)
        boots = self._get_boots_5(imap, nboots)
        zstat12 = abs(res[0]) / boots[0].std()
        zstat1 = abs(res[1]) / boots[1].std()
        zstat2 = abs(res[2]) / boots[2].std()
        return pd.Series(
            name="partitioned-D-statistic",
            dtype=float,
            data={
                "D12": res[0][0],
                "D1": res[0][1],
                "D2": res[0][2],
                "boot12std": res[1][0].std(),
                "boot1std": res[1][1].std(),
                "boot2std": res[1][2].std(),
                "Z12": zstat12,
                "Z1": zstat1,
                "Z2": zstat2,
                "ABBBA": res[0][3],
                "BABBA": res[0][4],
                "ABBAA": res[0][5],
                "BABAA": res[0][6],
                "ABABA": res[0][7],
                "BAABA": res[0][8],
                "nSNPs": self._ext.snpsmap.shape[0],
                "nloci": len(set(self._ext.snpsmap[:, 0])),
            },
        )

    def _get_pop_freqs(self, arr, imap):
        """Calculates frequencies of 'derived' alleles for each site x pop.

        This first chooses an 'ancestral' allele for each site based
        on the most frequent allele (0/1) in the outgroup (p4). All
        sites are already post-filtered, and thus have at least one
        allele per population (minmap minimum enforced).

        Note
        ----
        ma arrays do not support numba jit
        """
        barr = np.zeros((4, arr.shape[1]), dtype=np.float)

        # iterate over populations
        for pidx, pop in enumerate(("p1", "p2", "p3", "p4")):

            # get sample idxs
            sidx = [self._ext.names.index(i) for i in imap[pop]]

            # mask missing data
            marr = np.ma.array(data=arr[sidx, :], mask=arr[sidx, :] == 9)

            # proportion derived
            freq = (marr / 2).mean(axis=0).data
            barr[pidx] = freq

        # if p4 deriv freq is >50% then invert to make it 50% ancestor.
        flip = barr[-1] >= 0.5
        barr[:, flip] = 1 - barr[:, flip]
        return barr

    def _get_pop_freqs_5(self, arr, imap):
        """Calculates frequencies of 'derived' alleles for each site x pop.

        This first chooses an 'ancestral' allele for each site based
        on the most frequent allele (0/1) in the outgroup (p4). All
        sites are already post-filtered, and thus have at least one
        allele per population (minmap minimum enforced).

        Note
        ----
        ma arrays do not support numba jit
        """
        barr = np.zeros((5, arr.shape[1]), dtype=float)
        for pidx, pop in enumerate(["p1", "p2", "p3_1", "p3_2", "p4"]):
            sidx = [self._ext.names.index(i) for i in imap[pop]]
            marr = np.ma.array(data=arr[sidx, :], mask=arr[sidx, :] == 9)
            freq = (marr / 2).mean(axis=0).data
            barr[pidx] = freq
        flip = barr[-1] >= 0.5
        barr[:, flip] = 1 - barr[:, flip]
        return barr

    # @njit
    @staticmethod
    def _get_dstat(barr) -> (float, int, int):
        """Return D-stat and abba baba counts."""
        abba = (1 - barr[0]) * (barr[1]) * (barr[2]) * (1 - barr[3])
        baba = (barr[0]) * (1 - barr[1]) * (barr[2]) * (1 - barr[3])
        sabba = abba.sum()
        sbaba = baba.sum()
        dstat = (sabba - sbaba) / (sabba + sbaba)
        return dstat, sabba, sbaba

    @staticmethod
    def _get_partitioned_dstats(barr):
        """Return partitioned D-statistics and site counts."""
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

    def _get_boots(self, imap, nboots):
        """Return an array of bootstrap replicate D-stats."""
        boots = np.zeros(nboots)
        for idx in range(nboots):
            arr = self._ext.subsample_loci(log_level="INFO" if not idx else "DEBUG")
            barr = self._get_pop_freqs(arr, imap)
            dhat, _, _ = self._get_dstat(barr)
            boots[idx] = dhat
        return boots

    def _get_boots_5(self, imap, nboots):
        """Return an array of bootstrap replicate partitioned D-stats."""
        boots12 = np.zeros(nboots)
        boots1 = np.zeros(nboots)
        boots2 = np.zeros(nboots)
        for idx in range(nboots):
            arr = self._ext.subsample_loci(log_level="INFO" if not idx else "DEBUG")
            barr = self._get_pop_freqs_5(arr, imap)
            res = self._get_partitioned_dstats(barr)
            boots12[idx] = res[0]
            boots1[idx] = res[1]
            boots2[idx] = res[2]
        return boots12, boots1, boots2

    def run(
        self,
        imaps: List[Dict[str,List[str]]],
        minmaps: List[Union[Dict,float,int]],
        nboots: int=100,
        cores: int=4,
        ipyclient: Optional["ipyparallel.Client"]=None,
        ):
        """Run a batch of dstat tests in parallel.

        The imaps args takes as input a list of test dictionaries.
        The list of tests can either be set on the .tests attribute
        of the baba object, or auto-generated by calling
        .generate_tests_from_tree().

        Parameters
        ----------
        imaps: list of imap dictionaries
            This ...
        minmaps: 
            ...
        nboots: int
            Number of bootstrap replicates to run.
        cores: int
            The number of cores to parallelize the jobs on.
        ipyclient: ipyparallel.Client object
            An ipyparallel client object to distribute jobs to a
            cluster. This is an optional alternative to using an
            automatic cluster, which is done if left as None.
        """
        # distribute jobs in a wrapped cleaner function
        cluster = Cluster()
        cluster.start(cores=cores)
        lbview = ipyclient = cluster.ipyclient.load_balanced_view()

        # progress bar tracker
        prog = AssemblyProgressBar({}, "abba-baba", "ipa", True)

        # check and expand minmaps if None
        if minmaps is None:
            minmaps = [
                {i: 1 for i in ["p1", "p2", "p3", "p4"]}
                for j in range(len(imaps))
            ]
        # self.run_test(imaps[0], minmaps[0], nboots)

        # distribute jobs
        for idx, _ in enumerate(imaps):
            args = (imaps[idx], minmaps[idx], nboots)
            rasync = lbview.apply(self.run_test, *args)
            prog.jobs[idx] = rasync
        prog.update()
        prog.block()
        prog.check()        

        # concat results to df
        data = pd.concat(
            [prog.jobs[idx].get() for idx in sorted(prog.jobs)],
            axis=1,
            ignore_index=True,
        )
        self.results_table = data.T
        logger.info(f"{data.shape[0]} test results stored to `.results_table`")

        # concat sample names to df strings
        self.taxon_table = pd.DataFrame(imaps).applymap(lambda x: ",".join(x))
        cluster.cleanup_safely(None)


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



    def draw(
        self,
        tree,
        width=500,
        height=500,
        sort=False,
        prune=False,
        fade=False,
        zscoreTH=2.5,
        **kwargs,
        ):
        """Draw a multi-panel figure with tree, tests, and results

        Parameters
        ----------
        width: int
            Width in pixels
        height: int
            Height in pixels
        prune: bool
            Prune the tree to only draw tips that are involved in tests.
        sort: bool
            Sort tests

        fade: float
            Fade test blocks if the Z-score is not significant.
        """
        # make the plot
        drawing = Drawing(
            self.results_table,
            self.taxon_table,
            tree,
            width,
            height,
            sort=sort,
            prune=prune,
            fade=fade,
            zscoreTH=zscoreTH,
        )
        return drawing.canvas






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