#!/usr/bin/env python

"Sliding window (or sampling window) for phylo inference"

# py2/3 compat
from __future__ import print_function

# standard
import os
import sys
import time
import shutil
import tempfile

# third party
import pandas as pd
import numpy as np

# internal librries
from .raxml import Raxml as raxml
from .mrbayes import MrBayes as mrbayes
from .window_extracter import WindowExtracter as window_extracter
from .utils import ProgressBar
from ..core.Parallel import Parallel

# suppress the terrible h5 warning
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py

try:
    import toytree
except ImportError:
    pass

_MISSING_TOYTREE = """
This ipyrad.analysis tool requires the dependency 'toytree'. 
You can install it with the following command from a terminal:

conda install toytree -c eaton-lab
"""

"""
Infer whole scaffold if windowsize = 0, None
scaffold_idx = 0 is default.
extract PHY and trim for a given window entered...
"""


class TreeSlider(object):
    """
    Performs phylo inference across RAD data sampled in windows. Uses the
    hdf5 database output from ipyrad as input (".seqs.hdf5"). If no window 
    size is entered then entire scaffolds are used as windows. 

    Parameters:
    ------------
    name: name prefix for output files.
    workdir: directory for output files.
    data: .loci file
    imap: optional dictionary mapping of sample names to new names. 
    minmap: optional dictionary of minimum sampling per imap group.
    minsnps: minimum number of SNPs to include window in analysis.
    inference_method: 'raxml' or 'mb'
    """
    def __init__(
        self,
        data,
        name=None,
        workdir="./analysis-treeslider",
        window_size=None, 
        slide_size=None,
        scaffold_idxs=None,
        minsnps=1,
        mincov=0,
        imap={},
        inference_method="raxml",
        inference_args={},
        quiet=False,
        ):

        # check installations
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # store attributes
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))
        self.tmpdir = os.path.join(self.workdir, "tmpdir")

        # work
        self.scaffold_idxs = scaffold_idxs
        self.window_size = window_size
        self.slide_size = slide_size
        self.minsnps = minsnps
        self.imap = imap
        self.mincov = mincov
        self.inference_method = inference_method
        self.inference_args = inference_args
        self.quiet = quiet

        # use user name else create one
        if not self.name:
            self.name = "test"               

        # get outfile name
        self.tree_table_path = os.path.join(
            self.workdir,
            "{}.tree_table.csv".format(self.name))

        # parallelization
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

        # to-be parsed attributes
        self.tree_table = None
        self.scaffold_table = None
        self.phymap = None
        self._pnames = None

        # checks params and loads tree table if existing.
        self._parameter_check()

        # get scaffold names and lengths
        self._init_scaffold_table()

        # default to all scaffolds if none entered.
        if self.scaffold_idxs is None:
            self.scaffold_idxs = self.scaffold_table.index.tolist()
        elif isinstance(self.scaffold_idxs, (list, tuple, set)):
            self.scaffold_idxs = sorted(self.scaffold_idxs)
        elif isinstance(self.scaffold_idxs, int):
            self.scaffold_idxs = [self.scaffold_idxs]

        # build the tree table from the scaffolds, windows, and slides.
        if self.tree_table is None:
            self._init_tree_table()


    def _print(self, message):
        if not self.quiet:
            print(message)


    def show_inference_command(self, show_full=False):
        # show a warning if user entered threads
        # if "T" in self.inference_args:
        # print("\n".join([
        # "threads argument (T) must be set in ipcluster. ",
        # "e.g., tet.ipcluster['threads'] = 4"
        # ]))

        # debug inference args
        threads = {"T": max(1, self.ipcluster["threads"])}
        self.inference_args.update(threads)
        rax = raxml(
            data=self.data, 
            name=self.name,
            workdir=tempfile.gettempdir(),
            **self.inference_args
        )

        # return it
        if show_full:
            print(rax.command)

        # pretty print it
        else:
            printkwargs = {
                "s": "...", 
                "w": "...", 
                "n": "...",        
            }
            rax.params.update(printkwargs)
            print(rax.command)


    def _parameter_check(self):
        assert os.path.exists(self.data), "database file not found"
        assert self.data.endswith(".seqs.hdf5"), (
            "data must be '.seqs.hdf5' file.")

        # if not window then slide is set to window
        if (not self.window_size) or (not self.slide_size):
            self.slide_size = self.window_size

        # ensure workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)


    def _load_existing_tree_table(self):

        # if CSV exists then load it (user can overwrite with run(force))
        if os.path.exists(self.tree_table_path):

            # load the existing table
            self.tree_table = pd.read_csv(
                self.tree_table_path, sep=",", index_col=0)

            # is table finished or incomplete?
            completed = (0 if np.any(self.tree_table.tree == 0) else 1)               
            if not completed:
                # report message
                msg = "\n".join([
                    "Unfinished tree_table loaded from [workdir]/[name].",
                    "You can continue filling the table from this checkpoint",
                    "by calling .run without using the force flag."
                    "Path: {}"
                ])
            else:
                msg = "\n".join([
                    "Finished tree_table loaded from [workdir]/[name].",
                    "Call run with force=True to overwrite these results,", 
                    "or set a new name or workdir to use a new file path.",
                    "Path: {}"
                ])
            self._print(msg.format(self.tree_table_path))


    def _init_scaffold_table(self):
        "get chromosome lengths from the database"
        with h5py.File(self.data) as io5:

            # parse formatting from db
            self._pnames = np.array([
                i.decode() for i in io5["phymap"].attrs["phynames"]
            ])
            self._longname = 1 + max([len(i) for i in self._pnames])

            # parse names and lengths from db
            scafnames = [i.decode() for i in io5["scaffold_names"][:]]
            scaflens = io5["scaffold_lengths"][:]
            self.scaffold_table = pd.DataFrame(
                data={
                    "scaffold_name": scafnames,
                    "scaffold_length": scaflens,
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )


    def _init_tree_table(self):
        "Build DataFrame for storing results"
        dfs = []
        for scaffold in self.scaffold_idxs:

            # get the length of this scaffold
            chromlen = (
                self.scaffold_table.loc[scaffold, "scaffold_length"])
            
            # get start positions
            if self.window_size:
                starts = np.arange(0, chromlen - self.window_size, self.slide_size)
            else:
                starts = [0]

            # get end positions
            if self.window_size:
                ends = np.arange(self.window_size, chromlen, self.slide_size)
            else:
                ends = [chromlen]

            # build to table
            df = pd.DataFrame(data={
                "scaffold": scaffold,
                "start": starts,
                "end": ends,
                "snps": 0,  # np.nan,
                "sites": 0,
                "samples": 0,
                "missing": 0.0,  # np.nan,
                "tree": 0,
                }, 
                columns=[
                    "scaffold", "start", "end", 
                    "sites", "snps", "samples", "missing", "tree"
                ], 
            )
            dfs.append(df)

        # concat data from one or more scaffolds
        self.tree_table = pd.concat(dfs)
        self.tree_table.reset_index(drop=True, inplace=True)


    def _parse_scaffold_phymap(self, scaffold_idx):
        "scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table"
        with h5py.File(self.data) as io5:
            colnames = io5["phymap"].attrs["columns"]

            # mask to select this scaff
            mask = io5["phymap"][:, 0] == self.scaffold_idx + 1

            # load dataframe of this scaffold
            self.phymap = pd.DataFrame(
                data=io5["phymap"][mask, :],
                columns=[i.decode() for i in colnames],
            )
     

    def run(self, ipyclient=None, force=False, show_cluster=False, auto=False):
        """
        Distribute tree slider jobs in parallel. 

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

        # TODO: This is a little complicated b/c we need to store and load old
        # params also...
        # check for results
        # if not force:
        # self._load_existing_tree_table()

        # wrap analysis in parallel client
        pool = Parallel(
            tool=self, 
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            rkwargs={"force": force},
            )
        pool.wrap_run()


    def _run(self, force=False, ipyclient=None):
        """
        Hidden func to distribute jobs that is wrapped inside Parallel.
        """
        # do not overwrite tree table
        if os.path.exists(self.tree_table_path):
            if not force:
                print((
        "Existing tree table loaded from {}; Use force to instead overwrite."
        .format(self.tree_table_path)))
                return

        # THREADING set to match between ipcluster and raxml 
        if "T" in self.inference_args:
            self.ipcluster["threads"] = max(2, int(self.inference_args["T"]))
        self.inference_args["T"] = max(2, int(self.ipcluster["threads"]))
        threads = self.inference_args["T"]

        # load balance parallel jobs 2-threaded
        lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::threads])

        # initial progress ticker to run during job submission
        self._print(
            "building database: nwindows={}; minsnps={}"
            .format(
                self.tree_table.shape[0],
                self.minsnps,
            ))

        # submit jobs: (fname, scafidx, minpos, maxpos, minsnps, )
        jidx = 0
        finished = 0
        rasyncs = {}
        for idx in self.tree_table.index:

            # if continuing an existing job, skip if row already filled
            if self.tree_table.tree[idx] != 0:
                finished += 1
                continue

            # extract the phylip alignment for this window
            ext = window_extracter(
                data=self.data,
                workdir=os.path.join(self.workdir, "tmpdir"),
                scaffold_idx=self.tree_table.scaffold[idx],
                start=self.tree_table.start[idx],
                end=self.tree_table.end[idx],
                mincov=self.mincov,
                imap=self.imap,
                quiet=True,
            )

            # fill table stats
            self.tree_table.loc[idx, "snps"] = ext.stats.loc["postfilter", "snps"]
            self.tree_table.loc[idx, "sites"] = ext.stats.loc["postfilter", "sites"]
            self.tree_table.loc[idx, "missing"] = ext.stats.loc["postfilter", "missing"]
            self.tree_table.loc[idx, "samples"] = ext.stats.loc["postfilter", "samples"]

            # filter by SNPs
            if ext.stats.loc["postfilter", "snps"] < self.minsnps:
                self.tree_table.loc[idx, "tree"] = np.nan
                finished += 1

            else:
                # write phylip file to the tmpdir
                ext.run(force=True)

                # send remote tree inference job
                args = [ext.outfile, self.inference_args]
                rasyncs[jidx] = lbview.apply(remote_raxml, *args)
                jidx += 1

        # wait for jobs to finish
        prog = ProgressBar(self.tree_table.shape[0], None, "inferring trees")
        prog.finished = finished

        while 1:
            # check for completed
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                self.tree_table.iloc[idx, -1] = rasyncs[idx].get()
                self.tree_table.to_csv(self.tree_table_path)
                prog.finished += 1
                del rasyncs[idx]

            # show progress
            prog.update()
            time.sleep(0.9)
            if not rasyncs:
                self._print("")
                break
        shutil.rmtree(self.tmpdir)



def remote_raxml(phyfile, inference_args):
    """
    Call raxml on phy and returned parse tree result
    """
    # call raxml on the input phylip file with inference args
    rax = raxml(
        data=phyfile,
        name="temp_" + str(os.getpid()),
        workdir=tempfile.gettempdir(),
        **inference_args
    )
    rax.run(force=True, quiet=True, block=True)

    # get newick string from result
    if os.path.exists(rax.trees.bipartitions):
        tree = toytree.tree(rax.trees.bipartitions).newick
    else:
        tree = toytree.tree(rax.trees.bestTree).newick

    # remote tree files
    for tfile in rax.trees:
        tpath = getattr(rax.trees, tfile)
        if os.path.exists(tpath):
            os.remove(tpath)

    # return results
    return tree


# DEPRECATED 
def remote_tree_inference(
    database, 
    wmin, 
    wmax, 
    minsnps, 
    threads=2, 
    mincov=4,
    inference_method="raxml",
    inference_args={},
    ):
    "remote job to run phylogenetic inference."

    # get seqarr for phy indices    
    nsnps = np.nan
    tree = np.nan
    nsamplecov = np.nan

    # parse phylip string if there is sequence in this window
    if wmax > wmin:

        # PARSE hdf5 window into a dictionary {name+spacer: seq} and info
        nsnps, nsamplecov, phydict = parse_phydict(database, wmin, wmax)

        # infer a tree if there is variation in this window
        if (nsamplecov >= mincov) and (nsnps >= minsnps):

            # RAxML ML inference
            if inference_method == "raxml":

                # write a temporary phylip file
                fname = write_phydict_to_phy(phydict)

                # init raxml object and run with blocking
                initkwargs = {
                    "T": max(1, threads),
                }
                initkwargs.update(inference_args)
                rax = raxml(
                    data=fname,
                    name="temp_" + str(os.getpid()),
                    workdir=tempfile.gettempdir(),
                    **initkwargs
                )

                # run command
                rax.run(force=True, quiet=True, block=True)

                # get tree file result
                if os.path.exists(rax.trees.bipartitions):
                    tree = toytree.tree(rax.trees.bipartitions).newick
                else:
                    tree = toytree.tree(rax.trees.bestTree).newick

            else:
                # write a temporary NEXUS file
                fname = write_phydict_to_nex(phydict)

                # init raxml object and run with blocking
                mb = mrbayes(
                    data=fname,
                    name="temp_" + str(os.getpid()),
                    workdir=tempfile.gettempdir(),
                    clock_model={
                        "clockratepr": "lognorm(-7,0.6)",
                        "clockvarpr": "tk02",
                        "tk02varpr": "exp(1.0)",
                        "brlenspr": "clock:birthdeath",
                        "samplestrat": "diversity",
                        "sampleprob": "0.1",
                        "speciationpr": "exp(10)",
                        "extinctionpr": "beta(2, 200)",
                        "treeagepr": "offsetexp(1, 5)",
                        "ngen": "1000000",
                        "nruns": "2",
                        "nchains": "4",
                        "samplefreq": "1000",
                    },                   
                )
                mb.run(force=True, quiet=True, block=True)

                # get tree file result
                tree = toytree.tree(mb.trees.constre)

    return nsnps, nsamplecov, tree



def write_phydict_to_phy(phydict):

    # dress up as phylip format
    phylip = ["{} {}".format(len(phydict), len(next(iter(phydict.values()))))]
    for name, seq in phydict.items():
        phylip.append("{} {}".format(name, seq))
    phystring = ("\n".join(phylip))

    # write to temp file
    fname = os.path.join(
        tempfile.gettempdir(), str(os.getpid()) + ".tmp")
    with open(fname, 'w') as temp:
        temp.write(phystring)
    return fname


def write_phydict_to_nex(phydict):

    # dress up as phylip format
    ntax = len(phydict)
    nchar = len(next(iter(phydict.values())))
    matrix = ""
    for i in phydict.items():
        matrix += i[0] + i[1] + "\n"
    
    # format into NEXUS data string
    nex_string = NEX_MATRIX.format(
        **{"ntax": ntax, "nchar": nchar, "matrix": matrix})
    
    # write to temp file
    fname = os.path.join(
        tempfile.gettempdir(), str(os.getpid()) + ".tmp")
    with open(fname, 'w') as temp:
        temp.write(nex_string)
    return fname


def parse_phydict(database, wmin, wmax):
    "Returns phy data as a dict with stats on sub matrix."

    with h5py.File(database, 'r') as io5:

        # select all sequence data in window range
        seqarr = io5["phy"][:, wmin:wmax]  # cmap[:, 1].min():cmap[:, 2].max()]
        pnames = np.array([
            i.decode() for i in io5["phymap"].attrs["phynames"]
        ])
        longname = 1 + max([len(i) for i in pnames])

    # calculate stats on seqarr
    all_ns = np.all(seqarr == 78, axis=1)
    nsamplecov = seqarr.shape[0] - all_ns.sum()
    nsnps = count_snps(seqarr)

    # drop all N samples(?)
    pnames = pnames[~all_ns]
    pseqs = seqarr[~all_ns, :]

    # return dictionary
    phydict = {}
    for name, seq in zip(pnames, pseqs):
        name = name + " " * (longname - len(name))
        seq = b"".join(seq.view("S1")).decode()
        phydict[name] = seq
    return nsnps, nsamplecov, phydict



NEX_MATRIX = """
#NEXUS
begin data;
  dimensions ntax={ntax} nchar={nchar};
  format datatype=dna interleave=yes gap=- missing=N;
  matrix
{matrix}
    ;
end;
"""

# proc = subprocess.Popen([
#     self.raxml_binary, 
#     "--msa", fname, 
#     "--model", "JC", 
#     "--threads", "1", 
#     "--redo",
#     ], 
#     stderr=subprocess.PIPE, 
#     stdout=subprocess.PIPE,
# )
# out, _ = proc.communicate()
# if proc.returncode:
#     raise Exception("raxml error: {}".format(out.decode()))
# tre = toytree.tree(fname + ".raxml.bestTree")