#!/usr/bin/env python

"Sliding window (or sampling window) for phylo inference"

# py2/3 compat
from __future__ import print_function

# standard
import os
import sys
import time
import datetime
import tempfile

# third party
import pandas as pd
import numpy as np

# internal librries
from .raxml import Raxml as raxml
from .mrbayes import MrBayes as mrbayes
from .utils import count_snps
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


class TreeSlider():
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
        inference_method="raxml",
        inference_args={},
        # no_missing_taxa=True,
        # imap=None,
        # minmap=None,
        ):

        # check installations
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # store attributes
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))

        # work
        self.scaffold_idxs = scaffold_idxs
        self.window_size = window_size
        self.slide_size = slide_size
        self.minsnps = minsnps

        # use user name else create one
        if not self.name:
            self.name = "all-scaffolds"

        # get outfile name
        self.tree_table_path = os.path.join(
            self.workdir,
            "{}.tree_table.csv".format(self.name))

        # inference 
        self.inference_method = inference_method
        self.inference_args = inference_args

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

        # if not window then slide is set to window
        if (not self.window_size) or (not self.slide_size):
            self.slide_size = self.window_size

        # # parse attributes
        # self.imap = imap
        # self.minmap = minmap
        # self.rmap = {}
        # if self.imap:
        #     for k, v in self.imap.items():
        #         for i in v:
        #             self.rmap[i] = k

        # ensure workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # get outfile name
        if self.name:
            tree_table_path = os.path.join(
                self.workdir,
                "{}.tree_table.csv".format(self.name))
            if os.path.exists(tree_table_path):
                self.tree_table = pd.read_csv(tree_table_path, sep="\t")
                print("existing results loaded from {}".format(
                    tree_table_path))

        # # fill mindict
        # if not minmap:
        #     if imap:
        #         self.minmap = {i: 1 for i in self.imap}

        # parsed attributes
        self.scaffold_table = None
        self.tree_table = None
        self.phymap = None
        self._pnames = None
        self._parameter_check()

        # default to all scaffolds if none entered.
        self._parse_scaffolds()
        if not self.scaffold_idxs:
            self.scaffold_idxs = self.scaffold_table.index.tolist()
        if isinstance(self.scaffold_idxs, (list, tuple, set)):
            self.scaffold_idxs = sorted(self.scaffold_idxs)

        # build the tree table from the scaffolds, windows, and slides.
        self._parse_tree_table()


    def show_inference_command(self, show_full=False):
        # show a warning if user entered threads
        if "T" in self.inference_args:
            print("\n".join([
                "threads argument (T) must be set in ipcluster. ",
                "e.g., tet.ipcluster['threads'] = 4"
            ]))

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


    def _parse_scaffolds(self):
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


    def _parse_tree_table(self):
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
                "nsnps": np.nan,
                "nsamplecov": np.nan,
                "tree": np.nan,
                }, 
                columns=["scaffold", "start", "end", "nsnps", "nsamplecov", "tree"], 
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
        Run command for tree slider.
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
        """
        Hidden func to run parallel job.
        """

        # do not overwrite tree table
        if os.path.exists(self.tree_table_path):
            if not force:
                print((
        "Existing tree table loaded from {}; Use force to instead overwrite."
        .format(self.tree_table_path)))
                return

        # THREADING
        if "T" in self.inference_args:
            self.ipcluster["threads"] = int(self.inference_args["T"])

        # load balance parallel jobs 2-threaded
        threads = int(max(2, self.ipcluster["threads"]))
        lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::threads])

        # distribute jobs on client
        time0 = time.time()
        rasyncs = {}

        # initial progress ticker to run during job submission
        print(
            "building database: nwindows={}; minsnps={}"
            .format(
                self.tree_table.shape[0],
                self.minsnps,
            ))

        # subset phymap for this scaffold
        io5 = h5py.File(self.data, 'r')
        scaffs = io5["phymap"][:, 0]

        # submit jobs: (fname, scafidx, minpos, maxpos, minsnps, )
        jidx = 0
        for scaff in self.tree_table.scaffold.unique():

            # load phymap for each scaff at a time
            phymap = io5["phymap"][scaffs == scaff + 1, :]
            subtree = self.tree_table[self.tree_table.scaffold == scaff]

            # submit jobs each window at a time
            for idx in subtree.index:

                # get window margins
                start, stop = (
                    subtree.loc[idx, ["start", "end"]].astype(int))

                # subset phymap for this window range
                mask = (phymap[:, 3] > start) & (phymap[:, 4] < stop)
                cmap = phymap[mask, :]

                # only if there is RAD sequence data in this window send job
                if cmap.size:
                    wmin = cmap[:, 1].min()
                    wmax = cmap[:, 2].max()

                    # store async result
                    args = (
                        self.data, wmin, wmax, self.minsnps, 
                        threads, 4, 
                        self.inference_method, self.inference_args,
                    )
                    rasyncs[jidx] = lbview.apply(remote_tree_inference, *args)

                # advance fill counter, leaves NaN in rows with no data.
                jidx += 1

        # close database
        io5.close()

        # track progress and save result table
        message = "inferring raxml trees"
        self._track_progress_and_store_results(rasyncs, time0, message)
        self.tree_table.to_csv(self.tree_table_path)


    def _track_progress_and_store_results(self, rasyncs, time0, message):
        # track progress and collect results.
        nwindows = self.tree_table.shape[0]
        done = 0
        while 1:
            finished = [i for i in rasyncs if rasyncs[i].ready()]
            for idx in finished:
                if rasyncs[idx].successful():
                    self.tree_table.iloc[idx, 3:] = rasyncs[idx].get()
                    self.tree_table.to_csv(self.tree_table_path)
                    del rasyncs[idx]
                    done += 1
                else:
                    raise Exception(rasyncs[idx].get())

            # progress
            progressbar(done, nwindows, time0, message)
            time.sleep(0.5)
            if not rasyncs:
                print("")
                break



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
    tree = np.nan
    nsamplecov = np.nan
    nsnps = np.nan

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


def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
          .format(hashes + nohash, int(progress), elapsed, message),
          end="")
    sys.stdout.flush()





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