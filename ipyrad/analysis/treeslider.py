#!/usr/bin/env python

"Sliding window (or sampling window) for phylo inference"

# py2/3 compat
from __future__ import print_function

# standard
import os
import sys
import time
import glob
import shutil
import tempfile

# third party
import h5py
import pandas as pd
import numpy as np

# internal librries
from .raxml import Raxml as raxml
from .mrbayes import MrBayes as mrbayes
from .window_extracter import WindowExtracter as window_extracter
from .utils import ProgressBar
from ..core.Parallel import Parallel
from ..assemble.utils import IPyradError

# do not require
try:
    import toytree
except ImportError:
    pass

_MISSING_TOYTREE = """
This ipyrad.analysis tool requires the dependency 'toytree'. 
You can install it with the following command from a terminal:

conda install toytree -c eaton-lab
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
        imap=None,
        minmap=None,
        rmincov=0.0,
        consensus_reduce=False,
        inference_method="raxml",
        inference_args={},
        quiet=False,
        scaffold_minlen=0,
        keep_all_files=False,
        **kwargs
        ):

        # check installations
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # report bad arguments
        if kwargs:
            print(
                "Warning: Some arg names are not recognized and may have "
                "changed. Please check the documentation:\n"
                "{}".format(kwargs))

        # store attributes
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))
        self.keep_all_files = keep_all_files

        # work
        self.scaffold_idxs = scaffold_idxs
        self.window_size = (int(window_size) if window_size else None)
        self.slide_size = (int(slide_size) if slide_size else None)
        self.minsnps = minsnps
        self.imap = imap
        self.mincov = mincov
        self.minmap = minmap
        self.rmincov = float(rmincov if rmincov else 0.0)
        self.consensus_reduce = consensus_reduce
        self.inference_method = inference_method
        self.inference_args = inference_args
        self.quiet = quiet
        self.scaffold_minlen = scaffold_minlen
        self._nexus = bool("raxml" not in self.inference_method)

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

        # if entered then only grab idxs that are in the scaff table
        elif isinstance(self.scaffold_idxs, (list, tuple, set)):
            idxs = sorted(self.scaffold_idxs)
            idxs = [i for i in idxs if i in self.scaffold_table.index]
            self.scaffold_idxs = idxs
        # same ...
        elif isinstance(self.scaffold_idxs, int):
            self.scaffold_idxs = [self.scaffold_idxs]

        # do not allow indices beyond existing...
        self.scaffold_idxs = (
            self.scaffold_idxs[:self.scaffold_table.index.max() + 1])
        # if self.scaffold_minlen:
        # self.scaffold_idxs = np.array(self.scaffold_idxs)[self.mask_minlen]

        # build the tree table from the scaffolds, windows, and slides.
        if self.tree_table is None:
            self._init_tree_table()


    def _print(self, message):
        if not self.quiet:
            print(message)


    def show_inference_command(self, show_full=False):
        """
        Shows the inference command (and args if show_full=True).
        """

        # show raxml command and args
        if self.inference_method == "raxml":
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

        # show mrbayes command and args
        elif self.inference_method == "mb":

            mb = mrbayes(
                data=self.data + ".nex",
                name="temp_" + str(os.getpid()),
                workdir=tempfile.gettempdir(),
                **self.inference_args
            )
            mb.print_command()
            if show_full:
                mb.print_nexus_string()


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
        with h5py.File(self.data, 'r') as io5:

            # parse formatting from db
            try:
                self._pnames = np.array([
                    i.decode() for i in io5["phymap"].attrs["phynames"]
                ])
            except AttributeError:
                self._pnames = np.array([
                    i for i in io5["phymap"].attrs["phynames"]
                ])

            self._longname = 1 + max([len(i) for i in self._pnames])

            # parse names and lengths from db
            try:
                scafnames = [i.decode() for i in io5["scaffold_names"][:]]
            except AttributeError:
                scafnames = [i for i in io5["scaffold_names"][:]]

            scaflens = io5["scaffold_lengths"][:]

            # organize as a DF
            self.scaffold_table = pd.DataFrame(
                data={
                    "scaffold_name": scafnames,
                    "scaffold_length": scaflens,
                }, 
                columns=["scaffold_name", "scaffold_length"],
            )

            # mask for min scafflen
            if self.scaffold_minlen:
                self.scaffold_table = (
                    self.scaffold_table[self.scaffold_table.scaffold_length > self.scaffold_minlen]
                )
                # self.scaffold_table.reset_index(inplace=True)
                #self.scaffold_table.reset_index(inplace=True, drop=True)            
            # if self.scaffold_minlen:
            #     self.mask_minlen = np.array(scaflens) > self.scaffold_minlen
            #     scafnames = np.array(scafnames)[self.mask_minlen]
            #     scaflens = np.array(scaflens)[self.mask_minlen]


    def _init_tree_table(self):
        "Build DataFrame for storing results"
        dfs = []

        # TODO: minlen scaffs..., and make this faster... [not looped]?
        # faster by trimming scaffold table idxs already...
        # for scaffold in self.scaffold_table.index:
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
                "tree": "",
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
        """
        scaffs are 1-indexed in h5 phymap, 0-indexed in scaffold_table.
        I know, right?
        """
        with h5py.File(self.data, 'r') as io5:
            colnames = io5["phymap"].attrs["columns"]
            try:
                colnames = [i.decode() for i in colnames]
            except AttributeError:
                colnames = [i for i in colnames]

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
                    "Existing tree table loaded from {}; "
                    "Use force to instead overwrite.")
                    .format(self.tree_table_path))
                return

        # THREADING set to match between ipcluster and raxml 
        if self.inference_method == "raxml":
            if "T" in self.inference_args:
                self.ipcluster["threads"] = max(2, int(self.inference_args["T"]))
            self.inference_args["T"] = max(2, int(self.ipcluster["threads"]))
            threads = self.inference_args["T"]
        else:
            threads = 1           

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
        finished = []
        prog = ProgressBar(self.tree_table.shape[0], None, "inferring trees")
        prog.finished = 0
        prog.update()
        rasyncs = {}
        for idx in self.tree_table.index:

            # if continuing an existing job, skip if row already filled
            if self.tree_table.tree[idx] != "":
                prog.finished += 1
                continue

            # extract the alignment for this window (auto-generate name)
            keepdir = os.path.join(
                self.workdir, "{}-{}".format(self.name, "bootsdir"))
            ext = window_extracter(
                # name=str(np.random.randint(0, 1e15)),
                data=self.data,
                workdir=keepdir,
                scaffold_idxs=int(self.tree_table.scaffold[idx]),
                start=self.tree_table.start[idx],
                end=self.tree_table.end[idx],
                mincov=self.mincov,
                imap=self.imap,
                minmap=self.minmap,
                consensus_reduce=self.consensus_reduce,
                rmincov=self.rmincov,
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
                prog.finished += 1

            else:
                # write phylip (or nex) file to the tmpdir
                ext.run(force=True, nexus=self._nexus)

                # remote inference args
                args = [ext.outfile, self.inference_args, keepdir]

                # send remote tree inference job that will clean up itself
                if "raxml" in self.inference_method:
                    rasyncs[idx] = lbview.apply(remote_raxml, *args)
                elif "mb" in self.inference_method:
                    rasyncs[idx] = lbview.apply(remote_mrbayes, *args)
                elif self.inference_method is None:
                    pass
                else:
                    raise IPyradError(
                        "inference_method should be raxml, mb or None, you entered {}"
                        .format(self.inference_method))

            prog.update()

        # wait for jobs to finish
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

        # the tree table was written as CSV to the workdir so report it.
        self._print("tree_table written to {}".format(self.tree_table_path))    

        # if not keeping boot then remove bootsdir
        if not self.keep_all_files:
            if os.path.exists(keepdir):
                shutil.rmtree(keepdir)

        # or, write a boots file pointing to all bootsfiles
        if self.keep_all_files:
            if self.inference_args.get("N"):
                newbootsfile = os.path.join(
                    self.workdir, 
                    "{}.bootsfiles.txt".format(self.name)
                )
                blist = sorted(
                    glob.glob(os.path.join(keepdir, "RAxML_bootstrap*"))
                )
                with open(newbootsfile, 'w') as out:
                    out.write("\n".join(blist))




def remote_mrbayes(nexfile, inference_args, keepdir=None):
    """
    Call mb on phy and returned parse tree result
    """
    # convert phyfile to tmp nexus seqfile

    # if keep_all_files then use workdir as the workdir instead of tmp
    if keepdir:
        workdir = keepdir
    else:
        workdir = os.path.dirname(nexfile)

    # call mb on the input phylip file with inference args
    mb = mrbayes(
        data=nexfile,
        name="temp_" + str(os.getpid()),
        workdir=workdir,
        **inference_args
    )
    mb.run(force=True, quiet=True, block=True)

    # get newick string from result
    tree = toytree.tree(mb.trees.constre, tree_format=10).newick

    # cleanup remote tree files
    for tup in mb.trees:
        tpath = tup[1]
        if os.path.exists(tpath):
            os.remove(tpath)

    # remove the TEMP phyfile in workdir/tmpdir
    os.remove(nexfile)

    # return results
    return tree    



def remote_raxml(phyfile, inference_args, keepdir=None):
    """Call raxml on phy and returned parse tree result
    """
    # if keep_all_files then use workdir as the workdir instead of tmp
    if keepdir:
        workdir = keepdir
    else:
        workdir = os.path.dirname(phyfile)

    # call raxml on the input phylip file with inference args
    rax = raxml(
        data=phyfile,
        name=os.path.basename(phyfile).rsplit(".phy")[0],  # "temp_" + str(os.getpid()),
        workdir=workdir,
        **inference_args
    )
    rax.run(force=True, quiet=True, block=True)

    # get newick string from result
    if os.path.exists(rax.trees.bipartitions):
        tree = toytree.tree(rax.trees.bipartitions).write()
    else:
        tree = toytree.tree(rax.trees.bestTree).write()

    # remote tree files
    if keepdir is None:
        for tfile in rax.trees:
            tpath = getattr(rax.trees, tfile)
            if os.path.exists(tpath):
                os.remove(tpath)

    # remove the TEMP phyfile in workdir/tmpdir
    os.remove(phyfile)    

    # return results
    return tree



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
