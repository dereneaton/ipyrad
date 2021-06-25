#!/usr/bin/env python

"""
Sliding window (or sampling window) for phylo inference
"""

# standard
import os
import time
import glob
import shutil
import tempfile
import traceback
from typing import List, Dict, Optional, Union

# third party
import pandas as pd
import numpy as np
import toytree
from loguru import logger

# internal librries
from ipyrad.core.parallel import Cluster
from ipyrad.core.progress_bar import AssemblyProgressBar
from ipyrad.assemble.utils import IPyradError
from .raxml import Raxml as raxml
# from .mrbayes import MrBayes as mrbayes
from .window_extracter import window_extracter


class TreeSlider:
    """
    Performs phylogenetic inference on subsampled and filtered RADseq
    data in genomic windows. Uses the seqs.hdf5 database file from 
    ipyrad as input. Data must have been assembled by mapping to a
    reference genome. If no window size is entered then entire 
    scaffolds are used as windows. 

    Example:
    --------
    tool = ipa.tree_slider(
        data=data, 
        mincov=4, 
        minsnps=100, 
        window_size=1e5,
        slide_size=1e5,
        scaffold_idxs=[0, 1, 2, 4, 5],
    )
    tool.run(cores=8)

    Parameters:
    ------------
    data: (str)
        Database file (.seqs.hdf5) produced by ipyrad.
    name: (str)
        Name prefix for output files.
    workdir: (str)
        Directory for output files.
    imap: (dict) 
        Optional dictionary mapping population names to lists of 
        sample names. This can be used to apply minmap filters for
        sample coverage within groups, or for consensus_reduce 
        sampling of sites (see window_extracter docs).
    minmap: (dict) 
        Optional dictionary mapping population names to min sample
        coverage values as int or float (proportion).
    minsnps: int
        Minimum number of SNPs required to include window in analysis.
    inference_method: str
        'raxml' or 'mb'
    inference_args: dict
        A dictionary mapping method param names to their values. 
        See examples or contact developers to enable unsupported
        options.
    scaffold_minlen: int
        The minimum allowed length of post-filtered alignments from
        an extracted window after site filtering is applied. Windows
        that do not pass filtering are skipped (NA).
    keep_all_files: bool
        All alignment and tree inference files are saved in the 
        workdir instead of being retained only temporarily.
    quiet: bool
        Suppress progress info
    """
    def __init__(
        self,
        data: str,
        name: Optional[str]=None,
        workdir: str="./analysis-treeslider",
        window_size: Optional[int]=None, 
        slide_size: Optional[int]=None,
        scaffold_idxs:Optional[List[int]]=None,
        minsnps: int=1,
        mincov: int=0,
        imap: Dict[str, List[str]]=None,
        minmap: Dict[str, Union[float,int]]=None,
        rmincov: Union[float, int]=0.0,
        consensus_reduce: bool=False,
        inference_method: str="raxml",
        inference_args: Dict[str, str]=None,
        quiet: bool=False,
        minlen_alignment: int=0,
        keep_all_files: bool=False,
        **kwargs
        ):

        # report bad arguments
        if kwargs:
            print(
                "Warning: Some arg names are not recognized and may have "
                "changed. Please check the documentation:\n"
                "{}".format(kwargs))

        # store attributes
        self.name = name if name is not None else "tree_slider"
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.data = os.path.realpath(os.path.expanduser(data))
        self.keep_all_files = bool(keep_all_files)

        # work
        self.scaffold_idxs = scaffold_idxs
        self.window_size = int(window_size) if window_size else None
        self.slide_size = int(slide_size) if slide_size else None
        self.minsnps = minsnps
        self.imap = imap
        self.mincov = mincov
        self.minmap = minmap
        self.rmincov = float(rmincov if rmincov else 0.0)
        self.consensus_reduce = consensus_reduce
        self.inference_method = inference_method
        self.inference_args = inference_args if inference_args is not None else {}
        self.quiet = quiet
        self.minlen_alignment = minlen_alignment

        # get outfile name
        self.tree_table_path = os.path.join(
            self.workdir, f"{self.name}.tree_table.csv")           

        # to-be parsed attributes
        self.phymap = None

        # checks params and loads tree table if existing.
        self._parameter_check()

        # get scaffold names and lengths
        wex = window_extracter(self.data, scaffold_idxs=scaffold_idxs)
        self.scaffold_table = wex.scaffold_table
        self.scaffold_idxs = wex.scaffold_idxs

        # build the tree table from the scaffolds, windows, and slides.
        self.tree_table = None
        self._init_tree_table()


    def show_inference_command(self, show_full=False):
        """
        Shows the inference command (and args if show_full=True).
        """
        # show raxml command and args
        if self.inference_method == "raxml":
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
        """
        Checks data file, slide size, and window size
        """
        assert os.path.exists(self.data), "database file not found"
        assert self.data.endswith(".seqs.hdf5"), (
            "data must be '.seqs.hdf5' file.")
        assert self.inference_method in ['raxml', 'mb'], (
            "only raxml and mb are supported inference methods.")

        # if not window then slide is set to window
        if (not self.window_size) or (not self.slide_size):
            self.slide_size = self.window_size
            logger.info(
                "Setting slide_size = window_size because a "
                "value was not entered for both options.")

        # set a default threading
        if self.inference_method == "raxml":
            self.inference_args['T'] = self.inference_args.get("T", 4)


    def _init_tree_table(self):
        """
        Build DataFrame for storing results
        """
        dataframes = []
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
            data = pd.DataFrame(data={
                "scaffold": scaffold,
                "start": starts,
                "end": ends,
                "snps": 0,  # np.nan,
                "sites": 0,
                "samples": 0,
                "missing": 0.0,  # np.nan,
                "tree": "",
                }, 
            )
            dataframes.append(data)

        # concat data from one or more scaffolds
        self.tree_table = pd.concat(dataframes)
        self.tree_table.reset_index(drop=True, inplace=True)


    def run(self, cores=None, ipyclient=None, force=False):
        """
        Distribute tree slider jobs in parallel on N cores. The 
        number of threads per job will be detected from the inference
        args to start the appropriate number of concurrent jobs.

        Parameters:
        -----------
        cores: int
            Number of cores to use for parallelization. The tool will
            automatically manage a local cluster with N cores.

        ipyclient: (type=ipyparallel.Client); Default=None. 
            If you started an ipyclient manually then you can 
            connect to it and use it to distribute jobs here. This
            option may be used to setup a cluster over multiple
            machines/nodes using MPI.

        force: (type=bool); Default=False.
            Force overwrite of existing output with the same name.
        """
        # ensure outdir exists
        os.makedirs(self.workdir, exist_ok=True)

        # do not overwrite tree table
        if os.path.exists(self.tree_table_path) and (not force):
            logger.warning(
                "Results already exist for this name at "
                f"{self.tree_table_path}. Use force=True to overwrite.")
            return

        # init the ipyparallel cluster class wrapper
        cluster = Cluster(quiet=self.quiet)
        try:
            # establish connection to a new or running ipyclient
            cluster.start(cores=cores, ipyclient=ipyclient)
            self._run(ipyclient=cluster.ipyclient)

        except KeyboardInterrupt:
            logger.warning("keyboard interrupt by user, cleaning up.")

        # AssemblyProgressBar logs the traceback
        except IPyradError as inst:
            logger.error(f"An error occurred:\n{inst}")
            print("An error occurred, see logfile and below.")
            raise

        # logger.error logs the traceback
        except Exception as inst:
            logger.error(
                "An unexpected error occurred, see logfile "
                f"and trace:\n{traceback.format_exc()}")
            raise

        finally:
            cluster.cleanup_safely(None)


    def _run(self, ipyclient, quiet=False):
        """
        Hidden func to distribute jobs that is wrapped inside Parallel.
        """
        # THREADING set to match between ipcluster and raxml 
        nthreads = 1
        if self.inference_method == "raxml":
            if "T" in self.inference_args:
                nthreads = self.inference_args["T"]

        # load balance parallel jobs 2-threaded
        lbview = ipyclient.load_balanced_view(targets=ipyclient.ids[::nthreads])

        # initial progress ticker to run during job submission
        logger.info(
            f"building database: nwindows={self.tree_table.shape[0]}; "
            f"minsnps={self.minsnps}"
        )

        # iterate to submit jobs
        keepdir = os.path.join(self.workdir, f"{self.name}-tmpdir")
        jobs = {}        
        for idx in self.tree_table.index:

            # setup window extracter tool for this window
            wex = window_extracter(
                data=self.data,
                name=f"tree_slider_window-{idx}",
                workdir=keepdir,
                scaffold_idxs=int(self.tree_table.scaffold[idx]),
                start=self.tree_table.start[idx],
                end=self.tree_table.end[idx],
                mincov=self.mincov,
                imap=self.imap,
                minmap=self.minmap,
                consensus_reduce=self.consensus_reduce,
                rmincov=self.rmincov,
            )

            # extract the window to a file in tmpdir
            wex.run(nexus=bool(self.inference_method == "mb"))

            # fill table with filtered stats for window
            self.tree_table.loc[idx, "snps"] = wex.stats.snps_post[0]
            self.tree_table.loc[idx, "sites"] = wex.stats.sites_post[0]
            self.tree_table.loc[idx, "missing"] = wex.stats.missing_post[0]
            self.tree_table.loc[idx, "samples"] = wex.stats.samples_post[0]

            # if window has enough SNPS then submit inference job.
            cond1 = wex.stats.snps_post[0] >= self.minsnps
            cond2 = wex.stats.sites_post[0] >= self.minlen_alignment
            if cond1 & cond2:
                args = [wex.outfile, self.inference_args, keepdir]
                if "raxml" in self.inference_method:
                    jobs[idx] = lbview.apply(remote_raxml, *args)
                elif "mb" in self.inference_method:
                    jobs[idx] = lbview.apply(remote_mrbayes, *args)

        # submit jobs: (fname, scafidx, minpos, maxpos, minsnps, )
        msg1 = f"{self.inference_method} inference in sliding windows"
        msg2 = f"{self.tree_table.shape[0]} trees"
        prog = AssemblyProgressBar(jobs, msg1, msg2, quiet)
        prog.block()
        prog.check()

        # assemble results
        for idx in prog.results:
            newick = prog.results[idx]
            self.tree_table.loc[idx, 'tree'] = newick
        
        # write CSV to disk
        self.tree_table.to_csv(self.tree_table_path)
        logger.info(f"tree_table written to {self.tree_table_path}")

        # if not keeping then remove dir
        if not self.keep_all_files:
            if os.path.exists(keepdir):
                shutil.rmtree(keepdir)

        # if raxml "N" write a boots file pointing to all bootsfiles
        if self.keep_all_files:
            if self.inference_args.get("N"):
                bootsf = os.path.join(self.workdir, f"{self.name}.bootsfiles.txt")
                blist = sorted(
                    glob.glob(os.path.join(keepdir, "RAxML_bootstrap*"))
                )
                with open(bootsf, 'w') as out:
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
    """
    Call raxml on phy and returned parse tree result
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
        tree = toytree.tree(rax.trees.bipartitions).newick
    else:
        tree = toytree.tree(rax.trees.bestTree).newick

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