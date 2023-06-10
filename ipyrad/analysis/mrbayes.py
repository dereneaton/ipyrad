#!/usr/bin/python 

""" Wrapper to make automated and parallelized calls to mb.

Example
-------
>>> import ipyrad.analysis as ipa
>>> import toytree
>>> tool = ipa.mrbayes(name='test', data='...seqs_hdf5')
>>> # tool.params.[...]
>>> tool.run()
>>> tree = toytree.tree(tool.trees.contre)
>>> tree.draw()
"""

from typing import Union, TypeVar, Optional, List
from dataclasses import dataclass, field
from pathlib import Path
import os
import sys
import subprocess
import pandas as pd
import toytree
from ipyrad.analysis.mrbayes_utils import *
from ipyrad.assemble.utils import IPyradError

# first assumes mb will be in activated conda bin
BINPATH = Path(sys.prefix) / "bin"
MB_BINARY = BINPATH / "mb"
Tool = TypeVar("MrBayes")

@dataclass
class Trees:
    """Container class to hold trees output file paths."""
    tool: Tool = field(repr=False)
    constre: Path = None
    posttre: Path = None
    pstat: Path = None
    def __post_init__(self):
        self.constre = self.tool.workdir / (self.tool.name + ".nex.con.tre")
        self.posttre = self.tool.workdir / (self.tool.name + f".nex.{self.tool.runs}t")
        self.pstat = self.tool.workdir / (self.tool.name + ".nex.pstat")

@dataclass
class Params:
    ngen: int = 100000
    nruns: int = 1
    nchains: int = 4
    samplefreq: int = 1000

@dataclass
class MrBayes:
    """Convenience tool for running mrbayes analyses.
    
    Example
    -------
    >>> tool = ipa.mrbayes(name='test', data=DATA, ...)
    >>> tool.run()
    >>> tool....
    """
    name: str
    """: Name prefix used for output files."""
    data: Union[Path, str]
    """: Path to a NEXUS input data file."""
    workdir: Union[Path, str] = None
    """: Path to directory for output files."""
    clock_model: int = 0
    """: Int to select a class of relaxed molecular clock models."""
    binary: Union[Path, str] = None
    """: Path to a 'mb' binary."""
    runs: int = 2
    """: Param..."""
    force: bool = False
    """: Force overwriting of existing files."""

    params: Params = Params()
    """: Params..."""
    trees: Trees = None
    """: Params..."""

    def __post_init__(self):
        self.data = Path(self.data)
        if not self.data.exists():
            raise IOError(f"file does not exist at: {self.data}")
        self.workdir = Path(self.workdir if self.workdir else '.').expanduser().absolute()
        self.workdir.mkdir(exist_ok=True)
        self.trees = Trees(self)

    def _write_nexus(self) -> None:
        """Write a mrbayes block to a copy of the NEXUS file."""
        # get parameters for this model type
        cwargs = self.params.__dict__.copy()

        # always add I/O args
        cwargs["nexus"] = self.data
        cwargs["outname"] = self.nexus

        # is the tree topology fixed?
        if isinstance(self.constraints, toytree.ToyTree):
            self._fixed_tree()
            cwargs["treeblock"] = self.treeblock
            cwargs["topologypr"] = self.params.topologypr
            cwargs["constraints"] = ""

        elif isinstance(self.constraints, dict):          
            # set constraints from dict
            cwargs["constraints"] = self._get_constraint_str()
            cwargs["treeblock"] = ""
        else:
            cwargs["constraints"] = ""
            cwargs["treeblock"] = ""

        # expand NEXUS usign string formatting
        if self.clock_model == 1:
            self._nexstring = NEX_TEMPLATE_2.format(**cwargs)
        elif self.clock_model == 2:
            self._nexstring = NEX_TEMPLATE_3.format(**cwargs)
        else:
            self._nexstring = NEX_TEMPLATE_1.format(**cwargs)

        # write the NEX string
        with open(self.nexus, 'w', encoding="utf-8") as out:
            out.write(self._nexstring)

    def _check_binary(self) -> None:
        """Raises an error if full path to binary is missing."""
        bins = [MB_BINARY]
        if self.binary:
            bins = [Path(self.binary), MB_BINARY]
        for path in bins:
            if path.exists():
                self.binary = path
                return
        raise IOError(
            f"mb binary not found in {bins}. Enter a full path.")

    def run(self):
        """Distribute ...

        """
        # stop before trying in mrbayes
        if self.force:
            for key, oldfile in self.trees:
                if oldfile.exists():
                    oldfile.unlink()
        else:
            if self.trees.pstat.exists():
                print(
                    f"Results exist at {self.trees.pstat}. "
                    "Use force=True to overwrite.")
                return

        # rewrite nexus file in case params have been updated
        self._write_nexus_file()

        # submit it
        if not ipyclient:
            self.stdout = _call_mb([self.binary, self.nexus])

        else:
            # find all hosts and submit job to host with most available engines
            lbview = ipyclient.load_balanced_view()
            self.rasync = lbview.apply(
                _call_mb, [self.binary, self.nexus])

        # initiate random seed
        if not quiet:
            if not ipyclient:
                print("job {} finished successfully".format(self.name))

            else:               
                if block:
                    print("job {} running".format(self.name))
                    ipyclient.wait()
                    if self.rasync.successful():
                        print(
                            "job {} finished successfully"
                            .format(self.name))
                    else:
                        self.rasync.get()
                else:
                    print("job {} submitted to cluster".format(self.name))





class OLDMrBayes:
    """MrBayes analysis utility function for running simple commands. 

    Parameters
    ----------
    data: str
        The phylip formated sequence file (.phy from ipyrad).
    name: str
        The name for this run. An alias for '-n'.
    workdir: str
        The output directory for results. An alias for '-w'. 
    force: bool
        Overwrite/rm any existing mb results with this workdir/name prefix.

    Additional optional parameters
    -------------------------------
    ngen: int
        Number of MCMC generations to run.
    sfreq: int
        Frequency to sample from MCMC chain.
    burnin: int
        Number of generations to run before sampling starts.       
    clock_model: int
        0, 1, or 2 define a set of parameters for a relaxed molecular clock
        analysis that can then be further modified. 
    constraints: dict or ToyTree
        A dictionary mapping {constraint_names: [list of tips]}. To 
        constrain an entire tree you can enter a tree and a dict will 
        automatically be built to describe the tree structure.

    Attributes
    ----------
    params: dict
        parameters for this mb run
    cmd: 
        returns the command string to run mb

    Functions
    ---------
    run()
        submits a mrbayes job locally or on an ipyparallel client cluster. 
    """    

    # init object for params
    def __init__(
        self,
        data,
        name="test",
        workdir="analysis-mb", 
        clock_model=False,
        constraints=None,
        **kwargs,
        ):

        # path attributes
        self._kwargs = {}            
        self._kwargs.update(kwargs)

        # check workdir
        self.workdir = Path(workdir if workdir else '.').expanduser().absolute()
        self.workdir.mkdir(exist_ok=True)

        # entered args
        self.clock_model = clock_model
        self.constraints = constraints
        self.name = name
        self.workdir = workdir
        self.data = Path(data).expanduser().absolute()
        self.nexus = self.workdir / (self.name + ".nex")
        self.binary = ""
        self._get_binary(self._kwargs.get("binary"))

        self.params = Params()
        defaults = {
            "ngen": 100000,
            "nruns": "1",
            "nchains": 4,
            "samplefreq": 1000,
        }
        if self.clock_model == 1:
            defaults.update(TEMPLATE_2_DICT)
        elif self.clock_model == 2:
            defaults.update(TEMPLATE_3_DICT)

        for i, j in defaults.items():
            setattr(self.params, i, j)

        # set params (overrides defaults)
        for key in self._kwargs:
            setattr(self.params, key, self._kwargs[key])

        # attributes
        self.rasync = None
        self.stdout = None
        self.stderr = None

        # results files        
        self.trees = Params()
        runs = ("" if int(self.params.nruns) < 2 else "run1.")

        self._write_nexus_file()

        # check attribute for existing results at this name.
        if self.result_files:
            print(
                "Existing results loaded for run [{}], see .trees attribute."
                .format(len(self.result_files), self.name)
            )


    def _write_nexus_file(self):
        """Write a mrbayes block to a copy of the NEXUS file."""
        # get parameters for this model type
        cwargs = self.params.__dict__.copy()

        # always add I/O args
        cwargs["nexus"] = self.data
        cwargs["outname"] = self.nexus

        # is the tree topology fixed?
        if isinstance(self.constraints, toytree.ToyTree):
            self._fixed_tree()
            cwargs["treeblock"] = self.treeblock
            cwargs["topologypr"] = self.params.topologypr
            cwargs["constraints"] = ""

        elif isinstance(self.constraints, dict):          
            # set constraints from dict
            cwargs["constraints"] = self._get_constraint_str()
            cwargs["treeblock"] = ""
        else:
            cwargs["constraints"] = ""
            cwargs["treeblock"] = ""

        # expand NEXUS usign string formatting
        if self.clock_model == 1:
            self._nexstring = NEX_TEMPLATE_2.format(**cwargs)
        elif self.clock_model == 2:
            self._nexstring = NEX_TEMPLATE_3.format(**cwargs)
        else:
            self._nexstring = NEX_TEMPLATE_1.format(**cwargs)

        # write the NEX string
        with open(self.nexus, 'w', encoding="utf-8") as out:
            out.write(self._nexstring)


    def print_command(self):
        print("{} {}".format(self.binary, self.nexus))


    def print_nexus_string(self):
        "update nexus string and print"
        self._write_nexus_file(write=False)
        print(self._nexstring)


    @property 
    def command(self):
        return "{} {}".format(self.binary, self.nexus)


    @property
    def nexus_string(self):
        "update nexus string and return"        
        self._write_nexus_file(write=False)        
        return self._nexstring


    @property
    def result_files(self):
        "returns a list of files that have finished structure"
        resfiles = [i[1] for i in self.trees if os.path.exists(i[1])]
        return resfiles


    @property
    def convergence_stats(self):
        if not os.path.exists(self.trees.pstat):
            print("no stats available")

        else:
            stats = pd.read_csv(
                self.trees.pstat,
                sep="\t", 
                skiprows=1, 
                index_col=0,
            )
            return stats


    def run(
        self, 
        ipyclient=None, 
        quiet=False,
        force=False,
        block=False,
        ):
        """
        Submits mrbayes job to run. If no ipyclient object is provided then 
        the function will block until the mb run is finished. If an ipyclient
        is provided then the job is sent to a remote engine and an asynchronous 
        result object is returned which can be queried or awaited until it 
        finishes.

        Parameters
        -----------
        ipyclient:
            Not yet supported... 
        quiet: 
            suppress print statements
        force:
            overwrite existing results files with this job name. 
        block:
            will block progress in notebook until job finishes, even if job
            is running on a remote ipyclient.
        """

        # check for input data file
        if not os.path.exists(self.data):
            raise IPyradError("data file not found {}".format(self.data))

        # stop before trying in mrbayes
        if force:
            for key, oldfile in self.trees:
                if os.path.exists(oldfile):
                    os.remove(oldfile)
        if os.path.exists(self.trees.pstat):
            print("Error Files Exist: set a new name or use Force flag.\n{}"
                  .format(self.trees.pstat))
            return 

        # rewrite nexus file in case params have been updated
        self._write_nexus_file()

        # submit it
        if not ipyclient:
            self.stdout = _call_mb([self.binary, self.nexus])

        else:
            # find all hosts and submit job to host with most available engines
            lbview = ipyclient.load_balanced_view()
            self.rasync = lbview.apply(
                _call_mb, [self.binary, self.nexus])

        # initiate random seed
        if not quiet:
            if not ipyclient:
                print("job {} finished successfully".format(self.name))

            else:               
                if block:
                    print("job {} running".format(self.name))
                    ipyclient.wait()
                    if self.rasync.successful():
                        print(
                            "job {} finished successfully"
                            .format(self.name))
                    else:
                        self.rasync.get()
                else:
                    print("job {} submitted to cluster".format(self.name))



    def _get_binary(self, search_first=None):
        """ find binaries available"""

        # check for binary
        list_binaries = [
            search_first,
            os.path.join(sys.prefix, "bin", "mb"),
        ]
        list_binaries = [i for i in list_binaries if i]

        # check user binary first, then backups
        for binary in list_binaries:
            proc = subprocess.Popen(
                ["which", binary],
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT).communicate()
            # update the binary
            if proc[0]:
                self.binary = binary

        # if none then raise error
        if not proc[0]:
            raise Exception(
                "cannot find mb; "
                "run 'conda install mrbayes -c conda-forge -c bioconda'")



    def _get_constraints_from_tree(self):
        """
        Returns a dictionary mapping node idx to a list of names. This can 
        be entered to the ipa.mb object as a constraints argument to constrain
        an entire tree.
        """
        constraints = {}
        for node in self.constraints.treenode.traverse():

            # get all tips descendant
            lvs = [i.name for i in node.get_leaves()]

            if len(lvs) > 1:

                # constraint string
                constraints[node.idx] = lvs
        return constraints



    def _fixed_tree(self):
        """
        Based on a thread on GitHub it is recommended that for fixing a 
        tree you simply set the parameters affecting topoology to zero and 
        set the starting tree to your desired tree. This is better than 
        setting a constraint on every node.

        https://github.com/NBISweden/MrBayes/issues/38

        """
        # write a tree block in the NEXUS
        self.treeblock = (
            "begin trees;\n  tree fixedtree = {}\nend;".
            format(self.constraints.write(tree_format=9))
        )

        # set topologypr to use fixed tree
        self.params.topologypr = "fixed(fixedtree)"




    def _get_constraint_str(self):
        """
        constraint <constraint name> <probability-x> = <list of taxa>

        constraint example 100 = taxon_2 taxon_3;
        prset topologypr = example;
        """
        # self.constraints is a dictionary
        constraints = []
        for key in self.constraints:

            # const str
            constraint = "constraint idx-{} 100 = {};".format(
                key, " ".join(self.constraints[key])
            )

            # setting string
            setit = "prset topologypr = constraints (idx-{});".format(key)

            constraints.append(constraint)
            constraints.append(setit)

        # return empty string if no constraints
        return "\n".join(constraints)



def _call_mb(cmd: List[str]) -> str:
    """Run subprocess on mb command and return stdout."""
    kwargs = dict(stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    with subprocess.Popen(cmd, **kwargs) as proc:
        comm = proc.communicate()
    if proc.returncode:
        raise IPyradError(comm[0].decode())
    return comm[0].decode()


if __name__ == "__main__":

    mb = MrBayes(name='hi', workdir='/tmp', data="/home/deren/Downloads/T88270.nex")
    print(mb)