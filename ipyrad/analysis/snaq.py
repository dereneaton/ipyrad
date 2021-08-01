#!/usr/bin/env python

"""
Wrapper tool to conveniently call snaq from ipa.
"""

from typing import Optional
import os
import subprocess as sps
from loguru import logger
import numpy as np
import toytree
from ipyrad.assemble.utils import IPyradError
from toytree.utils.network import parse_network

# TODO: check pinky for previous examples.


class Snaq:
    """
    Wrapper to run simple snaq analyses on a list of gene trees.
    The input can be either a file with newick trees on separate 
    lines, or a list of newick strings, or a list of toytree 
    objects, or a DataFrame containing a column labeled .tree.

    It is assumed that Julia is installed in your $PATH.

    Parameters:
    -----------
    gtrees (str or list)
        A tree table inferred from ipa.tree_slider as a Dataframe
        or a filepath to a CSV file.
    name (str)
        Name prefix for output files.
    workdir (str)
        Directory for output files.
    bootsfile (str or None)
    
    ...
    """
    def __init__(
        self, 
        gtrees: str,
        netin: str,
        nedges: int,
        name: str="test",
        workdir: str="analysis-snaq",
        seed: Optional[int]=None,
        nruns: int=4,
        nproc: int=4,
        # cftable=None,
        force: bool=False,
        path_to_julia: Optional[str]=None
        ):

        # params
        self.name = name
        self.gtrees = gtrees
        self.netin = netin
        self.nedges = int(nedges)
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.nruns = nruns
        self.nproc = nproc
        self.seed = (np.random.randint(int(1e7)) if seed is None else seed)
        self.force = force
        self.binary = path_to_julia
        # self.cftable = cftable

        # i/o
        self.in_gt = os.path.realpath(os.path.expanduser(self.gtrees))
        self.in_net = os.path.realpath(os.path.expanduser(self.netin))
        self.io_table = os.path.join(self.workdir, self.name + ".CFs.csv")
        self.out_net = os.path.join(self.workdir, self.name)
        self.io_script = os.path.join(self.workdir, self.name + '.jl')

        # final result
        self.out_log = os.path.join(self.workdir, self.name + '.snaq.log')
        self.out_net = os.path.join(self.workdir, self.name + '.snaq')
        self.tree = None
        self.admix = None
        
        # prep
        self._check_binary()
        self._expand_script()


    def _check_binary(self):
        """
        Check that java is installed and get a tmp binary if needed.
        """
        # check for java
        if not self.binary:
            cmd = ["which", "julia"]
            proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
            comm = proc.communicate()
            if proc.returncode:
                raise IPyradError(f"julia not found: {comm[0].decode()}")
            if not comm[0]:
                raise IPyradError("julia must be installed and in your $PATH")
            self.binary = comm[0].decode()
            logger.info(f"julia path: {self.binary}")


    def _expand_script(self):
        """
        When the command is run we also save stderr to a log file.
        """
        os.makedirs(self.workdir, exist_ok=True)
        expand = {
            "nproc": self.nproc,
            "nruns": self.nruns,
            "io_table": self.io_table,
            "in_net": self.in_net,
            "nedges": self.nedges,
            "out_net": self.out_net,
            "seed": self.seed,
            "gtree_input": self.in_gt,
        }
        self._setup = SETUP.format(**expand)
        self._run = SCRIPT.format(**expand)

        # remove existing cf table
        if self.force:
            if os.path.exists(self.io_table):
                os.remove(self.io_table)

        # if table already exists then use it
        if os.path.exists(self.io_table):
            print("using existing CF table: {}".format(self.io_table))
            self._script = self._run

        else:
            self._script = self._setup + "\n" + self._run

        with open(self.io_script, 'w') as out:
            out.write(self._script)


    def _get_command(self):
        """ base command """
        cmd = [self.binary, self.io_script]
        return cmd


    def print_command(self):
        """
        Print the command line script
        """
        self._expand_script()
        print(self._script)


    def run(self):
        """
        Call SNAQ julia script 
        """
        print("[SNAQ v.x.y]")
        print("[nproc = {}]".format(self.nproc))
        print("julia {}".format(self.io_script))

        # setup the comamnd 
        proc = sps.Popen(
            self._get_command(), 
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        comm = proc.communicate()
        if proc.returncode:
            print("SNAQ Error:\n", comm[0].decode())
            raise IPyradError(
                "SNAQ Error: see .jl script and .err file in workdir")

        # try loading the tree result
        with open(self.out_log, 'r') as inlog:
            maxnet = inlog.read().split("MaxNet is ")[-1]
        self.tree, self.admix = parse_network(maxnet)

        # report result file
        print("inferred network written to ({})".format(self.out_net))



SCRIPT = """
#!/usr/bin/env julia

# check for required packages
using Pkg
Pkg.add("PhyloNetworks")
Pkg.add("CSV")

# parallelize
using Distributed
addprocs({nproc})

# load packages
using CSV, DataFrames
@everywhere using PhyloNetworks

# load quartet-CF object from table
# df_sp = CSV.read("{io_table}", categorical=false);
# d_sp = readTableCF!(df_sp);
d_sp = readTableCF("{io_table}")

# load starting network
netin = readTopology("{in_net}")

# infer the network
snaq!(netin, d_sp, hmax={nedges}, filename="{out_net}", seed={seed}, runs={nruns})
"""


SETUP = """
#!/usr/bin/env julia

# load required packages
using PhyloNetworks
using CSV

# load gene trees and starting tree
gtrees = readMultiTopology("{gtree_input}");

# count quartet CFs
q, t = countquartetsintrees(gtrees);

# reshape into dataframe
cfdf = writeTableCF(q, t);

# save table
CSV.write("{io_table}", cfdf);
"""
