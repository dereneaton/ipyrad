#!/usr/bin/env python

"""
Wrapper tool to conveniently call snaq from ipa
"""

# py2/3 compat
from __future__ import print_function
from builtins import range

import os
import sys
import subprocess as sps
import numpy as np
from ..assemble.utils import IPyradError

# import toytree
try:
    import toytree
except ImportError:
    pass
_MISSING_TOYTREE = """
You are missing required packages to use ipa.snaq().
First run the following conda install command:

conda install toytree -c conda-forge
"""


class Snaq:
    """
    Wrapper to run simple snaq analyses on a list of gene trees.
    The input can be either a file with newick trees on separate lines, 
    or a list of newick strings, or a list of toytree objects, or a DataFrame
    containing a column labeled .tree.

    It is assumed that Julia is installed in your $PATH.

    Parameters
    ===========
    data (str or list)

    name (str)

    workdir (str)

    bootsfile (str or None)

    ...
    """
    def __init__(
        self, 
        gtrees,
        netin,
        nedges,
        name="test",
        workdir="analysis-snaq",
        cftable=None,
        seed=None,
        nruns=4,
        nproc=4,
        **kwargs):

        # params
        self.name = name
        self.gtrees = gtrees
        self.netin = netin
        self.nedges = int(nedges)
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.nruns = nruns
        self.nproc = nproc
        self.seed = (np.random.randint(int(1e7)) if seed is None else seed)
        self.cftable = cftable

        # i/o
        self.in_gt = os.path.realpath(os.path.expanduser(self.gtrees))
        self.in_net = os.path.realpath(os.path.expanduser(self.netin))
        self.io_table = os.path.join(self.workdir, self.name + ".CFs.csv")
        self.out_net = os.path.join(self.workdir, self.name)
        self.io_script = os.path.join(self.workdir, self.name + '.jl')

        # final result
        self.out_log = os.path.join(self.workdir, self.name + '.log')
        self.out_net = os.path.join(self.workdir, self.name + '.networks')

        # prep
        self._check_binary()
        self._expand_script()



    def _check_binary(self):
        """
        Check that java is installed and get a tmp binary if needed.
        """
        # check for toytree
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # check for java
        cmd = ["which", "julia"]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        comm = proc.communicate()
        if proc.returncode:
            print(comm[0])
        if not comm[0]:
            raise IPyradError("julia must be installed and in your $PATH")




    def _expand_script(self):
        """
        When the command is run we also save stderr to a log file.
        """

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

        # if table already exists then use it
        if os.path.exists(self.io_table):
            print("using existing CF table: {}".format(self.io_table))
            self._script = self._run

        else:
            self._script = self._setup + "\n" + self._run

        with open(self.io_script, 'w') as out:
            out.write(self._script)


    def _get_command(self):
        # base command
        cmd = ["julia", self.io_script]
        return cmd


    def print_command(self):
        self._expand_script()
        print(self._script)



    def run(self):
        """
        Call SNAQ julia script 
        """
        print("[SNAQ v.x.y]")
        print("[nproc = {}]".format(self.nruns))
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
        self.tree, self.admix = toytree.utils.parse_network(maxnet)

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
using CSV
@everywhere using PhyloNetworks

# load quartet-CF object from table
df_sp = CSV.read("{io_table}", categorical=false);
d_sp = readTableCF!(df_sp);

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
