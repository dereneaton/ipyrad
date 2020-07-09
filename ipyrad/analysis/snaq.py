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
import pandas as pd
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
        **kwargs):

        # params
        self.name = name
        self.gtrees = gtrees
        self.netin = netin
        self.nedges = int(nedges)
        self.workdir = os.path.realpath(os.path.expanduser(workdir))

        # i/o
        self.gt_infile = os.path.join(self.workdir, self.name + ".gtin")
        self.net_infile = os.path.join(self.workdir, self.name + ".netin")
        self.out_table = os.path.join(self.workdir, self.name + ".CFs.csv")
        self.out_net = os.path.join(self.workdir, self.name + ".networks.nwk")

        # prep
        self._check_binary()
        self._write_treefiles()
        self._expand_script()


        self._check_args()
        self._parse_data_to_tmpfile()
        self._write_mapfile()



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
            "nedges": self.nedges,
            "in_net": self.netin,
            "out_net": self.netout,
            "gtree_input": self.gtrees,
            "out_table": self.out_table
        }
        self._script = SCRIPT.format(**expand)



    def null(self):
        # base command
        cmd = ["julia", "-jar", self.binary]

        # additional args
        for key, val in self._kwargs.items():
            if val is not None:
                if val is True:
                    cmd.extend([str(key)])
                elif val is False:
                    pass
                else:
                    cmd.extend([str(key), str(val)])
        return cmd



    def print_command(self):
        self._expand_script()
        print(self._script)



    def _parse_data_to_tmpfile(self):
        """
        Input can be a CSV file or DataFrame or list to 
        a bootfile that points the location of many bootstrap tree files.
        """
        # shorthand code
        data = self.data

        # data is a .tree series possiblly with NaN values.
        if isinstance(data, list):
            treelist = data

        # it is a filepath string
        elif isinstance(data, (str, bytes)):
            data = pd.read_csv(data)
            treelist = data[data.tree.notna()].tree.tolist()

        # assume this is the treeslider dataframe output with .tree column
        elif isinstance(data, pd.DataFrame):
            treelist = data[data.tree.notna()].tree.tolist()

        else:
            raise IPyradError("input format should be list or DataFrame.")

        # write to tmpfile
        self._tmptrees = os.path.join(self.workdir, "tmptrees.txt")
        with open(self._tmptrees, 'w') as out:
            out.write("\n".join(treelist))



    def run(self):
        """
        Call Astral command ()
        """
        print("[Julia SNAQ]")

        # setup the comamnd 
        proc = sps.Popen(
            self._get_command(), 
            stderr=sps.STDOUT, 
            stdout=sps.PIPE,
        )
        comm = proc.communicate()
        if proc.returncode:
            print("Astral Error:\n", comm[0].decode())
            raise IPyradError(
                "Astral Error: your command string was:\n{}"
                .format(" ".join(self._get_command())))

        # store stderr to logfile
        with open(self.logfile, 'w') as out:
            out.write(comm[0].decode())

        # cleanup 
        if os.path.exists(self._tmptrees):
            os.remove(self._tmptrees)

        # try loading the tree result
        self.tree = toytree.tree(self.treefile)

        # report result file
        print("inferred tree written to ({})".format(self.treefile))




SCRIPT = """

#!/usr/bin/env julia

# check for required packages
using Pkg
Pkg.add("PhyloNetworks")
Pkg.add("CSV")

# load required packages
using PhyloNetworks
using CSV

# i/o handles
io_gtrees = {gtree_input}
io_table = {out_table}
io_netin = {in_net}
io_netout = {out_net}

# load gene trees and starting tree
gtrees = readMultiTopology(io_gtrees);
netin = readTopology(io_netin)

# count quartet CFs and save table
q, t = countquartetsintrees(gtrees);
cfdf = writeTableCF(q, t);
qtcf = readTableCF(cfdf);
CSV.write(io_table, qtcf);

# infer the zero edge network
snaq!(netin, qtcf, hmax={nedges}, filename=netout, seed=1234)
"""
