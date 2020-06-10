#!/usr/bin/env python

"""
Wrapper tool to conveniently call astral from ipa
"""

# py2/3 compat
from __future__ import print_function
from builtins import range

import os
import sys
import tempfile
import requests
import subprocess as sps
import pandas as pd
from ..assemble.utils import IPyradError

# import toytree
try:
    import toytree
except ImportError:
    pass
_MISSING_TOYTREE = """
You are missing required packages to use ipa.bpp().
First run the following conda install command:

conda install toytree -c eaton-lab
"""


class Astral:
    """
    Wrapper to run simple astral analyses on a list of gene trees.
    """
    def __init__(
        self, 
        data, 
        name="test",
        workdir="analysis-astral", 
        bootsfile=None, 
        imap=None,
        annotation=1,
        gene_resampling=False,
        nboots=None,
        **kwargs):

        # i/o
        self.name = name
        self.data = data
        self._tmptrees = None
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.binary = None
        self.mapfile = None

        # results 
        self.tree = None
        self.quartet_score = None
        self.treefile = os.path.join(self.workdir, self.name + ".tre")
        self.logfile = os.path.join(self.workdir, self.name + ".log")

        # params
        self.bootsfile = bootsfile
        self.imap = imap
        self.gene_resampling = gene_resampling
        self.annotation = annotation
        self.nboots = nboots

        # prep
        self._check_binary()
        self._check_args()
        self._parse_data_to_tmpfile()
        self._write_mapfile()

        # the kwargs for astral
        self._kwargs = {
            "-i": self._tmptrees,
            "-o": self.treefile,
            "-a": self.mapfile,
            "-b": self.bootsfile,
            "-t": self.annotation,
            "-r": self.nboots,
            "-g": self.gene_resampling, 
        }


    def _write_mapfile(self):
        """
        species_name:individual_1,individual_2,...
        """
        if self.imap:
            self.mapfile = os.path.join(
                self.workdir, 
                "{}.imap.txt".format(self.name),
            )
            maplines = [
                "{}:{}".format(i, ",".join(j))
                for (i, j) in self.imap.items()
            ]
            with open(self.mapfile, 'w') as out:
                out.write("\n".join(maplines))



    def _check_args(self):
        """
        Check for bad or incompatible arguments
        """
        # check for workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # check bootstrap args
        if self.bootsfile:

            # cannot do bootstrapping unless bootsfile exists
            assert os.path.exists(self.bootsfile), "bootsfile not found"

        #     # cannot have more replicates than the length of boots
        #     with open(self.bootsfile, 'r') as infile:
        #         nboots = sum(1 for i in infile)
        #         assert nboots >= self.nboots, (
        #             "nboots ({}) cannot be > number of bootstrap trees ({})"
        #             .format(nboots, self.nboots))

        # # not bootstrapping
        # else:
        #     # cannot do gene-resampling unless bootstrapping
        #     assert self.gene_resampling is None, (
        #         "cannot do gene resampling unless using bootsfile")



    def _get_command(self):
        """
        When the command is run we also save stderr to a log file.
        """
        # base command
        cmd = ["java", "-jar", self.binary]

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
        print(" ".join(self._get_command()))


    def _check_binary(self):
        """
        Check that java is installed and get a tmp binary if needed.
        """
        if not sys.modules.get("toytree"):
            raise ImportError(_MISSING_TOYTREE)

        # check for java
        cmd = ["which", "java"]
        proc = sps.Popen(cmd, stderr=sps.STDOUT, stdout=sps.PIPE)
        comm = proc.communicate()
        if proc.returncode:
            print(comm[0])
        if not comm[0]:
            raise IPyradError(
                "java must be installed or loaded to use Astral.")

        # paths
        uzip = "https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.3.zip"
        dzip = "Astral.5.7.3.zip"
        adir = "Astral"
        jbin = "astral.5.7.3.jar"

        # look for existing binary in tmpdir
        self.binary = os.path.join(
            tempfile.gettempdir(), adir, jbin,
        )

        # check that bpp is installed and in path            
        cmd = ['which', self.binary]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if comm:
            return 

        # bpp not found in /tmp, download it.        
        res = requests.get(uzip, allow_redirects=True)
        tmpzip = os.path.join(tempfile.gettempdir(), dzip)
        with open(tmpzip, 'wb') as tz:
            tz.write(res.content)

        # decompress tar file 
        cmd = ["unzip", tmpzip, "-d", tempfile.gettempdir()]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()
        if proc.returncode:
            print(comm[0], comm[1])

        # check that binary now can be found
        cmd = ['which', self.binary]
        proc = sps.Popen(cmd, stderr=sps.PIPE, stdout=sps.PIPE)
        comm = proc.communicate()[0]
        if comm:
            return 
        raise IPyradError("astral binary not found.")



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
        if isinstance(data, (str, bytes)):
            data = pd.read_csv(data)
            treelist = data[data.tree.notna()].tree.tolist()

        # assume this is the treeslider dataframe output with .tree column
        if isinstance(data, pd.DataFrame):
            treelist = data[data.tree.notna()].tree.tolist()

        # write to tmpfile
        self._tmptrees = os.path.join(self.workdir, "tmptrees.txt")
        with open(self._tmptrees, 'w') as out:
            out.write("\n".join(treelist))



    def run(self):
        """
        Call Astral command ()
        """
        print("[astral.5.7.3.jar]")
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
