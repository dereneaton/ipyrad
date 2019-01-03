#!/usr/bin/env python

import os
import sys
import subprocess as sps
from ..assemble.utils import IPyradWarningExit


class Bins(object):
    """
    Get binaries that should be shipped with ipyrad
    """
    def __init__(self):
        self._getbins()
        self._testbins()

    def _getbins(self):
        "Find path to software versions"
    
        # Return error if system is 32-bit arch.
        if not sys.maxsize > 2 ** 32:
            raise IPyradWarningExit("ipyrad requires 64bit architecture")

        ## get binary directory
        ipyrad_path = os.path.dirname(os.path.dirname(
            os.path.abspath(os.path.dirname(__file__))))
        bin_path = os.path.join(ipyrad_path, "bin")

        ## get the correct binaries
        if 'linux' in sys.platform:
            self.vsearch = os.path.join(
                os.path.abspath(bin_path), "vsearch-linux-x86_64")
            self.muscle = os.path.join(
                os.path.abspath(bin_path), "muscle-linux-x86_64")
            self.smalt = os.path.join(
                os.path.abspath(bin_path), "smalt-linux-x86_64")
            self.bwa = os.path.join(
                os.path.abspath(bin_path), "bwa-linux-x86_64")
            self.samtools = os.path.join(
                os.path.abspath(bin_path), "samtools-linux-x86_64")
            self.bedtools = os.path.join(
                os.path.abspath(bin_path), "bedtools-linux-x86_64")
            self.qmc = os.path.join(
                os.path.abspath(bin_path), "QMC-linux-x86_64")
        else:
            self.vsearch = os.path.join(
                os.path.abspath(bin_path), "vsearch-osx-x86_64")
            self.muscle = os.path.join(
                os.path.abspath(bin_path), "muscle-osx-x86_64")
            self.smalt = os.path.join(
                os.path.abspath(bin_path), "smalt-osx-x86_64")
            self.bwa = os.path.join(
                os.path.abspath(bin_path), "bwa-osx-x86_64")
            self.samtools = os.path.join(
                os.path.abspath(bin_path), "samtools-osx-x86_64")
            self.bedtools = os.path.join(
                os.path.abspath(bin_path), "bedtools-osx-x86_64")
            ## only one compiled version available, works for all?
            self.qmc = os.path.join(
                os.path.abspath(bin_path), "QMC-osx-x86_64")


    def _testbins(self):
        # Test for existence of binaries
        for cmd in ['muscle', 'vsearch', 'smalt', 'bwa', 
                    'samtools', 'bedtools', 'qmc']:
            cmdcall = self.__getattribute__(cmd)
            assert self._cmd_exists(cmdcall), ("{} not found here: {}"
                .format(cmd, cmdcall))


    def _cmd_exists(self, cmd):
        """ check if dependency program is there """
        return sps.call("type " + cmd,
                        shell=True,
                        stdout=sps.PIPE,
                        stderr=sps.PIPE) == 0
