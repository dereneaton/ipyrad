#!/usr/bin/python 

""" 
Wrapper to make simple calls to raxml 
"""

from typing import Dict

import os
import sys
import glob
import subprocess
from loguru import logger
from ipyrad.analysis.utils import Params
from ipyrad.assemble.utils import IPyradError


class Raxml:
    """
    RAxML analysis utility function. This tool makes it easy to 
    build a raxml command line string, check it, run it, and collect
    the results. 

    Example:
    ---------
    rax = ipa.raxml(
        data="test.phy"
        name="test",
        workdir="./analysis-raxml",
        T=20,
        N=500,
    )
    print(rax.command)
    rax.run()
    tree = rax.trees.bipartitions

    Parameters:
    -----------
    data: str
        The phylip formated sequence file (.phy from ipyrad). 
        An alias for raxml command arg '-s'. 
    name: str
        The name for this run. An alias for '-n'.
    workdir: str
        The output directory for results. An alias for '-w'. 
    **kwargs
        Any parameter can be entered here and will be added to the
        raxml command string. Call self.command to see the command.
        To remove a pre-set default argument set the value to None.
        For example, N=None will remove the argument N from the 
        raxml command line.

    Common raxml command options
    ----------------------------
    f: str
        (-f a) The raxml function. Default is 'a'.
    T: str
        (-T 4) The number of threads. Default is 4. 
    m: str
        (-m GTRGAMMA) The model to use.
    N: str
        (-N 100) The number of distinct starting trees from which to 
        run full search, or number of bootstrap replicates to run if 
        using -f a.
    x: str
        (-x 12345) The bootstrap random seed.
    p: str
        (-p 54321) The parsimony random seed.
    n: str
        (-n test) The prefix name for output files
    w: str
        (-w outdir) The output directory
    s: str
        (-s seq.phy) The .phy formatted sequence file. 
    o: str or list
        (-o tax1,tax2) A list of outgroup sample names or a string. 

    Attributes:
    -----------
    kwargs: dict
        A dictionary with parameters for this raxml run
    command: 
        The command string that will be used to run raxml

    Functions:
    ----------
    run()
        starts the raxml job running. You can submit to an 
        ipyparallel cluster or run on the current kernel, and 
        optionally block or not until job is finished.
    """    
    def __init__(
        self,
        data: str,
        name: str="test",
        workdir: str="analysis-raxml", 
        **kwargs: Dict[str, str],
        ):

        # path attributes
        self.kwargs = {
            "f": "a", 
            "T": 4,
            "m": "GTRGAMMA",
            "N": 100,
            "x": 12345,
            "p": 54321,
            "o": None,
            "binary": "",
        }

        # update default kwargs with user kwargs
        self.kwargs.update(kwargs)

        # drop any kwargs that user set to None
        self.kwargs = {i: j for (i, j) in self.kwargs.items() if j is not None}

        # ensure workdir exists
        if workdir is None:
            workdir = os.path.abspath(os.path.curdir)
        workdir = os.path.abspath(os.path.expanduser(workdir))          
        os.makedirs(workdir, exist_ok=True)

        # check name and data
        name = "raxml" if name is None else name
        data = os.path.abspath(os.path.expanduser(data))

        # data, name and workdir override args
        self.kwargs['n'] = self.kwargs.get("n", name)
        self.kwargs['w'] = self.kwargs.get("w", workdir)
        self.kwargs['s'] = self.kwargs.get("s", data)

        # get binary
        binaries_list = _get_binary_paths()
        if self.kwargs["binary"]:
            binaries_list = [self.kwargs["binary"]] + binaries_list
        self.kwargs['binary'] = _check_binaries(binaries_list)

        # attributes
        self.rasync = None
        self.stdout = None
        self.stderr = None

        # results files        
        self.trees = Params()
        self._paths = {
            "info": os.path.join(workdir, "RAxML_info." + name),
            "best": os.path.join(workdir, "RAxML_bestTree." + name),
            "bipartitions": os.path.join(workdir, "RAxML_bipartitions." + name),
            "bipart_blens": os.path.join(workdir, "RAxML_bipartitionsBranchLabels." + name),
            "bootstrap": os.path.join(workdir, "RAxML_bootstrap." + name),
        }
        self.trees.info = self._paths['info']


    @property
    def _command_list(self):
        """ 
        build the command list 
        """
        cmd = [self.kwargs['binary']]
        for key in self.kwargs:
            if key not in ["binary", "o"]:
                cmd.extend([f"-{key}", str(self.kwargs[key])])

            # add ougroups from a list of taxa
            if key == "o":
                cmd.extend(["-o", ",".join(self.kwargs['o'])])
        return cmd


    @property
    def command(self):
        """ returns command as a string """
        return " ".join(self._command_list)


    def run(
        self, 
        ipyclient: 'ipyparallel.Client'=None, 
        quiet: bool=False,
        force: bool=False,
        block: bool=False,
        ):
        """
        Submits raxml job to run. If no ipyclient object is provided
        then the function will block until the raxml run is finished. 
        If an ipyclient is provided then the job is sent to a remote
        engine and an asynchronous result object is returned which
        can be queried or awaited until it finishes.

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
        # force removes old files, a bit risky here if names are subsets
        if force:
            opath = os.path.join(self.kwargs['w'], f"RAxML_*.{self.kwargs['n']}")
            oldfiles = glob.glob(opath)
            for oldfile in oldfiles:
                if os.path.exists(oldfile):
                    os.remove(oldfile)

        if os.path.exists(self.trees.info):
            logger.error(
                f"File Exists ({self.trees.info}): "
                "set a new name or use force=True."
            )
            return 

        # non-remote execution (most common)
        if not ipyclient:
            self.stdout = _call_raxml(self._command_list)
            if not quiet:
                if "Overall" not in self.stdout:
                    print("Error in raxml run\n" + self.stdout)
                else:
                    print(
                        f"job {self.kwargs['n']} finished successfully.")

        else:
            # find all hosts and submit job to the host with most available engines
            lbview = ipyclient.load_balanced_view()
            self.rasync = lbview.apply(_call_raxml, self._command_list)

            # initiate random seed
            if not quiet:
                if block:
                    print("job {} running".format(self.kwargs['n']))
                    ipyclient.wait()
                    if self.rasync.successful():
                        print(
                            "job {} finished successfully"
                            .format(self.kwargs['n']))
                    else:
                        raise IPyradError(self.rasync.get())
                else:
                    print("job {} submitted to cluster".format(self.kwargs['n']))

        # make outfiles accessible if they exist
        for key in self._paths:
            path = self._paths[key]
            if os.path.exists(path):
                self.trees[key] = path
            else:
                self.trees[key] = ""



def _get_binary_paths():
    """
    check for binary in PATH
    """
    list_binaries = [
        "raxmlHPC-PTHREADS-AVX2",
        "raxmlHPC-PTHREADS-AVX",            
        "raxmlHPC-PTHREADS-SSE3",
        "raxmlHPC-PTHREADS", 
        ]
    # expand for env path
    list_binaries = [os.path.join(sys.prefix, "bin", i) for i in list_binaries]
    return list_binaries



def _check_binaries(binaries):
    """ 
    find and return a working binary
    """
    # check user binary first, then backups
    for binary in binaries:

        # call which to find 
        proc = subprocess.Popen(
            ["which", binary],
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE, 
        ).communicate()

        # if it exists then update the binary
        if proc[0]:
            return binary

    # if you get here then no binaries were found
    raise NameError(BINARY_ERROR)


def _call_raxml(command_list):
    """ 
    Call the command in subprocess
    """
    proc = subprocess.Popen(
        command_list,
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
    )
    comm, _ = proc.communicate()
    if proc.returncode:
        raise IPyradError(f"RAxML error:\n{comm.decode()}")
    return comm.decode()



BINARY_ERROR = """
  RAxML binary not found. 

  Check that you have raxml installed. For example, with conda:
  'conda install raxml -c conda-forge -c bioconda'

  If you have a different binary installed you can select it using 
  the argument 'binary'. For example:

  rax = ipa.raxml(name='test', data='test.phy', binary='raxmlHPC')
"""
