#!/usr/bin/python 

""" wrapper to make simple calls to raxml """

import os
import subprocess
from ipyrad.assemble.utils import Params

## alias
OPJ = os.path.join


class Raxml(object):
    """
    RAxML analysis utility function for running simple commands. 

    Parameters:
    -----------
    data: str
        The phylip formated sequence file (.phy from ipyrad). An alias for '-s'. 
    name: str
        The name for this run. An alias for '-n'.
    workdir: str
        The output directory for results. An alias for '-w'. 

    Additional optional parameters
    -------------------------------
    f: str
        (-f a) The raxml function. Default is 'a'.
    T: str
        (-T 4) The number of threads. Default is 4. 
    m: str
        (-m GTRGAMMA) The model to use.
    N: str
        (-N 100) The number of bootstrap replicates to run.
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
    params: dict
        parameters for this raxml run
    cmd: 
        returns the command string to run raxml

    Functions:
    ----------
    run()
        submits a raxml job to locally or on an ipyparallel client cluster. 

    """    

    ## init object for params
    def __init__(self,
        data,
        name="test",
        workdir="analysis-raxml", 
        *args, 
        **kwargs):

        ## path attributes
        self._kwargs = {
            "f": "a", 
            "T": 4,  # <- change to zero !?
            "m": "GTRGAMMA",
            "N": 100,
            "x": 12345,
            "p": 54321,
            "o": None,
            "binary": "",
            }
        self._kwargs.update(kwargs)

        ## check workdir
        if workdir:
            workdir = os.path.abspath(os.path.expanduser(workdir))
        else:
            workdir = os.path.abspath(os.path.curdir)
        if not os.path.exists(workdir):
            os.makedirs(workdir)

        ## entered args
        self.params = Params()
        self.params.n = name
        self.params.w = workdir
        self.params.s = os.path.abspath(os.path.expanduser(data))

        ## find the binary
        if not self._kwargs["binary"]:
            self.params.binary = _find_binary()

        ## set params
        notparams = set(["workdir", "name", "data", "binary"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## check binary
        self._get_binary()

        ## attributesx
        self.rasync = None
        self.stdout = None
        self.stderr = None

        ## results files        
        self.trees = Params()
        self.trees.bestTree = OPJ(workdir, "RAxML_bestTree."+name)
        self.trees.bipartitionsBranchLabels = OPJ(workdir, "RAxML_bipartitionsBranchLabels."+name)
        self.trees.bipartitions = OPJ(workdir, "RAxML_bipartitions."+name)
        self.trees.bootstrap = OPJ(workdir, "RAxML_bootstrap."+name)
        self.trees.info = OPJ(workdir, "RAxML_info."+name)


    @property
    def _command_list(self):
        """ build the command list """
        cmd = [self.params.binary, 
                "-f", str(self.params.f), 
                "-T", str(self.params.T), 
                "-m", str(self.params.m),
                "-N", str(self.params.N),
                "-x", str(self.params.x),
                "-p", str(self.params.p),
                "-n", str(self.params.n),
                "-w", str(self.params.w),
                "-s", str(self.params.s),
               ]
        ## add ougroups
        if self.params.o:
            cmd += ["-o"]
            cmd += [",".join(self.params.o)]
        return cmd


    @property
    def command(self):
        """ returns command as a string """
        return " ".join(self._command_list)


    def run(self, 
        ipyclient=None, 
        quiet=False,
        force=False,
        block=False,
        ):
        """
        Submits raxml job to run. If no ipyclient object is provided then 
        the function will block until the raxml run is finished. If an ipyclient
        is provided then the job is sent to a remote engine and an asynchronous 
        result object is returned which can be queried or awaited until it finishes.

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

        ## stop before trying in raxml
        if force:
            for key, oldfile in self.trees:
                if os.path.exists(oldfile):
                    os.remove(oldfile)
        if os.path.exists(self.trees.info):
            print("Error: set a new name for this job or use Force flag.\nFile exists: {}"\
                  .format(self.trees.info))
            return 

        ## TODO: add a progress bar tracker here. It could even read it from
        ## the info file that is being written. 
        ## submit it
        if not ipyclient:
            proc = _call_raxml(self._command_list)
            self.stdout = proc[0]
            self.stderr = proc[1]
        else:
            ## find all hosts and submit job to the host with most available engines
            lbview = ipyclient.load_balanced_view()
            self.rasync = lbview.apply(_call_raxml, self._command_list)

        ## initiate random seed
        if not quiet:
            if not ipyclient:
                ## look for errors
                if "Overall execution time" not in self.stdout.decode():
                    print("Error in raxml run\n" + self.stdout.decode())
                else:
                    print("job {} finished successfully".format(self.params.n))
            else:
                print("job {} submitted to cluster".format(self.params.n))



    def _get_binary(self):
        """ find binaries available"""

        ## check for binary
        backup_binaries = ["raxmlHPC-PTHREADS", "raxmlHPC-PTHREADS-SSE3"]

        ## check user binary first, then backups
        for binary in [self.params.binary] + backup_binaries:
            proc = subprocess.Popen(["which", self.params.binary], 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT).communicate()
            ## update the binary
            if proc:
                self.params.binary = binary

        ## if none then raise error
        if not proc[0]:
            raise Exception(BINARY_ERROR.format(self.params.binary))



def _find_binary():
    ## check for binary
    list_binaries = [
        "raxmlHPC-PTHREADS-AVX2",
        "raxmlHPC-PTHREADS-AVX",            
        "raxmlHPC-PTHREADS-SSE3",
        "raxmlHPC-PTHREADS", 
        ]

    ## check user binary first, then backups
    for binary in list_binaries:
        proc = subprocess.Popen(["which", binary],
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT).communicate()
        ## update the binary
        if proc[0]:
            return binary

    ## if not binaries found
    raise Exception("cannot find raxml binary")


def _call_raxml(command_list):
    """ call the command as sps """
    proc = subprocess.Popen(
        command_list,
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE
        )
    comm = proc.communicate()
    return comm



BINARY_ERROR = """
  Binary {} not found. 

  Check that you have raxml installed. If you have a different binary
  installed you can select it using the argument 'binary'. 

  For example, 
    rax = ipa.raxml(name='test', data='test.phy', binary='raxmlHPC')

  or, you can set it after object creation with:
    rax.params.binary = "raxmlHPC-PTHREADS"

  If you install raxml with conda ipyrad should be able to find it:
  'conda install raxml -c bioconda'
"""
