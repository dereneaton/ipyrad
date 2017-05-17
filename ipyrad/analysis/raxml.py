#!/usr/bin/python 

""" wrapper to make simple calls to raxml """

import os
import subprocess


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
    submit_raxml_job()
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
            "T": 4,
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
            workdir = os.path.curdir        
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
        notparams = set(["workdir", "name", "data"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## attributes
        self.stdout = None
        self.stderr = None

        ## results files        
        self.trees = Params()
        self.trees.bestTree = None
        self.trees.bipartitionsBranchLabels = None
        self.trees.bipartitions = None
        self.trees.boostrap = None
        self.trees.info = None


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
        ):
        """
        Submits raxml job to run on the cluster. 

        Parameters
        -----------
        ipyclient:
            Not yet supported... 
        quiet: 
            suppress print statements
        force:
            remove existing results files with this job name. 
        """



        ## if none then raise error
        if not proc[0]:
            raise Exception(BINARY_ERROR.format(self.params.binary))

        ## attach tree paths to the results
        f1 = os.path.join(self.params.w, "RAxML_bestTree."+self.params.n)
        f2 = os.path.join(self.params.w, "RAxML_bipartitionsBranchLabels."+self.params.n)
        f3 = os.path.join(self.params.w, "RAxML_bipartitions."+self.params.n)
        f4 = os.path.join(self.params.w, "RAxML_bootstrap."+self.params.n)
        f5 = os.path.join(self.params.w, "RAxML_info."+self.params.n)

        ## stop before trying in raxml
        if force:
            for oldfile in [f1, f2, f3, f4, f5]:
                if os.path.exists(oldfile):
                    os.remove(oldfile)

        if os.path.exists(f5):
            print("Error: set a new name for this job:\nFile exists: {}".format(f5))
            return 

        ## TODO: add a progress bar tracker here. It could even read it from
        ## the info file that is being written. 
        ## submit it
        if not ipyclient:
            proc = _call_raxml(self._command_list)
            self.stdout = proc[0]
            self.stderr = proc[1]
        else:
            print("not yet supported")
            return 
            sync = ipyclient[0].apply(_call_raxml, self._command_list)
            sync.wait()

        ## TODO: error checking
        ## ...

        if os.path.exists(f1):
            self.trees.bestTree = f1
        if os.path.exists(f2):
            self.trees.bipartitionsBranchLabels = f2
        if os.path.exists(f3):
            self.trees.bipartitions = f3
        if os.path.exists(f4):
            self.trees.boostrap = f4
        if os.path.exists(f5):
            self.trees.info = f5

        ## initiate random seed
        if not quiet:
            print("job {} finished successfully".format(self.params.n))



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
        if proc:
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



class Params(object):
    """ 
    A dict-like object for storing params values with a custom repr
    """
    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __repr__(self):
        _repr = ""
        keys = sorted(self.__dict__.keys())      
        _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
        for key in keys:
            _val = str(self[key]).replace(os.path.expanduser("~"), "~")
            _repr += _printstr.format(key, _val)
        return _repr


BINARY_ERROR = """
  Binary {} not found. 

  Check that you have raxml installed. If you have a different binary
  installed you can select it using the argument 'binary'. 

  For example, 
    rax = ipa.raxml(name='test', data='test.phy', binary='raxmlHPC')

  or, you can set it after object creation with:
    rax.params.binary = "raxmlHPC-PTHREADS"
"""

  