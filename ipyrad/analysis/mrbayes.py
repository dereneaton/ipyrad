#!/usr/bin/python 

""" wrapper to make simple calls to mb """

import os
import subprocess
from ipyrad.assemble.utils import Params
from ipyrad.assemble.utils import IPyradError

OPJ = os.path.join


# template for standard tree inference
NEX_TEMPLATE_1 = """\
#NEXUS
execute {nexus};

begin mrbayes;
set autoclose=yes nowarn=yes;

lset nst=6 rates=gamma;

mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
mcmcp relburnin=yes burninfrac=0.25; 
mcmcp samplefreq={samplefreq} printfreq=10000;
mcmcp filename={outname}
mcmc;

sump filename={outname};
sumt filename={outname};
end;
"""


# template for clock model tree inference
# https://www.groundai.com/project/molecular-clock-dating-using-mrbayes/
NEX_TEMPLATE_2 = """\
#NEXUS
execute {nexus};

begin mrbayes;
set autoclose=yes nowarn=yes;

lset nst=6 rates=gamma;

prset clockratepr=lognorm(-7,0.6);
prset clockvarpr=tk02;
prset tk02varpr=exp(1.0);
prset brlenspr=clock:birthdeath;
prset samplestrat=diversity;
prset sampleprob=0.1;
prset speciationpr=exp(10);
prset extinctionpr=beta(2, 200);
prset treeagepr=offsetexp(1,5);

mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
mcmcp relburnin=yes burninfrac=0.25;
mcmcp samplefreq={samplefreq};
mcmcp printfreq=10000 diagnfr=5000;
mcmcp filename={outname};
mcmc;

sump filename={outname};
sumt filename={outname};
end;
"""



class MrBayes(object):
    """
    MrBayes analysis utility function for running simple commands. 

    Parameters:
    -----------
    data: str
        The phylip formated sequence file (.phy from ipyrad).
    name: str
        The name for this run. An alias for '-n'.
    workdir: str
        The output directory for results. An alias for '-w'. 

    Additional optional parameters
    -------------------------------
    ngen: int
        Number of MCMC generations to run.
    sfreq: int
        Frequency to sample from MCMC chain.
    burnin: int
        Number of generations to run before sampling starts.       

    Attributes:
    -----------
    params: dict
        parameters for this mb run
    cmd: 
        returns the command string to run mb

    Functions:
    ----------
    run()
        submits a mrbayes job locally or on an ipyparallel client cluster. 
    """    

    # init object for params
    def __init__(self,
        data,
        name="test",
        workdir="analysis-mb", 
        clock=False,        
        params={
            "clockratepr": "lognorm(-7,0.6)",
            "clockvarpr": "tk02",
            "tk02varpr": "exp(1.0)",
            "brlenspr": "clock:birthdeath",
            "samplestrat": "diversity",
            "sampleprob": "0.1",
            "speciationpr": "exp(10)",
            "extinctionpr": "beta(2, 200)",
            "treeagepr": "offsetexp(1, 5)",
            "ngen": "1000000",
            "nruns": "2",
            "nchains": 4,
            "samplefreq": 1000,
        },                   
        **kwargs):

        # path attributes
        self._kwargs = {}            
        self._kwargs.update(kwargs)

        # check workdir
        if workdir:
            workdir = os.path.abspath(os.path.expanduser(workdir))
        else:
            workdir = os.path.abspath(os.path.curdir)
        if not os.path.exists(workdir):
            os.makedirs(workdir)

        # entered args
        self.params = Params()
        self.params.name = name
        self.params.workdir = workdir
        self.params.data = os.path.abspath(os.path.expanduser(data))
        self.params.nexus = os.path.join(
            self.params.workdir, self.params.name + ".nex")
        self.params.binary = ""

        # sets self.params.binary
        self._get_binary(self._kwargs.get("binary"))

        # set params
        notparams = set(["workdir", "name", "data", "binary"])
        for key in set(self._kwargs.keys()) - notparams:
            self.params[key] = self._kwargs[key]

        ## attributesx
        self.rasync = None
        self.stdout = None
        self.stderr = None

        ## results files        
        self.trees = Params()
        self.trees.constre = self.params.data + ".con.tre"
        self.trees.info = self.params.data + ".lstat"

        # write a mrbayes block to a copy of the NEXUS file.
        cwargs = params
        cwargs.update({"nexus": self.params.data, "outname": self.params.nexus})
        if clock:
            self._nexstring = NEX_TEMPLATE_2.format(**cwargs)
        else:
            self._nexstring = NEX_TEMPLATE_1.format(**cwargs)
        with open(self.params.nexus, 'w') as out:
            out.write(self._nexstring)



    def print_command(self):
        print("mb {}".format(self.params.nexus))
    

    def print_nexus_string(self):
        print(self._nexstring)



    def run(self, 
        ipyclient=None, 
        quiet=False,
        force=False,
        block=False,
        ):
        """
        Submits mrbayes job to run. If no ipyclient object is provided then 
        the function will block until the mb run is finished. If an ipyclient
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

        # stop before trying in mrbayes
        if force:
            for key, oldfile in self.trees:
                if os.path.exists(oldfile):
                    os.remove(oldfile)
        if os.path.exists(self.trees.info):
            print("Error Files Exist: set a new name or use Force flag.\n{}"
                  .format(self.trees.info))
            return 

        ## TODO: add a progress bar tracker here. It could even read it from
        ## the info file that is being written. 
        ## submit it
        if not ipyclient:
            proc = _call_mb([self.params.binary, self.params.nexus])
            self.stdout = proc[0]
            self.stderr = proc[1]
        else:
            ## find all hosts and submit job to the host with most available engines
            lbview = ipyclient.load_balanced_view()
            self.rasync = lbview.apply(
                _call_mb, [self.params.binary, self.params.nexus])

        ## initiate random seed
        if not quiet:
            if not ipyclient:
                ## look for errors
                if "Overall execution time" not in self.stdout.decode():
                    print("Error in mb run\n" + self.stdout.decode())
                else:
                    print("job {} finished successfully".format(self.params.n))
            else:               
                if block:
                    print("job {} running".format(self.params.n))
                    ipyclient.wait()
                    if self.rasync.successful():
                        print("job {} finished successfully"
                            .format(self.params.n))
                    else:
                        raise IPyradError(self.rasync.get())
                else:
                    print("job {} submitted to cluster".format(self.params.n))



    def _get_binary(self, search_first=None):
        """ find binaries available"""

        # check for binary
        list_binaries = [
            search_first,
            "mb",
        ]
        list_binaries = [i for i in list_binaries if i]

        # check user binary first, then backups
        for binary in list_binaries:
            proc = subprocess.Popen(["which", binary],
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT).communicate()
            # update the binary
            if proc[0]:
                self.params.binary = binary

        ## if none then raise error
        if not proc[0]:
            raise Exception("cannot find mb; run 'conda install mrbayes -c bioconda'")



def _call_mb(command_list):
    """ call the command as sps """
    proc = subprocess.Popen(
        command_list,
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE
        )
    comm = proc.communicate()
    return comm
