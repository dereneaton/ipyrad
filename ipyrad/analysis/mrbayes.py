#!/usr/bin/python 

""" wrapper to make simple calls to mb """

import os
import sys
import subprocess
import pandas as pd
from ipyrad.assemble.utils import Params
from ipyrad.assemble.utils import IPyradError

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
mcmcp filename={outname};
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



# template for clock model tree inference
# https://www.groundai.com/project/molecular-clock-dating-using-mrbayes/
NEX_TEMPLATE_3 = """\
#NEXUS

log start filename={outname}.log replace;
execute {nexus};

begin mrbayes;
set autoclose=yes nowarn=yes;

lset nst=6 rates=gamma;

prset brlenspr=clock:uniform;
prset clockvarpr=igr;
prset igrvarpr=exp(10.0);
prset clockratepr=normal(0.01,0.005);

{other}

mcmcp ngen={ngen} nrun={nruns} nchains={nchains};
mcmcp relburnin=yes burninfrac=0.25;
mcmcp samplefreq={samplefreq};
mcmcp printfreq=10000 diagnfr=5000;
mcmcp filename={outname};
mcmc;

sump filename={outname};
sumt filename={outname} contype=allcompat;
log stop filename={outname}.log append;
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
    def __init__(
        self,
        data,
        name="test",
        workdir="analysis-mb", 
        clock_model=False,        
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
        self.clock_model = clock_model
        self.name = name
        self.workdir = workdir
        self.data = os.path.abspath(os.path.expanduser(data))
        self.nexus = os.path.join(self.workdir, self.name + ".nex")
        self.binary = ""
        self._get_binary(self._kwargs.get("binary"))

        self.params = Params()
        defaults = {
            "clockratepr": "lognorm(-7,0.6)",
            "clockvarpr": "tk02",
            "tk02varpr": "exp(1.0)",
            "brlenspr": "clock:birthdeath",
            "samplestrat": "diversity",
            "sampleprob": "0.1",
            "speciationpr": "exp(10)",
            "extinctionpr": "beta(2, 200)",
            "treeagepr": "offsetexp(1, 5)",
            "ngen": 100000,
            "nruns": "1",
            "nchains": 4,
            "samplefreq": 1000,
        }
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
        self.trees.constre = os.path.join(
            self.workdir, "{}.nex.con.tre".format(self.name))
        self.trees.posttre = os.path.join(
            self.workdir, "{}.nex.{}t".format(self.name, runs))
        self.trees.pstat = os.path.join(
            self.workdir, "{}.nex.pstat".format(self.name))
        self._write_nexus_file()

        # check attribute for existing results at this name.
        if self.result_files:
            print(
                "Existing results loaded for run [{}], see .trees attribute."
                .format(len(self.result_files), self.name)
            )


    def _write_nexus_file(self, write=True):
        "write a mrbayes block to a copy of the NEXUS file."
        cwargs = self.params.__dict__.copy()
        cwargs["nexus"] = self.data
        cwargs["outname"] = self.nexus
        cwargs["other"] = ""
        if self.clock_model == 1:
            self._nexstring = NEX_TEMPLATE_2.format(**cwargs)
        elif self.clock_model == 2:
            self._nexstring = NEX_TEMPLATE_3.format(**cwargs)
        else:
            self._nexstring = NEX_TEMPLATE_1.format(**cwargs)
        with open(self.nexus, 'w') as out:
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
                "run 'conda install mrbayes -c bioconda -c conda-forge'")



def _call_mb(command_list):
    """ call the command as sps """
    proc = subprocess.Popen(
        command_list,
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE
        )
    comm = proc.communicate()
    if proc.returncode:
        raise IPyradError(comm[0].decode())
    return comm[0].decode()
