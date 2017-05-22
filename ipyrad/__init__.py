#!/usr/bin/env python2


## dunders
__version__ = "0.6.22"
__author__ = "Deren Eaton & Isaac Overcast"

## Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
__debugflag__ = "./.debug"
__debugfile__ = "./ipyrad_log.txt"

## define state vars
__interactive__ = 1      ## CLI __main__ changes to 0

# pylint: disable=C0103
import os as _os
import atexit as _atexit

## main ip.functions
from . import load as _load
from . import assemble 
from ipyrad.core.parallel import get_client as _get_client

#from . import file_conversion     ## do not autoimport; import as fc
#from . import plotting            ## do not autoimport; import as ipp
#from . import analysis            ## do not autoimport; import as ipa

## bring nested functions to top for API
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from .load import save_json
from .load import load_json

####################################################################
## create logger for debugging
## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging as _logging
import logging.config as _lconfig

## debug is set based on whether the flag exists
if _os.path.exists(__debugflag__):
    __loglevel__ = "DEBUG"
else:
    __loglevel__ = "ERROR"#"INFO"

## ensures MPI is loaded at init on all engines
# actually it freezes if we use this...
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()

## check that all dependencies exist and are working
import subprocess as _subprocess
import sys as _sys
import socket as _socket


_LOGGER = _logging.getLogger(__name__)
if __loglevel__ == "DEBUG":
    _LOGGER.debug("Engine init")



def cluster_info(
    client=None, 
    cluster_id="", 
    profile="default",
    engines="Local", 
    timeout=60, 
    cores=0, 
    quiet=True,
    spacer=0,
    **kwargs):
    """ report info on the ipcluster instance """  

    ## get the client
    if not client:
        ipyclient = _get_client(cluster_id, profile, 
                                engines, timeout, cores, quiet)
    else:
        ipyclient = client

    ## report 
    hosts = ipyclient[:].apply_sync(_socket.gethostname)
    result = []
    if spacer:
        spacer = " " * spacer
    elif not __interactive__:
        spacer = "  "
    else:
        spacer = ""
    for hostname in set(hosts):
        result.append(spacer+"host compute node: [{} cores] on {}"\
                      .format(hosts.count(hostname), hostname))
    
    if not client:
        ipyclient.close()
    return "\n".join(result)



def _debug_on():
    """
    Turns on debugging by creating hidden tmp file
    This is only run by the __main__ engine.
    """
    ## make tmp file and set loglevel for top-level init
    with open(__debugflag__, 'w') as dfile:
        dfile.write("wat")
    __loglevel__ = "DEBUG"
    _LOGGER.info("debugging turned on and registered to be turned off at exit")
    _set_debug_dict(__loglevel__)



def _set_debug_dict(__loglevel__):
    """ set the debug dict """

    _lconfig.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,

    'formatters': {
        'standard': {
            'format': "%(asctime)s \t"\
                     +"pid=%(process)d \t"\
                     +"[%(filename)s]\t"\
                     +"%(levelname)s \t"\
                     +"%(message)s"
        },
    },
    'handlers': {
        __name__: {
            'level':__loglevel__,
            'class':'logging.FileHandler',
            'filename':__debugfile__,
            'formatter':"standard",
            'mode':'a+'
        }
    },
    'loggers':{
        __name__: {
            'handlers': [__name__],
            'level': __loglevel__,
            'propogate': True
        }
    }
    })

_set_debug_dict(__loglevel__)


def _debug_off():
    """ turns off debugging by removing hidden tmp file """
    if _os.path.exists(__debugflag__):
        _os.remove(__debugflag__)



def _cmd_exists(cmd):
    """ check if dependency program is there """
    return _subprocess.call("type " + cmd,
                           shell=True,
                           stdout=_subprocess.PIPE,
                           stderr=_subprocess.PIPE) == 0



def _getbins():
    """ gets the right version of vsearch, muscle, and smalt
    depending on linux vs osx """

    # Return error if system is 32-bit arch.
    # This is straight from the python docs:
    # https://docs.python.org/2/library/platform.html#cross-platform
    if not _sys.maxsize > 2**32:
        _sys.exit("ipyrad requires 64bit architecture")

    ## get platform mac or linux
    _platform = _sys.platform

    ## get current location
    path = _os.path.abspath(_os.path.dirname(__file__))

    ## find bin directory
    ipyrad_path = _os.path.dirname(path)
    bin_path = _os.path.join(ipyrad_path, "bin")

    ## get the correct binaries
    if 'linux' in _platform:
        vsearch = _os.path.join(
                       _os.path.abspath(bin_path),
                       "vsearch-linux-x86_64")
        muscle = _os.path.join(
                       _os.path.abspath(bin_path),
                       "muscle-linux-x86_64")
        smalt = _os.path.join(
                       _os.path.abspath(bin_path),
                       "smalt-linux-x86_64")
        bwa = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bwa-linux-x86_64")
        samtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "samtools-linux-x86_64")
        bedtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bedtools-linux-x86_64")
        qmc = _os.path.join(
                       _os.path.abspath(bin_path),
                       "QMC-linux-x86_64")
    else:
        vsearch = _os.path.join(
                       _os.path.abspath(bin_path),
                       "vsearch-osx-x86_64")
        muscle = _os.path.join(
                       _os.path.abspath(bin_path),
                       "muscle-osx-x86_64")
        smalt = _os.path.join(
                       _os.path.abspath(bin_path),
                       "smalt-osx-x86_64")
        bwa = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bwa-osx-x86_64")
        samtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "samtools-osx-x86_64")
        bedtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bedtools-osx-x86_64")
        ## only one compiled version available, works for all?
        qmc = _os.path.join(
                       _os.path.abspath(bin_path),
                       "QMC-osx-x86_64")

    # Test for existence of binaries
    assert _cmd_exists(muscle), "muscle not found here: "+muscle
    assert _cmd_exists(vsearch), "vsearch not found here: "+vsearch
    assert _cmd_exists(smalt), "smalt not found here: "+smalt
    assert _cmd_exists(bwa), "bwa not found here: "+bwa
    assert _cmd_exists(samtools), "samtools not found here: "+samtools
    assert _cmd_exists(bedtools), "bedtools not found here: "+bedtools
    #assert _cmd_exists(qmc), "wQMC not found here: "+qmc
    return vsearch, muscle, smalt, bwa, samtools, bedtools, qmc


## create globals for binaries that can be accessed as: ipyrad.bins.muscle
bins = assemble.util.ObjDict()
_binnames = ["vsearch", "muscle", "smalt", "bwa", "samtools", "bedtools", "qmc"]
for _binn, _binx in zip(_binnames, _getbins()):
    bins[_binn] = _binx
## clean up for the API
del _binnames



if __name__ == "__main__":
    pass
