#!/usr/bin/env ipython2

# pylint: disable=C0103

## define state vars
__interactive__ = 1      ## CLI __main__ changes to 0
__version__ = "0.1.45"

## Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
__loglevel__ = "ERROR"
__debugfile__ = "/tmp/ipyrad_debug.txt"

## main ip.functions
from . import load
from . import assemble 
from .load import save_json
from .load import load_json
#from . import plotting  ## do not autoimport plotting, import as ipp
#from . import analysis  ## do not autoimport analysis, import as ipa

## bring nested functions to ip.
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import paramsinfo


####################################################################
## create logger for debugging
## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging as _logging
import logging.config as _lconfig

## clear the logfile if it is too big
import os as _os
try:
    if _os.path.getsize(__debugfile__) > 5000000:
        with open(__debugfile__, 'w') as clear:
            clear.write("")
## in case system doesn't let you use /tmp            
except (OSError, IOError):
    __debugfile__ = _os.devnull
    _, __loglevel__ = "null", "ERROR"  ## hack for versioner


# set up logging to file 
_logging.basicConfig(level=__loglevel__,
                     format="%(asctime)s \t"\
                            +"pid=%(process)d \t"\
                            +"[%(filename)s]\t"\
                            +"%(levelname)s \t"\
                            +"%(message)s",
                     datefmt='%m-%d %H:%M',
                     filename='/tmp/ipyrad_debug.txt',
                     filemode='w')

# Define the logger and test
_LOGGER = _logging.getLogger(__name__)

## set globals
if __loglevel__ == "DEBUG":
    _LOGGER.debug("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
else:
    _LOGGER.info("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)

####################################################################

## check that all dependencies exist and are working
import subprocess as _subprocess
import sys as _sys

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
        samtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "samtools-linux-x86_64")
        bedtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bedtools-linux-x86_64")
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
        samtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "samtools-osx-x86_64")
        bedtools = _os.path.join(
                       _os.path.abspath(bin_path),
                       "bedtools-osx-x86_64")

    # Test for existence of binaries
    assert _cmd_exists(muscle), "muscle not found here: "+muscle
    assert _cmd_exists(vsearch), "vsearch not found here: "+vsearch
    assert _cmd_exists(smalt), "smalt not found here: "+smalt
    assert _cmd_exists(samtools), "samtools not found here: "+samtools
    assert _cmd_exists(bedtools), "bedtools not found here: "+bedtools
    return vsearch, muscle, smalt, samtools, bedtools


## create globals for binaries that can be accessed as: ipyrad.bins.muscle
bins = assemble.util.ObjDict()
binnames = ["vsearch", "muscle", "smalt", "samtools", "bedtools"]
for binn, binx in zip(binnames, _getbins()):
    bins[binn] = binx
## clean up for the API
del binnames, binn, binx



if __name__ == "__main__":
    pass
