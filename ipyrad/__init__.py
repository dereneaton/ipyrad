#!/usr/bin/env ipython2


## define state vars
__interactive__ = 1      ## CLI __main__ changes to 0
__version__ = "0.1.42"

## Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
__loglevel__ = "ERROR"
__debugfile__ = "/tmp/ipyrad_debug.txt"

## main ip.functions
from . import load
from . import assemble 
#from . import plotting  ## do not autoimport plotting, import as ipp
#from . import analysis  ## do not autoimport analysis, import as ipa

## bring nested functions to ip.
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import paramsinfo


## create logger for debugging
## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging as _logging
import logging.config as _lconfig

## clear the logfile if it is too big
import os
try:
    if os.path.getsize(__debugfile__) > 5000000:
        with open(__debugfile__, 'w') as clear:
            clear.write("")
## in case system doesn't let you use /tmp            
except IOError:
    __debugfile__ = os.devnull
__loglevel__ = "ERROR"

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
LOGGER = _logging.getLogger(__name__)

## set globals
if __loglevel__ == "DEBUG":
    LOGGER.debug("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
else:
    LOGGER.info("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)


if __name__ == "__main__":
    pass
