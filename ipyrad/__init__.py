#!/usr/bin/env ipython2

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

## define state vars
__version__ = "0.1.10"
__interactive__ = 1
__loglevel__ = "DEBUG"   ##  "DEBUG"  "INFO"  "WARN"  "ERROR" 


## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging as _logging
import logging.config as _lconfig

## clear the log file 
#with open("/tmp/ipyrad_debug.txt", 'w') as logfile: 
#    pass

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
            'filename':'/tmp/ipyrad_debug.txt',
            'formatter':"standard",
            'mode':'a'
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

_LOG = _logging.getLogger(__name__)
if __loglevel__ == "DEBUG":
    _LOG.debug("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
else:
    _LOG.info("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)


if __name__ == "__main__":
    __interactive__ = 0

