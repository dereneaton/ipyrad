#!/usr/bin/env ipython2

## main ip.functions
from . import assemble
from . import plotting
from . import dstats


## bring nested functions to ip.
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import get_params_info
from ipyrad.core.load_dataobj import load_assembly


## define state vars
__version__ = "0.0.66"
__interactive__ = 1
__loglevel__ = "DEBUG"   ##  "DEBUG"  "INFO"  "WARN"  "ERROR" 


## failed attempts at launching parallel code in __init__
#import ipyrad.core.parallel


## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging
import logging.config

## clear the log file 
#with open("/tmp/ipyrad_debug.txt", 'w') as logfile: 
#    pass

logging.config.dictConfig({
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

LOGGER = logging.getLogger(__name__)
if __loglevel__ == "DEBUG":
    LOGGER.debug("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
else:
    LOGGER.info("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)


