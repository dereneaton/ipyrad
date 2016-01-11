#!/usr/bin/env ipython2
#!/usr/bin/env ipython2


## main ip.functions
## main ip.functions
from . import load
from . import load
from . import assemble 
from . import assemble 
#from . import plotting  ## do not autoimport plotting, import as ipp
#from . import plotting  ## do not autoimport plotting, import as ipp
#from . import analysis  ## do not autoimport analysis, import as ipa
#from . import analysis  ## do not autoimport analysis, import as ipa


## bring nested functions to ip.
## bring nested functions to ip.
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import Assembly
from ipyrad.core.assembly import merge
from ipyrad.core.assembly import merge
from ipyrad.core.sample import Sample
from ipyrad.core.sample import Sample
from ipyrad.core.paramsinfo import paramsinfo
from ipyrad.core.paramsinfo import paramsinfo


## define state vars
## define state vars
__version__ = "0.1.14"
__version__ = "0.1.14"
__interactive__ = 1
__interactive__ = 1
__loglevel__ = "DEBUG"   ##  "DEBUG"  "INFO"  "WARN"  "ERROR" 
__loglevel__ = "ERROR"




## this needs to come after __loglevel__ definition
__loglevel__ = "ERROR"
## sets log config and prints warning if __loglevel__ is in hackers mode
__loglevel__ = "ERROR"
import logging as _logging
import logging as _logging
import logging.config as _lconfig
import logging.config as _lconfig


## clear the log file 
## clear the log file 
#with open("/tmp/ipyrad_debug.txt", 'w') as logfile: 
#with open("/tmp/ipyrad_debug.txt", 'w') as logfile: 
#    pass
#    pass


_lconfig.dictConfig({
_lconfig.dictConfig({
    'version': 1,              
    'version': 1,              
    'disable_existing_loggers': False,  
    'disable_existing_loggers': False,  


    'formatters': {
    'formatters': {
        'standard': {
        'standard': {
            'format': "%(asctime)s \t"\
            'format': "%(asctime)s \t"\
                     +"pid=%(process)d \t"\
                     +"pid=%(process)d \t"\
                     +"[%(filename)s]\t"\
                     +"[%(filename)s]\t"\
                     +"%(levelname)s \t"\
                     +"%(levelname)s \t"\
                     +"%(message)s"
                     +"%(message)s"
        },
        },
    },
    },
    'handlers': {
    'handlers': {
        __name__: {
        __name__: {
            'level':__loglevel__,
__loglevel__ = "ERROR"
            'class':'logging.FileHandler',
            'class':'logging.FileHandler',
            'filename':'/tmp/ipyrad_debug.txt',
            'filename':'/tmp/ipyrad_debug.txt',
            'formatter':"standard",
            'formatter':"standard",
            'mode':'a'
            'mode':'a'
        }
        }
    },
    },
    'loggers':{
    'loggers':{
        __name__: {
        __name__: {
            'handlers': [__name__],
            'handlers': [__name__],
            'level': __loglevel__,
__loglevel__ = "ERROR"
            'propogate': True
            'propogate': True
        }
        }
    }
    }
})
})


_LOG = _logging.getLogger(__name__)
_LOG = _logging.getLogger(__name__)
if __loglevel__ == "DEBUG":
__loglevel__ = "ERROR"
    _LOG.debug("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
__loglevel__ = "ERROR"
else:
else:
    _LOG.info("H4CKERZ-mode: __loglevel__ = %s", __loglevel__)
__loglevel__ = "ERROR"




if __name__ == "__main__":
if __name__ == "__main__":
    __interactive__ = 0
    __interactive__ = 0


