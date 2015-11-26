#!/usr/bin/env ipython2

""" A logger that works across engines and writes to /tmp/ipyrad_debug.txt """


import logging
import logging.config
from ipyrad import __loglevel__


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


## print init logging that should only show up during development
