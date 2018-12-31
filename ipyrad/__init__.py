#!/usr/bin/env python

## bring nested functions to top for API access
from .core.assembly import Assembly, merge
#from .core.sample import Sample
from .core.parallel import get_client as _get_client, cluster_info
from .core.startup import Bins as _Bins
from .load import save_json, load_json
import logging as _logging


# Dunders
__version__ = "0.8.0-dev"
__author__ = "Deren Eaton & Isaac Overcast"

# CLI __main__ changes to 0
__interactive__ = 1



# def set_logger_level(level=None):
#     #import os as _os
#     import logging.config as _lconfig

#     # set the log level
#     if level:
#         assert level in [0, 1, 2, 3, "DEBUG", "INFO", "WARN", "ERROR"]
#         __loglevel__ = level
#     else:
#         __loglevel__ = "ERROR"
#         # if _os.path.exists(__debugflag__):
#         #     __loglevel__ = "DEBUG"
#         # else:
#         #     __loglevel__ = "ERROR"

#     # set the debug dict
#     _lconfig.dictConfig({
#         'version': 1,
#         'disable_existing_loggers': False,

#         'formatters': {
#             'standard': {
#                 'format': "\t".join((
#                     "%(asctime)s", 
#                     "pid=%(process)d", 
#                     "[%(filename)s]", 
#                     "%(levelname)s", 
#                     "%(message)s"))
#             },
#         },
#         'handlers': {
#             __name__: {
#                 'level': __loglevel__,
#                 'class': 'logging.FileHandler',
#                 'filename': __debugfile__,
#                 'formatter': "standard",
#                 'mode': 'a+'
#             }
#         },
#         'loggers': {
#             __name__: {
#                 'handlers': [__name__],
#                 'level': __loglevel__,
#                 'propogate': True
#             }
#         }
#     })


# log file
__debugfile__ = "./ipyrad_log.txt"
__debugflag__ = "./.debug"
# logger = _logging.getLogger(__name__)
# set_logger_level()

# get binaries
bins = _Bins()
