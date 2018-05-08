#!/usr/bin/env python

import os
import logging.config


def _set_logger_config():

    # debug is set based on whether the flag exists
    # Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
    __debugflag__ = "./.debug"
    __debugfile__ = "./ipyrad_log.txt"

    # check for debugflag file
    if os.path.exists(__debugflag__):
        __loglevel__ = "DEBUG"
    else:
        __loglevel__ = "ERROR"

    # set the debug dict
    logging.config.dictConfig({
        'version': 1,
        'disable_existing_loggers': False,

        'formatters': {
            'standard': {
                'format': "\t".join((
                    "%(asctime)s", 
                    "pid=%(process)d", 
                    "[%(filename)s]", 
                    "%(levelname)s", 
                    "%(message)s"))
            },
        },
        'handlers': {
            __name__: {
                'level': __loglevel__,
                'class': 'logging.FileHandler',
                'filename': __debugfile__,
                'formatter': "standard",
                'mode': 'a+'
            }
        },
        'loggers': {
            __name__: {
                'handlers': [__name__],
                'level': __loglevel__,
                'propogate': True
            }
        }
    })



# logging.basicConfig(
#     filename="./ipyrad_log.txt",
#     filemode="a+",
#     level=logging.DEBUG,
#     format="\t".join((
#         "%(asctime)s", 
#         "pid=%(process)d", 
#         "[%(filename)s]", 
#         "%(levelname)s", 
#         "%(message)s"))
#     )


# # debug is set based on whether the flag exists
# # Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
# __debugflag__ = "./.debug"
# __debugfile__ = "./ipyrad_log.txt"


# if os.path.exists(__debugflag__):
#     __loglevel__ = "DEBUG"
# else:
#     __loglevel__ = "ERROR"


# _LOGGER = logging.getLogger(__name__)
# if __loglevel__ == "DEBUG":
#     _LOGGER.debug("Engine init")



# class Logger():
#     def __init__(self):
#         pass


# def _debug_on():
#     """
#     Turns on debugging by creating hidden tmp file
#     This is only run by the __main__ engine.
#     """
#     ## make tmp file and set loglevel for top-level init
#     with open(__debugflag__, 'w') as dfile:
#         dfile.write("wat")
#     __loglevel__ = "DEBUG"
#     _LOGGER.info("debugging turned on and registered to be turned off at exit")
#     _set_debug_dict(__loglevel__)




# def _set_debug_dict(__loglevel__):
#     """ set the debug dict """

#     config.dictConfig({
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

# _set_debug_dict(__loglevel__)




# def _debug_off():
#     """ turns off debugging by removing hidden tmp file """
#     if _os.path.exists(__debugflag__):
#         _os.remove(__debugflag__)
#     __loglevel__ = "ERROR"
#     _LOGGER.info("debugging turned off")
#     _set_debug_dict(__loglevel__)
  
