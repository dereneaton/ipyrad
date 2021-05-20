#!/usr/bin/env python

"""
Logger for stderr and optionally to file.
"""

import sys
from loguru import logger


STDFORMAT = (
    "<cyan>{time:hh:mm:ss}</cyan> <white>|</white> "
    "<level>{level: <7}</level> <white>|</white> "
    "<magenta>{file:<18}</magenta> <white>|</white> "
    "<level>{message}</level>"
)

LOGFORMAT = (
    "{time:YYYY/MM/DD} | {time:hh:mm:ss} | "
    "{level:<7} | "
    "{function:>18}:{line:<4} | "
    "{message}"
)


def colorize():
    """
    check whether terminal/tty supports color
    """
    # check if we're in IPython/jupyter
    try:
        import IPython
        tty1 = bool(IPython.get_ipython())
    except ImportError:
        tty1 = False
    # check if we're in a terminal
    tty2 = sys.stderr.isatty()
    return tty1 or tty2


def set_loglevel(loglevel="DEBUG", stderr=True, logfile=None):
    """
    Config and start the logger
    """
    config = {'handlers': []}

    if stderr:
        stderr_logger = dict(
            sink=sys.stderr, 
            format=STDFORMAT, 
            level=loglevel,
            colorize=colorize(),
        )
        config["handlers"].append(stderr_logger)

    if logfile:
        file_logger = dict(
            sink=logfile,
            format=LOGFORMAT, 
            level=loglevel,
            enqueue=True,
            rotation="50 MB",
            backtrace=True, 
            diagnose=True,
        )
        config["handlers"].append(file_logger)

    logger.configure(**config)
    logger.enable("ipyrad")
    logger.debug("")
    logger.debug(f"ipyrad logging enabled at loglevel={loglevel}")
