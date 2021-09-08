#!/usr/bin/env python

"""
Logger for stderr and optionally to file.
"""

import os
import sys
from loguru import logger
import IPython


ASSEMBLY_STDERR_LOGGER_FORMAT = (
    "<cyan>{time:hh:mm:ss}</cyan> <white>|</white> "
    "<level>{level: <7}</level> <white>|</white> "
    "<magenta>{file:<18}</magenta> <white>|</white> "
    "<level>{message}</level>"
)

ASSEMBLY_FILE_LOGGER_FORMAT = (
    "{time:YYYY/MM/DD} | {time:hh:mm:ss} | "
    "{level:<7} | "
    "{function:>18}:{line:<4} | "
    "{message}"
)

IP_LOGGERS = []


def colorize():
    """
    check whether terminal/tty supports color
    """
    # check if we're in IPython/jupyter
    tty1 = bool(IPython.get_ipython())
    # check if we're in a terminal
    tty2 = sys.stderr.isatty()
    return tty1 or tty2


def set_log_level(log_level="DEBUG", log_file=None):
    """Add logger for ipyrad to stderr and optionally to file.

    The IDs of the loggers are stored in a global variable
    IP_LOGGERS so they can be removed and updated, not duplicated.

    These loggers are bound to the 'extra' keyword 'ip'. Thus, any
    module in assembly that aims to use this formatted logger should
    put `logger = logger.bind(ip=True)` at the top of the module.
    """
    # remove any previously assigned loggers for ip and ipa. This uses
    # try/except in case run it multiple times in a row.
    for log_id in IP_LOGGERS:
        try:
            logger.remove(log_id)
        except ValueError:
            pass

    # add a logger for assembly
    log_id = logger.add(
        sink=sys.stderr,
        level=log_level,
        colorize=colorize(),
        format=ASSEMBLY_STDERR_LOGGER_FORMAT,
        filter=lambda x: x['extra'].get("ip"),
    )
    IP_LOGGERS.append(log_id)

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        log_id = logger.add(
            sink=log_file,
            level=log_level,
            colorize=False,
            format=ASSEMBLY_FILE_LOGGER_FORMAT,
            filter=lambda x: x['extra'].get("ipa"),
            enqueue=True,
            rotation="50 MB",
            backtrace=True,
            diagnose=True,
        )
        IP_LOGGERS.append(log_id)
