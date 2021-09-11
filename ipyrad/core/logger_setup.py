#!/usr/bin/env python

"""
Logger for stderr and optionally to file.
"""

import os
import sys
from loguru import logger
import IPython


def formatter(record):
    """Custom formatter that allows for progress bar."""
    end = record["extra"].get("end", "\n")
    fmessage = (
        "<level>{level:}</level> <white>|</white> "
        "<magenta>{file:<18}</magenta> <white>|</white> "
        "{message}"
    ) + end
    return fmessage

ASSEMBLY_FILE_LOGGER_FORMAT = (
    "{time:YYYY/MM/DD} | {time:hh:mm:ss} | "
    "{level:<7} | "
    "{function:>18}:{line:<4} | "
    "{message}"
)


def color_support():
    """Check for color support in stderr as a notebook or terminal/tty."""
    # check if we're in IPython/jupyter
    tty1 = bool(IPython.get_ipython())
    # check if we're in a terminal
    tty2 = sys.stderr.isatty()
    return tty1 or tty2


LOGGERS = [0]
def set_log_level(log_level="DEBUG", log_file=None):
    """Add logger for ipyrad to stderr and optionally to file.

    These loggers are bound to the 'extra' keyword 'ipyrad'. Thus, any
    module in assembly that aims to use this formatted logger should
    put `logger = logger.bind(name="ipyrad")` at the top of the module.
    """
    for idx in LOGGERS:
        try:
            logger.remove(idx)
        except ValueError:
            pass
    idx = logger.add(
        sink=sys.stderr,
        level=log_level,
        colorize=color_support(),
        format=formatter,#ASSEMBLY_STDERR_LOGGER_FORMAT,
        filter=lambda x: x['extra'].get("name") == "ipyrad",
    )
    LOGGERS.append(idx)

    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        idx = logger.add(
            sink=log_file,
            level=log_level,
            colorize=False,
            format=ASSEMBLY_FILE_LOGGER_FORMAT,
            filter=lambda x: x['extra'].get('name') == "ipyrad",
            enqueue=True,
            rotation="50 MB",
            backtrace=True,
            diagnose=True,
        )
        LOGGERS.append(idx)

    logger.enable("ipyrad")
