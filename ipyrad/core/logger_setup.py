#!/usr/bin/env python

"""Logger for ipyrad to STDERR and optionally also to a LOGFILE.

logging to STDERR
-----------------
DEBUG: used by developers to examine extra details.
INFO: info reported to users, including progress bars. (DEFAULT)
WARNING: warnings to users, if set to default then progress bars are not shown.
ERROR: sometimes printed along with raised errors.

logging to LOGFILE
------------------
DEBUG: developer stuff
INFO: same as above, w/ some extra info, but not progress bars. (DEFAULT)
same
same

Examples
--------
>>> import ipyrad as ip
>>> ip.set_log_level("DEBUG")
>>> ip.set_log_level("DEBUG", log_file="/tmp/ip-log.txt")
"""

import os
import sys
from loguru import logger
import IPython

def formatter(record):
    """Custom formatter that allows for progress bar."""
    end = record["extra"].get("end", "\n")
    fmessage = (
        "<level>{level:<8}</level> <white>|</white> "
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
    logger.bind(name='ipyrad').debug(f"ipyrad logging enabled: {log_level}")

if __name__ == "__main__":
    set_log_level("DEBUG")
