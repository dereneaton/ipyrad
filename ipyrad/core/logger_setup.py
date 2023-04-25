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

Note
----
Exceptions written to the logfile have color support, which 
can be viewed using `less -R logfile.txt`
"""

from typing import Optional
import sys
from pathlib import Path
from loguru import logger
import IPython

LOGGERS = [0]


def formatter(record):
    """Custom formatter that allows for progress bar."""
    end = record["extra"].get("end", "\n")
    fmessage = (
        "{time:hh:mm:ss} | "
        "<level>{level:<8}</level> <white>|</white> "
        "<magenta>{file:<18}</magenta> <white>|</white> "
        "{message}"
    ) + end
    return fmessage


def color_support():
    """Check for color support in stderr as a notebook or terminal/tty."""
    # check if we're in IPython/jupyter
    tty1 = bool(IPython.get_ipython())
    # check if we're in a terminal
    tty2 = sys.stderr.isatty()
    return tty1 or tty2


def set_log_level(log_level: str = "DEBUG", log_file: Optional[Path] = None):
    """Add logger for ipyrad to stderr and optionally to file.

    These loggers are bound to the 'extra' keyword 'ipyrad'. Thus, any
    module in assembly that aims to use this formatted logger should
    put `logger = logger.bind(name="ipyrad")` at the top of the module.

    The logger will use EITHER a STDERR or a LOGFILE, but not both.
    """
    # remove any previous loggers created by ipyrad
    for idx in LOGGERS:
        try:
            logger.remove(idx)
        except ValueError:
            pass

    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(exist_ok=True)
        log_file.touch(exist_ok=True)
        idx = logger.add(
            sink=log_file,
            level=log_level,
            colorize=False,
            format=formatter,
            filter=lambda x: x['extra'].get('name') == "ipyrad",
            enqueue=True,
            rotation="50 MB",
            # backtrace=True,
            # diagnose=True,
        )
    else:
        idx = logger.add(
            sink=sys.stderr,
            level=log_level,
            colorize=color_support(),
            format=formatter,  # ASSEMBLY_STDERR_LOGGER_FORMAT,
            filter=lambda x: x['extra'].get("name") == "ipyrad",
            enqueue=True,
        )
    LOGGERS.append(idx)

    # activate
    logger.enable("ipyrad")
    logger.bind(name='ipyrad').debug(f"ipyrad logging enabled: {log_level}")


def get_logger():
    return logger.bind(name="ipyrad")


if __name__ == "__main__":
    set_log_level("DEBUG")
    log = get_logger()
    log.info("HI")
