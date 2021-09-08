#!/usr/bin/env python

"""
ipyrad-analysis logger module.
"""

import sys
from loguru import logger
from ipyrad.core.logger_setup import colorize
import ipyrad.analysis as ipa


STDFORMAT = (
    "<level>ipa: {file}</level> <white>|</white> "
    "<black>{message}</black>"
)

def set_loglevel(loglevel: str="INFO"):
    """Config and start the ipa logger."""
    config = {'handlers': []}

    # stderr is always set to INFO if logfile is in use else it 
    # uses the specified loglevel
    stderr_logger = dict(
        sink=sys.stderr, 
        format=STDFORMAT, 
        level=loglevel,
        colorize=colorize(),
    )
    config["handlers"].append(stderr_logger)

    logger.configure(**config)
    logger.enable("ipyrad-analysis")

    # log message if set by user to DEBUG (default from init is INFO)
    logger.debug(f"ipyrad-analysis (v.{ipa.__version__}) logging={loglevel}")
