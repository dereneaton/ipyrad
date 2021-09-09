#!/usr/bin/env python

"""
ipyrad-analysis logger module.
"""

import sys
from loguru import logger
from ipyrad.core.logger_setup import colorize

ANALYSIS_STDERR_LOGGER_FORMAT = (
    "<level>ipa</level> <white>|</white> "
    "<level>{file}</level> <white>|</white> "
    "<black>{message}</black>"
)


LOGGERS = [0]
def set_log_level(log_level="DEBUG"):
    """Add logger for ipa to stderr.

    These loggers are bound to the 'extra' keyword 'ipa'. Thus, any
    module in assembly that aims to use this formatted logger should
    put `logger = logger.bind(name="ipa")` at the top of the module.
    """
    for idx in LOGGERS:
        try:
            logger.remove(idx)
        except ValueError:
            pass
    idx = logger.add(
        sink=sys.stderr,
        level=log_level,
        colorize=colorize(),
        format=ANALYSIS_STDERR_LOGGER_FORMAT,
        filter=lambda x: x['extra'].get("name") == "ipa",
    )
    LOGGERS.append(idx)
    logger.enable("ipa")
