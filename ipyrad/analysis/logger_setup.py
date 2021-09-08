#!/usr/bin/env python

"""
ipyrad-analysis logger module.
"""

import sys
from loguru import logger
from ipyrad.core.logger_setup import colorize

ANALYSIS_STDERR_LOGGER_FORMAT = (
    "<level>ipa: {file}</level> <white>|</white> "
    "<black>{message}</black>"
)

IPA_LOGGERS = []


def set_log_level(log_level="DEBUG"):
    """Add logger for ipyrad-analysis to stderr.

    The IDs of the loggers are stored in a global variable
    IP_LOGGERS so they can be removed and updated, not duplicated.
    """
    # remove any previously assigned loggers for ip and ipa. This uses
    # try/except in case run it multiple times in a row.
    for log_id in IPA_LOGGERS:
        try:
            logger.remove(log_id)
        except ValueError:
            pass

    log_id = logger.add(
        sink=sys.stderr,
        level=log_level,
        colorize=colorize(),
        format=ANALYSIS_STDERR_LOGGER_FORMAT,
        filter=lambda x: x['extra'].get("ipa"),
    )
    IPA_LOGGERS.append(log_id)
