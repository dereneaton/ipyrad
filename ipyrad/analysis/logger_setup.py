#!/usr/bin/env python

"""Starts separate logger for ipa.


"""

import sys
from loguru import logger
from ipyrad.core.logger_setup import color_support


# def formatter(record):
#     """Custom formatter that allows for progress bar."""
#     start = record["extra"].get("start", "")
#     end = record["extra"].get("end", "\n")
#     fmessage = start + (
#         "<level>ipa</level> <white>|</white> "
#         "<level>{file}</level> <white>|</white> "
#         "<black>{message}</black>"
#     ) + end
#     return fmessage


def formatter(record):
    """Custom formatter that allows for progress bar."""
    start = record["extra"].get("start", "")
    end = record["extra"].get("end", "\n")
    lev = record['level'].name[:4]
    fname = record['file'].name[:12]
    fmessage = start + (
        f"<level>{lev}</level> <white>|</white> "
        f"<fg #00005f>{fname: <12}</fg #00005f> <white>|</white> "
        "<level><normal>{message}</normal></level>"
    ) + end
    return fmessage

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
        colorize=color_support(),
        format=formatter,
        filter=lambda x: x['extra'].get("name") == "ipa",
    )
    LOGGERS.append(idx)
    logger.enable("ipa")


if __name__ == "__main__":
    pass