#!/usr/bin/env python

"""Globals used commonly.

"""

import string


# used in demux_raw.py and demux_sorted.py to fix sample names.
BADCHARS = (
    string.punctuation
    .replace("_", "")
    .replace("-", "")
    .replace(".", "") + " "
)

# used in demux_raw.py to resolve ambiguous cutters
AMBIGS = {
    "R": ("G", "A"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("T", "C"),
    "W": ("T", "A"),
    "M": ("C", "A"),
}


def comp(seq: str) -> str:
    """ returns a seq with complement. Preserves little n's for splitters."""
    return seq.replace("A", 't')\
              .replace('T', 'a')\
              .replace('C', 'g')\
              .replace('G', 'c')\
              .replace('n', 'Z')\
              .upper()\
              .replace("Z", "n")

