#!/usr/bin/env python

"""Concatenate reads from technical replicates into a new single file.

"""

from typing import List
from pathlib import Path
from subprocess import Popen, STDOUT
from ipyrad.core.exceptions import IPyradError


def concatenate_technical_replicates(
    outdir: Path,
    name: str,
    files: List[Path],
    pair: int = 1,
) -> Path:
    """Return Path to new concatenated file.

    """
    # cat command works for gzip or not.
    cmd = ["cat"] + [str(i) for i in files]

    # add '_concat' to name in case new name is same as an old one
    outfile = outdir / f"{name}.trimmed_R{pair}_concat.fastq.gz"

    # call and write to open outfile
    with open(outfile, 'w', encoding="utf-8") as out:
        with Popen(cmd, stderr=STDOUT, stdout=out, close_fds=True) as proc:
            res = proc.communicate()[0]
            if proc.returncode:
                raise IPyradError(f"error in: {cmd}, {res}")

    # remove the files
    for file in files:
        Path(file).unlink()

    # return the new handle
    return outfile


if __name__ == "__main__":

    OUTDIR = Path("/home/deren/Documents/ipyrad/pedtest/half-demuxed_trimmed_fastqs/")
    NAME = "kansuensis-merge1"
    FILES = [
        Path("/home/deren/Documents/ipyrad/pedtest/half-demuxed_trimmed_fastqs/kansuensis-DE662.trimmed_R2.fastq.gz"),
        Path("/home/deren/Documents/ipyrad/pedtest/half-demuxed_trimmed_fastqs/kansuensis-DE739.trimmed_R2.fastq.gz"),
    ]
    OUTFILE = concatenate_technical_replicates(OUTDIR, NAME, FILES, 2)
    print(OUTFILE)
