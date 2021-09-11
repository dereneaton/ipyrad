#!/usr/bin/env python

"""
Progress bars
"""

from typing import Dict
import sys
import time
from loguru import logger

logger = logger.bind(name="ipa")


class ProgressBar:
    """Print pretty progress bar to logger for ipa."""
    def __init__(self, jobs: Dict, message: str, delay=1):
        self.jobs = jobs
        self.message = message
        self.finished = 0
        self.results = {}
        self.delay = delay
        logger.bind(end="").info(f"[                    ] {self.message}")

    @property
    def progress(self):
        """returns the percent progress as a float"""
        if not self.jobs:
            return 0
        return 100 * (self.finished / float(len(self.jobs)))

    def update(self, final: bool=False):
        """Flushes stderr and logs progress bar at current state"""
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))
        sys.stderr.flush()
        logger.bind(start="\r", end="\n" if final else "").info(
            f"[{hashes + nohash}] {self.message}"
        )

    def block(self):
        """Block and show progress until all jobs in .jobs are complete."""
        while 1:
            # check for finished jobs and use stdout as a hack to log
            # and then clear messages from engines.
            self.finished = 0
            for job in self.jobs:
                if self.jobs[job].ready():
                    self.finished += 1
                # check each engine stdout for log messages
                self.engine_log(job)

            # flush progress and end
            if self.finished == len(self.jobs):               
                self.update(True)
                break
            self.update()
            time.sleep(self.delay)

    def check(self):
        """Gather results and raise exceptions if failed."""
        for job in self.jobs:
            self.results[job] = self.jobs[job].get()

    def engine_log(self, key):
        """Send engine stdout to logger if present and clear."""
        if self.jobs[key].stdout:
            logger.info(self.jobs[key].stdout.strip())
            self.jobs[key].stdout = ""


if __name__ == "__main__":


    import ipyrad.analysis as ipa

    prog = ProgressBar({3:4, 4:5, 5:6}, "hello world")
    logger.info(prog.update())
    prog.finished += 1
    logger.info(prog.update())    
    prog.finished += 1
    logger.info(prog.update())    
    prog.finished += 1
    logger.info(prog.update(final=True))
