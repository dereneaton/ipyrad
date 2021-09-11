#!/usr/bin/env python

"""
Progress bars
"""

import os
import sys
import time
import traceback
import datetime
import ipyparallel
import IPython
from loguru import logger
from ipyrad.assemble.utils import IPyradError


logger = logger.bind(name="ipyrad")


class AssemblyProgressBar:
    """
    Print pretty progress bar with printings specific to Assembly object
    """       
    def __init__(self, jobs, message, step, quiet=False, start=None):
        self.jobs = jobs
        self.start = (start if start else time.time())
        self.message = message
        self.step = step
        self.quiet = quiet
        self.finished = 0

        # filled by check_success
        self.results = {}

    @property
    def progress(self):
        "returns the percent progress as a float"
        if not self.jobs:
            return 0
        return 100 * (self.finished / float(len(self.jobs)))

    @property
    def elapsed(self):
        "returns the elapsed time in nice format"
        return datetime.timedelta(seconds=int(time.time() - self.start))

    def update(self, final=False):
        "flushes progress bar at current state to STDOUT"
        if self.quiet:
            return
        # build the bar
        hashes = '#' * int(self.progress / 5.)
        nohash = ' ' * int(20 - len(hashes))

        # print to stderr
        message = (
            f"\r[{hashes + nohash}] "
            f"{int(self.progress):>3}% {self.elapsed} | "
            f"{self.message.ljust(20)} | {self.step} |"
        )
        print(message, end="\n" if final else "")
        sys.stdout.flush()

    def block(self):
        """
        Tracks completion of asynchronous result objects in a while 
        loop until they are finished, checking every 0.5 seconds.
        Prints progress either continuously or occasionally, depending
        on TTY output type.
        """
        # get TTY output type of STDOUT (or os.isatty(1))
        isatty = sys.stdout.isatty() or bool(IPython.get_ipython())

        # store the current value
        cur_finished = 0

        # loop until all jobs are finished
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
                self.update(final=True)
                break

            # -- decide verbosity based on tty or change occurred ---
            # if value changed then print
            if cur_finished != self.finished:
                self.update()

            # else, only print every 1 minutes if not in tty
            elif not isatty:
                if not int(self.elapsed.seconds) % 30:
                    self.update()
            
            # normal tty print every second
            else:
                self.update()

            # update cur_finished to match finished
            cur_finished = self.finished
            time.sleep(1)            

    def check(self):
        """
        Will log and raise an error with traceback. Stores results 
        into self.results
        """
        # check for failures:
        for job in self.jobs:
            try:
                self.results[job] = self.jobs[job].get()
            except ipyparallel.RemoteError as inst:
                raise inst
                # logger.error(
                #     "An error occurred, see logfile "
                #     f"and trace:\n{traceback.format_exc()}"
                # )
                # logger.error(
                    # "An error occurred: SEE TRACE BELOW\n" + 
                    # '\n'.join(inst.render_traceback())
                # )
                raise IPyradError("Exception on remote engine.") from inst


    def engine_log(self, key):
        """
        Logs the current stdout from engine and then clears it.
        """
        if self.jobs[key].stdout:
            logger.info(self.jobs[key].stdout.strip())
            self.jobs[key].stdout = ""



if __name__ == "__main__":

    pass
