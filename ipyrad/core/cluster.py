#!/usr/bin/env python

"""ipyrad Cluster class for starting/stopping w/ ipyparallel.

This uses the new ipyparallel cluster API (v.>7.0) for starting and
stopping clusters using a context wrapper. Our custom superclass
of the ipyparallel.Cluster class suppresses some print statements
of that tool, and adds our own custom logging and exception handling.
"""

import re
import os
import sys
import time
import traceback
import platform
import subprocess
from datetime import timedelta
from loguru import logger
import IPython
import ipyparallel


logger = logger.bind(name="ipyrad")


class Cluster(ipyparallel.cluster.cluster.Cluster):
    """Custom superclass of ipyparallel cluster.

    This class is used to start an ipcluster with an optional set of 
    additional kwargs, return a connected Client instance, and 
    shutdown the ipcluster when the context manager closes.

    Compared to the ipyparallel parent class, this one suppresses 
    print statements and instead uses a logger, and the context manager
    exit function has custom exception handling and formatting for
    ipyrad specifically.
    """
    # suppress INFO calls from ipyparallel built-in logging.
    log_level = 30

    def __init__(
        self, 
        cores: int, 
        logger_name: str="ipyrad",
        **kwargs,
        ):

        # cores is an alias for .n, which is also stored for ipp parent.
        self.cores = self.n = cores if cores else get_num_cpus()

        # can limit logging to a name bound logger (see top of this module)
        self.logger_name = logger_name

        # init parent class with kwargs for ipcluster start (e.g., MPI)
        super().__init__(**kwargs)

        # hidden attributes
        self._client_start_time = None
        self._context_client = None

    def __enter__(self):
        """Starts a new cluster and connects a client."""
        logger.bind(name=self.logger_name, end="").info(
            "Establishing parallel ipcluster: ")
        self.start_controller_sync()
        self.start_engines_sync()
        client = self._context_client = self.connect_client_sync()
        if self.n:
            # wait for engine registration
            client.wait_for_engines(
                self.n,
                interactive=False,
                block=True,
                timeout=self.engine_timeout,
            )
            logger.bind(name=self.logger_name).opt(raw=True).info(
                f"{len(client)} engines.\n")
        self._client_start_time = time.time()
        return client

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Ensures shutdown of ipcluster and handling of exceptions."""
        if self._context_client:
            self._context_client.close()
            self._context_client = None
        self.stop_engines_sync()
        self.stop_controller_sync()
        elapsed = int(time.time() - self._client_start_time)
        elapsed = str(timedelta(seconds=elapsed))        
        logger.bind(name=self.logger_name).info(
            f"ipcluster stopped. Elapsed time: {elapsed}")

        # check for and raise any exceptions that occurred.
        if exc_value:

            # nice message for interrupt.
            if exc_type == KeyboardInterrupt:
                logger.error("keyboard interrupt by user, cleaning up.")

            # remote errors raise with IPython ansi colored traceback
            # which looks like garbage if color is not supported so we
            # strip this if needed.
            elif exc_type == ipyparallel.error.RemoteError:
                trace = "\n".join(exc_value.render_traceback())
                if not color_support():
                    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
                    trace = ansi_escape.sub('', trace)
                logger.error(f"An error occurred on an ipengine:\n{trace}")

            # error not on remote engine, occurred anywhere else.
            else:
                trace = traceback.format_exc()
                logger.error(f"An error occurred:\n{trace}")


def get_num_cpus():
    """Return the effective number of CPUs in the system.

    Returns an integer for either Unix or MacOSX (Darwin). This code
    is mostly copied from a similar implementation in IPython.
    If it can't find a sensible answer, it returns 1.
    """
    if platform.system() == "Linux":
        ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
    else:
        proc = subprocess.run(
            ['sysctl', '-n', 'hw.ncpu'], check=True, capture_output=True)
        ncpus = proc.stdout.decode().strip()
    try:
        ncpus = max(1, int(ncpus))
    except:
        ncpus = 1
    return ncpus


def color_support():
    """Check for color support in stderr as a notebook or terminal/tty."""
    # check if we're in IPython/jupyter
    tty1 = bool(IPython.get_ipython())
    # check if we're in a terminal
    tty2 = sys.stderr.isatty()
    return tty1 or tty2


if __name__ == "__main__":

    # import ipyrad as ip
    # ip.set_log_level("DEBUG")
    with Cluster(cores=0) as c:
        c
