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

# pylint: disable=invalid-name, abstract-method

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

    def __init__(self, cores: int, name: str = "ipyrad", **kwargs):
        # cores is an alias for .n, which is also stored for ipp parent.
        self.cores = self.n = cores if cores else get_num_cpus()

        # can limit logging to a name bound logger (see top of this module)
        self.logger_name = name

        # init parent class with kwargs for ipcluster start (e.g., MPI)
        super().__init__(**kwargs)

        # hidden attributes
        self._client_start_time = None
        # self._context_client = None

    def __enter__(self):
        """Starts a new cluster and connects a client."""
        # log engine starter to info
        logger.bind(name=self.logger_name, end="").info("Establishing parallel ipcluster: ")

        # start cluster and wait for n engines to start
        self.start_cluster_sync(n=self.n)
        client = self.connect_client_sync()
        client.wait_for_engines(n=self.n, block=True, interactive=False)

        # log engines started to same line of info
        logger.bind(name=self.logger_name).opt(raw=True).info(f"{len(client)} engines\n")

        # store client and start to log runtime at stoppage and close connection.
        self._context_client = client
        self._client_start_time = time.time()
        return client

    def __exit__(self, *args):
        """Ensures shutdown of ipcluster and handling of exceptions."""
        # abort outstanding (unstarted) jobs
        self._context_client.abort()

        # send SIGINT to any still running jobs
        self.signal_engines_sync(signum=2)

        # close client connection
        if self._context_client:
            self._context_client.close()
            self._context_client = None

        # stop engines and controller
        self.stop_cluster_sync()

        # log runtime to info
        elapsed = int(time.time() - self._client_start_time)
        elapsed = str(timedelta(seconds=elapsed))
        logger.bind(name=self.logger_name).info(
            f"ipcluster stopped. Elapsed time: {elapsed}")

        # raise traceback for any exceptions that occurred
        log_traceback(*args)


def log_traceback(exc_type, exc_value, exc_traceback):
    """Check for and raise any exceptions that occurred.

    """
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
            logger.error(f"An error occurred on an ipengine, see below:\n{trace}")

        # error not on remote engine, occurred anywhere else.
        else:
            trace = traceback.format_exc()
            logger.error(f"An error occurred, see below:\n{trace}")
            #logger.warning(exc_traceback)


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
    except (TypeError, ValueError):
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

    import ipyrad as ip
    ip.set_log_level("DEBUG")
    with Cluster(cores=0) as c:
        print('ipyclient', c)
