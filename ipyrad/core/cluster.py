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
import asyncio

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

    def __init__(self, cores: int, name: str="ipyrad", **kwargs):
        # cores is an alias for .n, which is also stored for ipp parent.
        self.cores = self.n = cores if cores else get_num_cpus()

        # can limit logging to a name bound logger (see top of this module)
        self.logger_name = name

        # init parent class with kwargs for ipcluster start (e.g., MPI)
        super().__init__(**kwargs)

        # hidden attributes
        self._client_start_time = None
        self._context_client = None

    def __enter__(self):
        """Starts a new cluster and connects a client."""
        logger.bind(name=self.logger_name, end="").info(
            "Establishing parallel ipcluster: ")
        self.start_cluster_sync(n=self.n)
        client = self.connect_client_sync()
        client.wait_for_engines(n=self.n, block=True, interactive=False)
        # self._context_client = self.start_and_connect_sync(n=self.n)
        logger.bind(name=self.logger_name).opt(raw=True).info(f"{len(client)} engines.\n")
        self._client_start_time = time.time()
        return client

    def __exit__(self, *args):
        """Ensures shutdown of ipcluster and handling of exceptions."""
        if self._context_client:
            self._context_client.close()
            self._context_client = None
        self.stop_cluster_sync()
        # self.stop_engines_sync()
        # self.stop_controller_sync()
        elapsed = int(time.time() - self._client_start_time)
        elapsed = str(timedelta(seconds=elapsed))
        logger.bind(name=self.logger_name).info(
            f"ipcluster stopped. Elapsed time: {elapsed}")

        # raise exception traceback
        log_traceback(*args)

    async def __aenter__(self):
        """Asynchronous start cluster.
        
        This is same as `client = self._context_client = await self.start_and_connect()`
        however, we need to call the more verbose version in order to 
        allow suppressing the default tqdm progress bar.
        """
        logger.bind(name=self.logger_name, end="").info("Establishing parallel ipcluster: ")
        await self.start_cluster(n=self.n)
        client = self._context_client = await self.connect_client()
        await asyncio.wrap_future(
            client.wait_for_engines(n=self.n, block=False, interactive=False)
        )
        logger.opt(raw=True).info(f"{len(client)} engines.\n")
        self._client_start_time = time.time()
        return client

    async def __aexit__(self, *args):
        """Asynchronous stop cluster."""
        if self._context_client is not None:
            self._context_client.close()
            self._context_client = None
        future = self.stop_cluster()

        elapsed = int(time.time() - self._client_start_time)
        elapsed = str(timedelta(seconds=elapsed))
        logger.info(f"ipcluster stopped. Elapsed time: {elapsed}")

        # raise exception traceback
        log_traceback(*args)
        await future

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
            logger.error(f"An error occurred on an ipengine:\n{trace}")

        # error not on remote engine, occurred anywhere else.
        else:
            trace = traceback.format_exc()
            logger.error(f"An error occurred:\n{trace}")
            #logger.warning(exc_traceback)


    # def __exit__(self, exc_type, exc_value, exc_traceback):
    #     """Ensures shutdown of ipcluster and handling of exceptions."""
    #     if self._context_client:
    #         self._context_client.close()
    #         self._context_client = None
    #     self.stop_engines_sync()
    #     self.stop_controller_sync()
    #     elapsed = int(time.time() - self._client_start_time)
    #     elapsed = str(timedelta(seconds=elapsed))        
    #     logger.bind(name=self.logger_name).info(
    #         f"ipcluster stopped. Elapsed time: {elapsed}")

    #     # check for and raise any exceptions that occurred.
    #     if exc_value:

    #         # nice message for interrupt.
    #         if exc_type == KeyboardInterrupt:
    #             logger.error("keyboard interrupt by user, cleaning up.")

    #         # remote errors raise with IPython ansi colored traceback
    #         # which looks like garbage if color is not supported so we
    #         # strip this if needed.
    #         elif exc_type == ipyparallel.error.RemoteError:
    #             trace = "\n".join(exc_value.render_traceback())
    #             if not color_support():
    #                 ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    #                 trace = ansi_escape.sub('', trace)
    #             logger.error(f"An error occurred on an ipengine:\n{trace}")

    #         # error not on remote engine, occurred anywhere else.
    #         else:
    #             trace = traceback.format_exc()
    #             logger.error(f"An error occurred:\n{trace}")


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
        print(c)
