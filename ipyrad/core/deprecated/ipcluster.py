#!/usr/bin/env python

"""
DEPRECATED BY RELEASE OF NEW IPYPARALLEL

A context wrapper around the ipcluster auto-start tool to ensure
start and cleanup.
"""

from typing import Dict
import re
import os
import sys
import time
import random
import signal
import traceback
import platform
import subprocess
from threading import Thread
from contextlib import AbstractContextManager
from loguru import logger
import IPython
import ipyparallel

#
# My old attempt before ipyparallel update.
#

class Ipcluster(AbstractContextManager):
    """Context Wrapper class for ipcluster management.

    Examples
    --------
    >>> with Ipcluster(cores=20, ipyclient=None) as ipyclient:
    >>>     lbview = ipyclient.load_balanced_view()
    >>>     for job in jobs:
    >>>         lbview.apply(func, args)
    """
    def __init__(self, cores: int=4, ipyclient: 'ipyparallel.client.Client'=None):
        self.cores=cores
        self.ipyclient: 'ipyparallel.client.Client'=ipyclient
        self.name: str=f"{os.getpid()}-{random.randint(0, 2**31)}"
        self.auto_started: bool=False
        self.engine_pids: dict={}

    def __repr__(self):
        return f"<ipcluster ncores={self.cores} auto-started={self.auto_started}>"

    def __enter__(self):
        """On entry ipyclient is stored and if empty one is started."""       
        # if ipyclient is already connected use it.
        logger.info(f"ipyclient={self.ipyclient}")
        if hasattr(self.ipyclient, "ids"):
            self.auto_started = False            
            self.cores = len(self.ipyclient)
            logger.info(f"Connected to cluster {self.ipyclient.profile}: {self.cores} engines.")

        # if no ipyclient, start a new one with number of requested cores.
        else:
            self.auto_started = True
            self.cores = self.cores if self.cores else get_num_cpus()
            self._start_ipcluster()
            self.ipyclient = self._auto_connect()

        # store engine pids for stopping later.
        self.engine_pids = self.ipyclient[:].apply(os.getpid).get_dict()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Ensures shutdown of ipcluster and handling of exceptions."""
        # send cleanup job
        self._cleanup_safely()

        # check for and raise any exceptions that occurred. Because 
        # ipyparallel returns 
        if exc_value:
            if exc_type == KeyboardInterrupt:
                logger.error("keyboard interrupt by user, cleaning up.")
            elif exc_type == ipyparallel.error.RemoteError:
                trace = "\n".join(exc_value.render_traceback())
                if not color_support():              
                    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
                    trace = ansi_escape.sub('', trace)
                logger.error(f"An error occurred on an ipengine:\n{trace}")               
            else:
                trace = traceback.format_exc()
                logger.error(f"An error occurred:\n{trace}")

    def _start_ipcluster(self):
        """Start ipcluster by calling CLI from subprocess. 

        We give it the argument profile_dir so that it will write 
        the pid of the job to a local file so we can easily shut it 
        down later.
        """
        # path to the ipcluster binary
        ipcluster_bin = os.path.join(sys.prefix, "bin", "ipcluster")

        # build the command
        cmd = [
            ipcluster_bin,
            "start",
            "--cluster-id={}".format(self.name),
            "--n={}".format(self.cores),
            "--quiet",
            "--daemonize",            
        ]

        # start binary running on fork
        subprocess.run(cmd, check=True)
        logger.debug(" ".join(cmd))
        time.sleep(1)

    def _auto_connect(self):
        """Return an ipyclient connected to the new auto-started client."""
        logger.bind(end="").info("Establishing parallel cluster: ")

        # It can sometimes take a long time for the ipcluster to startup
        # engines on some machines, especially HPC, so we check the 
        # number of engines it finds until all requested engines started.
        # Loop until all engines are found
        while 1:

            # if ipcluster file is not made *yet* this will raise an
            # OSError, in which case just want to wait a tiny bit for
            # the ipcluster CLI that has been started to create it.
            try:
                ipyclient = ipyparallel.Client(cluster_id=self.name)

                # if all engines are found then break
                if len(ipyclient) == self.cores:
                    logger.opt(raw=True).info(f"{len(ipyclient)} engines.\n")
                    break

                # close it again (prevents too many open files!)
                ipyclient.close()
                time.sleep(0.1)
            except OSError:
                pass
        return ipyclient

    def _shutdown(self):
        """Calls `ipcluster stop` to shutdown hub and engines.

        This seems to really really really be the only reliable way
        we could find to get everything to stop without zombies. 

        Note
        ----
        This function itself is called after engines are killed by
        SIGINT in `cleanup` and that function is run on a separate 
        Thread so that it cannot be interrupted.
        """
        # path to the ipcluster binary
        ipcluster_bin = os.path.join(sys.prefix, "bin", "ipcluster")

        # build the command
        cmd = [
            ipcluster_bin,
            "stop",
            "--cluster-id={}".format(self.name),
            "--quiet",
        ]

        # run command and block until finished
        logger.debug(" ".join(cmd))
        subprocess.run(cmd, check=True)
        time.sleep(1)

    def _cleanup(self, rasyncs: Dict=None):
        """Cancels any running jobs, kills engines and hub if auto-started."""
        if self.auto_started:
            # cancel future jobs
            self.ipyclient.abort()
            # interrupt current jobs
            if rasyncs:
                for job in rasyncs:
                    # remote jobs should print their pid to stdout
                    if rasyncs[job].stdout:
                        pid = int(rasyncs[job].stdout.strip())
                        os.kill(pid, signal.SIGINT)
            time.sleep(1)
            self.ipyclient.close()
            self._shutdown()
            logger.info("ipcluster stopped.")

        else:
            pass
            # self.ipyclient.abort()
            # for pid in self.engine_pids.values():
            #     os.kill(pid, signal.SIGINT)
            # time.sleep(1)
            # logger.info(self.ipyclient.outstanding)
            # if not self.ipyclient.outstanding:
            #     self.ipyclient.purge_everything()
            # logger.info("ipcluster reset.")

    def _cleanup_safely(self, rasyncs=None):
        """Starts cleanup on a separate thread so it cannot be interrupted."""
        try:
            waitjob = Thread(target=self._cleanup, args=(rasyncs,))
            waitjob.start()
            waitjob.join()
        except KeyboardInterrupt:
            logger.warning("too late sucker, already sent the job.")


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

    import ipyrad as ip
    import numpy as np
    ip.set_log_level("DEBUG", "/tmp/ipyrad.log")

    def test_func(x, arr):
        np.concatenate(arr)
        if not isinstance(x, int):
            raise ValueError("x must be an int")
        return x

    import time
    start = time.time()
    with Cluster2(cores=2) as ipyclient:
        print(ipyclient)
        # raise KeyboardInterrupt()
        ipyclient[0].apply(test_func, *(3, np.zeros(shape=(3,3,3)))).get()

    print(time.time() - start)

    # cluster = Cluster()
    # cluster.start()   
    # start = time.time()
    # with Ipcluster(cores=2) as cluster:
        # print(cluster)
    # print(time.time() - start)        
        # rasync = cluster.ipyclient[0].apply(test_func, '3')
    #     print(rasync.get())
    #     raise AttributeError("oy mate")
