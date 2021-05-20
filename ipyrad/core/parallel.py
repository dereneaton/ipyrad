#!/usr/bin/env python

"""
Starts an ipcluster (IPython cluster) and writes the pid/ipcluster.pid 
into a specified directory for killing cluster later.
"""

import os
import sys
import time
import platform
import signal
import random
import subprocess
from threading import Thread
import ipyparallel as ipp
from loguru import logger
from ipyrad.assemble.utils import IPyradError


class Cluster:
    """
    Stores cluster information.
    """
    def __init__(self, quiet=False):

        # get random name for this cluster to avoid conflicts
        # if one job stops and another starts right away the delay
        # in the kill thread can cause it to kill the new job if 
        # they don't differ by some random element in the name.
        self.name = f"{os.getpid()}-{random.randint(0, 9999999)}"
        self.quiet = quiet

        # to be filled: client, and pids for shutting down
        self.cores = 0
        self.ipyclient = None
        self.auto_started = False
        self.engine_pids = {}

        # futures should be stored in this dict
        self.rasyncs = {}


    def start(self, cores=0, ipyclient=None):
        """
        Create connection to an ipp.Client.
        """
        self.cores = (cores if cores else get_num_cpus())
        self.ipyclient = ipyclient
        self.auto_connect()
        self.store_pids()


    def store_pids(self):
        """
        Store process ids of ipcluster and engine processes
        """
        # store pids for interrupting engine jobs
        self.engine_pids = self.ipyclient[:].apply(os.getpid).get_dict()


    def auto_connect(self):
        """
        Start an ipcluster instance, or connect to a running instance.
        """
        # the user provided an already active ipclient through the API
        if hasattr(self.ipyclient, "ids"):
            return

        # start an ipcluster with a local PID stored in workdir/pids
        self.auto_started = True
        self.start_ipcluster()

        # It can sometimes take a long time for the ipcluster to startup
        # engines on some machines, especially HPC, so we check the 
        # number of engines it finds until all requested engines started.
        # Loop until all engines are found
        while 1:

            # if ipcluster file is not made yet this will raise an
            # OSError, in which case just want to wait a tiny bit.
            try: 
                # connect to client
                self.ipyclient = ipp.Client(cluster_id=self.name)

                # if all engines are found then break
                if len(self.ipyclient) == self.cores:
                    logger.info(
                        f"cluster established: {len(self.ipyclient)} engines")
                    break

                # close it again (prevents too many open files)
                if not self.quiet:
                    print("\rEstablishing parallel cluster ...", end="")                        
                self.ipyclient.close()
                sys.stdout.flush()
                time.sleep(0.2)

            except OSError:
                pass

        if not self.quiet:
            print(f"\rEstablishing parallel cluster ({len(self.ipyclient)} cores)")


    def start_ipcluster(self):
        """
        Start ipcluster by calling it from subprocess. We give it the
        argument profile_dir so that it will write the pid of the job
        to a local file so we can easily shut it down later. We also
        """
        # avoid starting a null cluster
        assert self.cores > 0, (
            'start ipcluster using the .run function, or set .cores > 0')

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


    def shutdown(self):
        """
        Calls stop from the ipcluster binary 
        This is really really really the only reliable way we could 
        find to get everything to stop without zombies.
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


    def cleanup(self, rasyncs=None):
        """
        If ipcluster was auto-started then it is shutdown, otherwise
        we simply cancell all jobs.
        """
        # protect from keyboard interrup while cleaning up
        try:
            # auto-started: stop ipcluster
            if self.auto_started:
                # stop future jobs
                self.ipyclient.abort()
                # stop current jobs
                if rasyncs:
                    for job in rasyncs:
                        if rasyncs[job].stdout:
                            pid = int(rasyncs[job].stdout.strip())
                            os.kill(pid, signal.SIGINT)
                # give it a second
                time.sleep(1)
                self.ipyclient.close()
                # send shutdown to cluster and controller
                self.shutdown()
                logger.info('ipcluster stopped')

            # user-entered ipclient: leave ipcluster and connected, just stop jobs
            else:
                self.ipyclient.abort()
                for pid in self.engine_pids.values():
                    os.kill(pid, signal.SIGINT)
                time.sleep(1)
                if not self.ipyclient.outstanding:
                    self.ipyclient.purge_everything()
                self.ipyclient.close()  # do this?

        except KeyboardInterrupt:
            self.ipyclient.close()
            logger.warning("cleaning up...")

            

    def cleanup_safely(self, rasyncs):
        """
        Start cleanup on a separate thread that cannot be interrupted
        by a keyboard interrupt.
        """
        try:
            waitjob = Thread(target=self.cleanup, args=(rasyncs,))
            waitjob.start()
            waitjob.join()
        except KeyboardInterrupt:
            logger.warning("too late sucker, already sent the job.")



def get_num_cpus():
    """
    Return the effective number of CPUs in the system as an integer for
    either Unix or MacOSX (Darwin). This code is mostly copied from a
    similar implementation in IPython.
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



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG")

    cluster = Cluster()
    try:
        cluster.start(cores=0, ipyclient=None)
        time.sleep(10)

    except KeyboardInterrupt:
        logger.warning("keyboard interrupt by user, cleaning up.")

    except IPyradError:
        logger.error("Error occurred, cleaning up.")
        raise

    except Exception:
        logger.error("An unexpected error occurred, cleaning up.")
        raise

    finally:
        cluster.cleanup_safely(None)
