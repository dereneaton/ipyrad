#!/usr/bin/env ipython2

## imports for running ipcluster
from __future__ import print_function


import ipyparallel as ipp
import ipyrad as ip
import subprocess
import cStringIO
import atexit
import shlex
import time
import sys
import os

# pylint: disable=W0212
# pylint: disable=W0142

import logging
LOGGER = logging.getLogger(__name__)


## start ipcluster
def start_ipcluster(data):
    """ Start ipcluster """

    ## if MPI argument then use --ip arg to view all sockets
    iparg = ""
    if "MPI" in data._ipcluster["engines"]:
        iparg = "--ip=*"

    ## make ipcluster arg call
    standard = """
        ipcluster start 
                  --daemonize 
                  --cluster-id={}
                  --engines={} 
                  --profile={}
                  --n={}
                  {}"""\
        .format(data._ipcluster["cluster_id"], 
                data._ipcluster["engines"], 
                data._ipcluster["profile"],
                data._ipcluster["cores"],
                iparg)
                   
    ## wrap ipcluster start
    try: 
        LOGGER.info(shlex.split(standard))
        subprocess.check_call(shlex.split(standard), 
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE)

    except subprocess.CalledProcessError as inst:
        LOGGER.debug("  ipcontroller already running.")
        raise

    except Exception as inst:
        sys.exit("  Error launching ipcluster for parallelization:\n({})\n".\
                 format(inst))




## decorated func for stopping. Does not need to be called?
def stop_ipcluster(profile, cluster_id):
    """ stop ipcluster at sys.exit """

    LOGGER.info("Shutting down [%s] remote Engines", cluster_id)
    stopcall = ["ipcluster", 
                "stop", 
                "--profile="+profile,
                "--cluster-id="+cluster_id]
    try:
        subprocess.check_call(stopcall, 
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE)
    except subprocess.CalledProcessError:
        pass



def register_ipcluster(data):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. 
    """
    ## check if this pid already has a running cluster
    data._ipcluster["cluster_id"] = "ipyrad-"+str(os.getpid())
    start_ipcluster(data)
    atexit.register(stop_ipcluster, 
                    data._ipcluster["profile"],
                    data._ipcluster["cluster_id"])
    return data



def get_client(cluster_id, profile, engines, timeout, quiet, **kwargs):
    """ 
    Creates a client to view ipcluster engines for a given profile and 
    returns it with at least one engine spun up and ready to go. If no 
    engines are found after nwait amount of time then an error is raised.
    If engines==MPI it waits a bit longer to find engines.
    """

    ## save stds for later, we're gonna hide them to prevent external printing 
    save_stdout = sys.stdout 
    save_stderr = sys.stderr
    sys.stdout = cStringIO.StringIO()
    sys.stderr = cStringIO.StringIO()

    try: 
        ## are we looking for a running ipcluster instance?
        if profile not in [None, "default"]:
            args = {'profile': profile, "timeout": timeout}
        else:
            clusterargs = [cluster_id, profile, timeout]
            argnames = ["cluster_id", "profile", "timeout"]
            args = {key:value for key, value in zip(argnames, clusterargs)}

        ## get connection within timeout window of wait time and hide messages
        ipyclient = ipp.Client(**args)
        sys.stdout = save_stdout
        sys.stderr = save_stderr

        ## check that all engines have connected            
        for _ in range(6000):
            initid = len(ipyclient)
            time.sleep(0.01)

            ## If MPI then wait for all engines to start so we can report
            ## how many cores are on each host. If Local then only wait for
            ## one engine to be ready and then just go.
            if not engines == "MPI":
                if initid:
                    break
            else:
                if not quiet:
                    print(\
            "\r  establishing MPI connection to remote hosts...", end="")
                    sys.stdout.flush()
                time.sleep(3)
                if initid:
                    if len(ipyclient) == initid:
                        print("")
                        break


    except KeyboardInterrupt as inst:
        ## ensure stdout is reset even if Exception was raised            
        sys.stdout = save_stdout
        sys.stderr = save_stderr
        raise inst

    ## This is raised if ipcluster is not running ------------
    except IOError as inst:
        ## ensure stdout is reset even if Exception was raised
        sys.stdout = save_stdout
        sys.stderr = save_stderr
        if ip.__interactive__:
            raise ip.assemble.util.IPyradWarningExit(NO_IPCLUSTER_API)
        else:
            raise ip.assemble.util.IPyradWarningExit(NO_IPCLUSTER_CLI)

    except (ipp.TimeoutError, ipp.NoEnginesRegistered) as inst:
        ## raised by ipp if no connection file is found for 'nwait' seconds
        sys.stdout = save_stdout
        sys.stderr = save_stderr
        raise inst

    except Exception as inst:
        ## if any other exceptions were missed...
        sys.stdout = save_stdout
        sys.stderr = save_stderr
        raise inst

    finally:
        ## ensure that no matter what we reset the stds
        sys.stdout = save_stdout
        sys.stderr = save_stderr

    return ipyclient


## GLOBALS AND EXCEPTIONS

NO_IPCLUSTER_CLI = """\
    No ipcluster instance found. This may be a problem with your installation
    setup. I would recommend that you contact the ipyrad developers. 
    """
NO_IPCLUSTER_API = """
    No ipcluster instance found. See documentation for the proper way to set 
    up an ipcluster instance when running the ipyrad Python API. In short, 
    you must run 'ipcluster start' to initiate a local or remote cluster. 
    """

