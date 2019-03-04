#!/usr/bin/env python

""" functions to auto-launch an ipcluster instance """

## imports for running ipcluster
from __future__ import print_function

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import ipyparallel as ipp
import subprocess
import random
import socket
import time
import sys
# import os

from ipyrad.assemble.utils import IPyradWarningExit


def cluster_info(ipyclient, spacer=""):
    """ reports host and engine info for an ipyclient """    
    ## get engine data, skips busy engines.    
    hosts = []
    for eid in ipyclient.ids:
        engine = ipyclient[eid]
        if not engine.outstanding:
            hosts.append(engine.apply(socket.gethostname))

    ## report it
    hosts = [i.get() for i in hosts]
    result = []
    for hostname in set(hosts):
        result.append(
            "{}host compute node: [{} cores] on {}"
            .format(spacer, hosts.count(hostname), hostname))
    print("\n".join(result))



## start ipcluster
def start_ipcluster(data):
    """ Start ipcluster """

    ## if MPI argument then use --ip arg to view all sockets
    iparg = ""
    if "MPI" in data._ipcluster["engines"]:
        iparg = "--ip=*"

    # make ipcluster arg call
    standard = [
        "ipcluster", "start",
        "--daemonize", 
        "--cluster-id={}".format(data._ipcluster["cluster_id"]),
        "--engines={}".format(data._ipcluster["engines"]),
        "--profile={}".format(data._ipcluster["profile"]),
        "--n={}".format(data._ipcluster["cores"]),
        "{}".format(iparg),
    ]
                   
    # wrap ipcluster start
    try:
        subprocess.check_call(
            standard,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE)

    # cluster with THIS ID is running then kill it and try again
    except subprocess.CalledProcessError as inst:
        subprocess.check_call([
            "ipcluster", "stop", 
            "--cluster-id", data._ipcluster["cluster_id"],
        ], 
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
        )
        time.sleep(3)

        try:
            # try again to start it
            subprocess.check_call(
                standard, 
                stderr=subprocess.STDOUT,
                stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as inst:
            print(inst)
            raise

    except Exception as inst:
        sys.exit("  Error launching ipcluster for parallelization:\n({})\n".\
                 format(inst))



def register_ipcluster(data):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. Run statements are wrapped
    so that ipcluster SHOULD be killed on exit.
    """
    ## check if this random/pid already has a running cluster
    rand = random.randint(0, 1000)
    data._ipcluster["cluster_id"] = "ipyrad-cli-{}".format(rand)
    start_ipcluster(data)
    return data



def get_client(data):
    """ 
    Creates a client to view ipcluster engines for a given profile and 
    returns it with at least one engine spun up and ready to go. If no 
    engines are found after nwait amount of time then an error is raised.
    If engines==MPI it waits a bit longer to find engines. If the number
    of engines is set then it waits even longer to try to find that number
    of engines.
    """
    # shorthand for ipcluster dict
    class ipclusterdict:
        def __init__(self, data):
            for key, val in data._ipcluster.items():
                self.__setattr__(key, val)
            self.spacer = data._spacer
    ipd = ipclusterdict(data)

    ## save stds for later, we're gonna hide them to prevent external printing 
    save_stdout = sys.stdout 
    save_stderr = sys.stderr
    sys.stdout = StringIO()
    sys.stderr = StringIO()

    ## get cluster_info print string
    connection_string = (
        "{}establishing parallel connection:".format(ipd.spacer))

    ## wrapped search for ipcluster
    try: 
        ## are we looking for a running ipcluster instance?
        if ipd.profile not in [None, "default"]:
            args = {'profile': ipd.profile, "timeout": ipd.timeout}
        else:
            clusterargs = [ipd.cluster_id, ipd.profile, ipd.timeout]
            argnames = ["cluster_id", "profile", "timeout"]
            args = {key: value for key, value in zip(argnames, clusterargs)}

        ## get connection within timeout window of wait time and hide messages
        ipyclient = ipp.Client(**args)
        sys.stdout = save_stdout
        sys.stderr = save_stderr

        ## check that all engines have connected            
        if (ipd.engines == "MPI") or ("ipyrad-cli-" in ipd.cluster_id):
            if not ipd.quiet:
                print(connection_string)

        for _ in range(6000):
            initid = len(ipyclient)
            time.sleep(0.01)
            ## If MPI then wait for all engines to start so we can report
            ## how many cores are on each host. If Local then only wait for
            ## one engine to be ready and then just go.
            if (ipd.engines == "MPI") or ("ipyrad-cli-" in ipd.cluster_id):
                ## wait for cores to be connected
                if ipd.cores:
                    time.sleep(0.1)
                    if initid == ipd.cores:
                        break
                if initid:
                    time.sleep(3)
                    if len(ipyclient) == initid:
                        break
            else:
                if ipd.cores:
                    if initid == ipd.cores:
                        break
                else:
                    if initid:
                        break


    except KeyboardInterrupt as inst:
        ## ensure stdout is reset even if Exception was raised            
        #sys.stdout = save_stdout
        #sys.stderr = save_stderr
        raise inst

    ## This is raised if ipcluster is not running ------------
    except IOError as inst:
        ## ensure stdout is reset even if Exception was raised
        # sys.stdout = save_stdout
        # sys.stderr = save_stderr
        if "ipyrad-cli-" in ipd.cluster_id:
            raise IPyradWarningExit(NO_IPCLUSTER_CLI)
        else:
            raise IPyradWarningExit(NO_IPCLUSTER_API)

    except (ipp.TimeoutError, ipp.NoEnginesRegistered) as inst:
        ## raised by ipp if no connection file is found for 'nwait' seconds
        raise
        # sys.stdout = save_stdout
        # sys.stderr = save_stderr
        # raise inst

    # except Exception as inst:
    #     ## if any other exceptions were missed...
    #     sys.stdout = save_stdout
    #     sys.stderr = save_stderr
    #     raise inst

    finally:
        ## ensure that no matter what we reset the stds
        sys.stdout = save_stdout
        sys.stderr = save_stderr

    return ipyclient


## GLOBALS AND EXCEPTIONS

NO_IPCLUSTER_CLI = """\
    No ipcluster instance found. This may be a problem with your installation
    setup, or it could be that the cluster instance isn't firing up fast enough.
    This most often happens on cluster nodes. One solution is to launch
    ipcluster by hand and then pass the `--ipcluster` flag to ipyrad. See
    the docs for more info: http://ipyrad.rtfd.io/HPC_script.html
    """
NO_IPCLUSTER_API = """
    No ipcluster instance found. See documentation for the proper way to set 
    up an ipcluster instance when running the ipyrad Python API. In short, 
    you must run 'ipcluster start' to initiate a local or remote cluster. 
    Also, if you changed the 'profile' or 'cluster_id' setting from their 
    default values you must enter these into the Assembly._ipcluster dictionary.
    """
