#!/usr/bin/env ipython2

## imports for running ipcluster
from __future__ import print_function

import subprocess
import atexit
import shlex
import sys
import os
from ipyrad.assemble.util import detect_cpus

# pylint: disable=W0212


import logging
LOGGER = logging.getLogger(__name__)


## start ipcluster
def start(data, quiet):
    """ Start ipcluster """

    ## open all ip views for MPI
    iparg = ""
    if "MPI" in data._ipcluster["engines"]:
        iparg = "--ip='*' "

    ## make ipcluster arg call
    standard = """
        ipcluster start 
                  --daemon 
                  --cluster-id={}
                  --engines={} 
                  --n={}
                  {}"""\
        .format(data._ipcluster["id"], 
                data._ipcluster["engines"], 
                data.cpus,
                iparg)
                   
    ## wrap ipcluster start
    try: 
        LOGGER.info(shlex.split(standard))
        subprocess.check_output(shlex.split(standard))

        print("  ipyparallel setup: {} connection to {} Engines\n"\
              .format(data._ipcluster["engines"], data.cpus))

    except subprocess.CalledProcessError as inst:
        LOGGER.debug("ipcontroller already running.")
        raise

    except Exception as inst:
        sys.exit("Error launching ipcluster for parallelization:\n({})\n".\
                 format(inst))



## decorated func for stopping. Does not need to be called?
def stop(cluster_id):
    """ stop ipcluster at sys.exit """
    #print("\nclosing remote Engines")
    LOGGER.info("Shutting down [%s] remote Engines", cluster_id)
    stopcall = ["ipcluster", "stop", 
                "--cluster-id="+cluster_id]
    try:
        subprocess.check_call(" ".join(stopcall), 
                              shell=True, 
                              stderr=subprocess.STDOUT,
                              stdout=subprocess.PIPE)
    except subprocess.CalledProcessError:
        pass



def ipcontroller_init(data, quiet=False):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. The controller option is 
    used to toggle between Local, MPI, PBS.
    """
    ## check if this pid already has a running cluster
    data._ipcluster["id"] = "ipyrad-"+str(os.getpid())

    start(data, quiet)
    atexit.register(stop, data._ipcluster["id"])

    return data
