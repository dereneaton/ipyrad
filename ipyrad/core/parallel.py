#!/usr/bin/env ipython2

## imports for running ipcluster
from __future__ import print_function

import subprocess
import psutil
import atexit
import sys
import os

# pylint: disable=W0212


import logging
LOGGER = logging.getLogger(__name__)


## start ipcluster
def start(data, quiet):
    """ Start ipcluster """

    ## select all cores if no entry for args.cores
    nproc = data._ipcluster["cores"]
    if not data._ipcluster["cores"]:
        nproc = str(psutil.cpu_count())

    ## open all ip views for MPI
    iparg = ""
    if "MPI" in data._ipcluster["engines"]:
        iparg = "ip='*' "

    ## make ipcluster arg call
    standard = ["ipcluster", "start", 
                "--daemon", 
                "--cluster-id="+data._ipcluster["id"],
                "--engines="+data._ipcluster["engines"],
                "--n="+str(nproc), 
                iparg]


    ## wrap ipcluster start
    try: 
        LOGGER.info(" ".join(standard))
        with open(os.devnull, 'w') as fnull:
            subprocess.Popen(standard, stderr=subprocess.STDOUT, 
                                       stdout=fnull).communicate()

        print("  ipyparallel setup: {} connection to {} Engines\n"\
              .format(data._ipcluster["engines"], nproc))

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



#global IPNAME
#IPNAME = "ipyrad-"+str(os.getpid())


# # def parallel(engines, controller="Local"):
# #     """
# #     The name is a unique id that keeps this __init__ of ipyrad distinct
# #     from interfering with other ipcontrollers. The controller option is 
# #     used to toggle between Local, MPI, PBS.
# #     """
# #     global __IPNAME__    
# #     print("Establishing {} connection.".format(controller))
# #     ipname = "ipyrad[id="+str(random.randint(1, 999))+"]"
# #     start(ipname, controller, delay="1.0")
# #     ## give engines time to connect... (longer?)    
# #     time.sleep(1)    
# #     atexit.register(stop, ipname)    
# #     __IPNAME__ = ipname



# def ipcontroller_set(controller="Local"):
#     """
#     The name is a unique id that keeps this __init__ of ipyrad distinct
#     from interfering with other ipcontrollers. The controller option is 
#     used to toggle between Local, MPI, PBS.
#     """
#     print("Establishing {} connection.".format(controller))
#     ipname = "ipyrad-"+str(os.getpid())
#     start(ipname, controller, delay="1.0")
#     ## give engines time to connect... (longer?)    
#     time.sleep(3)    
#     atexit.register(stop, ipname)    
#     return ipname  



# if __name__ == "__main__":

#     ## Start ipcluster and register exit call
#     NAME = "test"
#     start(NAME)
#     atexit.register(stop, NAME)

