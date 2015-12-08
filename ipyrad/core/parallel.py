#!/usr/bin/env ipython2

## imports for running ipcluster
from __future__ import print_function
import subprocess
import psutil
import atexit
#import time
import sys
import os

#import ipyparallel.nbextension.clustermanager as clustermanager
#from tornado import web

import logging
LOGGER = logging.getLogger(__name__)

# ## start and stop ipcluster
# CM = clustermanager.ClusterManager()

# ## check that u'default' profile exists
# PROFILES = [i['profile'] for i in CM.list_profiles()]
# assert u'default' in PROFILES, "ipyparallel is not properly configured."

# ## check if u'default ipcluster is already running
# RUNNING = CM.profile_info(u'default')['status']

# ## if it's stopped then start a new ipcluster
# if RUNNING == 'stopped':
#     print("in here")
#     print(CM.profile_info(u'default')['status'])    
#     CM.start_cluster(u'default')
#     atexit.register(LOGGER.debug('%s'), u'default')
#     print(CM.profile_info(u'default')['status'])
# else:
#     print("out here")

# ## launch ipcluster with profile=u'default'.
# try:
#     CM.check_profile(u'default')
# except web.HTTPError():
#     LOGGER.critical(u'jupyter profile not found')


## start ipcluster
def start(name, nproc, controller, quiet, delay):
    """ Start ipcluster """
    if nproc != None:
        nproc = str(psutil.cpu_count())
    standard = ["ipcluster", "start", 
                "--daemon", 
                "--delay="+delay,
                "--cluster-id="+name,
                "--controller="+controller,
                "-n", str(nproc)]

    try: 
        subprocess.check_call(" ".join(standard), 
         	                  shell=True, 
                              stderr=subprocess.STDOUT)
        LOGGER.info("%s connection to %s engines [%s]", controller, nproc, name)
        print("  ipyparallel setup: {} connection to {} engines\n".\
                 format(controller, nproc))

    except subprocess.CalledProcessError as inst:
        LOGGER.debug("ipcontroller already running.")
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



def ipcontroller_init(nproc=None, controller="Local", quiet=False):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. The controller option is 
    used to toggle between Local, MPI, PBS.
    """
    ## check if this pid already has a running cluster
    ipname = "ipyrad-"+str(os.getpid())
    start(ipname, nproc, controller, quiet, delay='1.0')
    #time.sleep(1)
    atexit.register(stop, ipname)
    return ipname


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

