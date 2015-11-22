#!/usr/bin/env ipython2

## imports for running ipcluster
from __future__ import print_function
import subprocess
import psutil
import atexit
import time
import os


## start ipcluster
def start(name, controller="Local", delay="1.0"):
    """ Start ipcluster """
    nproc = str(psutil.cpu_count())
    standard = ["ipcluster", "start", 
                "--daemon", 
                "--delay="+delay,
                "--cluster-id="+name,
                "--controller="+controller,
                "-n", str(nproc)]
                #"--IPClusterEngines.overwrite=True",
    try: 
        #print("[{}]".format(" ".join(standard)))        
        subprocess.check_call(" ".join(standard), 
    	                  shell=True, 
                          stderr=subprocess.STDOUT)
        print("ipyparallel setup: {} connection to {} engines.".\
              format(controller, nproc))

    except subprocess.CalledProcessError:
        print("error: ipcontroller could not connect to engines.")
        print(subprocess.STDOUT)



## decorated func for stopping. Does not need to be called?
def stop(cluster_id):
    """ stop ipcluster at sys.exit """
    print("Closing {} remote parallel engines:".format(cluster_id))
    stopcall = ["ipcluster", "stop", 
                "--cluster-id="+cluster_id]
    try:
        subprocess.check_call(" ".join(stopcall), shell=True)
    except subprocess.CalledProcessError:
        pass



def ipcontroller_init(controller="Local"):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. The controller option is 
    used to toggle between Local, MPI, PBS.
    """
    ipname = "ipyrad-"+str(os.getpid())
    start(ipname, controller, delay="1.0")
    ## give engines time to connect... TODO: make this smarter
    time.sleep(1)

    atexit.register(stop, ipname)
    return ipname    



# def parallel(engines, controller="Local"):
#     """
#     The name is a unique id that keeps this __init__ of ipyrad distinct
#     from interfering with other ipcontrollers. The controller option is 
#     used to toggle between Local, MPI, PBS.
#     """
#     global __IPNAME__    
#     print("Establishing {} connection.".format(controller))
#     ipname = "ipyrad[id="+str(random.randint(1, 999))+"]"
#     start(ipname, controller, delay="1.0")
#     ## give engines time to connect... (longer?)    
#     time.sleep(1)    
#     atexit.register(stop, ipname)    
#     __IPNAME__ = ipname



def ipcontroller_set(controller="Local"):
    """
    The name is a unique id that keeps this __init__ of ipyrad distinct
    from interfering with other ipcontrollers. The controller option is 
    used to toggle between Local, MPI, PBS.
    """
    print("Establishing {} connection.".format(controller))
    ipname = "ipyrad-"+str(os.getpid())
    start(ipname, controller, delay="1.0")
    ## give engines time to connect... (longer?)    
    time.sleep(3)    
    atexit.register(stop, ipname)    
    return ipname  



if __name__ == "__main__":

    ## Start ipcluster and register exit call
    NAME = "test"
    start(NAME)
    atexit.register(stop, NAME)

