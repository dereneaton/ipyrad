#!/usr/bin/env python

"utility functions for the analysis tools"

# py2/3 compat
from __future__ import print_function
#from builtins import range
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


# standard lib
import datetime
import time
import sys
import os

# third party
from numba import njit, prange
from ipyparallel import NoEnginesRegistered, TimeoutError, Client
from ..core.parallel import cluster_info, register_ipcluster
from ..assemble.utils import IPyradError, IPyradWarningExit



def get_parallel(self, ipyclient, show_cluster):
    "Start a parallel client with _ipcluster dict or start a new one."

    # connect to a running client or raise an error if not found.
    if not ipyclient:
        ipyclient = get_client(self)

    # print a message about the cluster status
    if show_cluster:
        cluster_info(ipyclient=ipyclient, spacer="")

    # store ipyclient engine pids so we can hard-interrupt
    # later if assembly is interrupted. Only stores pids of engines
    # that aren't busy at this moment, otherwise it would block.
    self._ipcluster["pids"] = {}
    for eid in ipyclient.ids:
        engine = ipyclient[eid]
        if not engine.outstanding:
            pid = engine.apply(os.getpid).get()
            self._ipcluster["pids"][eid] = pid

    # return new or existing ipyclient
    return ipyclient




def parallelize_run(self, run_kwargs, show_cluster=False, auto=False):
    """
    Takes an analysis tools object with an associated _ipcluster attribute
    dictionary and either launches an ipcluster instance or connects to a 
    running one. 
    """
    try:
        # launch ipcluster and get the parallel client with cluster-id
        if auto:
            if not self._ipcluster["cores"]:

                # set default to 4
                if not self._ipcluster["cores"]:
                    self._ipcluster["cores"] = 4

                # start ipcluster and attach ipyrad-cli cluster-id
                register_ipcluster(self)

        # find a running ipcluster instance and register client
        ipyclient = get_parallel(self, run_kwargs["ipyclient"], show_cluster)
        run_kwargs["ipyclient"] = ipyclient

        # call the run command
        self._run(**run_kwargs)


    # print the error and cleanup
    except KeyboardInterrupt:
        print("\nKeyboard Interrupt by user\n")

    except (IPyradWarningExit, IPyradError) as inst:
        print("\nEncountered an IPyradError:\n{}\n".format(inst))
        raise

    except Exception as inst:
        print("\nEncountered an unexpected error:\n{}\n".format(inst))
        raise

    # close client when done or interrupted
    finally:
        # unquiet 
        self.quiet = False

        try:
            # can't close client if it was never open
            if ipyclient:

                # send SIGINT (2) to all engines
                try:
                    ipyclient.abort()
                    time.sleep(1)
                    for engine_id, pid in self._ipcluster["pids"].items():
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                    time.sleep(1)
                except NoEnginesRegistered:
                    pass

                # shutdown the hub if it was auto-launched
                if auto:
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()

                # cleanup but keep alive
                else:
                    if not ipyclient.outstanding:
                        ipyclient.purge_everything()
                    else:
                        # nanny: kill everything, something bad happened
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        print("\nwarning: ipcluster was shutdown and must be restarted")
                
        # if exception is close and save, print and ignore
        except TypeError:
            # ipcluster was never started
            pass

        except Exception as inst2:
            print("warning: error during shutdown:\n{}".format(inst2))



def get_client(data):
    """ 
    Creates a client to view ipcluster engines for a given profile and 
    returns it with at least one engine spun up and ready to go. If no 
    engines are found after nwait amount of time then an error is raised.
    If engines==MPI it waits a bit longer to find engines. If the number
    of engines is set then it waits even longer to try to find that number
    of engines.
    """
    # save stds for later, we're gonna hide them to prevent external printing 
    save_stdout = sys.stdout 
    save_stderr = sys.stderr
    sys.stdout = StringIO()
    sys.stderr = StringIO()

    # wrapped search for ipcluster
    try: 
        args = {
            "profile": data._ipcluster["profile"],
            "timeout": data._ipcluster["timeout"],
            "cluster_id": data._ipcluster["cluster_id"],            
        }
        ipyclient = Client(**args)
        sys.stdout = save_stdout
        sys.stderr = save_stderr
        print("establishing parallel connection:")

        # allow time to find the connection; count cores to break
        for _ in range(6000):           
            initid = len(ipyclient)
            time.sleep(0.01)

            ## If MPI then wait for all engines to start so we can report
            ## how many cores are on each host. If Local then only wait for
            ## one engine to be ready and then just go.
            if (data._ipcluster["engines"] == "MPI") or ("ipyrad-cli-" in data._ipcluster["cluster_id"]):

                # wait for cores to be connected
                if data._ipcluster["cores"]:
                    time.sleep(0.1)
                    if initid == data._ipcluster["cores"]:
                        break

                if initid:
                    time.sleep(3)
                    if len(ipyclient) == initid:
                        break

            else:
                if data._ipcluster["cores"]:
                    if initid == data._ipcluster["cores"]:
                        break
                else:
                    if initid:
                        break

    except KeyboardInterrupt as inst:
        raise inst

    except (IOError, OSError):
        raise IPyradWarningExit("ipcluster not found, use 'auto=True'")

    except (TimeoutError, NoEnginesRegistered):
        raise IPyradWarningExit("ipcluster not found, use 'auto=True'")        

    finally:
        # no matter what we reset the stds
        sys.stdout = save_stdout
        sys.stderr = save_stderr

    return ipyclient








# parallel and prange give a xncores speedup?
njit(parallel=True)
def get_spans(maparr, spans):
    """ 
    Get span distance for each locus in original seqarray. This
    is used to create re-sampled arrays in each bootstrap to sample
    unlinked SNPs from. Used on snpsphy or str or ...
    """
    start = 0
    end = 0
    for idx in prange(1, spans.shape[0] + 1):
        lines = maparr[maparr[:, 0] == idx]
        if lines.size:
            end = lines[:, 3].max()
            spans[idx - 1] = [start, end]
        else: 
            spans[idx - 1] = [end, end]
        start = spans[idx - 1, 1]

    # drop rows with no span (invariant loci)
    spans = spans[spans[:, 0] != spans[:, 1]]
    return spans



def progressbar(finished, total, start, message):
    progress = 100 * (finished / float(total))
    hashes = '#' * int(progress / 5.)
    nohash = ' ' * int(20 - len(hashes))
    elapsed = datetime.timedelta(seconds=int(time.time() - start))
    print("\r[{}] {:>3}% {} | {:<12} "
        .format(hashes + nohash, int(progress), elapsed, message),
        end="")
    sys.stdout.flush()    


# New params class is iterable returning keys
class Params(object):
    "A dict-like object for storing params values with a custom repr"
    def __init__(self):
        self._i = 0

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __iter__(self):
        return self
    
    def __next__(self):
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if self._i > len(keys) - 1:
            self._i = 0
            raise StopIteration
        else:
            self._i += 1
            return keys[self._i - 1]
        
    def __repr__(self):
        "return simple representation of dict with ~ shortened for paths"
        _repr = ""
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if keys:
            _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
            for key in keys:
                _val = str(self[key]).replace(os.path.expanduser("~"), "~")
                _repr += _printstr.format(key, _val)
        return _repr

