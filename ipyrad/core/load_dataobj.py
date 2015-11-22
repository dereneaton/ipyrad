#!/usr/bin/env python2.7

""" loads an archived Assembly object. """

from __future__ import print_function

import os
import dill
from ipyrad.core.assembly import Assembly
from ipyrad.core.parallel import ipcontroller_init

def load_assembly(tryname, controller="Local"):
    """ loads an ipython pickled Assembly object """
    if ".assembly" not in tryname:
        tryname += ".assembly"

    if not os.path.exists(tryname):
        print("cannot find", tryname, "try entering the full path to file")

    else:
        ## load in the Assembly object
        with open(tryname, "rb") as pickin:
            data = dill.load(pickin)
        ## relaunch the ipcluster
        data.__ipname__ = ipcontroller_init(controller)


        return data


def save_dataobj():
	""" TODO: """
	pass





    ## iterate over Assembly objects





