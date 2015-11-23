#!/usr/bin/env python2.7

""" loads an archived Assembly object. """

from __future__ import print_function

import os
import dill
from ipyrad.core.assembly import Assembly
from ipyrad.core.parallel import ipcontroller_init


def load_assembly(name, controller="Local"):
    """ loads an ipython pickled Assembly object """
    ## flexible name entry
    if ".assembly" not in name:
        name += ".assembly"

    ## does Assembly save obj exist?
    if not os.path.exists(name):
        print("cannot find", name, "try entering the full path to file.")

    else:
        ## load in the Assembly object
        with open(name, "rb") as pickin:
            data = dill.load(pickin)
        ## relaunch ipcluster
        data.__ipname__ = ipcontroller_init(controller)

        return data


def save_dataobj():
	""" TODO: """
	pass





    ## iterate over Assembly objects





