#!/usr/bin/env python2.7

""" loads an archived Assembly object. """

from __future__ import print_function

import os
import dill
from ipyrad.core.parallel import ipcontroller_init


def load_assembly(name, controller="Local", quiet=False, launch=True):
    """ loads an ipython dill pickled Assembly object """
    ## flexible name entry
    locations = [name]
    locations.append(name+".assembly")

    ## does Assembly saved obj exist?
    for name in locations:
        try:
            ## load in the Assembly object
            with open(name, "rb") as pickin:
                data = dill.load(pickin)

            ## will raise Attribute error if not loaded
            fullcurdir = os.path.realpath(os.path.curdir)
            name = name.replace(fullcurdir, ".")
            if not quiet:
                print("  loading Assembly: {} [{}]".\
                      format(data.name, name))
    
            ## relaunch ipcluster
            if launch:
                data._ipclusterid = ipcontroller_init(nproc="",
                                                      controller=controller,
                                                      quiet=quiet)

        except (IOError, AttributeError):
            pass

    try:
        return data
    except UnboundLocalError:
        raise AssertionError(name)


def save_dataobj():
	""" TODO: """
	pass





    ## iterate over Assembly objects





