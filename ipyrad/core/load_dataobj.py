#!/usr/bin/env python2.7

""" loads an archived Assembly object. """

from __future__ import print_function
import dill
import os


def load_dataobj(tryname):
    """ loads an ipython pickled Assembly object """
    if ".dataobj" not in tryname:
        tryname += ".dataobj"
    if os.path.exists(tryname):
        with open(tryname, "rb") as pickin:
            data = dill.load(pickin)
        return data
    else:
        print("cannot find", tryname, "please enter full path to file")

