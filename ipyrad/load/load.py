#!/usr/bin/env python2.7

""" loads an archived Assembly object. """

from __future__ import print_function

import os
import sys
import dill
import time
from copy import deepcopy
from ipyrad.core.parallel import ipcontroller_init
from ipyrad.core.assembly import Assembly



def load_assembly(name, controller="Local", quiet=False, launch=False):
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

            ## Test if our assembly is currently up to date
            ## How to deal with assembly objects falling out of synch with the 
            ## currently assembly params in the code. If forceupdate is on
            ## then update the loaded assembly to the current version. If 
            ## forceupdate is not on, then test to check if the loaded
            ## assembly is current. If it is then fine, if not bail out with
            ## a hopefully useful error message.
            if test_assembly(data):
                print("  Attempting to update assembly to newest version.")
                data = update_assembly(data)

            ## relaunch ipcluster
            if launch:
                data._ipclusterid = ipcontroller_init(nproc="",
                                                      controller=controller,
                                                      quiet=quiet)
            else:
                data._ipclusterid = ""

        except (IOError, AttributeError):
            pass

    try:
        return data
    except UnboundLocalError:
        raise AssertionError("Attempting to load assembly. File not found: {}".format(name))


def test_assembly(data):
    """ Check to see if the assembly you're trying to load is concordant
        with the current assembly version. Basically it creates a new tmp
        assembly and tests whether the paramsdicts are the same. It also
        tests the _hackersonly dict."""

    new_assembly = Assembly(data.name, quiet=True)
    new_params = set(new_assembly.paramsdict.keys())

    my_params = set(data.paramsdict.keys())

    ## Find all params that are in the new paramsdict and not in the old one.
    params_diff = new_params.difference(my_params)

    result = False
    if params_diff:
        result = True

    ## Test hackersonly dict as well.
    my_hackerdict = set(data._hackersonly.keys())
    new_hackerdict = set(new_assembly._hackersonly.keys())
    hackerdict_diff = new_hackerdict.difference(my_hackerdict)

    if hackerdict_diff:
        result =  True

    return result


def update_assembly(data):
    """ Create a new Assembly() and convert as many of our old params to the new
        version as we can. Also report out any parameters that are removed
        and what their values are. 
    """

    print("##############################################################")
    print("Updating assembly to current version")
    ## New assembly object to update pdate from.
    new_assembly = Assembly("update", quiet=True)

    ## Hackersonly dict gets automatically overwritten
    ## Always use the current version for params in this dict.
    data._hackersonly = deepcopy(new_assembly._hackersonly)

    new_params = set(new_assembly.paramsdict.keys())

    my_params = set(data.paramsdict.keys())

    ## Find all params in loaded assembly that aren't in the new assembly.
    ## Make a new dict that doesn't include anything in removed_params
    removed_params = my_params.difference(new_params)
    for i in removed_params:
        print("Removing parameter: {} = {}".format(i, data.paramsdict[i]))
        
    ## Find all the params that are in the new paramsdict and not in the old one.
    ## If the set isn't emtpy then we create a new dictionary based on the new
    ## assembly parameters and populated with the currently loaded assembly values.
    ## Conditioning on not including any removed params. Magic.
    added_params = new_params.difference(my_params)
    for i in added_params:
        print("Adding parameter: {} = {}".format(i, new_assembly.paramsdict[i]))

    print("\nPlease take note of these changes. Every effort is made to\n"\
            +"ensure compatibility across versions of ipyrad. See online\n"\
            +"documentation for further details about new parameters.")
    time.sleep(5)
    print("##############################################################")
    
    if added_params:
        for i in data.paramsdict:
            if i not in removed_params:
                new_assembly.paramsdict[i] = data.paramsdict[i]
        data.paramsdict = deepcopy(new_assembly.paramsdict)

    data.save()
    return data

def save_dataobj():
	""" TODO: """
	pass





    ## iterate over Assembly objects





