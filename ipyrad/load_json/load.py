#!/usr/bin/env python

"loads an archived Assembly object."

from __future__ import print_function

import os
import json
import itertools
import pandas as pd
from ..core.sample import Sample
from ..core.assembly import Assembly
from ..assemble.utils import IPyradWarningExit, ObjDict, IPyradError



def load_json(path, quiet=False, cli=False):
    """ 
    Load a json serialized object and ensure it matches to the current 
    Assembly object format 
    """

    # load the JSON string and try with name+.json
    checkfor = [path + ".json", path]
    for inpath in checkfor:
        inpath = inpath.replace("~", os.path.expanduser("~"))
        try:
            with open(inpath, 'rb') as infile:
                fullj = json.loads(infile.read(), object_hook=tup_and_byte)
        except IOError:
            pass

    # create a new empty Assembly
    try:
        oldname = fullj["assembly"].pop("name")
        olddir = fullj["assembly"]["dirs"]["project"]
        oldpath = os.path.join(olddir, os.path.splitext(oldname)[0] + ".json")
        null = Assembly(oldname, quiet=True, cli=cli)

    except (UnboundLocalError, AttributeError) as inst:
        raise IPyradError("""
    Could not find saved Assembly file (.json) in expected location.
    Checks in: [project_dir]/[assembly_name].json
    Checked: {}
    """.format(inpath))

    # print msg with shortpath
    if not quiet:
        oldpath = oldpath.replace(os.path.expanduser("~"), "~")
        print("{}loading Assembly: {}".format(null._spacer, oldname))
        print("{}from saved path: {}".format(null._spacer, oldpath))

    # First get the samples. Create empty sample dict of correct length 
    samplekeys = fullj["assembly"].pop("samples")
    null.samples = {name: "" for name in samplekeys}

    # Next get paramsdict and use set_params to convert values back to 
    # the correct dtypes. Allow set_params to fail because the object will 
    # be subsequently updated by the params from the params file, which may
    # correct any errors/incompatibilities in the old params file
    oldparams = fullj["assembly"].pop("paramsdict")
    for param, val in oldparams.items():
        # backwards compatiblity with <0.8.0
        param = param.lstrip("_")
        if param not in ["assembly_name", "keys"]:
            try:
                null.set_params(param, val)
            except IPyradWarningExit as inst:
                pass

    # Import the hackersonly dict. In this case we don't have the nice
    # set_params so we're shooting from the hip to reset the values
    try:
        oldhackersonly = fullj["assembly"].pop("_hackersonly")
        for param, val in oldhackersonly.items():
            if val is None:
                null._hackersonly[param] = None
            else:
                null._hackersonly[param] = val

    except Exception as inst:
        pass
    #     LOGGER.warning("""
    # Load assembly error resetting hackersonly dict element. We will just use
    # the default value in the current assembly.""")

    # Check remaining attributes of Assembly and Raise warning if attributes
    # do not match up between old and new objects
    newkeys = list(null.__dict__.keys())
    oldkeys = list(fullj["assembly"].keys())

    ## find shared keys and deprecated keys
    sharedkeys = set(oldkeys).intersection(set(newkeys))
    # lostkeys = set(oldkeys).difference(set(newkeys))
    # # raise warning if there are lost/deprecated keys
    # if lostkeys:
    #     LOGGER.warning("""
    # load_json found {a} keys that are unique to the older Assembly.
    #     - assembly [{b}] v.[{c}] has: {d}
    #     - current assembly is v.[{e}]
    #     """.format(a=len(lostkeys), 
    #                b=oldname,
    #                c=fullj["assembly"]["_version"],
    #                d=lostkeys,
    #                e=null._version))

    # load in remaining shared Assembly attributes to null
    for key in sharedkeys:
        setattr(null, key, fullj["assembly"][key])
        null.__setattr__(key, fullj["assembly"][key])

    # ## load in svd results if they exist
    # try:
    #     if fullj["assembly"]["svd"]:
    #         null.__setattr__("svd", fullj["assembly"]["svd"])
    #         null.svd = ObjDict(null.svd)
    # except Exception:
    #     LOGGER.debug("skipping: no svd results present in old assembly")

    # Now, load in the Sample objects json dicts
    sample_names = list(fullj["samples"].keys())
    if not sample_names:
        raise IPyradWarningExit("""
    No samples found in saved assembly. If you are just starting a new
    assembly the file probably got saved erroneously, so it's safe to try 
    removing the assembly file (e.g., rm {}.json) and restarting.

    If you fully completed step 1 and you see this message you should probably
    contact the developers.
    """.format(inpath))
        
    sample_keys = list(fullj["samples"][sample_names[0]].keys())
    stats_keys = list(fullj["samples"][sample_names[0]]["stats"].keys())
    stats_dfs_keys = list(fullj["samples"][sample_names[0]]["stats_dfs"].keys())
    ind_statkeys = (
        [fullj["samples"][sample_names[0]]["stats_dfs"][i].keys() 
        for i in stats_dfs_keys])
    ind_statkeys = list(itertools.chain(*ind_statkeys))

    # check against a null sample
    nsamp = Sample()
    newkeys = list(nsamp.__dict__.keys())
    newstats = list(nsamp.__dict__["stats"].keys())
    newstatdfs = list(nsamp.__dict__["stats_dfs"].keys())
    newindstats = [
        nsamp.__dict__["stats_dfs"][i].keys() for i in newstatdfs]
    newindstats = list(itertools.chain(*[i.values for i in newindstats]))

    # different in attributes?
    # diffattr = set(sample_keys).difference(newkeys)
    # diffstats = set(stats_keys).difference(newstats)
    # diffindstats = set(ind_statkeys).difference(newindstats)
    # Raise warning if any oldstats were lost or deprecated
    # alldiffs = diffattr.union(diffstats).union(diffindstats)
    # if any(alldiffs):
    #     LOGGER.warning("""
    # load_json found {a} keys that are unique to the older Samples.
    #     - assembly [{b}] v.[{c}] has: {d}
    #     - current assembly is v.[{e}]
    #     """.format(a=len(alldiffs), 
    #                b=oldname,
    #                c=fullj["assembly"]["_version"],
    #                d=alldiffs,
    #                e=null._version))

    # save stats and statsfiles to Samples
    for sample in null.samples:
        # create a null Sample
        null.samples[sample] = Sample()

        # save stats
        sdat = fullj["samples"][sample]['stats']
        # Reorder the keys so they ascend by step, only include
        # stats that are actually in the sample. newstats is a
        # list of the new sample stat names, and stats_keys
        # are the names of the stats from the json file.
        newstats = [x for x in newstats if x in stats_keys]
        null.samples[sample].stats = pd.Series(sdat).reindex(newstats)

        # save stats_dfs
        for statskey in stats_dfs_keys:
            null.samples[sample].stats_dfs[statskey] = (
                pd.Series(fullj["samples"][sample]["stats_dfs"][statskey])
                .reindex(
                    list(nsamp.__dict__["stats_dfs"][statskey].keys())
                    )
                )

        # save Sample files
        for filehandle in fullj["samples"][sample]["files"].keys():
            null.samples[sample].files[filehandle] = (
                fullj["samples"][sample]["files"][filehandle])


    # build the Assembly object stats_dfs
    for statskey in stats_dfs_keys:
        indstat = null._build_stat(statskey)
        if not indstat.empty:
            null.stats_dfs[statskey] = indstat

    ## add remaning attributes to null Samples
    shared_keys = set(sample_keys).intersection(newkeys)
    shared_keys.discard("stats")
    shared_keys.discard("files")    
    shared_keys.discard("stats_files")
    shared_keys.discard("stats_dfs")

    for sample in null.samples:
        ## set the others
        for key in shared_keys:
            null.samples[sample].__setattr__(key, fullj["samples"][sample][key])

    ## ensure objects are object dicts
    null.dirs = ObjDict(null.dirs)
    null.stats_files = ObjDict(null.stats_files)
    null.stats_dfs = ObjDict(null.stats_dfs)    
    null.populations = ObjDict(null.populations)
    null.outfiles = ObjDict(null.outfiles)

    return null




def tup_and_byte(obj):
    """ this is used in loading """

    # convert all strings to bytes
    if isinstance(obj, (bytes)):
        return obj.decode()  # encode('utf-8')
        #return obj.encode('utf-8')

    # if this is a list of values, return list of byteified values
    if isinstance(obj, list):
        return [tup_and_byte(item) for item in obj]

    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(obj, dict):
        if "__tuple__" in obj:
            return tuple(tup_and_byte(item) for item in obj["items"])
        else:
            return {
                tup_and_byte(key): tup_and_byte(val) for
                key, val in obj.items()
                }

    # if it's anything else, return it in its original form
    return obj
