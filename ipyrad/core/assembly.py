#!/usr/bin/env python

""" ipyrad Assembly class object. """

from __future__ import print_function
import os
import glob
import sys
import copy
import time
import json
import string
import datetime
import numpy as np
import pandas as pd
import ipyrad as ip

from collections import OrderedDict
from ipyrad.assemble.utils import IPyradParamsError, IPyradError
from ipyrad.assemble.utils import ObjDict
from ipyrad.core.paramsinfo import paraminfo, paramname
from ipyrad.core.Parallel import Parallel
from ipyrad.core.params import Params, Hackers
# from ipyrad.core.parallel import register_ipcluster


    

class Assembly(object):
    """ 
    An ipyrad Assembly class object.

    The core object in ipyrad used to store and retrieve results, to
    call assembly functions, and to link to Sample objects.

    Parameters
    ----------
    name : str
        A name should be passed when creating a new Assembly object.
        This name will be used as a prefix for all files saved to disk
        associated with this Assembly. It is automatically set as the
        prefix name (parameter 14).


    Attributes
    ----------
    name : str
        A name for the Assembly object. Used for all saved files on disk.
    samples : dict
        Returns a dict with Sample names as keys and Sample objects as values.
    barcodes : dict
        Returns a dictionary with Sample names as keys and barcodes as values.
        The barcodes information is fetched from parameter 3
        '[Assembly].params.barcodes_path'

    Functions
    ---------
    run(step, force, ipyclient)
        runs step x of assembly
    write_params(filehandle, force)
        writes params dictionary to params.txt file format


    Returns
    -------
    object
         A new assembly object is returned.
    """
    def __init__(self, name, **kwargs):

        # If using CLI then "cli" is included in kwargs
        self._cli = False
        if kwargs.get("cli"):
            self._cli = True           

        # this is False and only updated during .run()
        self.quiet = False

        # No special characters in assembly name
        check_name(name)      
        self.name = name
        if (not kwargs.get("quiet")) and (not self._cli):
            self._print("New Assembly: {}".format(self.name))

        # Default ipcluster launch info
        self.ipcluster = {
            "cluster_id": "",
            "profile": "default",
            "engines": "Local",
            "quiet": 0,
            "timeout": 120,
            "cores": 0,  # detect_cpus(),
            "threads": 2,
            "pids": {},
        }
        # ipcluster settings can be set during init using kwargs
        for key, val in kwargs.items():
            if key in self.ipcluster:
                self.ipcluster[key] = val

        # statsfiles is a dict with file locations
        # stats_dfs is a dict with pandas dataframes
        self.stats_files = ObjDict({})
        self.stats_dfs = ObjDict({})

        # samples linked {sample-name: sample-object}
        self.samples = {}

        # populations {popname: poplist}
        self.populations = {}

        # multiplex files linked
        self.barcodes = {}

        # outfiles locations
        self.outfiles = ObjDict()
        self.outfiles.loci = ""

        # storing supercatg file
        self.clust_database = ""
        self.database = ""

        ## the default params dict
        self.params = Params(self)
        self.hackersonly = Hackers()

        ## Store data directories for this Assembly. Init with default project
        self.dirs = ObjDict({
            "project": os.path.realpath(self.params.project_dir),
            "fastqs": "",
            "edits": "",
            "clusts": "",
            "consens": "",
            "across": "",
            "outfiles": "",
        })


    def __str__(self):
        return "<ipyrad.Assembly object {}>".format(self.name)

    @property
    def _spacer(self):
        """ return print spacer for CLI versus API """
        if self._cli:
            return "  "
        return ""


    @property
    def stats(self):
        """ Returns a data frame with Sample data and state. """
        nameordered = list(self.samples.keys())
        nameordered.sort()

        ## Set pandas to display all samples instead of truncating
        pd.options.display.max_rows = len(self.samples)
        statdat = pd.DataFrame(
            data=[self.samples[i].stats for i in nameordered],
            index=nameordered,
        ).dropna(axis=1, how='all')
        # ensure non h,e columns print as ints
        for column in statdat:
            if column not in ["hetero_est", "error_est"]:
                statdat[column] = np.nan_to_num(statdat[column]).astype(int)

        # build step 6-7 stats from database
        # ...

        return statdat


    @property
    def files(self):
        """ Returns a data frame with Sample files. Not very readable... """
        nameordered = list(self.samples.keys())
        nameordered.sort()
        ## replace curdir with . for shorter printing
        sdf = pd.DataFrame(
            data=[self.samples[i].files for i in nameordered],
            index=nameordered,
        ).dropna(axis=1, how='all')
        return sdf


    def save(self):
        save_json(self)


    def _print(self, value):
        if not self.quiet:
            print("{}{}".format(self._spacer, value))


    def _progressbar(self, njobs, finished, start, msg):

        # bail
        if self.quiet:
            return
        # measure progress
        if njobs:
            progress = 100 * (finished / float(njobs))
        else:
            progress = 100

        # build the bar
        hashes = '#' * int(progress / 5.)
        nohash = ' ' * int(20 - len(hashes))

        # timestamp
        elapsed = datetime.timedelta(seconds=int(time.time() - start))

        # print to stderr
        if self._cli:
            print("\r{}[{}] {:>3}% {} | {:<12} ".format(
                self._spacer,
                hashes + nohash,
                int(progress),
                elapsed,
                msg[0],
            ), end="")
        else:
            print("\r{}[{}] {:>3}% {} | {:<12} | {} |".format(*[
                self._spacer,
                hashes + nohash,
                int(progress),
                elapsed,
                msg[0],
                msg[1],
            ]), end="")
        sys.stdout.flush()


    def _build_stat(self, idx):
        """ 
        Returns an Assembly stats data frame rom Sample stats data frames.
        e.g., data.stats_dfs.s1 = self.build_stats("s1")
        """
        # build steps 1-5 stats from samples
        nameordered = list(self.samples.keys())
        nameordered.sort()
        newdat = pd.DataFrame(
            (self.samples[i].stats_dfs[idx] for i in nameordered),
            index=nameordered,
        ).dropna(axis=1, how='all')
        return newdat


    def _link_barcodes(self):
        """
        Parses Sample barcodes to a dictionary from 'barcodes_path'. This 
        function is called whenever a barcode-ish param is changed. 
        """
        # find barcodefile
        barcodefile = glob.glob(self.params.barcodes_path)
        if not barcodefile:
            raise IPyradError(
                "Barcodes file not found. You entered: {}"
                .format(self.params.barcodes_path))

        # read in the file
        bdf = pd.read_csv(barcodefile[0], header=None, delim_whitespace=1)
        bdf = bdf.dropna()

        # make sure bars are upper case
        bdf[1] = bdf[1].str.upper()

        # if replicates are present then print a warning
        if bdf[0].value_counts().max() > 1:
            self._print("Warning: technical replicates (same name) present.")

            # adds -technical-replicate-N to replicate names (NON_DEFAULT)
            # if not self.hackersonly.merge_technical_replicates:                   
            repeated = (bdf[0].value_counts() > 1).index
            for rep in repeated:
                farr = bdf[bdf[0] == rep]
                for idx, index in enumerate(farr.index):
                    bdf.loc[index, 0] = (
                        "{}-technical-replicate-{}".format(rep, idx))
                            
        # make sure chars are all proper
        if not all(bdf[1].apply(set("RKSYWMCATG").issuperset)):
            raise IPyradError(BAD_BARCODE)

        # store barcodes as a dict
        self.barcodes = dict(zip(bdf[0], bdf[1]))

        # 3rad/seqcap use multiplexed barcodes
        if "3rad" in self.params.datatype:
            if not bdf.shape[1] == 3:
                raise IPyradError(
                    "pair3rad datatype should have two barcodes per sample.")
        
            # We'll concatenate them with a plus and split them later
            bdf[2] = bdf[2].str.upper()
            self.barcodes = dict(zip(bdf[0], bdf[1] + "+" + bdf[2]))               


    def _link_populations(self, popdict=None, popmins=None):
        """
        Creates self.populations dictionary to save mappings of individuals to
        populations/sites, and checks that individual names match with Samples.
        The self.populations dict keys are pop names and the values are lists
        of length 2. The first element is the min number of samples per pop
        for final filtering of loci, and the second element is the list of
        samples per pop.

        Population assigments are used for heirarchical clustering, for
        generating summary stats, and for outputing some file types (.treemix
        for example). Internally stored as a dictionary.

        Note
        ----
        By default a File is read in from `pop_assign_file` with one individual
        per line and space separated pairs of ind pop:

            ind1 pop1
            ind2 pop2
            ind3 pop3
            etc...

        Parameters
        ----------
        popdict : dict
            When using the API it may be easier to simply create a dictionary
            to pass in as an argument instead of reading from an input file.
            This can be done with the `popdict` argument like below:

            pops = {'pop1': ['ind1', 'ind2', 'ind3'], 'pop2': ['ind4', 'ind5']}
            [Assembly]._link_populations(popdict=pops).

        popmins : dict
            If you want to apply a minsamples filter based on populations
            you can add a popmins dictionary. This indicates the number of 
            samples in each population that must be present in a locus for 
            the locus to be retained. Example:

            popmins = {'pop1': 3, 'pop2': 2}

        """
        if not popdict:
            ## glob it in case of fuzzy matching
            popfile = glob.glob(self.params.pop_assign_file)[0]
            if not os.path.exists(popfile):
                raise IPyradError(
                    "Population assignment file not found: {}"
                    .format(self.params.pop_assign_file))

            try:
                ## parse populations file
                popdat = pd.read_csv(
                    popfile, header=None,
                    delim_whitespace=1,
                    names=["inds", "pops"], 
                    comment="#")

                popdict = {
                    key: group.inds.values.tolist() for key, group in
                    popdat.groupby("pops")}

                ## parse minsamples per population if present (line with #)
                mindat = [
                    i.lstrip("#").lstrip().rstrip() for i in 
                    open(popfile, 'r').readlines() if i.startswith("#")]

                if mindat:
                    popmins = {}
                    for i in range(len(mindat)):
                        minlist = mindat[i].replace(",", "").split()
                        popmins.update({i.split(':')[0]: int(i.split(':')[1])
                                        for i in minlist})
                else:
                    raise IPyradError(NO_MIN_SAMPLES_PER_POP)

            except (ValueError, IOError):
                raise IPyradError(
                    "  Populations file malformed - {}".format(popfile))

        else:
            ## pop dict is provided by user
            if not popmins:
                popmins = {i: 1 for i in popdict}

        ## check popdict. Filter for bad samples
        ## Warn user but don't bail out, could be setting the pops file
        ## on a new assembly w/o any linked samples.
        # badsamples = [
            # i for i in itertools.chain(*popdict.values())
            # if i not in self.samples.keys()]

        # if any(badsamples):
        #     ip.logger.warn(
        #         "Some names from population input do not match Sample "\
        #       + "names: ".format(", ".join(badsamples)))
        #     ip.logger.warn("If this is a new assembly this is normal.")

        ## check popmins
        ## cannot have higher min for a pop than there are samples in the pop
        popmax = {i: len(popdict[i]) for i in popdict}
        if not all([popmax[i] >= popmins[i] for i in popdict]):
            raise IPyradError(
                " minsample per pop value cannot be greater than the " +
                " number of samples in the pop. Modify the populations file.")

        ## return dict
        self.populations = {i: (popmins[i], popdict[i]) for i in popdict}


    def get_params(self, param=""):
        "Pretty prints parameter settings to stdout"
        return self.params


    def set_params(self, param, newvalue):
        """
        Set a parameter to a new value in a verbose way. 
        Raises error if newvalue is wrong type.
        Use .get_params() or .params to see current assembly parameters. 

        Parameters
        ----------
        param : int or str
            The index (e.g., 1) or string name (e.g., "project_dir")
            for the parameter that will be changed.

        newvalue : int, str, or tuple
            The new value for the parameter selected for `param`. Use
            `ipyrad.get_params_info()` to get further information about
            a given parameter. If the wrong type is entered for newvalue
            (e.g., a str when it should be an int), an error will be raised.
            Further information about each parameter is also available
            in the documentation.

        Examples
        --------
        ## param 'project_dir' takes only a str as input
        [Assembly].set_params('project_dir', 'new_directory')

        ## param 'restriction_overhang' must be a tuple or str, if str it is
        ## converted to a tuple with the second entry empty.
        [Assembly].set_params('restriction_overhang', ('CTGCAG', 'CCGG')

        ## param 'max_shared_Hs_locus' can be an int or a float:
        [Assembly].set_params('max_shared_Hs_locus', 0.25)

        ## Simpler alternative: set attribute directly
        [Assembly].params.max_shared_Hs_locus = 0.25
        """
        # check param in keys
        if "_" + param not in self.params._keys:
            raise IPyradParamsError(
                "Parameter key not recognized: {}".format(param))
        
        # set parameter newvalue
        setattr(self.params, param, newvalue)


    def write_params(self, outfile=None, force=False):
        """ 
        Write out the parameters of this assembly to a file properly
        formatted as input for `ipyrad -p <params.txt>`. A good and
        simple way to share/archive parameter settings for assemblies.
        This is also the function that's used by __main__ to
        generate default params.txt files for `ipyrad -n`
        """
        if outfile is None:
            outfile = "params-{}.txt".format(self.name)

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise IPyradError(PARAMS_EXISTS.format(outfile))

        with open(outfile, 'w') as paramsfile:
            ## Write the header. Format to 80 columns
            header = "------- ipyrad params file (v.{})".format(ip.__version__)
            header += ("-" * (80 - len(header)))
            paramsfile.write(header)

            ## Whip through the current params and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key in self.params._keys:
                val = getattr(self.params, key)

                # If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                padding = (" " * (30 - len(paramvalue)))
                paramkey = self.params._keys.index(key)
                paramindex = " ## [{}] ".format(paramkey)
                name = "[{}]: ".format(paramname(paramkey))
                description = paraminfo(paramkey, short=True)
                paramsfile.write(
                    "\n" + paramvalue + padding + paramindex + name + description)


    def branch(self, newname, subsamples=None, infile=None, force=False):
        """
        Returns a copy of the Assembly object. Does not allow Assembly
        object names to be replicated in namespace or path.
        """
        # subsample by removal or keeping.
        remove = 0

        # does json file already exist?
        exists = os.path.exists(os.path.join(
            self.params.project_dir, 
            "{}.json".format(newname)
        ))

        # print warning and return if name exists
        if (newname == self.name) or (exists and not force):
            self._print(
                "Assembly object named {} already exists".format(newname))
            return 

        # else do the branching
        else:
            # Make sure the new name doesn't have any wacky characters
            check_name(newname)

            # Bozo-check. Carve off 'params-' if it's in the new name.
            if newname.startswith("params-"):
                newname = newname.split("params-")[1]

            # create a copy of the Assembly obj
            newobj = copy.deepcopy(self)
            newobj.name = newname
            newobj.params._assembly_name = newname

            # warn user to only use one of these at a time
            if subsamples and infile:
                print(BRANCH_NAMES_AND_INPUT)

            # parse infile 
            if infile:
                if infile[0] == "-":
                    remove = 1
                    infile = infile[1:]
                if os.path.exists(infile):
                    subsamples = read_sample_names(infile)

            ## if remove then swap the samples
            if remove:
                subsamples = list(set(self.samples.keys()) - set(subsamples))

            ## create copies of each subsampled Sample obj
            if subsamples:
                for sname in subsamples:
                    if sname in self.samples:
                        newobj.samples[sname] = copy.deepcopy(self.samples[sname])
                    else:
                        if sname != "reference":
                            print("Sample name not found: {}".format(sname))

                ## reload sample dict w/o non subsamples
                newobj.samples = {
                    name: sample for name, sample in newobj.samples.items() 
                    if name in subsamples}

            ## create copies of each subsampled Sample obj
            else:
                for sample in self.samples:
                    newobj.samples[sample] = copy.deepcopy(self.samples[sample])

            ## save json of new obj and return object
            newobj.save()
            return newobj


    def _compatible_params_check(self):
        "check params that must be compatible at run time"

        # do not allow statistical < majrule
        val1 = self.params.mindepth_statistical
        val2 = self.params.mindepth_majrule
        if val1 < val2:
            raise IPyradParamsError(
                "mindepth_statistical cannot be < mindepth_majrule")
        # other params to check ...


    def run(self, 
        steps=None,
        force=False,
        ipyclient=None, 
        quiet=False,
        show_cluster=False, 
        auto=False):
        """

        """
        # save assembly at state of run start
        self.save()
        
        # hide all messages/progress bars       
        self.quiet = quiet

        # check that mindepth params are compatible, fix and report warning.
        self._compatible_params_check()

        # distribute filling jobs in parallel
        pool = Parallel(
            tool=self,
            rkwargs={"steps": steps, "force": force},
            ipyclient=ipyclient,
            show_cluster=show_cluster,
            auto=auto,
            )
        pool.wrap_run()



    def _run(
        self, 
        steps=None, 
        force=False, 
        ipyclient=None, 
        show_cluster=False,
        auto=False,
        ):
        """
        Run assembly steps (1-7) of an ipyrad analysis.
        
        Parameters:
        ===========
        steps: (str, default=None)
            The steps of assembly to run, e.g., "123", "1234567".
        force: (bool, default=False)
            Whether to overwrite an existing assembly with the same name.
        ipyclient: (obj, default=None)
            An ipyparallel.Client() object to tune parallelization. See
            docs for details. Or, use auto=True. 
        quiet: (bool, default=False)
            Print progress information to stdout.
        show_cluster: (bool, default=False)
            Print parallelization information to stdout.
        auto: (bool, default=False)
            Automatically launch an ipcluster instance for parallelization 
            of this run and shut it down when finished. 
        """
        # function dictionary
        stepdict = {
            "1": ip.assemble.demultiplex.Step1,
            "2": ip.assemble.rawedit.Step2, 
            "3": ip.assemble.clustmap.Step3,
            "4": ip.assemble.jointestimate.Step4, 
            "5": ip.assemble.consens_se.Step5, 
            "6": ip.assemble.clustmap_across.Step6, 
            "7": ip.assemble.write_outputs.Step7,
        }
          
        # run step fuctions and save and clear memory after each
        for step in steps:
            stepdict[step](self, force, ipyclient).run()
            self.save()
            ipyclient.purge_everything()





class Encoder(json.JSONEncoder):
    """ 
    Save JSON string with tuples embedded as described in stackoverflow
    thread. Modified here to include dictionary values as tuples.
    link: http://stackoverflow.com/questions/15721363/

    This Encoder Class is used as the 'cls' argument to json.dumps()
    """
    def encode(self, obj):
        """ function to encode json string"""
        def hint_tuples(item):
            """ embeds __tuple__ hinter in json strings """
            if isinstance(item, tuple):
                return {'__tuple__': True, 'items': item}
            if isinstance(item, list):
                return [hint_tuples(e) for e in item]
            if isinstance(item, dict):
                return {
                    key: hint_tuples(val) for key, val in item.items()
                    }
            else:
                return item
        return super(Encoder, self).encode(hint_tuples(obj))



def default(o):
    # https://stackoverflow.com/questions/11942364/
    # typeerror-integer-is-not-json-serializable-when-
    # serializing-json-in-python?utm_medium=organic&utm_
    # source=google_rich_qa&utm_campaign=google_rich_qa
    if isinstance(o, np.int64): 
        return int(o)  
    raise TypeError


def save_json(data):
    """ 
    Save assembly and samples as json 
    ## data as dict
    #### skip _ipcluster because it's made new
    #### statsfiles save only keys
    #### samples save only keys
    """
    # store params without the reference to Assembly object in params
    paramsdict = data.params.__dict__
    paramsdict = {i: j for (i, j) in paramsdict.items() if i != "_data"}

    # store all other dicts
    datadict = OrderedDict([
        ("name", data.__dict__["name"]), 
        ("dirs", data.__dict__["dirs"]),
        ("paramsdict", paramsdict),
        ("samples", list(data.__dict__["samples"].keys())),
        ("populations", data.__dict__["populations"]),
        ("database", data.__dict__["database"]),
        ("clust_database", data.__dict__["clust_database"]),        
        ("outfiles", data.__dict__["outfiles"]),
        ("barcodes", data.__dict__["barcodes"]),
        ("stats_files", data.__dict__["stats_files"]),
        ("hackersonly", data.hackersonly._data),
    ])

    ## sample dict
    sampledict = OrderedDict([])
    for key, sample in data.samples.items():
        sampledict[key] = sample._to_fulldict()

    ## json format it using cumstom Encoder class
    fulldumps = json.dumps({
        "assembly": datadict,
        "samples": sampledict
    },
        cls=Encoder,
        sort_keys=False, indent=4, separators=(",", ":"),
        default=default,
    )

    ## save to file
    assemblypath = os.path.join(data.dirs.project, data.name + ".json")
    if not os.path.exists(data.dirs.project):
        os.mkdir(data.dirs.project)
    
    ## protect save from interruption
    done = 0
    while not done:
        try:
            with open(assemblypath, 'w') as jout:
                jout.write(fulldumps)
            done = 1
        except (KeyboardInterrupt, SystemExit): 
            print('.')
            continue


def merge(name, assemblies, rename_dict=None):
    """
    Creates and returns a new Assembly object in which samples from two or more
    Assembly objects with matching names are 'merged'. Merging does not affect 
    the actual files written on disk, but rather creates new Samples that are 
    linked to multiple data files, and with stats summed.

    # merge two assemblies
    new = ip.merge('newname', (assembly1, assembly2))

    # merge two assemblies and rename samples
    rename = {"1A_0", "A", "1B_0", "A"}
    new = ip.merge('newname', (assembly1, assembly2), rename_dict=rename)
    """
    # null rename dict if empty
    if not rename_dict:
        rename_dict = {}

    # create new Assembly
    merged = Assembly(name)

    # one or multiple assemblies?
    try:
        _ = len(assemblies)
    except TypeError:
        assemblies = [assemblies]

    # inherit params setting from first assembly
    for key in assemblies[0].params._keys[5:]:
        value = getattr(assemblies[0].params, key)
        setattr(merged.params, key, value)

    # iterate over all sample names from all Assemblies
    for data in assemblies:

        # make a deepcopy
        ndata = copy.deepcopy(data)
        for sname, sample in ndata.samples.items():

            # rename sample if in rename dict
            if sname in rename_dict:
                sname = rename_dict[sname]
                sample.name = sname

            # is it in the merged assembly already
            if sname in merged.samples:
                msample = merged.samples[sname]

                # update stats
                msample.stats.reads_raw += sample.stats.reads_raw
                if sample.stats.reads_passed_filter:
                    msample.stats.reads_passed_filter += \
                    sample.stats.reads_passed_filter

                # append files
                if sample.files.fastqs:
                    msample.files.fastqs += sample.files.fastqs
                if sample.files.edits:
                    msample.files.edits += sample.files.edits

                # do not allow state >2 at merging (requires reclustering)
                # if merging WITHIN samples.
                msample.stats.state = min(sample.stats.state, 2)

            # merge its stats and files
            else:
                merged.samples[sname] = sample

    # Merged assembly inherits max of hackers values (max frag length)
    merged.hackersonly.max_fragment_length = max(
        [i.hackersonly.max_fragment_length for i in assemblies])

    # Set the values for some params that don't make sense inside mergers
    merged_names = ", ".join([i.name for i in assemblies])
    merged.params.raw_fastq_path = "Merged: " + merged_names
    merged.params.barcodes_path = "Merged: " + merged_names
    merged.params.sorted_fastq_path = "Merged: " + merged_names

    # return the new Assembly object
    merged.save()
    return merged



def check_name(name):
    invalid_chars = (
        string.punctuation.replace("_", "").replace("-", "") + " ")
    if any(char in invalid_chars for char in name):
        raise IPyradParamsError(BAD_ASSEMBLY_NAME.format(name))


def read_sample_names(fname):
    """ 
    Read in sample names from a text file, a convenience function for branching
    """
    try:
        with open(fname, 'r') as infile:
            subsamples = [x.split()[0] for x in infile.readlines() if x.strip()]

    except Exception as inst:
        print("Failed to read input file with sample names.\n{}".format(inst))
        raise inst

    return subsamples


### ERROR MESSAGES ###################################
UNKNOWN_EXCEPTION = """\
{}Encountered an unexpected error (see ./ipyrad_log.txt)"+\
{}Error message is below -------------------------------"+\
{}
"""

IPYRAD_EXCEPTION = """\

"""


MISSING_PAIRFILE_ERROR = """\
    Paired file names must be identical except for _R1_ and _R2_. 
    Example, there are not matching files for samples: \n{}
    """

PAIRED_FILENAMES_ERROR = """\
    Fastq filenames are malformed. R1 must contain the string _R1_ and
    R2 must be identical to R1, excepting the replacement of _R2_ for _R1_.
    """

REF_NOT_FOUND = """\
    "Warning: reference sequence file not found. This must be an absolute path
    (/home/wat/ipyrad/data/reference.gz) or relative to the directory where
    you're running ipyrad (./data/reference.gz). You entered:
    {}
    """


SORTED_NOT_FOUND = """\
    Error: fastq sequence files in sorted_fastq_path could not be found.
    Please check that the location was entered correctly and that a wild
    card selector (*) was used to select all or a subset of files.
    You entered: {}
    """

SORTED_ISDIR = """\
    Error: You entered the path to a directory for sorted_fastq_path. To
    ensure the correct files in the directory are selected, please use a
    wildcard selector to designate the desired files.
    Example: /home/user/data/*.fastq   ## selects all files ending in '.fastq'
    You entered: {}
    """




CANNOT_CHANGE_ASSEMBLY_NAME = """\
    Warning: Assembly name is set at Assembly creation time and is an immutable
    property: You may, however, branch the assembly which will create a copy
    with a new name, but retain a copy of the original Assembly. Here's how:

    Command Line Interface:
        ipyrad -p params-old-name.txt -b new-name

    API (Jupyter Notebook Users):
        new_assembly = my_assembly.branch("new_name")
    """

REQUIRE_ASSEMBLY_NAME = """\
    Assembly name _must_ be set. This is the first parameter in the params.txt
    file, and will be used as a prefix for output files. It should be a short
    string with no special characters, i.e., not a path (no \"/\" characters).
    If you need a suggestion, name it after the organism you're working on.
    """
REQUIRE_REFERENCE_PATH = """\
    Assembly method '{}' requires a 'reference_sequence' in parameter settings.
    """
BAD_ASSEMBLY_NAME = """\
    No spaces or special characters of any kind are allowed in the assembly 
    name. Special characters include all punctuation except dash '-' and 
    underscore '_'. A good practice is to replace spaces with underscores '_'.
    An example of a good assembly_name is: white_crowned_sparrows 
    
    Here's what you put:
    {}
    """
BAD_BARCODE = """\
    One or more barcodes contain invalid IUPAC nucleotide code characters.
    Barcodes must contain only characters from this list "RKSYWMCATG".
    Doublecheck your barcodes file is properly formatted.
    """
NO_MIN_SAMPLES_PER_POP = """\n\
    Population assignment file must include a line indicating the minimum
    number of samples per population. This line should come at the end
    of the file and should be preceded by a hash sign (#), e.g.:

    # pop1:3 pop2:3 pop3:3
    """

BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """
PARAMS_EXISTS = """
    Error: Params file already exists: {}
    Use force argument to overwrite.
    """
EDITS_EXIST = """\
    Skipping: All {} selected Samples already edited.
    (can overwrite with force argument)\
    """
CLUSTERS_EXIST = """\
    Skipping: All {} selected Samples already clustered.
    (can overwrite with force argument)\
    """
JOINTS_EXIST = """\
    Skipping: All {} selected Samples already joint estimated
    (can overwrite with force argument)\
    """
CONSENS_EXIST = """\
    Skipping: All {} selected Samples already consensus called
    (can overwrite with force argument)\
    """
DATABASE_EXISTS = """\
    Skipping: All {} selected Samples already clustered.
    (can overwrite with force argument)\
    """
NOT_CLUSTERED_YET = """\
    The following Samples do not appear to have been clustered in step6
    (i.e., they are not in {}).
    Check for typos in Sample names, or try running step6 including the
    selected samples.

    Missing: {}
    """
OUTPUT_EXISTS = """\
    Output files already created for this Assembly in:
    {}
    To overwrite, rerun using the force argument. 
    """
FIRST_RUN_1 = """\
    No Samples found. First run step 1 to load raw or demultiplexed fastq
    files from the raw_fastq_path or sorted_fastq_path, respectively.
    """
FIRST_RUN_2 = """\
    No Samples ready to be clustered. First run step 2.
    """
FIRST_RUN_3 = """\
    No Samples ready for estimation. First run step 3.
    """
FIRST_RUN_4 = """\
    No Samples ready for consensus calling. First run step 4.
    """
FIRST_RUN_5 = """\
    No Samples ready for clustering. First run step 5.
"""
FIRST_RUN_6 = """\
    Database file {} not found. First run step 6.
"""

BRANCH_NAMES_AND_INPUT = \
"""
    Attempting to branch passing in subsample names
    and an input file, ignoring 'subsamples' argument
"""

########################################################



if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.get_params()
    print(DATA.log)
