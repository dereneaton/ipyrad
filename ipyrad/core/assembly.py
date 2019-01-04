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
import ipyparallel as ipp

from collections import OrderedDict
from ipyrad.assemble.utils import IPyradParamsError, IPyradError
from ipyrad.assemble.utils import ObjDict, IPyradWarningExit
from ipyrad.core.paramsinfo import paraminfo, paramname


# GLOBALS
OUTPUT_FORMATS = {
    'l': 'loci',
    'p': 'phy',
    's': 'snps',
    'n': 'nex',
    'k': 'struct',
    'a': 'alleles',
    'g': 'geno',
    'G': "gphocs",
    'u': 'usnps',
    'v': 'vcf',
    't': 'treemix',
    'm': 'migrate-n',
    }
    

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
        self._ipcluster = {
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
            if key in self._ipcluster:
                self._ipcluster[key] = val

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


    def branch(self, newname, subsamples=None, infile=None):
        """
        Returns a copy of the Assembly object. Does not allow Assembly
        object names to be replicated in namespace or path.
        """
        ## subsample by removal or keeping.
        remove = 0

        ## is there a better way to ask if it already exists?
        if (newname == self.name or os.path.exists(
                os.path.join(
                    self.params.project_dir,
                    newname + ".assembly"))):
            self._print(
                "Assembly object named {} already exists".format(newname))

        else:
            # Make sure the new name doesn't have any wacky characters
            check_name(newname)

            ## Bozo-check. Carve off 'params-' if it's in the new name.
            if newname.startswith("params-"):
                newname = newname.split("params-")[1]

            ## create a copy of the Assembly obj
            newobj = copy.deepcopy(self)
            newobj.name = newname
            newobj.params._assembly_name = newname

            if subsamples and infile:
                print(BRANCH_NAMES_AND_INPUT)

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


    def _get_parallel(self, ipyclient, show_cluster):
        "Start a parallel client with _ipcluster dict or start a new one."

        # connect to a running client or raise an error if not found.
        if not ipyclient:
            ipyclient = ip.core.parallel.get_client(self)

        # print a message about the cluster status
        if self._cli or show_cluster:
            ip.cluster_info(ipyclient=ipyclient, spacer=self._spacer)

        # print an Assembly name header if inside API
        if not self._cli:
            self._print("Assembly: {}".format(self.name))

        # store ipyclient engine pids to the Assembly so we can hard-interrupt
        # later if assembly is interrupted. Only stores pids of engines that 
        # aren't busy at this moment, otherwise it would block.
        self._ipcluster["pids"] = {}
        for eid in ipyclient.ids:
            engine = ipyclient[eid]
            if not engine.outstanding:
                pid = engine.apply(os.getpid).get()
                self._ipcluster["pids"][eid] = pid
        return ipyclient


    def _cleanup_parallel(self, ipyclient):
        try:
            # save the Assembly... should we save here or not. I added save
            # to substeps in six so we don't need it here. By removing it
            # here there is less chance someone will bonk their json file
            # by executing a bad run() command.
            #self.save()  

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
                except ipp.NoEnginesRegistered:
                    pass

                # if CLI, stop jobs and shutdown. Don't use _cli here 
                # because you can have a CLI object but use the --ipcluster
                # flag, in which case we don't want to kill ipcluster.
                if self._cli:
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                else:
                    if not ipyclient.outstanding:
                        ipyclient.purge_everything()
                    else:
                        # nanny: kill everything, something bad happened
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        print("\nwarning: ipcluster shutdown and must be restarted")
                
        # if exception is close and save, print and ignore
        except Exception as inst2:
            print("warning: error during shutdown:\n{}".format(inst2))


    def run(self, steps=None, force=False, ipyclient=None, quiet=False, show_cluster=False):
        """
        Run assembly steps of an ipyrad analysis. Enter steps as a string,
        e.g., "1", "123", "12345". This step checks for an existing
        ipcluster instance otherwise it raises an exception. The ipyparallel
        connection is made using information from the _ipcluster dict of the
        Assembly class object.
        """
        # hide all messages/progress bars       
        self.quiet = quiet

        # check that mindepth params are compatible, fix and report warning.
        self._compatible_params_check()

        # wrap everything in a try statement to ensure that we save the
        # Assembly object if it is interrupted at any point, and also
        # to ensure proper cleanup of the ipyclient.
        try:
            # get a running ipcluster instance or start one 
            ipyclient = self._get_parallel(ipyclient, show_cluster)

            # get the list of steps to run
            if steps:
                if isinstance(steps, int):
                    steps = str(steps)
                steps = sorted(list(steps))
            else:
                print("No assembly steps selected (e.g., '123')")

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

        except KeyboardInterrupt as inst:
            print("\n{}Keyboard Interrupt by user\n".format(self._spacer))

        except (IPyradWarningExit, IPyradError) as inst:
            print("\n{}Encountered an IPyradError:\n{}{}\n".format(
                self._spacer, self._spacer, inst))
            raise

        except Exception as inst:
            print("\n{}Encountered an unexpected error:\n{}{}\n".format(
                self._spacer, self._spacer, inst))
            raise

        # close client when done or interrupted
        finally:
            # unquiet 
            self.quiet = False
            self._cleanup_parallel(ipyclient)



class Hackers:
    def __init__(self):

        # private dictionary so we can check dtypes before changing 
        self._data = dict([
            ("random_seed", 42),
            ("max_fragment_length", 50),
            ("max_inner_mate_distance", 60),
            ("p5_adapter", "AGATCGGAAGAGC"),
            ("p3_adapter", "AGATCGGAAGAGC"),
            ("p3_adapters_extra", []),
            ("p5_adapters_extra", []),
            ("query_cov", None),
            ("bwa_args", ""),
            ("demultiplex_on_i7_tags", False),
            ("declone_PCR_duplicates", False),
            ("merge_technical_replicates", False),
            # ("output_loci_name_buffer", 5),
        ])

    # pretty printing of object
    def __repr__(self):
        printstr = ""
        for idx, (key, val) in enumerate(self._data.items()):
            printstr += "{:<4}{:<28}{:<45}\n".format(idx, key, str(val))
        return printstr
    
    def __str__(self):
        return self.__repr__()        

    # setters
    @property
    def random_seed(self):
        return self._data["random_seed"]
    @random_seed.setter
    def random_seed(self, value):
        self._data["random_seed"] = int(value)

    @property
    def max_fragment_length(self):
        return self._data["max_fragment_length"]
    @max_fragment_length.setter
    def max_fragment_length(self, value):
        self._data["max_fragment_length"] = int(value)

    @property
    def max_inner_mate_distance(self):
        return self._data["max_inner_mate_distance"]
    @max_inner_mate_distance.setter
    def max_inner_mate_distance(self, value):
        self._data["max_inner_mate_distance"] = int(value)

    @property
    def p5_adapter(self):
        return self._data["p5_adapter"]
    @p5_adapter.setter
    def p5_adapter(self, value):
        self._data["p5_adapter"] = str(value)

    @property
    def p5_adapters_extra(self):
        return self._data["p5_adapters_extra"]
    @p5_adapters_extra.setter
    def p5_adapters_extra(self, value):
        if isinstance(value, str):
            self._data["p5_adapters_extra"] = [value]
        else:
            self._data["p5_adapters_extra"] = value

    @property
    def p3_adapter(self):
        return self._data["p3_adapter"]
    @p3_adapter.setter
    def p3_adapter(self, value):
        self._data["p3_adapter"] = str(value)

    @property
    def p3_adapters_extra(self):
        return self._data["p3_adapters_extra"]
    @p3_adapters_extra.setter
    def p3_adapters_extra(self, value):
        if isinstance(value, str):
            self._data["p3_adapters_extra"] = [value]
        else:
            self._data["p3_adapters_extra"] = value

    @property
    def query_cov(self):
        return self._data["query_cov"]
    @query_cov.setter
    def query_cov(self, value):
        self._data["query_cov"] = float(value)

    @property
    def bwa_args(self):
        return self._bwa_args
    @bwa_args.setter
    def bwa_args(self, value):
        self._data["bwa_args"] = str(value)

    @property
    def demultiplex_on_i7_tags(self):
        return self._data["demultiplex_on_i7_tags"]
    @demultiplex_on_i7_tags.setter
    def demultiplex_on_i7_tags(self, value):
        self._data["demultiplex_on_i7_tags"] = bool(value)

    @property
    def declone_PCR_duplicates(self):
        return self._data["declone_PCR_duplicates"]
    @declone_PCR_duplicates.setter
    def declone_PCR_duplicates(self, value):
        self._data["declone_PCR_duplicates"] = bool(value)

    @property
    def merge_technical_replicates(self):
        return self._data["merge_technical_replicates"]
    @merge_technical_replicates.setter
    def merge_technical_replicates(self, value):
        self._data["merge_technical_replicates"] = bool(value)

    

class Params:
    def __init__(self, data):

        # harder to 'update' values if data is here...
        self._data = data

        self._assembly_name = data.name
        self._project_dir = os.path.realpath("./analysis-ipyrad")
        self._raw_fastq_path = ""
        self._barcodes_path = ""
        self._sorted_fastq_path = ""
        self._assembly_method = "denovo"
        self._reference_sequence = ""
        self._datatype = "rad"
        self._restriction_overhang = ("TGCAG", "")
        self._max_low_qual_bases = 5
        self._phred_Qscore_offset = 33
        self._mindepth_statistical = 6
        self._mindepth_majrule = 6
        self._maxdepth = 10000
        self._clust_threshold = 0.85
        self._max_barcode_mismatch = 0
        self._filter_adapters = 0
        self._filter_min_trim_len = 35
        self._max_alleles_consens = 2
        self._max_Ns_consens = (5, 5)
        self._max_Hs_consens = (8, 8)
        self._min_samples_locus = 4
        self._max_SNPs_locus = (20, 20)
        self._max_Indels_locus = (8, 8)
        self._max_shared_Hs_locus = 0.5
        self._trim_reads = (0, 0, 0, 0)
        self._trim_loci = (0, 0, 0, 0)
        self._output_formats = list("psv")
        self._pop_assign_file = ""
        
        self._keys = [
            "_assembly_name",
            "_project_dir",
            "_raw_fastq_path",
            "_barcodes_path",
            "_sorted_fastq_path", 
            "_assembly_method",
            "_reference_sequence",
            "_datatype", 
            "_restriction_overhang",
            "_max_low_qual_bases", 
            "_phred_Qscore_offset", 
            "_mindepth_statistical", 
            "_mindepth_majrule", 
            "_maxdepth", 
            "_clust_threshold", 
            "_max_barcode_mismatch", 
            "_filter_adapters", 
            "_filter_min_trim_len",
            "_max_alleles_consens", 
            "_max_Ns_consens", 
            "_max_Hs_consens", 
            "_min_samples_locus", 
            "_max_SNPs_locus", 
            "_max_Indels_locus", 
            "_max_shared_Hs_locus", 
            "_trim_reads", 
            "_trim_loci", 
            "_output_formats", 
            "_pop_assign_file",            
        ]
                
        
    def __repr__(self):
        fullcurdir = os.path.realpath(os.path.curdir)
        printstr = ""
        for idx, key in enumerate(self._keys):
            value = self.__getattribute__(key)
            if isinstance(value, str):
                value = value.replace(fullcurdir + "/", "./")
                value = value.replace(os.path.expanduser("~"), "~")
            printstr += "{:<4}{:<28}{:<45}\n".format(idx, key[1:], str(value))
        return printstr
    
    def __str__(self):
        return self.__repr__()
        
    @property
    def assembly_name(self):
        return self._assembly_name
    @assembly_name.setter    
    def assembly_name(self, value):
        raise IPyradError(CANNOT_CHANGE_ASSEMBLY_NAME)


    @property
    def project_dir(self):
        return self._project_dir
    @project_dir.setter
    def project_dir(self, value):
        if " " in value:
            raise IPyradError(BAD_PROJDIR_NAME.format(value))
        self._project_dir = os.path.realpath(os.path.expanduser(value))
        self._data.dirs.project = self._project_dir


    @property
    def raw_fastq_path(self):
        return self._raw_fastq_path
    @raw_fastq_path.setter
    def raw_fastq_path(self, value):
        if value and ("Merged:" not in value):
            fullpath = os.path.realpath(os.path.expanduser(value))
            if os.path.isdir(fullpath):
                raise IPyradError(RAW_PATH_ISDIR.format(fullpath))
            elif glob.glob(fullpath):
                self._raw_fastq_path = fullpath
            else:
                raise IPyradError(NO_RAW_FILE.format(fullpath))
        # if 'Merged:' in value then set to ""
        else:
            self._raw_fastq_path = ""


    @property
    def barcodes_path(self):
        return self._barcodes_path
    @barcodes_path.setter
    def barcodes_path(self, value):
        if value and ("Merged:" not in value):

            # allow fuzzy name match
            fullbar = glob.glob(os.path.realpath(os.path.expanduser(value)))
            if not fullbar:
                raise IPyradError(BARCODE_NOT_FOUND.format(fullbar))

            # file must exist
            fullbar = fullbar[0]
            if not os.path.exists(fullbar):
                raise IPyradError(BARCODE_NOT_FOUND.format(fullbar))

            else:
                self._barcodes_path = fullbar
                self._data._link_barcodes()
        # if 'Merged:' in value then set to ""
        else:
            self._barcodes_path = ""


    @property
    def sorted_fastq_path(self):
        return self._sorted_fastq_path
    @sorted_fastq_path.setter
    def sorted_fastq_path(self, value):
        if value and ("Merged:" not in value):
            fullpath = os.path.realpath(os.path.expanduser(value))
            if os.path.isdir(fullpath):
                raise IPyradError(SORTED_ISDIR.format(fullpath))
            elif glob.glob(fullpath):
                self._sorted_fastq_path = fullpath
            else:
                raise IPyradError(SORTED_NOT_FOUND.format(fullpath))
        # if 'Merged:' in value then set to ""
        else:
            self._sorted_fastq_path = ""


    @property
    def assembly_method(self):
        return self._assembly_method
    @assembly_method.setter
    def assembly_method(self, value):
        allowed = ["denovo", "reference", "denovo+reference", "denovo-reference"]
        assert value in allowed, BAD_ASSEMBLY_METHOD.format(value)
        self._assembly_method = value


    @property
    def reference_sequence(self):
        return self._reference_sequence
    @reference_sequence.setter
    def reference_sequence(self, value):
        if value:
            fullpath = os.path.realpath(os.path.expanduser(value))
            if not os.path.exists(fullpath):
                raise IPyradError("reference sequence file not found")
            if fullpath.endswith(".gz"):
                raise IPyradError("reference sequence file must be decompressed.")
            self._reference_sequence = fullpath
        else:
            self._reference_sequence = ""


    @property
    def datatype(self):
        return self._datatype
    @datatype.setter
    def datatype(self, value):
        allowed = (
            'rad', 'gbs', 'ddrad', 'pairddrad', 
            'pairgbs', '2brad', 'pair3rad', 'merged'
        )
        assert value in allowed, (
            "datatype must be one of: {}".format(", ".join(allowed)))
        self._datatype = str(value)

        # update barcode dist for expected second bcode
        if "3rad" in value:
            if self.params.barcodes_path:
                self._data._link_barcodes()

            
    @property 
    def restriction_overhang(self):
        return self._restriction_overhang
    @restriction_overhang.setter
    def restriction_overhang(self, value):
        # returns string values as a tuple ("", "") or ("",)
        value = tuplecheck(value, str)
        
        # expand GBS for user if they set only one cutter 
        if (self.datatype == "GBS") & (len(value) == 1):
            value = (value[0], value[0])

        # Handle 3rad datatype with only 3 cutters
        elif len(value) == 3:
            value = (value[0], value[1], value[2], "")
            if self.params.barcodes_path:
                self._data._link_barcodes()

        assert len(value) <= 4, """
    most datasets require 1 or 2 cut sites, e.g., (TGCAG, '') or (TGCAG, CCGG).
    For 3rad/seqcap may be up to 4 cut sites."""
        self._restriction_overhang = value


    @property
    def max_low_qual_bases(self):
        return self._max_low_qual_bases
    @max_low_qual_bases.setter
    def max_low_qual_bases(self, value):
        try:
            value = int(value)
        except TypeError:
            raise IPyradParamsError("max_low_qual_bases must be an integer.")
        self._max_low_qual_bases = value


    @property
    def phred_Qscore_offset(self):
        return self._phred_Qscore_offset
    @phred_Qscore_offset.setter
    def phred_Qscore_offset(self, value):
        try:
            value = int(value)
        except TypeError:
            raise IPyradParamsError("phred_Qscore_offset must be an integer.")
        self._phred_Qscore_offset = int(value)


    @property
    def mindepth_statistical(self):
        return self._mindepth_statistical
    @mindepth_statistical.setter
    def mindepth_statistical(self, value):
        try:
            value = int(value)
        except TypeError:
            raise IPyradParamsError("mindepth_statistical must be an integer.")
        # do not allow values below 5
        assert int(value) >= 5, (
            "mindepth_statistical cannot be <5. Set mindepth_majrule instead.")
        self._mindepth_statistical = int(value)


    @property
    def mindepth_majrule(self):
        return self._mindepth_majrule
    @mindepth_majrule.setter
    def mindepth_majrule(self, value):
        try:
            value = int(value)
        except TypeError:
            raise IPyradParamsError("mindepth_majrule must be an integer.")
        self._mindepth_majrule = int(value)


    @property
    def maxdepth(self):
        return self._maxdepth
    @maxdepth.setter
    def maxdepth(self, value):
        self._maxdepth = int(value)


    @property
    def clust_threshold(self):
        return self._clust_threshold
    @clust_threshold.setter
    def clust_threshold(self, value):
        value = float(value)
        assert (value < 1) & (value > 0), (
            "clust_threshold must be a decimal value between 0 and 1.")
        self._clust_threshold = value


    @property
    def max_barcode_mismatch(self):
        return self._max_barcode_mismatch
    @max_barcode_mismatch.setter
    def max_barcode_mismatch(self, value):
        self._max_barcode_mismatch = int(value)


    @property
    def filter_adapters(self):
        return self._filter_adapters
    @filter_adapters.setter
    def filter_adapters(self, value):
        value = int(value)
        assert value in (0, 1, 2, 3), "filter_adapters must be 0, 1, 2, or 3"
        self._filter_adapters = value


    @property
    def filter_min_trim_len(self):
        return self._filter_min_trim_len
    @filter_min_trim_len.setter
    def filter_min_trim_len(self, value):
        self._filter_min_trim_len = int(value)


    @property
    def max_alleles_consens(self):
        return self._max_alleles_consens
    @max_alleles_consens.setter
    def max_alleles_consens(self, value):
        self._max_alleles_consens = int(value)


    @property
    def max_Ns_consens(self):
        return self._max_Ns_consens
    @max_Ns_consens.setter
    def max_Ns_consens(self, value):
        value = tuplecheck(value, int)
        assert isinstance(value, tuple), (
            "max_Ns_consens should be a tuple e.g., (5, 5)")
        self._max_Ns_consens = value


    @property
    def max_Hs_consens(self):
        return self._max_Hs_consens
    @max_Hs_consens.setter
    def max_Hs_consens(self, value):
        value = tuplecheck(value, int)
        assert isinstance(value, tuple), (
            "max_Hs_consens should be a tuple e.g., (5, 5)")
        self._max_Hs_consens = value


    @property
    def min_samples_locus(self):
        return self._min_samples_locus
    @min_samples_locus.setter
    def min_samples_locus(self, value):
        self._min_samples_locus = int(value)


    @property
    def max_shared_Hs_locus(self):
        return self._max_shared_Hs_locus
    @max_shared_Hs_locus.setter
    def max_shared_Hs_locus(self, value):
        if isinstance(value, str):
            if value.isdigit():
                value = int(value)
            else:
                try:
                    value = float(value)
                except Exception as inst:
                    raise IPyradParamsError("""
    max_shared_Hs_locus must be int or float, you put: {}""".format(alue))
        self._max_shared_Hs_locus = value


    @property
    def max_SNPs_locus(self):
        return self._max_SNPs_locus
    @max_SNPs_locus.setter
    def max_SNPs_locus(self, value):
        value = tuplecheck(value, int)
        assert isinstance(value, tuple), (
            "max_SNPs_locus should be a tuple e.g., (20, 20)")
        self._max_SNPs_locus = value


    @property
    def max_Indels_locus(self):
        return self._max_Indels_locus
    @max_Indels_locus.setter
    def max_Indels_locus(self, value):
        value = tuplecheck(value, int)
        assert isinstance(value, tuple), (
            "max_Indels_locus should be a tuple e.g., (5, 100)")
        self._max_Indels_locus = value


    @property
    def trim_reads(self):
        return self._trim_reads
    @trim_reads.setter
    def trim_reads(self, value):
        # cast to ints
        value = tuplecheck(value, int)

        # check that entries make sense 
        if value[1] > 0:
            if not value[1] > value[0]:
                raise IPyradError(BAD_TRIM_READS)
        if value[3] > 0:
            if not value[3] > value[2]:
                raise IPyradError(BAD_TRIM_READS)
        if (value[0] < 0) or (value[2] < 0):
            raise IPyradError(BAD_TRIM_READS)       
        self._trim_reads = value


    @property
    def trim_loci(self):
        return self._trim_loci
    @trim_loci.setter
    def trim_loci(self, value):
        value = tuplecheck(value, str)
        assert isinstance(value, tuple), (
            "trim_loci should be a tuple e.g., (0, -5, -5, 0)")
        self._trim_loci = tuple([int(i) for i in value])


    @property
    def output_formats(self):
        return self._output_formats
    @output_formats.setter
    def output_formats(self, value):
        # Handle the case where output formats is an empty string
        if isinstance(value, str):
            # strip commas and spaces from string so we have only letters
            value = value.replace(",", "").replace(" ", "")
            value = list(value)
            if not value:
                value = ["*"]
        if isinstance(value, tuple):
            value = list(value)

        if isinstance(value, list):
            # if more than letters, raise an warning
            if any([len(i) > 1 for i in value]):
                self._print("""
    'output_formats' params entry is malformed. Setting to * to avoid errors.""")
                value = "*"
        
        if "*" in value:
            value = list(OUTPUT_FORMATS.keys())

        # set the param
        self._output_formats = value


    @property
    def pop_assign_file(self):
        return self._pop_assign_file
    @pop_assign_file.setter
    def pop_assign_file(self, value):
        fullpath = os.path.realpath(os.path.expanduser(value))

        # if a path is entered, raise exception if not found
        if value:
            if not os.path.isfile(fullpath):
                raise IPyradError("""
    Warning: Population assignment file not found. This must be an
    absolute path (/home/wat/ipyrad/data/my_popfile.txt) or relative to
    the directory where you're running ipyrad (./data/my_popfile.txt)
    You entered: {}\n""".format(fullpath))
            self._pop_assign_file = fullpath
            self._link_populations()

        else:
            # Don't forget to possibly blank the populations dictionary
            self._pop_assign_file = ""
            self._data.populations = {}


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


def merge(name, assemblies):
    """
    Creates and returns a new Assembly object in which samples from two or more
    Assembly objects with matching names are 'merged'. Merging does not affect 
    the actual files written on disk, but rather creates new Samples that are 
    linked to multiple data files, and with stats summed.
    """
    # make sure assemblies is a list of more than one
    assert len(assemblies) > 1, (
        "You must enter a list of >1 Assemblies to merge")
    assemblies = list(assemblies)

    # create new Assembly as a branch (deepcopy) of the first in list
    merged = assemblies[0].branch(name)

    # get all sample names from all Assemblies
    allsamples = set(merged.samples.keys())
    for iterass in assemblies[1:]:
        allsamples.update(set(iterass.samples.keys()))

    # Make sure we have the max of all values for max frag length
    # from all merging assemblies.
    merged.hackersonly.max_fragment_length = max(
        [i.hackersonly.max_fragment_length for i in assemblies])

    # warning message?
    warning = 0

    # iterate over assembly objects, skip first already copied
    for iterass in assemblies[1:]:
        # iterate over allsamples, add if not in merged
        for sample in iterass.samples:
            # iterate over stats, skip 'state'
            if sample not in merged.samples:
                merged.samples[sample] = copy.deepcopy(iterass.samples[sample])
                # if barcodes data present then keep it
                if iterass.barcodes.get(sample):
                    merged.barcodes[sample] = iterass.barcodes[sample]
            else:
                # merge stats and files of the sample
                for stat in merged.stats.keys()[1:]:
                    merged.samples[sample].stats[stat] += (
                        iterass.samples[sample].stats[stat])
                # merge file references into a list
                for filetype in ['fastqs', 'edits']:
                    merged.samples[sample].files[filetype] += (
                        iterass.samples[sample].files[filetype])
                if iterass.samples[sample].files["clusters"]:
                    warning = 1

    # print warning if clusters or later was present in merged assembly
    if warning:
        print("""\
    Warning: the merged Assemblies contained Samples that are identically named,
    and so ipyrad has attempted to merge these Samples. This is perfectly fine to
    do up until step 3, but not after, because at step 3 all reads for a Sample
    should be included during clustering/mapping. Take note, you can merge Assemblies
    at any step *if they do not contain the same Samples*, however, here that is not
    the case. If you wish to proceed with this merged Assembly you will have to
    start from step 3, therefore the 'state' of the Samples in this new merged
    Assembly ({}) have been set to 2.
    """.format(name))
        for sample in merged.samples:
            merged.samples[sample].stats.state = 2
            ## clear stats
            for stat in ["refseq_mapped_reads", "refseq_unmapped_reads",
                         "clusters_total", "clusters_hidepth", "hetero_est",
                         "error_est", "reads_consens"]:
                merged.samples[sample].stats[stat] = 0
            ## clear files
            for ftype in ["mapped_reads", "unmapped_reads", "clusters",
                          "consens", "database"]:
                merged.samples[sample].files[ftype] = []

    # Set the values for some params that don't make sense inside
    # merged assemblies
    merged_names = ", ".join([i.name for i in assemblies])
    merged.params.raw_fastq_path = "Merged: " + merged_names
    merged.params.barcodes_path = "Merged: " + merged_names
    merged.params.sorted_fastq_path = "Merged: " + merged_names

    # return the new Assembly object
    merged.save()
    return merged


def tuplecheck(newvalue, dtype=str):
    """
    Takes a string argument and returns value as a tuple.
    Needed for paramfile conversion from CLI to set_params args
    """
    if isinstance(newvalue, list):
        newvalue = tuple(newvalue)

    if isinstance(newvalue, str):
        newvalue = newvalue.rstrip(")").strip("(")
        try:
            newvalue = tuple([dtype(i.strip()) for i in newvalue.split(",")])

        ## Type error is thrown by tuple if it's applied to a non-iterable.
        except TypeError:
            newvalue = tuple(dtype(newvalue))

        ## If dtype fails to cast any element of newvalue
        except ValueError:
            raise IPyradError(
                "Assembly.tuplecheck() failed to cast to {} - {}"
                .format(dtype, newvalue)
            )

        except Exception as inst:
            raise IPyradError(
                "\nError: Param`{}` is not formatted correctly.\n({})\n"
                .format(newvalue, inst)
            )
    return newvalue


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

BAD_ASSEMBLY_METHOD = """\
    The assembly_method parameter must be one of the following: denovo, reference,
    denovo+reference or denovo-reference. You entered:
    {}.
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

BARCODE_NOT_FOUND = """\
    Error: barcodes file not found. This must be an absolute path
    (/home/wat/ipyrad/data/data_barcodes.txt) or relative to the directory
    where you're running ipyrad (./data/data_barcodes.txt). You entered:
    {}
    """

RAW_PATH_ISDIR = """\
    Error: You entered the path to a directory for raw_fastq_path. To
    ensure the correct files in the directory are selected, please use a
    wildcard selector to designate the desired files.
    Example: /home/user/data/*.fastq  ## selects all files ending in '.fastq'
    You entered: {}
    """


BAD_PROJDIR_NAME = """\
    Error: Your project_dir contains a directory with a space in the name.
    This can cause all kinds of funny problems so please rename this
    directory and remove the space. Try replacing the space with an underscore.
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

NO_RAW_FILE = """\
    The value entered for the path to the raw fastq file is unrecognized.
    Please be sure this path is correct. Double check the file name and
    the file extension. If it is a relative path be sure the path is
    correct with respect to the directory you're running ipyrad from.
    You entered: {}
"""

BAD_TRIM_READS = """\
    Bad trim_reads entry. Think of these values like slice indices, but with 
    0 as a special character meaning no effect. So (0, 80, 0, 0) trims the 
    first read to 80 bp. (5, 80, 0, 0) trims R1 to keep only bp 5-80. 
    See documentation for details.
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
