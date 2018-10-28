#!/usr/bin/env python

""" ipyrad Assembly class object. """

from __future__ import print_function
import os
import glob
import sys
import copy
import time
import string
import datetime
import itertools
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
        `[Assembly].paramsdict['barcodes_path']`.
    dirs : dict
        Returns a dictionary with the location of directories that contain
        linked Sample object files and stats results.


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
    def __init__(self, name, quiet=False, **kwargs):

        #ip.logger.debug("new assembly: {}".format(name))
        self.quiet = quiet

        # record whether we're in the CLI or API
        self._cli = False
        if kwargs.get("cli"):
            self._cli = True           

        # Make sure assembly name is not empty
        if not name:
            raise IPyradParamsError(REQUIRE_ASSEMBLY_NAME)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        ## This will raise an error immediately if there are bad chars in name.
        self._check_name(name)

        ## print name, can't use run() spacer b/c no ipcluster settings yet.
        self.name = name
        self._print("{}New Assembly: {}".format(self._spacer, self.name))

        ## Store assembly version #
        self._version = ip.__version__

        ## stores default ipcluster launch info
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
        for key, val in kwargs.items():
            if key in self._ipcluster:
                self._ipcluster[key] = val

        ## print headers, this is used as a 'quiet' option
        ## or to control differences in printing between API and CLI
        self._headers = 0

        ## used to set checkpoints within step 6
        self._checkpoint = 0

        ## statsfiles is a dict with file locations
        ## stats_dfs is a dict with pandas dataframes
        self.stats_files = ObjDict({})
        self.stats_dfs = ObjDict({})

        ## samples linked
        self.samples = {}

        ## samples linked
        self.populations = OrderedDict()

        ## multiplex files linked
        self.barcodes = ObjDict()

        ## outfiles locations
        self.outfiles = ObjDict()
        self.outfiles.loci = ""

        ## storing supercatg file
        self.clust_database = ""
        self.database = ""

        ## the default params dict
        self.paramsdict = OrderedDict([
            ("assembly_name", name),
            ("project_dir", "./"),  # os.path.realpath(os.path.curdir)),
            ("raw_fastq_path", ""),
            ("barcodes_path", ""),
            ("sorted_fastq_path", ""),
            ("assembly_method", "denovo"),
            ("reference_sequence", ""),
            ("datatype", 'rad'),
            ("restriction_overhang", ("TGCAG", "")),
            ("max_low_qual_bases", 5),
            ("phred_Qscore_offset", 33),
            ("mindepth_statistical", 6),
            ("mindepth_majrule", 6),
            ("maxdepth", 10000),
            ("clust_threshold", 0.85),
            ("max_barcode_mismatch", 0),
            ("filter_adapters", 0),
            ("filter_min_trim_len", 35),
            ("max_alleles_consens", 2),
            ("max_Ns_consens", (5, 5)),
            ("max_Hs_consens", (8, 8)),
            ("min_samples_locus", 4),
            ("max_SNPs_locus", (20, 20)),
            ("max_Indels_locus", (8, 8)),
            ("max_shared_Hs_locus", 0.50),
            ("trim_reads", (0, 0, 0, 0)),
            ("trim_loci", (0, 0, 0, 0)),
            ("output_formats", ['p', 's', 'v']),
            ("pop_assign_file", ""),
        ])

        ## Store data directories for this Assembly. Init with default project
        self.dirs = ObjDict({
            "project": os.path.realpath(self.paramsdict["project_dir"]),
            "fastqs": "",
            "edits": "",
            "clusts": "",
            "consens": "",
            "across": "",
            "outfiles": "",
        })

        ## Default hackers only parameters dictionary
        ## None of the safeguards of the paramsdict, no nice accessor
        ## functions, so you have to be sure of what you're doing if
        ## you change these values.
        self._hackersonly = OrderedDict([
            ("random_seed", 42),
            ("max_fragment_length", 50),
            ("max_inner_mate_distance", 60),
            ("p5_adapter", "AGATCGGAAGAGC"),
            ("p3_adapter", "AGATCGGAAGAGC"),
            ("p3_adapters_extra", []),
            ("p5_adapters_extra", []),
            ("preview_step1", 4000000),
            ("preview_step2", 100000),
            ("output_loci_name_buffer", 5),
            ("query_cov", None),
            ("smalt_index_wordlen", 8),
            ("aligner", "bwa"),
            ("min_SE_refmap_overlap", 10),
            ("refmap_merge_PE", True),
            ("bwa_args", "")
        ])

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
        """ Returns a data frame with Sample stats for each step """
        nameordered = list(self.samples.keys())
        nameordered.sort()
        newdat = pd.DataFrame(
            (self.samples[i].stats_dfs[idx] for i in nameordered),
            index=nameordered,
        ).dropna(axis=1, how='all')
        return newdat


    def _check_name(self, name):
        "Test assembly name is valid and raise if any special characters"
        invalid_chars = (
            string.punctuation.replace("_", "").replace("-", "") + " ")
        if any(char in invalid_chars for char in name):
            raise IPyradParamsError(BAD_ASSEMBLY_NAME.format(name))


    def _link_barcodes(self):
        """
        Private function. Links Sample barcodes in a dictionary as
        [Assembly].barcodes, with barcodes parsed from the 'barcodes_path'
        parameter. This function is called during set_params() when setting
        the barcodes_path.
        """

        ## parse barcodefile
        try:
            ## allows fuzzy match to barcodefile name
            barcodefile = glob.glob(self.paramsdict["barcodes_path"])[0]

            ## read in the file
            bdf = pd.read_csv(barcodefile, header=None, delim_whitespace=1)
            bdf = bdf.dropna()

            ## make sure bars are upper case
            bdf[1] = bdf[1].str.upper()

            ## if replicates are present then print a warning
            reps = bdf[0].unique().shape[0] != bdf[0].shape[0]
            if reps:
                print("{spacer}Warning: technical replicates (same name) will be combined."\
                      .format(**{'spacer': self._spacer}))
                ## add -technical-replicate-N to replicate names
                reps = [i for i in bdf[0] if list(bdf[0]).count(i) > 1]
                ureps = list(set(reps))
                for name in ureps:
                    idxs = bdf[bdf[0] == ureps[0]].index.tolist()
                    for num, idx in enumerate(idxs):
                        bdf.ix[idx][0] = bdf.ix[idx][0] + "-technical-replicate-" + str(num+1)

            ## make sure chars are all proper
            if not all(bdf[1].apply(set("RKSYWMCATG").issuperset)):
                ip.logger.warn(BAD_BARCODE)
                raise IPyradError(BAD_BARCODE)

            ## 3rad/seqcap use multiplexed barcodes
            ## We'll concatenate them with a plus and split them later
            if "3rad" in self.paramsdict["datatype"]:
                try:
                    bdf[2] = bdf[2].str.upper()
                    self.barcodes = dict(zip(bdf[0], bdf[1] + "+" + bdf[2]))
                except KeyError as inst:
                    msg = "    3rad assumes multiplexed barcodes. Doublecheck your barcodes file."
                    ip.logger.error(msg)
                    raise IPyradError(msg)
            else:
                ## set attribute on Assembly object
                self.barcodes = dict(zip(bdf[0], bdf[1]))

        except (IOError, IndexError):
            raise IPyradError(
                "    Barcodes file not found. You entered: {}"
                .format(self.paramsdict["barcodes_path"]))

        except ValueError as inst:
            msg = "    Barcodes file format error."
            ip.logger.warn(msg)
            raise IPyradError(inst)


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
            popfile = glob.glob(self.paramsdict["pop_assign_file"])[0]
            if not os.path.exists(popfile):
                raise IPyradError("Population assignment file not found: {}"\
                                  .format(self.paramsdict["pop_assign_file"]))

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
                ip.logger.warn("Populations file may be malformed.")
                raise IPyradError("  Populations file malformed - {}"\
                                  .format(popfile))

        else:
            ## pop dict is provided by user
            if not popmins:
                popmins = {i: 1 for i in popdict}

        ## check popdict. Filter for bad samples
        ## Warn user but don't bail out, could be setting the pops file
        ## on a new assembly w/o any linked samples.
        badsamples = [
            i for i in itertools.chain(*popdict.values())
            if i not in self.samples.keys()]

        if any(badsamples):
            ip.logger.warn(
                "Some names from population input do not match Sample "\
              + "names: ".format(", ".join(badsamples)))
            ip.logger.warn("If this is a new assembly this is normal.")

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
        fullcurdir = os.path.realpath(os.path.curdir)
        if not param:
            for index, (key, value) in enumerate(self.paramsdict.items()):
                if isinstance(value, str):
                    value = value.replace(fullcurdir + "/", "./")
                sys.stdout.write(
                    "{}{:<4}{:<28}{:<45}\n"
                    .format(self._spacer, index, key, str(value)),
                )
        else:
            try:
                if int(param):
                    #sys.stdout.write(self.paramsdict.values()[int(param)-1])
                    return self.paramsdict.values()[int(param)]
            except (ValueError, TypeError, NameError, IndexError):
                try:
                    return self.paramsdict[param]
                except KeyError:
                    return 'key not recognized'


    def set_params(self, param, newvalue):
        """
        Set a parameter to a new value. Raises error if newvalue is wrong type.

        Note
        ----
        Use [Assembly].get_params() to see the parameter values currently
        linked to the Assembly object.

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

        """

        ## this includes current params and some legacy params for conversion
        legacy_params = ["edit_cutsites", "trim_overhang"]
        current_params = self.paramsdict.keys()
        allowed_params = list(current_params) + legacy_params

        ## require parameter recognition
        if param not in allowed_params:
            raise IPyradParamsError(
                "Parameter key not recognized: {}".format(param))

        ## make string
        param = str(param)
        ## get index if param is keyword arg (this index is now zero based!)
        if len(param) < 3:
            param = list(self.paramsdict.keys())[int(param)]

        ## run assertions on new param
        try:
            self = _paramschecker(self, param, newvalue)

        except Exception as inst:
            raise IPyradError(
                BAD_PARAMETER.format(param, inst, newvalue))


    def write_params(self, outfile=None, force=False):
        """ 
        Write out the parameters of this assembly to a file properly
        formatted as input for `ipyrad -p <params.txt>`. A good and
        simple way to share/archive parameter settings for assemblies.
        This is also the function that's used by __main__ to
        generate default params.txt files for `ipyrad -n`
        """
        if outfile is None:
            outfile = "params-" + self.name + ".txt"

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

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key, val in self.paramsdict.items():  
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                padding = (" " * (30 - len(paramvalue)))
                paramkey = list(self.paramsdict.keys()).index(key)
                paramindex = " ## [{}] ".format(paramkey)
                ip.logger.debug("\t".join([str(i) for i in
                    (key, val, paramindex)]))
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
                    self.paramsdict["project_dir"],
                    newname + ".assembly"))):
            print("{}Assembly object named {} already exists"\
                  .format(self._spacer, newname))

        else:
            ## Make sure the new name doesn't have any wacky characters
            self._check_name(newname)

            ## Bozo-check. Carve off 'params-' if it's in the new name.
            if newname.startswith("params-"):
                newname = newname.split("params-")[1]

            ## create a copy of the Assembly obj
            newobj = copy.deepcopy(self)
            newobj.name = newname
            newobj.paramsdict["assembly_name"] = newname

            if subsamples and infile:
                print(BRANCH_NAMES_AND_INPUT)

            if infile:
                if infile[0] == "-":
                    remove = 1
                    infile = infile[1:]
                if os.path.exists(infile):
                    subsamples = _read_sample_names(infile)

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


    def save(self):
        """
        Save Assembly object to disk as a JSON file. Used for checkpointing,
        ipyrad auto-saves after every assembly step. File is saved to:
        [project_dir]/[assembly_name].json
        """
        ip.save_json(self)


    def _compatible_params_check(self):
        "check params that must be compatible at run time"

        # do not allow statistical < majrule
        val1 = self.paramsdict["mindepth_statistical"]
        val2 = self.paramsdict['mindepth_majrule']
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
                            ip.logger.info(
                                'interrupted engine {} w/ SIGINT to {}'
                                .format(engine_id, pid))
                    time.sleep(1)
                except ipp.NoEnginesRegistered:
                    pass

                # if CLI, stop jobs and shutdown. Don't use _cli here 
                # because you can have a CLI object but use the --ipcluster
                # flag, in which case we don't want to kill ipcluster.
                if self._cli:
                    ip.logger.info("  shutting down engines")
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                    ip.logger.info("  finished shutdown")
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
            ip.logger.error("shutdown warning: %s", inst2)


    def run(self, steps=0, force=False, ipyclient=None, show_cluster=0):
        """
        Run assembly steps of an ipyrad analysis. Enter steps as a string,
        e.g., "1", "123", "12345". This step checks for an existing
        ipcluster instance otherwise it raises an exception. The ipyparallel
        connection is made using information from the _ipcluster dict of the
        Assembly class object.
        """
        ## check that mindepth params are compatible, fix and report warning.
        self._compatible_params_check()

        ## wrap everything in a try statement to ensure that we save the
        ## Assembly object if it is interrupted at any point, and also
        ## to ensure proper cleanup of the ipyclient.
        try:
            # get a running ipcluster instance or start one 
            ipyclient = self._get_parallel(ipyclient, show_cluster)

            # get the list of steps to run
            if isinstance(steps, int):
                steps = str(steps)
            steps = sorted(list(steps))

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

            # tell user if they forgot to enter steps
            if not steps:
                print("No assembly steps selected (e.g., '123')")

            # run step fuctions and save and clear memory after each
            for step in steps:
                stepdict[step](self, force, ipyclient).run()
                self.save()
                ipyclient.purge_everything()

        except KeyboardInterrupt as inst:
            print("\n{}Keyboard Interrupt by user".format(self._spacer))

        except (IPyradWarningExit, IPyradError) as inst:
            print("\n{}Encountered an IPyradError:\n{}{}".format(
                self._spacer, self._spacer, inst))
            raise

        except Exception as inst:
            print("\n{}Encountered an unexpected error:\n{}{}".format(
                self._spacer, self._spacer, inst))
            raise

        # close client when done or interrupted
        finally:
            self._cleanup_parallel(ipyclient)



def _read_sample_names(fname):
    """ Read in sample names from a plain text file. This is a convenience
    function for branching so if you have tons of sample names you can
    pass in a file rather than having to set all the names at the command
    line.
    """
    try:
        with open(fname, 'r') as infile:
            subsamples = [x.split()[0] for x in infile.readlines() if x.strip()]

    except Exception as inst:
        print("Failed to read input file with sample names.\n{}".format(inst))
        raise inst

    return subsamples



def _expander(namepath):
    """ expand ./ ~ and ../ designators in location names """
    if "~" in namepath:
        namepath = os.path.expanduser(namepath)
    else:
        namepath = os.path.abspath(namepath)
    return namepath



def merge(name, assemblies):
    """
    Creates and returns a new Assembly object in which samples from two or more
    Assembly objects with matching names are 'merged'. Merging does not affect 
    the actual files written on disk, but rather creates new Samples that are 
    linked to multiple data files, and with stats summed.
    """

    ## checks
    assemblies = list(assemblies)

    ## create new Assembly as a branch (deepcopy)
    merged = assemblies[0].branch(name)

    ## get all sample names from all Assemblies
    allsamples = set(merged.samples.keys())
    for iterass in assemblies[1:]:
        allsamples.update(set(iterass.samples.keys()))

    ## Make sure we have the max of all values for max frag length
    ## from all merging assemblies.
    merged._hackersonly["max_fragment_length"] =\
        max([x._hackersonly["max_fragment_length"] for x in assemblies])

    ## warning message?
    warning = 0

    ## iterate over assembly objects, skip first already copied
    for iterass in assemblies[1:]:
        ## iterate over allsamples, add if not in merged
        for sample in iterass.samples:
            ## iterate over stats, skip 'state'
            if sample not in merged.samples:
                merged.samples[sample] = copy.deepcopy(iterass.samples[sample])
                ## if barcodes data present then keep it
                if iterass.barcodes.get(sample):
                    merged.barcodes[sample] = iterass.barcodes[sample]
            else:
                ## merge stats and files of the sample
                for stat in merged.stats.keys()[1:]:
                    merged.samples[sample].stats[stat] += \
                                iterass.samples[sample].stats[stat]
                ## merge file references into a list
                for filetype in ['fastqs', 'edits']:
                    merged.samples[sample].files[filetype] += \
                                iterass.samples[sample].files[filetype]
                if iterass.samples[sample].files["clusters"]:
                    warning += 1

    ## print warning if clusters or later was present in merged assembly
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

    ## Set the values for some params that don't make sense inside
    ## merged assemblies
    merged_names = ", ".join([x.name for x in assemblies])
    merged.paramsdict["raw_fastq_path"] = "Merged: " + merged_names
    merged.paramsdict["barcodes_path"] = "Merged: " + merged_names
    merged.paramsdict["sorted_fastq_path"] = "Merged: " + merged_names

    ## return the new Assembly object
    merged.save()
    return merged


def _tuplecheck(newvalue, dtype=str):
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
            ip.logger.info("Assembly.tuplecheck() failed to cast to {} - {}"\
                        .format(dtype, newvalue))
            raise

        except Exception as inst:
            ip.logger.info(inst)
            raise SystemExit(\
            "\nError: Param`{}` is not formatted correctly.\n({})\n"\
                 .format(newvalue, inst))

    return newvalue


def _paramschecker(self, param, newvalue):
    """ Raises exceptions when params are set to values they should not be"""

    if param == 'assembly_name':
        ## Make sure somebody doesn't try to change their assembly_name, bad
        ## things would happen. Calling set_params on assembly_name only raises
        ## an informative error. Assembly_name is set at Assembly creation time
        ## and is immutable.
        raise IPyradError(CANNOT_CHANGE_ASSEMBLY_NAME)

    elif param == 'project_dir':
        expandpath = _expander(newvalue)
        if not expandpath.startswith("/"):
            if os.path.exists(expandpath):
                expandpath = _expander(expandpath)
        ## Forbid spaces in path names
        if " " in expandpath:
            raise IPyradError(BAD_PROJDIR_NAME.format(expandpath))
        self.paramsdict["project_dir"] = expandpath
        self.dirs["project"] = expandpath

    ## `Merged:` in newvalue for raw_fastq_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'raw_fastq_path':
        if newvalue and "Merged:" not in newvalue:
            fullrawpath = _expander(newvalue)
            if os.path.isdir(fullrawpath):
                raise IPyradError(RAW_PATH_ISDIR.format(fullrawpath))
            ## if something is found in path
            elif glob.glob(fullrawpath):
                self.paramsdict['raw_fastq_path'] = fullrawpath
            ## else allow empty, tho it can still raise an error in step1
            else:
                raise IPyradError(NO_RAW_FILE.format(fullrawpath))
        else:
            self.paramsdict['raw_fastq_path'] = ""


    ## `Merged:` in newvalue for barcodes_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'barcodes_path':
        ## if a value was entered check that it exists
        if newvalue and "Merged:" not in newvalue:
            ## also allow for fuzzy match in names using glob
            fullbarpath = glob.glob(_expander(newvalue))[0]
            ## raise error if file is not found
            if not os.path.exists(fullbarpath):
                raise IPyradError(BARCODE_NOT_FOUND.format(fullbarpath))
            else:
                self.paramsdict['barcodes_path'] = fullbarpath
                self._link_barcodes()

        ## if no path was entered then set barcodes path to empty.
        ## this is checked again during step 1 and will raise an error
        ## if you try demultiplexing without a barcodes file
        else:
            self.paramsdict['barcodes_path'] = newvalue


    ## `Merged:` in newvalue for sorted_fastq_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'sorted_fastq_path':
        if newvalue and not "Merged:" in newvalue:
            fullsortedpath = _expander(newvalue)

            if os.path.isdir(fullsortedpath):
                raise IPyradError(SORTED_ISDIR.format(fullsortedpath))
            elif glob.glob(fullsortedpath):
                self.paramsdict['sorted_fastq_path'] = fullsortedpath

            else:
                raise IPyradError(SORTED_NOT_FOUND.format(fullsortedpath))
        ## if no value was entered then set to "".
        else:
            self.paramsdict['sorted_fastq_path'] = ""


    elif param == 'assembly_method':
        methods = ["denovo", "reference", "denovo+reference", "denovo-reference"]
        assert newvalue in methods, BAD_ASSEMBLY_METHOD.format(newvalue)
        self.paramsdict['assembly_method'] = newvalue


    elif param == 'reference_sequence':
        if newvalue:
            fullrawpath = _expander(newvalue)
            if not os.path.isfile(fullrawpath):
                ip.logger.info("reference sequence file not found.")
                raise IPyradError(REF_NOT_FOUND.format(fullrawpath))
            self.paramsdict['reference_sequence'] = fullrawpath
        ## if no value was entered the set to "". Will be checked again
        ## at step3 if user tries to do refseq and raise error
        else:
            self.paramsdict['reference_sequence'] = ""


    elif param == 'datatype':
        ## list of allowed datatypes
        datatypes = ['rad', 'gbs', 'ddrad', 'pairddrad',
                     'pairgbs', 'merged', '2brad', 'pair3rad']
        ## raise error if something else
        if str(newvalue) not in datatypes:
            raise IPyradError("""
    datatype {} not recognized, must be one of: {}
    """.format(newvalue, datatypes))
        else:
            self.paramsdict['datatype'] = str(newvalue)
            ## link_barcodes is called before datatypes is set
            ## we need to know the datatype so we can read in
            ## the multiplexed barcodes for 3rad. This seems
            ## a little annoying, but it was better than any
            ## alternatives I could think of.
            if "3rad" in self.paramsdict['datatype'] and not \
            self.paramsdict['sorted_fastq_path'].strip():
                if "Merged:" not in self.paramsdict['barcodes_path']:
                    self._link_barcodes()

    elif param == 'restriction_overhang':
        newvalue = _tuplecheck(newvalue, str)
        assert isinstance(newvalue, tuple), """
    cut site must be a tuple, e.g., (TGCAG, '') or (TGCAG, CCGG)"""
        ## Handle the special case where the user has 1
        ## restriction overhang and does not include the trailing comma
        if len(newvalue) == 1:
            ## for gbs users might not know to enter the second cut site
            ## so we do it for them.
            if self.paramsdict["datatype"] == "gbs":
                newvalue += newvalue
            else:
                newvalue += ("",)

        ## Handle 3rad datatype with only 3 cutters
        if len(newvalue) == 3:
            newvalue = (newvalue[0], newvalue[1], newvalue[2], "")
        assert len(newvalue) <= 4, """
    most datasets require 1 or 2 cut sites, e.g., (TGCAG, '') or (TGCAG, CCGG).
    For 3rad/seqcap may be up to 4 cut sites."""
        self.paramsdict['restriction_overhang'] = newvalue

    elif param == 'max_low_qual_bases':
        assert isinstance(int(newvalue), int), """
    max_low_qual_bases must be an integer."""
        self.paramsdict['max_low_qual_bases'] = int(newvalue)

    elif param == 'phred_Qscore_offset':
        assert isinstance(int(newvalue), int), \
            "phred_Qscore_offset must be an integer."
        self.paramsdict['phred_Qscore_offset'] = int(newvalue)

    elif param == 'mindepth_statistical':
        assert isinstance(int(newvalue), int), \
            "mindepth_statistical must be an integer."
        ## do not allow values below 5
        if int(newvalue) < 5:
            raise IPyradError("""
    mindepth statistical cannot be set < 5. Use mindepth_majrule.""")
        else:
            self.paramsdict['mindepth_statistical'] = int(newvalue)

    elif param == 'mindepth_majrule':
        assert isinstance(int(newvalue), int), \
            "mindepth_majrule must be an integer."
        self.paramsdict['mindepth_majrule'] = int(newvalue)

    elif param == 'maxdepth':
        self.paramsdict['maxdepth'] = int(newvalue)

    elif param == 'clust_threshold':
        newvalue = float(newvalue)
        assert (newvalue < 1) & (newvalue > 0), \
        "clust_threshold must be a decimal value between 0 and 1."
        self.paramsdict['clust_threshold'] = newvalue

    elif param == 'max_barcode_mismatch':
        self.paramsdict['max_barcode_mismatch'] = int(newvalue)

    elif param == 'filter_adapters':
        self.paramsdict['filter_adapters'] = int(newvalue)

    elif param == 'filter_min_trim_len':
        self.paramsdict["filter_min_trim_len"] = int(newvalue)

    elif param == 'max_alleles_consens':
        self.paramsdict['max_alleles_consens'] = int(newvalue)

    elif param == 'max_Ns_consens':
        newvalue = _tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_Ns_consens should be a tuple e.g., (8, 8)"
        self.paramsdict['max_Ns_consens'] = newvalue

    elif param == 'max_Hs_consens':
        newvalue = _tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_Hs_consens should be a tuple e.g., (5, 5)"
        self.paramsdict['max_Hs_consens'] = newvalue

    elif param == 'min_samples_locus':
        self.paramsdict['min_samples_locus'] = int(newvalue)

    elif param == 'max_shared_Hs_locus':
        if isinstance(newvalue, str):
            if newvalue.isdigit():
                newvalue = int(newvalue)
            else:
                try:
                    newvalue = float(newvalue)
                except Exception as inst:
                    raise IPyradParamsError("""
    max_shared_Hs_locus must be int or float, you put: {}""".format(newvalue))
        self.paramsdict['max_shared_Hs_locus'] = newvalue

    elif param == 'max_SNPs_locus':
        newvalue = _tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_SNPs_locus should be a tuple e.g., (20, 20)"
        self.paramsdict['max_SNPs_locus'] = newvalue

    elif param == 'max_Indels_locus':
        newvalue = _tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_Indels_locus should be a tuple e.g., (5, 100)"
        self.paramsdict['max_Indels_locus'] = newvalue

    ## deprecated but retained for legacy, now uses trim_reads (below)
    elif param == 'edit_cutsites':
        ## Force into a string tuple
        newvalue = _tuplecheck(newvalue)
        ## try converting each tup element to ints
        newvalue = list(newvalue)
        for i in range(2):
            try:
                newvalue[i] = int(newvalue[i])
            except (ValueError, IndexError):
                newvalue.append(0)
                pass
        newvalue = tuple(newvalue)
        ## make sure we have a nice tuple
        if not isinstance(newvalue, tuple):
            raise IPyradError("""
    Error: edit_cutsites should be a tuple e.g., (0, 5) or ('TGCAG', 6),
    you entered {}
    """.format(newvalue))
        self.paramsdict['edit_cutsites'] = newvalue

    elif param == 'trim_reads':
        ## Force into a string tuple
        newvalue = _tuplecheck(newvalue)
        ## try converting each tup element to ints
        newvalue = list(newvalue)
        for i in range(4):
            try:
                newvalue[i] = int(newvalue[i])
            except (ValueError, IndexError):
                newvalue.append(0)
                pass
        newvalue = tuple(newvalue)
        ## make sure we have a nice tuple
        if not isinstance(newvalue, tuple):
            raise IPyradError("""
    Error: trim_reads should be a tuple e.g., (0, -5, -5, 0) 
    or (0, 90, 0, 90), or (0, 0, 0, 0). 
    You entered {}\n""".format(newvalue))
        self.paramsdict['trim_reads'] = newvalue        

    ## deprecated but retained for legacy, now named trim_loci 
    elif param == 'trim_overhang':
        newvalue = _tuplecheck(newvalue, str)
        assert isinstance(newvalue, tuple), \
        "trim_overhang should be a tuple e.g., (4, *, *, 4)"
        self.paramsdict['trim_overhang'] = tuple([int(i) for i in newvalue])

    elif param == 'trim_loci':
        newvalue = _tuplecheck(newvalue, str)
        assert isinstance(newvalue, tuple), \
        "trim_overhang should be a tuple e.g., (0, -5, -5, 0)"
        self.paramsdict['trim_loci'] = tuple([int(i) for i in newvalue])


    elif param == 'output_formats':

        ## Handle the case where output formats is an empty string
        if isinstance(newvalue, str):
            ## strip commas and spaces from string so we have only letters
            newvalue = newvalue.replace(",", "").replace(" ", "")
            newvalue = list(newvalue)
            if not newvalue:
                newvalue = ["*"]
        if isinstance(newvalue, tuple):
            newvalue = list(newvalue)

        if isinstance(newvalue, list):
            ## if more than letters, raise an warning
            if any([len(i) > 1 for i in newvalue]):
                ip.logger.warning("""
    'output_formats' params entry is malformed. Setting to * to avoid errors.""")
                newvalue = OUTPUT_FORMATS
            newvalue = tuple(newvalue)
            #newvalue = tuple([i for i in newvalue if i in allowed])
        if "*" in newvalue:
            newvalue = OUTPUT_FORMATS

        ## set the param
        self.paramsdict['output_formats'] = newvalue


    elif param == 'pop_assign_file':
        fullpoppath = _expander(newvalue)

        ## if a path is entered, raise exception if not found
        if newvalue:
            if not os.path.isfile(fullpoppath):
                ip.logger.warn("Population assignment file not found.")
                raise IPyradError("""
    Warning: Population assignment file not found. This must be an
    absolute path (/home/wat/ipyrad/data/my_popfile.txt) or relative to
    the directory where you're running ipyrad (./data/my_popfile.txt)
    You entered: {}\n""".format(fullpoppath))
        ## should we add a check here that all pop samples are in samples?

            self.paramsdict['pop_assign_file'] = fullpoppath
            self._link_populations()
        else:
            self.paramsdict['pop_assign_file'] = ""
            ## Don't forget to possibly blank the populations dictionary
            self.populations = {}

    return self


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
