#!/usr/bin/env ipython2

""" ipyrad Assembly class object. """

# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0142
# pylint: disable=W0212
# pylint: disable=C0301

from __future__ import print_function
import os
import glob
import sys
import gzip
import copy
import h5py
import string
import itertools
import numpy as np
import pandas as pd
import ipyrad as ip
import socket
import time
import datetime

from collections import OrderedDict
from ipyrad.assemble.util import *
from ipyrad.assemble.refmap import index_reference_sequence
from ipyrad.core.paramsinfo import paraminfo, paramname
from ipyrad.core.sample import Sample
from ipyrad import assemble

import logging
LOGGER = logging.getLogger(__name__)

## turn off traceback for the CLI
if not ip.__interactive__:
    sys.tracebacklimit = 0



class Assembly(object):
    """ An ipyrad Assembly class object.

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
    step[x]()
        runs step x of assembly
    toparamsfile(filehandle)
        writes params dictionary to params.txt file format


    Returns
    -------
    object
         A new assembly object is returned.

     """


    def __init__(self, name, quiet=False):

        ## Make sure assembly name is not empty
        if not name:
            raise IPyradParamsError(REQUIRE_ASSEMBLY_NAME)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        invalid_chars = string.punctuation.replace("_", "")\
                                          .replace("-", "")+ " "
        if any(char in invalid_chars for char in name):
            raise IPyradParamsError(BAD_ASSEMBLY_NAME.format(name))

        self.name = name
        if not quiet:
            print("  New Assembly: {}".format(self.name))

        ## Store assembly version #
        self._version = ip.__version__

        ## stores default ipcluster launch info
        self._ipcluster = {
            "cluster_id" : "",
            "profile" : "default",
            "engines" : "Local",
            "quiet" : 0,
            "timeout" : 120,
            "cores" : 0, #detect_cpus(),
            "threads" : 2
            }

        ## print headers, this is used as a 'quiet' option
        ## or to control differences in printing between API and CLI
        self._headers = 0

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

        ## storing supercatg file
        self.clust_database = ""
        self.database = ""

        ## the default params dict
        self.paramsdict = OrderedDict([
                       ("assembly_name", name),
                       ("project_dir", "./"),#os.path.realpath(os.path.curdir)),
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
                       ("edit_cutsites", (0, 0)),
                       ("trim_overhang", (0, 0, 0, 0)),
                       ("output_formats", ['l', 'p', 's', 'v']),
                       ("pop_assign_file", ""),
        ])

        ## Store data directories for this Assembly. Init with default project
        self.dirs = ObjDict({"project":
                              os.path.realpath(self.paramsdict["project_dir"]),
                             "fastqs": "",
                             "edits": "",
                             "clusts": "",
                             "consens": "",
                             "outfiles": ""})

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
                        ("aligner", "bwa")
        ])

    def __str__(self):
        return "<ipyrad.Assembly object {}>".format(self.name)

    @property
    def stats(self):
        """ Returns a data frame with Sample data and state. """
        nameordered = self.samples.keys()
        nameordered.sort()

        ## Set pandas to display all samples instead of truncating
        pd.options.display.max_rows = len(self.samples)
        statdat = pd.DataFrame([self.samples[i].stats for i in nameordered],
                      index=nameordered).dropna(axis=1, how='all')
        # ensure non h,e columns print as ints
        for column in statdat:
            if column not in ["hetero_est", "error_est"]:
                statdat[column] = np.nan_to_num(statdat[column]).astype(int)
        return statdat


    @property
    def files(self):
        """ Returns a data frame with Sample files. Not very readable... """
        nameordered = self.samples.keys()
        nameordered.sort()
        ## replace curdir with . for shorter printing
        #fullcurdir = os.path.realpath(os.path.curdir)
        return pd.DataFrame([self.samples[i].files for i in nameordered],
                      index=nameordered).dropna(axis=1, how='all')


    def build_stat(self, idx):
        """ Returns a data frame with Sample stats for each step """
        nameordered = self.samples.keys()
        nameordered.sort()
        return pd.DataFrame([self.samples[i].stats_dfs[idx] \
                      for i in nameordered],
                      index=nameordered)\
                      .dropna(axis=1, how='all')



    def link_fastqs(self, path=None, force=False, append=False, splitnames="_",
                    fields=None, ipyclient=None):
        """
        Create Sample objects from demultiplexed fastq files in sorted_fastq_path,
        or append additional fastq files to existing Samples. This provides
        more flexible file input through the API than available in step1 of the
        command line interface. If passed ipyclient it will run in parallel.

        Note
        ----
        This function is called during step 1 if files are specified in
        'sorted_fastq_path'.

        Parameters
        ----------
        path : str
            Path to the fastq files to be linked to Sample objects. The default
            location is to select all files in the 'sorted_fastq_path'.
            Alternatively a different path can be entered here.

        append : bool
            The default action is to overwrite fastq files linked to Samples if
            they already have linked files. Use append=True to instead append
            additional fastq files to a Sample (file names should be formatted
            the same as usual, e.g., [name]_R1_[optional].fastq.gz).

        splitnames : str
            A string character used to file names. In combination with the
            fields argument can be used to subselect filename fields names.

        fields : list
            A list of indices for the fields to be included in names after
            filnames are split on the splitnames character. Useful for appending
            sequence names which must match existing names. If the largest index
            is greater than the number of split strings in the name the index
            if ignored. e.g., [2,3,4] ## excludes 0, 1 and >4

        force : bool
            Overwrites existing Sample data and statistics.

        Returns
        -------
        str
            Prints the number of new Sample objects created and the number of
            fastq files linked to Sample objects in the Assembly object.
        """

        ## cannot both force and append at once
        if force and append:
            raise IPyradError("Cannot use force and append at the same time.")

        if self.samples and not (force or append):
            raise IPyradError("Files already linked to `{}`.".format(self.name)\
                +" Use force=True to replace all files, or append=True to add"
                +" additional files to existing Samples.")

        ## make sure there is a workdir and workdir/fastqdir
        self.dirs.fastqs = os.path.join(self.paramsdict["project_dir"],
                                        self.name+"_fastqs")
        if not os.path.exists(self.paramsdict["project_dir"]):
            os.mkdir(self.paramsdict["project_dir"])
        if not os.path.exists(self.dirs.fastqs):
            os.mkdir(self.dirs.fastqs)

        ## get path to data files
        if not path:
            path = self.paramsdict["sorted_fastq_path"]
        #print(LINKING_TO_MSG.format(path))

        ## but grab fastq/fq/gz, and then sort
        fastqs = glob.glob(path)
        fastqs = [i for i in fastqs if i.endswith(".gz") \
                                    or i.endswith(".fastq") \
                                    or i.endswith(".fq")]
        LOGGER.debug("linking these files:\n{}".format(fastqs))
        fastqs.sort()
        LOGGER.debug("Linking these fastq files:\n".format(fastqs))

        ## raise error if no files are found
        if not fastqs:
            raise IPyradError(NO_FILES_FOUND_PAIRS\
                        .format(self.paramsdict["sorted_fastq_path"]))

        ## link pairs into tuples
        if 'pair' in self.paramsdict["datatype"]:
            ## check that names fit the paired naming convention
            r1_files = [i for i in fastqs if "_R1_" in i]
            r2_files = [i.replace("_R1_", "_R2_") for i in r1_files]

            if r1_files:
                if not any(["_R1_" in i for i in fastqs]) or \
                       (len(r1_files) != len(r2_files)):
                    raise IPyradError("""
        Paired file names must be identical except for _R1_ and _R2_.""")

            ## Test R2 files actually exist
            if not all([os.path.exists(x) for x in r2_files]):
                raise IPyradError("""
        Paired file names must be identical except for _R1_ and _R2_.""")

            fastqs = [(i, j) for i, j in zip(r1_files, r2_files)]

        ## data are not paired, create empty tuple pair
        else:
            ## print warning if _R2_ is in names when not paired
            if any(["_R2_" in i for i in fastqs]):
                print(NAMES_LOOK_PAIRED_WARNING)
            fastqs = [(i, "") for i in fastqs]

        ## counters for the printed output
        created = 0
        linked = 0
        appended = 0

        ## clear samples if force
        if force:
            self.samples = {}

        ## track parallel jobs
        linkjobs = {}
        if ipyclient:
            lbview = ipyclient.load_balanced_view()

        ## iterate over input files
        for fastqtuple in list(fastqs):
            assert isinstance(fastqtuple, tuple), "fastqs not a tuple."
            ## local counters
            createdinc = 0
            linkedinc = 0
            appendinc = 0
            ## remove file extension from name
            sname = _name_from_file(fastqtuple[0], splitnames, fields)
            LOGGER.debug("New Sample name {}".format(sname))

            if sname not in self.samples:
                ## create new Sample
                LOGGER.debug("Creating new sample - ".format(sname))
                self.samples[sname] = Sample(sname)
                self.samples[sname].stats.state = 1
                self.samples[sname].barcode = None
                self.samples[sname].files.fastqs.append(fastqtuple)
                createdinc += 1
                linkedinc += 1
            else:
                ## if not forcing, shouldn't be here with existing Samples
                if append:
                    #if fastqtuple not in self.samples[sname].files.fastqs:
                    self.samples[sname].files.fastqs.append(fastqtuple)
                    appendinc += 1

                elif force:
                    ## overwrite/create new Sample
                    LOGGER.debug("Overwriting sample - ".format(sname))
                    self.samples[sname] = Sample(sname)
                    self.samples[sname].stats.state = 1
                    self.samples[sname].barcode = None
                    self.samples[sname].files.fastqs.append(fastqtuple)
                    createdinc += 1
                    linkedinc += 1
                else:
                    print("""
        The files {} are already in Sample. Use append=True to append additional
        files to a Sample or force=True to replace all existing Samples.
        """.format(sname))

            ## support serial execution w/o ipyclient
            if not ipyclient:
                if any([linkedinc, createdinc, appendinc]):
                    gzipped = bool(fastqtuple[0].endswith(".gz"))
                    nreads = 0
                    for alltuples in self.samples[sname].files.fastqs:
                        nreads += zbufcountlines(alltuples[0], gzipped)
                    self.samples[sname].stats.reads_raw = nreads/4
                    LOGGER.debug("Got reads for sample - {} {}".format(sname,\
                                    self.samples[sname].stats.reads_raw))
                    created += createdinc
                    linked += linkedinc
                    appended += appendinc

            ## do counting in parallel
            else:
                if any([linkedinc, createdinc, appendinc]):
                    gzipped = bool(fastqtuple[0].endswith(".gz"))
                    for sidx, tup in enumerate(self.samples[sname].files.fastqs):
                        key = sname+"_{}".format(sidx)
                        linkjobs[key] = lbview.apply(bufcountlines,
                                                    *(tup[0], gzipped))
                    LOGGER.debug("sent count job for {}".format(sname))
                    created += createdinc
                    linked += linkedinc
                    appended += appendinc

        ## wait for link jobs to finish if parallel
        if ipyclient:
            start = time.time()
            while 1:
                fin = [i.ready() for i in linkjobs.values()]
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(len(fin), sum(fin),
                    ' loading reads         | {} | s1 |'.format(elapsed))
                time.sleep(0.1)
                if len(fin) == sum(fin):
                    print("")
                    break

            ## collect link job results
            sampdict = {i:0 for i in self.samples}
            for result in linkjobs:
                sname = result.rsplit("_", 1)[0]
                nreads = linkjobs[result].result()
                sampdict[sname] += nreads

            for sname in sampdict:
                self.samples[sname].stats.reads_raw = sampdict[sname]/4

        ## print if data were linked
        #print("  {} new Samples created in '{}'.".format(created, self.name))
        if linked:
            ## double for paired data
            if 'pair' in self.paramsdict["datatype"]:
                linked = linked*2
            if self._headers:
                print("  {} fastq files loaded to {} Samples.".\
                      format(linked, len(self.samples)))
            ## save the location where these files are located
            self.dirs.fastqs = os.path.realpath(os.path.dirname(path))

        if appended:
            if self._headers:
                print("  {} fastq files appended to {} existing Samples.".\
                      format(appended, len(self.samples)))



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

            ## make sure chars are all proper
            if not all(bdf[1].apply(set("RKSYWMCATG").issuperset)):
                LOGGER.warn(BAD_BARCODE)
                raise IPyradError(BAD_BARCODE)

            ## 3rad/seqcap use multiplexed barcodes
            ## We'll concatenate them with a plus and split them later
            if "3rad" in self.paramsdict["datatype"]:
                try:
                    bdf[2] = bdf[2].str.upper()
                    self.barcodes = dict(zip(bdf[0], bdf[1] + "+" + bdf[2]))
                except KeyError as inst:
                    msg = "    3rad assumes multiplexed barcodes. Doublecheck your barcodes file."
                    LOGGER.error(msg)
                    raise IPyradError(msg)
            else:
                ## set attribute on Assembly object
                self.barcodes = dict(zip(bdf[0], bdf[1]))

        except (IOError, IndexError):
            raise IPyradWarningExit(\
                "    Barcodes file not found. You entered: {}"\
                .format(self.paramsdict["barcodes_path"]))

        except ValueError as inst:
            msg = "    Barcodes file format error."
            LOGGER.warn(msg)
            raise IPyradError(inst)


    def link_populations(self, popdict=None):
        """
        Creates self.populations dictionary to save mappings of individuals to
        populations/sites, and checks that individual names match with Samples.
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

            pops = {pop1: [ind1, ind2, ind3], pop2: [ind4, ind5]}
            [Assembly].link_populations(popdict=pops).


        """
        if not popdict:
            ## glob it in case of fuzzy matching
            popfile = glob.glob(self.paramsdict["pop_assign_file"])[0]
            if not os.path.exists(popfile):
                raise IPyradError("Population assignment file not found: {}"\
                                  .format(self.paramsdict["pop_assign_file"]))

            ## parse populations file
            try:
                popdat = pd.read_csv(popfile, header=None,
                                              delim_whitespace=1,
                                              names=["inds", "pops"])
                popdict = {key: group.inds.values.tolist() for key, group in \
                                                        popdat.groupby("pops")}
            except ValueError:
                LOGGER.warn("Populations file may be malformed.")
                raise IPyradError("Populations file malformed - {}"\
                                  .format(popfile))
        else:
            ## pop dict is provided by user
            pass

        ## filter for bad samples
        ## Warn user but don't bail out, could be setting the pops file
        ## on a new assembly w/o any linked samples.
        badsamples = [i for i in itertools.chain(*popdict.values()) \
                      if i not in self.samples.keys()]
        if any(badsamples):
            LOGGER.warn("Some names from population input do not match Sample "\
                + "names: ".format(", ".join(badsamples)))
            LOGGER.warn("If this is a new assembly this is normal.")

        ## return dict
        self.populations = popdict



    def get_params(self, param=""):
        """ pretty prints params if called as a function """
        fullcurdir = os.path.realpath(os.path.curdir)
        if not param:
            for index, (key, value) in enumerate(self.paramsdict.items()):
                if isinstance(value, str):
                    value = value.replace(fullcurdir, ".")
                sys.stdout.write("  {:<4}{:<28}{:<45}\n".format(index,
                           key, value))
        else:
            try:
                if int(param):
                    #sys.stdout.write(self.paramsdict.values()[int(param)-1])
                    return self.paramsdict.values()[int(param)]
            except (ValueError, TypeError, NameError, IndexError):
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
        ## param 1 takes only a str as input
        [Assembly].set_params(1, 'new_directory')
        [Assembly].set_params('project_dir', 'new_directory')

        ## param 8 must be a tuple or str, if str it is converted to a tuple
        ## with the second entry empty.
        [Assembly].set_params(8, 'TGCAG')
        [Assembly].set_params('restriction_overhang', ('CTGCAG', 'CCGG')

        ## param 24 can be an int or a float:
        [Assembly].set_params(24, 4)
        [Assembly].set_params('max_shared_Hs_locus', 0.25)

        """

        ## require parameter recognition
        if not ((param in range(50)) or \
                (param in [str(i) for i in range(50)]) or \
                (param in self.paramsdict.keys())):
            raise IPyradParamsError("Parameter key not recognized: {}"\
                                    .format(param))

        ## make string
        param = str(param)

        ## get index if param is keyword arg (this index is now zero based!)
        if len(param) < 3:
            param = self.paramsdict.keys()[int(param)]

        ## run assertions on new param
        try:
            self = paramschecker(self, param, newvalue)

        except Exception as inst:
            raise IPyradWarningExit(BAD_PARAMETER\
                                    .format(param, inst, newvalue))



    def write_params(self, outfile=None, force=False):
        """ Write out the parameters of this assembly to a file properly
        formatted as input for `ipyrad -p <params.txt>`. A good and
        simple way to share/archive parameter settings for assemblies.
        This is also the function that's used by __main__ to
        generate default params.txt files for `ipyrad -n`
        """
        if outfile is None:
            outfile = "params-"+self.name+".txt"

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise IPyradWarningExit(PARAMS_EXISTS.format(outfile))

        with open(outfile, 'w') as paramsfile:
            ## Write the header. Format to 80 columns
            header = "------- ipyrad params file (v.{})".format(ip.__version__)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key, val in self.paramsdict.iteritems():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)
                padding = (" "*(30-len(paramvalue)))
                paramkey = self.paramsdict.keys().index(key)
                paramindex = " ## [{}] "\
                             .format(paramkey)
                name = "[{}]: ".format(paramname(paramkey))
                description = paraminfo(paramkey, short=True)
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)



    def branch(self, newname, subsamples=None, infile=None):
        """
        Returns a copy of the Assembly object. Does not allow Assembly
        object names to be replicated in namespace or path.
        """
        ## subsample by removal or keeping.
        remove = 0

        ## is there a better way to ask if it already exists?
        if (newname == self.name or os.path.exists(
                                    os.path.join(self.paramsdict["project_dir"],
                                    newname+".assembly"))):
            print("    Assembly object named {} already exists".format(newname))

        else:
            ## create a copy of the Assembly obj
            newobj = copy.deepcopy(self)
            newobj.name = newname
            newobj.paramsdict["assembly_name"] = newname

            if subsamples and infile:
                print("Attempting to branch passing in subsample names "\
                        "and an input file, ignoring `subsamples` argument")

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
                        print("  Sample name not found: {}".format(sname))
                ## reload sample dict w/o non subsamples
                newobj.samples = {name:sample for name, sample in \
                           newobj.samples.items() if name in subsamples}

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
        #if self._headers:
        #    print("")
        ip.save_json(self)



    def _step1func(self, force, preview, ipyclient):
        """ hidden wrapped function to start step 1 """

        ## check input data files
        sfiles = self.paramsdict["sorted_fastq_path"]
        rfiles = self.paramsdict["raw_fastq_path"]

        ## do not allow both a sorted_fastq_path and a raw_fastq
        if sfiles and rfiles:
            raise IPyradWarningExit(NOT_TWO_PATHS)

        ## but also require that at least one exists
        if not (sfiles or rfiles):
            raise IPyradWarningExit(NO_SEQ_PATH_FOUND)

        ## print headers
        if self._headers:
            if sfiles:
                print("\n  Step 1: Loading sorted fastq data to Samples")
            else:
                print("\n  Step 1: Demultiplexing fastq data to Samples")
        #else:
        #    print("")

        ## if Samples already exist then no demultiplexing
        if self.samples:
            if not force:
                print(SAMPLES_EXIST.format(len(self.samples), self.name))
            else:
                ## overwrite existing data
                if glob.glob(sfiles):
                    if self._headers:
                        self.link_fastqs(ipyclient=ipyclient)

                ## otherwise do the demultiplexing
                else:
                    assemble.demultiplex.run(self, preview, ipyclient, force)

        ## Creating new Samples
        else:
            ## first check if demultiplexed files exist in sorted path
            if glob.glob(sfiles):
                self.link_fastqs(ipyclient=ipyclient)

            ## otherwise do the demultiplexing
            else:
                assemble.demultiplex.run(self, preview, ipyclient, force)



    def _step2func(self, samples, nreplace, force, preview, ipyclient):
        """ hidden wrapped function to start step 2"""

        ## print header
        if self._headers:
            print("\n  Step 2: Filtering reads ")
        #else:
        #    print("")

        ## If no samples in this assembly then it means you skipped step1,
        if not self.samples.keys():
            raise IPyradWarningExit(FIRST_RUN_1)

        ## Get sample objects from list of strings, if API.
        samples = _get_samples(self, samples)

        if not force:
            ## print warning and skip if all are finished
            if all([i.stats.state >= 2 for i in samples]):
                print(EDITS_EXIST.format(len(samples)))
                return

        ## Run samples through rawedit
        assemble.rawedit.run2(self, samples, force, ipyclient)
        #assemble.rawedit.run(self, samples, nreplace, force, preview, ipyclient)



    def _step3func(self, samples, noreverse, maxindels, force, preview, ipyclient):
        """ hidden wrapped function to start step 3 """
        ## print headers
        if self._headers:
            print("\n  Step 3: Clustering/Mapping reads")
        #else:
        #    print("")

        ## Require reference seq for reference-based methods
        if self.paramsdict['assembly_method'] != "denovo":
            if not self.paramsdict['reference_sequence']:
                raise IPyradError(REQUIRE_REFERENCE_PATH\
                            .format(self.paramsdict["assembly_method"]))
            else:
                ## index the reference sequence
                ## Allow force to reindex the reference sequence
                index_reference_sequence(self, force)

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 3, force):
            raise IPyradError(FIRST_RUN_2)

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 3 for i in samples]):
                print(CLUSTERS_EXIST.format(len(samples)))
                return

        ## run the step function
        assemble.cluster_within.run(self, samples, noreverse, maxindels,
                                    force, preview, ipyclient)



    def _step4func(self, samples, subsample, force, ipyclient):
        """ hidden wrapped function to start step 4 """

        if self._headers:
            print("\n  Step 4: Joint estimation of error rate and heterozygosity")
        #else:
        #    print("")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 4, force):
            raise IPyradError(FIRST_RUN_3)

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 4 for i in samples]):
                print(JOINTS_EXIST.format(len(samples)))
                return

        ## send to function
        assemble.jointestimate.run(self, samples, subsample, force, ipyclient)



    def _step5func(self, samples, force, ipyclient):
        """ hidden wrapped function to start step 5 """
        ## print header
        if self._headers:
            print("\n  Step 5: Consensus base calling ")
        #else:
        #    print("")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 5, force):
            raise IPyradError(FIRST_RUN_4)

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 5 for i in samples]):
                print(CONSENS_EXIST.format(len(samples)))
                return
        ## pass samples to rawedit
        assemble.consens_se.run(self, samples, force, ipyclient)



    def _step6func(self, samples, noreverse, force, randomseed, ipyclient):
        """ hidden function to start Step 6"""

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## remove samples that aren't ready
        csamples = self.samples_precheck(samples, 6, force)

        ## print CLI header
        if self._headers:
            print("\n  Step 6: Clustering at {} similarity across {} samples".\
                  format(self.paramsdict["clust_threshold"], len(csamples)))

        ## Check if all/none in the right state
        if not csamples:
            raise IPyradError(FIRST_RUN_5)

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 6 for i in csamples]):
                print(DATABASE_EXISTS.format(len(samples)))
                return

        ## run if this point is reached. We no longer check for existing
        ## h5 file, since checking Sample states should suffice.
        assemble.cluster_across.run(self, csamples, noreverse,
                                    force, randomseed, ipyclient)



    def _step7func(self, samples, force, ipyclient):
        """ Step 7: Filter and write output files """

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        if self._headers:
            print("\n  Step 7: Filter and write output files for {} Samples".\
                  format(len(samples)))

        ## Check if all/none of the samples are in the self.database
        try:
            with h5py.File(self.clust_database, 'r') as io5:
                dbset = set(io5["seqs"].attrs['samples'])
                iset = set([i.name for i in samples])

                ## TODO: Handle the case where dbdiff is not empty?
                ## This might arise if someone tries to branch and remove
                ## samples at step 7.
                dbdiff = dbset.difference(iset)
                idiff = iset.difference(dbset)
                if idiff:
                    print(NOT_CLUSTERED_YET\
                    .format(self.database, ", ".join(list(idiff))))

                    ## The the old way that failed unless all samples were
                    ## clustered successfully in step 6. Adding some flexibility
                    ## to allow writing output even if some samples failed.
                    ## raise IPyradWarningExit(msg)

                    ## Remove the samples that aren't ready for writing out
                    ## i.e. only proceed with the samples that are actually
                    ## present in the db
                    samples = [x for x in samples if x.name not in idiff]
        except (IOError, ValueError):
            raise IPyradError(FIRST_RUN_6.format(self.database))

        if not force:
            if os.path.exists(
                os.path.join(self.dirs.project, self.name+"_outfiles")):
                raise IPyradWarningExit(OUTPUT_EXISTS\
               .format(os.path.join(self.dirs.project, self.name+"_outfiles")))

        ## Run step7
        assemble.write_outfiles.run(self, samples, force, ipyclient)



    def samples_precheck(self, samples, mystep, force):
        """ Return a list of samples that are actually ready for the next step.
            Each step runs this prior to calling run, makes it easier to
            centralize and normalize how each step is checking sample states.

            mystep is the state produced by the current step.
        """
        subsample = []
        ## filter by state
        for sample in samples:
            if sample.stats.state < mystep - 1:
                LOGGER.debug("Sample {} not in proper state."\
                             .format(sample.name))
            else:
                subsample.append(sample)
        return subsample


    def compatible_params_check(self):
        """ check for mindepths after all params are set, b/c doing it while each
        is being set becomes complicated """

        ## do not allow statistical < majrule
        val1 = self.paramsdict["mindepth_statistical"]
        val2 = self.paramsdict['mindepth_majrule']
        if val1 < val2:
            msg = """
    Warning: mindepth_statistical cannot not be < mindepth_majrule.
    Forcing mindepth_majrule = mindepth_statistical = {}
    """.format(val1)
            LOGGER.warning(msg)
            print(msg)
            self.paramsdict["mindepth_majrule"] = val1



    def run(self, steps=0, force=False, preview=False, quiet=1):
        """
        Run assembly steps of an ipyrad analysis. Enter steps as a string,
        e.g., "1", "123", "12345". This step checks for an existing
        ipcluster instance otherwise it raises an exception. The ipyparallel
        connection is made using information from the _ipcluster dict of the
        Assembly class object.
        """
        ## check that mindepth params are compatible, fix and report warning.
        self.compatible_params_check()

        ## wrap everything in a try statement to ensure that we save the
        ## Assembly object if it is interrupted at any point, and also
        ## to ensure proper cleanup of the ipyclient.
        try:
            ## use an existing ipcluster instance
            ipyclient = ip.core.parallel.get_client(**self._ipcluster)

            ## print a message about the cluster status
            ## if MPI setup then we are going to wait until all engines are
            ## ready so that we can print how many cores started on each
            ## host machine exactly.
            if not quiet:
                if (self._ipcluster["profile"] != "default") or \
                   (self._ipcluster["engines"] == "MPI"):
                    hosts = ipyclient[:].apply_sync(socket.gethostname)
                    for hostname in set(hosts):
                        print("  host compute node: [{} cores] on {}"\
                              .format(hosts.count(hostname), hostname))

                ## if Local setup then we know that we can get all the cores for
                ## sure and we won't bother waiting for them to start, since
                ## they'll start grabbing jobs once they're started.
                else:
                    ## If `cores` is set then honor this request, else use all
                    ## available cores.
                    if self._ipcluster["cores"]:
                        _cpus = self._ipcluster["cores"]
                    else:
                        _cpus = detect_cpus()
                    print("  local compute node: [{} cores] on {}"\
                          .format(_cpus, socket.gethostname()))

            ## get the list of steps to run
            if isinstance(steps, int):
                steps = str(steps)
            steps = list(steps)
            steps.sort()

            ## print an Assembly name header if inside API
            if ip.__interactive__:
                print("\n  Assembly: {}".format(self.name))

            ## store pids in case we need to die hard (w/ a vengeance)
            #pids = ipyclient[:].apply_sync(os.getpid)

            ## has many fixed arguments right now, but we may add these to
            ## hackerz_only, or they may be accessed in the API.
            if '1' in steps:
                self._step1func(force, preview, ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '2' in steps:
                self._step2func(samples=None, nreplace=1, force=force,
                                preview=preview, ipyclient=ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '3' in steps:
                self._step3func(samples=None, noreverse=0, force=force,
                             maxindels=8, preview=preview, ipyclient=ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '4' in steps:
                self._step4func(samples=None, subsample=9999999, force=force,
                                ipyclient=ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '5' in steps:
                self._step5func(samples=None, force=force, ipyclient=ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '6' in steps:
                self._step6func(samples=None, noreverse=0, randomseed=12345,
                                force=force, ipyclient=ipyclient)
                self.save()
                ipyclient.purge_everything()

            if '7' in steps:
                self._step7func(samples=None, force=force, ipyclient=ipyclient)
                ipyclient.purge_everything()


        ## handle exceptions so they will be raised after we clean up below
        except KeyboardInterrupt as inst:
            LOGGER.info("assembly interrupted by user.")
            print("\n  Keyboard Interrupt by user. Cleaning up...")

        except IPyradWarningExit as inst:
            LOGGER.error("IPyradWarningExit: %s", inst)
            print("\n  Encountered an error, see ./ipyrad_log.txt. \n  {}"\
                  .format(inst))

        except Exception as inst:
            LOGGER.error(inst)
            print("\n  Encountered an unexpected error (see ./ipyrad_log.txt)"+\
                  "\n  Error message is below -------------------------------"+\
                  "\n{}".format(inst))

        ## close client when done or interrupted
        finally:
            try:
                ## save the Assembly
                self.save()

                ## can't close client if it was never open
                if ipyclient:

                    ## if CLI (has cluster_id), stop jobs and close
                    if self._ipcluster["cluster_id"]:
                        ## protect from KBD while killing jobs?
                        LOGGER.info("  shutting down engines")
                        #_cleanup_and_die(pids)
                        ipyclient.abort()
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        LOGGER.info("  finished shutdown")
                    ## elif API, stop jobs and clean queue but don't close
                    else:
                        ipyclient.abort()
                        ## purge fails if jobs are still running. We could do
                        ## a wait call, but what if its a really long running
                        ## job?. We could do os.kill, but that also kills the
                        ## engines. it seems there is a shutdown function in
                        ## the works that allows for engine restart, but its
                        ## not available yet. This is problem for the API.
                        if not ipyclient.outstanding:
                            ipyclient.purge_everything()
                        else:
                            ## nanny: kill the engines left running, report
                            ## that some engines were killed.
                            pass
                ## a final spacer
                if self._headers:
                    print("")

            ## if exception in close and save, print and ignore
            except Exception as inst2:
                LOGGER.warning("\
            error during ipcluster shutdown (%s)\
            some Python processes may have been orphaned and should be killed"
            , inst2)





def _cleanup_and_die(pids):
    """
    When engines are running external bins like vsearch they can't hear
    KeyboardInterrupts, and so they aren't killed properly. We can save the
    pids of Engines and kill them here instead.
    """
    ## get pids of engines
    for pid in pids:
        try:
            os.kill(pid, 9)
        except OSError:
            pass



def _get_samples(self, samples):
    """
    Internal function. Prelude for each step() to read in perhaps
    non empty list of samples to process. Input is a list of sample names,
    output is a list of sample objects."""
    ## if samples not entered use all samples
    if not samples:
        samples = self.samples.keys()

    ## Be nice and allow user to pass in only one sample as a string,
    ## rather than a one element list. When you make the string into a list
    ## you have to wrap it in square braces or else list makes a list of
    ## each character individually.
    if isinstance(samples, str):
        samples = list([samples])

    ## if sample keys, replace with sample obj
    assert isinstance(samples, list), \
        "to subselect samples enter as a list, e.g., [A, B]."
    newsamples = [self.samples.get(key) for key in samples \
                  if self.samples.get(key)]
    strnewsamples = [i.name for i in newsamples]

    ## are there any samples that did not make it into the dict?
    badsamples = set(samples).difference(set(strnewsamples))
    if badsamples:
        outstring = ", ".join(badsamples)
        raise IPyradError(\
        "Unrecognized Sample name(s) not linked to {}: {}"\
        .format(self.name, outstring))

    ## require Samples
    assert newsamples, \
           "No Samples passed in and none in assembly {}".format(self.name)

    return newsamples



def _name_from_file(fname, splitnames, fields):
    """ internal func: get the sample name from any pyrad file """
    ## allowed extensions
    file_extensions = [".gz", ".fastq", ".fq", ".fasta",
                       ".clustS", ".consens"]
    base, _ = os.path.splitext(os.path.basename(fname))

    ## remove read number from name
    base = base.replace("_R1_.", ".")\
               .replace("_R1_", "_")\
               .replace("_R1.", ".")

    ## remove extensions
    tmpb, tmpext = os.path.splitext(base)
    while tmpext in file_extensions:
        tmpb, tmpext = os.path.splitext(tmpb)
        base = tmpb

    if fields:
        namebits = base.split(splitnames)
        base = []
        for field in fields:
            try:
                base.append(namebits[field])
            except IndexError:
                pass
        base = splitnames.join(base)

    if not base:
        raise IPyradError("""
    Found invalid/empty filename in link_fastqs. Check splitnames argument.
    """)

    return base



def _read_sample_names(fname):
    """ Read in sample names from a plain text file. This is a convenience
    function for branching so if you have tons of sample names you can
    pass in a file rather than having to set all the names at the command
    line.
    """
    try:
        infile = open(fname, 'r')
        subsamples = [x.split()[0] for x in infile.readlines()]

    except Exception as inst:
        print("Failed to read input file with sample names.\n{}".format(inst))
        raise inst

    return subsamples



def expander(namepath):
    """ expand ./ ~ and ../ designators in location names """
    if "~" in namepath:
        namepath = os.path.expanduser(namepath)
    else:
        namepath = os.path.abspath(namepath)
    return namepath



def merge(name, assemblies):
    """
    Creates and returns a new Assembly object in which
    samples from two or more Assembly objects with matching names
    are 'merged'. Merging does not affect the actual files written
    on disk, but rather creates new Samples that are linked to
    multiple data files, and with stats summed.
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

    ## return the new Assembly object
    merged.save()
    return merged



def bufcountlines(filename, gzipped):
    """
    fast line counter. Used to quickly sum number of input reads when running
    link_fastqs to append files. """
    if gzipped:
        fin = gzip.open(filename)
    else:
        fin = open(filename)
    nlines = 0
    buf_size = 1024 * 1024
    read_f = fin.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        nlines += buf.count('\n')
        buf = read_f(buf_size)
    fin.close()
    return nlines


## Tried this out but it's slower than bufcountlines
def zbufcountlines(filename, gzipped):
    """ faster line counter """
    if gzipped:
        cmd1 = ["gunzip", "-c", filename]
    else:
        cmd1 = ["cat", filename]
    cmd2 = ["wc"]

    proc1 = sps.Popen(cmd1, stdout=sps.PIPE, stderr=sps.PIPE)
    proc2 = sps.Popen(cmd2, stdin=proc1.stdout, stdout=sps.PIPE, stderr=sps.PIPE)
    res = proc2.communicate()[0]
    if proc2.returncode:
        raise IPyradWarningExit("error zbufcountlines {}:".format(res))
    LOGGER.info(res)
    nlines = int(res.split()[0])
    return nlines



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
            LOGGER.info("Assembly.tuplecheck() failed to cast to {} - {}"\
                        .format(dtype, newvalue))
            raise

        except Exception as inst:
            LOGGER.info(inst)
            raise SystemExit(\
            "\nError: Param`{}` is not formatted correctly.\n({})\n"\
                 .format(newvalue, inst))

    return newvalue



def paramschecker(self, param, newvalue):
    """ Raises exceptions when params are set to values they should not be"""

    if param == 'assembly_name':
        ## Make sure somebody doesn't try to change their assembly_name, bad
        ## things would happen. Calling set_params on assembly_name only raises
        ## an informative error. Assembly_name is set at Assembly creation time
        ## and is immutable.
        raise IPyradWarningExit("""
    Warning: Assembly name is set at Assembly creation time and is an immutable
    property: You may, however, branch the assembly which will create a copy
    with a new name, a sort of roundabout name change. Here's how:

    Command Line Interface:
        ipyrad -p old_name-params.txt --branch new_name

    API (Jupyter Notebook Users):
        new_assembly = my_assembly.copy("new_name")
    """)

    elif param == 'project_dir':
        expandpath = expander(newvalue)
        if not expandpath.startswith("/"):
            if os.path.exists(expandpath):
                #expandpath = "./"+expandpath
                expandpath = expander(expandpath)
        ## Forbid spaces in path names
        if " " in expandpath:
            raise IPyradWarningExit("""
    Error: Your project_dir contains a directory with a space in the name.
    This can cause all kinds of funny problems so please rename this
    directory and remove the space. Try replacing the space with an underscore.
    You entered: {}
    """.format(expandpath))
        self.paramsdict["project_dir"] = expandpath
        self.dirs["project"] = expandpath

    ## `Merged:` in newvalue for raw_fastq_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'raw_fastq_path':
        if newvalue and not "Merged:" in newvalue:
            fullrawpath = expander(newvalue)
            if os.path.isdir(fullrawpath):
                raise IPyradWarningExit("""
    Error: You entered the path to a directory for raw_fastq_path. To
    ensure the correct files in the directory are selected, please use a
    wildcard selector to designate the desired files.
    Example: /home/user/data/*.fastq  ## selects all files ending in '.fastq'
    You entered: {}
    """.format(fullrawpath))
            ## if something is found in path
            elif glob.glob(fullrawpath):
                self.paramsdict['raw_fastq_path'] = fullrawpath
            ## else allow empty, tho it can still raise an error in step1
            else:
                raise IPyradWarningExit(NO_RAW_FILE.format(fullrawpath))
        else:
            self.paramsdict['raw_fastq_path'] = ""


    ## `Merged:` in newvalue for barcodes_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'barcodes_path':
        ## if a value was entered check that it exists
        if newvalue and not "Merged:" in newvalue:
            ## also allow for fuzzy match in names using glob
            fullbarpath = glob.glob(expander(newvalue))[0]
            ## raise error if file is not found
            if not os.path.exists(fullbarpath):
                raise IPyradWarningExit("""
    Error: barcodes file not found. This must be an absolute path
    (/home/wat/ipyrad/data/data_barcodes.txt) or relative to the directory
    where you're running ipyrad (./data/data_barcodes.txt). You entered:
    {}
    """.format(fullbarpath))
            else:
                self.paramsdict['barcodes_path'] = fullbarpath
                self._link_barcodes()

        ## if no path was entered then set barcodes path to empty.
        ## this is checked again during step 1 and will raise an error
        ## if you try demultiplexing without a barcodes file
        else:
            self.paramsdict['barcodes_path'] = ""


    ## `Merged:` in newvalue for sorted_fastq_path indicates that this
    ## assembly is a merge of several others, so this param has no
    ## value for this assembly
    elif param == 'sorted_fastq_path':
        if newvalue and not "Merged:" in newvalue:
            fullsortedpath = expander(newvalue)

            if os.path.isdir(fullsortedpath):
                raise IPyradWarningExit("""
    Error: You entered the path to a directory for sorted_fastq_path. To
    ensure the correct files in the directory are selected, please use a
    wildcard selector to designate the desired files.
    Example: /home/user/data/*.fastq   ## selects all files ending in '.fastq'
    You entered: {}
    """.format(fullsortedpath))
            elif glob.glob(fullsortedpath):
                self.paramsdict['sorted_fastq_path'] = fullsortedpath

            else:
                raise IPyradWarningExit("""
    Error: fastq sequence files in sorted_fastq_path could not be found.
    Please check that the location was entered correctly and that a wild
    card selector (*) was used to select all or a subset of files.
    You entered: {}
    """.format(fullsortedpath))
        ## if no value was entered then set to "".
        else:
            self.paramsdict['sorted_fastq_path'] = ""



    elif param == 'assembly_method':
        methods = ["denovo", "reference", "denovo+reference", "denovo-reference"]
        assert newvalue in methods, """
    The assembly_method parameter must be one of the following: denovo, reference,
    denovo+reference or denovo-reference. You entered: \n{}.
    """.format(newvalue)
        self.paramsdict['assembly_method'] = newvalue


    elif param == 'reference_sequence':
        if newvalue:
            fullrawpath = expander(newvalue)
            if not os.path.isfile(fullrawpath):
                LOGGER.info("reference sequence file not found.")
                raise IPyradWarningExit("""
    "Warning: reference sequence file not found. This must be an absolute path
    (/home/wat/ipyrad/data/reference.gz) or relative to the directory where
    you're running ipyrad (./data/reference.gz). You entered:
    {}
    """.format(fullrawpath))
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
                self._link_barcodes()

    elif param == 'restriction_overhang':
        newvalue = tuplecheck(newvalue, str)
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
        #=======
        #    newvalue = (newvalue[0], "")
        #>>>>>>> d40a5d5086a0d0aace04dd08338ec4ba5341d1f2

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
        self.paramsdict['clust_threshold'] = float(newvalue)

    elif param == 'max_barcode_mismatch':
        self.paramsdict['max_barcode_mismatch'] = int(newvalue)

    elif param == 'filter_adapters':
        self.paramsdict['filter_adapters'] = int(newvalue)

    elif param == 'filter_min_trim_len':
        self.paramsdict['filter_min_trim_len'] = int(newvalue)

    elif param == 'max_alleles_consens':
        self.paramsdict['max_alleles_consens'] = int(newvalue)

    elif param == 'max_Ns_consens':
        newvalue = tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_Ns_consens should be a tuple e.g., (8, 8)"
        self.paramsdict['max_Ns_consens'] = newvalue

    elif param == 'max_Hs_consens':
        newvalue = tuplecheck(newvalue, int)
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
        newvalue = tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_SNPs_locus should be a tuple e.g., (20, 20)"
        self.paramsdict['max_SNPs_locus'] = newvalue

    elif param == 'max_Indels_locus':
        newvalue = tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "max_Indels_locus should be a tuple e.g., (5, 100)"
        self.paramsdict['max_Indels_locus'] = newvalue

    elif param == 'edit_cutsites':
        ## Force into a string tuple
        newvalue = tuplecheck(newvalue)
        ## try converting each tup element to ints
        newvalue = list(newvalue)
        for i in range(2):
            try:
                newvalue[i] = int(newvalue[i])
            except (ValueError, IndexError):
                pass
        newvalue = tuple(newvalue)
        ## make sure we have a nice tuple
        if not isinstance(newvalue, tuple):
            raise IPyradWarningExit("""
    Error: edit_cutsites should be a tuple e.g., (0, 5) or ('TGCAG', 6),
    you entered {}
    """.format(newvalue))

        self.paramsdict['edit_cutsites'] = newvalue

    elif param == 'trim_overhang':
        newvalue = tuplecheck(newvalue, str)
        assert isinstance(newvalue, tuple), \
        "trim_overhang should be a tuple e.g., (4, *, *, 4)"
        self.paramsdict['trim_overhang'] = tuple([int(i) for i in newvalue])


    elif param == 'output_formats':
        ## let's get whatever the user entered as a tuple of letters
        allowed = assemble.write_outfiles.OUTPUT_FORMATS.keys()

        if "*" in newvalue:
            newvalue = allowed
        if isinstance(newvalue, tuple):
            newvalue = list(newvalue)
        if isinstance(newvalue, str):
            newvalue = [i.strip() for i in newvalue.split(",")]
        if isinstance(newvalue, list):
            ## if more than letters, raise an warning
            if any([len(i) > 1 for i in newvalue]):
                LOGGER.warning("""
    'output_formats' params entry is malformed. Setting to * to avoid errors.""")
                newvalue = allowed
            newvalue = tuple(newvalue)
            #newvalue = tuple([i for i in newvalue if i in allowed])

        ## set the param
        self.paramsdict['output_formats'] = newvalue


    elif param == 'pop_assign_file':
        fullpoppath = expander(newvalue)

        ## if a path is entered, raise exception if not found
        if newvalue:
            if not os.path.isfile(fullpoppath):
                LOGGER.warn("Population assignment file not found.")
                raise IPyradWarningExit("""
    Warning: Population assignment file not found. This must be an
    absolute path (/home/wat/ipyrad/data/my_popfile.txt) or relative to
    the directory where you're running ipyrad (./data/my_popfile.txt)
    You entered: {}\n""".format(fullpoppath))
        ## should we add a check here that all pop samples are in samples?

            self.paramsdict['pop_assign_file'] = fullpoppath
            self.link_populations()
        else:
            self.paramsdict['pop_assign_file'] = ""


    return self




### ERROR MESSAGES ###################################
REQUIRE_ASSEMBLY_NAME = """\
    Assembly name _must_ be set. This is the first parameter in the params.txt
    file, and will be used as a prefix for output files. It should be a short
    string with no special characters, i.e., not a path (no \"/\" characters).
    If you need a suggestion, name it after the organism you're working on.
    """
REQUIRE_REFERENCE_PATH = """\
    Assembly method {} requires that you enter a 'reference_sequence_path'.
    """
BAD_ASSEMBLY_NAME = """\
    No spaces or special characters are allowed in the assembly name. A good
    practice is to replace spaces with underscores '_'. An example of a good
    assembly_name is: white_crowned_sparrows. Here's what you put:
    {}
    """
NO_FILES_FOUND_PAIRS = """\
    No files found in 'sorted_fastq_path': {}
    Check that file names match the required convention for paired datatype
    i.e., paired file names should be identical save for _R1_ and _R2_
    (note the underscores before AND after R*).
    """
NAMES_LOOK_PAIRED_WARNING = """\
    Warning: '_R2_' was detected in a file name, which suggests the data may
    be paired-end. If so, you should set the parameter 'datatype' to a paired
    option (e.g., pairddrad or pairgbs) and re-run step 1, which will require
    using the force flag (-f) to overwrite existing data.
    """
BAD_BARCODE = """\
    One or more barcodes contain invalid IUPAC nucleotide code characters.
    Barcodes must contain only characters from this list "RKSYWMCATG".
    Doublecheck your barcodes file is properly formatted.
    """


BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """
NOT_TWO_PATHS = """\
    Error: Must enter a raw_fastq_path or sorted_fastq_path, but not both.
    """
NO_SEQ_PATH_FOUND = """\
    Error: Step 1 requires that you enter one of the following:
        (1) a sorted_fastq_path
        (2) a raw_fastq_path + barcodes_path
    """
PARAMS_EXISTS = """\
    Params file already exists: {}
    Use force argument to overwrite.
    """
SAMPLES_EXIST = """\
    Skipping: {} Samples already found in Assembly {}.
    (can overwrite with force argument)\
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
    You entered - {}
    """

LINKING_TO_MSG = """\
  Loading data from: {}"""

########################################################



if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.get_params()
    print(DATA.log)
