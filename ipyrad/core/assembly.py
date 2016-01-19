#!/usr/bin/env ipython2

""" ipyrad Assembly class object. """

# pylint: disable=E1101
# pylint: disable=E1103
# pylint: disable=W0142
# pylint: disable=W0212


from __future__ import print_function
import os
import time
import glob
import sys
import gzip
import dill
import copy
import itertools
import subprocess
import pandas as pd
import ipyparallel as ipp
import ipyrad as ip

from collections import OrderedDict
from ipyrad.assemble.util import *
from ipyrad.assemble.refmap import index_reference_sequence
from ipyrad.core.paramsinfo import paraminfo
from ipyrad.core.sample import Sample
from ipyrad import assemble

import logging
LOGGER = logging.getLogger(__name__)



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
    bins : dict
        Keys: values for the path to vsearch, muscle, smalt, samtools
        executables. 
    log : list
        A list of all modifications to the Assembly object and its Samples with
        time stamps. Use `print [Assembly].log` for easier viewing.
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

        ## obj name
        self.name = name
        if not quiet:
            print("  New Assembly: {}".format(self.name))

        ## Store assembly version #
        self._version = ip.__version__ 

        ## stores ipcluster launch info
        self._ipclusterid = ""
        self._ipprofile = ""

        ## print headers
        self._headers = 0

        ## get binaries of dependencies
        self.bins = ObjDict()
        binnames = ["vsearch", "muscle", "smalt", "samtools", "bedtools"]
        for binn, binx in zip(binnames, getbins()):
            self.bins[binn] = binx

        ## link a log history of executed workflow
        self.log = []
        self._stamp(self.name+" created")

        ## statsfiles is a dict where keys return a func... 
        ## can't get this to work with a @property func.
        self.statsfiles = ObjDict({})
        #                            {"s1": self.statsfile("s1"), 
        #                            "s2": self.statsfile("s2"), 
        #                            "s3": self.statsfile("s3"), 
        #                            "s4": self.statsfile("s4"), 
        #                            "s5": self.statsfile("s5"), 
        #                            "s6": self.statsfile("s6"), 
        #                            "s7": self.statsfile("s7")})

        ## samples linked 
        self.samples = ObjDict()

        ## samples linked 
        self.populations = ObjDict()

        ## multiplex files linked
        self.barcodes = ObjDict()

        ## outfiles locations
        self.outfiles = ObjDict()

        ## storing supercatg file
        self.database = ""

        ## the default params dict
        self.paramsdict = OrderedDict([
                       ("working_directory", os.path.realpath(
                                                os.path.curdir)),
                       ("raw_fastq_path", os.path.join(
                                            os.path.realpath(
                                                 os.path.curdir),
                                                 "*.fastq")),
                       ("barcodes_path", os.path.join(
                                            os.path.realpath(
                                                os.path.curdir),
                                                "*.barcodes.txt")),
                       ("sorted_fastq_path", ""),
                       ("assembly_method", "denovo"),
                       ("reference_sequence", ""), 
                       ("datatype", 'rad'), 
                       ("restriction_overhang", ("TGCAG", "")),
                       ("max_low_qual_bases", 5),
                       ("phred_Qscore_offset", 33),
                       ("mindepth_statistical", 6), 
                       ("mindepth_majrule", 6), 
                       ("maxdepth", 1000),
                       ("clust_threshold", .85),
                       ("max_barcode_mismatch", 1),
                       ("filter_adapters", 0), 
                       ("filter_min_trim_len", 35), 
                       ("max_alleles_consens", 2), 
                       ("max_Ns_consens", (5, 5)), 
                       ("max_Hs_consens", (8, 8)), 
                       ("min_samples_locus", 4), 
                       ("max_SNPs_locus", (100, 100)), 
                       ("max_Indels_locus", (5, 99)), 
                       ("max_shared_Hs_locus", .25), 
                       ("edit_cutsites", (0, 0)),
                       ("trim_overhang", (1, 2, 2, 1)),                        
                       ("output_formats", "*"),
                       ("pop_assign_file", ""),
                       ("excludes", ""),
                       ("outgroups", ""),
        ])

        ## Store data directories for this Assembly. Init with default working.
        self.dirs = ObjDict({"working": self.paramsdict["working_directory"]})

        ## Default hackers only parameters dictionary
        ## None of the safeguards of the paramsdict, no nice accessor
        ## functions, so you have to be sure of what you're doing if
        ## you change these values.
        self._hackersonly = OrderedDict([
                        ("random_seed", 42),
                        ("max_fragment_length", 150),
                        ("max_inner_mate_distance", 60),
                        ("preview_truncate_length", 500000),
                        ("output_loci_name_buffer", 5),
                        ("query_cov", None)
        ])

    def __str__(self):
        return "<ipyrad.Assembly object {}>".format(self.name)

    @property
    def stats(self):
        """ Returns a data frame with Sample data and state. """
        nameordered = self.samples.keys()
        nameordered.sort()
        return pd.DataFrame([self.samples[i].stats for i in nameordered], 
                      index=nameordered).dropna(axis=1, how='all')
                      #dtype=[int, int, int, int, int, float, float, int])

    #@property
    def statsfile(self, idx):
        """ Returns a data frame with Sample stats for each step """
        nameordered = self.samples.keys()
        nameordered.sort()
        return pd.DataFrame([self.samples[i].statsfiles[idx] \
                      for i in nameordered], 
                      index=nameordered).dropna(axis=1, how='all')

    @property
    def s1(self):
        """ Returns a data frame with Sample stats for step1 """
        nameordered = self.samples.keys()
        nameordered.sort()
        return pd.DataFrame([self.samples[i].statsfiles["s1"] \
                      for i in nameordered], 
                      index=nameordered).dropna(axis=1, how='all')

    @property
    def files(self):
        """ Returns a data frame with Sample files. Not very readable... """
        nameordered = self.samples.keys()
        nameordered.sort()
        ## replace curdir with . for shorter printing
        #fullcurdir = os.path.realpath(os.path.curdir)
        return pd.DataFrame([self.samples[i].files for i in nameordered], 
                      index=nameordered).dropna(axis=1, how='all')


    def cpus(self):
        """ View the connection  """



    def _stamp(self, event):
        """ Stamps an event into the log history. """
        tev = time.strftime("%m/%d/%y %H:%M:%S", time.gmtime())
        self.log.append((self.name, tev, event))



    def link_fastqs(self, path=None, merged=False, force=False, append=False):
        """ Create Sample objects from demultiplexed fastq files in 
        sorted_fastq_path, or append additional fastq files to 
        existing Samples.

        Note
        ----
        link_fastqs() is called automatically during step2() if no Samples
        are yet present in the Assembly object (data were not demultiplexed
        in step1(). It looks for demultiplexed data files located in the
        `sorted_fastq_path`.

        Parameters
        ----------
        path : str
            Path to the fastq files to be linked to Sample objects. The default
            location is to select all files in the 'sorted_fastq_path'. 
            Alternatively a different path can be entered here. 

        merged : bool
            Set to True if files represent first and second reads that were 
            merged using some external software such as `PEAR` or `VSEARCH`. 

        append : bool
            The default action is to overwrite fastq files linked to Samples if 
            they already have linked files. Use append=True to instead append 
            additional fastq files to a Sample (file names should be formatted 
            the same as usual, e.g., [name]_R1_[optional].fastq.gz).

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
        self.dirs.fastqs = os.path.join(self.paramsdict["working_directory"],
                                        self.name+"_fastqs")
        if not os.path.exists(self.paramsdict["working_directory"]):
            os.mkdir(self.paramsdict["working_directory"])
        if not os.path.exists(self.dirs.fastqs):
            os.mkdir(self.dirs.fastqs)

        ## get path to data files
        if not path:
            path = self.paramsdict["sorted_fastq_path"]

        ## does location exist, if no files selected, try selecting all
        if os.path.isdir(path):
            path += "*"

        ## but grab fastq/fq/gz, and then sort
        fastqs = glob.glob(path)
        fastqs = [i for i in fastqs if i.endswith(".gz") \
                                    or i.endswith(".fastq") \
                                    or i.endswith(".fq")]
        fastqs.sort()

        ## raise error if no files are found
        if not fastqs:
            raise IPyradError("No files found in `sorted_fastq_path`: {}".
                              format(self.paramsdict["sorted_fastq_path"]))

        ## link pairs into tuples        
        if 'pair' in self.paramsdict["datatype"]:
            ## check that names fit the paired naming convention
            r1_files = [i for i in fastqs if "_R1_" in i]
            r2_files = [i.replace("_R1_", "_R2_") for i in r1_files]

            if r1_files:
                if not any(["_R1_" in i for i in fastqs]) or \
                       (len(r1_files) != len(r2_files)):
                    raise IPyradError(\
               "Paired file names must be identical except for _R1_ and _R2_")
            fastqs = [(i, j) for i, j in zip(r1_files, r2_files)]

        ## data are not paired, create empty tuple pair
        else:
            ## print warning if _R2_ is in names when not paired
            if any(["_R2_" in i for i in fastqs]):
                print("""
        Warning: '_R2_' was detected in a file name, which suggests the data
        may be paired-end. If so, you should set the Assembly parameter
        `datatype` to a paired type (e.g., pairddrad or pairgbs) and run
        link_fastqs(force=True) to re-link fastq data.
        """)
            fastqs = [(i, ) for i in fastqs]

        ## counters for the printed output
        created = 0
        linked = 0
        appended = 0

        ## clear samples if force
        if force:
            self.samples = {}

        ## iterate over input files
        for fastqtuple in list(fastqs):
            assert isinstance(fastqtuple, tuple), "fastqs not a tuple."
            ## local counters
            createdinc = 0
            linkedinc = 0
            appendinc = 0
            ## remove file extension from name
            sname = _name_from_file(fastqtuple[0])

            if sname not in self.samples:
                ## create new Sample
                self.samples[sname] = Sample(sname)
                self.samples[sname].stats.state = 1
                self.samples[sname].barcode = None 
                self.samples[sname].files.fastqs.append(fastqtuple)
                createdinc += 1
                linkedinc += 1
            else:
                ## if not forcing, shouldn't be here with existing Samples
                if append:
                    if fastqtuple not in self.samples[sname].files.fastqs:
                        self.samples[sname].files.fastqs.append(fastqtuple)
                        appendinc += 1
                    else:
                        print("The files {} are already in Sample {}, "\
                              .format(fastqtuple, sname) \
                              +"cannot append duplicate files to a Sample.\n")
                elif force:
                    ## create new Sample
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

            ## record whether data were merged.
            if merged:
                self.samples[sname].merged = 1

            ## do not allow merged=False and .forward in file names
            if (merged == False) and ('forward' in fastqtuple[0]):
                print("""
        Warning: If R1 and R2 data are merged use link_fastqs(merge=True) to 
        indicate this. You may need force=True to overwrite existing files.
        """)

            ## if fastqs already demultiplexed, try to link stats
            if any([linkedinc, createdinc, appendinc]):
                gzipped = bool(fastqtuple[0].endswith(".gz"))
                nreads = 0
                ## iterate over files if there are multiple
                for alltuples in self.samples[sname].files.fastqs:
                    nreads += bufcountlines(alltuples[0], gzipped)
                self.samples[sname].stats.reads_raw = nreads/4
                created += createdinc
                linked += linkedinc
                appended += appendinc

        ## print if data were linked
        print("{} new Samples created in `{}`.".format(created, self.name))
        if linked:
            print("{} fastq files linked to {} new Samples.".\
                  format(linked, len(self.samples)))
            self.dirs.fastqs = os.path.realpath(os.path.dirname(path))
        if appended:
            print("{} fastq files appended to {} existing Samples.".\
                  format(appended, len(self.samples)))



    def link_barcodes(self):
        """ 
        Saves Sample barcodes in a dictionary as [Assembly].barcodes. Barcodes
        are parsed from the file `barcodes_path`.

        Note
        ----
        [Assembly].link_barcodes() is run automatically if set_params() is used
        to change barcodes_path.

        """
        ## in case fuzzy selected
        try: 
            barcodefile = glob.glob(self.paramsdict["barcodes_path"])[0]
        except IndexError: 
            print("Barcodes file not found:", self.paramsdict["barcodes_path"])

        ## parse barcodefile
        try:
            bdf = pd.read_csv(barcodefile, header=None, delim_whitespace=1)
            bdf = bdf.dropna()
            ## make sure upper case
            bdf[1] = bdf[1].str.upper()
            ## set attribute on Assembly object
            self.barcodes = dict(zip(bdf[0], bdf[1]))
        except ValueError:
            LOGGER.warn("Barcodes file not recognized.")



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
                sys.stdout.write("  {:<4}{:<28}{:<45}\n".format(index+1,
                           key, value))
        else:
            try:
                if int(param):
                    #sys.stdout.write(self.paramsdict.values()[int(param)-1])
                    return self.paramsdict.values()[int(param)-1]
            except (ValueError, TypeError, NameError, IndexError):
                return 'key not recognized'




    def set_params(self, param, newvalue):
        """ Set a parameter to a new value. Raises error if newvalue 
        is wrong type.

        Note
        ----
        Use [Assembly].get_params() to see the parameter values currently
        linked to the Assembly object.

        Parameters
        ----------
        param : int or str
            The index (e.g., 1) or string name (e.g., "working_directory")
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
        [Assembly].set_params('working_directory', 'new_directory')

        ## param 6 must be a tuple or str, if str it is converted to a tuple
        ## with the second entry empty.
        [Assembly].set_params(6, 'TGCAG')
        [Assembly].set_params('restriction_overhang', ('CTGCAG', 'CCGG')                            

        ## param 13 can be an int or a float:
        [Assembly].set_params(13, 4)
        [Assembly].set_params('max_shared_H_locus', 0.25)
            
        """

        ## require parameter recognition
        if not ((param in range(50)) or \
                (param in [str(i) for i in range(50)]) or \
                (param in self.paramsdict.keys())):
            raise IPyradParamsError("Parameter key not recognized: {}"\
                                    .format(param))

        ## make string
        param = str(param)

        ## get index if param is keyword arg (this index is not zero based!)
        if len(param) < 3:
            param = self.paramsdict.keys()[int(param)-1]

        ## run assertions on new param 
        try:
            self = paramschecker(self, param, newvalue)
        except Exception as inst:
            #print("\nError:", inst, "\n")
            raise IPyradParamsError(inst)


    def write_params(self, outfile=None, force=False):
        """ Write out the parameters of this assembly to a file properly
        formatted as input for `ipyrad -p <params.txt>`. A good and
        simple way to share/archive parameter settings for assemblies.
        This is also the function that's used by __main__ to 
        generate default params.txt files for `ipyrad -n`
        """
        if outfile is None:
            outfile = os.path.join(self.paramsdict["working_directory"],
                                self.name+"-params.txt")

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                LOGGER.error("Assembly.write_params() attempting to write"\
                            + " a file that already exists - {}".format(outfile))
                LOGGER.error("Use write_params(force=True) to override")
                raise IPyradError("File exists: {}. \nUse force=True to overwrite.".format(outfile))

        with open(outfile, 'w') as paramsfile:

            ## Write the header. Format to 80 columns
            header = "------ ipyrad params file (v.{})".format(ip.__version__)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number (paramsinfo is 1-based,
            ## so we have to increment the index we think it is. Also get the short 
            ## description from paramsinfo. Make it look pretty, pad nicely 
            ## if at all possible.
            for key, val in self.paramsdict.iteritems():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)
                padding = (" "*(30-len(paramvalue)))
                paramindex = " ## [{}] ".format(self.paramsdict.keys().index(key) + 1)
                description = paraminfo(self.paramsdict.keys().index(key) + 1, short=True)
                paramsfile.write("\n" + paramvalue + padding + paramindex + description)


    def copy(self, newname):
        """ Returns a copy of the Assembly object. Does not allow Assembly 
        object names to be replicated in namespace or path. """
        ## is there a better way to ask if it already exists?
        if (newname == self.name) or (os.path.exists(newname+".assembly")):
            print("Assembly object named {} already exists".format(newname))
        else:
            ## create a copy of the Assembly obj
            newobj = copy.deepcopy(self)
            newobj.name = newname

            ## create copies of each Sample obj
            for sample in self.samples:
                newobj.samples[sample] = copy.deepcopy(self.samples[sample])
            return newobj



    def filetree(self):
        """ prints the project data structure. TODO: this needs work.
        prints way too much other junk if [work] is home dir. """
        startpath = self.paramsdict["working_directory"]
        for root, _, files in os.walk(startpath):
            level = root.replace(startpath, '').count(os.sep)
            indent = ' ' * 4 * (level)
            print('{}{}/'.format(indent, os.path.basename(root)))
            subindent = ' ' * 4 * (level + 1)
            for fname in files:
                print('{}{}'.format(subindent, fname))



    def _save(self):
        """ Pickle the Assembly object. Could be used for checkpointing before
        and after assembly steps. Currently it is called after assembly steps.
        """
        dillout = open(os.path.join(
                          self.paramsdict["working_directory"],
                          self.name+".assembly"), "wb")
        dill.dump(self, dillout)
        dillout.close()



    def _launch(self, inittries):
        """ launch ipyclient.
        launch within try statement in case engines aren't ready yet
        and try 30 one second sleep/wait cycles before giving up on engines
        """
        tries = inittries
        while tries:
            try:
                ## launches ipcluster with arguments if present in self
                clusterargs = [self._ipclusterid, self._ipprofile]
                argnames = ["cluster_id", "profile"]
                args = {key:value for key, value in zip(argnames, clusterargs)}
                ipyclient = ipp.Client(**args)
                if tries > 1:
                    LOGGER.info('try %s: starting controller', tries)
                ## make sure all engines are connected
                try:
                    assert ipyclient.ids                    
                    if tries != inittries:                        
                        ## get initial number of ids
                        ## ugly hack to find all engines while they're spinng up
                        initid = ipyclient.ids
                        if len(initid) > 10:
                            LOGGER.warn("waiting 3 seconds to find Engines")
                            time.sleep(3)
                        else:
                            time.sleep(1)                                            
                        try:
                            ## make sure more engines aren't found
                            assert len(ipyclient.ids) == len(initid)
                            LOGGER.warn('OK! Connected to (%s) engines', 
                                        len(ipyclient.ids))
                            return ipyclient

                        except AssertionError as _: 
                            LOGGER.warn('finding engines (%s, %s)', 
                                         len(initid), len(ipyclient.ids))
                            raise
                    else:
                        LOGGER.debug('OK! Connected to (%s) engines', 
                                    len(ipyclient.ids))
                        return ipyclient

                except AssertionError as _: 
                    LOGGER.debug('connected to %s engines', len(ipyclient.ids))
                    raise
            except (IOError, ipp.NoEnginesRegistered, AssertionError) as _:
                time.sleep(1)
                tries -= 1
        raise ipp.NoEnginesRegistered



    def _clientwrapper(self, stepfunc, args, nwait):
        """ wraps a call with error messages for when ipyparallel fails"""
        try:
            ipyclient = self._launch(nwait)
            if ipyclient.ids:
                ## append client to args and call stepfunc
                args.append(ipyclient)
                stepfunc(*args)

        except (ipp.TimeoutError, ipp.NoEnginesRegistered):
            ## maybe different messages depending on whether it is CLI or API
            inst = """
        Check to that ipcluster is running. When using the API you must start
        ipcluster outside of IPython/Jupyter to launch parallel engines using
        either `ipcluster start`, or in the Clusters tab in a Jupyter notebook.
        (See Docs)
            """
            ## raise right away since there is no ipyclient to close
            raise IPyradError(inst)

        ## except user or system interrupt
        except KeyboardInterrupt as inst:
            ipyclient.abort()
            ipyclient.shutdown()
            ipyclient = self._launch(nwait)
            logging.error("assembly interrupted by user.")
            raise IPyradError("Keyboard Interrupt")

        except SystemExit as inst:
            logging.error("assembly interrupted by sys.exit.")
            raise IPyradError("SystemExit Interrupt")

        except AssertionError as inst:
            logging.error("Assertion: %s", inst)
            raise IPyradError(inst)

        except IPyradWarningExit as inst:
            print(inst)

        #except Exception as inst:
        #    print("Exception:", inst)
        #    raise inst

        ## close client when done or interrupted
        finally:
            ## can't close client if it was never open
            try:
                ipyclient.close()
                ## pickle the data obj
                self._save()                
            except UnboundLocalError:
                pass




    def _step1func(self, force, preview, ipyclient):
        """ testing"""
        msg1 = "  step1: Demultiplexing raw reads."
        msg2 = "  step1: Linking demultiplexed fastq data to Samples."

        ## if Samples already exist then no demultiplexing
        if self.samples:
            if not force:
                print("Skipping step1: {} ".format(len(self.samples)) \
                    +"Samples already found in `{}` ".format(self.name)\
                    +"(can overwrite with force option).")
            else:
                if self._headers:
                    print(msg1)
                assemble.demultiplex.run(self, preview, ipyclient)
                self._stamp("s1_demultiplexing:")

        ## Creating new Samples
        else:
            ## first check if demultiplexed files exist in path
            if os.path.exists(self.paramsdict["sorted_fastq_path"]):
                self.link_fastqs(path=os.path.join(
                    self.paramsdict["sorted_fastq_path"], "*"))
                if self._headers:
                    print(msg2, "\n  linking files from {}".\
                          format(self.paramsdict["sorted_fastq_path"]))
            ## otherwise do the demultiplexing
            else:
                if self._headers:
                    print(msg1)                        
                assemble.demultiplex.run(self, preview, ipyclient)
                self._stamp("s1_demultiplexing:")



    def _step2func(self, samples, nreplace, force, preview, ipyclient):
        if self._headers:
            print("  Step2: Filtering reads ")

        ## If no samples in this assembly then it means you skipped step1,
        ## so attempt to link existing demultiplexed fastq files
        if not self.samples.keys():
            self.link_fastqs()

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## pass samples to rawedit
        assemble.rawedit.run(self, samples, nreplace, force, preview, ipyclient)



    def _step3func(self, samples, noreverse, force, preview, ipyclient):
        """ hidden wrapped function to start step 3 """
        ## print headers
        if self._headers:
            print("  Step3: Clustering/Mapping reads")

        ## Require reference seq for reference-based methods
        if self.paramsdict['assembly_method'] != "denovo":
            if not self.paramsdict['reference_sequence']:
                raise IPyradError("Reference or hybrid assembly requires a "+\
                                  "value for reference_sequence_path paramter.")

            ## index the reference sequence
            index_reference_sequence(self)

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## skip if all are finished
        if not force:
            if all([i.stats.state >= 3 for i in samples]):
                print("  Skipping: All {} ".format(len(samples))\
                     +"selected Samples already clustered")
            else:
                assemble.cluster_within.run(self, samples, noreverse, 
                                            force, preview, ipyclient)
        else:
            assemble.cluster_within.run(self, samples, noreverse, 
                                        force, preview, ipyclient)



    def _step4func(self, samples, subsample, force, ipyclient):
        """ step 4: Joint estimation of error rate and heterozygosity. 
        If you want to overwrite data for a file, first set its state to 3:
        data.samples['sample'].stats['state'] = 3 """
        if self._headers:
            print("  Step4: Joint estimation of error rate and heterozygosity")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## if keys are in list
        #if any([isinstance(i, str) for i in samples]):
            ## make into a subsampled sample dict
        #    subsamples = {i: self.samples[i] for i in samples}

        ## send to function
        assemble.jointestimate.run(self, samples, subsample, force, ipyclient)


    def _step5func(self, samples, force, ipyclient):
        """ hidden wrapped function to start step 5 """
        ## print header
        if self._headers:
            print("  Step5: Consensus base calling ")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## pass samples to rawedit
        assemble.consens_se.run(self, samples, force, ipyclient)



    def _step6func(self, samples, noreverse, force, randomseed, ipyclient):
        """ hidden function to start Step 6"""

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        if self._headers:
            print("  Step 6: clustering across {} samples at {} similarity".\
            format(len(samples), self.paramsdict["clust_threshold"]))

        ## attach filename for all reads database
        self.database = os.path.join(self.dirs.consens, self.name+".hdf5")

        ## check for existing and force
        if not force:
            if os.path.exists(self.database):
                print("  Skipping step6: Clust file already exists:"\
                      +"{}\n".format(self.database))
            else:
                assemble.cluster_across.run(self, samples, noreverse,
                                            force, randomseed, ipyclient)
        else:
            assemble.cluster_across.run(self, samples, noreverse,
                                        force, randomseed, ipyclient)


    def _step7func(self, samples, force, ipyclient):
        """ Step 7: Filter and write output files """

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        if not force:
            try:
                if os.path.exists(self.dirs.outfiles):
                    print(\
    "  Step 7: Cowardly refusing to overwrite existing output directory {}".\
                    format(self.dirs.outfiles))
                    print(\
    "  Step 7: rerun with `force=True` to overwrite")
                    sys.exit()
            except AttributeError as _:
                ## If not force and directory doesn't exist then nbd.
                pass

        assemble.write_outfiles.run(self, samples, force, ipyclient)


    def step1(self, force=False, preview=False):
        """ test """
        self._clientwrapper(self._step1func, [force, preview], 10)

    def step2(self, samples=None, nreplace=True, force=False, preview=False):
        """ 
        Edit/Filter raw demultiplexed reads based on read quality scores and the
        presence of Illumina adapter sequences. 

        The following parameters are used in this step:
            - datatype
            - phred_Qscore_offset
            - max_low_qual_bases
            - filter_adapters
            - filter_min_trim_len
            - edit_cutsites
            - restriction_overhang
            ...

        Parameters
        ----------
        samples : list or str
            By default all Samples linked to an Assembly object are run. If a 
            subset of Sampled is entered as a list then only those Samples will 
            be run. 

        nreplace : bool
            If True (default) low quality base calls (Q < 20 given the 
            `phred_Qscore_offset`) are converted to Ns. If False, low quality 
            bases are not converted, but simply counted. Reads with > 
            `max_low_qual_bases` are excluded. 

        force : bool
            If force=True existing files are overwritten, otherwise Samples in 
            state 2 will return a warning that the Sample has already been run. 

        preview : bool
            ...
        """
        self._clientwrapper(self._step2func, 
                           [samples, nreplace, force, preview], 10)

    def step3(self, samples=None, noreverse=False, force=False, preview=False):
        """ 
        Demultiplex reads and then cluster/map denovo or with a reference 
        sequence file. 

        The following parameters are used in this step:
            - datatype
            - assembly_method
            - clust_threshold
            - 
            ...

        Parameters
        ----------
        samples : list or str
            By default all Samples linked to an Assembly object are run. If a 
            subset of Sampled is entered as a list then only those Samples will 
            be run. 

        noreverse : bool
            ...

        force : bool
            ...

        preview : bool
            ...
        """
        self._clientwrapper(self._step3func, 
                           [samples, noreverse, force, preview], 10)


    def step4(self, samples=None, subsample=None, force=False):
        """ test """
        self._clientwrapper(self._step4func, [samples, subsample, force], 10)



    def step5(self, samples=None, force=False):
        """ 
        Consensus base calling and filtering from within-sample clusters. 
        Samples must be in state 3 or 4 (passed step3 and/or step4).

        The following parameters are used in this step: 
            - max_Ns_consens
            - max_Hs_consens
            - maxdepth
            - mindepth_statistical
            - mindepth_majrule
            - ploidy

        If you want to overwrite data for a file, first set its state to 
        3 or 4. e.g., data.samples['sample'].stats['state'] = 3 

        Parameters
        ----------
        samples : list or str
            By default all Samples linked to an Assembly object are run. 
            If a subset of Samples is entered as a list then only those samples
            will be executed. 

        force : bool
            Force to overwrite existing files. By default files will not be 
            overwritten unless force=True. 
        """

        self._clientwrapper(self._step5func, [samples, force], 2)



    def step6(self, samples=None, noreverse=False, force=False, randomseed=123):
        """ 
        Cluster consensus reads across samples and align with muscle. 

        Parameters
        ----------
        samples : list or str
            By default all Samples linked to an Assembly object are clustered. 
            If a subset of Samples is entered as a list then only those samples
            will be clustered. It is recommended to create .copy() Assembly 
            objects if step6 is performed on different subsets of Samples. 

        noreverse : bool
            Reverse complement clustering is performed on gbs and pairgbs data
            types by default. If noreverse=True then reverse complement 
            clustering will not be performed. This can improve clustering speed.

        force : bool
            Force to overwrite existing files. By default files will not be 
            overwritten unless force=True. 

        randomseed : int
            Consensus reads are sorted by length and then randomized within 
            size classes prior to clustering. The order of sequences in this 
            list can (probably minimally) affect their clustering. The default
            randomseed is 123. Thus, unless it is changed results should be 
            reproducible. 
        """
        self._clientwrapper(self._step6func, [samples, noreverse, force,
                                              randomseed], 10)


    def step7(self, samples=None, force=False):
        """ 
        Create output files in a variety of formats. 

        Parameters
        ----------
        samples : list or str
            ...

        force : bool
            ...

        """
        self._clientwrapper(self._step7func, [samples, force], 10)




    def run(self, steps=0, force=False, preview=False):
        """ Select steps of an analysis. If no steps are entered then all
        steps are run. Enter steps as a string, e.g., "1", "123", "12345" """
        if not steps:
            steps = list("123457")
        else:
            steps = list(steps)

        if '1' in steps:
            self.step1(preview=preview)
        if '2' in steps:
            self.step2(force=force, preview=preview)
        if '3' in steps:
            self.step3(force=force, preview=preview)
        if '4' in steps:
            self.step4(force=force)
        if '5' in steps:
            self.step5(force=force)
        if '6' in steps:
            self.step6()            
        if '7' in steps:
            self.step7()


def _get_samples(self, samples):
    """ Internal function. Prelude for each step() to read in perhaps
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



def _name_from_file(fname):
    """ internal func: get the sample name from any pyrad file """
    file_extensions = [".gz", ".fastq", ".fq", ".fasta", 
                       ".clustS", ".consens"]
    base, ext = os.path.splitext(os.path.basename(fname))
    ## remove read number from name
    base = base.replace("_R1_.", ".")\
               .replace("_R1_", "_")\
               .replace("_R1.", ".")
    ## remove extensions
    while ext in file_extensions:
        base, ext = os.path.splitext(base)
    return base



def expander(namepath):
    """ expand ./ ~ and ../ designators in location names """        
    if "~" in namepath:
        namepath = namepath.replace("~", os.path.expanduser("~"))
    if "../" in namepath:
        _, post = namepath.split("../")
        namepath = os.path.abspath(
                    os.path.join(
                        os.path.dirname(""), '..', post))

    elif "./" in namepath:
        _, post = namepath.split("./")
        namepath = os.path.abspath("")+"/"+post
    return namepath



def cmd_exists(cmd):
    """ check if dependency program is there """
    return subprocess.call("type " + cmd,
                           shell=True, 
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) == 0



def getbins():
    """ gets the right version of vsearch, muscle, and smalt
    depending on linux vs osx """

    # Return error if system is 32-bit arch.
    # This is straight from the python docs:
    # https://docs.python.org/2/library/platform.html#cross-platform
    if not sys.maxsize > 2**32:
        sys.exit("iPyrad requires 64bit architecture") 

    ## get platform mac or linux
    _platform = sys.platform

    ## get current location
    path = os.path.abspath(os.path.dirname(__file__))

    ## find bin directory
    ipyrad_path = os.path.dirname(os.path.dirname(path))
    bin_path = os.path.join(ipyrad_path, "bin")

    ## get the correct binaries 
    if 'linux' in _platform:
        vsearch = os.path.join(
                       os.path.abspath(bin_path),
                       "vsearch-linux-x86_64")
        muscle = os.path.join(
                       os.path.abspath(bin_path),
                       "muscle-linux-x86_64")
        smalt = os.path.join(
                       os.path.abspath(bin_path),
                       "smalt-linux-x86_64")
        samtools = os.path.join(
                       os.path.abspath(bin_path),
                       "samtools-linux-x86_64")
        bedtools = os.path.join(
                       os.path.abspath(bin_path),
                       "bedtools-linux-x86_64")
    else:
        vsearch = os.path.join(
                       os.path.abspath(bin_path),
                       "vsearch-osx-x86_64")
        muscle = os.path.join(
                       os.path.abspath(bin_path),
                       "muscle-osx-x86_64")
        smalt = os.path.join(
                       os.path.abspath(bin_path),
                       "smalt-osx-x86_64")
        samtools = os.path.join(
                       os.path.abspath(bin_path),
                       "samtools-osx-x86_64")
        bedtools = os.path.join(
                       os.path.abspath(bin_path),
                       "bedtools-osx-x86_64")

    # Test for existence of binaries
    assert cmd_exists(muscle), "muscle not found here: "+muscle
    assert cmd_exists(vsearch), "vsearch not found here: "+vsearch
    assert cmd_exists(smalt), "smalt not found here: "+smalt
    assert cmd_exists(samtools), "samtools not found here: "+samtools
    assert cmd_exists(bedtools), "bedtools not found here: "+bedtools
    return vsearch, muscle, smalt, samtools, bedtools



def merge(name, assemblies):
    """ Creates and returns a new Assembly object in which 
    samples from two or more Assembly objects with matching names
    are 'merged'. Merging does not affect the actual files written
    on disk, but rather creates new Samples that are linked to 
    multiple data files, and with stats summed. """

    ## checks
    assemblies = list(assemblies)

    ## create new Assembly
    merged = assemblies[0].copy(name)

    ## get all sample names from all Assemblies
    allsamples = set(merged.samples.keys())
    for iterass in assemblies[1:]:
        allsamples.update(set(iterass.samples.keys()))

    ## iterate over assembly objects, skip first already copied
    for iterass in assemblies[1:]:
        ## iterate over stats, skip 'state'
        for stat in merged.stats.keys()[1:]:
            ## iterate over allsamples, add if not in merged
            for sample in iterass.samples:
                if sample not in merged.samples:
                    merged.samples[sample] = iterass.samples[sample]
                ## merge stats
                merged.samples[sample].stats[stat] += \
                                  iterass.samples[sample].stats[stat]
                ## merge file references
                for filetype in ["fastq", "edits", "clusters", "consens"]:
                    merged.samples[sample].files[filetype].append(
                                  iterass.samples[sample].files[filetype])

    ## return the new Assembly object
    return merged



def bufcountlines(filename, gzipped):
    """ fast line counter. Used to quickly sum number of input reads 
    when running link_fastqs to append files. """
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


def tuplecheck(newvalue, dtype=str):
    """ Takes a string argument and returns value as a tuple. 
    Needed for paramfile conversion from CLI to set_params args """
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
    if param == 'working_directory':
        expandpath = expander(newvalue)
        if not expandpath.startswith("/"): 
            if os.path.exists(expandpath):
                expandpath = "./"+expandpath
                expandpath = expander(expandpath)
        self._stamp("[{}] set to {}".format(param, newvalue))
        self.paramsdict["working_directory"] = expandpath
        self.dirs["working"] = expandpath

    elif param == 'raw_fastq_path':
        fullrawpath = expander(newvalue)
        if os.path.isdir(fullrawpath):
            fullrawpath = os.path.join(fullrawpath, "*.gz")
        self.paramsdict['raw_fastq_path'] = fullrawpath
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'barcodes_path':
        #assert type(newvalue) is StringType, "arg must be a string"
        fullbarpath = expander(newvalue)
        if glob.glob(fullbarpath):
            self.paramsdict['barcodes_path'] = fullbarpath
            self.link_barcodes()
            self._stamp("[{}] set to {}".format(param, newvalue))
        elif not fullbarpath:
            self.paramsdict['barcodes_path'] = fullbarpath                
            self._stamp("[{}] set to {}".format(param, newvalue))
        else:
            print(
        "Warning: barcodes file not found. This must be an absolute \n"\
       +"path (/home/wat/ipyrad/data/data_barcodes.txt) or relative to the \n"\
       +"directory where you're running ipyrad (./data/data_barcodes.txt).\n"\
       +"You entered: %s\n" % fullbarpath)

    elif param == 'sorted_fastq_path':
        assert isinstance(newvalue, str), \
        "sorted_fastq_path must be a string, e.g., /home/data/fastqs/*"
        newvalue = expander(newvalue)
        if os.path.isdir(newvalue):
            newvalue = os.path.join(newvalue, "*.gz")
        self.paramsdict['sorted_fastq_path'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'assembly_method':
        assert newvalue in ["denovo", "reference_only", "hybrid", "denovo_only"], \
            "The `assembly_method` parameter must be one of the following: "+\
            "denovo, reference, hybrid, or denovo_only. You entered: %s." % newvalue
        self.paramsdict['assembly_method'] = newvalue            
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'reference_sequence':
        fullrawpath = expander(newvalue)
        if not os.path.isfile(fullrawpath):
            LOGGER.info("reference sequence file not found.")
            print(\
        "Warning: reference sequence file not found. This must be an\n"\
       +"absolute path (/home/wat/ipyrad/data/reference.gz) or relative to \n"\
       +"the directory where you're running ipyrad (./data/reference.gz).\n"\
       +"You entered: %s\n" % fullrawpath)
        self.paramsdict['reference_sequence'] = fullrawpath
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'datatype':
        ## list of allowed datatypes
        datatypes = ['rad', 'gbs', 'ddrad', 'pairddrad',
                     'pairgbs', 'merged', '2brad']
        ## raise error if something else
        if str(newvalue) not in datatypes:
            sys.exit("error: datatype {} not recognized, must be one of: ".format( newvalue ), datatypes)
        else:
            self.paramsdict['datatype'] = str(newvalue)
            self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'restriction_overhang':
        newvalue = tuplecheck(newvalue, str)                        
        assert isinstance(newvalue, tuple), \
        "cut site must be a tuple, e.g., (TGCAG, '') or (TGCAG, CCGG)"
        if len(newvalue) == 1:
            newvalue = (newvalue, "")
        assert len(newvalue) == 2, \
        "must enter 1 or 2 cut sites, e.g., (TGCAG, '') or (TGCAG, CCGG)"
        self.paramsdict['restriction_overhang'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_low_qual_bases':
        assert isinstance(int(newvalue), int), \
            "max_low_qual_bases must be an integer."        
        self.paramsdict['max_low_qual_bases'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'phred_Qscore_offset':
        assert isinstance(int(newvalue), int), \
            "phred_Qscore_offset must be an integer."
        self.paramsdict['phred_Qscore_offset'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'mindepth_statistical':
        assert isinstance(int(newvalue), int), \
            "mindepth_statistical must be an integer."
        ## do not allow values below 5
        if int(newvalue) < 5:
            print(\
        "error: mindepth statistical cannot be set < 5. Use mindepth_majrule.")
        ## do not allow majrule to be > statistical
        elif int(newvalue) < self.paramsdict["mindepth_majrule"]:
            print(\
        "error: mindepth statistical cannot be less than mindepth_majrule")                
        else:
            self.paramsdict['mindepth_statistical'] = int(newvalue)
            ## TODO: calculate new clusters_hidepth if passed step3
            self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'mindepth_majrule':
        assert isinstance(int(newvalue), int), \
            "mindepth_majrule must be an integer."
        if int(newvalue) > self.paramsdict["mindepth_statistical"]:
            print(\
        "error: mindepth_majrule cannot be > mindepth_statistical")
        else:
            ## TODO: calculate new clusters_hidepth if passed step3
            self.paramsdict['mindepth_majrule'] = int(newvalue)
            self._stamp("[{}] set to {}".format(param, newvalue))

    ## TODO: not yet implemented
    elif param == 'maxdepth':
        self.paramsdict['maxdepth'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'clust_threshold':
        self.paramsdict['clust_threshold'] = float(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_barcode_mismatch':
        self.paramsdict['max_barcode_mismatch'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'filter_adapters':
        self.paramsdict['filter_adapters'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'filter_min_trim_len':
        self.paramsdict['filter_min_trim_len'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_alleles_consens':
        self.paramsdict['max_alleles_consens'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_Ns_consens':
        newvalue = tuplecheck(newvalue, int)                        
        assert isinstance(newvalue, tuple), \
        "max_Ns_consens should be a tuple e.g., (8, 8)"
        self.paramsdict['max_Ns_consens'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_Hs_consens':
        newvalue = tuplecheck(newvalue, int)                        
        assert isinstance(newvalue, tuple), \
        "max_Hs_consens should be a tuple e.g., (5, 5)"
        self.paramsdict['max_Hs_consens'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'min_samples_locus':
        self.paramsdict['min_samples_locus'] = int(newvalue)
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_shared_H_locus':
        self.paramsdict['max_shared_H_locus'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_SNPs_locus':
        newvalue = tuplecheck(newvalue, int)                        
        assert isinstance(newvalue, tuple), \
        "max_SNPs_locus should be a tuple e.g., (20, 20)"
        self.paramsdict['max_SNPs_locus'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'max_Indels_locus':
        newvalue = tuplecheck(newvalue, int)            
        assert isinstance(newvalue, tuple), \
        "max_Indels_locus should be a tuple e.g., (5, 100)" 
        self.paramsdict['max_Indels_locus'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))
 
    elif param == 'edit_cutsites':
        ## Check if edit_cutsites is int or string values.
        ## Try int, if it fails the fall back to str
        try:
            newvalue = tuplecheck(newvalue, int)
        except ValueError as e:
            print("edit_cutsites value error - {}".format(e))
            newvalue = tuplecheck(newvalue)
            assert isinstance(newvalue, tuple), \
                "edit_cutsites should be a tuple e.g., (0, 5), you entered {}"\
                .format(newvalue)

        ## If edit_cutsites params are ints, then cast the tuple values
        ## to ints. If they aren't ints then just leave them as strings.
        #try:
        #    newvalue = (int(newvalue[0]), int(newvalue[1]))
        #except ValueError as e:
        #    LOGGER.info("edit_cutsites values are strings - {} {}".format(\
        #                newvalue[0], newvalue[1]))

        self.paramsdict['edit_cutsites'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'trim_overhang':
        newvalue = tuplecheck(newvalue, int)
        assert isinstance(newvalue, tuple), \
        "trim_overhang should be a tuple e.g., (1, 2, 2, 1)"
        self.paramsdict['trim_overhang'] = newvalue
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'output_formats':
        ## Get all allowed file types from assembly.write_outfiles
        output_formats = assemble.write_outfiles.OUTPUT_FORMATS

        ## If wildcard, then just do them all
        if "*" in newvalue:
            requested_formats = output_formats
        else:
            ## output_formats should be comma separated, with optional spaces
            requested_formats = newvalue.replace(" ", "").split(',')

            ## Exit if requested formats are bad
            ## Only test here if no wildcard present
            for f in requested_formats:
                if f not in output_formats:
                    sys.exit("error: File format [ {} ] not recognized, must be one of: {}".format(f , output_formats))
        
        self.paramsdict['output_formats'] = requested_formats
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'pop_assign_file':
        fullpoppath = expander(newvalue)

        if not os.path.isfile(fullpoppath):
            LOGGER.warn("Population assignment file not found.")
            sys.exit(\
        "Warning: Population assignment file not found. This must be an\n"\
       +"absolute path (/home/wat/ipyrad/data/my_popfile.txt) or relative to \n"\
       +"the directory where you're running ipyrad (./data/my_popfile.txt).\n"\
       +"You entered: %s\n" % fullpoppath)

        self.paramsdict['pop_assign_file'] = fullpoppath
        self.link_populations( )
        self._stamp("[{}] set to {}".format(param,fullpoppath))

    elif param == 'excludes':
        excluded_individuals = newvalue.replace(" ", "").split(',')

        ## Test if the individuals requested for exclusion actually
        ## exist in sample list? I hate implicit failure, but it could
        ## be tricky to handle the case where people set excludes before
        ## they run step1?
        ## TODO: Maybe do this, maybe not.

        self.paramsdict['excludes'] = excluded_individuals
        self._stamp("[{}] set to {}".format(param, newvalue))

    elif param == 'outgroups':
        outgroup_individuals = newvalue.replace(" ", "").split(',')

        ## Test if the outgroup individuals actually
        ## exist in sample list? I hate implicit failure, but it could
        ## be tricky to handle the case where people set excludes before
        ## they run step1?
        ## TODO: Maybe do this, maybe not.

        self.paramsdict['excludes'] = outgroup_individuals
        self._stamp("[{}] set to {}".format(param, newvalue))

    return self




if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.get_params()
    print(DATA.log)
