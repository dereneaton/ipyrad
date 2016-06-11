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
import copy
import h5py
import string
import cStringIO
import itertools
import numpy as np
import pandas as pd
import ipyparallel as ipp
import ipyrad as ip

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
            msg = """\n
    Assembly name _must_ be set. This is the first parameter in the params.txt 
    file, and will be used as a prefix for output files. It should be a short 
    string with no special characters, i.e., not a path (no \"/\" characters). 
    If you need a suggestion, name it after the organism you're working on.\n"""
            raise IPyradParamsError(msg)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        invalid_chars = string.punctuation.replace("_", "")\
                                          .replace("-", "")+ " "
        if any(char in invalid_chars for char in name):
            msg = """\n
    No spaces or special characters are allowed in the assembly name. A good 
    practice is to replace spaces with underscores '_'. An example of a good 
    assembly_name is: white_crowned_sparrows. Here's what you put:
    {}""".format(name)
            raise IPyradParamsError(msg)

        self.name = name
        if not quiet:
            print("  New Assembly: {}".format(self.name))

        ## Store assembly version #
        self._version = ip.__version__ 

        ## stores ipcluster launch info
        self._ipcluster = {}
        for ipkey in ["id", "profile", "engines", "cores"]:
            self._ipcluster[ipkey] = None

        ## print headers
        self._headers = 0

        ## analysis save objects
        self.svd = ObjDict({})
        self.dstat = ObjDict({})

        ## all available cpus
        self.cpus = detect_cpus()

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
                       ("clust_threshold", .85),
                       ("max_barcode_mismatch", 1),
                       ("filter_adapters", 0), 
                       ("filter_min_trim_len", 35), 
                       ("max_alleles_consens", 2), 
                       ("max_Ns_consens", (5, 5)), 
                       ("max_Hs_consens", (8, 8)), 
                       ("min_samples_locus", 4), 
                       ("max_SNPs_locus", (50, 50)), 
                       ("max_Indels_locus", (8, 8)), 
                       ("max_shared_Hs_locus", .25), 
                       ("edit_cutsites", (0, 0)),
                       ("trim_overhang", (4, 4, 4, 4)),                        
                       ("output_formats", "*"),
                       ("pop_assign_file", ""),
        ])
        #               ("excludes", ""),
        #               ("outgroups", ""),

        ## Store data directories for this Assembly. Init with default project
        self.dirs = ObjDict({"project": self.paramsdict["project_dir"],
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
                        ("max_fragment_length", 150),
                        ("max_inner_mate_distance", 60),
                        ("preview_step1", 4000000),
                        ("preview_step2", 100000),
                        ("output_loci_name_buffer", 5),
                        ("query_cov", None),
                        ("smalt_index_wordlen", 8)
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



    def link_fastqs(self, path=None, merged=False, force=False, append=False,
                    splitnames="_", fields=None):
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

        print("""\
    Linking to demultiplexed fastq files in:
      {}""".format(path))

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
            inst = "No files found in `sorted_fastq_path`: {}".\
                   format(self.paramsdict["sorted_fastq_path"])
            ## check for simple naming error
            if any(["_R1." in i for i in fastqs]):
                inst += "\nNames should contain _R1_, not _R1."
            raise IPyradError(inst)

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
            fastqs = [(i, j) for i, j in zip(r1_files, r2_files)]

        ## data are not paired, create empty tuple pair
        else:
            ## print warning if _R2_ is in names when not paired
            if any(["_R2_" in i for i in fastqs]):
                print("""
        Warning: '_R2_' was detected in a file name, which suggests the data
        may be paired-end. If so, you should set the Assembly parameter
        `datatype` to a paired type (e.g., pairddrad or pairgbs) and run
        with the force argument to re-link fastq data.
        """)
            fastqs = [(i, "") for i in fastqs]

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
            sname = _name_from_file(fastqtuple[0], splitnames, fields)

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
                    #if fastqtuple not in self.samples[sname].files.fastqs:
                    self.samples[sname].files.fastqs.append(fastqtuple)
                    appendinc += 1
                    #else:
                    #    print("    The files {} are already in Sample {}, "\
                    #          .format(fastqtuple, sname) \
                    #          +"cannot append duplicate files to a Sample.\n")
                elif force:
                    ## overwrite/create new Sample
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
        print("    {} new Samples created in `{}`.".format(created, self.name))
        if linked:
            ## double for paired data
            if 'pair' in self.paramsdict["datatype"]:
                linked = linked*2
            print("    {} fastq files linked to {} new Samples.".\
                  format(linked, len(self.samples)))
            ## save the location where these files are located
            self.dirs.fastqs = os.path.realpath(os.path.dirname(path))

        if appended:
            print("    {} fastq files appended to {} existing Samples.".\
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
            raise IPyradWarningExit("""
    Error: Barcodes file not found: {}
    """.format(self.paramsdict["barcodes_path"]))

        ## parse barcodefile
        try:
            bdf = pd.read_csv(barcodefile, header=None, delim_whitespace=1)
            bdf = bdf.dropna()
            ## make sure upper case
            bdf[1] = bdf[1].str.upper()

            ## 3rad/seqcap use multiplexed barcodes
            ## We'll concatenate them with a plus and split them later
            if "3rad" in self.paramsdict["datatype"]:
                bdf[2] = bdf[2].str.upper()
                self.barcodes = dict(zip(bdf[0], bdf[1] + "+" + bdf[2]))
            else:
                ## set attribute on Assembly object
                self.barcodes = dict(zip(bdf[0], bdf[1]))

        except ValueError:
            msg = "Barcodes file not recognized."
            LOGGER.warn(msg)
            raise IPyradError(msg)



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
        """ Set a parameter to a new value. Raises error if newvalue 
        is wrong type.

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
            if ip.__interactive__:
                raise IPyradError("""
    Error setting parameter {}: 
    Exception: {}. 
    You entered: {}""".format(param, inst, newvalue))
            else:
                raise IPyradWarningExit("""
    Error setting parameter {}: {}. 
    You entered: {}""".format(param, inst, newvalue))



    def write_params(self, outfile=None, force=False):
        """ Write out the parameters of this assembly to a file properly
        formatted as input for `ipyrad -p <params.txt>`. A good and
        simple way to share/archive parameter settings for assemblies.
        This is also the function that's used by __main__ to 
        generate default params.txt files for `ipyrad -n`
        """
        if outfile is None:
            #usedir = self.paramsdict["project_dir"]
            #if not os.path.exists(usedir):
            #    usedir = "./"
            #outfile = os.path.join(
            outfile = "params-"+self.name+".txt"

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise IPyradWarningExit("""
  **Aborting** File exists: {}
  Use force argument to overwrite.
    """.format(outfile))

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



    def branch(self, newname, subsamples=[]):
        """ Returns a copy of the Assembly object. Does not allow Assembly 
        object names to be replicated in namespace or path. """
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

            ## create copies of each subsampled Sample obj
            if subsamples:
                for sname in subsamples:
                    if sname in self.samples:
                        newobj.samples[sname] = copy.deepcopy(self.samples[sname])
                    else:
                        print("  Sample name not found: {}".format(sname))
                ## reload sample dict w/o non subsamples
                newobj.samples = {name:sample for name, sample in \
                           self.samples.items() if name in subsamples}

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
        if self._headers:
            print("  Saving Assembly.\n")
        ip.save_json(self)



    def _launch2(self, nwait):
        """ 
        Creates a client for a given profile to connect to the running 
        clusters. If CLI, the cluster should have been started by __main__, 
        if API, you need to start the cluster on your own. 
        """

        #save_stdout = sys.stdout           
        try: 
            clusterargs = [self._ipcluster['id'], self._ipcluster["profile"]]
            argnames = ["cluster_id", "profile"]
            args = {key:value for key, value in zip(argnames, clusterargs)}

            ## wait for at least 1 engine to connect
            for _ in range(nwait):
                try:
                    ## using this wrap to avoid ipyparallel's ugly warnings
                    ## save orig stdout
                    save_stdout = sys.stdout 
                    save_stderr = sys.stderr
                    ## file-like obj to catch stdout
                    sys.stdout = cStringIO.StringIO()
                    sys.stderr = cStringIO.StringIO()                    
                    ## run func with stdout hidden
                    ipyclient = ipp.Client(**args)
                    ## resets stdout
                    sys.stdout = save_stdout
                    sys.stderr = save_stderr
                    break

                except IOError as inst:
                    time.sleep(0.1)

            ## check that all engines have connected            
            for _ in range(300):
                initid = len(ipyclient)
                time.sleep(0.1)
                if not ip.__interactive__:
                    ## no longer waiting for all engines, just one, and trust
                    ## that others will come. We let load_balance distribute
                    ## jobs to them.
                    if initid:
                        break
                else:
                    ## don't know how many to expect for interactive, but the
                    ## connection stays open, so just wait til no new engines
                    ## have been added for three seconds
                    time.sleep(3)
                    if initid:
                        if len(ipyclient) == initid:
                            break
                    else:
                        print("  connecting to Engines...")


        except KeyboardInterrupt as inst:
            ## ensure stdout is reset even if Exception was raised            
            sys.stdout = save_stdout
            raise inst

        except IOError as inst:
            ## ensure stdout is reset even if Exception was raised
            sys.stdout = save_stdout
            print(inst)
            raise inst

        return ipyclient



    def _clientwrapper(self, stepfunc, args, nwait):
        """ wraps a call with error messages for when ipyparallel fails"""
        ## emtpy error string
        inst = ""

        ## wrapper to ensure closure of ipyparallel
        try:
            ipyclient = ""
            ipyclient = self._launch2(nwait)
            args.append(ipyclient)
            stepfunc(*args)

        except (ipp.TimeoutError, ipp.NoEnginesRegistered) as inst:
            ## raise by ipyparallel if no connection file is found for 30 sec.
            msg = """
    No Engines found... ensure ipcluster is running (see API docs for details).
    When using the API you must start an ipyparallel instance using either 
    `ipcluster start` from a terminal, or the Clusters tab in a Jupyter notebook.
    """
            if not ip.__interactive__:
                msg = """
    There was a problem connecting to parallel engines. See Docs for advice.
            """
            ## raise right away since there is no ipyclient to close
            msg = "ipyrad error message - {}".format(inst) + "\n\n" + msg 
            raise IPyradError(msg)

        except IOError as inst:
            LOGGER.error("IOError: {}".format(inst))
            raise

        ## except user or system interrupt
        except KeyboardInterrupt as inst:
            ## abort and allow wrapper to save and close
            LOGGER.info("assembly interrupted by user.")
            print("  Keyboard Interrupt by user")
            sys.exit(2)

        except IPyradWarningExit as inst:
            ## save inst for raise error after finally statement
            LOGGER.info("IPyradWarningExit: %s", inst)
            print("  IPyradWarningExit: {}".format(inst))

        except SystemExit as inst:
            LOGGER.info("assembly interrupted by sys.exit.")
            print("  SystemExit Interrupt: {}".format(inst))

        ## An Engine Crashed. Raise a readable traceback message.
        except ipp.error.CompositeError as inst:
            ## print the trace if it's turned on, tho
            print(inst.print_traceback())

            ## find and print engine error for debugging
            for job in ipyclient.metadata:
                if ipyclient.metadata[job]['error']:
                    print(ipyclient.metadata[job]['error'])

        except IPyradError as inst:
            LOGGER.info(inst)
            print("  IPyradError: {}".format(inst))            

        except Exception as inst:
            ## Caught unhandled exception, print and reraise
            LOGGER.error(inst)
            print("  Caught unknown exception - {}".format(inst))
            raise  ## uncomment raise to get traceback


        ## close client when done or interrupted
        finally:
            try:
                ## pickle the data obj
                self.save()                
                ## can't close client if it was never open
                if ipyclient:
                    ## if CLI, stop jobs and shutdown
                    if not ip.__interactive__:
                        ipyclient.abort()                        
                    ## if API, stop jobs and clean queue
                    else:
                        ipyclient.abort()
                        #ipyclient.purge_everything()
                    ipyclient.close()
            ## if exception is close and save, print and ignore
            except Exception as inst2:
                LOGGER.error("shutdown warning: %s", inst2)

            if inst:
                IPyradWarningExit(inst)


    def _step1func(self, force, preview, ipyclient):
        """ hidden wrapped function to start step 1 """

        ## check input data files
        sfiles = self.paramsdict["sorted_fastq_path"]
        rfiles = self.paramsdict["raw_fastq_path"]

        ## do not allow both a sorted_fastq_path and a raw_fastq
        if sfiles and rfiles:
            raise IPyradWarningExit("""
    Error: Must enter a raw_fastq_path or sorted_fastq_path, but not both.""")

        ## but also require that at least one exists
        if not (sfiles or rfiles):
            raise IPyradWarningExit("""
    Error: Step 1 requires that you enter either:
        (1) a sorted_fastq_path or (2) a raw_fastq_path and barcodes_path
    """)

        ## print headers
        if self._headers:
            if sfiles:
                print("  Step1: Linking sorted fastq data to Samples")
            else:
                print("  Step1: Demultiplexing fastq data to Samples")                

        ## if Samples already exist then no demultiplexing
        if self.samples:
            if not force:
                print("""\
    Skipping: {} Samples already found in Assembly {}. 
    (can overwrite with force argument)\
    """.format(len(self.samples), self.name))
            else:
                ## overwrite existing data
                if glob.glob(sfiles):
                    if self._headers:
                        self.link_fastqs()

                ## otherwise do the demultiplexing
                else:
                    assemble.demultiplex.run(self, preview, ipyclient, force)

        ## Creating new Samples
        else:
            ## first check if demultiplexed files exist in sorted path
            if glob.glob(sfiles):
                self.link_fastqs()

            ## otherwise do the demultiplexing
            else:
                assemble.demultiplex.run(self, preview, ipyclient, force)



    def _step2func(self, samples, nreplace, force, preview, ipyclient):
        """ hidden wrapped function to start step 2"""

        ## print header
        if self._headers:
            print("  Step2: Filtering reads ")

        ## If no samples in this assembly then it means you skipped step1,
        ## so attempt to link existing demultiplexed fastq files
        if not self.samples.keys():
            raise IPyradWarningExit("""
    Error: No Samples found. First run step 1 to load raw or demultiplexed
    fastq data files from either the raw_fastq_path or sorted_fastq_path. """)

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        if not force:
            ## skip if all are finished
            if all([i.stats.state >= 2 for i in samples]):
                print("""\
    Skipping: All {} selected Samples already edited.
    (can overwrite with force argument)\
    """.format(len(samples)))
                return

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
                raise IPyradError("""
    {} assembly method requires a value for reference_sequence_path.
    """.format(self.paramsdict["assembly_method"]))
            else:
                ## index the reference sequence
                ## Allow force to reindex the reference sequence
                index_reference_sequence(self, force)

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 3, force):
            raise IPyradError("""
    No Samples ready to be clustered. First run step2().""")

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 3 for i in samples]):
                print("""\
    Skipping: All {} selected Samples already clustered.
    (can overwrite with force argument)\
    """.format(len(samples)))
                return

        ## run the step function
        assemble.cluster_within.run(self, samples, noreverse, 
                                    force, preview, ipyclient)


    def _step4func(self, samples, subsample, force, ipyclient):
        """ hidden wrapped function to start step 4 """

        if self._headers:
            print("  Step4: Joint estimation of error rate and heterozygosity")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 4, force):
            raise IPyradError("""
    No Samples ready for joint estimation. First run step3().""")

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 4 for i in samples]):
                print("""\
    Skipping: All {} selected Samples already joint estimated
    (can overwrite with force argument)\
    """.format(len(samples)))
                return

        ## send to function
        assemble.jointestimate.run(self, samples, subsample, force, ipyclient)


    def _step5func(self, samples, force, ipyclient):
        """ hidden wrapped function to start step 5 """
        ## print header
        if self._headers:
            print("  Step5: Consensus base calling ")

        ## Get sample objects from list of strings
        samples = _get_samples(self, samples)

        ## Check if all/none in the right state
        if not self.samples_precheck(samples, 5, force):
            raise IPyradError("""
    No Samples ready for consensus calling. First run step4().""")

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 5 for i in samples]):
                print("""\
    Skipping: All {} selected Samples already consensus called
    (can overwrite with force argument)\
    """.format(len(samples)))
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
            print("  Step6: Clustering across {} samples at {} similarity".\
                  format(len(csamples), self.paramsdict["clust_threshold"]))

        ## Check if all/none in the right state
        if not csamples:
            raise IPyradError("""
    No Samples ready for clustering. First run step5().""")

        elif not force:
            ## skip if all are finished
            if all([i.stats.state >= 6 for i in csamples]):
                print("""\
    Skipping: All {} selected Samples already clustered.
    (can overwrite with force argument)\
    """.format(len(samples)))
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
            print("  Step7: Filter and write output files for {} Samples".\
                  format(len(samples)))

        ## Check if all/none of the samples are in the self.database
        try:
            with h5py.File(self.clust_database, 'r') as ioh5:
                dbset = set(ioh5["seqs"].attrs['samples'])
                iset = set([i.name for i in samples])
                diff = iset.difference(dbset)
                if diff:
                    raise IPyradError("""
    The following Samples do not appear to have been clustered in step6
    (i.e., they are not in {}): 
    Missing: {}
    Check for typos in Sample names, or try running step6 including the 
    selected samples.
    """.format(self.database, ", ".join(list(diff))))

        except (IOError, ValueError):
            raise IPyradError("""
    Database file {} not found. First run step6
    """.format(self.database))

        if not force:
            if os.path.exists(
                os.path.join(self.dirs.project, self.name+"_outfiles")):
                raise IPyradWarningExit("""
    Output files already created for this Assembly in:
    {} 
    To overwrite, rerun using the force argument.""".format(self.dirs.outfiles))

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
        
        

    def step1(self, force=False, preview=False):
        """ docsting ... test """
        self._clientwrapper(self._step1func, [force, preview], 45)


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
                           [samples, nreplace, force, preview], 45)

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
                           [samples, noreverse, force, preview], 45)


    def step4(self, samples=None, subsample=None, force=False):
        """ test """
        self._clientwrapper(self._step4func, [samples, subsample, force], 45)



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

        self._clientwrapper(self._step5func, [samples, force], 45)



    def step6(self, samples=None, noreverse=False, force=False, randomseed=123):
        """ 
        Cluster consensus reads across samples and align with muscle. 

        Parameters
        ----------
        samples : list or str
            By default all Samples linked to an Assembly object are clustered. 
            If a subset of Samples is entered as a list then only those samples
            will be clustered. It is recommended to create .branch() Assembly 
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
                                              randomseed], 45)


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
        self._clientwrapper(self._step7func, [samples, force], 45)




    def run(self, steps=0, force=False, preview=False):
        """ Select steps of an analysis. If no steps are entered then all
        steps are run. Enter steps as a string, e.g., "1", "123", "12345" """
        if not steps:
            steps = list("1234567")
        else:
            if isinstance(steps, int):
                steps = str(steps)
            steps = list(steps)
        ## print a header in inside API
        if ip.__interactive__:
            print("\n  Assembly: {}".format(self.name))

        if '1' in steps:
            self.step1(force=force, preview=preview)
        if '2' in steps:
            self.step2(force=force, preview=preview)
        if '3' in steps:
            self.step3(force=force, preview=preview)
        if '4' in steps:
            self.step4(force=force)
        if '5' in steps:
            self.step5(force=force)
        if '6' in steps:
            self.step6(force=force)            
        if '7' in steps:
            self.step7(force=force)



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
    base, ext = os.path.splitext(os.path.basename(fname))

    ## remove read number from name
    base = base.replace("_R1_.", ".")\
               .replace("_R1_", "_")\
               .replace("_R1.", ".")

    ## remove extensions
    while ext in file_extensions:
        base, ext = os.path.splitext(base)

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


def expander(namepath):
    """ expand ./ ~ and ../ designators in location names """        
    if "~" in namepath:
        namepath = os.path.expanduser(namepath)
    else:
        namepath = os.path.abspath(namepath)
    return namepath



def merge(name, assemblies):
    """ Creates and returns a new Assembly object in which 
    samples from two or more Assembly objects with matching names
    are 'merged'. Merging does not affect the actual files written
    on disk, but rather creates new Samples that are linked to 
    multiple data files, and with stats summed. """

    ## checks
    assemblies = list(assemblies)

    ## create new Assembly
    merged = assemblies[0].branch(name)

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
        self.paramsdict["project_dir"] = expandpath
        self.dirs["project"] = expandpath

    elif param == 'raw_fastq_path':
        if newvalue:
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
                raise IPyradWarningExit("""
    Error: The value entered for the path to the raw fastq file is 
    unrecognized. Please be sure this path is correct. Double check the
    file name and the file extension. If it is a relative path be sure
    the path is correct with respect to the directory you're running ipyrad
    from. 
    You entered: {}
    """.format(fullrawpath))
        else:
            self.paramsdict['raw_fastq_path'] = ""


    elif param == 'barcodes_path':
        ## if a value was entered check that it exists
        if newvalue:
            fullbarpath = expander(newvalue)
        
            if not os.path.exists(fullbarpath):
                raise IPyradWarningExit("""
    Error: barcodes file not found. This must be an absolute path 
    (/home/wat/ipyrad/data/data_barcodes.txt) or relative to the directory 
    where you're running ipyrad (./data/data_barcodes.txt). You entered: 
    {}
    """.format(fullbarpath))
            else:
                self.paramsdict['barcodes_path'] = fullbarpath
                self.link_barcodes()

        ## if no path was entered then set barcodes path to empty. 
        ## this is checked again during step 1 and will raise an error 
        ## if you try demultiplexing without a barcodes file
        else:
            self.paramsdict['barcodes_path'] = ""


    elif param == 'sorted_fastq_path':
        if newvalue:
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
            if "3rad" in self.paramsdict['datatype']:
                self.link_barcodes()

    elif param == 'restriction_overhang':
        newvalue = tuplecheck(newvalue, str)                        
        assert isinstance(newvalue, tuple), """
    cut site must be a tuple, e.g., (TGCAG, '') or (TGCAG, CCGG)"""
        if len(newvalue) == 1:
            newvalue = (newvalue, "")
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
        ## do not allow majrule to be > statistical
        elif int(newvalue) < self.paramsdict["mindepth_majrule"]:
            raise IPyradError("""
    mindepth statistical cannot be less than mindepth_majrule""")
        else:
            self.paramsdict['mindepth_statistical'] = int(newvalue)
            ## TODO: calculate new clusters_hidepth if passed step3

    elif param == 'mindepth_majrule':
        assert isinstance(int(newvalue), int), \
            "mindepth_majrule must be an integer."
        if int(newvalue) > self.paramsdict["mindepth_statistical"]:
            print(\
        "error: mindepth_majrule cannot be > mindepth_statistical")
        else:
            ## TODO: calculate new clusters_hidepth if passed step3
            self.paramsdict['mindepth_majrule'] = int(newvalue)


    ## TODO: not yet implemented
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
        if newvalue.isdigit():
            self.paramsdict['max_shared_Hs_locus'] = int(newvalue)
        else:
            try:
                self.paramsdict['max_shared_Hs_locus'] = float(newvalue)
            except Exception as inst:
                sys.exit("max_shared_Hs_locs must be int or float, you put: "\
                        + newvalue)

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
            except ValueError:
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
        self.paramsdict['trim_overhang'] = newvalue


    elif param == 'output_formats':
        ## Get all allowed file types from assembly.write_outfiles
        output_formats = assemble.write_outfiles.OUTPUT_FORMATS

        if isinstance(newvalue, list):
            newvalue = ",".join(newvalue)

        ## If wildcard, then just do them all
        if "*" in newvalue:
            requested_formats = output_formats
            #"".join([i[0] for i in output_formats])
        else:
            #newvalue = newvalue.replace(",", "")
            #requested_formats = "".join([i[0] for i in newvalue.split()])
            requested_formats = [i for i in \
                                 newvalue.replace(" ", "").split(",")]
        
            ## Exit if requested formats are bad
            for form in requested_formats:
                if form not in [i for i in output_formats]:
                    raise IPyradWarningExit("""
    File format [{}] not recognized, must be one of: {}.
    """.format(form, output_formats))
        ## set the param
        self.paramsdict['output_formats'] = requested_formats


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




if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.get_params()
    print(DATA.log)
