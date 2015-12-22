#!/usr/bin/env ipython2

""" ipyrad Assembly class object. """

# pylint: disable=E1101
# pylint: disable=E1103


from __future__ import print_function
import os
import time
import glob
import sys
import gzip
import dill
import copy
import subprocess
import pandas as pd
import ipyparallel as ipp

from collections import OrderedDict
from ipyrad.assemble.worker import ObjDict
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
    vsearch : str
        The path to the default vsearch executable. If not found, this can be 
        changed by setting [Assembly].vsearch = [newpath].
    muscle : str
        The path to the default muscle executable. If not found, this can be 
        changed by setting `[Assembly].muscle = [newpath]`.
    smalt : str
        The path to the default smalt executable. If not found, this can be 
        changed by setting `[Assembly].smalt = [newpath]`.
    samtools : str
        The path to the default samtools executable. If not found, this can be 
        changed by setting `[Assembly].samtools = [newpath]`.
    bedtools : str
        The path to the default bedtools executable. If not found, this can be 
        changed by setting `[Assembly].bedtools = [newpath]`.
    log : list
        A list of all modifications to the Assembly object and its Samples with
        time stamps. Use `print [Assembly].log` for easier viewing.
    dirs : dict
        Returns a dictionary with the location of directories that contain 
        linked Sample object files and stats results.

        
    Returns
    -------
    object
         A new assembly object is returned.


     """


    def __init__(self, name):

        ## obj name
        self.name = name
        print("  New Assembly: {}".format(self.name))

        ## stores ipcluster launch info
        self._ipclusterid = ""
        self._ipprofile = ""

        ## get binaries of dependencies
        self.vsearch, self.muscle, self.smalt, self.samtools, self.bedtools = getbins()

        ## link a log history of executed workflow
        self.log = []
        self._stamp(self.name+" created")
        self.statsfiles = ObjDict()

        ## samples linked 
        self.samples = ObjDict()

        ## multiplex files linked
        self.barcodes = ObjDict()

        ## an object for storing data directories for this Assembly
        self.dirs = ObjDict()

        ## storing final results
        self.loci = ObjDict()

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
                       ("restriction_overhang", ("TGCAG", "")),
                       ("max_low_qual_bases", 5),
                       ("engines_per_job", 4),
                       ("mindepth_statistical", 6), 
                       ("mindepth_majrule", 6), 
                       ("datatype", 'rad'), 
                       ("clust_threshold", .85),
                       ("minsamp", 4), 
                       ("max_shared_heterozygosity", .25), 
                       ("prefix_outname", self.name),
                       ("phred_Qscore_offset", 33),
                       ("max_barcode_mismatch", 1),
                       ("filter_adapters", 0), 
                       ("filter_min_trim_len", 35), 
                       ("ploidy", 2), 
                       ("max_stack_size", 1000),
                       ("max_Ns_consens", (5, 5)), 
                       ("max_Hs_consens", (8, 8)), 
                       ("max_SNPs_locus", (100, 100)), 
                       ("max_Indels_locus", (5, 99)), 
                       ("trim_overhang", (1, 2, 2, 1)), 
                       ("hierarchical_clustering", 0),
                       ("assembly_method", "denovo"),
                       ("reference_sequence", "")
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

    @property
    def files(self):
        """ Returns a data frame with Sample files. Not very readable... """
        nameordered = self.samples.keys()
        nameordered.sort()
        ## replace curdir with . for shorter printing
        #fullcurdir = os.path.realpath(os.path.curdir)
        return pd.DataFrame([self.samples[i].files for i in nameordered], 
                      index=nameordered).dropna(axis=1, how='all')


                      
    def _stamp(self, event):
        """ Stamps an event into the log history. """
        tev = time.strftime("%m/%d/%y %H:%M:%S", time.gmtime())
        self.log.append((self.name, tev, event))



    def link_fastqs(self, path=None, merged=False, force=False, append=False):
        """ Create Sample objects for samples in sorted_fastq_path.

        Note
        ----
        link_fastqs() is called automatically during step2() if no Samples
        are yet present in the Assembly object (data were not demultiplexed
        in step1().) It looks for demultiplexed data files located in the
        [sorted_fastq_path].


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

        Returns
        -------
        str
            Prints the number of new Sample objects created and the number of 
            fastq files linked to Sample objects in the Assembly object. 
        
        """

        ## cannot both force and append at once
        if force and append:
            raise Exception("Cannot use force and append at the same time.")

        if self.samples and not (force or append):
            raise Exception("Files already linked to `{}`. ".format(self.name)\
                +"Use force=True to replace all files, or append=True to "
                +"add additional files to existing Samples.")

        ## get path to data files
        if not path:
            path = self.paramsdict["sorted_fastq_path"]

        ## make sure there is an out directory
        self.dirs.fastqs = os.path.join(self.paramsdict["working_directory"],
                                        self.name+"_fastqs")
        if not os.path.exists(self.paramsdict["working_directory"]):
            os.mkdir(self.paramsdict["working_directory"])
        if not os.path.exists(self.dirs.fastqs):
            os.mkdir(self.dirs.fastqs)

        ## does location exist, if no files selected, try selecting all
        if os.path.isdir(path):
            path += "*"

        ## grab fastqs/fq/gzip/all
        fastqs = glob.glob(path)
        fastqs = [i for i in fastqs if i.endswith(".gz") \
                                    or i.endswith(".fastq") \
                                    or i.endswith(".fq")]

        ## sort alphabetical
        fastqs.sort()

        ## link pairs into tuples        
        if 'pair' in self.paramsdict["datatype"]:
            ## check that names fit the paired naming convention
            r1_files = [i for i in fastqs if "_R1_" in i]
            r2_files = [i.replace("_R1_", "_R2_") for i in r1_files]

            if r1_files:
                if not any(["_R1_" in i for i in fastqs]) or \
                       (len(r1_files) != len(r2_files)):
                    raise Exception(\
                "File name format error: paired file names " \
                +"must be identical except for _R1_ and _R2_ in their names.")
            fastqs = [(i, j) for i, j in zip(r1_files, r2_files)]

        ## data are not paired, create empty tuple pair
        else:
            if any(["_R2_" in i for i in fastqs]):
                print("Given the presence of '_R2_' in file names, this "\
              +"is a warning that if your data are paired-end you should set "\
              +"the Assembly object datatype to a paired type (e.g., "\
              +"pairddrad or pairgbs) prior to running link_fastqs().")
            fastqs = [(i, ) for i in fastqs]

        ## counters for the printed output
        created = 0
        linked = 0
        appended = 0
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
                    print("The files {} are already in Sample.".format(sname) \
                    + " Use append=True to append additional files to a Sample"\
                    + " or force=True to replace all existing Samples.")

            ## record whether data were merged.
            if merged:
                self.samples[sname].merged = 1

            ## do not allow merged=False and .forward in file names
            if (merged == False) and ('forward' in fastqtuple[0]):
                print(\
                "If R1 and R2 data are merged (e.g., with PEAR) " \
              + "use link_fastqs(merge=True) to indicate this. You " \
              + "may need force=True to overwrite existing files.\n")

            ## if fastqs already demultiplexed, try to link stats
            if any([linkedinc, createdinc, appendinc]):
                gzipped = bool(fastqtuple[0].endswith(".gz"))
                nreads = 0
                ## iterate over files if there are multiple
                for alltuples in self.samples[sname].files.fastqs:
                    nreads += bufcount(alltuples[0], gzipped)
                self.samples[sname].stats.reads_raw = nreads/4
                created += createdinc
                linked += linkedinc
                appended += appendinc

        ## print if data were linked
        print("{} new Samples created in `{}`.".format(created, self.name))
        if linked:
            print("{} fastq files linked to {} new Samples.".\
                  format(linked, len(self.samples)))
        if appended:
            print("{} fastq files appended to {} existing Samples.".\
                  format(appended, len(self.samples)))



    def link_barcodes(self):
        """ creates a self.barcodes object to save barcodes info 
            as a dictionary, if there is a barcodes file in 
            self.paramsdict["barcodes_path"] """
        ## in case fuzzy selected
        try: 
            barcodefile = glob.glob(self.paramsdict["barcodes_path"])[0]
        except IndexError: 
            print("Barcodes file not found:", self.paramsdict["barcodes_path"])

        ## parse barcodefile
        bdf = pd.read_csv(barcodefile, header=None, delim_whitespace=1)
        bdf = bdf.dropna()
        ## make sure upper case
        bdf[1] = bdf[1].str.upper()
        ## set attribute on Assembly object
        self.barcodes = dict(zip(bdf[0], bdf[1]))



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
        [Assembly].set_params('max_shared_heterozygosity', 0.25)
            
        """

        ## require parameter recognition
        assert (param in range(50)) or \
               (param in [str(i) for i in range(50)]) or \
               (param in self.paramsdict.keys()), \
            "Parameter key not recognized: `{}`.".format(param)

        ## make string
        param = str(param)

        ## if matching
        if param in ['1', 'working_directory']:
            self.paramsdict['working_directory'] = expander(newvalue)
            self._stamp("[1] set to "+newvalue)
            self.dirs["working"] = self.paramsdict["working_directory"]


        elif param in ['2', 'raw_fastq_path']:
            fullrawpath = expander(newvalue)
            if os.path.isdir(fullrawpath):
                fullrawpath = os.path.join(fullrawpath, "*.gz")
            self.paramsdict['raw_fastq_path'] = fullrawpath
            self._stamp("[2] set to "+newvalue)


        elif param in ['3', 'barcodes_path']:
            #assert type(newvalue) is StringType, "arg must be a string"
            fullbarpath = expander(newvalue)
            if glob.glob(fullbarpath):
                self.paramsdict['barcodes_path'] = fullbarpath
                self.link_barcodes()
                self._stamp("[3] set to "+newvalue)
            elif not fullbarpath:
                self.paramsdict['barcodes_path'] = fullbarpath                
                self._stamp("[3] set to empty")
            else:
                print(
        "Warning: barcodes file not found. This must be an absolute \n"\
       +"path (/home/wat/ipyrad/data/data_barcodes.txt) or relative to the \n"\
       +"directory where you're running ipyrad (./data/data_barcodes.txt).\n"\
       +"You entered: %s\n" % fullbarpath)


        elif param in ['4', 'sorted_fastq_path']:
            assert isinstance(newvalue, str), \
            "sorted_fastq_path must be a string, e.g., /home/data/fastqs/*"
            newvalue = expander(newvalue)
            if os.path.isdir(newvalue):
                newvalue = os.path.join(newvalue, "*.gz")
            self.paramsdict['sorted_fastq_path'] = newvalue
            self._stamp("[4] set to "+newvalue)


        elif param in ['5', 'restriction_overhang']:
            newvalue = tuplecheck(newvalue, str)                        
            assert isinstance(newvalue, tuple), \
            "cut site must be a tuple, e.g., (TGCAG, '') or (TGCAG, CCGG)"
            if len(newvalue) == 1:
                newvalue = (newvalue, "")
            assert len(newvalue) == 2, \
            "must enter 1 or 2 cut sites, e.g., (TGCAG, '') or (TGCAG, CCGG)"
            self.paramsdict['restriction_overhang'] = newvalue
            self._stamp("[5] set to "+str(newvalue))


        elif param in ['6', 'max_low_qual_bases']:
            self.paramsdict['max_low_qual_bases'] = int(newvalue)
            self._stamp("[6] set to "+str(newvalue))


        elif param in ['7', "engines_per_job"]:
            self.paramsdict['engines_per_job'] = int(newvalue)
            self._stamp("[7] set to "+str(newvalue))


        elif param in ['8', 'mindepth_statistical']:
            ## do not allow values below 5
            if int(newvalue) < 5:
                print("error: mindepth statistical cannot be set < 5")
            ## do not allow majrule to be > statistical
            elif int(newvalue) < self.paramsdict["mindepth_majrule"]:
                print("error: mindepth statistical cannot be less than \
                       mindepth_majrule")                
            else:
                self.paramsdict['mindepth_statistical'] = int(newvalue)
                self._stamp("[8] set to "+str(newvalue))


        elif param in ['9', 'mindepth_majrule']:
            if int(newvalue) > self.paramsdict["mindepth_statistical"]:
                print("error: mindepth_majrule cannot be > \
                       mindepth_statistical")
            else:
                self.paramsdict['mindepth_majrule'] = int(newvalue)
                self._stamp("[9] set to "+str(newvalue))


        elif param in ['10', 'datatype']:
            ## list of allowed datatypes
            datatypes = ['rad', 'gbs', 'ddrad', 'pairddrad',
                         'pairgbs', 'merged', '2brad']
            ## raise error if something else
            if self.paramsdict['datatype'] not in datatypes:
                print("error: datatype not recognized")
            else:
                self.paramsdict['datatype'] = str(newvalue)
                self._stamp("[10] set to "+newvalue)


        elif param in ['11', 'clust_threshold']:
            self.paramsdict['clust_threshold'] = float(newvalue)
            self._stamp("[11] set to {}".format(newvalue))


        elif param in ['12', 'minsamp']:
            self.paramsdict['minsamp'] = int(newvalue)
            self._stamp("[12] set to {}".format(int(newvalue)))


        elif param in ['13', 'max_shared_heterozygosity']:
            self.paramsdict['max_shared_heterozygosity'] = newvalue
            self._stamp("[13] set to {}".format(newvalue))


        elif param in ['14', 'prefix_outname']:
            self.paramsdict['prefix_outname'] = newvalue
            self._stamp("[14] set to {}".format(newvalue))


        elif param in ['15', 'phred_Qscore_offset']:
            self.paramsdict['phred_Qscore_offset'] = int(newvalue)
            self._stamp("[15] set to {}".format(int(newvalue)))


        elif param in ['16', 'max_barcode_mismatch']:
            self.paramsdict['max_barcode_mismatch'] = int(newvalue)
            self._stamp("[16] set to {}".format(int(newvalue)))

        ### ....
        elif param in ['17', 'filter_adapters']:
            self.paramsdict['filter_adapters'] = int(newvalue)
            self._stamp("[17] set to "+str(newvalue))


        elif param in ['18', 'filter_min_trim_len']:
            self.paramsdict['filter_min_trim_len'] = int(newvalue)
            self._stamp("[18] set to {}".format(int(newvalue)))


        elif param in ['19', 'ploidy']:
            self.paramsdict['ploidy'] = int(newvalue)
            self._stamp("[19] set to {}".format(int(newvalue)))


        elif param in ['20', 'max_stack_size']:
            self.paramsdict['max_stack_size'] = int(newvalue)
            self._stamp("[20] set to {}".format(int(newvalue)))


        elif param in ['21', 'max_Ns_consens']:
            newvalue = tuplecheck(newvalue)                        
            assert isinstance(newvalue, tuple), \
            "max_Ns_consens should be a tuple e.g., (8,8)"
            self.paramsdict['max_Ns_consens'] = newvalue
            self._stamp("[21] set to {}".format(newvalue))


        elif param in ['22', 'max_Hs_consens']:
            newvalue = tuplecheck(newvalue)                        
            assert isinstance(newvalue, tuple), \
            "max_Hs_consens should be a tuple e.g., (1,2,2,1)"
            self.paramsdict['max_Hs_consens'] = newvalue
            self._stamp("[22] set to {}".format(newvalue))


        elif param in ['23', 'max_SNPs_locus']:
            newvalue = tuplecheck(newvalue)                        
            assert isinstance(newvalue, tuple), \
            "max_SNPs_locus should be a tuple e.g., (20,20)"
            self.paramsdict['max_SNPs_locus'] = newvalue
            self._stamp("[23] set to {}".format(newvalue))


        elif param in ['24', 'max_Indels_locus']:
            newvalue = tuplecheck(newvalue)            
            assert isinstance(newvalue, tuple), \
            "max_Indels_locus should be a tuple e.g., (5, 100)" 
            self.paramsdict['max_Indels_locus'] = newvalue
            self._stamp("[24] set to {}".format(newvalue))


        elif param in ['25', 'trim_overhang']:
            newvalue = tuplecheck(newvalue)
            assert isinstance(newvalue, tuple), \
            "trim_overhang should be a tuple e.g., (1,2,2,1)"
            self.paramsdict['trim_overhang'] = newvalue
            self._stamp("[25] set to {}".format(newvalue))


        elif param in ['27', 'assembly_method']:
            self.paramsdict['assembly_method'] = newvalue
            assert self.paramsdict['assembly_method'] in \
                              ["denovo", "reference", "hybrid"], \
                 "The assembly_method setting must be one of the following: "+\
                 "denovo, reference, or hybrid."
            self._stamp("[27] set to {}".format(newvalue))


        elif param in ['28', 'reference_sequence']:
            fullrawpath = expander(newvalue)
            if not os.path.isfile(fullrawpath):
                LOGGER.info("reference sequence file not found.")
                print(\
        "Warning: reference sequence file not found. This must be an\n"\
       +"absolute path (/home/wat/ipyrad/data/reference.gz) or relative to \n"\
       +"the directory where you're running ipyrad (./data/reference.gz).\n"\
       +"You entered: %s\n" % fullrawpath)
            self.paramsdict['reference_sequence'] = fullrawpath
            self._stamp("[28] set to "+fullrawpath)



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
            newobj.set_params('prefix_outname', newname)

            ## create copies of each Sample obj
            for sample in self.samples:
                newobj.samples[sample] = copy.deepcopy(self.samples[sample])
            return newobj



    def file_tree(self):
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



    def _launch(self, tries):
        """ launch ipyclient.
        launch within try statement in case engines aren't ready yet
        and try 30 one second sleep/wait cycles before giving up on engines
        """
        while tries:
            try:
                ## launches ipcluster with arguments if present in self
                clusterargs = [self._ipclusterid, self._ipprofile]
                argnames = ["--cluster-id", "--profile"]
                args = {key:value for key, value in zip(argnames, clusterargs)}
                ipyclient = ipp.Client(**args)

                if tries > 1:
                    LOGGER.debug('try %s: starting controller', tries)                  
                tries = 0
                ## make sure all engines are connected
                try:
                    assert ipyclient.ids
                    LOGGER.info('OK! Connected to (%s) engines', 
                                len(ipyclient.ids))
                except AssertionError as _: 
                    LOGGER.debug('connected to %s engines', len(ipyclient.ids))
                    raise
            except (IOError, ipp.NoEnginesRegistered, AssertionError):
                time.sleep(1)
                tries -= 1
        try:
            return ipyclient
        except UnboundLocalError as inst:
            LOGGER.error(inst)
            SystemExit("Error: failed to launch ipcluster")



    def step1(self, force=False, preview=0):
        """ step 1: demultiplex raw reads or link pre-demultiplexed files."""

        ## launch parallel client within guarded statement
        try: 
            ipyclient = self._launch(30) 
            msg1 = "  step1: Demultiplexing raw reads."
            msg2 = "  step1: Linking fastq data to Samples."
            ## if Samples already exist then no demultiplexing
            if self.samples:
                if not force:
                    print("Skipping step1: {} ".format(len(self.samples)) \
                         +"Samples already found in `{}` ".format(self.name) \
                         +"(see --force).")
                else:
                    assemble.demultiplex.run(self, preview, ipyclient)
                    self._stamp("s1_demultiplexing:")
                    print(msg1)
                    #print("Overwriting Samples with data files in {}".\
                    #      format(self.paramsdict["sorted_fastq_path"]))
                    #self.link_fastqs(force=True)
            ## Creating new Samples
            else:
                ## first check if demultiplexed files exist in path
                if os.path.exists(self.paramsdict["sorted_fastq_path"]):
                    try:
                        self.link_fastqs(path=os.path.join(
                            self.paramsdict["sorted_fastq_path"], "*"))
                        print("linking files from {}".\
                              format(self.paramsdict["sorted_fastq_path"]))
                    except AssertionError as inst:
                        SystemExit("inst")
                        raise
                ## otherwise do the demultiplexing
                else:
                    assemble.demultiplex.run(self, preview, ipyclient)
                    self._stamp("s1_demultiplexing:")
                    print(msg1)

        except (KeyboardInterrupt, SystemExit, AttributeError):
            logging.error("assembly step1 interrupted.")
            raise
        except Exception as inst:
            print("other error: %s" % inst)
            raise

        ## close client when done or if interrupted
        finally:
            try:
                ipyclient.close()
            except UnboundLocalError as inst:
                LOGGER.error("problem launching ipcluster: \n(%s)\n", inst)

        ## pickle the data obj
        self._save()



    def step2(self, samples="", preview=0, force=False):
        """ step 2: edit raw reads. Takes dictionary keys (sample names)
        either individually, or as a list, or it takes no argument to 
        select all samples in the Assembly object. Only samples in state
        =1 will be edited, all others are skipped. To overwrite data
        use the argument force=True. 

        """
        print("  Step2: Filtering reads ")
        ## launch parallel client within guarded statement
        try: 
            ipyclient = self._launch(30) 

            ## if samples not entered use all samples
            if not samples:
                samples = self.samples.keys()

            if not samples:
                self.link_fastqs()

            ## require Samples 
            assert samples, "No Samples in {}".format(self.name)

            ## if sample keys, replace with sample obj
            assert isinstance(samples, list), \
            "to subselect samples enter as a list, e.g., [A, B]."
            samples = [self.samples[key] for key in samples]

            ## pass samples to rawedit
            assemble.rawedit.run(self, samples, ipyclient, force)

        except (KeyboardInterrupt, SystemExit):
            LOGGER.error("assembly step2 interrupted by user.")
            raise
            
        except Exception as inst:
            LOGGER.error("crashed step2 error: %s", inst)
            print("error during step2 ...\n({})\n".format(inst))
            raise

        ## close parallel client if done or interrupted
        finally:
            try:
                ipyclient.close()
            except UnboundLocalError as inst:
                LOGGER.error("problem launching ipcluster: \n(%s)\n", inst)

        ## checkpoint the data obj
        self._save()



    def step3(self, samples=None, preview=0, noreverse=0, force=False):
        """ step 3: clustering within samples """
        print("  Step3: Clustering/Mapping reads")
        ## Require reference seq for reference-based methods
        if self.paramsdict['assembly_method'] != "denovo":
            assert self.paramsdict['reference_sequence'], \
            "Reference or hybrid assembly requires a value for "+\
            "reference_sequence_path paramter."

            ## index the reference sequence
            index_reference_sequence(self)

        try: 
            ipyclient = self._launch(30)

            ## if samples not entered use all samples
            if not samples:
                samples = self.samples.keys()

            ## if sample keys, replace with sample obj
            assert isinstance(samples, list), \
            "to subselect samples enter as a list, e.g., [A, B]."
            samples = [(key, self.samples[key]) for key in samples]

            ## require Samples 
            assert samples, "No Samples in {}".format(self.name)

            ## skip if all are finished
            if not force:
                if all([i[1].stats.state >= 3 for i in samples]):
                    print("Skipping step3: All {} ".format(len(self.samples))\
                         +"Samples already clustered in `{}`".format(self.name))
                else:
                    assemble.cluster_within.run(self, samples, ipyclient, 
                                                preview, noreverse, force)
            else:
                assemble.cluster_within.run(self, samples, ipyclient, 
                                            preview, noreverse, force)


        except (KeyboardInterrupt, SystemExit) as _:
            print("assembly step3 interrupted by user.")
            raise
        except Exception as inst:
            print("other error: %s" % inst)
            raise

        ## close parallel client if done or interrupted
        try:
            ipyclient.close()
        except UnboundLocalError as inst:
            LOGGER.error("problem launching ipcluster: \n(%s)\n", inst)

        ## pickle the data object
        self._save()



    def step4(self, samples=None, preview=0, force=False, subsample=None):
        """ step 4: Joint estimation of error rate and heterozygosity. 
        If you want to overwrite data for a file, first set its state to 3:
        data.samples['sample'].stats['state'] = 3 """
        print("  Step4: Joint estimation of error rate and heterozygosity")
        try: 
            ## launch parallel client            
            ipyclient = self._launch(30) 

            ## sampling
            if samples:
                ## make a list keys or samples
                if isinstance(samples, str):
                    samples = list([samples])
                else:
                    samples = list(samples)

                ## if keys are in list
                if any([isinstance(i, str) for i in samples]):
                    ## make into a subsampled sample dict
                    subsamples = {i: self.samples[i] for i in samples}

                ## send to function
                assemble.jointestimate.run(self, subsamples.values(), 
                                           ipyclient, force, subsample)
            else:
                ## if no sample, then do all samples
                if not self.samples:
                    print("Assembly object has no Samples.")
                ## run clustering for all samples
                assemble.jointestimate.run(self, self.samples.values(), 
                                           ipyclient, force, subsample)

        except (KeyboardInterrupt, SystemExit) as inst:
            LOGGER.error("assembly step4 interrupted: %s", inst)
            print("assembly step4 interrupted.")
            raise

        except Exception as inst:
            LOGGER.error("crashed step4 error: %s", inst)
            print("error during step4 ...\n({})\n".format(inst))
            raise

        ## close parallel client if done or interrupted
        finally:
            ipyclient.close()
            if preview:
                print(".")

        ## pickle the data object
        self._save()




    def step5(self, samples="", preview=0, force=False):
        """ step 5: Consensus base calling from clusters within samples.
        If you want to overwrite data for a file, first set its state to 
        3 or 4. e.g., data.samples['sample'].stats['state'] = 3 """

        ## launch parallel client within guarded statement
        try: 
            ipyclient = self._launch(30) 

            ## if samples not entered use all samples
            if not samples:
                samples = self.samples.keys()

            ## require Samples 
            assert samples, "No Samples in {}".format(self.name)

            ## if sample keys, replace with sample obj
            assert isinstance(samples, list), \
            "to subselect samples enter as a list, e.g., [A, B]."
            samples = [self.samples[key] for key in samples]

            ## pass samples to rawedit
            assemble.consens_se.run(self, samples, ipyclient, force)

        except (KeyboardInterrupt, SystemExit):
            LOGGER.error("assembly step5 interrupted by user.")
            raise
            
        except Exception as inst:
            LOGGER.error("crashed step5 error: %s", inst)
            print("error during step5 ...\n({})\n".format(inst))
            raise

        ## close parallel client if done or interrupted
        finally:
            try:
                ipyclient.close()
            except UnboundLocalError as inst:
                LOGGER.error("problem launching ipcluster: \n(%s)\n", inst)

        ## checkpoint the data obj
        self._save()


    def step6(self, samples="", noreverse=False, force=False):
        """ Step 6: Cluster consensus reads across samples and 
        align with muscle """
        try: 
            ipyclient = self._launch(30) 

            ## if samples not entered use all samples
            if not samples:
                samples = self.samples.keys()

            ## if sample keys, replace with sample obj
            assert isinstance(samples, list), \
            "to subselect samples enter as a list, e.g., [A, B]."
            samples = [(key, self.samples[key]) for key in samples]

            ## require Samples 
            assert samples, "No Samples in {}".format(self.name)

            ## file for all reads
            self.loci.clusters = os.path.join(self.dirs.consens, 
                                              "cat.clust.gz")
            if os.path.exists(self.loci.clusters):
                print("Skipping step6: Clust file already exists:"\
                     +"{}\n".format(self.loci.clusters))
            else:
                assemble.cluster_across.run(self, samples, ipyclient, 
                                            noreverse, force)

        except (KeyboardInterrupt, SystemExit) as _:
            print("assembly step3 interrupted by user.")
            raise

        ## close parallel client if done or interrupted
        try:
            ipyclient.close()
        except UnboundLocalError as inst:
            LOGGER.error("problem launching ipcluster: \n(%s)\n", inst)

        ## pickle the data object
        self._save()


    def step7():
        """ Apply filters on aligned clusters to get final data set and 
        provide output files in a number of formats """





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
            self.step4(force=force, preview=preview)            
        if '5' in steps:
            self.step5(force=force, preview=preview)            
        # if '6' in steps:
        #     self.step6()            
        # if '7' in steps:
        #     self.step7()



def _name_from_file(fname):
    """ internal func: get the sample name from any pyrad file """
    file_extensions = [".gz", ".fastq", ".fq", ".fasta", 
                       ".clustS", ".consens"]
    base, ext = os.path.splitext(os.path.basename(fname))
    ## remove read number from name
    base = base.replace("_R1_", "")\
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
    assert cmd_exists(bedtools), "samtools not found here: "+bedtools
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



def bufcount(filename, gzipped):
    """ fast line counter """
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



def index_reference_sequence(self):
    """ Attempt to index the reference sequence. This is a little naive
    in that it'll actually _try_ do to the reference every time, but it's
    quick about giving up if it detects the indices already exist. You could
    also test for existence of both index files, but i'm choosing to just let
    smalt do that for us ;) """

    refseq_file = self.paramsdict['reference_sequence']

    #TODO: Here test if the indices exist already
    # These are smalt specific index files. We don't ever reference
    # them directly except here to make sure they exist, so we don't need
    # to keep them around.
    index_sma = refseq_file+".sma"
    index_smi = refseq_file+".smi"

    if not os.path.isfile(index_sma) or not os.path.isfile(index_smi):
        cmd = self.smalt+\
            " index "\
            " -s 2 "+\
	    refseq_file+" "+\
	    refseq_file

        print(cmd)
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

def tuplecheck(newvalue, conv=int):
    """ Takes a string argument and returns value as a tuple. 
    Needed for paramfile conversion from CLI to set_params args """
    if isinstance(newvalue, str):
        newvalue = newvalue.rstrip(")").strip("(")
        try:
            newvalue = tuple([conv(i) for i in newvalue.split(",")])
        except TypeError:
            newvalue = tuple(conv(newvalue))
        except Exception as inst:
            LOGGER.info(inst)
            sys.exit("\nError: arg `{}` is not formatted correctly.\n({})\n"\
                     .format(newvalue, inst))
    return newvalue



if __name__ == "__main__":
    ## test...
    DATA = Assembly("test")
    DATA.get_params()
    DATA.set_params(1, "./")
    DATA.get_params()
    print(DATA.log)
