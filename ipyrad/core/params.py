


import os
import glob
from ..assemble.utils import IPyradError, IPyradParamsError


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



class Hackers(object):
    def __init__(self):

        # private dictionary so we can check dtypes before changing 
        self.__dict__["_data"] = dict([
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
            ("exclude_reference", False),
        ])

    # pretty printing of object
    def __repr__(self):
        printstr = ""
        for idx, (key, val) in enumerate(self._data.items()):
            printstr += "{:<4}{:<28}{:<45}\n".format(idx, key, str(val))
        return printstr
    
    def __str__(self):
        return self.__repr__()


    def __setattr__(self, key, value):
        "Checks keys during setattr's"
        if key == "_data":
            pass
        elif key in self._data:
            super(Hackers, self).__setattr__(key, value)
        else:
            print("'{}' is not a valid parameter".format(key))


    def __setitem__(self, key, value):
        "Checks keys during setattr's"
        if key == "_data":
            pass
        if key in self._data:
            super(Hackers, self).__setattr__(key, value)
        else:
            print("'{}' is not a valid parameter".format(key))


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
        return self._data["bwa_args"]
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

    @property
    def exclude_reference(self):
        return self._data["exclude_reference"]
    @exclude_reference.setter
    def exclude_reference(self, value):
        self._data["exclude_reference"] = bool(value)
    


class Params(object):
    def __init__(self, data):

        # harder to 'update' values if data is here...
        self._data = data

        # DEFAULT VALUES
        self._assembly_name = data.name
        self._project_dir = os.path.realpath("./")
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
        self._max_Ns_consens = 0.05
        self._max_Hs_consens = 0.05
        self._min_samples_locus = 4
        self._max_SNPs_locus = 0.20
        self._max_Indels_locus = 8
        self._max_shared_Hs_locus = 0.5
        self._trim_reads = (0, 0, 0, 0)
        self._trim_loci = (0, 0, 0, 0)
        self._output_formats = list("psl")
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
        

    # def __setattr__(self):

    # def update(self, dict):



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
            if self.barcodes_path:
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
            if self.barcodes_path:
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
        # warning if old style params
        if isinstance(value, (tuple, int)):
            print(
    "Warning: The format of max_Ns_consens should now be a float " + 
    "and was set on load to 0.05")
            value = 0.05
        self._max_Ns_consens = float(value)


    @property
    def max_Hs_consens(self):
        return self._max_Hs_consens
    @max_Hs_consens.setter
    def max_Hs_consens(self, value):
        # complain if old params format
        if isinstance(value, tuple):
            print(
    "Warning: The format of max_Hs_consens should now be a float " + 
    "and was set on load to 0.05")
            value = 0.05
        self._max_Hs_consens = float(value)


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
    max_shared_Hs_locus must be int or float, you put: {}""".format(value))
        self._max_shared_Hs_locus = value


    @property
    def max_SNPs_locus(self):
        return self._max_SNPs_locus
    @max_SNPs_locus.setter
    def max_SNPs_locus(self, value):
        #backwards compatible allow tuple and take first value
        if isinstance(value, tuple):
            value = value[0]
        if isinstance(value, str):
            if value.isdigit():
                value = int(value)
            else:
                try:
                    value = float(value)
                except Exception as inst:
                    raise IPyradParamsError("""
    max_SNPs_locus must be int or float, you put: {}""".format(value))
        self._max_SNPs_locus = value


    @property
    def max_Indels_locus(self):
        return self._max_Indels_locus
    @max_Indels_locus.setter
    def max_Indels_locus(self, value):
        if isinstance(value, tuple):
            value = value[0]
        try:
            value = int(value)
        except ValueError:
            raise IPyradError("max_Indels_locus should be an integer value e.g., 5. You entered: {}".format(value))
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



CANNOT_CHANGE_ASSEMBLY_NAME = """\
Warning: Assembly name is set at Assembly creation time and is an immutable
property: You may, however, branch the assembly which will create a copy
with a new name, but retain a copy of the original Assembly. Here's how:

Command Line Interface:
    ipyrad -p params-old-name.txt -b new-name

API (Jupyter Notebook Users):
    new_assembly = my_assembly.branch("new_name")
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

BAD_PROJDIR_NAME = """\
Error: Your project_dir contains a directory with a space in the name.
This can cause all kinds of funny problems so please rename this
directory and remove the space. Try replacing the space with an underscore.
You entered: {}
"""

BARCODE_NOT_FOUND = """\
Error: barcodes file not found. This must be an absolute path
(/home/wat/ipyrad/data/data_barcodes.txt) or relative to the directory
where you're running ipyrad (./data/data_barcodes.txt). You entered:
{}
"""

BAD_ASSEMBLY_METHOD = """\
The assembly_method parameter must be one of the following: denovo, reference,
denovo+reference or denovo-reference. You entered:
{}.
"""

RAW_PATH_ISDIR = """\
You entered the path to a directory for raw_fastq_path. To ensure the correct
files in the directory are selected, please use a wildcard selector to 
designate the desired files. 
Example: /home/user/data/*.fastq  ## selects all files ending in '.fastq'
You entered: {}
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
