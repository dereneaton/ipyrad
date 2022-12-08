#!/usr/bin/env python

"""Load a h5 database file in ipyrad.

This contains a simple function for loading an h5 file produced
by ipyrad. Importantly, it checks the version of the file, since the
format of the h5 database files changed after v.1.0. 

Example
-------
Instead of using the h5py contextmanager like this:
>>> with h5py.File(...) as io5:
>>>     ...

You should instead use the ipa context manager like this:
>>> with ipa.SnpsFile(...) as io5:
>>>     ...
"""

from typing import List
from pathlib import Path
from loguru import logger
import h5py

logger = logger.bind(name="ipyrad")


MSG_TO_CONVERT_DATABASE_FILES = """
  Your HDF5 file is from ipyrad version <1.0 and must be upgraded.
  Run the following to convert the file to the new format from the 
  command line. This will create a new file named `{name}_v1.snps_hdf5`.

  >>> ipyrad convert {your-current-file}
"""


class SnpsDatabase(h5py.File):
    """Child class of h5py.File for loading ipyrad h5 SNP files.

    Checks version for supported format and logs to debug. Also it
    converts all attrs names into strings, not bytes...
    """
    def __init__(self, name: Path, *args, **kwargs):
        self._path = Path(name)
        super().__init__(name, *args, **kwargs)

        self._check_file_type()
        self._check_version()
        self.summary("DEBUG")

    def _check_version(self) -> None:
        """if version is <1.0 then do not allow reading the file."""
        if "version" not in self.attrs:
            raise ValueError(MSG_TO_CONVERT_DATABASE_FILES)

    def _check_file_type(self) -> None:
        """if user entered seqs file raise exception asking for snps file."""
        if "hdf5" in self._path.suffix:
            if "seqs" in self._path.suffix:
                raise IOError("data should be .snps_hdf5 file not .seqs_hdf5.")
    
    def summary(self, log_level: str="INFO") -> None:
        """Prints a summary of the database attrs and keys to LOGGER INFO."""
        version = self.attrs['version']
        attrs = list(self.attrs.keys())
        keys = list(self.keys())        
        logger.log(log_level, f"loading h5 snps database v{version}; attrs={attrs}; keys={keys}")

    def get_sample_names(self) -> List[str]:
        """Return a list of sample names in database order."""
        return self.attrs["names"]


class SeqsDatabase(SnpsDatabase):

    def _check_file_type(self) -> None:
        """if user entered seqs file raise exception asking for snps file."""
        if "hdf5" in self._path.suffix:
            if ".snps" in self._path.suffix:
                raise IOError("data should be .snps_hdf5 file not .seqs_hdf5.")



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")

    FILE = "/home/deren/Documents/ipyrad/isolation/test-simrad_outfiles/test-simrad.seqs.hdf5"
    FILE = "/home/deren/Documents/ipyrad/sra-fastqs/cyatho_outfiles/cyatho.snps_hdf5"
    oseqs5 = SnpsDatabase(FILE, 'r')
