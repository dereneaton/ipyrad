#!/usr/bin/env python

"""Load SNPs array and SNPsMAP from ipyrad SNPs HDF5 database.

These data can be used to create subsampled SNP data sets, or for 
imputation, in other ipyrad.analysis tools.
"""


from typing import Optional, Dict, List, Union
from dataclasses import dataclass

import h5py
import numpy as np
import pandas as pd
from loguru import logger
from ipyrad.assemble.utils import IPyradError
from ipyrad.analysis.utils import jsubsample_snps, jsubsample_loci
from ipyrad.core.cluster import Cluster
from ipyrad.analysis.progress import ProgressBar


@dataclass
class SNPsExtracter:
    """...

    """
    data: Assembly
    imap: Optional[Dict[str, List[str]]] = None
    minmap: Optional[Dict[str, int]] = None
    mincov: Union[float, int] = 0.0
    minmaf: Union[float, int] = 0.0 
