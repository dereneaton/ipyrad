#!/usr/bin/env python2.7
""" Sample object """

#import os
import pandas as pd
import dill
from ipyrad.assemble.worker import ObjDict


class Sample(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self, name=""):
        ## a sample name
        self.name = name
        self.barcode = ""
        self.merged = 0

        ## stats dictionary
        self.stats = pd.Series(
            index=["state",
                   "reads_raw",
                   "reads_filtered",
                   "refseq_mapped_reads",
                   "refseq_unmapped_reads",
                   "clusters_total",
                   "clusters_kept",
                   "hetero_est",
                   "error_est",
                   "reads_consens",])

        ## link to files
        self.files = ObjDict({
              "fastqs": [],
              "edits": [],
              "mapped_reads": [],
              "unmapped_reads": [],
              "clusters": [],
              "depths": [],
              "consens": [],
              "database": [],
              "stats": []
              })

        ## store cluster depth information
        self.depths = ObjDict()
        self.depths.total = []
        self.depths.mjmin = []
        self.depths.statmin = []                

        ## assignments for hierarchical clustering
        self.group = []


    #def __getstate__(self):
    #    return self.__dict__

    #def __setstate__(self, dicto):
    #    self.__dict__.update(dicto)

    #def save(self):
    #    """ pickle the data object """
    #    dillout = open(self.name+".dataobj", "wb")
    #    dill.dump(self, dillout)
    #    dillout.close()


