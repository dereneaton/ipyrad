#!/usr/bin/env python2.7
""" Sample object """

#import os
import pandas as pd
from ipyrad.assemble.util import ObjDict, IPyradError

# pylint: disable=C0103
# pylint: disable=R0903


class Sample(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self, name=""):
        ## a sample name
        self.name = name
        self.barcode = ""
        self.merged = 0

        ## summary stats dictionary
        self.stats = pd.Series(
            index=["state",
                   "reads_raw",
                   "reads_filtered",
                   "reads_merged",
                   "refseq_mapped_reads",
                   "refseq_unmapped_reads",
                   "clusters_total",
                   "clusters_hidepth",
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
              "consens": [],
              "database": []
              })

        ## stats for each step
        self.statsfiles = ObjDict({
              "s1": pd.Series(index=["reads_raw"]),
              "s2": pd.Series(index=["reads_raw",
                                     "filter_qscore",
                                     "filter_adapter",
                                     "reads_passed"]),
              "s3": pd.Series(index=["reads_raw",
                                     "clusts_total",
                                     "clusts_hidepth",
                                     "avg.depth.tot",
                                     "avg.depth>mj",
                                     "avg.depth>stat"]),
              "s4": pd.Series(index=["hetero_est",
                                     "error_est"]),
              "s5": pd.Series(index=["nclusters",
                                     "depthfilter",
                                     "maxHfilter",
                                     "maxAllelefilter",
                                     "maxNfilter",
                                     "nconsens",
                                     "nsites",
                                     "nhetero",
                                     "heterozygosity"]),
              "s6": pd.Series(index=["null"]),
              "s7": pd.Series(index=["null"]),              
          })      

        ## store cluster depth information (biggest memory cost)
        self.depths = []


    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, dicto):
        if dicto == 'name':
            if self.stats.state:
                if self.stats.state > 1:
                    raise IPyradError("""
    Sample names cannot be changed after step2 without breaking file paths""")
                else:
                    self.__dict__.update(dicto)
        else:
            self.__dict__.update(dicto)

    #def save(self):
    #    """ pickle the data object """
    #    dillout = open(self.name+".dataobj", "wb")
    #    dill.dump(self, dillout)
    #    dillout.close()


