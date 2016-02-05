#!/usr/bin/env python2.7
""" Sample object """

#import os
import pandas as pd
import json
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
              "s1": pd.Series(index=["reads_raw",
                                     ]),
              "s2": pd.Series(index=["reads_raw",
                                     "filtered_by_qscore",
                                     "filtered_by_adapter",
                                     "reads_passed",
                                     ]),
              "s3": pd.Series(index=["merged_pairs",
                                     "clusters_total",
                                     "clusters_hidepth",
                                     "avg_depth_tot",
                                     "avg_depth_mj",
                                     "avg_depth_stat",
                                     ]),
              "s4": pd.Series(index=["hetero_est",
                                     "error_est",
                                     ]),
              "s5": pd.Series(index=["nclusters",
                                     "depthfilter",
                                     "maxHfilter",
                                     "maxAllelefilter",
                                     "maxNfilter",
                                     "nconsens",
                                     "nsites",
                                     "nhetero",
                                     "heterozygosity",
                                     ]),
              "s6": pd.Series(index=["null"]),
              "s7": pd.Series(index=["null"]),              
          })      

        ## store cluster depth information (biggest memory cost)
        self.depths = []


    def __str__(self):
        return "<ipyrad.Sample object {}>".format(self.name)

    def to_JSON(self):
        """ write to JSON to serialize object """
        return json.dumps(self, default=lambda o: o.__dict__, 
                          sort_keys=True, indent=4)

    #def save(self):
    #    """ pickle the data object """
    #    dillout = open(self.name+".dataobj", "wb")
    #    dill.dump(self, dillout)
    #    dillout.close()


