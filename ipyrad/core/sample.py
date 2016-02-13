#!/usr/bin/env python2.7
""" Sample object """

#import os
import pandas as pd
import numpy as np
import json
from collections import OrderedDict
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
                   "reads_consens",]).astype(np.object)

        ## stats for each step
        self.statsfiles = ObjDict({
              "s1": pd.Series(index=["reads_raw",
                                     ]).astype(np.object),

              "s2": pd.Series(index=["reads_raw",
                                     "filtered_by_qscore",
                                     "filtered_by_adapter",
                                     "reads_passed",
                                     ]).astype(np.object),

              "s3": pd.Series(index=["merged_pairs",
                                     "clusters_total",
                                     "clusters_hidepth",
                                     "avg_depth_tot",
                                     "avg_depth_mj",
                                     "avg_depth_stat",
                                     "sd_depth_tot", 
                                     "sd_depth_mj",
                                     "sd_depth_stat",
                                     "filtered_bad_align",
                                     ]).astype(np.object),

              "s4": pd.Series(index=["hetero_est",
                                     "error_est",
                                     ]).astype(np.object),

              "s5": pd.Series(index=["nclusters",
                                     "depthfilter",
                                     "maxHfilter",
                                     "maxAllelefilter",
                                     "maxNfilter",
                                     "nconsens",
                                     "nsites",
                                     "nhetero",
                                     "heterozygosity",
                                     ]).astype(np.object)
              })

        ## store cluster depth information (biggest memory cost), 
        ## should save this as a Counter dict, currently its a list
        self.depths = {}


    def __str__(self):
        return "<ipyrad.Sample object {}>".format(self.name)


    def _to_fulldict(self):
        """ Write to dict including data frames. All sample dicts 
        are combined in save() to dump JSON output """
        
        ## 
        returnjson = json.dumps({
            "name": self.name,
            "barcode": self.barcode,
            "files": self.files,
            "statsfiles": {
                "s1": self.statsfiles.s1.to_dict(),
                "s2": self.statsfiles.s2.to_dict(),                
                "s3": self.statsfiles.s3.to_dict(),
                "s4": self.statsfiles.s4.to_dict(),
                "s5": self.statsfiles.s5.to_dict(),
            },
            "stats": self.stats.to_dict()
            },
            default=lambda out: out.__dict__, 
            sort_keys=True,
            indent=4, 
            separators=(",", ";"))

        returndict = OrderedDict([
            ("name", self.name),
            ("barcode", self.barcode),
            ("files", self.files),
            ("statsfiles", {
                "s1": self.statsfiles.s1.to_dict(),
                "s2": self.statsfiles.s2.to_dict(),                
                "s3": self.statsfiles.s3.to_dict(),
                "s4": self.statsfiles.s4.to_dict(),
                "s5": self.statsfiles.s5.to_dict(),
            }),
            ("stats", self.stats.to_dict()),
            ])

        return returndict


    #def save(self):
    #    """ pickle the data object """
    #    dillout = open(self.name+".dataobj", "wb")
    #    dill.dump(self, dillout)
    #    dillout.close()


