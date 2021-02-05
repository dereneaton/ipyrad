#!/usr/bin/env python

""" Sample object """

from collections import OrderedDict
import pandas as pd
import numpy as np
from ipyrad.assemble.utils import ObjDict



class Sample(object):
    """ 
    ipyrad Sample object. Links to files associated with an individual 
    sample, used to combine samples into Assembly objects.
    """

    def __init__(self, name=""):
        self.name = name 
        self.barcode = ""

        # link to files
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
                   "reads_passed_filter",
                   "reads_merged",
                   "refseq_mapped_reads",
                   "refseq_unmapped_reads",
                   "clusters_total",
                   "clusters_hidepth",
                   "hetero_est",
                   "error_est",
                   "reads_consens",
                   ]).astype(np.object)

        ## stats for each step
        self.stats_dfs = ObjDict({
              "s1": pd.Series(index=["reads_raw",
                                     ]).astype(np.object),

              "s2": pd.Series(index=["reads_raw",
                                     # "trim_adapter_bp_read1",
                                     # "trim_adapter_bp_read2",                                     
                                     # "trim_quality_bp_read1",
                                     # "trim_quality_bp_read2",                                     
                                     "reads_filtered_by_Ns",
                                     "reads_filtered_by_low_quality",
                                     "reads_filtered_by_low_complexity",                                     
                                     "reads_filtered_by_minlen",
                                     "avg_len_R1_before_trimming",
                                     "avg_len_R1_after_trimming",                                     
                                     "avg_len_R2_before_trimming",
                                     "avg_len_R2_after_trimming",
                                     "reads_passed_filter",
                                     ]).astype(np.object),
                                     #"filtered_by_qscore",
                                     #"filtered_by_adapter",

              "s3": pd.Series(index=["merged_pairs",
                                     "clusters_total",
                                     "hidepth_min",
                                     "clusters_hidepth",
                                     "avg_depth_total",
                                     "avg_depth_mj",
                                     "avg_depth_stat",
                                     "sd_depth_total", 
                                     "sd_depth_mj",
                                     "sd_depth_stat",
                                     "filtered_bad_align",
                                     ]).astype(np.object),

              "s4": pd.Series(index=["hetero_est",
                                     "error_est",
                                     ]).astype(np.object),

              "s5": pd.Series(index=["clusters_total",
                                     "filtered_by_depth",
                                     "filtered_by_maxH",
                                     "filtered_by_maxAlleles",
                                     "filtered_by_maxN",
                                     "reads_consens",
                                     "nsites",
                                     "nhetero",
                                     "heterozygosity",
                                     ]).astype(np.object)
              })

        ## store cluster depth information (biggest memory cost), 
        self.depths = {}


    def __str__(self):
        return "<ipyrad.Sample object {}>".format(self.name)


    def _to_fulldict(self):
        """ 
        Write to dict including data frames. All sample dicts 
        are combined in save() to dump JSON output """
        
        ## 
        returndict = OrderedDict([
            ("name", self.name),
            ("barcode", self.barcode),
            ("files", self.files),
            ("stats_dfs", {
                "s1": self.stats_dfs.s1.to_dict(),
                "s2": self.stats_dfs.s2.to_dict(),                
                "s3": self.stats_dfs.s3.to_dict(),
                "s4": self.stats_dfs.s4.to_dict(),
                "s5": self.stats_dfs.s5.to_dict(),
            }),
            ("stats", self.stats.to_dict()),
            ("depths", self.depths)
            ])

        return returndict
