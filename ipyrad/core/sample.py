#!/usr/bin/env python2.7
""" Sample object """

#import os
import pandas as pd
from ipyrad.assemble.worker import ObjDict

# pylint: disable=C0103

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

        ## step stats files (builds into data.statsfiles.sX)
        self.s1 = ObjDict({"reads_raw": []})

        self.s2 = ObjDict({"reads_raw": [],
                           "filter_qscore": [],
                           "filter_adapter": [],
                           "reads_passed": []})

        self.s3 = ObjDict({"reads_raw": [],
                           "clusts_total": [],
                           "clusts_hidepth": [],
                           "avg.depth.tot": [],
                           "avg.depth>mj": [],
                           "avg.depth>stat": []}) 

        self.s4 = ObjDict({"hetero_est": [],
                           "error_est": []})

        self.s5 = ObjDict({"nclusters": [],
                           "depthfilter": [],
                           "maxHfilter": [],
                           "maxAllelefilter": [],
                           "maxNfilter": [],
                           "nconsens": [],
                           "nsites": [],
                           "nhetero": [],
                           "heterozygosity": []})

        self.s6 = ObjDict({})
        self.s7 = ObjDict({})        

        ## store cluster depth information (biggest memory cost)
        self.depths = []


    #def __getstate__(self):
    #    return self.__dict__

    #def __setstate__(self, dicto):
    #    self.__dict__.update(dicto)

    #def save(self):
    #    """ pickle the data object """
    #    dillout = open(self.name+".dataobj", "wb")
    #    dill.dump(self, dillout)
    #    dillout.close()


