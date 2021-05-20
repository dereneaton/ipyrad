#!/usr/bin/env python

"""
Assembly class object as the main API for calling assembly steps.
"""

import os
import sys
import shutil
from typing import List, Optional
from loguru import logger
import numpy as np
import pandas as pd
from ipyparallel import Client
from ipyrad.core.params_schema import ParamsSchema, HackersSchema
from ipyrad.core.schema import Project, SampleSchema
from ipyrad.core.parallel import Cluster
from ipyrad.assemble.utils import IPyradError
from ipyrad.assemble.s1_demux import Step1
from ipyrad.assemble.s2_trim_reads import Step2
from ipyrad.assemble.s3_clustmap_within import Step3
from ipyrad.assemble.s4_joint_estimate import Step4
from ipyrad.assemble.s5_consensus import Step5
from ipyrad.assemble.s6_clustmap_across import Step6
from ipyrad.assemble.s7_assemble import Step7


class Assembly:
    """
    Returns a new Assembly class instance with default parameter
    settings. The first steps of an API analysis typically involves
    initializing a new Assembly, setting params, and then calling
    self.run() to start running assembly steps.
    """
    def __init__(self, name:str):

        # core JSON file components
        self.params = ParamsSchema(assembly_name=name)
        self.hackers = HackersSchema()
        self.samples = {}

        # optional dict for setting cluster config.
        self.ipcluster = {
            "cores": 0,
            "threads": 2,
        }

    def __repr__(self):
        return "<ipyrad.Assembly object {}>".format(self.name)

    @property
    def name(self):
        """shortcut to return the assembly_name from params"""
        return self.params.assembly_name
    
    @property
    def json_file(self):
        """JSON file is project_dir/assembly_name.json"""
        return os.path.join(
            self.params.project_dir, 
            self.params.assembly_name + ".json"
        )

    @property
    def stats(self):
        """
        Returns a dataframe with *summarized stats* extracted from the 
        project JSON file given the current completed assembly steps.
        To see more detailed stats from specific steps see instead 
        the self.stats_dfs 
        """
        proj = Project.parse_file(self.json_file)
        stats = pd.DataFrame(
            index=sorted(proj.samples),
            columns=[
                'state', 
                'reads_raw', 
                'reads_passed_filter',
                'clusters_total',
                'clusters_hidepth',
                'mean_depth_total',
                'reads_mapped_to_ref_prop',
                'consensus_total',
                'heterozygosity',
            ],
        )
        for sname in stats.index:
            sample = proj.samples[sname]
            stats.loc[sname, 'state'] = sample.state
            stats.loc[sname, 'reads_raw'] = sample.stats_s1.reads_raw

            if sample.stats_s2:
                value = sample.stats_s2.reads_passed_filter
            else:
                value = pd.NA
            stats.loc[sname, 'reads_passed_filter'] = value

            if sample.stats_s3:
                value = sample.stats_s3.clusters_total
            else:
                value = pd.NA
            stats.loc[sname, 'clusters_total'] = value

            if sample.stats_s3:
                value = sample.stats_s3.clusters_hidepth
            else:
                value = pd.NA
            stats.loc[sname, 'clusters_hidepth'] = value

            if sample.stats_s3:
                value = sample.stats_s3.mean_depth_total
            else:
                value = pd.NA
            stats.loc[sname, 'mean_depth_total'] = value

            if sample.stats_s3:
                value = sample.stats_s3.reads_mapped_to_ref_prop
            else:
                value = pd.NA
            stats.loc[sname, 'reads_mapped_to_ref_prop'] = value

            if sample.stats_s5:
                value = sample.stats_s5.consensus_total
            else:
                value = pd.NA
            stats.loc[sname, 'consensus_total'] = value            

            if sample.stats_s5:
                value = sample.stats_s5.heterozygosity
            else:
                value = pd.NA
            stats.loc[sname, 'heterozygosity'] = value            

        # drop columns that are all NAN
        stats = stats.dropna(axis=1, how="all")

        # FIXME: keep adding summary stats
        return stats
    

    def branch(self, name:str, subsample:List[str]=None) -> 'Assembly':
        """
        Returns a new branched Assembly class instance.

        The new object will have the same parameter settings as the 
        current object, and inherits the same sample histories, but 
        will write to a different name prefix path
        going forward. 

        Creating a new branch does not write a JSON until you either 
        run an Assembly step by calling .run() or call .save_json()

        Examples:
        ---------
        data1 = ip.Assembly("data1")
        data1.params.sorted_fastq_path = "./data/*.gz"
        data1.run('1')
        data2 = data1.branch("data2", subsample=['1A_0', '1B_0'])
        data2.run('2')
        """
        # create new names Assembly and copy over all params except name
        branch = Assembly(name)
        params = self.params.dict()
        params['assembly_name'] = name
        branch.params = ParamsSchema(**params)
        branch.hackers = HackersSchema(**self.hackers.dict())

        # copy over all or just a subsamples of the samples.
        if subsample is None:
            branch.samples = {
                i: SampleSchema(**self.samples[i].dict()) for i in self.samples
            }
        else:
            branch.samples = {
                i: SampleSchema(**self.samples[i].dict()) 
                for i in self.samples if i in subsample
            }
            for i in subsample:
                if i not in self.samples:
                    logger.warning(f"sample name {i} does not exist.")
        return branch


    def write_params(self) -> None:
        """
        Write a CLI params file from the current params settings.
        """
        # TODO


    def save_json(self) -> None:
        """
        Save the current Assembly object to the project JSON file  
        (<project_dir>/<name>.json)
        """
        project = Project(
            params=ParamsSchema(**self.params.dict()),
            hackers=HackersSchema(**self.hackers.dict()),
            samples={sname: self.samples[sname] for sname in self.samples},
        )
        with open(self.json_file, 'w') as out:
            out.write(project.json(indent=4, exclude_none=True))
        logger.debug(f"wrote to {self.json_file}")


    def run(
        self, 
        steps:str, 
        force:bool=False, 
        quiet:bool=False,
        ipyclient:Optional[Client]=None,
        **kwargs,
        ) -> None:
        """
        Run one or more assembly steps (1-7) of an ipyrad assembly.

        Parameters
        ----------
        steps: str
            A string of steps to run, e.g., "1", or "123".
        force: bool
            Force overwrite of existing results for this step.
        quiet: bool
            Suppress printed outputs to stdout.
        ipyclient: Optional[ipyparallel.Client]
            Optional ipyparallel client to connect to for distributing
            jobs in parallel. This option is generally only useful if
            you start a Client using MPI to connect to multiple nodes
            of an HPC cluster. Otherwise, just configure the local 
            cluster parallelization using the .ipcluster attribute.
        """
        # save a tmp backup copy of the JSON file.
        # FIXME ...

        # ...
        step_map = {
            "1": Step1,
            "2": Step2, 
            "3": Step3,
            "4": Step4,
            "5": Step5, 
            "6": Step6,
            "7": Step7,
        }

        # could load the tool to check whether this job can be run 
        # before starting the ipcluster?...


        # start the cluster
        cluster = Cluster()
        try:
            # establish connection to a running ipyclient
            cores = kwargs.get('cores', self.ipcluster['cores'])
            cluster.start(cores=cores, ipyclient=ipyclient)

            # use client for all steps of assembly 
            for step in steps:
                tool = step_map[step](self, force, quiet, cluster.ipyclient)
                tool.run()
                shutil.rmtree(tool.tmpdir)

        except KeyboardInterrupt:
            logger.warning("keyboard interrupt by user, cleaning up.")

        # AssemblyProgressBar logs the traceback
        except IPyradError as inst:
            logger.error(f"An error occurred:\n{inst}")
            print("An error occurred, see logfile and below.")
            raise

        # logger.exception logs the traceback
        except Exception as inst:
            logger.error(f"An unexpected error occurred:\n{inst}")
            print("An unexpected error occurred, see logfile and below.")
            raise

        finally:
            cluster.cleanup_safely(None)



if __name__ == "__main__":
    

    import ipyrad as ip
    ip.set_loglevel("DEBUG", stderr=False, logfile="/tmp/test.log")

    # TEST = ip.Assembly("PEDIC")
    # TEST.params.sorted_fastq_path = "../../sra-fastqs/*.fastq"
    # TEST.params.project_dir = "/tmp"
    # TEST.run('12', force=True, quiet=True)
    # print(TEST.stats)

    TEST = ip.Assembly("TEST1")
    TEST.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_R1*.gz"    
    TEST.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes.txt"
    TEST.params.project_dir = "/tmp"
    TEST.params.max_barcode_mismatch = 1
    TEST.run('1', force=True, quiet=True)

    data = ip.Assembly('TEST')
    data.params.project_dir = "/tmp"
    data.params.raw_fastq_path = "../../tests/ipsimdata/rad_example_R1*.fastq.gz"
    data.params.barcodes_path = "../../tests/ipsimdata/rad_example_barcodes.txt"
    data.run("1", force=True, quiet=True)
    print(data.stats)

