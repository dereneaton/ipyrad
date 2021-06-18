#!/usr/bin/env python

"""
General Base class for Step class objects called by Assembly.run()
to call the basic setup funcs shared by all steps.
"""

import os
import shutil
from loguru import logger
from ipyrad.assemble.utils import IPyradError


class BaseStep:
    """
    General Base class for Step class objects called by Assembly.run()
    to call the basic setup funcs shared by all steps.
    """
    def __init__(self, data, step, quiet, force):

        # store step
        self.data = data
        self.quiet = quiet
        self.step = step
        self.force = force
        self.print_headers()

        # will be filled by .get_subsamples() 
        self.samples = {}
        self.get_subsamples()

        # dirs will be created for this step by .setup_dirs()
        self.stepdir = None
        self.tmpdir = None
        self.setup_dirs()

        # parse population file for step 7
        if step == 7:
            self.parse_populations()


    def print_headers(self) -> None:
        """print the CLI header for this step"""
        if self.quiet:
            return
        messages = {
            '1': "Step 1: Demultiplexing fastq data to samples",
            '1a': "Step 1: Loading demultiplexed fastq data files",
            '2': "Step 2: Adapter trimming and filtering reads",
            '3': "Step 3: Denovo clustering reads within samples",
            '3a': "Step 3: Mapping reads to reference within samples",
            '4': "Step 4: Jointly estimating error rate and heterozygosity",
            '5': "Step 5: Calling consensus alleles within samples",
            '6': "Step 6: Denovo clustering homologs across samples",
            '6a': "Step 6: Reference aligning homologs across samples", 
            '7': "Step 7: Filtering, trimming, and writing output files",
        }
        key = str(self.step)
        if key == '1':
            if self.data.params.sorted_fastq_path is not None:
                print(messages['1a'])
            else:
                print(messages['1'])
        elif key == '3':
            if self.data.params.assembly_method == "reference":
                print(messages['3a'])
            else:
                print(messages['3'])
        elif key == "6":
            if self.data.params.assembly_method == "reference":
                print(messages['6a'])
            else:
                print(messages['6'])
        else:
            print(messages[key])


    def get_subsamples(self) -> None:
        """
        Sets step.samples dictionary mapping {sname: Sample}

        Select samples ready for this step based on the .state 
        attribute. This is set at the end of each step if a sample
        completes successfully. If, for example, a sample recovered
        zero clusters at the end of step3 it will not advance to 
        state=3, and will be skipped when runnign step 4.
        """
        # no samples to get yet
        if self.step == 1:
            return

        # check samples for current states
        not_ready = {}
        already_done = {}
        todo_samples = {}
        for sname in self.data.samples:
            if self.data.samples[sname].state < self.step - 1:
                not_ready[sname] = self.data.samples[sname]
            elif self.data.samples[sname].state >= self.step:
                already_done[sname] = self.data.samples[sname]
            else:
                todo_samples[sname] = self.data.samples[sname]

        # warn about skipped samples that are not ready
        if not_ready:
            if len(not_ready) == len(self.data.samples):
                raise IPyradError(f"No samples ready for step {self.step}")
            logger.warning(
                f"Skipping samples not ready for step {self.step}. Create "
                "a new branch and drop these samples to suppress this "
                "warning message:\n"
                f"{', '.join(list(already_done))}"
            )

        # build list to run for samples being forced
        if self.force:
            self.samples = {**already_done, **todo_samples}
        else:
            if already_done:
                logger.warning(
                    f"skipping samples already finished step {self.step}:\n"
                    f"{', '.join(list(already_done))}")
            self.samples = todo_samples


    def setup_dirs(self) -> None:
        """
        Create the output dir and tmpdir for this step.
        """
        logger.debug(f"creating dirs for step {self.step}")
        suffix = {
            1: "fastqs",
            2: "edits",
            3: "clustmap",
            4: "estimates",
            5: "consensus",
            6: "across",
            7: "outfiles",
        }
        self.stepdir = os.path.join(
            self.data.params.project_dir, 
            self.data.name + f"_{suffix[self.step]}",
        )
        self.tmpdir = os.path.join(
            self.data.params.project_dir, 
            self.data.name + f"_tmp_{suffix[self.step]}",
        )

        # clear stepdir or raise an error depending on exists and force
        if os.path.exists(self.stepdir):            
            if self.force:
                logger.debug(
                    f"removing previous {suffix[self.step]} dir: {self.stepdir}"
                )
                shutil.rmtree(self.stepdir)
            else:
                raise IPyradError(
                    f"Directory {self.stepdir} already exists. "
                    " use force to overwrite it.")
            
        # always clear tmpdir
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)

        # always make both dirs
        os.makedirs(self.tmpdir, exist_ok=True)
        os.makedirs(self.stepdir, exist_ok=True)


    def parse_populations(self) -> None:
        """
        Parse the population file params input file. In the API
        a user _could_ alternatively add a pop dictionary to the 
        Assembly object as {popname: ([samps], minsamp), ...}
        """
        # TODO:



class TestStep(BaseStep):
    """
    Testing example setup for step subclasses
    """
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data=data, step=1, force=force, quiet=quiet)
        self.ipyclient = ipyclient



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_loglevel("DEBUG")
    DATA = ip.Assembly("TEST")
    DATA.params.project_dir = "/tmp"
    STEP = TestStep(DATA, force=True, quiet=False, ipyclient=None)
    print(STEP)
