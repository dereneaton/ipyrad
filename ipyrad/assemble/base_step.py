#!/usr/bin/env python

"""Abstract Base class for Step class objects called by Assembly.run().

Each step of the assembly creates a subclass, e.g., Step1, Step2, etc,
that inherit from the Base class. This base class prints a header, 
fetches the subsamples ready to process, and parses population data.
"""

from typing import TypeVar, Dict, Tuple, List
from pathlib import Path
import shutil
from abc import ABC
from loguru import logger
from ipyrad.assemble.utils import IPyradExit
from ipyrad.core.schema import SampleSchema

Assembly = TypeVar("Assembly")
Sample = TypeVar("Sample")
logger = logger.bind(name="ipyrad")


class BaseStep(ABC):
    """Abstract Base Class for Step class objects.

    This class must be subclassed to be used. Subclasses exist for each
    step of the ipyrad assembly process. This includes functions run
    by all steps, such as printing (logging) header info, getting 
    the samples, ensuring directories, and optionally loading popfile.
    """
    def __init__(self, data: Assembly, step: str, quiet: bool, force: bool):

        # store step
        self.data = data
        self.step = step
        self.quiet = quiet
        self.force = force

        # data to be filled, updates settings of Assembly
        self.samples: Dict[str, Sample] = {}
        """: Subsample of Assembly Samples ready for this step."""
        self.data.stepdir: Path = None
        """: Output dir used during this step."""
        self.data.tmpdir: Path = None
        """: Tmp dir used during this step."""
        self.data.populations: Dict[str, Tuple[List[str], int]] = {}
        """: Population mapping {popname: ([sample_names], minsamp)}."""
        self.data.outfiles: Dict[str, Path] = {}
        """: Output files."""  # reset's to empty on any run() step.

        # Initiate the step.
        self._print_headers()
        self._get_subsamples()
        self._setup_dirs()
        if self.step == 7:
            self._parse_populations()

    def _print_headers(self) -> None:
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
                logger.info(messages['1'])
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

    def _get_subsamples(self) -> None:
        """Sets step.samples dictionary mapping {sname: Sample}

        Select samples ready for this step based on the .state 
        attribute. This is set at the end of each step if a sample
        completes successfully. If, for example, a sample recovered
        zero clusters at the end of step3 it will not advance to 
        state=3, and will be skipped when running step 4.

        The 'reference' sample is created in step7 if not 
        hackers.exclude_reference. This sample is dropped if running
        any steps<7.
        """
        # no samples to get yet
        if self.step == 1:
            return
        if not self.data.samples:
            raise IPyradExit("Error: No samples found. You must first run step 1.")

        # check samples for current states
        not_ready = {}
        already_done = {}
        todo_samples = {}
        for sname, sample in self.data.samples.items():
            # drop the reference sample created in s7, if present
            if sname != "reference":
                if sample.state < self.step - 1:
                    not_ready[sname] = sample
                elif sample.state >= self.step:
                    already_done[sname] = sample
                else:
                    todo_samples[sname] = sample

        # warn about skipped samples that are not ready
        if not_ready:
            if len(not_ready) == len(self.data.samples):
                raise IPyradExit(f"Error: No samples ready for step {self.step}")
            logger.warning(
                f"Skipping samples not ready for step {self.step}. Create "
                "a new branch and drop these samples \nto suppress this "
                f"warning message: {', '.join(list(not_ready))}"
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

        # STEP7 only: add a 'reference' sample (even if we will exclude it later.)
        if (self.step == 7) and self.data.is_ref:
            self.samples['reference'] = SampleSchema(name="reference", state=7)

    def _setup_dirs(self) -> None:
        """Create the output dir and tmpdir for this step."""
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
        self.data.stepdir = (
            self.data.params.project_dir /
            f"{self.data.name}_{suffix[self.step]}")
        self.data.tmpdir = (
            self.data.params.project_dir /
            f"{self.data.name}_tmp_{suffix[self.step]}")

        # clear stepdir or raise an error depending on exists and force
        # for steps 2-7 we can go ahead and remove the entire directory
        # that we created, since it should not contain any other files
        # unless the user added them. *They did use the force flag*.
        if self.data.stepdir.exists():
            if self.force:
                msg = f"removing previous {suffix[self.step]} dir: {self.data.stepdir}"
                logger.info(msg)
                shutil.rmtree(self.data.stepdir)
            else:
                msg = f"Error: Directory {self.data.stepdir} exists.\nUse force (-f) to overwrite."
                raise IPyradExit(msg)

        # always clear tmpdir, and always make both new dirs.
        if self.data.tmpdir.exists():
            shutil.rmtree(self.data.tmpdir)
        self.data.tmpdir.mkdir(exist_ok=True)
        self.data.stepdir.mkdir(exist_ok=True)

    def _parse_populations(self) -> None:
        """Parse the population file params input file. 

        In the API a user can set the populations using a dictionary 
        on the Assembly object as {popname: ([samps], minsamp), ...}.

        In the CLI this information is parsed from a file.
        Default format is one sample per line, listing population name
        then sample name separated by whitespace. In addition, 'minsamp'
        population minimum coverage settings (analagous to the param
        `minsamples_locus` but applied to each population) must be set
        on an additional line starting with a hash character.

        >>> ind1 pop1
        >>> ind2 pop1
        >>> ind3 pop2
        >>> ind4 pop2
        >>> ...        
        >>> # pop1:2, pop2: 2, ...
        """
        # if Assembly already has population info saved to it (e.g., 
        # in its JSON) then if pop_assign_file is empty we do not 
        # want it to be removed... so hard to do both API and CLI...

        # TODO: just didn't get around to re-implementing and testing yet...
        if self.data.params.pop_assign_file:
            raise NotImplementedError(
                "populations assignments during assembly are not "
                "currently supported in this version. You can however "
                "filter by population sampling using ipyrad-analysis "
                "in the ipyrad API after finishing assembly."
            )




class TestStep(BaseStep):
    """Test BaseStep with a toy example subclass."""
    def __init__(self, data, force, quiet, ipyclient):
        super().__init__(data=data, step=1, force=force, quiet=quiet)
        self.ipyclient = ipyclient



if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")
    DATA = ip.Assembly("TEST")
    DATA.params.project_dir = "/tmp"
    STEP = TestStep(DATA, force=True, quiet=False, ipyclient=None)
    print(STEP.samples)
