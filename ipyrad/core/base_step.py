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
from ipyrad.schema import Sample
from ipyrad.schema.sample_schema import (
    Stats2, Stats3, Stats4, Stats5, Stats6, Stats7
)
from ipyrad.core.exceptions import IPyradError

Assembly = TypeVar("Assembly")
logger = logger.bind(name="ipyrad")


class BaseStep(ABC):
    """Abstract Base Class for Step class objects.

    This class must be subclassed to be used. Subclasses exist for each
    step of the ipyrad assembly process. This includes functions run
    by all steps, such as printing (logging) header info, getting
    the samples, ensuring directories, and optionally loading popfile.
    """
    def __init__(self, data: Assembly, step: str, force: bool):

        # store step
        self.data = data
        self.step = step
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
        self._log_headers()
        self._get_subsamples()
        self._add_sample_stats()
        self._setup_dirs()
        if self.step == 7:
            self._parse_populations()

    def _log_headers(self) -> None:
        """print the CLI header for this step"""
        messages = {
            '1': "Step 1: Loading/trimming demultiplexed fastq data files",
            '2': "Step 2: Clustering/Mapping reads within samples",
            '3': "Step 3: Building/Aligning clusters within samples",
            '4': "Step 4: Calling/Filtering consensus alleles within samples",
            '5': "Step 5: Clustering/Mapping orthologs across samples",
            '6': "Step 6: Building/Aligning loci across samples",
            '7': "Step 7: Trimming/Filtering and writing output files",
        }
        key = str(self.step)
        logger.info(messages[key])

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
            raise IPyradError("Error: No samples found. You must first run step 1.")

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
                raise IPyradError(f"Error: No samples ready for step {self.step}")
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
            self.samples['reference'] = Sample(name="reference", state=7)

    def _add_sample_stats(self) -> None:
        """Add/Reset stats_sX object for each sample."""
        for sname, sample in self.samples.items():
            # do not change any steps below current.
            # set current step to new Obj
            # set future steps to None
            for s in range(2, 8):
                if s == self.step:
                    setattr(sample, f"stats_s{s}", get_stat_obj(s)())
                if s > self.step:
                    setattr(sample, f"stats_s{s}", None)

    def _setup_dirs(self) -> None:
        """Create the output dir and tmpdir for this step."""
        suffix = {
            1: "trimmed_fastqs",
            2: "clustmaps_within",
            3: "clusters_within",
            4: "consensus",
            5: "clustmaps_across",
            6: "clusters_across",
            7: "outfiles",
        }
        self.data.stepdir = (
            self.data.params.project_dir / f"{self.data.name}_{suffix[self.step]}")
        self.data.tmpdir = (
            self.data.params.project_dir / f"{self.data.name}_tmp_{suffix[self.step]}")

        # clear stepdir or raise an error depending on exists and force
        # for steps 2-7 we can go ahead and remove the entire directory
        # that we created, since it should not contain any other files
        # unless the user added them. *They did use the force flag*.
        if self.data.stepdir.exists():
            if self.force:
                msg = f"removing old step {self.step} dir: {self.data.stepdir}"
                logger.info(msg)
                shutil.rmtree(self.data.stepdir)
            else:
                msg = f"Error: Directory {self.data.stepdir} exists.\nUse force (-f) to overwrite."
                logger.error(msg)
                raise IPyradError(msg)

        # always clear tmpdir, and always make both new dirs.
        logger.debug(f"creating step {self.step} dir: {self.data.stepdir}")
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


def get_stat_obj(step: int) -> "StatsX":
    if step == 2:
        return Stats2
    if step == 3:
        return Stats3
    if step == 4:
        return Stats4
    if step == 5:
        return Stats5
    if step == 6:
        return Stats6
    if step == 7:
        return Stats7


class TestStep(BaseStep):
    """Test BaseStep with a toy example subclass."""
    def __init__(self, data, force, ipyclient):
        super().__init__(data=data, step=1, force=force)
        self.ipyclient = ipyclient


if __name__ == "__main__":

    import ipyrad as ip
    ip.set_log_level("DEBUG")
    DATA = ip.Assembly("TEST")
    DATA.params.project_dir = "/tmp"
    STEP = TestStep(DATA, force=True, ipyclient=None)
    print(STEP.samples)