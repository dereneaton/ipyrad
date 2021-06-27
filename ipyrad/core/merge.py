#!/usr/bin/env python

"""
Merge assemblies
"""

from typing import List, Dict, Optional
from loguru import logger
from ipyrad.core.params_schema import ParamsSchema
from ipyrad.core.schema import SampleSchema
from ipyrad.core.assembly import Assembly
from ipyrad.assemble.utils import BADCHARS


def merge(
    name: str, 
    assemblies: List[Assembly], 
    rename_dict: Optional[Dict[str,str]] = None,
    ):
    """
    Creates and return a new Assembly containing the merged samples
    of all input Assemblies, and parameter settings set to the 
    first Assembly. Merging can be used to combine samples up to 
    step 5 (i.e., before running 6). Sample states will be set to 
    a maximum of 5 or their current state. Merging does not affect 
    the actual files that currently exist, but rather creates new 
    samples that reference multiple existing files.

    Examples:
    ---------
    # merge two assemblies
    new = ip.merge('newname', [assembly1, assembly2])

    # merge two assemblies and rename samples
    rename = {"1A_0", "A", "1B_0", "A"}
    new = ip.merge('newname', [assembly1, assembly2], rename_dict=rename)    
    """
    # update rename dict to avoid bad characters
    rename_dict = rename_dict if rename_dict is not None else {}
    for oldname in rename_dict:
        newname = rename_dict[oldname]
        newername = "".join([
            i.replace(i, "_") if i in BADCHARS else i for i in newname
        ])
        if any(i in newname for i in BADCHARS):
            logger.warning(
                "modifying {newname} to {newername} to avoid bad characters.")
        rename_dict[oldname] = newername

    # ensure Assemblies is multiple
    assert len(assemblies) >= 2, "must enter >1 assembly to be merged"

    # create new Merged assembly that inherits params from the 1st assembly
    merged = Assembly(name)
    params = assemblies[0].params.dict()
    params['assembly_name'] = merged.name
    merged.params = ParamsSchema(**params)

    # copy samples into new. For samples that are present in multiple
    # assemblies this requires merging their stats from previous steps.

    # A flag to set if there are technical replicates among merging
    # assemblies, so we can print a helpful message.
    any_replicates = False

    # iterate over all sample names from all Assemblies
    for data in assemblies:

        # make a deepcopy of the sample
        for sname in data.samples:
            sample = SampleSchema(**data.samples[sname].dict())

            # rename sample if in rename dict and update sname variable
            if sname in rename_dict:
                sname = rename_dict[sname]
                sample.name = sname

            # is it in the merged assembly already (technical replicate)
            if sname in merged.samples:
                msample = merged.samples[sname]

                # update stats for steps 1-2
                msample.stats_s1.reads_raw += sample.stats_s1.reads_raw
                # if both are not step >=2 then set to state=1
                if msample.stats_s2 and sample.stats_s2:
                    msample.stats_s2.reads_passed_filter += (
                        sample.stats_s2.reads_passed_filter)
                else:
                    msample.state = 1

                # append files
                if sample.files.fastqs:
                    msample.files.fastqs.extend(sample.files.fastqs)
                if sample.files.edits:
                    msample.files.edits.extend(sample.files.edits)

                # do not allow state >2 at merging (requires reclustering)
                # if merging WITHIN samples. Set the flag so we can inform
                # the user after all the samples have been handled
                if sample.state > 2:
                    msample.state = 2
                    any_replicates = True

            # merge its stats and files
            else:
                merged.samples[sname] = sample

    # set these to empty sinc they no longer relevant
    # merged_names = ", ".join([i.name for i in assemblies])
    merged.params.raw_fastq_path = None # "Merged: " + merged_names
    merged.params.barcodes_path = None # "Merged: " + merged_names
    merged.params.sorted_fastq_path = None # "Merged: " + merged_names

    if any_replicates:
        logger.warning(MERGED_TECHNICAL_REPLICATES)

    # return the new Assembly object (save to JSON also occurs in CLI command.)
    return merged


MERGED_TECHNICAL_REPLICATES = """\
    NB: One or more samples are present in one or more of the merged 
    assemblies, and are beyond step 3. Technical replicates need to be 
    clustered within samples so YOU MUST re-run these samples from at 
    least step 3. Sample states in the new Assembly are set to 2 to 
    enforce this.
    """


if __name__ == "__main__":

    import ipyrad as ip
    data1 = ip.Assembly("test1")
    data2 = ip.Assembly("test2")
    data3 = ip.merge('merged', assemblies=[data1, data2])
    data3.stats
