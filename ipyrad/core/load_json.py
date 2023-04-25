#!/usr/bin/env python

"""
Load an Assembly object from a project JSON file.
"""

from loguru import logger
from ipyrad.schema import Project
from ipyrad.core import Assembly


def load_json(json_file: str) -> Assembly:
    """Return an Assembly object loaded from a project JSON file.
    """
    proj = Project.parse_file(json_file)
    data = Assembly(proj.params.assembly_name)
    data.samples = {i: j for (i, j) in proj.samples.items() if i != "reference"}
    data.params = proj.params
    data.outfiles = proj.outfiles
    data.populations = proj.populations
    data.assembly_stats = proj.assembly_stats
    logger.info(f"loaded Assembly {data.name} from JSON file")
    return data


if __name__ == "__main__":

    DATA = load_json("/tmp/hi.json")
