#!/usr/bin/env python

"""
Install ipyrad with:
 `conda install ipyrad -c conda-forge -c bioconda`

Or, for developers, install dependencies and then ipyrad w/ pip local:
 `cd ipyrad/`
 `conda env create -f environment.yml`
 `pip install . --no-deps`
"""

import re
from setuptools import setup, find_packages


# Auto-update ipyrad version from git repo tag
# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
INITFILE = "ipyrad/__init__.py"
CUR_VERSION = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                        open(INITFILE, "r").read(),
                        re.M).group(1)

setup(
    name="ipyrad",
    version=CUR_VERSION,
    url="https://github.com/dereneaton/ipyrad",
    author="Deren Eaton and Isaac Overcast",
    author_email="de2356@columbia.edu",
    description="Interactive assembly and analysis of RAD-seq data sets",
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    packages=find_packages(),
    install_requires=[
        "future",
        "ipyparallel",
        "scipy",
        "numpy",
        "numba",
        "pandas",
        "pysam",
        "h5py",
        "requests",
        "loguru",
    ],
    entry_points={'console_scripts': ['ipyrad = ipyrad.__main__:CLI']},
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
