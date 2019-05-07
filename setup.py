#!/usr/bin/env python

from setuptools import setup, find_packages
import glob
import re


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
    packages=find_packages(),
    install_requires=[
        "future",
        "notebook",
        "ipyparallel>=6.0.2",
        "scipy>0.10",
        "numpy>=1.9",
        "numba>=0.37",
        "pandas>=0.16",
        "pysam>=0.10.0",  # ipyrad::
        "h5py",
        "mpi4py",
        "cutadapt",    # ipyrad::
        # "toytree",     # eaton-lab::
        "requests",
    ],
    entry_points={
        'console_scripts': [
            'ipyrad = ipyrad.__main__:main',
        ],
    },
    data_files=[('bin', glob.glob("./bin/*"))],
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
)
