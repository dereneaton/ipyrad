#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import glob
import re


def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()


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
    author="Deren Eaton",
    author_email="de2356@columbia.edu",
    description="Interactive assembly and analysis of RADseq data sets",
    long_description=open('README.rst').read(),
    packages=find_packages(),
    install_requires=requires(),
    # dependencies=dependency_links(),
    entry_points={
        'console_scripts': [
            'ipyrad = ipyrad.__main__:main',
            'tetrad = ipyrad.analysis.__tetrad_cli__:main',
        ],
    },
    data_files=[('bin', glob.glob("./bin/*"))],
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
