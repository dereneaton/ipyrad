#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import glob
import re



def requires():
    """ gets packages from requirements.txt """
    with open('requirements.txt') as infile:
        return infile.read().splitlines()


def dependency_links():
    """
    return: the package specifications
    """
    with open('constraints.txt') as infile:
        return infile.read().splitlines()


## Auto-update ipyrad version from git repo tag

# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
initfile = "ipyrad/__init__.py"

cur_version = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                    open(initfile, "r").read(),
                    re.M).group(1)

setup(
    name="ipyrad",

    version = cur_version,

    url="https://github.com/dereneaton/ipyrad",

    author="Deren Eaton",
    author_email="deren.eaton@yale.edu",

    description="Interactive assembly and analysis of RADseq data sets",
    long_description=open('README.rst').read(),

    packages=find_packages(),
    install_requires=requires(),
    dependencies=dependencies(),

    entry_points={
            'console_scripts': [
                'ipyrad = ipyrad.__main__:main',
            ],
    },

    data_files = [ ( 'bin', glob.glob("./bin/*") ) ],

    license='GPL',
    
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
