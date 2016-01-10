#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import subprocess
import fileinput
import glob
import os
import re

requirements = [
    'pip>7.0',
    'dill>0.2',
    'sphinx',
    'cython',
    'ipython>=4.0.0',
    'ipyparallel',
    'ipykernel>=4.1',
    'jupyter',
    'jupyter-client>=4.1',
    'nbconvert',
    'scipy>0.10',
    'numpy>=1.9',
    'pandas',
    'h5py'
    ]

## Auto-update ipyrad version from git repo tag

# Fetch version from git tags, and write to version.py.
# Also, when git is not available (PyPi package), use stored version.py.
initfile = "ipyrad/__init__.py"

cur_version = version_git = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                    open(initfile, "r").read(),
                    re.M).group(1)

# If git is available and working pull version from github, else read from ipyrad/__init__.py
try:
    version_git = subprocess.check_output(["git", "describe"]).rstrip().split("-")[0]
except:
    version_git = cur_version

if( not cur_version == version_git ):
    # Write version to ipyrad/__init__.py
    for line in fileinput.input(initfile, inplace=1):
            if "__version__" in line:
                line = "__version__ = \""+version_git+"\""    
            print(line.strip("\n"))
    subprocess.call(["git", "add", initfile])
    subprocess.call(["git", "commit", "-m \"Updating ipyrad/__init__.py to"+\
                        "version - {}".format(version_git)])
    subprocess.call(["git", "push"])

setup(
    name="ipyrad",

    version = version_git,
#    version=re.search(
#       r"^__version__ = ['\"]([^'\"]*)['\"]",
#           open(
#               "ipyrad/__init__.py",
#               "r").read(),
#       re.M).group(1),

    url="https://github.com/dereneaton/ipyrad",

    author="Deren Eaton",
    author_email="deren.eaton@yale.edu",

    description="Interactive assembly and analysis of RADseq data sets",
    long_description=open('README.rst').read(),

    packages=find_packages(),
    
    install_requires=[requirements],
    extras_require = {
        'plotting': ['toyplot>0.0.9'],
    },

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
