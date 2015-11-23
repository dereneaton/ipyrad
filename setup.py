#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import glob
import re

requirements = [
    'pip>7.0',
    'cython',
    'scipy>0.10',
    'numpy>=1.9',
    'pandas',
    'dill>0.2',
    'sphinx',
    'ipython>=4.0',
    'ipyparallel',
    'ipykernel>=4.1',
    'jupyter',
    'jupyter-client>=4.1'
    ]

#import ipyrad
#version=ipyrad.__version__,
    
setup(
    name="ipyrad",

    version=re.search(
       r"^__version__ = ['\"]([^'\"]*)['\"]",
           open(
               "ipyrad/__init__.py",
               "r").read(),
       re.M).group(1),

    url="https://github.com/dereneaton/ipyrad",

    author="Deren Eaton",
    author_email="deren.eaton@yale.edu",

    description="Interactive assembly and analysis of RADseq data sets",
    long_description=open('README.rst').read(),

    packages=find_packages(),
    
    install_requires=[requirements],
    extras_require = {
        'plotting': ['toyplot>0.0.8'],
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
