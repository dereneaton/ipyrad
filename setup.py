#!/usr/bin/env python2.7

from setuptools import setup, find_packages
import re

requirements = [
    'numpy>=1.9',
    'scipy>0.10',
    'h5py',
    'pyzmq>14.5',
    'dill>0.2',
    'sphinx',
    'ipyparallel'
    ]

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

    entry_points={
            'console_scripts': [
                'ipyrad = ipyrad.__main__:main',
            ],
    },

    license='GPL',

    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
)
