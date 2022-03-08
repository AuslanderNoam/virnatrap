import setuptools
from setuptools import setup
import os
import sys
import subprocess

long_description = "viRNAtrap is a package to generate predicted viral contigs from unmapped RNAseq reads"


requirements = [
    'numpy>=1.17.3',
    'tensorflow>=2.0.0',
    'argparse',
]


setup(
    name="virnatrap",
    version="1",
    author="Add list",
    author_email="aelbasir@wistar.org, nauslander@wistar.org",
    description="Extract viral contigs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=requirements,
    test_suite='nose.collector',
    tests_require=['nose'],
    entry_points = {
        'console_scripts': [
            'virnatrap-predict=virnatrap.command_line:virnatrap_predict',
        ],
    }
)
