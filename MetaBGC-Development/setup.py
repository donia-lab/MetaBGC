"""Setuptools magic to install MetaBGC."""
import glob
import os
from setuptools import setup, find_packages
import subprocess
import sys


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname),encoding="utf8").read()

long_description = read('README.md')

install_requires = [
    'click',
    'numpy',
    'matplotlib',
    'scipy',
    'biopython >= 1.72',
    'scikit-learn >= 0.20.1',
    'pandas >= 0.19.2',
    'rpy2 >= 2.9.1'
    'tzlocal'
]

def read_version():
    """Read the version from the appropriate place in the library."""
    for line in open(os.path.join("metabgc",'metabgc_cmds.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip('"')

setup(
    name="metabgc",
    python_requires='>=3.6',
    version=read_version(),
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_data={'': ['license.txt']},
    author='MetaBGC development team.',
    author_email='ab50@princeton.edu',
    description='Metagenomic identifier of Biosynthetic Gene Clusters.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['metabgc=metabgc.metabgc_cmds:main'],
    },
    url='https://github.com/donia-lab/MetaBGC',
    license='GNU General Public License v3',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
    ]
)