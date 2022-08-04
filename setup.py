"""Setuptools magic to install MetaBGC."""
import os
from setuptools import setup, find_packages


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname),encoding="utf8").read()

long_description = read('README.md')

install_requires = [
    'click',
    'numpy',
    'matplotlib',
    'scipy',
    'biopython >= 1.72, < 1.78',
    'scikit-learn >= 0.20.1',
    'pandas >= 0.19.2',
    'pytest'
]

def read_version():
    """Read the version from the appropriate place in the library."""
    for line in open(os.path.join("metabgc","src",'__main__.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip('"')

setup(
    name="metabgc",
    python_requires='>=3.6',
    version=read_version(),
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    author='MetaBGC development team.',
    author_email='ab50@princeton.edu',
    description='Metagenomic identifier of Biosynthetic Gene Clusters.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['metabgc=metabgc.src.__main__:main'],
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