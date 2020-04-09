#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import sys
from pkg_resources import VersionConflict, require
from setuptools import setup

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

try:
    require('setuptools>=38.3')
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)

setup(name='pydeseq2',
      version='0.1',
      description='python package to run DESeq2 in python',
      long_description=readme,
      author='Federica Gervasoni',
      author_email='federica.gervasoni@unimi.it',
      license='Apache Software License 2.0',
      packages=find_packages(include=['pydeseq2']),
      python_requires='>=3.7',
      include_package_data=True,
      zip_safe=False)
