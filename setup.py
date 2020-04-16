#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

from numpy.distutils.core import Extension

# try:
#     from setuptools import setup, find_packages
# except ImportError:
#     from distutils.core import setup

ext1 = Extension(name='tspack',
                 sources=['tspack.pyf', 'tspack.f'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='pytspack',
          description="Python interface to the fortran TSPACK",
          author="Barry Baker",
          author_email="barry.baker@noaa.gov",
          license='MIT',
          ext_modules=[ext1]
          )
