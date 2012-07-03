#/usr/bin/env python

#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup
from setuptools import Extension
#from Cython.Distutils import build_ext
import numpy


setup(
  name = 'pflexible',
  version = '0.9.0',
  author = 'John F. Burkhart',
  author_email = 'jfburkhart@gmail.com',
  url = 'http://niflheim.nilu.no/~burkhart/pflexible',
  description = 'A Python interface to FLEXPART data.',
  license = 'Creative Commons',
  #cmdclass = {'build_ext': build_ext},
  #ext_modules = [Extension('pflexcy', ['pflexcy.pyx'])],
  ext_modules = [Extension('pflexcy', ['pflexcy.c'], include_dirs=[numpy.get_include()])],
  py_modules = ['pflexible','mapping'],
)


