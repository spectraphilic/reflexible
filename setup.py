#/usr/bin/env python

#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup
from setuptools import Extension
#from Cython.Distutils import build_ext
import numpy

# pflexible version
VERSION = open('VERSION').read().strip()
# Create the version.py file
open('pflexible/version.py', 'w').write('__version__ = "%s"\n' % VERSION)


setup(
  name = 'pflexible',
  version = VERSION,
  author = 'John F. Burkhart',
  author_email = 'jfburkhart@gmail.com',
  url = 'http://niflheim.nilu.no/~burkhart/pflexible',
  description = 'A Python interface to FLEXPART data.',
  license = 'Creative Commons',
  ext_modules = [Extension('pflexible.pflexcy', ['pflexible/pflexcy.c'],
                           include_dirs=[numpy.get_include()])],
  packages = ['pflexible'],
  py_modules = ['mapping'],  # XXX that should be converted into a package too?
)
