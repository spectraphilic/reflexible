#/usr/bin/env python

from __future__ import print_function

from distutils.dep_util import newer
import os.path
from setuptools import setup
import subprocess


# pflexible version
VERSION = open('VERSION').read().strip()
# Create the version.py file
open('pflexible/version.py', 'w').write('__version__ = "%s"\n' % VERSION)

# Build the FortFlex extension if necessary
if (not os.path.exists("pflexible/FortFlex.so") or
    newer("fortflex/FortFlex.f", "pflexible/FortFlex.so")):
    try:
        print(subprocess.check_output(
            "cd fortflex; sh build_FortFlex.sh", shell=True))
    except:
        print("Problems compiling the FortFlex module.  "
              "Will continue using a slower fallback...")
    else:
        print("FortFlex.so extension has been created in pflexible/!")

setup(
  name = 'pflexible',
  version = VERSION,
  author = 'John F. Burkhart',
  author_email = 'jfburkhart@gmail.com',
  url = 'http://niflheim.nilu.no/~burkhart/pflexible',
  description = 'A Python interface to FLEXPART data.',
  license = 'Creative Commons',
  # ext_modules = [Extension('pflexible.pflexcy', ['pflexible/pflexcy.c'],
  #                          include_dirs=[numpy.get_include()])],
  packages = [
      'pflexible',
      #'pflexible.tests',   # not included because tests need data samples
      ],
)
