#/usr/bin/env python

from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
from setuptools import setup
import subprocess


# reflexible version
VERSION = open('VERSION').read().strip()
# Create the version.py file
open('reflexible/version.py', 'w').write('__version__ = "%s"\n' % VERSION)

# Build the FortFlex extension if necessary
if (not os.path.exists("reflexible/conv2netcdf4/FortFlex.so") or
    newer("reflexible/conv2netcdf4/fortflex/FortFlex.f",
          "reflexible/conv2netcdf4/FortFlex.so")):
    try:
        print(subprocess.check_output(
            "cd reflexible/conv2netcdf4/fortflex; "
            "sh build_FortFlex.sh", shell=True))
    except:
        print("Problems compiling the FortFlex module.  "
              "Will continue using a slower fallback...")
    else:
        print("FortFlex.so extension has been created in reflexible/conv2netcdf4/!")


def find_package_data(pdir):
    """Recursively get all the data files that are in pdir."""
    return [(d, [os.path.join(d, f) for f in files])
            for d,folders,files in os.walk(pdir)]


setup(
    name = 'reflexible',
    version = VERSION,
    author = 'John F. Burkhart',
    author_email = 'jfburkhart@gmail.com',
    url = 'http://niflheim.nilu.no/~burkhart/reflexible',
    description = 'A Python interface to FLEXPART data.',
    license = 'Creative Commons',
    packages = [
        'reflexible',
        'reflexible.scripts',
        'reflexible.conv2netcdf4',
        'reflexible.tests',
        'reflexible.legacy',
        ],
    data_files = [
        ('reflexible/conv2netcdf4', ['reflexible/conv2netcdf4/FortFlex.so'])] + \
        find_package_data('reflexible/uio_examples'),
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'create_ncfile = reflexible.scripts.create_ncfile:main',
            'readpartpositions = reflexible.scripts.readpartpositions:main',
        ]
    },

)
