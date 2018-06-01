#/usr/bin/env python

from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
import setuptools
import subprocess
import sysconfig


# reflexible version
VERSION = open('VERSION').read().strip()
# Create the version.py file
open('reflexible/version.py', 'w').write('__version__ = "%s"\n' % VERSION)

# Build the FortFlex extension if necessary
ext_suffix = sysconfig.get_config_var('EXT_SUFFIX') or '.so'
fortflex_f = os.path.join('reflexible', 'conv2netcdf4', 'fortflex', 'FortFlex.f')
fortflex_so = os.path.join('reflexible', 'conv2netcdf4', 'FortFlex' + ext_suffix)
if not os.path.exists(fortflex_so) or newer(fortflex_f, fortflex_so):
    try:
        print(subprocess.check_output(
            "cd reflexible/conv2netcdf4/fortflex; "
            "sh build_FortFlex.sh", shell=True))
    except subprocess.CalledProcessError as e:
        print(e.output)
        print("Problems compiling the FortFlex module.  "
              "Will continue using a slower fallback...")
    else:
        print("FortFlex extension has been created in {}".format(fortflex_so))


def find_package_data(pdir):
    """Recursively get all the data files that are in pdir."""
    return [(d, [os.path.join(d, f) for f in files])
            for d,folders,files in os.walk(pdir)]


data_files = [
    ('reflexible/conv2netcdf4', [fortflex_so]),
    ('reflexible', ['reflexible/mapping_db.yml']),
] + find_package_data('reflexible/uio_examples')

setuptools.setup(
    name = 'reflexible',
    version = VERSION,
    author = 'John F. Burkhart, Francesc Alted',
    author_email = 'jfburkhart@gmail.com',
    url = 'https://github.com/spectraphilic/reflexible',
    description = 'A Python interface to FLEXPART data.',
    license = 'Creative Commons',
    packages = [
        'reflexible',
        'reflexible.scripts',
        'reflexible.conv2netcdf4',
        'reflexible.tests',
        ],
    data_files = data_files,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'create_ncfile = reflexible.scripts.create_ncfile:main',
            'readpartpositions = reflexible.scripts.readpartpositions:main',
            'fprun = reflexible.scripts.fprun:main',
        ]
    },
)
