WELCOME = """
DESCRIPTION
===========

    reflexible: A python package for working with FLEXPART Output.


EXAMPLES
========

    #TODO:

    A lot! This is just a starting point. See
    http://reflexible.readthedocs.org for a gentle introduction to the package
    and the doc strings for more information about the various functions.

    Also, have a look at the different tests in reflexible/tests for more hints
    on how to use the package.

AUTHOR
======

    JFB: John F Burkhart <jfburkhart@gmail.com>

CONTRIBUTORS
============

    HSO: Harald Sodemann
    SEC: Sabine Eckhardt
    AST: Andreas Stohl

    Many functions are adaptations of Fortran / NCARG programs (AST)
    or Matlab functions (HSO/SEC).

LICENSE
=======

    This package follows creative commons usage.

"""
#print WELCOME

import os
import datetime as dt

from .version import __version__


this_dir = __path__[0]


def _get_hg_description(path_):
    """ Get the output of hg summary when executed in a given path. """
    import subprocess

    # make an absolute path if required, for example when running in a clone
    if not os.path.isabs(path_):
        path_ = os.path.join(os.environ['PWD'], path_)
    # look up the commit using subprocess and hg summary
    try:
        # redirect stderr to stdout to make sure the hg error message in case
        # we are not in a hg repo doesn't appear on the screen and confuse the
        # user.
        out = subprocess.check_output(
            ["hg", "summary"], cwd=path_, stderr=subprocess.STDOUT).strip()
        out = "".join(['\n  ' + l for l in out.split('\n')])
        return out
    except OSError:  # in case hg wasn't found
        pass
    except subprocess.CalledProcessError:  # not in hg repo
        pass

_hg_description = _get_hg_description(this_dir)


def print_versions():
    """Print all the versions for packages that reflexible relies on."""
    import numpy
    import sys
    print("-=" * 38)
    print("reflexible version: %s" % __version__)
    if _hg_description:
        print("reflexible hg summary:    %s" % _hg_description)
    print("NumPy version:     %s" % numpy.__version__)
    print("Python version:    %s" % sys.version)
    if os.name == "posix":
        (sysname, nodename, release, version, machine) = os.uname()
        print("Platform:          %s-%s" % (sys.platform, machine))
    print("Byte-ordering:     %s" % sys.byteorder)
    print("-=" * 38)


# Some data sources (for testing purposes mainly)
datasets = {
    'Fwd1_V10.0': os.path.join(this_dir, "uio_examples/Fwd1_V10.0"),
    'Bwd1_V9.02': os.path.join(this_dir, "uio_examples/Bwd1_V9.02/outputs"),
    'Bwd2_V9.2beta': os.path.join(this_dir, "uio_examples/Bwd2_V9.2beta/outputs"),
    'Fwd1_V9.02': os.path.join(this_dir, "uio_examples/Fwd1_V9.02/outputs"),
    'Fwd2_V9.02': os.path.join(this_dir, "uio_examples/Fwd2_V9.02/outputs"),
    'HelloWorld_V9.02': os.path.join(this_dir, "uio_examples/HelloWorld_V9.02/outputs"),
    }


# Import the public functions here
from .tests.all import test
from . import conv2netcdf4
from .scripts.create_ncfile import create_ncfile
from .data_structures import Header
from .data_structures import Command
from .data_structures import Release
from .data_structures import ReleasePoint
from .base_read import read_trajectories
