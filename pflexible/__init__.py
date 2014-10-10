WELCOME = """
SYNOPSIS
========

    pflexible [-h] [-v,--verbose] [--version]

DESCRIPTION
===========

    pflexible: A python module for working with FLEXPART Output.


EXAMPLES
========

    #TODO:

    A lot! This is just a starting point. See the doc strings
    for information about the various functions.


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

    This script follows creative commons usage.

"""
#print WELCOME

from .version import __version__
import subprocess
import os
import sys
import numpy


def _get_hg_description(path_):
    """ Get the output of hg summary when executed in a given path. """

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

_hg_description = _get_hg_description(__path__[0])


def print_versions():
    """Print all the versions for packages that pflexible relies on."""
    print("-=" * 38)
    print("pflexible version: %s" % __version__)
    if _hg_description:
        print("pflexible hg summary:    %s" % _hg_description)
    print("NumPy version:     %s" % numpy.__version__)
    print("Python version:    %s" % sys.version)
    if os.name == "posix":
        (sysname, nodename, release, version, machine) = os.uname()
        print("Platform:          %s-%s" % (sys.platform, machine))
    print("Byte-ordering:     %s" % sys.byteorder)
    print("-=" * 38)



# Import the public functions here
from .flexpart_read import read_header, read_trajectories
from .grid_read import read_grid, get_slabs, fill_grids
from .data_structures import BinaryFile, Structure, Header
from .tests.all import test
