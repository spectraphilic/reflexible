"""
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

from .version import __version__

# Import the public functions here
from .flexpart_read import read_header
from .grid_read import read_grid, get_slabs
