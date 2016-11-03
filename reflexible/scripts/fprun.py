#!/usr/bin/env python
from __future__ import print_function

"""Utility to print different data from a Flexpart run.

:Author: Francesc Alted
:Contact:  francesc@blosc.org, john.burkhart@geo.uio.no
:Created:  2016-11-03

Help to use this script can be get with::

  $ fprun.py -h

"""

import os
import sys

import reflexible as rf


def main():
    """Parses the passed command line arguments.

    The passed arguments will be used to create the netCDF4 file.
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--nested", action="store_true",
        help="Use a nested output.",
        )
    parser.add_argument(
        "-C", "--command", action="store_true",
        help="Print the COMMAND contents.",
        )
    parser.add_argument(
        "-R", "--releases", action="store_true",
        help="Print the RELEASES contents.",
        )
    parser.add_argument(
        "-S", "--species", action="store_true",
        help="Print the SPECIES contents.",
        )
    parser.add_argument(
        "-H", "--header-key",
        help="Print the contents of H[HEADER_KEY].",
        )
    parser.add_argument(
        "-K", "--header-keys", action="store_true",
        help="Print all the HEADER keys.",
        )
    parser.add_argument(
        "pathnames", nargs="?",
        help="The Flexpart pathnames file stating where options and output are."
        )

    args = parser.parse_args()

    if args.pathnames is None:
        # The FLEXPART pathnames file is mandatory
        parser.print_help()
        sys.exit(1)

    if not os.path.isfile(args.pathnames):
        raise IOError("PLEXPART pathnames not found in '%s'" % args.pathnames)

    fprun = rf.Flexpart(args.pathnames, nested=args.nested)
    print(fprun)

    if args.releases:
        print("Releases:", repr(fprun.Releases))
    if args.command:
        print("Command:", fprun.Command)
    if args.species:
        print("Species:", fprun.Species)
    if args.header_key:
        print("Header[%s]:\n" % args.header_key, fprun.Header[args.header_key])
    if args.header_keys:
        print("Header keys:", list(fprun.Header.keys()))


if __name__ == '__main__':
    main()
