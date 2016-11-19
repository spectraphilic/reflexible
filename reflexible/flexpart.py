"""
Class that can open a FLEXPART run and access to the configfiles and outputs.
"""

from __future__ import print_function
from __future__ import absolute_import

import glob
import warnings

from reflexible import Header
from reflexible.scripts import read_conffiles, read_species
import reflexible.conv2netcdf4 as conv


class Flexpart(object):

    def __init__(self, pathnames, nested=False):
        """Open a Flexpart run based on pathnames.

        Parameters
        ----------
        pathnames : pathname
            File that contains the paths for <options> and <output> dirs
        nested : bool
            Whether the nested version of the header would be read
        """
        self.pathnames = pathnames
        self.nested = nested
        self.fp_options, self.fp_output = conv.get_fpdirs(pathnames)
        self.ncfile = None

        # First try to use the NetCDF4 files, if they are in the output dir
        ncfiles = glob.glob(self.fp_output + '/*.nc')
        if ncfiles:
            ncfile = [f for f in ncfiles if ("nest" in f) == nested][0]
            self.Header = Header(ncfile, nested=nested, absolute_path=True)
            self.ncfile = ncfile
        else:
            # If not, fall back to the original Flexpart format
            warnings.warn(
                "NetCDF4 files not found in output directory '{}'.  "
                "You can always generate them from data there "
                "with the `create_ncfile` command line utility.".format(
                    self.fp_output))
            self.Header = conv.Header(self.fp_output, nested=nested)

        # Other config files
        self.Command = read_conffiles("COMMAND", self.fp_options, None)
        self.Releases = read_conffiles("RELEASES", self.fp_options, None)
        self.Species = read_species(self.fp_options, self.Header.nspec)

    def __str__(self):
        return "Flexpart('%s', nested=%s)" % (self.pathnames, self.nested)


if __name__ == "__main__":
    import sys

    fprun = Flexpart(sys.argv[1])
    print("Options dir:", fprun.fp_options)
    print("Output dir:", fprun.fp_output)
    print("Releases:", repr(fprun.Releases))
    print("Command:", fprun.Command)
    print("Species:", fprun.Species)
