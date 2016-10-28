from __future__ import print_function
from __future__ import absolute_import

from reflexible.scripts import create_ncfile, read_conffiles, read_species
import reflexible.conv2netcdf4 as conv


class Flexpart(object):

    def __init__(self, pathnames, nested=False):
        """Open a Flexpart run based on pathnames.

        Parameters
        ----------
        pathnames : pathname
            File that contains the paths for <options> and <output> dirs
        nested : Bool
            Whether the nested version of the header would be read
        """
        self.fp_options, self.fp_output = conv.get_fpdirs(pathnames)
        self.Header = conv.Header(self.fp_output, nested=nested)
        self.Command = read_conffiles("COMMAND", self.fp_options, None)
        self.Releases = read_conffiles("RELEASES", self.fp_options, None)
        self.Species = read_species(self.fp_options, self.Header.nspec)


if __name__ == "__main__":
    import sys

    fprun = Flexpart(sys.argv[1])
    print("Options dir:", fprun.fp_options)
    print("Output dir:", fprun.fp_output)
    print("Species:", fprun.Species)
