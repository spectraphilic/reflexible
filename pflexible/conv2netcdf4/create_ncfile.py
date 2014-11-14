'''
Script to convert a FLEXPART dataset into a NetCDF4 file.

:Author: Francesc Alted
:Contact:  francesc@blosc.org
:Created:  2014-11-14
:Acknowledgment: Funding for the development of this code is provided
     through the iiSPAC project (NSF-ARC-1023651)

'''

from __future__ import print_function

import sys
import os.path
import netCDF4 as nc

from pflexible.conv2netcdf4 import Header, read_grid


def nc_create(H, ncfile):
    FD = read_grid(H, time_ret=-1, nspec_ret=0)
#     fdkeys = sorted(FD.keys())
#     assert fdkeys == ['grid_dates', 'options', (0, '20070121100000')]
    print("grid_dates:", FD['grid_dates'], FD['options'])
#     fdkeys_ = sorted(FD[(0, '20070121100000')].keys())
#     assert fdkeys_ == self.fdkeys


def main():
    fddir = sys.argv[1]
    H = Header(fddir)    
    nc_create(H, os.path.join(fddir, ".nc"))


if __name__ == '__main__':
    main()