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


def write_metadata(H, ncid):
    # hes CF convention requires these attributes
    ncid.Conventions = 'CF-1.6'
    ncid.title = 'FLEXPART model output'
    # ncid.institution = institution.strip()
    ncid.source = H.version + ' model output'
    # ncid.history = date[:4] + '-' + date[4:6] + '-' + date[6:8] + ' ' + time[0:2] + ':' + time(3:4) + ' ' + zone + '  created by '+ trim(login_name) + ' on ' + trim(host_name)
    ncid.references = 'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200'

    # attributes describing model run
    ncid.outlon0 = H.outlon0
    ncid.outlat0 = H.outlat0
    ncid.dxout = H.dxout
    ncid.dyout = H. dyout

    # COMMAND file settings
    ncid.ibdate = H.ibdate
    ncid.ibtime = H.ibtime
    # ncid.iedate = H.iedate    # XXX not present in H header
    # ncid.ietime = H.ietime    # XXX not present in H header
    ncid.loutstep = H.loutstep
    ncid.loutaver = H.loutaver
    ncid.loutsample = H.loutsample
    ncid.lsubgrid = H.lsubgrid
    ncid.lconvection = H.lconvection
    ncid.ind_source = H.ind_source
    ncid.ind_receptor = H.ind_receptor

    # From here on, the COMMAND settings are not present in H header
    # ncid.ldirect = H.ldirect
    # ncid.itsplit = H.itsplit
    # ncid.linit_cond = H.linit_cond
    # ncid.lsynctime = H.lsynctime
    # ncid.ctl = H.ctl
    # ncid.ifine = H.ifine
    # ncid.iout = H.iout
    # ncid.ipout = H.ipout
    # ncid.lagespectra = H.lagespectra
    # ncid.ipin = H.ipin
    # ncid.ioutputforeachrelease = H.ioutputforeachrelease
    # ncid.iflux = H.iflux
    # ncid.mdomainfill = H.mdomainfill
    # ncid.mquasilag = H.mquasilag
    # ncid.nested_output = H.nested_output
    # ncid.surf_only = H.surf_only


def create_ncfile(H, ncfile):
    ncid = nc.Dataset(ncfile, 'w')
    write_metadata(H, ncid)
    FD = read_grid(H, time_ret=-1, nspec_ret=0)
    ncid.close()


def main():
    fddir = sys.argv[1]
    try:
        outname = sys.argv[2]
    except:
        outname = "output.nc"
    H = Header(fddir)
    create_ncfile(H, "output.nc")
    print("New netCDF4 files is available in: '%s'" % outname)


if __name__ == '__main__':
    main()
