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
    ncid.ldirect = 1 if H.direction == "forward" else -1
    ncid.ibdate = H.available_dates[0][:8]
    ncid.ibtime = H.available_dates[0][8:]
    ncid.iedate = H.available_dates[-1][:8]
    ncid.ietime = H.available_dates[-1][8:]
    ncid.loutstep = H.loutstep
    ncid.loutaver = H.loutaver
    ncid.loutsample = H.loutsample
    ncid.lsubgrid = H.lsubgrid
    ncid.lconvection = H.lconvection
    ncid.ind_source = H.ind_source
    ncid.ind_receptor = H.ind_receptor

    # From here on, the COMMAND settings are not present in H header
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


def write_header(H, ncid):
    pass


def create_ncfile(fddir, nested, outfile=None):
    H = Header(fddir, nested=nested)
    if H.direction == "forward":
        fprefix = 'grid_conc_'
    else:
        fprefix = 'grid_time_'
    path = os.path.dirname(fddir)
    fprefix = os.path.join(path, fprefix)
    if outfile is None:
        if H.nested:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + "_nest.nc"
        else:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + ".nc"
    else:
        ncfname = outfile
    ncid = nc.Dataset(ncfname, 'w')
    write_metadata(H, ncid)
    write_header(H, ncid)
    ncid.close()
    return ncfname


def main():
    try:
        fddir = sys.argv[1]
    except IndexError:
        print("USAGE: create_ncfile output_dir [nested]")
        sys.exit(1)
    if fddir.endswith('/'):
        # Remove the trailing '/'
        fddir = fddir[:-1]
    try:
        nested = bool(sys.argv[2])
    except IndexError:
        nested = False
    ncfname = create_ncfile(fddir, nested)
    print("New netCDF4 files is available in: '%s'" % ncfname)


if __name__ == '__main__':
    main()
