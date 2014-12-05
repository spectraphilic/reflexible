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
import platform
import getpass
import datetime
import os.path
import netCDF4 as nc

from reflexible.conv2netcdf4 import Header, read_grid, read_command, read_commandV9


def write_metadata(H, command, ncid):
    # hes CF convention requires these attributes
    ncid.Conventions = 'CF-1.6'
    ncid.title = 'FLEXPART model output'
    ncid.institution = 'NILU'
    ncid.source = H.version + ' model output'
    date = "%d-%d-%d %d:%d" % datetime.datetime.now().timetuple()[:5]
    zone = "NA"
    ncid.history = date + ' ' + zone + '  created by '+ getpass.getuser() + ' on ' + platform.node()
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

    # COMMAND settings
    ncid.itsplit = command['T_PARTSPLIT']
    ncid.linit_cond = command['LINIT_COND']
    ncid.lsynctime = command['SYNC']
    ncid.ctl = command['CTL']
    ncid.ifine = command['IFINE']
    ncid.iout = command['IOUT']
    ncid.ipout = command['IPOUT']
    ncid.lagespectra = command['LAGESPECTRA']
    ncid.ipin = command['IPIN']
    ncid.ioutputforeachrelease = command['OUTPUTFOREACHRELEASE']
    ncid.iflux = command['IFLUX']
    ncid.mdomainfill = command['MDOMAINFILL']
    ncid.mquasilag = command['MQUASILAG']
    ncid.nested_output = command['NESTED_OUTPUT']
    ncid.surf_only = 0  # XXX what's that?  a new option maybe?


def write_header(H, command, ncid):
    nnx = H.numxgrid
    nny = H.numygrid
    nnz = H.numzgrid

    # Create dimensions

    # time
    timeDimID = ncid.createDimension('time', None)
    adate, atime = str(H.ibdate), str(H.ibtime)
    timeunit = 'seconds since ' + adate[:4] + '-' + adate[4:6] + \
        '-'+ adate[6:8] + ' ' + atime[:2] + ':' + atime[2:4]

    # lon
    lonDimID = ncid.createDimension('longitude', nnx)
    # lat
    latDimID = ncid.createDimension('latitude', nny)
    # level
    levDimID = ncid.createDimension('height', nnz)
    # number of species
    nspecDimID = ncid.createDimension('numspec', H.nspec)
    # number of release points
    pointspecDimID = ncid.createDimension('pointspec', H.numpointspec)  # XXX or H.maxpoint?
    # number of age classes
    nageclassDimID = ncid.createDimension('nageclass', H.nageclass)
    # dimension for release point characters
    ncharDimID = ncid.createDimension('nchar', 45)
    # number of actual release points
    npointDimID = ncid.createDimension('numpoint', H.numpoint)

    # create variables

    # time
    tID = ncid.createVariable('time', 'i4', ('time',))
    tID.units = timeunit
    tID.calendar = 'proleptic_gregorian'

    # lon
    lonID = ncid.createVariable('longitude', 'f4', ('longitude',))
    lonID.long_name = 'longitude in degree east'
    lonID.axis = 'Lon'
    lonID.units = 'degrees_east'
    lonID.standard_name = 'grid_longitude'
    lonID.description = 'grid cell centers'

    # lat
    latID = ncid.createVariable('latitude', 'f4', ('latitude',))
    latID.long_name = 'latitude in degree north'
    latID.axis = 'Lat'
    latID.units = 'degrees_north'
    latID.standard_name = 'grid_latitude'
    latID.description = 'grid cell centers'

    # height
    levID = ncid.createVariable('height', 'f4', 'height')
    # levID.axis = 'Z'
    levID.units = 'meters'
    levID.positive = 'up'
    levID.standard_name = 'height'
    levID.long_name = 'height above ground'


def create_ncfile(fddir, nested, outfile=None):
    H = Header(fddir, nested=nested)
    if H.direction == "forward":
        fprefix = 'grid_conc_'
    else:
        fprefix = 'grid_time_'

    command_path = os.path.join(os.path.dirname(fddir), "options/COMMAND")
    # XXX This needs to be checked out, as I am not sure when the new format started
    command = read_commandV9(command_path)

    path = os.path.dirname(fddir)
    fprefix = os.path.join(path, fprefix)
    if outfile is None:
        if H.nested:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + "_nest.nc"
        else:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + ".nc"
    else:
        ncfname = outfile
    cache_size = 16 * H.numxgrid * H.numygrid * H.numzgrid

    ncid = nc.Dataset(ncfname, 'w', chunk_cache=cache_size)
    write_metadata(H, command, ncid)
    write_header(H, command, ncid)
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
