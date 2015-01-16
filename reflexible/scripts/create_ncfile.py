from __future__ import print_function

"""Script to convert a FLEXPART dataset into a NetCDF4 file.

:Author: Francesc Alted
:Contact:  francesc@blosc.org
:Created:  2014-11-14
:Acknowledgment: Funding for the development of this code is provided
     through the iiSPAC project (NSF-ARC-1023651)

This script demonstrates how to convert FLEXDATA output files into a
file with netCDF4 format.  Help to use this script can be get with::

  $ create_ncfile.py -h

"""

import sys
import warnings
import platform
import getpass
import datetime
import os.path

import netCDF4 as nc
import numpy as np

from reflexible.conv2netcdf4 import Header, read_command

UNITS = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
"""Used in combination with H.nested to determine the value of the
   ncid.iout attribute when the COMMAND file is not available.
"""

# Global variables for hosting different command line arguments.
# The default values here are not relevant.
COMPLEVEL = 9
MIN_SIZE = False


def read_releases(path):
    """Read metadata from a RELEASES path and return it as a dict.

    Only 'release_point_names' returned.
    """
    rpnames = []
    with open(path) as f:
        prev_line = None
        for line in f:
            if prev_line is not None and "comment" in line:
                rpnames.append(prev_line.strip())
            prev_line = line
    # Return just the release point names for now
    return {"release_point_names": np.array(rpnames, dtype="S45")}


def output_units(ncid):
    """The ncid.units attribute computation.

    This function computes the value of the ncid.units attribute required
    when concentration output, wet and dry deposition variables are created.

    Parameters
    ----------
    ncid : Python object
        the object associated to the netCDF4 file

    Return
    ------
    string
        the value of the ncid.units attributte
    """
    if ncid.ldirect == 1:
        # forward simulation
        if ncid.ind_source == 1:
            if ncid.ind_receptor == 1:
                units = 'ng m-3'
            else:
                units = 'ng kg-1'
        else:
            if ncid.ind_receptor == 1:
                units = 'ng m-3'
            else:
                units = 'ng kg-1'
    else:
        # backward simulation
        if ncid.ind_source == 1:
            if ncid.ind_receptor == 1:
                units = 's'
            else:
                units = 's m3 kg-1'
        else:
            if ncid.ind_receptor == 1:
                units = 's kg m-3'
            else:
                units = 's'

    return units


def write_metadata(H, command, ncid):
    """Write attributes to the netCDF4 file object (i.e. ncid).

    These function writes metadata (as attributes) to the netCDF4 file object.
    The metadata is read from the Python objects associated to the header and
    command files.

    Parameters
    ----------
    H: Python object
      the Header object.
    command: Python object
      the COMMAND file object.
    ncid: Python object
      the netCDF4 file object.
    """
    # the CF convention requires these attributes
    ncid.Conventions = 'CF-1.6'
    ncid.title = 'FLEXPART model output'
    ncid.institution = 'NILU'
    ncid.source = H.version + ' model output'
    date = "%d-%d-%d %d:%d" % datetime.datetime.now().timetuple()[:5]
    zone = "NA"
    ncid.history = (date + ' ' + zone + ' created by ' +
                    getpass.getuser() + ' on ' + platform.node())
    ncid.references = ("Stohl et al., Atmos. Chem. Phys., 2005, "
                       "doi:10.5194/acp-5-2461-200")

    # attributes describing model run
    ncid.outlon0 = H.outlon0
    ncid.outlat0 = H.outlat0
    ncid.dxout = H.dxout
    ncid.dyout = H.dyout

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

    # additional COMMAND settings
    if len(command) > 0:
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
        # This is a new option in V9.2
        ncid.surf_only = command.get('SURF_ONLY', 0)


def write_header(H, ncid, wetdep, drydep, write_releases):
    """Create netCDF4 dimensions and variables.

    Create the netCDF4 variables (and the required dimensions) that will be
    stored in the netCDF4 file. The created variables are available through
    the dictionary ncid.variables as::

      ncid.variables['var_name']

    Beware that variables are created but not filled with data.

    Parameters
    ----------
    H : Python object
      The Header object.
    ncid : Python object
      The netCDF4 file object.
    wetdep : boolean
      True if wet depositions have to be written into the netCDF4 file.
    drydep : boolean
      True if dry depositions have to be written into the netCDF4 file.

    Return
    ------
    int
      1 -> conc. output (ng/m3)
      2 -> mixing ratio (pptv)
      3 -> both
      4 -> plume traject
      5 -> 1 + 4
    """
    # global UNITS

    if hasattr(ncid, "iout"):
        iout = ncid.iout
    else:
        # If IOUT is not available (no COMMAND file), guess the value of IOUT
        # Add 1 because the count starts with 1
        unit_i = UNITS.index(H.unit) + 1
        iout = (unit_i) + (H.nested * 5)

    # Create dimensions

    # time
    ncid.createDimension('time', None)
    adate, atime = str(H.ibdate), str(H.ibtime)
    timeunit = 'seconds since ' + adate[:4] + '-' + adate[4:6] + \
        '-' + adate[6:8] + ' ' + atime[:2] + ':' + atime[2:4]
    # lon
    ncid.createDimension('longitude', H.numxgrid)
    # lat
    ncid.createDimension('latitude', H.numygrid)
    # level
    ncid.createDimension('height', H.numzgrid)
    # number of species
    ncid.createDimension('numspec', H.nspec)
    # number of release points
    ncid.createDimension('pointspec', H.numpointspec)
    # number of age classes
    ncid.createDimension('nageclass', H.nageclass)
    # dimension for release point characters
    ncid.createDimension('nchar', 45)
    # number of actual release points
    ncid.createDimension('numpoint', H.numpoint)

    # Create variables

    # time
    tID = ncid.createVariable('time', 'i4', ('time',),
                              zlib=True, complevel=COMPLEVEL)
    tID.units = timeunit
    tID.calendar = 'proleptic_gregorian'

    # lon
    lonID = ncid.createVariable('longitude', 'f4', ('longitude',),
                                zlib=True, complevel=COMPLEVEL)
    lonID.long_name = 'longitude in degree east'
    lonID.axis = 'Lon'
    lonID.units = 'degrees_east'
    lonID.standard_name = 'grid_longitude'
    lonID.description = 'grid cell centers'

    # lat
    latID = ncid.createVariable('latitude', 'f4', ('latitude',),
                                zlib=True, complevel=COMPLEVEL)
    latID.long_name = 'latitude in degree north'
    latID.axis = 'Lat'
    latID.units = 'degrees_north'
    latID.standard_name = 'grid_latitude'
    latID.description = 'grid cell centers'

    # height
    levID = ncid.createVariable('height', 'f4', ('height',),
                                zlib=True, complevel=COMPLEVEL)
    # levID.axis = 'Z'
    levID.units = 'meters'
    levID.positive = 'up'
    levID.standard_name = 'height'
    levID.long_name = 'height above ground'

    if write_releases:
        # RELCOM
        relcomID = ncid.createVariable('RELCOM', 'S45', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        # Fill RELCOM with default values ("NA" means Non-Available)
        relcomID[:] = np.array(["NA"] * H.numpoint, dtype="S45")
        relcomID.long_name = 'release point name'

        # RELLNG1
        rellng1ID = ncid.createVariable('RELLNG1', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellng1ID.units = 'degrees_east'
        rellng1ID.long_name = 'release longitude lower left corner'

        # RELLNG2
        rellng2ID = ncid.createVariable('RELLNG2', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellng2ID.units = 'degrees_east'
        rellng2ID.long_name = 'release longitude upper right corner'

        # RELLAT1
        rellat1ID = ncid.createVariable('RELLAT1', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellat1ID.units = 'degrees_north'
        rellat1ID.long_name = 'release latitude lower left corner'

        # RELLAT2
        rellat2ID = ncid.createVariable('RELLAT2', 'f4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        rellat2ID.units = 'degrees_north'
        rellat2ID.long_name = 'release latitude upper right corner'

        # RELZZ1
        relzz1ID = ncid.createVariable('RELZZ1', 'f4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relzz1ID.units = 'meters'
        relzz1ID.long_name = 'release height bottom'

        # RELZZ2
        relzz2ID = ncid.createVariable('RELZZ2', 'f4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relzz2ID.units = 'meters'
        relzz2ID.long_name = 'release height top'

        # RELKINDZ
        relkindzID = ncid.createVariable('RELKINDZ', 'i4', ('numpoint',),
                                         zlib=True, complevel=COMPLEVEL)
        relkindzID.long_name = 'release kind'

        # RELSTART
        relstartID = ncid.createVariable('RELSTART', 'i4', ('numpoint',),
                                         zlib=True, complevel=COMPLEVEL)
        relstartID.units = 'seconds'
        relstartID.long_name = 'release start relative to simulation start'

        # RELEND
        relendID = ncid.createVariable('RELEND', 'i4', ('numpoint',),
                                       zlib=True, complevel=COMPLEVEL)
        relendID.units = 'seconds'
        relendID.long_name = 'release end relative to simulation start'

        # RELPART
        relpartID = ncid.createVariable('RELPART', 'i4', ('numpoint',),
                                        zlib=True, complevel=COMPLEVEL)
        relpartID.long_name = 'number of release particles'

        # RELXMASS
        relxmassID = ncid.createVariable('RELXMASS', 'f4',
                                         ('numspec', 'numpoint'),
                                         zlib=True, complevel=COMPLEVEL)
        relxmassID.long_name = 'total release particles mass'

    # LAGE
    lageID = ncid.createVariable('LAGE', 'i4', ('nageclass',))
    lageID.units = 'seconds'
    lageID.long_name = 'age class'

    # ORO
    if not MIN_SIZE:
        oroID = ncid.createVariable('ORO', 'i4', ('longitude', 'latitude'),
                                    chunksizes=(H.numxgrid, H.numygrid),
                                    zlib=True, complevel=COMPLEVEL)
        oroID.standard_name = 'surface altitude'
        oroID.long_name = 'outgrid surface altitude'
        oroID.units = 'm'

    units = output_units(ncid)

    # Concentration output, wet and dry deposition variables (one per species)
    # Variables for concentration ouput have dimensions given by dIDs
    # Variables for wer and dry deposition have dimensions given by depdIDs
    dIDs = (
        'longitude', 'latitude', 'height', 'time', 'pointspec', 'nageclass')
    depdIDs = ('longitude', 'latitude', 'time', 'pointspec', 'nageclass')

    chunksizes = (H.numxgrid, H.numygrid, H.numzgrid, 1, 1, 1)
    dep_chunksizes = (H.numxgrid, H.numygrid, 1, 1, 1)
    for i in range(0, H.nspec):
        anspec = "%3.3d" % (i + 1)
        # iout: 1 conc. output (ng/m3), 2 mixing ratio (pptv), 3 both,
        # 4 plume traject, 5=1+4
        if iout in (1, 3, 5):
            var_name = "spec" + anspec + "_mr"
            sID = ncid.createVariable(var_name, 'f4', dIDs,
                                      chunksizes=chunksizes,
                                      zlib=True, complevel=COMPLEVEL)
            sID.units = units
            sID.long_name = H.species[i]
            # TODO: we still need to figure out how to get the next data
            # sID.decay = decay[i]
            # sID.weightmolar = weightmolar[i]
            # sID.ohreact = ohreact[i]
            # sID.kao = kao[i]
            # sID.vsetaver = vsetaver[i]
            # sID.spec_ass = spec_ass[i]
            ncid.variables[var_name] = sID
        if iout in (2, 3):
            var_name = "spec" + anspec + "_pptv"
            sID = ncid.createVariable(var_name, 'f4', dIDs,
                                      chunksizes=chunksizes,
                                      zlib=True, complevel=COMPLEVEL)
            sID.units = 'pptv'
            sID.long_name = H.species[i]
            # TODO: we still need to figure out how to get the next data
            # sID.decay = decay[i]
            # sID.weightmolar = weightmolar[i]
            # sID.ohreact = ohreact[i]
            # sID.kao = kao[i]
            # sID.vsetaver = vsetaver[i]
            # sID.spec_ass = spec_ass[i]
            ncid.variables[var_name] = sID

        if wetdep:
            var_name = "WD_spec" + anspec
            wdsID = ncid.createVariable(var_name, 'f4', depdIDs,
                                        chunksizes=dep_chunksizes,
                                        zlib=True, complevel=COMPLEVEL)
            wdsID.units = '1e-12 kg m-2'
            # TODO: we still need to figure out how to get the next data
            # wdsID.weta = weta[i]
            # wdsID.wetb = wetb[i]
            # wdsID.weta_in = weta_in[i]
            # wdsID.wetb_in = wetb_in[i]
            # wdsID.wetc_in = wetc_in[i]
            # wdsID.wetd_in = wetd_in[i]
            # wdsID.dquer = dquer[i]
            # wdsID.henry = henry[i]
            ncid.variables[var_name] = wdsID

        if drydep:
            var_name = "DD_spec" + anspec
            ddsID = ncid.createVariable(var_name, 'f4', depdIDs,
                                        chunksizes=dep_chunksizes,
                                        zlib=True, complevel=COMPLEVEL)
            ddsID.units = '1e-12 kg m-2'
            # TODO: we still need to figure out how to get the next data
            # dsID.dryvel = dryvel[i]
            # ddsID.reldiff = reldiff[i]
            # ddsID.henry = henry[i]
            # ddsID.f0 = f0[i]
            # ddsID.dquer = dquer[i]
            # ddsID.density = density[i]
            # ddsID.dsigma = dsigma[i]
            ncid.variables[var_name] = ddsID

    return iout


def write_variables(H, ncid, wetdep, drydep, iout, write_releases, releases):
    """Fill netCDF4 variables with data.

    The netCDF4 variables created in the ``write_header`` function are filled
    with data.

    Parameters
    ----------
    H : Python object
      the header object
    ncid : Python object
      the netCDF4 file object
    wetdep : boolean
      True if wet depositions have to be written into the netCDF4 file.
    drydep : boolean
      True if dry depositions have to be written into the netCDF4 file.
    iout : int
      1 -> conc. output (ng/m3)
      2 -> mixing ratio (pptv)
      3 -> both
      4 -> plume traject
      5 -> 1+4
    """
    # Fill variables with data.

    # time
    time_var = ncid.variables['time']
    time_var[:] = nc.date2num(H.available_dates_dt, units=time_var.units,
                              calendar=time_var.calendar)

    # longitudes (grid cell centers)
    ncid.variables['longitude'][:] = np.linspace(
        ncid.outlon0 + 0.5 * ncid.dxout,
        ncid.outlon0 + (H.numxgrid-0.5) * ncid.dxout,
        H.numxgrid)

    # latitudes (grid cell centers)
    ncid.variables['latitude'][:] = np.linspace(
        ncid.outlat0 + 0.5 * ncid.dyout,
        ncid.outlat0 + (H.numygrid-0.5) * ncid.dyout,
        H.numygrid)

    # levels
    ncid.variables['height'][:] = H.outheight

    if write_releases:
        # release point information
        ncid.variables['RELSTART'][:] = H.ireleasestart
        ncid.variables['RELEND'][:] = H.ireleaseend
        ncid.variables['RELKINDZ'][:] = H.kindz
        ncid.variables['RELLNG1'][:] = H.xp1
        ncid.variables['RELLNG2'][:] = H.xp2
        ncid.variables['RELLAT1'][:] = H.yp1
        ncid.variables['RELLAT2'][:] = H.yp2
        ncid.variables['RELZZ1'][:] = H.zpoint1
        ncid.variables['RELZZ2'][:] = H.zpoint2
        ncid.variables['RELPART'][:] = H.npart
        ncid.variables['RELXMASS'][:, :] = H.xmass
        relnames = releases["release_point_names"]
        ncid.variables['RELCOM'][:len(relnames)] = relnames

    # Age classes
    ncid.variables['LAGE'][:] = H.lage

    # Orography
    if not MIN_SIZE:
        ncid.variables['ORO'][:, :] = H.oro

    # Concentration output, wet and dry deposition variables (one per species)
    # Loop over all the species and dates
    for ispec in range(H.nspec):
        anspec = "%3.3d" % (ispec + 1)
        for idt, date in enumerate(H.available_dates):
            # read grid, as well as wet and dry depositions
            try:
                H.read_grid(nspec_ret=ispec, time_ret=idt,
                            getwet=wetdep, getdry=drydep)
                fd = H.FD[(ispec, date)]
            except IOError:
                # Oops, we have got an error while reading, so close the file
                ncid.close()
                # and re-raise the error
                raise

            # Fill concentration values
            if iout in (1, 3, 5):
                conc_name = "spec" + anspec + "_mr"
            if iout in (2, 3):
                conc_name = "spec" + anspec + "_pptv"
            conc = ncid.variables[conc_name]
            # (x, y, z, time, pointspec, nageclass) <-
            # (x, y, z, pointspec, nageclass)
            conc[:, :, :, idt, :, :] = fd.grid[:, :, :, np.newaxis, :, :]

            # wet and dry depositions
            # (x, y, time, pointspec, nageclass) <- (x, y, pointspec, nageclass)
            if wetdep:
                wet = ncid.variables["WD_spec" + anspec]
                wet[:, :, idt, :, :] = fd.wet[:, :, np.newaxis, :, :]
            if drydep:
                dry = ncid.variables["DD_spec" + anspec]
                dry[:, :, idt, :, :] = fd.dry[:, :, np.newaxis, :, :]


def read_conffiles(filename, fddir, path):
    """Read FLEXPART config files and return a dictionary of metadata."""
    if path is None:
        # Try in the options/ directory before giving up
        path = os.path.join(os.path.dirname(fddir), "options/%s" % filename)
    if not os.path.isfile(path):
        warnings.warn(
            "The %s file cannot be found.  Continuing without it!" % filename)
        return {}
    try:
        if filename == "COMMAND":
            return read_command(path)
        elif filename == "RELEASES":
            return read_releases(path)
    except IOError:
        warnings.warn(
            "The %s file format is not supported.  "
            "Continuing without it!" % filename)
    return {}


def create_ncfile(fddir, nested, wetdep=False, drydep=False,
                  command_path=None, releases_path=None,
                  write_releases=True, dirout=None, outfile=None):
    """Main function that create a netCDF4 file from a FLEXPART output.

    Parameters
    ----------
    fddir : string
      the directory where the FLEXDATA output files are stored.
    nested : bool
      use a nested output.
    wetdep : bool
      write wet deposition in the netCDF4 file.
    drydep : bool
      write dry deposition in the netCDF4 file.
    command_path : string
      path for the associated COMMAND file.
    releases_path : string
      path for the associated RELEASES file.
    write_releases : string
      whether output of release point information.
    dirout : string
      the dir where the netCDF4 file will be created.
    outfile : string
      the complete path of the output file (overrides the ``dirout`` argument)

    Return
    ------
    string
      the full path of the netCDF4 file to be created.
    """
    if fddir.endswith('/'):
        # Remove the trailing '/'
        fddir = fddir[:-1]

    H = Header(fddir, nested=nested)

    if H.direction == "forward":
        fprefix = 'grid_conc_'
    else:
        fprefix = 'grid_time_'

    command = read_conffiles("COMMAND", fddir, command_path)
    releases = read_conffiles("RELEASES", fddir, releases_path)

    if outfile:
        # outfile has priority over previous flags
        ncfname = outfile
    else:
        if dirout is None:
            path = os.path.dirname(fddir)
            fprefix = os.path.join(path, fprefix)
        else:
            fprefix = os.path.join(dirout, fprefix)
        if H.nested:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + "_nest.nc"
        else:
            ncfname = fprefix + "%s%s" % (H.ibdate, H.ibtime) + ".nc"

    cache_size = 16 * H.numxgrid * H.numygrid * H.numzgrid

    print("About to create new netCDF4 file: '%s'" % ncfname)
    ncid = nc.Dataset(ncfname, 'w', chunk_cache=cache_size)
    write_metadata(H, command, ncid)
    iout = write_header(H, ncid, wetdep, drydep, write_releases)
    write_variables(H, ncid, wetdep, drydep, iout, write_releases, releases)
    ncid.close()
    return ncfname


def main():
    """Parses the passed command line arguments.

    The passed arguments will be used to create the netCDF4 file.
    """
    import argparse
    global MIN_SIZE, COMPLEVEL

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--nested", action="store_true",
        help="Use a nested output.",
        )
    parser.add_argument(
        "-W", "--wetdep", action="store_true",
        help="Write wet depositions in the netCDF4 file.",
        )
    parser.add_argument(
        "-D", "--drydep", action="store_true",
        help="Write dry depositions into the netCDF4 file.",
        )
    parser.add_argument(
        "-d", "--dirout",
        help=("The dir where the netCDF4 file will be created."
              "If not specified, then the fddir/.. is used.")
        )
    parser.add_argument(
        "-o", "--outfile",
        help=("The complete path for the output file."
              "This overrides the --dirout flag.")
        )
    parser.add_argument(
        "-C", "--command-path",
        help=("The path for the associated COMMAND file."
              "If not specified, then the fddir/../options/COMMAND is used.")
        )
    parser.add_argument(
        "-R", "--releases-path",
        help=("The path for the associated RELEASES file."
              "If not specified, then the fddir/../options/RELEASES is used.")
        )
    parser.add_argument(
        "-r", "--dont-write-releases", action="store_true",
        help=("Don't write release point information.")
        )
    parser.add_argument(
        "fddir", nargs="?",
        help="The directory where the FLEXDATA output files are."
        )
    parser.add_argument(
        "--min-size", dest="min_size", action="store_true",
        help=("Do not write redundant fields (orographry) so as to reduce "
              "netCDF4 file size.")
        )
    parser.add_argument(
        "--complevel", type=int, default=9,
        help="Compression level for the netCDF4 file."
        )

    args = parser.parse_args()

    MIN_SIZE = args.min_size
    COMPLEVEL = args.complevel

    if args.fddir is None:
        # At least the FLEXDATA output dir is needed
        parser.print_help()
        sys.exit(1)

    ncfname = create_ncfile(args.fddir, args.nested, args.wetdep, args.drydep,
                            args.command_path, args.releases_path,
                            not args.dont_write_releases,
                            args.dirout, args.outfile)
    print("New netCDF4 file is available in: '%s'" % ncfname)


if __name__ == '__main__':
    main()
