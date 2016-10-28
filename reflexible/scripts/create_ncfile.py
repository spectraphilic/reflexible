#!/usr/bin/env python
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
from collections import defaultdict

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

    Only 'release_point_names' entry returned.
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


def read_species(options_dir, nspec):
    """Read metadata from SPECIES dir and return it as a dict.

    TODO: maybe the current version is specific of FLEXPART 9.
    """
    species_dir = os.path.join(options_dir, "SPECIES")
    if not os.path.isdir(species_dir):
        warnings.warn(
            "The SPECIES dir cannot be found.  Continuing without it!")
        return {}

    keymap = {
        # "Tracer name": "tname",
        "Species half life": "decay",
        "Wet deposition - A": "weta",
        "Wet deposition - B": "wetb",
        "Dry deposition (gases) - D": "reldiff",
        "Dry deposition (gases) - Henrys const.": "henry",
        "Dry deposition (gases) - f0 (reactivity)": "f0",
        # TODO: what's the equivalent to rho?
        # "Dry deposition (particles) - rho": "",
        "Dry deposition (particles) - dquer": "dquer",
        "Dry deposition (particles) - dsig": "dsigma",
        "Alternative: dry deposition velocity": "dryvel",
        "molweight": "weightmolar",
        "OH Reaction rate at 25 deg, [cm^3/sec]": "ohreact",
        "number of associated specias (neg. none)": "spec_ass",
        "KOA - organic matter air partitioning": "kao",
        }

    dspec = defaultdict(list)
    for spec in range(1, nspec+1):
        path = os.path.join(species_dir, "SPECIES_%03d" % spec)
        with open(path) as fspec:
            for line in fspec:
                for sstring in keymap:
                    if sstring in line:
                        try:
                            if keymap[sstring] == "spec_ass":
                                value = int(line.split()[0])
                            else:
                                value = float(line.split()[0])
                        except ValueError:
                            # Probably this line does not have a value
                            pass
                        dspec[keymap[sstring]].append(value)
    return dspec


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
        ncid.itsplit = command['ITSPLIT'] if 'ITSPLIT' in command else command['T_PARTSPLIT']
        ncid.linit_cond = command['LINIT_COND']
        ncid.lsynctime = command['LSYNCTIME'] if 'LSYNCTIME' in command else command['SYNC']
        ncid.ctl = command['CTL']
        ncid.ifine = command['IFINE']
        ncid.iout = command['IOUT']
        ncid.ipout = command['IPOUT']
        ncid.lagespectra = command['LAGESPECTRA']
        ncid.ipin = command['IPIN']
        ncid.ioutputforeachrelease = command['IOUTPUTFOREACHRELEASE'] \
            if 'IOUTPUTFOREACHRELEASE' in command else command['OUTPUTFOREACHRELEASE']
        ncid.iflux = command['IFLUX']
        ncid.mdomainfill = command['MDOMAINFILL']
        ncid.mquasilag = command['MQUASILAG']
        ncid.nested_output = command['NESTED_OUTPUT']
        # This is a new option in V9.2
        ncid.surf_only = command.get('SURF_ONLY', 0)
    else:
        # IOUT is not available (no COMMAND file) so guess the value:
        unit_i = UNITS.index(H.unit) + 1
        iout = (unit_i) + (H.nested * 5)
        ncid.iout = iout


def write_header(H, ncid, wetdep, drydep, write_releases, species):
    """Create netCDF4 dimensions and variables.

    Create the netCDF4 variables (and the required dimensions) that will be
    stored in the netCDF4 file.

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

    """
    iout = ncid.iout

    # Create dimensions

    # time
    ncid.createDimension('time', None)
    adate, atime = str(H.ibdate), str(H.ibtime).zfill(6)
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
        oroID = ncid.createVariable('ORO', 'i4', ('latitude', 'longitude'),
                                    chunksizes=(H.numygrid, H.numxgrid),
                                    zlib=True, complevel=COMPLEVEL)
        oroID.standard_name = 'surface altitude'
        oroID.long_name = 'outgrid surface altitude'
        oroID.units = 'm'

    units = output_units(ncid)

    # Concentration output, wet and dry deposition variables (one per species)
    # Variables for concentration ouput have dimensions given by dIDs
    # Variables for wet and dry deposition have dimensions given by depdIDs
    dIDs = (
        'longitude', 'latitude', 'height', 'time', 'pointspec', 'nageclass')
    depdIDs = ('longitude', 'latitude', 'time', 'pointspec', 'nageclass')
    # Reverse dims because we want Fortran order on-disk (FLEXPART convention)
    dIDs = dIDs[::-1]
    depdIDs = depdIDs[::-1]

    chunksizes = (H.numxgrid, H.numygrid, H.numzgrid, 1, 1, 1)[::-1]
    dep_chunksizes = (H.numxgrid, H.numygrid, 1, 1, 1)[::-1]
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
        if iout in (2,):   # XXX what to do with 3?
            var_name = "spec" + anspec + "_pptv"
            sID = ncid.createVariable(var_name, 'f4', dIDs,
                                      chunksizes=chunksizes,
                                      zlib=True, complevel=COMPLEVEL)
            sID.units = 'pptv'
        sID.long_name = H.species[i]
        sID.decay = species.get("decay", [])
        sID.weightmolar = species.get("weightmolar", [])
        sID.ohreact = species.get("ohreact", [])
        sID.kao = species.get("kao", [])
        sID.vsetaver = species.get("vsetaver", [])
        sID.spec_ass = species.get("spec_ass", np.array([], dtype=np.int_))

        if wetdep:
            var_name = "WD_spec" + anspec
            wdsID = ncid.createVariable(var_name, 'f4', depdIDs,
                                        chunksizes=dep_chunksizes,
                                        zlib=True, complevel=COMPLEVEL)
            wdsID.units = '1e-12 kg m-2'
            wdsID.weta = species.get("weta", [])
            wdsID.wetb = species.get("wetb", [])
            wdsID.weta_in = species.get("weta_in", [])
            wdsID.wetb_in = species.get("wetb_in", [])
            wdsID.wetc_in = species.get("wetc_in", [])
            wdsID.wetd_in = species.get("wetd_in", [])
            wdsID.dquer = species.get("dquer", [])
            wdsID.henry = species.get("henry", [])

        if drydep:
            var_name = "DD_spec" + anspec
            ddsID = ncid.createVariable(var_name, 'f4', depdIDs,
                                        chunksizes=dep_chunksizes,
                                        zlib=True, complevel=COMPLEVEL)
            ddsID.units = '1e-12 kg m-2'
            ddsID.dryvel = species.get("dryvel", [])
            ddsID.reldiff = species.get("reldiff", [])
            ddsID.henry = species.get("henry", [])
            ddsID.f0 = species.get("f0", [])
            ddsID.dquer = species.get("dquer", [])
            ddsID.density = species.get("density", [])
            ddsID.dsigma = species.get("dsigma", [])

    return iout


def write_variables(H, ncid, wetdep, drydep, write_releases, releases):
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
    """
    iout = ncid.iout

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
        ncid.variables['RELXMASS'][:, :] = H.xmass.T
        if "release_point_names" in releases:
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
                conc = ncid.variables[conc_name]
                # (x, y, z, time, pointspec, nageclass) <-
                # (x, y, z, pointspec, nageclass)
                conc[:, :, idt, :, :, :] = fd.grid[:, :, :, np.newaxis, :, :].T
            if iout in (2, 3):         # XXX what to do with the 3 case here?  John?
                conc_name = "spec" + anspec + "_pptv"
                # XXX fill the _pptv here...

            # wet and dry depositions
            # (x, y, time, pointspec, nageclass) <- (x, y, pointspec, nageclass)
            if wetdep:
                wet = ncid.variables["WD_spec" + anspec]
                wet[:, :, idt, :, :] = fd.wet[:, :, np.newaxis, :, :].T
            if drydep:
                dry = ncid.variables["DD_spec" + anspec]
                dry[:, :, idt, :, :] = fd.dry[:, :, np.newaxis, :, :].T


def read_conffiles(filename, fddir, path):
    """Read FLEXPART config files and return a dictionary of metadata."""
    if path is None:
        path = os.path.join(fddir, filename)
    if not os.path.isfile(path):
        warnings.warn(
            "The %s file cannot be found.  Continuing without it!" % path)
        return {}
    try:
        if filename == "COMMAND":
            return read_command(path)
        elif filename == "RELEASES":
            return read_releases(path)
    except IOError:
        warnings.warn(
            "The %s file format is not supported.  "
            "Continuing without it!" % path)
    return {}


def create_ncfile(pathnames, nested, wetdep=False, drydep=False,
                  command_path=None, releases_path=None,
                  write_releases=True,
                  dirout=None, outfile=None):
    """Main function that create a netCDF4 file from a FLEXPART output.

    Parameters
    ----------
    pathnames : string
      the file where the FLEXDATA <options> and <output> are specified.
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
    tuple
      (the path of the netCDF4 file, the options dir, the output dir).
    """

    def get_dir(dir, parent_dir):
        if dir.startswith('/'):
            # Absolute path.  Just keep the last level and append to parent.
            dir = dir[:-1] if dir.endswith('/') else dir
            dir = os.path.join(parent_dir, os.path.basename(dir))
        else:
            dir = os.path.join(parent_dir, dir)
        return dir

    if not os.path.isfile(pathnames):
        raise IOError("pathnames file is not found at '{}'".format(pathnames))
    # Get the <options> and <output> directories
    with open(pathnames) as f:
        options_dir = get_dir(f.readline().strip(), os.path.dirname(pathnames))
        output_dir = get_dir(f.readline().strip(), os.path.dirname(pathnames))
    H = Header(output_dir, nested=nested)

    if H.direction == "forward":
        fprefix = 'grid_conc_'
    else:
        fprefix = 'grid_time_'

    command = read_conffiles("COMMAND", options_dir, command_path)
    releases = read_conffiles("RELEASES", options_dir, releases_path)
    species = read_species(options_dir, H.nspec)

    if outfile:
        # outfile has priority over previous flags
        ncfname = outfile
    else:
        if dirout is None:
            path = os.path.dirname(output_dir)
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
    write_header(H, ncid, wetdep, drydep, write_releases, species)
    write_variables(H, ncid, wetdep, drydep, write_releases, releases)
    ncid.close()
    return (ncfname, options_dir, output_dir)


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
        help=("The dir where the netCDF4 file will be created.  "
              "If not specified, then the <output> dir is used.")
        )
    parser.add_argument(
        "-o", "--outfile",
        help=("The complete path for the output file."
              "This overrides the --dirout flag."))
    parser.add_argument(
        "-C", "--command-path",
        help=("The path for the associated COMMAND file.  "
              "If not specified, then the <options>/COMMAND is used.")
        )
    parser.add_argument(
        "-R", "--releases-path",
        help=("The path for the associated RELEASES file.  "
              "If not specified, then the <options>/RELEASES is used.")
        )
    parser.add_argument(
        "-r", "--dont-write-releases", action="store_true",
        help=("Don't write release point information.")
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
    parser.add_argument(
        "pathnames", nargs="?",
        help="The Flexpart pathnames file stating where options and output are."
        )

    args = parser.parse_args()

    MIN_SIZE = args.min_size
    COMPLEVEL = args.complevel

    if args.pathnames is None:
        # The FLEXDATA pathnames file is mandatory
        parser.print_help()
        sys.exit(1)

    ncfname, options_dir, output_dir = create_ncfile(
        args.pathnames, args.nested, args.wetdep, args.drydep,
        args.command_path, args.releases_path, not args.dont_write_releases,
        args.dirout, args.outfileN)

    print("New netCDF4 file is available in: '%s'" % ncfname)

if __name__ == '__main__':
    main()
