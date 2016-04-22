# !/usr/bin/env python

"""
Definition of the different data structures in reflexible.
"""
import os
import datetime as dt
import itertools
from collections import Iterable
import glob

from math import (pi, cos, sqrt)
import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr

import reflexible


class Header(object):
    """This is the primary starting point for processing FLEXPART output.

    It contains all the metadata from the simulation run and tries to
    fake the behaviour of the `Header` of former ``pflexible`` package
    (that still lives in the ``reflexible.conv2netcdf4`` subpackage).

    This version is using a netCDF4 file instead of a native FLEXPART
    format.

    Usage::

        > H = Header(inputpath)
        > H.keys()   # provides a list of available attributes

    Parameters
    -----------
      path : string
        The path of the netCDF4 file.

    """

    @property
    def alt_unit(self):
        # XXX this depends on H.kindz, which is not in netCDF4 file (I think)
        return 'unkn.'

    @property
    def area(self):
        return self._gridarea()

    def _gridarea(self):
        """returns an array of area corresponding to each nx,ny,nz

        Usage::

            > area = gridarea(H)


        Returns
            OUT = array area corresponding to nx,ny,nz

        Arguments
            H  = :class:`Header` object from readheader function.

        """

        pih = pi / 180.
        r_earth = 6.371e6
        cosfunc = lambda y: cos(y * pih) * r_earth

        nx = self.numxgrid
        ny = self.numygrid
        outlat0 = self.outlat0
        dyout = self.dyout
        dxout = self.dxout
        area = np.zeros((nx, ny))

        for iy in range(ny):
            ylata = outlat0 + (float(iy) + 0.5) * dyout  # NEED TO Check this, iy since arrays are 0-index
            ylatp = ylata + 0.5 * dyout
            ylatm = ylata - 0.5 * dyout
            if (ylatm < 0 and ylatp > 0):
                hzone = dyout * r_earth * pih
            else:
                # cosfact = cosfunc(ylata)
                cosfactp = cosfunc(ylatp)
                cosfactm = cosfunc(ylatm)
                if cosfactp < cosfactm:
                    hzone = sqrt(r_earth ** 2 - cosfactp ** 2) - sqrt(r_earth ** 2 - cosfactm ** 2)
                else:
                    hzone = sqrt(r_earth ** 2 - cosfactm ** 2) - sqrt(r_earth ** 2 - cosfactp ** 2)

            gridarea = 2. * pi * r_earth * hzone * dxout / 360.
            for ix in range(nx):
                area[ix, iy] = gridarea

        return area

    @property
    def outlon0(self):
        return self.nc.outlon0

    @property
    def outlat0(self):
        return self.nc.outlat0

    @property
    def dxout(self):
        return self.nc.dxout

    @property
    def dyout(self):
        return self.nc.dyout

    @property
    def ibdate(self):
        return self.nc.ibdate

    @property
    def ibtime(self):
        return self.nc.ibtime

    @property
    def iedate(self):
        return self.nc.iedate

    @property
    def ietime(self):
        return self.nc.ietime

    @property
    def loutstep(self):
        return self.nc.loutstep

    @property
    def loutaver(self):
        return self.nc.loutaver

    @property
    def loutsample(self):
        return self.nc.loutsample

    @property
    def lsubgrid(self):
        return self.nc.lsubgrid

    @property
    def lconvection(self):
        return self.nc.lconvection

    @property
    def ind_source(self):
        return self.nc.ind_source

    @property
    def ind_receptor(self):
        return self.nc.ind_receptor

    @property
    def ldirect(self):
        return self.nc.ldirect

    @property
    def iout(self):
        return self.nc.iout

    @property
    def direction(self):
        if self.nc.ldirect < 0:
            return "backward"
        else:
            return "forward"

    @property
    def nspec(self):
        return self.nc.dims['numspec']

    @property
    def species(self):
        l = []
        for i in range(self.nspec):
            if self.iout in (1, 3, 5):
                varname = "spec%03d_mr" % (i + 1)
            if self.iout in (2,):  # XXX what to do with 3?
                varname = "spec%03d_pptv" % (i + 1)
            ncvar = self.nc.variables[varname]
            l.append(ncvar.attrs['long_name'])
        return l

    @property
    def output_unit(self):
        if self.iout in (1, 3, 5):
            varname = "spec001_mr"
        if self.iout in (2,):  # XXX what to do with 3?
            varname = "spec001_pptv"
        ncvar = self.nc.variables[varname]
        return ncvar.attrs['units']

    @property
    def numpoint(self):
        return self.nc.dims['numpoint']

    @property
    def pointspec(self):
        return self.nc.dims['pointspec']

    @property
    def numpointspec(self):
        return self.pointspec

    @property
    def nageclass(self):
        return self.nc.dims['nageclass']

    @property
    def numageclasses(self):
        return self.nageclass

    @property
    def numxgrid(self):
        return self.nc.dims['longitude']

    @property
    def numygrid(self):
        return self.nc.dims['latitude']

    @property
    def numzgrid(self):
        return self.nc.dims['height']

    @property
    def longitude(self):
        return np.arange(self.outlon0,
                         self.outlon0 + (self.dxout * self.numxgrid),
                         self.dxout)

    @property
    def latitude(self):
        return np.arange(self.outlat0,
                         self.outlat0 + (self.dyout * self.numygrid),
                         self.dyout)

    @property
    def available_dates_dt(self):
        loutstep = self.nc.loutstep
        nsteps = self.nc.dims['time']
        if self.nc.ldirect < 0:
            # backward direction
            d = dt.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                     "%Y%m%d%H%M%S")
            return [(d + dt.timedelta(seconds=t))
                    for t in range(loutstep * (nsteps - 1), -loutstep, -loutstep)]
        else:
            # forward direction
            d = dt.datetime.strptime(self.nc.ibdate + self.nc.ibtime,
                                     "%Y%m%d%H%M%S")
            return [(d + dt.timedelta(seconds=t))
                    for t in range(0, loutstep * nsteps, loutstep)]

    @property
    def available_dates(self):
        return [d.strftime("%Y%m%d%H%M%S") for d in self.available_dates_dt]

    @property
    def ireleasestart(self):
        return self.nc.variables['RELSTART'][:]

    @property
    def ireleaseend(self):
        return self.nc.variables['RELEND'][:]

    @property
    def releasestart(self):
        if self.nc.ldirect < 0:
            rel_start = self.ireleasestart[::-1]
            d = dt.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                     "%Y%m%d%H%M%S")
            # note xarray converts netcdf file times to timedelta x64 [ns]
            return [(d + dt.timedelta(seconds=int(t) * 10e-9)) for t in rel_start]
        else:
            rel_start = self.ireleasestart[:]
            d = dt.datetime.strptime(self.nc.ibdate + self.nc.ibtime,
                                     "%Y%m%d%H%M%S")
            return [(d + dt.timedelta(seconds=int(t) * 10e-9)) for t in rel_start]

    @property
    def releaseend(self):
        if self.nc.ldirect < 0:
            rel_end = self.ireleaseend[::-1]
            d = dt.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                     "%Y%m%d%H%M%S")
            # note xarray converts netcdf file times to timedelta x64 [ns]
            return [(d + dt.timedelta(seconds=int(t) * 10e-9)) for t in rel_end]
        else:
            rel_end = self.ireleaseend[:]
            d = dt.datetime.strptime(self.nc.ibdate + self.nc.ibtime,
                                     "%Y%m%d%H%M%S")
            return [(d + dt.timedelta(seconds=int(t) * 10e-9)) for t in rel_end]

    @property
    def releasetimes(self):
        return [b - ((b - a) / 2)
                for a, b in zip(self.releasestart, self.releaseend)]

    @property
    def ORO(self):
        if 'ORO' in self.nc.variables:
            return self.nc.variables['ORO'][:].T
        else:
            return None

    @property
    def outheight(self):
        return self.nc.variables['height'][:].T

    @property
    def Heightnn(self):

        outheight = self.outheight[:]
        if self.ORO is not None:
            oro = self.ORO[:]
            Heightnn = outheight + oro
        else:
            Heightnn = outheight
        return Heightnn

    @property
    def zpoint1(self):
        return self.nc.variables['RELZZ1'][:].T

    @property
    def zpoint2(self):
        return self.nc.variables['RELZZ2'][:].T

    def __getitem__(self, key):
        return getattr(self, key)

    def keys(self):
        not_listed = ["keys", "fill_backward", "add_trajectory"]
        return [k for k in dir(self)
                if not k.startswith('_') and k not in not_listed]

    def fill_grids(self):
        return self.C

    def add_trajectory(self):
        """ see :func:`read_trajectories` """
        self.trajectory = reflexible.read_trajectories(self)

    @property
    def options(self):
        return {'readp': None,
                'nested': self.nested,
                'absolute_path': self.absolute_path}

    @property
    def FD(self):
        return FD(self.nc, self.nspec, self.species, self.available_dates,
                  self.direction, self.iout)

    @property
    def C(self):
        return C(self.nc, self.releasetimes, self.species, self.available_dates,
                 self.direction, self.iout, self.Heightnn, self.FD)

    def __init__(self, path=None, nested=False, absolute_path=True):
        if absolute_path:
            files = [path]
        else:
            # print("Warning assuming files have .nc extension")
            files = glob.glob(os.path.join(path, '*.nc'))

        self.nested = nested
        self.absolute_path = absolute_path

        # check for nested or not, assumes only two nc files in 
        # output directory
        if nested:
            ncfile = [d for d in files if 'nest' in d][0]
        else:
            ncfile = [d for d in files if '_nest' not in d][0]

        self.ncfile = ncfile
        self.fp_path = os.path.split(ncfile)[0]
        self.nc = xr.open_dataset(ncfile)


class FD(object):
    """Class that contains FD data indexed with (spec, date)."""

    def __init__(self, nc, nspec, species, available_dates, direction, iout):
        self.nc = nc
        self.nspec = nspec
        self.species = species
        self.available_dates = available_dates
        self.grid_dates = available_dates
        self.direction = direction
        self.iout = iout
        self._keys = [(s, k) for s, k in itertools.product(
            range(nspec), available_dates)]

    @property
    def keys(self):
        return self._keys

    def __getitem__(self, item):
        nspec, date = item
        idate = self.available_dates.index(date)
        if self.iout in (1, 3, 5):
            varname = "spec%03d_mr" % (nspec + 1)
        if self.iout in (2,):  # XXX what to do with the 3 case?
            varname = "spec%03d_pptv" % (nspec + 1)
        fdc = FDC()
        fdc.data_cube = self.nc.variables[varname][:, :, idate, :, :, :].T
        fdc.itime = self.nc.variables['time'][idate]
        fdc.timestamp = dt.datetime.strptime(
            self.available_dates[idate], "%Y%m%d%H%M%S")
        fdc.spec_i = nspec
        if self.direction == "forward":
            fdc.rel_i = 0
        else:
            fdc.rel_i = 'k'
        fdc.species = self.species
        # fdc.wet  # TODO
        # fdc.dry  # TODO
        return fdc


class C(object):
    """Class that contains C data indexed with (spec, date)."""

    def __init__(self, nc, releasetimes, species, available_dates,
                 direction, iout, Heightnn, FD):
        self.nc = nc
        self.nspec = nc.dims['numspec']
        self.pointspec = nc.dims['pointspec']
        self.releasetimes = releasetimes
        self.species = species
        self.available_dates = available_dates
        self.direction = direction
        self.iout = iout
        self.Heightnn = Heightnn
        self._FD = FD
        self._keys = [(s, k) for s, k in itertools.product(
            range(self.nspec), range(self.pointspec))]

    @property
    def keys(self):
        return self._keys()

    def __dir__(self):
        """ necessary for Ipython tab-completion """
        return self._keys

    def __iter__(self):
        return iter(self._keys)

    def __getitem__(self, item):
        """
        Calculates the 20-day sensitivity at each release point.

        This will cycle through all available_dates and create the filled
        array for each k in pointspec.

        Parameters
        ----------
        item : tuple
            A 2-element tuple specifying (nspec, pointspec)

        Return
        ------
        FDC instance
            An instance with grid, timestamp, species and other properties.

        Each element in the dictionary is a 3D array (x,y,z) for each species,k

        """
        assert type(item) is tuple and len(item) == 2
        nspec, pointspec = item
        assert type(nspec) is int and type(pointspec) is int

        if self.direction == 'backward':
            c = FDC()
            c.itime = None
            c.timestamp = self.releasetimes[pointspec]
            c.species = self.species[nspec]
            c.gridfile = 'multiple'
            c.rel_i = pointspec
            c.spec_i = nspec

            # read data grids and attribute/sum sensitivity
            if self.iout in (1, 3, 5):
                varname = "spec%03d_mr" % (nspec + 1)
            if self.iout in (2,):  # XXX what to do with the 3 case?
                varname = "spec%03d_pptv" % (nspec + 1)
            specvar = self.nc.variables[varname][0, pointspec, :, :, :, :]
            if True:
                # changed to use direct summation
                c.data_cube = specvar
                c.time_integrated = np.sum(c.data_cube, axis=0).T
                c.total_column = np.sum(c.time_integrated, axis=2)
                c.foot_print = c.time_integrated[:, :, 0]

            else:
                # Same than the above, but it comsumes more memory
                # Just let it here for future reference
                c.data_cube = specvar[:, :, :, :, pointspec, :].sum(axis=(-2, -1))
            c.slabs = get_slabs(self.Heightnn, c.time_integrated)
        else:
            # forward direction
            FD = self._FD
            d = FD.grid_dates[pointspec]
            c = FD[(nspec, d)]
            c.slabs = get_slabs(self.Heightnn, c.data_cube)

        return c


def get_slabs(Heightnn, grid):
    """Preps grid for plotting.

    Arguments
    ---------
    Heightnn : numpy array
      Height (outheight + topography).
    grid : numpy array
      A grid from the FLEXPARTDATA.

    Returns
    -------
    dictionary
      dictionary of rank-2 arrays corresponding to vertical levels.

    """
    normAreaHeight = True

    slabs = {}
    for i in range(grid.shape[2]):
        if normAreaHeight:
            data = grid[:, :, i] / Heightnn[:, :, i]
        else:
            data = grid[:, :, i]
        slabs[i + 1] = data.T  # XXX why?  something to do with http://en.wikipedia.org/wiki/ISO_6709 ?

    slabs[0] = np.sum(grid, axis=2).T  # XXX why?  something to do with http://en.wikipedia.org/wiki/ISO_6709 ?
    return slabs


class FDC(object):
    """Data container for FD and C grids."""

    def __init__(self):
        self._keys = [
            'grid', 'gridfile', 'itime', 'timestamp', 'species', 'rel_i',
            'spec_i', 'dry', 'wet', 'slabs', 'shape', 'max', 'min']
        for key in self._keys:
            setattr(self, "_" + key, None)

    def keys(self):
        return self._keys

    @property
    def grid(self):
        return self._grid

    @grid.setter
    def grid(self, value):
        self._grid = value
        self._shape = value.shape
        self._max = value.max()
        self._min = value.min()

    @property
    def gridfile(self):
        return self._gridfile

    @gridfile.setter
    def gridfile(self, value):
        self._gridfile = value

    @property
    def itime(self):
        return self._itime

    @itime.setter
    def itime(self, value):
        self._itime = value

    @property
    def timestamp(self):
        return self._timestamp

    @timestamp.setter
    def timestamp(self, value):
        self._timestamp = value

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, value):
        self._species = value

    @property
    def rel_i(self):
        return self._rel_i

    @rel_i.setter
    def rel_i(self, value):
        self._rel_i = value

    @property
    def spec_i(self):
        return self._spec_i

    @spec_i.setter
    def spec_i(self, value):
        self._spec_i = value

    @property
    def wet(self):
        """I'm the 'wet' property."""
        return self._wet

    @wet.setter
    def wet(self, value):
        self._wet = value

    @property
    def dry(self):
        return self._dry

    @dry.setter
    def dry(self, value):
        self._dry = value

    @property
    def slabs(self):
        return self._slabs

    @slabs.setter
    def slabs(self, value):
        self._slabs = value

    # Read-only properties
    @property
    def shape(self):
        return self._shape

    @property
    def max(self):
        return self.time_integrated.max()

    @property
    def min(self):
        return self.time_integrated.min()


class Command(object):
    """ General COMMAND input for Flexpart 

    #TODO: use properties ??
    """

    def __init__(self, **options):

        self._OPTIONS = {

            'IBDATE': [None, '''String simulation date start'''],
            'IBTIME': [None, '''string, simulation time start'''],
            'IEDATE': [None, '''string, simulation date end'''],
            'IETIME': [None, '''string, simulation time end'''],
            'LDIRECT': [1, '''Simulation direction, 1 for forward, -1 for backward in time'''],
            'LOUTSTEP': [10800, '''Average concentrations are calculated every SSSSS seconds.'''],
            'LOUTAVER': [10800, '''The average concentrations are time averages of SSSSS seconds
                duration. If SSSSS is 0, instantaneous concentrations are outputted.'''],
            'LOUTSAMPLE': [900,
                           '''The concentrations are sampled every SSSSS seconds to calculate  the time average concentration. This period must be shorter than the averaging time.'''],
            'ITSPLIT': [999999999,
                        '''Time constant for particle splitting. Particles are split into two after S SSSS seconds, 2xSSSSS seconds, 4xSSSSS seconds, and so on.'''],
            'LSYNCTIME': [900, '''All processes are synchronized with this time interval (lsynctime). Therefore, all other time constants must be multiples of this value.
                Output interval and time average of output must be at least twice lsynctime.'''],
            'CTL': [-5.0,
                    '''CTL must be >1 for time steps shorter than the Lagrangian time scale. If CTL<0, a purely random walk simulation is done'''],
            'IFINE': [4, '''IFINE=Reduction factor for time step used for vertical wind'''],
            'IOUT': [5, '''IOUT determines how the output shall be made: concentration
                  (ng/m3, Bq/m3), mixing ratio (pptv), or both, or plume trajectory mode,
                  or concentration + plume trajectory mode.
                  In plume trajectory mode, output is in the form of average trajectories.'''],
            'IPOUT': [0, '''IPOUT determines whether particle positions are outputted (in addition to the gridded concentrations or mixing ratios) or not.
                0=no output, 1 output every output interval, 2 only at end of the
                simulation'''],
            'LSUBGRID': [1, '''Switch on/off subgridscale terrain parameterization (increase of
                mixing heights due to subgridscale orographic variations)'''],
            'LCONVECTION': [1, '''Switch on/off the convection parameterization'''],
            'LAGESPECTRA': [1,
                            '''Switch on/off the calculation of age spectra: if yes, the file AGECLASSES must be available'''],
            'IPIN': [0,
                     '''If IPIN=1, a file "partposit_end" from a previous run must be available in the output directory. Particle positions are read in and previous simulation is continued. If IPIN=0, no particles from a previous run are used'''],
            'IOUTPUTFOREACHRELEASE': [0, '''Switch on/off writing out each release.'''],
            'IFLUX': [0, '''If IFLUX is set to 1, fluxes of each species through each of the output boxes are calculated. Six fluxes, corresponding to northward, southward,
                eastward, westward, upward and downward are calculated for each grid cell of
                the output grid. The control surfaces are placed in the middle of each
                output grid cell. If IFLUX is set to 0, no fluxes are determined.'''],
            'MDOMAINFILL': [0,
                            '''If MDOMAINFILL is set to 1, the first box specified in file   RELEASES is used as the domain where domain-filling trajectory calculations are to be done. Particles are initialized uniformly distributed (according to the air mass distribution) in that domain at the beginning of the simulation, and are created at the boundaries throughout the simulation perio'''],
            'IND_SOURCE': [1, '''IND_SOURCE switches between different units for concentrations at  the source NOTE that in backward simulations the release of computational particles takes place at   the "receptor" and the sampling of particles at the "source".
                1=mass units (for bwd-runs = concentration)   
                2=mass mixing ratio units'''],
            'IND_RECEPTOR': [1, '''IND_RECEPTOR switches between different units for concentrations at the receptor
                          1=mas s units (concentrations) 
                          2=mas s mixing ratio units'''],
            'MQUASILAG': [0,
                          '''MQUASILAG indicates whether particles shall be numbered consecutively (1) or with their release location number (0). The first option allows tracking of individual particles using the partposit output files'''],
            'NESTED_OUTPUT': [0,
                              '''NESTED_OUTPUT decides whether model output shall be made also   for a nested output field (normally with higher resolution)'''],
            'LINIT_COND': [0, '''For Backward Runs, sets initial conditions: [0]=No, 1=Mass Unit, 2=Mass Mixing'''],
            'SURF_ONLY': [0, '''SURF_ONLY: When set to 1, concentration/emission sensitivity'''],
            'CBLFLAG': [0,
                        '''CBLFLAG: When set to 1, a skewed rather than Gaussian turbulence in  the convective PBL is used.'''],

            ## below here, not actual COMMAND input
            'HEADER': """**********************************************\n\n\n  Input file for   FLEXPART\n\n*********************************************\n\n""",
            'FLEXPART_VER': [10, '''FLEXPART VERSION Used to define format of COM   MAND File'''],
            'SIM_START': [dt.datetime(2000, 01, 01, 00, 00, 00),
                          '''Beginning date and    time of   simulation. Must be given in format YYYYMMDD HHMISS, where YYYY is YEAR, MM  is MONTH, DD is DAY, HH is HOUR, MI is MINUTE and SS is SECOND. Current  version utilizes UTC.'''],
            'SIM_END': [dt.datetime(2000, 02, 01, 00, 00, 00),
                        '''Ending date and time of simulation. Same format as 2'''],
            'AGECLASSES': [[86400 * 30], '''list of ageclasses (seconds) in the simulation'''],
            'RELEASE_SECONDS': [86400, '''duration of the releases in seconds''']
        }

        self._overrides = options

        # set the default options as attributes
        for key, value in self._OPTIONS.iteritems():
            setattr(self, key.lower(), value[0])

        # override the attributes with options
        for key, value in options.iteritems():
            setattr(self, key.lower(), value)

        if self.ibdate is None:
            self.ibdate = self.sim_start.strftime('%Y%m%d')

        if self.iedate is None:
            self.iedate = self.sim_end.strftime('%Y%m%d')

        if self.ibtime is None:
            self.ibtime = self.sim_start.strftime('%H%M%S')

        if self.ietime is None:
            self.ietime = self.sim_end.strftime('%H%M%S')

        self.timedelta = dt.timedelta(seconds=max(self.ageclasses))  # 50 days, time offset with start/end time
        self.release_seconds = dt.timedelta(seconds=self.release_seconds)

    def help(self, key):
        if key.upper() in self._OPTIONS:
            return self._OPTIONS[key.upper()][1]
        else:
            return 'no help available'

    def to_file(self, cfile):
        """ write out the command file """

        if self.ldirect == -1:
            # backward run
            tstart = self.sim_start - self.timedelta
            tend = self.sim_end + self.release_seconds
        elif self.ldirect == 1:
            tstart = self.sim_start
            tend = self.sim_end + self.timedelta + self.release_seconds

        with open(cfile, 'wb') as outf:

            outf.write('&COMMAND\n')
            outf.write(' LDIRECT={0},\n'.format(self.ldirect))
            outf.write(' IBDATE=  {0},\n'.format(tstart.strftime('%Y%m%d')))
            outf.write(' IBTIME=  {0},\n'.format(tstart.strftime('%H%M%S')))
            outf.write(' IEDATE=  {0},\n'.format(tend.strftime('%Y%m%d')))
            outf.write(' IETIME=  {0},\n'.format(tend.strftime('%H%M%S')))
            outf.write(' LOUTSTEP=             {0},\n'.format(self.loutstep))
            outf.write(' LOUTAVER=             {0},\n'.format(self.loutaver))
            outf.write(' LOUTSAMPLE=           {0},\n'.format(self.loutsample))
            outf.write(' ITSPLIT=              {0},\n'.format(self.itsplit))
            outf.write(' LSYNCTIME=            {0},\n'.format(self.lsynctime))
            outf.write(' CTL=                  {0},\n'.format(self.ctl))
            outf.write(' IFINE=                {0},\n'.format(self.ifine))
            outf.write(' IOUT=                 {0},\n'.format(self.iout))
            outf.write(' IPOUT=                {0},\n'.format(self.ipout))
            outf.write(' LSUBGRID=             {0},\n'.format(self.lsubgrid))
            outf.write(' LCONVECTION=          {0},\n'.format(self.lconvection))
            outf.write(' LAGESPECTRA=          {0},\n'.format(self.lagespectra))
            outf.write(' IPIN=                 {0},\n'.format(self.ipin))
            outf.write(' IOUTPUTFOREACHRELEASE={0},\n'.format(self.ioutputforeachrelease))
            outf.write(' IFLUX=                {0},\n'.format(self.iflux))
            outf.write(' MDOMAINFILL=          {0},\n'.format(self.mdomainfill))
            outf.write(' IND_SOURCE=           {0},\n'.format(self.ind_source))
            outf.write(' IND_RECEPTOR=         {0},\n'.format(self.ind_receptor))
            outf.write(' MQUASILAG=            {0},\n'.format(self.mquasilag))
            outf.write(' NESTED_OUTPUT=        {0},\n'.format(self.nested_output))
            outf.write(' LINIT_COND=           {0},\n'.format(self.linit_cond))
            outf.write(' SURF_ONLY=            {0},\n'.format(self.surf_only))
            outf.write(' CBLFLAG=              {0},\n/\n'.format(self.cblflag))
            outf.close()


class Ageclass(object):
    """ General COMMAND input for Flexpart 

    """

    def __init__(self, ageclasses=[86400 * 50]):
        self._keys = ['ageclasses']
        for key in self._keys:
            setattr(self, "_" + key, ageclasses)

    def keys(self):
        return self._keys

    def __dir__(self):
        """ necessary for Ipython tab-completion """
        return self._keys

    def __iter__(self):
        return iter(self._keys)

    @property
    def ageclasses(self):
        return self._ageclasses

    @ageclasses.setter
    def ageclasses(self, value):
        self._ageclasses = value

    def to_file(self, acfile):
        """ write out an ageclasses files """
        # get number of AGECLASSES
        assert isinstance(self.ageclasses, Iterable), 'ageclasses argument must be an iterable of seconds'
        nageclass = len(self.ageclasses)

        with open(acfile, 'w') as outf:
            # WRITE TO FILE namelist style
            outf.write('&AGECLASS\n')
            outf.write(' NAGECLASS={0},\n'.format(nageclass))
            outf.write(' LAGE=')
            for i in range(len(self.ageclasses)):
                outf.write(' {0},'.format(self.ageclasses[i]))

            outf.write('\n/\n')
            outf.close()

            # print('WRITE AGECLASSES: wrote: {0} \n'.format(acfile))


class Release():
    """ subclass of a pandas dataframe to allow for some special properties
    and methods.

    The pandas object is presently set up only to handle mass of all nspec equal
    """

    def __init__(self, data):

        self.releases = data

    def to_file(self, rfile):
        """ write out all the releases """

        with open(rfile, 'w') as outf:

            self.rel_file = outf
            outf.write('&RELEASES_CTRL\n')
            outf.write(' NSPEC=        {0},\n'.format(self.releases.attrs['nspec']))
            outf.write(' SPECNUM_REL=')
            idx = range(self.releases.attrs['nspec'])
            for i in idx:
                outf.write(' {0},'.format(self.releases.attrs['specnum_rel']))
            outf.write('\n /\n')

            for row in self.releases.iterrows():
                # for some reason itertuples would be better?
                self._write_single_release(row, self.releases.attrs)

    def _write_single_release(self, row, attrs):
        """ a nice exercise would be to create a custom formatter from the pandas
        class types, but requires cython magic. """

        # write out the release to file,
        # assumes it is appending
        outf = self.rel_file
        # during OBuoy processing, lat, lon in index caused problems
        # set only time as index, adjusted below
        # t = row[0][2] # seems row returns the MultiIndex as a tuple, 'time' is [2]
        # lat = row[0][0]
        # lon = row[0][1]
        # import pdb
        # pdb.set_trace()
        t = row[0]
        d = row[1]
        lat = d['lat']
        lon = d['lon']

        # get the attrs
        seconds = attrs['release_seconds']
        nspec = attrs['nspec']
        dx = attrs['dx']
        dy = attrs['dy']
        name = attrs['name']

        # vars we can calculate
        lon1 = lon - dx
        lon2 = lon + dx
        lat1 = lat - dy
        lat2 = lat + dy
        rel_ident = '{0}_{1}_{2}|{3}'.format(name, t.strftime('%Y%j'), lat, lon)
        # print(t.strftime('%Y%m%d'), d.lat1, d.lon1)
        t2 = t + dt.timedelta(seconds=seconds)

        outf.write('&RELEASE\n')
        outf.write(t.strftime(' IDATE1=  %Y%m%d,\n'))
        outf.write(t.strftime(' ITIME1=  %H%M%S,\n'))
        outf.write(t2.strftime(' IDATE2=  %Y%m%d,\n'))
        outf.write(t2.strftime(' ITIME2=  %H%M%S,\n'))
        outf.write(' LON1=    {0:3.4f},\n'.format(lon1))  # LON values -180 180
        outf.write(' LON2=    {0:3.4f},\n'.format(lon2))
        outf.write(' LAT1=    {0:3.4f},\n'.format(lat1))  # LAT values -90 90
        outf.write(' LAT2=    {0:3.4f},\n'.format(lat2))
        outf.write(' Z1=      {0:f},\n'.format(d.z1))  # altitude in meters
        outf.write(' Z2=      {0:f},\n'.format(d.z2))
        outf.write(' ZKIND=   {0:d},\n'.format(int(d.zkind)))  # M)ASL= MAG=
        outf.write(' MASS=')

        for i in range(nspec):
            outf.write('    {:8.4f},'.format(d.mass))

        outf.write('\n PARTS=   {0:d},\n'.format(int(d.parts)));
        outf.write(' COMMENT= "{0}"\n /\n'.format(rel_ident))


class ReleasePoint(object):
    """ An individual release entity (point, line, or area)

    """

    def __init__(self, idt1=None, idt2=None, **options):
        """ A release point is a single entity within a Release.

        Each point has the attributes of the release point.

        idt1 and idt2 are datetime objects that will override anything provided
        in idate1, itime1, idate2, or itime2
        """

        self._OPTIONS = {
            'idate1': ['20010101', '''YYYYMMDD begin date of release '''],
            'itime1': ['000000', '''YYYYMMDD begin time of release '''],
            'idate2': ['20010201', '''YYYYMMDD end date of release '''],
            'itime2': ['000000', '''YYYYMMDD end time of release '''],
            'lon1': [120, '''lowerleft Longitude'''],
            'lon2': [130, '''upperright Longitude'''],
            'lat1': [55, '''lowerleft Latitude'''],
            'lat2': [60, '''upperright Latitude'''],
            'z1': [20, '''lower boundary of release point (m)'''],
            'z2': [100, '''upper z-level of release point (m)'''],
            'zkind': [3, ''' 1 for m above ground, 2 for m above sea level, 3 for pressure in hPa'''],
            'mass': [[1.0], '''mass of species'''],
            'nspec': [1, '''number of species'''],
            'parts': [50000, '''total number of particles in release'''],
            'specnum_rel': [(22,), '''tuple of species number id'''],
            'run_ident': ['comment', '''character*40 comment''']
        }

        self._overrides = options

        for key, value in self._OPTIONS.iteritems():
            setattr(self, key.lower(), value[0])

        for key, value in options.iteritems():
            setattr(self, key.lower(), value)

        if idt1:
            assert isinstance(idt1, dt.datetime)
            setattr(self, 'idate1', idt1.strftime('%Y%m%d'))
            setattr(self, 'itime1', idt1.strftime('%H%M%S'))

        if idt2:
            assert isinstance(idt2, dt.datetime)
            setattr(self, 'idate2', idt2.strftime('%Y%m%d'))
            setattr(self, 'itime2', idt2.strftime('%H%M%S'))

    def help(self, key):
        if key in self._OPTIONS:
            return self._OPTIONS[key.lower()][1]
        else:
            return 'no help available'

    def _write_single_release(self, rfile):
        """ write out the release to file, assumes it is appending """

        if isinstance(rfile, str):

            outf = open(rfile, 'w')
        else:
            outf = rfile

        outf.write('&RELEASE\n')
        outf.write(' IDATE1=  {0},\n'.format(self.idate1))
        outf.write(' ITIME1=  {0},\n'.format(self.itime1))
        outf.write(' IDATE2=  {0},\n'.format(self.idate2))
        outf.write(' ITIME2=  {0},\n'.format(self.itime2))
        outf.write(' LON1=    {0},\n'.format(self.lon1))  # LON values -180 180
        outf.write(' LON2=    {0},\n'.format(self.lon2))
        outf.write(' LAT1=    {0},\n'.format(self.lat2))  # LAT values -90 90
        outf.write(' LAT2=    {0},\n'.format(self.lat2))
        outf.write(' Z1=      {0},\n'.format(self.z1))  # altitude in meters
        outf.write(' Z2=      {0},\n'.format(self.z2))
        outf.write(' ZKIND=   {0},\n'.format(self.zkind))  # M)ASL= MAG=
        outf.write(' MASS=')
        for j in range(self.nspec):
            outf.write('    {%8.2d},'.format(self.mass))

        outf.write('\n PARTS=   {0},\n'.format(self.parts));
        outf.write(' COMMENT= "{0}"\n /\n'.format(self.run_ident))


class Releases(object):
    """ class for a group of releases """

    def __init__(self, release_points):
        ''' takes an iterable of `ReleasePoint` classes '''
        self.releases = release_points
        # hack, take the first release to get some general info
        self.nspec = release_points[0].nspec
        self.specnum_rel = release_points[0].specnum_rel
        self.nreleases = len(release_points)

    def to_file(self, rfile):
        """ write out all the releases """

        with open(rfile, 'w') as outf:

            outf.write('&RELEASES_CTRL\n')
            outf.write(' NSPEC=        {0},\n'.format(self.nspec))
            outf.write(' SPECNUM_REL=')
            for i in range(self.nspec):
                outf.write(' {0},   '.format(self.specnum_rel[i]))
            outf.write('\n /\n')

            for r in self.releases:
                r._write_single_release(outf)

            outf.close()


class Trajectory(dict, object):
    """ a dictionary-like container for the trajectories.
    This should become a pandas dataframe in the future. """

    FutureWarning('to be a pandas.DataFrame in the future')

    def __getattr__(self, attr):
        # Fake a __getstate__ method that returns None
        if attr == "__getstate__":
            return lambda: None
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def __dir__(self):
        """ necessary for Ipython tab-completion """
        return self.keys()

    def set_with_dict(self, D):
        """ set attributes with a dict """
        for k in D.keys():
            self.__setattr__(k, D[k])
