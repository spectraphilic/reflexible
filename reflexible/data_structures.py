"""
Definition of the different data structures in reflexible.
"""

import datetime
import itertools

import numpy as np
import netCDF4 as nc


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
    def direction(self):
        if self.nc.ldirect < 0:
            return "backward"
        else:
            return "forward"

    @property
    def nspec(self):
        return len(self.nc.dimensions['numspec'])

    @property
    def species(self):
        l = []
        for i in range(self.nspec):
            varname = "spec%03d_pptv" % (i + 1)   # XXX check this with IOUT
            ncvar = self.nc.variables[varname]
            l.append(ncvar.long_name)
        return l

    @property
    def numpoint(self):
        return len(self.nc.dimensions['numpoint'])

    @property
    def numpointspec(self):
        return len(self.nc.dimensions['pointspec'])

    @property
    def numageclasses(self):
        return len(self.nc.dimensions['nageclass'])

    @property
    def numxgrid(self):
        return len(self.nc.dimensions['longitude'])

    @property
    def numygrid(self):
        return len(self.nc.dimensions['latitude'])

    @property
    def numzgrid(self):
        return len(self.nc.dimensions['height'])

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
    def available_dates(self):
        loutstep = self.nc.loutstep
        nsteps = len(self.nc.dimensions['time'])
        d = datetime.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                       "%Y%m%d%H%M%S")
        l = [(d + datetime.timedelta(seconds=t)).strftime("%Y%m%d%H%M%S")
             for t in range(loutstep * (nsteps - 1), -loutstep, -loutstep)]
        return l

    @property
    def ireleasestart(self):
        return self.nc.variables['RELSTART']

    @property
    def ireleaseend(self):
        return self.nc.variables['RELEND']

    @property
    def releasestart(self):
        rel_start = self.nc.variables['RELSTART'][::-1]
        d = datetime.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                       "%Y%m%d%H%M%S")
        return [(d - datetime.timedelta(seconds=int(t))) for t in rel_start]

    @property
    def releaseend(self):
        rel_end = (self.nc.variables['RELEND'][::-1] +
                   self.nc.variables['RELSTART'][-1] * 2)  # XXX ugly workaround
        d = datetime.datetime.strptime(self.nc.iedate + self.nc.ietime,
                                       "%Y%m%d%H%M%S")
        return [(d - datetime.timedelta(seconds=int(t))) for t in rel_end]

    @property
    def releasetimes(self):
        return [b - ((b - a) / 2)
                for a, b in zip(self.releasestart, self.releaseend)]

    @property
    def ORO(self):
        if 'ORO' in self.nc.variables:
            return self.nc.variables['ORO']
        else:
            return None

    @property
    def outheight(self):
        return self.nc.variables['height']

    @property
    def Heightnn(self):
        nx, ny, nz = (self.numxgrid, self.numygrid, self.numzgrid)
        Heightnn = np.zeros((nx, ny, nz), np.float)
        if self.ORO is not None:
            oro = self.ORO[:].T
        else:
            oro = None
        for ix in range(nx):
            if oro is not None:
                Heightnn[ix, :, 0] = oro[ix, :]  # XXX this value is overwritten later on.  Have a look John.
                for iz in range(nz):
                    Heightnn[ix, :, iz] = self.outheight[iz] + oro[ix, :]
            else:
                Heightnn[ix, :, :] = self.outheight[:]
        return Heightnn

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
        self.trajectory = reflexible.conv2netcdf4.read_trajectories(self)

    @property
    def FD(self):
        return FD(self.nc, self.nspec, self.species, self.available_dates,
                  self.direction)

    @property
    def C(self):
        return C(self.nc, self.releasetimes, self.species,
                 self.available_dates, self.direction, self.Heightnn)

    def __init__(self, path=None):
        self.nc = nc.Dataset(path, 'r')


class FD(object):
    """Class that contains FD data indexed with (spec, date)."""

    def __init__(self, nc, nspec, species, available_dates, direction):
        self.nc = nc
        self.nspec = nspec
        self.species = species
        self.available_dates = available_dates
        self.direction = direction
        self._keys = [(s, k) for s, k in itertools.product(
            range(nspec), available_dates)]

    @property
    def keys(self):
        return self._keys()

    def __getitem__(self, item):
        nspec, date = item
        idate = self.available_dates.index(date)
        varname = "spec%03d_pptv" % (nspec + 1)   # XXX check this with IOUT
        fdc = FDC()
        fdc.grid = self.nc.variables[varname][:, :, idate, :, :, :].T
        fdc.itime = self.nc.variables['time'][idate]
        fdc.timestamp = datetime.datetime.strptime(
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
                 direction, Heightnn):
        self.nc = nc
        # self._FD = FD
        self.nspec = len(nc.dimensions['numspec'])
        self.pointspec = len(nc.dimensions['pointspec'])
        self.releasetimes = releasetimes
        self.species = species
        self.available_dates = available_dates
        self.direction = direction
        self.Heightnn = Heightnn
        self._keys = [(s, k) for s, k in itertools.product(
            range(self.nspec), range(self.pointspec))]

    @property
    def keys(self):
        return self._keys()

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
            varname = "spec%03d_pptv" % (nspec + 1)   # XXX check this with IOUT
            specvar = self.nc.variables[varname]
            if False:
                c.grid = np.zeros((
                    len(self.nc.dimensions['longitude']),
                    len(self.nc.dimensions['latitude']),
                    len(self.nc.dimensions['height'])))
                for date in self.available_dates:
                    idate = self.available_dates.index(date)
                    # cycle through all the date grids
                    c.grid += specvar[:, pointspec, idate, :, :, :].T.sum(axis=-1)
            else:
                # Same than the above, but it probably comsumes more memory
                c.grid = specvar[:, pointspec, :, :, :, :].T.sum(axis=(3, 4))
        else:
            # forward direction
            FD = self._FD
            d = FD.grid_dates[pointspec]
            c = FD[(nspec, d)]

        # add total column
        c.slabs = get_slabs(self.Heightnn, c.grid)

        return c


def get_slabs(Heightnn, G):
    """Preps grid for plotting.

    Accepts an 3D or 5D GRID from readgrid with optional index.

    Arguments
    ---------
    Heightnn : numpy array
      Height (outheight + topography)
    G : numpy array
      A grid from the FLEXPARTDATA.

    Returns
    -------
    dictionary
      dictionary of rank-2 arrays corresponding to vertical levels.

    """
    index = 0
    nageclass = 0   # XXX assume nageclass equal to 0
    normAreaHeight = True

    Slabs = {}
    grid_shape = G.shape
    if len(grid_shape) is 5:
        g = G[:, :, :, index, nageclass]
    elif len(grid_shape) is 3:
        g = G
    else:
        raise(ValueError, "len(grid_shape) cannot be 4")

    for i in range(g.shape[2]):
        # first time sum to create Total Column
        if i == 0:
            TC = np.sum(g, axis=2).T
        if normAreaHeight:
            data = g[:, :, i] / Heightnn[:, :, i]
        else:
            data = g[:, :, i]
        Slabs[i + 1] = data.T

    Slabs[0] = TC
    return Slabs


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
        return self._max

    @property
    def min(self):
        return self._min
