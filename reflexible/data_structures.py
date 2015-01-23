"""
Definition of the different data structures in reflexible.
"""

import datetime

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
        > dir(H)   # provides a list of available attributes

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
    def numxgrid(self):
        return len(self.nc.dimensions['longitude'])

    @property
    def numygrid(self):
        return len(self.nc.dimensions['latitude'])

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
            for t in range(loutstep * (nsteps - 1),
                           -loutstep, -loutstep)]
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

    def __init__(self, path=None):
        self.nc = nc.Dataset(path, 'r')

    def __getitem__(self, key):
        return getattr(self, key)

    def keys(self):
        not_listed = ["keys", "fill_backward", "add_trajectory"]
        return [k for k in dir(self)
                if not k.startswith('_') and k not in not_listed]

    def _read_grid(self, **kwargs):
        """ see :func:`read_grid` """
        self.FD = FD = FDC()
        self.FD = reflexible.conv2netcdf4.read_grid(self, **kwargs)

    def fill_backward(self, **kwargs):
        """ see :func:`fill_backward` """
        reflexible.conv2netcdf4.fill_grids(self, **kwargs)

    def add_trajectory(self):
        """ see :func:`read_trajectories` """
        self.trajectory = reflexible.conv2netcdf4.read_trajectories(self)


class FDC(object):
    """ Brief explanation of what represents this class.
    """
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
