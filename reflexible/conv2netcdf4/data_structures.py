"""
Definition of the different data structures in reflexible.
"""

import datetime
import struct
import traceback

import numpy as np

import reflexible.conv2netcdf4
from .helpers import closest


class BinaryFile(object):

    """

    BinaryFile: A class for accessing data to/from large binary files

    The data is meant to be read/write sequentially from/to a binary file.
    One can request to read a piece of data with a specific type and shape
    from it.  Also, it supports the notion of Fortran and C ordered data,
    so that the returned data is always well-behaved (C-contiguous and
    aligned).

    This class is seeking capable.

    :Author:   Francesc Alted
    :Contact:  faltet@pytables.org
    :Created:  2010-03-18
    :Acknowledgment: Funding for the development of this code is provided
         through the Norwegian Research Council VAUUAV project #184724, 2010

    """

    # Common types whose conversion can be accelerated via the struct
    # module
    structtypes = {
        'i1': 'b', 'i2': 'h', 'i4': 'i',
        'f4': 'f', 'f8': 'd',
    }

    def __init__(self, filename, mode="r", order="fortran"):
        """Open the `filename` for write/read binary files.

        The `mode` can be 'r', 'w' or 'a' for reading (default),
        writing or appending.  The file will be created if it doesn't
        exist when opened for writing or appending; it will be
        truncated when opened for writing.  Add a '+' to the mode to
        allow simultaneous reading and writing.

        `order` specifies whether the file is is written in 'fortran'
        or 'c' order.
        """
        self.mode = mode + "b"
        self.file = open(filename, mode=self.mode, buffering=1)
        """The file handler."""
        if order not in ['fortran', 'c']:
            raise ValueError("order should be either 'fortran' or 'c'.")
        self.order = order
        """The order for file ('c' or 'fortran')."""

    def read(self, dtype, shape=(1,)):
        """Read an array of `dtype` and `shape` from current position.

        `shape` must be any tuple made of integers or even () for scalars.

        The current position will be updated to point to the end of
        read data.
        """
        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)
        if isinstance(shape, int):
            shape = (shape,)
        if not isinstance(shape, tuple):
            raise ValueError("shape must be a tuple")
        length = dtype.itemsize
        rank = len(shape)
        if rank == 1:
            length *= shape[0]
        elif rank > 1:
            length *= np.array(shape).prod()

        # Correct the shape in case dtype is multi-dimensional
        if shape != (1,):
            shape = shape + dtype.shape
        else:
            shape = dtype.shape
        rank = len(shape)

        if shape in (1, (1,)):
            order = "c"
        else:
            order = self.order

        # Read the data from file
        data = self.file.read(length)
        if len(data) < length:
            raise EOFError("Asking for more data than available in file.")

        # Convert read string into a regular array, or scalar
        dts = dtype.base.str[1:]
        if rank == 0:
            if dts[1] == "S":
                data = str(data)
            elif dts in self.structtypes:
                data = struct.unpack(self.structtypes[dts], data)[0]
        else:
            data = np.ndarray(shape=shape, buffer=data, dtype=dtype.base)
            if rank == 0:
                # Retrieve the scalar out of the 0-dim array
                data = data[()]

        if rank > 1:
            # If original data file is in fortran mode, reverse the
            # shape first
            if order == "fortran":
                shape = [i for i in shape[::-1]]
            data = data.reshape(shape)
            # If original data file is in fortran mode, do a transpose.
            # As the shape was reversed previously, we get the original
            # shape again.
            if self.order == "fortran":
                data = data.transpose().copy()
            # Do an additional copy just in case the array is not
            # well-behaved (i.e., it is not aligned or not contiguous).
            elif not data.flags.behaved:
                data = data.copy()
        return data

    def write(self, arr):
        """Write an `arr` to current position.

        The current position will be updated to point to the end of
        written data.
        """
        # Transpose data if case we need to
        if (self.order == "fortran") != (arr.flags.fortran):
            arr = arr.transpose().copy()
        # Write the data to file
        self.file.write(arr.data)

    def seek(self, offset, whence=0):
        """Move to new file position.

        Argument offset is a byte count.  Optional argument whence
        defaults to 0 (offset from start of file, offset should be >=
        0); other values are 1 (move relative to current position,
        positive or negative), and 2 (move relative to end of file,
        usually negative, although many platforms allow seeking beyond
        the end of a file).  If the file is opened in text mode, only
        offsets returned by tell() are legal.  Use of other offsets
        causes undefined behavior.
        """
        self.file.seek(offset, whence)

    def tell(self):
        "Returns current file position, an integer (may be a long integer)."
        return self.file.tell()

    def flush(self):
        "Flush buffers to file."
        self.file.flush()

    def close(self):
        "End access to file."
        self.file.close()


class Structure(dict, object):

    """ A 'fancy' dictionary that provides 'MatLab' structure-like
    referencing.

    .. warning::
        may be replaced with a pure dict in future release.

    """

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

    def __dir__(self):
        return self.keys()


class Header(Structure):

    """This is the primary starting point for processing FLEXPART output.
    The Header class ( :class:`Structure` ) behaves like a dictionary.
    It contains all the metadata from the simulation run as read from the
    "header" or "header_nest" binary files from the model output.

    This version is using the BinaryFile class rather than FortFlex.

    Usage::

        > H = Header(inputpath)
        > H.keys() #provides a list of keys available

    Returns a dictionary

        H = dictionary like object with all the run metadata. TODO: Fill in keys.


    Arguments

      .. tabularcolumns::  |l|L|

      ==============        ========================================
      keyword               Description [default]
      ==============        ========================================
      path                  path to the run directory
      headerfile            name of the header file if non standard
      readheader_ops        optional dictionary to pass readheader
      ==============        ========================================

    Arguments for readheader_ops

      .. tabularcolumns::  |l|L|

      =============       ========================================
      keyword             Description [default]
      =============       ========================================
      pathname            FLEXPART run output directory
      readp               read release points 0=no, [1]=y
      nested              nested output True or [False]
      version             version of FLEXPART, default = 'V8'
      =============       ========================================


    .. note::
        **This function is in development**

        This function is being developed so that there is no dependence on
        using f2py to compile the FortFlex module. It is working using the
        :class:`BinaryFile`, but is notably slower than :class:`FortFlex`.
        Please report any bugs found.


    """

    def __init__(self, path=None, headerfile=None, version='V8',
                 **readheader_ops):
        """


        """
        try:
            h = reflexible.conv2netcdf4.read_header(path, **readheader_ops)
            self.set_with_dict(h)
            self.lonlat()
            self.version = 'V8'
        except:
            traceback.print_exc()
            raise IOError('''
            Could not set header variables.
            Does the `header` or `header_nest` file exist in path?\n{0}'''.format(path))

    def lonlat(self):
        """ Add longitude and latitude attributes using data from header """
        lons = np.arange(self.outlon0,
                         self.outlon0 + (self.dxout * self.numxgrid),
                         self.dxout)
        lats = np.arange(self.outlat0,
                         self.outlat0 + (self.dyout * self.numygrid),
                         self.dyout)
        self.longitude = lons
        self.latitude = lats

    def read_grid(self, **kwargs):
        """ see :func:`read_grid` """
        self.FD = reflexible.conv2netcdf4.read_grid(self, **kwargs)

    def fill_backward(self, **kwargs):
        """ see :func:`fill_backward` """
        reflexible.conv2netcdf4.fill_grids(self, **kwargs)

    def add_trajectory(self):
        """ see :func:`read_trajectories` """
        self.trajectory = reflexible.conv2netcdf4.read_trajectories(self)

    def closest_dates(self, dateval, fmt=None, take_set=False):
        """ given an iterable of datetimes, finds the closest dates.
            if passed a list, assumes it is a list of datetimes

            if take_set=True, then a set of unique values will be returned.
            This can be used with H.read_grid as the value for time_ret to
            return only the grids matching the array of times.
            See (e.g. `extract_curtain`).
        """

        try:
            vals = [closest(d, self['available_dates_dt']) for d in dateval]
            if take_set:
                return list(set(vals))
            else:
                return vals

        except IOError:
            print('If dateval is iterable, must contain datetimes')

    def closest_date(self, dateval, fmt=None):
        """ given a datestring or datetime, tries to find the closest date.
            if passed a list, assumes it is a list of datetimes

        """

        if isinstance(dateval, str):
            if not fmt:
                if len(dateval) == 8:
                    fmt = '%Y%m%d'
                if len(dateval) == 14:
                    fmt = '%Y%m%d%H%M%S'
                else:
                    raise IOError("no format provided for datestring")
            print("Assuming date format: {0}".format(fmt))
            dateval = datetime.datetime.strptime(dateval, fmt)

        return closest(dateval, self['available_dates_dt'])


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
        """I'm the 'grid' property."""
        return self._grid

    @grid.setter
    def grid(self, value):
        """Example of setter.  Add additional code here if desired."""
        self._grid = value
        self._shape = value.shape
        self._max = value.max()
        self._min = value.min()

    @property
    def gridfile(self):
        """I'm the 'gridfile' property."""
        return self._gridfile

    @gridfile.setter
    def gridfile(self, value):
        self._gridfile = value

    @property
    def itime(self):
        """I'm the 'itime' property."""
        return self._itime

    @itime.setter
    def itime(self, value):
        self._itime = value

    @property
    def timestamp(self):
        """I'm the 'timestamp' property."""
        return self._timestamp

    @timestamp.setter
    def timestamp(self, value):
        self._timestamp = value

    @property
    def species(self):
        """I'm the 'species' property."""
        return self._species

    @species.setter
    def species(self, value):
        self._species = value

    @property
    def rel_i(self):
        """I'm the 'rel_i' property."""
        return self._rel_i

    @rel_i.setter
    def rel_i(self, value):
        self._rel_i = value

    @property
    def spec_i(self):
        """I'm the 'spec_i' property."""
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
        """I'm the 'dry' property."""
        return self._dry

    @dry.setter
    def dry(self, value):
        self._dry = value

    @property
    def slabs(self):
        """I'm the 'slabs' property."""
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
