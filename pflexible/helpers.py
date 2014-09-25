########### HELPER FUNCTIONS ##########

import datetime
import struct

from matplotlib.dates import date2num
import matplotlib.image as image
from matplotlib.patches import Ellipse


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


def closest(num, numlist):
    """ returns the index of the *closest* value in a list """
    # check if we're using datetimes
    dates = False
    if isinstance(num, datetime.datetime):
        dates = True
    if dates:
        num = date2num(num)
        assert isinstance(numlist[0], datetime.datetime), \
            "num is date, numlist must be a list of dates"
        numlist = date2num(numlist)

    return (np.abs(numlist - num)).argmin()


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

        > H = pf.Header(inputpath)
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
            h = read_header(path, **readheader_ops)
            self.set_with_dict(h)
            self.lonlat()
            self.version = 'V8'
        except:
            traceback.print_exc()
            raise IOError('''
            Could not set header variables.
            Does the `header` file exist in path?\n{0}'''.format(path))

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
        self.FD = read_grid(self, **kwargs)

    def fill_backward(self, **kwargs):
        """ see :func:`fill_backward` """
        fill_grids(self, add_attributes=True, **kwargs)

    def add_trajectory(self, **kwargs):
        """ see :func:`read_trajectories` """
        self.trajectory = read_trajectories(self)

    def add_fires(self, **kwargs):
        """ uses the :mod:`emissions` module to read the MODIS hotspot data and
        add it to the header class as a 'fires' attribute.

        **This function is only available within UIO.**

        """

        from jfb.pflexpart import emissions as em
        self.fires = None
        for day in self.available_dates_dt:
            # day = day[:8]
            firedata = em.MODIS_hotspot(day)
            daily = firedata.daily
            # pdb.set_trace()
            if daily is None:
                continue
            else:
                if self.fires is None:
                    self.fires = daily
                else:

                    self.fires = np.hstack(
                        (self.fires, daily)).view(
                        np.recarray)

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


def data_range(data, min='median'):
    """
    return a data range for flexpart data

    optional keyword min = ['median', 'mean', 'min']
    """
    dmax = np.nanmax(data)
    if np.isnan(dmax):
        dmax = 1e5

    if min == 'mean':
        dmin = np.mean(data[data.nonzero()])
    elif min == 'median':
        dmin = np.median(data[data.nonzero()])
    else:
        dmin = np.nanmin(data[data.nonzero()])

    if np.isnan(dmin):
        dmin = 1e-5

    return [dmin, dmax]


def _normdatetime(Y, M, D, h, m, s):
    return datetime.datetime(
        Y, M, D) + datetime.timedelta(hours=h, minutes=m, seconds=s)


def _log_clevs(dat_min, dat_max):
    """
    ## create logorithmic color scale
    ## method uses strings to get the magnitude of variable
    #ss = "%1.2e" % (dat_max)
    #dmx = int('%s%s' % (ss[-3],int(ss[-2:])))
    #dmx=float(len(str(int(np.ceil(dat_max)))))
    #if data_range is None:
    #     dmn = int(np.round(np.log((1e-10*(dat_max-dat_min)))))
    #    ss = "%1.2e" % (1e-10*(dat_max-dat_min))
    #else:
    #    ss = "%1.2e" % (dat_min)
    #dmn = int('%s%s' % (ss[-3],int(ss[-2:])))
    """

    if dat_max > 0:
        dmx = int(np.round(np.log10(dat_max))) + 1
    else:
        # print 'dat_max not positive'
        # print dat_max
        dmx = 1

    if dat_min > 0:
        dmn = int(np.round(np.log10(dat_min)))
    elif dat_min == 0. or np.isnan(dat_min):
        # print 'dat_min not positive'
        # print dat_min
        dmn = dmx - 3

    # create equally spaced range
    if dmx == dmn:
        dmx = dmn + 1
    clevs = np.logspace(dmn, dmx, 100)

    return clevs


def add_nilu_logo(fig, xo=10, yo=520, origin='upper'):
    """ Adds NILU logo to background of plot """
    im = image.imread('/xnilu_wrk/flex_wrk/WEB_APS/imgs/logo_nilu.png')
    im[:, :, -1] = 0.5
    fig.figimage(im, xo, yo, origin=origin)
    return fig


def _cum_spec(inspectra, cum=True):
    """ helper function for the plot_spectra routines """

    if cum == 'norm':
        # Normalize the data so it fills to 100%
        # spectra = np.zeros(inspectra.shape)
        spectra = (
            inspectra.transpose() /
            np.sum(
                inspectra,
                axis=1)).transpose()
        spectra = np.cumsum(spectra[:, :], axis=1)
        # sums = np.sum(inspectra,axis=1)
        # for i,elem in enumerate(inspectra):
        #    spectra[i,:] = elem/sums[i]
    elif cum is True and cum != 'norm':
        spectra = np.cumsum(inspectra[:, :], axis=1)
    else:
        spectra = inspectra

    return spectra


def _getfile_lines(infile):
    """ returns all lines from a file or file string
    reverts to beginning of file."""

    if isinstance(infile, str):
        return file(infile, 'r').readlines()
    else:
        infile.seek(0)
        return infile.readlines()


def _gen_MapPar_fromHeader(H):
    """
    Define some default map parameters from the Header File.
    """

    MapPar = Structure()
    MapPar.llcrnrlat = H.outlat0
    MapPar.llcrnrlon = H.outlon0
    MapPar.urcrnrlat = H.outlat0 + (H.dyout * H.numygrid)
    MapPar.urcrnrlon = H.outlon0 + (H.dxout * H.numxgrid - 1)
    for k in MapPar.keys():
        print k, MapPar[k]

    return MapPar


def _gen_altitude_color(p, altlim=(0, 5000), numcon=1, cmap_id='jet'):
    """
    Generates color based on altlim = (min,max)
    """
    norm = mpl.colors.Normalize(vmin=altlim[0], vmax=altlim[1])
    cmap = plt.cm.get_cmap(cmap_id)
    # p = p/altmax
    cm = cmap(norm(p))
    # print p, cm
    return cm


def _gen_daylabels(P, H=None, dt=86400):
    """
    Uses H.loutstep to calculate 'days back' for plotting on clusters and trajectories.
    """
    if isinstance(P, int):

        if H:
            dt = abs(H.loutstep)
        return str(1 + int(abs(P)) / dt)
    else:
        return [str(1 + int(abs(p)) / dt) for p in P]


def _datarange(H, G, index=None):
    """ returns max and min for footprint and total column
        for a given range of releases [default is all] """
    Heightnn = H['Heightnn']
    fpmax = -999.
    if index is None:
        seek = range(G.shape[-1])
    else:
        seek = [index]

    for i in seek:
        zpmax = np.max(G[:, :, 0, i] / Heightnn[:, :, 0])
        if zpmax > fpmax:
            fpmax = zpmax
        # print fpmax
    tcmax = np.max(G)
    return ((0, fpmax), (0, tfcmax))


def _genEllipse(data, m, sizescale=20000,
                altlim=None):
    """

    Generates ellipses based on input array 'data'. Requires basemap instance, 'm'. NOTE::

        data = [itime,x,y,z,[size]]

        r,c = data.shape

        if c == 5:
            size function of data[:,4]
        else:
            size = 1*sizescale

    uses functions:

        * _gen_altitude_color
        * _gen_daylabels

    for labeling/color of ellipses

    """
    r, c = data.shape
    if altlim is None:
        altlim = (np.min(data[:, 3]), np.max(data[:, 3]))

    if c == 5:

        ell = [(Ellipse(xy=np.array(m(p[1], p[2])),
                        width=p[4] * sizescale, height=p[4] * sizescale,
                        angle=360,
                        facecolor=_gen_altitude_color(
                            p[3],
                            altlim=altlim,
                            cmap_id='gray'),
                        label=_gen_daylabels(p[0])),
                np.array(m(p[1], p[2]))) for p in data]
        # np.array( m(p[1],p[2])) ) for p in data]
    else:
        ell = [(Ellipse(xy=np.array(m(p[1], p[2])),
                        width=1e4 * sizescale, height=1e4 * sizescale,
                        angle=360,
                        facecolor=_gen_altitude_color(
                            p[3],
                            altlim=altlim,
                            cmap_id='gray'),
                        label=_gen_daylabels(p[0])),
                np.array(m(p[1], p[2]))) for p in data]

    return ell


def _shout(string, verbose=True):
    """ write a string to stdout if 'verbose' flag given """
    w = sys.stdout.write
    if string[-1] != '\n':
        string = string + '\n'
    if verbose:
        w(string)


def _gen_flexpart_colormap(ctbfile=None, colors=None):
    """
    Generate the ast colormap for FLEXPART
    """
    from matplotlib.colors import ListedColormap
    if ctbfile:
        try:
            colors = np.loadtxt(ctbfile)
        except:
            print "WARNING: cannot load ctbfile. using colors"
    if colors:
        name = 'user_colormap'
    if not colors:
        # AST Colorset for FLEXPART
        colors = [
            1.0000000e+00, 1.0000000e+00, 1.0000000e+00,
            9.9607843e-01, 9.1372549e-01, 1.0000000e+00,
            9.8431373e-01, 8.2352941e-01, 1.0000000e+00,
            9.6470588e-01, 7.1764706e-01, 1.0000000e+00,
            9.3333333e-01, 6.0000000e-01, 1.0000000e+00,
            8.9019608e-01, 4.4705882e-01, 1.0000000e+00,
            8.3137255e-01, 2.0000000e-01, 1.0000000e+00,
            7.5686275e-01, 0.0000000e+00, 1.0000000e+00,
            6.6274510e-01, 0.0000000e+00, 1.0000000e+00,
            5.4901961e-01, 0.0000000e+00, 1.0000000e+00,
            4.0784314e-01, 0.0000000e+00, 1.0000000e+00,
            2.4705882e-01, 0.0000000e+00, 1.0000000e+00,
            7.4509804e-02, 0.0000000e+00, 1.0000000e+00,
            0.0000000e+00, 2.8235294e-01, 1.0000000e+00,
            0.0000000e+00, 4.8627451e-01, 1.0000000e+00,
            0.0000000e+00, 6.3137255e-01, 1.0000000e+00,
            0.0000000e+00, 7.4509804e-01, 1.0000000e+00,
            0.0000000e+00, 8.4705882e-01, 1.0000000e+00,
            0.0000000e+00, 9.3725490e-01, 1.0000000e+00,
            0.0000000e+00, 1.0000000e+00, 9.7647059e-01,
            0.0000000e+00, 1.0000000e+00, 8.9411765e-01,
            0.0000000e+00, 1.0000000e+00, 8.0000000e-01,
            0.0000000e+00, 1.0000000e+00, 6.9019608e-01,
            0.0000000e+00, 1.0000000e+00, 5.6470588e-01,
            0.0000000e+00, 1.0000000e+00, 4.0000000e-01,
            0.0000000e+00, 1.0000000e+00, 0.0000000e+00,
            3.9607843e-01, 1.0000000e+00, 0.0000000e+00,
            5.6470588e-01, 1.0000000e+00, 0.0000000e+00,
            6.9019608e-01, 1.0000000e+00, 0.0000000e+00,
            7.9607843e-01, 1.0000000e+00, 0.0000000e+00,
            8.9411765e-01, 1.0000000e+00, 0.0000000e+00,
            9.7647059e-01, 1.0000000e+00, 0.0000000e+00,
            1.0000000e+00, 9.4509804e-01, 0.0000000e+00,
            1.0000000e+00, 8.7450980e-01, 0.0000000e+00,
            1.0000000e+00, 7.9215686e-01, 0.0000000e+00,
            1.0000000e+00, 7.0588235e-01, 0.0000000e+00,
            1.0000000e+00, 6.0392157e-01, 0.0000000e+00,
            1.0000000e+00, 4.8235294e-01, 0.0000000e+00,
            1.0000000e+00, 3.1372549e-01, 0.0000000e+00,
            1.0000000e+00, 0.0000000e+00, 1.4901961e-01,
            1.0000000e+00, 0.0000000e+00, 3.3333333e-01,
            1.0000000e+00, 0.0000000e+00, 4.4705882e-01,
            1.0000000e+00, 0.0000000e+00, 5.3725490e-01,
            1.0000000e+00, 0.0000000e+00, 6.1176471e-01,
            9.7647059e-01, 0.0000000e+00, 6.6666667e-01,
            8.9411765e-01, 0.0000000e+00, 6.6666667e-01,
            7.9607843e-01, 0.0000000e+00, 6.3921569e-01,
            6.9019608e-01, 0.0000000e+00, 5.9215686e-01,
            5.6470588e-01, 0.0000000e+00, 5.0980392e-01,
            3.9607843e-01, 0.0000000e+00, 3.8039216e-01]
        colors = np.reshape(colors, (-1, 3))
        name = 'flexpart_cmap'
    cmap = ListedColormap(colors, name)
    return cmap
