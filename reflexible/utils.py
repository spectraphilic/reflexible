# General utilities live here

from datetime import datetime

from matplotlib.dates import date2num
import numpy as np
import xarray as xr


class Structure(dict):
    """Basically a dictionary whose keys are attributes too."""

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def set_with_dict(self, d):
        """ set attributes with a dict """
        for k in d.keys():
            self.__setattr__(k, d[k])


class CacheDict(dict):
    """A dictionary that prevents itself from growing too much."""

    def __init__(self, maxentries):
        self.maxentries = maxentries
        super(CacheDict, self).__init__(self)

    def __setitem__(self, key, value):
        # Protection against growing the cache too much
        if len(self) > self.maxentries:
            # Remove a 10% of (arbitrary) elements from the cache
            entries_to_remove = self.maxentries // 10
            entries_to_remove += 1  # to avoid removing no entries
            for k in self.keys()[:entries_to_remove]:
                super(CacheDict, self).__delitem__(k)
        super(CacheDict, self).__setitem__(key, value)


def closest(num, numlist):
    """ returns the index of the *closest* value in a list """
    # check if we're using datetimes
    if isinstance(num, datetime):
        num = date2num(num)
        assert isinstance(numlist[0], datetime), \
               "num is date, numlist must be a list of dates"
        numlist = date2num(numlist)

    return (np.abs(numlist - num)).argmin()


def data_range(data, min='median'):
    """
    return a data range for flexpart data

    optional keyword min = ['median', 'mean', 'min']
    """
    if isinstance(data, xr.Variable):
        data = data.values

    dmax = np.nanmax(data)
    if np.isnan(dmax):
        dmax = 1e5

    data = data[data.nonzero()]
    if min == 'mean':
        dmin = np.mean(data)
    elif min == 'median':
        dmin = np.median(data)
    else:
        dmin = np.nanmin(data)

    if np.isnan(dmin):
        dmin = 1e-5

    return [dmin, dmax]


#
# Curtain related functions
#

def __groundlevel_for_line(H, X, Y, coords, index=0):
    """
    extracts ground level from H.heightnn along a track of lon, lat

    input:  H or H.heightnn (takes the lowest level)

            X, Y ,= H.longitude, H.latitude

            coords = zip(x, y)

    output: groundlevel (a 1-d array with shape (len(flighttrack)

    """
    try:
        hgt = H.Heightnn[:, :, 0]
    except Exception:
        hgt = H

    # fix for hgt offset
    hgt = hgt - hgt.min()

    grndlvl = np.zeros(len(coords))
    for i, (x, y) in enumerate(coords):
        I = closest(x, X)
        J = closest(y, Y)
        grndlvl[i] = hgt[I, J]
        grndlvl = np.nan_to_num(grndlvl)

    return grndlvl - grndlvl.min()


def curtain_agltoasl(H, curtain_agl, coords, below_gl=0.0):
    """ converts the agl curtain to asl

        adds .asl_axis attribute to H
    """

    gl = __groundlevel_for_line(H, H.longitude, H.latitude, coords)
    H.asl_axis = np.linspace(0, H.outheight[-1])
    xp = H.outheight - H.outheight[0]
    casl = np.zeros((len(H.asl_axis), len(coords)))

    for i in range(len(coords)):
        casl[:, i] = np.interp(H.asl_axis, xp + gl[i], \
                              curtain_agl[:, i], left=below_gl)

    return casl
