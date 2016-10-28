########### HELPER FUNCTIONS ##########

from __future__ import print_function

import os
import datetime

import numpy as np
import matplotlib.image as image

# Matplotlib
from matplotlib.dates import date2num


def get_fpdirs(pathnames):
    """Return the <options> and <output> dirs from a `pathnames` file."""
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
    return options_dir, output_dir


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


def add_nilu_logo(fig, xo=10, yo=520, origin='upper'):
    """ Adds NILU logo to background of plot """
    im = image.imread('/xnilu_wrk/flex_wrk/WEB_APS/imgs/logo_nilu.png')
    im[:, :, -1] = 0.5
    fig.figimage(im, xo, yo, origin=origin)
    return fig


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
        # print(fpmax)
    tcmax = np.max(G)
    return ((0, fpmax), (0, tcmax))


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


def _shout(string, verbose=True):
    """ write a string to stdout if 'verbose' flag given """
    import sys
    w = sys.stdout.write
    if string[-1] != '\n':
        string += '\n'
    if verbose:
        w(string)

