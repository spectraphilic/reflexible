from __future__ import division
from builtins import range
from past.utils import old_div
########### HELPER FUNCTIONS ##########

import datetime

import numpy as np
import matplotlib.image as image

# Matplotlib
from matplotlib.dates import date2num



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
        seek = list(range(G.shape[-1]))
    else:
        seek = [index]

    for i in seek:
        zpmax = np.max(old_div(G[:, :, 0, i], Heightnn[:, :, 0]))
        if zpmax > fpmax:
            fpmax = zpmax
        # print fpmax
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

