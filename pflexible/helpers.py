########### HELPER FUNCTIONS ##########

import datetime

import numpy as np
import matplotlib.image as image
from matplotlib.patches import Ellipse

# Matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

import pflexible as pf


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


def _gen_MapPar_fromHeader(H):
    """
    Define some default map parameters from the Header File.
    """

    MapPar = pf.Structure()
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
    return ((0, fpmax), (0, tcmax))


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
