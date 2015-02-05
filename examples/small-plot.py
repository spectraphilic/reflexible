#!/usr/bin/env python

"""
This is an example providing a very brief overview of some of the
pflexible functionality.

Test Data:
    Note, FLEXPART output is large. I have made available a complete
    run for one day from the POLARCAT NOAA-ICEALOT Cruise here:
    http://niflheim.nilu.no/~burkhart/sharing/pflexpart_testdata.tgz
"""

#builtins
import os

import numpy as np
import mapping as mp
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import shiftgrid, addcyclic

#import pflexible as pf
import reflexible as rf
import reflexible.conv2netcdf4 as conv


def plot_backward(SOURCE_FILE, OUTPUT_DIR):
    #read the header of a FLEXPART output
    # H = conv.Header(SOURCE_FILE)
    H = rf.Header(SOURCE_FILE)

    # Below creates a D object as an attribute of H, nothing is returned
    # H.C is keyed by (s,k) where s = speciesID, and k is rel_i
    # NOTE: speciesID is 0-indexed
    # H.fill_backward(nspec=(0,1)) #get the first species, all releases

    # # Get the trajectories, these are overlayed when we plot the
    # # sensitivities below.
    # T = pf.read_trajectories(H)

    # Plot everything up
    TC = None #TC Empty figure object, the idea is to reuse this
    FP = None #FP Empty figure object

    # iterate over every species and every timestep (k)
    for s,k in H.C:
        data = H.C[(s,k)]
        # total column
        TC = plot_totalcolumn(H, data, datainfo_str='NORTHSEA', FIGURE=TC)
        # TC = plot_trajectory(H, T, k, FIGURE=TC)
        filename = '%s_tc_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        TC.fig.savefig(ofilename)

        # footprint
        FP = plot_at_level(H, data, datainfo_str='NORTHSEA', FIGURE=FP)
        # FP = plot_trajectory(H, T, k, FIGURE=FP)
        filename = '%s_fp_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        FP.fig.savefig(ofilename)


# Plotting funtions follows

def plot_at_level(H, C, level=1, \
                   ID=' ', \
                   map_region=5, projection='lcc', \
                   overlay=False,
                   datainfo_str=None, log=True,
                   data_range=None, coords=None, FIGURE=None,
                   plot_title=None,
                   units=None,
                   **kwargs
                   ):
    """
    TODO:
    -make units a function of H['species']
    """
    if units is None:
        units = H.output_unit
    rel_i = C.rel_i
    species = C.species
    timestamp = C.timestamp
    data = C.slabs[level]

    if level == 1:
        level_desc = 'Footprint'
    else:
        level_desc = H.outheight[level - 1]

    dmax = data.max()
    dmin = data.min()
    # print dmax,dmin

    if data_range == None:
        data_range = [dmin, dmax]

    if H.direction == 'backward' and H.options['readp']:
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            # need to find a way to set m.a.s.l. or hPa here
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f, Z2: %.2f (%s)\n""" \
                           % (dmax, units, zp1, zp2, H.alt_unit)
        if plot_title is None:
            plot_title = """
        %s Sensitivity at %s %s: %s\n
        Release Start: %s, Release End: %s""" % \
            (ID, level_desc, H['alt_unit'], species, H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s """ % (dmax, units)
        if plot_title is None:
            plot_title = """ %s Sensitivity at %s %s: %s \n %s """ % (ID, level_desc, H['alt_unit'], species, timestamp)

    FIGURE = plot_sensitivity(H, data, \
                           data_range=data_range, \
                           rel_i=rel_i, log=log,
                           map_region=map_region, projection=projection,
                           units=units, datainfo_str=datainfo_str,
                           overlay=overlay,
                           coords=coords, FIGURE=FIGURE, **kwargs)


    FIGURE.ax.set_title(plot_title, fontsize=10)
    return FIGURE


def plot_totalcolumn(H, C=None, \
                    ID=' ', \
                    map_region=5, projection='lcc', \
                    data_range=None, coords=None,
                    FIGURE=None, overlay=False,
                    datainfo_str=None, **kwargs):

    if C is None:
        C = H.C[(0, 0)]

    if 'units' in kwargs:
        units = kwargs.pop('units')
    else:
        units = 'ns m kg-1'

    rel_i = C.rel_i
    species = C.species
    timestamp = C.timestamp
    data = C.slabs[0]


    if data_range == None:
        dmax = data.max()
        dmin = data.min()
        data_range = [dmin, dmax]
    else:
        dmin, dmax , = data_range
        # print dmin,dmax


    if H.direction == 'backward':
        rel_i = C.rel_i
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f, Z2: %.2f (%s)\n""" % \
                               (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n
        Release Start: %s, Release End: %s""" % \
                               (ID, species, H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s""" % (dmax, units)
        plot_title = """
        %s Total Column Sensitivity: %s\n %s """ % (ID, species, timestamp)


    FIGURE = plot_sensitivity(H, data, \
                           data_range=data_range, \
                           rel_i=rel_i, map_region=map_region, \
                           projection=projection, units=units, \
                           datainfo_str=datainfo_str, coords=coords, \
                           FIGURE=FIGURE, overlay=overlay, **kwargs)

    FIGURE.ax.set_title(plot_title, fontsize=10)

    return FIGURE


def plot_sensitivity(H, data, \
             data_range=None, \
             units='ns m^2 / kg', \
             datainfo_str=None, \
             plottitle=None, \
             rel_i=None,
             map_region=None, projection=None,
             dropm=None, coords=None,
             overlay=False,
             transform=True,
             log=True,
             FIGURE=None,
             MapPar=None,
             FigPar=None,
             cax_title=None,
             autofit=False, method='contourf', lsmask=False):
    # autofit=False,method='imshow'):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::

        > FIG = plot_sensitivity(H,data,*kwargs)

    This returns the FIGURE object, and plots the sensitivity from the data contained in the "C"
    array.

    Inputs
      H = a :class:`Header` instance for a FLEXPART run.
      data = a 2d data array containing the sensitivity values to plot, this can be extracted from a
      grid instance (see :func:`readgridV8` and :func:`get_slabs`)

    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ================================================
      keyword               Description [default]
      =============         ================================================
      data_range            range of data for scale bars, if None will
                            be taken from min/max of data
      cax_title             string to be passed as colorbar title (units will
                            be passed to the format argument)
      units                 units for the scale bar
      datainfo_str          A string for labeling the scale bar.
      plottitle             Title for the plot.
      rel_i                 Release index to plot from the data array
      map_region                A map_region specified in mapping.py
      projection            [deprecated] use pre-defined map_regions.
      dropm                 Force creation of a new basemap instance
      coords                Used with autofit option. An array of lat,lon
                            values for the mapping module to autofit a
                            basemap instance to.
      autofit               Try to generate map_region automagically (flakey)
      overlay               Force removal of previous figure elements.
      transform             For use with imshow method, if your data is not
                            in same coordinates as projection, try to transform
                            the data to the basemap projection.
      log                   Create a logarithmic color scale.
      FIGURE                A FIGURE instance from mapping module get_FIGURE
      MapPar                A Structure of paramters to be passed to the
                            basemap class when creating an instance.
      method                The method to use for plotting array data. May be
                            one of: [pcolormesh], imshow, or contourf
      lsmask                set to True to draw a grey landseamask [False]
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.
        lsmask = True seems to only work with certain projections (POLARCAT)

    .. note::
        This is the primary main function for creating plots of flexpart output.
        Most the other routines are simply wrappers to this function, passing
        arguments in with some predefined settings. For more information on the
        mechanics of this function, see the mapping.py module and the matplotlib
        basemap toolkit.


    """



    methods = ['imshow', 'pcolormesh', 'contourf', 'contour', 'None']
    assert method in methods, "method keyword must be one of: %s" % methods

    if FIGURE is None:
        if map_region is None:
            if MapPar is None:
                if autofit:
                    MapPar = _gen_MapPar_fromHeader(H)

        FIGURE = mp.get_FIGURE(map_region=map_region,
                               projection=projection, coords=coords,
                               MapPar=MapPar, FigPar=FigPar)
    else:
        if FIGURE.m is None:
            FIGURE = mp.get_FIGURE(fig=FIGURE.fig, ax=FIGURE.ax, map_region=map_region,
                               projection=projection, coords=coords,
                               MapPar=MapPar, FigPar=FigPar)

    if overlay is False:
        del FIGURE.ax.images[FIGURE.indices.images:]
        del FIGURE.ax.collections[FIGURE.indices.collections:]
        del FIGURE.ax.lines[FIGURE.indices.lines:]

    if dropm != None:
        try:
            del m
            plt.close('all')
        except:
            print 'could not drop m'

    # # make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax

    # # make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    # # set up transformations for the data array
    if method == 'imshow':
        if m.projection not in ['cyl', 'merc', 'mill']:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)[:-1]
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)[:-1]
            data = data[:-1, :-1]
        else:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)

        # # transform to nx x ny regularly spaced native projection grid
        if transform:
            dx = 2.*np.pi * m.rmajor / len(lons)
            nx = int((m.xmax - m.xmin) / dx) + 1; ny = int((m.ymax - m.ymin) / dx) + 1
            if nx is 1:
                topodat = data
            else:
                topodat = m.transform_scalar(data, lons, lats, nx, ny)
        else:
            topodat = data

    if method != 'imshow':
        # # Check to see if a cyclic wraparound is required

        lons = np.arange(H.outlon0, H.outlon0 + (H.dxout * H.numxgrid), H.dxout)
        lats = np.arange(H.outlat0, H.outlat0 + (H.dyout * H.numygrid), H.dyout)
        # # if required add polar coordinates
        if 'npstere' in m.projection:
            if lats[-1] != 90.:
                npole = np.ones(len(lons)).T * data[0, :]
                npole = np.reshape(npole, (1, -1))
                data = np.vstack((npole, data))
                lats = np.hstack((lats, [lats[-1] + H.dyout]))

        if m.projection != 'merc':
            if lons[-1] - lons[0] < 360.:
                print "data.shape:", data.shape, len(lons)
                topodat, lons = addcyclic(data, lons)

        if m.projection == 'merc':
            topodat = data

    # # get min/max range
    if data_range != None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = data_range(data)


    if log:
        clevs = _log_clevs(dat_min, dat_max)

    else:
        clevs = [i for i in np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    # # draw land sea mask
    # m.fillcontinents(zorder=0)
    if lsmask:
        m.drawlsmask(ocean_color='grey', zorder=-10)

    # # Plot Release Location if points were read
    if H.options['readp']:
        if rel_i:
            releaselocation = (H.xpoint[rel_i], H.ypoint[rel_i])
            xpt, ypt = m(releaselocation[0], releaselocation[1])
            # # Remove prior location point
            try:
                del ax.lines[-1]
            except:
                pass
            location, = m.plot([xpt], [ypt], 'bx', linewidth=6, markersize=20, zorder=1000)

    # # Plot the footprint

    # # Set up the IMAGE
    # # cmapnames = ['jet', 'hsv', 'gist_ncar', 'gist_rainbow', 'cool', 'spectral']
    # colmap = plt.get_cmap('jet')
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    # # Plotting METHODS (pcolormesh now default, imshow is smoother)
    # print topodat.max(), topodat.min(), topodat.shape
    if method == 'imshow':
        im = m.imshow(topodat, cmap=colmap, zorder=-1,
                      norm=mpl.colors.LogNorm(vmin=clevs[0],
                                              vmax=clevs[-1]))

    if method == 'pcolormesh':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.pcolormesh(nx, ny, topodat, cmap=colmap,
                          norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                  vmax=clevs[-1]))
    if method == 'contourf':
        # # Trying some fancier scaling
        # cnts,bins = np.histogram(topodat,bins=100)
        # topodat = np.ma.masked_where(topodat< .05* np.average((0,bins[1])),topodat)
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    if method == 'contour':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contour(nx, ny, topodat, cmap=colmap,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    # # Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    # # CREATE COLORBAR
    # # Note, with upgrades to matplotlib and basemap had to make some
    # # changes here... no more 'ghost' axes
    # # does a colorbar already exist?
    try:
        cb = FIGURE.cb
        cax = FIGURE.cax
        cb.update_normal(im)
    except:
    # # make a copy of the image object, change
    # # colormap to linear version of the precip colormap.
    # pdb.set_trace()
    # im2 = copy.copy(im)
    # im2.set_cmap(colmap)
    # # create new axis for colorbar.
        h = 0.5 * h
        l = l + w + .03
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
    # # using im2, not im (hack to prevent colors from being
    # # too compressed at the low end on the colorbar - results
    # # from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax)  # , format='%3.2g') # draw colorbar
        FIGURE.cax = cax
        FIGURE.cb = cb
    # cb.update_normal(im2)




    # # set colorbar label and ticks
    # pdb.set_trace()
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10])  # #clevs, by 10 steps
    clabels.append(clevs[-1])  # # add the last label
    # cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])
                        # fontproperties=p_cax)
    if H.direction == 'forward':
        cax.set_title('%s' % units,
                      fontproperties=p_cax)
    else:
        if cax_title:
            cax.set_title(cax_title.format(units), fontproperties=p_cax)
        else:
            cax.set_title('sensitivity\n({0})'.format(units),
                          fontproperties=p_cax)

    # # make the original axes current again
    plt.axes(ax)

    # # write text information block on plot
    # # first try to remove prior text by accessing
    # # the last text element added to the axes from
    # # the prior iteration.
    # # This is tricky when using together with plot_clusters...
    # # need to figure out how to resolve the indexing
    # # of what texts, collections, etc to delete, when iterating.
    try:
        del ax.texts[FIGURE.indices.texts:]
        del ax.artists[FIGURE.indices.artists:]
    except:
        pass
    if datainfo_str:
        plt.text(l, b + 1000,
             datainfo_str,
             fontsize=10,
             bbox=dict(boxstyle="round",
                     ec=(1., 0.5, 0.5),
                     fc=(1., 0.8, 0.8),
                     alpha=0.8
                     )
             )

    FIGURE.ax = ax
    FIGURE.m = m
    FIGURE.fig = fig

    if plottitle != None:
        # plt.title(plottitle,fontproperties=p_cax)
    # plt = plt
        FIGURE.ax.set_title(plottitle, fontsize=10)
    return FIGURE


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
        # # AST Colorset for FLEXPART
        colors = [1.0000000e+00, 1.0000000e+00, 1.0000000e+00
                , 9.9607843e-01, 9.1372549e-01, 1.0000000e+00
                , 9.8431373e-01, 8.2352941e-01, 1.0000000e+00
                , 9.6470588e-01, 7.1764706e-01, 1.0000000e+00
                , 9.3333333e-01, 6.0000000e-01, 1.0000000e+00
                , 8.9019608e-01, 4.4705882e-01, 1.0000000e+00
                , 8.3137255e-01, 2.0000000e-01, 1.0000000e+00
                , 7.5686275e-01, 0.0000000e+00, 1.0000000e+00
                , 6.6274510e-01, 0.0000000e+00, 1.0000000e+00
                , 5.4901961e-01, 0.0000000e+00, 1.0000000e+00
                , 4.0784314e-01, 0.0000000e+00, 1.0000000e+00
                , 2.4705882e-01, 0.0000000e+00, 1.0000000e+00
                , 7.4509804e-02, 0.0000000e+00, 1.0000000e+00
                , 0.0000000e+00, 2.8235294e-01, 1.0000000e+00
                , 0.0000000e+00, 4.8627451e-01, 1.0000000e+00
                , 0.0000000e+00, 6.3137255e-01, 1.0000000e+00
                , 0.0000000e+00, 7.4509804e-01, 1.0000000e+00
                , 0.0000000e+00, 8.4705882e-01, 1.0000000e+00
                , 0.0000000e+00, 9.3725490e-01, 1.0000000e+00
                , 0.0000000e+00, 1.0000000e+00, 9.7647059e-01
                , 0.0000000e+00, 1.0000000e+00, 8.9411765e-01
                , 0.0000000e+00, 1.0000000e+00, 8.0000000e-01
                , 0.0000000e+00, 1.0000000e+00, 6.9019608e-01
                , 0.0000000e+00, 1.0000000e+00, 5.6470588e-01
                , 0.0000000e+00, 1.0000000e+00, 4.0000000e-01
                , 0.0000000e+00, 1.0000000e+00, 0.0000000e+00
                , 3.9607843e-01, 1.0000000e+00, 0.0000000e+00
                , 5.6470588e-01, 1.0000000e+00, 0.0000000e+00
                , 6.9019608e-01, 1.0000000e+00, 0.0000000e+00
                , 7.9607843e-01, 1.0000000e+00, 0.0000000e+00
                , 8.9411765e-01, 1.0000000e+00, 0.0000000e+00
                , 9.7647059e-01, 1.0000000e+00, 0.0000000e+00
                , 1.0000000e+00, 9.4509804e-01, 0.0000000e+00
                , 1.0000000e+00, 8.7450980e-01, 0.0000000e+00
                , 1.0000000e+00, 7.9215686e-01, 0.0000000e+00
                , 1.0000000e+00, 7.0588235e-01, 0.0000000e+00
                , 1.0000000e+00, 6.0392157e-01, 0.0000000e+00
                , 1.0000000e+00, 4.8235294e-01, 0.0000000e+00
                , 1.0000000e+00, 3.1372549e-01, 0.0000000e+00
                , 1.0000000e+00, 0.0000000e+00, 1.4901961e-01
                , 1.0000000e+00, 0.0000000e+00, 3.3333333e-01
                , 1.0000000e+00, 0.0000000e+00, 4.4705882e-01
                , 1.0000000e+00, 0.0000000e+00, 5.3725490e-01
                , 1.0000000e+00, 0.0000000e+00, 6.1176471e-01
                , 9.7647059e-01, 0.0000000e+00, 6.6666667e-01
                , 8.9411765e-01, 0.0000000e+00, 6.6666667e-01
                , 7.9607843e-01, 0.0000000e+00, 6.3921569e-01
                , 6.9019608e-01, 0.0000000e+00, 5.9215686e-01
                , 5.6470588e-01, 0.0000000e+00, 5.0980392e-01
                , 3.9607843e-01, 0.0000000e+00, 3.8039216e-01]
        colors = np.reshape(colors, (-1, 3))
        name = 'flexpart_cmap'
    cmap = ListedColormap(colors, name)
    return cmap


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


    # # create equally spaced range
    if dmx == dmn:
        dmx = dmn + 1
    clevs = np.logspace(dmn, dmx, 100)

    return clevs


if __name__ == "__main__":
    """ run through plotting routines """

    #define paths
    SOURCE_FILE="/tmp/grid_time_20100601210000.nc"
    OUTPUT_DIR='/tmp/plots2'

    plot_backward(SOURCE_FILE, OUTPUT_DIR)
