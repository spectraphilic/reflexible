#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib as mpl
from matplotlib import font_manager, colors
import matplotlib.pyplot as plt
from mpl_toolkits import basemap

# local imports
import reflexible.mapping as mp
from reflexible import Structure, CacheDict

# Dictionary for caching the figures
_figure_cache = CacheDict(10)



def plot_at_level(H, data, level=1,
                  ID=' ', rel_i=None, species=None,
                  timestamp=None,
                  map_region=5,
                  overlay=False,
                  datainfo_str=None, log=True,
                  data_range=None,
                  plot_title=None,
                  units=None,
                  **kwargs):
    """
    TODO: make units a function of H['species']
    """
    if units is None:
        units = H.output_unit

    if level == 1:
        level_desc = 'Footprint'
    else:
        level_desc = H.outheight[level - 1]

    if data_range is None:
        dmax = data.max()
        dmin = data.min()
        data_range = [dmin, dmax]
    else:
        dmin, dmax, = data_range

    if H.direction == 'backward' and H.options['readp']:
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            # need to find a way to set m.a.s.l. or hPa here
            datainfo_str = (
                " Max Value: %.2g %s\n Release Z1: %.2f, Z2: %.2f (%s)\n" % (
                    dmax, units, zp1, zp2, H.alt_unit))
        if plot_title is None:
            plot_title = """
        %s Sensitivity at %s %s: %s\n
        Release Start: %s, Release End: %s""" % (
                ID, level_desc, H['alt_unit'], species,
                H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s """ % (dmax, units)
        if plot_title is None:
            plot_title = """ %s Sensitivity at %s %s: %s \n %s """ % (
                ID, level_desc, H['alt_unit'], species, timestamp)

    figure = plot_sensitivity(H, data,
                              data_range=data_range,
                              rel_i=rel_i, log=log,
                              map_region=map_region,
                              units=units, datainfo_str=datainfo_str,
                              overlay=overlay,
                              **kwargs)

    figure.ax.set_title(plot_title, fontsize=10)
    return figure


def plot_totalcolumn(H, data=None,
                     ID=' ', rel_i=None, species=None,
                     timestamp=None,
                     map_region=5,
                     data_range=None,
                     overlay=False,
                     datainfo_str=None, **kwargs):
    if 'units' in kwargs:
        units = kwargs.pop('units')
    else:
        units = 'ns m kg-1'

    if data_range is None:
        dmax = data.max()
        dmin = data.min()
        data_range = [dmin, dmax]
    else:
        dmin, dmax, = data_range

    if H.direction == 'backward':
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            datainfo_str = (
                " Max Value: %.2g %s\n Release Z1: %.2f, Z2: %.2f (%s)\n" % (
                    dmax, units, zp1, zp2, H.alt_unit))
        plot_title = """
        %s Total Column Sensitivity: %s\n
        Release Start: %s, Release End: %s""" % (
            ID, species, H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s""" % (dmax, units)
        plot_title = """
        %s Total Column Sensitivity: %s\n %s """ % (ID, species, timestamp)

    figure = plot_sensitivity(H, data,
                              data_range=data_range,
                              rel_i=rel_i, map_region=map_region,
                              units=units,
                              datainfo_str=datainfo_str,
                              overlay=overlay, **kwargs)

    figure.ax.set_title(plot_title, fontsize=10)

    return figure


def _get_figure(fig=None, ax=None, m=None, map_region=None,
               map_par=None, fig_par=None, image=None):
    """Returns a matplotlib figure based on the parameters.

    The idea is that I create a :class:`Structure` that contains the figure,
    ax, and m instance. I also add a field for "indices".
    This is so as to be able to reuse figures. This saves a huge amount of
    time, as creating then basemap instance can be time consuming.

    USAGE::
        >>> fig = _get_figure()
        or
        >>> fig = _get_figure(map_region='polarcat')

    Returns
       This will return the "FIG" object, which has attributes: `fig`, `ax`,
       `m`, and `indices`. The indices are used for deleting lines, texts,
       collections, etc. if and when we are reusing the figure instance. The
       indices basically give us a reference to the *empty* map, so we can
       delete lines without losing meridians or parallels for example.

         ============      ======================================
         keys              description
         ============      ======================================
         fig               A pyplot.fig instance, use
                           plt.figure(FIG.fig.number) to make the
                           fig active (for example to use
                           plt.savefig('filename.png')
         m                 The basemap instance so you can do:
                           x,y = FIG.m(lon,lat)
         ax                The axes
         indices           with index for texts, images, lines,
                           and collections
         ============      ======================================

    """
    figure = Structure()

    if m is None:
        if image:
            fig, m = mp.get_base_image(image, map_region=map_region,
                                       map_par=map_par, fig_par=fig_par)
        else:
            fig, m = mp.get_base1(map_region=map_region, map_par=map_par,
                                  fig_par=fig_par, fig=fig)
            figure.fig = fig
            figure.m = m
            figure.ax = fig.gca()
    else:
        figure.m = m

    if fig is None:
        figure.fig = plt.figure()
        fig = figure.fig
    else:
        figure.fig = fig

    figure.ax = ax if ax is not None else fig.gca()

    figure.indices = Structure()
    figure.indices.texts = len(figure.ax.texts)
    figure.indices.images = len(figure.ax.images)
    figure.indices.collections = len(figure.ax.collections)
    figure.indices.lines = len(figure.ax.lines)

    print("Using figure: %s" % figure.fig.number)
    return figure


def plot_sensitivity(H, data,
                     data_range=None,
                     units='ns m^2 / kg',
                     datainfo_str=None,
                     plottitle=None,
                     rel_i=None,
                     map_region=None,
                     overlay=False,
                     transform=True,
                     log=True,
                     map_par=None,
                     fig_par=None,
                     cax_title=None,
                     method='contourf',
                     lsmask=False):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::
        FIG = plot_sensitivity(H,data,*kwargs)

    This returns the figure object, and plots the sensitivity from the data
    contained in the "D" array.

    Inputs
        H : a :class:`Header` instance for a FLEXPART run.
        data : a 2d data array containing the sensitivity values to plot,
            this can be extracted from a grid instance
            (see :func:`readgridV8` and :func:`get_slabs`)

    Returns
        A "mapping.py" ``figure`` object.

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
      map_region            A map_region specified in mapping.py
      overlay               Force removal of previous figure elements.
      transform             For use with imshow method, if your data is not
                            in same coordinates as projection, try to transform
                            the data to the basemap projection.
      log                   Create a logarithmic color scale.
      map_par               A Structure of parameters to be passed to the
                            basemap class when creating an instance.
      fig_par               A Structure of parameters to be passed to the
                            figure when creating an instance.
      method                The method to use for plotting array data. May be
                            one of: [pcolormesh], imshow, contour or contourf.
      lsmask                set to True to draw a grey landseamask [False]
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.
        lsmask = True seems to only work with certain projections (polarcat).

    .. note::
        This is the primary main function for creating plots of flexpart
        output.  Most the other routines are simply wrappers to this function,
        passing arguments in with some predefined settings. For more
        information on the mechanics of this function, see the mapping.py
        module and the matplotlib basemap toolkit.


    """
    data = data.T

    # Look if the figure is in cache already, and if not, cache it
    figure_key = map_region + str(map_par) + str(fig_par)
    try:
        figure = _figure_cache[figure_key]
    except KeyError:
        figure = _figure_cache[figure_key] = _get_figure(
            map_region=map_region, map_par=map_par, fig_par=fig_par)

    if overlay is False:
        del figure.ax.images[figure.indices.images:]
        del figure.ax.collections[figure.indices.collections:]
        del figure.ax.lines[figure.indices.lines:]

    # make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    fig = figure.fig
    m = figure.m
    ax = figure.ax

    # make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    # set up transformations for the data array
    if method == 'imshow':
        if m.projection not in ['cyl', 'merc', 'mill']:
            lats = np.arange(
                H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)[:-1]
            lons = np.arange(
                H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)[:-1]
            data = data[:-1, :-1]
        else:
            lats = np.arange(
                H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)
            lons = np.arange(
                H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)

        # transform to nx x ny regularly spaced native projection grid
        if transform:
            dx = 2. * np.pi * m.rmajor / len(lons)
            nx = int((m.xmax - m.xmin) / dx) + 1
            ny = int((m.ymax - m.ymin) / dx) + 1
            if nx is 1:
                topodat = data
            else:
                topodat = m.transform_scalar(data, lons, lats, nx, ny)
        else:
            topodat = data
    else:
        # Check to see if a cyclic wraparound is required
        lons = np.arange(
            H.outlon0, H.outlon0 + (H.dxout * H.numxgrid), H.dxout)
        lats = np.arange(
            H.outlat0, H.outlat0 + (H.dyout * H.numygrid), H.dyout)
        # if required add polar coordinates
        if 'npstere' in m.projection:
            if lats[-1] != 90.:
                npole = np.ones(len(lons)).T * data[0, :]
                npole = np.reshape(npole, (1, -1))
                data = np.vstack((npole, data))
                lats = np.hstack((lats, [lats[-1] + H.dyout]))

        if m.projection != 'merc':
            if lons[-1] - lons[0] < 360.:
                topodat, lons = basemap.addcyclic(data, lons)
            else:
                topodat = data

    # get min/max range
    if data_range is not None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = data_range(data)

    if log:
        clevs = _gen_log_clevs(dat_min, dat_max)
    else:
        clevs = [i for i in
                 np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    # draw land sea mask
    # m.fillcontinents(zorder=0)
    if lsmask:
        m.drawlsmask(ocean_color='grey', zorder=-10)

    # Plot Release Location if points were read
    if H.options['readp']:
        if rel_i:
            releaselocation = (H.xpoint[rel_i], H.ypoint[rel_i])
            xpt, ypt = m(releaselocation[0], releaselocation[1])
            # Remove prior location point
            try:
                del ax.lines[-1]
            except:
                pass
            location, = m.plot([xpt], [ypt], 'bx', linewidth=6, markersize=20,
                               zorder=1000)

    # custom flexpart colormap
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    # Plotting METHODS (pcolormesh now default, imshow is smoother)
    # print(topodat.max(), topodat.min(), topodat.shape)
    # ERROR: tododat might be used uninitialized
    if method == 'imshow':
        im = m.imshow(topodat, cmap=colmap, zorder=-1,
                      norm=mpl.colors.LogNorm(vmin=clevs[0],
                                              vmax=clevs[-1]))
    elif method == 'pcolormesh':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.pcolormesh(nx, ny, topodat, cmap=colmap,
                          norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                  vmax=clevs[-1]))
    elif method == 'contourf':
        # Trying some fancier scaling
        # cnts,bins = np.histogram(topodat,bins=100)
        # topodat = np.ma.masked_where(topodat< .05* np.average((0,bins[1])),
        #                              topodat)
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))
    elif method == 'contour':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contour(nx, ny, topodat, cmap=colmap,
                       norm=mpl.colors.LogNorm(vmin=clevs[0],
                                               vmax=clevs[-1]))
    else:
        raise ValueError("`method` param '%s' is not a valid one." % method)

    # Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    # CREATE COLORBAR
    # Note, with upgrades to matplotlib and basemap had to make some
    # changes here... no more 'ghost' axes
    # does a colorbar already exist?
    try:
        cb = figure.cb
        cax = figure.cax
        cb.update_normal(im)
    except:
        # make a copy of the image object, change
        # colormap to linear version of the precip colormap.
        # pdb.set_trace()
        # im2 = copy.copy(im)
        # im2.set_cmap(colmap)
        # create new axis for colorbar.
        h = 0.7 * h
        l = l + w  # + .01
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
        # using im2, not im (hack to prevent colors from being
        # too compressed at the low end on the colorbar - results
        # from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax)
        figure.cax = cax
        figure.cb = cb
    # cb.update_normal(im2)

    # set colorbar label and ticks
    # pdb.set_trace()
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10])  # #clevs, by 10 steps
    clabels.append(clevs[-1])  # add the last label
    # cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])

    if H.direction == 'forward':
        cax.set_title('%s' % units,
                      fontproperties=p_cax)
    else:
        if cax_title:
            cax.set_title(cax_title.format(units), fontproperties=p_cax)
        else:
            cax.set_title('sensitivity\n({0})'.format(units),
                          fontproperties=p_cax)

    # make the original axes current again
    plt.axes(ax)

    # write text information block on plot
    # first try to remove prior text by accessing
    # the last text element added to the axes from
    # the prior iteration.
    # This is tricky when using together with plot_clusters...
    # need to figure out how to resolve the indexing
    # of what texts, collections, etc to delete, when iterating.
    try:
        del ax.texts[figure.indices.texts:]
        del ax.artists[figure.indices.artists:]
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

    figure.ax = ax
    figure.m = m
    figure.fig = fig

    if plottitle is not None:
        figure.ax.set_title(plottitle, fontsize=10)

    return figure


def plot_trajectory(H, T, rel_i,
                    map_region=None,
                    overlay=True,
                    draw_circles=True,
                    draw_labels=True, days_back=20,
                    cbar2=True,
                    cbar2_title=None,
                    map_par=None,
                    fig_par=None):
    """Plot the center trajectory of the plume on the map

    Usage::

        > FIG = plot_trajectory(H,RT,rel_id,*kwargs)

    .. note::
        You can substitude "None" for the :class:`Header` instance if you
        haven't created one.


    Returns
      A :mod:`mapping` "figure" object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         =============================================
      keyword               Description [default]
      =============         =============================================
      rel_i                 **required** release index
      overlay               [True] will overlay the trajectory
                            on top of another map instance.
      draw_circles          [True] will mark the trajectory with
                            circles. If [False] nothing is drawn.
      draw_labels           [True] will draw the 'day' numbers on the
                            trajectories.
      days_back             For how many days back should the labels be
                            shown? [20]
      cbar2                 [True] draws the scale bar as a second axis.
      cbar2_title           Optional argument to overide the cbar title.
      map_par               A Structure of mapping parameters to pass
                            to the basemap instance if desired.
      fig_par               A Structure of parameters to be passed to the
                            figure when creating an instance.
      =============         =============================================

    .. todo::
        Who knows?! Could probably handle the overlaying better.

    .. note::
        Just set H to "None" if you haven't already created a "Header"
        instance.

    """

    if H:
        pass

    # Set up the figure
    figure_key = map_region + str(map_par) + str(fig_par)
    try:
        figure = _figure_cache[figure_key]
    except KeyError:
        figure = _figure_cache[figure_key] = _get_figure(
            map_region=map_region, map_par=map_par, fig_par=fig_par)

    # Get fig info and make active
    fig = figure.fig
    m = figure.m
    ax = figure.fig.axes[0]
    plt.figure(fig.number)
    plt.axes(ax)

    # prepare the data
    trjs = T['Trajectories']
    rel = rel_i + 1  # account for zero indexing

    # extract only releases of interest
    t = trjs[np.where(trjs[:, 0] == rel), :][0]

    # Get the data for the days_back we're interested in
    lon = t[:days_back, 2]
    lat = t[:days_back, 3]
    zlevel = t[:days_back, 4]
    zsize = np.ones(len(lon)) * 50
    marker = 'o'

    # clear the previous track
    if 'circles' in figure.keys():
        del figure['circles']
    if overlay is False:
        del ax.collections[figure.indices.collections:]
        del ax.texts[figure.indices.texts:]

    # plot the track
    cx, cy = m(lon, lat)
    if draw_circles:
        cmap = plt.get_cmap('gist_gray')
        circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                            marker=marker, edgecolor=None,
                            zorder=10, alpha=0.85)

        try:
            figure.circles = circles
        except:
            pass
        # make the figure active again
        plt.figure(fig.number)
        # draw the legend and title
        # CREATE COLORBAR
        # does a colorbar already exist?
        # Get the current axes, and properties for use later
        ax0 = fig.axes[1]
        pos = ax0.get_position()
        l, b, w, h = pos.bounds
        if cbar2:
            try:
                cb2 = figure.cb2
                cax2 = figure.cax2
                cb2.update_normal(circles)
            except:
                cax2 = plt.axes([l + w + 0.03, b, 0.02, h])
                # using im2, not im (hack to prevent colors from being
                # too compressed at the low end on the colorbar - results
                # from highly nonuniform colormap)
                cb2 = fig.colorbar(circles, cax=cax2)  # draw colorbar
                figure.cax2 = cax2
                figure.cb2 = cb2
        # check if a colorbar legend exists (this will cause
        # problems for subplot routines!)
        else:
            # ERROR: cax2 is unitialized here
            cax2 = plt.axes([l + w + 0.03, b, 0.025, h - 0.2])
            # ERROR: What's jnkmap?
            cb2 = fig.colorbar(jnkmap, cax=cax2)  # draw colorbar

        p_leg = mpl.font_manager.FontProperties(size='6')
        if cbar2_title:
            # ERROR: What is cax?
            cax.set_title(cbar2_title, fontproperties=p_leg)
        else:
            cax2.set_title('altitude\n(m)', fontproperties=p_leg)
            # delete the ghost instance
            # plt.close(jnkfig.number)
            # del jnkax, jnkfig, jnkmap

    plt.axes(ax)
    plt.setp(ax, xticks=[], yticks=[])

    if draw_labels:
        day_labels = _gen_daylabels(t[:days_back, 1], dt=3600)
        p_cax = mpl.font_manager.FontProperties(size='10',
                                                style='italic',
                                                weight='bold',
                                                )

        for i, (p, x, y) in enumerate(zip(day_labels, cx, cy)):
            if x > m.llcrnrx and x < m.urcrnrx:
                if y > m.llcrnry and y < m.urcrnry:
                    ax.text(x, y + y * .02, '{0}'.format(p), va='bottom',
                            ha='left',
                            fontproperties=p_cax, zorder=11,
                            color='white',
                            bbox=dict(facecolor='green', alpha=0.5)
                            )

    figure.fig = fig
    figure.m = m
    figure.ax = ax
    if draw_circles:
        # ERROR: cax2 can be uninitialized here
        figure.cax2 = cax2
        # ERROR: cb2 can be uninitialized here
        figure.cb2 = cb2
    return figure


def _gen_daylabels(P, H=None, dt=None):
    """
    Uses H.loutstep to calculate time steps back for plotting on clusters
    and trajectories. If dt = 86400, then labels will show 'days back'.
    """
    if dt is None and H is None:
        dt = 86400

    elif H:
        dt = abs(H.loutstep)

    if isinstance(P, int):
        return str(1 + int(abs(P)) / dt)

    else:
        # return [str(int(abs(p))) for p in P]
        return [str(1 + int(abs(p)) / dt) for p in P]


def _gen_log_clevs(dat_min, dat_max):
    """Creates a logarithmic color scale."""

    if dat_max > 0:
        dmx = int(np.round(np.log10(dat_max))) + 1
    else:
        dmx = 1

    # TODO: What's the default value of dmn?
    if dat_min > 0:
        dmn = int(np.round(np.log10(dat_min)))
    elif dat_min == 0. or np.isnan(dat_min):
        dmn = dmx - 3

    # create equally spaced range
    # ERROR: dmn could be uninitialized
    if dmx == dmn:
        dmx = dmn + 1
    clevs = np.logspace(dmn, dmx, 100)

    return clevs


def _gen_flexpart_colormap(ctbfile=None, colors=None):
    """Generate the ast colormap for FLEXPART."""

    from matplotlib.colors import ListedColormap
    if ctbfile:
        try:
            colors = np.loadtxt(ctbfile)
        except:
            print("WARNING: cannot load ctbfile. using colors")
    if colors:
        name = 'user_colormap'
    else:
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
