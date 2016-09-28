#!/usr/bin/env python
# -*- coding: utf-8 -*-
# John F Burkhart - 2010
# Licence : this code is released under the matplotlib license
""" Matplotlib Basemap Tool Suite """

from __future__ import print_function

import sys
import os
import os.path
import re

import yaml
import numpy as np
from PIL import Image

import matplotlib as mpl
import matplotlib.image
import matplotlib.pyplot as plt
from mpl_toolkits import basemap

# mpl.use("Agg")
# mp.interactive(False)
# mpl.use('Agg')

__author__ = "John F Burkhart <jfburkhart@gmail.com>"
__version__ = "0.03"


class Structure(dict):

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def set_with_dict(self, D):
        """ set attributes with a dict """
        for k in D.keys():
            self.__setattr__(k, D[k])


def map_regions(map_region='default', map_par=None, fig_par=None):
    """Given a `map_region`, return the associated parameters in mapping DB.

    USAGE::
        map_par, fig_par = map_regions(map_region="polarcat")

    The list of regions know by reflexible by default is located in the
    package file named 'mapping_db.yml'.  If you want to create a new region,
    or override an existing one in reflexible itself, you can create it in
    your own YAML file.  For example, suppose that you have the next
    'myregions.yml' file::

      northern_hemisphere:
        descr: Northern Hemisphere
        alias: my_own_nh
        map_par:
            projection: cyl
            llcrnrlat: 0.5
            urcrnrlat: 90
            llcrnrlon: -180.7
            urcrnrlon: 180
            resolution: c
            # anchor: W
        fig_par:
            figsize: [8, 3]  # w,h tuple
            axlocs: [0.1, 0.1, .7, .8]

    For informing reflexible on where your mapping file is located, you just
    create the REFLEXIBLE_MAPDB environment variable with its path::

      $ export REFLEXIBLE_MAPDB = $HOME/my_analysis/myregions.yml

    Since this moment on, the definitions in your file will be used (with
    highest priority) for finding the regions specified in `map_region`.

    Returns
      Two dictionaries, first a map_par dictionary with keywords that are the
      same as what is need to create a basemap instance of matplotlib:

      ============      ==========================
      keys              description
      ============      ==========================
      llcrnrlat         lower left latitude
      llcrnrlon         lower left longitude
      urcrnrlat         upper right latitude
      urcrnrlon         upper right longitude
      area_thresh       area threshold
      resolution        resolution
      projection        projection
      lat_1             lat_1
      lon_0             lon_0
      rsphere           (6378137.00,6356752.3142)
      m                 you can pass an m object
                        which is needed for some
                        regions
      ============      ==========================

      Second, a fig_par dictionary that contains options that may be passed to
      the :mod:`matplotlib.pyplot` :func:`figure` function:

      ============      ==========================
      keys              description
      ============      ==========================
      figsize           size of the figure
      axlocs            locations of the axes
      ============      ==========================

      .. note::
          You can override the returned map_par and fig_par dicts by passing
          dicts through the optional `map_par` and `fig_par` parameters.
    """
    # Set some default values
    map_par_ = Structure()
    fig_par_ = Structure()
    map_par_.anchor = 'C'
    fig_par_.figsize = [8, 7]   # w,h tuple
    fig_par_.axlocs = [0.05, 0.01, .8, .9]   # rect = l,b,w,h

    # Get the database out of the system YAML file
    mapdb_file = os.path.join(os.path.dirname(__file__), 'mapping_db.yml')
    with open(mapdb_file) as mapdb:
        mapping_db = yaml.load(mapdb)

    # and merge it with a possible one pointed by REFLEXIBLE_MAPDB env var
    if 'REFLEXIBLE_MAPDB' in os.environ:
        user_mapdb_file = os.environ['REFLEXIBLE_MAPDB']
        with open(user_mapdb_file) as mapdb:
            mapping_db.update(yaml.load(mapdb))

    # Lookup the region and its aliases
    try:
        region = mapping_db[map_region]
    except KeyError:
        # Lookup aliases
        for key in mapping_db:
            if 'alias' in mapping_db[key]:
                alias = mapping_db[key]['alias']
                if map_region in re.split(',\s*', alias):
                    region = mapping_db[key]
                    break
        else:
            raise KeyError("region {} not found".format(map_region))

    # Get the mapping params (should be always there)
    map_region_par = region['map_par']
    map_par_.set_with_dict(map_region_par)

    # Get the figure params (not always present)
    if 'fig_par' in region:
        fig_par_.set_with_dict(region['fig_par'])

    # Override params if `map_par` or `fig_par` are passed
    if map_par is not None:
        map_par_.set_with_dict(map_par)
    if fig_par is not None:
        fig_par_.set_with_dict(fig_par)

    # print(map_region, map_par_, fig_par_)
    try:
        del map_par_['m']  # in case
    except:
        pass
    return map_par_, fig_par_


def draw_grid(m, xdiv=10., ydiv=5., location=[1, 0, 0, 1],
              linewidth=0.5, color='k'):
    """
    draw parallels and meridians on a map.

    Currently draws labels on left, right and bottom of map.

    Keyword arguments:

        =========       ===========================
        keyword         description
        =========       ===========================
        xdiv            The division between lon
        ydiv            The division between lat
        location        Where to draw labels
        linewidth       Width of the grid lines.
        color           color of the grid lines
        =========       ===========================


    """
    # Label properties
    p_leg = mpl.font_manager.FontProperties(size='8')

    pd_options = [0.0001, 0.001, 0.01, 0.1, 0.2, 0.5, 1, 2.5, 5, 10, 20]
    md_options = [0.0001, 0.001, 0.01, 0.1, 0.2, 0.5, 1, 2.5, 5, 10, 20]

    xdiff = np.abs(m.urcrnrlon - m.llcrnrlon)
    ydiff = np.abs(m.urcrnrlat - m.llcrnrlat)

    if m.projection in ['npstere', 'spstere']:
        ydiff = 90. - np.abs(m.urcrnrlat)
        maxlat = 90.
    else:
        maxlat = None

    # setup map to have meridians:
    if xdiff > xdiv:
        md = np.round((xdiff / xdiv))
    else:
        md = np.round((xdiff / xdiv), 1)

    md = md_options[(np.abs(np.array(md_options) - md)).argmin()]
    # print(md)
    if m.projection in ['npstere', 'spstere']:
        meridians = np.arange(-180, 181, 20)
    else:
        meridians = np.arange(m.llcrnrlon, m.urcrnrlon, md)
    m_m = m.drawmeridians(meridians, labels=location,
                          linewidth=linewidth, color=color,
                          fontproperties=p_leg)

    # setup map to have parallels.
    if ydiff > ydiv:
        pd = np.round((ydiff / ydiv))
    else:
        pd = np.round((ydiff / ydiv), 2)
    if not maxlat:
        maxlat = m.urcrnrlat

    pd = pd_options[(np.abs(np.array(pd_options) - pd)).argmin()]
    if m.projection in ['npstere', 'spstere']:
        parallels = np.arange(30, 91, 10)
    else:
        if pd > 0.01:
            parallels = np.arange(np.round(m.llcrnrlat), maxlat, pd)
        else:
            parallels = np.arange(m.llcrnrlat, maxlat, pd)

    print(parallels, m.llcrnrlat, maxlat, pd)
    m_p = m.drawparallels(parallels, labels=location,
                          linewidth=linewidth, color=color,
                          fontproperties=p_leg)
    return m_p, m_m


def get_FIGURE(fig=None, ax=None, m=None, map_region=None,
               getm=True, map_par=None, fig_par=None,
               image=None):
    """
    This is a core function, used throughout this module. It is called
    by various functions. The idea is that I create a :class:`Structure` that
    contains the figure, ax, and m instance. I also add a field for "indices".
    I'm not sure this all makes the most sense, but it is what I came up with
    in order to be able to reuse figures. This saves a huge amount of time, as
    creating then basemap instance can be time consuming.
    This whole concept overall really needs to be reviewed!

    .. note::
        Generally you won't use this function direction.

    USAGE::

        > FIG = get_FIGURE()
        or
        > FIG = get_FIGURE(map_region='polarcat')

    Returns
       This will return the "FIG" object, which has attributes: `fig`, `ax`,
       `m`, and `indices`. The indices are used for deleting lines, texts,
       collections, etc. if and when we are reusing the figure instance. The
       indices basically give us a reference to the *empty* map, so we can
       delete lines without losing meridians or parallels for example.


    ============      ======================================
    keys              description
    ============      ======================================
    fig               a pyplot.fig instance, use
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
    FIGURE = Structure()

    if getm:
        if m is None:
            if image:
                fig, m = get_base_image(image, map_region=map_region,
                                        map_par=map_par,
                                        fig_par=fig_par,
                                        )
            else:

                fig, m = get_base1(map_region=map_region,
                                   map_par=map_par,
                                   fig_par=fig_par,
                                   fig=fig,
                                   )

            FIGURE.fig = fig
            FIGURE.m = m
            FIGURE.ax = fig.gca()
    elif m is None:
        FIGURE.m = None
    else:
        FIGURE.m = m

    if fig is None:
        FIGURE.fig = plt.figure()
        fig = FIGURE.fig
    else:
        FIGURE.fig = fig

    if ax is None:
        FIGURE.ax = fig.gca()
        ax = FIGURE.ax
    else:
        FIGURE.ax = ax

    FIGURE.indices = Structure()
    FIGURE.indices.texts = len(FIGURE.ax.texts)
    FIGURE.indices.images = len(FIGURE.ax.images)
    FIGURE.indices.collections = len(FIGURE.ax.collections)
    FIGURE.indices.lines = len(FIGURE.ax.lines)

    print("Using figure: %s" % FIGURE.fig.number)
    return FIGURE


def get_base1(map_region=1,
              figname=None,
              fig=None,
              drawlsmask=False,
              map_par=None,
              fig_par=None):
    """
    Primarily an internally used function, creates a
    basemap for plotting. Returns a fig object and
    a basemap instance.

    Usage::

      >fig, m=get_base1(map_region="region_name")
    """

    # Use map_regions function to define
    # input parameters for Basemap
    map_par_sd, fig_par_sd = map_regions(map_region=map_region, map_par=map_par)
    if map_par:
        map_par_sd.set_with_dict(map_par)
    if fig_par:
        fig_par_sd.set_with_dict(fig_par)

    # create the figure
    if fig is None:
        axlocs = fig_par_sd.pop('axlocs')
        fig = plt.figure(**fig_par_sd)
        ax = fig.add_axes(axlocs)
    else:
        ax = fig.gca()

    m = basemap.Basemap(**map_par_sd)
    print('getting base1')
    print(m)
    plt.axes(ax)   # make sure axes ax are current
    # draw coastlines and political boundaries.
    m.drawcoastlines(linewidth=0.8)
    m.drawcountries(linewidth=0.2)
    m.drawstates(linewidth=0.2)
    if drawlsmask:
        m.drawlsmask(ocean_color='#008EBA', zorder=0)
    # m.fillcontinents(zorder=0.)
    # draw parallels and meridians
    # use draw_grid function
    m_p, m_m = draw_grid(m)

    if figname is not None:
        plt.savefig(figname)
    return fig, m


def get_base_image(imagefile, **kwargs):
    """Warps NASA Blue Marble Image version.

    Create basemap figure for plotting on top of
    returns a figure and a basemap instance.

    Usage::

        >fig,m=get_base1(map_region="myregion")


    """
    # Get Keyword Arguments
    reg = kwargs.get('map_region', 'default')
    figname = kwargs.get('figname', None)
    fig = kwargs.get('figure', None)

    # shows how to warp an image from one map projection to another.
    # image from http://visibleearth.nasa.gov/

    # read in jpeg image to rgba array of normalized floats.
    pilImage = Image.open(imagefile)
    rgba = mpl.image.pil_to_array(pilImage)
    rgba = rgba.astype(np.float32) / 255.  # convert to normalized floats.

    # define lat/lon grid that image spans (projection='cyl').
    nlons = rgba.shape[1]
    nlats = rgba.shape[0]
    delta = 360. / float(nlons)
    lons = np.arange(-180. + 0.5 * delta, 180., delta)
    lats = np.arange(-90. + 0.5 * delta, 90., delta)

    # create new figure
    fig = plt.figure(1, figsize=(8, 6))

    # define Lambert Conformal basemap for North America.
    mr = {'map_region': reg}
    mp, fig_par = map_regions(map_region=reg)

    m = basemap.Basemap(**mp)
    ax = fig.add_axes(
        [0.1, 0.1, 0.7, 0.7])  # need to change back to [0.1,0.1,0.7,0.7]
    plt.axes(ax)  # make the original axes current again
    # transform to nx x ny regularly spaced native projection grid
    # nx and ny chosen to have roughly the same horizontal res as original
    # image.
    dx = 2. * np.pi * m.rmajor / float(nlons)
    nx = int((m.xmax - m.xmin) / dx) + 1
    ny = int((m.ymax - m.ymin) / dx) + 1
    rgba_warped = np.zeros((ny, nx, 4), np.float64)
    # interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
    try:
        for k in range(4):
            rgba_warped[:, :, k] = m.transform_scalar(
                rgba[:, :, k], lons, lats, nx, ny)
    except:
        rgba_warped = rgba
        print('problem with transform_scalar')
    # plot warped rgba image.
    im = m.imshow(rgba_warped)
    # draw coastlines.
    m.drawcoastlines(linewidth=0.5, color='0.5')
    # draw parallels and meridians.
    draw_grid(m, linewidth=0.5, color='0.5')
    # draw()
    return fig, m


if __name__ == '__main__':
    # Some tests for reading the mapping YAML database
    map_regions(map_region=sys.argv[1])
