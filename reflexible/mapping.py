#!/usr/bin/env python 
""" Matplotlib Basemap Tool Suite """

# -*- coding: utf-8 -*-
# John F Burkhart - 2010
# Licence : this code is released under the matplotlib license

from __future__ import print_function

__author__ = "John F Burkhart <jfburkhart@gmail.com>"
__version__ = "0.03"

import numpy as np
import math
import os.path
import yaml
import re

import datetime as dt
import pdb
# from matplotlib import interactive, use
# use("Agg")
# interactive(False)
import matplotlib as mpl
from matplotlib import colors, cm
from matplotlib.collections import LineCollection
# mpl.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset as NetCDFFile

TEX = False
# !!NEED TO FIX!! #
from matplotlib.ticker import NullFormatter


class Structure(dict):

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value

    def set_with_dict(self, D):
        """ set attributes with a dict """
        for k in D.keys():
            self.__setattr__(k, D[k])


class KML_File:
    "For creating KML files used for Google Earth"
    __author__ = ("Jon Goodall <jon.goodall@gmail.com> "
                  "http://www.duke.edu/~jgl34 "
                  "modifications, John Burkhart <jfburkhart@gmail.com")
    __version__ = "0.0.2"
    __license__ = ""
    __copyright__ = ""

    def __init__(self, filepath):
        self.filepath = filepath
        "adds the kml header to a file (includes a default style)"
        file = open(filepath, "w")
        file.write(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            "<kml xmlns=\"http://earth.google.com/kml/2.0\">\n"
            "<Document>\n"
            "<Style id='normalPlaceMarker'>\n"
            "  <IconStyle>\n"
            "    <Icon>\n"
            "      <href>root://icons/palette-3.png</href>\n"
            "      <x>96</x>\n"
            "      <y>160</y>\n"
            "      <w>32</w>\n"
            "      <h>32</h>\n"
            "    </Icon>\n"
            "  </IconStyle>\n"
            "</Style>\n")
        file.close()

    def close(self):
        file = open(self.filepath, "a")
        file.write(
            "</Document>\n"
            "</kml>")
        file.close()

    def open_folder(self, name):
        file = open(self.filepath, "a")
        file.write(
            "<Folder>\n" \
            "  <name>" + name + "</name>\n")
        file.close()

    def close_folder(self):
        file = open(self.filepath, "a")
        file.write(
            "</Folder>\n")
        file.close()

    def add_placemarker(self, latitude, longitude, altitude=0.0,
                        description=" ", name=" ",
                        range=6000, tilt=45, heading=0):
        "adds the point to a kml file"

        file = open(self.filepath, "a")
        if not np.iterable(longitude) and not np.iterable(latitude):
            point = True
            print('no lines')
        else:
            point = False
            print('lines!')

        if point:
            file.write(
                "<Placemark>\n"
                "  <description>" + description +
                "</description>\n"
                "  <name>" + name + "</name>\n"
                "  <styleUrl>#normalPlaceMarker</styleUrl>"
                "  <LookAt>\n"
                "    <longitude>" + str(longitude) + "</longitude>\n"
                "    <latitude>" + str(latitude) + "</latitude>\n"
                "    <range>" + str(range) + "</range>\n"
                "    <tilt>" + str(tilt) + "    </tilt>\n"
                "    <heading>" + str(heading) + "</heading>\n"
                "  </LookAt>\n"
                "  <visibility>0</visibility>\n"
                "   <Point>\n"
                "    <extrude>1</extrude>\n"
                "    <altitudeMode>relativeToGround</altitudeMode>\n"
                "    <coordinates>" + str(longitude) + "," + str(latitude) +
                ", " + str(altitude) + "</coordinates>\n"
                "   </Point>\n"
                " </Placemark>\n")
            file.close()

        else:
            ## let's create line strings

            file.write(
                "<Placemark>\n" \
                "  <description>" + description + "</description>\n"
                "  <name>" + name + "</name>\n"
                "  <styleUrl>#normalPlaceMarker</styleUrl>"
                "  <visibility>1</visibility>\n"
                "       <LineString>\n"
                "         <extrude>0</extrude>\n"
                "           <tessellate>0</tessellate>\n"
                "              <coordinates>\n")
            for i, lat in enumerate(latitude):
                file.write("%s,%s,0 " % (longitude[i], lat))
            file.write("\n")
            file.write(
                "              </coordinates>\n"
                "       </LineString>\n"
                " </Placemark>\n")
            file.close()


def map_regions(map_region='default', map_par=None, fig_par=None):
    """
    This is a core function, used throughout this module. It is called
    by the :func:`get_FIGURE` function. If you want to create a new
    region, this is where it should be added. Follow the protocol for the
    existing regions below. The idea is to move these regions out of here
    eventually, and have them be called from a file or database.

    .. note::
        Generally, I just use region names anymore, and not projection. It is
        easiest to name your new region uniquely, then simply define
        region="myregion" when you call one of the plotting routines. For the most
        part they take "region" as a keyword.

    USAGE::

        > map_par, fig_par = map_regions(map_region="POLARCAT")

    Returns
      Two dictionaries, first  a map_par dictionary with keywords that are the
      same as what is need to create a
      `basemap`<http://matplotlib.sourceforge.net/basemap/doc/html/users/mapsetup.html>_ instance.

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

      Second, a fig_par dictionary that contains options that may be passed to the
      :mod:`matplotlib.pyplot` :func:`figure` function.

      ============      ==========================
      keys              description
      ============      ==========================
      figsize           size of the figure
      axlocs            locations of the axes
      ============      ==========================

      .. note::
          You can override the fig_par.figsize in your region definition.



    """
    # Set some default values
    map_par_ = Structure()
    if map_par is not None:
        map_par_.set_with_dict(map_par)
    fig_par_ = Structure()
    if fig_par is not None:
        fig_par_.set_with_dict(fig_par)
    map_par, fig_par = map_par_, fig_par_

    map_par.anchor = 'C'
    fig_par.figsize = [8, 7]  # w,h tuple
    # rect = l,b,w,h
    fig_par.axlocs = [0.05, 0.01, .8, .9]

    # Get the database out of the YAML file
    mapdb_file = os.path.join(os.path.dirname(__file__), 'mapping_db.yml')
    with open(mapdb_file) as mapdb:
        mapping_db = yaml.load(mapdb)

    # Lookup the region and its aliases
    try:
        region = mapping_db[map_region]
    except KeyError:
        # Lookup aliases
        for key in mapping_db:
            if 'alias' in mapping_db[key]:
                alias = mapping_db[key]['alias']
                if map_region in re.split(',\s*', alias):
                    map_region = key
                    region = mapping_db[key]
                    break
        else:
            raise KeyError("region {} not found".format(map_region))

    # Get the mapping params (should be always there)
    map_region_par = region['map_par']
    map_par.set_with_dict(map_region_par)

    # Get the figure params (not always present)
    if 'fig_par' in region:
        fig_par.set_with_dict(region['fig_par'])

    #print(map_region, map_par, fig_par)
    # map_par.boundinglat=50.
    try:
        del map_par['m']  # in case
    except:
        pass
    return map_par, fig_par


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
    MD = m.drawmeridians(meridians, labels=location,
                         linewidth=linewidth, color=color, fontproperties=p_leg)

    # setup map to have parallels.
    if ydiff > ydiv:
        pd = np.round((ydiff / ydiv))
    else:
        pd = np.round((ydiff / ydiv), 2)
    if not maxlat:
        maxlat = m.urcrnrlat

    pd = pd_options[(np.abs(np.array(pd_options) - pd)).argmin()]
    # print(pd)
    if m.projection in ['npstere', 'spstere']:
        parallels = np.arange(30, 91, 10)
    else:
        if pd > 0.01:
            parallels = np.arange(np.round(m.llcrnrlat), maxlat, pd)
        else:
            parallels = np.arange(m.llcrnrlat, maxlat, pd)

    print(parallels, m.llcrnrlat, maxlat, pd)
    MP = m.drawparallels(parallels, labels=location,
                         linewidth=linewidth, color=color, fontproperties=p_leg)
    return MP, MD


def plot_ECMWF(nc, variable, time, level, map_region='NorthAmerica'):
    """ uses :function:`plot_grid`_ to plot a netcdf file downloaded
    from ECMWF MARS data store.
    """

    if isinstance(nc, str):
        nc = NetCDFFile(nc)

    alldata = nc.variables[variable][:]
    # read the data desired, but also reverse lats
    # we need lats increasing for the transform
    lats = nc.variables['latitude'][:]
    lat2 = lats[::-1]
    lons = nc.variables['longitude'][:]
    lon1 = np.arange(-180, 180, 360. / len(lons))
    data = alldata[time, level, :, :]
    data = data[::-1, :]  # flip lats
    # shift the grid to go from -180 to 180
    d2, lon2 = shiftgrid(-180, data, lon1)
    plot_grid((lon2, lat2, d2), map_region=map_region)


def default_netcdf(nco_filename,
                   lon0=-179.5, lat0=-89.5,
                   nx=720, ny=360,
                   dx=0.5, dy=0.5,
                   title="",
                   institution="Norwegian Institute for Air Research",
                   source="FLEXPART Model Data",
                   author="John F. Burkhart",
                   contact="john.burkhart@nilu.no",
                   comment="",
                   references="",
                   history=""
                   ):
    """ A function to add default information to a netcdf file.
    Adds default attributes to the ncfile using keyword arguments.

    Usage::
        nco = mp.default_netcdf('test.nc',lon0=H.outlon0,lat0=H.outlat0,
                                nx=H.nxmax,ny=H.nymax,dy=H.dyout,dx=H.dxout)

        then, nco.variables['data'][:] = np.random.rand(360,180)

    Keyword arguments:

        =============       ===========================
        keyword             description
        =============       ===========================
        lon0                The division between lon
        lat0                The division between lat
        dx                  lon spacing
        dy                  lat spacing
        title               A data name title
        institution         institute of data origin
        source              source of data origin
        author              author of data origin
        contact             contact of data origin
        comment             comment of data origin
        references          references of data origin
        history             history of data origin
        =============       ===========================


    """
    nco = NetCDFFile(nco_filename, 'w')
    nco.author = author
    nco.createdate = dt.datetime.now().ctime()
    nco.contact = contact
    nco.Conventions = "CF-1.4"
    nco.institution = institution
    nco.title = title
    nco.source = source
    nco.comment = comment
    nco.references = references
    nco.history = history

    nco.createDimension('lon', nx)
    nco.createDimension('lat', ny)
    nco.createVariable('lat', 'd', ('lat',))
    nco.createVariable('lon', 'd', ('lon',))
    lon = np.arange(lon0, lon0 + (nx * dx), dx)
    lat = np.arange(lat0, lat0 + (ny * dy), dy)
    nco.variables['lat'][:] = lat
    nco.variables['lat'].long_name = 'latitude'
    nco.variables['lat'].units = 'degrees_north'

    nco.variables['lon'][:] = lon
    nco.variables['lon'].long_name = 'longitude'
    nco.variables['lon'].units = 'degrees_east'

    # nco.createVariable('data','d',('lon','lat'))

    return nco


def grid_to_netcdf(H, D, nco_filename, spc, time):
    """INCOMPLETE, NOT YETWORKING
    Dump a grid, returned from :module:`rf.get_grid` into a netcdf file.

    USAGE::
        nco = grid_to_netcdf(H,D)
        nco.sync()
        nco.close()
    """
    print("Not working")

    return


"""


    if H.direction == 'forward':
        data = 'concentration'
    else:
        data = 'sensitivity'

    for level in D:
        if level == 0:
            units = 'mg/m2'
        else:
            units = 'ug/m2'

    nco = mp.default_netcdf(nco_filename,lon0=H.outlon0,lat0=H.outlat0,
                                    nx=H.nxmax,ny=H.nymax,dy=H.dyout,dx=H.dxout)

    # Define attributes and time variable
    nco.createDimension('lev',len(H.outheight)
    nco.createVariable('lev','d',('lev',))
    nco.variables['lev'][:] = H.outheight
    nco.variables['lev'].units = 'hPa'
    nco.variables['lev'].long_name = 'pressure level'
    nco.createDimension('time',1)
    nco.createVariable('time','d',('time',))
    nco.variables['time'][:] = d

    nco.createVariable(data,'d',('lev','lat','lon'))
    nco.variables[data].units = units
    nco.variables[data].standard_name = "SO2 concentration"
    nco.standard_name = "FLEXPART product: %s" % spc
    nco.coordinates = 'latitude,longitude'

    return nco

"""


def get_FIGURE(fig=None, ax=None, m=None, map_region=None,
               getm=True, map_par=None, fig_par=None,
               image=None):
    """
    This is a core function, used throughout this module. It is called
    by various functions. The idea is that I create a :class:`Structure` that 
    contains the figure, ax, and m instance. I also add a field for "indices".
    I'm not sure this all makes the most sense, but it is what I came up with in
    order to be able to reuse figures. This saves a huge amount of time, as
    creating then basemap instance can be time consuming.
    This whole concept overall really needs to be reviewed!

    .. note::
        Generally you won't use this function direction.

    USAGE::

        > FIG = get_FIGURE()
        or
        > FIG = get_FIGURE(map_region='POLARCAT')

    Returns
       This will return the "FIG" object, which has attributes: `fig`, `ax`, 
       `m`, and `indices`. The indices are used for deleting lines, texts,
       collections, etc. if and when we are reusing the figure instance. The 
       indices basically give us a reference to the *empty* map, so we can delete
       lines without losing meridians or parallels for example.


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
        if m == None:
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
    elif m == None:
        FIGURE.m = None
    else:
        FIGURE.m = m

    if fig == None:
        FIGURE.fig = plt.figure()
        fig = FIGURE.fig
    else:
        FIGURE.fig = fig

    if ax == None:
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

      >fig,m=get_base1(map_region="regionname")
    """

    ## Use map_regions function to define
    ## input paramters for Basemap
    map_par_sd, fig_par_sd = map_regions(map_region=map_region, map_par=map_par)
    if map_par:
        map_par_sd.set_with_dict(map_par)
    if fig_par:
        fig_par_sd.set_with_dict(fig_par)

    # create the figure
    if fig == None:
        axlocs = fig_par_sd.pop('axlocs')
        fig = plt.figure(**fig_par_sd)
        ax = fig.add_axes(axlocs)
    else:
        ax = fig.gca()

    print(Basemap)
    m = Basemap(**map_par_sd)
    print('getting base1')
    print(m)
    plt.axes(ax)  ## make sure axes ax are current
    ## draw coastlines and political boundaries.
    m.drawcoastlines(linewidth=0.8)
    m.drawcountries(linewidth=0.2)
    m.drawstates(linewidth=0.2)
    if drawlsmask:
        m.drawlsmask(ocean_color='#008EBA', zorder=0)
    # m.fillcontinents(zorder=0.)
    ## draw parallels and meridians.
    ## use draw_grid function
    MP, MD = draw_grid(m)

    if figname != None:
        plt.savefig(figname)
    return fig, m


def get_base2(**kwargs):
    """ Warps NASA Blue Marble Image version
    Create basemap figure for plotting on top of
    returns a figure and a basemap instance:
    usage: fig,m=get_base1(map_region=2) """
    ## Get Keyword Arguments
    if 'map_region' in kwargs.keys():
        reg = kwargs['map_region']
    else:
        reg = 1
    if 'figname' in kwargs.keys():
        figname = kwargs['figname']
    else:
        figname = None
    if 'figure' in kwargs.keys():
        fig = kwargs['figure']
    else:
        fig = None
    ##IMPORTS
    import pylab as P
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        from matplotlib.toolkits.basemap import Basemap
    from matplotlib.image import pil_to_array
    from PIL import Image

    # shows how to warp an image from one map projection to another.
    # image from http://visibleearth.nasa.gov/

    # read in jpeg image to rgba array of normalized floats.
    pilImage = Image.open('land_shallow_topo_2048.jpg')
    rgba = pil_to_array(pilImage)
    rgba = rgba.astype(P.float32) / 255.  # convert to normalized floats.

    # define lat/lon grid that image spans (projection='cyl').
    nlons = rgba.shape[1]
    nlats = rgba.shape[0]
    delta = 360. / float(nlons)
    lons = P.arange(-180. + 0.5 * delta, 180., delta)
    lats = P.arange(-90. + 0.5 * delta, 90., delta)

    # create new figure
    fig = P.figure(1, figsize=(8, 6))

    # define Lambert Conformal basemap for North America.
    mr = {'map_region': reg}
    mp, fig_par = map_regions(map_region=reg)

    m = Basemap(**mp)
    ax = fig.add_axes(
        [0.1, 0.1, 0.7, 0.7])  # need to change back to [0.1,0.1,0.7,0.7]
    P.axes(ax)  # make the original axes current again
    # transform to nx x ny regularly spaced native projection grid
    # nx and ny chosen to have roughly the same horizontal res as original image.
    dx = 2. * P.pi * m.rmajor / float(nlons)
    nx = int((m.xmax - m.xmin) / dx) + 1;
    ny = int((m.ymax - m.ymin) / dx) + 1
    rgba_warped = P.zeros((ny, nx, 4), P.Float64)
    # interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
    for k in range(4):
        rgba_warped[:, :, k] = m.transform_scalar(rgba[:, :, k], lons, lats, nx,
                                                  ny)
    # plot warped rgba image.
    im = m.imshow(rgba_warped)
    # draw coastlines.
    m.drawcoastlines(linewidth=0.5, color='0.5')
    # draw parallels and meridians.
    draw_grid(m, linewidth=0.5, color='0.5')
    # draw()
    return fig, m


def get_base3(**kwargs):
    """ Custom version for sea ice. """

    ## Get Keyword Arguments
    if 'map_region' in kwargs.keys():
        map_region = kwargs['map_region']
    else:
        map_region = 1
    if 'figname' in kwargs.keys():
        figname = kwargs['figname']
    else:
        figname = None
    if 'figure' in kwargs.keys():
        fig = kwargs['figure']
    else:
        fig = None
    if 'basefile' in kwargs.keys():
        basefile = kwargs['basefile']
    else:
        basefile = None
    # the data is interpolated to the native projection grid.
    from PIL import Image
    from matplotlib import interactive
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        from matplotlib.toolkits.basemap import Basemap
    from matplotlib.image import pil_to_array
    from pylab import title, colorbar, show, close, \
        axes, cm, load, arange, figure, \
        text, savefig, setp, draw, clf, cla

    # read in jpeg image to rgba array of normalized floats.
    if basefile != None:
        pilImage = Image.open(basefile)
        rgba = pil_to_array(pilImage)
        rgba = rgba.astype(np.float32) / 255.  # convert to normalized floats.

    interactive(False)
    ## create the figure.
    if fig == None:
        fig = figure(1, \
                     dpi=100,
                     frameon=False,
                     facecolor=None,
                     edgecolor=None,
                     figsize=(12, 9.11)
                     )
        # fig=figure(1,figsize=(57.625,43.75))
    # Use map_regions function to define input paramters for Basemap
    mp, fig_par = map_regions(map_region)
    m = Basemap(**mp)
    ax = fig.add_axes([0, 0, 1, 1],
                      frameon=False)  # need to change back to [0.1,0.1,0.7,0.7]
    axes(ax)  # make the original axes current again
    ## draw coastlines and political boundaries.
    m.drawcoastlines(linewidth=.5)
    m.drawcountries(linewidth=.5)
    m.drawstates(linewidth=.5)
    ## draw parallels and meridians.
    draw_grid(m, linewidth=0.5)

    if figname != None:
        savefig(figname)
    return fig, m


def get_base_image(imagefile, **kwargs):
    """ Warps NASA Blue Marble Image version
    Create basemap figure for plotting on top of
    returns a figure and a basemap instance.

    Usage::

        >fig,m=get_base1(map_region="myregion")


    """
    ## Get Keyword Arguments
    if 'map_region' in kwargs.keys():
        reg = kwargs['map_region']
    else:
        reg = 1
    if 'figname' in kwargs.keys():
        figname = kwargs['figname']
    else:
        figname = None
    if 'figure' in kwargs.keys():
        fig = kwargs['figure']
    else:
        fig = None
    ##IMPORTS
    import matplotlib.pyplot as plt
    try:
        from mpl_toolkits.basemap import Basemap
    except:
        from matplotlib.toolkits.basemap import Basemap
    from matplotlib.image import pil_to_array
    from PIL import Image

    # shows how to warp an image from one map projection to another.
    # image from http://visibleearth.nasa.gov/

    # read in jpeg image to rgba array of normalized floats.
    pilImage = Image.open(imagefile)
    rgba = pil_to_array(pilImage)
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

    m = Basemap(**mp)
    ax = fig.add_axes(
        [0.1, 0.1, 0.7, 0.7])  # need to change back to [0.1,0.1,0.7,0.7]
    plt.axes(ax)  # make the original axes current again
    # transform to nx x ny regularly spaced native projection grid
    # nx and ny chosen to have roughly the same horizontal res as original image.
    dx = 2. * np.pi * m.rmajor / float(nlons)
    nx = int((m.xmax - m.xmin) / dx) + 1
    ny = int((m.ymax - m.ymin) / dx) + 1
    rgba_warped = np.zeros((ny, nx, 4), np.float64)
    # interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
    try:
        for k in range(4):
            rgba_warped[:, :, k] = m.transform_scalar(rgba[:, :, k], lons, lats,
                                                      nx, ny)
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


def set_plotkwargs(kwargs, plot_kwargs=None):
    """ Internal function to check for valid plotting keywords
    and extracts them into a new kwarg dict for the plotting routine """

    valid_kwargs = ['alpha', 'animated', 'antialiased', 'axes', 'clip_box',
                    'clip_on ', 'clip_path ', 'color', 'contains',
                    'dash_capstyle',
                    'dash_joinstyle', 'dashes', 'data', 'drawstyle', 'figure ',
                    'label', 'linestyle', 'linewidth ', 'lod', 'marker',
                    'markeredgecolor', 'markeredgewidth ', 'markerfacecolor ',
                    'markersize', 'picker ', 'snap ', 'solid_capstyle',
                    'solid_joinstyle ', 'transform', 'url ',
                    'visible', 'xdata', 'ydata', 'zorder']
    if plot_kwargs is None:
        plot_kwargs = {}
    for arg in kwargs.keys():
        if arg in valid_kwargs:
            plot_kwargs[arg] = kwargs[arg]
    return plot_kwargs


def plot_track(lon, lat,
               figname=None, FIGURE=None, overlay=False, cbar2=False,
               zlevel=None, zsize=None, base=1, marker='o',
               plotargs=None,
               map_region=None,
               scatter=False,
               units=''):
    """ Plot a longitude,latitude course over a basemap. Accepts several
    keyword arguments and makes use of the FIG.

    Usage::

        > FIG = plot_track(lon,lat,**kwargs)

    ===============     ========================================
    keyword             description
    ===============     ========================================
    figname             A figurename can be passed. If
                        given the figure will be saved
                        with the given name.
    FIGURE              a FIG object from :func:`get_FIGURE`
    overlay             overlay the track
    cbar2               create a second color bar
    zlevel              if given, the track will be colored
                        accordingly. Must be equal to
                        length (lon)
    zsize               Size of points. Can be int or len(lat)
    base                which get_base function to use [1]
    marker              marker style
    plotargs            dictionary of plot arguments
    map_region              a region name
    projection          a projection
    units               a string than can be passed.
    ===============     ========================================



    """

    ## make tick lables smaller and set parameters for colorbar
    mpl.rcParams['xtick.labelsize'] = 8.5
    mpl.rcParams['ytick.labelsize'] = 8.5
    p_cax = mpl.font_manager.FontProperties(size='7')

    # check that lon/lat are iterable
    if not np.iterable(lon):
        lon = [lon]
    if not np.iterable(lat):
        lat = [lat]
        scatter = True

    if zsize == None:
        zsize = np.ones(len(lon)) * 1000

    ## Default plot_kwargs
    plot_kwargs = {}
    if plotargs != None:
        plot_kwargs = set_plotkwargs(plotargs, plot_kwargs)
    else:
        plotargs = {'zorder': 10,
                    'alpha': 0.35,
                    'edgecolor': None,
                    # 'marker':marker  #should get rid of marker
                    }
        plot_kwargs = set_plotkwargs(plotargs, plot_kwargs)

    if FIGURE == None:
        FIGURE = get_FIGURE(map_region=map_region,)

    ##Get fig if exists
    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax
    nullfmt = NullFormatter()
    if ax == None:
        az = plt.gca()
        pos = az.get_position()
        l, b, w, h = getattr(pos, 'bounds', pos)
        ax = fig.add_axes([l, b, w, h], frameon=False)
    if overlay == False:
        del ax.collections[FIGURE.indices.collections:]

    # Make sure the axes and figure are current
    plt.figure(fig.number)
    plt.axes(ax)
    # PRINT CRUISE TRACK
    cx, cy = m(lon, lat)
    if zlevel != None:
        if scatter:
            ## Set plotting kwargs
            # if plotargs == None:
            #    plotargs = {'zorder':10,
            #                'alpha':0.35,
            #                'edgecolor':None
            #                }
            plot_kwargs = set_plotkwargs(plotargs, plot_kwargs)
            # from pylab import linspace, cm, colorbar,scatter, cla, axes
            # m.scatter(cx,cy,25*linspace(0,1,(len(cx))),zlevel,cmap=cm.spectral,marker='o',faceted=False,zorder=10)
            cmap = plt.get_cmap('jet')
            # c=m.scatter(cx,cy,zsize,zlevel,cmap=cmap,marker=marker,
            #            **plot_kwargs)
            c = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                          **plot_kwargs)
            # pos = az.get_position()
            # l, b, w, h = getattr(pos, 'bounds', pos)
            # cax = axes([l+w+0.075, b, 0.05, h])
            # c.set_alpha(0.01)
            # m.scatter has no color bar, so create a ghost 'scatter' instance
            pos = ax.get_position()
            l, b, w, h = getattr(pos, 'bounds', pos)
            jnkfig = plt.figure()
            jnkax = jnkfig.add_axes([l, b, w, h], frameon=False)
            jnkax.scatter(cx, cy, zsize, zlevel, cmap=cmap, marker=marker,
                          **plot_kwargs)
            plt.figure(fig.number)
            if cbar2 is True:
                cax = plt.axes([l + w + 0.12, b, 0.02, h - 0.035])
            else:
                cax = plt.axes([l + w + 0.03, b, 0.025, h - 0.035])
            # cax = plt.axes([l+w+0.03, b, 0.02, h])
            plt.colorbar(cax=cax)  # draw colorbar
            if TEX:
                cax.set_title(r'textit{%s}' % units)
            else:
                cax.set_title('%s' % units,
                              fontproperties=p_cax)

            # delete the ghost instance
            plt.close(jnkfig.number)
            plt.axes(ax)
        else:
            # raise IOError('Vidit is going to add this method')
            print('trying now')
            if np.rank(cx) == 1:
                points = np.array([cx, cy]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
                                    norm=plt.Normalize(zlevel.min(),
                                                       zlevel.max()))
                lc.set_array(zlevel.flatten())
                lc.set_linewidth(3)
                plt.gca().add_collection(lc)
            elif np.rank(cx) == 2:
                for i in range(cx.shape[0]):
                    points = np.array([cx[i, :], cy[i, :]]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)

                    lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
                                        norm=plt.Normalize(zlevel.min(),
                                                           zlevel.max()))
                    lc.set_array(zlevel.flatten())
                    lc.set_linewidth(3)
                    plt.gca().add_collection(lc)
            else:
                raise IOError(
                    "input array shape zlevel cannot be greater than rank-2")


    else:
        if scatter:
            m.scatter(cx, cy, **plot_kwargs)
        else:
            m.plot(cx, cy, **plot_kwargs)
        print('boring')
    plt.axes(ax)
    ax.xaxis.set_major_formatter(nullfmt)
    ax.yaxis.set_major_formatter(nullfmt)
    plt.setp(ax, xticks=[], yticks=[])
    # ax.axesPatch.set_alpha(0.0)
    if figname:
        plt.savefig(figname)
    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.ax = ax
    return FIGURE


def plot_grid(D, map_region='POLARCAT', dres=0.5,
              transform=True, figname=None, fillcontinents=False,
              points=False):
    """ plot an array over a basemap. The required argument "D" is
    either:
    * a tuple of (x,y,z)
    * a tuple with (lon,lat,grid)
    * a tuple with (grid,)

    Usage::

        > FIG = plot_grid(D,**kwargs)

    ===============     ========================================
    keyword             description
    ===============     ========================================
    dres                resolution of the grid
    fillcontinents      fills continents if True
    plotargs            dictionary of plot arguments
    figname             A figurename, if passed will save the
                        figure with figname
    points              set to True if passing a x,y,z matrix of
                        points
    map_region              A region from :func:`get_base1`
    ===============     ========================================


    """
    print("length of D: %s" % len(D))

    if isinstance(D, np.ndarray):
        # pdb.set_trace()
        assert len(D.shape) == 2, "D grid must be 2d"
        # if len(D)==1:
        # assume full earth grid with dres
        print('received grid of shape:', D.shape)
        if D.shape[0] == 720:
            lons = np.arange(-180, 180, dres)
        elif D.shape[0] == 721:
            lons = np.arange(-180, 180.01, dres)
        if D.shape[1] == 360:
            lats = np.arange(-90, 90, dres)
            lons = np.arange(-180, 180.01, dres)
        elif D.shape[1] == 361:
            lats = np.arange(-90, 90.01, dres)
            lons = np.arange(-180, 180.01, dres)
        points = False
        z = D.T

    elif len(D) == 3:
        points = True
        x = D[0]
        y = D[1]
        z = D[2]

        if len(z.shape) > 1:
            points = False
            lons = x
            lats = y

            ## CHANGED THIS HERE FOR THE CODEX ##
            ## Be sure to set map_region=None                          ##
    if isinstance(map_region, str):
        print("getting basemap with map_region: %s" % map_region)
        fig, m = get_base1(map_region=map_region)
    else:
        fig = plt.figure()
        ax = fig.add_axes()
        a = 6378.273e3
        ec = 0.081816153
        b = a * np.sqrt(1. - ec ** 2)
        m = Basemap(projection='stere', lat_0=90, lon_0=-45, lat_ts=70, \
                    llcrnrlat=33.92, llcrnrlon=279.96, \
                    urcrnrlon=102.34, urcrnrlat=31.37, \
                    rsphere=(a, b))
        #        m = Basemap(width=12000000,height=8000000,
        #                resolution='l',projection='npstere',\
        #                lat_ts=50,lat_0=50,lon_0=-107.)#Set up a basemap
        m.drawcoastlines()

    if points == True:
        # Plot individual data points
        norm = colors.normalize(z.min(), z.max())
        xpt, ypt = m(x, y)
        cmap = [cm.jet(norm(i)) for i in z]
        m.scatter(xpt, ypt, s=z / 100., edgecolors='none', color=cmap)
        # for i in range(len(y)):
        #     xpt,ypt = m(x[i],y[i])
        #     cmap = cm.jet(norm(z[i]))
        #     #cmap = 'k'
        #     m.plot([xpt],[ypt],'.',color=cmap,markersize=2)
        #     #plt.plot(x[i],y[i],'.',color=cmap,markersize=2)

    if points == False:
        # transform Z data into projection
        # transform to nx x ny regularly spaced native projection grid
        dx = 2. * np.pi * m.rmajor / len(lons)
        nx = int((m.xmax - m.xmin) / dx) + 1;
        ny = int((m.ymax - m.ymin) / dx) + 1
        if transform:
            # Need this if we choose lon,lat approach
            Zt, xx, yy = m.transform_scalar(z, lons, lats, nx, ny,
                                            returnxy=True)
        else:
            Zt = z
        print(m.projection)

        if 'moll' in m.projection:
            x, y = m(xx, yy)
            m.contourf(x, y, Zt)
        else:
            m.imshow(Zt)
        if fillcontinents:
            m.fillcontinents()
        plt.colorbar()
        # plt.imshow(Z)
    if figname != None:
        # plt.ylim([40,90])
        # plt.title('data locations on mercator grid')
        plt.savefig(figname)
    else:
        plt.show()

    return fig, m


#### LateX #####
if TEX:
    mpl.rc('font', **{'family': 'sans-serif',
                      'sans-serif': ['Helvetica']
                      }
           )

    mpl.rc('text', usetex=True)


def read_satdata(c, dat, D, start_time, stop_time=None, satellite='Calipso'):
    """ Reads text files for a few defined format satellite groundtrack data

    .. note::
        Not ready 
    """
    import datetime
    if stop_time == None:
        stop_time = start_time + datetime.timedelta(1)
    if satellite == 'Calipso':
        for l in c.readlines():
            l = l.strip().split()
            if len(l) == 7:
                dat.append(l)
        for d in dat:
            t = ' '.join(d[:4])
            t = t[:-4]
            t = datetime.datetime.strptime(t, '%d %b %Y %H:%M:%S')
            if t > start_time and t < stop_time:
                lat, lon, rng = d[4:]
                D.append((t, lon, lat, rng))
    # pdb.set_trace()
    if satellite in ['Terra', 'Aqua']:
        for l in c.readlines():
            l = l.strip().split()
            if l != []:
                if l[0] != 'GMT':
                    if len(l) == 3:
                        y, m, d = l[0].split('/')
                    if len(l) > 3:
                        hr, mn, sc = l[0].split(':')
                        t = datetime.datetime(int(y), int(m), int(d), int(hr),
                                              int(mn), int(sc))
                        # print(start_time, t, stop_time)
                        if t > start_time and t < stop_time:
                            # print('yes!')
                            try:
                                lat, lon, rng = l[1:]
                                D.append((t, lon, lat, rng))
                            except:
                                lat, lon, hdg, ltlat, ltlon,\
                                rtlat, rtlon = l[1:]
                                D.append((t, lon, lat, hdg))

    t0 = datetime.datetime(2009, 4, 15)
    Sets = {}
    i = 0

    for t in D:
        if t[0] - t0 < datetime.timedelta(minutes=60):
            Sets[i].append(t)
        else:
            t0 = t[0]
            i += 1
            Sets[i] = [t]

    return Sets


def read_LARC_predict(filepath):
    """ read path prediction file from the LARC website """
    lines = open(filepath).readlines()
    out = []
    for l in lines:
        s = l.strip().split()
        if len(s) == 3:
            day = dt.datetime.strptime(s[0], '%Y/%m/%d')
            continue
        if s:
            if 'GMT' in s[0]:
                pass
            else:
                H, M, S = [int(i) for i in s[0].split(':')]
                time = dt.datetime(day.year, day.month, day.day, H, M, S)
                lat = float(s[1])
                lon = -1 * float(s[2])
                hdg = float(s[3])
                out.append([time, lon, lat, hdg])

    out = np.array(out)
    Sets = {}
    i = 0

    for t in out:
        if not Sets.keys():
            Sets[i] = [t]
            t0 = t[0]
        else:
            if t[0] - t0 < dt.timedelta(minutes=60):
                Sets[i].append(t)
            else:
                t0 = t[0]
                i += 1
                Sets[i] = [t]

    return Sets


def Sat_tracks(map, start_time=None, stop_time=None,
               satellite='Calipso', ifile=None,
               plotargs=None,
               get_Sets=False, Sets=None,
               print_times=True, plot_frq=30):
    """ Plots CALIPSO/Terra Ground tracks.

    .. note::
        Not ready.

    Requires:
    map
    start_time
    stop_time
    input_data_file
    """
    import datetime, sys
    from pylab import text
    from matplotlib import colors, cm
    if ifile is None:
        if satellite == 'Calipso':
            if 'lin' in sys.platform:
                # SATELLITEDAT='/mnt/win/07_jfb/jfbin/pybin/data/CALIPSO_SpecialProduct_NorthAtlantic_20080408.txt'
                SATELLITEDAT = '/xnilu_wrk/jfb/DATASETS_TEMPLATES/CALIPSO/CALIPSO_SpecialProduct_PAM-ARCMIP_20110324.txt'
            else:
                SATELLITEDAT = 'C:\07_jfb\jfbin\pybin\data\CALIPSO_SpecialProduct_NorthAtlantic_20080408.txt'
        elif satellite == 'Terra':
            if 'lin' in sys.platform:
                SATELLITEDAT = '/mnt/win/07_jfb/jfbin/pybin/data/terra.dat'
            else:
                SATELLITEDAT = 'C:\07_jfb\jfbin\pybin\data\terra.dat'
        else:
            SATELLITEDAT = None
            print('Need satellite data file!!!')
            AttributeError
    else:
        SATELLITEDAT = ifile
        print('Using: ' + ifile)

    if plotargs == None:
        plot_kwargs = {'linewidth': 1,
                       'color': 'b',
                       'linestyle': '--'}
    else:
        plot_kwargs = set_plotkwargs(plotargs)

    if start_time == None:
        start_time = datetime.datetime(2010, 5, 11)
    if stop_time == None:
        stop_time = start_time + datetime.timedelta(30)
    if Sets == None:
        c = open(SATELLITEDAT, 'r')
        dat, D = [], []
        # doy=start_time.strftime('%Y%m%d')
        Sets = read_satdata(c, dat, D, start_time, stop_time, satellite)

    if get_Sets:
        return Sets

    # PRINT TRACK
    tinit = Sets[1][0][0]
    tend = Sets[24][0][0]
    norm = colors.normalize(0, len(Sets))
    for k, D in enumerate(Sets.values()):
        cmap = cm.jet(norm(k))
        clon = [-1 * float(i[1]) for i in D]
        clat = [float(i[2]) for i in D]
        cx, cy = map(clon, clat)
        plot_kwargs.update({'color': cmap})
        map.plot(cx, cy, **plot_kwargs)
        # Print legends
        if 1:
            label = D[0][0].strftime('%Y%m%d')

        # Print TRACK TIMES
        if print_times:
            frq = plot_frq
            tlon = [-1 * float(i[1]) for i in D[::frq]]
            tlat = [float(i[2]) for i in D[::frq]]
            tx, ty = map(tlon, tlat)
            map.plot(tx, ty, 'r+', label=label, **plot_kwargs)
            ty = [i + 1000 for i in ty]  # reposition for text placement
            tx = [i + 1000 for i in tx]
            t_txt = [i[0].strftime('%H:%M') for i in D[::frq]]

            for i in range(len(tx)):
                if (tx[i] < map.urcrnrx and
                    ty[i] < map.urcrnry and
                    tx[i] > map.llcrnrx and
                    ty[i] > map.llcrnry):
                    text(tx[i], ty[i], t_txt[i], size=8)
                    # title('%s Track: %s' % (satellite,doy))

    return map, [tinit, tend]


def greatCircleDistance(RLAT1, RLON1, RLAT2, RLON2):
    """
    C
    C PROGRAM HISTORY LOG:
    C   96-04-10  IREDELL
    C
    C USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
    C
    C   INPUT ARGUMENT LIST:
    C     RLAT1    - REAL LATITUDE OF POINT 1 IN DEGREES
    C     RLON1    - REAL LONGITUDE OF POINT 1 IN DEGREES
    C     RLAT2    - REAL LATITUDE OF POINT 2 IN DEGREES
    C     RLON2    - REAL LONGITUDE OF POINT 2 IN DEGREES
    C
    C   OUTPUT ARGUMENT LIST:
    C     DISTANCE - REAL GREAT CIRCLE DISTANCE IN KILOMETERS
    C
    C
    """
    RERTH = 6.3712E6
    PI = 3.14159265358979
    DPR = 180. / PI
    if abs(RLAT1 - RLAT2) < 0.03 and abs(RLON1 - RLON2) < 0.03:
        DISTANCE = 0.
    else:
        CLAT1 = math.cos(RLAT1 / DPR)
        SLAT1 = math.sin(RLAT1 / DPR)
        CLAT2 = math.cos(RLAT2 / DPR)
        SLAT2 = math.sin(RLAT2 / DPR)
        CDLON = math.cos((RLON1 - RLON2) / DPR)
        CRD = SLAT1 * SLAT2 + CLAT1 * CLAT2 * CDLON
        DISTANCE = RERTH * math.acos(CRD) / 1000.

    return DISTANCE


gcd = greatCircleDistance


def gridarea(xl, yl, xr, yr, method='gcd'):
    """ calculates grid area as trapezoid using great circle distance """
    # yres = (yr-yl)
    # lat = yl + yres
    # R = 6371.009
    # if method == 'badc':
    # A = (R**2)*(np.radians(xr)-np.radians(xl))*(np.sin(yr)-np.sin(yl))
    # return A
    # if method == 'drmath':
    # A = ((np.pi/180)*R**2)*abs((np.sin(yr)-np.sin(yl)))*abs(xr-xl)
    # return A
    # if method == 'radians':
    ##Area is calculated using the latitude at the upper bound of the grid-cell, using Radians
    ##lat = latitude of center of the grid cell.
    # radians = (90.0 - (lat+yres/2))*3.141593/180.0
    ##Calculate cosines:
    # cosines = np.cos(radians)-np.cos(radians + (yres*3.141593)/180.0)
    ##Calculate area in square kilometers:
    # area = (6371221.3*6371221.3*3.141593*cosines/360.0)*1.0e-6
    # return area
    if method == 'gcd':
        dxa = gcd(yr, xl, yr, xr)
        dxb = gcd(yl, xl, yl, xr)
        dy = gcd(yl, xl, yr, xl)
        # areafract = (dlon*dlat)/((xr-xl)/(yr-yl))
        # print(dxa, dxb, dy)
        area = 0.5 * dy * (dxa + dxb)
        return area


__version__ = '1.0.1'


class GreatCircle(object):
    """
   formula for perfect sphere from Ed Williams' 'Aviation Formulary'
   (http://williams.best.vwh.net/avform.htm)

   code for ellipsoid posted to GMT mailing list by Jim Leven in Dec 1999

   Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
   """

    def __init__(self, rmajor, rminor, lon1, lat1, lon2, lat2):
        """
        Define a great circle by specifying:
        rmajor - radius of major axis of ellipsoid
        rminor - radius of minor axis of ellipsoid.
        lon1 - starting longitude of great circle
        lat1 - starting latitude
        lon2 - ending longitude
        lat2 - ending latitude
        All must be given in degrees.

        Instance variables:
        distance - distance along great circle in radians.
        lon1,lat1,lon2,lat2 - start and end points (in radians).
        """
        if rmajor == 0 and rminor == 0:
            # WGS84
            a = 6378137.0
            b = 6356752.3142
            f = (a - b) / a
            rmajor = (2 * a + b) / 3.
            rminor = (2 * a + b) / 3.
        # convert to radians from degrees.
        lat1 = math.radians(lat1)
        lon1 = math.radians(lon1)
        lat2 = math.radians(lat2)
        lon2 = math.radians(lon2)
        self.a = rmajor
        self.f = (rmajor - rminor) / rmajor
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        # distance along geodesic in meters.
        d, a12, a21 = vinc_dist(self.f, self.a, lat1, lon1, lat2, lon2)
        self.distance = d
        self.azimuth12 = a12
        self.azimuth21 = a21
        # great circle arc-length distance (in radians).
        self.gcarclen = 2. * math.asin(
            math.sqrt((math.sin((lat1 - lat2) / 2)) ** 2 + \
                      math.cos(lat1) * math.cos(lat2) * (
                      math.sin((lon1 - lon2) / 2)) ** 2))
        # check to see if points are antipodal (if so, route is undefined).
        if self.gcarclen == math.pi:
            self.antipodal = True
        else:
            self.antipodal = False

    def points(self, npoints):
        """
        compute arrays of npoints equally spaced
        intermediate points along the great circle.

        input parameter npoints is the number of points
        to compute.

        Returns lons, lats (lists with longitudes and latitudes
        of intermediate points in degrees).

        For example npoints=10 will return arrays lons,lats of 10
        equally spaced points along the great circle.
        """
        # must ask for at least 2 points.
        if npoints <= 1:
            raise ValueError('npoints must be greater than 1')
        elif npoints == 2:
            return [math.degrees(self.lon1), math.degrees(self.lon2)],\
                   [math.degrees(self.lat1), math.degrees(self.lat2)]
        # can't do it if endpoints are antipodal, since
        # route is undefined.
        if self.antipodal:
            raise ValueError(
                'cannot compute intermediate points on a great circle '
                'whose endpoints are antipodal')
        d = self.gcarclen
        delta = 1.0 / (npoints - 1)
        f = delta * np.arange(npoints)  # f=0 is point 1, f=1 is point 2.
        incdist = self.distance / (npoints - 1)
        lat1 = self.lat1
        lat2 = self.lat2
        lon1 = self.lon1
        lon2 = self.lon2
        # perfect sphere, use great circle formula
        if self.f == 0.:
            A = np.sin((1 - f) * d) / math.sin(d)
            B = np.sin(f * d) / math.sin(d)
            x = A * math.cos(lat1) * math.cos(lon1) + B * math.cos(lat2) * math.cos(lon2)
            y = A * math.cos(lat1) * math.sin(lon1) + B * math.cos(lat2) * math.sin(lon2)
            z = A * math.sin(lat1) + B * math.sin(lat2)
            lats = np.arctan2(z, np.sqrt(x ** 2 + y ** 2))
            lons = np.arctan2(y, x)
            lons = map(math.degrees, lons.tolist())
            lats = map(math.degrees, lats.tolist())
        # use ellipsoid formulas
        else:
            latpt = self.lat1
            lonpt = self.lon1
            azimuth = self.azimuth12
            lons = [math.degrees(lonpt)]
            lats = [math.degrees(latpt)]
            for n in range(npoints - 2):
                latptnew, lonptnew, alpha21 = vinc_pt(self.f, self.a, latpt,
                                                      lonpt, azimuth, incdist)
                d, azimuth, a21 = vinc_dist(self.f, self.a, latptnew, lonptnew,
                                            lat2, lon2)
                lats.append(math.degrees(latptnew))
                lons.append(math.degrees(lonptnew))
                latpt = latptnew
                lonpt = lonptnew
                lons.append(math.degrees(self.lon2))
                lats.append(math.degrees(self.lat2))
        return lons, lats


#
# --------------------------------------------------------------------- 
# |                                                                    |
# |     geodetic.py -  a collection of geodetic functions              |
# |                                                                    |
# --------------------------------------------------------------------- 
# 
# 
# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual      |
# |                                                                     |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html                  |
# |                                                                     |
# | This page last updated 11 May 1999                                  |
# |                                                                     |
# | Computations on the Ellipsoid                                       |
# |                                                                     |
# | There are a number of formulae that are available                   |
# | to calculate accurate geodetic positions,                           |
# | azimuths and distances on the ellipsoid.                            |
# |                                                                     |
# | Vincenty's formulae (Vincenty, 1975) may be used                    |
# | for lines ranging from a few cm to nearly 20,000 km,                |
# | with millimetre accuracy.                                           |
# | The formulae have been extensively tested                           |
# | for the Australian region, by comparison with results               |
# | from other formulae (Rainsford, 1955 & Sodano, 1965).               |
# |                                                                     |
# | * Inverse problem: azimuth and distance from known                  |
# |                     latitudes and longitudes                        |
# | * Direct problem: Latitude and longitude from known                 |
# |                     position, azimuth and distance.                 |
# | * Sample data                                                       |
# | * Excel spreadsheet                                                 |
# |                                                                     |
# | Vincenty's Inverse formulae                                         |
# | Given: latitude and longitude of two points                         |
# |                     (phi1, lembda1 and phi2, lembda2),              |
# | Calculate: the ellipsoidal distance (s) and                         |
# | forward and reverse azimuths between the points (alpha12, alpha21). |
# |                                                                     |
# ---------------------------------------------------------------------- 

def vinc_dist(f, a, phi1, lembda1, phi2, lembda2):
    """ 

      Returns the distance between two geographic points on the ellipsoid
      and the forward and reverse azimuths between these points.
      lats, longs and azimuths are in radians, distance in metres 

      Returns ( s, alpha12,  alpha21 ) as a tuple

      """

    if (abs(phi2 - phi1) < 1e-8) and (abs(lembda2 - lembda1) < 1e-8):
        return 0.0, 0.0, 0.0

    two_pi = 2.0 * math.pi

    b = a * (1.0 - f)

    TanU1 = (1 - f) * math.tan(phi1)
    TanU2 = (1 - f) * math.tan(phi2)

    U1 = math.atan(TanU1)
    U2 = math.atan(TanU2)

    lembda = lembda2 - lembda1
    last_lembda = -4000000.0  # an impossibe value
    omega = lembda

    # Iterate the following equations, 
    #  until there is no significant change in lembda 

    while (last_lembda < -3000000.0 or
        lembda != 0 and
        abs((last_lembda - lembda) / lembda) > 1.0e-9):
        sqr_sin_sigma = pow(math.cos(U2) * math.sin(lembda), 2) +\
                        pow((math.cos(U1) * math.sin(U2) -
                        math.sin(U1) * math.cos(U2) * math.cos(lembda)), 2)

        Sin_sigma = math.sqrt(sqr_sin_sigma)

        Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lembda)

        sigma = math.atan2(Sin_sigma, Cos_sigma)

        Sin_alpha = math.cos(U1) * math.cos(U2) * math.sin(lembda) / math.sin(sigma)
        alpha = math.asin(Sin_alpha)

        Cos2sigma_m = math.cos(sigma) - (
        2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha), 2))

        C = (f / 16) * pow(math.cos(alpha), 2) * (
        4 + f * (4 - 3 * pow(math.cos(alpha), 2)))

        last_lembda = lembda

        lembda = omega + (1 - C) * f * math.sin(alpha) * (
        sigma + C * math.sin(sigma) *
        (Cos2sigma_m + C * math.cos(sigma) * (
            -1 + 2 * pow(Cos2sigma_m, 2))))

    u2 = pow(math.cos(alpha), 2) * (a * a - b * b) / (b * b)

    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B / 4) *
                                   (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) -
                                    (B / 6) * Cos2sigma_m * (
                                    -3 + 4 * sqr_sin_sigma) *
                                    (-3 + 4 * pow(Cos2sigma_m, 2))))

    s = b * A * (sigma - delta_sigma)

    alpha12 = math.atan2((math.cos(U2) * math.sin(lembda)),
                         (math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(
                             U2) * math.cos(lembda)))

    alpha21 = math.atan2((math.cos(U1) * math.sin(lembda)),
                         (
                         -math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(
                             U2) * math.cos(lembda)))

    if alpha12 < 0.0:
        alpha12 += two_pi
    if alpha12 > two_pi:
        alpha12 -= two_pi

    alpha21 = alpha21 + two_pi / 2.0
    if alpha21 < 0.0:
        alpha21 += two_pi
    if alpha21 > two_pi:
        alpha21 -= two_pi

    return s, alpha12, alpha21

    # END of Vincenty's Inverse formulae 


# ----------------------------------------------------------------------------
# Vincenty's Direct formulae                                                |
# Given: latitude and longitude of a point (phi1, lembda1) and              |
# the geodetic azimuth (alpha12)                                            |
# and ellipsoidal distance in metres (s) to a second point,                 |
#                                                                           |
# Calculate: the latitude and longitude of the second point (phi2, lembda2) |
# and the reverse azimuth (alpha21).                                        |
#                                                                           |
# ----------------------------------------------------------------------------

def vinc_pt(f, a, phi1, lembda1, alpha12, s):
    """

      Returns the lat and long of projected point and reverse azimuth
      given a reference point and a distance and azimuth to project.
      lats, longs and azimuths are passed in decimal degrees

      Returns ( phi2,  lambda2,  alpha21 ) as a tuple 

      """

    two_pi = 2.0 * math.pi

    if alpha12 < 0.0:
        alpha12 += two_pi
    if alpha12 > two_pi:
        alpha12 -= two_pi

    b = a * (1.0 - f)

    TanU1 = (1 - f) * math.tan(phi1)
    U1 = math.atan(TanU1)
    sigma1 = math.atan2(TanU1, math.cos(alpha12))
    Sinalpha = math.cos(U1) * math.sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha * Sinalpha

    u2 = cosalpha_sq * (a * a - b * b) / (b * b)
    A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 *
                                           (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    # Starting with the approximation
    sigma = (s / (b * A))

    last_sigma = 2.0 * sigma + 2.0  # something impossible

    # Iterate the following three equations 
    # until there is no significant change in sigma 

    # two_sigma_m , delta_sigma

    while abs((last_sigma - sigma) / sigma) > 1.0e-9:
        two_sigma_m = 2 * sigma1 + sigma

        delta_sigma = B * math.sin(sigma) * (
            math.cos(two_sigma_m) +
            (B / 4) * (math.cos(sigma) *
            (-1 + 2 * math.pow(math.cos(two_sigma_m), 2) -
            (B / 6) * math.cos(two_sigma_m) *
            (-3 + 4 * math.pow(math.sin(sigma), 2)) *
            (-3 + 4 * math.pow(math.cos(two_sigma_m), 2)))))
        last_sigma = sigma
        sigma = (s / (b * A)) + delta_sigma

    phi2 = math.atan2((math.sin(U1) * math.cos(sigma) +
                       math.cos(U1) * math.sin(
                       sigma) * math.cos(alpha12)),
                       ((1 - f) * math.sqrt(math.pow(Sinalpha, 2) +
                       pow(math.sin(U1) * math.sin(sigma) -
                           math.cos(U1) * math.cos(sigma) *
                           math.cos(alpha12), 2))))

    lembda = math.atan2((math.sin(sigma) * math.sin(alpha12)),
                        (math.cos(U1) * math.cos(sigma) -
                         math.sin(U1) * math.sin(sigma) * math.cos(alpha12)))

    C = (f / 16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq))

    omega = lembda - (1 - C) * f * Sinalpha * \
                     (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) +
                                                     C * math.cos(sigma) * (
                                                         -1 + 2 * math.pow(
                                                             math.cos(
                                                                 two_sigma_m),
                                                             2))))

    lembda2 = lembda1 + omega

    alpha21 = math.atan2(Sinalpha, (-math.sin(U1) * math.sin(sigma) +
                                    math.cos(U1) * math.cos(sigma) * math.cos(
                                        alpha12)))

    alpha21 += two_pi / 2.0
    if alpha21 < 0.0:
        alpha21 += two_pi
    if alpha21 > two_pi:
        alpha21 -= two_pi

    return phi2, lembda2, alpha21

    # END of Vincenty's Direct formulae

    ##---------------------------------------------------------------------------
    # Notes:
    #
    # * "The inverse formulae may give no solution over a line
    #       between two nearly antipodal points. This will occur when
    #       lembda ... is greater than pi in absolute value". (Vincenty, 1975)
    #
    # * In Vincenty (1975) L is used for the difference in longitude,
    #       however for consistency with other formulae in this Manual,
    #       omega is used here.
    #
    # * Variables specific to Vincenty's formulae are shown below,
    #       others common throughout the manual are shown in the Glossary.
    #
    #
    # alpha = Azimuth of the geodesic at the equator
    # U = Reduced latitude
    # lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2
    #               are the geodetic longitudes of points 1 & 2)
    # sigma = Angular distance on a sphere, from point 1 to point 2
    # sigma1 = Angular distance on a sphere, from the equator to point 1
    # sigma2 = Angular distance on a sphere, from the equator to point 2
    # sigma_m = Angular distance on a sphere, from the equator to the
    #               midpoint of the line from point 1 to point 2
    # u, A, B, C = Internal variables
    #
    #
    # Sample Data
    #
    # Flinders Peak
    # -37o57'03.72030"
    # 144o25'29.52440"
    # Buninyong
    # -37o39'10.15610"
    # 143o55'35.38390"
    # Ellipsoidal Distance
    # 54,972.271 m
    #
    # Forward Azimuth
    # 306o52'05.37"
    #
    # Reverse Azimuth
    # 127o10'25.07"
    #
    #
    ##*******************************************************************

    # Test driver

    # if __name__ == "__main__" :

    ## WGS84

    # a = 6378137.0
    # b = 6356752.3142
    # f = (a-b)/a

    # print("\n Ellipsoidal major axis =  %12.3f metres\n" % ( a ))
    # print("\n Inverse flattening     =  %15.9f\n" % ( 1.0/f ))

    # print("\n Test Flinders Peak to Buninyon")
    # print("\n ****************************** \n")
    # phi1 = -(( 3.7203 / 60. + 57) / 60. + 37 )
    # lembda1 = ( 29.5244 / 60. + 25) / 60. + 144
    # print("\n Flinders Peak = %12.6f, %13.6f \n" % ( phi1, lembda1 ))
    # deg = int(phi1)
    # minn = int(abs( ( phi1 - deg) * 60.0 ))
    # sec = abs(phi1 * 3600 - deg * 3600) - minn * 60
    # print(" Flinders Peak =   %3i\xF8%3i\' %6.3f\", % ( deg, minn, sec ),)
    # deg = int(lembda1)
    # minn = int(abs( ( lembda1 - deg) * 60.0 ))
    # sec = abs(lembda1 * 3600 - deg * 3600) - minn * 60
    # print(" %3i\xF8%3i\' %6.3f\" \n" % ( deg, minn, sec ))

    # phi2 = -(( 10.1561 / 60. + 39) / 60. + 37 )
    # lembda2 = ( 35.3839 / 60. + 55) / 60. + 143
    # print("\n Buninyon      = %12.6f, %13.6f \n" % ( phi2, lembda2 ))

    # deg = int(phi2)
    # minn = int(abs( ( phi2 - deg) * 60.0 ))
    # sec = abs(phi2 * 3600 - deg * 3600) - minn * 60
    # print(" Buninyon      =   %3i\xF8%3i\' %6.3f\",  " % ( deg, minn, sec ))
    # deg = int(lembda2)
    # minn = int(abs( ( lembda2 - deg) * 60.0 ))
    # sec = abs(lembda2 * 3600 - deg * 3600) - minn * 60
    # print(" %3i\xF8%3i\' %6.3f\" \n" % ( deg, minn, sec ))

    # dist, alpha12, alpha21   = vinc_dist  ( f, a, math.radians(phi1), math.radians(lembda1), math.radians(phi2),  math.radians(lembda2) )

    # alpha12 = math.degrees(alpha12)
    # alpha21 = math.degrees(alpha21)

    # print("\n Ellipsoidal Distance = %15.3f metres\n            should be         54972.271 m\n" % ( dist ))
    # print("\n Forward and back azimuths = %15.6f, %15.6f \n" % ( alpha12, alpha21 ))
    # deg = int(alpha12)
    # minn =int( abs(( alpha12 - deg) * 60.0 ) )
    # sec = abs(alpha12 * 3600 - deg * 3600) - minn * 60
    # print(" Forward azimuth = %3i\xF8%3i\' %6.3f\"\n" % ( deg, minn, sec ))
    # deg = int(alpha21)
    # minn =int(abs( ( alpha21 - deg) * 60.0 ))
    # sec = abs(alpha21 * 3600 - deg * 3600) - minn * 60
    # print(" Reverse azimuth = %3i\xF8%3i\' %6.3f\"\n" % ( deg, minn, sec ))


    ## Test the direct function */
    # phi1 = -(( 3.7203 / 60. + 57) / 60. + 37 )
    # lembda1 = ( 29.5244 / 60. + 25) / 60. + 144
    # dist = 54972.271
    # alpha12 = ( 5.37 / 60. + 52) / 60. + 306
    # phi2 = lembda2 = 0.0
    # alpha21 = 0.0

    ephi2, lembda2, alpha21 = vinc_pt(f, a,
                                      math.radians(phi1),
                                      math.radians(lembda1),
                                      math.radians(alpha12),
                                      dist)

    # phi2 = math.degrees(phi2)
    # lembda2 = math.degrees(lembda2)
    # alpha21 = math.degrees(alpha21)

    # print("\n Projected point =%11.6f, %13.6f \n" % ( phi2, lembda2 ))
    # deg = int(phi2)
    # minn =int(abs( ( phi2 - deg) * 60.0 ))
    # sec = abs( phi2 * 3600 - deg * 3600) - minn * 60
    # print(" Projected Point = %3i\xF8%3i\' %6.3f\", " % ( deg, minn, sec ),)
    # deg = int(lembda2)
    # minn =int(abs( ( lembda2 - deg) * 60.0 ))
    # sec = abs(lembda2 * 3600 - deg * 3600) - minn * 60
    # print("  %3i\xF8%3i\' %6.3f\"\n" % ( deg, minn, sec ))
    # print(" Should be Buninyon \n")
    # print("\n Reverse azimuth = %10.6f \n" % ( alpha21 ))
    # deg = int(alpha21)
    # minn =int(abs( ( alpha21 - deg) * 60.0 ))
    # sec = abs(alpha21 * 3600 - deg * 3600) - minn * 60
    # print(" Reverse azimuth = %3i\xF8%3i\' %6.3f\"\n\n" % ( deg, minn, sec ))

    ## lat/lon of New York
    # lat1 = 40.78
    # lon1 = -73.98
    ## lat/lon of London.
    # lat2 = 51.53
    # lon2 = 0.08
    # print('New York to London:')
    # gc = GreatCircle((2*a+b)/3.,(2*a+b)/3.,lon1,lat1,lon2,lat2)
    # print('geodesic distance using a sphere with WGS84 mean radius = ',gc.distance)
    # print('lon/lat for 10 equally spaced points along geodesic:')
    # lons,lats = gc.points(10)
    # for lon,lat in zip(lons,lats):
    # print(lon, lat)
    # gc = GreatCircle(a,b,lon1,lat1,lon2,lat2)
    # print('geodesic distance using WGS84 ellipsoid = ', gc.distance)
    # print('lon/lat for 10 equally spaced points along geodesic:')
    # lons,lats = gc.points(10)
    # for lon,lat in zip(lons,lats):
    #     print(lon,lat)


if __name__ == '__main__':
    import sys
    # Some tests for reading the mapping YAML database
    map_regions(map_region=sys.argv[1])
