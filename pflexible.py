#!/usr/bin/env python
"""
SYNOPSIS
========

    pflexible [-h] [-v,--verbose] [--version]

DESCRIPTION
===========

    pflexible: A python module for working with FLEXPART Output.


EXAMPLES
========

    #TODO:

    A lot! This is just a starting point. See the doc strings
    for information about the various functions.


AUTHOR
======

    JFB: John F Burkhart <jfburkhart@gmail.com>

CONTRIBUTORS
============

    HSO: Harald Sodemann
    SEC: Sabine Eckhardt
    AST: Andreas Stohl

    Many functions are adaptations of Fortran / NCARG programs (AST)
    or Matlab functions (HSO/SEC).

LICENSE
=======

    This script follows creative commons usage.

VERSION
=======

    ID: $Id$: $Rev$ 
"""
#builtin imports
import pdb
import sys
import os
import struct
import re
import traceback
import optparse
import time
import datetime
import itertools
from math import pi, sqrt, cos

#Dependencies:
# Numpy
import numpy as np
# Matplotlib
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import matplotlib.image as image
from matplotlib.patches import Ellipse
# Basemap
try:
    from mpl_toolkits.basemap import shiftgrid, addcyclic
except ImportError:
    from matplotlib.toolkits.basemap import shiftgrid, addcyclic


#local imports
import mapping as mp

__version__ = '0.8.8'
__path__ = os.path.abspath(os.curdir)

#### Functions for reading FLEXPART output ##### 


def read_command(path, headerrows=7):
    """
    Reads a FLEXPART COMMAND file.

    .. note::
        Assumes releases file has a header of 7 lines. Use
        option "headerrows" to override this.

    USAGE::

        > A = read_command(filepath)

    where filepath is either a file object or a path.


    Returns
      a dictionary with the following keys

      ============      =================================================================
      fields            description
      ============      =================================================================
      SIM_DIR           Simulation direction
      SIM_START         Text str of YYYYMMDD HHMMSS
      SIM_END           Text str of YYYYMMDD HHMMSS
      AVG_CNC_INT       Average concentrations are calculated every SSSSS seconds
      AVG_CNC_TAVG      The average concentrations are time averages of SSSSS sec
      CNC_SAMP_TIME     The concentrations are sampled every SSSSS seconds to
                        calculate the time average concentration.
      T_PARTSPLIT       Time constant for particle splitting.
      SYNC              Alll processes are synchronized with this time interval
      CTL
      IFINE             IFINE=Reduction factor for time step used for vertical wind
      IOUT              IOUT determines how the output shall be made: concentration
                        (ng/m3, Bq/m3), mixing ratio (pptv), or both, or plume
                        trajectory mode, or concentration + plume trajectory mode.
      IPOUT             IPOUT determines whether particle positions are outputted
                        (in addition to the gridded concentrations or mixing ratios)
                        or not. 0=no output, 1 output every output interval, 2 only
                        at end of the simulation
      LSUBGRID          Switch on/off subgridscale terrain parameterization
                        (increase of mixing heights due to subgridscale orog. var
      LCONVECTION       Switch on/off the convection parameterization
      LAGESPECTRA       Switch on/off the calculation of age spectra: if yes, the
                        file AGECLASSES must be available
      IPIN              If IPIN=1, a file "partposit_end" from a previous run must
                        be available in the output directory. Particle positions
                        are read in and previous simulation is continued. If
                        IPIN=0, no particles from a previous run are used
      OUTPUTFOREACHRELEASE Switch on/off writing out each release.
      IFLUX             If IFLUX is set to 1, fluxes of each species through each
                        of the output boxes are calculated. Six fluxes,
                        corresponding to northward, southward, eastward, westward,
                        upward and downward are calculated for each grid cell of
                        the output grid. The control surfaces are placed in the
                        middle of each output grid cell. If IFLUX is set to 0,
                        no fluxes are determined.
      MDOMAINFILL       If MDOMAINFILL is set to 1, the first box specified in file
                        RELEASES is used as the domain where domain-filling
                        trajectory calculations are to be done. Particles are
                        initialized uniformly distributed (according to the air mass
                        distribution) in that domain at the beginning of the
                        simulation, and are created at the boundaries throughout
                        the simulation period
      IND_SOURCE        IND_SOURCE switches between different units for
                        concentrations at the source. NOTE that in backward
                        simulations the release of computational particles takes
                        place at the "receptor" and the sampling of particles at
                        the "source".  1=mass units (for bwd-runs = concentration)
                                       2=mass mixing ratio units'''],

      IND_RECEPTOR      IND_RECEPTOR switches between different units for
                        concentrations at the receptor 1=mass units (concentrations)
                                                       2=mass mixing ratio units

      MQUASILAG         MQUASILAG indicates whether particles shall be numbered
                        consecutively (1) or with their release location number (0).
                        The first option allows tracking of individual particles
                        using the partposit output files

      NESTED_OUTPUT     NESTED_OUTPUT decides whether model output shall be made
                        also for a nested output field (normally with higher resolution)
      LINIT_COND        For Backward Runs, sets initial conditions:
                        [0]=No, 1=Mass Unit, 2=Mass Mixing

      ===============   =================================================================

    """
    lines = file(path, 'r').readlines()
    command_vals = [i.strip() for i in lines[headerrows:]] #clean line ends
    COMMAND_KEYS = (
               'SIM_DIR',
               'SIM_START',
               'SIM_END',
               'AVG_CNC_INT',
               'AVG_CNC_TAVG',
               'CNC_SAMP_TIME',
               'T_PARTSPLIT',
               'SYNC',
               'CTL',
               'IFINE',
               'IOUT',
               'IPOUT',
               'LSUBGRID',
               'LCONVECTION',
               'LAGESPECTRA',
               'IPIN',
               'OUTPUTFOREACHRELEASE',
               'IFLUX',
               'MDOMAINFILL',
               'IND_SOURCE',
               'IND_RECEPTOR',
               'MQUASILAG',
               'NESTED_OUTPUT',
               'LINIT_COND')
    float_keys=['CTL']
    date_keys=['SIM_START', 'SIM_END']
    
    COMMAND = {}
    
    for i, key in enumerate(COMMAND_KEYS):
        val = command_vals[i].split()
        if key in date_keys:
            val = val[:2]
        elif key in float_keys:
            val = float(val[0])
        else:
            val = int(val[0])
        
        COMMAND[key] = val
    return COMMAND




def read_releases(path, headerrows=11):
    """
    Reads a FLEXPART releases file.

    .. note::
        Assumes releases file has a header of 11 lines. Use
        option "headerrows" to override this.

    USAGE::

        > A = read_releases(filepath)

    where filepath is either a file object or a path.


    Returns
      a record array with fields:

      ============      ==========================
      fields            description
      ============      ==========================
      start_time        datetime start
      end_time          datetime end
      lllon             lower left longitude
      llat              lower left latitude
      urlon             upper right longitude
      urlat             upper right latitude
      altunit           1=magl, 2=masl, 3=hPa
      elv1              lower z level
      elv2              upper z level
      numpart           numparticles
      mass              mass for each spec, so the
                        array will actually have
                        fields: mass0, mass1,..
      id                ID for each release
      ============      ==========================

    """
    lines = _getfile_lines(path)
    lines = [i.strip() for i in lines] #clean line ends

    # we know nspec is at line 11
    nspec = int(lines[headerrows])
    blocklength = headerrows + nspec
    spec = []
    for i in range(nspec):
        spec.append(int(lines[headerrows + 1 + i]))
    indx = headerrows + 1 + (i + 2)
    blocks = []
    for i in range(indx, len(lines), blocklength + 1):
        blocks.append(lines[i:i + blocklength])
    for b in blocks:
        b[0] = datetime.datetime.strptime(b[0], '%Y%m%d %H%M%S')
        b[1] = datetime.datetime.strptime(b[1], '%Y%m%d %H%M%S')
        b[2:6] = [np.float(i) for i in b[2:6]]
        b[6] = int(b[6])
        b[7] = float(b[7])
        b[8] = float(b[8])
        b[9] = int(b[9])
        for i in range(nspec):
            b[10 + i] = float(b[10 + i])
        b = tuple(b)


    names = ['start_time', 'end_time', 'lllon', 'lllat', 'urlon', 'urlat', \
                          'altunit', 'elv1', 'elv2', 'numpart']
    #formats = [object, object, np.float, np.float, np.float, np.float,\
    #                      int, np.float, np.float, int]
    for i in range(nspec):
            names.append('mass%s' % i)
            #formats.append(np.float)
    names.append('id')
    #formats.append('S30')

    #dtype = {'names':names, 'formats':formats}
    #RELEASES = np.rec.array(blocks,dtype=dtype)
    return np.rec.fromrecords(blocks, names=names)


def read_trajectories(H, trajfile='trajectories.txt', \
                     ncluster=5, \
                     ageclasses=20):
    """
    Reads the trajectories.txt file in a FLEXPART run output directory.

    Based on output from plumetraj.f::

        ********************************************************************************
        *                                                                              *
        * Determines a plume centroid trajectory for each release site, and manages    *
        * clustering of particle locations. Certain parameters (average PV,            *
        * tropopause height, etc., are provided along the plume trajectories.          *
        * At the end, output is written to file 'trajectories.txt'.                    *
        *                                                                              *
        *     Author: A. Stohl                                                         *
        *                                                                              *
        *     24 January 2002                                                          *
        *                                                                              *
        * Variables:                                                                   *
        * fclust          fraction of particles belonging to each cluster              *
        * hmixcenter      mean mixing height for all particles                         *
        * ncluster        number of clusters to be used                                *
        * pvcenter        mean PV for all particles                                    *
        * pvfract         fraction of particles with PV<2pvu                           *
        * rms             total horizontal rms distance after clustering               *
        * rmsdist         total horizontal rms distance before clustering              *
        * rmsclust        horizontal rms distance for each individual cluster          *
        * topocenter      mean topography underlying all particles                     *
        * tropocenter     mean tropopause height at the positions of particles         *
        * tropofract      fraction of particles within the troposphere                 *
        * zrms            total vertical rms distance after clustering                 *
        * zrmsdist        total vertical rms distance before clustering                *
        * xclust,yclust,  Cluster centroid positions                                   *
        * zclust                                                                       *
        *                                                                              *
        ********************************************************************************

    USAGE::

        > T = read_trajectories(H_OR_path_to_directory, **kwargs)

    .. note::
        The first argument is either a path to a trajectory file, or simply a :class:`Header`
        instance.

    Returns a dictionary of the trajectories for each release.

      .. tabularcolumns::  |l|L|

      =============       ========================================
      Keys                Description
      =============       ========================================
      Trajectories        array_of_floats(j,it1,xi,yi,zi,topoi,hmixi,tropoi,pvi,
                          rmsdisti,rmsi,zrmsdisti,zrmsi,hfri,pvfri,trfri,
                          (xclusti(k),yclusti(k),zclusti(k),fclusti(k),rmsclusti(k),k=1,5))
      RELEASE_ID          1array(i1,i2,xp1,yp1,xp2,yp2,zp1,zp2,k,npart)
      numpspec            number of species
      numageclass         number of ageclasses
      =============       ========================================


    Arguments

      .. tabularcolumns::  |l|L|

      =============       ========================================
      keyword             Description [default]
      =============       ========================================
      trajfile            over the name of the input 
                          file ['trajectories.txt']
      ncluster            number of clusters [5]
      ageclasses          number of ageclasses [20]
      =============       ========================================

    """

    if isinstance(H, str):
        try:
            alltraj = file(H, 'r').readlines()
        except:
            raise IOError('Could not open file: %s' % H)
    else:
        path = H.path
        alltraj = file(os.path.join(path, trajfile), 'r').readlines()


    try:
        ibdate, ibtime, model, version = alltraj[0].strip().split()[:4]
    except:
        ibdate, ibtime = alltraj[0].strip().split()[:2]
        model = 'Flexpart'
        version = 'V.x'
    dt = datetime.datetime.strptime(ibdate + ibtime.zfill(6), '%Y%m%d%H%M%S')
    numpoint = int(alltraj[2].strip())
    # Fill a dictionary with the Release points and information keyed by name
    # RelTraj['RelTraj_ID'] = (i1,i2,xp1,yp1,xp2,yp2,zp1,zp2,k,npart)
    RelTraj = Structure()
    Trajectories = []

    for i in range(3, 3 + (numpoint * 2), 2):
        i1, i2, xp1, yp1, xp2, yp2, zp1, zp2, k, npart , = \
          tuple([float(j) for j in alltraj[i].strip().split()])
        itimerel1 = dt + datetime.timedelta(seconds=i1)
        itimerel2 = dt + datetime.timedelta(seconds=i2)
        Xp = (xp1 + xp2) / 2
        Yp = (yp1 + yp2) / 2
        Zp = (zp1 + zp2) / 2
        RelTraj[alltraj[i + 1].strip()] = np.array((itimerel1, itimerel2, Xp, Yp, Zp, k, npart))

    for i in range(3 + (numpoint * 2), len(alltraj)):
        raw = alltraj[i]
        format = [0, 5, 8, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6] + ncluster * [8, 8, 7, 6, 8]

        data = [raw[sum(format[:ii]):sum(format[:ii + 1])] for ii in range(1, len(format) - 1)] + \
                     [raw[sum(format[:-1]):]]
        ### FIX ###
        ### To get rid of '******' that is now in trajectories.txt
        data = [float(r.replace('********', 'NaN')) for r in data]

        Trajectories.append(data)

    data = np.array(Trajectories)
    RelTraj['version'] = model + ' ' + version
    RelTraj['date'] = dt
    RelTraj['Trajectories'] = data
    RelTraj['labels'] = \
            ['release number', 'seconds prior to release', 'lon', 'lat', 'height', 'mean topography', \
             'mean mixing height', 'mean tropopause height', 'mean PV index', \
             'rms distance', 'rms', 'zrms distance', 'zrms', \
             'fraction height??', 'fraction PV<2pvu', 'fraction in troposphere'] + \
            ncluster * ['xcluster', 'ycluster', 'zcluster', 'fcluster', 'rmscluster']
    RelTraj['info'] = \
    """
    Returns a dictionary:
        R['Trajectories'] = array_of_floats(
            releasenum,it1,xi,yi,zi,topoi,hmixi,tropoi,pvi,
            rmsdisti,rmsi,zrmsdisti,zrmsi,hfri,pvfri,trfri,
            (xclusti(k),yclusti(k),zclusti(k),fclusti(k),rmsclusti(k),k=1,5))

        R['RELEASE_ID'] = (dt_i1,dt_i2,xp1,yp1,xp2,yp2,zp1,zp2,k,npart)
        R['info'] = this message

    To plot a trajectory timeseries:
        RT = read_trajectories(H)
        T = RT['Trajectories']
        rel = 1
        t = T[np.where(T[:,0]==rel),:][0]
        plt.plot(t[:,1],t[:,14])
        plt.savefig('trajectories.png')
    """
    return RelTraj


def curtain_for_flightrack(H, flighttrack, nspec=0, npspec_int=0, index=0, get_track=False):
    """
    extracts curtain data along a given flight track

    input:  H (a pf.Header with a FD data object from H.read_grid)
            if the read_grid method has not been called, we'll call
            it below.

            flighttrack
            (a numpy array with datetime, lon, lat, elv values)
            
            nspec: is optional for reference to the FD dict.
            
            index: is optional in case numpointspec > 1
            
            get_track: if True extracts the points along the flight
            track rather than the curtain.

    output: curtain (a 2-d array with shape (len(flighttrack), len(H.outheight)
            
    TODO::
        add interpolation, extract a track, not curtain (e.g. z points)

    """

    #assert(isinstance(H, Header), 'ERROR: H Wrong Data Type')
    if H.direction == 'forward':
        if 'FD' in H.keys():
            pass
        else:
            # need to call H.read_grid to fill read the grids
            flighttimes = flighttrack[:, 0]
            time_ret = H.closest_dates(flighttimes, take_set=True)
            H.read_grid(nspec_ret=nspec, time_ret=time_ret,
                        npspec_int=npspec_int)

    elif H.direction == 'backward':
        if 'C' in H.keys():
            pass
        else:
            H.fill_backward()

    if get_track:
        curtain = np.zeros(len(flighttrack))
    else:
        curtain = np.zeros((len(flighttrack), len(H.outheight)))

    for i, (t, lon, lat, elv) in enumerate(flighttrack):

        
        #loop over FD tuples using t to extract the right time.
        # assuming t is a datetime
        idx = H.closest_date(t)
        datestring = H.available_dates[idx]
        if H.direction == 'forward':
            grid = H.FD[(nspec, datestring)]['grid']
            grid = grid[:, :, :, index]
            
        elif H.direction == 'backward':
            grid = H.C[(0, idx)]['grid']

        I = closest(lon, H.longitude)
        J = closest(lat, H.latitude)
        if get_track:
            K = closest(elv, H.outheight)
        #print lon, I, H.longitude[I], '===========' , lat, J, H.latitude[J]
            curtain[i] = grid[I, J, K]
        else:
            curtain[i] = grid[I, J, :]

    return curtain.T

def curtain_for_line(grid, X, Y, coords, index=0):
    """
    extracts curtain data from a grid given pairs of lon, lat

    input:  H (a pf.Header with a FD data object from H.read_grid)
            if the read_grid method has not been called, we'll call
            it below.

            flighttrack
            (a numpy array with datetime, lon, lat, elv values)
            
            nspec: is optional for reference to the FD dict.
            
            index: is optional in case numpointspec > 1
            
            get_track: if True extracts the points along the flight
            track rather than the curtain.

    output: curtain (a 2-d array with shape (len(flighttrack), len(H.outheight)
            
    TODO::
        add interpolation, extract a track, not curtain (e.g. z points)

    """
    if len(grid.shape) == 4:
        grid = grid[:, :, :, index]
    else:
        grid = grid
        
    curtain = np.zeros((len(coords), grid.shape[2]))
    for i, (x,y) in enumerate(coords):
       
        I = closest(x, X)
        J = closest(y, Y)
        curtain[i] = grid[I, J, :]
        curtain = np.nan_to_num(curtain)
        #curtain = np.ma.fix_invalid(np.array(curtain,dtype=float)) 
    return curtain.T


def groundlevel_for_line(H, X, Y, coords, index=0):
    """
    extracts ground level from H.heightnn along a track of lon, lat

    input:  H or H.heightnn (takes the lowest level)
    
            X, Y ,= H.longitude, H.latitude

            coords = zip(x, y) 
            
    output: groundlevel (a 1-d array with shape (len(flighttrack)
          

    """
    try:
        hgt = H.Heightnn[:,:,0]
    except:
        hgt = H
                
    # fix for hgt offset
    hgt = hgt - hgt.min()

    grndlvl = np.zeros(len(coords))
    
    for i, (x,y) in enumerate(coords):
       
        I = closest(x, X)
        J = closest(y, Y)
        grndlvl[i] = hgt[I, J]
        grndlvl = np.nan_to_num(grndlvl)
         
    return grndlvl - grndlvl.min()


def curtain_agltoasl(H, curtain_agl, coords, below_gl=0.0):
    """ converts the agl curtain to asl
    
        adds .asl_axis attribute to H
    """
    
    gl = groundlevel_for_line(H, H.longitude, H.latitude, coords)
    
    H.asl_axis = np.linspace(0, H.outheight[-1])
    
    xp = H.outheight - H.outheight[0]
    
    casl = np.zeros((len(H.asl_axis), len(coords)))
    
    for i in xrange(len(coords)): 
        casl[:,i] = np.interp(H.asl_axis, xp + gl[i], \
                              curtain_agl[:,i], left=below_gl)
        
    return casl


def read_agespectrum(filename, part=False, ndays=20):
    """
    Reads the spectrum.txt files generated from the "make_webpages" scripts.

    .. note::
        This functionality depends on having run the make_webpages scripts
        in-house. It is not supported in the public API.

    USAGE::

        > A = read_agespectrum(filename)

    Returns a dictionary containing ageclass information

      .. tabularcolumns::  |l|L|

      =============       ========================================
      Keys                Description
      =============       ========================================
      agespectrum         2d-array
                          (dt, lat, lon, alt, tracer[:numageclass]
                          for each release)
      ageclasses          1d-array(ageclass ages)
      numpspec            number of species
      numageclass         number of ageclassesA
      part                if it is from make webpages without the
                          header information.
      ndays               The number of ageclasses
      =============       ========================================

    """

    f = file(filename, 'r').readlines()

    line1 = f[0].strip().split()
    if part:
        numageclass = ndays
        ageclasses = np.arange(1, ndays * 1) * 60 * 60 * 24
        numspec = 1
    else:
        numageclass = int(line1[0])
        ageclasses = np.array([int(i) for i in line1[1:]])
        numspec = int(f[1].strip().split()[0])
    A = Structure()
    D = []
    for line in f[2:]:
        line = line.strip().split()
        # convert ibdate, ibtime to datetime
        dt = datetime.datetime.strptime(line[0] + line[1].zfill(6), '%Y%m%d%H%M%S')
        data = [float(i) for i in line[2:]]
        # reordering the array, dt object will be '0'
        data = [data[1]] + [data[0]] + data[2:]
        D.append([dt] + data)
    D = np.array(D)
    A['agespectrum'] = D
    A['numageclass'] = numageclass
    A['ageclasses'] = ageclasses
    A['numspec'] = numspec
    A['filename'] = filename
    A['info'] = \
    """
    A dictionary containing ageclass information:
    keys:
         agespectrum = 2d-array(dt, lon, lat, alt, tracer[:numageclass] for each release)
         ageclasses = 1d-array(ageclass ages)
         numpspec = number of species
         numageclass = number of ageclasses
         info = this message
    """
    return A


def save_spectrum(outf, H, agespectra, spectype='agespec',
                  header='## AGECLASS File'):
    """ Save an ageclass or continents spectrum to an outfile. """

    if spectype == 'agespec':
        ftype = 'AGECLASS'
        try:
            T = H.releasetimes
            spectrum = agespectra
            nClasses = spectrum.shape[1]
        except:
            #assume H is a list or None
            if H:
                nClasses = H[0]
                T = H[1]
                spectrum = agespectra
            else:
                T = agespectra[:, 0]
                nClasses = agespectra.shape[1] - 1
                spectrum = agespectra[:, 1:]
    elif spectype == 'contspec':
        #assume H is a list or None
        ftype = 'SPECTRUM'
        try:
            T = H.releasetimes
            spectrum = agespectra
            nClasses = spectrum.shape[1]
            header = '## Continental Spectrum File'
        except:
            #assume H is a list or None
            if H:
                nClasses = H[0]
                T = H[1]
                spectrum = agespectra
            else:
                T = agespectra[:, 0]
                nClasses = agespectra.shape[1] - 1
                spectrum = agespectra[:, 1:]

    T = np.array(T)
    T = np.reshape(T, (len(T), 1))
    xout = np.hstack((T, spectrum))
    fmt = '%s ' + nClasses * '%10.5f '
    outf.write("%s \n" % (header))
    np.savetxt(outf, xout, fmt=fmt)
    outf.close()




def gridarea(H):
    """returns an array of area corresponding to each nx,ny,nz

    Usage::

        > area = gridarea(H)


    Returns
        OUT = array area corresponding to nx,ny,nz

    Arguments
        H  = :class:`Header` object from readheader function.

    """


    pih = pi / 180.
    r_earth = 6.371e6
    cosfunc = lambda y : cos(y * pih) * r_earth
    nz = H['numzgrid']
    nx = H['numxgrid']
    ny = H['numygrid']
    outlat0 = H['outlat0']
    dyout = H['dyout']
    dxout = H['dxout']
    area = np.zeros((nx, ny))

    for iy in range(ny):
        ylata = outlat0 + (float(iy) + 0.5) * dyout #NEED TO Check this, iy since arrays are 0-index
        ylatp = ylata + 0.5 * dyout
        ylatm = ylata - 0.5 * dyout
        if (ylatm < 0 and ylatp > 0): hzone = dyout * r_earth * pih
        else:
            cosfact = cosfunc(ylata)
            cosfactp = cosfunc(ylatp)
            cosfactm = cosfunc(ylatm)
            if cosfactp < cosfactm:
                hzone = sqrt(r_earth ** 2 - cosfactp ** 2) - sqrt(r_earth ** 2 - cosfactm ** 2)
            else:
                hzone = sqrt(r_earth ** 2 - cosfactm ** 2) - sqrt(r_earth ** 2 - cosfactp ** 2)

        gridarea = 2.*pi * r_earth * hzone * dxout / 360.
        for ix in range(nx):
            area[ix, iy] = gridarea

    return area



def read_header(pathname, **kwargs):
    """
    The readheader function returns a special class (Structure) which behaves
    like a dictionary. It contains all the metadata from the simulation which
    is contained in the "header" or "header_nest" binary files from the model 
    output.
    
    .. warning::
        It is recommended to use the :class:`Header` class: H = pf.Header(path)

    This version is using the BinaryFile class rather than FortFlex. 

    Usage::

        > H = readheader(inputpath)

    Returns a dictionary

        H = dictionary like object with all the run metadata. TODO: Fill in keys.

    Arguments

      .. tabularcolumns::  |l|L|

      =============       ========================================
      keyword             Description [default]
      =============       ========================================
      pathname            FLEXPART run output directory
      readp               read release points 0=no, [1]=y
      readp_ff            readp_ff (read releases using Fortran [False]
      nest                nested output [0]=no, 1=yes
      version             version of FLEXPART, default = 'V8'
      =============       ========================================

    .. note::
        **This function is in development**

        This function is being developed so that there is no dependence on
        using f2py to compile the FortFlex module. So far it seems to work, but is
        notably slower than FortFlex. Please report any bugs found.

    .. todo::

        probably a lot of things.

    """

    OPS = Structure()
    OPS.readp = True
    OPS.readp_ff = False
    OPS.nest = False
    OPS.ltopo = 1 # 1 for AGL, 0 for ASL
    OPS.version = 'V8'
    OPS.headerfile = None
    OPS.update(kwargs)
#    for o in OPS:
#        print "Reading Header with:"
#        print "%s ==> %s" % (o, OPS[o])


    # Utility functions
    skip = lambda n = 8 : bf.seek(n, 1)
    getbin = lambda dtype, n = 1 : bf.read(dtype, (n,))

    #H={} #create dictionary for header
    h = Structure()

    if OPS.nest is True:
        filename = os.path.join(pathname, 'header_nest')
        h['nested'] = 1;
    else:
        filename = os.path.join(pathname, 'header')
        h['nested'] = 0;
    

    if OPS.headerfile:
        filename = os.path.join(pathname, OPS.headerfile);

    # Open header file in binary format
    bf = BinaryFile(filename, order="fortran")
    #Get available_dates from dates file in same directory as header
    datefile = os.path.join(pathname, 'dates')
    fd = file(datefile, 'r').readlines()
    #get rid of any duplicats (a fix for the forecast system)
    fd = sorted(list(set(fd)))
    h['available_dates'] = [d.strip('\n') for d in fd]

    #Define Header format and create Dictionary Keys
    I = {0:'rl0', 1:'ibdate', 2:'ibtime', 3:'version', \
         4:'rl1', 5:'loutstep', 6:'loutaver', 7:'loutsample', \
         8:'rl2', 9:'outlon0', 10:'outlat0', 11:'numxgrid', \
         12:'numygrid', 13:'dxout', 14:'dyout', 15:'rl3', 16:'numzgrid', \
         }
    #format for binary reading first part of the header file
    Dfmt = ['i', 'i', 'i', '13S', '2i', 'i', 'i', 'i', '2i', 'f', 'f', 'i', 'i', 'f', 'f', '2i', 'i']
    if bf:
        a = [bf.read(fmt) for fmt in Dfmt]
        for j in range(len(a)):
            h[I[j]] = a[j]
        #add items to the dictionary
        ss_start = datetime.datetime.strptime(str(h['ibdate']) + str(h['ibtime']).zfill(6), \
                                           '%Y%m%d%H%M%S')

        h['simulationstart'] = ss_start
        h['pathname'] = pathname
        h['decayconstant'] = 0
        h['outheight'] = np.array([getbin('f') for i in range(h['numzgrid'])])
        skip()
        h['jjjjmmdd'] = getbin('i')
        h['hhmmss'] = getbin('i')
        skip()
        h['nspec'] = getbin('i') / 3
        h['numpointspec'] = getbin('i')
        skip()
        #Read in the species names and levels for each nspec
        h['numzgrid'] = []
        h['species'] = []
        #temp dictionaries
        for i in range(h['nspec']):
            one = getbin('i'); # input skipped ??
            wd = getbin('c', 10);  # input skipped ??
            skip();
            one = getbin('i'); # input skipped ??
            dd = getbin('c', 10);  # input skipped ??
            skip();
            h['numzgrid'].append(getbin('i'))
            h['species'].append(''.join([getbin('c') for i in range(10)]).strip())
            skip();
        h['numpoint'] = getbin('i')

        # read release info if requested
        # if OPS.readp pass has changed, we cycle through once,
        # then break the loop if OPS.readp is false. This is done
        # in order to get some of the information into the header.
        before_readp = bf.tell()

        # initialise release fields
        I = {2:'kindz', 3:'xp1', 4:'yp1', 5:'xp2', \
           6:'yp2', 7:'zpoint1', 8:'zpoint2', 9:'npart', 10:'mpart'}
        h['ireleasestart'] = []
        h['ireleaseend'] = []
        for k, v in I.iteritems():
            h[v] = np.zeros(h['numpoint']) #create zero-filled lists in H dict
        h['compoint'] = []
        h['xmass'] = np.zeros((h['numpoint'], h['nspec']))
        r1 = getbin('i')
        for i in range(h['numpoint']):
            r2 = getbin('i')
            i1 = getbin('i')
            i2 = getbin('i')
            #h['ireleasestart'].append( ss_start + datetime.timedelta(seconds=float(i1)) )
            #h['ireleaseend'].append( ss_start + datetime.timedelta(seconds=float(i2)) )
            h['ireleasestart'].append(i1)
            h['ireleaseend'].append(i2)
            h['kindz'][i] = getbin('h') # This is an int16, might need to to change something
            skip() #get xp, yp,...

            h['xp1'][i] = getbin('f')
            h['yp1'][i] = getbin('f')
            h['xp2'][i] = getbin('f')
            h['yp2'][i] = getbin('f')
            h['zpoint1'][i] = getbin('f')
            h['zpoint2'][i] = getbin('f')

            skip() #get n/mpart
            h['npart'][i] = getbin('i')
            h['mpart'][i] = getbin('i')

            r3 = getbin('i')
            # initialise release fields

            l = getbin('i')#get compoint length?
            gt = bf.tell() + l #create 'goto' point
            sp = ''
            while re.search("\w", getbin('c')): #collect the characters for the compoint
                bf.seek(-1, 1)
                sp = sp + getbin('c')

            bf.seek(gt) #skip ahead to gt point

            h['compoint'].append(sp) #species names in dictionary for each nspec
            #h['compoint'].append(''.join([getbin('c') for i in range(45)]))
            r1 = getbin('i')
            #now loop for nspec to get xmass
            for v in range(h['nspec']):
                Dfmt = ['i', 'f', '2i', 'f', '2i', 'f', 'i']
                a = [bf.read(fmt) for fmt in Dfmt]
                h['xmass'][i, v] = a[1]

            if OPS.readp is False:
                """
                We get the first set of release points here in order
                to get some information, but skip reading the rest
                """
                bf.seek(before_readp)
                bf.seek(119 * h['numpoint'] + (h['nspec'] * 36) * h['numpoint'] + 4, 1)
                break


        rl = getbin('i')
        jnk_method = getbin('i')
        h['lsubgrid'] = getbin('i')
        h['lconvection'] = getbin('i')
        if rl == 20:
            h['ind_source'] = getbin('i')
            h['ind_receptor'] = getbin('i')
        skip()
        h['nageclass'] = getbin('i')
        Lage_fmt = ['i'] * h.nageclass
        jnk_lage = [bf.read(fmt) for fmt in Lage_fmt]
        nx = h['numxgrid']
        ny = h['numygrid']
        Dfmt = ['f'] * nx
        h['oro'] = np.zeros((nx, ny), np.float)
        skip()
        for ix in range(nx):
            #h['oro'][ix]=[getbin('f') for jx in range(ny)]
            # The next is *much* faster!
            h['oro'][ix] = getbin('f', ny)
            skip()
        if h['loutstep'] < 0:
            h['nspec'] = h['numpoint']

        bf.close

        # Calculate Height (outheight + topography)
        # There is an offset issue here related to the 0-indexing. Be careful.
        Z = h['oro'] #z is a numpy array
        nz = h['numzgrid'][0]
        Heightnn = np.zeros((nx, ny, nz), np.float)
        for ix in range(nx):
            if OPS.ltopo == 1:
                Heightnn[ix, :, 0] = Z[ix, :]
            else:
                Heightnn[ix, :, 0] = np.zeros(ny)

            for iz in range(nz):
                if OPS.ltopo == 1:
                    Heightnn[ix, :, iz] = [h['outheight'][iz] + Z[ix, y] for y in range(ny)]
                else:
                    Heightnn[ix, :, iz] = h['outheight'][iz]

        h['Area'] = gridarea(h)
        h['Heightnn'] = Heightnn


    #############  A FEW ADDITIONS ###########
    # add a few default attributes
    # optionally, use fortran routine to read
    # release points (deprecated)
    h.path = pathname
    if OPS.readp_ff:
        h = _read_headerFF(filename, h,
                          nxmax=h.numxgrid, nymax=h.numygrid, nzmax=h.numzgrid,
                          maxspec=h.nspec, maxageclass=h.nageclass,
                          maxpoint=h.numpoint)

    # Convert ireleasestart and ireleaseend to datetimes
    if OPS.readp:
        releasestart, releaseend = [], []
        for i in range(h.numpointspec):
            releasestart.append(h.simulationstart + \
                                 datetime.timedelta(seconds=int(h.ireleasestart[i])))
            releaseend.append(h.simulationstart + \
                               datetime.timedelta(seconds=int(h.ireleaseend[i])))
        h.releasestart = releasestart
        h.releaseend = releaseend[:h.numpointspec]
        h.releasetimes = [b - ((b - a) / 2) for a, b in zip(h.releasestart, h.releaseend)]

    # Add datetime objects for dates
    available_dates_dt = []
    for i in h.available_dates:
        available_dates_dt.append(datetime.datetime(
            int(i[:4]), int(i[4:6]), int(i[6:8]), int(i[8:10]), int(i[10:12]), int(i[12:])))
    h.available_dates_dt = available_dates_dt
    h.first_date = available_dates_dt[0]
    h.last_date = available_dates_dt[-1]
    h.ageclasses = np.array([act - h.simulationstart for act in h.available_dates_dt])
    h.numageclasses = len(h.ageclasses)

    # Add other helpful attributes
    h.nxmax = h.numxgrid
    h.nymax = h.numygrid
    h.nzmax = h.numzgrid
    h.maxspec = h.nspec
    h.maxpoint = h.numpoint
    h.area = h.Area

    if OPS.readp:
        h.xpoint = h.xp1
        h.ypoint = h.yp1

    # Add release unit derived from kindz
    if 'kindz' not in h.keys():
        h.kindz = [0]
        h.alt_unit = 'unkn.'
    if 3 in h.kindz:
        h.alt_unit = 'hPa'
    elif 2 in h.kindz:
        h.alt_unit = 'm.a.s.l.'
    elif 1 in h.kindz:
        h.alt_unit = 'm.a.g.l.'

    # 
    if h.loutstep > 0:
        h.direction = 'forward'
        h.unit = 'conc' #could be pptv
        h.plot_unit = 'ppb'
    else:
        h.direction = 'backward'
        h.unit = 'time'
        h.plot_unit = 'ns / kg' #Not sure about this

    # Units based on Table 1, ACP 2005
    if h.direction == 'forward':
        if h.ind_source == 1:
            if h.ind_receptor == 1:
                h.output_unit = 'ng m-3'
            if h.ind_receptor == 2:
                h.output_unit = 'pptm'
        if h.ind_source == 2:
            if h.ind_receptor == 1:
                h.output_unit = 'ng m-3'
            if h.ind_receptor == 2:
                h.output_unit = 'pptm'
    if h.direction == 'backward':
        if h.ind_source == 1:
            if h.ind_receptor == 1:
                h.output_unit = 's'
            if h.ind_receptor == 2:
                h.output_unit = 's m^3 kg-1'
        if h.ind_source == 2:
            if h.ind_receptor == 1:
                h.output_unit = 's kg m-3'
            if h.ind_receptor == 2:
                h.output_unit = 's'
    
                

    # Add layer thickness
    layerthickness = [h.outheight[0]]
    for i, lh in enumerate(h.outheight[1:]):
        layerthickness.append(lh - h.outheight[i])
    h.layerthickness = layerthickness

    h.options = OPS



    print 'Header read: %s' % filename

    return h

readheader = read_header

def readheaderV6(pathname, **kwargs):
    """
    see :func:`read_header`
    
    This is the 'V6' function. It is less tested than the 'V8' version.

    """
    warning = """
    This is in development, and many variables are not yet correctly registered.
    USE WITH CAUTION.
    """

    fail = -1
    OPS = Structure()
    OPS.readp = True
    OPS.readp_ff = False
    OPS.nest = False
    OPS.ltopo = 1 # 1 for AGL, 0 for ASL
    OPS.version = 'V6'
    OPS.update(kwargs)
#    for o in OPS:
#        print "Reading Header with:"
#        print "%s ==> %s" % (o, OPS[o])


    # Utility functions
    skip = lambda n = 8 : f2.seek(n, 1)
    getbin = lambda dtype, n = 1 : f2.read(dtype, (n,))

    #H={} #create dictionary for header
    h = Structure()

    if OPS.nest == False:
        filename = os.path.join(pathname, 'header'); h['nested'] = 0;
    else:
        filename = os.path.join(pathname, 'header_nest'); h['nested'] = 1;

    # Open header file in binary format
    f2 = BinaryFile(filename, order="fortran")
    #Get available_dates from dates file in same directory as header
    datefile = os.path.join(pathname, 'dates')
    fd = file(datefile, 'r').readlines()
    h['available_dates'] = [d.strip('\n') for d in fd]

    #Define Header format and create Dictionary Keys
    I = {0:'rl0', 1:'ibdate', 2:'ibtime', 3:'version', \
         4:'rl1', 5:'loutstep', 6:'loutaver', 7:'loutsample', \
         8:'rl2', 9:'outlon0', 10:'outlat0', 11:'numxgrid', \
         12:'numygrid', 13:'dxout', 14:'dyout', 15:'rl3', 16:'numzgrid', \
         }
    #format for binary reading first part of the header file
    Dfmt = ['i', 'i', 'i', '13S', '2i', 'i', 'i', 'i', '2i', 'f', 'f', 'i', 'i', 'f', 'f', '2i', 'i']
    if f2:
        a = [f2.read(fmt) for fmt in Dfmt]
        for j in range(len(a)):
            h[I[j]] = a[j]
        #add items to the dictionary
        ss_start = datetime.datetime.strptime(str(h['ibdate']) + str(h['ibtime']).zfill(6), \
                                           '%Y%m%d%H%M%S')

        h['simulationstart'] = ss_start
        h['pathname'] = pathname
        h['decayconstant'] = 0
        h['outheight'] = [getbin('f') for i in range(h['numzgrid'])]
        skip()
        h['jjjjmmdd'] = getbin('i')
        h['hhmmss'] = getbin('i')
        skip()
        h['nspec'] = getbin('i') / 3
        print 'nspec:', h['nspec']
        h['numpointspec'] = getbin('i')
        skip()
        #Read in the species names and levels for each nspec
        h['numzgrid'] = []
        h['species'] = []
        #temp dictionaries
        for i in range(h['nspec']):
            #pdb.set_trace()
            #one=getbin('i'); # input skipped ??
            wd = getbin('c', 10);  # input skipped ??
            skip();
            one = getbin('i'); # input skipped ??
            dd = getbin('c', 10);  # input skipped ??
            #print wd,dd
            skip();
            h['numzgrid'].append(getbin('i'))
            h['species'].append(''.join([getbin('c') for i in range(10)]).strip())
            #skip();
        skip()
        #print f2.read('i')

        h['numpoint'] = getbin('i')
        skip()

        # read release info if requested, this is deprecated. Use FortFlex routine.
        if OPS.readp is True:
            # initialise release fields
            I = {2:'kindz', 3:'xp1', 4:'yp1', 5:'xp2', \
               6:'yp2', 7:'zpoint1', 8:'zpoint2', 9:'npart', 10:'mpart'}
            h['ireleasestart'] = []
            h['ireleaseend'] = []
            for k, v in I.iteritems(): h[v] = np.zeros(h['numpoint']) #create zero-filled lists in H dict
            h['compoint'] = []
            h['xmass'] = np.zeros((h['numpoint'], h['nspec']))
            #r1=getbin('i')
            for i in range(h['numpoint']):
                #r2=getbin('i')
                i1 = getbin('i')
                i2 = getbin('i')
                #h['ireleasestart'].append( ss_start + datetime.timedelta(seconds=float(i1)) )
                #h['ireleaseend'].append( ss_start + datetime.timedelta(seconds=float(i2)) )
                h['ireleasestart'].append(i1)
                h['ireleaseend'].append(i2)
                h['kindz'][i] = getbin('i') # This is an int16, might need to to change something
                skip() #get xp, yp,...

                h['xp1'][i] = getbin('f')
                h['yp1'][i] = getbin('f')
                h['xp2'][i] = getbin('f')
                h['yp2'][i] = getbin('f')
                h['zpoint1'][i] = getbin('f')
                h['zpoint2'][i] = getbin('f')

                skip() #get n/mpart
                h['npart'][i] = getbin('i')
                h['mpart'][i] = getbin('i')

                r3 = getbin('i')
                # initialise release fields

                l = getbin('i')#get compoint length?
                gt = f2.tell() + l #create 'goto' point
                sp = ''
                sp = sp.join(getbin('c', l))
                f2.seek(gt) #skip ahead to gt point

                h['compoint'].append(sp) #species names in dictionary for each nspec
                #h['compoint'].append(''.join([getbin('c') for i in range(45)]))
                skip()
                #r1=getbin('i')

                #now loop for nspec to get xmass
                for v in range(h['nspec']):
                    Dfmt = ['i', 'f', '2i', 'f', '2i', 'f', 'i']
                    a = [f2.read(fmt) for fmt in Dfmt]
                    h['xmass'][i, v] = a[1]

        else:
            #skip reading points, fill data structure with zeroes
            # different for different species!!
            #f2.seek(119*h['numpoint']+(h['nspec']*36)*h['numpoint']+4,1);
            f2.seek(157 * h['numpoint'], 1);
        #rl = getbin('i')
        h['method'] = getbin('i')
        h['lsubgrid'] = getbin('i')
        h['lconvection'] = getbin('i')
        skip()
        #h['ind_source'] = getbin('i')
        #h['ind_receptor'] = getbin('i')
        h['nageclass'] = getbin('i')
        Lage_fmt = ['i'] * h.nageclass
        jnk_lage = [f2.read(fmt) for fmt in Lage_fmt]
        nx = h['numxgrid']
        ny = h['numygrid']
        Dfmt = ['f'] * nx
        h['oro'] = np.zeros((nx, ny), np.float)
        skip()
        for ix in range(nx):
            #h['oro'][ix]=[getbin('f') for jx in range(ny)]
            # The next is *much* faster!
            h['oro'][ix] = getbin('f', ny)
            skip()
        if h['loutstep'] < 0:
            h['nspec'] = h['numpoint']

        f2.close
        fail = 0

        # Calculate Height (outheight + topography)
        # There is an offset issue here related to the 0-indexing. Be careful.
        Z = h['oro'] #z is a numpy array
        nz = h['numzgrid'][0]
        Heightnn = np.zeros((nx, ny, nz), np.float)
        for ix in range(nx):
            if OPS.ltopo == 1: Heightnn[ix, :, 0] = Z[ix, :]
            else: Heightnn[ix, :, 0] = np.zeros(ny)

            for iz in range(nz):
                if OPS.ltopo == 1: Heightnn[ix, :, iz] = [h['outheight'][iz] + Z[ix, y] for y in range(ny)]
                else: Heightnn[ix, :, iz] = h['outheight'][iz]

        h['Area'] = gridarea(h)
        h['Heightnn'] = Heightnn


    #############  A FEW ADDITIONS ###########
    # add a few default attributes

    h.path = pathname
#    if OPS.readp_ff:
#        h = _read_headerFF(filename,h,
#                          nxmax=h.numxgrid,nymax=h.numygrid,nzmax=h.numzgrid,
#                          maxspec=h.nspec,maxageclass=h.nageclass,
#                          maxpoint=h.numpoint)
#
    # Convert ireleasestart and ireleaseend to datetimes

    releasestart, releaseend = [], []
    h.numpointspec = h.numpoint
    for i in range(h.numpointspec):
        releasestart.append(h.simulationstart + \
                             datetime.timedelta(seconds=int(h.ireleasestart[i])))
        releaseend.append(h.simulationstart + \
                           datetime.timedelta(seconds=int(h.ireleaseend[i])))
    h.releasestart = releasestart
    h.releaseend = releaseend[:h.numpointspec]
    h.releasetimes = [b - ((b - a) / 2) for a, b in zip(h.releasestart, h.releaseend)]
    # Add datetime objects for dates
    available_dates_dt = []
    for i in h.available_dates:
        available_dates_dt.append(datetime.datetime(
            int(i[:4]), int(i[4:6]), int(i[6:8]), int(i[8:10]), int(i[10:12]), int(i[12:])))
    h.available_dates_dt = available_dates_dt
    h.ageclasses = np.array([act - h.simulationstart for act in h.available_dates_dt])
    h.numageclasses = len(h.ageclasses)

    # Add other helpful attributes
    h.xpoint = h.xp1
    h.ypoint = h.yp1
    h.nxmax = h.numxgrid
    h.nymax = h.numygrid
    h.nzmax = h.numzgrid
    h.maxspec = h.nspec
    h.maxpoint = h.numpoint
    h.area = h.Area
    if h.loutstep > 0:
        h.direction = 'forward'
        h.unit = 'conc' #could be pptv
        h.plot_unit = 'ppb'
    else:
        h.direction = 'backward'
        h.unit = 'time'
        h.plot_unit = 'ns / kg' #Not sure about this

    # Add release unit derived from kindz
    if 3 in h.kindz:
        h.alt_unit = 'hPa'
    elif 2 in h.kindz:
        h.alt_unit = 'masl'
    else:
        h.alt_unit = 'magl'

    if h.loutstep > 0:
        h.direction = 'forward'
        h.unit = 'conc' #could be pptv
        h.plot_unit = 'ppb'
    else:
        h.direction = 'backward'
        h.unit = 'time'
        h.plot_unit = 'ns / kg' #Not sure about this


    print 'Header read: %s,\n%s' % (filename, warning)

    return h

readheaderV8 = read_header


def _readgrid_noFF(H, **kwargs):
    """ accepts a header dictionary as input, returns dictionary of Grid values
    %===========================================
    %
    %-------------------------------------------
    % input
    %   - H: required header dictionary
    %
    % optional (most is obtained from header dict)
    %   - date: which yyyymmddhhmmss from dates
    %   - unit: 'conc', 'pptv', ['time'], 'footprint'
    %   - nspec_ret: numspecies 
    %   - pspec_ret: 
    %   - age_ret:
    %   - time_ret: 
    %   - nest: nested = 1
    %
    % 
    % output
    %   - grid dictionary object
    %   - fail indicator (-1=fail, 0=success)
    %-------------------------------------------
    % FLEXPART python import routines
    %-------------------------------------------
    % last changes: JFB, 10.10.2008
    %===========================================
    """
    #if os.sys.platform != 'win32':
    #    try:
    #        from FortFlex import readgrid, sumgrid
    #        useFortFlex = 1
    #        print 'using FortFlex'
    #    except:
    #        useFortFlex = 0
    #        print 'Cannot find FortFlex.so'
    useFortFlex = 0

    if 'date' in kwargs.keys():
        date = kwargs['date']
    #else: date = H['ibtime']        
    else: date = None

    if 'unit' in kwargs.keys():
        unitname = kwargs['unit']
    else: unitname = 'time'

    units = ['conc', 'pptv', 'time', 'footprint']
    unit = units.index(unitname)

    if 'nspec_ret' in kwargs.keys():
        nspec_ret = kwargs['nspec_ret']
    else: nspec_ret = 1

    if 'pspec_ret' in kwargs.keys():
        pspec_ret = kwargs['ppsec_ret']
    else: pspec_ret = 1

    if 'age_ret' in kwargs.keys():
        age_ret = kwargs['age_ret']
    else: age_ret = 1

    if 'nest' in kwargs.keys():
        nest = kwargs['nest']
    else: nest = 0

    if 'time_ret' in kwargs.keys():
        time_ret = kwargs['time_ret']
    else:
        time_ret = range(len(H['available_dates']))

    if 'scaledepo' in kwargs.keys():
        scaledepo = kwargs['scaledepo']
    else:
        scaledepo = 1.0

    if 'scaleconc' in kwargs.keys():
        scaleconc = kwargs['scaleconc']
    else:
        scaleconc = 1.e9

    if 'decaycons' in kwargs.keys():
        decaycons = kwargs['decaycons']
    else:
        decaycons = 99999999990
    fail = -1

    #set filenames
    prefix = ['grid_conc_', 'grid_pptv_', \
              'grid_time_', 'footprint_total', \
              'grid_conc_nest_', 'grid_pptv_nest_', \
              'grid_time_nest_', 'footprint_total_nest']

    # local functions for reading binary data and creating grid dump
#    skip = lambda n=8 : f.read(n)
#    getbin = lambda fmt,n=1 : struct.unpack(fmt*n,f.read(struct.calcsize(fmt)))[0]

    # Utility functions
    skip = lambda n = 8 : f2.seek(n, 1)
    getbin = lambda dtype, n = 1 : f2.read(dtype, (n,))


    def getdump(n, fmt='f'):
        """ function to get the dump values for the sparse format """
        #print n
        skip()
        #Dfmt=[fmt]*n
#        a=[struct.unpack(ft,f.read(struct.calcsize(ft))) for ft in Dfmt]
        a = f2.read(fmt, n)
#        dumplist=[a[j][0] for j in range(len(a))]
        #dumplist=[a[j] for j in range(len(a))]
        return a #dumplist

    def key2var(D, key):
        cmd = "global %s; %s = D['%s'];" % (key, key, key)
        exec(cmd)

    def dumpgrid(dmp_i, cnt_r, dmp_r, grd, k, nage):
        """ function to dump sparse elements into grid """
        ii = 0
        fact = 1
        for ir in range(cnt_r):
            if dmp_r[ir] * fact > 0:
                n = dmp_i[ii]; ii = ii + 1; fact = fact * -1.
            else:
                n = n + 1

            kz = n / (numxgrid * numygrid)
            jy = (n - kz * numxgrid * numygrid) / numxgrid
            ix = n - numxgrid * numygrid * kz - numxgrid * jy
            #print "n  ==> ix,jy,kz,k,nage"
            #print "%s ==> %s,%s,%s,%s,%s" % (n,ix,jy,kz,k,nage)
            #print grd.shape
            #print grd[0,0,0,0,0]
            grd[ix, jy, kz - 1, k, nage] = abs(dmp_r[ir])
        return grd #flipud(grd.transpose())

    G = {}
    #get values from header file
    H['nageclass'] = H['numageclasses']
    headervars = ['nested', 'nspec', 'numxgrid', 'numygrid', 'numzgrid', 'nageclass', \
                  'available_dates', 'pathname', 'decayconstant', 'numpoint', 'Area', 'Heightnn']

    for k in headervars:
        key2var(H, k)

    #create zero arrays for datagrid and zplot
    datagrid = np.zeros((numxgrid, numygrid, numzgrid[0], numpoint), dtype=np.float32)
    zplot = np.empty((numxgrid, numygrid, numzgrid[0], numpoint))

    #--------------------------------------------------
    # Loop over all times, given in field H['dates']
    #--------------------------------------------------
    for ks in range(nspec_ret, nspec_ret + 1):
            #print 'processing: ' + str(H['dates'][date_i]) + '   ' + str(ks)
            # print 'processing: ' + str(date) + '   ' + str(ks)
        specs = '_' + str(ks).zfill(3)
        fpmax = -999.
        if date == None and time_ret == None:
            get_dates = [available_dates[0]]
        elif time_ret == None:
            get_dates = []
            date = date.strip().split(',')
            for d in date:
                try:
                    get_dates.append(available_dates.index(d))
                except:
                    shout("Cannot find date: %s in H['dates']\n" % d)
        else:
            get_dates = available_dates

        print 'getting grid for: ', get_dates

        for date_i in range(len(get_dates)):

            if unit != 4:
                datestring = H['available_dates'][date_i]
                #datestring = date
                filename = os.path.join(H['pathname'], \
                            prefix[(unit) + (nest * 4)] + datestring + specs)
            else:
                filename = os.path.join(H['pathname'], \
                            prefix[(unit) + (nest * 4)])
            
            #print 'reading: ' + filename

            if os.path.exists(filename):
                #print nspec,numxgrid,numygrid,numzgrid,nageclass,scaledepo,scaleconc,decaycons
                #print ks
                #print date

                if useFortFlex == 1:
                    #print 'Using FortFLEX'
                    #####################################################################################
                    ## FORTRAN WRAPPER CODE - only works on linux
                    ##
                    ## USAGE: grid = readgrid(filegrid,numxgrid,numygrid,numzgrid,\
                    ##                        nspec,nageclass,scaleconc,decayconstant)
                    ##
                    ##        zplot = sumgrid(zplot,grid,\
                    ##                        area,heightnn,\
                    ##                        [numxgrid,numygrid,numzgrid,numpoint,nageclass])
                    ##
                    ##
                    ## NOTE: numpoint = number of releases, ageclass always 1 for backward
                    ##
                    ## RETURN: grid(numxgrid,numygrid,numzgrid,numpoint,nageclass)
                    #####################################################################################

                    concgrid = readgrid(filename, numxgrid, numygrid, numzgrid[0], \
                                    numpoint, nageclass, \
                                    scaleconc, decayconstant)

                    #contribution[:,:,:] = concgrid[:,:,:,:,0]
                    print np.min(concgrid)
                    print np.max(concgrid)

                    #altitude = 50000
                    zplot = sumgrid(zplot, concgrid, \
                                    H.area, H.Heightnn)


                else:
                    dat_cnt = 0
                    nage = 0
                    #read data:
                    ##datagrid=np.zeros((numxgrid,numygrid,numzgrid[nspec-1],nspec,nageclass),np.float)
                    datagrid = np.zeros((numxgrid, numygrid, numzgrid[0], 1, 1), np.float)
                    #f = file(filename, 'rb')
                    #print filename
                    f2 = BinaryFile(filename, order='fortran')
                    skip(4)
                    G['itime'] = getbin('i');
                    print H['available_dates'][date_i]

                    #Read Wet Depostion
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
                    #wet=dumpgrid(dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage)
                    #Read Dry Deposition
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
                    #dry=dumpgrid(dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage)

                    #Read Concentrations
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
            #        print dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage
                    concgrid = dumpgrid(dmp_i, cnt_r, dmp_r, datagrid, ks - 1, nage)

                    #G[H['ibtime']].append(concgrid)
                    G[H['ibtime']] = concgrid
                    f2.close()
                fail = 0
            else:
                print "\n\n INPUT ERROR: Could not find file: %s" % filename
                raise IOError('No file: %s' % filename)

        if useFortFlex == 1: G = zplot
    return G

readgridBF = _readgrid_noFF

def _readgridBF(H, filename):
    """ Read grid using BinaryFile class"""
    # Utility functions
    skip = lambda n = 8 : f2.seek(n, 1)
    getbin = lambda dtype, n = 1 : f2.read(dtype, (n,))


    def getdump(n, fmt='f'):
        """ function to get the dump values for the sparse format """
        skip()
        #Dfmt=[fmt]*n
#        a=[struct.unpack(ft,f.read(struct.calcsize(ft))) for ft in Dfmt]
        a = f2.read(fmt, n)
#        dumplist=[a[j][0] for j in range(len(a))]
        #dumplist=[a[j] for j in range(len(a))]
        return a #dumplist

    def key2var(D, key):
        cmd = "global %s; %s = D['%s'];" % (key, key, key)
        exec(cmd)

    def _dumpgrid(dmp_i, cnt_r, dmp_r, grd, k, nage, nx, ny):
        """ function to dump sparse elements into grid, fall back method if 
        pflexcy.so module (cython) not available -- it is SLOW """
        conc = False
        if len(grd.shape) == 5:
            conc = True
        ii = 0
        fact = 1
        pos = 0
        for ir in range(cnt_r):

            if conc:
                if dmp_r[ir] * fact > 0:
                    n = dmp_i[ii]
                    ii = ii + 1
                    fact = fact * -1.
                else:
                    n = n + 1

                kz = n / (H.numxgrid * H.numygrid)
                jy = (n - kz * H.numxgrid * H.numygrid) / H.numxgrid
                ix = n - H.numxgrid * H.numygrid * kz - H.numxgrid * jy
                grd[ix, jy, kz - 1, k, nage] = abs(dmp_r[ir])

#
#                print "n  ==> ix,jy,kz,k,nage"
#                print "%s ==> %s,%s,%s,%s,%s" % (n,ix,jy,kz,k,nage)
#                print grd.shape
#                print grd[0,0,0,0,0]

            else:
                if dmp_r[ir] * fact > 0:
                    n = dmp_i[ii]
                    ii = ii + 1
                    fact = fact * -1.
                else:
                    pos = pos + 1
                jy = n / H.numxgrid
                ix = n - H.numxgrid * jy
                grd[ix, jy, k, nage] = abs(dmp_r[ir])

        return grd #flipud(grd.transpose())

    ## Import pflexcy.so (cython compiled version of dumpgrid)
    try:
        from pflexcy import dumpdatagrid, dumpdepogrid
        #print 'using pflexcy'
    except:
        print """WARNING: Using PURE Python to readgrid, execution will be slow.
         Try compiling the FortFlex module or the pflexcy module
         for your machine. For more information see the
         pflexible/f2py_build directory or use cython with pflexcy.pyx
        """
        dumpdatagrid = _dumpgrid
        dumpdepogrid = _dumpgrid

    dat_cnt = 0
    nage = 1
    ##datagrid=np.zeros((numxgrid,numygrid,numzgrid[nspec-1],nspec,nageclass),np.float)
    wetgrid = np.zeros((H.numxgrid, H.numygrid, H.numpointspec, 1), np.float)
    drygrid = np.zeros((H.numxgrid, H.numygrid, H.numpointspec, 1), np.float)
    datagrid = np.zeros((H.numxgrid, H.numygrid, H.numzgrid[0], H.numpointspec, nage), np.float)
    #f = file(filename,'rb')
    #print filename
    f2 = BinaryFile(filename, order='fortran')
    #read data:
    skip(4)
    itime = getbin('i')

    for na in range(nage):

        for ks in range(H.numpointspec):
            #Read Wet Depostion
            skip()
            cnt_i = getbin('i')
            dmp_i = getdump(cnt_i, 'i')
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            if dmp_r:
                print dmp_r, dmp_i
                wetgrid = dumpdepogrid(dmp_i, cnt_r, dmp_r, wetgrid, ks, na, H.numxgrid, H.numygrid)

            #Read Dry Deposition
            skip()
            cnt_i = getbin('i')
            dmp_i = getdump(cnt_i, 'i')
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            if dmp_r:
                print dmp_r, dmp_i
                drygrid = dumpdepogrid(dmp_i, cnt_r, dmp_r, drygrid, ks, na, H.numxgrid, H.numygrid)

            #Read Concentrations
            skip()
            cnt_i = getbin('i')
            dmp_i = getdump(cnt_i, 'i')
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            #print len(dmp_r),len(dmp_i)
            #print cnt_r,cnt_i
            #print dmp_i
            #print type(dmp_i),type(cnt_r),type(dmp_r),type(datagrid),type(ks),type(na)
            #print type(H.numxgrid),type(H.numygrid)
            datagrid = dumpdatagrid(dmp_i, cnt_r, dmp_r, datagrid, ks, na, H.numxgrid, H.numygrid)

    #G[H['ibtime']].append(concgrid)
    f2.close()

    return  datagrid, wetgrid, drygrid, itime



def _read_headerFF(pathname, h=None,
               maxpoint=800000, maxspec=4, maxageclass=20,
               nxmax=722, nymax=362, nzmax=40, verbose=True):
    """Called from read_header if readp_ff is True, uses FortFlex.readheader 
    to read all releasepoints.

    This function is dependant on the FortFlex.so module
    see FortFlex.f and the f2py directory
    """
    try:
        from FortFlex import readheader
    except:
        print "Error with FortFlex.readheader, use read_header"

    headervars = ['numxgrid', 'numygrid', 'numzgrid', 'outlon0', 'outlat0', 'compoint', \
                  'dxout', 'dyout', 'outheight', 'ibdate', 'ibtime', 'loutstep', \
                  'nspec', 'nageclass', 'lage', 'ireleasestart', 'ireleaseend', \
                  'numpoint', 'xpoint', 'ypoint', 'zpoint1', 'zpoint2', 'heightnn', 'area', \
                  'maxpoint', 'maxspec', 'maxageclass', 'nxmax', 'nymax', 'nzmax', \
                  'npart', 'kind', 'lage', 'loutaver', 'loutsample', 'yyyymmdd', 'hhmmss', 'method']
    if h == None:
        h = Structure()

    if verbose:
        print """Reading Header with:
                    maxpoint : %s
                    maxspec : %s
                    maxageclass : %s
                    nxmax : %s
                    nymax : %s
                    nzmax : %s
                    """ % (maxpoint, maxspec, maxageclass, nxmax, nymax, nzmax)


    numxgrid, numygrid, numzgrid, outlon0, outlat0, dxout, dyout, outheight, \
            ibdate, ibtime, loutstep, nspec, nageclass, lage, ireleasestart, ireleaseend, \
            numpoint, xpoint, ypoint, zpoint1, zpoint2, heightnn, area, compoint, species_f, \
            npart, kind, loutaver, loutsample, yyyymmdd, hhmmss, method = \
            readheader(pathname, maxpoint, maxspec, maxageclass, nxmax, nymax, nzmax)

    for v in headervars:
        if v not in h.keys():
            exec("h.%s = %s" % (v, v))
    return h

def read_grid(H, **kwargs):
    """
    Accepts a header object as input, returns dictionary of Grid values
    keyed by species and datestring from H['available_dates']. A grid instance from this
    dictionary may be passed to the get_slabs function to return a dictionary of slabs 
    suitable for plotting with the plotting routines. See below for more information.

    **DEPENDENCY**
        Requires FortFlex.so module compiled using f2py. See FortFlex.f for more details.

    Usage::

        > FLEXDATA = read_grid(H,**kwargs)

    Returns:

        A grid dictionary key by tuples with (species,available_dates), where species is
        and integer. See grid.keys()

        FLEXDATA[(s,datestring)]['grid']
        FLEXDATA[(s,datestring)]['itime']
        FLEXDATA[(s,datestring)]['shape']
        FLEXDATA[(s,datestring)]['max']
        FLEXDATA[(s,datestring)]['min']
        FLEXDATA[(s,datestring)]['timestamp']
        FLEXDATA[(s,datestring)]['species']
        FLEXDATA[(s,datestring)]['gridfile']
        FLEXDATA[(s,datestring)]['rel_i']
        FLEXDATA[(s,datestring)]['spec_i']

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      date                  which yyyymmddhhmmss from available_dates
                            or use (time_ret)
      time_ret              index to time
      unit                  'conc', 'pptv', ['time'], 'footprint'
      nspec_ret             numspecies
      pspec_ret             index to ???
      age_ret               index to ageclass
      nest                  obtained from H['nested']
      BinaryFile            Use BinaryFile vs. FortFlex [False]
      getwet                True, [False]
      getdry                True, [False]
      scaledepo             A float value to scale deposition [1.0]
      scaleconc             A float value to scale conc [1.0]
      decayconstant         A float for decay const. [9e6]
      calcfoot              Will cause footprint to be calculated
                            [False]
      verbose               more output
      version               A string 'V8' or 'V6', ['V8']
      =============         ========================================


    .. note::
        most arguments are able to be extracted fro the header "H"
    
    """
    if H.version == 'V8':
        return readgridV8(H, **kwargs)
    if H.version == 'V6':
        return readgridV6(H, **kwargs)
    else:
        raise IOError("No version attribute defined for Header.")


def readgridV8(H, **kwargs):

    """
    Accepts a header object as input, and selects appropriate readgrid function
    to use for reading in data from the flexpart binary Fortran files.

    See the :func:`read_grid` for more information on keyword arguments
    
    This is the 'V8' version of the function.

    """
    ## OPS is the options Structure, sets defaults, then update w/ kwargs
    OPS = Structure()
    OPS.unit = H.unit
    OPS.getwet = False
    OPS.getdry = False
    OPS.nspec_ret = 0
    OPS.npspec_int = False  # allows to select an index of npsec when calling readgrid
    OPS.pspec_ret = 0
    OPS.age_ret = 0
    OPS.time_ret = 0
    OPS.scaledepo = 1.0
    OPS.scaleconc = 1.0
    OPS.decayconstant = 9e6
    OPS.date = None
    OPS.calcfoot = False
    OPS.verbose = False
    OPS.BinaryFile = False
    OPS.version = 'V8'
    ## add keyword overides and options to header
    OPS.update(kwargs)
    #H.update(OPS)

    ## set up the return dictionary (FLEXDATA updates fd, fd is returned)
    FLEXDATA = {}
    fd = Structure()
    fd.options = Structure()

    ## What direction is the run?
    unit = OPS.unit
    
    if H['loutstep'] > 0:
        forward = True
        if unit == 'time':
            ## default forward unit
            unit = 'conc'
            OPS.unit = unit
    else:
        forward = False

    ## What species to return?
    nspec_ret = OPS.nspec_ret
    if isinstance(nspec_ret, int):
        nspec_ret = [nspec_ret]
    assert iter(nspec_ret), "nspec_ret must be iterable."

    ## get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, int) == True:
            time_ret = [time_ret]

        if time_ret[0] < 0:
            if forward == False:
                ## get all dates for calculating footprint.
                time_ret = np.arange(len(H.available_dates))
            else:
                raise ValueError("Must enter a positive time_ret for forward runs")

        for t in time_ret:
            get_dates.append(H.available_dates[t])


    ## define what dates to extract if user has explicitly defined a 'date'
    if OPS.date != None:
        date = OPS.date
        if time_ret is not None:
            Warning("overwriting time_ret variable, date was requested")
        get_dates = []
        if not isinstance(date, list):
            date = date.strip().split(',')
        for d in date:
            try:
                get_dates.append(H.available_dates[H.available_dates.index(d)])
                time_ret = None
            except:
                _shout("Cannot find date: %s in H['available_dates']\n" % d)

    if get_dates is None:
        raise ValueError("Must provide either time_ret or date value.")
    else:
        ## assign grid dates for indexing fd
        fd.grid_dates = get_dates[:]

    print 'getting grid for: ', get_dates
    # Some predifinitions
    fail = 0
    # set filename prefix
    prefix = ['grid_conc_', 'grid_pptv_', \
              'grid_time_', 'footprint_', 'footprint_total', \
              'grid_conc_nest_', 'grid_pptv_nest_', \
              'grid_time_nest_', 'footprint_nest_', 'footprint_total_nest'
              ]

    units = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
    unit_i = units.index(unit)

    # Determine what module to read, try to use FortFlex, then dumpgrid, lastly pure Python
    # import the FortFlex / Fortran module
    try:
        if OPS.version == 'V6':
            from FortFlex import readgrid_v6 as readgrid
            from FortFlex import sumgrid
            useFortFlex = True
            print 'using FortFlex VERSION 6'
        else:
            print('Assumed V8 Flexpart')
            from FortFlex import readgrid, sumgrid
            useFortFlex = True
    except:
        # get the original module (no memory allocation)
        try:
            from nilu.pflexpart.FortFlex import readgrid, sumgrid
            useFortFlex = True
            print 'using nilu.pflexpart FortFlex'
        except:
            useFortFlex = False
            print('Cannot load FortFlex, reverting to BinaryFile.')
    if not useFortFlex:
        readgrid = _readgridBF
        OPS.BinaryFile = True


   # reserve output fields
    print H.numxgrid, H.numygrid, H.numzgrid, OPS.nspec_ret, OPS.pspec_ret, OPS.age_ret, len(get_dates), H.numpoint

    # -------------------------------------------------

    ## add the requests to the fd object to be returned
    OPS.unit = unit
    fd.options.update(OPS)

    #--------------------------------------------------
    # Loop over all times, given in field H['available_dates']
    #--------------------------------------------------

    for date_i in range(len(get_dates)):
        datestring = get_dates[date_i]
        print datestring
        for s in nspec_ret: #range(OPS.nspec_ret,OPS.nspec_ret+1):A
            
            FLEXDATA[(s, datestring)] = Structure()
            spec_fid = '_' + str(s + 1).zfill(3)

            if unit_i != 4:
                filename = os.path.join(H['pathname'], \
                            prefix[(unit_i) + (H.nested * 5)] + datestring + spec_fid)
                H.zdims = H.numzgrid[0]

            else:
                #grid total footprint
                print "Total footprint"
                filename = os.path.join(H['pathname'], \
                            prefix[(unit_i) + (H.nested * 5)] + spec_fid)
                H.zdims = 1

            if os.path.exists(filename):
                H.filename = filename
                #print 'reading: ' + filename
                if OPS.verbose:
                    print 'with values:'
                    inputvars = ['filename', 'numxgrid', 'numygrid',
                                 'zdims', 'numpoint', 'nageclass', \
                                 'scaledepo', 'scaleconc',
                                 'decayconstant', 'numpointspec']
                    for v in inputvars:
                        print v, " ==> ", H[v]


                if OPS.BinaryFile:
                    gridT, wetgrid, drygrid, itime = _readgridBF(H, filename)
                else:
                    ## Quick fix for Sabine's Ship releases, added nspec_int so that only one
                    ## field of the nspec dimension is actually read
                    if OPS.npspec_int is not False:
                        npspec_int=OPS.npspec_int 
                        numpointspec = 1
                    else:
                        npspec_int=0
                        numpointspec = H.numpointspec
                    
                    gridT, wetgrid, drygrid, itime = readgrid(filename, \
                                                  H.numxgrid, H.numygrid,
                                                  H.zdims, numpointspec, H.nageclass, \
                                                  OPS.scaledepo, OPS.scaleconc, H.decayconstant, npspec_int)
                
                if OPS.getwet:
                    return wetgrid
                if OPS.getdry:
                    return drygrid
                if forward:
                    zplot = gridT[:, :, :, :, 0]
                else:
                    zplot = gridT[:, :, :, :, 0]

                if OPS.calcfoot:
                
                    zplot = sumgrid(zplot, gridT, \
                                    H.area, H.Heightnn)


                ## get the total column and prep the grid
                if H.direction == 'forward':
                    #not trying to do anything here... must be done
                    #after retrieving the grid
                    #D = get_slabs(H,np.squeeze(zplot))
                    rel_i = H.available_dates.index(datestring)
                    D = zplot
                    
                else:
                    D = zplot
                    rel_i = 'k'

                ## NOTE:
                ## If you're changing things here, you might want to change
                ## them in fill_backward as well, yes I know... something is
                ## poorly designed ;(
                
                FLEXDATA[(s, datestring)]['grid'] = D #zplot
                
                FLEXDATA[(s, datestring)]['itime'] = itime
                
                FLEXDATA[(s, datestring)]['shape'] = zplot.shape
                
                FLEXDATA[(s, datestring)]['max'] = zplot.max()
                
                FLEXDATA[(s, datestring)]['min'] = zplot.min()
                FLEXDATA[(s, datestring)]['timestamp'] = datetime.datetime.strptime(datestring, '%Y%m%d%H%M%S')
                FLEXDATA[(s, datestring)]['species'] = H['species'][s]
                FLEXDATA[(s, datestring)]['gridfile'] = filename
                FLEXDATA[(s, datestring)]['rel_i'] = rel_i
                FLEXDATA[(s, datestring)]['spec_i'] = s


            else:
                _shout('***ERROR: file %s not found! \n' % filename)
                fail = 1
        
        fd.set_with_dict(FLEXDATA)
        try:
            # just for testing, set the first available grid as a shortcut
            # this will be removed.
            qind = (nspec_ret[0], fd.grid_dates[0])
            fd.grid = fd[qind][fd[qind].keys()[0]].grid
        except:
            pass


    return fd

def readgridV6(H, **kwargs):

    """
    Accepts a header object as input, and selects appropriate readgrid function
    to use for reading in data from the flexpart binary Fortran files.

    See the :func:`read_grid` for more information on keyword arguments
    
    This is the 'V6' version of the function.

    """
    ## OPS is the options Structure, sets defaults, then update w/ kwargs
    OPS = Structure()
    OPS.unit = 'time'
    OPS.nspec_ret = 0
    OPS.pspec_ret = 0
    OPS.age_ret = 0
    OPS.time_ret = 0
    OPS.scaledepo = 1.0
    OPS.scaleconc = 1.0
    OPS.decayconstant = 9e6
    OPS.date = None
    OPS.calcfoot = False
    OPS.verbose = False
    OPS.BinaryFile = False
    ## add keyword overides and options to header
    OPS.update(kwargs)
    #H.update(OPS)

    ## set up the return dictionary (FLEXDATA updates fd, fd is returned)
    FLEXDATA = {}
    fd = Structure()
    fd.options = Structure()
    ## add the requests to the fd object to be returned
    fd.options.update(OPS)

    ## What direction is the run?
    unit = OPS.unit
    if H['loutstep'] > 0:
        forward = True
        if unit == 'time':
            ## default forward unit
            unit = 'conc'
    else:
        forward = False

    ## What species to return?
    nspec_ret = OPS.nspec_ret
    if isinstance(nspec_ret, int):
        nspec_ret = [nspec_ret]
    assert iter(nspec_ret), "nspec_ret must be iterable."

    ## get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, int) == True:
            time_ret = [time_ret]

        if time_ret[0] < 0:
            if forward == False:
                ## get all dates for calculating footprint.
                time_ret = np.arange(len(H.available_dates))
            else:
                raise ValueError("Must enter a positive time_ret for forward runs")

        for t in time_ret:
            get_dates.append(H.available_dates[t])


    ## define what dates to extract if user has explicitly defined a 'date'
    if OPS.date != None:
        date = OPS.date
        if time_ret is not None:
            Warning("overwriting time_ret variable, date was requested")
        get_dates = []
        if not isinstance(date, list):
            date = date.strip().split(',')
        for d in date:
            try:
                get_dates.append(H.available_dates[H.available_dates.index(d)])
                time_ret = None
            except:
                _shout("Cannot find date: %s in H['available_dates']\n" % d)

    if get_dates is None:
        raise ValueError("Must provide either time_ret or date value.")
    else:
        ## assign grid dates for indexing fd
        fd.grid_dates = get_dates[:]

    print 'getting grid for: ', get_dates
    # Some predifinitions
    fail = 0
    # set filename prefix
    prefix = ['grid_conc_', 'grid_pptv_', \
              'grid_time_', 'footprint_', 'footprint_total', \
              'grid_conc_nest_', 'grid_pptv_nest_', \
              'grid_time_nest_', 'footprint_nest_', 'footprint_total_nest'
              ]

    units = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
    unit_i = units.index(unit)
    # Determine what module to read, try to use FortFlex, then dumpgrid, lastly pure Python
    # import the FortFlex / Fortran module
    try:
        from FortFlex import readgrid_v6 as readgrid
        from FortFlex import sumgrid
        useFortFlex = True
        OPS.useFortFlex = useFortFlex
        print 'using FortFlex VERSION 6'
    except:
        # get the original module (no memory allocation)
        useFortFlex = False
        raise Warning('Cannot import FortFlex trying to use BinaryFile')
    if not useFortFlex:
        readgrid = _readgridBF
        OPS.BinaryFile = True

    #cheat: use key2var function to get values from header dict, H
    # need to change this to using H.xxxxx
    #headervars = ['nested','nspec','numxgrid','numygrid','numzgrid','nageclass',\
    #              'dates','pathname','decayconstant','numpoint','numpointspec',
    #              'area','Heightnn','lage']
    #for k in headervars:
    #    key2var(H,k)



   # reserve output fields
    print H.numxgrid, H.numygrid, H.numzgrid, OPS.nspec_ret, OPS.pspec_ret, OPS.age_ret, len(get_dates), H.numpoint

    # -------------------------------------------------

    #if forward:
    #    numpointspec = 1
    #else:
    #    numpointspec = numpoint
    #numpointspec = H['numpointspec']
    #if unit_i == 4:
        #grid=np.empty((numxgrid,numygrid,1,nspec_ret,pspec_ret,age_ret,len(get_dates)))
    #else:
        #grid=np.empty((numxgrid,numygrid,numzgrid[0],nspec_ret,pspec_ret,age_ret,len(get_dates)))
    #zplot=np.empty((numxgrid,numygrid,numzgrid[0],numpointspec))

    #--------------------------------------------------
    # Loop over all times, given in field H['available_dates']
    #--------------------------------------------------

    for date_i in range(len(get_dates)):
        datestring = get_dates[date_i]
        print datestring
        FLEXDATA[datestring] = {}
        for s in nspec_ret: #range(OPS.nspec_ret,OPS.nspec_ret+1):
            total_footprint = False
            FLEXDATA[(s, datestring)] = Structure()
            #spec_fid = '_'+str(s+1).zfill(3)

            if unit_i != 4:
                filename = os.path.join(H['pathname'], \
                            prefix[(unit_i) + (H.nested * 5)] + datestring)
                H.zdims = H.numzgrid[0]

            else:
                #grid total footprint
                print "Total footprint"
                total_footprint = True
                filename = os.path.join(H['pathname'], \
                            prefix[(unit_i) + (H.nested * 5)])
                H.zdims = 1

            if os.path.exists(filename):
                H.filename = filename
                #print 'reading: ' + filename
                if OPS.verbose:
                    print 'with values:'
                    inputvars = ['filename', 'numxgrid', 'numygrid',
                                 'zdims', 'numpoint', 'nageclass', \
                                 'scaledepo', 'scaleconc',
                                 'decayconstant', 'numpointspec']
                    for v in inputvars:
                        print v, " ==> ", H[v]


                if OPS.BinaryFile:
                    #print 'Using BinaryFile'
                    gridT, wetgrid, drygrid, itime = _readgridBF(H, filename)
                else:
                    gridT, wetgrid, drygrid, itime = readgrid(filename, \
                                                  H.numxgrid, H.numygrid,
                                                  H.zdims, H.numpointspec, H.nageclass, \
                                                  OPS.scaledepo, OPS.scaleconc, H.decayconstant)


                if forward:
                    zplot = gridT[:, :, :, :, 0]
                else:
                    zplot = gridT[:, :, :, :, 0]

                #if total_footprint:
                #    zplot = np.squeeze(gridT)

                if OPS.calcfoot:
                    zplot = sumgrid(zplot, gridT, \
                                    H.area, H.Heightnn)


                ## get the total column and prep the grid
                if H.direction == 'forward':
                    #not trying to do anything here... must be done
                    #after retrieving the grid
                    #D = get_slabs(H,np.squeeze(zplot))
                    #rel_i = H.available_dates.index(datestring)
                    D = zplot
                    rel_i = 'k'
                else:
                    D = zplot
                    rel_i = 'k'

                ## NOTE:
                ## If you're changing things here, you might want to change
                ## them in fill_backward as well, yes I know... something is
                ## poorly designed ;(
                FLEXDATA[(s, datestring)]['grid'] = D #zplot
                FLEXDATA[(s, datestring)]['itime'] = itime
                FLEXDATA[(s, datestring)]['shape'] = zplot.shape
                FLEXDATA[(s, datestring)]['max'] = zplot.max()
                FLEXDATA[(s, datestring)]['min'] = zplot.min()
                FLEXDATA[(s, datestring)]['timestamp'] = datetime.datetime.strptime(datestring, '%Y%m%d%H%M%S')
                FLEXDATA[(s, datestring)]['species'] = H['species'][s]
                FLEXDATA[(s, datestring)]['gridfile'] = filename
                FLEXDATA[(s, datestring)]['rel_i'] = rel_i
                FLEXDATA[(s, datestring)]['spec_i'] = s


            else:
                _shout('***ERROR: file %s not found! \n' % filename)
                fail = 1
        fd.set_with_dict(FLEXDATA)

        try:
            # just for testing, set the first available grid as a shortcut
            # this will be removed.
            qind = (0, fd.grid_dates[0])
            fd.grid = fd[qind][fd[qind].keys()[0]].grid
        except:
            pass


    return fd

def fill_backward(H, nspec=0, FD=None, add_attributes=False):
    """ for backward runs, calculates the 20-day sensitivity at each release point.

    Usage::

        > C = fill_backward(H,nspec=(0))


    This will cycle through all available_dates and create the filled backward array
    for each k in H.numpointspec.

    Returns

        A dictionary keyed by a (species,k) tuple.
        OR
        C & FD attributes on H

    Each element in the dictionary is a 3D array (x,y,z) for each species,k
    .. note::
        USE with *caution*, this is **MEMORY** intensive!

    Arguments

      .. tabularcolumns::  |l|L|

      ==============        ========================================
      keyword               Description [default]
      ==============        ========================================
      nspec                 the specied ID or a tuple of species IDs
      FD                    FD can be passed if it is already read
      add_attributes        will add C and FD as attributes to H,
                            rather than returning just C
      ==============        ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes, maybe there is a
        better way to handle this.

    """

    assert H.direction == 'backward', "fill_backward is only valid for backward runs"
    ## make sure npsec is iterable
    if isinstance(nspec, int):
        species = [nspec]
    else:
        species = nspec
    assert iter(species), 'nspec must be iterable, or you can pass an int'
    ## initialize variables
    C = Structure()
    for s, k in itertools.product(species, range(H.numpointspec)):
        C[(s, k)] = Structure()
        C[(s, k)].grid = np.zeros((H.numxgrid, H.numygrid, H.numzgrid[0]))
        C[(s, k)]['itime'] = None
        C[(s, k)]['timestamp'] = H.releasetimes[k]
        C[(s, k)]['species'] = H['species'][s]
        C[(s, k)]['gridfile'] = 'multiple'
        C[(s, k)]['rel_i'] = k
        C[(s, k)]['spec_i'] = s

    if FD is None:
        ## then we need to read the grids
        FD = read_grid(H, time_ret= -1, nspec_ret=species)

    ## read data grids and attribute/sum sensitivity
    print species
    for s in species:
        ## cycle through all the date grids (20days back)
        for d in FD.grid_dates:
            ## cycle through each release point
            for k in range(H.numpointspec):
                contribution = FD[(s, d)].grid[:, :, :, k]
                C[(s, k)].grid = C[(s, k)].grid + contribution

    for s, k in C:
        ## add total column
        C[(s, k)].slabs = get_slabs(H, C[(s, k)].grid)
        ## shape, min, max based on total column
        C[(s, k)]['shape'] = C[(s, k)].grid[0].shape
        C[(s, k)]['max'] = C[(s, k)].grid[0].max()
        C[(s, k)]['min'] = C[(s, k)].grid[0].min()


    if add_attributes:
        H.C = C
        H.FD = FD
    else:
        return C


def read_emissions(emissionsfile, E=None, maxemissions=1):
    """ Use reademissions.so module to read emissions file """
    try:
        from FortFlex import reademissions
    except:
        raise ImportError("Cannot find FortFlex module or missing reademissions")
    if not E:
        E = Structure()
        # set defaults for global 0.5 degree emissions
        defaults = {'nxmax':720, 'nymax':360, 'outlon0':-180, 'outlat0':-90, \
                    'numxgrid':720, 'numygrid':360, 'dxout':0.5, 'dyout':0.5}
        for k, v in defaults.iteritems():
            exec("E.%s = %s" % (k, v))

    emissions = reademissions(emissionsfile, maxemissions, E.nxmax, E.nymax, \
                       E.outlon0, E.outlat0, E.numxgrid, E.numygrid, E.dxout, E.dyout)
    E.grid = emissions
    E.filename = emissionsfile
    return E


def get_slabs(H, G, index=None, normAreaHeight=True, scale=1.0):
    """ Preps grid from readgridV8 for plotting.

    Accepts an 3D or 4D GRID from readgrid with optional index.
    shape ~= (360, 180, 3, 88, 1)

    Usage::

        plotgrid = get_slabs(H,G)

    Inputs

        H = A Header instance from a FLEXPART run
        G = A grid from the FLEXPARTDATA or FD dictionary returned from read_grid

    Returns

        slabs : a dictionary of rank-2 arrays corresponding to vertical levels.

        **Note**: transposes the grid for plotting.
        Total Column as D[0]
        Level 1 as D[1]
        Level 2 as D[2]
        ...

    Arguments

      .. tabularcolumns::  |l|L|

      ===============         ========================================
      keyword                 Description [default]
      ===============         ========================================
      index                   a release point index (k)
      normAreaHeight          [True], normalizes by the area/height
      ===============         ========================================

    .. todo::
        Need to check the normalization and indexing.

    """
    Heightnn = H['Heightnn']
    area = H['area']
    Slabs = Structure()
    grid_shape = G.shape
    if len(grid_shape) is 4:
        if index is None:
            if H.direction is 'forward':
                index = 0
            else:
                index = 0
            g = G[:, :, :, index]
        else:
            try:
                g = G[:, :, :, index]
            except:
                print '######### ERROR: Which Release Point to get? ########## '
                raise IOError
    if len(grid_shape) is 3:
        if index is None:
            g = G
        else:
            g = G


    for i in range(g.shape[2]):
        # first time sum to create Total Column
        if i == 0:
            TC = np.sum(g, axis=2).T
            TC = TC * scale
        if normAreaHeight:
            data = g[:, :, i] / Heightnn[:, :, i]
        else:
            data = g[:, :, i]

        Slabs[i + 1] = data.T * scale

    Slabs[0] = TC
    return Slabs


#### FLEXPART Plotting Functions  ########

def _plot_dropm():
    try:
        del FIGURE.m
        plt.close('all')
        del m
        "m dropped"
        return 1
    except:
        print 'Cannot delete m, does not exist.'
        return 0




def plot_trajchar(H, R, numrelease=1, varindex=[4, 15, 14], \
              FIGURE=None, map_region=None, projection=None, coords=None):
    """ plot FLEXPART trajectory (or set of trajectories) characteristics
        R is a dictionary returned from pflexible.read_trajectories

        varindex is a list of the variables to plot from the trajectories
    """
    if FIGURE == None:
        FIGURE = mp.get_FIGURE(getm=False)

    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    try:
        del fig.axes[0]
    except:
        print 'ax not removed'


    T = R['Trajectories']
    labels = R['labels']
    t = T[np.where(T[:, 0] == numrelease), :][0]
    numplot = len(varindex)

    for i in range(numplot):
        sp = fig.add_subplot(numplot, 1, i + 1)
        sp.plot(t[:, 1], t[:, varindex[i]])
        ax = plt.gca()
        if i + 1 != numplot: plt.setp(ax, xticklabels=[])
        plt.ylabel(labels[varindex[i]])
        plt.grid('on')

    return FIGURE


def plot_releases(R, FIGURE=None, threedim=False,
                    map_region=None, projection=None, coords=None,
                    overlay=True,
                    draw_circles=True,
                    draw_labels=True,
                    cbar2=True,
                    MapPar=None):
    """ plot the llon,llat,elv1 of the releases.

    Usage::
        >F = plot_releases(R)

    See the :func:`read_releases` function for information on how to generate "R"
    """

    if threedim:
        try:
            from mpl_toolkits.mplot3d import Axes3D
        except:
            raise ImportError('mplot3d not available from mpl_toolkits!')

        fig = plt.figure()
        ax = Axes3D(fig)

        ax.plot(R.lllon, R.lllat, R.elv1)

        return fig
    else:
        ## Set up the FIGURE
        if FIGURE == None:
            FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)
        ##Get fig info and make active
        fig = FIGURE.fig
        m = FIGURE.m
        ax = FIGURE.ax
        plt.figure(fig.number)
        plt.axes(ax)


        ## prepare the data
        lon = R.lllon
        lat = R.lllat
        zlevel = R.elv1
        zsize = np.ones(len(lon)) * 30
        marker = 'o'


        ## clear the previous track
        if 'circles' in FIGURE.keys():
            del FIGURE['circles']
        if overlay is False:
            del ax.collections[FIGURE.indices.collections:]
            del ax.texts[FIGURE.indices.texts:]

        ## plot the track
        cx, cy = m(lon, lat)
        if draw_circles:
            cmap = plt.get_cmap('gist_gray')
            circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                              marker=marker, edgecolor=None,
                              zorder=10, alpha=0.85)
            ## m.scatter has no color bar,
            ## so create a ghost 'scatter' instance:
            pos = ax.get_position()
            l, b, w, h = getattr(pos, 'bounds', pos)
            jnkfig = plt.figure()
            jnkax = jnkfig.add_axes([l, b, w, h], frameon=False)
            jnkmap = jnkax.scatter(cx, cy, zsize, zlevel,
                                   cmap=cmap, marker=marker,
                                   edgecolor=None)

            try:
                FIGURE.circles = circles
            except:
                pass
            ## make the figure active again
            plt.figure(fig.number);
            ## draw the legend and title
            ## check if a colorbar legend exists (this will cause
            ## problems for subplot routines!)
            if cbar2 is True:
                cax = plt.axes([l + w + 0.12, b, 0.02, h - 0.035])
            else:
                cax = plt.axes([l + w + 0.03, b, 0.025, h - 0.035])
            cb = fig.colorbar(jnkmap, cax=cax) # draw colorbar
            p_leg = mpl.font_manager.FontProperties(size='6')
            cax.set_title('altitude\n(m)', fontproperties=p_leg)

            ## delete the ghost instance
            plt.close(jnkfig.number)
            del jnkax, jnkfig, jnkmap

        plt.axes(ax)
        plt.setp(ax, xticks=[], yticks=[])

        if draw_labels:
            p_cax = mpl.font_manager.FontProperties(size='8',
                                                    style='italic',
                                                    weight='bold',
                                                   )

        FIGURE.fig = fig
        FIGURE.m = m
        FIGURE.ax = ax
        return FIGURE


def plot_spectra(inspectra,
                    plt_title='', fig_title='',
                    y_label=None,
                    spectra_label='bins',
                    FIGURE=None, y_datarange=None,
                    cum=False, labels=[],
                    bars=False, debug=False):
    """ plot an spectra

    Usage::

        > FIG = plot_agespectra(spectra,*kwargs)


    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      runid                 an ID for the title
      yunits                units for the y-axis [None]
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      data_range            y-axis data range
      web_version           Set to True if using the make_webpages
                            function. see: :func:`read_agespectrum`
      continental           Set to True if continental_spectrum
      cum                   Set to true if data should be cumulative
      bars                  Plot using stacked bars rather than a
                            filled line plot
      =============         ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes, maybe there is a
        better way to handle this.

    .. note::
        Required attributes of 'H' if you want to make a dummy H. Or just set
        H to "None" and it will be taken from the agespectra input.
        H.numageclasses
        H.releasetimes


    """
    ## make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    if FIGURE == None:
        FIGURE = Structure()
        fig = plt.figure(figsize=(8, 6))
        FIGURE.fig = fig
        ax = fig.add_subplot(111)#,pos=[0.1,0.2,.8,.7])
        FIGURE.ax = ax
    else:
        fig = FIGURE.fig
        ax = FIGURE.ax


    try:
        numageclasses = inspectra.numageclass
        releasetimes = inspectra.releasetimes
        inspectra = inspectra[:]
    except:
        numageclasses = inspectra.shape[1] - 1
        releasetimes = inspectra[:, 0]
        inspectra = inspectra[:, 1:]

    if cum == 'norm':
        ## Normalize the data so it fills to 100%
        #spectra = np.zeros(inspectra.shape)
        spectra = (inspectra.transpose() / np.sum(inspectra, axis=1)).transpose()
        spectra = np.cumsum(spectra[:, :], axis=1)
        #sums = np.sum(inspectra,axis=1)
        #for i,elem in enumerate(inspectra):
        #    spectra[i,:] = elem/sums[i]
    elif cum is True and cum != 'norm':
        spectra = np.cumsum(inspectra[:, :], axis=1)
    else:
        spectra = inspectra

    # Set up plotting environment colors
    Nc = np.array([float(i) / numageclasses for i in range(numageclasses)])
    norm = mpl.colors.normalize(Nc.min(), Nc.max())
    #jet = plt.cm.get_cmap('jet')
    jet = _gen_flexpart_colormap()
    plt.hold('on')
    facecolors = []
    #for i in range(0,H.numageclasses-1):
    for i in range(numageclasses):
        facecolors.append(jet(norm(Nc[i])))
        if labels:
            lbl = labels[i]
        else:
            lbl = i
        #ax.plot(ageclass[:,i])
        if i == 0:
            # create the baseline
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1])
            else:
                ax.fill_between(releasetimes, np.zeros(len(spectra[:, i])), spectra[:, i],
                                color=facecolors[-1], label='%s' % lbl)
        else:
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1], bottom=spectra[:, i - 1])
            else:
                ax.fill_between(releasetimes, spectra[:, i - 1], spectra[:, i],
                                color=facecolors[-1], label='%s' % (lbl))
                #facecolors.append(jet(norm(Nc[i+1])))


    #ax.set_yscale('log')
    if y_datarange:
        print 'setting data range'
        ax.set_ylim(y_datarange)
    else:
        ax.set_ylim((0, spectra.max()))
    ax.grid(True)
    if y_label:
        ax.set_ylabel('%s' % y_label, rotation=90, ha='center', size='small')

    if plt_title:
        plt.title('%s' % (plt_title), size='small')
    if fig_title:
        fig.suptitle('Emissions by %s' % (spectra_type), size='small')

    fig.autofmt_xdate()
    #for xl in ax.get_xticklabels():
    #    plt.setp(xl,size='x-small')
    ## ListedColormap
    pos = ax.get_position()
    l, b, w, h = getattr(pos, 'bounds', pos)
    #polygons = ax.collections
    #for i,p in enumerate(polygons): p.set_color(facecolors[i])
    # BoundaryNorm, and extended ends to show the "over" and "under"
    # value co
    ax2 = fig.add_axes([l, .08, w, 0.03])
    cmap = mpl.colors.ListedColormap(facecolors)
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')

    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.
    bounds = range(numageclasses + 1)
    #norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
    #                                norm=norm,
    #                                # to use 'extend', you must
    #                                # specify two extra boundaries:
                                    boundaries=bounds,
    #                                extend='both',
                                    #values=range(H.numageclasses+1),
                                    ticks=bounds, # optional
                                    spacing='proportional',
                                    orientation='horizontal')
    if spectra_label:
        cb2.set_label(spectra_label, size='x-small')
    if labels:
        ax2.set_xticklabels(labels, va='center', ha='left', rotation=0, size='xx-small')

    if debug:
        plt.show()

    # make main axes current
    fig.add_subplot(ax)

    FIGURE.ax = ax
    FIGURE.fig = fig

    return FIGURE


def plot_agespectra(H, agespectra,
                    plt_title='', y_label=None,
                    cbar_labels=None,
                    FIGURE=None, y_datarange=None,
                    web_version=False, continental=False, cum=False,
                    bars=False, debug=False):
    """ plot an agespectra

    Usage::

        > FIG = plot_agespectra(H,agespectra,*kwargs)

    This creates an filled in between agespectra plot using information from
    the header "H".

    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ========================================
      keyword               Description [default]
      =============         ========================================
      runid                 an ID for the title
      yunits                units for the y-axis [None]
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      data_range            y-axis data range
      web_version           Set to True if using the make_webpages
                            function. see: :func:`read_agespectrum`
      continental           Set to True if continental_spectrum
      cum                   Set to true if data should be cumulative
      bars                  Plot using stacked bars rather than a
                            filled line plot
      =============         ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes, maybe there is a
        better way to handle this.

    .. note::
        Required attributes of 'H' if you want to make a dummy H. Or just set
        H to "None" and it will be taken from the agespectra input.
        H.numageclasses
        H.releasetimes


    """
    ## make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    if FIGURE == None:
        FIGURE = Structure()
        fig = plt.figure()
        FIGURE.fig = fig
        ax = fig.add_subplot(111)
        FIGURE.ax = ax
    else:
        fig = FIGURE.fig
        ax = FIGURE.ax

    if web_version:
        inspectra = agespectra.agespectrum[:, 4:]
    else:
        inspectra = agespectra[:]

    if continental:
        from nilu.pflexpart import emissions as em
        spectra_label = "Continents"
        spectra_type = "Continent"
        conts = em.Continents()
        cbar_labels = [i[1] for i in conts.flexpart_continents]
    else:
        spectra_label = "Ageclasses (Days)"
        spectra_type = "Ageclass"

    if not H:
        H = Structure()
        try:
            numageclasses = inspectra.numageclass
            releasetimes = inspectra.agespectrum[:, 0]
        except:
            numageclasses = inspectra.shape[1] - 1
            releasetimes = inspectra[:, 0]
            inspectra = inspectra[:, 1:]
    else:
        numageclasses = H.numageclasses
        releasetimes = H.releasetimes

    #if cum == 'norm':
    #    ## Normalize the data so it fills to 100%
    #    #spectra = np.zeros(inspectra.shape)
    #    spectra = (inspectra.transpose()/np.sum(inspectra, axis=1)).transpose()
    #    spectra = np.cumsum(spectra[:,:],axis=1)
    #    #sums = np.sum(inspectra,axis=1)
    #    #for i,elem in enumerate(inspectra):
    #    #    spectra[i,:] = elem/sums[i]
    #elif cum is True and cum != 'norm':
    #    spectra = np.cumsum(inspectra[:,:],axis=1)
    #else:
    #    spectra = inspectra
    if cum is not None:
        spectra = _cum_spec(inspectra, cum=cum)

    # Set up plotting environment colors
    Nc = np.array([float(i) / numageclasses for i in range(numageclasses)])
    norm = mpl.colors.normalize(Nc.min(), Nc.max())
    #jet = plt.cm.get_cmap('jet')
    jet = _gen_flexpart_colormap()
    plt.hold('on')
    facecolors = []
    #for i in range(0,H.numageclasses-1):
    for i in range(numageclasses):
        facecolors.append(jet(norm(Nc[i])))
        #ax.plot(ageclass[:,i])
        if i == 0:
            # create the baseline
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1])
            else:
                ax.fill_between(releasetimes, np.zeros(len(spectra[:, i])), spectra[:, i],
                                color=facecolors[-1], label='%s' % i)
        else:
            if bars:
                ax.bar(releasetimes, spectra[:, i], 0.03, color=facecolors[-1], bottom=spectra[:, i - 1])
            else:
                ax.fill_between(releasetimes, spectra[:, i - 1], spectra[:, i],
                                color=facecolors[-1], label='%s' % (i))
                #facecolors.append(jet(norm(Nc[i+1])))


    #ax.set_yscale('log')
    if y_datarange:
        #print 'setting data range'
        ax.set_ylim(y_range)
    else:
        ax.set_ylim((0, spectra.max()))
    ax.grid(True)
    if y_label:
        ax.set_ylabel('%s' % y_label, rotation=90, ha='center', size='small')

    plt.title('%s' % (plt_title), size='small')
    fig.suptitle('Emissions by %s' % (spectra_type), size='small')
    fig.autofmt_xdate()
    #for xl in ax.get_xticklabels():
    #    plt.setp(xl,size='x-small')
    ## ListedColormap
    pos = ax.get_position()
    l, b, w, h = getattr(pos, 'bounds', pos)
    #polygons = ax.collections
    #for i,p in enumerate(polygons): p.set_color(facecolors[i])
    # BoundaryNorm, and extended ends to show the "over" and "under"
    # value co
    ax2 = fig.add_axes([l, .08, w, 0.03])
    cmap = mpl.colors.ListedColormap(facecolors)
    #cmap.set_over('0.25')
    #cmap.set_under('0.75')

    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.
    bounds = range(numageclasses + 1)
    #norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(ax2,
                                    cmap=cmap,
    #                                norm=norm,
    #                                # to use 'extend', you must
    #                                # specify two extra boundaries:
                                    boundaries=bounds,
    #                                extend='both',
                                    #values=range(H.numageclasses+1),
                                    ticks=bounds, # optional
                                    spacing='proportional',
                                    orientation='horizontal')
    cb2.set_label(spectra_label, size='x-small')
    if cbar_labels:
        cb2.set_ticklabels(cbar_labels)

    if debug:
        plt.show()

    return FIGURE


def plot_clusters(H, T, rel_i=0,
                  ncluster=0, sizescale=10000,
                  map_region=None, projection=None, coords=None,
                  FIGURE=None, overlay=False, opacity=0.7,
                  draw_circles=True, draw_labels=True,
                  MapPar=None):
    """ plots the clusters """
    trjs = T['Trajectories']
    labels = T['labels']

    #Set legend properties:
    p_legend = mpl.font_manager.FontProperties(size='8')



    #extract only releases of interest
    rel += 1  #account for zero indexing
    t = trjs[np.where(trjs[:, 0] == rel), :][0]
    if FIGURE == None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)


    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    ## Remove prior retroplume
    if not overlay:
        if len(ax.artists) > 1:
            try:
                art_ind = len(ax.artists)
                del ax.artists[:]
                # del previous texts, but not lat/lon labels
                del ax.texts[-art_ind:]
            except:
                'could not delete artists'
                pass

    #get cluster from data
    #set up cluster index, start column of T array
    cindx = [16, 20, 24, 28, 32, 36]
    if ncluster == 'all':
        #plot all clusters
        ncluster = range(5)
    elif isinstance(ncluster, int):
        ncluster = [ncluster]

    for nc in ncluster:
        #clstrindx = 16 + (5*(nc))
        clstrindx = cindx[nc]
        indx = [1] + range(clstrindx, clstrindx + 4)
        data = t[:, indx]
        #use ellipses
        ells = _genEllipse(data, m, sizescale=sizescale)
        texts = []
        for item in ells:
            e = item[0]
            xy = item[1]
            texts.append([xy, e.get_label()])
            e.set_clip_box(ax.bbox)
            e.set_alpha(opacity)
            if draw_circles:
                ax.add_artist(e)
        for txt in texts:
            x = txt[0][0]
            y = txt[0][1]
            #print x,y, t[1]
            if draw_labels:
                ax.text(x, y, txt[1])
            #plt.colorbar(ax)
    #ax.legend(ax.artists,(o.get_label() for o in ax.artists), prop=p_legend )
    FIGURE.ax = ax
    FIGURE.fig = fig
    FIGURE.m = m

    return FIGURE


def plot_trajectory_ellipses(H, T, rel_i=0,
                  ncluster=0, sizescale=10000,
                  map_region=None, projection=None, coords=None,
                  FIGURE=None, overlay=False, opacity=0.7,
                  draw_circles=True, draw_labels=True,
                  MapPar=None):
    """
    Plots trajectories from FLEXPART output.

    """

    trjs = T['Trajectories']
    labels = T['labels']

    #Set legend properties:
    p_legend = mpl.font_manager.FontProperties(size='8')

    #extract only releases of interest according to rel_i
    rel += 1  #account for zero indexing
    t = trjs[np.where(trjs[:, 0] == rel), :][0]
    if FIGURE == None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)

    try:
        ax = FIGURE.ax
        fig = FIGURE.fig
        m = FIGURE.m
    except:
        print 'problem getting ax,m, or fig.'

    ## Remove prior retroplume
    if overlay:
        if len(ax.artists) > 1:
            try:
                art_ind = len(ax.artists)
                del ax.artists[:]
                # del previous texts, but not lat/lon labels
                del ax.texts[-art_ind:]
            except:
                'could not delete artists'
                pass

    #get trajectory from data
    indx = [1, 2, 3, 4] #return, j,x,y,z
    data = t[:, indx]
    #use ellipses
    ells = _genEllipse(data, m, sizescale=sizescale)
    texts = []
    for item in ells:
        e = item[0]
        xy = item[1]
        texts.append([xy, e.get_label()])
        e.set_clip_box(ax.bbox)
        e.set_alpha(opacity)
        if draw_circles:
            ax.add_artist(e)
    for txt in texts:
        x = txt[0][0]
        y = txt[0][1]
        #print x,y, t[1]
        if draw_labels:
            ax.text(x, y, txt[1])
    plt.draw()
    FIGURE.ax = ax
    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.circles = ells

    return FIGURE

def plot_markers(H, lon, lat, zsize=None, zlevel=None,
                    FIGURE=None,
                    map_region=None, projection=None, coords=None,
                    overlay=True,
                    draw_circles=True,
                    draw_labels=False,
                    cbar2=True,
                    cbar2_title=None,
                    MapPar=None,
                    color='blue', edgecolor='none',
                    zorder=10, alpha=0.85):
    
    """Plot a group of x,y pairs, with optional size and color information.

    Usage::

        > FIG = plot_trajectory(H,RT,rel_id,*kwargs)

    .. note::
        You can substitude "None" for the :class:`Header` instance if you haven't
        created one


    Returns
      A :mod:`mapping` "FIGURE" object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         =============================================
      keyword               Description [default]
      =============         =============================================
      rel_i                 **required** release index
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      projection            A projection pre-defined in :mod:`mapping`
      coords                A set of lon,lat coords for autosetting
                            map map_region (not really working).
      overlay               [True] will overlay the trajectory 
                            on top of another map instance.
      draw_circles          [True] will mark the trajectory with
                            circles. If [False] nothing is drawn.
      draw_labels           [False] unimplemented
      cbar2                 [True] draws the scale bar as a second axis.
      cbar2_title            Optional argument to overide the cbar title.
      MapPar                A Structure of mapping parameters to pass
                            to the basemap instance if desired.
      =============         =============================================

    .. todo::
        Who knows?! Could probably handle the overlaying better.

    .. note::
        Just set H to "None" if you haven't already created a "Header" instance.


    """

    if H:
        pass

    ## Set up the FIGURE
    if FIGURE == None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)
    ##Get fig info and make active
    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.fig.axes[0]
    plt.figure(fig.number)
    plt.axes(ax)

    ## prepare the data
    cx, cy = m(lon, lat)
    if zlevel is None:
        zlevel = np.ones(len(lon)) * 10.
    if zsize is None:
        zsize = np.ones(len(lon)) * 1
    marker = 'o'


    ## clear the previous track
    if 'circles' in FIGURE.keys():
        del FIGURE['circles']
    if overlay is False:
        del ax.collections[FIGURE.indices.collections:]
        del ax.texts[FIGURE.indices.texts:]

    ## plot the track
    if draw_circles:
        cmap = plt.get_cmap('jet')
        circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                          marker=marker, edgecolor=edgecolor,
                          zorder=zorder, alpha=alpha)

        try:
            FIGURE.circles = circles
        except:
            pass
        ## make the figure active again
        plt.figure(fig.number)
        ## draw the legend and title
        ## CREATE COLORBAR
        ## does a colorbar already exist?
        ## Get the current axes, and properties for use later
        ax0 = fig.axes[0]
        pos = ax0.get_position()
        l, b, w, h = pos.bounds
        if cbar2:
            try:
                cb2 = FIGURE.cb2
                cax2 = FIGURE.cax2
                cb2.update_normal(circles)
            except:
                cax2 = plt.axes([l + w + 0.08, b, 0.02, h])
        ## using im2, not im (hack to prevent colors from being
        ## too compressed at the low end on the colorbar - results
        ## from highly nonuniform colormap)
                cb2 = fig.colorbar(circles, cax=cax2) #, format='%3.2g') # draw colorbar
                FIGURE.cax2 = cax2
                FIGURE.cb2 = cb2
        ## check if a colorbar legend exists (this will cause
        ## problems for subplot routines!)
            p_leg = mpl.font_manager.FontProperties(size='6')
            if cbar2_title:
                cax2.set_title(cbar2_title, fontproperties=p_leg)
            else:
                cax2.set_title('altitude\n(m)', fontproperties=p_leg)
        else:
            pass
        #cax2 = plt.axes([l+w+0.03, b, 0.025, h-0.2])
        #    cb2 = fig.colorbar(jnkmap,cax=cax2) # draw colorbar

        ## delete the ghost instance
        #plt.close(jnkfig.number)
        #del jnkax, jnkfig, jnkmap

    plt.axes(ax)
    plt.setp(ax, xticks=[], yticks=[])


    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.ax = ax
    try:
        FIGURE.cax2 = cax2
        FIGURE.cb2 = cb2
    except:
        pass
    return FIGURE

def plot_trajectory(H, T, rel_i, FIGURE=None,
                    map_region=None, projection=None, coords=None,
                    overlay=True,
                    draw_circles=True,
                    draw_labels=True, days_back=20,
                    cbar2=True,
                    cbar2_title=None,
                    MapPar=None):
    """Plot the center trajectory of the plume on the map

    Usage::

        > FIG = plot_trajectory(H,RT,rel_id,*kwargs)

    .. note::
        You can substitude "None" for the :class:`Header` instance if you haven't
        created one


    Returns
      A :mod:`mapping` "FIGURE" object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         =============================================
      keyword               Description [default]
      =============         =============================================
      rel_i                 **required** release index
      FIGURE                A "FIGURE" object[None] (see mapping.py)
      projection            A projection pre-defined in :mod:`mapping`
      coords                A set of lon,lat coords for autosetting
                            map map_region (not really working).
      overlay               [True] will overlay the trajectory 
                            on top of another map instance.
      draw_circles          [True] will mark the trajectory with
                            circles. If [False] nothing is drawn.
      draw_labels           [True] will draw the 'day' numbers on the
                            trajectories.
      days_back             For how many days back should the labels be
                            shown? [20]
      cbar2                 [True] draws the scale bar as a second axis.
      cbar2_title            Optional argument to overide the cbar title.
      MapPar                A Structure of mapping parameters to pass
                            to the basemap instance if desired.
      =============         =============================================

    .. todo::
        Who knows?! Could probably handle the overlaying better.

    .. note::
        Just set H to "None" if you haven't already created a "Header" instance.


    """

    if H:
        pass

    ## Set up the FIGURE
    if FIGURE == None:
        FIGURE = mp.get_FIGURE(map_region=map_region, projection=projection,
                               coords=coords, MapPar=MapPar)
    ##Get fig info and make active
    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.fig.axes[0]
    plt.figure(fig.number)
    plt.axes(ax)

    ## prepare the data
    trjs = T['Trajectories']
    rel = rel_i + 1  ## account for zero indexing
    
    ##extract only releases of interest
    t = trjs[np.where(trjs[:, 0] == rel), :][0]
    
    ## Get the data for the days_back we're interested in
    day_labels = _gen_daylabels(t[:days_back, 1])
    lon = t[:days_back, 2]
    lat = t[:days_back, 3]
    zlevel = t[:days_back, 4]
    zsize = np.ones(len(lon)) * 50
    marker = 'o'


    ## clear the previous track
    if 'circles' in FIGURE.keys():
        del FIGURE['circles']
    if overlay is False:
        del ax.collections[FIGURE.indices.collections:]
        del ax.texts[FIGURE.indices.texts:]

    ## plot the track
    cx, cy = m(lon, lat)
    if draw_circles:
        cmap = plt.get_cmap('gist_gray')
        circles = m.scatter(cx, cy, zsize, zlevel, cmap=cmap,
                          marker=marker, edgecolor=None,
                          zorder=10, alpha=0.85)

        try:
            FIGURE.circles = circles
        except:
            pass
        ## make the figure active again
        plt.figure(fig.number)
        ## draw the legend and title
        ## CREATE COLORBAR
        ## does a colorbar already exist?
        ## Get the current axes, and properties for use later
        ax0 = fig.axes[1]
        pos = ax0.get_position()
        l, b, w, h = pos.bounds
        if cbar2:
            try:
                cb2 = FIGURE.cb2
                cax2 = FIGURE.cax2
                cb2.update_normal(circles)
            except:
                cax2 = plt.axes([l + w + 0.08, b, 0.02, h])
        ## using im2, not im (hack to prevent colors from being
        ## too compressed at the low end on the colorbar - results
        ## from highly nonuniform colormap)
                cb2 = fig.colorbar(circles, cax=cax2) #, format='%3.2g') # draw colorbar
                FIGURE.cax2 = cax2
                FIGURE.cb2 = cb2
        ## check if a colorbar legend exists (this will cause
        ## problems for subplot routines!)
        else:
            cax2 = plt.axes([l + w + 0.03, b, 0.025, h - 0.2])
            cb2 = fig.colorbar(jnkmap, cax=cax2) # draw colorbar

        p_leg = mpl.font_manager.FontProperties(size='6')
        if cbar2_title:
            cax.set_title(cbar2_title, fontproperties=p_leg)
        else:
            cax2.set_title('altitude\n(m)', fontproperties=p_leg)
        ## delete the ghost instance
        #plt.close(jnkfig.number)
        #del jnkax, jnkfig, jnkmap

    plt.axes(ax)
    plt.setp(ax, xticks=[], yticks=[])

    if draw_labels:
        p_cax = mpl.font_manager.FontProperties(size='10',
                                                style='italic',
                                                weight='bold',
                                               )
        
        
        
        for i, (p, x, y) in enumerate(zip(day_labels, cx, cy)):
            if x > m.llcrnrx and x < m.urcrnrx:
                if y > m.llcrnry and y < m.urcrnry:
                    ax.text(x, y+y*.02, '{0}'.format(p), va='bottom', ha='left',
                            fontproperties=p_cax, zorder=11,
                            color='white', bbox=dict(facecolor='green', alpha=0.5)
                            )
                                
    FIGURE.fig = fig
    FIGURE.m = m
    FIGURE.ax = ax
    FIGURE.cax2 = cax2
    FIGURE.cb2 = cb2
    return FIGURE


def plot_at_level(H, D, level=1, \
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
    if D is None:
        D = H.D
    rel_i = D.rel_i
    species = D.species
    timestamp = D.timestamp
    data = D.slabs[level]

    if level == 1:
        level_desc = 'Footprint'
    else:
        level_desc = H.outheight[level - 1]

    dmax = data.max()
    dmin = data.min()
    #print dmax,dmin

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

plot_footprint = plot_at_level

def plot_totalcolumn(H, D=None, \
                   ID=' ', \
                   map_region=5, projection='lcc', \
                   data_range=None, coords=None,
                   FIGURE=None, overlay=False,
                   datainfo_str=None, **kwargs):



    if D is None:
        D = H.D[(0, 0)]
    
    if 'units' in kwargs:
        units = kwargs.pop('units')
    else:
        units = H.output_unit

    rel_i = D.rel_i
    species = D.species
    timestamp = D.timestamp
    data = D.slabs[0]


    if data_range == None:
        dmax = data.max()
        dmin = data.min()
        data_range = [dmin, dmax]
    else:
        dmin, dmax ,= data_range
        #print dmin,dmax
        

    if H.direction == 'backward':
        rel_i = D.rel_i
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


def plot_sourcecontribution(H, D,
                            rel_i=None, s=None,
                            ID=' ',
                            map_region='POLARCAT',
                            data_range=None, coords=None,
                            datainfo_str=None,
                            units=None,
                            FIGURE=None,
                            overlay=False,
                            cax_title=None,
                            **kwargs):
    """ plot_sourcecontribution: a wrapper around :func:`plot_sensitivity` for
    plotting a source contribution field.

    .. warning:
        This function is still quite finicky, I'm working on resolving it.

    Usage::

        > FIG = plot_sourcecontribution(H,D,*kwargs)

    Inputs
      H = a :class:`Header` instance from a FLEXPART run
      D = a grid instance from a FLEXDATA dictionary. See :func:`readgridV8`

      .. note::
          I am trying to create this so D can be a number of different types
          from a 2-d array to the grid instance type.


    Returns
      A "mapping.py" ``FIGURE`` object.

    Arguments

      .. tabularcolumns::  |l|L|

      =============         ================================================
      keyword               Description [default]
      =============         ================================================
      ID                    Can be passed to prefix the title.
      data_range            range of data for scale bars, if None will
                            be taken from min/max of data
      datainfo_str          A string for labeling the scale bar.
      cax_title             A string to pass to plot_sensitivity for colorbar
                            title (units will be passed as format argument)
      k                     Release index to plot from the data array (aka rel_i)
      s                     Species index
      map_region                A map_region specified in mapping.py
      projection            [deprecated] use pre-defined map_regions.
      overlay               Force removal of previous figure elements.
      FIGURE                A FIGURE instance from mapping module get_FIGURE
      =============         ================================================

    .. todo::
        A lot!! There are some problems here and it is sensitive to options.

    .. note::
        k and rel_i are used throughout pflexible, should be more consistent in future.


    """
    if D is None:
        D = H.D[(0, 0)]
        data = D.slabs[1]
    elif isinstance(D, list):
        data = D[0]
    elif isinstance(D, np.ndarray):
        data = D
    else:
        try:
            D = D[(s, rel_i)]
            data = D.slabs[0] #take total column
        except:
            raise IOError("Unrecognized 'grid' type passed")

    if s is None:
        try:
            species = D.species
        except:
            species = 'unknown'
    else:
        species = H.species[s]

    try:
        timestamp = D.timestamp
    except:
        if rel_i is None:
            try:
                timestamp = D.timestamp
                rel_i = D.rel_i
            except:
                raise IOError("Either a release index (rel_i) must be passed, \
                           or the grid must have a timestamp attribute")
        else:
            timestamp = H.releasetimes[rel_i]

    if units is None:
        units = 'emission units'

    dmax = data.max()
    dmin = data.min()

    if data_range == None:
        data_range = [dmin, dmax]

    if H.direction == 'backward':
        zp1 = H['zpoint1'][rel_i]
        zp2 = H['zpoint2'][rel_i]
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f , Z2: %.2f (%s)\n""" % \
                                (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n
        Release Start: %s, Release End: %s""" % \
                       (ID, species, H['releasestart'][rel_i], H['releaseend'][rel_i])
    else:
        if datainfo_str is None:
            datainfo_str = """ Max Value: %.2g %s\n Release Z1: %.2f , Z2: %.2f (%s)\n""" % \
                                (dmax, units, zp1, zp2, H.alt_unit)
        plot_title = """
        %s Total Column Sensitivity: %s\n %s""" % (ID, species, timestamp)

    FIGURE = plot_sensitivity(H, data, \
                           data_range=data_range, \
                           rel_i=rel_i, map_region=map_region, \
                           units=units, \
                           datainfo_str=datainfo_str, coords=coords, \
                           FIGURE=FIGURE,
                           cax_title=cax_title,
                           overlay=overlay, **kwargs)
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
    #autofit=False,method='imshow'):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::

        > FIG = plot_sensitivity(H,data,*kwargs)

    This returns the FIGURE object, and plots the sensitivity from the data contained in the "D"
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

    print FIGURE
    if dropm != None:
        try:
            del m
            plt.close('all')
        except:
            print 'could not drop m'

    ## make tick lables smaller
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax

    ## make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    ## set up transformations for the data array
    if method == 'imshow':
        if m.projection not in ['cyl', 'merc', 'mill']:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)[:-1]
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)[:-1]
            data = data[:-1, :-1]
        else:
            lats = np.arange(H.outlat0, (H.outlat0 + (H.numygrid * H.dyout)), H.dyout)
            lons = np.arange(H.outlon0, (H.outlon0 + (H.numxgrid * H.dxout)), H.dxout)

        ## transform to nx x ny regularly spaced native projection grid
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
        ## Check to see if a cyclic wraparound is required

        lons = np.arange(H.outlon0, H.outlon0 + (H.dxout * H.numxgrid), H.dxout)
        lats = np.arange(H.outlat0, H.outlat0 + (H.dyout * H.numygrid), H.dyout)
        ## if required add polar coordinates
        if 'npstere' in m.projection:
            if lats[-1] != 90.:
                npole = np.ones(len(lons)).T * data[0, :]
                npole = np.reshape(npole, (1, -1))
                data = np.vstack((npole, data))
                lats = np.hstack((lats, [lats[-1] + H.dyout]))

        if m.projection != 'merc':
            if lons[-1] - lons[0] < 360.:
                topodat, lons = addcyclic(data, lons)

        if m.projection == 'merc':
            topodat = data

    ## get min/max range
    if data_range != None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = _data_range(data)
        
        
    if log:
        clevs = _log_clevs(dat_min, dat_max)
        
    else:
        clevs = [i for i in np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    ## draw land sea mask
    #m.fillcontinents(zorder=0)
    if lsmask:
        m.drawlsmask(ocean_color='grey', zorder= -10)

    ## Plot Release Location if points were read
    if H.options['readp']:
        if rel_i:
            releaselocation = (H.xpoint[rel_i], H.ypoint[rel_i])
            xpt, ypt = m(releaselocation[0], releaselocation[1])
            ## Remove prior location point
            try:
                del ax.lines[-1]
            except:
                pass
            location, = m.plot([xpt], [ypt], 'bx', linewidth=6, markersize=20, zorder=1000)

    ## Plot the footprint

    ## Set up the IMAGE
    ## cmapnames = ['jet', 'hsv', 'gist_ncar', 'gist_rainbow', 'cool', 'spectral']
    #colmap = plt.get_cmap('jet')
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    ## Plotting METHODS (pcolormesh now default, imshow is smoother)
    #print topodat.max(), topodat.min(), topodat.shape
    if method == 'imshow':
        im = m.imshow(topodat, cmap=colmap, zorder= -1,
                      norm=mpl.colors.LogNorm(vmin=clevs[0],
                                              vmax=clevs[-1]))

    if method == 'pcolormesh':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.pcolormesh(nx, ny, topodat, cmap=colmap,
                          norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                  vmax=clevs[-1]))
    if method == 'contourf':
        ## Trying some fancier scaling
        #cnts,bins = np.histogram(topodat,bins=100)
        #topodat = np.ma.masked_where(topodat< .05* np.average((0,bins[1])),topodat)
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    if method == 'contour':
        nx, ny = m(*np.meshgrid(lons, lats))
        im = m.contour(nx, ny, topodat, cmap=colmap,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    ## Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    ## CREATE COLORBAR
    ## Note, with upgrades to matplotlib and basemap had to make some
    ## changes here... no more 'ghost' axes
    ## does a colorbar already exist?
    try:
        cb = FIGURE.cb
        cax = FIGURE.cax
        cb.update_normal(im)
    except:
    ## make a copy of the image object, change
    ## colormap to linear version of the precip colormap.
    #pdb.set_trace()
    #im2 = copy.copy(im)
    #im2.set_cmap(colmap)
    ## create new axis for colorbar.
        h = 0.5 * h
        l = l + w + .03
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
    ## using im2, not im (hack to prevent colors from being
    ## too compressed at the low end on the colorbar - results
    ## from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax) #, format='%3.2g') # draw colorbar
        FIGURE.cax = cax
        FIGURE.cb = cb
    #cb.update_normal(im2)




    ## set colorbar label and ticks
    #pdb.set_trace()
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10]) ##clevs, by 10 steps
    clabels.append(clevs[-1]) ## add the last label
    #cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])
                        #fontproperties=p_cax)
    if H.direction == 'forward':
        cax.set_title('%s' % units,
                      fontproperties=p_cax)
    else:
        if cax_title:
            cax.set_title(cax_title.format(units), fontproperties=p_cax)
        else:
            cax.set_title('sensitivity\n({0})'.format(units),
                          fontproperties=p_cax)

    ## make the original axes current again
    plt.axes(ax)

    ## write text information block on plot
    ## first try to remove prior text by accessing
    ## the last text element added to the axes from
    ## the prior iteration.
    ## This is tricky when using together with plot_clusters...
    ## need to figure out how to resolve the indexing
    ## of what texts, collections, etc to delete, when iterating.
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
        #plt.title(plottitle,fontproperties=p_cax)
    #plt = plt
        FIGURE.ax.set_title(plottitle, fontsize=10)
    return FIGURE


def plot_curtain(H, data, \
                 nx=None,
                 ny=None,
                 data_range=None, \
                 units='ppbv', \
                 datainfo_str=None, \
                 asl = True,
                 plottitle=None, \
                 log=True,
                 FIGURE=None,
                 cax_title=None,
                 method='contourf',
                 figPar=None):
    """ plot_sensitivity: core function for plotting FLEXPART output.

    Usage::

        > FIG = plot_sensitivity(H,data,*kwargs)

    This returns the FIGURE object, and plots the sensitivity from the data contained in the "D"
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
      asl                   [True] plot the units in m.a.s.l
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
        FIGURE = Structure()
        
            
        fig= plt.figure(**figPar)
        ax = fig.add_subplot(111)
        
        FIGURE['fig'] = fig
        FIGURE['ax'] = ax
    
    
    fig = FIGURE.fig
    ax = FIGURE.ax
    
    ## make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    ## get min/max range
    if data_range != None:
        dat_min = data_range[0]
        dat_max = data_range[1]
    else:
        dat_min, dat_max = _data_range(data)
        
    if log:
        clevs = _log_clevs(dat_min, dat_max)
    else:
        clevs = [i for i in np.arange(dat_min, dat_max, (dat_max - dat_min) / 100)]

    ## Set up the IMAGE
    ## cmapnames = ['jet', 'hsv', 'gist_ncar', 'gist_rainbow', 'cool', 'spectral']
    colmap = _gen_flexpart_colormap()
    colmap.set_over(color='k', alpha=0.8)
    topodat = data
    
     
    if method == 'imshow':
        im = plt.imshow(np.flipud(topodat), cmap=colmap, zorder= -1,
                      norm=mpl.colors.LogNorm(vmin=clevs[0],
                                              vmax=clevs[-1]))

    if method == 'pcolormesh':
        im = plt.pcolormesh(nx, ny, topodat, cmap=colmap,
                          norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                  vmax=clevs[-1]))
    if method == 'contourf':
        im = plt.contourf(nx, ny, topodat, cmap=colmap, levels=clevs,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    if method == 'contour':
        im = plt.contour(nx, ny, topodat, cmap=colmap,
                        norm=mpl.colors.LogNorm(vmin=clevs[0],
                                                vmax=clevs[-1]))

    ## Get the current axes, and properties for use later
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    ## CREATE COLORBAR
    ## Note, with upgrades to matplotlib and basemap had to make some
    ## changes here... no more 'ghost' axes
    ## does a colorbar already exist?
    try:
        cb = FIGURE.cb
        cax = FIGURE.cax
        cb.update_normal(im)
    except:
    ## make a copy of the image object, change
    ## colormap to linear version of the precip colormap.
    ## create new axis for colorbar.
        h = 0.8 * h
        l = l + w + .02
        b = 0.5 - (h / 2)
        w = 0.025
        cax = plt.axes([l, b, w, h])
    ## using im2, not im (hack to prevent colors from being
    ## too compressed at the low end on the colorbar - results
    ## from highly nonuniform colormap)
        cb = fig.colorbar(im, cax=cax) #, format='%3.2g') # draw colorbar
        FIGURE.cax = cax
        FIGURE.cb = cb


    ## set colorbar label and ticks
    p_cax = mpl.font_manager.FontProperties(size='6')
    clabels = list(clevs[::10]) ##clevs, by 10 steps
    clabels.append(clevs[-1]) ## add the last label
    #cax.set_yticks(np.linspace(clabels[0],clabels[-1],len(clabels)))
    cax.set_yticks(np.linspace(0, 1, len(clabels)))
    cax.set_yticklabels(['%3.2g' % cl for cl in clabels])
                        #fontproperties=p_cax)
    
    if cax_title:
        cax.set_title(cax_title.format(units), fontproperties=p_cax)
    else:
        cax.set_title('sensitivity\n({0})'.format(units),
                          fontproperties=p_cax)

    ## make the original axes current again
    plt.axes(ax)
    plt.grid(True)
    FIGURE.ax = ax
    FIGURE.fig = fig

    if plottitle != None:
        FIGURE.ax.set_title(plottitle, fontsize=10)
        
    return FIGURE


def plot_METDATA(METDATA, FIGURE, date=None, level=None):
    """ plots met data returned from :module:`mp.get_OpenDAP`

    """


    fig = FIGURE.fig
    m = FIGURE.m
    ax = FIGURE.ax

    ## make the figure current
    plt.figure(fig.number)
    plt.axes(ax)

    datelabels = METDATA['extracted_dates']
    slpin = METDATA['prmslmsl']['data'] * .01
    uin = METDATA['ugrdprs']['data']
    vin = METDATA['vgrdprs']['data']
    longitudes = METDATA['longitudes']
    latitudes = METDATA['latitudes']

    if date:
        cnt = datelabels.index(date)
    else:
        cnt = 0
    ## FOR OPENDAP OVERLAY
    slp_ = slpin[cnt, :, :]
    u_ = uin[cnt, :, :]
    v_ = vin[cnt, :, :]

    # plot wind vectors on projection grid (looks better).
    # first, shift grid so it goes from -180 to 180 (instead of 0 to 360
    # in longitude).  Otherwise, interpolation is messed uplt.

    if longitudes[-1] - longitudes[0] < 360.:
        longarray = np.array(longitudes)
        slp, cyclon = addcyclic(slp_, longarray)
        u, cyclon = addcyclic(u_, longarray)
        v, cyclon = addcyclic(v_, longarray)
        ugrid, newlons = shiftgrid(180., u, cyclon, start=False)
        vgrid, newlons = shiftgrid(180., v, cyclon, start=False)
        slpgrid, newlons = shiftgrid(180., slp, cyclon, start=False)
    # transform vectors to projection grid.

    lons, lats = np.meshgrid(newlons, latitudes)
    # set desired contour levels.
    clevs = np.arange(960, 1061, 5)
    # compute native x,y coordinates of grid.
    x, y = m(lons, lats)
    #pos = FIG.ax.get_position()
    #l, b, w, h = pos.bounds

    CS = m.contour(x, y, slpgrid, clevs, linewidths=1, colors='k', animated=True)
    #CS = FIG.m.contourf(x,y,slpgrid,clevs,cmap=plt.cm.RdBu_r,animated=True,alpha=0.7)
   # CS = FIG.m.contour(x,y,slp[0,:,:],clevs,linewidths=0.5,colors='k',animated=True)
#            CS = FIG.m.contourf(x,y,slp[0,:,:],clevs,cmap=plt.cm.RdBu_r,animated=True,alpha=0.7)
    #plt.clabel(CS,inline=1,fontsize=10)
    # plot wind vectors over maplt.
    urot, vrot, xx, yy = m.transform_vector(ugrid, vgrid, newlons, latitudes, 51, 51,
                                         returnxy=True, masked=True)
    Q = m.quiver(xx, yy, urot, vrot, scale=500)
#            # make quiver key.
    qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')


    FIGURE.ax = ax
    FIGURE.m = m
    FIGURE.fig = fig

    return FIGURE












########### HELPER FUNCTIONS ##########

class BinaryFile(object):

    """
    BinaryFile: A class for accessing data to/from large binary files
    =================================================================

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
            raise ValueError, "order should be either 'fortran' or 'c'."
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
        if type(shape) is int:
            shape = (shape,)
        if type(shape) is not tuple:
            raise ValueError, "shape must be a tuple"
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
            raise EOFError, "Asking for more data than available in file."

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



class Structure(dict):
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

    def set_with_dict(self, D):
        """ set attributes with a dict """
        for k in D.keys():
            self.__setattr__(k, D[k])



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
      readp_ff            readp_ff (read releases using Fortran [False]
      nest                nested output [0]=no, 1=yes
      version             version of FLEXPART, default = 'V8'
      =============       ========================================
  
    
    .. note::
        **This function is in development**

        This function is being developed so that there is no dependence on
        using f2py to compile the FortFlex module. It is working using the
        :class:`BinaryFile`, but is notably slower than :class:`FortFlex`. 
        Please report any bugs found.

  
    """


    def __init__(self, path=None, headerfile=None, **readheader_ops):
        """

        
        """

        try:
            h = readheader(path, **readheader_ops)
            self.set_with_dict(h)
            self.lonlat()
            self.version = 'V8'
        except:
            try:
                h = readheaderV6(path, **readheader_ops)
                self.set_with_dict(h)
                self.lonlat()
                self.version = 'V6'
            except:
                print path
                raise IOError('Could not set header variables. Does the path exist?')

    def lonlat(self):
        """ Add longitude and latitude attributes using data from header """
        lons = np.arange(self.outlon0, self.outlon0 + (self.dxout * self.numxgrid), self.dxout)
        lats = np.arange(self.outlat0, self.outlat0 + (self.dyout * self.numygrid), self.dyout)
        self.longitude = lons
        self.latitude = lats

    def read_grid(self, **kwargs):
        """ see :func:`read_grid` """
        self.FD = read_grid(self, **kwargs)

    def fill_backward(self, **kwargs):
        """ see :func:`fill_backward` """
        fill_backward(self, add_attributes=True, **kwargs)

    def add_trajectory(self, **kwargs):
        """ see :func:`read_trajectories` """
        self.trajectory = read_trajectories(self)

    def add_fires(self, **kwargs):
        """ uses the :mod:`emissions` module to read the MODIS hotspot data and
        add it to the header class as a 'fires' attribute.
        """

        from nilu.pflexpart import emissions as em
        self.fires = None
        for day in self.available_dates_dt:
            # day = day[:8]
            firedata = em.MODIS_hotspot(day)
            daily = firedata.daily
            if self.fires == None:
                self.fires = daily
            else:
                self.fires = np.hstack((self.fires, daily)).view(np.recarray)

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

def _data_range(data, min='median'):
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
        print 'dat_max not positive'
        print dat_max
        dmx = 1
            
    if dat_min > 0:
        dmn = int(np.round(np.log10(dat_min)))
    elif dat_min == 0. or np.isnan(dat_min):
        print 'dat_min not positive'
        print dat_min
        dmn = dmx - 3
        
        
    ## create equally spaced range
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
        ## Normalize the data so it fills to 100%
        #spectra = np.zeros(inspectra.shape)
        spectra = (inspectra.transpose() / np.sum(inspectra, axis=1)).transpose()
        spectra = np.cumsum(spectra[:, :], axis=1)
        #sums = np.sum(inspectra,axis=1)
        #for i,elem in enumerate(inspectra):
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
    #p = p/altmax
    cm = cmap(norm(p))
    #print p, cm
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
    if index == None:
        seek = range(G.shape[-1])
    else:
        seek = [index]

    for i in seek:
        zpmax = np.max(G[:, :, 0, i] / Heightnn[:, :, 0])
        if zpmax > fpmax:
            fpmax = zpmax
        #print fpmax
    tcmax = np.max(G)
    return ((0, fpmax), (0, tfcmax))


def _genEllipse(data, m, sizescale=20000,
                altlim=None):
    """
    Generates ellipses based on input array 'data'. Requires basemap instance, 'm'. NOTE:

    data = [itime,x,y,z,[size]]

    r,c = data.shape
    if c == 5:
        size function of data[:,4]
    else:
        size = 1*sizescale

    uses function:
    _gen_altitude_color
    _gen_daylabels
    for labeling/color of ellipses

    """
    r, c = data.shape
    if altlim == None:
        altlim = (np.min(data[:, 3]), np.max(data[:, 3]))

    if c == 5:

        ell = [(Ellipse(xy=np.array(m(p[1], p[2])),
                width=p[4] * sizescale, height=p[4] * sizescale,
                angle=360,
                facecolor=_gen_altitude_color(p[3], altlim=altlim, cmap_id='gray'), \
                label=_gen_daylabels(p[0])), \
                np.array(m(p[1], p[2]))) for p in data]
                     #np.array( m(p[1],p[2])) ) for p in data]
    else:
        ell = [(Ellipse(xy=np.array(m(p[1], p[2])),
                width=1e4 * sizescale, height=1e4 * sizescale,
                angle=360, \
                facecolor=_gen_altitude_color(p[3], altlim=altlim, cmap_id='gray'), \
                label=_gen_daylabels(p[0])), \
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
        ## AST Colorset for FLEXPART
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



#### Below for use from command line  ##########################################

def main ():

    global options, args
    here = os.path.curdir
    dataDir = options.pathname
    H = read_header(dataDir, **{'readp':1}) #readheader and release pointspiu
    G = readgridV8(H, **{'unit':'time'})
    D = get_slabs(H, G)

    plt = plot_map(D[0], H)
    plot_name = 'test'
    plt.title(plot_name, fontsize=10)
    plt.savefig(os.path.join(here, plot_name + '.png'))
    print plot_name


if __name__ == '__main__':
    print 'testing'
    try:
        start_time = time.time()
        parser = optparse.OptionParser(
            formatter=optparse.TitledHelpFormatter(),
            usage=globals()['__doc__'],
            version='$Id: py.tpl 327 2008-10-05 23:38:32Z root $')
        parser.add_option ('-v', '--verbose', action='store_true',
                           default=False, help='verbose output')
        parser.add_option ('-T', '--latex', action='store_true',
                           default=False, help='latex usage flag')
        parser.add_option ('-d', '--directory', dest='pathname',
                           default=False, help='pathname of directory containg FLEXPART Output')
        (options, args) = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
        if options.verbose: print time.asctime()

        #### LateX #####
        if options.latex:
            mpl.rc('font', **{'family':'sans-serif',
                              'sans-serif':['Helvetica']}
                   )

            mpl.rc('text', usetex=True)

        exit_code = main()
        if exit_code is None:
            exit_code = 0
        if options.verbose: print time.asctime()
        if options.verbose: print 'TOTAL TIME IN MINUTES:',
        if options.verbose: print (time.time() - start_time) / 60.0
        sys.exit(exit_code)
    except KeyboardInterrupt, e: # Ctrl-C
        raise e
    except SystemExit, e: # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

# vim:set sr et ts=4 sw=4 ft=python fenc=utf-8: // See Vim, :help 'modeline'

