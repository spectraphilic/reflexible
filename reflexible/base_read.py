"""
functions for reading fp data in reflexible.

SYNOPSIS
========

    reflexible.base_read

DESCRIPTION
===========

    base_read: module for reading FLEXPART Output.


AUTHOR
======

    JFB: John F Burkhart <jfburkhart@gmail.com>

CONTRIBUTORS
============

    Francesc Alted <falted@gmail.com>
LICENSE
=======

    This script follows creative commons usage.

VERSION
=======

    ID: $Id$: $Rev$
"""
# builtin imports
import os
import datetime as dt

# Dependencies:
# Numpy
import numpy as np

from .data_structures import Trajectory


def read_trajectories(H, trajfile='trajectories.txt', ncluster=5,
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
            alltraj = open(H, 'r').readlines()
        except:
            raise IOError('Could not open file: %s' % H)
    else:
        path = H.fp_path
        alltraj = open(os.path.join(path, trajfile), 'r').readlines()

    try:
        ibdate, ibtime, model, version = alltraj[0].strip().split()[:4]
    except:
        ibdate, ibtime = alltraj[0].strip().split()[:2]
        model = 'Flexpart'
        version = 'V.x'
    tdelta = dt.datetime.strptime(ibdate + ibtime.zfill(6), '%Y%m%d%H%M%S')
    numpoint = int(alltraj[2].strip())

    RelTraj = Trajectory()
    Trajectories = []
    # format change?? Try block below because some trajectories.txt
    # files have linebreaks, but seems more recent FP v10. does not.
    try:
        for i in range(3, 3 + (numpoint * 2), 3):
            i1, i2, xp1, yp1, xp2, = tuple((float(j) for j in
                                            alltraj[i].strip().split()))
            yp2, zp1, zp2, k, npart, = tuple((float(j) for j in
                                              alltraj[i + 1].strip().split()))
    except:
        for i in range(3, 3 + (numpoint * 2), 2):
            i1, i2, xp1, yp1, xp2, yp2, zp1, zp2, k, npart, = \
                tuple([float(j) for j in alltraj[i].strip().split()])

        itimerel1 = tdelta + dt.timedelta(seconds=i1)
        itimerel2 = tdelta + dt.timedelta(seconds=i2)
        Xp = (xp1 + xp2) / 2
        Yp = (yp1 + yp2) / 2
        Zp = (zp1 + zp2) / 2
        RelTraj[alltraj[i + 1].strip()] = np.array(
            (itimerel1, itimerel2, Xp, Yp, Zp, k, npart))

    for i in range(3 + (numpoint * 3), len(alltraj)):
        raw = alltraj[i]
        FMT = [0, 5, 8, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6] + \
              ncluster * [8, 8, 7, 6, 8]

        data = [raw[sum(FMT[:ii]):sum(FMT[:ii + 1])]
                for ii in range(1, len(FMT) - 1)] + [raw[sum(FMT[:-1]):]]
        ### FIX ###
        # ## To get rid of '******' that is now in trajectories.txt
        data = [float(r.replace('********', 'NaN')) for r in data]

        Trajectories.append(data)

    data = np.array(Trajectories)
    RelTraj['version'] = model + ' ' + version
    RelTraj['date'] = tdelta
    RelTraj['Trajectories'] = data
    RelTraj['labels'] = \
        ['release number', 'seconds prior to release', 'lon', 'lat',
         'height', 'mean topography',
         'mean mixing height', 'mean tropopause height', 'mean PV index',
         'rms distance', 'rms', 'zrms distance', 'zrms',
         'fraction mixing layer', 'fraction PV<2pvu',
         'fraction in troposphere'] + \
        ncluster * ['xcluster', 'ycluster', 'zcluster', 'fcluster',
                    'rmscluster']
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
