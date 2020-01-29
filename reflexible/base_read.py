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
import os
import datetime as dt

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


# The amount of records that are read from disk in a single shot.  Keep this as
# small as possible.
CHUNKSIZE = 10 * 1000


def get_quantized_ctable(dtype, cparams, quantize=None, expectedlen=None):
    """Return a ctable with the quantize filter enabled for floating point cols.
    """
    import bcolz

    columns, names = [], []
    for fname, ftype in dtype.descr:
        names.append(fname)
        if 'f' in ftype:
            cparams2 = bcolz.cparams(clevel=cparams.clevel, cname=cparams.cname, quantize=quantize)
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams2, expectedlen=expectedlen))
        else:
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams, expectedlen=expectedlen))
    return bcolz.ctable(columns=columns, names=names)


def read_partpositions(filename, nspec, ctable=True, clevel=5, cname="lz4", quantize=None):
    """Read the particle positions in `filename`.

    This function strives to use as less memory as possible; for this, a
    bcolz ctable container is used for holding the data.  Besides to be compressed
    in-memory, its chunked nature makes a natural fit for data that needs to
    be appended because it does not need expensive memory resize operations.

    NOTE: This code reads directly from un UNFORMATTED SEQUENTIAL data Fortran
    file so care has been taken to skip the record length at the beginning and
    the end of every record.  See:
    http://stackoverflow.com/questions/8751185/fortran-unformatted-file-format

    Parameters
    ----------
    filename : string
        The file name of the particle raw data
    nspec : int
        number of species in particle raw data
    ctable : bool
        Return a bcolz ctable container.  If not, a numpy structured array is returned instead.
    clevel : int
        Compression level for the ctable container
    cname : string
        Codec name for the ctable container.  Can be 'blosclz', 'lz4', 'zlib' or 'zstd'.
    quantize : int
        Quantize data to improve (lossy) compression.  Data is quantized using
        np.around(scale*data)/scale, where scale is 2**bits, and bits is
        determined from the quantize value.  For example, if quantize=1, bits
        will be 4.  0 means that the quantization is disabled.

    Returns
    -------
    ctable object OR structured_numpy_array

    Returning a ctable is preferred because it is used internally so it does not require to be
    converted to other formats, so it is faster and uses less memory.

    Note: Passing a `quantize` param > 0 can increase the compression ratio of the ctable
    container, but it may also slow down the reading speed significantly.
    """
    import bcolz

    xmass_dtype = [('xmass_%d' % (i + 1), 'f4') for i in range(nspec)]
    # note age is calculated from itramem by adding itimein
    out_fields = [
                     ('npoint', 'i4'), ('xtra1', 'f4'), ('ytra1', 'f4'), ('ztra1', 'f4'),
                     ('itramem', 'i4'), ('topo', 'f4'), ('pvi', 'f4'), ('qvi', 'f4'),
                     ('rhoi', 'f4'), ('hmixi', 'f4'), ('tri', 'f4'), ('tti', 'f4')] + xmass_dtype
    raw_fields = [('begin_recsize', 'i4')] + out_fields + [('end_recsize', 'i4')]
    raw_rectype = np.dtype(raw_fields)
    recsize = raw_rectype.itemsize

    cparams = bcolz.cparams(clevel=clevel, cname=cname)
    if quantize is not None and quantize > 0:
        out = get_quantized_ctable(raw_rectype, cparams=cparams, quantize=quantize, expectedlen=int(1e6))
    else:
        out = bcolz.zeros(0, dtype=raw_rectype, cparams=cparams, expectedlen=int(1e6))

    with open(filename, "rb", buffering=1) as f:
        # The timein value is at the beginning of the file
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        itimein = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        nrec = 0
        while True:
            # Try to read a complete chunk
            data = f.read(CHUNKSIZE * recsize)
            read_records = int(len(data) / recsize)  # the actual number of records read
            chunk = np.ndarray(shape=(read_records,), buffer=data, dtype=raw_rectype)
            # Add the chunk to the out array
            out.append(chunk[:read_records])
            nrec += read_records
            if read_records < CHUNKSIZE:
                # We reached the end of the file
                break

    # Truncate at the max length (last row is always a sentinel, so remove it)
    out.trim(1)
    # Remove the first and last columns
    out.delcol("begin_recsize")
    out.delcol("end_recsize")

    if ctable:
        return out
    else:
        return out[:]


def read_shortpositions(filename, nspec, xlon0, ylat0, dx, dy, dataframe=False, header=None):
    """Read the particle positions in filename.

    Parameters
    ----------
    filename : string
      the file name of the particle raw files
    nspec : int
      number of species
    xlon0 : float
      the longitude
    xlat0 : float
      the latitude
    dx : float
      the cell size (x)
    dy : float
      the cell size (u)
    dataframe : bool
      Return a pandas DataFrame

    Returns
    -------
    structured_numpy_array OR pandas DataFrame
    """
    #TODO: read_shortposit see partoutput_short.f90
    pass

