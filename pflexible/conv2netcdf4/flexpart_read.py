#### Functions for reading FLEXPART output #####

import os
import re
import datetime
from math import pi, sqrt, cos

import numpy as np

import pflexible.conv2netcdf4
from .helpers import closest


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

    ================    =================================================================
    fields              description
    ================    =================================================================
    SIM_DIR             Simulation direction
    SIM_START           Text str of YYYYMMDD HHMMSS
    SIM_END             Text str of YYYYMMDD HHMMSS
    AVG_CNC_INT         Average concentrations are calculated every SSSSS seconds
    AVG_CNC_TAVG        The average concentrations are time averages of SSSSS sec
    CNC_SAMP_TIME       The concentrations are sampled every SSSSS seconds to
                        calculate the time average concentration.
    T_PARTSPLIT         Time constant for particle splitting.
    SYNC                All processes are synchronized with this time interval
    CTL                 --
    IFINE               IFINE=Reduction factor for time step used for vertical wind
    IOUT                IOUT determines how the output shall be made: concentration
                        (ng/m3, Bq/m3), mixing ratio (pptv), or both, or plume
                        trajectory mode, or concentration + plume trajectory mode.
    IPOUT               IPOUT determines whether particle positions are outputted
                        (in addition to the gridded concentrations or mixing ratios)
                        or not. 0=no output, 1 output every output interval, 2 only
                        at end of the simulation
    LSUBGRID            Switch on/off subgridscale terrain parameterization
                        (increase of mixing heights due to subgridscale orog. var
    LCONVECTION         Switch on/off the convection parameterization
    LAGESPECTRA         Switch on/off the calculation of age spectra: if yes, the
                        file AGECLASSES must be available
    IPIN                If IPIN=1, a file "partposit_end" from a previous run must
                        be available in the output directory. Particle positions
                        are read in and previous simulation is continued. If
                        IPIN=0, no particles from a previous run are used
    IOFR                 Switch on/off writing out each release.
    IFLUX               If IFLUX is set to 1, fluxes of each species through each
                        of the output boxes are calculated. Six fluxes,
                        corresponding to northward, southward, eastward, westward,
                        upward and downward are calculated for each grid cell of
                        the output grid. The control surfaces are placed in the
                        middle of each output grid cell. If IFLUX is set to 0,
                        no fluxes are determined.
    MDOMAINFILL         If MDOMAINFILL is set to 1, the first box specified in file
                        RELEASES is used as the domain where domain-filling
                        trajectory calculations are to be done. Particles are
                        initialized uniformly distributed (according to the air mass
                        distribution) in that domain at the beginning of the
                        simulation, and are created at the boundaries throughout
                        the simulation period
    IND_SOURCE          IND_SOURCE switches between different units for
                        concentrations at the source. NOTE that in backward
                        simulations the release of computational particles takes
                        place at the "receptor" and the sampling of particles at
                        the "source".  1=mass units (for bwd-runs = concentration)
                        2=mass mixing ratio units'''],

    IND_RECEPTOR        IND_RECEPTOR switches between different units for
                        concentrations at the receptor 1=mass units (concentrations)
                        2=mass mixing ratio units

    MQUASILAG           MQUASILAG indicates whether particles shall be numbered
                        consecutively (1) or with their release location number (0).
                        The first option allows tracking of individual particles
                        using the partposit output files

    NESTED_OUTPUT       NESTED_OUTPUT decides whether model output shall be made
                        also for a nested output field (normally with higher resolution)
    LINIT_COND          For Backward Runs, sets initial conditions:
                        [0]=No, 1=Mass Unit, 2=Mass Mixing

    ================    =================================================================

    """
    lines = file(path, 'r').readlines()
    command_vals = [i.strip() for i in lines[headerrows:]]  # clean line ends
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
    float_keys = ['CTL']
    date_keys = ['SIM_START', 'SIM_END']

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

    def getfile_lines(infile):
        """ returns all lines from a file or file string
        reverts to beginning of file."""

        if isinstance(infile, str):
            return file(infile, 'r').readlines()
        else:
            infile.seek(0)
            return infile.readlines()

    lines = getfile_lines(path)
    lines = [i.strip() for i in lines]  # clean line ends

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

    names = ['start_time', 'end_time', 'lllon', 'lllat', 'urlon', 'urlat',
             'altunit', 'elv1', 'elv2', 'numpart']
    # formats = [object, object, np.float, np.float, np.float, np.float,\
    #                      int, np.float, np.float, int]
    for i in range(nspec):
        names.append('mass%s' % i)
        # formats.append(np.float)
    names.append('id')
    # formats.append('S30')

    # dtype = {'names':names, 'formats':formats}
    # RELEASES = np.rec.array(blocks,dtype=dtype)
    return np.rec.fromrecords(blocks, names=names)


def read_trajectories(H, trajfile='trajectories.txt',
                      ncluster=5,
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
    RelTraj = pflexible.conv2netcdf4.Structure()
    Trajectories = []

    for i in range(3, 3 + (numpoint * 2), 2):
        i1, i2, xp1, yp1, xp2, yp2, zp1, zp2, k, npart, = \
            tuple([float(j) for j in alltraj[i].strip().split()])
        itimerel1 = dt + datetime.timedelta(seconds=i1)
        itimerel2 = dt + datetime.timedelta(seconds=i2)
        Xp = (xp1 + xp2) / 2
        Yp = (yp1 + yp2) / 2
        Zp = (zp1 + zp2) / 2
        RelTraj[
            alltraj[
                i +
                1].strip()] = np.array(
            (itimerel1,
             itimerel2,
             Xp,
             Yp,
             Zp,
             k,
             npart))

    for i in range(3 + (numpoint * 2), len(alltraj)):
        raw = alltraj[i]
        FMT = [0, 5, 8, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6] + \
              ncluster * [8, 8, 7, 6, 8]

        data = [raw[sum(FMT[:ii]):sum(FMT[:ii + 1])]
                for ii in range(1, len(FMT) - 1)] + [raw[sum(FMT[:-1]):]]
        ### FIX ###
        # To get rid of '******' that is now in trajectories.txt
        data = [float(r.replace('********', 'NaN')) for r in data]

        Trajectories.append(data)

    data = np.array(Trajectories)
    RelTraj['version'] = model + ' ' + version
    RelTraj['date'] = dt
    RelTraj['Trajectories'] = data
    RelTraj['labels'] = \
        ['release number', 'seconds prior to release', 'lon', 'lat', 'height',
         'mean topography', 'mean mixing height', 'mean tropopause height',
         'mean PV index', 'rms distance', 'rms', 'zrms distance', 'zrms',
         'fraction height??', 'fraction PV<2pvu',
         'fraction in troposphere'] + ncluster * [
            'xcluster', 'ycluster', 'zcluster', 'fcluster', 'rmscluster']
    RelTraj['info'] = """
    Returns a dictionary:
        R['Trajectories'] = array_of_floats(
            releasenum, it1, xi, yi, zi, topoi, hmixi, tropoi, pvi,
            rmsdisti, rmsi, zrmsdisti, zrmsi, hfri, pvfri, trfri,
            (xclusti(k), yclusti(k), zclusti(k), fclusti(k), rmsclusti(k),
             k=1,5))

        R['RELEASE_ID'] = (dt_i1, dt_i2, xp1, yp1, xp2, yp2, zp1, zp2, k, npart)
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


def groundlevel_for_line(H, X, Y, coords, index=0):
    """
    extracts ground level from H.heightnn along a track of lon, lat

    input:  H or H.heightnn (takes the lowest level)

            X, Y ,= H.longitude, H.latitude

            coords = zip(x, y)

    output: groundlevel (a 1-d array with shape (len(flighttrack)

    """
    try:
        hgt = H.Heightnn[:, :, 0]
    except:
        hgt = H

    # fix for hgt offset
    hgt = hgt - hgt.min()

    grndlvl = np.zeros(len(coords))

    for i, (x, y) in enumerate(coords):

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
        casl[:, i] = np.interp(H.asl_axis, xp + gl[i],
                               curtain_agl[:, i], left=below_gl)
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
    A = pflexible.conv2netcdf4.Structure()
    D = []
    for line in f[2:]:
        line = line.strip().split()
        # convert ibdate, ibtime to datetime
        dt = datetime.datetime.strptime(
            line[0] +
            line[1].zfill(6),
            '%Y%m%d%H%M%S')
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
        # ftype = 'AGECLASS'
        try:
            T = H.releasetimes
            spectrum = agespectra
            nClasses = spectrum.shape[1]
        except:
            # assume H is a list or None
            if H:
                nClasses = H[0]
                T = H[1]
                spectrum = agespectra
            else:
                T = agespectra[:, 0]
                nClasses = agespectra.shape[1] - 1
                spectrum = agespectra[:, 1:]
    elif spectype == 'contspec':
        # assume H is a list or None
        # ftype = 'SPECTRUM'
        try:
            T = H.releasetimes
            spectrum = agespectra
            nClasses = spectrum.shape[1]
            header = '## Continental Spectrum File'
        except:
            # assume H is a list or None
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
    cosfunc = lambda y: cos(y * pih) * r_earth

    nx = H['numxgrid']
    ny = H['numygrid']
    outlat0 = H['outlat0']
    dyout = H['dyout']
    dxout = H['dxout']
    area = np.zeros((nx, ny))

    for iy in range(ny):
        # NEED TO Check this, iy since arrays are 0-index
        ylata = outlat0 + (float(iy) + 0.5) * dyout
        ylatp = ylata + 0.5 * dyout
        ylatm = ylata - 0.5 * dyout
        if (ylatm < 0 and ylatp > 0):
            hzone = dyout * r_earth * pih
        else:
            # cosfact = cosfunc(ylata)
            cosfactp = cosfunc(ylatp)
            cosfactm = cosfunc(ylatm)
            if cosfactp < cosfactm:
                hzone = sqrt(r_earth ** 2 - cosfactp ** 2) - \
                    sqrt(r_earth ** 2 - cosfactm ** 2)
            else:
                hzone = sqrt(r_earth ** 2 - cosfactm ** 2) - \
                    sqrt(r_earth ** 2 - cosfactp ** 2)

        gridarea = 2. * pi * r_earth * hzone * dxout / 360.
        for ix in range(nx):
            area[ix, iy] = gridarea

    return area


def _get_header_version(bf):
    """
    Open and read the binary file (bf) header only to the point of
    the Flexpart version string, return the string, seek to the
    start of the file.

    """
    try:
        bf = pflexible.conv2netcdf4.BinaryFile(bf)
    except:
        bf = bf

    ret = bf.tell()
    bf.seek(12)  # start of version string
    version = bf.read('13S')
    bf.seek(ret)

    return version


def read_header(pathname, **kwargs):
    """
    The readheader function returns a special class (Structure) which behaves
    like a dictionary. It contains all the metadata from the simulation which
    is contained in the "header" or "header_nest" binary files from the model
    output.

    .. warning::
        It is recommended to use the :class:`Header` class: H = Header(path)

    This version is using the BinaryFile class rather than FortFlex.

    Usage::

        > H = read_header(pathname) #Don't include header filename

    Returns a dictionary

        H = dictionary like object with all the run metadata.
    Arguments

      .. tabularcolumns::  |l|L|

      =============       ========================================
      keyword             Description [default]
      =============       ========================================
      pathname            FLEXPART run output directory
      readp               read release points 0=no, [1]=y
      nested              nested output [False] or True
      headerfile          provide custom name for header file
      datefile            provide a custom name for the date file
      verbose             print information while loading header
      =============       ========================================

    .. note::
        **This function is in development**

        Please report any bugs found.

    .. TODO::

        probably a lot of things, among which ...

       [] choose skip/getbin or direct seek/read
       [] define output keys in docstring


    .. note::
        The user is no longer required to indicate which version of FLEXPART
        the header is from. Checks are made, and the header is read accordingly.
        Please report any problems...

    """
    h = pflexible.conv2netcdf4.Structure()

    h.options = OPS = pflexible.conv2netcdf4.Structure()
    OPS.readp = True
    OPS.nested = False
    OPS.ltopo = 1  # 1 for AGL, 0 for ASL
    OPS.verbose = False
    OPS.headerfile = None
    OPS.datefile = None

    # add keyword overides and options to header
    for k in kwargs.keys():
        if k not in OPS.keys():
            print("WARNING: {0} not a valid input option.".format(k))

    # BW compat fixes
    if 'nest' in kwargs.keys():
        raise IOError(
            "nest is no longer a valid keyword, see docs. \n "
            "Now use nested=True or nested=False")

    if 'nested' in kwargs.keys():
        # Force the use of true boolean values
        kwargs['nested'] = bool(kwargs['nested'])

    OPS.update(kwargs)

    if OPS.verbose:
        print "Reading Header with:\n"

        for o in OPS:
            print "%s ==> %s" % (o, OPS[o])

    # Define utility functions for reading binary file
    skip = lambda n = 8: bf.seek(n, 1)
    getbin = lambda dtype, n = 1: bf.read(dtype, (n,))

    if OPS.headerfile:
        filename = os.path.join(pathname, OPS.headerfile)

    elif OPS.nested is True:
        filename = os.path.join(pathname, 'header_nest')
        h['nested'] = True
    else:
        filename = os.path.join(pathname, 'header')
        h['nested'] = False

    # Open header file in binary format
    if not os.path.exists(filename):
        raise IOError("No such file: {0}".format(filename))
    else:
        try:
            bf = pflexible.conv2netcdf4.BinaryFile(filename, order="fortran")
        except:
            raise IOError(
                "Error opening: {0} with BinaryFile class".format(filename))

    # Get available_dates from dates file in same directory as header
    if OPS.datefile:
        datefile = os.path.join(pathname, OPS.datefile)
    else:
        datefile = os.path.join(pathname, 'dates')

    if not os.path.exists(datefile):
        raise IOError("No DATEFILE: {0}".format(datefile))
    else:
        try:
            fd = file(datefile, 'r').readlines()
        except:
            raise IOError("Could not read datefile: {0}".format(datefile))

    # get rid of any duplicate dates (a fix for the forecast system)
    fd = sorted(list(set(fd)))
    h['available_dates'] = [d.strip('\n') for d in fd]

    # which version format is header file:
    version = _get_header_version(bf)

    # required containers
    junk = []  # for catching unused output
    h['nz_list'] = []
    h['species'] = []
    h['wetdep'] = []
    h['drydep'] = []
    h['ireleasestart'] = []
    h['ireleaseend'] = []
    h['compoint'] = []

    # Define Header format and create Dictionary Keys
    I = {0: '_0', 1: 'ibdate', 2: 'ibtime', 3: 'flexpart',
         4: '_1', 5: 'loutstep', 6: 'loutaver', 7: 'loutsample',
         8: '_2', 9: 'outlon0', 10: 'outlat0', 11: 'numxgrid',
         12: 'numygrid', 13: 'dxout', 14: 'dyout', 15: '_3', 16: 'numzgrid',
         }
    # format for binary reading first part of the header file
    Dfmt = ['i', 'i', 'i', '13S', '2i', 'i', 'i', 'i', '2i', 'f', 'f', 'i',
            'i', 'f', 'f', '2i', 'i']
    if bf:
        a = [bf.read(fmt) for fmt in Dfmt]
        for j in range(len(a)):
            h[I[j]] = a[j]

        h['outheight'] = np.array([bf.read('f') for i in range(h['numzgrid'])])
        junk.append(bf.read('2i'))
        h['jjjjmmdd'] = bf.read('i')
        h['hhmmss'] = bf.read('i')
        junk.append(bf.read('2i'))
        h['nspec'] = bf.read('i') / 3  # why!?
        h['numpointspec'] = bf.read('i')
        junk.append(bf.read('2i'))

        # Read in the species names and levels for each nspec
        # temp dictionaries
        for i in range(h['nspec']):
            if 'V6' in version:
                h['wetdep'].append(
                    ''.join([bf.read('c') for i in range(10)]).strip())
                junk.append(bf.read('2i'))
                junk.append(bf.read('i'))
                h['drydep'].append(
                    ''.join([bf.read('c') for i in range(10)]).strip())
                junk.append(bf.read('2i'))
                h['nz_list'].append(getbin('i'))
                h['species'].append(
                    ''.join([getbin('c') for i in range(10)]).strip())

            else:
                junk.append(bf.read('i'))
                h['wetdep'].append(
                    ''.join([bf.read('c') for i in range(10)]).strip())
                junk.append(bf.read('2i'))
                junk.append(bf.read('i'))
                h['drydep'].append(
                    ''.join([bf.read('c') for i in range(10)]).strip())
                junk.append(bf.read('2i'))
                h['nz_list'].append(bf.read('i'))
                h['species'].append(
                    ''.join([bf.read('c') for i in range(10)]).strip())
                junk.append(bf.read('2i'))

        if 'V6' in version:
            bf.seek(8, 1)
        # pdb.set_trace()
        h['numpoint'] = bf.read('i')

        # read release info if requested
        # if OPS.readp pass has changed, we cycle through once,
        # then break the loop if OPS.readp is false. This is done
        # in order to get some of the information into the header.
        before_readp = bf.tell()

        # initialise release fields
        I = {2: 'kindz', 3: 'xp1', 4: 'yp1', 5: 'xp2',
             6: 'yp2', 7: 'zpoint1', 8: 'zpoint2', 9: 'npart', 10: 'mpart'}

        for k, v in I.iteritems():
            # create zero-filled lists in H dict
            h[v] = np.zeros(h['numpoint'])

        h['xmass'] = np.zeros((h['numpoint'], h['nspec']))

        if 'V6' in version:
            skip()
            _readV6(bf, h)
        else:
            junk.append(bf.read('i'))
            for i in range(h['numpoint']):
                junk.append(bf.read('i'))

                h['ireleasestart'].append(bf.read('i'))
                h['ireleaseend'].append(bf.read('i'))
                # This is an int16, might need to to change something
                h['kindz'][i] = bf.read('h')
                junk.append(bf.read('2i'))

                h['xp1'][i] = bf.read('f')
                h['yp1'][i] = bf.read('f')
                h['xp2'][i] = bf.read('f')
                h['yp2'][i] = bf.read('f')
                h['zpoint1'][i] = bf.read('f')
                h['zpoint2'][i] = bf.read('f')

                junk.append(bf.read('2i'))
                h['npart'][i] = bf.read('i')
                h['mpart'][i] = bf.read('i')

                junk.append(bf.read('i'))
                # initialise release fields

                l = bf.read('i')  # get compoint length?
                gt = bf.tell() + l  # create 'goto' point
                sp = ''
                # collect the characters for the compoint
                while re.search("\w", bf.read('c')):
                    bf.seek(-1, 1)
                    sp = sp + bf.read('c')

                bf.seek(gt)  # skip ahead to gt point

                # species names in dictionary for each nspec
                h['compoint'].append(sp)
                # h['compoint'].append(''.join([bf.read('c') for i in range(45)]))

                junk.append(bf.read('i'))
                # now loop for nspec to get xmass
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
                    if 'V6' in version:
                        bf.seek(157 * h['numpoint'], 1)
                    else:
                        bf.seek(119 * h['numpoint'] + (
                            h['nspec'] * 36) * h['numpoint'] + 4, 1)
                    break

        if 'V6' in version:
            h['method'] = bf.read('i')

        else:
            junk.append(bf.read('i'))
            junk.append(bf.read('i'))

        h['lsubgrid'] = bf.read('i')
        h['lconvection'] = bf.read('i')
        h['ind_source'] = bf.read('i')
        h['ind_receptor'] = bf.read('i')

        if 'V6' in version:
            pass
        else:
            junk.append(bf.read('2i'))

        h['nageclass'] = bf.read('i')

        Lage_fmt = ['i'] * h.nageclass
        h['lage'] = [bf.read(fmt) for fmt in Lage_fmt]

        # Orography
        nx = h['numxgrid']
        ny = h['numygrid']
        Dfmt = ['f'] * nx
        h['oro'] = np.zeros((nx, ny), np.float)
        junk.append(bf.read('2i'))

        for ix in range(nx):
            h['oro'][ix] = bf.read('f', ny)
            bf.seek(8, 1)

        # Why was this? / deprecated.
        # if h['loutstep'] < 0:
            # h['nspec'] = h['numpoint']

        bf.close()

    #############  ADDITIONS ###########
    # attributes to header that can be added after reading below here
    h['pathname'] = pathname
    h['decayconstant'] = 0
    h.path = pathname

    # Calculate Height (outheight + topography)
    # There is an offset issue here related to the 0-indexing. Be careful.
    oro = h['oro']  # z is a numpy array
    nx = h['numxgrid']
    ny = h['numygrid']
    nz = h['numzgrid']

    Heightnn = np.zeros((nx, ny, nz), np.float)
    for ix in range(nx):
        if OPS.ltopo == 1:
            Heightnn[ix, :, 0] = oro[ix, :]
        else:
            Heightnn[ix, :, 0] = np.zeros(ny)

        for iz in range(nz):
            if OPS.ltopo == 1:
                Heightnn[ix, :, iz] = [h['outheight'][iz] + oro[ix, y]
                                       for y in range(ny)]
            else:
                Heightnn[ix, :, iz] = h['outheight'][iz]

    h['Area'] = gridarea(h)
    h['Heightnn'] = Heightnn
    h['nx'] = nx
    h['ny'] = ny
    h['nz'] = nz

    # Convert ireleasestart and ireleaseend to datetimes, fix issue #10
    start_day = datetime.datetime.strptime(str(h['ibdate']), '%Y%m%d')
    H, M, S, = [int(str(h['ibtime']).zfill(6)[i:i + 2])
                for i in range(0, 6, 2)]
    start_time = datetime.timedelta(hours=H, minutes=M, seconds=S)
    h['simulationstart'] = start_day + start_time

    if OPS.readp:
        releasestart, releaseend = [], []
        for i in range(h.numpointspec):
            releasestart.append(h.simulationstart +
                                datetime.timedelta(seconds=int(h.ireleasestart[i])))
            releaseend.append(h.simulationstart +
                              datetime.timedelta(seconds=int(h.ireleaseend[i])))
        h.releasestart = releasestart
        h.releaseend = releaseend[:h.numpointspec]
        h.releasetimes = [b - ((b - a) / 2)
                          for a, b in zip(h.releasestart, h.releaseend)]

    # Add datetime objects for dates
    available_dates_dt = []
    for i in h.available_dates:
        available_dates_dt.append(datetime.datetime(
            int(i[:4]), int(i[4:6]), int(i[6:8]), int(i[8:10]), int(i[10:12]), int(i[12:])))
    h.available_dates_dt = available_dates_dt
    h.first_date = available_dates_dt[0]
    h.last_date = available_dates_dt[-1]
    h.ageclasses = np.array(
        [act - h.simulationstart for act in h.available_dates_dt])
    h.numageclasses = len(h.ageclasses)

    # Add other helpful attributes
    h.nxmax = h.numxgrid
    h.nymax = h.numygrid
    h.nzmax = h.numzgrid
    h.maxspec = h.nspec
    h.maxpoint = h.numpoint
    h.maxageclass = h.numageclasses

    h.area = h.Area  # fix an annoyance

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
        h.unit = 'conc'  # could be pptv
        h.plot_unit = 'ppb'
    else:
        h.direction = 'backward'
        h.unit = 'time'
        h.plot_unit = 'ns / kg'  # Not sure about this

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

    h.fp_version = h.flexpart
    h.junk = junk  # the extra bits that were read... more for debugging

    print('Read {0} Header: {1}'.format(h.fp_version, filename))

    return h

# BW Compatability
readheader = read_header
readheaderV8 = read_header
readheaderV6 = read_header


def _readV6(bf, h):
    """
    Version 6, 7X read routines...

    """
    # Utility functions
    skip = lambda n = 8: bf.seek(n, 1)
    getbin = lambda dtype, n = 1: bf.read(dtype, (n,))

    if bf:
        # bf.read('i')
        print h['numpoint']
        for i in range(h['numpoint']):
            # r2=getbin('i')
            i1 = getbin('i')
            i2 = getbin('i')
            # h['ireleasestart'].append( ss_start + datetime.timedelta(seconds=float(i1)) )
            # h['ireleaseend'].append( ss_start + datetime.timedelta(seconds=float(i2)) )
            h['ireleasestart'].append(i1)
            h['ireleaseend'].append(i2)
            # This is an int16, might need to to change something
            h['kindz'][i] = getbin('i')
            skip()  # get xp, yp,...

            h['xp1'][i] = getbin('f')
            h['yp1'][i] = getbin('f')
            h['xp2'][i] = getbin('f')
            h['yp2'][i] = getbin('f')
            h['zpoint1'][i] = getbin('f')
            h['zpoint2'][i] = getbin('f')

            skip()  # get n/mpart
            h['npart'][i] = getbin('i')
            h['mpart'][i] = getbin('i')

            getbin('i')
            # initialise release fields

            l = getbin('i')  # get compoint length?
            gt = bf.tell() + l  # create 'goto' point
            sp = ''
            sp = sp.join(getbin('c', l))
            bf.seek(gt)  # skip ahead to gt point

            # species names in dictionary for each nspec
            h['compoint'].append(sp)
            # h['compoint'].append(''.join([getbin('c') for i in range(45)]))
            skip()
            # r1=getbin('i')

            # now loop for nspec to get xmass
            for v in range(h['nspec']):
                Dfmt = ['i', 'f', '2i', 'f', '2i', 'f', 'i']
                a = [bf.read(fmt) for fmt in Dfmt]
                h['xmass'][i, v] = a[1]

    return
