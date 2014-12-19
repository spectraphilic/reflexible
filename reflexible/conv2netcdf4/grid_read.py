########### Grid Reading Routines ###############

import itertools
import datetime
import os

import numpy as np

import reflexible.conv2netcdf4
from .FortFlex import sumgrid
from .helpers import _shout



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
    %   - nested: nested = True
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
    # if os.sys.platform != 'win32':
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
    # else: date = H['ibtime']
    else:
        date = None

    if 'unit' in kwargs.keys():
        unitname = kwargs['unit']
    else:
        unitname = 'time'

    units = ['conc', 'pptv', 'time', 'footprint']
    unit = units.index(unitname)

    if 'nspec_ret' in kwargs.keys():
        nspec_ret = kwargs['nspec_ret']
    else:
        nspec_ret = 1

    if 'pspec_ret' in kwargs.keys():
        pspec_ret = kwargs['ppsec_ret']
    else:
        pspec_ret = 1

    if 'age_ret' in kwargs.keys():
        age_ret = kwargs['age_ret']
    else:
        age_ret = 1

    if 'nested' in kwargs.keys():
        nested = kwargs['nested']
    else:
        nested = False

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
        scaleconc = 1.0

    if 'decaycons' in kwargs.keys():
        decaycons = kwargs['decaycons']
    else:
        decaycons = 99999999990
    fail = -1

    # set filenames
    prefix = ['grid_conc_', 'grid_pptv_',
              'grid_time_', 'footprint_total',
              'grid_conc_nest_', 'grid_pptv_nest_',
              'grid_time_nest_', 'footprint_total_nest']

    # local functions for reading binary data and creating grid dump
#    skip = lambda n=8 : f.read(n)
#    getbin = lambda fmt,n=1 : struct.unpack(fmt*n,f.read(struct.calcsize(fmt)))[0]

    # Utility functions
    skip = lambda n = 8: f2.seek(n, 1)
    getbin = lambda dtype, n = 1: f2.read(dtype, (n,))

    def getdump(n, fmt='f'):
        """ function to get the dump values for the sparse format """
        # print n
        skip()
        # Dfmt=[fmt]*n
#        a=[struct.unpack(ft,f.read(struct.calcsize(ft))) for ft in Dfmt]
        a = f2.read(fmt, n)
#        dumplist=[a[j][0] for j in range(len(a))]
        # dumplist=[a[j] for j in range(len(a))]
        return a  # dumplist

    def key2var(D, key):
        cmd = "global %s; %s = D['%s'];" % (key, key, key)
        exec(cmd)

    def dumpgrid(dmp_i, cnt_r, dmp_r, grd, k, nage):
        """ function to dump sparse elements into grid """
        ii = 0
        fact = 1
        for ir in range(cnt_r):
            if dmp_r[ir] * fact > 0:
                n = dmp_i[ii]
                ii = ii + 1
                fact = fact * -1.
            else:
                n = n + 1  # XXX n is actually set before this?

            kz = n / (numxgrid * numygrid)
            jy = (n - kz * numxgrid * numygrid) / numxgrid
            ix = n - numxgrid * numygrid * kz - numxgrid * jy
            # print "n  ==> ix,jy,kz,k,nage"
            # print "%s ==> %s,%s,%s,%s,%s" % (n,ix,jy,kz,k,nage)
            # print grd.shape
            # print grd[0,0,0,0,0]
            grd[ix, jy, kz - 1, k, nage] = abs(dmp_r[ir])
        return grd  # flipud(grd.transpose())

    G = {}
    # get values from header file
    H['nageclass'] = H['numageclasses']
    headervars = ['nested', 'nspec', 'numxgrid', 'numygrid', 'numzgrid',
                  'nageclass', 'available_dates', 'pathname', 'decayconstant',
                  'numpoint', 'Area', 'Heightnn']

    for k in headervars:
        key2var(H, k)

    # create zero arrays for datagrid and zplot
    datagrid = np.zeros((numxgrid, numygrid, numzgrid, numpoint),
                        dtype=np.float32)
    zplot = np.empty((numxgrid, numygrid, numzgrid, numpoint))

    #--------------------------------------------------
    # Loop over all times, given in field H['dates']
    #--------------------------------------------------
    for ks in range(nspec_ret, nspec_ret + 1):
            # print 'processing: ' + str(H['dates'][date_i]) + '   ' + str(ks)
            # print 'processing: ' + str(date) + '   ' + str(ks)
        specs = '_' + str(ks).zfill(3)
        fpmax = -999.
        if date is None and time_ret is None:
            get_dates = [available_dates[0]]
        elif time_ret is None:
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
                # datestring = date
                filename = os.path.join(H['pathname'],
                                        prefix[(unit) + (nested * 4)] + datestring + specs)
            else:
                filename = os.path.join(H['pathname'],
                                        prefix[(unit) + (nested * 4)])

            # print 'reading: ' + filename

            if os.path.exists(filename):
                # print nspec,numxgrid,numygrid,numzgrid,nageclass,scaledepo,scaleconc,decaycons
                # print ks
                # print date

                if useFortFlex == 1:
                    # print 'Using FortFLEX'
                    ###########################################################
                    # FORTRAN WRAPPER CODE - only works on linux
                    # #
                    # USAGE: grid = readgrid(filegrid,numxgrid,numygrid,numzgrid,\
                    # nspec,nageclass,scaleconc,decayconstant)
                    # #
                    # zplot = sumgrid(zplot, grid, area, heightnn,
                    #         [numxgrid,numygrid,numzgrid,numpoint,nageclass])
                    # #
                    # #
                    # NOTE: numpoint = number of releases, ageclass always 1 for backward
                    # #
                    # RETURN: grid(numxgrid,numygrid,numzgrid,numpoint,nageclass)
                    ###########################################################

                    concgrid = readgrid(filename, numxgrid, numygrid, numzgrid,
                                        numpoint, nageclass,
                                        scaleconc, decayconstant)

                    # contribution[:,:,:] = concgrid[:,:,:,:,0]
                    print np.min(concgrid)
                    print np.max(concgrid)

                    # altitude = 50000
                    zplot = sumgrid(zplot, concgrid,
                                    H.area, H.Heightnn)

                else:
                    dat_cnt = 0
                    nage = 0
                    # read data:
                    # datagrid=np.zeros((numxgrid,numygrid,numzgrid[nspec-1],
                    #                    nspec,nageclass),np.float)
                    datagrid = np.zeros(
                        (numxgrid, numygrid, numzgrid, 1, 1), np.float)
                    # f = file(filename, 'rb')
                    # print filename
                    f2 = reflexible.conv2netcdf4.BinaryFile(filename, order='fortran')
                    skip(4)
                    G['itime'] = getbin('i')
                    print H['available_dates'][date_i]

                    # Read Wet Depostion
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
                    # wet=dumpgrid(dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage)
                    # Read Dry Deposition
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
                    # dry=dumpgrid(dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage)

                    # Read Concentrations
                    skip()
                    cnt_i = getbin('i')
                    dmp_i = getdump(cnt_i, 'i')
                    skip()
                    cnt_r = getbin('i')
                    dmp_r = getdump(cnt_r)
                    # print dmp_i, cnt_r, dmp_r, datagrid, ks-1, nage
                    concgrid = dumpgrid(dmp_i, cnt_r, dmp_r, datagrid,
                                        ks - 1, nage)
                    # G[H['ibtime']].append(concgrid) G[H['ibtime']] = concgrid
                    f2.close()
                fail = 0
            else:
                print "\n\n INPUT ERROR: Could not find file: %s" % filename
                raise IOError('No file: %s' % filename)

        if useFortFlex == 1:
            G = zplot
    return G

readgridBF = _readgrid_noFF


def _readgridBF(H, filename):
    """ Read grid using BinaryFile class"""
    # Utility functions
    skip = lambda n = 8: f2.seek(n, 1)
    getbin = lambda dtype, n = 1: f2.read(dtype, (n,))

    def getdump(n, fmt='f'):
        """ function to get the dump values for the sparse format """
        skip()
        # Dfmt=[fmt]*n
#        a=[struct.unpack(ft,f.read(struct.calcsize(ft))) for ft in Dfmt]
        a = f2.read(fmt, n)
#        dumplist=[a[j][0] for j in range(len(a))]
        # dumplist=[a[j] for j in range(len(a))]
        return a  # dumplist

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
                    n = n + 1
                jy = n / H.numxgrid
                ix = n - H.numxgrid * jy
                grd[ix, jy, k, nage] = abs(dmp_r[ir])

        return grd  # flipud(grd.transpose())

    # Import pflexcy.so (cython compiled version of dumpgrid)
    try:
        # Using Pyximport for compiling on-the flight.
        # See: http://docs.cython.org/src/userguide/source_files_and_compilation.html#pyximport
        import pyximport; pyximport.install()
        from pflexcy import dumpdatagrid, dumpdepogrid
        print 'using pflexcy'
    except:
        print """WARNING: Using PURE Python to readgrid, execution will be slow.
         Try compiling the FortFlex module or the pflexcy module
         for your machine. For more information see the
         reflexible/f2py_build directory or use cython with pflexcy.pyx
        """
        dumpdatagrid = _dumpgrid
        dumpdepogrid = _dumpgrid

    dat_cnt = 0
    nage = 1
    # datagrid=np.zeros((numxgrid,numygrid,numzgrid[nspec-1],nspec,nageclass),np.float)
    wetgrid = np.zeros((H.numxgrid, H.numygrid, H.numpointspec, 1), np.float)
    drygrid = np.zeros((H.numxgrid, H.numygrid, H.numpointspec, 1), np.float)
    datagrid = np.zeros(
        (H.numxgrid, H.numygrid, H.numzgrid, H.numpointspec, nage), np.float)
    # f = file(filename,'rb')
    # print filename
    f2 = reflexible.conv2netcdf4.BinaryFile(filename, order='fortran')
    # read data:
    skip(4)
    itime = getbin('i')

    for na in range(nage):

        for ks in range(H.numpointspec):
            # Read Wet Depostion
            skip()
            cnt_i = getbin('i')
            dmp_i = getdump(cnt_i, 'i')
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            if dmp_r.any():
                # print dmp_r, dmp_i
                wetgrid = dumpdepogrid(
                    dmp_i, cnt_r, dmp_r, wetgrid, ks, na, H.numxgrid,
                    H.numygrid)

            # Read Dry Deposition
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            if dmp_r.any():
                # print dmp_r, dmp_i
                drygrid = dumpdepogrid(
                    dmp_i, cnt_r, dmp_r, drygrid, ks, na, H.numxgrid,
                    H.numygrid)

            # Read Concentrations
            skip()
            cnt_i = getbin('i')
            dmp_i = getdump(cnt_i, 'i')
            skip()
            cnt_r = getbin('i')
            dmp_r = getdump(cnt_r)
            # print len(dmp_r),len(dmp_i)
            # print cnt_r,cnt_i
            # print dmp_i
            # print type(dmp_i),type(cnt_r),type(dmp_r),type(datagrid),type(ks),type(na)
            # print type(H.numxgrid),type(H.numygrid)
            datagrid = dumpdatagrid(
                dmp_i, cnt_r, dmp_r, datagrid, ks, na, H.numxgrid, H.numygrid)

    # G[H['ibtime']].append(concgrid)
    f2.close()

    return datagrid, wetgrid, drygrid, itime


def _read_headerFF(pathname, h=None,
                   maxpoint=800000, maxspec=4, maxageclass=20,
                   nxmax=722, nymax=362, nzmax=40, verbose=True):
    """
    DEPRECATED

    Called from read_header if readp_ff is True, uses FortFlex.readheader
    to read all releasepoints.

    This function is dependant on the FortFlex.so module
    see FortFlex.f and the f2py directory
    """
    try:
        from .FortFlex import readheader
    except:
        print "Error with FortFlex.readheader, use read_header"

    headervars = ['numxgrid', 'numygrid', 'numzgrid', 'outlon0', 'outlat0',
                  'compoint', 'dxout', 'dyout', 'outheight', 'ibdate',
                  'ibtime', 'loutstep', 'nspec', 'nageclass', 'lage',
                  'ireleasestart', 'ireleaseend', 'numpoint',
                  'xpoint', 'ypoint', 'zpoint1', 'zpoint2', 'heightnn', 'area',
                  'maxpoint', 'maxspec', 'maxageclass', 'nxmax', 'nymax',
                  'nzmax', 'npart', 'kind', 'lage', 'loutaver', 'loutsample',
                  'yyyymmdd', 'hhmmss', 'method']
    if h is None:
        h = reflexible.conv2netcdf4.Structure()

    if verbose:
        print """Reading Header with:
                    maxpoint : %s
                    maxspec : %s
                    maxageclass : %s
                    nxmax : %s
                    nymax : %s
                    nzmax : %s
                    """ % (maxpoint, maxspec, maxageclass, nxmax, nymax, nzmax)

    (numxgrid, numygrid, numzgrid, outlon0, outlat0, dxout, dyout, outheight,
     ibdate, ibtime, loutstep, nspec, nageclass, lage, ireleasestart,
     ireleaseend, numpoint, xpoint, ypoint, zpoint1, zpoint2, heightnn, area,
     compoint, species_f, npart, kind, loutaver, loutsample, yyyymmdd, hhmmss,
     method) = readheader(pathname, maxpoint, maxspec, maxageclass,
                          nxmax, nymax, nzmax)

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
        an integer. See grid.keys().

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
      nested                obtained from H['nested']
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
        most arguments are able to be extracted from the header "H"

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
    # set up the return dictionary (FLEXDATA updates fd, fd is returned)
    FLEXDATA = {}
    fd = reflexible.conv2netcdf4.Structure()

    # OPS is the options Structure, sets defaults, then update w/ kwargs
    fd.options = OPS = reflexible.conv2netcdf4.Structure()
    OPS.unit = H.unit
    OPS.getwet = False
    OPS.getdry = False
    OPS.nspec_ret = 0
    # allows to select an index of npsec when calling readgrid
    OPS.npspec_int = False
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
    # add keyword overrides and options to header
    OPS.update(kwargs)
    # H.update(OPS)

    # What direction is the run?
    unit = OPS.unit

    if H['loutstep'] > 0:
        forward = True
        if unit == 'time':
            # default forward unit
            unit = 'conc'
            OPS.unit = unit
    else:
        forward = False

    # What species to return?
    nspec_ret = OPS.nspec_ret
    if isinstance(nspec_ret, int):
        nspec_ret = [nspec_ret]
    assert iter(nspec_ret), "nspec_ret must be iterable."

    # get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, int) == True:
            time_ret = [time_ret]

        if time_ret[0] < 0:
            # if forward == False:
                # get all dates for calculating footprint.
            time_ret = np.arange(len(H.available_dates))
            # else:
            #    raise ValueError("Must enter a positive time_ret for forward runs")

        for t in time_ret:
            get_dates.append(H.available_dates[t])

    # define what dates to extract if user has explicitly defined a 'date'
    if OPS.date is not None:
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
        # assign grid dates for indexing fd
        fd.grid_dates = get_dates[:]

    print 'getting grid for: ', get_dates
    # Some pre-definitions
    fail = 0
    # set filename prefix
    prefix = ['grid_conc_', 'grid_pptv_',
              'grid_time_', 'footprint_', 'footprint_total',
              'grid_conc_nest_', 'grid_pptv_nest_',
              'grid_time_nest_', 'footprint_nest_', 'footprint_total_nest'
              ]

    units = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
    unit_i = units.index(unit)

    # Determine what module to read, try to use FortFlex, then dumpgrid, lastly pure Python
    # import the FortFlex / Fortran module
    try:
        print('Assumed V8 Flexpart')
        from .FortFlex import readgrid, sumgrid
        useFortFlex = True
    except:
        # get the original module (no memory allocation)
        try:
            from nilu.pflexpart.FortFlex import readgrid, sumgrid
            useFortFlex = True
            # print 'using nilu.pflexpart FortFlex'
        except:
            useFortFlex = False
            print('Cannot load FortFlex, reverting to BinaryFile.')
    if not useFortFlex:
        readgrid = _readgridBF
        OPS.BinaryFile = True

    # reserve output fields
    print H.numxgrid, H.numygrid, H.numzgrid, OPS.nspec_ret, OPS.pspec_ret, OPS.age_ret, len(get_dates), H.numpoint

    # -------------------------------------------------

    # add the requests to the fd object to be returned
    OPS.unit = unit

    #--------------------------------------------------
    # Loop over all times, given in field H['available_dates']
    #--------------------------------------------------

    for date_i in range(len(get_dates)):
        datestring = get_dates[date_i]
        print datestring
        for s in nspec_ret:  # range(OPS.nspec_ret,OPS.nspec_ret+1):A
            FLEXDATA[(s, datestring)] = fdc = reflexible.conv2netcdf4.FDC()
            spec_fid = '_' + str(s + 1).zfill(3)

            if unit_i != 4:
                filename = os.path.join(
                    H['pathname'],
                    prefix[(unit_i) + (H.nested * 5)] + datestring + spec_fid)
                H.zdims = H.numzgrid

            else:
                # grid total footprint
                print "Total footprint"
                filename = os.path.join(
                    H['pathname'],
                    prefix[(unit_i) + (H.nested * 5)] + spec_fid)
                H.zdims = 1

            if os.path.exists(filename):
                H.filename = filename
                # print 'reading: ' + filename
                if OPS.verbose:
                    print 'with values:'
                    inputvars = ['filename', 'numxgrid', 'numygrid',
                                 'zdims', 'numpoint', 'nageclass',
                                 'scaledepo', 'scaleconc',
                                 'decayconstant', 'numpointspec']
                    for v in inputvars:
                        print v, " ==> ", H[v]

                if OPS.BinaryFile:
                    print("Reading {0} with BinaryFile".format(filename))
                    gridT, wetgrid, drygrid, itime = _readgridBF(H, filename)
                else:
                    # Quick fix for Sabine's Ship releases, added nspec_int so that only one
                    # field of the nspec dimension is actually read
                    if OPS.npspec_int is not False:
                        npspec_int = OPS.npspec_int
                        numpointspec = 1
                    else:
                        npspec_int = 0
                        numpointspec = H.numpointspec

                    gridT, wetgrid, drygrid, itime = readgrid(
                        filename, H.numxgrid, H.numygrid, H.zdims,
                        numpointspec, H.nageclass, OPS.scaledepo,
                        OPS.scaleconc, H.decayconstant, npspec_int)

                if OPS.getwet:
                    wet = wetgrid
                if OPS.getdry:
                    dry = drygrid
                if forward:
                    zplot = gridT[:, :, :, :, 0]
                else:
                    zplot = gridT[:, :, :, :, 0]

                if OPS.calcfoot:

                    zplot = sumgrid(zplot, gridT, H.area, H.Heightnn)

                # get the total column and prep the grid
                if H.direction == 'forward':
                    # not trying to do anything here... must be done
                    # after retrieving the grid
                    # D = get_slabs(H,np.squeeze(zplot))
                    rel_i = 0  # H.available_dates.index(datestring)
                    D = zplot

                else:
                    D = zplot
                    rel_i = 'k'

                # NOTE:
                # If you're changing things here, you might want to change
                # them in fill_backward as well, yes I know... something is
                # poorly designed ;(

                fdc.grid = D  # zplot

                fdc.itime = itime

                fdc.timestamp = \
                    datetime.datetime.strptime(datestring, '%Y%m%d%H%M%S')
                fdc.species = H['species'][s]
                fdc.gridfile = filename
                fdc.rel_i = rel_i
                fdc.spec_i = s
                if OPS.getwet:
                    fdc.wet = wet
                else:
                    fdc.wet = None

                if OPS.getdry:
                    fdc.dry = dry
                else:
                    fdc.dry = None

            else:
                _shout('***ERROR: file %s not found! \n' % filename)
                fail = 1

        fd.set_with_dict(FLEXDATA)  # XXX keys are tuples, is this really intended?
        try:
            # just for testing, set the first available grid as a shortcut
            # this will be removed.
            # TODO: this can be removed now?
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
    # set up the return dictionary (FLEXDATA updates fd, fd is returned)
    FLEXDATA = {}
    fd = reflexible.conv2netcdf4.Structure()

    # OPS is the options Structure, sets defaults, then update w/ kwargs
    fd.options = OPS = reflexible.conv2netcdf4.Structure()
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
    # add keyword overrides and options to header
    OPS.update(kwargs)

    # What direction is the run?
    unit = OPS.unit
    if H['loutstep'] > 0:
        forward = True
        if unit == 'time':
            # default forward unit
            unit = 'conc'
    else:
        forward = False

    # What species to return?
    nspec_ret = OPS.nspec_ret
    if isinstance(nspec_ret, int):
        nspec_ret = [nspec_ret]
    assert iter(nspec_ret), "nspec_ret must be iterable."

    # get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, int) == True:
            time_ret = [time_ret]

        if time_ret[0] < 0:
            if forward == False:
                # get all dates for calculating footprint.
                time_ret = np.arange(len(H.available_dates))
            else:
                raise ValueError(
                    "Must enter a positive time_ret for forward runs")

        for t in time_ret:
            get_dates.append(H.available_dates[t])

    # define what dates to extract if user has explicitly defined a 'date'
    if OPS.date is not None:
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
        # assign grid dates for indexing fd
        fd.grid_dates = get_dates[:]

    print 'getting grid for: ', get_dates
    # Some pre-definitions
    fail = 0
    # set filename prefix
    prefix = ['grid_conc_', 'grid_pptv_',
              'grid_time_', 'footprint_', 'footprint_total',
              'grid_conc_nest_', 'grid_pptv_nest_',
              'grid_time_nest_', 'footprint_nest_', 'footprint_total_nest'
              ]

    units = ['conc', 'pptv', 'time', 'footprint', 'footprint_total']
    unit_i = units.index(unit)
    # Determine what module to read, try to use FortFlex, then dumpgrid, lastly pure Python
    # import the FortFlex / Fortran module
    try:
        from .FortFlex import readgrid_v6 as readgrid
        from .FortFlex import sumgrid
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

    # cheat: use key2var function to get values from header dict, H
    # need to change this to using H.xxxxx
    # headervars = ['nested','nspec','numxgrid','numygrid','numzgrid','nageclass',\
    #              'dates','pathname','decayconstant','numpoint','numpointspec',
    #              'area','Heightnn','lage']
    # for k in headervars:
    #    key2var(H,k)

   # reserve output fields
    print H.numxgrid, H.numygrid, H.numzgrid, OPS.nspec_ret, OPS.pspec_ret, OPS.age_ret, len(get_dates), H.numpoint

    # -------------------------------------------------

    # if forward:
    #    numpointspec = 1
    # else:
    #    numpointspec = numpoint
    # numpointspec = H['numpointspec']
    # if unit_i == 4:
    # grid=np.empty((numxgrid,numygrid,1,nspec_ret,pspec_ret,age_ret,len(get_dates)))
    # else:
    # grid=np.empty((numxgrid,numygrid,numzgrid,nspec_ret,pspec_ret,age_ret,len(get_dates)))
    # zplot=np.empty((numxgrid,numygrid,numzgrid,numpointspec))

    #--------------------------------------------------
    # Loop over all times, given in field H['available_dates']
    #--------------------------------------------------

    for date_i in range(len(get_dates)):
        datestring = get_dates[date_i]
        print datestring
        for s in nspec_ret:  # range(OPS.nspec_ret,OPS.nspec_ret+1):
            FLEXDATA[(s, datestring)] = fdc = reflexible.conv2netcdf4.FDC()
            # spec_fid = '_'+str(s+1).zfill(3)

            if unit_i != 4:
                filename = os.path.join(H['pathname'],
                                        prefix[(unit_i) + (H.nested * 5)] + datestring)
                H.zdims = H.numzgrid

            else:
                # grid total footprint
                print "Total footprint"
                filename = os.path.join(H['pathname'],
                                        prefix[(unit_i) + (H.nested * 5)])
                H.zdims = 1

            if os.path.exists(filename):
                H.filename = filename
                # print 'reading: ' + filename
                if OPS.verbose:
                    print 'with values:'
                    inputvars = ['filename', 'numxgrid', 'numygrid',
                                 'zdims', 'numpoint', 'nageclass',
                                 'scaledepo', 'scaleconc',
                                 'decayconstant', 'numpointspec']
                    for v in inputvars:
                        print v, " ==> ", H[v]

                if OPS.BinaryFile:
                    # print 'Using BinaryFile'
                    gridT, wetgrid, drygrid, itime = _readgridBF(H, filename)
                else:
                    gridT, wetgrid, drygrid, itime = readgrid(
                        filename, H.numxgrid, H.numygrid, H.zdims,
                        H.numpointspec, H.nageclass, OPS.scaledepo,
                        OPS.scaleconc, H.decayconstant)

                if forward:
                    zplot = gridT[:, :, :, :, 0]
                else:
                    zplot = gridT[:, :, :, :, 0]

                if OPS.calcfoot:
                    zplot = sumgrid(zplot, gridT,
                                    H.area, H.Heightnn)

                # get the total column and prep the grid
                if H.direction == 'forward':
                    # not trying to do anything here... must be done
                    # after retrieving the grid
                    # D = get_slabs(H,np.squeeze(zplot))
                    # rel_i = H.available_dates.index(datestring)
                    D = zplot
                    rel_i = 'k'
                else:
                    D = zplot
                    rel_i = 'k'

                # NOTE:
                # If you're changing things here, you might want to change
                # them in fill_backward as well, yes I know... something is
                # poorly designed ;(
                fdc.grid = D  # zplot
                fdc.itime = itime
                fdc.timestamp = datetime.datetime.strptime(
                    datestring, '%Y%m%d%H%M%S')
                fdc.species = H['species'][s]
                fdc.gridfile = filename
                fdc.rel_i = rel_i
                fdc.spec_i = s

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


def monthly_footprints(H):

    footprints = np.zeros((H.ny, H.nx, len(H.C.keys())))
    for i, key in enumerate(H.C):
        footprints[:, :, i] = H.C[key].slabs[0]

    #footprints = np.average(footprints, axis=2)

    return footprints


def fill_grids(H, nspec=0, FD=None):
    """ for backward runs, calculates the 20-day sensitivity at each release point.

    Usage::

        > FDC = fill_backward(H, nspec=(0))


    This will cycle through all available_dates and create the filled backward array
    for each k in H.numpointspec.

    Returns

        A dictionary keyed by a (species,k) tuple.

    Each element in the dictionary is a 3D array (x,y,z) for each species,k

    .. note::
        USE with *caution*, this is **MEMORY** intensive!

    Arguments


    ==============        ========================================
    keyword               Description [default]
    ==============        ========================================
    nspec                 the specied ID or a tuple of species IDs
    FD                    FD can be passed if it is already read
    ==============        ========================================

    .. todo::
        There's a lot of redundancy in the storage of attributes, maybe there is a
        better way to handle this.

    """

#    assert H.direction == 'backward', "fill_backward is only valid for backward runs"
    # make sure npsec is iterable
    if isinstance(nspec, int):
        species = [nspec]
    else:
        species = nspec
    assert iter(species), 'nspec must be iterable, or you can pass an int'
    # initialize variables
    if FD is None:
        # then we need to read the grids
        FD = read_grid(H, time_ret=-1, nspec_ret=species)

    C = reflexible.conv2netcdf4.Structure()

    if H.direction == 'backward':

        for s, k in itertools.product(species, range(H.numpointspec)):
            C[(s, k)] = c = reflexible.conv2netcdf4.FDC()
            c.grid = np.zeros((H.numxgrid, H.numygrid, H.numzgrid))
            c.itime = None
            c.timestamp = H.releasetimes[k]
            c.species = H['species'][s]
            c.gridfile = 'multiple'
            c.rel_i = k
            c.spec_i = s

        # read data grids and attribute/sum sensitivity
        print species
        for s in species:

            for d in FD.grid_dates:
                # cycle through all the date grids (20days back)
                for k in range(H.numpointspec):
                    # cycle through each release point
                    contribution = FD[(s, d)].grid[:, :, :, k]
                    C[(s, k)].grid += contribution

    # added this, not sure if it makes sense to have this for 'forward'
    # however, adding 'C' does not add to memory, as it points to same
    # memory location as FD
    if H.direction == 'forward':
        for s in species:
            for k in range(len(FD.grid_dates)):
                d = FD.grid_dates[k]
                C[(s, k)] = FD[(s, d)]

    for s, k in C:
        # add total column
        C[(s, k)].slabs = get_slabs(H, C[(s, k)].grid)

    H.C = C
    H.FD = FD
    return C


def read_emissions(emissionsfile, E=None, maxemissions=1):
    """ Use reademissions.so module to read emissions file """
    try:
        from .FortFlex import reademissions
    except:
        raise ImportError(
            "Cannot find FortFlex module or missing reademissions")
    if not E:
        E = reflexible.conv2netcdf4.Structure()
        # set defaults for global 0.5 degree emissions
        defaults = {'nxmax': 720, 'nymax': 360, 'outlon0': -180,
                    'outlat0': -90, 'numxgrid': 720, 'numygrid': 360,
                    'dxout': 0.5, 'dyout': 0.5}
        for k, v in defaults.iteritems():
            exec("E.%s = %s" % (k, v))

    emissions = reademissions(emissionsfile, maxemissions, E.nxmax, E.nymax,
                              E.outlon0, E.outlat0, E.numxgrid, E.numygrid,
                              E.dxout, E.dyout)
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
    Slabs = reflexible.conv2netcdf4.Structure()
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
                raise IOError(
                    '######### ERROR: Which Release Point to get? ########## ')

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
