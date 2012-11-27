'''
Created on Jul 10, 2012

@author: jfb

partoutput module for working with FLEXPART partoutput data
'''

import os

import numpy as np
import datetime
from collections import OrderedDict

import pflexible as pf


def read_particlepositions(H, **kwargs):
    """
    Reads the partposit files when using FLEXPART with IOFR set to 1
    
    The particle positions are dumped to a unformatted binary fortran
    output file. The routine reads the binary grids and creates a dictionary
    keyed by the datestring of the output file.
    
    Each dictionary has two attributes:
    
    
    partposit
    ##########
    
    An array of values, each row corresponding to a particle and columns
    as described in the table below.
        
    ----------    ---------    ----------------------------------
    variable      column        description
    ----------    ---------    ----------------------------------
    xlonin        0            longitude
    ylatin        1            latitude
    ztraj         2            elevation (includes topography)
    itramem       3            memorized particle release time
    topo          4            topography
    pvi           5            potential vorticity
    qvi           6            specific humidity
    rhoi          7            density, converted to pressure (TODO: check units)
    hmixi         8            PBL height
    tri           9           Tropopause height 
    tti           10           temperature [K]
    ----------    ---------    -----------------------------------
    
    
    xmass
    ########
    
    An array of size (numparticles, numspecies) giving the mass of the species
    along the trajectory.
    
    (xmass(i,j),j=1,nspec)
    
    npoint
    ########
    
    An array of size (numpoarticles) giving the npoint (int identifier) of each
    particle (if enabled in FLEXPART, otherwise '1')
    
    See the source code for conversion conventions and units.
    
    """
    from FortFlex import readparticles
    
    ## OPS is the options Structure, sets defaults, then update w/ kwargs
    OPS = pf.Structure()
    OPS.npspec_int = False  # allows to select an index of npsec when calling readgrid
    OPS.time_ret = 0
    OPS.date = None
    OPS.BinaryFile = False
    ## add keyword overides and options to header
    OPS.update(kwargs)
    
    ## set up the return dictionary (FLEXDATA updates fd, fd is returned)
    PP = OrderedDict()
    
    ## get times to return
    get_dates = None
    if OPS.time_ret is not None:
        get_dates = []
        time_ret = OPS.time_ret
        if isinstance(time_ret, int) == True:
            time_ret = [time_ret]
        
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
                raise IOError("Cannot find date: %s in H['available_dates']\n" % d)

    if get_dates is None:
        raise ValueError("Must provide either time_ret or date value.")
    else:
        ## assign grid dates for indexing fd
        pass
    #PP.grid_dates = get_dates[:]

    # Some predifinitions
    prefix = ['partposit_', 'partposit_nest_']
    join = os.path.join
    
    for date_i in range(len(get_dates)):
        datestring = get_dates[date_i]
        print datestring
        dt_key = datetime.datetime.strptime(datestring, '%Y%m%d%H%M%S')
            
        PP[dt_key] = pf.Structure()
        
        if H.nested is True:
            filename = join(H['pathname'], \
                        prefix[1] + datestring)
            

        else:
            filename = join(H['pathname'], \
                        prefix[0] + datestring)
            
        if os.path.exists(filename):

            if OPS.BinaryFile:
                raise IOError("BinaryFile not yet implemented")
                #pp, xm = _readpartBF(H, filename)
            else:
                pp, xm, ip = readparticles(filename, H.npart[0], H.nspec, H.dxout, H.dyout, H.outlon0, H.outlat0)
            """ pp = (xlonin, ylatin, ztraj, itramem, topo, pvi, qvi, rhoi, hmixi, tri, tti) """
            """ xm = (xmass(i,j),j=1,nspec) """
            """ np = (npoint(i)) """
            
        else:
            print("Cannot find partposit file: {0}".format(filename))    
        

        
        ## Some unit conversions
        #pp[:,7] = (pp[:,7] /1200000.0+0.025) * 1000
        #tconv = pp[:,10]/180.0 + 100.0 
        #pp[:,7] = pp[:,7] / 10000.0 * 287.0 * tconv / 100.0 ## convert density to pressure, assuming dry air
        pp[:,7] = pp[:,7]*287.058*pp[:,10] # method of vidit.
        
        
        ## Drop the -9999.99 values
        #itramem = pp[:,4]
        ppc = np.ma.masked_where(np.isnan(pp), pp, copy=False)
        pp = np.ma.masked_where(pp <= -9999., ppc, copy=False)
        #ppc[:,4] = itramem
        
        xmc = np.ma.masked_where(np.isnan(xm), xm, copy=False)
        xm = np.ma.masked_where(xm <= -9999., xmc, copy=False) 
        
        PP[dt_key].partposit = pp
        PP[dt_key].xmass = xm
        PP[dt_key].npoint = ip
        
            
    return PP

def dictarrays_as_array(PARTS, array_key='partposit'):
    """ creates one single array of all arrays in a dictionary.
    Assumes a dictionary of dictionaries, with the array keyed by
    "array_key".
    
    Written to work with defaults for dictionary returned from
    :func:`read_particlepositions` """

    return np.dstack((PARTS[k][array_key] for k in sorted(PARTS.keys())))


def eminusp(H, allparticles, qindx=6, per_hr=True):
    """ returns E-P for all particles and for all output timesteps.
    Needs the Header to determine which direction to take the difference.
    """
    
    if H.direction == 'forward':
        emp = np.diff(allparticles[:,qindx,:])
    elif H.direction == 'backward':
        emp = np.diff(allparticles[:,qindx,-1::-1]) #calculate starting from the end
    
    if per_hr:
        emp = (emp / np.abs(H.loutstep) ) / 3600.
        
    return emp
    
    