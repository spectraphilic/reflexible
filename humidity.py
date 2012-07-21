'''
Humidity module, converted from watersip
Created on Jul 22, 2012

@author: zvj

'''
import numpy as np

#returns the saturation water
#vapour pressure over water for
#absolute temperature T
#after Flatau et al. 1992 (range +/-50\C)

def esw(T):
    Ttr=273.15
    a1=6.11176750
    a2=0.443986062
    a3=0.143053301e-01
    a4=0.265027242e-03
    a5=0.302246994e-05
    a6=0.203886313e-07
    a7=0.638780966e-10
    esw = (a1 + a2*(T-Ttr) + a3*((T-Ttr)**2) + a4*((T-Ttr)**3) + a5*((T-Ttr)**4) + a6*((T-Ttr)**5) + a7*((T-Ttr)**6 ))*100
    return esw


#Reurns the saturation water
#vapour pressure over ice for
#absolute temperature T
#after Flatau et al. 1992 (range +/-50C)

def esi(T):
    Ttr=273.15
    a1=6.10952665
    a2=0.501948366
    a3=0.186288989e-01
    a4=0.403488906e-03
    a5=0.539797852e-05
    a6=0.420713632e-07
    a7=0.147271071e-09
    esi = (a1 + a2*(T-Ttr) + a3*((T-Ttr)**2) + a4*((T-Ttr)**3) + a5*((T-Ttr)**4) + a6*((T-Ttr)**5) + a7*((T-Ttr)**6 ))*100
    return esi


#calculates the relative humidity
#with respect to water
#from specific humidity given a
#certain tempertature and pressure

def relative_humidity_water(qv,T,p):
   
    wv = qv / (1-qv)             # convert to mixing ratio
    e  = p*wv / (wv+0.622)       # calculate vapour pressure
    rhw = e / esw(T) * 100.00    #calculate rh from saturation over water
    return rhw


#calculates the relative humidity
#with respect to ice
#from specific humidity given a
#certain tempertature and pressure

def relative_humidity_ice(qv,T,p):

    wv = qv / (1-qv)          #convert to mixing ratio
    e  = p*wv / (wv+0.622)    #calculate vapour pressure
    rhi = e/esi(T) * 100.0    #calculate rh from saturation over water
    return rhi