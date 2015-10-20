import sys
sys.path.append('/home/johnbur/jfb/git/reflexible')
import reflexible as ref
import datetime as dt

from dateutil import relativedelta

def gen_releases():
    # Set up the period of simulations
    start_month = dt.datetime(2005, 1, 1, 0, 0, 0)
    end_month = dt.datetime(2015, 1, 1, 0, 0, 0)
    delta_run = relativedelta.relativedelta(months=1)
    
    
    # Set up the COMMAND options
    command_defaults = {'sim_start' : dt.datetime(2008, 1, 1, 0, 0, 0),
            'sim_end' : dt.datetime(2008, 2, 1, 0, 0, 0),
            'loutstep' : 216000,
            'loutaver' : 216000,
            'loutsample' : 3600,
            'itsplit' : 9999999,
            'lsynctime' : 900,
            'ctl' : -5,
            'ifine' : 4,
            'iout' : 1,
            'ipout' : 0,
            'lsubgrid' :             1,
            'lconvection' : 1,
            'lagespectra' : 1,
            'ipin' : 0,
            'ioutputforeachrelease' : 0,
            'iflux' : 0,
            'mdomainfill' : 0,
            'ind_source' : 1,
            'ind_receptor' : 1,
            'mquasilag' : 0,
            'nested_output' : 0,
            'linit_cond' : 1,
            'surf_only' : 0,
            'cblflag' : 0,
                }
    
    
    release_defaults = {
            'idate1' : ['20010101', '''YYYYMMDD begin date of release '''],
            'itime1' : ['000000', '''YYYYMMDD begin time of release '''],
            'idate2' : ['20010201', '''YYYYMMDD end date of release '''],
            'itime2' : ['000000', '''YYYYMMDD end time of release '''],
            'lon1'  : [ 120, '''lowerleft Longitude'''],
            'lon2'  : [ 130, '''upperright Longitude'''],
            'lat1'  : [ 55, '''lowerleft Latitude'''],
            'lat2'  : [ 60, '''upperright Latitude'''],
            'z1'    : [ 20, '''lower boundary of release point (m)'''],
            'z2'    : [ 100, '''upper z-level of release point (m)'''],
            'zkind' : [ 3, ''' 1 for m above ground, 2 for m above sea level, 3 for pressure in hPa'''],
            'mass'  : [ [1.0], '''mass of species'''],
            'nspec' : [ 1, '''number of species'''],
            'parts' : [50000, '''total number of particles in release'''],
            'specnum_rel': [(22,), '''tuple of species number id'''],
            'run_ident': ['comment', '''character*40 comment''']
            }


    












    run_month = start_month
    while run_month < end_month:
        overrides = {'sim_start' : run_month,
                    'sim_end' : run_month + delta_run,
                    'ageclasses' : [86400*20]}
        command.update( overrides )
    
        C = ref.Command(**command)
        cfile = run_month.strftime('RUN_TEST_%Y%m.command')
        C.write_command(cfile)
    
        run_month += delta_run
        
    
    
    