#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Post processing of Multiple FLEXPART Simulations.

"""

import os
import inspect
import datetime as dt
import pandas as pd
import xarray as xr

_relative_home = os.path.join(os.environ['HOME'], 'hycamp/')


def merge_gains_gfed_runs(years=(2003, 2005)):
    retrievals = {'gains': gains_files,
                  'gfed': gfed_files
                  }

    t_axis = pd.date_range(start=dt.datetime(years[0], 1, 1),
                           end=dt.datetime(years[1], 1, 1),
                           freq='D')

    files = []
    for retrieval in retrievals:
        files.extend(retrievals[retrieval](t_axis))

    return files


def gfed_files(time_axis):
    """ return gfed files """
    filepath = os.path.join(_relative_home,
                            "team/jfb/FLEXPART/CASES/Felix/BCRUNS/OUTPUT/ECMWF/GFED_{0}/GFED_{0}.001/grid_conc_{0}01000000.nc")
    dates = pd.date_range(start=time_axis[0], end=time_axis[-1], freq='M')

    return [filepath.format((date.strftime('%Y%m'))) for date in dates]


def gains_files(time_axis):
    """ return gains filenames """

    filepath = os.path.join(_relative_home,
                            "team/jfb/FLEXPART/CASES/Felix/AnthroBC/OUTPUT/ECMWF/GAINS_{0}01/GAINS_{0}01.001/grid_conc_{0}0101000000.nc")
    dates = list(set([t.year for t in time_axis]))[:-1]
    return [filepath.format(str(year)) for year in dates]


def concat_simulations(fp_files, variables=['spec001_mr', 'WD_spec001', 'DD_spec001']):
    """ concatentate a series of flexpart simulations

    works with a list of netcdf input files from FP > v10 and a list of variable
    names. Only those variables will be returned in the concatenated dataset.

    Inputs:
        fp_files = list of flexpart netcdf files
        variables = list of variable names to concat

    Assumptions:
        Uses a function `reduce_dataset` which extracts the variables, it also
        takes a finite difference for all the deposition variables (cumulative)
        Requires all runs have the same:
            loutstep, dxout, dyout and dimesions

    """
    # The following could work, but is too memory intensive
    # Need to do some profiling, to see where speedup/opt is possible
    # with one year of runs (1 yearly GAINS, 12 monthly GFED) this method
    # is 3 minutes, vs. 4 minutes as implemented, but memory cons. higher
    if False:
        dsets = [xr.open_dataset(f) for f in fp_files]
        nds = [reduce_dataset(d, variables) for d in dsets]
        nds = concat_runs(nds)

        clsd = [d.close() for d in dsets]
        ds = nds

    for i, f in enumerate(fp_files):
        print(f)

        with xr.open_dataset(f) as dset:
            dset = xr.open_dataset(f)
            nds = reduce_dataset(dset, variables)

        if i == 0:
            ds = nds
            continue

        ds = concat_runs([ds, nds])
        dset.close()

        _end_datetime = {'iedate': dset.attrs['iedate'],
                         'ietime': dset.attrs['ietime']}

    # provide some further attribute information
    ds.attrs.update(_end_datetime)
    ds.attrs['history'] = """ Concatentation of:
     {0}
     using: {1}""".format(fp_files, inspect.stack()[0][3])
    ds.attrs['institution'] = "UiO"
    ds.attrs['contact'] = "John F. Burkhart <john.burkhart@geo.uio.no>"

    return ds.squeeze()  # we return and get rid of nageclass/pointspec dims


def reduce_dataset(dset, variables):
    """ extract variables and attrs and return new dataset

    assumes deposition variables should be differenced to provide daily
    values (required in order to sum multiple run simulations)
    """

    _cum_vars = ['WD_spec001', 'DD_spec001']

    nds = {}
    for v in variables:

        if v in _cum_vars:

            # dset[v].loc[str(dset[v].time[1].values) : str(dset[v].time[-1].values)] = dset[v].diff(dim='time')
            dset[v].loc[dict(time=slice(str(dset[v].time[1].values), str(dset[v].time[-1].values)))] = dset[v].diff(
                dim='time')

            nds[v] = dset[v]
        else:
            nds[v] = dset[v]

    nds = xr.Dataset(nds)
    nds.attrs = dset.attrs

    return nds


def concat_runs(dsets):
    """ simple concatentation of runs

    sums across 'new' dimension to deal with overlapping dates
    """

    cc = xr.concat(dsets, dim='new')

    return cc.sum(dim='new', keep_attrs=True)


def main():
    period = [2003, 2004]

    fp_files = merge_gains_gfed_runs(years=period)
    dset = concat_simulations(fp_files)
    dset.to_netcdf('fp_global-bc_{0}-{1}.nc'.format(*period))
    return dset


if __name__ == '__main__':
    dset = main()
