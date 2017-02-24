#!/usr/bin/env python

"""
This is an example providing a very brief overview of some of the
pflexible functionality.

Test Data:
    Note, FLEXPART output is large. I have made available a complete
    run for one day from the POLARCAT NOAA-ICEALOT Cruise here:
    http://niflheim.nilu.no/~burkhart/sharing/pflexpart_testdata.tgz
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import addcyclic

import reflexible as rf
import reflexible.mapping as mp
from reflexible.plotting import plot_sensitivity, plot_totalcolumn, plot_at_level


def plot_backward(SOURCE_FILE, OUTPUT_DIR):
    #read the header of a FLEXPART output
    H = rf.Header(SOURCE_FILE, nested=True)

    # # Get the trajectories, these are overlayed when we plot the
    # # sensitivities below.
    # XXX reading trajectories does not work yet
    # T = pf.read_trajectories(H)

    # Plot everything up
    TC = None #TC Empty figure object, the idea is to reuse this
    FP = None #FP Empty figure object

    # iterate over every species and every timestep (k)

    # The region depends on the file, this should be a parameter.
    map_region = 'north_sea'
    datainfo_str = 'NORTHSEA'

    for s,k in H.C:
        data = H.C[(s,k)]
        # total column
        TC = plot_totalcolumn(H, data, datainfo_str=datainfo_str, map_region=map_region)
        # TC = plot_trajectory(H, T, k, FIGURE=TC)
        filename = '%s_tc_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        TC.fig.savefig(ofilename)

        # footprint
        FP = plot_at_level(H, data, datainfo_str=datainfo_str, map_region=map_region)
        # FP = plot_trajectory(H, T, k, FIGURE=FP)
        filename = '%s_fp_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        FP.fig.savefig(ofilename)


if __name__ == "__main__":
    """
    Run through plotting routines.

    Expected two params:
    - The netCDF4 file
    - The output dir

    This example is meant to be run with the data in the file
    http://folk.uio.no/johnbur/sharing/stads2_V10.tar
    Use the nested file as, at the time of this writing, the non-nested file is
    broken.
    """

    src_file, out_dir = sys.argv[1:]
    plot_backward(src_file, out_dir)
