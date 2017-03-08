#!/usr/bin/env python

"""
This is an example providing a very brief overview of some of the
pflexible functionality.

Test Data:
    Note, FLEXPART output is large. I have made available a complete
    run for one day from the POLARCAT NOAA-ICEALOT Cruise here:
    http://niflheim.nilu.no/~burkhart/sharing/pflexpart_testdata.tgz
"""


from argparse import ArgumentParser
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import addcyclic

import reflexible as rf
import reflexible.mapping as mp
from reflexible.plotting import plot_sensitivity, plot_totalcolumn, plot_at_level


def plot_backward(SOURCE_FILE, OUTPUT_DIR, region, nested=False, label=''):
    #read the header of a FLEXPART output
    H = rf.Header(SOURCE_FILE, nested=nested)

    # # Get the trajectories, these are overlayed when we plot the
    # # sensitivities below.
    # XXX reading trajectories does not work yet
    # T = pf.read_trajectories(H)

    # Plot everything up
    TC = None #TC Empty figure object, the idea is to reuse this
    FP = None #FP Empty figure object

    # iterate over every species and every timestep (k)
    for s,k in H.C:
        data = H.C[(s,k)]
        # total column
        TC = plot_totalcolumn(H, data, map_region=region, datainfo_str=label)
        # TC = plot_trajectory(H, T, k, FIGURE=TC)
        filename = '%s_tc_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        TC.fig.savefig(ofilename)

        # footprint
        FP = plot_at_level(H, data, map_region=region, datainfo_str=label)
        # FP = plot_trajectory(H, T, k, FIGURE=FP)
        filename = '%s_fp_%s.png' % (data.species,
                                     data.timestamp.strftime('%Y%m%dT%H:%M:%S'))
        ofilename = os.path.join(OUTPUT_DIR, filename)
        FP.fig.savefig(ofilename)


if __name__ == "__main__":
    """
    Run through plotting routines.

    Expected three params:
    - The netCDF4 file
    - The output dir
    - The map region

    Optional args:
    - nested (default False)
    - label (default empty)

    This example is meant to be run with the data in the file
    http://folk.uio.no/johnbur/sharing/stads2_V10.tar
    Use the nested file as, at the time of this writing, the non-nested file is
    broken.
    """

    parser = ArgumentParser(description='Plot the given netCDF4 file')
    parser.add_argument('src', help='netCDF4 file')
    parser.add_argument('out', help='output directory')
    parser.add_argument('region', default='north_sea',
                        help='region name (e.g. north_sea)')
    parser.add_argument('--label',
                        help='label to draw over the plot (e.g. NORTHSEA)')
    parser.add_argument('--nested', dest='nested', action='store_true',
                        help='whether the source file is nested or not')

    args = parser.parse_args()
    plot_backward(args.src, args.out, args.region, nested=args.nested,
                  label=args.label)
