import sys
import os.path

import reflexible as rf
import reflexible.conv2netcdf4 as conv

H = conv.Header("/tmp/test_data")
print "******* Original Header ***********"
print "keys:", H.keys()
print "available_dates:", H.available_dates
print "releasestart:", H.releasestart
print "releaseend:", H.releaseend
print "releasetimes", H.releasetimes
print ""

Hnc = rf.Header("/tmp/grid_time_20100601210000.nc")
print "******* netCDF4 Header ***********"
print "keys:", Hnc.keys()
print "available_dates:", Hnc.available_dates
print "releasestart:", Hnc.releasestart
print "releaseend:", Hnc.releaseend
print "releasetimes:", Hnc.releasetimes
print "species:", Hnc.species
print "FD:", Hnc.FD
print "FD[(0,'20100518210000')].species", Hnc.FD[(0,'20100518210000')].species
#print "FD[(0, '20100518210000')].grid:", Hnc.FD[(0, '20100518210000')].grid
# print "FD[(0, '20100518210000')].itime:", Hnc.FD[(0, '20100518210000')].itime
# print "FD[(0, '20100518210000')].timestamp:", Hnc.FD[(0, '20100518210000')].timestamp
# print "FD[(0, '20100518210000')].rel_i:", Hnc.FD[(0, '20100518210000')].rel_i
C = Hnc.C
print "C:", C
print "C[(0,1)].grid[1:3]:", C[(0,1)].grid[1:3]
print "C[(0,1)].slabs[1]:", C[(0,1)].slabs[1]
print ""

# print "Calling fill_backward():"
# H.fill_backward(nspec=(0,))

# # H.C should be available by now
# print "H.C.keys():"
# print H.C[(0,2)].keys()

# # H.FD should be available by now
# print "H.FD.keys():"
# print H.FD[(0, '20070121220000')].keys()

# print "Calling read_grid():"
# FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

# T = conv.read_trajectories(H)
# print "T.keys():"
# print T.keys()

# # Backward Runs
# #
# H = conv.Header(rf.datasets['Bwd1_V9.02'])
# print H.keys()

# print "Calling fill_backward():"
# H.fill_backward(nspec=(0,))

# # H.C should be available by now
# print "H.C.keys():"
# print H.C[(0,0)].keys()

# # H.FD should be available by now
# print "H.FD.keys():"
# print H.FD[(0, H.available_dates[0])].keys()

# print "Calling read_grid():"
# FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

