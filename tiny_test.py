import reflexible as pf
import reflexible.conv2netcdf4 as conv

H = conv.Header(pf.datasets['Fwd1_V9.02'])
print H.keys()

print "Calling fill_backward():"
H.fill_backward(nspec=(0,))

# H.C should be available by now
print "H.C.keys():"
print H.C[(0,2)].keys()

# H.FD should be available by now
print "H.FD.keys():"
print H.FD[(0, '20070121220000')].keys()

print "Calling read_grid():"
FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

T = conv.read_trajectories(H)
print "T.keys():"
print T.keys()

# Backward Runs
#
H = conv.Header(pf.datasets['Bwd1_V9.02'])
print H.keys()

print "Calling fill_backward():"
H.fill_backward(nspec=(0,))

# H.C should be available by now
print "H.C.keys():"
print H.C[(0,0)].keys()

# H.FD should be available by now
print "H.FD.keys():"
print H.FD[(0, H.available_dates[0])].keys()

print "Calling read_grid():"
FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

