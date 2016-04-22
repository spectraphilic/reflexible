from __future__ import print_function
import reflexible as rf
import reflexible.conv2netcdf4 as conv

H = conv.Header(rf.datasets['Fwd1_V9.02'])
print(list(H.keys()))

print("Calling fill_backward():")
H.fill_backward(nspec=(0,))

# H.C should be available by now
print("H.C.keys():")
print(list(H.C[(0,2)].keys()))

# H.FD should be available by now
print("H.FD.keys():")
print(list(H.FD[(0, '20070121220000')].keys()))

print("Calling read_grid():")
FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

T = conv.read_trajectories(H)
print("T.keys():")
print(list(T.keys()))

# Backward Runs
#
H = conv.Header(rf.datasets['Bwd1_V9.02'])
print(list(H.keys()))

print("Calling fill_backward():")
H.fill_backward(nspec=(0,))

# H.C should be available by now
print("H.C.keys():")
print(list(H.C[(0,0)].keys()))

# H.FD should be available by now
print("H.FD.keys():")
print(list(H.FD[(0, H.available_dates[0])].keys()))

print("Calling read_grid():")
FD = conv.read_grid(H, time_ret=0,nspec_ret=0)

