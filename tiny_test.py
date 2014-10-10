import pflexible as pf

H = pf.Header('../test_data/')
print H.keys()

print "Calling fill_backward():"
H.fill_backward(nspec=(0,1))

# H.C should be available by now
print "H.C.keys():"
print H.C[(0,1)].keys()  # this does not work yet

# H.FD should be available by now
print "H.FD.keys():"
print H.FD[(0, '20100527210000')].keys()  # this does not work yet

print "Calling read_grid():"
FD = pf.read_grid(H, time_ret=0,nspec_ret=0)

T = pf.read_trajectories(H)
print "T.keys():"
print T.keys()
