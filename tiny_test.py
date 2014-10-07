import pflexible as pf

H = pf.Header('../test_data/')
print H.keys()

print "Calling fill_backward():"
H.fill_backward()

# H.C should be available by now
print "H.C.keys():"
print H.C.keys()  # this does not work yet

# H.FD should be available by now
print "H.FD.keys():"
print H.FD.keys()  # this does not work yet

print "Calling read_grid():"
FD = pf.read_grid(H, time_ret=0,nspec_ret=0)
