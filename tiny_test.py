import pflexible as pf

H = pf.Header('../test_data/')
print H.keys()

#print "Calling fill_backward():"
#H.fill_backward()

print "Calling read_grid():"
FD = pf.read_grid(H, time_ret=0,nspec_ret=0)
