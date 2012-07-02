import numpy as np
cimport numpy as np
def dumpgrid( np.ndarray[int, ndim=1] dmp_i, int cnt_r, np.ndarray[float, ndim=1] dmp_r, np.ndarray[np.float64_t, ndim=5] grd, int k, int nage,int nx,int ny):
    """ function to dump sparse elementVs into grid """

    cdef int ir,ii,pos,n,kz,jy,ix
#    cdef double fact
#    cdef np.ndarray[np.float64_t, ndim=5]
    conc = True
    #if len(grd.shape) == 5:
    #    conc = True
    ii=0
    fact=1
    pos = 0
    for ir in range(cnt_r):

        if conc:
            if dmp_r[ir]*fact>0:
                n=dmp_i[ii]
                ii=ii+1
                fact=fact*-1.
            else:
                n=n+1

            kz=n/(nx*ny)
            jy=(n-kz*nx*ny)/nx
            ix=n-nx*ny*kz-nx*jy
            grd[ix,jy,kz-1,k,nage]=abs(dmp_r[ir])

#
#                print "n  ==> ix,jy,kz,k,nage"
#                print "%s ==> %s,%s,%s,%s,%s" % (n,ix,jy,kz,k,nage)
#                print grd.shape
#                print grd[0,0,0,0,0]

        else:
            if dmp_r[ir]*fact>0:
                n=dmp_i[ii]
                ii=ii+1
                fact=fact*-1.
            else:
                pos=pos+1
            jy=n/nx
            ix=n-nx*jy
            grd[ix,jy,k,nage]=abs(dmp_r[ir])

    return grd #flipud(grd.transpose())

