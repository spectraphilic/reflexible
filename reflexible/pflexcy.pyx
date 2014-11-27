import numpy as np
cimport numpy as np

def dumpdatagrid( np.ndarray[int, ndim=1] dmp_i, int cnt_r, np.ndarray[float, ndim=1] dmp_r, np.ndarray[np.float64_t, ndim=5] grd, int k, int nage,int nx,int ny):
    """ function to dump sparse elementVs into grid """

    cdef int ir,ii,pos,n,kz,jy,ix
    cdef float fact
    ii=0
    fact=1
    pos = 0
    for ir in range(cnt_r):

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


    return grd


def dumpdepogrid( np.ndarray[int, ndim=1] dmp_i, int cnt_r, np.ndarray[float, ndim=1] dmp_r, np.ndarray[np.float64_t, ndim=4] grd, int k, int nage,int nx,int ny):
    """ function to dump sparse elementVs into grid """

    cdef int ir,ii,pos,n,kz,jy,ix
    cdef float fact
    ii=0
    fact=1
    pos = 0
    for ir in range(cnt_r):

        if dmp_r[ir]*fact>0:
            n=dmp_i[ii]
            ii=ii+1
            fact=fact*-1.
        else:
            n=n+1
        jy=n/nx
        ix=n-nx*jy
        grd[ix,jy,k,nage]=abs(dmp_r[ir])

    return grd #flipud(grd.transpose())

