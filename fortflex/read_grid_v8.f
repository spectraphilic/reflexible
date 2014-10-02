      subroutine readgrid(filegrid,numpointspec,numxgrid,
     +numygrid,numzgrid,
     +nageclass,wetgrid,drygrid,grid,scaledepo,scaleconc,
     +decayconstant)

      parameter(smallnum=1.e-38)

      real decayconstant
      integer lage(0:nageclass),itime,numpointspec
      real grid(numxgrid,numygrid,numzgrid,1,nageclass)
Cf2py intent(out) grid
      real wetgrid(numxgrid,numygrid,1,nageclass)
Cf2py intent(out) wetgrid
      real drygrid(numxgrid,numygrid,1,nageclass)
Cf2py intent(out) drygrid

***************** new sparse output
      integer fact,ii,ir,npspec,np
      real sparse_dump_r(numxgrid*numygrid*numzgrid)
      integer sparse_dump_i(numxgrid*numygrid*numzgrid)
      integer sp_count_i,sp_count_r
***************** new sparse output

      character*500 filegrid

c allocate memory for sparse matrix
c      allocate(sparse_dump_r(numxgrid*numygrid*numzgrid),STAT=stat)
c      if (stat.ne.0) print*,'***error allocating sparse_dump_r'
c      allocate(sparse_dump_i(numxgrid*numygrid*numzgrid),STAT=stat)
c      if (stat.ne.0) print*,'***error allocating sparse_dump_i'
c
c     print*,numxgrid,numygrid,numzgrid,nageclass,numpointspec

      open(10,file=filegrid,form='unformatted',status='old')
      read(10,end=99) itime 
c     write(*,*) itime

C Initialize age class fields
*****************************

C      wetgrid=0.
C      drygrid=0.
C      grid=0.

C Loop about all age classes
****************************

      do 38 nage=1,nageclass

      do 38 np=1,numpointspec
      npspec=1

c     print*,npspec

c       do 63 ix=1,numxgrid
c         do 63 jy=1,numygrid
c           wetgrid(ix,jy,npspec,nage)=0.
c           drygrid(ix,jy,npspec,nage)=0.
c           do 63 kz=1,numzgrid
c63           grid(ix,jy,kz,npspec,nage)=0.

C Read wet deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 32 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             pos=pos+1
          endif

          jy=n/numxgrid
          ix=n-numxgrid*jy
          wetgrid(ix+1,jy+1,npspec,nage)=wetgrid(ix+1,jy+1,npspec,nage)+
     >      sparse_dump_r(ir)

32      continue

       print*,'wet grid read'

C Read dry deposition
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 36 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             pos=pos+1
          endif

          jy=n/numxgrid
          ix=n-numxgrid*jy
          drygrid(ix+1,jy+1,npspec,nage)=drygrid(ix+1,jy+1,npspec,nage)+
     >      sparse_dump_r(ir)

36      continue

       print*,'dry grid read'

C Read concentrations
*********************

        fact=1
        read(10) sp_count_i
        read(10) (sparse_dump_i(ix),ix=1,sp_count_i)
        read(10) sp_count_r
        read(10) (sparse_dump_r(ix),ix=1,sp_count_r)
c       write(*,*) sp_count_i,sp_count_r
c       write(*,*) (sparse_dump_i(ix),ix=1,sp_count_i)
c       write(*,*) (sparse_dump_r(ix),ix=1,sp_count_r)

        ii=0
        do 47 ir=1,sp_count_r
          if ((sparse_dump_r(ir)*fact).gt.smallnum) then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1.)
          else
             n=n+1
          endif

          kz=n/(numxgrid*numygrid)
          jy=(n-kz*numxgrid*numygrid)/numxgrid
          ix=n-numxgrid*numygrid*kz-numxgrid*jy
          grid(ix+1,jy+1,kz,npspec,nage)=grid(ix+1,jy+1,kz,npspec,nage)+
     >      abs(sparse_dump_r(ir))

47      continue

       print*,'conc grid read'
       print*,kz,n,numxgrid*numygrid
       print*,'scaling factors: fac decay scaleconc scaledep'

C Scale values as defined by user
*********************************

          fac=exp(-1.*float(lage(nage)+lage(nage-1))/2./decayconstant)
          if (lage(nage).gt.1e8) fac=1.0
c         print*,fac
          do 41 ix=1,numxgrid
            do 41 jy=1,numygrid
              drygrid(ix,jy,npspec,nage)=
     >          drygrid(ix,jy,npspec,nage)*scaledepo
              wetgrid(ix,jy,npspec,nage)=
     >          wetgrid(ix,jy,npspec,nage)*scaledepo
              do 41 kz=1,numzgrid
               grid(ix,jy,kz,npspec,nage)=
     >           grid(ix,jy,kz,npspec,nage)*fac*scaleconc
41             continue
          print*,fac,decayconstant,scaleconc,scaledepo

C End species loop, end age class loop
**************************************

138       continue
38      continue

      close(10)

99    continue

c      if (allocated(sparse_dump_r)) deallocate(sparse_dump_r)
c      if (allocated(sparse_dump_i)) deallocate(sparse_dump_i)

      return
      end
